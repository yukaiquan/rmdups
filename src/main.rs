use anyhow::{Context, Result};
use bstr::BStr;
use clap::Parser;
use noodles::bam;
use noodles::bgzf::io::Writer as BgzfWriter;
use rayon::prelude::*;
use roaring::RoaringBitmap;
use std::collections::{BinaryHeap, HashMap, HashSet};
use std::fs::File;
use std::io::{BufReader, Write};
use std::sync::Arc;
use std::time::Instant;
use tempfile::Builder;

#[cfg(not(windows))]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

mod args;
mod metadata;
mod algorithm;
mod io;
mod utils;

use args::{Args, effective_threads};
use metadata::{Metadata, MergeItem};
use algorithm::{get_5p_pos, get_score, identify_dups};
use io::{write_header, record_to_bytes, toggle_duplicate_flag, open_chunk_reader};
use utils::format_duration;

fn main() -> Result<()> {
    let args = Args::parse();

    // Determine effective thread count
    let threads = effective_threads(&args);

    // Set rayon thread pool size (only affects parallel operations)
    if threads > 1 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .ok();
    }

    let total_start = Instant::now();
    let tmp_dir = Builder::new()
        .prefix("markdup_rust")
        .tempdir_in(args.tmp_dir.clone().unwrap_or_else(std::env::temp_dir))?;

    eprintln!("rmduprs: using {} threads{}", threads, if args.single_threaded { " (single-threaded mode)" } else { "" });

    let mut reader = bam::io::reader::Builder::default().build_from_path(&args.input)?;
    let header = Arc::new(reader.read_header()?);

    // Build library map
    let mut lib_map = HashMap::new();
    for (_id, rg) in header.read_groups() {
        let lib_name = rg
            .other_fields()
            .get(noodles::sam::alignment::record::data::field::Tag::LIBRARY.as_ref())
            .map(|v| v.to_string())
            .unwrap_or_else(|| "unknown".to_string());
        let next_id = lib_map.len() as i32;
        lib_map.entry(lib_name).or_insert(next_id);
    }

    let header_clone = header.clone();
    let get_lib_id = move |rec: &bam::Record| -> i32 {
        rec.data()
            .get(noodles::sam::alignment::record::data::field::Tag::READ_GROUP.as_ref())
            .and_then(|v| v.ok())
            .and_then(|v| {
                if let noodles::sam::alignment::record::data::field::Value::String(s) = v {
                    header_clone
                        .read_groups()
                        .get::<BStr>(s.as_ref())
                        .and_then(|rg| {
                            let lib_name = rg
                                .other_fields()
                                .get(noodles::sam::alignment::record::data::field::Tag::LIBRARY.as_ref())
                                .map(|v| v.to_string())
                                .unwrap_or_else(|| "unknown".to_string());
                            lib_map.get(&lib_name).cloned()
                        })
                } else {
                    None
                }
            })
            .unwrap_or(0)
    };

    let find_start = Instant::now();
    let mut pe_count: u64 = 0;
    let mut se_count: u64 = 0;
    let mut unmatched_pairs_count: u64 = 0;

    eprintln!("finding positions of the duplicate reads in the file...");

    let mut pending_pairs: HashMap<Vec<u8>, (i32, i32, i32, bool, u32, u64)> = HashMap::new();
    let mut chunk = Vec::with_capacity(args.batch_size);
    let mut tmp_files = Vec::new();

    // Also collect PE second-end positions during first pass
    let mut pe_second_ends: HashSet<(i32, i32, i32, u8)> = HashSet::new();

    for (index, result) in reader.records().enumerate() {
        let record = result?;
        let flags = record.flags();
        if flags.is_unmapped() || flags.is_secondary() || flags.is_supplementary() {
            continue;
        }

        let lib_id = get_lib_id(&record);
        let pos = get_5p_pos(&record)?;
        let score = get_score(&record);
        let ref_id = record
            .reference_sequence_id()
            .transpose()?
            .map(|i| i as i32)
            .unwrap_or(-1);
        let rev = flags.is_reverse_complemented();

        if flags.is_segmented() && !flags.is_mate_unmapped() {
            let name = record.name().context("no name")?.to_vec();
            if let Some((m_lib, m_ref, m_pos, m_rev, m_score, m_idx)) = pending_pairs.remove(&name)
            {
                let (r1, p1, rv1, i1, r2, p2, rv2, i2) =
                    if (ref_id < m_ref) || (ref_id == m_ref && pos < m_pos) {
                        (ref_id, pos, rev, index as u64, m_ref, m_pos, m_rev, m_idx)
                    } else {
                        (m_ref, m_pos, m_rev, m_idx, ref_id, pos, rev, index as u64)
                    };

                pe_second_ends.insert((m_lib, r2, p2, rv2 as u8));

                chunk.push(Metadata {
                    lib_id: m_lib,
                    ref_id1: r1,
                    pos1: p1,
                    rev1: rv1 as u8,
                    ref_id2: r2,
                    pos2: p2,
                    rev2: rv2 as u8,
                    score: score + m_score,
                    idx1: i1,
                    idx2: i2,
                    paired_end: 1,
                });
                pe_count += 1;
            } else {
                pending_pairs.insert(name, (lib_id, ref_id, pos, rev, score, index as u64));
            }
        } else if flags.is_mate_unmapped() {
            chunk.push(Metadata {
                lib_id,
                ref_id1: ref_id,
                pos1: pos,
                rev1: rev as u8,
                ref_id2: -1,
                pos2: 0,
                rev2: 0,
                score,
                idx1: index as u64,
                idx2: 0,
                paired_end: 0,
            });
            se_count += 1;
        } else {
            chunk.push(Metadata {
                lib_id,
                ref_id1: ref_id,
                pos1: pos,
                rev1: rev as u8,
                ref_id2: -1,
                pos2: 0,
                rev2: 0,
                score,
                idx1: index as u64,
                idx2: 0,
                paired_end: 0,
            });
            se_count += 1;
        }

        if chunk.len() >= args.batch_size {
            let chunk_to_save = std::mem::replace(&mut chunk, Vec::with_capacity(args.batch_size));
            tmp_files.push(io::save_chunk_parallel(chunk_to_save, tmp_dir.path())?);
        }
    }

    // Handle remaining pending pairs
    for (_, (lib, r, p, rv, s, idx)) in pending_pairs {
        chunk.push(Metadata {
            lib_id: lib,
            ref_id1: r,
            pos1: p,
            rev1: rv as u8,
            ref_id2: -1,
            pos2: 0,
            rev2: 0,
            score: s,
            idx1: idx,
            idx2: 0,
            paired_end: 1,
        });
        se_count += 1;
        unmatched_pairs_count += 1;
    }
    if !chunk.is_empty() {
        tmp_files.push(io::save_chunk_parallel(chunk, tmp_dir.path())?);
    }

    eprintln!("  sorted {} end pairs", pe_count);
    eprintln!(
        "     and {} single ends (among them {} unmatched pairs)",
        se_count, unmatched_pairs_count
    );

    // Single pass merge and dedup
    eprint!("  collecting indices of duplicate reads... ");
    let collect_start = Instant::now();
    let mut dup_mask = RoaringBitmap::new();

    let mut heap = BinaryHeap::new();
    let mut readers: Vec<_> = tmp_files
        .iter()
        .map(|p| open_chunk_reader(p))
        .collect();

    for (i, r) in readers.iter_mut().enumerate() {
        if let Some(m) = Metadata::read_from(r)? {
            heap.push(MergeItem { data: m, f_idx: i });
        }
    }

    let mut group: Vec<Metadata> = Vec::with_capacity(1000);
    let mut total_orphan = 0usize;
    let mut total_pe = 0usize;
    let mut total_se_only = 0usize;

    while let Some(item) = heap.pop() {
        if let Some(first) = group.first() {
            let d = &item.data;
            if d.lib_id != first.lib_id
                || d.ref_id1 != first.ref_id1
                || d.pos1 != first.pos1
                || d.rev1 != first.rev1
            {
                let (o, p, s) = identify_dups(&group, &mut dup_mask, &pe_second_ends);
                total_orphan += o;
                total_pe += p;
                total_se_only += s;
                group.clear();
            }
        }
        group.push(item.data);
        if let Some(m) = Metadata::read_from(&mut readers[item.f_idx])? {
            heap.push(MergeItem {
                data: m,
                f_idx: item.f_idx,
            });
        }
    }
    let (o, p, s) = identify_dups(&group, &mut dup_mask, &pe_second_ends);
    total_orphan += o;
    total_pe += p;
    total_se_only += s;

    let collect_dur = collect_start.elapsed();
    eprintln!("done in {} ms", collect_dur.as_millis());
    eprintln!("  found {} duplicates", dup_mask.len());
    eprintln!(
        "  (orphan={}, pe={}, se_only={})",
        total_orphan, total_pe, total_se_only
    );

    let find_dur = find_start.elapsed();
    let (find_m, find_s) = format_duration(find_dur);
    eprintln!(
        "collected list of positions in {} min {} sec",
        find_m, find_s
    );

    // Write output - direct bytes modification
    eprintln!("marking duplicates...");
    let write_start = Instant::now();

    let out_file = File::create(&args.output)?;
    let mut bgzf_writer = BgzfWriter::new(out_file);

    let mut reader = bam::io::reader::Builder::default().build_from_path(&args.input)?;
    reader.read_header()?;

    // Write header using BGZF compression
    write_header(&mut bgzf_writer, &header)?;

    // Read records, modify flag, and write directly
    let mut record_count = 0usize;
    for (idx, result) in reader.records().enumerate() {
        let record = result?;

        // Get raw bytes from record
        let mut data = record_to_bytes(&header, &record)?;

        // Modify flag directly in bytes if not special
        if !record.flags().is_secondary() && !record.flags().is_supplementary() {
            let is_dup = dup_mask.contains(idx as u32);
            toggle_duplicate_flag(&mut data, is_dup);
        }

        bgzf_writer.write_all(&data)?;
        record_count += 1;
    }
    bgzf_writer.finish()?;

    let write_dur = write_start.elapsed();
    eprintln!("wrote output in {:.1} sec", write_dur.as_secs_f64());
    eprintln!("  processed {} records", record_count);

    let total_dur = total_start.elapsed();
    let (total_m, total_s) = format_duration(total_dur);
    eprintln!("done in {} min {} sec", total_m, total_s);

    Ok(())
}
