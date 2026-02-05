//! BAM file I/O utilities
//!
//! This module provides utilities for reading and writing BAM files,
//! including header handling and flag modification.

use anyhow::Result;
use lz4_flex::frame::{FrameDecoder, FrameEncoder};
use noodles::bam;
use noodles::bgzf::io::Writer as BgzfWriter;
use noodles::sam::alignment::io::Write as SamWrite;
use noodles::sam::header::Header as SamHeader;
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::path::Path;

/// Offset of the flag field in a BAM record's binary format
///
/// The flag field is at bytes 12-13 (after ref_id=4 + pos=4 + bin_mq_nl=4)
pub const FLAG_OFFSET: usize = 12;

/// The DUPLICATE flag bit in BAM format
pub const DUPLICATE_FLAG: u16 = 0x400;

/// Modify the DUPLICATE flag in raw BAM record bytes
///
/// Returns the modified flag value.
#[inline]
pub fn toggle_duplicate_flag(data: &mut [u8], is_duplicate: bool) -> Option<u16> {
    if data.len() < FLAG_OFFSET + 2 {
        return None;
    }

    let flag = u16::from_le_bytes([data[FLAG_OFFSET], data[FLAG_OFFSET + 1]]);
    let new_flag = if is_duplicate {
        flag | DUPLICATE_FLAG
    } else {
        flag & !DUPLICATE_FLAG
    };

    data[FLAG_OFFSET] = new_flag as u8;
    data[FLAG_OFFSET + 1] = (new_flag >> 8) as u8;

    Some(new_flag)
}

/// Check if a record is a duplicate based on index
#[inline]
pub fn is_duplicate(idx: usize, dup_mask: &roaring::RoaringBitmap) -> bool {
    dup_mask.contains(idx as u32)
}

/// Write header to BGZF-compressed BAM file
pub fn write_header(
    writer: &mut BgzfWriter<File>,
    header: &SamHeader,
) -> Result<()> {
    let mut header_buf = Vec::new();
    {
        let mut writer = bam::io::Writer::from(&mut header_buf);
        writer.write_header(header)?;
    }
    writer.write_all(&header_buf)?;
    writer.flush()?;
    Ok(())
}

/// Serialize a BAM record to raw bytes
pub fn record_to_bytes(
    header: &SamHeader,
    record: &bam::Record,
) -> Result<Vec<u8>> {
    let mut data = Vec::new();
    {
        let mut writer = bam::io::Writer::from(&mut data);
        writer.write_alignment_record(header, record)?;
    }
    Ok(data)
}

/// Write a BAM record with optional duplicate flag modification
#[allow(dead_code)]
pub fn write_record_with_dup_flag(
    writer: &mut BgzfWriter<File>,
    header: &SamHeader,
    record: &bam::Record,
    idx: usize,
    dup_mask: &roaring::RoaringBitmap,
) -> Result<()> {
    // Skip secondary/supplementary records for duplicate marking
    if record.flags().is_secondary() || record.flags().is_supplementary() {
        let data = record_to_bytes(header, record)?;
        writer.write_all(&data)?;
        return Ok(());
    }

    let mut data = record_to_bytes(header, record)?;
    let is_dup = is_duplicate(idx, dup_mask);
    toggle_duplicate_flag(&mut data, is_dup);
    writer.write_all(&data)?;
    Ok(())
}

/// Parallel chunk saving with LZ4 compression
///
/// Sorts the chunk in parallel before saving.
pub fn save_chunk_parallel(
    mut chunk: Vec<super::metadata::Metadata>,
    dir: &Path,
) -> Result<std::path::PathBuf> {
    chunk.par_sort_unstable();
    let path = dir.join(format!("{}.lz4", fastrand::u64(..)));
    let mut enc = FrameEncoder::new(BufWriter::with_capacity(1 << 20, File::create(&path)?));
    for m in chunk {
        m.write_to(&mut enc)?;
    }
    enc.finish()?;
    Ok(path)
}

/// Open a chunk file for reading
pub fn open_chunk_reader(path: &Path) -> BufReader<FrameDecoder<File>> {
    BufReader::with_capacity(1 << 18, FrameDecoder::new(File::open(path).unwrap()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_flag_offset_constant() {
        // Verify flag offset matches BAM spec
        // ref_id (4 bytes) + pos (4 bytes) + bin_mq_nl (4 bytes) = 12
        assert_eq!(FLAG_OFFSET, 12);
    }

    #[test]
    fn test_duplicate_flag_constant() {
        // 0x400 = 1024 = bit 10 (DUPLICATE flag in SAM/BAM)
        assert_eq!(DUPLICATE_FLAG, 0x400);
    }

    #[test]
    fn test_toggle_duplicate_flag_set() {
        let mut data = [0u8; 20];
        data[12] = 0x00; // flag = 0
        data[13] = 0x00;

        let result = toggle_duplicate_flag(&mut data, true);
        assert_eq!(result, Some(0x400));
        assert_eq!(u16::from_le_bytes([data[12], data[13]]), 0x400);
    }

    #[test]
    fn test_toggle_duplicate_flag_clear() {
        let mut data = [0u8; 20];
        data[12] = 0x00;
        data[13] = 0x04; // flag = 0x400

        let result = toggle_duplicate_flag(&mut data, false);
        assert_eq!(result, Some(0x000));
        assert_eq!(u16::from_le_bytes([data[12], data[13]]), 0x000);
    }

    #[test]
    fn test_toggle_duplicate_flag_preserve_other_bits() {
        let mut data = [0u8; 20];
        data[12] = 0x02; // flag = 0x402 (PAIRED | DUPLICATE)
        data[13] = 0x00;

        let result = toggle_duplicate_flag(&mut data, true);
        assert_eq!(result, Some(0x402)); // PAIRED | DUPLICATE
        assert_eq!(data[12], 0x02);
    }

    #[test]
    fn test_toggle_duplicate_flag_insufficient_data() {
        let mut data = [0u8; 12]; // Too short
        let result = toggle_duplicate_flag(&mut data, true);
        assert!(result.is_none());
    }
}
