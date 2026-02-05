//! Sambamba-consistent duplicate detection algorithm
//!
//! This module implements the core duplicate detection logic that matches
// Sambamba's markdup algorithm behavior.

use anyhow::Result;
use noodles::bam;
use noodles::sam::alignment::record::cigar::op::Kind;
use roaring::RoaringBitmap;
use std::collections::HashSet;

/// Calculate the 5' position of a read
///
/// For forward reads, this is the alignment start minus soft-clipped bases.
/// For reverse reads, this is the alignment end plus soft-clipped bases.
pub fn get_5p_pos(record: &bam::Record) -> Result<i32> {
    let start = record
        .alignment_start()
        .transpose()?
        .map(|p| p.get() as i32 - 1)
        .unwrap_or(-1);
    if start < 0 {
        return Ok(-1);
    }
    let cigar = record.cigar();

    if !record.flags().is_reverse_complemented() {
        let mut clipped = 0;
        for op in cigar.iter() {
            let op = op?;
            match op.kind() {
                Kind::SoftClip | Kind::HardClip => clipped += op.len() as i32,
                _ => break,
            }
        }
        Ok(start - clipped)
    } else {
        let mut ref_span = 0;
        for op in cigar.iter() {
            let op = op?;
            if op.kind().consumes_reference() {
                ref_span += op.len() as i32;
            }
        }
        let ops: Vec<_> = cigar.iter().collect::<Result<_, _>>()?;
        let mut clipped_end = 0;
        for op in ops.iter().rev() {
            match op.kind() {
                Kind::SoftClip | Kind::HardClip => clipped_end += op.len() as i32,
                _ => break,
            }
        }
        Ok(start + ref_span + clipped_end)
    }
}

/// Calculate the duplicate scoring metric
///
/// Sum of quality scores >= 15. This is used to select the best copy
/// when multiple duplicates exist.
#[inline]
pub fn get_score(record: &bam::Record) -> u32 {
    record
        .quality_scores()
        .as_ref()
        .iter()
        .map(|&q| u8::from(q))
        .filter(|&q| q >= 15)
        .map(|q| q as u32)
        .sum()
}

/// Identify duplicates within a group of reads with the same position
///
/// Returns a tuple of (orphan_count, pe_count, se_only_count) for the group.
///
/// - **orphan**: SE read in a group that also has PE reads
/// - **pe**: PE read where both reads have another duplicate pair
/// - **se_only**: SE read where no PE reads exist in the group
pub fn identify_dups(
    group: &[super::metadata::Metadata],
    mask: &mut RoaringBitmap,
    pe_second_ends: &HashSet<(i32, i32, i32, u8)>,
) -> (usize, usize, usize) {
    if group.is_empty() {
        return (0, 0, 0);
    }

    let mut orphan_marked = 0;
    let mut pe_marked = 0;
    let mut se_only_marked = 0;

    let (pes, ses): (Vec<_>, Vec<_>) = group.iter().partition(|m| m.ref_id2 != -1);

    // paired_end == 0: fragment (read with unmapped mate or SE read)
    // paired_end == 1: PE second end (mate is also in this group)
    let paired_0: Vec<_> = ses.iter().filter(|se| se.paired_end == 0).collect();
    let paired_1: Vec<_> = ses.iter().filter(|se| se.paired_end == 1).collect();

    let k_pe = pes.len();
    let group_pos = (
        group[0].lib_id,
        group[0].ref_id1,
        group[0].pos1,
        group[0].rev1,
    );
    let k_pos = if pe_second_ends.contains(&group_pos) {
        1
    } else {
        0
    };

    let k = paired_0.len() + paired_1.len();
    let total = k + k_pe + k_pos;

    let seen_fragment = !paired_0.is_empty();
    let seen_paired_read = !paired_1.is_empty() || k_pe > 0 || k_pos > 0;

    // SE-only deduplication logic
    if total >= 2 && seen_fragment {
        if seen_paired_read {
            // Orphan handling: mark fragments when paired reads exist
            for se in &paired_0 {
                mask.insert(se.idx1 as u32);
                orphan_marked += 1;
            }
        } else if paired_0.len() >= 2 {
            // Fragment deduplication: keep highest scoring
            let mut best_idx = 0;
            for i in 1..paired_0.len() {
                if paired_0[i].score > paired_0[best_idx].score {
                    best_idx = i;
                }
            }
            for i in 0..paired_0.len() {
                if i != best_idx {
                    mask.insert(paired_0[i].idx1 as u32);
                    se_only_marked += 1;
                }
            }
        }
    }

    // PE internal deduplication
    if pes.len() >= 2 {
        let mut i = 0;
        while i < pes.len() {
            let mut j = i + 1;
            let mut best_idx = i;
            // Find reads with same mate position/orientation
            while j < pes.len()
                && pes[i].rev2 == pes[j].rev2
                && pes[i].ref_id2 == pes[j].ref_id2
                && pes[i].pos2 == pes[j].pos2
            {
                if pes[j].score >= pes[best_idx].score {
                    best_idx = j;
                }
                j += 1;
            }
            // Mark all except best scoring
            for k in i..j {
                if k != best_idx {
                    mask.insert(pes[k].idx1 as u32);
                    mask.insert(pes[k].idx2 as u32);
                    pe_marked += 2;
                }
            }
            i = j;
        }
    }

    (orphan_marked, pe_marked, se_only_marked)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::metadata::Metadata;
    use roaring::RoaringBitmap;
    use std::collections::HashSet;

    /// Create a test metadata for SE read
    fn make_se(
        lib_id: i32,
        ref_id: i32,
        pos: i32,
        rev: u8,
        score: u32,
        idx: u64,
        paired_end: u8,
    ) -> Metadata {
        Metadata {
            lib_id,
            ref_id1: ref_id,
            pos1: pos,
            rev1: rev,
            rev2: 0,
            ref_id2: -1,
            pos2: 0,
            score,
            idx1: idx,
            idx2: 0,
            paired_end,
        }
    }

    /// Create a test metadata for PE read
    fn make_pe(
        lib_id: i32,
        ref_id1: i32,
        pos1: i32,
        rev1: u8,
        ref_id2: i32,
        pos2: i32,
        rev2: u8,
        score: u32,
        idx1: u64,
        idx2: u64,
    ) -> Metadata {
        Metadata {
            lib_id,
            ref_id1,
            pos1,
            rev1,
            rev2,
            ref_id2,
            pos2,
            score,
            idx1,
            idx2,
            paired_end: 1,
        }
    }

    #[test]
    fn test_empty_group() {
        let mask = &mut RoaringBitmap::new();
        let pe_second_ends: HashSet<(i32, i32, i32, u8)> = HashSet::new();
        let (orphan, pe, se_only) = identify_dups(&[], mask, &pe_second_ends);
        assert_eq!((orphan, pe, se_only), (0, 0, 0));
    }

    #[test]
    fn test_single_read_not_marked() {
        let group = vec![make_se(0, 0, 100, 0, 50, 0, 0)];
        let mask = &mut RoaringBitmap::new();
        let pe_second_ends: HashSet<(i32, i32, i32, u8)> = HashSet::new();
        let (orphan, pe, se_only) = identify_dups(&group, mask, &pe_second_ends);
        assert_eq!((orphan, pe, se_only), (0, 0, 0));
        assert!(mask.is_empty());
    }

    #[test]
    fn test_fragment_deduplication() {
        // Two SE reads at same position, no PE reads
        let group = vec![
            make_se(0, 0, 100, 0, 50, 0, 0),
            make_se(0, 0, 100, 0, 70, 1, 0), // higher score, should be kept
            make_se(0, 0, 100, 0, 40, 2, 0),
        ];
        let mask = &mut RoaringBitmap::new();
        let pe_second_ends: HashSet<(i32, i32, i32, u8)> = HashSet::new();
        let (orphan, pe, se_only) = identify_dups(&group, mask, &pe_second_ends);
        assert_eq!((orphan, pe, se_only), (0, 0, 2));
        assert_eq!(mask.len(), 2);
        assert!(!mask.contains(1)); // best one not marked
    }

    #[test]
    fn test_orphan_handling() {
        // SE fragment + PE reads = orphan
        let group = vec![
            make_se(0, 0, 100, 0, 50, 0, 0),            // fragment (paired_end=0)
            make_pe(0, 0, 100, 0, 1, 200, 1, 60, 1, 2), // PE
        ];
        let mask = &mut RoaringBitmap::new();
        let pe_second_ends: HashSet<(i32, i32, i32, u8)> = HashSet::new();
        let (orphan, pe, se_only) = identify_dups(&group, mask, &pe_second_ends);
        assert_eq!((orphan, pe, se_only), (1, 0, 0));
        assert!(mask.contains(0));
    }

    #[test]
    fn test_pe_deduplication() {
        // Two PE pairs at same position
        let group = vec![
            make_pe(0, 0, 100, 0, 1, 200, 1, 70, 0, 1), // higher score, kept
            make_pe(0, 0, 100, 0, 1, 200, 1, 50, 2, 3), // lower score, marked
        ];
        let mask = &mut RoaringBitmap::new();
        let pe_second_ends: HashSet<(i32, i32, i32, u8)> = HashSet::new();
        let (orphan, pe, se_only) = identify_dups(&group, mask, &pe_second_ends);
        assert_eq!((orphan, pe, se_only), (0, 2, 0));
        assert_eq!(mask.len(), 2);
        assert!(!mask.contains(0));
        assert!(!mask.contains(1));
    }

    #[test]
    fn test_different_library_separate() {
        // Note: identify_dups assumes all reads in the group have the same
        // (lib_id, ref_id1, pos1, rev1). Reads with different lib_ids are
        // processed in separate groups by the main algorithm.
        // This test verifies same-lib behavior.
        let group = vec![
            make_se(0, 0, 100, 0, 50, 0, 0),
            make_se(0, 0, 100, 0, 60, 1, 0), // same lib_id
        ];
        let mask = &mut RoaringBitmap::new();
        let pe_second_ends: HashSet<(i32, i32, i32, u8)> = HashSet::new();
        let (orphan, pe, se_only) = identify_dups(&group, mask, &pe_second_ends);
        // Same position, same library, should mark one as duplicate
        assert_eq!((orphan, pe, se_only), (0, 0, 1));
        assert_eq!(mask.len(), 1);
    }

    #[test]
    fn test_pe_second_ends_contains_check() {
        // PE second end check via pe_second_ends set
        let group = vec![
            make_se(0, 0, 100, 0, 50, 0, 0), // SE fragment
        ];
        let mut pe_second_ends = HashSet::new();
        pe_second_ends.insert((0, 0, 100, 0)); // This read IS a PE second end

        let mask = &mut RoaringBitmap::new();
        let (orphan, pe, se_only) = identify_dups(&group, mask, &pe_second_ends);
        assert_eq!((orphan, pe, se_only), (1, 0, 0));
    }
}
