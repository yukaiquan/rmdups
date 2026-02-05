//! Metadata for duplicate detection - matches Sambamba's comparator order
//!
//! This module defines the Metadata struct that stores read information
//! for duplicate detection, with serialization support for temporary files.

use anyhow::Result;
use std::io::{Read, Write};

/// Metadata for a read or read pair used in duplicate detection
///
/// The ordering of fields matches Sambamba's markdup comparator:
/// lib_id -> ref_id1 -> pos1 -> rev1 -> ref_id2 -> pos2 -> rev2 -> score
#[derive(Debug, PartialEq, Eq, PartialOrd, Ord, Clone)]
pub struct Metadata {
    pub lib_id: i32,
    pub ref_id1: i32,
    pub pos1: i32,
    pub rev1: u8,
    pub rev2: u8,
    pub ref_id2: i32,
    pub pos2: i32,
    pub score: u32,
    pub idx1: u64,
    pub idx2: u64,
    pub paired_end: u8, // 0 = SE/fragment, 1 = PE/second end
}

impl Metadata {
    /// Create new metadata for a single-end read
    #[inline]
    pub fn new_se(
        lib_id: i32,
        ref_id1: i32,
        pos1: i32,
        rev1: u8,
        score: u32,
        idx1: u64,
    ) -> Self {
        Self {
            lib_id,
            ref_id1,
            pos1,
            rev1,
            rev2: 0,
            ref_id2: -1,
            pos2: 0,
            score,
            idx1,
            idx2: 0,
            paired_end: 0,
        }
    }

    /// Create new metadata for a paired-end read
    #[inline]
    pub fn new_pe(
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
    ) -> Self {
        Self {
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

    /// Serialize metadata to binary format (little-endian)
    pub fn write_to<W: Write>(&self, w: &mut W) -> Result<()> {
        w.write_all(&self.lib_id.to_le_bytes())?;
        w.write_all(&self.ref_id1.to_le_bytes())?;
        w.write_all(&self.pos1.to_le_bytes())?;
        w.write_all(&[self.rev1, self.rev2])?;
        w.write_all(&self.ref_id2.to_le_bytes())?;
        w.write_all(&self.pos2.to_le_bytes())?;
        w.write_all(&self.score.to_le_bytes())?;
        w.write_all(&self.idx1.to_le_bytes())?;
        w.write_all(&self.idx2.to_le_bytes())?;
        w.write_all(&[self.paired_end])?;
        Ok(())
    }

    /// Deserialize metadata from binary format (little-endian)
    ///
    /// Returns `Ok(None)` if end of stream is reached.
    pub fn read_from<R: Read>(r: &mut R) -> Result<Option<Self>> {
        let mut buf4 = [0u8; 4];
        if r.read_exact(&mut buf4).is_err() {
            return Ok(None);
        }
        let lib_id = i32::from_le_bytes(buf4);

        r.read_exact(&mut buf4)?;
        let ref_id1 = i32::from_le_bytes(buf4);
        r.read_exact(&mut buf4)?;
        let pos1 = i32::from_le_bytes(buf4);

        let mut buf2 = [0u8; 2];
        r.read_exact(&mut buf2)?;
        let (rev1, rev2) = (buf2[0], buf2[1]);

        r.read_exact(&mut buf4)?;
        let ref_id2 = i32::from_le_bytes(buf4);
        r.read_exact(&mut buf4)?;
        let pos2 = i32::from_le_bytes(buf4);
        r.read_exact(&mut buf4)?;
        let score = u32::from_le_bytes(buf4);

        let mut buf8 = [0u8; 8];
        r.read_exact(&mut buf8)?;
        let idx1 = u64::from_le_bytes(buf8);
        r.read_exact(&mut buf8)?;
        let idx2 = u64::from_le_bytes(buf8);

        let mut buf1 = [0u8; 1];
        r.read_exact(&mut buf1)?;
        let paired_end = buf1[0];

        Ok(Some(Self {
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
            paired_end,
        }))
    }

    /// Get the binary size of metadata
    pub fn binary_size() -> usize {
        4 + 4 + 4 + 2 + 4 + 4 + 4 + 8 + 8 + 1 // 43 bytes
    }
}

/// Merge item for heap-based multi-way merge
#[derive(Eq, PartialEq)]
pub struct MergeItem {
    pub data: Metadata,
    pub f_idx: usize,
}

impl Ord for MergeItem {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        other.data.cmp(&self.data)
    }
}

impl PartialOrd for MergeItem {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_metadata_se_roundtrip() {
        let meta = Metadata::new_se(1, 0, 1000, 0, 50, 42);
        let mut buf = Vec::new();
        meta.write_to(&mut buf).unwrap();

        let mut cursor = Cursor::new(buf);
        let read_meta = Metadata::read_from(&mut cursor).unwrap().unwrap();

        assert_eq!(meta, read_meta);
        assert_eq!(meta.paired_end, 0);
        assert_eq!(meta.idx2, 0);
    }

    #[test]
    fn test_metadata_pe_roundtrip() {
        let meta = Metadata::new_pe(
            1,   // lib_id
            0,   // ref_id1
            1000, // pos1
            0,   // rev1
            0,   // ref_id2
            1050, // pos2
            1,   // rev2
            50,  // score
            42,  // idx1
            43,  // idx2
        );
        let mut buf = Vec::new();
        meta.write_to(&mut buf).unwrap();

        let mut cursor = Cursor::new(buf);
        let read_meta = Metadata::read_from(&mut cursor).unwrap().unwrap();

        assert_eq!(meta, read_meta);
        assert_eq!(meta.paired_end, 1);
        assert_eq!(meta.idx2, 43);
    }

    #[test]
    fn test_metadata_binary_size() {
        assert_eq!(Metadata::binary_size(), 43);
    }

    #[test]
    fn test_metadata_ordering() {
        // Test that ordering matches Sambamba's comparator
        // Ordering: lib_id -> ref_id1 -> pos1 -> rev1 -> ref_id2 -> pos2 -> rev2 -> score
        let m1 = Metadata::new_se(0, 0, 100, 0, 50, 1);
        let m2 = Metadata::new_se(0, 0, 200, 0, 50, 2);

        // m1.pos1=100 < m2.pos1=200, so m1 < m2 (for min-heap)
        assert!(m1 < m2);

        // m3 has same pos as m1 but rev1=1 (forward < reverse)
        let m3 = Metadata::new_se(0, 0, 100, 1, 50, 3);
        assert!(m1 < m3); // same pos, but rev1: 0 < 1
        assert!(m3 < m2); // m3.pos1=100 < m2.pos1=200
    }

    #[test]
    fn test_metadata_read_from_empty() {
        let mut cursor = Cursor::new(Vec::new());
        let result = Metadata::read_from(&mut cursor).unwrap();
        assert!(result.is_none());
    }
}
