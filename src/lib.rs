//! rmduprs - A Sambamba-consistent MarkDuplicates implementation in Rust
//!
//! This library provides efficient duplicate detection and marking for BAM files.
//! The algorithm is designed to produce results identical to Sambamba's markdup.
//!
//! # Example
//!
//! ```ignore
//! use rmduprs::{Args, run_markdup};
//!
//! let args = Args {
//!     input: "input.bam".to_string(),
//!     output: "output.bam".to_string(),
//!     remove_duplicates: false,
//!     threads: 8,
//!     batch_size: 2_000_000,
//!     tmp_dir: None,
//! };
//!
//! run_markdup(&args)?;

pub mod algorithm;
pub mod args;
pub mod io;
pub mod metadata;
pub mod utils;

// Re-export commonly used items
pub use algorithm::{get_5p_pos, get_score, identify_dups};
pub use args::Args;
pub use io::{DUPLICATE_FLAG, FLAG_OFFSET, toggle_duplicate_flag};
pub use metadata::Metadata;
