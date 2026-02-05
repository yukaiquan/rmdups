// Command-line argument parsing
use clap::Parser;

#[derive(Parser, Debug)]
#[command(name = "rmduprs", about = "Sambamba-consistent MarkDuplicates (Rust)")]
pub struct Args {
    #[arg(short, long)]
    pub input: String,
    #[arg(short, long)]
    pub output: String,
    #[arg(short = 'r', long)]
    pub remove_duplicates: bool,
    #[arg(short = 't', long, default_value_t = num_cpus())]
    pub threads: usize,
    #[arg(long, default_value_t = 2_000_000)]
    pub batch_size: usize,
    #[arg(long)]
    pub tmp_dir: Option<std::path::PathBuf>,
    /// Force single-threaded mode (useful for Windows or I/O-bound workloads)
    #[arg(long)]
    pub single_threaded: bool,
}

pub fn num_cpus() -> usize {
    std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(4)
}

/// Get effective thread count based on args and platform
#[inline]
pub fn effective_threads(args: &Args) -> usize {
    if args.single_threaded {
        1
    } else {
        args.threads
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_args_default_threads() {
        let args = Args {
            input: "test.bam".to_string(),
            output: "out.bam".to_string(),
            remove_duplicates: false,
            threads: 4,
            batch_size: 2_000_000,
            tmp_dir: None,
            single_threaded: false,
        };
        assert_eq!(args.input, "test.bam");
        assert_eq!(effective_threads(&args), 4);
    }

    #[test]
    fn test_single_threaded_flag() {
        let args = Args {
            input: "test.bam".to_string(),
            output: "out.bam".to_string(),
            remove_duplicates: false,
            threads: 8,
            batch_size: 2_000_000,
            tmp_dir: None,
            single_threaded: true,
        };
        assert_eq!(effective_threads(&args), 1);
    }
}
