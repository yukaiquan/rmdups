//! Utility functions
//!
//! Common helper functions used throughout the project.

use std::time::Duration;

/// Format a duration into (minutes, seconds) tuple
///
/// Useful for printing elapsed time in human-readable format.
#[inline]
pub fn format_duration(dur: Duration) -> (u64, u64) {
    let secs = dur.as_secs();
    (secs / 60, secs % 60)
}

/// Format duration as a human-readable string
#[inline]
pub fn format_duration_verbose(dur: Duration) -> String {
    let secs = dur.as_secs();
    if secs >= 60 {
        format!("{} min {} sec", secs / 60, secs % 60)
    } else {
        format!("{:.1} sec", dur.as_secs_f64())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_format_duration_seconds() {
        let dur = Duration::from_secs(45);
        assert_eq!(format_duration(dur), (0, 45));
    }

    #[test]
    fn test_format_duration_minutes() {
        let dur = Duration::from_secs(125); // 2 min 5 sec
        assert_eq!(format_duration(dur), (2, 5));
    }

    #[test]
    fn test_format_duration_verbose_seconds() {
        let dur = Duration::from_millis(500);
        let result = format_duration_verbose(dur);
        assert!(result.contains("sec"));
    }

    #[test]
    fn test_format_duration_verbose_minutes() {
        let dur = Duration::from_secs(125);
        let result = format_duration_verbose(dur);
        assert_eq!(result, "2 min 5 sec");
    }
}
