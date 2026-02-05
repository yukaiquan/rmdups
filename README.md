# rmduprs

A high-performance, Sambamba-consistent MarkDuplicates implementation in Rust.

## Features

- **Sambamba-compatible**: Produces identical duplicate marking results
- **Multi-threaded**: Parallel sorting with configurable thread count
- **Cross-platform**: Supports Linux, macOS, and Windows
- **Memory efficient**: Uses RoaringBitmap for duplicate tracking
- **Direct bytes modification**: Avoids RecordBuf clone overhead

## Performance

| Platform | Time | Compared to Sambamba |
|----------|------|---------------------|
| Linux/macOS | ~9-10s | ~30% faster |
| Windows | ~10-12s | ~20% faster (with `--single-threaded`) |

*Test data: 2.2M reads, 27,479 duplicates*

## Installation

```bash
# From source
cargo build --release

# Or install via cargo
cargo install --git https://github.com/yourusername/rmduprs
```

## Usage

```bash
# Basic usage with default threads
rmduprs -i input.bam -o output.bam

# Specify thread count
rmduprs -t 8 -i input.bam -o output.bam

# Single-threaded mode (recommended for Windows)
rmduprs -t 4 --single-threaded -i input.bam -o output.bam

# Use custom temp directory
rmduprs --tmp-dir /path/to/tmp -i input.bam -o output.bam

# Help
rmduprs --help
```

## Options

| Option | Description | Default |
|--------|-------------|---------|
| `-i, --input` | Input BAM file | Required |
| `-o, --output` | Output BAM file | Required |
| `-t, --threads` | Number of threads | CPU count |
| `-r, --remove-duplicates` | Remove duplicates instead of marking | false |
| `--batch-size` | Batch size for sorting | 2,000,000 |
| `--tmp-dir` | Temp directory for intermediate files | System temp |
| `--single-threaded` | Force single-threaded mode | false |
| `-h, --help` | Show help | - |
| `-V, --version` | Show version | - |

## How It Works

### Algorithm Overview

1. **First Pass**: Collect read positions and mate information
   - Parse BAM records to extract key metadata
   - Match paired-end reads by name
   - Save chunks to temporary LZ4-compressed files

2. **Merge & Deduplicate**: Multi-way merge with heap
   - Sort all chunks by (library, position, orientation)
   - Group reads at the same position
   - Identify duplicates using Sambamba's algorithm:
     - **PE duplicates**: Pairs with same mate position/orientation
     - **Orphan fragments**: SE reads in groups with PE reads
     - **SE duplicates**: Multiple SE reads at same position

3. **Mark Duplicates**: Direct bytes modification
   - Serialize records to BAM format
   - Set/clear `0x400` (DUPLICATE) flag at byte offset 12-13
   - Write BGZF-compressed output

### Key Implementation Details

- **Metadata Structure** (43 bytes):
  ```
  lib_id (4) | ref_id1 (4) | pos1 (4) | rev1 rev2 (2)
  ref_id2 (4) | pos2 (4) | score (4) | idx1 (8) | idx2 (8) | paired_end (1)
  ```

- **Duplicate Flag**: Bit 10 in BAM flag (0x400)
- **5' Position Calculation**:
  - Forward reads: `alignment_start - soft-clipped bases`
  - Reverse reads: `alignment_end + soft-clipped bases`

## Architecture

```
src/
├── lib.rs              # Library entry point
├── main.rs             # CLI entry point
├── args.rs             # Command-line arguments
├── metadata.rs         # Metadata struct & serialization
├── algorithm.rs        # Core duplicate detection
├── utils.rs            # Helper functions
└── io/
    └── mod.rs          # BAM I/O utilities
```

## Building for Different Platforms

### Linux/macOS
```bash
cargo build --release
```

### Windows
```bash
# Cross-compile from Linux
cargo build --release --target x86_64-pc-windows-gnu

# Or build directly on Windows
cargo build --release
```

Note: On Windows, using `--single-threaded` often performs better due to reduced thread overhead.

## Library Usage

```rust
use rmduprs::{Args, identify_dups, Metadata};
use roaring::RoaringBitmap;
use std::collections::HashSet;

fn main() {
    let args = Args {
        input: "input.bam".to_string(),
        output: "output.bam".to_string(),
        remove_duplicates: false,
        threads: 8,
        batch_size: 2_000_000,
        tmp_dir: None,
        single_threaded: false,
    };

    // Use the library functions
    let mask = RoaringBitmap::new();
    let pe_second_ends: HashSet<(i32, i32, i32, u8)> = HashSet::new();
    let group = vec![/* Metadata items */];

    let (orphan, pe, se_only) = identify_dups(&group, &mut mask, &pe_second_ends);
}
```

## Testing

```bash
# Run all tests
cargo test

# Run specific module tests
cargo test --lib

# Run with output
cargo test -- --nocapture
```

## Dependencies

- **noodles**: BAM/SAM parsing and writing
- **rayon**: Parallel iteration
- **roaring**: Efficient bitmap operations
- **lz4_flex**: Fast temporary file compression
- **clap**: Command-line argument parsing
- **tempfile**: Secure temporary file handling

## License

MIT License - see LICENSE file for details.

## Acknowledgments

- [Sambamba](https://github.com/biod/Sambamba) - Original algorithm reference
- [noodles](https://github.com/zaeleus/noodles-rs) - Excellent BAM/SAM library

---

# rmduprs

高性能、与 Sambamba 兼容的 Rust 实现版 MarkDuplicates 工具。

## 功能特性

- **Sambamba 兼容**: 产生完全一致的重复标记结果
- **多线程**: 支持可配置线程数的并行排序
- **跨平台**: 支持 Linux、macOS 和 Windows
- **内存高效**: 使用 RoaringBitmap 进行重复标记
- **直接字节修改**: 避免 RecordBuf clone 开销

## 性能表现

| 平台 | 时间 | 与 Sambamba 对比 |
|------|------|------------------|
| Linux/macOS | ~9-10s | 快约 30% |
| Windows | ~10-12s | 快约 20%（使用 `--single-threaded`） |

*测试数据: 220万 reads，27,479 个重复*

## 安装方法

```bash
# 从源码编译
cargo build --release

# 或通过 cargo 安装
cargo install --git https://github.com/yourusername/rmduprs
```

## 使用方法

```bash
# 使用默认线程数
rmduprs -i input.bam -o output.bam

# 指定线程数
rmduprs -t 8 -i input.bam -o output.bam

# 单线程模式（推荐 Windows 使用）
rmduprs -t 4 --single-threaded -i input.bam -o output.bam

# 使用自定义临时目录
rmduprs --tmp-dir /path/to/tmp -i input.bam -o output.bam

# 查看帮助
rmduprs --help
```

## 命令行参数

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `-i, --input` | 输入 BAM 文件 | 必填 |
| `-o, --output` | 输出 BAM 文件 | 必填 |
| `-t, --threads` | 线程数 | CPU 核心数 |
| `-r, --remove-duplicates` | 删除重复而非标记 | false |
| `--batch-size` | 排序批次大小 | 2,000,000 |
| `--tmp-dir` | 中间文件临时目录 | 系统临时目录 |
| `--single-threaded` | 强制单线程模式 | false |
| `-h, --help` | 显示帮助 | - |
| `-V, --version` | 显示版本 | - |

## 工作原理

### 算法概述

1. **第一遍扫描**: 收集 reads 位置和配对信息
   - 解析 BAM 记录提取关键元数据
   - 通过名称匹配双端 reads
   - 将数据块保存到 LZ4 压缩的临时文件

2. **归并去重**: 使用堆的多路归并
   - 按（文库、位置、方向）排序所有数据块
   - 将相同位置的 reads 分组
   - 使用 Sambamba 算法识别重复:
     - **PE 重复**: 配对位置/方向相同的 pairs
     - **Orphan 片段**: 在有 PE reads 的组中的单端 reads
     - **SE 重复**: 相同位置的多个单端 reads

3. **标记重复**: 直接字节修改
   - 将记录序列化为 BAM 格式
   - 在字节偏移 12-13 处设置/清除 `0x400` (DUPLICATE) 标志
   - 写入 BGZF 压缩的输出

### 关键实现细节

- **元数据结构** (43 字节):
  ```
  lib_id (4) | ref_id1 (4) | pos1 (4) | rev1 rev2 (2)
  ref_id2 (4) | pos2 (4) | score (4) | idx1 (8) | idx2 (8) | paired_end (1)
  ```

- **重复标志位**: BAM flag 的第 10 位 (0x400)
- **5' 位置计算**:
  - 正向 reads: `alignment_start - soft-clipped bases`
  - 反向 reads: `alignment_end + soft-clipped bases`

## 项目结构

```
src/
├── lib.rs              # 库入口
├── main.rs             # CLI 入口
├── args.rs             # 命令行参数
├── metadata.rs         # 元数据结构与序列化
├── algorithm.rs        # 核心去重算法
├── utils.rs            # 辅助函数
└── io/
    └── mod.rs          # BAM I/O 工具
```

## 不同平台编译

### Linux/macOS
```bash
cargo build --release
```

### Windows
```bash
# 从 Linux 交叉编译
cargo build --release --target x86_64-pc-windows-gnu

# 或在 Windows 上直接编译
cargo build --release
```

注意: 在 Windows 上，使用 `--single-threaded` 通常表现更好，因为减少了线程开销。

## 作为库使用

```rust
use rmduprs::{Args, identify_dups, Metadata};
use roaring::RoaringBitmap;
use std::collections::HashSet;

fn main() {
    let args = Args {
        input: "input.bam".to_string(),
        output: "output.bam".to_string(),
        remove_duplicates: false,
        threads: 8,
        batch_size: 2_000_000,
        tmp_dir: None,
        single_threaded: false,
    };

    // 使用库函数
    let mask = RoaringBitmap::new();
    let pe_second_ends: HashSet<(i32, i32, i32, u8)> = HashSet::new();
    let group = vec![/* 元数据项 */];

    let (orphan, pe, se_only) = identify_dups(&group, &mut mask, &pe_second_ends);
}
```

## 测试

```bash
# 运行所有测试
cargo test

# 运行特定模块测试
cargo test --lib

# 显示详细输出
cargo test -- --nocapture
```

## 依赖库

- **noodles**: BAM/SAM 解析和写入
- **rayon**: 并行迭代
- **roaring**: 高效位图操作
- **lz4_flex**: 快速临时文件压缩
- **clap**: 命令行参数解析
- **tempfile**: 安全临时文件处理

## 开源协议

MIT License - 详见 LICENSE 文件。

## 致谢

- [Sambamba](https://github.com/biod/Sambamba) - 原始算法参考
- [noodles](https://github.com/zaeleus/noodles-rs) - 优秀的 BAM/SAM 库
