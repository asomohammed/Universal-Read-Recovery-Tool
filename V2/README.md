# Universal Read Recovery Tool v2.0 🧬

## 🎯 What's New in v2.0

- **Modular Architecture**: Completely restructured codebase for better maintainability
- **Beautiful Reports**: Redesigned HTML reports with interactive visualizations

## 🚀 Features

### Core Features

- **Smart Barcode Discovery**: Automatic identification of sample barcodes
- **Advanced Fuzzy Matching**: Sophisticated algorithm for barcode error correction
- **Multi-Index Support**: Handles single, dual, and custom indexing schemes
- **Quality-Aware Processing**: Considers quality scores in recovery decisions
- **Comprehensive Reporting**: Detailed HTML and CSV reports with visualizations


## 📋 Requirements

- Python 3.8 or higher
- 4GB RAM minimum (8GB recommended)
- Operating System:
  - Linux (Recommended)
  - macOS 10.14+
  - Windows 10/11

### Dependencies

```bash
numpy>=1.20.0
pandas>=1.3.0
matplotlib>=3.4.0
seaborn>=0.11.0
```

## 🔧 Installation

### Using pip (Recommended)

```bash
pip install universal-read-recovery
```

### From source

```bash
git clone https://github.com/yourusername/Universal-Read-Recovery-Tool.git
cd Universal-Read-Recovery-Tool
pip install -e .
```

### Docker

```bash
docker pull yourusername/universal-read-recovery:2.0
docker run -v /path/to/data:/data yourusername/universal-read-recovery:2.0
```

## 📖 Usage

### Quick Start

```bash
universal-read-recovery -d /path/to/fastq/files
```

### Advanced Usage

```bash
universal-read-recovery \
    -d /path/to/fastq/files \
    -m 2 \                    # Maximum mismatches
    -s 1000000 \             # Sample size
    -r 10000 \               # Discovery reads
    -t 4 \                   # Number of threads
    --quality-threshold 30 \  # Minimum quality score
    -o custom_report.csv     # Report name
```

### Configuration File

```yaml
# config.yaml
input:
  directory: /path/to/fastq
  pattern: "*_R{1,2}*.fastq.gz"

parameters:
  max_mismatches: 2
  sample_reads: 1000000
  discovery_reads: 10000
  threads: 4
  quality_threshold: 30

output:
  report_name: recovery_report
  format: [html, csv]
  save_intermediates: false
```

## 📊 Output Structure

```
output_directory/
├── recovered/
│   ├── Sample1_R1.fastq.gz
│   ├── Sample1_R2.fastq.gz
│   └── ...
├── reports/
│   ├── recovery_report.html
│   ├── recovery_report.csv
│   └── recovery_metrics.json
└── logs/
    ├── recovery.log
    └── performance_metrics.log
```

## 🔍 Advanced Features

### Quality Score Integration

```bash
universal-read-recovery --quality-aware \
    --min-quality 30 \
    --quality-window 5
```

### Custom Barcode Patterns

```bash
universal-read-recovery --custom-patterns patterns.txt \
    --pattern-format illumina
```

### Performance Tuning

```bash
universal-read-recovery --threads 8 \
    --chunk-size 100000 \
    --buffer-size 1GB
```

## 🛠 Development

### Project Structure

```
Universal-Read-Recovery-Tool/
├── src/
│   └── universal_read_recovery/
│       ├── __init__.py
│       ├── __main__.py
│       ├── barcode_utils.py
│       ├── file_utils.py
│       ├── analysis_utils.py
│       └── report_utils.py
├── tests/
├── docs/
└── examples/
```

### Running Tests

```bash
pytest tests/
pytest tests/ --cov=universal_read_recovery
```

### Building Documentation

```bash
cd docs
make html
```

## 📈 Performance Metrics

| Dataset Size | Processing Time | Memory Usage | Recovery Rate |
| ------------ | --------------- | ------------ | ------------- |
| 1M reads     | 2 min           | 0.5 GB       | 15-20%        |
| 10M reads    | 15 min          | 2 GB         | 12-18%        |
| 100M reads   | 2 hours         | 8 GB         | 10-15%        |


### Development Setup

```bash
# Create virtual environment
python -m venv venv
source venv/bin/activate  # or `venv\Scripts\activate` on Windows

# Install development dependencies
pip install -r requirements-dev.txt

# Run pre-commit hooks
pre-commit install
```

## 📝 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- The Bioinformatics Community
- Contributors and Testers
- Supporting Institutions

### Common Issues

1. **Memory Errors**

   ```bash
   universal-read-recovery --low-memory \
       --chunk-size 50000
   ```

2. **Performance Issues**
   ```bash
   universal-read-recovery --optimize-speed \
       --threads 4
   ```

---

```
