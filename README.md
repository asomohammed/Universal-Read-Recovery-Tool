# Universal Read Recovery 🧬

A powerful Python tool for recovering misassigned reads from Illumina sequencing demultiplexing, specifically designed to rescue reads that were incorrectly placed in the "Unassigned" or "Undetermined" files due to sequencing errors in barcode sequences.

## 🎯 What it does

- **Automatically discovers** barcode patterns from your assigned sample files
- **Fuzzy matches** unassigned reads against known samples using configurable mismatch tolerance
- **Recovers and combines** rescued reads with original sample files
- **Generates professional reports** in both HTML and CSV formats
- **Preserves original filenames** in the recovered dataset

## 🚀 Quick Start

```bash
# Basic usage - analyze current directory with default settings
python3 universal_read_recovery.py

# Custom parameters
python3 universal_read_recovery.py -m 3 -s 2000000 -r 20000

# Specify directory
python3 universal_read_recovery.py -d /path/to/fastq/files
```

## 📋 Requirements

- Python 3.6+
- Standard libraries only (no additional installations required!)
- Gzipped FASTQ files from Illumina sequencing or anyother sequencing platforms 

## 🔧 Installation

```bash
# Clone the repository
git clone https://github.com/yasomohammed/Universal-Read-Recovery-Tool.git
cd universal-read-recovery

# Make executable
chmod +x universal_read_recovery.py

# Run directly
python3 universal_read_recovery.py --help
```

## 📊 Input Files Expected

Your directory should contain:

```
your_data/
├── Sample1_R1_001.fastq.gz    # Assigned samples
├── Sample1_R2_001.fastq.gz
├── Sample2_R1_001.fastq.gz
├── Sample2_R2_001.fastq.gz
├── ...
├── Unassigned_R1_001.fastq.gz # Unassigned reads
└── Unassigned_R2_001.fastq.gz
```

**Note**: Also works with "Undetermined" files from different Illumina software versions.

## 🎛️ Command Line Options

| Option | Default | Description |
|--------|---------|-------------|
| `-m, --max-mismatches` | 2 | Maximum barcode mismatches allowed |
| `-s, --sample-reads` | 1,000,000 | Number of unassigned reads to analyze |
| `-r, --discovery-reads` | 10,000 | Reads per sample for barcode discovery |
| `-d, --directory` | `.` | Directory containing FASTQ files |
| `-o, --output` | `recovery_report.csv` | CSV report filename |

## 📈 Output

### Files Created:
- `recovered/` - Directory with combined original + recovered reads
- `recovery_report.html` - Professional HTML report with tables and statistics
- `recovery_report.csv` - Machine-readable CSV summary

### Example Output Structure:
```
your_data/
├── recovered/
│   ├── Sample1_R1_001.fastq.gz  # Original + recovered reads
│   ├── Sample1_R2_001.fastq.gz
│   └── ...
├── recovery_report.html
└── recovery_report.csv
```

## 🧪 How It Works

1. **Discovery Phase**: Samples reads from each assigned file to learn barcode patterns
2. **Analysis Phase**: Examines unassigned reads to identify recoverable barcodes
3. **Matching Phase**: Uses fuzzy matching with configurable mismatch tolerance
4. **Recovery Phase**: Extracts matching read pairs from unassigned files
5. **Combination Phase**: Merges recovered reads with original sample files

## 📊 Example Results

```
Recovery Summary:
  Total sampled reads: 1,000,000
  Recoverable reads: 45,231 (4.5%)
  Unique recoverable barcodes: 127

Sample Results:
  Sample1: 12,445 recovered reads (2.1% increase)
  Sample2: 8,932 recovered reads (1.8% increase)
  Sample3: 15,678 recovered reads (3.2% increase)
```

## 🎯 Use Cases

- **Low-quality sequencing runs** with high error rates
- **Dual-index libraries** with sequencing errors in one or both indices  
- **Salvaging data** from runs with poor demultiplexing performance
- **Maximizing sample yields** from expensive sequencing experiments
- **Quality control** analysis of demultiplexing efficiency

## ⚡ Performance Tips

- **Large datasets**: Reduce `-s` to 500,000 for faster analysis
- **High accuracy**: Increase `-r` to 50,000 for better barcode discovery
- **Noisy data**: Start with `-m 1`, then increase if recovery is low
- **Clean data**: Use `-m 3` for maximum recovery

## 🔍 Troubleshooting

### No assigned samples found
- Ensure FASTQ files follow Illumina naming convention
- Check that files contain `R1` and don't start with "Unassigned"

### Low recovery rates
- Try increasing `--max-mismatches`
- Check barcode quality in your sequencing run
- Verify sample sheet was correct during demultiplexing

### Memory issues
- Reduce `--sample-reads` parameter
- Process smaller batches of samples

## 🤝 Contributing

Contributions welcome! Please feel free to submit issues, feature requests, or pull requests.

### Development Setup
```bash
git clone https://github.com/asomohammed/Universal-Read-Recovery-Tool.git
cd universal-read-recovery
# Make your changes
python3 universal_read_recovery.py  # Test locally
```

## 🙏 Citation

If you use this tool in your research, please cite:

```
Universal Read Recovery Tool
(https://github.com/asomohammed/Universal-Read-Recovery-Tool)

```

## 📞 Support

- 🐛 **Bug reports**: [Open an issue](https://github.com/asomohammed/Universal-Read-Recovery-Tool/)
- 💡 **Feature requests**: [Open an issue](https://github.com/asomohammed/Universal-Read-Recovery-Tool/issues)
- 📧 **Questions**: [Discussions](https://github.com/asomohammed/Universal-Read-Recovery-Tool/discussions)

## 🏆 Acknowledgments

- Inspired by the need to maximize data recovery from expensive sequencing runs
- Built for the genomics community to reduce data waste
- Designed with bioinformatics workflows in mind

---

**⭐ If this tool helped rescue your data, please star the repository!**


## Authors

**Aso Omer Mohammed**  
3DBM, Neurosurgery, University Hospital Freiburg  
GitHub: [@asomohammed](https://github.com/asomohammed)  
Email: aso.mohammed@uniklinik-freiburg.de

**Dr. Ing. Kevin Joseph**  
Neurosurgery, University Hospital Freiburg  
GitHub: [@kevinj24fr](https://github.com/kevinj24fr)  
Email: kevin.joseph@uniklinik-freiburg.de

