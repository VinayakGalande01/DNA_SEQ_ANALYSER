# 🧬 DNA Sequence Analyzer

A comprehensive DNA sequence analysis tool with both command-line and web interface options. Analyze DNA sequences for GC content, reverse complement, protein translation, and motif search.

## 🚀 Features

### Core Analysis
- **GC Content Calculation**: Percentage of G and C bases in sequences
- **Reverse Complement**: Generate reverse complement of DNA sequences
- **Protein Translation**: Translate DNA to protein (stops at first stop codon)
- **Motif Search**: Find all occurrences of specified motifs (case-insensitive)
- **Sequence Validation**: Check for valid DNA bases

### Command Line Interface
- **Flexible Input**: Accept any FASTA file
- **Customizable Output**: Specify output location and format
- **CLI Options**: Command-line arguments for motif and output file
- **Error Handling**: Robust error handling for malformed sequences

### Web Interface (Streamlit)
- **Interactive Upload**: Drag-and-drop FASTA file upload
- **Real-time Analysis**: Instant results with progress indicators
- **Visualizations**: Interactive charts for GC content, sequence lengths, and motif counts
- **Customizable Display**: Toggle sequence details, reverse complements, and protein translations
- **Download Results**: Export analysis results as CSV

## 📦 Installation

1. **Clone or download the project files**

2. **Install dependencies**:
   ```bash
   pip install -r requirements.txt
   ```

## 🖥️ Command Line Usage

### Basic Usage
```bash
python dna_sequence_analyzer.py example.fasta
```

### Advanced Usage
```bash
# Specify custom motif
python dna_sequence_analyzer.py sequences.fasta --motif TATA

# Custom output location
python dna_sequence_analyzer.py data.fasta --output results/analysis.csv

# Skip console preview
python dna_sequence_analyzer.py example.fasta --no-preview

# Full example
python dna_sequence_analyzer.py sequences.fasta --motif ATG --output results/output.csv
```

### Command Line Options
- `fasta_file`: Input FASTA file (required)
- `--motif`: Motif to search for (default: ATG)
- `--output`: Output CSV file path (default: results/output.csv)
- `--no-preview`: Skip printing results to console

## 🌐 Web Interface Usage

### Start the Streamlit App
```bash
streamlit run streamlit_app.py
```

### Using the Web Interface
1. **Upload File**: Use the sidebar to upload a FASTA file
2. **Configure Settings**: Set the motif to search for
3. **View Results**: See interactive charts and detailed analysis
4. **Download**: Export results as CSV

## 📁 File Structure

```
DNA_ANALYZER/
├── dna_sequence_analyzer.py  # Command-line tool
├── streamlit_app.py          # Web interface
├── example.fasta             # Sample DNA sequences
├── requirements.txt          # Python dependencies
└── README.md                 # This file
```

## 📊 Output Format

The analysis generates a CSV file with the following columns:

| Column | Description |
|--------|-------------|
| `Sequence_ID` | Identifier from FASTA file |
| `GC_Content` | GC content percentage |
| `Reverse_Complement` | Reverse complement sequence |
| `Protein_Translation` | Translated protein sequence |
| `Motif_Positions` | Comma-separated positions of motif |

## 🧪 Example FASTA Format

```
>sequence1
ATGCGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGC
>sequence2
ATGCGCGGCTTAACTGACTGACTGACTGACGATCGAT
>sequence3
ATGAAATTTGGGCCCTTTAAAGGGCCCTTTAAAGGGC
```

## 🔧 Technical Details

### Dependencies
- **Biopython**: FASTA parsing and sequence manipulation
- **Streamlit**: Web interface framework
- **Pandas**: Data manipulation and CSV export
- **Plotly**: Interactive visualizations
- **NumPy**: Numerical operations

### Analysis Methods
- **GC Content**: Uses Biopython's GC() function
- **Reverse Complement**: Biopython's reverse_complement() method
- **Translation**: Biopython's translate() with to_stop=True
- **Motif Search**: Case-insensitive string search with 1-based indexing

## 🐛 Error Handling

The tool includes comprehensive error handling for:
- Missing or malformed FASTA files
- Invalid DNA sequences
- Translation errors
- File I/O issues

## 📈 Visualizations (Web Interface)

- **GC Content Chart**: Bar chart showing GC content for each sequence
- **Sequence Length Chart**: Comparison of sequence lengths
- **Motif Count Chart**: Frequency of motif occurrences
- **Interactive Features**: Hover tooltips, zoom, and pan

## 🤝 Contributing

Feel free to submit issues, feature requests, or pull requests to improve the tool.

## 📄 License

This project is open source and available under the MIT License.

---
