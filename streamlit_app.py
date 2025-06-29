import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from io import StringIO
import base64
from typing import List, Tuple
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
import tempfile
import os

# Page configuration
st.set_page_config(
    page_title="DNA Sequence Analyzer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS for better styling
st.markdown("""
<style>
    .main-header {
        font-size: 3rem;
        color: #1f77b4;
        text-align: center;
        margin-bottom: 2rem;
    }
    .metric-card {
        background-color: #f0f2f6;
        padding: 1rem;
        border-radius: 0.5rem;
        border-left: 4px solid #1f77b4;
    }
    .sequence-box {
        background-color: #f8f9fa;
        padding: 1rem;
        border-radius: 0.5rem;
        font-family: 'Courier New', monospace;
        font-size: 0.9rem;
        overflow-x: auto;
    }
</style>
""", unsafe_allow_html=True)

def calculate_gc_content(sequence: str) -> float:
    """Calculate GC content percentage of a DNA sequence."""
    return round(gc_fraction(sequence) * 100, 2)

def get_reverse_complement(sequence: str) -> str:
    """Get the reverse complement of a DNA sequence."""
    seq_obj = Seq(sequence)
    return str(seq_obj.reverse_complement())

def translate_protein(sequence: str) -> str:
    """Translate DNA sequence to protein, stopping at first stop codon."""
    seq_obj = Seq(sequence)
    protein = seq_obj.translate(to_stop=True)
    return str(protein)

def find_motif_positions(sequence: str, motif: str) -> List[int]:
    """Find all positions of a motif in a DNA sequence (case-insensitive)."""
    sequence_upper = sequence.upper()
    motif_upper = motif.upper()
    positions = []
    start = 0
    
    while True:
        pos = sequence_upper.find(motif_upper, start)
        if pos == -1:
            break
        positions.append(pos + 1)  # 1-based indexing
        start = pos + 1
    
    return positions

def analyze_sequences(sequences: List[Tuple[str, str]], motif: str = "ATG") -> List[dict]:
    """Analyze a list of DNA sequences."""
    results = []
    
    for seq_id, sequence in sequences:
        try:
            # Validate sequence contains only valid DNA bases
            valid_bases = set('ATCGN')
            if not all(base in valid_bases for base in sequence.upper()):
                st.warning(f"Sequence {seq_id} contains invalid DNA bases")
            
            # Perform analysis
            gc_content = calculate_gc_content(sequence)
            reverse_complement = get_reverse_complement(sequence)
            protein_translation = translate_protein(sequence)
            motif_positions = find_motif_positions(sequence, motif)
            
            results.append({
                'Sequence_ID': seq_id,
                'Length': len(sequence),
                'GC_Content': gc_content,
                'Reverse_Complement': reverse_complement,
                'Protein_Translation': protein_translation,
                'Motif_Positions': motif_positions,
                'Motif_Count': len(motif_positions),
                'Original_Sequence': sequence
            })
            
        except Exception as e:
            st.error(f"Error analyzing sequence {seq_id}: {e}")
            results.append({
                'Sequence_ID': seq_id,
                'Length': 0,
                'GC_Content': 'ERROR',
                'Reverse_Complement': 'ERROR',
                'Protein_Translation': 'ERROR',
                'Motif_Positions': 'ERROR',
                'Motif_Count': 0,
                'Original_Sequence': sequence
            })
    
    return results

def create_gc_content_chart(results: List[dict]):
    """Create a bar chart of GC content for all sequences."""
    df = pd.DataFrame(results)
    df_filtered = df[df['GC_Content'] != 'ERROR']
    
    if not df_filtered.empty:
        fig = px.bar(
            df_filtered, 
            x='Sequence_ID', 
            y='GC_Content',
            title='GC Content by Sequence',
            labels={'GC_Content': 'GC Content (%)', 'Sequence_ID': 'Sequence ID'},
            color='GC_Content',
            color_continuous_scale='viridis'
        )
        fig.update_layout(height=400)
        return fig
    return None

def create_sequence_length_chart(results: List[dict]):
    """Create a bar chart of sequence lengths."""
    df = pd.DataFrame(results)
    df_filtered = df[df['Length'] > 0]
    
    if not df_filtered.empty:
        fig = px.bar(
            df_filtered,
            x='Sequence_ID',
            y='Length',
            title='Sequence Lengths',
            labels={'Length': 'Length (bp)', 'Sequence_ID': 'Sequence ID'},
            color='Length',
            color_continuous_scale='plasma'
        )
        fig.update_layout(height=400)
        return fig
    return None

def create_motif_count_chart(results: List[dict], motif: str):
    """Create a bar chart of motif counts."""
    df = pd.DataFrame(results)
    df_filtered = df[df['Motif_Count'] >= 0]
    
    if not df_filtered.empty:
        fig = px.bar(
            df_filtered,
            x='Sequence_ID',
            y='Motif_Count',
            title=f'Count of "{motif}" Motif by Sequence',
            labels={'Motif_Count': f'Count of {motif}', 'Sequence_ID': 'Sequence ID'},
            color='Motif_Count',
            color_continuous_scale='inferno'
        )
        fig.update_layout(height=400)
        return fig
    return None

def download_csv(results: List[dict]) -> str:
    """Generate CSV download link."""
    df = pd.DataFrame(results)
    # Remove the original sequence column for CSV download
    df_clean = df.drop('Original_Sequence', axis=1, errors='ignore')
    
    csv = df_clean.to_csv(index=False)
    b64 = base64.b64encode(csv.encode()).decode()
    href = f'<a href="data:file/csv;base64,{b64}" download="dna_analysis_results.csv">Download CSV</a>'
    return href

def main():
    # Header
    st.markdown('<h1 class="main-header">🧬 DNA Sequence Analyzer</h1>', unsafe_allow_html=True)
    
    # Sidebar
    st.sidebar.title("Settings")
    
    # File upload (allow multiple files)
    st.sidebar.header("📁 Upload FASTA File")
    uploaded_files = st.sidebar.file_uploader(
        "Choose one or more FASTA files",
        type=['fasta', 'fa', 'txt'],
        accept_multiple_files=True,
        help="Upload one or more FASTA files containing DNA sequences"
    )
    
    # Motif input
    st.sidebar.header("🔍 Motif Search")
    motif = st.sidebar.text_input(
        "Motif to search for:",
        value="ATG",
        help="Enter the DNA motif you want to search for (case-insensitive)"
    )
    
    # Analysis options
    st.sidebar.header("⚙️ Analysis Options")
    show_sequences = st.sidebar.checkbox("Show original sequences", value=False)
    show_reverse_complement = st.sidebar.checkbox("Show reverse complements", value=False)
    show_proteins = st.sidebar.checkbox("Show protein translations", value=True)
    
    # Main content
    if uploaded_files:
        sequences = []
        for uploaded_file in uploaded_files:
            try:
                content = uploaded_file.read().decode('utf-8')
                fasta_io = StringIO(content)
                for record in SeqIO.parse(fasta_io, "fasta"):
                    sequences.append((record.id, str(record.seq), uploaded_file.name))
            except Exception as e:
                st.error(f"Error processing file {uploaded_file.name}: {e}")
        if sequences:
            st.success(f"✅ Successfully loaded {len(sequences)} sequences from {len(uploaded_files)} file(s)")
            # Analyze sequences, passing file info
            results = []
            for seq_id, sequence, file_name in sequences:
                try:
                    valid_bases = set('ATCGN')
                    if not all(base in valid_bases for base in sequence.upper()):
                        st.warning(f"Sequence {seq_id} in {file_name} contains invalid DNA bases")
                    gc_content = calculate_gc_content(sequence)
                    reverse_complement = get_reverse_complement(sequence)
                    protein_translation = translate_protein(sequence)
                    motif_positions = find_motif_positions(sequence, motif)
                    results.append({
                        'File': file_name,
                        'Sequence_ID': seq_id,
                        'Length': len(sequence),
                        'GC_Content': gc_content,
                        'Reverse_Complement': reverse_complement,
                        'Protein_Translation': protein_translation,
                        'Motif_Positions': motif_positions,
                        'Motif_Count': len(motif_positions),
                        'Original_Sequence': sequence
                    })
                except Exception as e:
                    st.error(f"Error analyzing sequence {seq_id} in {file_name}: {e}")
                    results.append({
                        'File': file_name,
                        'Sequence_ID': seq_id,
                        'Length': 0,
                        'GC_Content': 'ERROR',
                        'Reverse_Complement': 'ERROR',
                        'Protein_Translation': 'ERROR',
                        'Motif_Positions': 'ERROR',
                        'Motif_Count': 0,
                        'Original_Sequence': sequence
                    })
            # The rest of the app (charts, tables, etc.) should use the 'File' column for grouping/comparison
            # Update all DataFrame and chart code to use 'File' as a color/grouping variable
            
            # Display summary metrics
            st.header("📊 Analysis Summary")
            col1, col2, col3, col4 = st.columns(4)
            
            with col1:
                st.metric("Total Sequences", len(results))
            
            with col2:
                avg_gc = sum(r['GC_Content'] for r in results if r['GC_Content'] != 'ERROR') / len([r for r in results if r['GC_Content'] != 'ERROR'])
                st.metric("Average GC Content", f"{avg_gc:.1f}%")
            
            with col3:
                total_motifs = sum(r['Motif_Count'] for r in results)
                st.metric(f"Total '{motif}' Motifs", total_motifs)
            
            with col4:
                avg_length = sum(r['Length'] for r in results if r['Length'] > 0) / len([r for r in results if r['Length'] > 0])
                st.metric("Average Length", f"{avg_length:.0f} bp")
            
            # Charts
            st.header("📈 Visualizations")
            col1, col2 = st.columns(2)
            
            with col1:
                gc_chart = create_gc_content_chart(results)
                if gc_chart:
                    st.plotly_chart(gc_chart, use_container_width=True)
            
            with col2:
                length_chart = create_sequence_length_chart(results)
                if length_chart:
                    st.plotly_chart(length_chart, use_container_width=True)
            
            # Motif count chart
            motif_chart = create_motif_count_chart(results, motif)
            if motif_chart:
                st.plotly_chart(motif_chart, use_container_width=True)
            
            # Detailed results
            st.header("🔬 Detailed Analysis")
            
            # Create DataFrame for display
            df = pd.DataFrame(results)
            df_display = df.copy()
            
            # Format motif positions for display
            df_display['Motif_Positions'] = df_display['Motif_Positions'].apply(
                lambda x: ', '.join(map(str, x)) if isinstance(x, list) and x else 'None found'
            )
            
            # Remove original sequence from display unless requested
            if not show_sequences:
                df_display = df_display.drop('Original_Sequence', axis=1, errors='ignore')
            
            # Remove reverse complement unless requested
            if not show_reverse_complement:
                df_display = df_display.drop('Reverse_Complement', axis=1, errors='ignore')
            
            # Remove protein translation unless requested
            if not show_proteins:
                df_display = df_display.drop('Protein_Translation', axis=1, errors='ignore')
            
            # Display table
            st.dataframe(df_display, use_container_width=True)
            
            # Download section
            st.header("💾 Download Results")
            st.markdown(download_csv(results), unsafe_allow_html=True)
            
            # Individual sequence details
            if show_sequences or show_reverse_complement or show_proteins:
                st.header("🧬 Sequence Details")
                
                for result in results:
                    with st.expander(f"Sequence: {result['Sequence_ID']} ({result['Length']} bp)"):
                        col1, col2 = st.columns(2)
                        
                        with col1:
                            st.subheader("Basic Information")
                            st.metric("GC Content", f"{result['GC_Content']}%")
                            st.metric(f"'{motif}' Count", result['Motif_Count'])
                            
                            if result['Motif_Positions'] and result['Motif_Positions'] != 'ERROR':
                                st.write(f"**'{motif}' Positions:** {', '.join(map(str, result['Motif_Positions']))}")
                        
                        with col2:
                            if show_sequences:
                                st.subheader("Original Sequence")
                                st.code(result['Original_Sequence'], language=None)
                            
                            if show_reverse_complement and result['Reverse_Complement'] != 'ERROR':
                                st.subheader("Reverse Complement")
                                st.code(result['Reverse_Complement'], language=None)
                            
                            if show_proteins and result['Protein_Translation'] != 'ERROR':
                                st.subheader("Protein Translation")
                                st.code(result['Protein_Translation'], language=None)
            
    else:
        # Welcome message and instructions
        st.info("👆 Please upload a FASTA file using the sidebar to begin analysis.")
        
        # Example section
        st.header("📋 Example FASTA Format")
        st.code(""">sequence1
ATGCGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGC
>sequence2
ATGCGCGGCTTAACTGACTGACTGACTGACGATCGAT
>sequence3
ATGAAATTTGGGCCCTTTAAAGGGCCCTTTAAAGGGC""", language=None)
        
        # Features section
        st.header("✨ Features")
        col1, col2 = st.columns(2)
        
        with col1:
            st.markdown("""
            **🔬 Analysis Features:**
            - GC content calculation
            - Reverse complement generation
            - Protein translation (stops at first stop codon)
            - Motif search (case-insensitive)
            - Sequence length analysis
            """)
        
        with col2:
            st.markdown("""
            **📊 Visualization Features:**
            - Interactive charts with Plotly
            - GC content distribution
            - Sequence length comparison
            - Motif frequency analysis
            - Downloadable results
            """)

if __name__ == "__main__":
    main() 