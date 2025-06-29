#!/usr/bin/env python3
"""
DNA Sequence Analyzer

A Python script that analyzes DNA sequences from FASTA files using Biopython.
Performs GC content analysis, reverse complement, protein translation, and motif search.
"""

import argparse
import csv
import os
import sys
from typing import List, Tuple, Optional
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction


def calculate_gc_content(sequence: str) -> float:
    """
    Calculate GC content percentage of a DNA sequence.
    
    Args:
        sequence (str): DNA sequence
        
    Returns:
        float: GC content percentage
    """
    return round(gc_fraction(sequence) * 100, 2)


def get_reverse_complement(sequence: str) -> str:
    """
    Get the reverse complement of a DNA sequence.
    
    Args:
        sequence (str): DNA sequence
        
    Returns:
        str: Reverse complement sequence
    """
    seq_obj = Seq(sequence)
    return str(seq_obj.reverse_complement())


def translate_protein(sequence: str) -> str:
    """
    Translate DNA sequence to protein, stopping at first stop codon.
    
    Args:
        sequence (str): DNA sequence
        
    Returns:
        str: Translated protein sequence
    """
    seq_obj = Seq(sequence)
    protein = seq_obj.translate(to_stop=True)
    return str(protein)


def find_motif_positions(sequence: str, motif: str) -> List[int]:
    """
    Find all positions of a motif in a DNA sequence (case-insensitive).
    
    Args:
        sequence (str): DNA sequence
        motif (str): Motif to search for
        
    Returns:
        List[int]: 1-based positions where motif occurs
    """
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


def read_fasta_file(filename: str) -> List[Tuple[str, str]]:
    """
    Read sequences from a FASTA file using Biopython.
    
    Args:
        filename (str): Path to FASTA file
        
    Returns:
        List[Tuple[str, str]]: List of (sequence_id, sequence) tuples
        
    Raises:
        FileNotFoundError: If file doesn't exist
        ValueError: If file is malformed
    """
    if not os.path.exists(filename):
        raise FileNotFoundError(f"FASTA file '{filename}' not found")
    
    sequences = []
    try:
        for record in SeqIO.parse(filename, "fasta"):
            sequences.append((record.id, str(record.seq)))
    except Exception as e:
        raise ValueError(f"Error parsing FASTA file: {e}")
    
    if not sequences:
        raise ValueError("No sequences found in FASTA file")
    
    return sequences


def analyze_sequences(sequences: List[Tuple[str, str]], motif: str = "ATG") -> List[dict]:
    """
    Analyze a list of DNA sequences.
    
    Args:
        sequences (List[Tuple[str, str]]): List of (sequence_id, sequence) tuples
        motif (str): Motif to search for
        
    Returns:
        List[dict]: List of analysis results
    """
    results = []
    
    for seq_id, sequence in sequences:
        try:
            # Validate sequence contains only valid DNA bases
            valid_bases = set('ATCGN')
            if not all(base in valid_bases for base in sequence.upper()):
                print(f"Warning: Sequence {seq_id} contains invalid DNA bases")
            
            # Perform analysis
            gc_content = calculate_gc_content(sequence)
            reverse_complement = get_reverse_complement(sequence)
            protein_translation = translate_protein(sequence)
            motif_positions = find_motif_positions(sequence, motif)
            
            results.append({
                'Sequence_ID': seq_id,
                'GC_Content': gc_content,
                'Reverse_Complement': reverse_complement,
                'Protein_Translation': protein_translation,
                'Motif_Positions': motif_positions
            })
            
        except Exception as e:
            print(f"Error analyzing sequence {seq_id}: {e}")
            # Add error entry
            results.append({
                'Sequence_ID': seq_id,
                'GC_Content': 'ERROR',
                'Reverse_Complement': 'ERROR',
                'Protein_Translation': 'ERROR',
                'Motif_Positions': 'ERROR'
            })
    
    return results


def save_results_to_csv(results: List[dict], output_file: str):
    """
    Save analysis results to CSV file.
    
    Args:
        results (List[dict]): Analysis results
        output_file (str): Output CSV file path
    """
    # Create results directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    with open(output_file, 'w', newline='', encoding='utf-8') as csvfile:
        fieldnames = ['Sequence_ID', 'GC_Content', 'Reverse_Complement', 
                     'Protein_Translation', 'Motif_Positions']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        
        writer.writeheader()
        for result in results:
            # Convert motif positions list to string
            result_copy = result.copy()
            if isinstance(result_copy['Motif_Positions'], list):
                result_copy['Motif_Positions'] = ','.join(map(str, result_copy['Motif_Positions']))
            writer.writerow(result_copy)


def print_results_table(results: List[dict]):
    """
    Print analysis results as a formatted table to console.
    
    Args:
        results (List[dict]): Analysis results
    """
    print("\n" + "="*100)
    print("DNA SEQUENCE ANALYSIS RESULTS")
    print("="*100)
    
    for result in results:
        print(f"\nSequence ID: {result['Sequence_ID']}")
        print(f"GC Content: {result['GC_Content']}%")
        print(f"Reverse Complement: {result['Reverse_Complement']}")
        print(f"Protein Translation: {result['Protein_Translation']}")
        
        if isinstance(result['Motif_Positions'], list):
            positions_str = ', '.join(map(str, result['Motif_Positions'])) if result['Motif_Positions'] else 'None found'
        else:
            positions_str = result['Motif_Positions']
        print(f"Motif Positions: {positions_str}")
        print("-" * 50)


def main():
    """Main function to run the DNA sequence analyzer."""
    parser = argparse.ArgumentParser(
        description="Analyze DNA sequences from a FASTA file",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python dna_sequence_analyzer.py example.fasta
  python dna_sequence_analyzer.py sequences.fasta --motif ATG
  python dna_sequence_analyzer.py data.fasta --motif TATA --output results/analysis.csv
        """
    )
    
    parser.add_argument(
        'fasta_file',
        help='Input FASTA file containing DNA sequences'
    )
    
    parser.add_argument(
        '--motif',
        default='ATG',
        help='Motif to search for (default: ATG)'
    )
    
    parser.add_argument(
        '--output',
        default='results/output.csv',
        help='Output CSV file (default: results/output.csv)'
    )
    
    parser.add_argument(
        '--no-preview',
        action='store_true',
        help='Skip printing results to console'
    )
    
    args = parser.parse_args()
    
    try:
        # Read sequences from FASTA file
        print(f"Reading sequences from {args.fasta_file}...")
        sequences = read_fasta_file(args.fasta_file)
        print(f"Found {len(sequences)} sequences")
        
        # Analyze sequences
        print(f"Analyzing sequences with motif '{args.motif}'...")
        results = analyze_sequences(sequences, args.motif)
        
        # Save results to CSV
        print(f"Saving results to {args.output}...")
        save_results_to_csv(results, args.output)
        print(f"Results saved successfully!")
        
        # Print preview to console
        if not args.no_preview:
            print_results_table(results)
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main() 