#!/usr/bin/env python3
"""
Pairwise Distance Matrix Generator for Haplotype Sequences

This script generates multiple types of pairwise distance matrices:
1. Raw mutation counts (total mutations between sequences)
2. Gene-normalized raw mutations (mutations per gene length)
3. Hotspot mutation counts (mutations only at hotspot positions)
4. Gene-normalized hotspot mutations
5. Combined matrix (weighted combination of raw and hotspot metrics)

Author: GitHub Copilot
Date: 2025-07-11
"""

import pandas as pd
import numpy as np
from Bio import SeqIO
from itertools import combinations
import os
from typing import Dict, List, Tuple, Set, Optional
import argparse

class PairwiseMatrixGenerator:
    def __init__(self, fasta_file: str, gene_coords_file: str, 
                 hotspot_files: Dict[str, str], output_dir: str = "."):
        """
        Initialize the matrix generator.
        
        Args:
            fasta_file: Path to aligned amino acid FASTA file
            gene_coords_file: Path to gene coordinates CSV
            hotspot_files: Dictionary mapping serotype to hotspot file path
            output_dir: Directory to save output files
        """
        self.fasta_file = fasta_file
        self.gene_coords_file = gene_coords_file
        self.hotspot_files = hotspot_files
        self.output_dir = output_dir
        
        # Data containers
        self.sequences = {}
        self.gene_coords = {}
        self.hotspots = {}
        self.sequence_ids = []
        
        # Load all data
        self._load_sequences()
        self._load_gene_coordinates()
        self._load_hotspots()
        
    def _load_sequences(self):
        """Load aligned sequences from FASTA file."""
        print("Loading sequences...")
        temp_sequences = {}
        for record in SeqIO.parse(self.fasta_file, "fasta"):
            temp_sequences[record.id] = str(record.seq)
        
        # Check sequence lengths and handle mismatched lengths
        lengths = [len(seq) for seq in temp_sequences.values()]
        if len(set(lengths)) > 1:
            print(f"Warning: Found sequences of different lengths: {set(lengths)}")
            # Find the most common length (excluding potential reference sequences)
            from collections import Counter
            length_counts = Counter(lengths)
            most_common_length = length_counts.most_common(1)[0][0]
            
            # If there's a significantly longer sequence (likely reference), use the second most common
            if len(length_counts) > 1:
                sorted_lengths = length_counts.most_common()
                if sorted_lengths[0][0] > sorted_lengths[1][0] * 1.5:  # Reference is >1.5x longer
                    most_common_length = sorted_lengths[1][0]
                    print(f"Using length {most_common_length} (excluding reference sequence)")
            
            # Filter sequences to only those of the target length
            filtered_sequences = {}
            for seq_id, seq in temp_sequences.items():
                if len(seq) == most_common_length:
                    filtered_sequences[seq_id] = seq
                else:
                    print(f"Excluding sequence {seq_id} (length {len(seq)})")
            
            self.sequences = filtered_sequences
        else:
            self.sequences = temp_sequences
            
        self.sequence_ids = list(self.sequences.keys())
        print(f"Loaded {len(self.sequences)} sequences of length {len(next(iter(self.sequences.values())))}")
        
    def _load_gene_coordinates(self):
        """Load gene coordinates from CSV file."""
        print("Loading gene coordinates...")
        df = pd.read_csv(self.gene_coords_file, header=None, 
                        names=['serotype', 'gene', 'start', 'end'])
        
        for _, row in df.iterrows():
            serotype = row['serotype']
            if serotype not in self.gene_coords:
                self.gene_coords[serotype] = {}
            self.gene_coords[serotype][row['gene']] = (row['start'], row['end'])
        
        print(f"Loaded gene coordinates for {len(self.gene_coords)} serotypes")
        
    def _load_hotspots(self):
        """Load hotspot positions for each serotype."""
        print("Loading hotspot positions...")
        for serotype, file_path in self.hotspot_files.items():
            if os.path.exists(file_path):
                df = pd.read_csv(file_path, header=None, names=['position'])
                self.hotspots[serotype] = set(df['position'].values)
                print(f"Loaded {len(self.hotspots[serotype])} hotspots for {serotype}")
            else:
                print(f"Warning: Hotspot file not found for {serotype}: {file_path}")
                self.hotspots[serotype] = set()
    
    def _get_serotype_from_id(self, seq_id: str) -> str:
        """Extract serotype from sequence ID."""
        if '_D1_' in seq_id:
            return 'D1'
        elif '_D2_' in seq_id:
            return 'D2'
        elif '_D3_' in seq_id:
            return 'D3'
        elif '_D4_' in seq_id:
            return 'D4'
        else:
            # Default to D2 if can't determine (as reference seems to be D2)
            return 'D2'
    
    def _count_mutations(self, seq1: str, seq2: str, positions: Optional[Set[int]] = None) -> int:
        """
        Count mutations between two sequences.
        
        Args:
            seq1, seq2: Aligned sequences to compare
            positions: Set of positions to check (1-based). If None, check all positions.
            
        Returns:
            Number of mutations
        """
        if len(seq1) != len(seq2):
            raise ValueError("Sequences must be of equal length")
        
        mutations = 0
        for i, (aa1, aa2) in enumerate(zip(seq1, seq2)):
            pos = i + 1  # Convert to 1-based indexing
            
            # Skip gaps and ambiguous characters
            if aa1 in '-X' or aa2 in '-X':
                continue
                
            # If specific positions are provided, only check those
            if positions is not None and pos not in positions:
                continue
                
            if aa1 != aa2:
                mutations += 1
                
        return mutations
    
    def _count_gene_mutations(self, seq1: str, seq2: str, serotype: str, 
                            gene: str, hotspot_only: bool = False) -> Tuple[int, int]:
        """
        Count mutations in a specific gene region.
        
        Args:
            seq1, seq2: Aligned sequences
            serotype: Serotype to get gene coordinates
            gene: Gene name
            hotspot_only: If True, only count hotspot mutations
            
        Returns:
            Tuple of (mutation_count, gene_length)
        """
        if serotype not in self.gene_coords or gene not in self.gene_coords[serotype]:
            return 0, 0
            
        start, end = self.gene_coords[serotype][gene]
        gene_seq1 = seq1[start-1:end]  # Convert to 0-based indexing
        gene_seq2 = seq2[start-1:end]
        gene_length = end - start + 1
        
        if hotspot_only:
            # Get hotspot positions within this gene
            gene_hotspots = {pos - start + 1 for pos in self.hotspots.get(serotype, set()) 
                           if start <= pos <= end}
            mutations = self._count_mutations(gene_seq1, gene_seq2, gene_hotspots)
        else:
            mutations = self._count_mutations(gene_seq1, gene_seq2)
            
        return mutations, gene_length
    
    def generate_raw_mutation_matrix(self) -> pd.DataFrame:
        """Generate matrix of raw mutation counts."""
        print("Generating raw mutation matrix...")
        n_seqs = len(self.sequence_ids)
        matrix = np.zeros((n_seqs, n_seqs))
        
        for i, seq_id1 in enumerate(self.sequence_ids):
            for j, seq_id2 in enumerate(self.sequence_ids):
                if i <= j:  # Only calculate upper triangle
                    seq1 = self.sequences[seq_id1]
                    seq2 = self.sequences[seq_id2]
                    mutations = self._count_mutations(seq1, seq2)
                    matrix[i, j] = mutations
                    matrix[j, i] = mutations  # Symmetric matrix
                    
        return pd.DataFrame(matrix, index=self.sequence_ids, columns=self.sequence_ids)
    
    def generate_hotspot_mutation_matrix(self) -> pd.DataFrame:
        """Generate matrix of hotspot mutation counts."""
        print("Generating hotspot mutation matrix...")
        n_seqs = len(self.sequence_ids)
        matrix = np.zeros((n_seqs, n_seqs))
        
        for i, seq_id1 in enumerate(self.sequence_ids):
            for j, seq_id2 in enumerate(self.sequence_ids):
                if i <= j:  # Only calculate upper triangle
                    seq1 = self.sequences[seq_id1]
                    seq2 = self.sequences[seq_id2]
                    
                    # Determine serotype (use serotype of first sequence)
                    serotype = self._get_serotype_from_id(seq_id1)
                    hotspot_positions = self.hotspots.get(serotype, set())
                    
                    mutations = self._count_mutations(seq1, seq2, hotspot_positions)
                    matrix[i, j] = mutations
                    matrix[j, i] = mutations  # Symmetric matrix
                    
        return pd.DataFrame(matrix, index=self.sequence_ids, columns=self.sequence_ids)
    
    def generate_gene_normalized_matrices(self) -> Dict[str, pd.DataFrame]:
        """Generate gene-normalized mutation matrices."""
        print("Generating gene-normalized matrices...")
        
        # Get all unique genes across serotypes
        all_genes = set()
        for serotype_genes in self.gene_coords.values():
            all_genes.update(serotype_genes.keys())
        
        results = {}
        
        for mutation_type in ['raw', 'hotspot']:
            print(f"Processing {mutation_type} mutations by gene...")
            
            for gene in all_genes:
                print(f"  Processing gene: {gene}")
                n_seqs = len(self.sequence_ids)
                matrix = np.zeros((n_seqs, n_seqs))
                
                for i, seq_id1 in enumerate(self.sequence_ids):
                    for j, seq_id2 in enumerate(self.sequence_ids):
                        if i <= j:  # Only calculate upper triangle
                            seq1 = self.sequences[seq_id1]
                            seq2 = self.sequences[seq_id2]
                            
                            # Determine serotype
                            serotype = self._get_serotype_from_id(seq_id1)
                            
                            # Count mutations in this gene
                            hotspot_only = (mutation_type == 'hotspot')
                            mutations, gene_length = self._count_gene_mutations(
                                seq1, seq2, serotype, gene, hotspot_only)
                            
                            # Normalize by gene length
                            normalized_mutations = mutations / gene_length if gene_length > 0 else 0
                            
                            matrix[i, j] = normalized_mutations
                            matrix[j, i] = normalized_mutations  # Symmetric matrix
                
                df = pd.DataFrame(matrix, index=self.sequence_ids, columns=self.sequence_ids)
                results[f'{mutation_type}_{gene}_normalized'] = df
        
        return results
    
    def generate_total_normalized_matrices(self) -> Dict[str, pd.DataFrame]:
        """Generate matrices normalized by total genome length."""
        print("Generating total genome normalized matrices...")
        
        results = {}
        
        for mutation_type in ['raw', 'hotspot']:
            print(f"Processing total normalized {mutation_type} mutations...")
            n_seqs = len(self.sequence_ids)
            matrix = np.zeros((n_seqs, n_seqs))
            
            for i, seq_id1 in enumerate(self.sequence_ids):
                for j, seq_id2 in enumerate(self.sequence_ids):
                    if i <= j:  # Only calculate upper triangle
                        seq1 = self.sequences[seq_id1]
                        seq2 = self.sequences[seq_id2]
                        
                        # Determine serotype
                        serotype = self._get_serotype_from_id(seq_id1)
                        
                        if mutation_type == 'hotspot':
                            hotspot_positions = self.hotspots.get(serotype, set())
                            mutations = self._count_mutations(seq1, seq2, hotspot_positions)
                            # Normalize by number of hotspot positions
                            total_length = len(hotspot_positions) if hotspot_positions else 1
                        else:
                            mutations = self._count_mutations(seq1, seq2)
                            # Normalize by total sequence length (excluding gaps)
                            total_length = len([aa for aa in seq1 if aa not in '-X'])
                        
                        normalized_mutations = mutations / total_length if total_length > 0 else 0
                        
                        matrix[i, j] = normalized_mutations
                        matrix[j, i] = normalized_mutations  # Symmetric matrix
            
            df = pd.DataFrame(matrix, index=self.sequence_ids, columns=self.sequence_ids)
            results[f'{mutation_type}_total_normalized'] = df
        
        return results
    
    def generate_combined_matrix(self, raw_matrix: pd.DataFrame, 
                               hotspot_matrix: pd.DataFrame,
                               alpha: float = 0.5, beta: float = 0.5) -> pd.DataFrame:
        """
        Generate combined matrix using weighted sum of raw and hotspot matrices.
        
        Args:
            raw_matrix: Raw mutation matrix
            hotspot_matrix: Hotspot mutation matrix
            alpha: Weight for raw mutations
            beta: Weight for hotspot mutations
            
        Returns:
            Combined matrix
        """
        print(f"Generating combined matrix (α={alpha}, β={beta})...")
        
        # Normalize matrices to [0,1] range for fair combination
        raw_norm = (raw_matrix - raw_matrix.min().min()) / (raw_matrix.max().max() - raw_matrix.min().min())
        hotspot_norm = (hotspot_matrix - hotspot_matrix.min().min()) / (hotspot_matrix.max().max() - hotspot_matrix.min().min())
        
        combined = alpha * raw_norm + beta * hotspot_norm
        return combined
    
    def save_matrices(self, matrices: Dict[str, pd.DataFrame]):
        """Save all matrices to CSV files."""
        print("Saving matrices...")
        os.makedirs(self.output_dir, exist_ok=True)
        
        for name, matrix in matrices.items():
            filename = f"pairwise_matrix_{name}.csv"
            filepath = os.path.join(self.output_dir, filename)
            matrix.to_csv(filepath)
            print(f"Saved {name} matrix to {filepath}")
    
    def generate_all_matrices(self, alpha: float = 0.5, beta: float = 0.5) -> Dict[str, pd.DataFrame]:
        """Generate all types of matrices."""
        matrices = {}
        
        # Basic matrices
        matrices['raw_mutations'] = self.generate_raw_mutation_matrix()
        matrices['hotspot_mutations'] = self.generate_hotspot_mutation_matrix()
        
        # Total normalized matrices
        total_normalized = self.generate_total_normalized_matrices()
        matrices.update(total_normalized)
        
        # Gene-specific normalized matrices
        gene_normalized = self.generate_gene_normalized_matrices()
        matrices.update(gene_normalized)
        
        # Combined matrix
        matrices['combined'] = self.generate_combined_matrix(
            matrices['raw_total_normalized'], 
            matrices['hotspot_total_normalized'],
            alpha, beta
        )
        
        return matrices

def main():
    """Main function to run the matrix generation."""
    parser = argparse.ArgumentParser(description='Generate pairwise distance matrices')
    parser.add_argument('--fasta', default='aln_std_nt.fasta',
                       help='Aligned nucleotide FASTA file (output from filter_fasta_keep_ids.py)')
    parser.add_argument('--coords', default='Genome_Coordinates.csv',
                       help='Gene coordinates CSV file')
    parser.add_argument('--output', default='.',
                       help='Output directory')
    parser.add_argument('--alpha', type=float, default=0.5,
                       help='Weight for raw mutations in combined matrix')
    parser.add_argument('--beta', type=float, default=0.5,
                       help='Weight for hotspot mutations in combined matrix')
    
    args = parser.parse_args()
    
    # Define hotspot files
    hotspot_files = {
        'D1': 'DENV1.csv',
        'D2': 'DENV2.csv',
        'D3': 'DENV3.csv'
    }
    
    # Initialize generator
    generator = PairwiseMatrixGenerator(
        fasta_file=args.fasta,
        gene_coords_file=args.coords,
        hotspot_files=hotspot_files,
        output_dir=args.output
    )
    
    # Generate all matrices
    matrices = generator.generate_all_matrices(alpha=args.alpha, beta=args.beta)
    
    # Save matrices
    generator.save_matrices(matrices)
    
    # Print summary
    print("\nMatrix generation completed!")
    print(f"Generated {len(matrices)} matrices:")
    for name in matrices.keys():
        print(f"  - {name}")
    
    print(f"\nMatrix dimensions: {matrices['raw_mutations'].shape}")
    print(f"Output directory: {args.output}")

if __name__ == "__main__":
    main()
