#!/usr/bin/env python3
"""
Patristic Distance Calculator
This script calculates patristic distances from a phylogenetic tree file.
"""

import os
import sys
import subprocess
import shutil
from pathlib import Path
import argparse
from datetime import datetime

# Try to import optional dependencies
try:
    import pandas as pd
    import numpy as np
    PANDAS_AVAILABLE = True
except ImportError:
    PANDAS_AVAILABLE = False
    print("Warning: pandas/numpy not available. CSV output will be limited.")

try:
    from Bio import Phylo
    from Bio.Phylo.TreeConstruction import DistanceMatrix
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False
    print("Warning: BioPython not available. Patristic distance calculation will be limited.")

class PhylogeneticAnalyzer:
    def __init__(self, fasta_file, output_prefix="phylo_analysis"):
        self.fasta_file = fasta_file
        self.output_prefix = output_prefix
        self.work_dir = Path.cwd()
        
        # Check if input file exists
        if not Path(fasta_file).exists():
            raise FileNotFoundError(f"Input file '{fasta_file}' not found!")
        
        # Check if IQ-TREE2 is available
        if not self.check_iqtree():
            raise RuntimeError("IQ-TREE2 not found! Please install IQ-TREE2 and ensure it's in your PATH.")
    
    def check_iqtree(self):
        """Check if IQ-TREE2 is available in the system"""
        try:
            result = subprocess.run(['iqtree2', '--version'], 
                                  capture_output=True, text=True, timeout=30)
            if result.returncode == 0:
                print(f"✓ Found IQ-TREE2: {result.stdout.split()[0]} {result.stdout.split()[1]}")
                return True
        except (subprocess.TimeoutExpired, FileNotFoundError):
            pass
        
        # Try alternative command names
        for cmd in ['iqtree', 'iqtree-2']:
            try:
                result = subprocess.run([cmd, '--version'], 
                                      capture_output=True, text=True, timeout=30)
                if result.returncode == 0:
                    print(f"✓ Found IQ-TREE: {result.stdout.split()[0]} {result.stdout.split()[1]}")
                    self.iqtree_cmd = cmd
                    return True
            except (subprocess.TimeoutExpired, FileNotFoundError):
                continue
        
        return False
    
    def run_iqtree_ml(self):
        """Run IQ-TREE2 to build Maximum Likelihood tree"""
        print("=== Step 1: Building ML tree with IQ-TREE2 ===")
        
        iqtree_cmd = getattr(self, 'iqtree_cmd', 'iqtree2')
        output_prefix = f"{self.output_prefix}_ml"
        
        # IQ-TREE2 command for ML tree construction
        cmd = [
            iqtree_cmd,
            '-s', self.fasta_file,           # Input alignment
            '--prefix', output_prefix,        # Output prefix
            '-m', 'MFP',                     # Model finder plus (finds best model)
            '-bb', '1000',                   # Bootstrap replicates
            '-alrt', '1000',                 # SH-aLRT branch support
            '-T', 'AUTO',                    # Automatic number of threads
            '--keep-ident',                  # Keep identical sequences
            '-v'                             # Verbose output
        ]
        
        print(f"Running command: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=9600)
            
            if result.returncode != 0:
                print(f"Error running IQ-TREE2: {result.stderr}")
                return False
            
            print("✓ ML tree construction completed successfully!")
            print(f"✓ Tree file: {output_prefix}.treefile")
            print(f"✓ Log file: {output_prefix}.log")
            print(f"✓ Model info: {output_prefix}.iqtree")
            
            return True
            
        except subprocess.TimeoutExpired:
            print("Error: IQ-TREE2 analysis timed out (>1 hour)")
            return False
        except Exception as e:
            print(f"Error running IQ-TREE2: {e}")
            return False
    
    def calculate_patristic_distances(self):
        """Calculate patristic distances from the ML tree"""
        print("=== Step 2: Calculating patristic distance matrix ===")
        
        tree_file = f"{self.output_prefix}_ml.treefile"
        
        if not Path(tree_file).exists():
            print(f"Error: Tree file '{tree_file}' not found!")
            return False
        
        if not BIOPYTHON_AVAILABLE:
            print("Warning: BioPython not available. Skipping patristic distance calculation.")
            print("To enable this feature, install BioPython: pip install biopython")
            return False
        
        try:
            # Read the tree using Bio.Phylo
            tree = Phylo.read(tree_file, 'newick')
            
            # Get all terminal nodes (leaf names)
            terminals = [term.name for term in tree.get_terminals()]
            n_taxa = len(terminals)
            
            print(f"Calculating patristic distances for {n_taxa} taxa...")
            
            # Calculate distance matrix manually since numpy might not be available
            distance_data = []
            
            # Create header row
            header = ['taxon'] + terminals
            distance_data.append(header)
            
            for i, taxon1 in enumerate(terminals):
                row = [taxon1]
                for j, taxon2 in enumerate(terminals):
                    if i == j:
                        row.append('0.000000')
                    else:
                        dist = tree.distance(taxon1, taxon2)
                        row.append(f"{dist:.6f}")
                distance_data.append(row)
            
            # Save as CSV manually
            patristic_file = f"{self.output_prefix}_patristic_distances.csv"
            with open(patristic_file, 'w') as f:
                for row in distance_data:
                    f.write(','.join(row) + '\n')
            
            print(f"✓ Patristic distance matrix saved: {patristic_file}")
            
            # Also save in phylip format for compatibility
            phylip_file = f"{self.output_prefix}_patristic_distances.phylip"
            with open(phylip_file, 'w') as f:
                f.write(f"    {n_taxa}\n")
                for i, taxon in enumerate(terminals):
                    # Truncate taxon name to 10 characters for phylip format
                    taxon_short = taxon[:10].ljust(10)
                    distances = []
                    for j, taxon2 in enumerate(terminals):
                        if i == j:
                            distances.append("0.000000")
                        else:
                            dist = tree.distance(taxon, taxon2)
                            distances.append(f"{dist:.6f}")
                    f.write(f"{taxon_short}  {'  '.join(distances)}\n")
            
            print(f"✓ Patristic distance matrix (PHYLIP format): {phylip_file}")
            return True
            
        except Exception as e:
            print(f"Error calculating patristic distances: {e}")
            return False
    
    def generate_summary_report(self):
        """Generate a summary report of the analysis"""
        print("=== Generating Summary Report ===")
        
        report_file = f"{self.output_prefix}_analysis_report.txt"
        
        with open(report_file, 'w') as f:
            f.write("=" * 60 + "\n")
            f.write("PHYLOGENETIC ANALYSIS REPORT\n")
            f.write("=" * 60 + "\n")
            f.write(f"Analysis date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Input file: {self.fasta_file}\n")
            f.write(f"Output prefix: {self.output_prefix}\n\n")
            
            f.write("ANALYSIS STEPS PERFORMED:\n")
            f.write("-" * 30 + "\n")
            f.write("1. Maximum Likelihood tree construction (IQ-TREE2)\n")
            f.write("   - Model selection: MFP (ModelFinder Plus)\n")
            f.write("   - Bootstrap support: 1000 replicates\n")
            f.write("   - SH-aLRT support: 1000 replicates\n\n")
            
            f.write("2. Patristic distance matrix calculation\n")
            f.write("   - Based on ML tree branch lengths\n")
            f.write("   - Output formats: CSV and PHYLIP\n\n")
            
            f.write("OUTPUT FILES:\n")
            f.write("-" * 15 + "\n")
            
            # List all output files
            output_files = [
                f"{self.output_prefix}_ml.treefile",
                f"{self.output_prefix}_ml.iqtree",
                f"{self.output_prefix}_ml.log",
                f"{self.output_prefix}_patristic_distances.csv",
                f"{self.output_prefix}_patristic_distances.phylip"
            ]
            
            for file in output_files:
                if Path(file).exists():
                    file_size = Path(file).stat().st_size
                    f.write(f"✓ {file} ({file_size:,} bytes)\n")
                else:
                    f.write(f"✗ {file} (not found)\n")
            
            f.write(f"\nReport file: {report_file}\n")
        
        print(f"✓ Summary report saved: {report_file}")
    
    def run_full_analysis(self):
        """Run the complete phylogenetic analysis pipeline"""
        print("=" * 60)
        print("PHYLOGENETIC ANALYSIS PIPELINE")
        print("=" * 60)
        print(f"Input file: {self.fasta_file}")
        print(f"Output prefix: {self.output_prefix}")
        print()
        
        success_steps = 0
        total_steps = 2
        
        # Step 1: Build ML tree
        if self.run_iqtree_ml():
            success_steps += 1
        else:
            print("✗ Failed to build ML tree. Stopping analysis.")
            return False
        
        print()
        
        # Step 2: Calculate patristic distances
        if self.calculate_patristic_distances():
            success_steps += 1
        else:
            print("✗ Failed to calculate patristic distances.")
        
        print()
        
        # Generate summary report
        self.generate_summary_report()
        
        print("=" * 60)
        print(f"ANALYSIS COMPLETED: {success_steps}/{total_steps} steps successful")
        print("=" * 60)
        
        if success_steps == total_steps:
            print("✓ All analysis steps completed successfully!")
            return True
        else:
            print(f"⚠ Analysis completed with {total_steps - success_steps} failed step(s)")
            return False

def main():
    parser = argparse.ArgumentParser(
        description="Phylogenetic analysis using IQ-TREE2",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python phylo_analysis.py aligned1075.fasta
  python phylo_analysis.py aligned1075.fasta --prefix my_analysis
        """
    )
    
    parser.add_argument('fasta_file', 
                       help='Input aligned FASTA file')
    parser.add_argument('--prefix', '-p', 
                       default='phylo_analysis',
                       help='Output prefix for all files (default: phylo_analysis)')
    
    args = parser.parse_args()
    
    try:
        analyzer = PhylogeneticAnalyzer(args.fasta_file, args.prefix)
        analyzer.run_full_analysis()
        
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
