#!/usr/bin/env python3
"""
Build patristic distance matrix from a FASTA alignment using IQ-TREE2.

Steps:
- Run IQ-TREE2 with modelFinder and 1000 bootstrap replicate (-m MFP -bb 1000)
- Read the resulting ML tree (.treefile)
- Compute pairwise patristic distances from the tree
- Save a CSV distance matrix in the output directory

Requirements:
- IQ-TREE2 available in PATH (executable name `iqtree2` or `iqtree`)
- Python packages: biopython, pandas
  Install with: pip install biopython pandas

Usage example (PowerShell):
python .\build_patristic_matrix_iqtree.py --fasta aln_std_nt.fasta --output results --iqtree iqtree2

Author: GitHub Copilot
Date: 2025-08-15
"""

import os
import argparse
import subprocess
import shutil
from Bio import SeqIO, Phylo
import pandas as pd


def run_iqtree(fasta, iqtree_cmd, output_dir):
    """Run IQ-TREE2 on the input FASTA using model finder and 1000 bootstraps.
    Returns path to the produced treefile.
    """
    basename = os.path.splitext(os.path.basename(fasta))[0]
    cwd = os.path.abspath(output_dir)
    os.makedirs(cwd, exist_ok=True)

    # Build command
    # -s: input alignment
    # -m MFP: model finder
    # -bb 1000: ultrafast bootstrap 1000 replicates
    # -nt AUTO: use available cores
    # -pre: output prefix so IQ-TREE writes files with this prefix
    cmd = [iqtree_cmd, '-s', os.path.abspath(fasta), '-m', 'MFP', '-bb', '1000', '-nt', 'AUTO', '-pre', basename]

    print('Running IQ-TREE: ' + ' '.join(cmd))
    try:
        proc = subprocess.run(cmd, cwd=cwd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        print(proc.stdout)
    except subprocess.CalledProcessError as e:
        # Print stderr and log file contents to help debugging
        print('IQ-TREE failed. Stderr:\n', e.stderr)
        # Check for log file in both current working directory and output directory
        for logdir in [cwd, os.getcwd()]:
            logpath = os.path.join(logdir, basename + '.log')
            if os.path.exists(logpath):
                print(f'\n===== IQ-TREE LOG ({logpath}) =====')
                with open(logpath, 'r', encoding='utf-8', errors='ignore') as f:
                    print(f.read())
                break
        raise RuntimeError(f'IQ-TREE failed: {e}')

    # IQ-TREE writes files into the working directory. Look for output files in both cwd and current directory
    search_dirs = [cwd, os.getcwd()]
    possible_extensions = ['.treefile', '.contree', '.tree', '.bestTree', '.iqtree']
    treefile = None
    
    for search_dir in search_dirs:
        for ext in possible_extensions:
            candidate = os.path.join(search_dir, basename + ext)
            if os.path.exists(candidate):
                treefile = candidate
                break
        if treefile:
            break
    
    if treefile is None:
        # Try searching both directories for files that end with .treefile or .contree
        for search_dir in search_dirs:
            if os.path.exists(search_dir):
                for fname in os.listdir(search_dir):
                    if fname.endswith('.treefile') or fname.endswith('.contree'):
                        treefile = os.path.join(search_dir, fname)
                        break
                if treefile:
                    break
    
    if treefile is None:
        all_searched = []
        for search_dir in search_dirs:
            for ext in possible_extensions:
                all_searched.append(os.path.join(search_dir, basename + ext))
        raise FileNotFoundError(f'Expected treefile not found after IQ-TREE run. Searched: {all_searched}')

    print(f'IQ-TREE finished. Treefile: {treefile}')
    return treefile


def compute_patristic_matrix(treefile, fasta):
    """Compute pairwise patristic distances using Bio.Phylo from the given treefile.

    Returns a pandas DataFrame (index and columns are sequence IDs from FASTA).
    """
    # Read tree
    try:
        tree = Phylo.read(treefile, 'newick')
    except Exception:
        # Fallback: import reader from Bio.Phylo._io (static checkers may require this)
        try:
            from Bio.Phylo._io import read as _phylo_read
            tree = _phylo_read(treefile, 'newick')
        except Exception as e:
            raise RuntimeError(f'Unable to read treefile with Bio.Phylo: {e}')

    # Read sequence IDs in order from FASTA
    seq_ids = [record.id for record in SeqIO.parse(fasta, 'fasta')]

    # Map IDs to clade objects in the tree
    clade_map = {}
    for sid in seq_ids:
        # Bio.Phylo: find_clades can be used to locate leaf with given name
        clade = next(tree.find_clades(name=sid), None)
        if clade is None:
            # Try matching by stripping after first whitespace (IQ-TREE sometimes changes names)
            short = sid.split()[0]
            clade = next(tree.find_clades(name=short), None)
        if clade is None:
            print(f'Warning: sequence id "{sid}" not found in tree leaves')
        clade_map[sid] = clade

    # Compute distances
    n = len(seq_ids)
    mat = [[0.0]*n for _ in range(n)]

    for i, id1 in enumerate(seq_ids):
        c1 = clade_map.get(id1)
        for j, id2 in enumerate(seq_ids):
            if j < i:
                mat[i][j] = mat[j][i]
                continue
            c2 = clade_map.get(id2)
            if c1 is None or c2 is None:
                dist = float('nan')
            else:
                try:
                    dist = tree.distance(c1, c2)
                except Exception:
                    # Fallback: if distance fails, set NaN
                    dist = float('nan')
            mat[i][j] = dist
            mat[j][i] = dist

    df = pd.DataFrame(mat, index=seq_ids, columns=seq_ids)
    return df


def main():
    parser = argparse.ArgumentParser(description='Build patristic distance matrix using IQ-TREE2')
    parser.add_argument('--fasta', default='aln_std_nt.fasta', help='Aligned FASTA (nucleotide)')
    parser.add_argument('--output', default='.', help='Output directory')
    parser.add_argument('--iqtree', default='iqtree2', help='IQ-TREE2 executable name or path (default: iqtree2)')

    args = parser.parse_args()

    fasta = args.fasta
    outdir = args.output
    iqtree_cmd = args.iqtree

    # Check FASTA exists
    if not os.path.exists(fasta):
        raise FileNotFoundError(f'Fasta file not found: {fasta}')

    # Check IQ-TREE
    if shutil.which(iqtree_cmd) is None:
        alt = 'iqtree'
        if shutil.which(alt) is not None:
            print(f'IQ-TREE executable "{iqtree_cmd}" not found, using "{alt}"')
            iqtree_cmd = alt
        else:
            raise FileNotFoundError(f'IQ-TREE executable not found in PATH: {iqtree_cmd}. Install IQ-TREE2 and ensure it is in PATH.')

    # Run IQ-TREE
    treefile = run_iqtree(fasta, iqtree_cmd, outdir)

    # Compute patristic distances
    df = compute_patristic_matrix(treefile, fasta)

    # Save matrix
    out_csv = os.path.join(outdir, 'patristic_distance_matrix.csv')
    df.to_csv(out_csv)
    print(f'Patristic distance matrix saved to: {out_csv}')


if __name__ == '__main__':
    main()
