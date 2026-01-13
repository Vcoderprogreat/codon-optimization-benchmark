#!/usr/bin/env python3
import argparse
import time
from dnachisel import DnaOptimizationProblem, CodonOptimize, EnforceTranslation
from cai2 import CAI  # Updated to use the working cai2 library

HUMAN_WEIGHTS = {
    'TTT': 0.45, 'TTC': 1.0, 'TTA': 0.08, 'TTG': 0.13, 'CTT': 0.13, 'CTC': 0.2,
    'CTA': 0.07, 'CTG': 1.0, 'ATT': 0.36, 'ATC': 0.47, 'ATA': 0.17, 'ATG': 1.0,
    'GTT': 0.18, 'GTC': 0.24, 'GTA': 0.12, 'GTG': 1.0, 'TCT': 0.19, 'TCC': 0.22,
    'TCA': 0.15, 'TCG': 0.05, 'AGT': 0.15, 'AGC': 1.0, 'CCT': 0.28, 'CCC': 0.32,
    'CCA': 0.28, 'CCG': 0.11, 'ACT': 0.25, 'ACC': 0.36, 'ACA': 0.28, 'ACG': 0.11,
    'GCT': 0.27, 'GCC': 0.4, 'GCA': 0.23, 'GCG': 0.11, 'TAT': 0.44, 'TAC': 1.0,
    'CAT': 0.42, 'CAC': 1.0, 'CAA': 0.27, 'CAG': 1.0, 'AAT': 0.47, 'AAC': 1.0,
    'AAA': 0.43, 'AAG': 1.0, 'GAT': 0.46, 'GAC': 1.0, 'GAA': 0.42, 'GAG': 1.0,
    'TGT': 0.45, 'TGC': 1.0, 'TGG': 1.0, 'CGT': 0.08, 'CGC': 0.18, 'CGA': 0.11,
    'CGG': 0.2, 'AGA': 0.21, 'AGG': 0.2, 'GGT': 0.16, 'GGC': 0.34, 'GGA': 0.25,
    'GGG': 0.25, 'TAA': 1.0, 'TAG': 1.0, 'TGA': 1.0
}

def run_dnachisel(dna_seq):
    """
    Optimizes a DNA sequence using DNA Chisel for Human (H. sapiens) expression.
    """
    start_time = time.time()
    
    problem = DnaOptimizationProblem(
        sequence=dna_seq,
        constraints=[EnforceTranslation()],
        objectives=[CodonOptimize(species="h_sapiens")]
    )
    
    problem.optimize()
    optimized_seq = problem.sequence
    runtime = time.time() - start_time

    # Calculate CAI scores (Quality Metric)
    orig_cai = CAI(dna_seq, weights=HUMAN_WEIGHTS)
    opt_cai = CAI(optimized_seq, weights=HUMAN_WEIGHTS)

    print(f"Original: {dna_seq}")
    print(f"Optimized: {optimized_seq}")
    print(f"Original CAI: {orig_cai:.4f}")
    print(f"Optimized CAI: {opt_cai:.4f}")
    print(f"Runtime: {runtime:.6f}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Professional Codon Optimization Tool (DNA Chisel)")
    parser.add_argument("-s", "--sequence", type=str, required=True, help="Input DNA sequence")
    args = parser.parse_args()
    
    if len(args.sequence) % 3 != 0:
        print("Error: Sequence length must be a multiple of 3.")
    else:
        run_dnachisel(args.sequence.upper())
