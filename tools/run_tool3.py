#!/usr/bin/env python3
import argparse
import time
from cai2 import CAI

# Top codons for Humans (The 'Greedy' approach)
GREEDY_MAP = {
    'A': 'GCC', 'C': 'TGC', 'D': 'GAC', 'E': 'GAG', 'F': 'TTC', 'G': 'GGC',
    'H': 'CAC', 'I': 'ATC', 'L': 'CTG', 'K': 'AAG', 'M': 'ATG', 'N': 'AAC',
    'P': 'CCC', 'Q': 'CAG', 'R': 'CGC', 'S': 'AGC', 'T': 'ACC', 'V': 'GTG',
    'W': 'TGG', 'Y': 'TAC', '*': 'TGA'
}

# Genetic Code for translation
TABLE = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*', 'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
}

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

def run_tool_3(dna_seq):
    start_time = time.time()
    
    # Translate DNA to Protein, then rebuild using the #1 human codons
    aas = [TABLE[dna_seq[i:i+3]] for i in range(0, len(dna_seq), 3)]
    optimized_seq = "".join([GREEDY_MAP[aa] for aa in aas])
    
    runtime = time.time() - start_time
    orig_cai = CAI(dna_seq, weights=HUMAN_WEIGHTS)
    opt_cai = CAI(optimized_seq, weights=HUMAN_WEIGHTS)

    print(f"Original: {dna_seq}")
    print(f"Optimized: {optimized_seq}")
    print(f"Original CAI: {orig_cai:.4f}")
    print(f"Optimized CAI: {opt_cai:.4f}")
    print(f"Runtime: {runtime:.6f}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--sequence", type=str, required=True)
    args = parser.parse_args()
    run_tool_3(args.sequence.upper())
