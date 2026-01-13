#!/usr/bin/env python3
import csv
import math
import random
import argparse
from Bio.Seq import Seq

# Load codon usage table

def load_codon_usage(csv_path):
    codon_usage = {}
    with open(csv_path, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            aa = row['Amino Acid']
            codon = row['Codon'].replace('T','U')  # make RNA-like if needed
            freq = float(row['Frequency'])
            if aa not in codon_usage:
                codon_usage[aa] = {}
            codon_usage[aa][codon] = freq
    return codon_usage


# CAI Calculation

def calculate_cai(dna_seq, codon_usage):
    seq = str(dna_seq)
    trim_len = len(seq) - len(seq) % 3
    seq = seq[:trim_len]
    codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
    weights = []
    for codon in codons:
        aa = str(Seq(codon).translate())
        if aa == '*' or aa not in codon_usage:
            continue
        codon_freqs = codon_usage[aa]
        max_freq = max(codon_freqs.values())
        weight = codon_freqs.get(codon, 0) / max_freq
        # Avoid zero weight
        if weight > 0:
            weights.append(weight)
    if not weights:
        return 0
    return math.exp(sum(math.log(w) for w in weights) / len(weights))


# Rare codons

def rare_codon_count(dna_seq, codon_usage, threshold=0.1):
    seq = str(dna_seq)
    trim_len = len(seq) - len(seq) % 3
    codons = [seq[i:i+3] for i in range(0, trim_len, 3)]
    count = 0
    for codon in codons:
        aa = str(Seq(codon).translate())
        if aa == '*' or aa not in codon_usage:
            continue
        if codon_usage[aa].get(codon, 0) < threshold:
            count += 1
    return count


# Homopolymers

def homopolymer_count(dna_seq, run_length=5):
    import re
    seq = str(dna_seq)
    matches = re.findall(r'(A{%d,}|T{%d,}|G{%d,}|C{%d,})' % (run_length, run_length, run_length, run_length), seq)
    return len(matches)


# Total scoring

def score_sequence(dna_seq, codon_usage, target_gc=0.5):
    cai_score = calculate_cai(dna_seq, codon_usage)
    seq = str(dna_seq)
    gc_frac = (seq.count('G') + seq.count('C')) / max(len(seq), 1)
    gc_penalty = abs(gc_frac - target_gc)
    rare_penalty = rare_codon_count(dna_seq, codon_usage)
    homo_penalty = homopolymer_count(dna_seq)
    total_score = cai_score - 0.5*gc_penalty - 0.05*rare_penalty - 0.2*homo_penalty
    return total_score


# Random synonymous swap

def random_synonymous_swap(dna_seq, codon_usage):
    seq = list(str(dna_seq))
    codon_index = random.randint(0, (len(seq)//3)-1)
    start = codon_index*3
    codon = ''.join(seq[start:start+3])
    aa = str(Seq(codon).translate())
    if aa == '*' or aa not in codon_usage:
        return ''.join(seq)
    choices = list(codon_usage[aa].keys())
    new_codon = random.choice(choices)
    seq[start:start+3] = list(new_codon)
    return ''.join(seq)


# Greedy optimization

def greedy_optimize(dna_seq, codon_usage):
    seq = str(dna_seq)
    trim_len = len(seq) - len(seq) % 3
    seq = seq[:trim_len]
    codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
    optimized = []
    for codon in codons:
        aa = str(Seq(codon).translate())
        if aa == '*' or aa not in codon_usage:
            optimized.append(codon)
            continue
        codon_freqs = codon_usage[aa]
        best_codon = max(codon_freqs, key=codon_freqs.get)
        optimized.append(best_codon)
    return ''.join(optimized)


# Stochastic optimization

def stochastic_optimize(dna_seq, codon_usage, iterations=2000, temperature=1.0, cooling_rate=0.995):
    current_seq = greedy_optimize(dna_seq, codon_usage)
    current_score = score_sequence(current_seq, codon_usage)
    best_seq = current_seq
    best_score = current_score
    for _ in range(iterations):
        new_seq = random_synonymous_swap(current_seq, codon_usage)
        new_score = score_sequence(new_seq, codon_usage)
        diff = new_score - current_score
        if diff > 0 or random.random() < math.exp(diff / max(temperature, 1e-8)):
            current_seq = new_seq
            current_score = new_score
        if current_score > best_score:
            best_seq = current_seq
            best_score = current_score
        temperature *= cooling_rate
    return best_seq


# CLI

def main():
    parser = argparse.ArgumentParser(description="Codon Optimizer CLI")
    parser.add_argument("-s", "--sequence", type=str, required=True, help="Input DNA sequence")
    parser.add_argument("-c", "--codon_csv", type=str, required=True, help="Path to codon usage CSV")
    parser.add_argument("--method", type=str, default="stochastic", choices=["greedy","stochastic"], help="Optimization method")
    args = parser.parse_args()

    codon_usage = load_codon_usage(args.codon_csv)
    if args.method == "greedy":
        optimized = greedy_optimize(args.sequence, codon_usage)
    else:
        optimized = stochastic_optimize(args.sequence, codon_usage)

    print("Original sequence:", args.sequence)
    print("Optimized sequence:", optimized)
    print("Original CAI:", calculate_cai(args.sequence, codon_usage))
    print("Optimized CAI:", calculate_cai(optimized, codon_usage))

if __name__ == "__main__":
    main()
