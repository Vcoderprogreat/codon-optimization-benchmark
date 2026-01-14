# mRNA Codon Optimization Benchmark 

## 1. Project Overview
I built a custom human mRNA optimizer (MyTool) and wanted to see how it actually stacks up against other tools. I focused on three things: speed, determinism, and usability.

## 2. Tools I Tested

* **MyTool (Stochastic Search):** A custom tool built by me that explores the sequence space using probabilistic selection to maximize the Codon Adaptation Index (CAI).
* **DNAChisel (Constraint-Solver):** An industry-standard tool. I used this to see how a constraint-based solver compares to my stochastic approach.
* **Greedy Baseline (Lookup Table):** A deterministic frequency-maximizer. 

## 3. Benchmarking Setup

I wanted to make sure the comparison was as fair as possible, so I set up a standardized pipeline in Google Colab:

* **Standardized Inputs:** I didn't just use random sequences. I created three specific FASTA files—short, medium, and long—to see how each tool handles different scales of data. Every tool had to process the exact same amino acids.
* **Fixing the Math (DNA vs RNA):** One of the biggest hurdles was that my tool was outputting RNA (`U`), but the CAI calculator and other tools use DNA (`T`). I wrote a normalization step to convert everything to DNA notation so the CAI scores would actually be comparable.
* **Testing from Scratch:** To make sure this code is reproducible for others, I tested the installation on a fresh Python 3.12 (Ubuntu) environment. This helped me catch a few "broken" library dependencies that I had to fix in the `requirements.txt` file.

## 4. Key Results

This summarizes the performance on the most complex test case (protein_long.faa):
MyTool, 1.00 Optimized CAI, 8.258 Runtime (seconds)
DNAChisel, 0.42 Optimized CAI, 2.561 Runtime (seconds)
Greedy Baseline, 0.40 Optimized CAI, 0.980 Runtime (seconds)

<img width="1584" height="614" alt="benchmark_summary" src="https://github.com/user-attachments/assets/55a1dfc9-1ea4-4fc9-b8fd-4c0304f60e90" />

## 5. Technical Discussion

- The "U vs T" Standardization: I ran into a major issue early on where my tool was outputting mRNA sequences (using Uracil), while the benchmarking tools and the CAI calculator were expecting DNA (using Thymine). This was making my tool's CAI scores look like 0.0. I fixed this by adding a normalization step in run_my_tool.py to swap all 'U's to 'T's so the comparison was scientifically fair.

- Dependency Issues (Python 3.12): When I tried to do a fresh install in the Colab environment, I realized the original CAI library is essentially "dead" and wouldn't install on modern Python versions. I had to pivot and migrate the whole project to use cai2, which is much more stable and reproducible.

- Proving Randomness: Since my algorithm uses a stochastic search, I needed to know if it produced the same result every time. I built check_determinism.py to run the same sequence three times and compare the outputs. It confirmed that my tool is stochastic (it explores different paths), while DNAChisel and the Greedy tool are deterministic (they give the exact same answer every time).

- Handling Nested Paths: Cloning the repo into Colab created some weird folder structures. I updated the scripts to use dynamic paths so that the benchmark doesn't break just because a folder name changed.

## 6. How to Reproduce
Clone
```bash
   git clone [https://github.com/](https://github.com/) Vcoderprogreat/codon-optimization-benchmark.git cd codon-optimization-benchmark
```
Install Dependencies
```bash
pip install -r requirements.txt
```
Execute Benchmark
```bash
python master_benchmark.py
