import subprocess
import time
import csv

sequences = [
    "ATGTTTCCCGGG",
    "ATGCGTACGTAG",
    "ATGAAACCCGGGTTT"
]

output_file = "benchmark_results.csv"

with open(output_file, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["Sequence", "Tool", "Optimized_Sequence", "Original_CAI", "Optimized_CAI", "Runtime_sec"])

for seq in sequences:
    start = time.time()
    # Call your CLI tool
    result = subprocess.run(
        ["python", "codon_optimizer_cli.py", "-s", seq, "-c", "human_codon_usage.csv", "--method", "stochastic"],
        capture_output=True,
        text=True
    )
    end = time.time()
    runtime = end - start

  
    lines = result.stdout.splitlines()
    orig_seq = lines[0].split(": ")[1]
    opt_seq = lines[1].split(": ")[1]
    orig_cai = float(lines[2].split(": ")[1])
    opt_cai = float(lines[3].split(": ")[1])

   
    with open(output_file, "a", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([seq, "MyTool", opt_seq, orig_cai, opt_cai, runtime])

print(f"Benchmark finished. all the  results saved to {output_file}")
