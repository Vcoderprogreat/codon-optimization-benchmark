import os
import subprocess
import csv
import time
from Bio import SeqIO


INPUT_DIR = "inputs"
TOOLS_DIR = "tools"
OUTPUT_FILE = "benchmark_results.csv"


TOOLS = [
    {"name": "MyTool", "script": "run_my_tool.py"},
    {"name": "DNAChisel", "script": "run_tool2.py"},
    {"name": "GreedyBaseline", "script": "run_tool3.py"}
]

def run_benchmark():
 
    with open(OUTPUT_FILE, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Sequence_Name", "Tool", "Optimized_Seq", "Original_CAI", "Optimized_CAI", "Runtime_Sec"])

    print("--- Starting mRNA Codon Optimization Benchmark ---")


    for filename in os.listdir(INPUT_DIR):
        if filename.endswith(".faa"):
            path = os.path.join(INPUT_DIR, filename)
            
     
            for record in SeqIO.parse(path, "fasta"):
              
                protein_seq = str(record.seq)
                dna_input = "".join(["ATG" if aa=='M' else "GCT" for aa in protein_seq]) 

                for tool in TOOLS:
                    print(f"Running {tool['name']} on {filename}...")
                    
                    script_path = os.path.join(TOOLS_DIR, tool['script'])
                    
                
                    start = time.time()
            
                    proc = subprocess.run(
                        ["python", script_path, "-s", dna_input],
                        capture_output=True, text=True
                    )
                    runtime = time.time() - start

                    if proc.returncode == 0:
                        lines = proc.stdout.strip().split('\n')
                    
                        try:
                            opt_seq = lines[1].split(": ")[1]
                            orig_cai = lines[2].split(": ")[1]
                            opt_cai = lines[3].split(": ")[1]

                            with open(OUTPUT_FILE, "a", newline="") as f:
                                writer = csv.writer(f)
                                writer.writerow([record.id, tool['name'], opt_seq, orig_cai, opt_cai, f"{runtime:.4f}"])
                        except:
                            print(f"Error parsing output from {tool['name']}")
                    else:
                        print(f"Tool Error: {proc.stderr}")

    print(f"--- Benchmark Complete! Results in {OUTPUT_FILE} ---")

if __name__ == "__main__":
    run_benchmark()
