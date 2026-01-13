import subprocess
import time
import sys

def run_custom_optimizer(dna_seq):
    """
    Wrapper for MyTool (Stochastic/Greedy optimization).
    """
    start_time = time.time()
    
    try:
        result = subprocess.run(
            ["python", "codon_optimizer_cli.py", "-s", dna_seq, "-c", "human_codon_usage.csv", "--method", "stochastic"],
            capture_output=True,
            text=True,
            check=True
        )
        runtime = time.time() - start_time
        
        lines = result.stdout.strip().split('\n')
        
        print(lines[0]) # Original
        print(lines[1]) # Optimized
        print(lines[2]) # Original CAI
        print(lines[3]) # Optimized CAI
        print(f"Runtime: {runtime:.6f}")

    except Exception as e:
        print(f"Error running MyTool: {e}")

if __name__ == "__main__":
    if len(sys.argv) < 3 or sys.argv[1] != "-s":
        print("Usage: python run_my_tool.py -s <SEQUENCE>")
    else:
        run_custom_optimizer(sys.argv[2].upper())
