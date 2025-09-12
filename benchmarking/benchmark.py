import argparse
import importlib.util
import sys
import os
import time
import matplotlib.pyplot as plt
from tqdm import tqdm
import platform
import psutil


# I/O counters
io_counts = {'read': 0, 'write': 0}

# Save original os.read and os.write
original_read = os.read
original_write = os.write

def get_hardware_info():
    info = {}

    # OS info
    info['platform'] = platform.platform()
    info['system'] = platform.system()
    info['node'] = platform.node()
    info['release'] = platform.release()
    info['version'] = platform.version()
    info['machine'] = platform.machine()
    info['processor'] = platform.processor()

    # CPU info
    info['cpu_count_physical'] = psutil.cpu_count(logical=False)
    info['cpu_count_logical'] = psutil.cpu_count(logical=True)
    # info['cpu_freq'] = psutil.cpu_freq().max if psutil.cpu_freq() else None

    # Memory info
    mem = psutil.virtual_memory()
    info['total_memory_GB'] = mem.total / (1024 ** 3)

    return info

def counting_read(fd, n):
    io_counts['read'] += 1
    return original_read(fd, n)

def counting_write(fd, b):
    io_counts['write'] += 1
    return original_write(fd, b)

# Decorator to mark benchmark target function
def benchmark_target(func):
    func._is_benchmark_target = True
    return func

def find_benchmark_targets(module):
    funcs = []
    for attr_name in dir(module):
        attr = getattr(module, attr_name)
        if callable(attr) and getattr(attr, '_is_benchmark_target', False):
            funcs.append((attr_name, attr))
    return funcs

def main():
    parser = argparse.ArgumentParser(description='Benchmark a python function using CPU time and I/O counts')
    parser.add_argument('--file', '-f', help='Python file to benchmark')
    parser.add_argument('--n_runs', '-n', type=int, help='Number of runs for benchmarking')
    args = parser.parse_args()

    # Dynamically load the module from the given file
    spec = importlib.util.spec_from_file_location("benchmod", args.file)
    module = importlib.util.module_from_spec(spec)
    sys.modules["benchmod"] = module
    spec.loader.exec_module(module)

    # Find all functions decorated with @benchmark_target
    target_funcs = find_benchmark_targets(module)
    if not target_funcs:
        print("Error: No functions decorated with @benchmark_target found in the file.")
        sys.exit(1)

    all_results = {}
    
    func_count = 0
    for fname, func in target_funcs:
        func_count += 1
        print(f"\n=== Benchmarking function {func_count}/{len(target_funcs)}: {fname} ===")

        # Patch os.read and os.write
        os.read = counting_read
        os.write = counting_write

        # Cold start
        global io_counts
        io_counts = {'read': 0, 'write': 0}
        start_cold = time.process_time()
        func()
        end_cold = time.process_time()
        print(f"Cold start CPU time: {end_cold - start_cold:.6f} s")

        # Benchmark runs
        cpu_times = []
        io_runs = []
        for _ in tqdm(range(args.n_runs), total=args.n_runs):
            io_counts = {'read': 0, 'write': 0}
            start = time.process_time()
            func()
            end = time.process_time()
            cpu_times.append(end - start)
            io_runs.append(io_counts.copy())

        all_results[fname] = cpu_times

    # Restore original os functions
    os.read = original_read
    os.write = original_write

    # Hardware info
    hi = get_hardware_info()
    python_info = f"Python [Compiler]: {platform.python_implementation()} {sys.version.split()[0]} [{platform.python_compiler().strip()}]"
    device_info = f"Platform: {hi['platform']}\n{hi['system']} {hi['release']}, {hi['cpu_count_physical']}-core CPU {hi['machine']}, {int(hi['total_memory_GB'])} GB RAM\n{python_info}"

    # Plot boxplot: one box per function
    plt.boxplot(all_results.values(), labels=all_results.keys())
    plt.title(f'CPU Time Elapsed (n={args.n_runs} runs per experiment)\n{device_info}')
    plt.ylabel('CPU time (seconds)')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    main()
