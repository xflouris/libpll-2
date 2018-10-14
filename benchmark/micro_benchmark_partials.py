import subprocess
import os
import os.path
import sys


def run_benchmark(output_dir, n_sites):
    output_file = os.path.join("benchmark_%d.txt" % n_sites)

    output_file_path = os.path.join(output_dir, output_file)

    command = ['./obj/partials-micro-benchmark',
               '-n-benchmark-repeats=10',
               '-n-itr=10',
               '-n-sites=%d' % n_sites]
    print("Executing benchmark for -n-sites=%d" % n_sites)

    with open(output_file_path, 'w') as out_file:
        subprocess.call(command,
                        stdout=out_file,
                        shell=False)

if len(sys.argv) < 5:
    print("Usage: script <start> <end> <step> <out_dir>")
    sys.exit()

start = int(sys.argv[1])
end = int(sys.argv[2])
step = int(sys.argv[3])
output_dir = sys.argv[4]


if not os.path.exists(output_dir):
    os.makedirs(output_dir)

for n_site in range(start, end, step):
    run_benchmark(output_dir, n_site)