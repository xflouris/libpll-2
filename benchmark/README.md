This directory contains basic benchmark programs and scripts for performance benchmarks

# Basic micro kernel benchmark program 

Execute the `make` command to build primitive benchmark programs.

The benchmark program `partials-micro-benchmark` in particular is useful to measure the execution time of individual kernels. It has the following command line parameters:

|Argument (With Default)   | Meaning|
|---------------------|----|
alpha_values = [0.1] | List of used alpha values
-n-ites = 1000 | number of MSA sites
-n-states = 20 | number of model states
-p-invar = 0.0 | percentage of invariant sites
-n-categories = 4 | number of model categories
-n-itr = 1 | number of repeated executions within a benchmark run
-n-benchmark-repeat = 1 | Number of repeated benchmark runs (will produce a CSV output if it is greater 1)
-seed = <TIME_STAMP> | Random seed for generating a MSA 
-print-seq = 0| if 1, it will print the generated MSA to the console
-n-pmatrix-itr = 1 | Number of times the pmatrix kernel is repeatedly executed within a benchmark run


## Benchmark scripts

The python script `micro_benchmark_partials.py` is useful for generating a series of benchmark reports on a machine. 

For instance, the command

```
nohup python micro_benchmark_partials.py 10 500 50 out_dir > cmd.out 2>&1 &
```
Will execute `partials-micro-benchmark` with `-n-sites` ranging from 10 to 500 in increments by 50. It will write the resulting CSV file from each run into the folder `out_dir`.

### Commands used for plotting (per test host / server)

```
nohup python micro_benchmark_partials.py 10 500 10 partials_10_500_10 > partials_10_500_10.out 2>&1 &

nohup python micro_benchmark_partials.py 100 2000 100 
partials_100_2000_100 > partials_100_2000_100.out 2>&1 &

nohup python micro_benchmark_partials.py 100 100000 1000 partials_100_100000_1000 > partials_100_100000_1000.out 2>&1 &
```

## Visualizations

We used an R Markdown to generate the line charts. The Markdown file is called `PartialsMicroBenchmark.Rmd`