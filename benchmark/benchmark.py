#!/usr/bin/python
import subprocess
import re
import sys
import shlex

real_pattern = re.compile('real\s+(\d+)m([\d.]+)s')
user_pattern = re.compile('user\s+(\d+)m([\d.]+)s')
sys_pattern =  re.compile('sys\s+(\d+)m([\d.]+)s')

def mean(numbers):
    return float(sum(numbers)) / max(len(numbers), 1)

def benchmark_time(command_line):

    proc = subprocess.Popen(shlex.split('bash -c "time %s"' % " ".join(command_line)),
                            stderr=subprocess.PIPE,
                            stdout=subprocess.PIPE,
                            shell=False)

    (outdata, errdata) = proc.communicate()

    if proc.returncode != 0:
        exit("Benchmark target exited with error code: %d" % proc.returncode)

    real_time = 0
    user_time = 0
    sys_time = 0

    for line in errdata.rsplit("\n"):
        real_m = real_pattern.match(line.rstrip())
        if real_m:
            real_time = float(real_m.group(1)) * 60 + float(real_m.group(2))
        user_m = user_pattern.match(line.rstrip())
        if user_m:
            user_time = float(user_m.group(1)) * 60 + float(user_m.group(2))
        sys_m = sys_pattern.match(line.rstrip())
        if sys_m:
            sys_time = float(sys_m.group(1)) * 60 + float(sys_m.group(2))

    return {'real': real_time,
            'user': user_time,
            'sys': sys_time}

print "Itr,real,user,sys"

results = [(x, benchmark_time(sys.argv[1:])) for x in range(10)]
for entry in results:
    print "%d,%.4f,%.4f,%.4f" % (entry[0], entry[1]['real'], entry[1]['user'], entry[1]['sys'])

print "MEAN,%.4f,%.4f,%.4f" % (mean([x[1]['real'] for x in results]),
                               mean([x[1]['user'] for x in results]),
                               mean([x[1]['sys'] for x in results]))
