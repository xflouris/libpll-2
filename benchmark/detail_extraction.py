#!/lrz/sys/tools/python/intelpython27_u4/intelpython2/bin/python
import os
import sys
import string
import re
import numpy as np
import pandas as pd

detail_pattern = re.compile('([^:]+):\s+(\d+\.\d+)')

def extract_tuples(file_name):
    with open(file_name) as f:
        lines = (line.rstrip() for line in f.readlines()) 
        lines = (line for line in lines if line)
        return (string.replace(os.path.basename(file_name), '.txt', ''),
                {key: value for (key, value) in [extract_tuple(x) for x in lines]})

def extract_tuple(line):
    detail_m = detail_pattern.match(line.rstrip())
    if detail_m:
        return (detail_m.group(1), float(detail_m.group(2)))

def extract_df(file_names):
    raw_data = {key: value for (key, value) in 
                    [extract_tuples(file_name) for file_name in file_names]}
    return pd.DataFrame(raw_data)

print extract_df(sys.argv[1:]).to_csv(sep=";",decimal=',')