#!/usr/bin/python
import os
import sys
import string
import re
import numpy as np
import pandas as pd

def extract_benchmark(file_name):
    prefix = string.replace(os.path.basename(file_name), '.txt', '')
    df = pd.read_csv(file_name)
    del df['real']
    del df['sys']
    df.columns = ['Itr'] + ['%s (%s)' % (prefix, col) for col in df.columns if col != 'Itr']
    return df

def extract_df(file_names):
    raw_data = [extract_benchmark(file_name) for file_name in file_names]
    df = reduce(lambda x, y: pd.merge(x, y, on = 'Itr'), raw_data)
    df.index = df.Itr
    del df['Itr']
    return df

print extract_df(sys.argv[1:]).to_csv(sep=";",decimal=',')