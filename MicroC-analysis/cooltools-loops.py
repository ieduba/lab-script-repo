import numpy as np
import pandas as pd
import cooler
import bioframe
import cooltools
from cooltools.lib.numutils import fill_diag
from optparse import OptionParser

opts = OptionParser()
opts.add_option('-m', help = '<.mcool> path to mcool file')
opts.add_option('-r', default = 10000, help = '<res> resolution in bp of mcool for loop calling')
options,arguments = opts.parse_args()

mc_file = options.m
res = options.r
clr = cooler.Cooler(f'{mc_file}::/resolutions/{res}')
expected = cooltools.expected_cis(clr)

dots_df = cooltools.dots(
    clr,
    expected=expected,
    max_loci_separation=10_000_000,
)
dots_df.to_csv(f'{mc_file}-{res}-loops.csv', index=False) 
