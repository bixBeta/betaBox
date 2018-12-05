#!/usr/bin/python

import pandas as pd
import os;

os.getcwd( )
os.chdir('/Users/bix/Desktop/subset.0021/0036')
os.listdir()
#pd.read_csv("cpGs.genes.csv")
barCode = pd.read_csv("key.csv")

betas = pd.read_csv("bb.csv")

betas.set_index('cpGs')

final = betas.rename(columns=dict(zip(barCode['centrix'], barCode['Sample_Name'])))

final.to_csv("TEST.csv")
