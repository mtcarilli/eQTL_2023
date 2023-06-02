import pandas as pd
import numpy as np
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('--name', type=str)
args = parser.parse_args()

name = args.name

eigenvecs = pd.read_table(f'../genotypeData/{name}.eigenvec',header=None,sep=' ')
eigenvecs[1] = [ '0_' + eigenvecs[1][i] for i in range(len(eigenvecs))]
eigenvecs = eigenvecs.T.drop(0)
eigenvecs.to_csv(f'../data/{name}.eigenvec.T',header=None,sep = '\t',index=True)
print('DONE with eigen reformatting')
