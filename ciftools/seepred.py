from ast import arg, parse
import os,sys
import json
from pymol import cmd
import argparse



parser = argparse.ArgumentParser()
parser.add_argument('-L','--ligand', type=str)
parser.add_argument('-s','--source', type=str)
parser.add_argument('-t','--target', type=str)

args = parser.parse_args()
LIG  = args.ligand.upper()
SRC  = args.source.upper()
TGT  = args.target.upper()


cmd.load('/home/rxz/dev/ribetl/static/{}/{}.cif'.format(TGT,TGT))

with  open('/home/rxz/dev/ribetl/static/{}/PREDICTION_{}_{}_{}.json'.format(TGT,LIG,SRC,TGT), 'rb') as infile:
	data = json.load(infile)

cmd.color('gray','all')
for chain in data:
	# print(data[ chain ]['source']['strand'])
	tgt_strand = data[ chain ]['target']['strand']
	# src_resids =  data[ chain ]['source']['src_ids']
	tgt_resids =  data[ chain ]['target']['tgt_ids']

	for resid in tgt_resids:
		cmd.color('cyan', f'c. {tgt_strand} and resi {resid}')
		