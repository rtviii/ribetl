import os,sys
import json
from pymol import cmd


# pdbid  = sys.argv[1].upper()
# ligand = sys.argv[2].upper()


with  open(f'/home/rxz/dev/ribetl/static/4UG0/PREDICTION_PAR_6AZ1_4UG0.json', 'rb') as infile:
	data = json.load(infile)

# print(data)
cmd.color('gray','all')
for chain in data:

	# print(data[ chain ]['source']['strand'])
	tgt_strand = data[ chain ]['target']['strand']

	# src_resids =  data[ chain ]['source']['src_ids']
	tgt_resids =  data[ chain ]['target']['tgt_ids']

	for resid in tgt_resids:
		cmd.color('cyan', f'c. {tgt_strand} and resi {resid}')
		