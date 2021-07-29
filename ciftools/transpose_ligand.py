import dataclasses
import operator
import json
from pprint import pprint
import re
import struct
from neo4j import GraphDatabase, Result
from dotenv import load_dotenv
import os
from typing import Dict, List, Tuple, TypedDict, Union, Callable
import operator
import sys
from Bio.PDB.MMCIFParser import FastMMCIFParser, MMCIFParser
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure
from Bio.Align import MultipleSeqAlignment
from Bio import pairwise2
import Bio.AlignIO
import argparse
import itertools
from dataclasses import dataclass,field
from asyncio import run
import itertools
import neo4j
import numpy as np
import render_binding_site as bsite
flatten = itertools.chain.from_iterable
n1      = np.array



# def seqformat(seq:str):
# 	print("Seqlen is ", len(seq))
# 	if len(seq)<=80:
# 		return seq
# 	else:
# 		subs = len(seq)//80
# 		print("Found ", subs,"Subdivisions")
# 		for g in range(1,subs+1):
# 			seq =  seq[:g*80] + "\n" + seq[g*80+1:] 
# 	print(seq)
# 	return seq
		
		

	




#! Include sequence into the ligand profiles
#! exclude the ligand itself from the transposition
#* match matchable chains from target to prediction




#3j7z  ERY ---> to 5hl7, 3j9z, 6otl
#6AZ1  PAR ---> 5t2a l.donovani

origin     = str( sys.argv[1] ).upper()
project_to = str( sys.argv[2] ).upper()

with open(f'../static/{origin}/LIGAND_PAR.json', 'rb') as infile:
	data = json.load(infile)


bs = bsite.BindingSite(data)

#? Protocol: 
#? for every chain in origin, grab same nomenclature in tgt
#? track ids back, apply forth
#? apply ids 


origin_chains = {

}

target_chains = {

}


for chain in bs.data:
	if len( bs.data[chain]['nomenclature'] ) <1:
		continue
	else:
		res_align_mapping :Dict[int,Tuple[int,int]] = {
			resid:(-1,-1) for  resid in [*map(lambda x : x['residue_id'], bs.data[chain]['residues'])]
		}
		origin_chains[bs.data[chain]['nomenclature'][0]] = {
			'strand':chain,
			'seq': bs.data[chain]['sequence'],
			# 'ids': [*map(lambda x : x['residue_id'], bs.data[chain]['residues'])]
			'ids': res_align_mapping
		}


for nom in origin_chains:
	name_matches = []
	cypher       = f"""match (n:RibosomeStructure {{rcsb_id:"{project_to}"}})-[]-(c)-[]-(r {{class_id:"{nom}"}}) return c.entity_poly_seq_one_letter_code, c.entity_poly_strand_id"""
	response     = bsite._neoget(cypher)
	if len( response )  < 1:
		print(f"No chain-class matches for {nom} in {project_to} in the database.")
		continue
	else:
		match = response[0]

	strand              = match[1]
	seq                 = match[0]
	target_chains[nom] ={
		'strand': strand,
		'seq'   : seq
	}

# ?The driving data strucutre for matching is Dict[int,Tuple(int,int)]
# ?where the key index is the index in the original_unaligned, tup[0] is the index in ori

class SeqMatch:

	def __init__(self,sourceseq:str,targetseq:str) -> None:
		"""S1 and S2"""
		self.seq1    = sourceseq
		self.seq2    = targetseq
		_            = pairwise2.align.globalxx(origin_seq,target_seq, one_alignment_only=True)
		self.src_aln = _[0].seqA
		self.tgt_aln = _[0].seqB


	def backwards_match(self, alntgt:str, resid:int):
		"""Returns the target-sequence  index of a residue in the (aligned) target sequence"""
		counter_proper = 0
		for i,char in enumerate(alntgt):
			if i == resid:
				if char == "-":
					return None
				return counter_proper
			if char =='-':
				continue
			else: 
				counter_proper  +=1

	def forwards_match(self,alnsrc:str, resid:int):
		"""Returns the index of a source-sequence residue in the aligned source sequence."""
		count_proper = 0
		for alignment_indx,char in enumerate( alnsrc ):
			if count_proper == resid:
				return alignment_indx
			if char =='-':
				continue
			else: 
				count_proper  +=1

	@staticmethod
	def hl_subseq(sequence:str, subsequence:str, index:int=None):
		"""Highlight subsequence"""
		CRED = '\033[91m'
		CEND = '\033[0m'
		_ = [ ]
		if index != None:
			return sequence[:index-1] + CRED + sequence[index] + CEND +sequence[index+1:]
		for item in re.split(re.compile(f'({subsequence})'),sequence):
			if item == subsequence:
				_.append(CRED + item + CEND)
			else:
				_.append(item)
		return ''.join(_)

	@staticmethod
	def hl_ixs(sequence:str,  ixs:List[int]):
		"""Highlight indices"""
		CRED = '\033[91m'
		CEND = '\033[0m'
		_ = ''
		for i,v in enumerate(sequence):
			if i in ixs: _ += CRED + v +CEND
			else: 	 	 _ += v
		return _


		
for  name in origin_chains:
	if name not in target_chains:
		continue
	source_ids = origin_chains[name]['ids']

	origin_seq     = origin_chains[name]['seq']
	target_seq     = target_chains[name]['seq']

	aligned        = pairwise2.align.globalxx(origin_seq,target_seq, one_alignment_only=True)

	source_aligned = aligned[0].seqA
	target_aligned = aligned[0].seqB

	
	aligned_ids = []
	target_ids  = []

	for orig in source_ids:
		aligned_ids.append(forwards_match(source_aligned,orig))

	for algn in aligned_ids:
		target_ids.append(backwards_match(target_aligned,algn))
	

	# print("Aligned ids", source_ids.keys())
	# print("To ------->", aligned_ids)
	# print("To targets :", target_ids)

	print(f"Protein {name}")

	print("ORG   : \t",hl_ixs(origin_seq,source_ids),"\n")
	print("ORG AL: \t",hl_ixs(source_aligned,aligned_ids),"\n")
	print("TGT AL: \t",hl_ixs(target_aligned,aligned_ids),"\n")
	print("TGT   : \t",hl_ixs(target_seq, target_ids),"\n")


