from ast import expr_context
from pydoc import resolve
from neoget import _neoget
import dataclasses
import json
from pprint import pprint
from neo4j import GraphDatabase, Result
from dotenv import load_dotenv
import os
from typing import Dict, List, Tuple, TypedDict, Union, Callable
import operator
import sys
import pandas as pd
from Bio.PDB.MMCIFParser import FastMMCIFParser, MMCIFParser
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure
from Bio.PDB.Chain import Chain
import argparse
import itertools
from dataclasses import dataclass,field
from asyncio import run
import itertools
import numpy as np
flatten = itertools.chain.from_iterable
n1      = np.array

load_dotenv(dotenv_path='/home/rxz/dev/ribetl/.env'      )
STATIC_ROOT = os  .getenv   ('STATIC_ROOT')
#? Charges
AMINO_ACIDS = {
"ALA":0,
'ARG':1,
'ASN':0,
'ASP':-1,
'CYS':0,
'GLN':0,
'GLU':-1,
'GLY':0,
'HIS':0,
'ILE':0,
'LEU':0,
'LYS':1,
'MET':0,
'PHE':0,
'PRO':0,
'SER':0,
'THR':0,
'TRP':0,
'TYR':0,
'VAL':0,
'SEC':0,
'PYL':0}
NUCLEOTIDES = ['A','T','C','G','U']



# * Ecoli structs :  3j7z, 7k00, 6q97, 5j30
# ! yeast : 6z6n, 5mrc, 3j6b,6xir, 4u4n

#? PAR:
# 6az1 , 5tcu, 5iqr, 5el7,4wsd,4l71

#? KIR:
# 5afi, 4v8q, 4v5s,4v5g, 


@dataclass(unsafe_hash=True, order=True)
class ResidueLite     : 

      residue_name     : str               = field(              hash=True ,compare=False)
      residue_id       : int               = field(              hash=True ,compare=True)
      parent_strand_id : str               = field(              hash=True ,compare=False)
      
      @staticmethod
      def res2reslite(r:Residue):
          biopy_id_tuple = r.get_full_id()
          parent_chain   = biopy_id_tuple[2]
          resname        = r.resname
          resid          = r.id[1]
          return ResidueLite(resname,resid,parent_chain)

@dataclass
class BindingSiteChain: 
      sequence        : str
      nomenclature    : List[str]
      residues        : List[ResidueLite]

class BindingSite:
    def __init__(self,data:Dict[str,BindingSiteChain]) -> None:
        self.data : Dict[str,BindingSiteChain] = data

    def transpose_site(self, target_struct)->None:
        tgt = openStructutre(target_struct)
        print("Target:", tgt)
        ...

    def to_json(self,pathtofile:str)->None:
        with open(pathtofile, 'w') as outf: 
            serialized= {}
            for x in self.data.items():
                serialized.update({x[0]: dataclasses.asdict(x[1])})  
            json.dump(serialized,outf)
            print("Saved successfuly.")

    def to_csv(self,pathtofile:str)->None:

        k = [
        "chainname",
        "nomenclature",
        "residue_id",
        "residue_name"
        ]

        serialized = dict.fromkeys(k,[])
        
        pprint(serialized)

        # for chain in self.data:
        # # print(serialized['chainname'])
        # print("nomenclature", serialized['nomenclature'])
        # pd.DataFrame(serialized,columns=k).to_csv(pathtofile, index=False)
        # print("Saved as csv: ",pathtofile)

def getLigandResIds(ligchemid:str, struct: Structure)->List[Residue]:
    """Returns a list of dictionaries specifying each _ligand_ of type @ligchemid as a biopython-residue inside a given @struct."""
    """*ligchemids are of type https://www.rcsb.org/ligand/IDS"""

    ligandResidues: List[Residue] = list(filter(lambda x: x.get_resname() == ligchemid, list( struct.get_residues() )))
    return ligandResidues

async def matchStrandToClass(pdbid:str, strand_id:str)->List[str]:
    """Request Ban nomenclature classes from the db given a protein's entity_poly_strand_id."""

    RPCYPHER="""
    match (r:RibosomeStructure{{rcsb_id: "{}"}})-[]-(rp{{entity_poly_strand_id:"{}"}})-[]-(n:RPClass)
    return n.class_id""".format(pdbid.upper(), strand_id)

    RNACYPHER="""
    match (r:RibosomeStructure{{rcsb_id: "{}"}})-[]-(rp{{entity_poly_strand_id:"{}"}})-[]-(n:RNAClass)
    return n.class_id""".format(pdbid.upper(), strand_id)

    rpresp  = _neoget(RPCYPHER)
    rnaresp = _neoget(RNACYPHER)
    resp    = [*rpresp,*rnaresp]
    return resp

def openStructutre(pdbid:str, cifpath: str = None)->Structure:
    pdbid = pdbid.upper()
    if cifpath == None: cifpath = os.path.join(STATIC_ROOT, pdbid, f'{pdbid}.cif')
    return FastMMCIFParser(QUIET=True).get_structure(pdbid,cifpath)

def openLigand(pdbid:str, ligid:str, ligpath: str = None):
    pdbid = pdbid.upper()
    ligid = ligid.upper()

    if ligpath == None: ligpath = os.path.join(STATIC_ROOT, pdbid, f'LIGAND_{ligid}.json')
    with open(ligpath, 'rb') as infile: 
            data = json.load(infile)
    return  data

def get_lig_ids_struct(pdbid:str)->List[str]:
    pdbid=pdbid.upper() 
    db_response     = _neoget("""
    match (l:Ligand)-[]-(r:RibosomeStructure{{rcsb_id:"{pdbid}"}}) 
    return l.chemicalId, l.chemicalName """.format_map({ "pdbid":pdbid }))
    return db_response

def get_ligand_nbrs(
      lig_chemid     : str,
      ligand_residues: List[Residue],
      struct         : Structure,
    )-> BindingSite  : 
    """KDTree search the neighbors of a given list of residues(which constitue a ligand) 
    and return unique having tagged them with a ban identifier proteins within 5 angstrom of these residues. """

    pdbid = struct.get_id()
    ns           = NeighborSearch(list(struct.get_atoms()))
    nbr_residues = []

    # print(f"Ligand has {len(ligand_residues)} residue(s).")
    #? Searching Phase
    # a ligand consists of residues
    for lig_res in ligand_residues:
        for atom in lig_res.child_list:
                nbr_residues.extend(ns.search(atom.get_coord(), 10,level='R'))

    #? Filtering phase
    nbr_residues = list(set([* map(ResidueLite.res2reslite, nbr_residues) ]))
    nbr_residues = list(filter(lambda resl:resl.residue_name  in [*AMINO_ACIDS.keys(),  *NUCLEOTIDES], nbr_residues))

    # x: ResidueLite
    # lambda _: x.residue_name == "H_HOH"
    # print(nbr_residues)
    nbr_dict     = {}
    chain_names  = list(set(map(lambda _:  _.parent_strand_id, nbr_residues)))

    for c in chain_names:
        seq = _neoget(f"""match (n:RibosomeStructure{{rcsb_id:"{struct.get_id().upper()}"}})-[]-(r {{entity_poly_strand_id:"{c}"}}) 
        return r.entity_poly_seq_one_letter_code""")
        seq = [*flatten(seq)]
        nbr_dict[c]= BindingSiteChain(
             seq[0] if len(seq) >0 else '',
             list  ( flatten    (run       (matchStrandToClass(pdbid, c))) ),
            [* map   ( dataclasses.asdict,sorted( list       (filter    (lambda _ : _.parent_strand_id == c, nbr_residues)),key= operator   .attrgetter('residue_id')) )]
        )

    return BindingSite(nbr_dict)

def dropions(s:str): return False if "ion" in s[1].lower()else  True

def ParseLigand(ligid:str, rcsbid:str, force:bool=False):
    # try:
    outfile_json = os.path.join(STATIC_ROOT,rcsbid.upper(), f'LIGAND_{ligid}.json')
    # outfile_csv  = os.path.join(STATIC_ROOT,rcsbid.upper(), f'LIGAND_{ligid}.csv')

    if os.path.exists(outfile_json) and not force:
        print(outfile_json, "exists. Skipping.")
        exit(1)

    struct                     = openStructutre           (rcsbid                                            )
    residues   : List[Residue] = getLigandResIds          (ligid       , struct                              )
    bs:BindingSite = get_ligand_nbrs          (ligid,residues   , struct                               )
    pprint(bs.data)
    bs.to_json(outfile_json)
    # bs.to_csv(outfile_csv)

    # except Exception as e:
    #     print("Error:", e)

if __name__ =="__main__" :

    parser      = argparse.ArgumentParser(                                             )
    parser .add_argument ('-l','--ligand'   , type  = str                )
    parser .add_argument ('-s','--structure', type  = str ,required =True)
    parser .add_argument ('-V','--verbose'  , action='store_true'        )
    parser .add_argument ('-f','--force'  , action='store_true'        )
    args = parser.parse_args()

    VERBOSE     = args.verbose
    PDBID       = args.structure.upper()
    LIGID       = args.ligand
    FORCE       = args.force

    if LIGID == None:
        ligids = n1([ *filter(dropions, get_lig_ids_struct(PDBID) ) ],dtype=object)
        for l in ligids:
            ParseLigand(l[0].upper(), PDBID,force=True)
    else:
        print("Parsing ligand ", LIGID.upper(), PDBID)
        ParseLigand(LIGID.upper(),PDBID, force=FORCE)


#?  Goals:
# - encapsulate ligand neighborhood into a class
# - RC to mafft to match two sequencce
# - wrap the call into a class method
# - match back to the original residue ids
# - match forth to the destination residue ids