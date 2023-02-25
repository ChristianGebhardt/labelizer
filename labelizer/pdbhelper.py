# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 16:44:03 2020

@author: gebha
"""
import warnings

from Bio import PDB
from Bio.PDB import PDBIO


def remove_hetatoms(model):
    for chain in model:
        het_atoms = []
        for residue in chain:
            res_Id = residue.get_id()
            if res_Id[0] != " ":
                # Optional differntiation between water, other HETATM (e.g. ligands, ATP, heavy metal ions)
                if res_Id[0] == "W":
                    pass
                if res_Id[0][:2] == "H_":
                    pass
                het_atoms.append((res_Id[1], res_Id))
        het_atoms.sort(key=lambda x: x[0])
        het_atoms.reverse()
        for het_id in het_atoms:
            chain.detach_child(het_id[1])
    return model


def load_pdb(identifier, path, remove_hetatm=True):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        parser = PDB.PDBParser(PERMISSIVE=1)
        structure = parser.get_structure(identifier, path)
        model = structure[0]
        if remove_hetatm:
            model = remove_hetatoms(model)
    return structure, model


def save_pdb(structure, file_path):
    with open(file_path, "w") as f:
        io = PDBIO()
        io.set_structure(structure)
        io.save(f)


# class PDBOperation():

# #     def set_chain(self,chainName):
# #         for chain in self.model:
# #             print("Chain: <"+str(chain)+">")
# # #            chain.id = chainName
# #             self.rename_chain(chain.id,chainName)

# #     def rename_chain(self,oldChain, newChain):
# #         chain = self.model[oldChain]
# #         try:
# #             chain.id = newChain
# #         except:
# #             print("New chain <"+newChain+"> already exists. Chains will be merged.")
# #             chainGoal = self.model[newChain]
# #             for residue in chain:
# #                 chainGoal.add(residue)
# #             self.model.detach_child(oldChain)
#     pass

# TBD MOVE TO PDBOperation

#     """
#     Set and alter chain identifier
#     """
# #    from AuxiliaryFunction import PDBOperation
# #    import os
#     def set_chain(self, prot_nbr, chain_letter):
# #    for protein in proteins:
#         newChain = PDBOperation()

#         path1 = os.path.join(self.folder, self.proteins[prot_nbr] +'.pdb')
#         newChain.load_pdb(path1,self.proteins[prot_nbr])

#         #set chain name or rename chain name
#         newChain.set_chain(chain_letter)
#     #    newChain.rename_chain("D","A")

#         newChain.file_tag = chain_letter
#         #save modified pdb file
#         newChain.save_pdb()

#     def rename_chain(self, prot_nbr, chain_letter_old, chain_letter_new):
# #    for protein in proteins:
#         newChain = PDBOperation()

#         path1 = os.path.join(self.folder, self.proteins[prot_nbr] +'.pdb')
#         newChain.load_pdb(path1,self.proteins[prot_nbr])

#         #set chain name or rename chain name
# #        newChain.set_chain(chain_letter)
#         newChain.rename_chain(chain_letter_old,chain_letter_new)

#         newChain.file_tag = chain_letter_new
#         #save modified pdb file
#         newChain.save_pdb()
