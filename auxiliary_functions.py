# -*- coding: utf-8 -*-
"""
Created on Thu Dec 14 13:52:03 2017

@author: LAdmin
"""
import platform
import shutil
from pathlib import Path

from .labeling_parameter import LabelingParameter


class ConservationCopy(LabelingParameter):
    # constants / parameter
    file_tag = "With_Conservation_Scores"
    # class variables

    def copy_conservation_score(self, chainOrigin, chainCopy):
        chainO = self.model[chainOrigin]
        chainC = self.model[chainCopy]
        for residue in chainO:
            resId = residue.get_id()
            try:
                bFact = residue.get_atoms().next().get_bfactor()
                for atom in chainC[resId]:
                    atom.set_bfactor(bFact)
            except:
                pass


class PDBOperation(LabelingParameter):
    # constants / parameter
    #    file_tag = 'With_Conservation_Scores'
    # class variables

    def set_chain(self, chainName):
        for chain in self.model:
            # print("Chain: <"+str(chain)+">")
            #            chain.id = chainName
            self.rename_chain(chain.id, chainName)

    def rename_chain(self, oldChain, newChain):
        chain = self.model[oldChain]
        try:
            chain.id = newChain
        except:
            # print("New chain <"+newChain+"> already exists. Chains will be merged.")
            chainGoal = self.model[newChain]
            for residue in chain:
                chainGoal.add(residue)
            self.model.detach_child(oldChain)


def find_msms(default_location: str) -> str:
    """In docker, we install it properly, but in development, we can just use the binary from tree"""
    # If the default `msms` is in PATH, use that
    if shutil.which(default_location):
        return default_location

    if platform.system() == "Windows":
        msms = r"msms.exe"
    elif platform.system() == "Linux":
        processor = platform.processor()
        if processor == "x86_64":
            msms = "msms.x86_64Linux2.2.6.1"
        else:
            raise NotImplementedError(f"Unsupported processor {processor}")
    elif platform.system() == "Darwin":
        msms = "msms.MacOSX.2.6.1"
    else:
        raise NotImplementedError(f"msms for {platform.system()}")

    # This is suboptimal because it breaks when moving this file, but as long as this isn't a package there's no
    # better solution
    msms = str(Path(__file__).parent.joinpath("msms").joinpath(msms))
    if not shutil.which(msms):
        raise FileNotFoundError(
            f"Neither `{default_location}` nor `{msms}` were executable"
        )

    return msms
