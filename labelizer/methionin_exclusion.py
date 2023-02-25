import logging
from pathlib import Path
from typing import Optional, Dict, Tuple

from Bio.PDB import NeighborSearch, Selection
from Bio.PDB import PDBParser, Model
from Bio.PDB.ResidueDepth import get_surface, residue_depth

from . import Labelizer
from .auxiliary_functions import find_msms
from .labeling_parameter import LabelingParameter

logger = logging.getLogger(__name__)


class AAExclusion(LabelingParameter):
    """Gives a very bad score to residues that a too close to a specific bad amino acid.
    The main usage is to ban MET to close to the labeled position as this can lead to
    quenching (TODO: reference once online). The other AA must be solvent exposed since
    we assume this is how contact and quenching happens."""

    # This was originally Methionine Exclusion we just keep the tag
    file_tag = "me"
    # Uppercase three letter code of the bad AA we don't want to have close
    aa: str = "MET"
    # To be excluded, a bad aa must be in the range of this many Ångström
    max_bad_aa_distance = 6
    # If the bad aa is buried (solvent exposure is lower than this value) we ignore it,
    # given as 1 / residue depth
    min_bad_aa_solvent_exposure = 0.4
    # MSMS location, the default will be checked before the custom logic
    msms: str = "msms"
    # We compute the solvent exposure values for each bad aa once and then story it here
    bad_aa_se_values: Dict[Tuple[str, int, str, Tuple[str, int, str]], float]
    model: Model
    neighbor_search: NeighborSearch

    def __init__(self) -> None:
        super().__init__()
        self.msms = find_msms(self.msms)

    def set_up(
        self,
        labelizer: Labelizer,
        protein: str,
        sensitivity: str = "m",
        para_model: Optional[str] = None,
    ):
        super().set_up(labelizer, protein, sensitivity, para_model)
        filepath = Path(self.file_path).parent.joinpath(self.identifier + ".pdb")
        # PDB warnings aren't useful for us
        parser = PDBParser(QUIET=True)
        model = parser.get_structure(protein, filepath)[0]
        self.set_up_freestanding(model)

    def set_up_freestanding(self, model: Model):
        """For use outside the labelizer codebase"""
        self.model = model
        surface = get_surface(self.model, MSMS=self.msms)
        # A is for atoms (not for the chain)
        atoms = Selection.unfold_entities(self.model, "A")
        self.neighbor_search = NeighborSearch(atoms)
        self.bad_aa_se_values = {}
        for chain in self.model:
            for target_aa in chain:
                if target_aa.resname != self.bad_aa:
                    continue

                # residue depth and surface exposure are inversely correlated
                se_value = 1 / max(residue_depth(target_aa, surface), 10**-10)
                self.bad_aa_se_values[target_aa.full_id] = se_value

    def is_too_close_to_met(self, chain_id: str, res_id_int: int) -> bool:
        residue = self.model[chain_id][res_id_int]
        target_ca = residue["CA"]
        close_atoms = self.neighbor_search.search(
            target_ca.coord, self.max_bad_aa_distance
        )
        is_too_close_to_met = False
        for atom in close_atoms:
            # For simplicity, we identify nearby residues by Cɑ-Cɑ distance only
            if atom.get_id() != "CA" or atom.parent.resname != self.bad_aa:
                continue

            # If residue is a BAD_AA residue, we'd mutate this out again anyway
            if atom.parent == residue:
                continue

            close_res_se_value = self.bad_aa_se_values[atom.parent.full_id]
            if close_res_se_value > self.min_bad_aa_solvent_exposure:
                is_too_close_to_met = True
                logger.debug(
                    f"{chain_id} {res_id_int} is too close to "
                    f"{atom.parent.parent.id} {atom.parent.id[1]} "
                    f"with solvent exposure {close_res_se_value}"
                )
                break
        return is_too_close_to_met

    def calc_parameter_score(self, chain_id: str, res_id_int: int) -> float:
        if self.is_too_close_to_met(chain_id, res_id_int):
            return 0.001
        else:
            return 0.999
