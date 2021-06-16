"""Class that is responsible for building and assessing proposed
   bonding patterns."""

from typing import Dict, List, Optional, Tuple

import numpy as np

from smu import dataset_pb2
from smu.parser import smu_utils_lib


class SmuMolecule:
  """Holds information about partially built molecules"""

  def __init__(self, hydrogens_attached: dataset_pb2.BondTopology,
               bonds_to_scores: Dict[Tuple[int, int], np.array]):
    """Class to perform bonding assessments.
    Args:
      hydrogens_attached: a BondTopology that has all atoms, and the bonds
        associated with the Hydrogen atoms.
      bonds_to_scores: A dict that maps tuples of pairs of atoms, to a
        numpy array of scores [0,3], for each possible bond type.
    Returns:
    """
    self._starting_bond_topology = hydrogens_attached
    natoms = len(hydrogens_attached.atoms)

    # For each atom, the maximum number of bonds that can be attached.
    self._max_bonds = np.zeros(natoms, dtype=np.int32)
    for i in range(0, natoms):
      self._max_bonds[i] = smu_utils_lib.ATOM_TYPE_TO_MAX_BONDS[hydrogens_attached.atoms[i]]

    # With the Hydrogens attached, the number of bonds to each atom.
    self._bonds_with_hydrogens_attached = np.zeros((natoms), dtype=np.int32)
    for bond in hydrogens_attached.bonds:
      self._bonds_with_hydrogens_attached[bond.atom_a] += 1
      self._bonds_with_hydrogens_attached[bond.atom_b] += 1

    self._current_bonds_attached = np.zeros((natoms), dtype=np.int32)

    # We turn bonds_to_scores into two arrays. So they can be iterated
    # via itertools.

    self._bonds = list(bonds_to_scores.keys())
    self._scores = list(bonds_to_scores.values())

  def _initialize(self):
    """Make the molecule reading for adding bonds between heavy atoms.
    """
    self._current_bonds_attached = self._bonds_with_hydrogens_attached

  def _place_bond(self, a1: int, a2: int, btype: int) -> bool:
    """Possibly add a new bond to the current config.

    If the bond can be placed, updates self._current_bonds_attached for 
    both `a`` and `a2`.
      Args:
        a1:
        a2:
        btype:
      Returns:
    print(f"Trying to place bond {btype} current {self._current_bonds_attached[a1]} and {self._current_bonds_attached[a2]}")
    """
    if self._current_bonds_attached[a1] + btype > self._max_bonds[a1]:
      return False
    if self._current_bonds_attached[a2] + btype > self._max_bonds[a1]:
      return False

    self._current_bonds_attached[a1] += btype
    self._current_bonds_attached[a2] += btype
    return True

  def generate_search_state(self) -> List[List[int]]:
    """For each pair of atoms, return a list of plausible bond types.

    This will be passed to itertools.product, which thereby enumerates all
    possible bonding combinations.
    Args:
    Returns:
      List of lists - one for each atom pair.
    """
    result: List[List[int]] = []
    for ndx, bond in enumerate(self._bonds):
      # For each pair of atoms, the plausible bond types - non zero score.
      plausible_types: List[int] = []
      for i, score in enumerate(self._scores[ndx]):
        if self._scores[ndx][i] > 0.0:
          plausible_types.append(i)

      result.append(plausible_types)

    return result

  def place_bonds(self, state: List[int]) -> Optional[dataset_pb2.BondTopology]:
    """Place bonds corresponding to `state`.

    Args:
      state: for each pair of atoms, the kind of bond to be placed.
    Returns:
      If successful, the score.
    """
    self._current_bonds_attached = np.copy(self._bonds_with_hydrogens_attached)

    result = dataset_pb2.BondTopology()
    result.CopyFrom(self._starting_bond_topology)    # only Hydrogens attached.
    #   result = self._starting_bond_topology  # Only the Hydrogens attached.
    result.score = 0.0

    for i in range(0, len(state)):
      a1 = self._bonds[i][0]
      a2 = self._bonds[i][1]
      btype = state[i]
      if not self._place_bond(a1, a2, btype):
        return None

      result.score += self._scores[i][btype]
      if btype > 0:
        result.bonds.append(
            dataset_pb2.BondTopology.Bond(atom_a=a1,
                                          atom_b=a2,
                                          bond_type=smu_utils_lib.INTEGER_TO_BOND_TYPE[btype]))

    return result
