#"""Functions related to discerning the BondTopology from the geometry."""
from typing import Optional

import numpy as np

import apache_beam as beam

import utilities

from smu import dataset_pb2
from smu.geometry import bond_length_distribution
from smu.parser import smu_utils_lib

# The longest distance considered.
THRESHOLD = 2.0

class PairOfAtoms:
  def __init__(self, a1: int, a2: int, scores: np.array):
    self._a1 = a1
    self._a2 = a2
    self._scores = scores

  def place_bond(self, bond_topology: molecule_pb2.BondTopology):
    """
    Args:
      bond_topology:
      
    Returns:
    """


def hydrogen_to_nearest_atom(bond_topology: dataset_pb2.BondTopology,
                              distances: np.array) -> Optional[dataset_pb2.BondTopology]:
  """Generate a BondTopology that joins each Hydrogen atom to its nearest
      heavy atom.
  Args:
    bond_topology:
    distances:
  Returns:
  """
  result = dataset_pb2.BondTopology()
  result.atoms[:] = bond_topology.atoms
  natoms = len(bond_topology.atoms)
  for a1 in range(0, natoms):
    if bond_topology.atoms[a1] != dataset_pb2.BondTopology.AtomType.ATOM_H:
      continue

    shortest_distance = 1.0e+30
    closest_heavy_atom = -1
    for a2 in range(0, natoms):
      if bond_topology.atoms[a2] == dataset_pb2.BondTopology.AtomType.ATOM_H:
        continue

      if distances[a1, a2] >= THRESHOLD:
        continue

      if distances[a1, a2] < shortest_distance:
        shortest_distance = distances[a1, a2]
        closest_heavy_atom = a2

    if closest_heavy_atom < 0:
      return None

    bond = dataset_pb2.BondTopology.Bond(atom_a=a1,
                                         atom_b=a2,
                                         bond_type=dataset_pb2.BondTopology.BondType.BOND_SINGLE)
    result.bonds.append(bond)

  return result


class TopologyFromGeom(beam.DoFn):
  """Beam class for extracting BondTopology from Conformer protos."""

  def __init__(self, bond_lengths: bond_length_distribution.AllAtomPairLengthDistributions):
    super().__init__()
    self.bond_lengths = bond_lengths

  def process(self, conformer: dataset_pb2.Conformer):
    """Called by Beam.
    Args:
      conformer:
    Yields:
      dataset_pb2.TopologyMatches
    """
    yield self._topology_from_geom(conformer.bond_topologies[0], conformer.optimized_geometry)

  def _topology_from_geom(self, bond_topology: dataset_pb2.BondTopology,
                          geometry: dataset_pb2.Geometry) -> dataset_pb2.TopologyMatches:
    """Discern the topology from `geometry`.
    Args:
      bond_topology: BondTopology
      geometry: atomic positions.
    Returns:
      dataset_pb2.TopologyMatches
    """

    # This will be yield'd.
    result = dataset_pb2.TopologyMatches()
    natoms = len(bond_topology.atoms)
    if natoms == 1:
      return result

    bonded = utilities.bonded(bond_topology)
    distances = utilities.distances(geometry)

    # First join each Hydrogen to its nearest heavy atom, thereby
    # creating a starting BondTopology from which all others can grow
    starting_bond_topology = hydrogen_to_nearest_atom(bond_topology, distances)

    heavy_atoms = [a for a in bond_topology.atoms if a != dataset_pb2.BondTopology.AtomType.ATOM_H]
    number_heavy_atoms = len(heavy_atoms)

    # For each atom pair, a list of possible bond types.
    # Key is a tuple of the two atom numbers, value is an np.array
    # with the score for each bond type.

    possible = {}
    for i in range(0, number_heavy_atoms):
      itype = heavy_atoms[i]
      for j in range(i + 1, number_heavy_atoms):
        dist = distances[i, j]
        if dist > THRESHOLD:
          continue
        jtype = heavy_atoms[j]
        btypes = np.full((4), 0.0, np.float32)
        print(f"Looking for pdfs of #{itype} {jtype}")
        for btype in range(0, 4):
          btypes[btype] = self.bond_lengths.pdf_length_given_type(itype, jtype, btype, dist)

        if np.count_nonzero(btypes) > 0:
          possible[(i, j)] = btypes

    # For each possibly bonded pair of atoms, we now have a list of scores for each of
    # the plausible bond types.

    atoms = []
    btypes = []
    scores = []
    for key,value in possible.items():
      atoms.append(key)

      
      

    return result
