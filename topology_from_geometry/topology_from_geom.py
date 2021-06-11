from typing import Optional

import numpy as np

import apache_beam as beam

from smu import dataset_pb2
from smu.geometry import bond_length_distribution
from smu.parser import smu_utils_lib
import utilities

# The longest distance considered.
THRESHOLD = 2.0

class TopologyFromGeom(beam.DoFn):
  """Beam class for extracting BondTopology from Conformer protos."""
  def __init__(self, bond_lengths: bond_length_distribution.AllAtomPairLengthDistributions):
    self.bond_lengths = bond_lengths

  def _hydrogen_to_nearest_atom(self, bond_topology: dataset_pb2.BondTopology,
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

      bond = dataset_pb2.BondTopology.Bond(atom_a=a1, atom_b=a2, bond_type=dataset_pb2.BondTopology.BondType.BOND_SINGLE)
      result.bonds.append(bond)

    return result
  def process(self, conformer:dataset_pb2.Conformer):
    """
    """
    yield self._topology_from_geom(conformer.bond_topologies[0], conformer.optimized_geometry)

  def _topology_from_geom(self, bond_topology: dataset_pb2.BondTopology,
                         geometry: dataset_pb2.Geometry):
    """
    Args:
    Returns:
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
    starting_bond_topology = self._hydrogen_to_nearest_atom(bond_topology, distances)

    heavy_atoms = [a for a in bond_topology.atoms if a != dataset_pb2.BondTopology.AtomType.ATOM_H]
    number_heavy_atoms = len(heavy_atoms)

    # For each atom pair, a list of possible bond types.

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
        possible[(i, j)] = btypes

    return result

