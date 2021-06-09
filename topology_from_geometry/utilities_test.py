# Tester for SMU utilities functions.

from parameterized import parameterized, parameterized_class

from rdkit import Chem

from smu import dataset_pb2
from smu.parser import smu_utils_lib
import utilities

import unittest
from google.protobuf import text_format

def zero2() -> dataset_pb2.Geometry:
  """Return a Geometry with two points at the origin."""
  return text_format.Parse("""
        atom_positions: {
          x:0.0,
          y:0.0,
          z:0.0
        },
        atom_positions: {
          x:0.0,
          y:0.0,
          z:0.0
        }
      
""", dataset_pb2.Geometry())

class TestUtilities(unittest.TestCase):

  def test_zero_distance(self):
    coords = zero2()
    self.assertEqual(utilities.distance_between_atoms(coords, 0, 1), 0.0)

  def test_unit_x(self):
    coords = zero2()
    coords.atom_positions[1].x = 1.0 / smu_utils_lib.BOHR_TO_ANGSTROMS
    self.assertAlmostEqual(utilities.distance_between_atoms(coords, 0, 1), 1.0)

  def test_unit_y(self):
    coords = zero2()
    coords.atom_positions[1].y = 1.0 / smu_utils_lib.BOHR_TO_ANGSTROMS
    self.assertAlmostEqual(utilities.distance_between_atoms(coords, 0, 1), 1.0)

  def test_unit_z(self):
    coords = zero2()
    coords.atom_positions[1].z = 1.0 / smu_utils_lib.BOHR_TO_ANGSTROMS
    self.assertAlmostEqual(utilities.distance_between_atoms(coords, 0, 1), 1.0)

  def test_connected(self):
    pass

  @parameterized.expand(
  [
    ["[H]", 0, dataset_pb2.BondTopology.ATOM_H],
    ["C", 0, dataset_pb2.BondTopology.ATOM_C],
    ["N", 0, dataset_pb2.BondTopology.ATOM_N],
    ["[N+]", 1, dataset_pb2.BondTopology.ATOM_NPOS],
    ["O", 0, dataset_pb2.BondTopology.ATOM_O],
    ["[O-]", -1, dataset_pb2.BondTopology.ATOM_ONEG],
    ["F", 0, dataset_pb2.BondTopology.ATOM_F]
  ]
  )
  def test_molecule_to_bond_topology_geom(self, smiles, charge, expected):
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    bt,geom = utilities.molecule_to_bond_topology_geom(mol)
    self.assertEqual(len(bt.atoms), mol.GetNumAtoms())
    self.assertEqual(bt.atoms[0], expected)

  @parameterized.expand(
  [
    ["CC", dataset_pb2.BondTopology.BOND_SINGLE],
    ["C=C", dataset_pb2.BondTopology.BOND_DOUBLE],
    ["C#C", dataset_pb2.BondTopology.BOND_TRIPLE]
  ]
  )
  def test_bonds(self, smiles, expected):
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    bt,geom = utilities.molecule_to_bond_topology_geom(mol)
    self.assertEqual(len(bt.atoms), mol.GetNumAtoms())

if __name__ == "__main__":
  unittest.main()
