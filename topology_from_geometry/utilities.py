# Utility functions for SMU

import math

import numpy as np

from rdkit import Chem

from smu import dataset_pb2
from smu.parser import smu_utils_lib

BOHR2ANSTROM = 0.529177

DISTANCE_BINS = 10000

def distance_between_atoms(geom: dataset_pb2.Geometry,
                a1: int, a2: int) -> float:
  """Return the distance between atoms `a1` and `a2` in `geom`.
  Args:
    geom:
    a1:
    a2:
  Returns:
    Distance in Angstroms.
  """
  return smu_utils_lib.bohr_to_angstroms(math.sqrt((geom.atom_positions[a1].x - geom.atom_positions[a2].x) *
                             (geom.atom_positions[a1].x - geom.atom_positions[a2].x) +
                             (geom.atom_positions[a1].y - geom.atom_positions[a2].y) *
                             (geom.atom_positions[a1].y - geom.atom_positions[a2].y) +
                             (geom.atom_positions[a1].z - geom.atom_positions[a2].z) *
                             (geom.atom_positions[a1].z - geom.atom_positions[a2].z)))

def bonded(bond_topology: dataset_pb2.BondTopology) -> np.array:
  """Return an int array of the bonded atoms in `bond_topology`.
  Args:
  Returns:
    a numpy array of BondType's
  """
  natoms = len(bond_topology.atoms)
  connected = np.full((natoms, natoms), 0, dtype=np.int32)
  for bond in bond_topology.bonds:
    a1 = bond.atom_a
    a2 = bond.atom_b
    connected[a1, a2] = connected[a2, a1] = bond.bond_type
  return connected

def distances(geometry: dataset_pb2.Geometry) -> np.array:
  """Return a float array of the interatomic distances in `geometry`.
  Args:
    geometry:
  Returns:
    a numpy array of distances
  """
  natoms = len(geometry.atom_positions)
  distances = np.full((natoms, natoms), 0.0, dtype=np.float32)
  for i in range(0, natoms):
    for j in range(i+1, natoms):
      distances[i,j] = distances[j,i] = distance_between_atoms(geometry, i, j)
  return distances

def rdkit_atom_to_atom_type(atom: Chem.Atom) -> dataset_pb2.BondTopology.AtomType:
  """
    Args:
      atom: 
    Returns:
      AtpmType
  """
  if atom.GetAtomicNum() == 1:
    return dataset_pb2.BondTopology.ATOM_H
  if atom.GetAtomicNum() == 6:
    return dataset_pb2.BondTopology.ATOM_C
  if atom.GetAtomicNum() == 7:
    if atom.GetFormalCharge() == 0:
      return dataset_pb2.BondTopology.ATOM_N
    else:
      return dataset_pb2.BondTopology.ATOM_NPOS
  if atom.GetAtomicNum() == 8:
    if atom.GetFormalCharge() == 0:
      return dataset_pb2.BondTopology.ATOM_O
    else:
      return dataset_pb2.BondTopology.ATOM_ONEG
  if atom.GetAtomicNum() == 9:
     return dataset_pb2.BondTopology.ATOM_F

  raise ValueError(f"Unrecognized atom type {atom.GetAtomicNum()}")
def rdkit_bond_type_to_btype(bond_type: Chem.BondType) -> dataset_pb2.BondTopology.BondType:
  """
    Args:
    Returns:
  """
  if bond_type == Chem.rdchem.BondType.SINGLE:
    return dataset_pb2.BondTopology.BOND_SINGLE
  if bond_type == Chem.rdchem.BondType.DOUBLE:
    return dataset_pb2.BondTopology.BOND_DOUBLE
  if bond_type == Chem.rdchem.BondType.TRIPLE:
    return dataset_pb2.BondTopology.BOND_TRIPLE

  raise ValueError(f"Unrecognized bond type #{bond_type}")

def molecule_to_bond_topology(mol: Chem.RWMol) -> dataset_pb2.BondTopology:
  """
  """
  bond_topology = dataset_pb2.BondTopology()
  for atom in mol.GetAtoms():
    bond_topology.atoms.append(rdkit_atom_to_atom_type(atom))

  for bond in mol.GetBonds():
    btype = rdkit_bond_type_to_btype(bond.GetBondType())
    bt_bond = dataset_pb2.BondTopology.Bond()
    bt_bond.atom_a = bond.GetBeginAtom().GetIdx()
    bt_bond.atom_b = bond.GetEndAtom().GetIdx()
    bt_bond.bond_type = btype
    bond_topology.bonds.append(bt_bond)

  return bond_topology
