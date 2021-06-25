"""Conversion from smiles to 3D TFDatarecord."""
import tensorflow as tf

from absl import app
from absl import flags
from absl import logging

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

import utilities
from smu import dataset_pb2
from smu.parser import smu_utils_lib

FLAGS = flags.FLAGS

flags.DEFINE_string("smiles", None, "Smiles input file")
flags.DEFINE_string("output", None, "TFRecord output file")

def contains_aromatic(mol: Chem.RWMol) -> bool:
  """Returns True of `mol` contains any aromatic atoms."""
  for atom in mol.GetAtoms():
    if atom.GetIsAromatic():
      return True
  return False

def smi23d(unused_argv):
  """Converts a smiles file to 3D TFDatarecord proto Conformer"""
  del unused_argv

  df = pd.read_csv(FLAGS.smiles)
  mols = [Chem.MolFromSmiles(smi) for smi in df.SMILES]
  mols = [m for m in mols if not contains_aromatic(m)]
  logging.info("Read %d molecules", len(mols))
  hmols = [Chem.AddHs(m) for m in mols]
  for mol  in hmols:
    AllChem.EmbedMolecule(mol,AllChem.ETKDG())
  smiles = list(df.SMILES)
  sid = list(df.ID)

  with tf.io.TFRecordWriter(FLAGS.output) as file_writer:

    for n, hmol in enumerate(hmols):
      hmol.SetProp("_Name","%s"%sid[n])
      hmol.SetProp("_ID","%s"%sid[n])
      hmol.SetProp("_SMILES","%s"%smiles[n])
      if hmol.GetNumConformers() == 0:
        continue
      conf = hmol.GetConformer(0)
      natoms = hmol.GetNumAtoms()
      geom = dataset_pb2.Geometry()
      for i in range(0, natoms):
        atom = dataset_pb2.Geometry.AtomPos()
        atom.x = conf.GetAtomPosition(i).x / smu_utils_lib.BOHR_TO_ANGSTROMS
        atom.y = conf.GetAtomPosition(i).y / smu_utils_lib.BOHR_TO_ANGSTROMS
        atom.z = conf.GetAtomPosition(i).z / smu_utils_lib.BOHR_TO_ANGSTROMS
        geom.atom_positions.append(atom)

      conformer = dataset_pb2.Conformer()
      conformer.bond_topologies.append(utilities.molecule_to_bond_topology(hmol))
      conformer.bond_topologies[-1].smiles = smiles[n]
      conformer.optimized_geometry.CopyFrom(geom)
      file_writer.write(conformer.SerializeToString())


if __name__ == "__main__":
  flags.mark_flag_as_required("smiles")
  app.run(smi23d)
