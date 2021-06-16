import tensorflow as tf

from absl import app
from absl import flags
from absl import logging

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

from smu import dataset_pb2
import utilities

FLAGS = flags.FLAGS

flags.DEFINE_string("smiles", None, "Smiles input file")
flags.DEFINE_string("output", None, "TFRecord output file")

def contains_aromatic(mol: Chem.RWMol) -> bool:
  """
  """
  for atom in mol.GetAtoms():
    if atom.GetIsAromatic():
      return True
  return False

def smi23d(unused_argv):
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

    for n in range(len(hmols)):
      hmols[n].SetProp("_Name","%s"%sid[n])
      hmols[n].SetProp("_ID","%s"%sid[n])
      hmols[n].SetProp("_SMILES","%s"%smiles[n])
      if hmols[n].GetNumConformers() == 0:
        continue
      conf = hmols[n].GetConformer(0)
      natoms = hmols[n].GetNumAtoms()
      geom = dataset_pb2.Geometry()
      for i in range(0, natoms):
        atom = dataset_pb2.Geometry.AtomPos()
        atom.x = conf.GetAtomPosition(i).x
        atom.y = conf.GetAtomPosition(i).y
        atom.z = conf.GetAtomPosition(i).z
        geom.atom_positions.append(atom)

      conformer = dataset_pb2.Conformer()
      conformer.bond_topologies.append(utilities.molecule_to_bond_topology(hmols[n]))
      conformer.optimized_geometry.CopyFrom(geom)
      file_writer.write(conformer.SerializeToString())


if __name__ == "__main__":
  flags.mark_flag_as_required("smiles")
  app.run(smi23d)
