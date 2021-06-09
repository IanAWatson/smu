import apache_beam as beam
from apache_beam.options.pipeline_options import PipelineOptions
from apache_beam.io.tfrecordio import ReadFromTFRecord

from smu import dataset_pb2

from absl import app
from absl import flags
from absl import logging

from smu import dataset_pb2

import bond_lengths

FLAGS = flags.FLAGS

flags.DEFINE_string("input", None, "TFDataRecord file containg Conformer protos")
flags.DEFINE_string("output", None, "Output file")

class BondDistToString(beam.DoFn):
  def process(self, bond_dist):
    key, value = bond_dist
    print(f"BondDistToString: key {key} value {value}")
    yield f"{key[0]}.{key[1]}.{key[2]}.{key[3]}.{value}"
#   yield f"{key[0]}.{key[1]}.{key[2]}.{value[0]} {value[1]}"

class GroupBondTypes(beam.DoFn):
  def process(self, bond_dist):
    key, value = bond_dist
    print(f"GroupBondTypes: key #{key} value {value}")
    yield (key[0], key[1], key[2]), (key[3], value)

def get_bond_length_distribution_inner(input: str, output: str):
  """
  """
  print("Reading from {input} output to {output}")
  with beam.Pipeline(options=PipelineOptions()) as p:
    protos = (p
      | beam.io.tfrecordio.ReadFromTFRecord(input, coder=beam.coders.ProtoCoder(dataset_pb2.Conformer().__class__))
      | beam.ParDo(bond_lengths.GetBondLengthDistribution())
      | beam.CombinePerKey(sum)
#     | beam.ParDo(GroupBondTypes())
#     | beam.GroupByKey()
      | beam.ParDo(BondDistToString())
      | beam.io.WriteToText(output)
    )
    print(protos)

def get_bond_length_distribution(unused_argv):
  """Scan Conformer protos to extract bond length distributions."""
  del unused_argv

  get_bond_length_distribution_inner(FLAGS.input, FLAGS.output)


if __name__ == "__main__":
  app.run(get_bond_length_distribution)
