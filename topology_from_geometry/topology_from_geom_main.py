from absl import app
from absl import flags
from absl import logging

import apache_beam as beam
from apache_beam.options.pipeline_options import PipelineOptions
from apache_beam.io.tfrecordio import ReadFromTFRecord

from smu import dataset_pb2

import topology_from_geom

FLAGS = flags.FLAGS

flags.DEFINE_string("input", None, "TFDataRecord file containg Conformer protos")
flags.DEFINE_string("output", None, "Output file")

def ReadConFormer(input:str, output:str):
  """
  """
  class GetAtoms(beam.DoFn):
    def process(self, item):
      yield item.optimized_geometry.atom_positions[0].x
  with beam.Pipeline(options=PipelineOptions()) as p:
    protos = (p
      | beam.io.tfrecordio.ReadFromTFRecord(input, coder=beam.coders.ProtoCoder(dataset_pb2.Conformer().__class__))
      | beam.ParDo(topology_from_geom.TopologyFromGeom())
      | beam.io.textio.WriteToText(output)
    )

    return protos

def topology_from_geometry_main(unused_argv):
  del unused_argv

  protos = ReadConFormer(FLAGS.input, FLAGS.output)
  print(protos)



if __name__ == "__main__":
  flags.mark_flag_as_required("input")
  app.run(topology_from_geometry_main)
