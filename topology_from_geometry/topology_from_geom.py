import apache_beam as beam

from smu import dataset_pb2

class TopologyFromGeom(beam.DoFn):
  """Beam class for extracting BondTopology from Conformer protos."""
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

    natoms = len(bond_topology.atoms)
    if natoms == 1:
      return geometry.atom_positions[0].x
    return geometry.atom_positions[0].x

