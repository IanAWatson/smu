# Discern connectivity from the coordinates in a Geometry

using ArgMacros
using Logging
using ProtoBuf
using TFRecord

append!(LOAD_PATH, ["smu/jlout/", "smu/topology_from_geometry/", "smu/geometry"])

using BondLengthDistributions
using dataset_pb2
using TopologyFromGeometry

function get_topology_from_geometry(bond_length_distributions::AllBondLengthDistributions,
                bond_topology::BondTopology,
                geometry::Geometry,
                output_stream::IOStream)::Bool
  matching_parameters = MatchingParameters()
  result = bond_topologies_from_geom(bond_length_distributions, bond_topology,
                                     geometry, matching_parameters)
  hasproperty(result, :bond_topology) || return false
  println(output_stream, length(result.bond_topology))
  return true
end

function main()
  @inlinearguments begin
    @argumentrequired String input_fname "-i" "--input"
    @argumentrequired String output_fname "-o" "--output"
    @argumentrequired String bonds "-b" "--bonds"
    @argumentflag debug "-debug" "--debug"
    @argumentflag verbose "-v" "--verbose"
  end
  println("input $(input_fname) bonds $(bonds) output $(output_fname)")
  debug && Logging.configure(level=Logging.Debug)

  bond_length_distributions = AllBondLengthDistributions()

  add_from_files!(bonds, bond_length_distributions) || @error("Cannot build bond length distribution $bonds")

  molecules_read = 0
  molecules_processed = 0

  open(output_fname, "w") do output_stream
    for conformer in TFRecord.read([input_fname], record_type=Conformer)
      molecules_read += 1
#     println(conformer)
      if ! hasproperty(conformer, :optimized_geometry)
        @warn("No optimized geometry, skipping")
        continue
      end

      if ! hasproperty(conformer, :bond_topologies)
        @warn("No bond_topology, skipping")
        continue
      end
      if get_topology_from_geometry(bond_length_distributions, conformer.bond_topologies[1],
                                    conformer.optimized_geometry, output_stream) 
        molecules_processed += 1
      end
      molecules_read > 10 && break
    end
  end
  verbose && println("Read $(molecules_read) molecules, processed $(molecules_processed)")
  0
end


main()
