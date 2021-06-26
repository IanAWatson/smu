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
  number_bond_topologies = length(result.bond_topology)
  to_write = string(number_bond_topologies)
  for bt in result.bond_topology
    to_write *= ",$(round(bt.score, digits=3))"
    number_bond_topologies == 1 && break
    to_write *= bt.is_starting_topology ? ",t" : ",f"
  end
  println(output_stream, to_write)
  true
end

function main()
  @inlinearguments begin
    @argumentrequired String input_fname "-i" "--input"
    @argumentrequired String output_fname "-o" "--output"
    @argumentrequired String bonds "-b" "--bonds"
    @argumentoptional Int64 nprocess "-N" "--nprocess"
    @argumentflag exclude_non_bonded "-x" "--xcldnonbond"
    @argumentflag debug "-debug" "--debug"
    @argumentflag verbose "-v" "--verbose"
  end
  if debug
    logger=Logging.SimpleLogger(stderr,Logging.Debug)
    global_logger(logger)
  end
# global_logger(Logging.SimpleLogger(stderr, Logging.Error))

  @info("input $(input_fname) bonds $(bonds) output $(output_fname) $nprocess")
  flush(stdout)
  nprocess === nothing && (nprocess = typemax(Int64))

  bond_length_distributions = AllBondLengthDistributions()
  exclude_non_bonded !== nothing && (bond_length_distributions.include_non_bonded = false)

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
      @info("processing $(conformer.conformer_id)")
      if get_topology_from_geometry(bond_length_distributions, conformer.bond_topologies[1],
                                    conformer.optimized_geometry, output_stream) 
        molecules_processed += 1
      end
      molecules_read > nprocess && break
    end
  end
  verbose && println("Read $(molecules_read) molecules, processed $(molecules_processed)")
  0
end


main()
