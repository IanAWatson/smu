# Discern connectivity from the coordinates in a Geometry

using ArgMacros
using ProtoBuf
using TFRecord

append!(LOAD_PATH, ["smu/jlout/", "smu/topology_from_geometry/", "smu/geometry"])

using dataset_pb2
using TopologyFromGeometry

function main()
@inlinearguments begin
    @argumentrequired String input "-i" "--input"
    @argumentrequired String output "-o" "--output"
    @argumentrequired String bonds "-b" "--bonds"
    @argumentflag verbose "-v"
  end
  println("input $(input) bonds $(bonds) output $(output)")

  for rawdata in TFRecord.read(input)
    conf = dataset_pb2.Conformer()
    ProtoBuf.readproto(rawdata, conf)
    println(conf)
  end
end


main()
