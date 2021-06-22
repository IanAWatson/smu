module UtilitiesTest

export test_distance_between

using Test

include("../jlout/dataset_pb.jl")

include("utilities.jl")

function test_bohr_angstroms()
  @test isapprox(Utilities.angstroms_to_bohr(Utilities.bohr_to_angstroms(1.0)), 1.0)
end

function test_distance_between()
  geom = Geometry()
  atoms = Vector{Geometry_AtomPos}()
  a1 = Geometry_AtomPos(x=0.0, y=0.0, z=0.0)
  push!(atoms, a1)
  a2 = Geometry_AtomPos(x=0.0, y=0.0, z=Utilities.angstroms_to_bohr(1.0))
  push!(atoms, a2)
  setproperty!(geom, :atom_positions, atoms)
  @test isapprox(Utilities.distance_between_atoms(geom, 1, 2), 1.0)
  setproperty!(geom.atom_positions[2], :x, Utilities.angstroms_to_bohr(1.0))
  setproperty!(geom.atom_positions[2], :y, Utilities.angstroms_to_bohr(1.0))
  setproperty!(geom.atom_positions[2], :z, Utilities.angstroms_to_bohr(1.0))
  @test isapprox(Utilities.distance_between_atoms(geom, 1, 2), sqrt(3.0))
end

function test_bonded()
  bond_topology = BondTopology()
  atoms = Vector{Int32}()
  bonds = Vector{BondTopology_Bond}()
  push!(atoms, BondTopology_AtomType.ATOM_C)
  push!(atoms, BondTopology_AtomType.ATOM_C)
  push!(atoms, BondTopology_AtomType.ATOM_C)
  push!(bonds, BondTopology_Bond(atom_a=0, atom_b=1, bond_type=BondTopology_BondType.BOND_SINGLE))
  push!(bonds, BondTopology_Bond(atom_a=1, atom_b=2, bond_type=BondTopology_BondType.BOND_SINGLE))
  setproperty!(bond_topology, :atoms, atoms)
  setproperty!(bond_topology, :bonds, bonds)

  bonded = Utilities.bonded(bond_topology)
  @test length(atoms) == size(bonded, 1)
  @test bonded[1,2] > 0
  @test bonded[2,3] > 0
  @test count(!iszero, bonded) == 4  # Bonds are placed twice.
  push!(bond_topology.bonds, BondTopology_Bond(atom_a=0, atom_b=2, bond_type=BondTopology_BondType.BOND_DOUBLE))
  bonded = Utilities.bonded(bond_topology)
  @test count(!iszero, bonded) == 6  # Bonds are placed twice.
  @test bonded[1,2] == bonded[2,3]
  @test bonded[1,2] != bonded[3,1]
end

function test_connections()
  bond_topology = BondTopology()
  atoms = Vector{Int32}()
  bonds = Vector{BondTopology_Bond}()
  push!(atoms, BondTopology_AtomType.ATOM_C)
  push!(atoms, BondTopology_AtomType.ATOM_C)
  push!(atoms, BondTopology_AtomType.ATOM_C)
  push!(bonds, BondTopology_Bond(atom_a=0, atom_b=1, bond_type=BondTopology_BondType.BOND_SINGLE))
  push!(bonds, BondTopology_Bond(atom_a=1, atom_b=2, bond_type=BondTopology_BondType.BOND_SINGLE))
  setproperty!(bond_topology, :atoms, atoms)
  setproperty!(bond_topology, :bonds, bonds)

  connections = Utilities.connections(Utilities.bonded(bond_topology))
  @test length(connections) == length(atoms)
  @test length(connections[1]) == 1
  @test length(connections[2]) == 2
  @test length(connections[3]) == 1
  @test connections[1][1] == 2
  @test connections[2][1] == 1
  @test connections[2][2] == 3
  @test connections[3][1] == 2
end

end  # module UtilitiesTest
