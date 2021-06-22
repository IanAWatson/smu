module UtilitiesTest

export test_distance_between

using ProtoBuf
using Test

using dataset_pb2

using TopologyFromGeom

function test_bohr_angstroms()
  @test isapprox(TopologyFromGeom.angstroms_to_bohr(TopologyFromGeom.bohr_to_angstroms(1.0)), 1.0)
end

function test_add_atom()
  bond_topology = BondTopology()
  add_atom!(bond_topology, 6)
  @test length(bond_topology.atoms) == 1
  @test bond_topology.atoms[1] == BondTopology_AtomType.ATOM_C

  add_atom!(bond_topology, 7)
  @test length(bond_topology.atoms) == 2
  @test bond_topology.atoms[2] == BondTopology_AtomType.ATOM_N

  add_atom!(bond_topology, 7, 1)
  @test length(bond_topology.atoms) == 3
  @test bond_topology.atoms[3] == BondTopology_AtomType.ATOM_NPOS

  # Do it systematically
  atype = [[1, 0], [6, 0], [7, 0], [7, 1], [8, 0], [8, 1], [9, 0]]
  expected = [
    BondTopology_AtomType.ATOM_H,
    BondTopology_AtomType.ATOM_C,
    BondTopology_AtomType.ATOM_N,
    BondTopology_AtomType.ATOM_NPOS,
    BondTopology_AtomType.ATOM_O,
    BondTopology_AtomType.ATOM_ONEG,
    BondTopology_AtomType.ATOM_F,
  ]
  @test length(atype) == length(expected)

  bt = BondTopology()
  for i in zip(atype, expected)
    add_atom!(bt, i[1][1], i[1][2])
    @test bt.atoms[end] == i[2]
  end
  @test length(bt.atoms) == length(expected)
  true
end

function test_add_single_bond()
  bond_topology = BondTopology()
  atoms = Vector{Int32}()
  push!(atoms, BondTopology_AtomType.ATOM_C)
  push!(atoms, BondTopology_AtomType.ATOM_C)
  setproperty!(bond_topology, :atoms, atoms)

  add_bond!(0, 1, 1, bond_topology)
  @test length(bond_topology.bonds) == 1
  @test bond_topology.bonds[1].bond_type == BondTopology_BondType.BOND_SINGLE
  true
end

function test_add_double_bond()
  bond_topology = BondTopology()
  atoms = Vector{Int32}()
  push!(atoms, BondTopology_AtomType.ATOM_C)
  push!(atoms, BondTopology_AtomType.ATOM_C)
  setproperty!(bond_topology, :atoms, atoms)

  add_bond!(0, 1, 2, bond_topology)
  @test length(bond_topology.bonds) == 1
  @test bond_topology.bonds[1].bond_type == BondTopology_BondType.BOND_DOUBLE
  true
end

function test_add_triple_bond()
  bond_topology = BondTopology()
  atoms = Vector{Int32}()
  push!(atoms, BondTopology_AtomType.ATOM_C)
  push!(atoms, BondTopology_AtomType.ATOM_C)
  setproperty!(bond_topology, :atoms, atoms)

  add_bond!(0, 1, 3, bond_topology)
  @test length(bond_topology.bonds) == 1
  @test bond_topology.bonds[1].bond_type == BondTopology_BondType.BOND_TRIPLE
  true
end

function test_add_multiple_bonds()
  bond_topology = BondTopology()
  atoms = Vector{Int32}()
  push!(atoms, BondTopology_AtomType.ATOM_C)
  push!(atoms, BondTopology_AtomType.ATOM_C)
  push!(atoms, BondTopology_AtomType.ATOM_C)
  setproperty!(bond_topology, :atoms, atoms)

  add_bond!(0, 1, 1, bond_topology)
  add_bond!(1, 2, 1, bond_topology)
  @test length(bond_topology.bonds) == 2
  @test bond_topology.bonds[1].bond_type == BondTopology_BondType.BOND_SINGLE
  @test bond_topology.bonds[2].bond_type == BondTopology_BondType.BOND_SINGLE
  true
end

function test_distance_between()
  geom = Geometry()
  atoms = Vector{Geometry_AtomPos}()
  a1 = Geometry_AtomPos(x=0.0, y=0.0, z=0.0)
  push!(atoms, a1)
  a2 = Geometry_AtomPos(x=0.0, y=0.0, z=TopologyFromGeom.angstroms_to_bohr(1.0))
  push!(atoms, a2)
  setproperty!(geom, :atom_positions, atoms)
  @test isapprox(TopologyFromGeom.distance_between_atoms(geom, 1, 2), 1.0)
  setproperty!(geom.atom_positions[2], :x, TopologyFromGeom.angstroms_to_bohr(1.0))
  setproperty!(geom.atom_positions[2], :y, TopologyFromGeom.angstroms_to_bohr(1.0))
  setproperty!(geom.atom_positions[2], :z, TopologyFromGeom.angstroms_to_bohr(1.0))
  @test isapprox(TopologyFromGeom.distance_between_atoms(geom, 1, 2), sqrt(3.0))

  distances = TopologyFromGeom.distances(geom)
  @test distances[1,1] == 0.0
  @test distances[2,2] == 0.0
  @test isapprox(distances[1,2], sqrt(3.0))
  @test isapprox(distances[2,1], sqrt(3.0))
  true
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

  bonded = TopologyFromGeom.bonded(bond_topology)
  @test length(atoms) == size(bonded, 1)
  @test bonded[1,2] > 0
  @test bonded[2,3] > 0
  @test count(!iszero, bonded) == 4  # Bonds are placed twice.
  push!(bond_topology.bonds, BondTopology_Bond(atom_a=0, atom_b=2, bond_type=BondTopology_BondType.BOND_DOUBLE))
  bonded = TopologyFromGeom.bonded(bond_topology)
  @test count(!iszero, bonded) == 6  # Bonds are placed twice.
  @test bonded[1,2] == bonded[2,3]
  @test bonded[1,2] != bonded[3,1]
  true
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

  connections = TopologyFromGeom.connections(TopologyFromGeom.bonded(bond_topology))
  @test length(connections) == length(atoms)
  @test length(connections[1]) == 1
  @test length(connections[2]) == 2
  @test length(connections[3]) == 1
  @test connections[1][1] == 2
  @test connections[2][1] == 1
  @test connections[2][2] == 3
  @test connections[3][1] == 2
  true
end

function test_canonical_bond_topology()
  bond_topology = BondTopology()
  atoms = Vector{Int32}()
  bonds = Vector{BondTopology_Bond}()
  push!(atoms, BondTopology_AtomType.ATOM_C)
  push!(atoms, BondTopology_AtomType.ATOM_C)
  push!(atoms, BondTopology_AtomType.ATOM_C)
  push!(bonds, BondTopology_Bond(atom_a=2, atom_b=1, bond_type=BondTopology_BondType.BOND_SINGLE))
  push!(bonds, BondTopology_Bond(atom_a=1, atom_b=0, bond_type=BondTopology_BondType.BOND_SINGLE))
  setproperty!(bond_topology, :atoms, atoms)
  setproperty!(bond_topology, :bonds, bonds)

  TopologyFromGeom.canonical_bond_topology(bond_topology)

  expected = BondTopology()
  atoms = Vector{Int32}()
  bonds = Vector{BondTopology_Bond}()
  push!(atoms, BondTopology_AtomType.ATOM_C)
  push!(atoms, BondTopology_AtomType.ATOM_C)
  push!(atoms, BondTopology_AtomType.ATOM_C)
  push!(bonds, BondTopology_Bond(atom_a=0, atom_b=1, bond_type=BondTopology_BondType.BOND_SINGLE))
  push!(bonds, BondTopology_Bond(atom_a=1, atom_b=2, bond_type=BondTopology_BondType.BOND_SINGLE))
  setproperty!(expected, :atoms, atoms)
  setproperty!(expected, :bonds, bonds)
  @test bond_topology == expected
  @test TopologyFromGeom.same_bond_topology(bond_topology, expected)

  true
end

function test_single_fragment1()
  bond_topology = BondTopology()
  atoms = Vector{Int32}()
  bonds = Vector{BondTopology_Bond}()
  push!(atoms, BondTopology_AtomType.ATOM_C)  # atom 0
  push!(atoms, BondTopology_AtomType.ATOM_C)  # atom 1
  push!(atoms, BondTopology_AtomType.ATOM_C)  # atom 2
  push!(bonds, BondTopology_Bond(atom_a=0, atom_b=1, bond_type=BondTopology_BondType.BOND_SINGLE))
  setproperty!(bond_topology, :atoms, atoms)
  setproperty!(bond_topology, :bonds, bonds)

  @test !TopologyFromGeom.is_single_fragment(bond_topology)

  push!(bond_topology.bonds, BondTopology_Bond(atom_a=1, atom_b=2, bond_type=BondTopology_BondType.BOND_SINGLE))
  @test TopologyFromGeom.is_single_fragment(bond_topology)
  
  push!(bond_topology.atoms, BondTopology_AtomType.ATOM_C)  # atom 3
  @test !TopologyFromGeom.is_single_fragment(bond_topology)

  push!(bond_topology.atoms, BondTopology_AtomType.ATOM_C)  # atom 4
  @test !TopologyFromGeom.is_single_fragment(bond_topology)
  push!(bond_topology.bonds, BondTopology_Bond(atom_a=3, atom_b=4, bond_type=BondTopology_BondType.BOND_SINGLE))
  @test !TopologyFromGeom.is_single_fragment(bond_topology)

  push!(bond_topology.bonds, BondTopology_Bond(atom_a=2, atom_b=3, bond_type=BondTopology_BondType.BOND_SINGLE))
  @test TopologyFromGeom.is_single_fragment(bond_topology)
  push!(bond_topology.bonds, BondTopology_Bond(atom_a=0, atom_b=4, bond_type=BondTopology_BondType.BOND_SINGLE))
  @test TopologyFromGeom.is_single_fragment(bond_topology)

  true
end

function all_tests()
  @test test_add_atom()
  @test test_add_single_bond()
  @test test_add_double_bond()
  @test test_add_triple_bond()
  @test test_distance_between()
  @test test_bonded()
  @test test_connections()
  @test test_canonical_bond_topology()
  @test test_single_fragment1()
  print("All tests complete\n")
end


end  # module UtilitiesTest
