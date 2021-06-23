"""Testing the SmuMolecule"""

module TestSmuMolecule

using ProtoBuf
using Test

using SmuMolecules
using SmuUtilities

using dataset_pb2

export all_tests

"""Make a BondTopology consisting of `n` disconnected Carbon atoms.
"""
function carbon_atoms(n::T)::BondTopology where {T<:Integer}
  result = BondTopology()
  for i in 1:n
    add_atom_by_atomic_number!(result, 6)
  end
  return result
end

"""The simplest multi-atom molecule, CC"""
function test_ethane()
  matching_parameters = MatchingParameters()
  ethane = carbon_atoms(2)
  scores = [0.1, 1.1, 2.1, 3.1]
  bonds_to_scores = Dict{Tuple{Int32, Int32}, Vector{Float32}}((0, 1) => scores)
  mol = SmuMolecule(ethane, bonds_to_scores, matching_parameters)
  state = generate_search_state(mol)
  @test length(state) == 1
  @test state == [0, 1, 2, 3]

  for s in product(state...)
    res = place_bonds!(s, mol)
    @test res !== nothing
    # self.assertAlmostEqual(res.score, scores[i])
  end
end

function test_ethane_all_btypes()
  matching_parameters = MatchingParameters()

  test_cases = [
    [0, dataset_pb2.BondTopology.BOND_UNDEFINED],
    [1, dataset_pb2.BondTopology.BOND_SINGLE],
    [2, dataset_pb2.BondTopology.BOND_DOUBLE],
    [3, dataset_pb2.BondTopology.BOND_TRIPLE],
  ]

  for t in test_cases
    numeric_btype = t[1]
    smu_btype = t[2]
    cc = carbon_atoms(2)
    bonds_to_scores = Dict{Tuple{Int32,Int32}}((0, 1) => zeros(Float32, 4))
    bonds_to_scores[(0, 1)][numeric_btype] = 1.0
    mol = SmuMolecule(cc, bonds_to_scores, matching_parameters)
    state = mol.generate_search_state()
    for s in Iterators.product(state...)
      res = place_bonds!(s, mol)
      @test res !== nothing
      if numeric_btype == 0
        @test length(res.bonds) == 0
      else
        @test length(res.bonds) == 1
        @test res.bonds[1].bond_type == smu_btype
      end
    end
  end
  true
end

function test_propane_all()
  matching_parameters = MatchingParameters()

  test_cases = [
    [0, 0, 0, 2.0],
    [0, 1, 1, 2.0],
    [0, 2, 1, 2.0],
    [0, 3, 1, 2.0],
    [1, 1, 2, 2.0],
    [1, 2, 2, 2.0],
    [1, 3, 2, 2.0],
    [2, 2, 2, 2.0],
    [2, 3, 0, nothing],
    [3, 3, 0, nothing]
  ]
  for t in test_cases
    btype1 = t[1]
    btype2 = t[2]
    expected_bonds = t[3]
    expected_score = t[4]
    cc = carbon_atoms(3)
    bonds_to_scores = Dict{Tuple{Int32,Int32}}((0, 1) => zeros(Float32, 4),
                       (1, 2) => zeros(Float32, 4))
    bonds_to_scores[(0, 1)][btype1] = 1.0
    bonds_to_scores[(1, 2)][btype2] = 1.0
    mol = smu_molecule.SmuMolecule(cc, bonds_to_scores, matching_parameters)
    state = mol.generate_search_state()
    for s in itertools.product(state...)
      res = place_bonds!(s, mol)
      if expected_score !== nothing
        @test res !== nothing
        @test length(res.bonds) == expected_bonds
        @test isapprox(res.score, expected_score)
        if btype1 == 0
          if btype2 > 0
            @test res.bonds[1].bond_type == btype2
          end
        else
          @test res.bonds[1].bond_type == btype1
          @test res.bonds[2].bond_type == btype2
        end
      else
        @test res === nothing
      end
    end
  end
  true
end

function test_operators()
  matching_parameters = MatchingParameters()
  cc = carbon_atoms(3)
  bonds_to_scores::Dict{Tuple{Int32,Int32}}((0, 1) => zeros(Float32, 4),
                                            (1, 2) => zeros(Float32, 4))
  scores = [1.0, 3.0]
  bonds_to_scores[(0, 1)][1] = scores[1]
  bonds_to_scores[(1, 2)][1] = scores[2]
  mol = smu_molecule.SmuMolecule(cc, bonds_to_scores, matching_parameters)
  set_initial_score_and_incrementer!(1.0, Base.:*, mol)
  state = generate_search_state(mol)
  for s in Iterators.product(state...)
    res = place_bonds!(s, mol)
    @test isapprox(res.score, cumprod(scores))
  end
  true
end


function all_tests()
  @test test_ethane()
  @test test_ethane_all_btypes()
  @test test_propane_all()
  @test test_operators()
end



end  # module TestSmuMolecule
