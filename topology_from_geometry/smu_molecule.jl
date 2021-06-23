using ProtoBuf

using BondLengthDistributions
using SmuUtilities
using dataset_pb2

export MatchingParameters
export SmuMolecule
export set_initial_score_and_incrementer!
export generate_search_state
export place_bonds!

struct MatchingParameters
  must_match_all_bonds::Bool
end
MatchingParameters() = MatchingParameters(true)

mutable struct SmuMolecule
  starting_bond_topology::BondTopology
  natoms::Int32
  max_bonds::Vector{Int32}
  bonds_with_hydrogens_attached::Vector{Int32}
  current_bonds_attached::Vector{Int32}

  initial_score::Float32
  accumulate_score::Function

  bonds::Vector{Tuple{Int32, Int32}}
  scores::Vector{Float32}

  must_match_all_bonds::Bool

  SmuMolecule(hydrogens_attached::BondTopology, 
              bonds_to_scores::Dict{Tuple{Int32,Int32}, Vector{Float32}},
              matching_parameters::MatchingParameters) =
    (
    mol = new();
#   copy!(mol.starting_bond_topology, hydrogens_attached);
    mol.starting_bond_topology = hydrogens_attached;
    natoms = length(hydrogens_attached.atoms);
    mol.max_bonds = zeros(Int32, natoms);
    for i in 1:natoms
      mol.max_bonds[i] = SmuUtilities.smu_atom_type_to_max_con(hydrogens_attached.atoms[i])
    end;
    # The starting configuration
    mol.bonds_with_hydrogens_attached = zeros(Int32, natoms);
    for bond in hydrogens_attached.bonds
      mol.bonds_with_hydrogens_attached[bond.atom_a + 1] += 1
      mol.bonds_with_hydrogens_attached[bond.atom_b + 1] += 1
    end;
    mol.current_bonds_attached = zeros(Int32, natoms);

    # separately extract the bonds and scores
    mol.bonds = keys(bonds_to_scores);
    mol.scores = values(bonds_to_scores)
    )
end

function set_initial_score_and_incrementer!(score::T, op, smu_molecule::SmuMolecule) where{T<:Number}
  smu_molecule.initial_score = score
  smu_molecule.accumulate_score = op
end

function accumulate_score(mol::SmuMolecule, existing_score::T, increment::T)::Float32 where {T<:Number}
  return mol.op(existing_score, increment)
end

"""Possibly add a new bond to the current config.

If the bond can be placed, updates self._current_bonds_attached for
both `a`` and `a2`.
  Args:
    a1:
    a2:
    btype:
  Returns:
print(f"Trying to place bond {btype} current {self._current_bonds_attached[a1]} and {self._current_bonds_attached[a2]}")
"""
function _place_bond!(a1::Int32, a2::Int32, btype::Int32, mol::SmuMolecule)::Bool
  mol.current_bonds_attached[a1] + btype > mol.max_bonds[a1] && return false
  mol.current_bonds_attached[a2] + btype > mol.max_bonds[a2] && return false

  mol.current_bonds_attached[a1] += btype
  mol.current_bonds_attached[a2] += btype
  return true
end


"""For each pair of atoms, return a list of plausible bond types.

The result will be used to enumerate all the possible bonding forms.
The resulting Vector has the same length as mol.bonds, and each item
in that vector are the indices of the plausible bond types.
Args:
  mol:
Returns:
  Vector of Vectors - one for each atom pair.
"""
function generate_search_state(mol::SmuMolecule)::Vector{Vector{Int32}}
  result = Vector{Vector{Int32}}()
  for ndx in 1:length(mol.bonds)
    # For each pair of atoms, the plausible bond types - non zero score.
    scores = mol.scores[ndx]
    push!(result, findall(scores->scores > 0.0, scores))
  end

  return result
end

"""Place bonds corresponding to `state`.

Args:
  state: for each pair of atoms, the kind of bond to be placed.
Returns:
  If successful, the score.
"""
function place_bonds!(state::Vector{Int32},
                      mol::SmuMolecule)::Union{dataset_pb2.BondTopology, nothing}
  result = dataset_pb2.BondTopology()  # To be returned.
  copy!(result, mol.starting_bond_topology)  # Only Hydrogens attached.
  result.score = mol.initial_score

  # Initialize state in mol to Hydrogens attached.
  mol.current_bonds_attached = copy!(mol.bonds_with_hydrogens_attached)

  for i in 1:length(state)
    a1 = mol.bonds[i][0]
    a2 = mol.bonds[i][1]
    btype = state[i]
    _place_bond!(a1, a2, btype, mol) || return nothing

    result.score = accumulate_score(mol, result.score, mol.scores[i][btype])
    # If the bond is anything other than BOND_UNKNOWN, add it to result.
    btype == BondTopology_BondType.BOND_UNKNOWN || add_bond!(a1, a2, btype, result)
  end

  # Optionally check whether all bonds have been matched
  mol.must_match_all_bonds || return result

  mol.current_bonds_attached == mol.max_bonds || return nothing

  return result
end
