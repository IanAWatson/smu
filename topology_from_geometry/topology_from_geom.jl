
using dataset_pb2

using BondLengthDistributions
using SmuMolecules
using SmuUtilities

"""Generate a BondTopology that joins each Hydrogen atom to its nearest
    heavy atom.
Args:
  bond_topology:
  distances:
Returns:
"""
function hydrogen_to_nearest_atom(bond_topology::dataset_pb2.BondTopology,
                             distances::Array{Float32,2})::Union{dataset_pb2.BondTopology, Nothing}
  result = dataset_pb2.BondTopology()
  # Not sure if this is safe or not.
  setproperty!(result, :atoms, bond_topology.atoms)

  single_bond = 1

  natoms = length(bond_topology.atoms)
  for a1 in 1:natoms
    bond_topology.atoms[a1] == dataset_pb2.BondTopology_AtomType.ATOM_H || continue

    shortest_distance = 1.0e+30
    closest_heavy_atom = -1
    for a2 in 1:natoms
      bond_topology.atoms[a2] == dataset_pb2.BondTopology.AtomType.ATOM_H && continue

      distances[a1, a2] >= THRESHOLD && continue

      if distances[a1, a2] < shortest_distance
        shortest_distance = distances[a1, a2]
        closest_heavy_atom = a2
      end
    end

    closest_heavy_atom < 0 && return nothing

    add_bond!(a1, closest_heavy_atom, single_bond, result)
  end

  return result
end

"""Return the indices of the heavy atoms in `bond_topology`.
  Args:
  Returns:
"""
function indices_of_heavy_atoms(bond_topology::dataset_pb2.BondTopology)::Vector{Int32}
  atoms = bond_topology.atoms
  findall(atoms->atoms!= dataset_pb2.BondTopology_AtomType.ATOM_H, atoms)
end


"""Return all BondTopology's that are plausible.

  Given a molecule described by `bond_topology` and `geometry`, return all possible
  BondTopology that are consistent with that.
  Note that `bond_topology` will be put in a canonical form.
  Args:
    bond_length_distribution:
    bond_topology:
    geometry:
  Returns:
    TopologyMatches
"""
function bond_topologies_from_geom(
    bond_lengths::AllBondLengthDistributions,
    bond_topology::dataset_pb2.BondTopology,
    geometry::dataset_pb2.Geometry,
    matching_parameters::MatchingParameters)::dataset_pb2.TopologyMatches
  result = dataset_pb2.TopologyMatches()    # To be returned.
  len(bond_topology.atoms) == 1 && return result  # return an empty result

  canonical_bond_topology(bond_topology)
  distances = distances(geometry)

  # First join each Hydrogen to its nearest heavy atom, thereby
  # creating a starting BondTopology from which all others can grow
  starting_bond_topology = hydrogen_to_nearest_atom(bond_topology, distances)
  starting_bond_topology === nothing && return result

  atoms = bond_topology.atoms
  heavy_atom_indices = findall(atoms->atoms != dataset_pb2.BondTopology_AtomType.ATOM_H)

  # For each atom pair, a list of possible bond types.
  # Key is a tuple of the two atom numbers, value is a Vector
  # with the score for each bond type.

  bonds_to_scores = Dict{Tuple{Int32,Int32}, Vector{Float32}}()

  found_topologies = Vector{BondTopology}()

  for c in combinations(heavy_atom_indices, 2)  # All pairs
     i = c[1]
     j = c[2]
     dist = distances[i, j]
     dist > THRESHOLD && continue

    btypes = [bond_lengths.pdf_length_given_type(bond_topology.atoms[i],
                                                 bond_topology.atoms[j], btype, dist) for btype in 0:4]

#   Should work does not compile.
#   any(nonzero, btypes) && bonds_to_scores[(i, j)] = btypes
    if any(nonzero, btypes)
      bonds_to_scores[(i, j)] = btypes
    end
  end

  isempty(bonds_to_scores) && return result

  mol = smu_molecule.SmuMolecule(starting_bond_topology, bonds_to_scores, matching_parameters)

  search_space = mol.generate_search_state()
  for s in product(search_space...)
    bt = mol.place_bonds(s...)
    bt || continue

    canonical_bond_topology(bt)
#   should work, fails to compile...
#   same_bond_topology(bond_topology, bt) && bt.is_starting_topology = true
    if same_bond_topology(bond_topology, bt)
      bt.is_starting_topology = true
    end
    # smiles not set
    push!(found_topologies, bt)
  end

  if length(found_topologies) > 1
    sort!(found_topologies, by = b -> b.score, reverse=true)
  end

  setproperty!(result, :bond_topology, found_topologies)

  return result
end
