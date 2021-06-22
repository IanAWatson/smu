# utility functions needed for discerning topology from geometry

module Utilities

using Combinatorics

include("../jlout/dataset_pb.jl")

export bohr_to_angstroms, angstroms_to_bohr, distance_between_atoms, bonded, distances

function bohr_to_angstroms(distance)
  return distance * 0.529177249
end

function angstroms_to_bohr(distance)
  return distance / 0.529177249
end

function distance_between_atoms(geom, a1, a2)::Float32
  return bohr_to_angstroms(sqrt(
    (geom.atom_positions[a1].x - geom.atom_positions[a2].x) *
    (geom.atom_positions[a1].x - geom.atom_positions[a2].x) +
    (geom.atom_positions[a1].y - geom.atom_positions[a2].y) *
    (geom.atom_positions[a1].y - geom.atom_positions[a2].y) +
    (geom.atom_positions[a1].z - geom.atom_positions[a2].z) *
    (geom.atom_positions[a1].z - geom.atom_positions[a2].z)
  ))
end

"""Return an int array of the bonded atoms in `bond_topology`.
Args:
Returns:
  a numpy array of BondType's
"""
function bonded(bond_topology)::Array{Int32, 2}
  natoms = length(bond_topology.atoms)
  connected = zeros(Int32, (natoms, natoms))  # to be returned
  # Note need to convert to 1 indexing for Julia
  for bond in bond_topology.bonds
    a1 = bond.atom_a + 1
    a2 = bond.atom_b + 1
    connected[a1, a2] = connected[a2, a1] = bond.bond_type
  end
  return connected
end

"""Given a connection matrix `bonded`, return a vector of
vectors, one for each atom, containing the indices of the
bonded atoms.
Args:
  bonded:
Returns:
"""
function connections(bonded::Array{Int32, 2})::Vector{Vector{Int32}}
  natoms = size(bonded, 1)
  result = Vector{Vector{Int32}}(undef, (natoms))
  for i in 1:natoms
    nbrs = Vector{Int32}()
    for j in 1:natoms
      if bonded[i,j] > 0
        push!(nbrs, j)
      end
    end
    result[i] = nbrs
  end
  return result
end

"""Return a float array of the interatomic distances in `geometry`.
Args:
  geometry:
Returns:
  a numpy array of distances
"""
function distances(geometry::Geometry)::Array{Float, 2}
  natoms = length(geometry.atom_positions)
  result = zeros(Float32, (natoms, natoms))
  for atoms in combinations(1:natoms, 2)
    result[i, j] = result[j, i] = distance_between_atoms(geometry, atoms[1], atoms[2])
  end
  return result
end


"""Transform the bonds attribute of `bond_topology` to a canonical form.

Args:
  bond_topology:
Returns:
  BondTopology
"""
function canonical_bond_topology(bond_topology::BondTopology)
  if length(bond_topology.bonds) < 2
    return
  end

  for bond in bond_topology.bonds
    if bond.atom_a > bond.atom_b
      bond.atom_a, bond.atom_b = bond.atom_b, bond.atom_a
    end
  end

  sort!(bond_topology.bonds, b -> (b.atom_a, b.atom_b))
end

Base.:(==)(b1::BondTopology_Bond, b2::BondTopology_Bond)::Bool = b1.atom_a == b2.atom_a &&
                b1.atom_b == b2.atom_b && b1.bond_type == b2.bond_type

"""Return True if bt1 == bt2.
Note that there is no attempt to canonialise the protos.
Args:
Returns:
"""
function same_bond_topology(bt1::BondTopology, bt2::BondTopology)::Bool
  if length(bt2.atoms) != length(bt2.atoms)
    return false
  end

  if length(bt2.bonds) != length(bt2.bonds)
    return false
  end

  # https://stackoverflow.com/questions/47564825/check-if-all-the-elements-of-a-julia-array-are-equal
  # all(y->y==x[1], x);

  atoms1 = bt1.atoms
  if ! all(atoms1->a1==bt2.atoms[1], bt2.atoms)
    return false
  end

  bonds1 = bt1.bonds
  return all(bonds1->b1==bt2.bonds[1], bt2.bonds)
end

"""Recusrively visit nodes in the graph defined by `nbrs`.
Args:
  nbrs:
  atom:
  visited:
Returns:
  The number of nodes visited - including `atom`.
"""
function visit(nbrs::Vector, atom::Int, visited::Vector{Bool})::Int32
  visited[atom] = 1
  result = 1    # To be returned.
  for nbr in nbrs[atom]
    if visited[nbr] > 0
      continue
    end
    result += visit(nbrs, nbr, visited)
  end

  return result
end

"""Return True if `bond_topology` is a single fragment.
Args:
  bond_topology:
Returns:
  True if `bond_topology` is a single fragment.
"""
function is_single_fragment(bond_topology::BondTopology)::Bool
  natoms = length(bond_topology.atoms)
  nbonds = length(bond_topology.bonds)
  # Some special cases are easy.
  if natoms == 1
    return true
  end
  if natoms == 2 && nbonds == 1
    return true
  end
  if natoms == 3 && nbonds == 2
    return true
  end
  if natoms == nbonds && natoms <= 4
    return true
  end

  attached = connections(bonded(bond_topology))
  # Any atom with zero neighbors means a detached atom.
  if any(n->length(n) == 0, attached)
    return false
  end

  visited = zeros(Int32, natoms)
  # Mark anything with a single connection as visited.
  # Record the index of an atom that has multiple connections.
  a_multiply_connected_atom = -1
  for i in range(1, natoms)
    if bond_topology.atoms[i] == BondTopology.AtomType.ATOM_H
      visited[i] = 1
      continue
    end

    if length(attached[i]) > 1
      a_multiply_connected_atom = i
      continue
    end

    # A singly connected heavy atom. Mark visited if not of a two atom fragment.
    if length(attached[attached[i][0]]) > 1
      visited[i] = 1
    end
  end

  if a_multiply_connected_atom < 0     # Cannot happen
    return false
  end

  number_visited = sum(visited) + visit(attached, a_multiply_connected_atom, visited)
  return number_visited == natoms
end

end  # module Utilities
