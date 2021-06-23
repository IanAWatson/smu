# Implementation empirical bond length distributions, supporting pdf() functions.

using Combinatorics
using CSV
using DataFrames
using LinearAlgebra

export BondLengthDistribution, AllBondLengthDistributions
export from_file, from_arrays, pdf

"""A single bond length distribution, holding a relationship between
  distance and value.
  It is assumed that the X values are uniformly distributed.
"""
mutable struct BondLengthDistribution
    bucket_size::Float32  # Distance between X values.
    min_distance::Float32
    max_distance::Float32
    pdf::Vector{Float32}
end

BondLengthDistribution() = BondLengthDistribution(0.0, 0.0, 0.0, [])

function from_file(fname::AbstractString, distribution::BondLengthDistribution)::Bool
  df = CSV.File(fname) |> DataFrame
  return from_arrays(df[!,1], df[!,2], distribution)
end

function from_arrays(lengths::Vector, counts::Vector, distribution::BondLengthDistribution)::Bool
    min_diff, max_diff = extrema(diff(lengths))
    if abs(max_diff - min_diff) > 1.0e-05
      @warn("non uniform buckets $(min_diff) $(max_diff)")
      return false
    end
    distribution.bucket_size = (max_diff + min_diff) / 2.0  # For accuracy.
    distribution.min_distance = lengths[1]
    distribution.max_distance =  lengths[end] + distribution.bucket_size
    distribution.pdf = Float32.(counts)
    normalize!(distribution.pdf, 1)
    return true
end

function pdf(distribution::BondLengthDistribution, distance::Real)::Float32
  if distance <= distribution.min_distance || distance >= distribution.max_distance
    return 0.0
  end
  idx = round(Int32, (distance - distribution.min_distance) / distribution.bucket_size) + 1
  return distribution.pdf[idx]
end

"""The AllBondLengthDistributions struct keeps a map from this type
  to the corresponding BondLengthDistribution.
"""
struct AtypesBtype
  type1::Int
  type2::Int
  btype::Int
end 

function Base.isequal(t1::AtypesBtype, t2::AtypesBtype)::Bool
  t1.type1 == t1.type1 && t1.type2 == t2.type2 && t1.btype == t2.btype
end

function Base.hash(atypes::AtypesBtype)
  return 200 * atypes.t1 + 100 * atypes.t2 + atypes.btype
end

mutable struct AllBondLengthDistributions
  distribution::Dict{AtypesBtype, BondLengthDistribution}
end

function add_file(fname::String, distributions::AllBondLengthDistributions, key::AtypesBtype)::Bool
  dist = BondLengthDistribution()
  if ! from_file(fname, dist)
    @warn("cannot read distribution from $(fname)")
    return false
  end

  distributions.distribution[key] = dist

  return true
end

function add_from_files(setm::String, distributions::AllBondLengthDistributions)::Bool
  atom_types = [1, 6, 7, 8, 9]
  bond_types = [0, 1, 2, 3]
  for atypes in combinations(atom_types, 2)
    t1 = atypes[1]
    t2 = atypes[2]
    for btype in bond_types
      fname = joinpath(sten, ".$(t1).$(btype).$(t2)")
      if ! isfile(fname)
        @warn("skipping non existent file $(fname)")
        continue
      end
      key = AtypesBtype(t1, t2, btype)
      if ! add_file(fname, distributions, key)
        @error("error adding $(fname)")
        return false
      end
    end
  end
end

function pdf(distributions::AllBondLengthDistributions, type1::Int, type2::Int, btype::Int, distance::Real)::Float32
  key = AtypesBtype(type1, type2, btype)
  value = get(distributions.distribution, key, nothing)
  if value !== nothing
    return 0.0
  end

  return value.pdf(distance)
end
