module SimpleFiniteElements

import StaticArrays: SVector

struct GeometryModel
    name            :: String
    entities        :: Vector{Tuple{Int32,Int32}}
    physical_groups :: Dict{String,Tuple{Int32,Int32}}
    entities_in     :: Dict{String,Vector{Int32}}
end

struct FEMesh
    gmodel             :: GeometryModel
    coord              :: Matrix{Float64}
    elt_type_in        :: Dict{String,Int32}
    elt_tags_in        :: Dict{String,Vector{Int64}}
    elt_node_tags_in   :: Dict{String,Matrix{Int64}}
    elt_properties     :: Dict{Int32,NamedTuple}
end

const FuncOrConst = Union{Function,Float64, SVector}

struct DegreesOfFreedom
    mesh            :: FEMesh
    essential_bc    :: Vector{Tuple}
    node_tag        :: Vector{Int64}
    num_free        :: Int64
    num_fixed       :: Int64
    node_tag_to_dof :: Vector{Int64}
    elt_dof         :: Dict{Int32,Vector{Int64}}
end

mutable struct Edge
    midpt          :: Int64
    left_phys_grp  :: Int64 # tag number of physical group 
    left_triangle  :: Int64 # elt_tag of triangle
    left_opposite  :: Int64 # index from {1,2,3} of vertex opposite the edge
    right_phys_grp :: Int64
    right_triangle :: Int64
    right_opposite :: Int64
end

struct EdgeEnumeration
    edge          :: Vector{Edge}
    midpt_to_edge :: Vector{Int64}
end

const SUCC = Int64[ 2, 3, 1 ]
const PRED = Int64[ 3, 1, 2 ]

include("submodules/Utils.jl")
include("submodules/MeshGen.jl")
include("submodules/FEM.jl")
include("submodules/Poisson.jl")
include("submodules/NonConformingPoisson.jl")
include("submodules/Elasticity.jl")
include("submodules/NonConformingElasticity.jl")

export GeometryModel, FEMesh, DegreesOfFreedom, EdgeEnumeration

import .MeshGen: get_nodal_values, reorder_dof_steklov!, max_elt_diameter
export get_nodal_values, reorder_dof_steklov!, max_elt_diameter
import .FEM: assemble_matrix, assemble_vector
export assemble_matrix, assemble_vector
import .Utils: SchurComplement
export SchurComplement

end # module SimpleFiniteElements
