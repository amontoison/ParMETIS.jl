const PARMETIS_MAJOR_VERSION = 4
const PARMETIS_MINOR_VERSION = 0
const PARMETIS_SUBMINOR_VERSION = 3

const idx_t = Cint
const real_t = Cfloat

# Operation type codes
const PARMETIS_OP_KMETIS = Cint(0)
const PARMETIS_OP_GKMETIS = Cint(1)
const PARMETIS_OP_GMETIS = Cint(2)
const PARMETIS_OP_RMETIS = Cint(3)
const PARMETIS_OP_AMETIS = Cint(4)
const PARMETIS_OP_OMETIS = Cint(5)
const PARMETIS_OP_M2DUAL = Cint(6)
const PARMETIS_OP_MKMETIS = Cint(7)

# Return codes
const METIS_OK = Cint(1)   # Returned normally
const METIS_ERROR_INPUT = Cint(-2)  # Returned due to erroneous inputs and/or options
const METIS_ERROR_MEMORY = Cint(-3)  # Returned due to insufficient memory
const METIS_ERROR = Cint(-4)  # Some other errors

# Matching types
const PARMETIS_MTYPE_LOCAL = Cint(1)  # Restrict matching to within processor vertices
const PARMETIS_MTYPE_GLOBAL = Cint(2)  # Remote vertices can be matched

# Separator refinement types
const PARMETIS_SRTYPE_GREEDY = Cint(1)  # Vertices are visted from highest to lowest gain
const PARMETIS_SRTYPE_2PHASE = Cint(2)  # Separators are refined in a two-phase fashion using PARMETIS_SRTYPE_GREEDY for the 2nd phase

# Coupling types for ParMETIS_V3_RefineKway & ParMETIS_V3_AdaptiveRepart
const PARMETIS_PSR_COUPLED = Cint(1)  # number of partitions == number of processors
const PARMETIS_PSR_UNCOUPLED = Cint(2)  # number of partitions != number of processors

# Debug levels
const PARMETIS_DBGLVL_TIME = Cint(1)   # Perform timing analysis
const PARMETIS_DBGLVL_INFO = Cint(2)   # Perform timing analysis
const PARMETIS_DBGLVL_PROGRESS = Cint(4)   # Show the coarsening progress
const PARMETIS_DBGLVL_REFINEINFO = Cint(8)   # Show info on communication during folding
const PARMETIS_DBGLVL_MATCHINFO = Cint(16)  # Show info on matching
const PARMETIS_DBGLVL_RMOVEINFO = Cint(32)  # Show info on communication during folding
const PARMETIS_DBGLVL_REMAP = Cint(64)  # Determines if remapping will take place

struct MetisError <: Exception
  code::Cint
end

function Base.showerror(io::IO, me::MetisError)
  print(io, "MetisError: ")
  if me.code == Metis.METIS_ERROR_INPUT
    print(io, "input error")
  elseif me.code == Metis.METIS_ERROR_MEMORY
    print(io, "could not allocate the required memory")
  else
    print(io, "unknown error")
  end
  print(io, " (error code $(me.code)).")
end

# Partition a graph
function ParMETIS_V3_PartKway(
  vtxdist,
  xadj,
  adjncy,
  vwgt,
  adjwgt,
  wgtflag,
  numflag,
  ncon,
  nparts,
  tpwgts,
  ubvec,
  options,
  edgecut,
  part,
  comm,
)
  r = ccall(
    (:ParMETIS_V3_PartKway, libparmetis),
    Cint,
    (
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{real_t},
      Ptr{real_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
    ),
    vtxdist,
    xadj,
    adjncy,
    vwgt,
    adjwgt,
    wgtflag,
    numflag,
    ncon,
    nparts,
    tpwgts,
    ubvec,
    options,
    edgecut,
    part,
    comm,
  )
  r == METIS_OK || throw(MetisError(r))
  return nothing
end

function ParMETIS_V3_PartGeomKway(
  vtxdist,
  xadj,
  adjncy,
  vwgt,
  adjwgt,
  wgtflag,
  numflag,
  ndims,
  xyz,
  ncon,
  nparts,
  tpwgts,
  ubvec,
  options,
  edgecut,
  part,
  comm,
)
  r = ccall(
    (:ParMETIS_V3_PartGeomKway, libparmetis),
    Cint,
    (
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{real_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{real_t},
      Ptr{real_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
    ),
    vtxdist,
    xadj,
    adjncy,
    vwgt,
    adjwgt,
    wgtflag,
    numflag,
    ndims,
    xyz,
    ncon,
    nparts,
    tpwgts,
    ubvec,
    options,
    edgecut,
    part,
    comm,
  )
  r == METIS_OK || throw(MetisError(r))
  return nothing
end

function ParMETIS_V3_PartGeom(vtxdist, ndims, xyz, part, comm)
  r = ccall(
    (:ParMETIS_V3_PartGeom, libparmetis),
    Cint,
    (Ptr{idx_t}, Ptr{idx_t}, Ptr{real_t}, Ptr{idx_t}, Ptr{idx_t}),
    vtxdist,
    ndims,
    xyz,
    part,
    comm,
  )
  r == METIS_OK || throw(MetisError(r))
  return nothing
end

# Partition a mesh
function ParMETIS_V3_PartMeshKway(
  elmdist,
  eptr,
  eind,
  elmwgt,
  wgtflag,
  numflag,
  ncon,
  ncommonnodes,
  nparts,
  tpwgts,
  ubvec,
  options,
  edgecut,
  part,
  comm,
)
  r = ccall(
    (:ParMETIS_V3_PartMeshKway, libparmetis),
    Cint,
    (
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{real_t},
      Ptr{real_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
    ),
    elmdist,
    eptr,
    eind,
    elmwgt,
    wgtflag,
    numflag,
    ncon,
    ncommonnodes,
    nparts,
    tpwgts,
    ubvec,
    options,
    edgecut,
    part,
    comm,
  )
  r == METIS_OK || throw(MetisError(r))
  return nothing
end

# Construct a graph from a mesh
function ParMETIS_V3_Mesh2Dual(elmdist, eptr, eind, numflag, ncommonnodes, xadj, adjncy, comm)
  r = ccall(
    (:ParMETIS_V3_Mesh2Dual, libparmetis),
    Cint,
    (
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{Ptr{idx_t}},
      Ptr{Ptr{idx_t}},
      Ptr{idx_t},
    ),
    elmdist,
    eptr,
    eind,
    numflag,
    ncommonnodes,
    xadj,
    adjncy,
    comm,
  )
  r == METIS_OK || throw(MetisError(r))
  return nothing
end

# Repartition a graph corresponding to an adaptively refined mesh
function ParMETIS_V3_AdaptiveRepart(
  vtxdist,
  xadj,
  adjncy,
  vwgt,
  vsize,
  adjwgt,
  wgtflag,
  numflag,
  ncon,
  nparts,
  tpwgts,
  ubvec,
  ipc2redist,
  options,
  edgecut,
  part,
  comm,
)
  r = ccall(
    (:ParMETIS_V3_AdaptiveRepart, libparmetis),
    Cint,
    (
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{real_t},
      Ptr{real_t},
      Ptr{real_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
    ),
    vtxdist,
    xadj,
    adjncy,
    vwgt,
    vsize,
    adjwgt,
    wgtflag,
    numflag,
    ncon,
    nparts,
    tpwgts,
    ubvec,
    ipc2redist,
    options,
    edgecut,
    part,
    comm,
  )
  r == METIS_OK || throw(MetisError(r))
  return nothing
end

# Refine the quality of a partitioning
function ParMETIS_V3_RefineKway(
  vtxdist,
  xadj,
  adjncy,
  vwgt,
  adjwgt,
  wgtflag,
  numflag,
  ncon,
  nparts,
  tpwgts,
  ubvec,
  options,
  edgecut,
  part,
  comm,
)
  r = ccall(
    (:ParMETIS_V3_RefineKway, libparmetis),
    Cint,
    (
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{real_t},
      Ptr{real_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
    ),
    vtxdist,
    xadj,
    adjncy,
    vwgt,
    adjwgt,
    wgtflag,
    numflag,
    ncon,
    nparts,
    tpwgts,
    ubvec,
    options,
    edgecut,
    part,
    comm,
  )
  r == METIS_OK || throw(MetisError(r))
  return nothing
end

# Compute a fill−reducing ordering
function ParMETIS_V3_NodeND(vtxdist, xadj, adjncy, numflag, options, order, sizes, comm)
  r = ccall(
    (:ParMETIS_V3_NodeND, libparmetis),
    Cint,
    (
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
    ),
    vtxdist,
    xadj,
    adjncy,
    numflag,
    options,
    order,
    sizes,
    comm,
  )
  r == METIS_OK || throw(MetisError(r))
  return nothing
end

function ParMETIS_V32_NodeND(
  vtxdist,
  xadj,
  adjncy,
  vwgt,
  numflag,
  mtype,
  rtype,
  p_nseps,
  s_nseps,
  ubfrac,
  seed,
  dbglvl,
  order,
  sizes,
  comm,
)
  r = ccall(
    (:ParMETIS_V32_NodeND, libparmetis),
    Cint,
    (
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{real_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
    ),
    vtxdist,
    xadj,
    adjncy,
    vwgt,
    numflag,
    mtype,
    rtype,
    p_nseps,
    s_nseps,
    ubfrac,
    seed,
    dbglvl,
    order,
    sizes,
    comm,
  )
  r == METIS_OK || throw(MetisError(r))
  return nothing
end

function ParMETIS_SerialNodeND(vtxdist, xadj, adjncy, numflag, options, order, sizes, comm)
  r = ccall(
    (:ParMETIS_SerialNodeND, libparmetis),
    Cint,
    (
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
      Ptr{idx_t},
    ),
    vtxdist,
    xadj,
    adjncy,
    numflag,
    options,
    order,
    sizes,
    comm,
  )
  r == METIS_OK || throw(MetisError(r))
  return nothing
end
