module ReaderIterTools

abstract type AbstractRecordIterator{T} end

function eachrecord end

mutable struct RecordIteratorState{T}
    # Machine state
    state::Int
    # Line number
    linenum::Int
    # Is record filled?
    filled::Bool
    # Placeholder
    record::T
end

Base.iteratorsize(::Type{<:AbstractRecordIterator}) = Base.SizeUnknown()

function Base.eltype(::Type{AbstractRecordIterator{T}}) where T
    return T
end

function Base.start(iter::AbstractRecordIterator{T}) where T
    return ReaderIterTools.RecordIteratorState(1, 1, false, T())
end

function Base.next(iter::AbstractRecordIterator{T}, state) where T
    @assert state.filled
    record = copy(state.record)
    state.filled = false
    return record, state
end