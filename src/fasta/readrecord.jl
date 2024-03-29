# Automa.jl generated readrecord! function
# ========================================

# The current implementation of the machine has the following debatable choices
# * You can put anything except \r and \n in the description, including trailing
#   whitespace.
# * You can put anything in the sequence lines except \n, \r and, and cannot begin with >.
#   The sequence must include at least one newline, i.e. ">A>A" is not valid FASTA,
#   but ">A\n>A\n" is. The newlines are not considered part of the sequence lines themselves.
#   This implies all whitespace except newlines, including trailing whitespace, is part
#   of the sequence.
machine = let
    hspace = re"[ \t\v]"
    newline = let
        lf = onenter!(re"\n", :countline)
        Re.opt('\r') * lf
    end
    space = hspace | newline
    
    # Identifier: Leading non-space
    identifier = onexit!(onenter!(Re.rep(Re.any() \ Re.space()), :mark), :identifier)
    
    # Description here include trailing whitespace.
    # This is needed for the FSM, since the description can contain arbitrary
    # whitespace, the only way to know the description ends is to encounter a newline.
    # NB: Make sure to also change the Index machine to match this is you change it.
    description = onexit!(identifier * Re.opt(hspace * re"[^\r\n]*"), :description)

    # Header: '>' then description
    header = re">" * description
    
    # Sequence line: Anything except \r, \n and >
    # Note: Must be consistent with the disallowed bytes in reader.jl used for seeking
    sequence_line = onexit!(onenter!(re"[^\n\r>]+", :mark), :seqline)

    # Sequence: This is intentionally very liberal with whitespace.
    # Any trailing whitespace is simply considered part of the sequence.
    # Is this bad? Maybe.
    sequence = Re.rep1(Re.opt(sequence_line) * Re.rep1(newline))
    
    # We have sequence_eof to allow the final sequence to not end in whitespace
    sequence_eof = Re.opt(sequence_line) * Re.rep(Re.rep1(newline) * Re.opt(sequence_line))

    record = onexit!(header * newline * sequence, :record)
    record_eof = onexit!(header * newline * sequence_eof, :record)

    fasta = Re.rep(space) * Re.rep(record) * Re.opt(record_eof)
    
    Automa.compile(fasta)
end

actions = Dict(
    :mark => :(@mark),
    :countline => :(linenum += 1),
    :identifier => :(record.identifier_len = Int32(@relpos(p-1))),
    # Append entire header line to buffer from pos 1
    :description => quote
        let n = @relpos(p-1)
            appendfrom!(record.data, 1, data, @markpos, n)
            filled += n
            record.description_len = Int32(n)
        end
    end,
    # Append sequence line to buffer
    :seqline => quote
        let n = @relpos(p-1)
            appendfrom!(record.data, filled + 1, data, @markpos, n)
            filled += n
            record.sequence_len += n
        end
    end,
    :record => quote
        found = true
        @escape
    end
)

initcode = quote
    pos = 0
    filled = 0
    found = false
    empty!(record)
    cs, linenum = state
end

loopcode = quote
    if cs < 0
        throw_parser_error(data, p, linenum < 0 ? nothing : linenum)
    end
    found && @goto __return__
end

returncode = :(return cs, linenum, found)

Automa.generate_reader(
    :readrecord!,
    machine,
    arguments = (:(record::Record), :(state::Tuple{Int,Int})),
    actions = actions,
    context = CONTEXT,
    initcode = initcode,
    loopcode = loopcode,
    returncode = returncode
) |> eval

validator_actions = Dict(k => quote nothing end for k in keys(actions))
validator_actions[:countline] = :(linenum += 1)

Automa.generate_reader(
    :validate_fasta,
    machine,
    arguments = (),
    actions= validator_actions,
    context = CONTEXT,
    initcode = :(linenum = 1),
    loopcode = :(cs < 0 && return linenum),
    returncode = :(iszero(cs) ? nothing : linenum)
) |> eval

# Currently returns linenumber if it is not, but we might remove
# this from the readers, since this state cannot be kept when seeking.
"""
    validate_fasta(io::IO) >: Nothing

Check if `io` is a valid FASTA file.
Return `nothing` if it is, and an instance of another type if not.

# Examples
```jldoctest
julia> validate_fasta(IOBuffer(">a bc\\nTAG\\nTA")) === nothing
true

julia> validate_fasta(IOBuffer(">a bc\\nT>G\\nTA")) === nothing
false
```
"""
validate_fasta(io::IO) = validate_fasta(NoopStream(io))