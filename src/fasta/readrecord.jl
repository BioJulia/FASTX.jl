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
    isinteractive() && @info "Compiling FASTA FSM..." 

    re = Automa.RegExp
    
    hspace = re"[ \t\v]"
    newline = let
        lf = re"\n"
        lf.actions[:enter] = [:countline]
        re.opt('\r') * lf
    end
    space = hspace | newline
    
    # Identifier: Leading non-space
    identifier = re.rep(re.any() \ re.space())
    identifier.actions[:enter] = [:mark]
    # Action: Store length of identifier
    identifier.actions[:exit]  = [:identifier]
    
    # Description here include trailing whitespace.
    # This is needed for the FSM, since the description can contain arbitrary
    # whitespace, the only way to know the description ends is to encounter a newline.
    # NB: Make sure to also change the Index machine to match this is you change it.
    description = identifier * re.opt(hspace * re"[^\r\n]*")

    # Action: Store length of description and append description to record.data
    description.actions[:exit]  = [:description]

    # Header: '>' then description
    header = re">" * description
    
    # Sequence line: Anything except \r, \n and >
    sequence_line = re"[^\n\r>]+"
    sequence_line.actions[:enter] = [:mark]
    # Action: Append letters to sequence_line
    sequence_line.actions[:exit]  = [:seqline]

    # Sequence: This is intentionally very liberal with whitespace.
    # Any trailing whitespace is simply considered part of the sequence.
    # Is this bad? Maybe.
    sequence = re.rep1(re.opt(sequence_line) * re.rep1(newline))
    
    # We have sequence_eof to allow the final sequence to not end in whitespace
    sequence_eof = re.opt(sequence_line) * re.rep(re.rep1(newline) * re.opt(sequence_line))

    record = header * newline * sequence
    record.actions[:exit] = [:record]
    record_eof = header * newline * sequence_eof
    record_eof.actions[:exit] = [:record]

    fasta = re.rep(space) * re.rep(record) * re.opt(record_eof)
    
    Automa.compile(fasta)
end

#write("fasta.dot", Automa.machine2dot(machine))
#run(`dot -Tsvg -o fasta.svg fasta.dot`)

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

function appendfrom!(dst, dpos, src, spos, n)
    if length(dst) < dpos + n - 1
        resize!(dst, dpos + n - 1)
    end
    copyto!(dst, dpos, src, spos, n)
    return dst
end

initcode = quote
    pos = 0
    filled = 0
    found = false
    empty!(record)
    cs, linenum = state
end

loopcode = quote
    if cs < 0
        throw(ArgumentError("malformed FASTA file at line $(linenum)"))
    end
    found && @goto __return__
end

returncode = :(return cs, linenum, found)

context = Automa.CodeGenContext(generator=:goto)

isinteractive() && @info "Generating FASTA parsing code..."

Automa.Stream.generate_reader(
    :readrecord!,
    machine,
    arguments = (:(record::Record), :(state::Tuple{Int,Int})),
    actions = actions,
    context = context,
    initcode = initcode,
    loopcode = loopcode,
    returncode = returncode
) |> eval
