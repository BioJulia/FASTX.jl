# Automa.jl generated readrecord! function
# ========================================

machine = (function ()
    isinteractive() && @info "Compiling FASTA FSM..." 

    re = Automa.RegExp
    
    hspace = re"[ \t\v]"
    
    identifier = re.rep(re.any() \ re.space())
    identifier.actions[:enter] = [:pos]
    identifier.actions[:exit]  = [:identifier]
    
    description = re.cat(re.any() \ re.space(), re"[^\n\r]*")
    description.actions[:enter] = [:pos]
    description.actions[:exit]  = [:description]

    header = re">" * identifier * re.opt(re.rep1(hspace) * re.opt(description))
    header.actions[:enter] = [:mark]
    header.actions[:exit]  = [:header]
    
    # '*': terminal, `-': gap
    letters = re"[A-Za-z*\-]+"
    letters.actions[:enter] = [:mark]
    letters.actions[:exit]  = [:letters]
    
    newline = let
        lf = re"\n"
        lf.actions[:enter] = [:countline]
        re.cat(re.opt('\r'), lf)
    end

    # Sequence: A sequence can be any free combination of newline, letters and hspace
    # It cannot be re.rep(letters | hspace | newline), because back-to-back repeated
    # letters would cause an FSM ambiguity between nothing and [:letters, :mark]
    sequence = re.rep(re.opt(letters) * (newline | hspace)) * re.opt(letters)
    
    record = re.cat(header, re.opt(re.cat(newline, sequence)), re.rep1(newline))
    record.actions[:exit] = [:record]
    
    record_trailing = re.cat(header, re.rep1(newline), sequence)
    record_trailing.actions[:exit] = [:record]
    
    fasta = re.cat(re.rep(newline), re.rep(record), re.opt(record_trailing))
    
    Automa.compile(fasta)
end)()

#write("fasta.dot", Automa.machine2dot(machine))
#run(`dot -Tsvg -o fasta.svg fasta.dot`)

actions = Dict(
    :mark => :(@mark),
    :pos => :(pos = @relpos(p)),
    :countline => :(linenum += 1),
    :identifier => :(record.identifier = pos:@relpos(p-1)),
    :description => :(record.description = pos:@relpos(p-1)),
    :header => quote
        let n = p - @markpos
            appendfrom!(record.data, 1, data, @markpos, n)
            filled += n
            appendfrom!(record.data, filled + 1, b"\n", 1, 1)
            filled += 1
        end
    end,
    :letters => quote
        let markpos = @markpos(), n = @relpos(p-1) - @relpos(markpos) + 1
            appendfrom!(record.data, filled + 1, data, markpos, n)
            if isempty(record.sequence)
                record.sequence = filled+1:filled+n
            else
                record.sequence = first(record.sequence):last(record.sequence)+n
            end
            filled += n
        end
    end,
    :record => quote
        record.filled = 1:filled
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
    initialize!(record)
    cs, linenum = state
end

loopcode = quote
    if cs < 0
        throw(ArgumentError("malformed FASTA file at line $(linenum)"))
    end
    found && @goto __return__
end

returncode = :(return cs, linenum, found)

context = Automa.CodeGenContext(generator = :goto, checkbounds = false, loopunroll = 8)

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
