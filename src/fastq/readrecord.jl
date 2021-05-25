machine = (function ()
    isinteractive() && @info "Compiling FASTQ FSM..." 

    re = Automa.RegExp
    
    hspace = re"[ \t\v]"
    
    header1 = let
        identifier = re.rep(re.any() \ re.space())
        identifier.actions[:enter] = [:pos]
        identifier.actions[:exit]  = [:header1_identifier]
        
        description = re.cat(re.any() \ re.space(), re"[^\r\n]*")
        description.actions[:enter] = [:pos]
        description.actions[:exit]  = [:header1_description]
        
        re.cat('@', identifier, re.opt(re.cat(re.rep1(hspace), re.opt(description))))
    end
    
    sequence = re"[A-z]*"
    sequence.actions[:enter] = [:pos]
    sequence.actions[:exit]  = [:sequence]
    
    header2 = let
        identifier = re.rep1(re.any() \ re.space())
        description = re.cat(re.any() \ hspace, re"[^\r\n]*")
        re.cat('+', re.opt(re.cat(identifier, re.opt(re.cat(re.rep1(hspace), description)))))
    end
    header2.actions[:enter] = [:pos]
    header2.actions[:exit]  = [:header2]
    
    quality = re"[!-~]*"
    quality.actions[:enter] = [:pos]
    quality.actions[:exit]  = [:quality]
    
    newline = let
        lf = re"\n"
        lf.actions[:enter] = [:countline]
        re.cat(re.opt('\r'), lf)
    end

    sep = re.opt('\r') * re"\n"
    sep.actions[:exit] = [:countline, :sep]
    
    record = re.cat(header1, newline, sequence, newline, header2, newline, quality)
    record.actions[:enter] = [:mark]
    record.actions[:exit] = [:record]
    
    fastq = re.rep(re.cat(record, sep)) * re.opt(record) * re.rep(sep)
    
    Automa.compile(fastq)
end)()

#write("fastq.dot", Automa.machine2dot(machine))
#run(`dot -Tsvg -o fastq.svg fastq.dot`)

function appendfrom!(dst, dpos, src, spos, n)
    if length(dst) < dpos + n - 1
        resize!(dst, dpos + n - 1)
    end
    copyto!(dst, dpos, src, spos, n)
    return dst
end

function is_same_mem(data, pos1, pos2, len)
    checkbounds(data, 1:pos1+len-1)
    checkbounds(data, 1:pos2+len-1)
    return ccall(:memcmp, Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Csize_t), pointer(data, pos1), pointer(data, pos2), len) == 0
end

actions = Dict(
    :mark => :(@mark),
    :pos => :(pos = @relpos(p)),
    :countline => :(linenum += 1),
    :header1_identifier => :(record.identifier = pos:@relpos(p-1)),
    :header1_description => :(record.description = pos:@relpos(p-1)),
    :sequence => :(record.sequence = pos:@relpos(p-1)),
    :header2 => :(second_header_pos = pos+1; second_header_len = @relpos(p)-(pos+1)),
    :quality => :(record.quality = pos:@relpos(p-1)),
    :sep => :(@escape),
    :record => quote
        appendfrom!(record.data, 1, data, @markpos, p-@markpos)
        record.filled = 1:(p-@markpos)
        found = true
    end,
)

initcode = quote
    pos = 0
    second_header_pos = 0
    second_header_len = 0
    found = false
    initialize!(record)
    cs, linenum = state
end

loopcode = quote
    if cs < 0
        throw(ArgumentError("malformed FASTQ file at line $(linenum)"))
    elseif found && length(record.sequence) != length(record.quality)
        throw(ArgumentError("mismatched sequence and quality length"))
    elseif found && second_header_len > 0
        first_header_pos = first(record.identifier)
        first_header_len = max(last(record.identifier), last(record.description)) - first_header_pos + 1
        if first_header_len != second_header_len || !is_same_mem(record.data, first_header_pos, second_header_pos, first_header_len)
            throw(ArgumentError("mismatched headers"))
        end
    elseif found && transform != nothing
        transform(record.data, record.sequence)
    end
    found && @goto __return__
end

returncode = :(return cs, linenum, found)

context = Automa.CodeGenContext(generator=:simd, checkbounds=false)

isinteractive() && @info "Generating FASTQ parsing code..."

Automa.Stream.generate_reader(
    :readrecord!,
    machine,
    arguments = (:(record::Record), :(state::Tuple{Int,Int}), :(transform)),
    actions = actions,
    context = context,
    initcode = initcode,
    loopcode = loopcode,
    returncode = returncode
) |> eval
