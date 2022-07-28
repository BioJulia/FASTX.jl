machine = let
    isinteractive() && @info "Compiling FASTQ FSM..." 

    re = Automa.RegExp
    
    hspace = re"[ \t\v]"
    
    header1 = let
        identifier = re.rep(re.any() \ re.space())
        identifier.actions[:enter] = [:mark]
        identifier.actions[:exit]  = [:header1_identifier]
        
        # Description here means "after whitespace", not whole line
        description = re.cat(re.any() \ re.space(), re"[^\r\n]*")
        re.cat('@', identifier, re.opt(re.cat(re.rep1(hspace), re.opt(description))))
    end
    header1.actions[:exit] = [:header1_description]
    
    sequence = re"[A-z]*"
    sequence.actions[:enter] = [:mark]
    sequence.actions[:exit]  = [:sequence]
    
    # The pattern recognized by header2 should be identical to header1
    # with the only difference being that h1 is split into identifier
    # and description
    header2 = let
        description2 = re"[^\r\n]+"
        description2.actions[:enter] = [:mark]
        description2.actions[:exit]  = [:header2_description]
        re.cat('+', re.opt(description2))
    end
    
    quality = re"[!-~]*"
    quality.actions[:enter] = [:mark]
    quality.actions[:exit]  = [:quality]
    
    newline = let
        lf = re"\n"
        lf.actions[:enter] = [:countline]
        re.cat(re.opt('\r'), lf)
    end
    
    record = re.cat(header1, newline, sequence, newline, header2, newline, quality)
    record.actions[:enter] = [:mark]
    record.actions[:exit] = [:record]
    
    fastq = re.opt(record) * re.rep(newline * record) * re.opt(newline)
    
    Automa.compile(fastq)
end

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
    :countline => :(linenum += 1),
    # Since the identifier is contained in the description, we just need
    # to store the length of the identifier. The bytes are copied in the description action
    :header1_identifier => :(record.identifier_len = Int32(@relpos(p-1))),

    # Copy description bytes and keep track of how many bytes copied
    :header1_description => quote
        let n = @relpos(p-1)
            appendfrom!(record.data, 1, data, @markpos, n)
            filled += n
            record.description_len = Int32(n)
        end
    end,

    # Copy sequence bytes and keep track of how many bytes copied
    :sequence => quote
        let n = @relpos(p-1)
            appendfrom!(record.data, filled + 1, data, @markpos, n)
            filled += n
            record.has_description_seq_len = UInt(n)
        end
    end,

    # Verify the second description is identical to the first,
    # and set the top bit in record.has_description_seq_len 
    :header2_description => quote
        let n = @relpos(p-1)
            recdata = record.data
            # This might look horribly unsafe, and it is.
            # But Automa should guarantee that data is under GC preservation, and that the
            # pointer p-n is valid. SHOULD, at least.
            GC.@preserve recdata begin
                if n != record.description_len || !iszero(memcmp(pointer(recdata), pointer(data, p-n), n%UInt))
                    error("First and second description line not identical")
                end
            end
            record.has_description_seq_len |= (0x80_00_00_00_00_00_00_00 % UInt)
        end
    end,

    # Verify the length of quality and sequence is identical, then copy bytes over
    :quality => quote
        let n = @relpos(p-1)
            n == seqlen(record) || error("Length of quality must be identical to length of sequence")
            appendfrom!(record.data, filled + 1, data, @markpos, n)
        end
    end,

    # Break from the loop
    :record => quote
        found = true
        @escape
    end,
)

initcode = quote
    pos = 0
    found = false
    filled = 0
    empty!(record)
    cs, linenum = state
end

loopcode = quote
    if cs < 0
        throw(ArgumentError("malformed FASTQ file at line $(linenum)"))
    end
    found && @goto __return__
end

returncode = :(return cs, linenum, found)

context = Automa.CodeGenContext(generator=:goto)

isinteractive() && @info "Generating FASTQ parsing code..."

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
