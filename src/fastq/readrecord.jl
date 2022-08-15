machine = let
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
        throw_parser_error(data, p, linenum)
    end
    found && @goto __return__
end

returncode = :(return cs, linenum, found)

Automa.Stream.generate_reader(
    :readrecord!,
    machine,
    arguments = (:(record::Record), :(state::Tuple{Int,Int})),
    actions = actions,
    context = CONTEXT,
    initcode = initcode,
    loopcode = loopcode,
    returncode = returncode
) |> eval

validator_actions = Dict(
    :mark => :(@mark),
    :countline => :(linenum += 1),
    :header1_identifier => quote nothing end,

    # Copy description to buffer to check if second description is same
    :header1_description => quote
        let n = @relpos(p-1)
            appendfrom!(headerbuffer, 1, data, @markpos, n)
            description_len = n
        end
    end,

    # Copy sequence bytes and keep track of how many bytes copied
    :sequence => :(sequence_len = @relpos(p-1)),

    # Verify the second description is identical to the first,
    :header2_description => quote
        let n = @relpos(p-1)
            n == description_len || return linenum
            GC.@preserve headerbuffer begin
                iszero(memcmp(pointer(headerbuffer), pointer(data, p-n), n%UInt)) || return linenum
            end
        end
    end,

    # Verify the length of quality and sequence is identical
    :quality => quote
        let n = @relpos(p-1)
            n == sequence_len || return linenum
        end
    end,
    :record => quote nothing end
)

initcode = quote
    linenum = 1
    description_len = 0
    sequence_len = 0
    headerbuffer = Vector{UInt8}(undef, 1024)
end

Automa.Stream.generate_reader(
    :validate_fastq,
    machine,
    arguments = (),
    actions= validator_actions,
    context = CONTEXT,
    initcode = initcode,
    loopcode = :(cs < 0 && return linenum),
    returncode = :(iszero(cs) ? nothing : linenum)
) |> eval

# Currently returns linenumber if it is not, but we might remove
# this from the readers, since this state cannot be kept when seeking.
"""
    validate_fastq(io::IO) >: Nothing

Check if `io` is a valid FASTQ file.
Return `nothing` if it is, and an instance of another type if not.

# Examples
```jldoctest
julia> validate_fastq(IOBuffer("@i1 r1\\nuuag\\n+\\nHJKI")) === nothing
true

julia> validate_fastq(IOBuffer("@i1 r1\\nu;ag\\n+\\nHJKI")) === nothing
false
```
"""
validate_fastq(io::IO) = validate_fastq(NoopStream(io))
