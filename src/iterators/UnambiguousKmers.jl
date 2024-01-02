"""
    UnambiguousKmers{A <: TwoBit, K, S}

Iterator of 2-bit nucleic acid kmers. This differs from `FwKmers` in that any kmers
containing ambiguous nucleotides are skipped, whereas using `FwKmers`, they result
in an error.

Can be constructed more conventiently with the constructors `UnambiguousDNAMers{K}(s)`
and `UnambiguousRNAMers{K}(s)`.

!!! note
    To obtain canonical unambiguous kmers, simply call `canonical` on each kmer output
by `UnambiguousKmers`.

# Examples:
```jldoctest
julia> it = UnambiguousRNAMers{4}(dna"TGAGCWKCATC");

julia> collect(it)
3-element Vector{Tuple{Kmer{RNAAlphabet{2}, 4, 1}, Int64}}:
 (UGAG, 1)
 (GAGC, 2)
 (CAUC, 8)
```
"""
struct UnambiguousKmers{A <: TwoBit, K, S}
    it::FwKmers{A, K, S}
end

Base.IteratorSize(::Type{<:UnambiguousKmers}) = Base.SizeUnknown()
function Base.eltype(::Type{<:UnambiguousKmers{A, K}}) where {A, K}
    Tuple{derive_type(Kmer{A, K}), Int}
end

source_type(::Type{UnambiguousKmers{A, K, S}}) where {A, K, S} = S

# Constructors
function UnambiguousKmers{A, K}(s::S) where {S, A <: TwoBit, K}
    UnambiguousKmers{A, K, S}(FwKmers{A, K}(s))
end
function UnambiguousKmers{A, K, S}(s::S) where {S, A <: TwoBit, K}
    UnambiguousKmers{A, K, S}(FwKmers{A, K}(s))
end

"`UnambiguousDNAMers{K, S}`: Alias for `UnambiguousKmers{DNAAlphabet{2}, K, S}`"
const UnambiguousDNAMers{K, S} = UnambiguousKmers{DNAAlphabet{2}, K, S}

"`UnambiguousRNAMers{K, S}`: Alias for `UnambiguousKmers{RNAAlphabet{2}, K, S}`"
const UnambiguousRNAMers{K, S} = UnambiguousKmers{RNAAlphabet{2}, K, S}

@inline function Base.iterate(it::UnambiguousKmers{A, K, S}) where {A, K, S}
    T = derive_type(Kmer{A, K})
    state = (T(unsafe, zero_tuple(T)), K, 1)
    iterate_kmer(RecodingScheme(A(), S), it, state)
end

@inline function Base.iterate(it::UnambiguousKmers{A, K, S}, state) where {A, K, S}
    iterate_kmer(RecodingScheme(A(), S), it, state)
end

@inline function iterate_kmer(
    ::RecodingScheme,
    it::UnambiguousKmers,
    state::Tuple{Kmer, Int, Int},
)
    (kmer, remaining, index) = state
    K = ksize(kmer)
    while !iszero(remaining)
        index > lastindex(it.it.seq) && return nothing
        symbol = convert(typeof(kmer), it.it.seq[index])
        index += 1
        if isambiguous(symbol)
            remaining = K
        else
            remaining -= 1
            kmer = shift(kmer, symbol)
        end
    end
    ((kmer, index - K), (kmer, 1, index))
end

# Here, we can forward directly to FwKmers
@inline function iterate_kmer(
    ::Copyable,
    it::UnambiguousKmers,
    state::Tuple{Kmer, Int, Int},
)
    (kmer, _, index) = state
    itval = iterate(it.it, (kmer, index))
    if itval === nothing
        nothing
    else
        (kmer, s) = itval
        ((kmer, index - ksize(kmer)), s)
    end
end

@inline function iterate_kmer(
    ::AsciiEncode,
    it::UnambiguousKmers{A, K, S},
    state::Tuple{Kmer, Int, Int},
) where {A <: TwoBit, K, S}
    src = used_source(RecodingScheme(A(), S), it.it.seq)
    Base.require_one_based_indexing(src)
    (kmer, remaining, index) = state
    while !iszero(remaining)
        index > lastindex(src) && return nothing
        byte = @inbounds src[index]
        index += 1
        encoding = @inbounds ASCII_SKIPPING_LUT[(byte + 0x01) % Int]
        if encoding == 0xff
            throw(BioSequences.EncodeError(Alphabet(kmer), repr(byte)))
        elseif encoding == 0xf0
            remaining = K
        else
            remaining -= 1
            kmer = shift_encoding(kmer, encoding % UInt)
        end
    end
    ((kmer, index - K), (kmer, 1, index))
end

@inline function iterate_kmer(
    ::FourToTwo,
    it::UnambiguousKmers{A, K, S},
    state::Tuple{Kmer, Int, Int},
) where {A <: TwoBit, K, S}
    (kmer, remaining, index) = state
    while !iszero(remaining)
        index > lastindex(it.it.seq) && return nothing
        encoding = UInt(BioSequences.extract_encoded_element(it.it.seq, index))::UInt
        kmer = shift_encoding(kmer, (trailing_zeros(encoding)) % UInt)
        index += 1
        remaining = isone(count_ones(encoding)) ? remaining - 1 : K
    end
    ((kmer, index - K), (kmer, 1, index))
end
