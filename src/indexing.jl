@inline encoded_data(x::Kmer) = x.data

@inline bitindex(seq::Kmer, i::Integer) = bitindex(BioSequences.BitsPerSymbol(seq), encoded_data_eltype(seq), i + n_unused(seq))