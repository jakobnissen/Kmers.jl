```@meta
CurrentModule = Kmers
DocTestSetup = quote
    using BioSequences
    using Test
    using Kmers
    using FASTX
    using MinHash
end
```
## MinHash
The MinHash algorithm is used in tools such as [Mash](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0997-x) and [sourmash](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6720031/) to quickly compute approximate similarities of genomes, collections of genomes, or collections of reads.

```jldoctest; filter = r"^\d\d\d? MB/s$" => s"***"
using BioSequences, MinHash, FASTX, Kmers

# Write 100 MB of DNA in 50 genomes to buffer
buffer = IOBuffer()
writer = FASTAWriter(buffer)
n_bytes = sum(1:50) do genome
    rec = FASTARecord("seq_$(genome)", randdnaseq(2_000_000))
    write(writer, rec)
end
flush(writer)

# Time minhashing the 50 genomes
timing = @timed FASTAReader(seekstart(buffer); copy=false) do reader
    map(reader) do record
        seq = codeunits(sequence(record))
        minhash(fx_hash, CanonicalDNAMers{16}(sequence(record)), 1000)
    end
end
println(round(Int, n_bytes / (timing.time * 1e6)), " MB/s")

# output

200 MB/s
```
