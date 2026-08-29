using Test
using Kmers
using BioSequences
using BitIntegers

@testset "Construction" begin
    @testset "Same alphabet, and ASCII alphabet" begin
        for s in [
                dna"ATGTCGTTAGT",
                dna"",
                dna"YATGC-ATGwTCTDV",
                rna"AUGUCGAGUGUAUGC",
                rna"AUGCWSYGANN--CA",
                aa"KWOP",
                aa"",
                aa"TYAPC",
            ]
            for i in (s, string(s), (i for i in s))
                m = Oligomer{typeof(Alphabet(s)), UInt64}(i)
                @test length(m) == length(i)
                @test string(m) == string(s)
            end
        end

        # Test macro equivalence with constructor
        @test dmer"ATGTCGTTAGT"d == DNAOligomer{UInt64}(dna"ATGTCGTTAGT")
        @test dmer"AUGUCGAGUGUAUGC"r == RNAOligomer{UInt64}(rna"AUGUCGAGUGUAUGC")
        @test dmer"KWOP"a == AAOligomer{UInt64}(aa"KWOP")

        d = DNAOligomer{UInt}(dna"ATGTCGTTAGT")
        @test RNAOligomer(d) == d
        @test RNAOligomer{UInt64}(d) == d

        # From a large kmer
        m = mer"TAGTGCTGTAGTAGTGCTGTATGATGTCTGCATGC"d
        dm = DNAOligomer{UInt256}(m)
        @test LongSequence(m) == LongSequence(dm)

        # Using large bitintegers - 240 2-bit nt = 480 bits should fit in
        # 512 bits with plenty room for the runtime length
        seq = randdnaseq(240)
        dm = DNAOligomer{UInt512}(seq)
        @test string(seq) == string(dm)
    end

    @testset "Two to four bit alphabet" begin
        for s in Any[
                dna"ATGCTGTGACCA",
                dna"ATGTCGA",
                dna"",
            ]
            for A in Any[DNAAlphabet, RNAAlphabet]
                for (srcB, dstB) in [(2, 4), (4, 2)]
                    src = LongSequence{A{srcB}}(s)
                    dst = Oligomer{A{dstB}, UInt64}(src)

                    @test string(src) == string(dst)
                end
            end
        end
    end

    @testset "Four to two bit alphabet" begin
        for s in [dna"TAGCTGAC", dna"ATGCTA", dna""]
            for A in [DNAAlphabet{2}, RNAAlphabet{2}]
                m = Oligomer{A, UInt64}(LongSequence{DNAAlphabet{4}}(s))
                @test string(m) == string(LongSequence{A}(s))
            end
        end
    end

    @testset "Generic alphabet" begin
        for s in ["HE", "", "中Å!"]
            m = Oligomer{CharAlphabet, UInt512}(s)
            @test length(m) == length(s)
            @test string(m) == s
        end
    end

    @testset "From Kmer" begin
        for s in [dna"TAGCTA", rna"UGCUGA", aa"PLKWM"]
            kmer = Kmer{typeof(Alphabet(s)), length(s)}(s)
            dkmer = Oligomer{typeof(Alphabet(s)), UInt64}(kmer)
            @test string(dkmer) == string(s)
            @test length(dkmer) == length(kmer)
        end
    end

    @testset "From Kmer into wider backing types" begin
        for U in [UInt128, UInt256, UInt512]
            # One symbol, a partially filled word, and a full source word.
            for s in [dna"T", dna"TAG", dna"T"^32]
                kmer = DNAKmer{length(s)}(s)
                oligomer = DNAOligomer{U}(kmer)
                @test string(oligomer) == string(kmer)
                @test DNAKmer{length(s)}(oligomer) == kmer
            end

            # Keep a multiword source case covered too.
            s = dna"T"^33
            kmer = DNAKmer{33}(s)
            oligomer = DNAOligomer{U}(kmer)
            @test string(oligomer) == string(kmer)
            @test DNAKmer{33}(oligomer) == kmer
        end
    end

    @testset "From Oligomer" begin
        # Same backing type, same alphabet
        s = dna"TAGCTGA"
        d1 = dmer"TAGCTGA"d
        d2 = DNAOligomer{UInt64}(d1)
        @test d2 == d1
        @test d2 === d1  # Should be identical
        @test length(d2) == length(d1)

        # Different backing type, same alphabet - widening
        d32 = DNAOligomer{UInt32}(dna"TAGC")
        d64 = DNAOligomer{UInt64}(d32)
        @test_throws MethodError d64 == d32
        @test_throws MethodError d32 == d64
        @test length(d64) == length(d32)
        @test string(d64) == "TAGC"

        # Different backing type, same alphabet - narrowing
        d128 = DNAOligomer{UInt128}(dna"TAGC")
        d32_narrow = DNAOligomer{UInt32}(d128)
        @test_throws MethodError d32_narrow == d128
        @test_throws MethodError d128 == d32_narrow
        @test length(d32_narrow) == length(d128)
        @test string(d32_narrow) == "TAGC"

        # Same alphabet family (DNA/RNA), different alphabets
        d_dna = dmer"ATGTCGTTAGT"d
        d_rna = RNAOligomer{UInt64}(d_dna)
        @test d_rna == d_dna
        @test length(d_rna) == length(d_dna)

        # Different backing type AND different alphabet
        d_dna32 = DNAOligomer{UInt32}(dna"TAGC")
        d_rna64 = RNAOligomer{UInt64}(d_dna32)
        @test_throws MethodError d_rna64 == d_dna32
        @test_throws MethodError d_dna32 == d_rna64
        @test length(d_rna64) == length(d_dna32)

        # Test with amino acids
        aa_d64 = dmer"KWOP"a
        aa_d128 = AAOligomer{UInt128}(aa_d64)
        @test_throws MethodError aa_d128 == aa_d64
        @test length(aa_d128) == length(aa_d64)
        aa_d256 = AAOligomer{UInt256}(aa_d64)
        @test_throws MethodError aa_d256 == aa_d64
        @test length(aa_d256) == length(aa_d64)

        # Test with empty kmers
        empty_d32 = DNAOligomer{UInt32}(dna"")
        empty_d64 = DNAOligomer{UInt64}(empty_d32)
        @test_throws MethodError empty_d64 == empty_d32
        @test length(empty_d64) == 0

        # Test with 4-bit alphabets
        s_4bit = LongSequence{DNAAlphabet{4}}(dna"TAGCN")
        d_4bit_32 = Oligomer{DNAAlphabet{4}, UInt32}(s_4bit)
        d_4bit_64 = Oligomer{DNAAlphabet{4}, UInt64}(d_4bit_32)
        @test_throws MethodError d_4bit_64 == d_4bit_32
        @test length(d_4bit_64) == length(d_4bit_32)

        # Narrowing must compare symbol capacities, not coding bits with symbols.
        for (source, target, fitting, overflowing) in [
                (DNAOligomer{UInt16}, DNAOligomer{UInt8}, dna"TAG", dna"TAGC"),
                (RNAOligomer{UInt16}, RNAOligomer{UInt8}, rna"UAG", rna"UAGC"),
            ]
            original = source(fitting)
            narrowed = target(original)
            @test string(narrowed) == string(original)
            @test length(narrowed) == length(original)
            @test_throws MethodError narrowed == original
            @test_throws Exception target(source(overflowing))
        end

        # AAOligomer{UInt8} has capacity zero, exercising the smallest target
        # backing type too.
        empty_aa = AAOligomer{UInt16}(aa"")
        @test_throws MethodError AAOligomer{UInt8}(empty_aa) == empty_aa
        @test_throws Exception AAOligomer{UInt8}(AAOligomer{UInt16}(aa"K"))
    end

    @testset "To Kmer" begin
        for s in [dna"TAGCTA", rna"UGCUGA", aa"PLKWM"]
            dkmer = Oligomer{typeof(Alphabet(s)), UInt64}(s)
            kmer = Kmer{typeof(Alphabet(s)), length(s)}(dkmer)
            @test length(kmer) == length(dkmer)
            @test string(dkmer) == string(kmer)
            @test_throws MethodError kmer == dkmer
            @test_throws MethodError dkmer == kmer
        end
    end

    @testset "To Kmer with N=2 (>64 coding bits)" begin
        # For 2-bit DNA: need >32 bases for >64 coding bits
        s_dna = dna"TAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT"  # 33 bases = 66 bits
        dkmer_dna = Oligomer{DNAAlphabet{2}, UInt256}(s_dna)
        kmer_dna = Kmer{DNAAlphabet{2}, 33}(dkmer_dna)
        @test length(kmer_dna) == length(dkmer_dna)
        @test string(dkmer_dna) == string(kmer_dna)

        # For 2-bit RNA
        s_rna = rna"UGCUGAUGCUGAUGCUGAUGCUGAUGCUGAUGA"  # 33 bases = 66 bits
        dkmer_rna = Oligomer{RNAAlphabet{2}, UInt256}(s_rna)
        kmer_rna = Kmer{RNAAlphabet{2}, 33}(dkmer_rna)
        @test length(kmer_rna) == length(dkmer_rna)
        @test string(dkmer_rna) == string(kmer_rna)

        # For 8-bit amino acids: need >8 bases for >64 coding bits
        s_aa = aa"KWOPPLKWM"  # 9 bases = 72 bits
        dkmer_aa = Oligomer{AminoAcidAlphabet, UInt256}(s_aa)
        kmer_aa = Kmer{AminoAcidAlphabet, 9}(dkmer_aa)
        @test length(kmer_aa) == length(dkmer_aa)
        @test string(dkmer_aa) == string(kmer_aa)

        # Test error on length mismatch
        dkmer = dmer"TAG"d
        @test_throws Exception Kmer{DNAAlphabet{2}, 5}(dkmer)
    end

    @testset "Capacity limits" begin
        # Test that exceeding capacity throws
        @test_throws Exception DNAOligomer{UInt32}(dna"T"^30)
        @test_throws Exception AAOligomer{UInt32}(aa"A"^8)
    end
end

@testset "Indexing and iteration" begin
    @testset "Scalar indexing" begin
        m = dmer"TAGCTGAC"d
        @test m[1] == DNA_T
        @test m[3] == DNA_G
        @test m[8] == DNA_C
        @test first(m) == DNA_T
        @test last(m) == DNA_C

        @test_throws BoundsError m[0]
        @test_throws BoundsError m[9]
    end

    @testset "Range indexing" begin
        m = dmer"TAGCTGAC"d
        @test m[1:3] == dmer"TAG"d
        @test m[2:5] == dmer"AGCT"d
        @test m[6:8] == dmer"GAC"d
        @test m[1:0] == dmer""d

        @test_throws BoundsError m[0:3]
        @test_throws BoundsError m[6:9]
    end

    @testset "Iteration" begin
        s = dna"TAGCTGAC"
        m = dmer"TAGCTGAC"d
        @test collect(m) == collect(s)
    end
end

@testset "Comparison and equality" begin
    @testset "Equality" begin
        @test dmer"TAG"d == dmer"TAG"d
        @test dmer"TAG"d != dmer"TAC"d
        @test dmer""d == dmer""d
    end

    @testset "Ordering" begin
        @test dmer"TAG"d < dmer"TGA"d
        @test dmer"AAA"d < dmer"AAC"d
        @test dmer"TAG"d > dmer"TAC"d
    end

    @testset "Comparison of different lengths" begin
        @test dmer"TAG"d < dmer"TAGA"d
        @test dmer"TAGA"d > dmer"TAG"d
    end

    @testset "cmp function" begin
        # Test cmp returns -1, 0, or 1
        @test cmp(dmer"TAG"d, dmer"TAG"d) == 0
        @test cmp(dmer"TAG"d, dmer"TGA"d) < 0
        @test cmp(dmer"TGA"d, dmer"TAG"d) > 0

        # Test with RNA
        @test cmp(dmer"UAG"r, dmer"UAG"r) == 0
        @test cmp(dmer"UAG"r, dmer"UGA"r) < 0

        # Test with different lengths
        @test cmp(dmer"TAG"d, dmer"TAGA"d) < 0
        @test cmp(dmer"TAGA"d, dmer"TAG"d) > 0
    end

    @testset "Comparison rejects different backing types" begin
        s1 = dna"TAG"
        s2 = dna"TAC"
        s3 = dna"TAGA"

        m32_1 = DNAOligomer{UInt32}(s1)
        m64_1 = DNAOligomer{UInt64}(s1)
        m32_2 = DNAOligomer{UInt32}(s2)
        m64_2 = DNAOligomer{UInt64}(s2)

        for (a, b) in [(m32_1, m64_1), (m32_2, m64_2)]
            @test_throws MethodError a == b
            @test_throws MethodError b == a
            @test_throws MethodError isequal(a, b)
            @test_throws MethodError a < b
            @test_throws MethodError cmp(a, b)
        end

        # Different lengths do not make the backing types comparable.
        m32_3 = DNAOligomer{UInt32}(s3)
        @test m32_1 < m32_3
        @test cmp(m32_1, m32_3) < 0
        @test_throws MethodError m64_1 < m32_3
        @test_throws MethodError cmp(m64_1, m32_3)

        # Also for AA types
        s1 = aa"KPRCRLF"
        s2 = aa"KPRCRLFAAA"

        m64_1 = AAOligomer{UInt64}(s1)
        m128_1 = AAOligomer{UInt128}(s1)
        m128_2 = AAOligomer{UInt128}(s2)
        m256_1 = AAOligomer{UInt256}(s1)

        @test_throws MethodError m64_1 == m128_1
        @test_throws MethodError cmp(m64_1, m128_1)
        @test m128_1 < m128_2
        @test_throws MethodError m64_1 == m256_1
        @test_throws MethodError cmp(m256_1, m128_1)
    end
end

@testset "Integer conversion" begin
    @testset "as_integer" begin
        m1 = dmer"TAG"d
        u1 = as_integer(m1)
        @test u1 isa Unsigned

        m2 = dmer"TAC"d
        u2 = as_integer(m2)
        @test u1 != u2

        # Empty kmer
        @test as_integer(dmer""d) == 0
    end

    @testset "from_integer" begin
        for s in [dna"TAG", dna"TAGCTGA", dna"ATGCTAGC"]
            m = DNAOligomer{UInt64}(s)
            u = as_integer(m)
            m2 = from_integer(typeof(m), u, length(m))
            @test m == m2
        end

        # Error on exceeding capacity
        @test_throws Exception from_integer(DNAOligomer{UInt32}, UInt32(0), 30)

        # A zero-length value must discard every input bit, including inputs
        # which would otherwise wrap a full-width shift to a zero-bit shift.
        for U in [UInt8, UInt16, UInt32, UInt64, UInt128, UInt256, UInt512]
            T = DNAOligomer{U}
            for value in [zero(U), one(U), typemax(U)]
                oligomer = from_integer(T, value, 0)
                @test length(oligomer) == 0
                @test isempty(oligomer)
                @test oligomer == empty(T)
                @test oligomer.x == zero(U)
            end
        end
    end

    @testset "Round-trip conversion" begin
        for s in [dna"TAGCTGA", rna"UGCUGA", aa"PLKWM"]
            m = Oligomer{typeof(Alphabet(s)), UInt64}(s)
            u = as_integer(m)
            m2 = from_integer(typeof(m), u, length(m))
            @test m === m2
        end
    end
end

@testset "Hashing" begin
    m1 = dmer"TAG"d
    m2 = dmer"TAG"d
    m3 = dmer"TAC"d

    # Same kmers hash to same value
    @test hash(m1) == hash(m2)

    # Different kmers likely hash to different values (not guaranteed but likely)
    @test hash(m1) != hash(m3)

    # Hash with seed
    h1 = hash(m1, UInt(123))
    h2 = hash(m1, UInt(456))
    @test h1 != h2

    # Length must be part of hash: TCA and TC have identical coding bits
    # (since A encodes to 00, which looks like padding), but different lengths.
    # They must hash differently despite having the same as_integer representation.
    m1 = dmer"TCA"d
    m2 = dmer"TC"d
    @test hash(m1) != hash(m2)

    @testset "fx_hash" begin
        m1 = dmer"TAG"d
        m2 = dmer"TAG"d
        m3 = dmer"TAC"d

        # Same kmers should produce same fx_hash
        @test Kmers.fx_hash(m1, UInt(0)) == Kmers.fx_hash(m2, UInt(0))

        # Different kmers should produce different fx_hash values
        @test Kmers.fx_hash(m1, UInt(0)) != Kmers.fx_hash(m3, UInt(0))

        # Different seeds should produce different hashes
        @test Kmers.fx_hash(m1, UInt(123)) != Kmers.fx_hash(m1, UInt(456))

        # Length must be part of hash: TCA and TC have identical coding bits
        # (since A encodes to 00, which looks like padding), but different lengths.
        # They must hash differently despite having the same as_integer representation.
        m1 = dmer"TCA"d
        m2 = dmer"TC"d
        @test Kmers.fx_hash(m1, UInt(0)) != Kmers.fx_hash(m2, UInt(0))
    end
end

@testset "Biological operations" begin
    @testset "Reverse" begin
        # Test with 2-bit DNA
        s = dna"TAGCTGA"
        m = dmer"TAGCTGA"d
        @test reverse(m) == DNAOligomer{UInt64}(reverse(s))

        # Test empty sequence
        @test reverse(dmer""d) == dmer""d

        # Test with 2-bit RNA
        s_rna = rna"UAGCUGA"
        m_rna = dmer"UAGCUGA"r
        @test reverse(m_rna) == RNAOligomer{UInt64}(reverse(s_rna))

        # Test with 4-bit DNA
        s_4bit = dna"TAGCTGA"
        m_4bit = Oligomer{DNAAlphabet{4}, UInt64}(s_4bit)
        @test reverse(m_4bit) == Oligomer{DNAAlphabet{4}, UInt64}(reverse(s_4bit))

        # Test with amino acids
        s_aa = aa"KWOP"
        m_aa = dmer"KWOP"a
        @test reverse(m_aa) == AAOligomer{UInt64}(reverse(s_aa))

        # Test with larger bit integers
        s_large = dna"TAGCTAGCTAGCTAGC"
        m_large = DNAOligomer{UInt256}(s_large)
        @test reverse(m_large) == DNAOligomer{UInt256}(reverse(s_large))
    end

    @testset "Complement" begin
        # Test with 2-bit DNA
        s = dna"TAGCTGA"
        m = dmer"TAGCTGA"d
        @test complement(m) == DNAOligomer{UInt64}(complement(s))

        # Test with 2-bit RNA
        s_rna = rna"UAGCUGA"
        m_rna = dmer"UAGCUGA"r
        @test complement(m_rna) == RNAOligomer{UInt64}(complement(s_rna))

        # Test with 4-bit DNA (includes ambiguous bases)
        s_4bit = LongSequence{DNAAlphabet{4}}(dna"TAGCTGAW")  # W = A or T
        m_4bit = Oligomer{DNAAlphabet{4}, UInt64}(s_4bit)
        @test complement(m_4bit) == Oligomer{DNAAlphabet{4}, UInt64}(complement(s_4bit))

        # Test with 4-bit RNA
        s_rna_4bit = LongSequence{RNAAlphabet{4}}(rna"UAGCUGAW")  # W = A or U
        m_rna_4bit = Oligomer{RNAAlphabet{4}, UInt64}(s_rna_4bit)
        @test complement(m_rna_4bit) == Oligomer{RNAAlphabet{4}, UInt64}(complement(s_rna_4bit))
    end

    @testset "Reverse complement" begin
        # Test with 2-bit DNA
        s = dna"TAGCTGA"
        m = dmer"TAGCTGA"d
        @test reverse_complement(m) == DNAOligomer{UInt64}(reverse_complement(s))

        # Test with 2-bit RNA
        s_rna = rna"UAGCUGA"
        m_rna = dmer"UAGCUGA"r
        @test reverse_complement(m_rna) == RNAOligomer{UInt64}(reverse_complement(s_rna))

        # Test with 4-bit DNA
        s_4bit = LongSequence{DNAAlphabet{4}}(dna"TAGCTGA")
        m_4bit = Oligomer{DNAAlphabet{4}, UInt64}(s_4bit)
        @test reverse_complement(m_4bit) == Oligomer{DNAAlphabet{4}, UInt64}(reverse_complement(s_4bit))

        # Test with 4-bit RNA
        s_rna_4bit = LongSequence{RNAAlphabet{4}}(rna"UAGCUGA")
        m_rna_4bit = Oligomer{RNAAlphabet{4}, UInt64}(s_rna_4bit)
        @test reverse_complement(m_rna_4bit) == Oligomer{RNAAlphabet{4}, UInt64}(reverse_complement(s_rna_4bit))

        # Test with larger bit integers
        s_large = dna"TAGCTAGCTAGCTAGCTAGC"
        m_large = DNAOligomer{UInt512}(s_large)
        @test reverse_complement(m_large) == DNAOligomer{UInt512}(reverse_complement(s_large))
    end

    @testset "Canonical" begin
        m1 = dmer"TAGCTGA"d
        m2 = dmer"TCAGCTA"d
        @test canonical(m1) == canonical(m2) == m1
        @test iscanonical(dmer"AATT"d)
        @test iscanonical(dmer"TTAA"d)
        @test iscanonical(empty(m1))
        @test !iscanonical(dmer"TGGA"d)

        # Test with larger bit integers
        m_large = DNAOligomer{UInt256}(dna"TAGCTAGCTAGC")
        m_rc = reverse_complement(m_large)
        @test canonical(m_large) == min(m_large, m_rc)
    end
end

@testset "Counting" begin
    @testset "Count GC" begin
        @test count(isGC, dmer"TATCGGAGA"d) == 4
        @test count(isGC, dmer"TATATATAAAAA"d) == 0
        @test count(isGC, dmer""d) == 0

        @test count(isGC, dmer"AUGUCGUAG"r) == 4
    end

    @testset "Count symbols" begin
        m = dmer"TAGCTGA"d
        @test count(==(DNA_A), m) == 2
        @test count(==(DNA_T), m) == 2
        @test count(==(DNA_G), m) == 2
        @test count(==(DNA_C), m) == 1

        m = AAOligomer{UInt256}(aa"WLAKWVMARQKW")
        @test count(==(AA_W), m) == 3
        @test count(==(AA_Q), m) == 1
        @test count(==(AA_A), m) == 2
        @test count(==(AA_C), m) == 0

        # Test with even larger integers
        m_large = DNAOligomer{UInt512}(dna"TAGCTAGCTAGCTAGC")
        @test count(==(DNA_T), m_large) == 4
        @test count(==(DNA_A), m_large) == 4

        # Test all standard alphabets and backing widths against explicit
        # iteration. This covers zero encodings, absent symbols, empty values,
        # and values that fill their backing type to capacity.
        for A in (
                    DNAAlphabet{2},
                    RNAAlphabet{2},
                    DNAAlphabet{4},
                    RNAAlphabet{4},
                    AminoAcidAlphabet,
                ), U in (UInt8, UInt16, UInt32, UInt64, UInt128, UInt256, UInt512)
            T = Oligomer{A, U}
            alphabet_symbols = collect(symbols(A()))
            sequence = [
                alphabet_symbols[mod1(i, length(alphabet_symbols))] for
                    i in 1:capacity(T)
            ]
            oligomer = T(sequence)
            @test length(oligomer) == capacity(T)

            for symbol in alphabet_symbols
                @test count(==(symbol), oligomer) == count(==(symbol), collect(oligomer))
                @test count(==(symbol), empty(T)) == 0
            end
        end
    end
end

@testset "shift_encoding" begin
    m = DNAOligomer{UInt32}(dna"TAGA")
    enc = UInt32(BioSequences.encode(DNAAlphabet{2}(), DNA_C))
    m2 = Kmers.shift_encoding(m, enc)
    @test m2 == DNAOligomer{UInt32}(dna"AGAC")
    @test length(m2) == length(m)
end

@testset "Mixed integer types" begin
    for U in [UInt8, UInt16, UInt32, UInt64, UInt128, UInt256, UInt512]
        s = dna"TAG"
        m = Oligomer{DNAAlphabet{2}, U}(s)
        @test string(m) == string(s)
        @test_throws MethodError m == s
        @test_throws MethodError s == m
        @test length(m) == 3
    end
end

@testset "push and push_first" begin
    for (dkmer_fn, longseq_fn) in Any[[push, push!], [push_first, pushfirst!]]
        for (dkmer, symbol) in Any[
                (dmer"TAGC"d, DNA_A),
                (dmer"AUGC"r, RNA_U),
                (dmer""d, DNA_T),
                (dmer""r, RNA_G),
                (DNAOligomer{UInt32}(dna"TAG"), DNA_G),
                (RNAOligomer{UInt32}(rna"AUG"), RNA_C),
                (AAOligomer{UInt64}(aa"KWOP"), AA_L),
                (AAOligomer{UInt256}(aa"MWP"), AA_K),
            ]
            # Apply operation to Oligomer
            result_dkmer = dkmer_fn(dkmer, symbol)

            # Apply operation to LongSequence for comparison
            longseq = LongSequence{typeof(Alphabet(dkmer))}(dkmer)
            longseq_fn(longseq, symbol)

            # Test content against the LongSequence reference.
            @test string(result_dkmer) == string(longseq)
            @test length(result_dkmer) == length(longseq)
            @test length(result_dkmer) == length(dkmer) + 1

            # Test that result has correct type
            @test result_dkmer isa typeof(dkmer)

            # Test that original is unchanged (immutability)
            @test length(dkmer) == length(result_dkmer) - 1
        end
    end

    # Test type conversion (pushing DNA to RNA kmer and vice versa)
    @test push(dmer"ATGC"d, RNA_U) == dmer"ATGCT"d
    @test push(dmer"AUGC"r, DNA_T) == dmer"AUGCU"r
    @test push_first(dmer"ATGC"d, RNA_U) == dmer"TATGC"d
    @test push_first(dmer"AUGC"r, DNA_T) == dmer"UAUGC"r

    # Test character conversion
    @test push(dmer"TAG"d, 'C') == dmer"TAGC"d
    @test push_first(dmer"TAG"d, 'A') == dmer"ATAG"d

    # Test chaining operations
    @test push(push(dmer"TAG"d, DNA_C), DNA_A) == dmer"TAGCA"d
    @test push_first(push_first(dmer"TAG"d, DNA_C), DNA_A) == dmer"ACTAG"d

    # Test error when kmer is at max capacity
    # For UInt32 with 2-bit DNA: capacity = 14
    d_full = DNAOligomer{UInt32}(dna"T"^14)
    @test_throws BoundsError push(d_full, DNA_A)
    @test_throws BoundsError push_first(d_full, DNA_A)

    # For UInt64 with 2-bit DNA: capacity = 29
    d64_full = DNAOligomer{UInt64}(dna"A"^29)
    @test_throws BoundsError push(d64_full, DNA_T)
    @test_throws BoundsError push_first(d64_full, DNA_G)

    # For UInt32 with 8-bit AA: capacity = 3
    aa32_full = AAOligomer{UInt32}(aa"WPK")
    @test_throws BoundsError push(aa32_full, AA_L)
    @test_throws BoundsError push_first(aa32_full, AA_M)

    # UInt8 has just enough length bits for a full 2-bit oligomer. Incrementing
    # its packed representation first would wrap the length and corrupt it.
    for (T, sequence, symbol) in [
            (DNAOligomer{UInt8}, dna"TAG", DNA_A),
            (RNAOligomer{UInt8}, rna"UAG", RNA_U),
        ]
        full = T(sequence)
        original = full.x
        @test length(full) == capacity(T)
        @test_throws BoundsError push(full, symbol)
        @test_throws BoundsError push_first(full, symbol)
        @test full.x == original
    end

    # Verify we can push to capacity-1 without error
    d_almost_full = DNAOligomer{UInt64}(dna"T"^28)
    @test length(push(d_almost_full, DNA_A)) == 29
    @test length(push_first(d_almost_full, DNA_G)) == 29

    # Test with larger bit integers
    d256 = DNAOligomer{UInt256}(dna"TAGC")
    @test push(d256, DNA_G) == DNAOligomer{UInt256}(dna"TAGCG")
    @test push_first(d256, DNA_A) == DNAOligomer{UInt256}(dna"ATAGC")
end

@testset "pop and pop_first" begin
    for (dkmer_fn, longseq_fn) in Any[[pop, pop!], [pop_first, popfirst!]]
        for dkmer in Any[
                dmer"TAGCA"d,
                dmer"AUGCU"r,
                dmer"T"d,
                dmer"A"r,
                DNAOligomer{UInt32}(dna"TAGG"),
                RNAOligomer{UInt32}(rna"AUGC"),
                AAOligomer{UInt64}(aa"KWOPL"),
                AAOligomer{UInt128}(aa"MWPK"),
            ]
            # Apply operation to Oligomer
            result_dkmer = dkmer_fn(dkmer)

            # Apply operation to LongSequence for comparison
            longseq = LongSequence{typeof(Alphabet(dkmer))}(dkmer)
            longseq_fn(longseq)

            # Test content against the LongSequence reference.
            @test string(result_dkmer) == string(longseq)
            @test length(result_dkmer) == length(longseq)
            @test length(result_dkmer) == length(dkmer) - 1

            # Test that result has correct type
            @test result_dkmer isa typeof(dkmer)

            # Test that original is unchanged (immutability)
            @test length(dkmer) == length(result_dkmer) + 1
        end
    end

    # Test specific sequences to verify correctness
    @test pop(dmer"TAGCA"d) == dmer"TAGC"d
    @test pop(dmer"AUGCU"r) == dmer"AUGC"r
    @test pop_first(dmer"TAGCA"d) == dmer"AGCA"d
    @test pop_first(dmer"AUGCU"r) == dmer"UGCU"r

    # Test chaining operations
    @test pop(pop(dmer"TAGCA"d)) == dmer"TAG"d
    @test pop_first(pop_first(dmer"TAGCA"d)) == dmer"GCA"d

    # Test popping down to empty
    @test pop(dmer"T"d) == dmer""d
    @test pop_first(dmer"A"r) == dmer""r

    # Test error when popping empty kmer
    @test_throws BoundsError pop(dmer""d)
    @test_throws BoundsError pop_first(dmer""d)
    @test_throws BoundsError pop(AAOligomer{UInt64}(aa""))
    @test_throws BoundsError pop_first(AAOligomer{UInt128}(aa""))

    # Test with larger bit integers
    d512 = DNAOligomer{UInt512}(dna"TAGCA")
    @test pop(d512) == DNAOligomer{UInt512}(dna"TAGC")
    @test pop_first(d512) == DNAOligomer{UInt512}(dna"AGCA")
end

@testset "setindex" begin
    # Basic functionality
    d = dmer"TAGCTGA"d
    @test Base.setindex(d, DNA_C, 1) == dmer"CAGCTGA"d
    @test Base.setindex(d, DNA_A, 4) == dmer"TAGATGA"d
    @test Base.setindex(d, DNA_T, 7) == dmer"TAGCTGT"d
    @test d == dmer"TAGCTGA"d  # Original unchanged

    # Different alphabets and backing types
    @test Base.setindex(dmer"AUGC"r, RNA_G, 2) == dmer"AGGC"r
    @test Base.setindex(dmer"KWOP"a, AA_L, 3) == dmer"KWLP"a
    @test Base.setindex(DNAOligomer{UInt32}(dna"ATGC"), DNA_G, 2) == DNAOligomer{UInt32}(dna"AGGC")

    # Type conversion
    @test Base.setindex(dmer"TAG"d, 'C', 2) == dmer"TCG"d

    # Bounds checking
    @test_throws BoundsError Base.setindex(dmer"TAGC"d, DNA_A, 0)
    @test_throws BoundsError Base.setindex(dmer"TAGC"d, DNA_A, 5)
    @test_throws BoundsError Base.setindex(dmer""d, DNA_A, 1)

    # Test with larger bit integers
    d256 = DNAOligomer{UInt256}(dna"TAGCAT")
    @test Base.setindex(d256, DNA_G, 3) == DNAOligomer{UInt256}(dna"TAGCAT")
    @test Base.setindex(d256, DNA_C, 1) == DNAOligomer{UInt256}(dna"CAGCAT")
end

@testset "capacity" begin
    # Create a zero-BPS alphabet for testing
    struct ZeroBPSAlphabet <: Alphabet end
    Base.eltype(::Type{ZeroBPSAlphabet}) = DNA
    BioSequences.BitsPerSymbol(::ZeroBPSAlphabet) = BioSequences.BitsPerSymbol{0}()

    # Test zero BPS: capacity should be clamped to typemax(Int)
    @test capacity(Oligomer{ZeroBPSAlphabet, UInt8}) == clamp(typemax(UInt8), Int)
    @test capacity(Oligomer{ZeroBPSAlphabet, UInt32}) == clamp(typemax(UInt32), Int)
    @test capacity(Oligomer{ZeroBPSAlphabet, UInt128}) == typemax(Int)

    # Zero-BPS values have no coding bits, so from_integer must retain only the
    # requested length regardless of its integer input.
    zero_bps_type = Oligomer{ZeroBPSAlphabet, UInt8}
    for value in [zero(UInt8), one(UInt8), typemax(UInt8)]
        @test from_integer(zero_bps_type, value, 0).x == 0x00
    end
    @test from_integer(zero_bps_type, typemax(UInt8), 1).x == 0x01

    # Test non-zero BPS: capacity should be in range 0:div(8 * sizeof(U), B)
    for (A, bps) in [(DNAAlphabet{2}, 2), (DNAAlphabet{4}, 4), (AminoAcidAlphabet, 8)]
        for U in [UInt8, UInt32, UInt64]
            cap = capacity(Oligomer{A, U})
            max_possible = div(8 * sizeof(U), bps)
            @test 0 <= cap <= max_possible
        end
    end

    # Test that larger backing type gives larger or equal capacity
    @test capacity(DNAOligomer{UInt32}) <= capacity(DNAOligomer{UInt64})
    @test capacity(AAOligomer{UInt32}) <= capacity(AAOligomer{UInt64})
    @test capacity(DNAOligomer{UInt64}) <= capacity(DNAOligomer{UInt256})
    @test capacity(AAOligomer{UInt128}) <= capacity(AAOligomer{UInt512})
end

@testset "translate" begin
    # Basic translation - 2-bit alphabets
    @testset "Basic 2-bit" begin
        # DNA 2-bit with different backing types
        @test string(translate(dmer"ATGTAA"d)) == "M*"
        @test string(translate(DNAOligomer{UInt32}(dna"ATGTAA"))) == "M*"
        @test string(translate(DNAOligomer{UInt64}(dna"ATGTAA"))) == "M*"

        # RNA 2-bit with different backing types
        @test string(translate(dmer"AUGUAA"r)) == "M*"
        @test string(translate(RNAOligomer{UInt32}(rna"AUGUAA"))) == "M*"
        @test string(translate(RNAOligomer{UInt64}(rna"AUGUAA"))) == "M*"

        # Longer sequences
        @test string(translate(dmer"TCTACACCCTAG"d)) == "STP*"
        @test string(translate(dmer"UCUACACCCUAG"r)) == "STP*"
        s = randdnaseq(249)
        t = translate(DNAOligomer{UInt512}(s))
        @test LongSequence(t) == translate(s)
    end

    # Basic translation - 4-bit alphabets
    @testset "Basic 4-bit" begin
        # DNA 4-bit with different backing types
        d = Oligomer{DNAAlphabet{4}, UInt64}(dna"TGGCCCGATTGA")
        @test string(translate(d)) == "WPD*"

        # UInt128 works with 4-bit since capacity is smaller
        d128 = Oligomer{DNAAlphabet{4}, UInt128}(dna"ATGTAA")
        @test string(translate(d128)) == "M*"

        # RNA 4-bit with different backing types
        r = Oligomer{RNAAlphabet{4}, UInt64}(rna"UGGCCCGAUUGA")
        @test string(translate(r)) == "WPD*"

        r32 = Oligomer{RNAAlphabet{4}, UInt32}(rna"AUGUAA")
        @test string(translate(r32)) == "M*"

        r128 = Oligomer{RNAAlphabet{4}, UInt128}(rna"AUGUAA")
        @test string(translate(r128)) == "M*"

        # Long sequence
        s = randdnaseq(252)
        t = translate(Oligomer{DNAAlphabet{4}, UInt1024}(s))
        @test LongSequence(t) == translate(s)
    end

    # Different genetic codes
    @testset "Different genetic codes" begin
        # Vertebrate mitochondrial code: AGA and AGG are stop codons
        vert_mito = BioSequences.vertebrate_mitochondrial_genetic_code
        @test string(translate(dmer"ATGAGA"d; code = vert_mito)) == "M*"  # AGA is stop

        # Standard code: AGA and AGG code for R (Arginine)
        @test string(translate(dmer"ATGAGA"d)) == "MR"
        @test string(translate(dmer"ATGAGG"d)) == "MR"

        # Test with different backing types
        @test string(translate(DNAOligomer{UInt64}(dna"ATGAGA"); code = vert_mito)) == "M*"
        @test string(translate(RNAOligomer{UInt32}(rna"AUGAGA"); code = vert_mito)) == "M*"

        # Test with 4-bit alphabets
        @test string(translate(Oligomer{DNAAlphabet{4}, UInt64}(dna"ATGAGA"); code = vert_mito)) == "M*"
    end

    # alternative_start flag
    @testset "alternative_start" begin
        # Without alternative_start: TTG codes for L (Leucine)
        @test string(translate(dmer"TTGCCC"d; alternative_start = false)) == "LP"
        @test string(translate(dmer"UUGCCC"r; alternative_start = false)) == "LP"

        # With alternative_start: first codon becomes M regardless
        @test string(translate(dmer"TTGCCC"d; alternative_start = true)) == "MP"
        @test string(translate(dmer"UUGCCC"r; alternative_start = true)) == "MP"
    end

    # Ambiguous codons (only for 4-bit alphabets)
    @testset "Ambiguous codons" begin
        # With allow_ambiguous_codons=true (default), ambiguous codons translate
        d_ambig = Oligomer{DNAAlphabet{4}, UInt128}("AAAACWGCSWTARACADA")
        @test string(translate(d_ambig)) == "KTAJBX"

        # With allow_ambiguous_codons=false, ambiguous codons throw
        @test_throws Exception translate(d_ambig; allow_ambiguous_codons = false)

        # Test various ambiguous nucleotides
        # W = A or T, so TWG could be AAG (K) or TAG (*)
        # With allow_ambiguous, should give X (ambiguous)
        d_w = Oligomer{DNAAlphabet{4}, UInt64}(dna"ATGTWG")
        result_w = translate(d_w; allow_ambiguous_codons = true)
        @test length(result_w) == 2  # M and something
    end

    # Error: Length not divisible by 3
    @testset "Length not divisible by 3" begin
        @test_throws ArgumentError translate(dmer"A"d)
        @test_throws ArgumentError translate(dmer"UG"r)
        @test_throws ArgumentError translate(dmer"CUGUAGUUGUCGC"r)
        @test_throws ArgumentError translate(Oligomer{RNAAlphabet{4}, UInt32}(rna"AUGC"))
    end

    # Error: Sequences with gaps (only for 4-bit alphabets)
    @testset "Sequences with gaps" begin
        @test_throws Exception translate(Oligomer{RNAAlphabet{4}, UInt64}(rna"-UGAUG"))
        @test_throws Exception translate(Oligomer{DNAAlphabet{4}, UInt64}(dna"AT-ATG"))
        @test_throws Exception translate(Oligomer{RNAAlphabet{4}, UInt64}(rna"AUGAU-"))
        @test_throws Exception translate(Oligomer{DNAAlphabet{4}, UInt64}(dna"A--"))
    end

    # Error: Input type capacity too large for output
    @testset "Capacity overflow" begin
        # These tests have loaded BitIntegers, so the maximum backing integer
        # is currently UInt1024.
        # This means 1024-bit 2-bit alphabets overflow, but not 4-bit alphabets.
        @test_throws ErrorException translate(DNAOligomer{UInt1024}(dna"ATG"))
        @test_throws ErrorException translate(RNAOligomer{UInt1024}(rna"AUG"))

        # Even empty sequences should error due to type-based capacity check
        @test_throws ErrorException translate(DNAOligomer{UInt1024}(dna""))
        @test_throws ErrorException translate(RNAOligomer{UInt1024}(rna""))

        # Not an overflow - 512 bits fit
        @test string(translate(RNAOligomer{UInt512}(rna""))) == ""

        # 4-bit nucleotides do not overflow even when 1024 bits
        @test string(translate(Oligomer{RNAAlphabet{4}, UInt1024}(rna"AUG"))) == "M"
    end

    # Edge cases
    @testset "Edge cases" begin
        # Empty sequence (length 0 is divisible by 3, produces 0 AA)
        for A in (DNAAlphabet{2}, RNAAlphabet{2}, DNAAlphabet{4}, RNAAlphabet{4}),
                U in (UInt8, UInt32), alternative_start in (false, true)
            input = Oligomer{A, U}("")
            result = translate(input; alternative_start)
            empty_result = empty(typeof(result))
            @test length(result) == 0
            @test isempty(result)
            @test result == empty_result
            @test hash(result) == hash(empty_result)
            @test string(result) == ""
            @test result.x == zero(result.x)
        end

        # Very long sequence (near capacity) - 2-bit
        # For UInt64 with 2-bit DNA, capacity is 31, so 30 nt = 10 AA
        long_dna = dmer"ATGATGATGATGATGATGATGATGATG"d
        @test length(translate(long_dna)) == 9
        @test LongSequence(translate(long_dna)) == aa"MMMMMMMMM"
    end
end

@testset "Misc" begin
    d = AAOligomer{UInt32}("WPK")
    @test only([d]') === d
end
