### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ 58f1747c-e2a1-11eb-2c26-95ef11565850
using GaloisFields,  LinearAlgebraX,  SHA

# ╔═╡ c7472820-0a79-4156-8d0e-e4367868c57e
using Base.Iterators, Base.Threads, Printf

# ╔═╡ 8ad81282-7a8e-4769-81d4-bb7fd3b7b1a5
md"
# List hash as matrices over finite fields

In my [previous post](https://blog.infogulch.com/merklist/2021/07/07/Merklist.html), we explored a definition for the hash of a list: hash each element, interpret the hash digest as a matrix with 8-bit unsigned int elements, and define the hash of the whole list to be the result of multiplying the matrices together. We found that this has some neat properties and works for a small number of elements, but  we also found that this definition is deficient because after enough entries are multiplied the list hash degenerates into the zero matrix.

If there is a primary reason this construction of the idea doesn't work, it's because the element hashes interpreted as matrices aren't invertible. The chance that a random matrix with integer mod 256 elements is invertible is equivalent to the chance of every element being odd; with 64 elements and a 50% chance per element, that's an overall $\frac{1}{2^{64}}$ chance. Meaning it's *possible* but very unlikely, and definitely not computationally efficient to find.

However, if we choose a different element type this degeneracy problem can be solved. Here we experiment with using [finite field](https://en.wikipedia.org/wiki/Finite_field) elements instead -- specifically $GF(256)$ which provides a $\frac{255}{256}$ chance of being invertible (much better!). At the end we also observe what in my opinion is the biggest benefit of this kind of definition for the hash of a list.
"

# ╔═╡ d8ab88dc-2278-42a7-b6e0-2dd67aad9c0f
md"
## Galois Fields

Finite fields aka Galois Fields have existed for a long time and binary fields in particular (like $GF(256)$ we use here) are used extensively by [cryptographic routines like AES](https://www.cs.uaf.edu/2015/spring/cs463/lecture/03_23_AES.html), and some CPUs even have [dedicated instructions](https://en.wikipedia.org/wiki/CLMUL_instruction_set) to enable fast math operations in the field.

For a number-theory description of galois fields, this Cryptography StackOverflow [answer to 'Galois fields in cryptography'](https://crypto.stackexchange.com/a/2718/6556) gives a great overview.

For an on-the-metal description of how math operations $+-\div×$ work in $GF(256)$, this page seems to cover it nicely: [samiam.org/galois.html](https://samiam.org/galois.html)


"

# ╔═╡ d573e0aa-4110-4d83-9ef3-68addd3ff9b8
md"
## Code

The code in this post is written in [Julia](https://julialang.org/) using a computational notebook [Pluto.jl](https://plutojl.org/).

The main packages that we'll use are [GaloisFields.jl](https://github.com/tkluck/GaloisFields.jl) which supports defining a much wider range of fields than what we use here, and [LinearAlgebraX.jl](https://github.com/scheinerman/LinearAlgebraX.jl) which defines functions that operate on exact/finite matrices (as opposed to Real/float).
"

# ╔═╡ 0ad5b059-9539-437b-97f5-e73ec1c06419
md"
First lets define the field we'll be using today:
"

# ╔═╡ b5fc120c-7727-4b53-bdc5-5352cd96dac1
const G = @GaloisField! 2^8 β

# ╔═╡ ef00d3c8-6ab4-4f86-ae6b-9fca3d71d5cb
md"## Quick tour of the field 𝔽₂₅₆"

# ╔═╡ e7621cea-c449-4b43-bd50-ba4646e2b48a
md"
Typical use of fields offered by this library includes exponentiating $β$ and performing operations on the results.
"

# ╔═╡ 9b964dbc-c78b-4246-bf7a-0761b9c97bc7
β^10

# ╔═╡ 52fcf820-cddd-4e36-86e9-afcb32f09403
β^10 * (β^10 + β^6)

# ╔═╡ 62e4fa8b-7fb7-4904-963f-a1b9c4176dc7
md"
The most notable feature of finite fields is that the results of all operations are also members of the field.
"

# ╔═╡ 710eef31-f07d-41f0-b410-dc6727f3589e
β^256

# ╔═╡ d9855f4f-c6d9-4bcf-ac25-e7f43b8b443b
md"
This is fun to play with but we need to convert bytes into field values and back. So we define some (ill-advised) conversion functions:
"

# ╔═╡ cec1b298-bf51-4500-97ae-5189a2049286
begin # conversions between G <-> UInt8 and Matrix{G} <-> Matrix{UInt8}

Base.convert(::Type{UInt8}, gf::G) = gf.n
Base.convert(::Type{G}, i::UInt8) = G(GaloisFields.Bits(), i)

display(x::G) = convert(UInt8, x)
display(x::Matrix{G}) = convert(Matrix{UInt8}, x)
	
end;

# ╔═╡ 89c9ca94-a428-442e-8dee-3a7b1b4f7338
md"
> Note, the author of the GaloisFields package [recommended against](https://github.com/tkluck/GaloisFields.jl/issues/17#issuecomment-879399545) implementing these conversion definitions as overrides on `Base.convert` because they can cause problems if you're not careful. We use them here in lieu of a better built-in solution.
"

# ╔═╡ f3e8221a-8d34-4442-a58f-9708952cb7b3
md"Byte values successfully round-trip through these conversions:"

# ╔═╡ 5bb70b21-5cb4-4451-8229-3d806e611408
display(convert(G, 0x11))

# ╔═╡ fb8b6da2-c515-4afa-b31c-8bd3259cb0f7
md"
By the way, you can also call functions that take a single parameter like this, which I use throughout this post:
"

# ╔═╡ 95a21d02-68d2-477f-b98d-1b83d0da1cdb
convert(G, 0x11) |> display

# ╔═╡ 0c0f5ee4-a3e0-49df-aafd-903a4eaca2d0
md"Addition is actually XOR:"

# ╔═╡ 523979df-d971-48cf-a868-174a4ebc644a
convert(G, 0x18) + convert(G, 0x08) |> display

# ╔═╡ a76ee3eb-a025-46f9-b4d8-1a9f47ab9c18
md"And multiplication is... lets say not what you learned in high school algebra. That's ok though, it still offers the properties of multiplication that we need."

# ╔═╡ c80af06b-3564-4e49-ae18-fcaf70cbc811
convert(G, 0x20) * convert(G, 0x08) |> display

# ╔═╡ fec99716-09a8-4cd5-80de-b57e7eb90540
md"Another notable property of finite fields is that every (nonzero) element has an inverse, which is not true of the ring of integers mod 256:"

# ╔═╡ 40c3a591-1b86-4abe-b7e8-5467f7f173a1
inv(convert(G, 0x08)) |> display

# ╔═╡ 3312903b-af7e-4088-bab9-f0e4c42c25fc
md"
## Deriving matrix entries from bytes

We have a new requirement compared to the last version: the matrix returned by `mk_entry` must be [invertible](https://en.wikipedia.org/wiki/Invertible_matrix).
"

# ╔═╡ 1d7a824a-73c5-4e21-9b2b-6b4923ea49c7
md"
`mk_entry` uses the `sha2_512` hash of `data` to derive an invertible 8×8 Matrix{𝔽₂₅₆}. Every candidate matrix hash has a 1/256 chance of not being invertible, so if the invertibility check fails a new matrix hash based on `data` is calculated and attempted again. Each round is calculated by prefixing `data` with a byte indicating the attempt number, as in:

```julia
# note: * is concatenation in julia, because monoid
digest_round1 = sha2_512([UInt8(1)] * data)
digest_round2 = sha2_512([UInt8(2)] * data)
digest_round3 = sha2_512([UInt8(3)] * data)
...
```

The first invertible digest in this sequence is taken to be the matrix hash of `data`. It is very rare that more than a couple rounds are needed to find an invertible matrix.

`mk_entry` is defined below.
"

# ╔═╡ 786df805-b24b-482d-a842-1756a1c38510
begin
	
# keep count of rejected candidates by round number because I'm curious
mk_entry_rejections = zeros(Int, 64)

function mk_entry(data::SHA.AbstractBytes)
	for i in 1:64
		ctx = SHA2_512_CTX()
		update!(ctx, [UInt8(i)])
		update!(ctx, data)
		digest = digest!(ctx)
		mat = convert(Matrix{G},reshape(digest, 8, 8))
		if detx(mat) != 0
			return mat
		end
		mk_entry_rejections[i] += 1 # note this rejection in the bucket for this round
	end
end

# Future: Use a different rejection design based on blake3 extended hash output?
# https://github.com/pixel27/Blake3Hash.jl

mk_entry(data::AbstractString) = mk_entry(String(data))
mk_entry(data::String) = mk_entry(codeunits(data))

end;

# ╔═╡ 533b55c8-cbfc-4ab8-bb3c-3813e5b58bea
md"
So what do entries look like?
"

# ╔═╡ 38f2b3fc-5078-4ee2-86f2-19d99c254166
mk_entry("Hello world!")

# ╔═╡ 7b2c1def-cdc7-4102-a66d-5a914a1b1eb0
md"
Hm, that's not very easy to read or compare, so we can convert entries back to UInt8 for display:
"

# ╔═╡ 8b251bfd-f61a-464f-a1b6-458dde096e14
mk_entry("Hello world!") |> display

# ╔═╡ c58360f7-010d-459e-a9ab-20469f0a98fc
md"
## A tour of matrix multiplication properties

We are using the matrix multiplication operation because it exhibits two desirable properties in particular:

* Associativity
* Non-Commutativity

This is an unusual combination of properties for an operation:

|     | associative | commutative |
| --- | ---         | ---         |
| $a+b$   | ✅ | ✅ |
| $a×b$   | ✅ | ✅ |
| $a-b$   | ❌ | ❌ |
| $a\div b$   | ❌ | ❌ |
| $a^b$ | ❌ | ❌ |
| $M×N$ | ✅ | ❌ |

Lets explore this with some simple examples for Merklist matrices. Here are a couple entries we can use to demonstrate.
"

# ╔═╡ 3f6fefd6-82c7-4f1f-b4c2-406faf39d0d2
a = mk_entry("A"); b = mk_entry("B"); c = mk_entry("C");

# ╔═╡ 7feb5e3c-1a71-432e-896a-263c8fff1815
md"
### Associativity

If matrix multiplication is associative, then calculating the result of $a*b*c$ as $(a*b)*c$ vs $a*(b*c)$ should result in the same matrix, which you can see below:
"

# ╔═╡ 80028718-f19c-492a-9a0e-09c5819dd633
[ (a*b)*c, a*(b*c) ] .|> display

# ╔═╡ 35bb2711-b8a9-4561-9e38-dc070b35d5dc
md"
**This is the key property that will unlock use-cases that have gone unaddressed due to the design of hash functions up to this point.**
"

# ╔═╡ a82c8cc1-e69c-4c19-bf40-e6efeadb00de
md"
### Non-Commutativity

If matrix multiplication is not commutative, then multiplying a matrix from the left ($a*b$) vs from the right ($b*a$) *should* cause the result to be different, which you can see is the case below:
"

# ╔═╡ 40a13273-b40f-4287-8b73-e1a8bcf119e2
[ a*b, b*a ] .|> display

# ╔═╡ 57f46d87-ca7a-4618-90f3-5a9ea8655bb2
md"
We *want* non-commutativity, because otherwise the result wouldn't describe a list.  Elements in a list have a defined order, that's what makes it a list and not something else like a set.
"

# ╔═╡ 201e4088-4481-4bdc-8c79-4a4fa5f23d61
md"
## Finding the hash of a large list
"

# ╔═╡ 87bedc95-8423-4be9-b25c-08e775ffb540
count = 1000000;

# ╔═╡ d9c3ebda-3678-4421-980e-0efd24f8a8f9
md"
Now, invertible matrices over finite fields form a [group](https://en.wikipedia.org/wiki/Group_(mathematics)), meaning that we should never run into any degenerate cases like a zero matrix. But the minimum we might do to validate this emperically is test it with a list that has a lot of elements. Lets try $count.
"

# ╔═╡ a8a7559f-2f9f-4f30-b8b8-f964b59da354
md"#### Helper definitions"

# ╔═╡ 85d19b91-c0a2-43be-abbb-28a7387195a6
# format_int is shorthand for formatting the digits of an an int into a string
format_int(x::Int) = @sprintf("%d", x);

# ╔═╡ 258da454-caf5-46c9-a45a-4df2fd44bf02
# chunk_range partitions a UnitRange into `n` ranges, with leftover elements spread 
# among the first chunks
function chunk_range(r::UnitRange{}, n::Int)
	# this is why I prefer 1-based indexing, c'est la vie
	chunk_size = (r.stop-r.start+1) ÷ n
	if chunk_size * n < (r.stop-r.start+1)
		chunk_size += 1
	end
	collect(partition(r, chunk_size))
end;

# ╔═╡ f4ca59b9-ca2f-4458-950a-b87dee7e3115
# zero is the zero-matrix over G, and is just used for initialization
zero = zeros(G, 8, 8);

# ╔═╡ 3d7c33df-74fb-49a8-98e2-c32456c1bf76
md"
### Single thread 
"

# ╔═╡ f8bc7a62-3b4a-47f7-936f-f5370a6efe97
# Normal left-to-right reduction
digest = reduce(*, 1:count .|> format_int .|> mk_entry); digest |> display

# ╔═╡ 67ee30f8-240d-4cfd-b03f-ca620790e0b0
md"
You might read this as:

* `1:count` - take the range of numbers from `1` to `count` (1-$count)
* `.|> format_int` - for each number convert it into a string
* `.|> mk_entry` - for each string element generate its matrix hash
* `reduce(*, ...)` - reduce all the matrix elements by `*` / matrix multiplication from left to right
* `digest = ...;` - assign the result to the variable `digest`
* `digest |> display` - display `digest` in a readable format
"

# ╔═╡ 786a6537-296e-40d6-830f-93ac232e2cf7
md"### Parallel reduction"

# ╔═╡ ceb7e3f7-031f-4b9e-9fa8-7678924db12c
md"
If the reduction is associative, then we should be able to break this calculation up into several chunks, calculate the hash of each chunk separately, then combine the hashes and get the same result as the single-threaded case.
"

# ╔═╡ 8e3a5919-d825-48bf-bde7-0a1cffca46b6
# Chunk up the range and reduce in multiple threads, then join the results
begin
	pdigest_partial = repeat([zero], Threads.nthreads())
	
	@threads for chunk in chunk_range(1:count, Threads.nthreads())
		rd = reduce(*, chunk .|> format_int .|> mk_entry)
		pdigest_partial[threadid()] = rd
	end
	
	pdigest = reduce(*, pdigest_partial)
	pdigest |> display
end

# ╔═╡ 1167b981-6baa-40f6-934a-928eea39ee1a
md"
You can't see from rendering this notebook as HTML, but the second one completes 2-3x faster on my machine. Neat! And yes they are identical:
"

# ╔═╡ 8dca48af-89c2-4528-a976-fd3c9992e2c2
digest == pdigest

# ╔═╡ 57adf406-ae46-41d9-a9ad-0c5c66489c12
md"
#### Analyzing rejected candidate matrices
By inspecting the rejections array we can see that we got a fair number of rejected candidate matrices in `mk_entry()`, and even some that were rejected more than once. This isn't exploitable, since the chance of finding a series of non-invertible matrices 64 times in a row is low enough as to be impossible to find. Note that the ratio of any pair of adjacent elements in this array should be approx. $1/256$. This rejection sampling approach may leak some pesky timing information, but that's well beyond the scope of this post.
"

# ╔═╡ be3d561c-75ab-4e17-91b6-a23921ad3aa4
mk_entry_rejections

# ╔═╡ a9283164-7b8d-40e0-9a3a-ff26fcc477d0
md"
## What's next?

First, I'd like to validate that the list hashes aren't susceptible to preimage attacks, collisions, or other various vulnerabilities. The statement I'd like to prove cryptographically/probabalistically is something along the lines of 'it's impossible to find any two sequences of byte buffers that have the same list hash'.

There's a couple applications that I think would be a good use-case:

* Integrate list hashes into some ordered kv-store to demonstrate that it's possible to maintain the hash of the whole database without any additional IO and minimal impact on the data structures.
* Replace [merkle-DAG data structure used by IPFS](https://docs.ipfs.io/how-to/work-with-blocks/) to store large files in chunks.
* Use list hashes to develop a system to syncronize two data stores, like rsync. (Say, a bloom filter of the list hashes of all sub-sequences?)
* Go head-to-head with merkle trees as the defacto primitive for crypto data structures. (Ok maybe it's a little soon for this one 😅.)
"

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
GaloisFields = "8d0d7f98-d412-5cd4-8397-071c807280aa"
LinearAlgebraX = "9b3f67b0-2d00-526e-9884-9e4938f8fb88"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"
SHA = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[compat]
GaloisFields = "~1.1.0"
LinearAlgebraX = "~0.1.7"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "84918055d15b3114ede17ac6a7182f68870c16f7"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.1"

[[ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[ExprTools]]
git-tree-sha1 = "b7e3d17636b348f005f11040025ae8c6f645fe92"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.6"

[[Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[GaloisFields]]
deps = ["HTTP", "JSON", "LinearAlgebra", "Polynomials", "Primes", "Random", "Requires", "Serialization"]
git-tree-sha1 = "4de10dbfd6c51fca43664e5efb03557c88a80841"
uuid = "8d0d7f98-d412-5cd4-8397-071c807280aa"
version = "1.1.0"

[[HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "c6a1fff2fd4b1da29d3dccaffb1e1001244d844e"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.12"

[[IniFile]]
deps = ["Test"]
git-tree-sha1 = "098e4d2c533924c921f9f9847274f2ad89e018b8"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.0"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[Intervals]]
deps = ["Dates", "Printf", "RecipesBase", "Serialization", "TimeZones"]
git-tree-sha1 = "323a38ed1952d30586d0fe03412cde9399d3618b"
uuid = "d8418881-c3e1-53bb-8760-2df7ec849ed5"
version = "1.5.0"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "81690084b6198a2e1da36fcfda16eeca9f9f24e4"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.1"

[[LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[LinearAlgebra]]
deps = ["Libdl"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[LinearAlgebraX]]
deps = ["LinearAlgebra", "Mods", "SimplePolynomials"]
git-tree-sha1 = "0941dac2304b2c8354f75974421f910d880741cd"
uuid = "9b3f67b0-2d00-526e-9884-9e4938f8fb88"
version = "0.1.7"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "Random", "Sockets"]
git-tree-sha1 = "1c38e51c3d08ef2278062ebceade0e46cefc96fe"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.0.3"

[[MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[Mocking]]
deps = ["ExprTools"]
git-tree-sha1 = "916b850daad0d46b8c71f65f719c49957e9513ed"
uuid = "78c3b35d-d492-501b-9361-3d52fe80e533"
version = "0.7.1"

[[Mods]]
git-tree-sha1 = "3747d335e703ecea870ff473c7677d457fd575d0"
uuid = "7475f97c-0381-53b1-977b-4c60186c8d62"
version = "1.3.0"

[[MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[Multisets]]
git-tree-sha1 = "3c478ef38e8d858aed1aeba3a2043be72154e3c7"
uuid = "3b2b4ff1-bcff-5658-a3ee-dbcf1ce5ac09"
version = "0.4.2"

[[NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "2bf78c5fd7fa56d2bbf1efbadd45c1b8789e6f57"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.10.2"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "c8abc88faa3f7a3950832ac5d6e690881590d6dc"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "1.1.0"

[[Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[Polynomials]]
deps = ["Intervals", "LinearAlgebra", "OffsetArrays", "RecipesBase"]
git-tree-sha1 = "0b15f3597b01eb76764dd03c3c23d6679a3c32c8"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "1.2.1"

[[Primes]]
git-tree-sha1 = "afccf037da52fa596223e5a0e331ff752e0e845c"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.0"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[RecipesBase]]
git-tree-sha1 = "b3fb709f3c97bfc6e948be68beeecb55a0b340ae"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.1.1"

[[Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "4036a3bd08ac7e968e27c203d45f5fff15020621"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.1.3"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[SimplePolynomials]]
deps = ["Mods", "Multisets", "Polynomials", "Primes"]
git-tree-sha1 = "12e65afc86adb18265de033311b442b2575747ad"
uuid = "cc47b68c-3164-5771-a705-2bc0097375a0"
version = "0.2.7"

[[Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[TimeZones]]
deps = ["Dates", "Future", "LazyArtifacts", "Mocking", "Pkg", "Printf", "RecipesBase", "Serialization", "Unicode"]
git-tree-sha1 = "81753f400872e5074768c9a77d4c44e70d409ef0"
uuid = "f269a46b-ccf7-5d73-abea-4c690281aa53"
version = "1.5.6"

[[URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╟─8ad81282-7a8e-4769-81d4-bb7fd3b7b1a5
# ╟─d8ab88dc-2278-42a7-b6e0-2dd67aad9c0f
# ╟─d573e0aa-4110-4d83-9ef3-68addd3ff9b8
# ╠═58f1747c-e2a1-11eb-2c26-95ef11565850
# ╠═c7472820-0a79-4156-8d0e-e4367868c57e
# ╟─0ad5b059-9539-437b-97f5-e73ec1c06419
# ╠═b5fc120c-7727-4b53-bdc5-5352cd96dac1
# ╟─ef00d3c8-6ab4-4f86-ae6b-9fca3d71d5cb
# ╟─e7621cea-c449-4b43-bd50-ba4646e2b48a
# ╠═9b964dbc-c78b-4246-bf7a-0761b9c97bc7
# ╠═52fcf820-cddd-4e36-86e9-afcb32f09403
# ╟─62e4fa8b-7fb7-4904-963f-a1b9c4176dc7
# ╠═710eef31-f07d-41f0-b410-dc6727f3589e
# ╟─d9855f4f-c6d9-4bcf-ac25-e7f43b8b443b
# ╠═cec1b298-bf51-4500-97ae-5189a2049286
# ╟─89c9ca94-a428-442e-8dee-3a7b1b4f7338
# ╟─f3e8221a-8d34-4442-a58f-9708952cb7b3
# ╠═5bb70b21-5cb4-4451-8229-3d806e611408
# ╟─fb8b6da2-c515-4afa-b31c-8bd3259cb0f7
# ╠═95a21d02-68d2-477f-b98d-1b83d0da1cdb
# ╟─0c0f5ee4-a3e0-49df-aafd-903a4eaca2d0
# ╠═523979df-d971-48cf-a868-174a4ebc644a
# ╟─a76ee3eb-a025-46f9-b4d8-1a9f47ab9c18
# ╠═c80af06b-3564-4e49-ae18-fcaf70cbc811
# ╟─fec99716-09a8-4cd5-80de-b57e7eb90540
# ╠═40c3a591-1b86-4abe-b7e8-5467f7f173a1
# ╟─3312903b-af7e-4088-bab9-f0e4c42c25fc
# ╟─1d7a824a-73c5-4e21-9b2b-6b4923ea49c7
# ╠═786df805-b24b-482d-a842-1756a1c38510
# ╟─533b55c8-cbfc-4ab8-bb3c-3813e5b58bea
# ╠═38f2b3fc-5078-4ee2-86f2-19d99c254166
# ╟─7b2c1def-cdc7-4102-a66d-5a914a1b1eb0
# ╠═8b251bfd-f61a-464f-a1b6-458dde096e14
# ╟─c58360f7-010d-459e-a9ab-20469f0a98fc
# ╠═3f6fefd6-82c7-4f1f-b4c2-406faf39d0d2
# ╟─7feb5e3c-1a71-432e-896a-263c8fff1815
# ╠═80028718-f19c-492a-9a0e-09c5819dd633
# ╟─35bb2711-b8a9-4561-9e38-dc070b35d5dc
# ╟─a82c8cc1-e69c-4c19-bf40-e6efeadb00de
# ╠═40a13273-b40f-4287-8b73-e1a8bcf119e2
# ╟─57f46d87-ca7a-4618-90f3-5a9ea8655bb2
# ╟─201e4088-4481-4bdc-8c79-4a4fa5f23d61
# ╟─d9c3ebda-3678-4421-980e-0efd24f8a8f9
# ╠═87bedc95-8423-4be9-b25c-08e775ffb540
# ╟─a8a7559f-2f9f-4f30-b8b8-f964b59da354
# ╠═85d19b91-c0a2-43be-abbb-28a7387195a6
# ╠═258da454-caf5-46c9-a45a-4df2fd44bf02
# ╠═f4ca59b9-ca2f-4458-950a-b87dee7e3115
# ╟─3d7c33df-74fb-49a8-98e2-c32456c1bf76
# ╠═f8bc7a62-3b4a-47f7-936f-f5370a6efe97
# ╟─67ee30f8-240d-4cfd-b03f-ca620790e0b0
# ╟─786a6537-296e-40d6-830f-93ac232e2cf7
# ╟─ceb7e3f7-031f-4b9e-9fa8-7678924db12c
# ╠═8e3a5919-d825-48bf-bde7-0a1cffca46b6
# ╟─1167b981-6baa-40f6-934a-928eea39ee1a
# ╠═8dca48af-89c2-4528-a976-fd3c9992e2c2
# ╟─57adf406-ae46-41d9-a9ad-0c5c66489c12
# ╠═be3d561c-75ab-4e17-91b6-a23921ad3aa4
# ╟─a9283164-7b8d-40e0-9a3a-ff26fcc477d0
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
