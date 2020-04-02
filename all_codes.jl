#HMMBase

"""
Hidden Markov Models for Julia.

[Documentation](https://maxmouchet.github.io/HMMBase.jl/stable/).
[Issues](https://github.com/maxmouchet/HMMBase.jl/issues).
"""
module HMMBase

using ArgCheck
using Clustering
using Distributions
using Hungarian
using LinearAlgebra

import Base: ==, copy, rand, size, OneTo
import Distributions: fit_mle, loglikelihood
import Random: AbstractRNG, GLOBAL_RNG

export
    # hmm.jl
    AbstractHMM,
    HMM,
    copy,
    rand,
    size,
    nparams,
    permute,
    statdists,
    istransmat,
    likelihoods,
    # messages_api.jl
    forward,
    backward,
    posteriors,
    loglikelihood,
    # mle_api.jl
    fit_mle,
    # viterbi_api.jl
    viterbi,
    # utilities.jl,
    gettransmat,
    randtransmat,
    remapseq

include("hmm.jl")
include("mle.jl")
include("mle_api.jl")
include("mle_init.jl")
include("messages.jl")
include("messages_api.jl")
include("viterbi.jl")
include("viterbi_api.jl")
include("likelihoods.jl")
include("likelihoods_api.jl")
include("utilities.jl")

# To be removed in a future version
# ---------------------------------

export
    n_parameters,
    log_likelihoods,
    forward_backward,
    messages_backwards,
    messages_backwards_log,
    messages_forwards,
    messages_forwards_log,
    compute_transition_matrix,
    rand_transition_matrix

@deprecate n_parameters(hmm) nparams(hmm)
@deprecate log_likelihoods(hmm, observations) likelihoods(hmm, observations, logl = true)

@deprecate forward_backward(init_distn, trans_matrix, log_likelihoods) posteriors(init_distn, trans_matrix, log_likelihoods, logl = true)
@deprecate messages_forwards(init_distn, trans_matrix, log_likelihoods) forward(init_distn, trans_matrix, log_likelihoods, logl = true)
@deprecate messages_backwards(init_distn, trans_matrix, log_likelihoods) backward(init_distn, trans_matrix, log_likelihoods, logl = true)

@deprecate forward_backward(hmm, observations) posteriors(hmm, observations, logl = true)
@deprecate messages_forwards(hmm, observations) forward(hmm, observations, logl = true)
@deprecate messages_backwards(hmm, observations) backward(hmm, observations, logl = true)

@deprecate messages_forwards_log(init_distn, trans_matrix, log_likelihoods) log.(forward(init_distn, trans_matrix, log_likelihoods, logl = true)[1])
@deprecate messages_backwards_log(trans_matrix, log_likelihoods) log.(backward(init_distn, trans_matrix, log_likelihoods, logl = true)[1])

@deprecate compute_transition_matrix(seq) gettransmat(seq, relabel = true)
@deprecate rand_transition_matrix(K, α = 1.0) randtransmat(K, α)

end


#hmm.jl
"""
    AbstractHMM{F<:VariateForm}
A custom HMM type must at-least implement the following interface:
```julia
struct CustomHMM{F,T} <: AbstractHMM{F}
    a::AbstractVector{T}               # Initial state distribution
    A::AbstractMatrix{T}               # Transition matrix
    B::AbstractVector{Distribution{F}} # Observations distributions
    # Optional, custom, fields ....
end
```
"""
abstract type AbstractHMM{F<:VariateForm} end

"""
    HMM([a, ]A, B) -> HMM
Build an HMM with transition matrix `A` and observation distributions `B`.
If the initial state distribution `a` is not specified, a uniform distribution is assumed.
Observations distributions can be of different types (for example `Normal` and `Exponential`),
but they must be of the same dimension.
**Arguments**
- `a::AbstractVector{T}`: initial probabilities vector.
- `A::AbstractMatrix{T}`: transition matrix.
- `B::AbstractVector{<:Distribution{F}}`: observations distributions.
**Example**
```julia
using Distributions, HMMBase
hmm = HMM([0.9 0.1; 0.1 0.9], [Normal(0,1), Normal(10,1)])
```
"""
struct HMM{F,T} <: AbstractHMM{F}
    a::Vector{T}
    A::Matrix{T}
    B::Vector{Distribution{F}}
    HMM{F,T}(a, A, B) where {F,T} = assert_hmm(a, A, B) && new(a, A, B)
end

HMM(a::AbstractVector{T}, A::AbstractMatrix{T}, B::AbstractVector{<:Distribution{F}}) where {F,T} = HMM{F,T}(a, A, B)
HMM(A::AbstractMatrix{T}, B::AbstractVector{<:Distribution{F}}) where {F,T} = HMM{F,T}(ones(size(A)[1])/size(A)[1], A, B)

"""
    assert_hmm(a, A, B)
Throw an `ArgumentError` if the initial state distribution and the transition matrix rows does not sum to 1,
and if the observation distributions do not have the same dimensions.
"""
function assert_hmm(a::AbstractVector,
                    A::AbstractMatrix,
                    B::AbstractVector{<:Distribution})
    @argcheck isprobvec(a)
    @argcheck istransmat(A)
    @argcheck all(length.(B) .== length(B[1])) ArgumentError("All distributions must have the same dimensions")
    @argcheck length(a) == size(A,1) == length(B)
    return true
end

"""
    issquare(A) -> Bool
Return true if `A` is a square matrix.
"""
issquare(A::AbstractMatrix) = size(A,1) == size(A,2)

"""
    istransmat(A) -> Bool
Return true if `A` is square and its rows sums to 1.
"""
istransmat(A::AbstractMatrix) = issquare(A) && all([isprobvec(A[i,:]) for i in 1:size(A,1)])

==(h1::AbstractHMM, h2::AbstractHMM) = (h1.a == h2.a) && (h1.A == h2.A) && (h1.B == h2.B)

"""
    rand([rng, ]hmm, T; init, seq) -> Array | (Vector, Array)
Sample a trajectory of `T` timesteps from `hmm`.
**Keyword Arguments**
- `init::Integer = rand(Categorical(hmm.a))`: initial state.
- `seq::Bool = false`: whether to return the hidden state sequence or not.
**Output**
- `Vector{Int}` (if `seq == true`): hidden state sequence.
- `Vector{Float64}` (for `Univariate` HMMs): observations (`T`).
- `Matrix{Float64}` (for `Multivariate` HMMs): observations (`T x dim(obs)`).
**Examples**
```julia
using Distributions, HMMBase
hmm = HMM([0.9 0.1; 0.1 0.9], [Normal(0,1), Normal(10,1)])
y = rand(hmm, 1000) # or
z, y = rand(hmm, 1000, seq = true)
size(y) # (1000,)
```
```julia
using Distributions, HMMBase
hmm = HMM([0.9 0.1; 0.1 0.9], [MvNormal(ones(2)), MvNormal(ones(2))])
y = rand(hmm, 1000) # or
z, y = rand(hmm, 1000, seq = true)
size(y) # (1000, 2)
```
"""
function rand(rng::AbstractRNG, hmm::AbstractHMM, T::Integer; init = rand(rng, Categorical(hmm.a)), seq = false)
    z = Vector{Int}(undef, T)
    (T >= 1) && (z[1] = init)
    for t in 2:T
        z[t] = rand(rng, Categorical(hmm.A[z[t-1],:]))
    end
    y = rand(rng, hmm, z)
    seq ? (z, y) : y
end

"""
    rand([rng, ]hmm, z) -> Array
Sample observations from `hmm` according to trajectory `z`.
**Output**
- `Vector{Float64}` (for `Univariate` HMMs): observations (`T`).
- `Matrix{Float64}` (for `Multivariate` HMMs): observations (`T x dim(obs)`).
**Example**
```julia
using Distributions, HMMBase
hmm = HMM([0.9 0.1; 0.1 0.9], [Normal(0,1), Normal(10,1)])
y = rand(hmm, [1, 1, 2, 2, 1])
```
"""
function rand(rng::AbstractRNG, hmm::AbstractHMM{Univariate}, z::AbstractVector{<:Integer})
    y = Vector{Float64}(undef, length(z))
    for t in eachindex(z)
        y[t] = rand(rng, hmm.B[z[t]])
    end
    y
end

function rand(rng::AbstractRNG, hmm::AbstractHMM{Multivariate}, z::AbstractVector{<:Integer})
    y = Matrix{Float64}(undef, length(z), size(hmm, 2))
    for t in eachindex(z)
        y[t,:] = rand(rng, hmm.B[z[t]])
    end
    y
end

rand(hmm::AbstractHMM, T::Integer; kwargs...) = rand(GLOBAL_RNG, hmm, T; kwargs...)

rand(hmm::AbstractHMM, z::AbstractVector{<:Integer}) = rand(GLOBAL_RNG, hmm, z)

"""
    size(hmm, [dim]) -> Int | Tuple
Return the number of states in `hmm` and the dimension of the observations.
**Example**
```jldoctest
using Distributions, HMMBase
hmm = HMM([0.9 0.1; 0.1 0.9], [Normal(0,1), Normal(10,1)])
size(hmm)
# output
(2, 1)
```
"""
size(hmm::AbstractHMM, dim=:) = (length(hmm.B), length(hmm.B[1]))[dim]

"""
    copy(hmm) -> HMM
Return a copy of `hmm`.
"""
copy(hmm::HMM) = HMM(copy(hmm.a), copy(hmm.A), copy(hmm.B))

"""
    permute(hmm, perm) -> HMM
Permute the states of `hmm` according to `perm`.
**Arguments**
- `perm::Vector{<:Integer}`: permutation of the states.
**Example**
```julia
using Distributions, HMMBase
hmm = HMM([0.8 0.2; 0.1 0.9], [Normal(0,1), Normal(10,1)])
hmm = permute(hmm, [2, 1])
hmm.A # [0.9 0.1; 0.2 0.8]
hmm.B # [Normal(10,1), Normal(0,1)]
```
"""
function permute(hmm::AbstractHMM, perm::Vector{<:Integer})
    a = hmm.a[perm]
    B = hmm.B[perm]
    A = zeros(size(hmm.A))
    for i in 1:size(A,1), j in 1:size(A,2)
        A[i,j] = hmm.A[perm[i],perm[j]]
    end
    HMM(a, A, B)
end

"""
    nparams(hmm) -> Int
Return the number of _free_ parameters in `hmm`, without counting the
observation distributions parameters.
**Example**
```jldoctest
using Distributions, HMMBase
hmm = HMM([0.9 0.1; 0.1 0.9], [Normal(0,1), Normal(10,1)])
nparams(hmm)
# output
3
```
"""
function nparams(hmm::AbstractHMM)
    (length(hmm.a) - 1) + (length(hmm.A) - size(hmm.A, 1))
end

"""
    statdists(hmm) -> Vector{Vector}
Return the stationnary distribution(s) of `hmm`.
That is, the eigenvectors of transpose(hmm.A) with eigenvalues 1.
"""
function statdists(hmm::AbstractHMM)
    eig = eigen(collect(transpose(hmm.A)))
    dists = []
    for (i, val) in enumerate(eig.values)
        if val == 1.0
            dist = eig.vectors[:,i]
            dist /= sum(dist)
            push!(dists, dist)
        end
    end
    dists
end


## likelihoods
function likelihoods!(L::AbstractMatrix, hmm::AbstractHMM{Univariate}, observations)
    T, K = size(observations, 1), size(hmm, 1)
    @argcheck size(L) == (T, K)
    @inbounds for i in OneTo(K), t in OneTo(T)
        L[t,i] = pdf(hmm.B[i], observations[t])
    end
end

function loglikelihoods!(LL::AbstractMatrix, hmm::AbstractHMM{Univariate}, observations)
    T, K = size(observations, 1), size(hmm, 1)
    @argcheck size(LL) == (T, K)
    @inbounds for i in OneTo(K), t in OneTo(T)
        LL[t,i] = logpdf(hmm.B[i], observations[t])
    end
end

function likelihoods!(L::AbstractMatrix, hmm::AbstractHMM{Multivariate}, observations)
    T, K = size(observations, 1), size(hmm, 1)
    @argcheck size(L) == (T, K)
    @inbounds for i in OneTo(K), t in OneTo(T)
        L[t,i] = pdf(hmm.B[i], view(observations, t, :))
    end
end

function loglikelihoods!(LL::AbstractMatrix, hmm::AbstractHMM{Multivariate}, observations)
    T, K = size(observations, 1), size(hmm, 1)
    @argcheck size(LL) == (T, K)
    @inbounds for i in OneTo(K), t in OneTo(T)
        LL[t,i] = logpdf(hmm.B[i], view(observations, t, :))
    end
end



#likelihoods_api
"""
    likelihoods(hmm, observations; logl) -> Matrix
Return the likelihood per-state and per-observation.
**Output**
- `Matrix{Float64}`: likelihoods matrix (`T x K`).
**Example**
```julia
using Distributions, HMMBase
hmm = HMM([0.9 0.1; 0.1 0.9], [Normal(0,1), Normal(10,1)])
y = rand(hmm, 1000)
L = likelihoods(hmm, y)
LL = likelihoods(hmm, y, logl = true)
```
"""
function likelihoods(hmm::AbstractHMM, observations; logl = false, robust = false)
    T, K = size(observations, 1), size(hmm, 1)
    L = Matrix{Float64}(undef, T, K)

    if logl
        loglikelihoods!(L, hmm, observations)
        robust && replace!(L, -Inf => nextfloat(-Inf), Inf => log(prevfloat(Inf)))
    else
        likelihoods!(L, hmm, observations)
        robust && replace!(L, -Inf => nextfloat(-Inf), Inf => prevfloat(Inf))
    end

    L
end


#message_jl
# Original implementations by @nantonel
# https://github.com/maxmouchet/HMMBase.jl/pull/6

# *log! methods use the samples log-likelihood instead of the likelihood.

# In-place forward pass, where α and c are allocated beforehand.
function forward!(α::AbstractMatrix, c::AbstractVector, a::AbstractVector, A::AbstractMatrix, L::AbstractMatrix)
    @argcheck size(α, 1) == size(L, 1) == size(c, 1)
    @argcheck size(α, 2) == size(L, 2) == size(a, 1) == size(A, 1) == size(A, 2)

    T, K = size(L)
    (T == 0) && return

    fill!(α, 0.0)
    fill!(c, 0.0)

    for j in OneTo(K)
        α[1,j] = a[j] * L[1,j]
        c[1] += α[1,j]
    end

    for j in OneTo(K)
        α[1,j] /= c[1]
    end

    @inbounds for t in 2:T
        for j in OneTo(K)
            for i in OneTo(K)
                α[t,j] += α[t-1,i] * A[i,j]
            end
            α[t,j] *= L[t,j]
            c[t] += α[t,j]
        end

        for j in OneTo(K)
            α[t,j] /= c[t]
        end
    end
end

# In-place backward pass, where β and c are allocated beforehand.
function backward!(β::AbstractMatrix, c::AbstractVector, a::AbstractVector, A::AbstractMatrix, L::AbstractMatrix)
    @argcheck size(β, 1) == size(L, 1) == size(c, 1)
    @argcheck size(β, 2) == size(L, 2) == size(a, 1) == size(A, 1) == size(A, 2)

    T, K = size(L)
    (T == 0) && return

    fill!(β, 0.0)
    fill!(c, 0.0)

    for j in OneTo(K)
        β[end,j] = 1.0
    end

    @inbounds for t in T-1:-1:1
        for j in OneTo(K)
            for i in OneTo(K)
                β[t,j] += β[t+1,i] * A[j,i] * L[t+1,i]
            end
            c[t+1] += β[t,j]
        end

        for j in OneTo(K)
            β[t,j] /= c[t+1]
        end
    end

    for j in OneTo(K)
        c[1] += a[j] * L[1,j] * β[1,j]
    end
end

# In-place forward pass, where α and c are allocated beforehand.
function forwardlog!(α::AbstractMatrix, c::AbstractVector, a::AbstractVector, A::AbstractMatrix, LL::AbstractMatrix)
    @argcheck size(α, 1) == size(LL, 1) == size(c, 1)
    @argcheck size(α, 2) == size(LL, 2) == size(a, 1) == size(A, 1) == size(A, 2)

    T, K = size(LL)
    (T == 0) && return

    fill!(α, 0.0)
    fill!(c, 0.0)

    m = vec_maximum(view(LL, 1, :))

    for j in OneTo(K)
        α[1,j] = a[j] * exp(LL[1,j] - m)
        c[1] += α[1,j]
    end

    for j in OneTo(K)
        α[1,j] /= c[1]
    end

    c[1] *= exp(m) + eps()

    @inbounds for t in 2:T
        m = vec_maximum(view(LL, t, :))

        for j in OneTo(K)
            for i in OneTo(K)
                α[t,j] += α[t-1,i] * A[i,j]
            end
            α[t,j] *= exp(LL[t,j] - m)
            c[t] += α[t,j]
        end

        for j in OneTo(K)
            α[t,j] /= c[t]
        end

        c[t] *= exp(m) + eps()
    end
end

# In-place backward pass, where β and c are allocated beforehand.
function backwardlog!(β::AbstractMatrix, c::AbstractVector, a::AbstractVector, A::AbstractMatrix, LL::AbstractMatrix)
    @argcheck size(β, 1) == size(LL, 1) == size(c, 1)
    @argcheck size(β, 2) == size(LL, 2) == size(a, 1) == size(A, 1) == size(A, 2)

    T, K = size(LL)
    L = zeros(K)
    (T == 0) && return

    fill!(β, 0.0)
    fill!(c, 0.0)

    for j in OneTo(K)
        β[end,j] = 1.0
    end

    @inbounds for t in T-1:-1:1
        m = vec_maximum(view(LL, t+1, :))

        for i in OneTo(K)
            L[i] = exp(LL[t+1,i] - m)
        end

        for j in OneTo(K)
            for i in OneTo(K)
                β[t,j] += β[t+1,i] * A[j,i] * L[i]
            end
            c[t+1] += β[t,j]
        end

        for j in OneTo(K)
            β[t,j] /= c[t+1]
        end

        c[t+1] *= exp(m) + eps()
    end

    m = vec_maximum(view(LL, 1,:))

    for j in OneTo(K)
        c[1] += a[j] * exp(LL[1,j] - m) * β[1,j]
    end

    c[1] *= exp(m) + eps();
end

# In-place posterior computation, where γ is allocated beforehand.
function posteriors!(γ::AbstractMatrix, α::AbstractMatrix, β::AbstractMatrix)
    @argcheck size(γ) == size(α) == size(β)
    T, K = size(α)
    for t in OneTo(T)
        c = 0.0
        for i = OneTo(K)
            γ[t,i] = α[t,i] * β[t,i]
            c += γ[t,i]
        end

        for i in OneTo(K)
            γ[t,i] /= c
        end
    end
end


#messages_api
# Forward/Backward

# {forward,backward}(a, A, L)
for f in (:forward, :backward)
    f!  = Symbol("$(f)!")    # forward!
    fl! = Symbol("$(f)log!") # forwardlog!

    @eval begin
        """
            $($f)(a, A, L; logl) -> (Vector, Float)
        Compute $($f) probabilities using samples likelihoods.
        See [Forward-backward algorithm](https://en.wikipedia.org/wiki/Forward–backward_algorithm).

        **Output**
        - `Vector{Float64}`: $($f) probabilities.
        - `Float64`: log-likelihood of the observed sequence.
        """
        function $(f)(a::AbstractVector, A::AbstractMatrix, L::AbstractMatrix; logl = false)
            m = Matrix{Float64}(undef, size(L))
            c = Vector{Float64}(undef, size(L)[1])
            if logl
                $(fl!)(m, c, a, A, L)
            else
                warn_logl(L)
                $(f!)(m, c, a, A, L)
            end
            m, sum(log.(c))
        end
    end
end

# {forward,backward}(hmm, observations)
for f in (:forward, :backward)
    f!  = Symbol("$(f)!")    # forward!
    fl! = Symbol("$(f)log!") # forwardlog!

    @eval begin
        """
            $($f)(hmm, observations; logl, robust) -> (Vector, Float)
        Compute $($f) probabilities of the `observations` given the `hmm` model.
        **Output**
        - `Vector{Float64}`: $($f) probabilities.
        - `Float64`: log-likelihood of the observed sequence.
        **Example**
        ```julia
        using Distributions, HMMBase
        hmm = HMM([0.9 0.1; 0.1 0.9], [Normal(0,1), Normal(10,1)])
        y = rand(hmm, 1000)
        probs, tot = $($f)(hmm, y)
        ```
        """
        function $(f)(hmm::AbstractHMM, observations; robust = false, kwargs...)
            L = likelihoods(hmm, observations; robust = robust, kwargs...)
            $(f)(hmm.a, hmm.A, L; kwargs...)
        end
    end
end


# Posteriors

"""
    posteriors(α, β) -> Vector
Compute posterior probabilities from `α` and `β`.
**Arguments**
- `α::AbstractVector`: forward probabilities.
- `β::AbstractVector`: backward probabilities.
"""
function posteriors(α::AbstractMatrix, β::AbstractMatrix)
    γ = Matrix{Float64}(undef, size(α))
    posteriors!(γ, α, β)
    γ
end

"""
    posteriors(a, A, L; logl) -> Vector
Compute posterior probabilities using samples likelihoods.
"""
function posteriors(a::AbstractVector, A::AbstractMatrix, L::AbstractMatrix; kwargs...)
    α, _ = forward(a, A, L; kwargs...)
    β, _ = backward(a, A, L; kwargs...)
    posteriors(α, β)
end

"""
    posteriors(hmm, observations; logl, robust) -> Vector
Compute posterior probabilities using samples likelihoods.
**Example**
```julia
using Distributions, HMMBase
hmm = HMM([0.9 0.1; 0.1 0.9], [Normal(0,1), Normal(10,1)])
y = rand(hmm, 1000)
γ = posteriors(hmm, y)
```
"""
function posteriors(hmm::AbstractHMM, observations; robust = false, kwargs...)
    L = likelihoods(hmm, observations; robust = robust, kwargs...)
    posteriors(hmm.a, hmm.A, L; kwargs...)
end

# Likelihood

"""
    loglikelihood(hmm, observations; logl, robust) -> Float64
Compute the log-likelihood of the observations under the model.
This is defined as the sum of the log of the normalization coefficients in the forward filter.
**Output**
- `Float64`: log-likelihood of the observations sequence under the model.
**Example**
```jldoctest
using Distributions, HMMBase
hmm = HMM([0.9 0.1; 0.1 0.9], [Normal(0,1), Normal(10,1)])
loglikelihood(hmm, [0.15, 0.10, 1.35])
# output
-4.588183811489616
```
"""
function loglikelihood(hmm::AbstractHMM, observations; logl = false, robust = false)
    forward(hmm, observations, logl = logl, robust = robust)[2]
end


#MLE
# In-place update of the initial state distribution.
function update_a!(a::AbstractVector, α::AbstractMatrix, β::AbstractMatrix)
    @argcheck size(α, 1) == size(β, 1)
    @argcheck size(α, 2) == size(β, 2) == size(a, 1)

    K = length(a)
    c = 0.0

    for i in OneTo(K)
        a[i] = α[1,i] * β[1,i]
        c += a[i]
    end

    for i in OneTo(K)
        a[i] /= c
    end
end

# In-place update of the transition matrix.
function update_A!(A::AbstractMatrix, ξ::AbstractArray, α::AbstractMatrix, β::AbstractMatrix, LL::AbstractMatrix)
    @argcheck size(α, 1) == size(β, 1) == size(LL, 1) == size(ξ, 1)
    @argcheck size(α, 2) == size(β, 2) == size(LL, 2) == size(A, 1) == size(A, 2) == size(ξ, 2) == size(ξ, 3)

    T, K = size(LL)

    @inbounds for t in OneTo(T - 1)
        m = vec_maximum(view(LL, t+1, :))
        c = 0.0

        for i in OneTo(K), j in OneTo(K)
            ξ[t,i,j] = α[t,i] * A[i,j] * exp(LL[t + 1,j] - m) * β[t + 1,j]
            c += ξ[t,i,j]
        end

        for i in OneTo(K), j in OneTo(K)
            ξ[t,i,j] /= c
        end
    end

    fill!(A, 0.0)

    @inbounds for i in OneTo(K)
        c = 0.0

        for j in OneTo(K)
            for t in OneTo(T - 1)
                A[i,j] += ξ[t,i,j]
            end
            c += A[i,j]
        end

        for j in OneTo(K)
            A[i,j] /= c
        end
    end
end

# In-place update of the observations distributions.
function update_B!(B::AbstractVector, γ::AbstractMatrix, observations, estimator)
    @argcheck size(γ, 1) == size(observations, 1)
    @argcheck size(γ, 2) == size(B, 1)
    K = length(B)
    for i in OneTo(K)
        if sum(γ[:,i]) > 0
            B[i] = estimator(typeof(B[i]), permutedims(observations), γ[:,i])
        end
    end
end

function fit_mle!(hmm::AbstractHMM, observations; display = :none, maxiter = 100, tol=1e-3, robust = false, estimator=fit_mle)
    @argcheck display in [:none, :iter, :final]
    @argcheck maxiter >= 0

    T, K = size(observations, 1), size(hmm, 1)
    history = EMHistory(false, 0, [])

    # Allocate memory for in-place updates
    c = zeros(T)
    α = zeros(T, K)
    β = zeros(T, K)
    γ = zeros(T, K)
    ξ = zeros(T, K, K)
    LL = zeros(T, K)

    loglikelihoods!(LL, hmm, observations)
    robust && replace!(LL, -Inf => nextfloat(-Inf), Inf => log(prevfloat(Inf)))

    forwardlog!(α, c, hmm.a, hmm.A, LL)
    backwardlog!(β, c, hmm.a, hmm.A, LL)
    posteriors!(γ, α, β)

    logtot = sum(log.(c))
    (display == :iter) && println("Iteration 0: logtot = $logtot")

    for it in 1:maxiter
        update_a!(hmm.a, α, β)
        update_A!(hmm.A, ξ, α, β, LL)
        update_B!(hmm.B, γ, observations, estimator)

        # Ensure the "connected-ness" of the states,
        # this prevents case where there is no transitions
        # between two extremely likely observations.
        robust && (hmm.A .+= eps())

        @check isprobvec(hmm.a)
        @check istransmat(hmm.A)

        loglikelihoods!(LL, hmm, observations)
        robust && replace!(LL, -Inf => nextfloat(-Inf), Inf => log(prevfloat(Inf)))

        forwardlog!(α, c, hmm.a, hmm.A, LL)
        backwardlog!(β, c, hmm.a, hmm.A, LL)
        posteriors!(γ, α, β)

        logtotp = sum(log.(c))
        (display == :iter) && println("Iteration $it: logtot = $logtotp")

        push!(history.logtots, logtotp)
        history.iterations += 1

        if abs(logtotp - logtot) < tol
            (display in [:iter, :final]) && println("EM converged in $it iterations, logtot = $logtotp")
            history.converged = true
            break
        end

        logtot = logtotp
    end

    if !history.converged
        if display in [:iter, :final]
            println("EM has not converged after $(history.iterations) iterations, logtot = $logtot")
        end
    end

    history
end

mutable struct EMHistory
    converged::Bool
    iterations::Int
    logtots::Vector{Float64}
end


#mle_api
"""
    fit_mle(hmm, observations; ...) -> AbstractHMM
Estimate the HMM parameters using the EM (Baum-Welch) algorithm, with `hmm` as the initial state.
**Keyword Arguments**
- `display::Symbol = :none`: when to display convergence logs, can be set to `:iter` or `:final`.
- `init::Symbol = :none`: if set to `:kmeans` the HMM parameters will be initialized using a K-means clustering.
- `maxiter::Integer = 100`: maximum number of iterations to perform.
- `tol::Integer = 1e-3`: stop the algorithm when the improvement in the log-likelihood is less than `tol`.
**Output**
- `<:AbstractHMM`: a copy of the original HMM with the updated parameters.
"""
function fit_mle(hmm::AbstractHMM, observations; init = :none, kwargs...)
    hmm = copy(hmm)

    if init == :kmeans
        kmeans_init!(hmm, observations, display = get(kwargs, :display, :none))
    end

    history = fit_mle!(hmm, observations; kwargs...)
    hmm, history
end


#mle_init.jl
# Functions for initializing HMM parameters from the observations

function kmeans_init!(hmm::AbstractHMM, observations; kwargs...)
    K = size(hmm, 1)

    res = kmeans(permutedims(observations), size(hmm, 1); kwargs...)
    seq = res.assignments

    # Initialize A
    copyto!(hmm.A, gettransmat(seq)[2])
    @check istransmat(hmm.A)

    # Initialize B
    for i in OneTo(K)
        observations_ = view(observations, seq .== i, :)
        if length(observations_) > 0
            hmm.B[i] = fit_mle(typeof(hmm.B[i]), permutedims(observations_))
        end
    end
end




#utilities
"""
    gettransmat(seq; relabel = false) -> (Dict, Matrix)
Return the transition matrix associated to the label sequence `seq`.
The labels must be positive integer.
**Arguments**
- `seq::Vector{<:Integer}`: positive label sequence.
**Keyword Arguments**
- `relabel::Bool = false`: if set to true the sequence
  will be made contiguous. E.g. `[7,7,9,9,1,1]` will become `[2,2,3,3,1,1]`.
**Output**
- `Dict{Integer,Integer}`: the mapping between the original and the new labels.
- `Matrix{Float64}`: the transition matrix.
"""
function gettransmat(seq::Vector{<:Integer}; relabel = false)
    @argcheck all(seq .>= 0)

    if relabel
        # /!\ Sort is important here, so that we don't relabel already contiguous states.
        mapping = Dict([(x[2], x[1]) for x in enumerate(sort(unique(seq)))])
    else
        mapping = Dict([(x, x) for x in unique(seq)])
    end

    (length(mapping) == 0) && return mapping, Float64[]
    K = maximum(values(mapping))

    transmat = zeros(K, K)
    for i in 1:length(seq)-1
        transmat[mapping[seq[i]], mapping[seq[i+1]]] += 1
    end
    transmat = transmat ./ sum(transmat, dims=2)
    transmat[isnan.(transmat)] .= 0.0

    mapping, transmat
end

"""
    randtransmat([rng,] prior) -> Matrix{Float64}
Generate a transition matrix where each row is sampled from `prior`.
The prior must be a multivariate probability distribution, such as a
Dirichlet distribution.
**Arguments**
- `prior::MultivariateDistribution`: distribution over the transition matrix rows.
**Example**
```julia
A = randtransmat(Dirichlet([0.1, 0.1, 0.1]))
```
"""
function randtransmat(rng::AbstractRNG, prior::MultivariateDistribution)
    K = length(prior)
    A = Matrix{Float64}(undef, K, K)
    for i in OneTo(K)
        A[i,:] = rand(rng, prior)
    end
    @check istransmat(A); A
end

randtransmat(prior::MultivariateDistribution) = randtransmat(GLOBAL_RNG, prior)

"""
    randtransmat([rng, ]K, α = 1.0) -> Matrix{Float64}
Generate a transition matrix where each row is sampled from
a Dirichlet distribution of dimension `K` and concentration
parameter `α`.
**Arguments**
- `K::Integer`: number of states.
- `α::Float64 = 1.0`: concentration parameter of the Dirichlet distribution.
**Example**
```julia
A = randtransmat(4)
```
"""
randtransmat(rng::AbstractRNG, K::Integer, α = 1.0) = randtransmat(rng, Dirichlet(K, α))

randtransmat(K::Integer, args...) = randtransmat(GLOBAL_RNG, K, args...)

"""
    remapseq(seq, ref) -> Vector{Integer}
Find the permutations of `seq` indices that maximize the overlap with `ref`.
**Arguments**
- `seq::Vector{Integer}`: sequence to be remapped.
- `ref::Vector{Integer}`: reference sequence.
**Example**
```julia
ref = [1,1,2,2,3,3]
seq = [2,2,3,3,1,1]
remapseq(seq, ref)
# [1,1,2,2,3,3]
```
"""
function remapseq(seq::Vector{<:Integer}, ref::Vector{<:Integer})
    seqlabels, reflabels = unique(seq), unique(ref)
    @argcheck all(seqlabels .> 0) && all(reflabels .> 0)

    # C[i,j]: cost of assigning seq. label `i` to ref. label `j`
    C = zeros(maximum(seqlabels), maximum(reflabels))
    for i in seqlabels, j in reflabels
        C[i,j] = - sum((seq .== i) .& (ref .== j))
    end

    # TODO: Own implementation of the hungarian alg.,
    # to avoid pulling another dependency ?
    assignment, _ = hungarian(C)
    [assignment[x] for x in seq]
end


function warn_logl(L::AbstractMatrix)
    if any(L .< 0)
        @warn "Negative likelihoods values, use the `logl = true` option if you are using log-likelihoods."
    end
end

# ~2x times faster than Base.maximum
# v = rand(25)
# @btime maximum(v)
# @btime vec_maximum(v)
#   63.909 ns (1 allocation: 16 bytes)
#   30.307 ns (1 allocation: 16 bytes)
function vec_maximum(v::AbstractVector)
    m = v[1]
    @inbounds for i = OneTo(length(v))
        if v[i] > m
            m = v[i]
        end
    end
    m
end




##viterbi
# Original implementations by @nantonel
# https://github.com/maxmouchet/HMMBase.jl/pull/6

# *log! methods use the samples log-likelihood instead of the likelihood.

function viterbi!(T1::AbstractMatrix, T2::AbstractMatrix, z::AbstractVector, a::AbstractVector, A::AbstractMatrix, L::AbstractMatrix)
    @argcheck size(T1, 1) == size(T2, 1) == size(L, 1) == size(z, 1)
    @argcheck size(T1, 2) == size(T2, 2) == size(L, 2) == size(a, 1) == size(A, 1) == size(A, 2)

    T, K = size(L)
    (T == 0) && return

    fill!(T1, 0.0)
    fill!(T2, 0)

    c = 0.0

    for i in OneTo(K)
        T1[1,i] = a[i] * L[1,i]
        c += T1[1,i]
    end

    for i in OneTo(K)
        T1[1,i] /= c
    end

    @inbounds for t in 2:T
        c = 0.0

        for j in OneTo(K)
            # TODO: If there is NaNs in T1 this may
            # stay to 0 (NaN > -Inf == false).
            # Hence it will crash when computing z[t].
            # Maybe we should check for NaNs beforehand ?
            amax = 0
            vmax = -Inf

            for i in OneTo(K)
                v = T1[t-1,i] * A[i,j]
                if v > vmax
                    amax = i
                    vmax = v
                end
            end

            T1[t,j] = vmax * L[t,j]
            T2[t,j] = amax
            c += T1[t,j]
        end

        for i in OneTo(K)
            T1[t,i] /= c
        end
    end

    z[T] = argmax(T1[T,:])
    for t in T-1:-1:1
        z[t] = T2[t+1,z[t+1]]
    end
end

function viterbilog!(T1::AbstractMatrix, T2::AbstractMatrix, z::AbstractVector, a::AbstractVector, A::AbstractMatrix, LL::AbstractMatrix)
    T, K = size(LL)
    (T == 0) && return

    fill!(T1, 0.0)
    fill!(T2, 0)

    al = log.(a)
    Al = log.(A)

    for i in OneTo(K)
        T1[1,i] = al[i] + LL[1,i]
    end

    @inbounds for t in 2:T
        for j in OneTo(K)
            amax = 0
            vmax = -Inf

            for i in OneTo(K)
                v = T1[t-1,i] + Al[i,j]
                if v > vmax
                    amax = i
                    vmax = v
                end
            end

            T1[t,j] = vmax + LL[t,j]
            T2[t,j] = amax
        end
    end

    z[T] = argmax(T1[T,:])
    for t in T-1:-1:1
        z[t] = T2[t+1,z[t+1]]
    end
end


## Convenience functions

# The following methods are defined:
# viterbi(a, A, L)              -> z
# viterbi(hmm, observations)    -> z

"""
    viterbi(a, A, L; logl) -> Vector
Find the most likely hidden state sequence, see [Viterbi algorithm](https://en.wikipedia.org/wiki/Viterbi_algorithm).
"""
function viterbi(a::AbstractVector, A::AbstractMatrix, L::AbstractMatrix; logl = false)
    T1 = Matrix{Float64}(undef, size(L))
    T2 = Matrix{Int}(undef, size(L))
    z = Vector{Int}(undef, size(L,1))
    if logl
        viterbilog!(T1, T2, z, a, A, L)
    else
        warn_logl(L)
        viterbi!(T1, T2, z, a, A, L)
    end
    z
end

"""
    viterbi(hmm, observations; logl, robust) -> Vector
**Example**
```julia
using Distributions, HMMBase
hmm = HMM([0.9 0.1; 0.1 0.9], [Normal(0,1), Normal(10,1)])
y = rand(hmm, 1000)
zv = viterbi(hmm, y)
```
"""
function viterbi(hmm::AbstractHMM, observations; robust = false, kwargs...)
    L = likelihoods(hmm, observations; robust = robust, kwargs...)
    viterbi(hmm.a, hmm.A, L; kwargs...)
end
