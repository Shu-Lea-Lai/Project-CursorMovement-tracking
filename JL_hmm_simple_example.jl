

using Distributions
using HMMBase


hmm = HMM([0.9 0.1; 0.1 0.9], [Normal(0,1), Normal(10,1)])
y = rand(hmm, 1000)

using ArgCheck
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

function warn_logl(L::AbstractMatrix)
    if any(L .< 0)
        @warn "Negative likelihoods values, use the `logl = true` option if you are using log-likelihoods."
    end
end
function viterbi_aa!(T1::AbstractMatrix, T2::AbstractMatrix, z::AbstractVector, a::AbstractVector, A::AbstractMatrix, L::AbstractMatrix)

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

function viterbi_aa(a::AbstractVector, A::AbstractMatrix, L::AbstractMatrix; logl = false)
    T1 = Matrix{Float64}(undef, size(L))
    T2 = Matrix{Int}(undef, size(L))
    z = Vector{Int}(undef, size(L,1))
    if logl
        viterbilog!(T1, T2, z, a, A, L)
    else
        warn_logl(L)
        viterbi_aa!(T1, T2, z, a, A, L)
    end
    z
end

function viterbi_aa(hmm::AbstractHMM, observations; robust = false, kwargs...)
    L = likelihoods(hmm, observations; robust = robust, kwargs...)
    viterbi_aa(hmm.a, hmm.A, L; kwargs...)
end
zv = viterbi_aa(hmm, y)
