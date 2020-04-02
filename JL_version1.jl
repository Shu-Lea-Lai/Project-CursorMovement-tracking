ddd# # Basic Usage

using Distributions
using HMMBase
using PyPlot
using Seaborn

rc("axes", xmargin = 0) # hide
set_style("whitegrid")  # hide

# ### Model Specification

# B can contains any probability distribution from the `Distributions` package

a = [0.6, 0.4]
A = [0.9 0.1; 0.1 0.9]
B = [MvNormal([0.0, 5.0], ones(2)*1), MvNormal([0.0, 5.0], ones(2)*3)]
hmm = HMM(a, A, B)
size(hmm) # (number of states, observations dimension)

# ### Sampling

z, y = rand(hmm, 500, seq = true)

# Let's plot the observations and the hidden state sequence:


# ### Parameters Estimation

hmm = HMM(randtransmat(2), [MvNormal(rand(2), ones(2)), MvNormal(rand(2), ones(2))])
hmm, hist = fit_mle(hmm, y, display = :iter, init = :kmeans)
hmm
#-

figure(figsize = (4,3)) # hide
plot(hist.logtots)
gcf() # hide
