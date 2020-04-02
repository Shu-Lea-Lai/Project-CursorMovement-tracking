pwd()
cd("HMMBase.jl-master/src")
push!(LOAD_PATH, pwd())
using Multi_HMMBase
using CSV
using PyPlot
using Distributions
using JuMP,NLopt
using DataFrames
import Dates
include("JL_UsefulFunctin.jl")


zv = CSV.read("zv_prd.csv") |> DataFrame
Nm(a,b) = Normal(a,b)

trans_matrix = Matrix{Float64}(undef, 7, 7)
trans_matrix[4,:] = [1/7 1/7 1/7 1/7 1/7 1/7 1/7]
for i in 1:7 if i!=4 trans_matrix[i,4] = 1 end end
distr_matrix = [Nm(-50,10) Nm(-20,8) Nm(-10,5) Nm(5,5) Nm(10,5) Nm(18,10) Nm(40,20);
                Nm(-50,10) Nm(-20,8) Nm(-10,5) Nm(5,5) Nm(10,5) Nm(18,10) Nm(40,20)]


""" ############# Objective and Constraints #############"""
tot_p = 14
function obj(all_vars::Vector)
    # μ_1, σ_1, μ_2, σ_2, μ_3, σ_3, μ_4, σ_4, μ_5, σ_5,
    # μ_6, σ_6, μ_7, σ_7, μ_8, σ_8, μ_9, σ_9, μ_10, σ_10,
    # μ_11, σ_11, μ_12, σ_12, μ_13, σ_13, μ_14, σ_14,
    # t1, t2, t3, t4, t5, t6, t7 = all_vars
    μ_1, σ_1, μ_2, σ_2, μ_3, σ_3, μ_4, σ_4, μ_5, σ_5,
    μ_6, σ_6, μ_7, σ_7= all_vars

    # trans_matrix = Matrix{Float64}(undef, 7, 7)
    trans_matrix = zeros(7,7)
    # trans_matrix[4,:] = [t1 t2 t3 t4 t5 t6 t7]
    trans_matrix[4,:] = [1/7 1/7 1/7 1/7 1/7 1/7 1/7]
    for i in 1:7 if i!=4 trans_matrix[i,4] = 1 end end
    distr_matrix = [Nm(μ_1,σ_1) Nm(μ_2,σ_2) Nm(μ_3,σ_3) Nm(μ_4,σ_4) Nm(μ_3,σ_3) Nm(μ_2,σ_2) Nm(μ_1,σ_1);
                    Nm(μ_4,σ_4) Nm(μ_5,σ_5) Nm(μ_6,σ_6) Nm(μ_7,σ_7) Nm(μ_6,σ_6) Nm(μ_5,σ_5) Nm(μ_4,σ_4)]

    curs_HMM = HMM(trans_matrix, distr_matrix)
    obs = [zv.s_curs;zv.a_curs];obs = reshape(obs, (2,962))
    zv.state = viterbi(curs_HMM, obs)
    zv.stateT = zv.state .- 4
    zv.curPrd = zv.curs .+ zv.stateT*100
    temp = convert(Array{Float64},zv.curPrd)
    @inbounds for i in 1:(size(zv.curPrd,1)-1)

        if i == 1 continue end
        if i+1 > size(zv.curPrd,1) break end

        if i%2 == 0
            temp[i] = (zv.curPrd[i-1] + zv.curPrd[i+1])/2
        end
    end
    zv.Prdflat = temp

    sum((zv.Prdflat - zv.dot).^2)
end
function nlopt_opt(x::Vector, grad::Vector)
    obj(x)
end


bdd = ((-200.0, 200.0),(1.0, 50.0)) # c_ANpure
toInt(x) = convert(Int64,x)
bdd_lower = [i[1] for i in bdd]; bdd_lower = repeat(bdd_lower,toInt(tot_p/2))
bdd_upper = [i[2] for i in bdd]; bdd_upper = repeat(bdd_upper,toInt(tot_p/2))
# append!(bdd_lower,[0.001 0.001 0.001 0.001 0.001 0.001 0.001])
# append!(bdd_upper,[0.999 0.999 0.999 0.999 0.999 0.999 0.999])
function random_start()
    temp = [rand(Uniform(bdd_lower[i],bdd_upper[i])) for i in 1:tot_p]|>Array
    # append!(temp, repeat([1/7],7))
end
function myconstraint1(x::Vector, grad::Vector)
    return(sum(x[35-7:35])-1)
end

opt = Opt(:LN_COBYLA, tot_p)
# opt = Opt(:LD_SLSQP, tot_p)
opt.lower_bounds = bdd_lower
opt.upper_bounds = bdd_upper

opt.xtol_rel = 1e-10
# opt.maxtime = 60*2
# opt.maxeval = 10000
opt.min_objective = nlopt_opt


# equality_constraint!(opt,(x,g) -> myconstraint1(x,g), 1e-8)


time1=Dates.now()
(minf,minx,ret) = optimize(opt, random_start())
println(Dates.now()-time1)


numevals = opt.numevals
println("got $minf at $minx after $numevals iterations (returned $ret)")

figure(figsize=(10,5))
rg = 100:500
plot(zv.curs[rg])
# plot(zv.curPrd[rg])
plot(zv.Prdflat[rg])
plot(zv.dot[rg])
gcf()

# save1 = minx
