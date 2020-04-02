pwd()
cd("HMMBase.jl-master/src")
push!(LOAD_PATH, pwd())
using Multi_HMMBase
using Distributions
using CSV
using DataFrames
using PyPlot
using RecursiveArrayTools
include("JL_UsefulFunctin.jl")

curdic = "/home/shulai/Insync/shulai@iu.edu/Google Drive/IUB/Project-Cursor Movement"#dirname(@__FILE__)
path = string(curdic,"/data")
df = CSV.read(joinpath(path,"labpc_2mins.csv"))|> DataFrame!
df = df[1:2000,:]
s_curs, a_curs = calc_speed(df.Xpos_cursor)

"""############################## Finish Loading ##############################"""
Nobs = size(s_curs,1)
trans_matrix =zeros(7,7)
trans_matrix[4,:] = [1.75/14, 2/14, 2/14, 2.5/14, 2/14, 2/14, 1.75/14]#[1/7 1/7 1/7 1/7 1/7 1/7 1/7]
# trans_matrix[4,:] = [1/7 1/7 1/7 1/7 1/7 1/7 1/7]
for i in 1:7 if i!=4 trans_matrix[i,4] = 1 end end
trans_matrix



cs_μ, cs_σ, cs_Nm, cs_sdNm = get_dist(s_curs)

~,cuts_σ = tot_rev_split(cs_σ,cs_sdNm,0.2)
cuts_μ = x_split(cs_μ, cs_Nm, 0.1)


cs_Nms = Nm.(cuts_μ,cuts_σ./5)
plot_dist(cs_Nms,cs_Nm)
plot_dist(cs_Nms)


distr_matrix = [Nm(-10,10) Nm(-8,8) Nm(-5,5) Nm(0,5) Nm(5,5) Nm(8,10) Nm(10,10);
                Nm(-10,10) Nm(-8,8) Nm(-5,5) Nm(0,5) Nm(5,5) Nm(8,10) Nm(10,10)]
distr_matrix = vcat(reshape(cs_Nms,(1,7)),reshape(cs_Nms,(1,7)))
# distr_matrix = [Nm(-4,1) Nm(-2,.8) Nm(-1,.5) Nm(0,5) Nm(1,.5) Nm(2,.8) Nm(4,1);
# Nm(-4,1) Nm(-2,.8) Nm(-1,.5) Nm(0,5) Nm(1,.5) Nm(2,.8) Nm(4,1)]
curs_HMM = HMM(trans_matrix, distr_matrix)
obs = [s_curs;a_curs];obs = reshape(obs, (2,Nobs))
zv = viterbi(curs_HMM, obs)
myplot(zv)

zv = DataFrame(state = zv)
finish_zv!(zv,100)
myplot(zv.stateT)

# tot = cs_σ
# x = cs_sdNm
# lim = 0.01
# num_p = 7
# rg = (lim/2,1-lim/2)


#
# figure()
# hist(s_curs,bins = 20, label="Cursor speed")
# hist(a_curs,bins = 20, label = "Cursor acceleration")
# legend()
# gcf()
plot_predflat(zv,200:500)
# plot_predorg(zv,200:500)
# figure(figsize = (20,6))




# CSV.write("zv_prd.csv",zv)
