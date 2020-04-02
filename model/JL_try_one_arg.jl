
using HMMBase
using Distributions
using CSV
using DataFrames
using PyPlot

# hmm = HMM([0.9 0.1; 0.1 0.9], [Normal(0,1) Normal(2,1); Normal(2,1) Normal(0,1)])
# y = rand(hmm, 1000)
# likelihoods(hmm,y)
# zv = viterbi(hmm, y)
# myplot(zv)

curdic = "/home/shulai/Insync/shulai@iu.edu/Google Drive/IUB/Project-Cursor Movement"#dirname(@__FILE__)
path = string(curdic,"/data")
df = CSV.read(joinpath(path,"simple1_jacob.csv"))|> DataFrame!

function calc_speed(Xpos_cursor)
    s = Vector{Float64}(undef,size(Xpos_cursor,1))
    a = Vector{Float64}(undef,size(Xpos_cursor,1))
    s[1] = 0; a[1] = 0
    for i in 2:size(s,1)
        s[i] = Xpos_cursor[i] - Xpos_cursor[i-1]
        if i==2
            a[2] = 0
        else
            a[i] = s[i]-s[i-2]
        end
    end
    s,a
end
s_curs, a_curs = calc_speed(df.Xpos_cursor)

trans_matrix = [0.30 0.20 0.15 0.15 0.10 0.05 0.05;
                0.18 0.25 0.18 0.15 0.10 0.07 0.07;
                0.10 0.20 0.25 0.20 0.10 0.10 0.05;
                0.05 0.10 0.20 0.30 0.20 0.10 0.05;
                0.05 0.10 0.10 0.20 0.25 0.20 0.10;
                0.07 0.07 0.10 0.15 0.18 0.25 0.18;
                0.05 0.05 0.10 0.15 0.15 0.20 0.30]

Nm(a,b) = Normal(a,b)
distr_matrix = [Nm(-50,10) Nm(-20,8) Nm(-10,5) Nm(5,5) Nm(10,5) Nm(18,10) Nm(40,20)]
distr_matrix = reshape(distr_matrix,(7))
curs_HMM = HMM(trans_matrix, distr_matrix)
obs = s_curs
zv = viterbi(curs_HMM, obs)

function myplot(zv)
    figure()
    plot(zv)
    gcf()
end
myplot(zv)
# figure() # hide
# plot(df.time,df.Xpos_cursor)
# plot(df.time,df.Xpos_dot)
# gcf() # h
#
# figure() # hide
# plot(s_curs[1:100])
# plot(a_curs[1:100])
# plot(df.Xpos_cursor[1:100])
# gcf() # h
#
#
# figure()
# hist(s_curs[200:300],bins = 20)
# hist(a_curs[200:300],bins = 20)
# gcf()
