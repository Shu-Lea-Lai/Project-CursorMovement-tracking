
using CSV
using PyPlot
include("JL_UsefulFunctin.jl")

zv = CSV.read("zv.csv")
curdic = "/home/shulai/Insync/shulai@iu.edu/Google Drive/IUB/Project-Cursor Movement"#dirname(@__FILE__)
path = string(curdic,"/data")
df = CSV.read(joinpath(path,"simple1_jacob.csv"))|> DataFrame!



zv.curs = df.Xpos_cursor;
zv.stateT = zv.state .- 4
zv.dot = df.Xpos_dot
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

figure(figsize = (20,6))
# figure()
rg = 1:962
plot(zv.curs[rg],label="Cursor")
# plot(zv.curPrd[rg],ls = "--",label = "Prediction")
plot(zv.Prdflat[rg],label = "Prediction",ls = "--")
plot(zv.dot[rg],label = "Target")
legend()
gcf()

myplot(zv.stateT)

CSV.write("zv_prd.csv",zv)
