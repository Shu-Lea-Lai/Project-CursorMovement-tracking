using PyPlot
function myplot(zv)
    figure()
    plot(zv)
    gcf()
end
DeleteAll(plot_distri) = (Base.delete_method).(methods(plot_distri))
function get_space(a::AbstractArray{Float64,1})
    space = zeros(size(a,1)-1)
    for i in eachindex(a)
        if i!=1
            space[i-1] = a[i]-a[i-1]
        end
    end
    space
end
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
function calc_predflat(curPrd::AbstractArray)
    temp = convert(Array{Float64},curPrd)
    @inbounds for i in 1:(size(curPrd,1)-1)

        if i == 1 continue end
        if i+1 > size(curPrd,1) break end

        if i%2 == 0
            temp[i] = (curPrd[i-1] + curPrd[i+1])/2
        end
    end
    temp
end
function plot_predflat(zv::AbstractDataFrame, rg::AbstractUnitRange)
    Nobs = size(zv.curPrd,1)
    figlen = Int64(round(size(rg,1)*60/Nobs))
    figure(figsize=(figlen,4))
    plot(zv.curs[rg],label="Cursor")
    # plot(zv.curPrd[rg],ls = "--",label = "Prediction")
    plot(zv.Prdflat[rg],label = "Prediction",ls = "--")
    plot(zv.dot[rg],label = "Target")
    legend()
    gcf()
end
function plot_predorg(zv::AbstractDataFrame, rg::AbstractUnitRange)
    Nobs = size(zv.curPrd,1)
    figlen = Int64(round(size(rg,1)*60/Nobs))
    figure(figsize=(figlen,4))
    plot(zv.curs[rg],label="Cursor")
    plot(zv.curPrd[rg],ls = "--",label = "Prediction")
    # plot(zv.Prdflat[rg],label = "Prediction",ls = "--")
    plot(zv.dot[rg],label = "Target")
    legend()
    gcf()
end

function finish_zv!(zv::AbstractDataFrame, space::Int64=100)

    zv.s_curs = s_curs
    zv.a_curs = a_curs
    zv.curs = df.Xpos_cursor;
    zv.stateT = zv.state .- 4
    zv.dot = df.Xpos_dot
    zv.curPrd = zv.curs .+ zv.stateT*space
    zv.Prdflat = calc_predflat(zv.curPrd)
    return nothing
end

"""
Plot distribution and the Sum of it, NOT including the data distr
"""
function plot_dist(distr::AbstractArray{}, rg::StepRangeLen = -150.0:0.05:150.0)
    narg = size(distr,1)
    pdfi = [pdf.(distr[i],rg) for i in 1:narg]
    va = VectorOfArray(pdfi)
    pdfi = convert(Array, va)
    ∑pdf = sum(pdfi,dims=2)
    # ∑pdf ./= sum(∑pdf)
    figure()
    if narg>1
        for i in 1:narg plot(rg,pdfi[:,i]) end
    end
    plot(rg,∑pdf)
    gcf()
end


"""
Plot a signle distribution
"""
function plot_dist(distr::Distribution, rg::StepRangeLen = -150.0:0.05:150.0)

    pdfi = pdf.(distr,rg)
    figure()
    plot(rg,pdfi)
    gcf()
end

"""
Plot distribution and the Sum of it, Including the data distr
"""
function plot_dist(distr::AbstractArray{},datadist::Distribution, rg::StepRangeLen = -150.0:0.05:150.0)
    narg = size(distr,1)
    pdfi = [pdf.(distr[i],rg) for i in 1:narg]
    va = VectorOfArray(pdfi)
    pdfi = convert(Array, va)
    ∑pdf = sum(pdfi,dims=2)
    ∑pdf ./= sum(∑pdf)
    println(sum(∑pdf))
    pdf_data = pdf.(datadist,rg)
    figure()
    if narg>1
        for i in 1:narg plot(rg,pdfi[:,i]) end
    end
    plot(rg,∑pdf,label="sum")
    plot(rg,pdf_data,ls = "-+",label="Data")
    legend()
    gcf()
end

"""
Proportionaly split a number tot in to sperate value given distribution x
The sum of proportion(cdf) euqals to one
This is especially helpful in parameter search for transition matrix

# Examples
```julia-repl
julia> tot = cs_σ;
julia> x = cs_sdNm;
julia> lim = 0.01;
julia> num_p = 7;
julia> tot_split(tot,x,lim,num_p=num_p)
([0.05642617059447854, 0.18982776105109617, 0.2904741897759673, ...
```
"""
function tot_split(tot::AbstractFloat,x::Distribution,lim::AbstractFloat; num_p=7)
    rg = (lim/2,1-lim/2)
    xval_rg = quantile.(x,rg)
    xval_cuts = Cuts(xval_rg,num_p+1)
    cdf_cuts = cdf.(x,xval_cuts)
    cdf_each = get_space(cdf_cuts)
    cdf_each .+= lim/num_p
    x_each = cdf_each .* tot
    cdf_each,x_each
end


function x_split(tot::AbstractFloat,x::Distribution,lim::AbstractFloat; num_p=7)
    rg = (lim/2,1-lim/2)
    xval_rg = quantile.(x,rg)
    xval_cuts = Cuts(xval_rg,num_p+1)
    cdf_cuts = cdf.(x,xval_cuts)
    cdf_each = get_space(cdf_cuts)
    cdf_each .+= lim/num_p
    x_each = cdf_each .* tot
    cdf_each,x_each
end

"""
Proportionaly split a number tot in to sperate value given distribution x
Though,the sum of proportion euqals to one, it take the reverse of proportions
"""
function tot_rev_split(tot::AbstractFloat,x::Distribution,lim::AbstractFloat, num_p=7)
    rg = (lim/2,1-lim/2)
    xval_rg = quantile.(x,rg)
    xval_cuts = Cuts(xval_rg,num_p+1)
    cdf_cuts = cdf.(x,xval_cuts)
    cdf_each = get_space(cdf_cuts)
    cdf_each .+= lim/num_p
    x_each = 1 ./ cdf_each .* tot
    cdf_each,x_each
end

plot_dist(cs_Nms)
cs_save = cs_Nms

Cuts(rg::Tuple,lt::Int) = range(rg[1],rg[2],length=lt)

plot_dist([Nm(-10,10),Nm(-7,10),Nm(-3,10),Nm(-0,10),Nm(3,10)],Nm(-3,2.2))
