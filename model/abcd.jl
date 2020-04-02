#=====
LTM global
======#

function assign_LTM_global(Filecondi,item_condi, walk, item)

    global F,L,Fnow,Lnow
    if Filecondi != "VMpure"
        if item_condi == "CM"

            if walk * "i" * item in ["oldiold","oldinew"]

                Fnow = F[Filecondi*"_CM_oldiold_oldinew"]
            else Fnow = 0
            end

            if walk * "i" * item in ["oldiold","newinew"]

                Lnow = L[Filecondi*"_CM_oldiold_newinew"]
            elseif walk * "i" * item in ["oldinew","newiold"]

                Lnow = L[Filecondi*"_CM_oldinew_newiold"]
            else Lnow=0
            end

        elseif item_condi == "AN"

            if walk * "i" * item in ["oldiold","oldinew"]

                Fnow = F[Filecondi*"_AN_oldiold_oldinew"]
            else Fnow=0
            end

            Lnow = 0

        else
            println("wrong condi 1",Filecondi,item_condi, walk*"i"*item)
        end
    elseif Filecondi == "VMpure"

        if item_condi == "VM"
            if walk*"i"*item in ["oldiold","oldinew"]
                Fnow = F["VMpure_VM_oldiold_oldinew"]
            else
                Fnow = 0
            end
        else print("wrong condi 2")
        end

        Lnow = L["VMpure"]
    else
        print("wrong filecondi")
    end
    return(Fnow + Lnow)
end


function calcA(df, all_vars)

    alpha, beta, s, c ,Old_crit ,New_crit , t0 ,kao , F, L = all_vars

    dnow = df
    # a = Array{Union{Nothing,Float64}}(nothing,size(dnow)[1], 8)  #activation
    a =   zeros(Float64, size(dnow)[1], 8)
    betanow=beta["all"]
    alphanow=alpha["all"]

    m = [(j^(-betanow) + alphanow) for j in 1:8]
    for j in 1:8

        if vary_ss==0
            a[(dnow.Lag .== j )|> Array,j] .= m[j]
            a[(dnow.Lag .!= j )|> Array,j] .= m[j] .* s["all"]
        elseif vary_ss==1
            a[((dnow.Lag .!= j) .& (dnow.Setsize .== 2))  |> Array,j] .= m[j] .* s["ss2"]
            a[((dnow.Lag .!= j) .& (dnow.Setsize .== 4))  |> Array,j] .= m[j] .* s["ss4"]
            a[((dnow.Lag .!= j) .& (dnow.Setsize .== 8))  |> Array,j] .= m[j] .* s["ss8"]
            a[(dnow.Lag .== j) |> Array,j] .= m[j]

        end
    end

    for i in 1:size(a)[1]

        a[i,dnow.Setsize[i]+1:end] .= .0 #a_ij suit for the correct amount of setsize
    end

    debug = 0
    if debug==1
        for i in 1:size(a)[1]
            println(i, "\nbeginn","\na is",a[i,:],
                  "\n m is", m ,"n Probtype is",
                  "\n Lag is ",dnow.Lag[i],
                  "\n Setsize is", dnow.Setsize[i],
                  "\n Probtype:",
                  dnow.Probtype[i],
                  "\n Oldnew: ",dnow.Oldnew[i],
                  "\n Ai is", sum(a[i,:], dims = 2),
                  "\n snow", s["ss2"]," ",s["ss4"]," ",s["ss8"],
                 "\n------------------------------------------------" )
        end
    end

    A = sum(a, dims = 2)

    return(A)
end


# F["AN_oldiold_oldinew"] = 0.099
function calcp(df,all_vars)

    alpha, beta, s, c ,Old_crit ,New_crit , t0 ,kao , F, L = all_vars

    A = calcA(df, all_vars)
    dnow = df
    p = Iterators.repeated(3.0, size(A)[1]) |> collect
    filecondis = CategoricalArray(dnow.FileCondi) |> levels; #get categories name

    for icondi in filecondis

        #only ieterate through probs in that condition
        probs = filter( t -> t[:FileCondi] == icondi, dnow).Probtype |> CategoricalArray |> levels
        for iprob in probs

            for ion in ["old","new"]

                tf_ion = (dnow.Oldnew .== ion)
                tf_iprob = (dnow.Probtype .== iprob)
                tf_icondi = (dnow.FileCondi .== icondi)

                tf_combi1 = (tf_ion .& tf_iprob) |> Array
                tf_all = (tf_combi1 .& tf_icondi) |> Array

                IR_old_current = assign_LTM_global(icondi,iprob,"old",ion)
                IR_new_current = assign_LTM_global(icondi, iprob,"new",ion)

                p[tf_all] = (A[tf_all] .+ IR_old_current)./
                (A[tf_all] .+ IR_old_current .+ c[icondi] .+ IR_new_current)
            end
        end
    end

    return(p)

end


function calc_theoretical_RW(df, which_return, all_vars)

    alpha, beta, s, c ,Old_crit ,New_crit , t0 ,kao , F, L = all_vars

    A = calcA(df, all_vars)
    p = calcp(df, all_vars)
    q = 1 .- p

    dnow = df
    filecondis = CategoricalArray(dnow.FileCondi) |> levels; #get categories name

    p_resp_old = zeros(size(p)[1])
    pred_correct = zeros(size(p)[1])
    pred_rt = zeros(size(p)[1])

    theta1 = zeros(size(p)[1])
    theta2 = zeros(size(q)[1])
    theta11 = zeros(size(p)[1])
    theta22 = zeros(size(q)[1])
    exp_nstep = zeros(size(p)[1])


    for ifile in filecondis

        gen_tf = (dnow.FileCondi .== ifile)|> Array
        AA = Old_crit[ifile]
        BB = -New_crit[ifile] # bb is a postive number

        if which_return=="crt"

            pq = (p[gen_tf]./q[gen_tf])
            theta1[gen_tf] = (pq .^(AA+BB) .+1)./(pq .^(AA+BB) .-1)
            theta2[gen_tf] = (pq .^BB .+1)./(pq .^BB .-1)

            tf = ((p .!= q) .& ((dnow.Oldnew .== "old").& gen_tf )) |> Array
            exp_nstep[tf] = (1 ./ (p[tf] .-q[tf])).*(theta1[tf].*(AA+BB) .- theta2[tf].*BB)

            tf = ((p .== q) .& ((dnow.Oldnew .== "old").& gen_tf)) |> Array
            exp_nstep[tf] .= (AA/3)*(2*BB+AA)

            theta11[gen_tf] = (pq .^ (-(AA+BB)) .+1) ./ (pq.^(-(AA+BB)) .-1)
            theta22[gen_tf] = (pq .^ -AA .+1)./(pq .^ -AA .-1)

            tf = ((p .!= q) .& (dnow.Oldnew .== "new") .& gen_tf) |> Array
            exp_nstep[tf] = (1 ./(q[tf] .- p[tf])).*(theta11[tf].*(AA+BB) .- theta22[tf].*AA)

            tf = ((p .== q) .& (dnow.Oldnew .== "new") .& gen_tf) |> Array
            exp_nstep[tf] .= (BB/3)*(2*AA+BB);

        #------------------- correct response
        elseif which_return == "err"

            qp = (q[gen_tf]./p[gen_tf])
            qptfn = ((p .!= q) .& gen_tf) |> Array
            qptf = ((p .== q) .& gen_tf) |> Array

            p_resp_old[qptfn] = ((1 .- qp.^BB) ./ (1 .-qp.^(AA+BB)))[Array(p[gen_tf] .!= q[gen_tf])];
            p_resp_old[qptf] .= BB/(AA +BB);

        end
    end
    #------------------


    if which_return=="crt"
        suprise = (((dnow.FileCondi.=="MIX") .| (dnow.FileCondi.=="MIX2")) .&
                (dnow.Probtype.=="AN") .& (dnow.Oldnew.=="new") )|> Array

        notsuprise = .!suprise

        pred_rt[suprise] = t0["ann"] .+ kao .* exp_nstep[suprise]
        pred_rt[notsuprise] = t0["all"] .+ kao .* exp_nstep[notsuprise]
    end

    #------------------

    if which_return == "err"
        pred_correct[Array(dnow.Oldnew.=="old")] = p_resp_old[Array(dnow.Oldnew.=="old")]
        pred_correct[Array(dnow.Oldnew.=="new")] = 1 .- p_resp_old[Array(dnow.Oldnew.=="new")]

    end


    if which_return == "crt"
        return(pred_rt)
    elseif which_return == "err"
        return(pred_correct)
    end

end


function calc_wssd(df_err,df_crt, all_vars)

    w=Dict()

    w["new_rt"] = 4*2
    w["old_rt"] = 1*2
    w["new_err"] = 4
    w["old_err"] = 1

    alpha, beta, s, c ,Old_crit ,New_crit , t0 ,kao , F, L = all_vars


    df_sub_err = df_err |> copy
    df_sub_crt = df_crt |> copy

    if search_MIX1 == 1
        df_sub_err = filter(t -> t[:FileCondi]=="MIX", df_sub_err)
        df_sub_crt = filter(t -> t[:FileCondi]=="MIX", df_sub_crt)
    end

    pred_correct = calc_theoretical_RW(df_sub_err,"err", all_vars)
    pred_crt = calc_theoretical_RW(df_sub_crt,"crt", all_vars)

    df_sub_crt.pred_crt = pred_crt ./ 1000 #translate crt to seconds.
    df_sub_err.pred_error = 1 .- pred_correct

    df_sub_crt.RT = df_sub_crt.RT ./1000

    df_sub_crt.SSD_RT = ((df_sub_crt.RT .- df_sub_crt.pred_crt) .^2)
    df_sub_err.SSD_err = ((df_sub_err.Error .- df_sub_err.pred_error).^2)

    df_sub_crt.wSSD_RT = [df_sub_crt.SSD_RT[i] * w[df_sub_crt.Oldnew[i] * "_rt"]
        for i in 1:size(df_sub_crt)[1]]

    df_sub_err.wSSD_err = [df_sub_err.SSD_err[i] * w[df_sub_err.Oldnew[i] * "_err"]
        for i in 1:size(df_sub_err)[1]]

    WSSD = (df_sub_err.wSSD_err |> sum)  .+ (df_sub_crt.wSSD_RT |> sum)

    return(WSSD)
end
