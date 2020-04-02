
# using Distributions
# import Plots

xs = -5.0:0.01:5.0
Plots.plot(xs, pdf.(Normal(), xs)+pdf.(Normal(), xs), legend=nothing)
Plots.ylabel!("\$f_X(x)\$")
Plots.xlabel!("\$x\$")


xs = -5.0:0.01:5.0
Plots.plot(xs, pdf.(Normal(), xs), legend=nothing)
Plots.ylabel!("\$f_X(x)\$")
Plots.xlabel!("\$x\$")
# Plots.title!("Gaussian mixture PDF")



plot_distri(Normal(0,10),Normal(10,10))
