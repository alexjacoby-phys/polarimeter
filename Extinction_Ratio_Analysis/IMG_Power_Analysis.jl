import Images, FileIO, LinearAlgebra

home = "/Users/alexjacoby/Documents/Research_Code/polarimeter/Extinction_Ratio_Analysis"
dat = "/Volumes/AJACOBY/extinction/"



function array_mean(M::Array{Float64})
    return sum(M)/*(size(M)...)
end

fn_vec = [string(dat,i,"_exp=100ms.png") for i in -15:15]

pre_dat = Images.load.(fn_vec)


dat = [Float64.(Images.Gray.(img)) for img in pre_dat]

plotdat = array_mean.(dat)



import Plots

Plots.plot(-15:15,plotdat,seriestype=:scatter)



import DelimitedFiles

cd(home)

pm_dat = DelimitedFiles.readdlm("Power_Meter_Low.txt")


plotdat[1:2:15]

Plots.plot(pm_dat[:,1],pm_dat[:,2],seriestype = :scatter)



mult_rng = 0.:0.001:10.



optim = [sum((pm_dat[1:8,2]- val*plotdat[16:2:30]) .^2) for val in mult_rng]

findmin(optim)





Plots.plot(pm_dat[1:8, 2], label = "Power Meter");
Plots.plot!(0.225*plotdat[16:2:30] .+ 0.055,label = "Camera")


Plots.plot(-15:15,0.225*plotdat .+0.055,label = "Camera")
Plots.plot!(pm_dat[:,1],pm_dat[:,2],label = "pwr mtr")

Plots.plot(pm_dat[1:8, 2] - plotdat[16:2:30])



(0.225 * plotdat .+ 0.055)[16]


