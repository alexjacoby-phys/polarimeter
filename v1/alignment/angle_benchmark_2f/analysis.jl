import DelimitedFiles, Plots
using LaTeXStrings

dat = DelimitedFiles.readdlm("angle_data_2f.txt")


Plots.plot(dat[:, 1] * 0.06, dat[:, 2], xlabel="Relative Angle from Micrometer (Degrees)", ylabel="Calculated Absolute Angle (Degrees)", legend=:none, linewidth=2.0)
Plots.savefig("2f_angle_benchmark1.pdf")

Plots.plot(dat[:, 1] * 0.06, dat[:, 2], xlabel="Relative Angle from Micrometer (Degrees)", ylabel="Calculated Absolute Angle (Degrees)", legend=:none, seriestype=:scatter)
Plots.savefig("angle_benchmark2.pdf")



datv2 = [dat[i+1, 2] - dat[i, 2] for i in 1:(length(dat[:, 2])-1)]


datv2 ./ 0.06

Plots.plot(Vector(-9.5:1.0:9.5), datv2 / 0.06, xlabel=L"Increment \ Number", ylabel=L"\frac{\delta Numerical}{\delta Micrometer = 0.06^{\circ}}", ylims=(0.5, 1.5), fontfamily="Times")
Plots.savefig("angle_benchmark3.pdf")