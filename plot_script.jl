using CSV
using PyPlot

filename = "result_400_100_3000_10.csv"
file = CSV.read(filename,header=false,transpose=true)
weights = Array(file.Column1)
capacities = Array(file.Column2)

fig = figure("")
ax = gca()
ax.grid()
ax.plot(weights,-capacities,"b*")
ax.set(xlabel="Total Weight", ylabel="-Electrical Capacity", xticklabels=[], yticklabels=[])
savefig("test.png",dpi=300)