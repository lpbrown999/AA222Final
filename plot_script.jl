using CSV
using PyPlot



# Load the file data
filename = "results/test_result_400_10000_10.csv"
DF = CSV.read(filename)
weights = Array(DF.Weight)
capacities = Array(DF.Capacity)
SF = Array(DF.SafetyFactor)
δ = Array(DF.Deflection)

#Grab the objectives we are trading off
y1s = weights 
y2s = -capacities
ys = []

#Sepearate into pareto optimal, non pareto optimal design points 
for (y1,y2) in zip(y1s,y2s)
	push!(ys,[y1,y2])
end

pareto_ys = []
dominated_ys = []

dominates(y, y′) = all(y .<= y′) && any(y .< y′)	#returns true if y dominates y'
for y in ys
	if !any(dominates(y′,y) for y′ in ys) 	#If no other point dominates it, it is pareto optimal
		push!(pareto_ys,y)
	else
		push!(dominated_ys,y)
	end
end

fig = figure("")
ax = gca()
ax.grid()

for y in pareto_ys
	ax.plot(y[1],y[2],"k.")
end
for y in dominated_ys
	ax.plot(y[1],y[2],"r.")
end

ax.set(xlabel="Total Weight", ylabel="-Electrical Capacity", xticklabels=[], yticklabels=[])
savefig(string(splitext(filename)[1],".png"),dpi=300)


