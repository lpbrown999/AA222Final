using CSV
using PyPlot



#Load the file data
filename = "results/result_400_100_3000_10.csv"
file = CSV.read(filename,header=false,transpose=true)
weights = Array(file.Column1)
capacities = Array(file.Column2)
SF = Array(file.Column3)
δ = Array(file.Column4)

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

dominates(y, y′) = all(y′ - y .≥ 0) && any(y′ - y .> 0)
for y in ys
	if !any(dominates(y,y′) for y′ in ys) 
		push!(pareto_ys,y)
	else
		push!(dominated_ys,y)
	end
end

fig = figure("")
ax = gca()
ax.grid()

for y in pareto_ys
	ax.plot(y[1],y[2],"k*")
end
for y in dominated_ys
	ax.plot(y[1],y[2],"r*")
end
ax.set(xlabel="Total Weight", ylabel="-Electrical Capacity", xticklabels=[], yticklabels=[])
# savefig("test.png",dpi=300)


