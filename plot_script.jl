using CSV
using PyPlot
using Colors

function main()

	Ps = 0:100:500
	fig = figure("")
	ax = gca()
	ax.grid()

	#Generate maximally distinguishable colors for this plot
	cols = distinguishable_colors(length(Ps)+1, [RGB(1,1,1)])[2:end]
	pcols = map(col -> (red(col), green(col), blue(col)), cols)

	for (P,pcol) in zip(Ps,pcols)

		# Load the file data
		filename = "results/test_result_$(P)_20000_20.csv"
		DF = CSV.read(filename)

		weights = Array(DF.Weight)
		capacities = Array(DF.Capacity)
		SF = Array(DF.SafetyFactor)
		δ = Array(DF.Deflection)
		height = Array(DF.Height)
		width = Array(DF.Width)
		battery_width = Array(DF.BatteryWidth)

		#Grab the objectives we are trading off
		y1s = weights 
		y2s = -capacities
		ys = []

		#Sepearate into pareto optimal, non pareto optimal design points 
		for (y1,y2) in zip(y1s,y2s)
			push!(ys,[y1,y2])
		end

		pareto_ys = []
		infeasible_pareto_ys = []
		dominated_ys = []
		infeasible_dominated_ys = []

		dominates(y, y′) = all(y .<= y′) && any(y .< y′)	#returns true if y dominates y'
		for y in ys
			if !any(dominates(y′,y) for y′ in ys) 	#If no other point dominates it, it is pareto optimal
				push!(pareto_ys,y)
			else
				push!(dominated_ys,y)
			end
		end

		for y in pareto_ys
			ax.plot(y[1],y[2],c=pcol,marker=".")
		end
		# for y in dominated_ys
		# 	ax.plot(y[1],y[2],c="gray",marker=".")
		# end

		# ax.set(xlabel="Total Weight", ylabel="-Electrical Capacity", xticklabels=[], yticklabels=[])
		# savefig(string(splitext(filename)[1],".png"),dpi=300)
	end
	savefig("test.png",dpi=300)
end

function is_feasible(SF,δ,height,width,battery_width)
	hmax = 2.0
	hmin = 1.0
	wmax = 1.0
	wmin = 0.5
	SFmin = 1.1
	SFmax = Inf
	δmax = .005

end







