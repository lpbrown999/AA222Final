using CSV
using PyPlot
using Colors

function main()

	Ps = [0,50,300]


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
		SFs = Array(DF.SafetyFactor)
		δs = Array(DF.Deflection)
		heights = Array(DF.Height)
		widths = Array(DF.Width)
		battery_widths = Array(DF.BatteryWidth)
		replace!(SFs, NaN=>Inf)

		#Sepearate into pareto optimal, non pareto optimal design points 
		ys = []
		for (weight,capacity,SF,δ,height,width,battery_width) in zip(weights,capacities,SFs,δs,heights,widths,battery_widths)
			#Here is where we can check to remove the non-feasible points if we want
			if is_feasible(weight,capacity,SF,δ,height,width,battery_width)
				push!(ys,[weight,-capacity])
			end
		end

		println(P," ",length(ys))
		pareto_ys = []
		dominated_ys = []

		for y in ys
			#If no other point dominates it, it is pareto optimal
			if !any(dominates(y′,y) for y′ in ys)
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

function dominates(y, y′) #returns true if y dominates y'
	return all(y .<= y′) && any(y .< y′)	
end

function is_feasible(weight,capacity,SF,δ,height,width,battery_width)
	#Returns true if all constraints satisfied

	#Constraint values
	hmax = 2.0
	hmin = 1.0
	wmax = 1.0
	wmin = 0.5
	SFmin = 1.1
	SFmax = Inf
	δmax = .005

	#Constraint vector -> Allow slight deviation from constraints
	constraint_vec = [height<=hmax*1.01,
				      height>=hmin*.99,
				      width<=wmax*1.01,
				      width>=wmin*.99,
				      battery_width<=width,
				      SF>=1.0,
				      SF<=SFmax
				      δ<=δmax*1.01]

	return all(constraint_vec)
end

main()






