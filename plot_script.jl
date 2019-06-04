using CSV
using PyPlot
using Colors

function main()
	fontname="Computer Modern"
	rc("text", usetex=true)
	Ps = [100,200,300]
	main_plot_labels = [L"\mathrm{Lowest \ Load}",L"\mathrm{Medium \ Load}",L"\mathrm{Highest \ Load}",""]
	figures = []
	axs = []

	for i in 1:length(Ps)+1
		push!(figures,figure(i,figsize=(16/2,9/2)))
		push!(axs,gca())
		axs[i].set(xlabel=L"\mathrm{Total \ Weight}", ylabel=L"\mathrm{-Electrical \ Capacity}", title=main_plot_labels[i], xticklabels=[], yticklabels=[])
	end
	axs[end].set(title=L"\mathrm{All \ Pareto \ Frontiers}")

	#Generate maximally distinguishable colors for this plot
	cols = distinguishable_colors(length(Ps)+1, [RGB(1,1,1)])[2:end]
	pcols = map(col -> (red(col), green(col), blue(col)), cols)

	for (P,pcol,i) in zip(Ps,pcols,1:length(Ps))

		# Load the file data
		filename = "results_joined_final/result_$(P).csv"
		DF = CSV.read(filename)

		weights = Array(DF.Weight)
		capacities = Array(DF.Capacity)
		SFs = Array(DF.SafetyFactor)
		δs = Array(DF.Deflection)
		heights = Array(DF.Height)
		widths = Array(DF.Width)
		battery_widths = Array(DF.BatteryWidth)
		replace!(SFs, NaN=>Inf)
		println(P," ", length(weights))
		#Sepearate into pareto optimal, non pareto optimal design points 
		ys = []
		pareto_ys = []
		dominated_ys = []
		for (weight,capacity,SF,δ,height,width,battery_width) in zip(weights,capacities,SFs,δs,heights,widths,battery_widths)
			if is_feasible(weight,capacity,SF,δ,height,width,battery_width)
				push!(ys,[weight,-capacity])
			end
		end

		for y in ys
			#If no other point dominates it, it is pareto optimal
			if !any(dominates(y′,y) for y′ in ys)
				push!(pareto_ys,y)
			else
				push!(dominated_ys,y)
			end
		end

		#Plot on main plot and individual plots. Plot one more to get legend.
		println(P," ",length(ys))
		println(P," ",length(pareto_ys))
		for y in dominated_ys
			axs[i].plot(y[1],y[2],c="gray",marker=".",alpha=0.5)
		end
		axs[i].plot(dominated_ys[end][1],dominated_ys[end][2],c="gray",marker=".",alpha=0.5, label=L"\mathrm{Dominated \ Points}")

		for y in pareto_ys
			axs[end].plot(y[1],y[2],c=pcol,marker=".")
			axs[i].plot(y[1],y[2],c=pcol,marker=".")
		end
		axs[end].plot(pareto_ys[end][1],pareto_ys[end][2],c=pcol,marker=".", label=main_plot_labels[i])
		axs[i].plot(pareto_ys[end][1],pareto_ys[end][2],c=pcol,marker=".", label=L"\mathrm{Pareto \ Frontier}")

	end

	for i in 1:length(Ps)+1
		figure(i)
		ax = gca()
		legend()
		# ax.set_aspect(0.0021)
		savefig("figure_$i.png",dpi=500)
	end

	return true
end

function dominates(y, y′) #returns true if y dominates y"
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
				      SF<=SFmax,
				      δ<=δmax*1.2]

	return all(constraint_vec)
end

main()






