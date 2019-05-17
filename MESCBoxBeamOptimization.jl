using AA222Final
using CSV, DataFrames
# using PyPlot

# all units in lb - in

## Notes
# 100,000 = n_eval_outer*n_eval_nner takes about 20 seconds for 1 capacity.
# Maximum feasible capacity to reach is ~ 91.326 (hmax*wmax*l*ρe = 2*1*9.3*4.91 = 91.326)
# If we set higher than this no way it will be able to reach it

P = 400 					#Tip load (integer)
n_eval_outer = 100			#Number of evaluations used to determine number of plies, doesnt need as many
n_eval_inner = 3000 		#Number of evaluations used to evaluate the design for number of plies, needs more
num_growths = 10			#Number of growths
Capacitymins = 0:1:90	    #capacities to iterate over
csv_name = "result_$(P)_$(n_eval_outer)_$(n_eval_inner)_$(num_growths).csv"

function mesc_box_beam_objective(input_vector::Vector; num_ply_vec::Vector, constraint_vec::Vector, allowable_ply_angles::Vector, p1::Float64, p2::Float64, mode::String)
	global P

	##################### DESCRIPTION / extract inputs #####################
	#num_ply_vec is a vector of the number of plies in the top flange, web, and bottom flange
	n_top_flange = num_ply_vec[1]
	n_webs       = num_ply_vec[2]
	n_bot_flange = num_ply_vec[3]
	
	## input_vector is a vector of beam geometry and layups, which is what we are optimizing over
		#Element 1: beam total height (enforce feasibility)
		#Element 2: beam total width  (enforce feasibility)
		#Element 3: battery width     (enforce feasibility)
		#Elements 4 to end: layups 
	h = max(input_vector[1],.05)
	w = max(input_vector[2],.05)
	wb = max(input_vector[3],.05)

	layups = input_vector[4:end]
	top_layup = layups[1:n_top_flange]
	web_layup = layups[1+n_top_flange:n_top_flange+n_webs]
	bot_layup = layups[1+n_top_flange+n_webs:n_top_flange+n_webs+n_bot_flange]

	## constraint_vec is a vector of constraints
		#Element1 : maximum height
		#Element2 : minimum height
		#Element3 : maximum width 
		#Element4 : minimum width
		#Element5 : minimum safety factor
		#Element6 : maximum safety factor -> over designed
		#Element7 : maximum defleciton
		#Element8 : minCapacity, the minimum level of the the negative capacity since we are doing pareto curve, maximizing capactiy
	hmax = constraint_vec[1]
	hmin = constraint_vec[2]
	wmax = constraint_vec[3]
	wmin = constraint_vec[4]
	SFmin = constraint_vec[5]
	SFmax = constraint_vec[6]
	δmax = constraint_vec[7]
	Capacitymin = constraint_vec[8]

	#Symmetry is also imposed

	## p1, p2 are penalty factors

	## Mode is what we are doing really.
	# options: "optimize" returns weight + constraint penalties 
	#		   "evaluate" returns height, width, weight, capacity, SF


	##################### Start main function #######################
	######Fixed Stuff#########
	l = 9.3		#Fixed length of 30 inches, can change

	## Battery properties, assume isotropic, homogenized behavior. Does not include out of plane effects
	E_bat = .78e6		#In plane
	G_bat = .304e6 		#In plane 
	ρe_bat = 4.91		#Energy density (300 W-h/L = 4.91 W-h/in^3)
	ρ_bat = .015		#Realize that also given as 130 W-H / kg -> 130/300 kg/L -> lighter than carbon prepreg

	battery_properties = BatteryProperties(E_bat,G_bat,ρ_bat,ρe_bat)

	##Composite properties.
	# https://www.rockwestcomposites.com/materials-tools/fabrics-pre-pregs-tow/prepregs/14033-d-group
	E1 = 8.26e6
	E2 = 7.83e6
	ν12 = .1
	G12 = .71e6
	tply = .01
	ρ = .056
	composite_properties = MaterialProperties(E1,E2,ν12,G12,tply,ρ)

	F1t = 70e3
	F2t = 66e3
	F1c = 55e3
	F2c = 55e3
	F12 = 10e3 
	composite_strengths = MaterialStrengths(F1t,F2t,F1c,F2c,F12)

	##Form laminates, then plates
	#Top and bot laminate are as wide as the whole boxbeam
	#The web plates are as wide as the height of the beam minus the height of the top and bottom laminates
	top_laminate  = laminate_analyzer(top_layup,composite_properties,composite_strengths)
	top_plate = composite_plate(top_laminate, w, l)		

	bot_laminate  = laminate_analyzer(bot_layup,composite_properties,composite_strengths)
	bot_plate = composite_plate(bot_laminate, w, l)		

	web_laminate = laminate_analyzer(web_layup,composite_properties,composite_strengths)
	web_plate = composite_plate(web_laminate, h-top_plate.h-bot_plate.h, l)


	#generate MESC box beam and can evaluate its performance in the loading case	
	#Enforce physical meaning of battery width. Cannot be be larger than the beam.	
	wb = min(wb, w-2*web_plate.h)
	box_beam = MESC_box_beam(top_plate,web_plate,bot_plate,battery_properties,h,w,l,wb)
	bending_case = cantilever_bending(box_beam,float(P))

	##################### Value return #####################
	# all constrants of form <= 0
	# constraint_vec1 -> non layup constraints
	# constraint_vec2 -> layup orientation constraints within allowable angles, and symmetry penalties
	SF = bending_case.SF
	δ = bending_case.δ
	capacity = box_beam.capacity
	weight   = box_beam.weight
	
	cvec1 = [h-hmax, 						#h <= hmax 
			 w-wmax, 						#w <= wmax 
			 wb-w+2*box_beam.webs.h,		#wbattery <= w - 2*webthickness	#Battery must fit inside!

			 hmin-h, 						#h >= hmin
			 wmin-w,						#w >= wmin  

			 (SFmin-bending_case.SF), 		#SF >= SFmin 
			 (bending_case.SF-SFmax),		#SF <= SFmax
			 δ-δmax, 						#δ <= δmax
			 Capacitymin-capacity]			#Capacity >= Capacitymmin
	cvec1_weights = [1.,1.,1.,1.,1.,100.,100.,1.,1.]	#Constraint weighting vector


	cvec2 = []
	for layer in layups	
		push!(cvec2,minimum(abs.(layer.-allowable_ply_angles)))		# c evaluates to the minimum distance from an allowed angle
	end
	cvec2 = [cvec2; symmetry_penalty(top_layup); symmetry_penalty(web_layup); symmetry_penalty(bot_layup)]

	#Objective + penalties on the constraints
	#no count penalty on the layup orientation since we will never be exact
	if mode == "optimize"
		return weight + combined_penalty(cvec1.*cvec1_weights,p1,p2) + combined_penalty(cvec2,p1/10,0)

	elseif mode == "evaluate"
		return [weight,capacity,SF,δ]
	
	elseif mode == "check_constraints"
		return cvec1,cvec2
	elseif mode == "debug"
		return bending_case
	else
		return Nothing
	end
end

function symmetry_penalty(layup::Vector)
	#Returns a vector penalizing the difference between the ply angles
	nplies = length(layup)
	if nplies == 1
		return 0.

	elseif isodd(nplies)
		idx = Int((nplies-1)/2)
		front_half = layup[1:idx]
		back_half = layup[idx+2:end]

	elseif iseven(nplies)
		idx = Int(nplies/2)
		front_half = layup[1:idx]
		back_half = layup[idx+1:end]

	end
	return abs.(front_half .- reverse(back_half))
end

function fix_num_ply_vec(num_ply_vec::Vector)
	return floor.(Int,max.(num_ply_vec,1.0)) 
end

function optimize(f, ndim; neval=100, a=-1., b=1., p1=1., p1g=2., p2=1., p2g=2., num_growths=1)
	# Inputs
	#f -> function to optimize. Should include constraints already, and needs to be fed p1, p2.
	#ndim -> problem dimension
	#neval -> number of evals total
	#a -> vector or float64 that is lower bound for simplex generation
	#b -> vector or float64 that is upper bound for simplex generation 
	#p1,p2 -> qudratic, count penalty factors. p1g,p2g is growth factor on penalty
	#num_growths -> number of times to increase the penalty and randomly restart while preserving the current best
	
	#Generate initial simplex
	S = generate_simplex(ndim,a,b)
	xbest = zeros(ndim)

	#Run nelder mead, random restart, increase penalty num_growths number of times
	for i in 1:num_growths
		f_penalized(input_vector) = f(input_vector,p1,p2)
		_, xbest = nelder_mead(f_penalized, S, floor(neval/num_growths))

		#Generate new random simplex, inserting best point into it. Gives random restart effect while preserving progress.
		S = generate_simplex(ndim,a,b)
		S[1] = xbest
		p1 *= p1g
		p2 *= p2g
	end
	return xbest
end

#This function will return an optimal design value for the given number of plies and constraint penalties
#This works by optimizing the geometry, layups for the given number of plies.
function design_for_num_ply(num_ply_vec::Vector; constraint_vec::Vector, allowable_ply_angles::Vector, p1::Float64, p2::Float64, n_eval::Int64, num_growths::Int64)
	
	#Inputs: 
	#num_ply_vec: Number of plies in the top flange, webs, bottom flange 
	#constraint vec which is fed into the MESC box beam solver 
	#allowable ply angles which is alos fed into the MESC box beam solver 
	#p1, p2 are penalties fed into the box beam solver when evaluating design.
	
	#Begin
	num_ply_vec = fix_num_ply_vec(num_ply_vec)  #Make this integers

	#Set up inner optimization objective
	a_inner = [constraint_vec[2]; constraint_vec[4]; constraint_vec[4]; -90. *ones(sum(num_ply_vec))]	#LB on h,w,wb, ply angles for simplex generation
	b_inner = [constraint_vec[1]; constraint_vec[3]; constraint_vec[3]; 135. *ones(sum(num_ply_vec))]	#UB
	ndim_inner = length(a_inner)

	#Optimize the beam for the given number of plies
	geometry_layup_objective(input_vector,p1,p2) = mesc_box_beam_objective(input_vector, p1=p1, p2=p2, num_ply_vec=num_ply_vec, constraint_vec=constraint_vec, allowable_ply_angles=allowable_ply_angles, mode="optimize")
	optimal_design = optimize(geometry_layup_objective, ndim_inner, neval=n_eval, a=a_inner, b=b_inner, p1=1.0, p1g=1.3,  p2=10., p2g=1.3, num_growths=num_growths)	
	
	design_value = geometry_layup_objective(optimal_design,p1,p2)
	return design_value
end

#Creates pareto frontier of weight and capacity by iteratively constraining capacity.
function pareto_weight_capacity(Capcitymins, n_eval_outer::Int64, n_eval_inner::Int64, num_growths::Int64)
	## Inputs
	#Capacitymins -> vector of minimum capacities to run
	#n_eval_outer -> number of evaluations used to determine number of plies
	#n_eval_inner
	#num_growths

	#Set up constraints
	hmax = 2.0
	hmin = 1.0
	wmax = 1.0
	wmin = 0.5
	SFmin = 1.1
	SFmax = 1.5
	δmax = 1.0
	allowable_ply_angles = [-45.,0.,45.,90.]


	#Storage
	best_num_ply_vecs = []	
	best_designs = []
	best_designs_evaluations = []
	ys = []

	#Within each loop, find the optimal design for this minimum capacity
	for Capacitymin in Capacitymins
		println(Capacitymin)
		constraint_vec = [hmax,hmin,wmax,wmin,SFmin,SFmax,δmax,Capacitymin]

		#Set up the outer optimization objective that will give us the best number of ply distribution
		#Get the best ply distribution for the given constraints
		ndim_outer = 3
		a_outer = 1.
		b_outer = 5.
		num_ply_objective(num_ply_vec,p1,p2) = design_for_num_ply(num_ply_vec, p1=p1, p2=p2, constraint_vec=constraint_vec, allowable_ply_angles=allowable_ply_angles, n_eval = n_eval_inner, num_growths=num_growths)
		num_ply_vec = fix_num_ply_vec( optimize(num_ply_objective, ndim_outer, neval=n_eval_outer, a=a_outer, b=b_outer, p1=1.0, p1g=1.3, p2=10., p2g=1.3, num_growths=num_growths) )

		#Set up inner optimization objective, copied from the outer objective function
		#Optimize the beam for the given number of plies
		a_inner = [hmin; wmin; wmin; -90. *ones(sum(num_ply_vec))]	#LB on h,w,wb, ply angles for simplex generation
		b_inner = [hmax; wmax; wmax; 135. *ones(sum(num_ply_vec))]	#UB
		ndim_inner = length(a_inner)
		geometry_layup_objective(input_vector,p1,p2) = mesc_box_beam_objective(input_vector, p1=p1, p2=p2, num_ply_vec=num_ply_vec, constraint_vec=constraint_vec, allowable_ply_angles=allowable_ply_angles, mode="optimize")
		optimal_design = optimize(geometry_layup_objective, ndim_inner, neval=10000, a=a_inner, b=b_inner, p1=1.0, p1g=1.3, p2=10., p2g=1.3, num_growths=100)	
		

		#Evalutate the design, store the weight and capacity of it in design_objectives
		design_evaluation = mesc_box_beam_objective(optimal_design, p1=1.0, p2=1.0, num_ply_vec=num_ply_vec, constraint_vec=constraint_vec, allowable_ply_angles=allowable_ply_angles, mode="evaluate")
		push!(best_num_ply_vecs, num_ply_vec)
		push!(best_designs, optimal_design)

		#Ys is weight, -capacity since capacity is maximized and this works for minimizing.
		push!(ys, [design_evaluation[1], -design_evaluation[2]] )
		push!(best_designs_evaluations, design_evaluation)				#For making sure designs are feasible
	end

	#Get only non-dominated points
	pareto_ys = []
	pareto_ply_vec = []
	pareto_design = []
	pareto_design_evaluations = []

	for (y, num_ply_vec, design, design_evaluation) in zip(ys,best_num_ply_vecs,best_designs, best_designs_evaluations)
		if !any(all(y′ - y .≥ 0) && any(y′ - y .> 0) for y′ in ys) #If it is not dominated by any point, add it to list
			push!(pareto_ys,y)
			push!(pareto_ply_vec,num_ply_vec)
			push!(pareto_design,design)
			push!(pareto_design_evaluations, design_evaluation)
		end
	end

	return pareto_ys, pareto_ply_vec, pareto_design, pareto_design_evaluations
end

pareto_ys, pareto_ply_vec, pareto_design, pareto_design_evaluations = pareto_weight_capacity(Capacitymins, n_eval_outer, n_eval_inner, num_growths)
CSV.write(csv_name,  DataFrame(pareto_design_evaluations), writeheader=false)



