using AA222Final
using PyPlot

# all units in lb - in

#Main objective + constraint function
function mesc_box_beam_objective(input_vector::Vector; num_ply_vec::Vector, constraint_vec::Vector, allowable_ply_anlges::Vector, p1::Float64, p2::Float64, mode::String)
	
	##################### DESCRIPTION / extract inputs #####################
	#num_ply_vec is a vector of the number of plies in the top flange, web, and bottom flange
	n_top_flange = num_ply_vec[1]
	n_webs       = num_ply_vec[2]
	n_bot_flange = num_ply_vec[3]
	
	## input_vector is a vector of beam geometry and layups, which is what we are optimizing over
		#Element 1: beam total height (enforce feasibility)
		#Element 2: beam total width  (enforce feasibility)
		#Element 3: battery width (enforce feasibility)
		#Elements 4 to end: layups 
	h = max(input_vector[1],.1)
	w = max(input_vector[2],.1)
	wb = max(input_vector[3],.1)

	layups = input_vector[4:end]
	top_layup = layups[1:n_top_flange]
	web_layup = layups[1+n_top_flange:n_top_flange+n_webs]
	bot_layup = layups[1+n_top_flange+n_webs:n_top_flange+n_webs+n_bot_flange]

	## constraint_vec is a vector of constraints
		#Element1 : maximum height
		#Element2 : maximum width
		#Element3 : minimum height
		#Element4 : minimum width
		#Element5 : minimum safety factor
		#Element6 : maximum defleciton
		#Element7 : minCapacity, the minimum level of the the negative capacity since we are doing pareto curve, maximizing capactiy
	hmax = constraint_vec[1]
	wmax = constraint_vec[2]
	hmin = constraint_vec[3]
	wmin = constraint_vec[4]
	SFmin = constraint_vec[5]
	δmax = constraint_vec[6]
	Capacitymin = constraint_vec[7]
	
	## p1, p2 are penalty factors

	## Mode is what we are doing really.
	# options: "optimize" returns weight + constraint penalties 
	#		   "evaluate" returns height, width, weight, capacity, SF


	##################### Start main function #######################
	######Fixed Stuff#########
	l = 30.					#Fixed length of 30 inches, can change
	P = 0.1				    #Fixed tip load

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
	box_beam = MESC_box_beam(top_plate,web_plate,bot_plate,battery_properties,h,w,l,wb)
	bending_case = cantilever_bending(box_beam,P)

	##################### Value return #####################
	# all constrants of form <= 0
	# constraint_vec1 -> non layup constraints
	# constraint_vec2 -> layup orientation constraints
	SF = bending_case.SF
	δ = bending_case.δ
	capacity = box_beam.capacity
	weight   = box_beam.weight
	
	cvec1 = [h-hmax, 						#h <= hmax 
			 w-wmax, 						#w <= wmax 
			 wb-w+2*box_beam.webs.h,		#wbattery <= w - 2*webthickness	#Battery must fit inside

			 hmin-h, 						#h >= hmin
			 wmin-w,						#w >= wmin  

			 SFmin-bending_case.SF, 		#SF >= SFmin 
			 δ-δmax, 						#δ <= δmax
			 Capacitymin-capacity]			#Capacity >= Capacitymmin


	cvec2 = []
	for layer in layups	
		push!(cvec2,minimum(abs.(layer.-allowable_ply_anlges)))		# c evaluates to the minimum distance from an allowed angle
	end

	#Objective + penalties on the constraints
	#no count penalty on the layup orientation since we will never be exact
	
	if mode == "optimize"
		return weight + combined_penalty(cvec1,p1,p2) + combined_penalty(cvec2,p1,0)

	elseif mode == "evaluate"
		return [h,w,weight,capacity,SF,δ]
	
	else
		return Nothing
	end
end

function optimize(f, ndim; neval=100, a=-1, b=1, p1=1., p1g=2., p2=1., p2g=2., num_growths=1)
	## Inputs
	#f -> function to optimize. Should include constraints already except for p1, p2. 
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

# function main()

#Constraints 
hmax = 2.0
wmax = 1.0
hmin = 1.0
wmin = 0.5
SFmin = 1.1
δmax = 1.0
Capacitymin = 10.
constraint_vec = [hmax,wmax,hmin,wmin,SFmin,δmax,Capacitymin]
allowable_ply_angles = [-45.,0.,45.,90.]

num_ply_vec =[4,3,2]       #Number of plys in top flange, webs, bot flange
a = [zeros(3); -180. *ones(sum(num_ply_vec))]	#LB on h,w,wb, ply angles for simplex generation
b = [5*ones(3); 180. *ones(sum(num_ply_vec))]	#UB
ndim = length(a)

f_optimize(input_vector,p1,p2) = mesc_box_beam_objective(input_vector, p1=p1, p2=p2, num_ply_vec=num_ply_vec, constraint_vec=constraint_vec, allowable_ply_anlges=allowable_ply_angles, mode="optimize")
f_evaluate(input_vector)       = mesc_box_beam_objective(input_vector, p1=0., p2=0., num_ply_vec=num_ply_vec, constraint_vec=constraint_vec, allowable_ply_anlges=allowable_ply_angles, mode="evaluate")

design = optimize(f_optimize, ndim, neval=100000, a=a, b=b, p1=1.0, p1g=1.3,  p2=10., p2g=1.3, num_growths=100)
design_eval = f_evaluate(design)




