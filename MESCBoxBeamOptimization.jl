using AA222Final
using CSV, DataFrames
# using PyPlot

# all units in lb - in

## Notes
# New method

function mesc_box_beam_objective(input_vector::Vector; P, num_dummy_plies::Int64, constraint_vec::Vector, allowable_ply_angles::Vector, p1::Float64, p2::Float64, mode::String)
	#Input_vector: now merged with num_ply_vec
	n_top_flange = Int( floor(max(min(input_vector[1],num_dummy_plies*2), 1.0)) ) 
	n_webs       = Int( floor(max(min(input_vector[2],num_dummy_plies*2), 1.0)) )
	n_bot_flange = Int( floor(max(min(input_vector[3],num_dummy_plies*2), 1.0)) )

	h = max(input_vector[4],.05)
	w = max(input_vector[5],.05)
	wb = max(input_vector[6],.05)

	layups = input_vector[7:end]
	top_layup_dummy_plies = layups[1:num_dummy_plies]
	web_layup_dummy_plies = layups[1+num_dummy_plies:2*num_dummy_plies]
	bot_layup_dummy_plies = layups[1+2*num_dummy_plies:3*num_dummy_plies]

	#Construct symmetric layup from the dummy plies
	top_layup = construct_symm_layup(top_layup_dummy_plies, n_top_flange)
	web_layup = construct_symm_layup(web_layup_dummy_plies, n_webs)
	bot_layup = construct_symm_layup(bot_layup_dummy_plies, n_bot_flange)
	layups = [top_layup; web_layup; bot_layup]
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

	## p1, p2 are penalty factors

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

	#Objective + penalties on the constraints
	#no count penalty on the layup orientation since we will never be exact
	if mode == "optimize"
		return weight + combined_penalty(cvec1.*cvec1_weights,p1,p2) + combined_penalty(cvec2,p1/10,0)
	elseif mode == "evaluate"
		return [weight,capacity,SF,δ,h,w,wb]		#For plotting and to check if ok!
	elseif mode == "check_constraints"
		return cvec1,cvec2
	elseif mode == "debug"
		return bending_case
	else
		return Nothing
	end
end

function construct_symm_layup(dummy_plies::Vector,nplies::Int64)
	if nplies == 1
		return [dummy_plies[1]]
	elseif isodd(nplies)
		front_half = dummy_plies[1:Int((nplies+1)/2)]
		back_half = reverse(front_half[1:end-1])
		return [front_half;back_half]
	elseif iseven(nplies)
		front_half = dummy_plies[1:Int(nplies/2)]
		return [front_half;reverse(front_half)]
	end
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

function weight_capacity_tradeoff(P,Capacitymins,neval,num_growths)
	
	#Set up constraints
	hmax = 2.0
	hmin = 1.0
	wmax = 1.0
	wmin = 0.5
	SFmin = 1.1
	SFmax = Inf
	δmax = .005
	allowable_ply_angles = [-45.,0.,45.,90.]

	#Set up simplex generation
	num_dummy_plies = 10
	ndim = 6 + num_dummy_plies*3
	a = [zeros(3);                  hmin; wmin; wmin; -90. *ones(num_dummy_plies*3)] #LB on h,w,wb, ply angles for simplex generation
	b = [2*num_dummy_plies*ones(3); hmax; wmax; wmax; 135. *ones(num_dummy_plies*3)] #UB
	ndim = length(a)

	ys = []
	
	for Capacitymin in Capacitymins
		println("P: $P, Capacitymins: $Capacitymin")
		constraint_vec = [hmax,hmin,wmax,wmin,SFmin,SFmax,δmax,Capacitymin]
		
		function_to_optimize(input_vector,p1,p2) = mesc_box_beam_objective(input_vector,P=P,p1=p1,p2=p2,num_dummy_plies=num_dummy_plies,constraint_vec=constraint_vec,allowable_ply_angles=allowable_ply_angles, mode="optimize")
		function_evaluation(input_vector,p1,p2)  = mesc_box_beam_objective(input_vector,P=P,p1=p1,p2=p2,num_dummy_plies=num_dummy_plies,constraint_vec=constraint_vec,allowable_ply_angles=allowable_ply_angles, mode="evaluate")
		
		@time x = optimize(function_to_optimize, ndim; neval=neval, a=a, b=b, p1=2.0, p1g=1.5, p2=2.0, p2g=1.5, num_growths=num_growths)
		y = function_evaluation(x,0.,0.)
		push!(ys, y)	
	end

	return ys
end

Ps = [100,200,300]
for P in Ps
	neval = 20001
	num_growths = 20
	Capacitymins = repeat(0:2:90,100)	   								#capacities to iterate over
	ys = weight_capacity_tradeoff(P, Capacitymins, neval, num_growths)
	csv_name = "results_may28/test_result_$(P)_$(neval)_$(num_growths).csv"
	CSV.write(csv_name,  DataFrame(hcat(ys...)'), header=["Weight","Capacity","SafetyFactor","Deflection","Height","Width","BatteryWidth"], transpose=true)
end