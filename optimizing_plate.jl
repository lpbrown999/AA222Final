# using AA222Final

# function objective(layup::Vector)
# 	#Material Properties - Assuming Carbon (AGP370) / Epoxy (3501-6S)
# 	#All units lbs, in
# 	E1 = 11.2e6;		
# 	E2 = 10.94e6;      
# 	ν12 = .06			
# 	G12 = .94e6			
# 	tply = .01
# 	ρ = .056		#lb/in^3			

# 	#Failure Properties of composite
# 	F1t = 140e3			
# 	F2t = 124e3	
# 	F1c = 130e3
# 	F2c = 20000e3
# 	F12 = 10.3e3
# 	composite_properties = MaterialProperties(E1,E2,ν12,G12,tply,ρ)
# 	composite_strengths = MaterialStrengths(F1t,F2t,F1c,F2c,F12)

# 	#Loading
# 	N = [-1000,0,0]
# 	M = [0,0,0]

# 	return -sheet_SF(layup, composite_properties, composite_strengths, N, M)
# end

# function penalty(layup::Vector)
# 	allowable_angles=[-90., -45., -30., 0., 30., 45., 90.]

# 	return ply_orientation_penalty(layup, allowable_angles)
# end

# function unconstrained_objective(layup::Vector, t::Float64)
# 	return objective(layup::Vector) + t*penalty(layup::Vector)
# end


# function optimize()
# 	#Number of plies in the layup
# 	nplies = 3
# 	t = 1
# 	tgrowth = 2
# 	n = 100
# 	num_iteraitons = 10

# 	dist = Uniform(-180, 180)
# 	S = [rand(dist,nplies) for i in 1:nplies+1]
# 	for i in 1:num_iteraitons
# 		f(layup) = unconstrained_objective(layup,t)
# 		S, layup = nelder_mead(f,S)
# 		t*=tgrow
# 	end

# end


