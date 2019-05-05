using MFCOpt
using Distributions

#Global Stuff
E1 = 11.2e6;					
E2 = 10.9e6
ν12 = .06
G12 = .94e6
tply = .01	
F1t = 140e3	
F2t = 124e3	
F1c = 130e3
F2c = 20000e3
F12 = 10.3e3
composite_properties = MaterialProperties(E1,E2,ν12,G12,tply)
composite_strengths = MaterialStrengths(F1t,F2t,F1c,F2c,F12)
N = [-1000,0,0]
M = [0,0,0]

allowable_angles=[-90., -45., -30., 0., 30., 45., 90.]


function objective(layup::Vector)
	return -sheet_SF(layup, composite_properties, composite_strengths, N, M)
end

function penalty(layup::Vector)
	return ply_orientation_penalty(layup, allowable_angles)
end

function unconstrained_objective(layup::Vector, t::Float64)
	return objective(layup::Vector) + t*penalty(layup::Vector)
end

#Number of plies in the layup
# nplies = 3
# t = .0001
# tgrow = 2

# dist = Uniform(-180, 180)
# S = [rand(dist,nplies) for i in 1:nplies+1]
# while t <= 1000
# 	global S
# 	global t
# 	global layup
# 	f(layup::Vector) = unconstrained_objective(layup, t)
# 	S, layup = nelder_mead(f, S)
# 	t *= tgrow
# end

layups = []
for ntries in 1:100
	global layups
	S = [rand(dist,nplies) for i in 1:nplies+1]
	f(layup::Vector) = unconstrained_objective(layup, 1.)
	S, layup = nelder_mead(f, S)
	push!(layups,layup)
end