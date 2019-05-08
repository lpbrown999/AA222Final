using AA222Final

#Runs the sandwiwch bending problem
function sandwich_bending_case(inpvec)
	##This input
	#Fist n-1 entries are the facesheet ply orientations.
	#last entry is the thickness of the core!!!
	face_sheet_layup = inpvec[1:end-1]
	tcore  = max(1e-10,inpvec[end])

	#Composite 
	composite_properties = MaterialProperties(11.2e6,0.94e6,.06,.94e6,.01,.056)
	composite_strengths = MaterialStrengths(140e3,124e3,130e3,130e3,10.3e3)
	face_sheet_laminate = laminate_analyzer(face_sheet_layup,composite_properties,composite_strengths)

	#Core
	core_properties = MaterialProperties(1e6,1e6,.3,1e6/(2*1.3),tcore,.01)
	core_strengths = MaterialStrengths(Inf,Inf,Inf,Inf,Inf)
	core_layup = [0.]
	core_laminate = laminate_analyzer(core_layup,core_properties,core_strengths)

	#Geometry
	w = 10.
	l = 20.
	P = 60.

	sandwich = sandwich_plate(face_sheet_laminate,core_laminate,face_sheet_laminate,w,l)
	bending_sandwich = three_point_bending(sandwich,P)
	three_point_bending_case = three_point_bending(sandwich,P)
	return three_point_bending_case
end

#Gets the weight
function sandwich_weight(inpvec)
	three_point_bending_case = sandwich_bending_case(inpvec)
	comp_struct = three_point_bending_case.composite_structure
	weight = comp_struct.top_sheet.weight + comp_struct.core.weight + comp_struct.bot_sheet.weight
	return weight
end

#Gets the deflection
function sandwich_defl(inpvec)
	three_point_bending_case = sandwich_bending_case(inpvec)
	return three_point_bending_case.δ
end

#Deflection constraint
function defl_constraint(inpvec,δmax::Float64)
	# δ <= δmax, c<=0
	δ = sandwich_defl(inpvec)
	return δ - δmax
end

function optimize(f,c,x0; n=100, a=-1,b=1, p1=1,p1g=2,p2=1,p2g=2)
	#f -> function to optimize
	#c -> constraint function
	#x0 -> initial design point
	#n -> number of evals
	f_penalized(x) = f(x) + combined_penalty(c(x),p1,p2)
	ndim = length(x0)
	S = generate_simplex(ndim,a,b)
	# S[1] = x0
	_, xbest = nelder_mead(f_penalized, S, n)
	return xbest
end

#Constraint function
function c(inpvec)
	cret = []
	push!(cret, defl_constraint(inpvec,.5))
	return cret 
end

ndim = 5
x0 = [0.,0.,0.,0.]
a = [-180*ones(ndim-1);0.]		
b = [180*ones(ndim-1); 5.]
best_layup = optimize(sandwich_weight, c, x0, n=10000, a=a, b=b, p1=1000, p2=1000)
# Now want to optimize ndim to provide lowest weight 
