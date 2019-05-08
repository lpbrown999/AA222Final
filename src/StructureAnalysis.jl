
abstract type CompositeStructure end

struct CompositePlate <: CompositeStructure
	laminate::LaminateProperties		
	w::Float64						#Width
	l::Float64						#length
	h::Float64						#height
	weight::Float64					#plate weight
end
function composite_plate(laminate::LaminateProperties, w::Float64, l::Float64)
	return CompositePlate(laminate,w,l,laminate.t,w*l*laminate.t*laminate.ρ)
end


struct SandwichPlate <: CompositeStructure
	top_sheet::CompositePlate
	core::CompositePlate
	bot_sheet::CompositePlate
end
function sandwich_plate(top_laminate::LaminateProperties, core_laminate::LaminateProperties, bot_laminate::LaminateProperties, w::Float64, l::Float64)
	top_sheet = composite_plate(top_laminate,w,l)
	core = composite_plate(core_laminate,w,l)
	bot_sheet = composite_plate(bot_laminate,w,l)
	return SandwichPlate(top_sheet,core,bot_sheet)
end


##Different load cases. Each load case can be analyzed for a different structure type!
abstract type LoadCase end

struct ThreePointBending <: LoadCase
	composite_structure::CompositeStructure
	P::Float64
	δ::Float64
	SF::Float64
end

#Three point bending for a plate
function three_point_bending(composite_structure::CompositePlate,P::Float64)
	E = composite_structure.laminate.Ex
	L = composite_structure.l
	h = composite_structure.h
	w = composite_structure.w
	Iyy = (1/12)*(w*h^3)
	δ = (P*L^3)/(48*E*Iyy)

	N = [0.,0.,0.]
	M = [(P*L/4)/w, 0., 0.]
	SF = safety_factor(composite_structure.laminate, N, M)
	return ThreePointBending(composite_structure,P,δ,SF)
end

#Three point bending for a sandwich plate
function three_point_bending(composite_structure::SandwichPlate,P::Float64)
	#Grab all of the component properties
	E_top  = composite_structure.top_sheet.laminate.Ex
	E_core = composite_structure.core.laminate.Ex
	E_bot = composite_structure.bot_sheet.laminate.Ex

	w_top = composite_structure.top_sheet.w
	w_core = composite_structure.core.w
	w_bot = composite_structure.bot_sheet.w

	h_top = composite_structure.top_sheet.h
	h_core = composite_structure.core.h
	h_bot = composite_structure.bot_sheet.h

	A_top = h_top*w_top
	A_core = h_core*w_core 
	A_bot = h_bot*w_bot

	Z_bot = h_bot/2
	Z_core = h_bot + h_core/2
	Z_top = h_bot + h_core + h_top/2

	L = composite_structure.top_sheet.l 

	#Compute modulus weighted area, centroid 
	EA = E_top*A_top + E_core*A_core + E_bot*A_bot
	Zbar = (Z_top*E_top*A_top + Z_core*E_core*A_core + Z_bot*E_bot*A_bot)/EA

	#Compute bending stifness
	Iyy_top  = (1/12)*w_top*(h_top^3)  + A_top*(Z_top-Zbar)^2
	Iyy_core = (1/12)*w_core*(h_core^3) + A_core*(Z_core-Zbar)^2
	Iyy_bot  = (1/12)*w_bot*(h_bot^3) + A_bot*(Z_bot-Zbar)^2
	EIyy = E_top*Iyy_top + E_core*Iyy_core + E_bot*Iyy_bot

	#Deflection
	δ = (P*L^3)/(48*EIyy)

	#Saftey Factor
	#Critical stresses -> Since doubly symmetric reduces to this.
	My = -P*L/4
	sxx_top = E_top*(Z_top-Zbar)*My/EIyy
	sxx_bot = E_bot*(Z_bot-Zbar)*My/EIyy
	sxx_core1 = E_core*((Z_top-h_top/2) - Zbar)*My/EIyy 
	sxx_core2 = E_core*((Z_bot+h_bot/2) - Zbar)*My/EIyy

	M = [0.,0.,0.]
	N_top = [sxx_top*h_top,0.,0.]
	N_bot = [sxx_bot*h_bot,0.,0.]
	N_core1 = [sxx_core1*h_core,0.,0.]
	N_core2 = [sxx_core2*h_core,0.,0.]

	SF_top   = safety_factor(composite_structure.top_sheet.laminate, N_top, M)
	SF_bot   = safety_factor(composite_structure.bot_sheet.laminate, N_bot, M)
	SF_core1 = safety_factor(composite_structure.core.laminate, N_core1, M)
	SF_core2 = safety_factor(composite_structure.core.laminate, N_core2, M)

	SF = minimum([SF_top,SF_bot,SF_core1,SF_core2])
	return ThreePointBending(composite_structure,P,δ,SF)
end