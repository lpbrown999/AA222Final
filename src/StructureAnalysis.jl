abstract type CompositeStructure end
abstract type CoreStructure end


struct CompositePlate <: CompositeStructure
	laminate::LaminateProperties		
	h::Float64						#height
	w::Float64						#Width
	l::Float64						#length
	weight::Float64					#plate weight
end
function composite_plate(laminate::LaminateProperties, w::Float64, l::Float64)
	return CompositePlate(laminate,laminate.t,w,l,w*l*laminate.t*laminate.ρ)
end

struct BatteryCore <: CoreStructure
	E::Float64	#Modulus 
	G::Float64	#Shear modulus 

	h::Float64	#height 
	w::Float64	#width 
	l::Float64	#length 

	ρ::Float64	#density
	ρe::Float64 #energy density 
	weight::Float64 #weight 
	capacity::Float64 #energy capacity
end
function battery_core(battery::BatteryProperties,h::Float64,w::Float64,l::Float64)
	E = battery.E
	G = battery.G
	ρ = battery.ρ
	ρe = battery.ρe
	V = w*h*l
	weight = ρ*V 
	capacity = ρe*V
	return BatteryCore(E,G,h,w,l,ρ,ρe,weight,capacity)
end	

struct MESCBoxBeam <: CompositeStructure
	#Geometry of the MESC box beam
	h::Float64
	w::Float64
	l::Float64

	#Materials of the box beam
	top_flange::CompositePlate
	webs::CompositePlate
	bot_flange::CompositePlate
	core::BatteryCore

	#Properties of the box beam
	weight::Float64 
	capacity::Float64
end
function MESC_box_beam(top_flange::CompositePlate,webs::CompositePlate,bot_flange::CompositePlate,battery_material::BatteryProperties,h::Float64,w::Float64,l::Float64,wb::Float64)
	# battery_w = w - 2*webs.h

	#Assume battery is as tall as it can be but scales width wise. 
	#This makes the most sense structurally since we do not have to define a location.

	battery_w = wb
	battery_h = h - (top_flange.h + bot_flange.h)
	core = battery_core(battery_material,battery_h,battery_w,l)
	weight = top_flange.weight + 2*webs.weight + bot_flange.weight + core.weight
	capacity = core.capacity
	return MESCBoxBeam(h,w,l,top_flange,webs,bot_flange,core,weight,capacity)
end


##Different load cases. Each load case can be analyzed for a different structure type
abstract type LoadCase end

struct ThreePointBending <: LoadCase
	composite_structure::CompositeStructure
	P::Float64
	δ::Float64
	SF::Float64
end

struct CantileverBending <: LoadCase
	composite_structure::CompositeStructure
	P::Float64 
	δ::Float64
	SF::Float64
end
function cantilever_bending(box_beam::MESCBoxBeam,P::Float64)
	L = box_beam.l
	#MESC Box beam in cantiliver loading. P defined positive upwards
	#In YZ frame (webs rotated 90)
	E_top = box_beam.top_flange.laminate.Ex
	E_web = box_beam.webs.laminate.Ex
	E_bot = box_beam.bot_flange.laminate.Ex
	E_core = box_beam.core.E

	h_top = box_beam.top_flange.h
	h_web = box_beam.webs.w
	h_bot = box_beam.bot_flange.h
	h_core = box_beam.core.h

	w_top = box_beam.top_flange.w
	w_web = box_beam.webs.h 
	w_bot = box_beam.bot_flange.w
	w_core = box_beam.core.w

	A_top = h_top*w_top
	A_web = h_web*w_web
	A_bot = h_bot*w_bot 
	A_core = h_core*w_core 

	Z_bot = h_bot/2
	Z_core = h_bot+h_core/2
	Z_web = h_bot+h_web/2 
	Z_top = h_bot+h_web+h_top/2

	EA = E_top*A_top + E_web*A_web + E_bot*A_bot + E_core*A_core
	Zbar = (E_top*A_top*Z_top + E_web*A_web*Z_web + E_bot*A_bot*Z_bot + E_core*A_core*Z_core)/EA

	Iyy_top  = (1/12)*w_top*(h_top^3)   + A_top*(Z_top-Zbar)^2
	Iyy_web  = (1/12)*w_web*(h_web^3)   + A_web*(Z_web-Zbar)^2
	Iyy_bot  = (1/12)*w_bot*(h_bot^3)   + A_bot*(Z_bot-Zbar)^2
	Iyy_core = (1/12)*w_core*(h_core^3) + A_core*(Z_core-Zbar)^2

	EIyy = E_top*Iyy_top + 2*E_web*Iyy_web + E_bot*Iyy_bot + E_core*Iyy_core

	#Compute defleciton. Ignoring shear deformation for now
	My = P*L
	Vz = P

	δ = (P*L^3)/(48*EIyy)

	##Compute safety factors
	##ASSUMPTIONS FOR NOW: 
	##Webs carry all of shear stress, uniformly distributed
	##Battery will not fail 
	##Critical locations:
	#1 - Top flange (bending stress)
	#2 - Bot flagne (bending stress)
	#3 - Highest point on web (bending + shear)
	#4 - Lowest point on web (bending + shear)
	#5 - Middle of web (peak shear, ignore for now since do not have expression)
	τ_ave = Vz/(2*A_web)
	sxx_top     = -E_top*(Z_top-Zbar)*My/EIyy
	sxx_bot     = -E_top*(Z_bot-Zbar)*My/EIyy
	sxx_web_top = -E_top*((Z_web+h_web/2)-Zbar)*My/EIyy
	sxx_web_bot = -E_top*((Z_web-h_web/2)-Zbar)*My/EIyy

	#Laminate level stress resultants
	#N = [Nxx,Nyy,Nxy]
	N_loc1 = [sxx_top*h_top,0.,0.]
	N_loc2 = [sxx_bot*h_bot,0.,0.]
	N_loc3 = [sxx_web_top*w_web,0.,τ_ave*w_web]
	N_loc4 = [sxx_web_bot*w_web,0.,τ_ave*w_web]

	SF_1   = safety_factor(box_beam.top_flange.laminate, N_loc1, [0.,0.,0.])
	SF_2   = safety_factor(box_beam.bot_flange.laminate, N_loc2, [0.,0.,0.])
	SF_3   = safety_factor(box_beam.webs.laminate, 		N_loc3, [0.,0.,0.])
	SF_4   = safety_factor(box_beam.webs.laminate, 		N_loc4, [0.,0.,0.])

	SF = minimum([SF_1,SF_2,SF_3,SF_4])
	return CantileverBending(box_beam,P,δ,SF)
end
