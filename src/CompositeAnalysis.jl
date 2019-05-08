## Classes
abstract type CompositeInformation end

struct MaterialProperties <: CompositeInformation
	E1::Float64
	E2::Float64 
	ν12::Float64
	G12::Float64
	tply::Float64
	ρ::Float64			#Density
end

struct MaterialStrengths <: CompositeInformation
	F1t::Float64
	F2t::Float64
	F1c::Float64
	F2c::Float64
	F12::Float64
end

struct LaminateProperties <: CompositeInformation
	layup::Vector 		 
	Qbar::Array
	A::Array
	B::Array
	D::Array
	Ex::Float64
	Ey::Float64
	Gxy::Float64
	z::Vector
	t::Float64
	ρ::Float64     			#Density 
	strengths               #This can be either a Vector of material strengths for each ply 
							#or a single material strength that is for all plies
end

function laminate_analyzer(layup::Vector, material::MaterialProperties,strengths::MaterialStrengths)
	## For only one material throughout
	# layup:  bottom up of layup orientations
	# material: material properties stored in a MaterialProperties struct
	# OUTPUTS
	# laminate

	n_layup = length(layup)
	tply = material.tply
	ρ    = material.ρ

	## Q, Qbar for each ply
	Q = zeros(3,3,n_layup)
	Qbar = zeros(3,3,n_layup)
	for i in 1:n_layup
		Q[:,:,i] = Q_calc(material)
		Qbar[:,:,i] = Qbar_calc(Q[:,:,i], layup[i])
	end

	## z vector based on tply -> accomodates non uniform ply thickness
	#	According to Eng Mechanics of Composite Materials 2ed (Daniel, Ishai) pg 
	#	159, z is always defined from mid thickness. not modulus weighted.
	t_tot = tply*n_layup
	z = zeros(n_layup+1)
	z[1] = -t_tot/2
	for i in 1:n_layup
	    z[i+1] = z[i]+tply
	end

	## ABD, Equivalent properties
	A,B,D = ABD_calc(Qbar, z)
	Ex,Ey,Gxy = effective_moduli(A, t_tot)

	return LaminateProperties(layup,Qbar,A,B,D,Ex,Ey,Gxy,z,t_tot,ρ,strengths)
end

function laminate_analyzer(layup::Vector, materials::Vector, strengths::Vector)
	## For different materials throughout
	## Pass in a vector of MaterialProperties
	# layup:  bottom up of layup orientations
	# materials: vector material properties stored in a MaterialProperties struct
	# OUTPUTS
	# ABD matrixes, Qbar, z, Exx, Gxy equivalent for the material
	n_layup = length(layup)
	if n_layup != length(materials)
		println("Wrong number of materials for layup.")
		return
	end

	## Q, Qbar for each ply
	## As well as total thickness
	## and density calc
	Q = zeros(3,3,n_layup)
	Qbar = zeros(3,3,n_layup)
	t_tot = 0
	ρt = zeros(n_layup)

	for i in 1:n_layup
		Q[:,:,i] = Q_calc(materials[i])
		Qbar[:,:,i] = Qbar_calc(Q[:,:,i], layup[i])
		t_tot += materials[i].tply
		ρt[i] = materials[i].tply * materials[i].ρ
	end

	ρ = sum(ρt)/t_tot

	## z vector based on tply -> accomodates non uniform ply thickness
	#	According to Eng Mechanics of Composite Materials 2ed (Daniel, Ishai) pg 
	#	159, z is always defined from mid thickness. not modulus weighted.
	z = zeros(n_layup+1)
	z[1] = -t_tot/2
	for i in 1:n_layup
	    z[i+1] = z[i]+materials[i].tply
	end

	## ABD, Equivalent properties
	A,B,D = ABD_calc(Qbar, z)
	Ex,Ey,Gxy = effective_moduli(A, t_tot)

	return LaminateProperties(layup,Qbar,A,B,D,Ex,Ey,Gxy,z,t_tot,ρ,strengths)
end


## Helper functions
function T_eps(theta::Float64)
	# Calculates the engineering strain transformation matrix T_epsilon
	# theta should be given in degrees not radians
	m = cosd(theta);
	n = sind(theta);

	T_e = [m^2    n^2     m*n;
		   n^2    m^2     -m*n;
		   -2*m*n 2*m*n   (m^2 - n^2)];   
	return T_e
end

function T_sig(theta::Float64)
	# Calculates the engineering stress transformation matrix T_sigma
	# theta should be given in degrees not radians
	m = cosd(theta);
	n = sind(theta);

	T_s = [m^2    n^2     2*m*n;
		   n^2    m^2     -2*m*n;
		   -m*n   m*n   (m^2 - n^2)];   
	return T_s
end

## Property calculations
function Q_calc(material::MaterialProperties)
	#computes reduced stiffness Q of a single orthotropic layer 
	#refered to the principal material axes
	E1 = material.E1
	E2 = material.E2
	ν12 = material.ν12
	G12 = material.G12

	ν21 = E2*ν12/E1				#(4.49)
	Q11 = E1/(1-ν12*ν21)		#(4.56)
	Q22 = E2/(1-ν12*ν21)
	Q12 = ν21*E1/(1-ν12*ν21)
	Q66 = G12

	Q = [Q11 Q12 0;
		 Q12 Q22 0; 
		 0   0   Q66]
	return Q
end

function Qbar_calc(Q::Array,theta::Float64)
	#Compute reduced stiffness Qbar of a single layer
	# in global coordinates 
	T_e = T_eps(theta)
	Qbar = T_e'*Q*T_e
	return Qbar
end

function ABD_calc(Qbar::Array,z::Vector)
	# Computes Force / deformation coupling ABD
	# Qbars 3x3xn array of global reduced stiffnesses of each ply
	# z: top, bottom levels of each ply

	n = size(Qbar,3)
	A = zeros(3,3); B = zeros(3,3); D = zeros(3,3);
	for k = 1:n
		A +=       Qbar[:,:,k]*(z[k+1]   - z[k]  );
		B += (1/2)*Qbar[:,:,k]*(z[k+1]^2 - z[k]^2);
		D += (1/3)*Qbar[:,:,k]*(z[k+1]^3 - z[k]^3);
	end
	return A,B,D
end

function effective_moduli(A::Array,t_tot::Float64)
	#Only valid for symmetric layups
	# A -> Global stiffness 
	# t -> Laminate total thickness

	Ap = inv(A)
	Ex  = 1/(t_tot*Ap[1,1])
	Ey  = 1/(t_tot*Ap[2,2])
	Gxy = 1/(t_tot*Ap[3,3])
	return Ex,Ey,Gxy
end


## Load Deformation functions
function load_deformation_solver(laminate::LaminateProperties,norm_in::Vector,bend_in::Vector,mode::String)
	#INPUTS
	# laminate properties
	# inplane_input: either N (F/L) or mid plane straine (3x1)
	# bending_input: either M (F*L/L) or curvature K (3x1)
	# mode: string, either "load" (input N,M) or "strain" (input e_o, K)
	#OUTPUTS
	# norm_out: either N or e_o
	# bend_out: either M or K

	ABD = [laminate.A laminate.B; laminate.B laminate.D]
	if mode == "load"
		N = norm_in
		M = bend_in
		strains = inv(ABD)*[N;M]
		eps_o = strains[1:3]
		K = strains[4:6]

		norm_out = eps_o
		bend_out = K

	elseif mode == "strain"
		eps_o = norm_in 
		K = bend_in 
		forces = ABD*[eps_o;K]
		N = forces[1:3]
		M = forces[4:6]

		norm_out = N
		bend_out = M
	end
	return norm_out, bend_out
end

function ply_local_stress_strain(laminate::LaminateProperties,ϵo::Vector,K::Vector)
	#INPUTS
	# layup -> vector of ply angles
	# laminate propeties
	# eo -> mid plain strain components
	# K -> curvature components

	layup = laminate.layup
	z = laminate.z
	Qbar = laminate.Qbar

	#computes the local stress at top and botom of each ply
	n_layup = length(layup)
	ϵxy_bot = zeros(3,n_layup)
	ϵxy_top = zeros(3,n_layup)
	ϵ12_bot = zeros(3,n_layup)
	ϵ12_top = zeros(3,n_layup)
	σxy_bot = zeros(3,n_layup)
	σxy_top = zeros(3,n_layup)
	σ12_bot = zeros(3,n_layup)
	σ12_top = zeros(3,n_layup)

	T_e = zeros(3,3,n_layup)
	T_s = zeros(3,3,n_layup)

	#Rotation matrices
	for i in 1:n_layup
		T_e[:,:,i] = T_eps(layup[i])
		T_s[:,:,i] = T_sig(layup[i])
	end

	for i in 1:n_layup
		#Global strain, stress. (7.8) (4.31) (lecture pack 4)
		ϵxy_bot[:,i] = ϵo + z[i]*K
		ϵxy_top[:,i] = ϵo + z[i+1]*K

		σxy_bot[:,i] = Qbar[:,:,i]*ϵxy_bot[:,i]
		σxy_top[:,i] = Qbar[:,:,i]*ϵxy_top[:,i]

		#Local strain, stress
		ϵ12_bot[:,i] = T_e[:,:,i]*ϵxy_bot[:,i]
		ϵ12_top[:,i] = T_e[:,:,i]*ϵxy_top[:,i]

		σ12_bot[:,i] = T_s[:,:,i]*σxy_bot[:,i]
		σ12_top[:,i] = T_s[:,:,i]*σxy_top[:,i]
	end
	return σ12_bot, σ12_top, ϵ12_bot, ϵ12_top
end

## Safety factor stuff -> tsai wu
#Failure criteria - Tsai Wu
#Laminates with one material
function TsaiWuSF(σ12::Array,strengths::MaterialStrengths)
	# for laminates with all same material
	
	#Inputs
	#σ12:    3xn LOCAL stress state 
	#Rest: lamina strengths. Either single input or vector for dif matl.
	#	Returns: safety factor of each ply. <1 is failure. Multiply SF onto loading to get failure load
	#	Computes laminate and layer saftery factors (9.7,9.8,9.11)
	#	See Ch6 for description of lamina strengths
	
	F1t = strengths.F1t
	F2t = strengths.F2t
	F1c = strengths.F1c
	F2c = strengths.F2c
	F12 = strengths.F12
	f1  = 1/F1t - 1/F1c
	f11 = 1/(F1t*F1c)
	f2  = 1/F2t - 1/F2c
	f22 = 1/(F2t*F2c)
	f66 = 1/(F12^2)
	f12 = -.5*sqrt(f11*f22)

	n_layup = size(σ12,2)
	Sfa = zeros(n_layup)			#layer safety factor		

	for i in 1:n_layup
		a = f11*σ12[1,i]^2 + f22*σ12[2,i]^2 + f66*σ12[3,i]^2 + 2*f12*σ12[1,i]*σ12[2,i]
		b = f1*σ12[1,i] + f2*σ12[2,i]
		Sfa[i] = (-b+sqrt(b^2 + 4*a))/(2*a)
	end

	return Sfa
end

#Laminates with many materials -> supply vector of material strenths
function TsaiWuSF(σ12::Array,strengths::Vector)
	# for laminates with different materials
	#Inputs
	#σ12:    3xn LOCAL stress state 
	#Rest: lamina strengths. Either single input or vector for dif matl.
	#	Returns: safety factor of each ply. <1 is failure. Multiply SF onto loading to get failure load
	#	Computes laminate and layer saftery factors (9.7,9.8,9.11)
	#	See Ch6 for description of lamina strengths
	
	n_layup = size(σ12,2)
	Sfa = zeros(n_layup)			#layer safety factor		

	for i in 1:n_layup

		F1t = strengths[i].F1t
		F2t = strengths[i].F2t
		F1c = strengths[i].F1c
		F2c = strengths[i].F2c
		F12 = strengths[i].F12
		f1  = 1/F1t - 1/F1c
		f11 = 1/(F1t*F1c)
		f2  = 1/F2t - 1/F2c
		f22 = 1/(F2t*F2c)
		f66 = 1/(F12^2)
		f12 = -.5*sqrt(f11*f22)

		a = f11*σ12[1,i]^2 + f22*σ12[2,i]^2 + f66*σ12[3,i]^2 + 2*f12*σ12[1,i]*σ12[2,i]
		b = f1*σ12[1,i] + f2*σ12[2,i]
		Sfa[i] = (-b+sqrt(b^2 + 4*a))/(2*a)
	end
	
	return Sfa
end

#Wrapper for the otehr 
function safety_factor(laminate::LaminateProperties, N::Vector, M::Vector)
	ϵo, K = load_deformation_solver(laminate,N,M,"load")
	σ12b, σ12t, _, _ = ply_local_stress_strain(laminate,ϵo,K)
	SF1 = TsaiWuSF(σ12b, laminate.strengths)		#bottom of each ply
	SF2 = TsaiWuSF(σ12t, laminate.strengths)		#laminate top of each ply
	SF = minimum([SF1;SF2])
	return SF
end

