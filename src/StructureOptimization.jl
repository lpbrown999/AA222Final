
function sheet_SF(layup::Vector, material_properties::MaterialProperties, material_strengths::MaterialStrengths, N::Vector, M::Vector)

	#Composite properties
	plate_laminate = laminate_analyzer(layup, material_properties)
	ϵo, K = load_deformation_solver(plate_laminate, N, M, "load")
	σ12_bot, σ12_top, _, _ = ply_local_stress_strain(plate_laminate, ϵo, K)
	SFbot = safety_factor(σ12_bot, material_strengths)
	SFtop = safety_factor(σ12_top, material_strengths)		
	SF = minimum([SFbot;SFtop])
	return SF
end

function ply_orientation_penalty(layup::Vector, allowable_angles::Vector)
	min_deviation = zeros(length(layup))
	penalty_term = 0
	for i in 1:length(layup)
		min_deviation[i] = minimum(abs.(layup[i].-allowable_angles))
	end
	penalty_term = quad_penalty_equality(min_deviation)
	return penalty_term
end