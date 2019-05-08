module AA222Final

#Stuff from composite analysis. functions then things
export laminate_analyzer, safety_factor
export MaterialProperties, MaterialStrengths, LaminateProperties

#stuff from optimization methods
export nelder_mead, generate_simplex, combined_penalty

#stuff from structure analysis
export CompositePlate, SandwichPlate,  ThreePointBending
export composite_plate, sandwich_plate, three_point_bending

include("CompositeAnalysis.jl")
include("OptimizationMethods.jl")
include("StructureAnalysis.jl")

end