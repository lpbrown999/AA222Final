module AA222Final

#Stuff from composite analysis. functions then things
export MaterialProperties, MaterialStrengths, LaminateProperties, BatteryProperties
export laminate_analyzer, safety_factor

#stuff from optimization methods
export nelder_mead, generate_simplex, combined_penalty

#stuff from structure analysis
export CompositePlate, MESCBoxBeam, CantileverBending
export composite_plate, MESC_box_beam, cantilever_bending


include("CompositeAnalysis.jl")
include("OptimizationMethods.jl")
include("StructureAnalysis.jl")

end