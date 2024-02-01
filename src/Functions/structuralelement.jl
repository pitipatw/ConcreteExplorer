#Structural Element
abstract type StructuralElement end 
abstract type ReinforcedConcreteElement <: StructuralElement  end 
abstract type PixelFrameElement <: StructuralElement end 

#Structural Section
abstract type StructuralSection end 
abstract type ReinforcedConcreteSection <: StructuralSection end 
abstract type PixelFrameSection <: StructuralSection end 

mutable struct PixelFrameSection
    section

end 

function 


mutable struct PixelFrameBeam 