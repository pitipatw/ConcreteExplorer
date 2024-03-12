using CSV
using DataFrames
using AsapToolkit
"""
function get_C, input area A and out put C from a csv file. 
csv file dataformat is A,C
"""
function get_C(A::Float64; test::Bool = false)
    if test
    #read the csv file
    filename = joinpath(@__DIR__,"table_AtoC.csv") #Make sure you got csv in the name
    data = CSV.read(filename,DataFrame)
    #get the column of A
    A_data = data[:,1]
    #get the column of C
    C_data = data[:,2]
    #get the index of the closest value of A
    index = findmin(abs.(A_data .- A))[2]
    #get the C value by interpolate between the two closest values of A
    Ai = A_data[index]
    Ai_1 = A_data[index+1] #the next value of A
    C = C_data[index] + (C_data[index+1] - C_data[index])/(Ai_1 - Ai)*(A - Ai)
    return C
    else
        println("to do, linked to AsapSection function")
        return nothing
    end

end

""" 
    function get_C(pixelframeelement::PixelFrameElement, araa::Float64)
Get compression depth of the section, given a compression area.
"""
function get_C(pixelframeelement::PixelFrameElement, area::Float64)
    c = depth_from_area(pixelframeelement.compoundsection, area)
    return c
end 

"""
    function get_Icrack input depth C output Icrack from a csv file. 
csv file dataformat is C,Icrack
"""
function get_Icrack(C::Float64; test::Bool = false)
    if test
    #read the csv file
    filename = joinpath(@__DIR__,"table_CtoIcrack.csv") #Make sure you got csv in the name
    data = CSV.read(filename, DataFrame)
    #get the column of C
    C_data = data[:,1]
    #get the column of Icrack
    Icrack_data = data[:,2]
    #get the index of the closest value of C
    index = findmin(abs.(C_data .- C))[2]
    #get the Icrack value by interpolate between the two closest values of C
    Ci = C_data[index]
    Ci_1 = C_data[index+1] #the next calue of C
    #get the Icrack value
    Icrack = Icrack_data[index]+(Icrack_data[index+1] - Icrack_data[index])/(Ci_1 - Ci)*(C - Ci)
    return Icrack
    else 
        println("Link to AsapSections' function")
        return nothing
    end
end

function get_Icrack(pixelframeelement::PixelFrameElement, c_depth::Float64)
    compoundsection = pixelframeelement.compoundsection

    c_depth_global = compoundsection.ymax - c_depth #global coordinate
    new_sections = Vector{SolidSection}()
    for sub_s in compoundsection.solids
        sub_s_ymax = sub_s.ymax #global coordinate
        sub_s_ymin = sub_s.ymin 
        c_depth_local = sub_s_ymax - c_depth_global
        if c_depth_local > 0
            c_depth_local = clamp(sub_s_ymax - c_depth_global, 0, sub_s_ymax - sub_s_ymin)
            push!(new_sections, sutherland_hodgman(sub_s, c_depth_local, return_section = true))
        end
    end
    return new_sections.Ix
end 
