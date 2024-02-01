function pixelframe_properties!(pixelframesection::PixelFrameSection)::PixelFrameSection

    pixelframesection.fc′::Float64 = 0.0
    pixelframesection.fR1::Float64 = 0.0
    pixelframesection.fR3::Float64 = 0.0

    pixelframesection.pt_area::Vector{Float64} = [0.0]
    pixelframesection.pt_force::Vector{Float64} = [0.0]
    pixelframesection.pt_pos::Matrix{Float64} = [0.0 0.0]

    #Material Properties
    pixelframesection.Ec::Float64  = 0.0
    pixelframesection.Eps::Float64 = 0.0


    return pixelframesection
end

