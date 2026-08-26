using LinearAlgebra
using CSV
using DataFrames
df = CSV.File("earthsunmoonposdata.csv") |> DataFrame

# Getting the positions into a list
earthposraw = df[!, "Earth Position"]
earthpos = [
    parse.(Float64, split(strip(position, ['{', '}']), ",")) 
    for position in earthposraw
]
sunposraw = df[!, "Sun Position"]
sunpos = [
    parse.(Float64, split(strip(position, ['{', '}']), ",")) 
    for position in sunposraw
]
moonposraw = df[!, "Moon Position"]
moonpos = [
    parse.(Float64, split(strip(position, ['{', '}']), ",")) 
    for position in moonposraw
]

# Get the distances between bodies
earthmoondists=Float64[]
sunmoondists=Float64[]
earthsundists=Float64[]
thetas=Float64[]
for i in 1:length(earthpos)
    earthmoondistsi=sqrt((earthpos[i][1]-moonpos[i][1])^2+(earthpos[i][2]-moonpos[i][2])^2+(earthpos[i][3]-moonpos[i][3])^2)
    push!(earthmoondists,earthmoondistsi)
    sunmoondistsi=sqrt((sunpos[i][1]-moonpos[i][1])^2+(sunpos[i][2]-moonpos[i][2])^2+(sunpos[i][3]-moonpos[i][3])^2)
    push!(sunmoondists,sunmoondistsi)
    earthsundistsi=sqrt((earthpos[i][1]-sunpos[i][1])^2+(earthpos[i][2]-sunpos[i][2])^2+(earthpos[i][3]-sunpos[i][3])^2)
    push!(earthsundists,earthsundistsi)
    # Calculate angle theta between sunmoon and earthmoon
    thetai=acos((sunmoondistsi^2+earthmoondistsi^2-earthsundistsi^2)/(2*sunmoondistsi*earthmoondistsi)) # Law of Cosines
    push!(thetas,thetai)
end

earthposrel=Vector{Vector{Float64}}()
sunposrel=Vector{Vector{Float64}}()
# Scaling the distances (set average distance so that radius of the moon can be 1)
for t in 1:length(earthpos)
    earthposrelt=(earthpos[t][1:2] - moonpos[t][1:2])/(1737400)
    push!(earthposrel, earthposrelt)
    sunposrelt=(sunpos[t][1:2] - moonpos[t][1:2])/(1737400)
    push!(sunposrel, sunposrelt)
end
