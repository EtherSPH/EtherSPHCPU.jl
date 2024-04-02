#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com>
  @ date: 2024/03/10 22:44:55
  @ license: MIT
  @ description:
 =#

abstract type AbstractDataIO end

const kWallTimeFormat = "yyyy/mm/dd__HH:MM:SS.SS"
const kPositionStrings = ["X", "Y", "Z"]
const kVelocityString = "Velocity"
const kVelocityStrings = ["VelocityX", "VelocityY", "VelocityZ"]
const kTypeString = "Type"
const kMassString = "Mass"
const kDensityString = "Density"
const kStringToTypeDict = Dict(
    "Float32" => Float32,
    "Float64" => Float64,
    "UInt8" => UInt8,
    "UInt16" => UInt16,
    "UInt32" => UInt32,
    "UInt64" => UInt64,
    "Int8" => Int8,
    "Int16" => Int16,
    "Int32" => Int32,
    "Int64" => Int64,
    "Int128" => Int128,
    "String" => String,
)

include("./VTPIO.jl")
include("./VTPToDataFrame.jl")

function assureDirPathExist(data_io::DataIOType where {DataIOType <: AbstractDataIO})::Nothing
    if !isdir(data_io.dir_path_)
        mkdir(data_io.dir_path_)
    end
    return nothing
end
