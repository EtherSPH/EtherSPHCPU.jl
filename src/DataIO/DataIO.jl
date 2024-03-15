#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com>
  @ date: 2024/03/10 22:44:55
  @ license: MIT
  @ description:
 =#

abstract type AbstractDataIO end;

const kWallTimeFormat = "yyyy/mm/dd__HH:MM:SS.SS";

include("./VTPIO.jl");

function assureDirPathExist(data_io::DataIOType where {DataIOType <: AbstractDataIO})::Nothing
    if !isdir(data_io.dir_path_)
        mkdir(data_io.dir_path_)
    end
    return nothing
end
