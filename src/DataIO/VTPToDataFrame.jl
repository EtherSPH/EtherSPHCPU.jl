#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com>
  @ date: 2024/03/26 18:14:41
  @ license: MIT
  @ description:
 =#

function getPointsDataFrameFromVTPFile(vtp_file::VTKFile)::Tuple{GroupedDataFrame, Int}
    points = get_points(vtp_file)
    velocitys = get_point_data(vtp_file)[kVelocityString] |> get_data
    dim = size(velocitys, 1)
    dataframe = DataFrame(points[1:dim, :]', Symbol.(kPositionStrings)[1:dim])
    for key in get_point_data(vtp_file) |> keys
        if key != kVelocityString
            dataframe[!, Symbol(key)] .= get_point_data(vtp_file)[key] |> get_data
        else
            for i_dim in 1:dim
                dataframe[!, Symbol(kVelocityStrings[i_dim])] .= velocitys[i_dim, :]
            end
        end
    end
    return (groupby(dataframe, :Type), dim)
end

function getFieldDataFrameFromVTPFile(vtp_file::VTKFile)::DataFrame
    field_data = LightXML.root(vtp_file.xml_file)[vtp_file.file_type][1]["FieldData"][1]
    names = String[]
    data_arrays = XMLElement[]
    numbers = Int[]
    types = []
    offsets = Int[]
    for xml_element in child_elements(field_data)
        if LightXML.name(xml_element) in ["DataArray"]
            push!(names, LightXML.attribute(xml_element, "Name", required = true))
            push!(data_arrays, xml_element)
            push!(numbers, parse(Int, LightXML.attribute(xml_element, "NumberOfComponents")))
            push!(types, kStringToTypeDict[LightXML.attribute(xml_element, "type")])
            push!(offsets, parse(Int, LightXML.attribute(xml_element, "offset")))
        else
            println("Unknown element: ", LightXML.name(xml_element))
            continue
        end
    end
    vtp_data_arrays =
        [VTKDataArray{types[i], numbers[i], ReadVTK.FormatAppended}(names[i], offsets[i], data_arrays[i], vtp_file) for i in eachindex(names)]
    return DataFrame([get_data(fd) for fd in vtp_data_arrays], [fd.name for fd in vtp_data_arrays])
end

mutable struct VTPIODataFrameView{IntType <: Int} <: AbstractDataIO
    vtp_file_::VTKFile
    dir_path_::String
    dim_::IntType
    grouped_points_data_frame_::GroupedDataFrame
    field_data_frame_::DataFrame
end

@inline function VTPIODataFrameView(vtp_file::VTKFile)::VTPIODataFrameView
    grouped_points_data_frame, dim = getPointsDataFrameFromVTPFile(vtp_file)
    field_data_frame = getFieldDataFrameFromVTPFile(vtp_file)
    file_path = vtp_file.filename
    dir_path = dirname(file_path)
    return VTPIODataFrameView{typeof(dim)}(vtp_file, dir_path, dim, grouped_points_data_frame, field_data_frame)
end

VTPIODataFrameView(vtp_file_name::String)::VTPIODataFrameView = VTPIODataFrameView(VTKFile(vtp_file_name))

@inline function insertMassColumn!(data_frame::DataFrameType where {DataFrameType <: AbstractDataFrame}; mass = NaN)::Nothing
    insertcols!(data_frame, Symbol(kDensityString), Symbol(kMassString) => mass, after = true)
    return nothing
end

@inline function insertDensityColumn!(data_frame::DataFrameType where {DataFrameType <: AbstractDataFrame}; density = NaN)::Nothing
    insertcols!(data_frame, Symbol(kTypeString), Symbol(kDensityString) => density, after = true)
    return nothing
end

function kernelValueInterpolationFromDataFrame!(
    interpolated_points::Matrix,
    dataframe::DataFrameType,
    field::Symbol,
    dim::IntType,
    smooth_kernel::SmoothKernel,
    kernel_weights::AtomicArrayType,
    weighted_values::AtomicArrayType,
)::Nothing where {
    IntType <: Integer,
    RealType <: AbstractFloat,
    AtomicArrayType <: AbstractVector{Threads.Atomic{RealType}},
    DataFrameType <: AbstractDataFrame,
}
    @assert dim == size(interpolated_points, 1) == smooth_kernel.dim_
    @assert kPositionStrings[1:dim] == names(dataframe)[1:dim]
    @assert hasproperty(dataframe, field) == true
    @assert hasproperty(dataframe, kDensityString) == true
    @assert hasproperty(dataframe, kMassString) == true
    @views points = dataframe[!, 1:dim] |> Matrix |> transpose |> Matrix
    neighbours = neighborlist(interpolated_points, points, smooth_kernel.influence_radius_, parallel = true)
    @floop @simd for neighbour in neighbours
        i, j, r = neighbour
        kernel_value = kernelValue(r, smooth_kernel)
        @inbounds mass_j = dataframe[j, kMassString]
        @inbounds rho_j = dataframe[j, kDensityString]
        @inbounds field_j = dataframe[j, field]
        @inbounds Threads.atomic_add!(kernel_weights[i], kernel_value * mass_j / rho_j)
        @inbounds Threads.atomic_add!(weighted_values[i], kernel_value * mass_j / rho_j * field_j)
    end
    return nothing
end

function kernelValueInterpolation(
    interpolated_points::Matrix,
    vtpio_data_frame_view::VTPIODataFrameView,
    field::Symbol,
    smooth_kernel::SmoothKernel,
    data_frame_index_list::IntArrayType;
    none_value::RealType where {RealType <: AbstractFloat} = NaN,
)::AbstractVector where {IntType <: Integer, IntArrayType <: AbstractVector{IntType}}
    RealType = vtpio_data_frame_view.grouped_points_data_frame_[1][1, field] |> typeof
    n_interpolated_points = size(interpolated_points, 2)
    kernel_weights = [Threads.Atomic{RealType}(0.0) for _ in 1:n_interpolated_points]
    weighted_values = [Threads.Atomic{RealType}(0.0) for _ in 1:n_interpolated_points]
    for data_frame_index in data_frame_index_list
        kernelValueInterpolationFromDataFrame!(
            interpolated_points,
            vtpio_data_frame_view.grouped_points_data_frame_[data_frame_index],
            field,
            vtpio_data_frame_view.dim_,
            smooth_kernel,
            kernel_weights,
            weighted_values,
        )
    end
    interpolated_values = zeros(RealType, n_interpolated_points)
    @floop @simd for i in 1:n_interpolated_points
        if kernel_weights[i][] > 0
            @inbounds interpolated_values[i] = weighted_values[i][] / kernel_weights[i][]
        else
            @inbounds interpolated_values[i] = none_value
        end
    end
    return interpolated_values
end
