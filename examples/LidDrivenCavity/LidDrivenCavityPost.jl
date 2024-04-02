### A Pluto.jl notebook ###
# v0.19.40

#> [frontmatter]
#> title = "CollapseDryPost"
#> description = "collapse dry SPH demo  post-process"
#> 
#>     [[frontmatter.author]]
#>     name = "bcynuaa"

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try
            Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value
        catch
            b -> missing
        end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 91c675d8-e105-11ee-23a2-e30fe3863c8b
begin
    using ReadVTK
    using WriteVTK
    using GLMakie
    using Makie
    Makie.inline!(true)
    using PlutoUI
    using Pkg
    Pkg.activate("../..")
    using EtherSPHCPU
end

# ╔═╡ 37768d7e-5e33-41cf-85ec-7d90d24cefb6
begin
    const file_dir = "./LidDrivenCavityData/"
    const file_lists = readdir(file_dir) |> sort!
    const total_time = 1.0
    vtpFile(step::Int64) = file_dir * file_lists[step + 1]
end

# ╔═╡ 7ffedf25-d515-4bf1-80eb-5f9a5331e347
begin
    const dr = 0.01
    const rho = 1.0
    const dim = 2
    const influence_radius = 3 * dr
    const smooth_kernel = SmoothKernel(influence_radius, dim, WendlandC2)
    middle_y = LinRange(0.0, 1.0, 101)
    middle_x = zeros(length(middle_y)) .+ 0.5
    middle_points = hcat(middle_x, middle_y) |> transpose |> Matrix
    ref_y = [1.0, 0.9766, 0.9688, 0.9609, 0.9531, 0.8516, 0.7344, 0.6172, 0.5, 0.4531, 0.2813, 0.1719, 0.1016, 0.0703, 0.0625, 0.0547, 0.0]
    ref_u = [
        1.0,
        0.84123,
        0.78871,
        0.73722,
        0.68717,
        0.23151,
        0.00332,
        -0.13641,
        -0.20581,
        -0.2109,
        -0.15662,
        -0.1015,
        -0.06434,
        -0.04775,
        -0.04192,
        -0.03717,
        0.0,
    ]
end

# ╔═╡ 5481ff11-75f4-47f6-974e-2ef8c43c7d77
function getData(step::Int64)
    vtpio_dataframe_view = VTPIODataFrameView(vtpFile(step))
    vtpio_dataframe_view.grouped_points_data_frame_[1][:, kMassString] .= dr^2 * 1
    vtpio_dataframe_view.grouped_points_data_frame_[2][:, kDensityString] .= 1.0
    vtpio_dataframe_view.grouped_points_data_frame_[2][:, kMassString] .= dr^2 * 1
    middle_u = kernelValueInterpolation(middle_points, vtpio_dataframe_view, :VelocityX, smooth_kernel, [1, 2])
    @views points = vtpio_dataframe_view.vtp_file_ |> get_points
    @views x = points[1, :]
    @views y = points[2, :]
    @views vs = get_point_data(vtpio_dataframe_view.vtp_file_)["Velocity"] |> get_data
    v = sum(vs .^ 2, dims = 1)[1, :] .|> sqrt
    t = vtpio_dataframe_view.field_data_frame_.Time[1]
    return x, y, v, middle_u, t
    return
end

# ╔═╡ d4959301-7222-439f-b56b-951900032908
function plot(step::Int64)
    x, y, v, middle_u, t = getData(step)
    fig = Figure(size = (1200, 400))
    ax1 = Axis(fig[1, 1], title = "Step=$(string(step))  Time=$(string(round(t, digits=4)))  Velocity", aspect = DataAspect())
    particles = meshscatter!(
        x,
        y,
        color = v,
        colormap = :turbo,
        colorrange = (0, 1),
        nan_color = :lightgrey,
        markersize = 0.007,
        shading = FastShading,
        disffuse = Vec3f(0.1),
        specular = Vec3f(0.1),
        shininess = 20,
        backlight = 0.0f0,
        ssao = true,
    )
    cbar = Colorbar(fig[1, 2], particles, label = "Velocity Norm", size = 20)
    cbar.ticks = 0:0.1:1
    cbar.limits = (0, 1)
    ax2 = Axis(fig[1, 3], title = "Step=$(string(step))  Time=$(string(round(t, digits=4)))  Middle U", xlabel = L"$x$", ylabel = L"$y$")
    lines!(ax2, middle_y, middle_u, color = :blue, label = "SPH Code-Development")
    scatter!(ax2, ref_y, ref_u, color = :red, label = "Reference-Steady")
    axislegend(ax2, position = :lt)
    return fig
end

# ╔═╡ 580dd940-4f4e-4d3f-a9eb-58af93af9e9e
@bind step PlutoUI.Slider(0:(length(file_lists) - 1), default = 0)

# ╔═╡ 8d72d6e3-33a3-4e79-bd77-94ed18e56f86
plot(step)

# ╔═╡ Cell order:
# ╠═91c675d8-e105-11ee-23a2-e30fe3863c8b
# ╠═37768d7e-5e33-41cf-85ec-7d90d24cefb6
# ╠═7ffedf25-d515-4bf1-80eb-5f9a5331e347
# ╠═5481ff11-75f4-47f6-974e-2ef8c43c7d77
# ╠═d4959301-7222-439f-b56b-951900032908
# ╠═580dd940-4f4e-4d3f-a9eb-58af93af9e9e
# ╠═8d72d6e3-33a3-4e79-bd77-94ed18e56f86
