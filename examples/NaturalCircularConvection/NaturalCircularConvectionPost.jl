### A Pluto.jl notebook ###
# v0.19.40

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

# ╔═╡ b25ae9a0-e119-11ee-3e61-c764c01ad4b3
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

# ╔═╡ 56e1d827-c33a-479f-b8de-67d0c7933e3a
begin
    const file_dir = "./NaturalCircularConvectionData/"
    const file_lists = readdir(file_dir) |> sort!
    const total_time = 30.0
    vtpFile(step::Int64) = file_dir * file_lists[step + 1]
end

# ╔═╡ 09174b63-cece-4bed-862e-92996cff491b
begin
    const dr = 0.02
    const dim = 2
    const influence_radius = 3 * dr
    const rho_0 = 1.1538824720077858
    const smooth_kernel = SmoothKernel(influence_radius, dim, WendlandC2)
    const n_points_along_radius = 21
    thetas = [0, 30, 60, 90, 120, 150, 180]
    theta_s = Float64.(thetas)
    theta_s .+= 90
    theta_s .*= pi / 180
    const n_theta = length(theta_s)
    points_along_radius = zeros(2, n_theta * n_points_along_radius)
    r_inner = 1 / 2.6
    r_outer = 1.0
    r_s = LinRange(r_inner, r_outer, n_points_along_radius)
    for i_theta in eachindex(theta_s)
        for i_point_theta in eachindex(r_s)
            x, y = r_s[i_point_theta] * [cos(theta_s[i_theta]), sin(theta_s[i_theta])]
            index = (i_theta - 1) * n_points_along_radius + i_point_theta
            points_along_radius[:, index] .= [x, y]
        end
    end
    const ref_radius = LinRange(0, 1, n_points_along_radius)
end

# ╔═╡ 0c877a7a-597f-47aa-854c-a26126bf5434
function getData(step::Int64)
    vtpio_dataframe_view = VTPIODataFrameView(vtpFile(step))
    @views points = vtpio_dataframe_view.vtp_file_ |> get_points
    @views x = points[1, :]
    @views y = points[2, :]
    @views temp = get_point_data(vtpio_dataframe_view.vtp_file_)["Temperature"] |> get_data
    t = vtpio_dataframe_view.field_data_frame_.Time[1]
    vtpio_dataframe_view.grouped_points_data_frame_[1][:, kMassString] .= dr^2 * vtpio_dataframe_view.grouped_points_data_frame_[1][:, kDensityString]
    vtpio_dataframe_view.grouped_points_data_frame_[2][:, kDensityString] .= rho_0
    vtpio_dataframe_view.grouped_points_data_frame_[2][:, kMassString] .= dr^2 * rho_0
    t_along_radius = kernelValueInterpolation(points_along_radius, vtpio_dataframe_view, :Temperature, smooth_kernel, [1, 2])
    t_along_radius = reshape(t_along_radius, (n_points_along_radius, n_theta))
    t_along_radius = Matrix(t_along_radius')
    return x, y, temp, t, t_along_radius
end

# ╔═╡ e889ccf2-fb7c-4da2-8a72-6f305657079f
function plot(step::Int64)
    x, y, temp, t, t_along_radius = getData(step)
    fig = Figure(size = (1400, 400))
    ax1 = Axis(fig[1, 1], title = "Step=$(string(step))  Time=$(string(round(t, digits=4)))  Temperature", aspect = DataAspect())
    particles = meshscatter!(
        x,
        y,
        color = temp,
        colormap = :balance,
        nan_color = :lightgrey,
        markersize = 0.02,
        shading = FastShading,
        disffuse = Vec3f(0.1),
        specular = Vec3f(0.1),
        shininess = 20,
        backlight = 0.0f0,
        ssao = true,
    )
    cbar = Colorbar(fig[1, 2], particles, label = "Temperature")
    cbar.ticks = 0:0.1:1
    ax2 = Axis(fig[1, 3], title = "Temperature Along Radius")
    for i_theta in 1:n_theta
        lines!(ref_radius, t_along_radius[i_theta, :], label = "SPH-$(thetas[i_theta])")
    end
    axislegend(ax2, position = :rt)
    return fig
end

# ╔═╡ 92864fce-9e09-4329-a1cd-677c16892289
@bind step PlutoUI.Slider(0:(length(file_lists) - 1), default = 0)

# ╔═╡ 0ce27dcf-fe29-47c1-8105-4d81d33308ad
plot(step)

# ╔═╡ Cell order:
# ╠═b25ae9a0-e119-11ee-3e61-c764c01ad4b3
# ╠═56e1d827-c33a-479f-b8de-67d0c7933e3a
# ╠═09174b63-cece-4bed-862e-92996cff491b
# ╠═0c877a7a-597f-47aa-854c-a26126bf5434
# ╠═e889ccf2-fb7c-4da2-8a72-6f305657079f
# ╠═92864fce-9e09-4329-a1cd-677c16892289
# ╠═0ce27dcf-fe29-47c1-8105-4d81d33308ad
