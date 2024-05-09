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

# ╔═╡ eacfb796-0542-11ef-38c1-759ed83492b3
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

# ╔═╡ 011ca230-bffb-4371-ba42-f821322e4bc4
begin
    const file_dir = "./PoiseuilleFlow2DData/"
    const file_lists = readdir(file_dir) |> sort!
    vtpFile(step::Int64) = file_dir * file_lists[step + 1]
end

# ╔═╡ b670dc62-cad8-4256-b86f-e17658dcf623
begin
    const dr = 2e-4
    const dim = 2
    const influence_radius = 3 * dr
    const rho_0 = 1000.0
    const smooth_kernel = SmoothKernel(influence_radius, dim, WendlandC2)
    const mass = rho_0 * dr^2
    const l = 0.01
    const monitor_y = LinRange(0, l, 51)
    const monitor_x = zeros(length(monitor_y)) .+ l * 2
    const monitor_points = hcat(monitor_x, monitor_y) |> transpose |> Matrix
    const ref_y = LinRange(0, l, 11)
end

# ╔═╡ 615fc739-4081-4e26-bf35-cc246f7c3396
function uxTheory(y::Float64, t::Float64)::Float64
    nu = 1e-6
    f = 8e-5
    poiseuille_flow_theory_solutio_order = 10
    ux = f / 2 / nu * y * (l - y)
    for n in 0:poiseuille_flow_theory_solutio_order
        ux -= 4 * f * l^2 / nu / pi^3 / (2 * n + 1)^3 * sin(pi * y / l * (2 * n + 1)) * exp(-(2 * n + 1)^2 * pi^2 * nu * t / l^2)
    end
    return ux
end

# ╔═╡ 8982db5d-0e49-4f98-b0c8-62fe32e9686b
function getData(step::Int64)
    vtpio_dataframe_view = VTPIODataFrameView(vtpFile(step))
    @views points = vtpio_dataframe_view.vtp_file_ |> get_points
    @views x = points[1, :]
    @views y = points[2, :]
    @views vx = (get_point_data(vtpio_dataframe_view.vtp_file_)["Velocity"] |> get_data)[1, :]
    t = vtpio_dataframe_view.field_data_frame_.Time[1]
    vtpio_dataframe_view.grouped_points_data_frame_[1][:, kMassString] .= dr^2 * vtpio_dataframe_view.grouped_points_data_frame_[1][:, kDensityString]
    vtpio_dataframe_view.grouped_points_data_frame_[4][:, kDensityString] .= rho_0
    vtpio_dataframe_view.grouped_points_data_frame_[4][:, kMassString] .= mass
    vtpio_dataframe_view.grouped_points_data_frame_[4][:, "VelocityX"] .= 0.0
    monitor_vx = kernelValueInterpolation(monitor_points, vtpio_dataframe_view, :VelocityX, smooth_kernel, [1, 2, 3, 4])
    ref_vx = uxTheory.(ref_y, t)
    return x, y, vx, t, monitor_vx, ref_vx
end

# ╔═╡ c71cf81e-ae20-4377-983a-71791b1bc0ce
function plot(step::Int64)
    x, y, vx, t, monitor_vx, ref_vx = getData(step)
    fig = Figure(size = (1200, 400))
    ax1 = Axis(fig[1, 1], title = L"Step=%$(string(step))  Time=%$(string(round(t, digits=4)))  $V_x$ ", aspect = DataAspect())
    particles = meshscatter!(
        x,
        y,
        color = vx,
        colormap = :turbo,
        nan_color = :lightgrey,
        markersize = 0.0002,
        shading = FastShading,
        disffuse = Vec3f(0.1),
        specular = Vec3f(0.1),
        shininess = 20,
        backlight = 0.0f0,
        ssao = true,
    )
    cbar = Colorbar(fig[1, 2], particles, label = L"$V_x$")
    cbar.ticks = 0:0.1:1
    ax2 = Axis(fig[1, 3], title = L"$V_x$ along y", xlabel = L"$V_x$", ylabel = L"$y$", limits = (0.0, 0.0011, 0.0, 0.01))
    lines!(monitor_vx, monitor_y, label = L"SPH-$V_x$", color = :blue)
    scatter!(ref_vx, ref_y, label = L"stable ref-$V_x$", color = :red, markersize = 20)
    axislegend(ax2, position = :rt)
    return fig
end

# ╔═╡ 3b5efb16-206b-4c6f-ae26-f2b433bf4e02
@bind step PlutoUI.Slider(0:(length(file_lists) - 1), default = 0)

# ╔═╡ a1e5635f-3384-4a55-8027-a5db3d1c666a
plot(step)

# ╔═╡ Cell order:
# ╠═eacfb796-0542-11ef-38c1-759ed83492b3
# ╠═011ca230-bffb-4371-ba42-f821322e4bc4
# ╠═b670dc62-cad8-4256-b86f-e17658dcf623
# ╠═615fc739-4081-4e26-bf35-cc246f7c3396
# ╠═8982db5d-0e49-4f98-b0c8-62fe32e9686b
# ╠═c71cf81e-ae20-4377-983a-71791b1bc0ce
# ╠═3b5efb16-206b-4c6f-ae26-f2b433bf4e02
# ╠═a1e5635f-3384-4a55-8027-a5db3d1c666a
