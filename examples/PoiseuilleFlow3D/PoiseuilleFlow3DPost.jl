### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ c91e55b4-0df6-11ef-2ae1-f96233ba3d33
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

# ╔═╡ b05fc988-d804-465d-a80b-4fe797918f36
begin
    const file_dir = "./PoiseuilleFlow3DData/"
    const file_lists = readdir(file_dir) |> sort!
    vtpFile(step::Int64) = file_dir * file_lists[step + 1]
end

# ╔═╡ d22c2440-52d6-425a-81ed-4cc12cf0c1bb
begin
	const ll = 0.01
    const dr = 0.01 / 40
    const dim = 3
    const influence_radius = 3 * dr
    const rho_0 = 1000.0
    const smooth_kernel = SmoothKernel(influence_radius, dim, WendlandC2)
    const mass = rho_0 * dr^3
    const l = 0.01
    const r0 = l / 2;
	const monitor_y = -r0:r0/40:r0
	const monitor_z = deepcopy(monitor_y)
	const monitor_x = l
	monitor_points = zeros(length(monitor_y) * length(monitor_z), 3)
	for j in 1: length(monitor_y)
		for k in 1: length(monitor_z)
			monitor_points[(j-1) * length(monitor_y) + k, :] .= [monitor_x, monitor_y[j], monitor_z[k]]
		end
	end
	monitor_points = monitor_points' |> Matrix
end

# ╔═╡ 26de3af6-fe61-4f1e-a241-c9ac07f06ece
function getData(step::Int64)
    vtpio_dataframe_view = VTPIODataFrameView(vtpFile(step))
    @views points = vtpio_dataframe_view.vtp_file_ |> get_points
    @views x = points[1, :]
    @views y = points[2, :]
	@views z = points[3, :]
    @views vx = (get_point_data(vtpio_dataframe_view.vtp_file_)["Velocity"] |> get_data)[1, :]
    t = vtpio_dataframe_view.field_data_frame_.Time[1]
    vtpio_dataframe_view.grouped_points_data_frame_[1][:, kMassString] .= dr^3 * vtpio_dataframe_view.grouped_points_data_frame_[1][:, kDensityString]
    vtpio_dataframe_view.grouped_points_data_frame_[2][:, kDensityString] .= rho_0
    vtpio_dataframe_view.grouped_points_data_frame_[2][:, kMassString] .= mass
    vtpio_dataframe_view.grouped_points_data_frame_[2][:, "VelocityX"] .= 0.0
    monitor_vx = kernelValueInterpolation(monitor_points, vtpio_dataframe_view, :VelocityX, smooth_kernel, [1, 2]; none_value=NaN)
	monitor_vx = reshape(monitor_vx, (length(monitor_y), length(monitor_z)))
	for j in 1: length(monitor_y)
		for k in 1: length(monitor_z)
			yy = monitor_y[j]
			zz = monitor_z[k]
			if yy^2+zz^2>r0^2
				monitor_vx[j, k] = NaN
			end
		end
	end
    return x, y, z, vx, t, monitor_vx
end

# ╔═╡ 56044fe7-c755-482c-96d5-7595f7980c18
function plot(step::Int64, view::Int64)
    x, y, z, vx, t, monitor_vx = getData(step)
    fig = Figure(size = (1200, 400))
    ax1 = Axis3(fig[1, 1], title = L"Step=%$(string(step))  Time=%$(string(round(t, digits=4)))  $V_x$ ", aspect = (1, 1, 1), azimuth = view * pi / 20)
    particles = meshscatter!(
        x,
        y,
		z,
        color = vx,
        colormap = :turbo,
        # nan_color = :lightgrey,
        markersize = 0.0001,
        shading = FastShading,
        disffuse = Vec3f(0.1),
        specular = Vec3f(0.1),
        shininess = 20,
        backlight = 0.0f0,
        ssao = true,
		alpha = 0.2,
    )
    cbar = Colorbar(fig[1, 2], particles, label = L"$V_x$")
	cbar.ticks = 0:0.00005:0.001
	ax2 = Axis(fig[1, 3], aspect = DataAspect())
	contourf!(ax2, monitor_y, monitor_z, monitor_vx, colormap=:turbo, levels=21)
	cbar2 = Colorbar(fig[1, 4], particles, label = L"$V_x$")
	cbar2.ticks = 0:0.00005:0.001
    return fig
end

# ╔═╡ 7bc1cba7-958c-4dcb-94d0-f9d8facb84a8
@bind step PlutoUI.Slider(0:(length(file_lists) - 1), default = 0)

# ╔═╡ 7e5eca22-067a-41e5-8f8d-647bd1141956
@bind view PlutoUI.Slider(-10:10, default = 0)

# ╔═╡ e8952ba5-6022-4629-8220-aa3d8f3eb6f7
plot(step, view)

# ╔═╡ Cell order:
# ╠═c91e55b4-0df6-11ef-2ae1-f96233ba3d33
# ╠═b05fc988-d804-465d-a80b-4fe797918f36
# ╠═d22c2440-52d6-425a-81ed-4cc12cf0c1bb
# ╠═26de3af6-fe61-4f1e-a241-c9ac07f06ece
# ╠═56044fe7-c755-482c-96d5-7595f7980c18
# ╠═7bc1cba7-958c-4dcb-94d0-f9d8facb84a8
# ╠═7e5eca22-067a-41e5-8f8d-647bd1141956
# ╠═e8952ba5-6022-4629-8220-aa3d8f3eb6f7
