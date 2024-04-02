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
    const file_dir = "./NaturalConvectionCavityNormalizationData/"
    const file_lists = readdir(file_dir) |> sort!
    const total_time = 30.0
    vtpFile(step::Int64) = file_dir * file_lists[step + 1]
end

# ╔═╡ 0c877a7a-597f-47aa-854c-a26126bf5434
function getData(step::Int64)
    vtpio_dataframe_view = VTPIODataFrameView(vtpFile(step))
    @views points = vtpio_dataframe_view.vtp_file_ |> get_points
    @views x = points[1, :]
    @views y = points[2, :]
    @views temp = get_point_data(vtpio_dataframe_view.vtp_file_)["Temperature"] |> get_data
    t = vtpio_dataframe_view.field_data_frame_.Time[1]
    return x, y, temp, t
    return
end

# ╔═╡ e889ccf2-fb7c-4da2-8a72-6f305657079f
function plot(step::Int64)
    x, y, temp, t = getData(step)
    fig = Figure(size = (600, 400))
    ax1 = Axis(fig[1, 1], title = "Step=$(string(step))  Time=$(string(round(t, digits=4)))  Temperature", aspect = DataAspect())
    particles = meshscatter!(
        x,
        y,
        color = temp,
        colormap = :balance,
        nan_color = :lightgrey,
        markersize = 0.01,
        shading = FastShading,
        disffuse = Vec3f(0.1),
        specular = Vec3f(0.1),
        shininess = 20,
        backlight = 0.0f0,
        ssao = true,
    )
    cbar = Colorbar(fig[1, 2], particles, label = "Temperature")
    cbar.ticks = 0:0.1:1
    return fig
end

# ╔═╡ 92864fce-9e09-4329-a1cd-677c16892289
@bind step PlutoUI.Slider(0:(length(file_lists) - 1), default = 0)

# ╔═╡ 0ce27dcf-fe29-47c1-8105-4d81d33308ad
plot(step)

# ╔═╡ Cell order:
# ╠═b25ae9a0-e119-11ee-3e61-c764c01ad4b3
# ╠═56e1d827-c33a-479f-b8de-67d0c7933e3a
# ╠═0c877a7a-597f-47aa-854c-a26126bf5434
# ╠═e889ccf2-fb7c-4da2-8a72-6f305657079f
# ╠═92864fce-9e09-4329-a1cd-677c16892289
# ╠═0ce27dcf-fe29-47c1-8105-4d81d33308ad
