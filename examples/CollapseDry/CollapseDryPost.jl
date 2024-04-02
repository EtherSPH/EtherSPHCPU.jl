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
    const file_dir = "./CollapseDryData/"
    const file_lists = readdir(file_dir) |> sort!
    const total_time = 3.0
    vtpFile(step::Int64) = file_dir * file_lists[step + 1]
end

# ╔═╡ 5481ff11-75f4-47f6-974e-2ef8c43c7d77
function getData(step::Int64)
    vtpio_dataframe_view = VTPIODataFrameView(vtpFile(step))
    @views points = vtpio_dataframe_view.vtp_file_ |> get_points
    @views x = points[1, :]
    @views y = points[2, :]
    @views vs = get_point_data(vtpio_dataframe_view.vtp_file_)["Velocity"] |> get_data
    v = sum(vs .^ 2, dims = 1)[1, :] .|> sqrt
    t = vtpio_dataframe_view.field_data_frame_.Time[1]
    return x, y, v, t
    return
end

# ╔═╡ d4959301-7222-439f-b56b-951900032908
function plot(step::Int64)
    x, y, v, t = getData(step)
    fig = Figure(size = (600, 400))
    ax1 = Axis(fig[1, 1], title = "Step=$(string(step))  Time=$(string(round(t, digits=4)))  Velocity", aspect = DataAspect())
    particles = meshscatter!(
        x,
        y,
        color = v,
        colormap = :turbo,
        nan_color = :lightgrey,
        markersize = 0.01,
        shading = FastShading,
        disffuse = Vec3f(0.1),
        specular = Vec3f(0.1),
        shininess = 20,
        backlight = 0.0f0,
        ssao = true,
    )
    cbar = Colorbar(fig[1, 2], particles, label = "Velocity Norm")
    cbar.ticks = 0:1:20
    return fig
end

# ╔═╡ 580dd940-4f4e-4d3f-a9eb-58af93af9e9e
@bind step PlutoUI.Slider(0:(length(file_lists) - 1), default = 0)

# ╔═╡ 8d72d6e3-33a3-4e79-bd77-94ed18e56f86
plot(step)

# ╔═╡ Cell order:
# ╠═91c675d8-e105-11ee-23a2-e30fe3863c8b
# ╠═37768d7e-5e33-41cf-85ec-7d90d24cefb6
# ╠═5481ff11-75f4-47f6-974e-2ef8c43c7d77
# ╠═d4959301-7222-439f-b56b-951900032908
# ╠═580dd940-4f4e-4d3f-a9eb-58af93af9e9e
# ╠═8d72d6e3-33a3-4e79-bd77-94ed18e56f86
