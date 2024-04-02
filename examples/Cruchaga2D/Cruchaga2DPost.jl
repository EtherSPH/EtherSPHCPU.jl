### A Pluto.jl notebook ###
# v0.19.40

#> [frontmatter]
#> title = "Cruchaga2DPost"
#> description = "Cruchaga 2D SPH post-processing"
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

# ╔═╡ 6a25a422-e117-11ee-2f34-b57460ca2983
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

# ╔═╡ 2506e8ec-d5ce-4404-b797-bea5dd6f3234
begin
    const file_dir = "./Cruchaga2DData/"
    const file_lists = readdir(file_dir) |> sort!
    const total_time = 1.0
    vtpFile(step::Int64) = file_dir * file_lists[step + 1]
end

# ╔═╡ 837e3eb4-9d7e-4616-a79b-baf65fe2f637
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

# ╔═╡ 10f88558-54dd-4a82-9c99-39083bc22618
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
        markersize = 0.003,
        shading = FastShading,
        disffuse = Vec3f(0.1),
        specular = Vec3f(0.1),
        shininess = 20,
        backlight = 0.0f0,
        ssao = true,
    )
    cbar = Colorbar(fig[1, 2], particles, label = "Velocity Norm")
    cbar.ticks = 0:0.1:2
    return fig
end

# ╔═╡ ace83ba7-8026-4f0c-a843-7398b48ed9ac
@bind step PlutoUI.Slider(0:(length(file_lists) - 1), default = 0)

# ╔═╡ b58fea1b-ac92-4045-8e6c-49e6447a8609
plot(step)

# ╔═╡ Cell order:
# ╠═6a25a422-e117-11ee-2f34-b57460ca2983
# ╠═2506e8ec-d5ce-4404-b797-bea5dd6f3234
# ╠═837e3eb4-9d7e-4616-a79b-baf65fe2f637
# ╠═10f88558-54dd-4a82-9c99-39083bc22618
# ╠═ace83ba7-8026-4f0c-a843-7398b48ed9ac
# ╠═b58fea1b-ac92-4045-8e6c-49e6447a8609
