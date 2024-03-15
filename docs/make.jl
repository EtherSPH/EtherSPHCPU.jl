#=
  @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com>
  @ date: 2024/03/15 13:38:39
  @ license: MIT
  @ description:
 =#

using Documenter
using EtherSPHCPU

makedocs(sitename = "EtherSPHCPU", format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"))

deploydocs(root = "./", repo = "github.com:EtherSPH/EtherSPHCPU.git", target = "build", branch = "gh-pages")
