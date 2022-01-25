using ParMETIS
using Documenter

DocMeta.setdocmeta!(ParMETIS, :DocTestSetup, :(using ParMETIS); recursive = true)

makedocs(;
  modules = [ParMETIS],
  doctest = true,
  linkcheck = false,
  strict = false,
  authors = "Abel Soares Siqueira <abel.s.siqueira@gmail.com> and contributors",
  repo = "https://github.com/JuliaSmoothOptimizers/ParMETIS.jl/blob/{commit}{path}#{line}",
  sitename = "ParMETIS.jl",
  format = Documenter.HTML(;
    prettyurls = get(ENV, "CI", "false") == "true",
    canonical = "https://JuliaSmoothOptimizers.github.io/ParMETIS.jl",
    assets = ["assets/style.css"],
  ),
  pages = ["Home" => "index.md", "Reference" => "reference.md"],
)

deploydocs(;
  repo = "github.com/JuliaSmoothOptimizers/ParMETIS.jl",
  push_preview = true,
  devbranch = "main",
)
