using TimeseriesBase
using Unitful
using Documenter

DocMeta.setdocmeta!(TimeseriesBase, :DocTestSetup,
                    :(using Unitful, TimeseriesBase);
                    recursive = true)

makedocs(;
         modules = [TimeseriesBase],
         authors = "brendanjohnharris <brendanjohnharris@gmail.com> and contributors",
         repo = "https://github.com/brendanjohnharris/TimeseriesBase.jl/blob/{commit}{path}#{line}",
         sitename = "TimeseriesBase.jl",
         format = Documenter.HTML(;
                                  prettyurls = get(ENV, "CI", "false") == "true",
                                  canonical = "https://brendanjohnharris.github.io/TimeseriesBase.jl",
                                  edit_link = "main",
                                  assets = String[],),
         pages = ["Home" => "index.md",
             "Types" => "types.md",
             "Utils" => "utils.md"],)

deploydocs(;
           repo = "github.com/brendanjohnharris/TimeseriesBase.jl",
           devbranch = "main",)
