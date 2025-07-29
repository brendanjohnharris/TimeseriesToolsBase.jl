using TimeseriesToolsBase
using Unitful
using Documenter

DocMeta.setdocmeta!(TimeseriesToolsBase, :DocTestSetup,
                    :(using Unitful, TimeseriesToolsBase);
                    recursive = true)

makedocs(;
         modules = [TimeseriesToolsBase],
         authors = "brendanjohnharris <brendanjohnharris@gmail.com> and contributors",
         repo = "https://github.com/brendanjohnharris/TimeseriesToolsBase.jl/blob/{commit}{path}#{line}",
         sitename = "TimeseriesToolsBase.jl",
         format = Documenter.HTML(;
                                  prettyurls = get(ENV, "CI", "false") == "true",
                                  canonical = "https://brendanjohnharris.github.io/TimeseriesToolsBase.jl",
                                  edit_link = "main",
                                  assets = String[],),
         pages = ["Home" => "index.md",
             "Types" => "toolsarrays.md",
             "Utils" => "utils.md"],)

deploydocs(;
           repo = "github.com/brendanjohnharris/TimeseriesToolsBase.jl",
           devbranch = "main",)
