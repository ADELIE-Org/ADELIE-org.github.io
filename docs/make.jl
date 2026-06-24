using Documenter

makedocs(
    authors = "ADELIE-org contributors",
    sitename = "ADELIE",
    format = Documenter.HTML(
        canonical = "https://ADELIE-Org.github.io",
        repolink = "https://github.com/ADELIE-Org",
        collapselevel = 1,
        prettyurls = true,
        inventory_version = "0",
    ),
    pages = [
        "Home" => "index.md",
        "Packages" => "packages.md",
    ],
    pagesonly = true,
    warnonly = true,
    remotes = nothing,
)

if get(ENV, "CI", "") == "true"
    deploydocs(
        repo = "github.com/ADELIE-Org/ADELIE-org.github.io",
        push_preview = false,
    )
end
