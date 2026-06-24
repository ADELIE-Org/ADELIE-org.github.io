using Documenter

makedocs(
    authors = "ADELIE-org contributors",
    sitename = "ADELIE",
    format = Documenter.HTML(
        canonical = "https://ADELIE-org.github.io",
        repolink = "https://github.com/ADELIE-org",
        collapselevel = 1,
        prettyurls = true,
        inventory_version = "0",
    ),
    pages = [
        "Home" => "index.md",
        "Packages" => "packages.md",
    ],
    warnonly = true,
    remotes = nothing,
)

if get(ENV, "CI", "") == "true"
    deploydocs(
        repo = "github.com/ADELIE-org/ADELIE-org.github.io",
        branch = "gh-pages",
        push_preview = false,
    )
end
