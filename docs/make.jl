using Documenter, SpeechBox

makedocs(
    sitename = "My Documentation",
    Pages = Any[
        "Library" => "index.md",
        "signal_objects.md",
        "speechframing.md",
        "spectralanalyses.md"
    ],
)

deploydocs(
    repo = "github.com/unsw-edu-au/SpeechBox.jl.git",
)