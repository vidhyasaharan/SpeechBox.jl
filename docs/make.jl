using Documenter, SpeechBox

makedocs(
    sitename = "SpeechBox.jl",
    authors = "Vidhyasaharan Sethu",
    Pages = [
        "Home" => "index.md",
        "Library" => [
            "signal_objects.md",
        "speechframing.md",
        "spectralanalyses.md"
        ],
    ],
)

deploydocs(
    repo = "github.com/unsw-edu-au/SpeechBox.jl.git",
)