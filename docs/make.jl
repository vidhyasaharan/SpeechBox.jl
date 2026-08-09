using Documenter, SpeechBox

DocMeta.setdocmeta!(SpeechBox, :DocTestSetup, :(using SpeechBox); recursive=true)

makedocs(
    sitename = "SpeechBox.jl",
    authors = "Vidhyasaharan Sethu",
    modules = [SpeechBox],
    checkdocs = :exports,
    pages = [
        "Home" => "index.md",
        "Contents" => [
            "signal_objects.md",
            "speechframing.md",
            "vad.md",
            "spectralanalyses.md",
            "correlations.md",
            "lpc.md",
            "pitch.md",
            "mfcc.md",
            "kmeans.md",
            "utilities.md"
        ],
    ],
)

deploydocs(
    repo = "github.com/unsw-edu-au/SpeechBox.jl.git",
)
