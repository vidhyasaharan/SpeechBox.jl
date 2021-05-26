using Documenter, SpeechBox

makedocs(
    sitename = "SpeechBox.jl",
    authors = "Vidhyasaharan Sethu",
    pages = [
        "Home" => "index.md",
        "Contents" => [
            "signal_objects.md",
            "speechframing.md",
            "vad.md",
            "spectralanalyses.md",
            "lpc.md",
            "mfcc.md"
        ],
    ],
)

deploydocs(
    repo = "github.com/unsw-edu-au/SpeechBox.jl.git",
)