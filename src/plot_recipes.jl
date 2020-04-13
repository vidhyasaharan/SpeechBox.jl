#Plot recipe for plotting spectra
@recipe function f(spec::SpeechBox.spectrum; nticks = 8)
    frqs = spec.frqs
    tindx = Int.(round.(range(1, length(frqs), length = nticks)))
    xtks = (tindx,string.(round.(frqs[tindx],digits=1)))

    legend := false

    @series begin
        seriestype := :line
        xticks := xtks
        xlabel := "Frequency(Hz)"
        title := spec.title
        spec.components
    end
end
