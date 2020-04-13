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


#Plot recipe for plotting time domain waveform
@recipe function f(sig::speech_waveform)
    x = sig.x
    fs = sig.fs
    t = (0:length(x)-1)/fs

    legend := false

    @series begin
        seriestype := :line
        xlabel := "Time (secs)"
        t,x
    end
end
