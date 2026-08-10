function generate_ticks(label_values::AbstractVector{<:Real}, nticks::Int)
    indx = Int.(round.(range(1, length(label_values), length = nticks)))
    tks = (indx,string.(round.(label_values[indx],digits=1)))
    return tks
end

#Plot recipe for plotting pitch
@recipe function f(p::pitch_timefreq; nxticks = 8, nyticks = 8)
    msp = p.msp
    frqs = msp.frqs
    time = msp.time
    xtks = generate_ticks(time,nxticks)
    ytks = generate_ticks(frqs,nyticks)

    @series begin
        seriestype := :heatmap
        xticks := xtks
        yticks := ytks
        msp.components
    end

    @series begin
        seriestype := :line
        xticks := xtks
        seriescolor := :green
        legend := false
        xguide := "Time (sec)"
        yguide := "Frequency (Hz)"
        if(~isnothing(msp.title))
            title := msp.title
        end
        p.pindx
    end
end



#Plot recipe for plotting spectrogram
@recipe function f(msp::SpeechBox.timefreq; nxticks = 8, nyticks = 8)
    frqs = msp.frqs
    time = msp.time
    xtks = generate_ticks(time,nxticks)
    ytks = generate_ticks(frqs,nyticks)

    @series begin
        seriestype := :heatmap
        xticks := xtks
        yticks := ytks
        xguide := "Time (sec)"
        yguide := "Frequency (Hz)"
        if(~isnothing(msp.title))
            title := msp.title
        end
        msp.components
    end
end


#Plot recipe for plotting spectra
@recipe function f(spec::SpeechBox.spectrum; nticks = 8)
    frqs = spec.frqs
    xtks = generate_ticks(frqs,nticks)

    legend := false

    @series begin
        seriestype := :line
        xticks := xtks
        xguide := "Frequency(Hz)"
        if(~isnothing(spec.title))
            title := spec.title
        end
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
        xguide := "Time (secs)"
        t,x
    end
end
