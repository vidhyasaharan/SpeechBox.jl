#Plot recipes for the pitch overlays. Recipes for the generic containers
#(signal/spectrum/timefreq) live in TimeFrequencyAnalysis.jl's RecipesBase package
#extension; generate_ticks is imported from TimeFrequencyAnalysis in the module file.

#Plot recipe for plotting pitch: spectrogram heatmap with the pitch track overlaid
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
