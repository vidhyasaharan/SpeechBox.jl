@testset "speech_waveform" begin
    x = rand(Float,3,10000)
    println("Testing multichannel input, message about input being Matrix is expected")
    signal = speech_waveform(x,1000)
    @test typeof(signal.x) == Array{Float,1}
    @test length(signal.x) == 10000
    @test signal.fs == 1000.0
    @test typeof(signal.fs) == Float
end


@testset "framed_signal" begin
    frames = framed_signal(x,fs) #for Array{Float} input
    frames_alt = framed_signal(signal)
    @test frames == frames_alt
    @test typeof(frames) == framed_signal
    @test typeof(frames.signal.x) <: Array{Float}
    @test length(frames.signal.x) > 0
    @test typeof(frames.frame_length) <: Int
    @test typeof(frames.frame_shift) <: Int
    @test typeof(frames.num_signal_frames) <: Int
    @test typeof(frames.num_frames) <: Int
    @test frames.frame_length > 0
    @test frames.frame_shift > 0
    @test frames.num_signal_frames > 0
    @test frames.num_frames > 0
    @test frames.num_frames >= frames.num_signal_frames
    @test frames.num_signal_frames >= (length(frames.signal.x)-frames.frame_length)/frames.frame_shift
end



@testset "extract_frame" begin
    frames = framed_signal(x,fs)
    @test typeof(extract_frame(frames,1))<:Array{<:AbstractFloat,1}
    @test length(extract_frame(frames,1)) > 0
    @test length(extract_frame(frames,frames.num_frames)) > 0
    @test maximum(abs.(extract_frame(frames,frames.num_frames))) > 0
end


@testset "spectrum" begin
    x = rand(Float,220)
    fs = 220
    c = rand(Complex{Float},111)
    f = convert.(Float,collect(1:length(c)))
    sp = spectrum(speech_waveform(x,fs),c,f)
    @test typeof(sp.signal) == speech_waveform
    @test typeof(sp.components) == Array{Complex{Float},1}
    @test typeof(sp.frqs) == Array{Float,1}
    @test typeof(sp.title) <: AbstractString
end


@testset "timefreq" begin
    sw = speech_waveform(rand(Float,1000),1000)
    comps = rand(Float,5,4)
    frqs = collect(1.0:5.0)
    t = collect(0.1:0.1:0.4)
    tf = SpeechBox.timefreq(sw, comps, frqs, t, "test") #construct directly from a waveform (no frames)
    @test tf.signal == sw
    @test tf.frames === nothing
    @test tf.components == comps
    @test tf.frqs == frqs
    @test tf.time == t
    @test tf.title == "test"
end
