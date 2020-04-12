@testset "framed_signal" begin
    frames = framed_signal(x,fs) #for Array{Float} input
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
