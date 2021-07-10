#Itakura Saito Distortion

#Convert inverse filter coefficients to autocorrelation coefficients
function lpcar2ra(a::Vector{Float})
    na = 1/dotavx(a)
    p = length(a)
    b = Vector{Float}(undef,p)
    b[1] = 1
    for i ∈ 2:p
        b[i] = na*dotavx(a[1:p-i+1],a[i:p])
    end
    return b
end