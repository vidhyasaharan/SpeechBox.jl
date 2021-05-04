#Dot product using LoopVectorization (real ⋅ real)
function dotavx(a::AbstractVector{T}, b::AbstractVector{T}) where {T}
    s = zero(T)
    @avx unroll=8 for i ∈ eachindex(a,b)
        s += a[i] * b[i]
    end
    return s
end

#Dot product using LoopVectorization (real ⋅ complex)
function dotavx(a::AbstractVector{T}, cb::AbstractVector{Complex{T}}) where {T}
    re = zero(T)
    im = zero(T)
    b = reinterpret(reshape, T, cb)
    @avx for i ∈ eachindex(a)
        re += a[i] * b[1,i]
        im += a[i] * b[2,i]
    end
    return Complex(re, im)
end

#Dot product using LoopVectorization (real ⋅ complex)
function dotavx(ca::AbstractVector{Complex{T}}, b::AbstractVector{T}) where {T}
    re = zero(T)
    im = zero(T)
    a = reinterpret(reshape, T, ca)
    @avx for i ∈ eachindex(b)
        re += a[1,i] * b[i]
        im += - (a[2,i] * b[i])
    end
    return Complex(re, im)
end

#Dot product using LoopVectorization (complex ⋅ complex)
function dotavx(ca::AbstractVector{Complex{T}}, cb::AbstractVector{Complex{T}}) where {T}
    re = zero(T)
    im = zero(T)
    a = reinterpret(reshape, T, ca)
    b = reinterpret(reshape, T, cb)
    @avx for i ∈ axes(a,2) #Conjugate(a) × b
        re += (a[1,i] * b[1,i]) + (a[2,i] * b[2,i])
        im += (a[1,i] * b[2,i]) - (a[2,i] * b[1,i])
    end
    return Complex(re, im)
end

#Cross corrleation with zero padding (and using dotavx)
"""
    xcorr(x, h[, z=1])

Computes the cross correlation between `x` and `h`, with the optional `z` indicating the position of the zero index of the array `h`. The output is of the same length as `x` and the cross correlation is computed with zero padding.
"""
function xcorr(x::Vector{Float},h::Vector{Float},z::Int=1)
    padded_x = zeros(length(x)+length(h)-1)
    padded_x[z:z+length(x)-1] = x
    y = Vector{Float}(undef,length(x))
    @inbounds for i ∈ eachindex(y)
        y[i] = dotavx(padded_x[i:i+length(h)-1],h)
    end
    return y
end


#Matrix multiplcation (in place for resultnant matrix)
function A_mul_B!(C::AbstractMatrix{T}, A::AbstractMatrix{T}, B::AbstractMatrix{T}) where {T}
    @avx for n ∈ indices((C,B), 2), m ∈ indices((C,A), 1)
        Cmn = zero(eltype(C))
        for k ∈ indices((A,B), (2,1))
            Cmn += A[m,k] * B[k,n]
        end
        C[m,n] = Cmn
    end
end

function A_mul_B!(cC::AbstractMatrix{Complex{T}}, cA::AbstractMatrix{Complex{T}}, B::AbstractMatrix{T}) where {T}
    A = reinterpret(reshape, T, cA)
    @avx for n ∈ indices((cC,B), 2), m ∈ indices((cC,cA), 1)
        re = zero(T)
        im = zero(T)
        for k ∈ indices((cA,B), (2,1))
            re += A[1,m,k] * B[k,n]
            im += A[2,m,k] * B[k,n]
        end
        cC[m,n] = Complex(re,im)
    end
end