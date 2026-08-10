
#Levinson-Durbin Recursion to estimate LPC from
function levinson_durbin(rxx::AbstractVector{T}) where {T<:AbstractFloat}
    p = length(rxx) - 1
    α = Vector{T}(undef,p)
    k = Vector{T}(undef,p)
    E = Vector{T}(undef,p+1)
    @views E[1] = rxx[1]
    @views k[1] = rxx[2]/rxx[1]
    @views α[p] = k[1]
    update_E!(E,k,1)
    for i ∈ 2:p
        update_k!(k, α, rxx, E, i)
        update_α!(α, k, i)
        update_E!(E, k, i)
    end
    return α, E[end], k
end


function update_E!(E::AbstractVector{T}, k::AbstractVector{T}, iter::Int) where {T<:AbstractFloat}
    @views E[iter+1] = (1-(k[iter]^2))*E[iter]
end

function update_k!(k::AbstractVector{T}, α::AbstractVector{T}, rxx::AbstractVector{T}, E::AbstractVector{T}, iter::Int) where {T<:AbstractFloat}
    @views a = α[end-iter+2:end]
    @views k[iter] = (rxx[iter+1] - dot(a,rxx[2:iter]))/E[iter]
end

function update_α!(α::AbstractVector{T}, k::AbstractVector{T}, iter::Int) where {T<:AbstractFloat}
    @views α[end-iter+1] = k[iter]
    a = α[end-iter+2:end]
    @views for i ∈ 1:iter-1
        α[end-i+1] -= k[iter]*a[i]
    end
end
