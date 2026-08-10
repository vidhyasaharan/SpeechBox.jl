#Peak finding in sequences (used by the pitch estimators). Generic frequency/sampling
#utilities that used to live here have moved to TimeFrequencyAnalysis.jl.

#Generate first order difference of a sequence y[i] = x[i+1] - x[i] (output sequence length is 1 less than input sequence length)
function Δ(x::AbstractVector)
    len = length(x)
    Δx = zeros(typeof(x[1]),len-1)
    if(len>1)
        @views for i=1:len-1
            Δx[i] = x[i+1]-x[i]
        end
    end
    return Δx
end


#Find local peaks/maximas in a sequence, keeping peaks at least min_dist samples apart
function findpeaks(x::Vector; min_dist::Int = 2)
    ind = Int[]
    mag = eltype(x)[]
    if x[1]>x[2]
        push!(ind,1)
        push!(mag,x[1])
    end
    for i=2:length(x)-1
        if(x[i-1]<x[i]>x[i+1])
            push!(ind,i)
            push!(mag,x[i])
        end
    end
    if(x[end]>x[end-1])
        push!(ind,length(x))
        push!(mag,x[end])
    end
    if(length(ind)>1)
        while(minimum(Δ(ind))<min_dist)
            remove_nearest_peak!(ind,mag)
        end
    end
    return ind,mag
end


#Local peaks/maximas sorted by magnitude (optionally only the largest num_peaks)
function findpeaks_sorted(x::Vector; min_dist::Int = 2, num_peaks::Int = 0)
    ind, mag = findpeaks(x;min_dist)
    si = sortperm(mag, rev=true)
    if((num_peaks>0)&&(num_peaks<length(mag)))
        si = si[1:num_peaks]
    end
    return ind[si], mag[si]
end


#Support function for findpeaks() - removes the smaller of the two closest peaks in a set of local peaks
function remove_nearest_peak!(ind::Vector,mag::Vector)
    npks = length(ind)
    if(npks>1)
        dist = Δ(ind)
        m_i = argmin(dist)
        if(mag[m_i+1]<mag[m_i])
            m_i += 1
        end
        deleteat!(ind,m_i)
        deleteat!(mag,m_i)
    end
end
