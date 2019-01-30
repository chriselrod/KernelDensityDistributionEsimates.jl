module KernelDensityDistributionEsimates

using   KernelDensity,
        Interpolations,
        Distributions,
        StatsBase, Statistics,
        Gadfly

export KDE, kdemean

struct KDE{T,ITP <: ScaledInterpolation} <: ContinuousUnivariateDistribution
    x::StepRangeLen{T,Base.TwicePrecision{T},Base.TwicePrecision{T}}
    density::Vector{T}
    cumulative_density::Vector{T}
    pdf::ITP
    cdf::ITP
end
function KDE(kde::UnivariateKDE)
    x = kde.x
    # if length(x) == 1
    #     @show x
    #     return KDE(x, kde.density, [1.0], [1.0], [1.0])
    # end
    density = kde.density
    cumulative_density = cumsum(density)
    cumulative_density ./= cumulative_density[end]

    pdf = Interpolations.scale(interpolate(density, BSpline(Quadratic(Line(OnGrid())))), x)
    cdf = Interpolations.scale(interpolate(cumulative_density, BSpline(Quadratic(Line(OnGrid())))), x)

    KDE(x, density, cumulative_density, pdf, cdf)
end
KDE(distances::AbstractVector) = KDE(kde(distances))
KDE(distances::AbstractVector, weights::AbstractVector) = KDE(kde(distances, weights = weights))

kdemean(kde::KDE) = mean(kde)
function kdemean(kde::UnivariateKDE)
    if length(kde.x) == 1
        return kde.x[1]
    else
        return mean(KDE(kde))
    end
end

# Parameters are summarized by the full (x, density) set
StatsBase.params(kde::KDE) = (kde.x, kde.density)
Statistics.mean(kde::KDE{T}) where T = kde.x' * kde.density * T(kde.x.step)
function Statistics.var(kde::KDE{T}) where T
    μ = zero(T)
    σ² = zero(T)
    x = kde.x; density = kde.density
    @inbounds @fastmath for i ∈ eachindex(x)
        xd = x[i] * density[i]
        μ += xd
        σ² += x[i] * xd
    end
    σ² * T(kde.x.step) - (μ * T(kde.x.step))^2
    # σ², μ, T(kde.x.step)
end
Statistics.std(kde::KDE) = sqrt(var(kde))
function Statistics.median(kde::KDE)
    kde.x[findfirst(p -> p > 0.5, kde.cumulative_density)]
end
function StatsBase.mode(kde::KDE)
    kde.x[argmax(kde.density)]
end
function StatsBase.entropy(kde::KDE{T}) where T
    out = zero(T)
    density = kde.density
    @inbounds for i ∈ eachindex(density)
        out -= density[i] * log(density[i])
    end
    out * T(kde.x.step)
end
function Distributions.pdf(kde::KDE, x::T) where T <: Real
    if (x < minimum(kde.x)) || (x > maximum(kde.x))
        zero(T)
    else
        max(zero(T), kde.pdf(x))
    end
end
Distributions.logpdf(kde::KDE, x::Real) = log(pdf(kde, x))
function Distributions.cdf(kde::KDE{T}, x::Real) where T
    if x < minimum(kde.x)
        p = zero(T)
    elseif x > maximum(kde.x)
        p = one(T)
    else
        # the min and max should ideally be unnecessary, but
        # special cases have come up in practice.
        # until we have a better way
        p = min(max(zero(T), kde.cdf(x)),one(T))
    end
    p
end

function Gadfly.layer(kde::KDE, args::Vararg{Union{Function, Gadfly.Element, Theme, Type},N} where N)
    layer(x = kde.x, y = kde.density, Geom.line, args...)
end
Gadfly.plot(kde::KDE, args::Vararg{Union{Function, Gadfly.Element, Theme, Type},N} where N) = plot(layer(kde, args...))

end # module
