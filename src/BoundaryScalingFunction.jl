abstract type AbstractBoundaryScalingFunction <: AbstractScalingFunction end

struct LeftScalingFunction <: AbstractBoundaryScalingFunction
    values::OffsetArrays.OffsetVector{Float64, Vector{Float64}}
    vanishing_moments::Int64
    index::Int64
    scale::Int64

    #= function LeftScalingFunction(values, p, index, scale, side) =#
    #=     if p < 0 =#
			#= throw(DomainError(p, "Not a valid number of vanishing moments")) =#
    #=     end =#

    #=     if scale < 0 =#
			#= throw(DomainError(scale, "Not a valid scale")) =#
    #=     end =#

    #=     if !(0 <= index < p) =#
			#= throw(DomainError(index, "Index should be between 0 and " * p - 1)) =#
    #=     end =#

    #=     new(values, p, index, scale) =#
    #= end =#
end


struct RightScalingFunction <: AbstractBoundaryScalingFunction
    values::OffsetArrays.OffsetVector{Float64, Vector{Float64}}
    vanishing_moments::Int64
    index::Int64
    scale::Int64
end


index(phi::AbstractBoundaryScalingFunction) = phi.index
support_boundaries(phi::LeftScalingFunction) = 0, vanishing_moments(phi) + index(phi)
support_boundaries(phi::RightScalingFunction) = -(vanishing_moments(phi) + index(phi)), 0

# It is not possible to have a functor for an abstract type
function value(phi::AbstractBoundaryScalingFunction, x::DyadicRational)
    phi[x]
end
(phi::LeftScalingFunction)(x::DyadicRational) = value(phi, x)
(phi::RightScalingFunction)(x::DyadicRational) = value(phi, x)


struct BoundaryScalingFunctions
    functions::Vector{AbstractBoundaryScalingFunction}
    filters::BoundaryFilters
    phi::InteriorScalingFunction
    # TODO: vanishing_moments is not needed when we have filter
    vanishing_moments::Int64
    scale::Int64
    side::Sides

    #= function BoundaryScalingFunctions(functions, filter, phi, p, scale, sides) =#
    #=     if p < 0 =#
			#= throw(DomainError(p, "Not a valid number of vanishing moments")) =#
    #=     end =#

    #=     if scale < 0 =#
			#= throw(DomainError(scale, "Not a valid scale")) =#
    #=     end =#

    #=     if support_boundaries(phi) != (-p + 1, p) =#
    #=         ErrorException("Interior scaling function has wrong domain") =#
    #=     end =#

    #=     new(functions, filter, phi, p, scale, sides) =#
    #= end =#
end


functions(Phi::BoundaryScalingFunctions) = Phi.functions
filters(Phi::BoundaryScalingFunctions) = Phi.filters
vanishing_moments(Phi::BoundaryScalingFunctions) = Phi.vanishing_moments
side(Phi::BoundaryScalingFunctions) = Phi.side
scale(Phi::BoundaryScalingFunctions) = Phi.scale
Base.getindex(Phi::BoundaryScalingFunctions, idx::Integer) = functions(Phi)[idx + 1]


function initialize_boundary_scaling_functions(b::BoundaryFilters, phi::InteriorScalingFunction)
    p = vanishing_moments(b)

    if side(b) == LEFT
    Y = [LeftScalingFunction(
            OffsetArrays.OffsetVector{Float64}(undef, 0:p + k), p, k, 0
        ) for k = 0:p - 1]
    elseif side(b) == RIGHT
    Y = [RightScalingFunction(
            OffsetArrays.OffsetVector{Float64}(undef, -(p + k):0), p, k, 0
        ) for k = 0:p - 1]
    end

    BoundaryScalingFunctions(Y, b, phi, p, 0, side(b))
end


function support_union(Phi::BoundaryScalingFunctions)
    p = vanishing_moments(Phi)
    support(Phi[p - 1])
end


function supports(Phi::BoundaryScalingFunctions)::Vector{Vector{DyadicRational}}
    p = vanishing_moments(Phi)
    [support(Phi[k]) for k = 0:p - 1]
end


function boundary_scaling_functions(b::BoundaryFilters, phi::InteriorScalingFunction)
    Phi = initialize_boundary_scaling_functions(b, phi)

    p = vanishing_moments(Phi)

    # TODO: Fix for right side
    support_values = support_union(Phi) |> reverse

    for x in support_values
        for k in p-1:-1:0
            if x ∉ support(Phi[k])
                continue
            end

            phi_val = 0.0

            # Boundary contribution
            for l in 0:p - 1
                phi_val += sqrt2 * filters(Phi)[k][l] * Phi[l](2*x)
            end
            
            # Interior contribution
            for m in p:p + 2k
                phi_val += sqrt2 * filters(Phi)[k][m] * phi(2*x - m)
            end

            Phi[k][x] = phi_val
        end
    end

    for k in 0:p - 1
        # TODO: getindex/setindex! with integers?
        # TODO: I think it would make more sense to set value to be
        # "missing". This requires a Union in BoundaryScalingFunction
        Phi[k][DyadicRational(0, 0)] = 0.0
    end

    return Phi
end

