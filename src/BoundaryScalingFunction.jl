struct BoundaryScalingFunction <: AbstractScalingFunction
    values::OffsetArrays.OffsetVector{Float64, Vector{Float64}}
    vanishing_moments::Int64
    index::Int64
    scale::Int64
    side::Sides

    function BoundaryScalingFunction(values, p, index, scale, side)
        if p < 0
			throw(DomainError(p, "Not a valid number of vanishing moments"))
        end

        if scale < 0
			throw(DomainError(scale, "Not a valid scale"))
        end

        if !(0 <= index < p)
			throw(DomainError(index, "Index should be between 0 and " * p))
        end

        new(values, p, index, scale, side)
    end
end


(phi::BoundaryScalingFunction)(x::DyadicRational) = phi[x]
index(phi::BoundaryScalingFunction) = phi.index
support_boundaries(phi::BoundaryScalingFunction) = 0, vanishing_moments(phi) + index(phi)


struct BoundaryScalingFunctions
    #= values::OffsetArrays.OffsetVector{OffsetArrays.OffsetVector{Float64, Vector{Float64}}} =#
    functions::Vector{BoundaryScalingFunction}
    filters::BoundaryFilters
    phi::InteriorScalingFunction
    # TODO: vanishing_moments is not needed when we have filter
    vanishing_moments::Int64
    scale::Int64
    side::Sides

    function BoundaryScalingFunctions(functions, filter, phi, p, scale, sides)
        if p < 0
			throw(DomainError(p, "Not a valid number of vanishing moments"))
        end

        if scale < 0
			throw(DomainError(scale, "Not a valid scale"))
        end

        if support_boundaries(phi) != (-p + 1, p)
            ErrorException("Interior scaling function has wrong domain")
        end

        new(functions, filter, phi, p, scale, sides)
    end
end


functions(Phi::BoundaryScalingFunctions) = Phi.functions
filters(Phi::BoundaryScalingFunctions) = Phi.filters
vanishing_moments(Phi::BoundaryScalingFunctions) = Phi.vanishing_moments
side(Phi::BoundaryScalingFunctions) = Phi.side
scale(Phi::BoundaryScalingFunctions) = Phi.scale
Base.getindex(Phi::BoundaryScalingFunctions, idx::Integer) = functions(Phi)[idx + 1]


function initialize_boundary_scaling_functions(b::BoundaryFilters, phi::InteriorScalingFunction)
    p = vanishing_moments(b)

    Y = Vector{BoundaryScalingFunction}(undef, p)

    for k in 0:p - 1
        # TODO: Handle both sides
        # TODO: Handle different scales
        Y[k + 1] = BoundaryScalingFunction(OffsetArrays.OffsetVector{Float64}(undef, 0:p + k), p, k, 0, side(b))
    end

    BoundaryScalingFunctions(Y, b, phi, p, 0, side(b))
end


function support_union(Phi::BoundaryScalingFunctions)
    R = scale(Phi)
    p = vanishing_moments(Phi)

    # TODO: Right boundary function
    DyadicRational.(0:2^R*(2*p - 1), R)
end


#= function supports(Phi::BoundaryScalingFunctions)#::Vector{UnitRange} =#
#=     p = vanishing_moments(Phi) =#

#=     return map(x -> DyadicRational.(axes(x)[1].indices, scale(Phi)), values(Phi)) =#

#=     if side(Phi) == LEFT =#
#=         #1= return [0:p + k for k in 0:p - 1] =1# =#
#=         return [0:p + k for k in 0:p - 1] =#
#=     elseif side(Phi) == RIGHT =#
#=         #1= return [-(p + k):0 for k in 0:p - 1] =1# =#
#=         return [-(p + k):0 for k in 0:p - 1] =#
#=     end =#
#= end =#


function boundary_scaling_functions(b::BoundaryFilters, phi::InteriorScalingFunction)
    Phi = initialize_boundary_scaling_functions(b, phi)

    p = vanishing_moments(Phi)

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
                phi_val += sqrt2 * filters(Phi)[k][m] * Phi.phi(2*x - m)
            end

            Phi[k][x] = phi_val
        end
    end

    for k in 0:p - 1
        # TODO: getindex/setindex! with integers?
        Phi[k][DyadicRational(0, 0)] = 0.0
        #= Phi[k][DyadicRational(0, 0)] = 1 - sum(values(Phi[k])) =#
    end

    return Phi
end

