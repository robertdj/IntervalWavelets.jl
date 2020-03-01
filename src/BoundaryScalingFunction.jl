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


struct BoundaryScalingFunctions{T <: AbstractBoundaryScalingFunction}
    functions::Vector{T}
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
interior(Phi::BoundaryScalingFunctions) = Phi.phi
side(Phi::BoundaryScalingFunctions) = Phi.side
resolution(Phi::BoundaryScalingFunctions) = Phi.scale
Base.getindex(Phi::BoundaryScalingFunctions, idx::Integer) = functions(Phi)[idx + 1]

Base.iterate(Phi::BoundaryScalingFunctions) = Phi[0], 0
Base.iterate(Phi::BoundaryScalingFunctions, state) = state > length(Phi) ? nothing : (Phi[state], state + 1)
Base.length(Phi::BoundaryScalingFunctions) = length(functions(Phi)) - 1


function initialize_boundary_scaling_functions(b::BoundaryFilters, phi::InteriorScalingFunction, R::Integer)
    p = vanishing_moments(b)

    if side(b) == LEFT
    Y = [LeftScalingFunction(
            OffsetArrays.OffsetVector{Float64}(undef, 0:2^R*(p + k)), p, k, R
        ) for k = 0:p - 1]
    elseif side(b) == RIGHT
    Y = [RightScalingFunction(
            OffsetArrays.OffsetVector{Float64}(undef, -2^R*(p + k):0), p, k, R
        ) for k = 0:p - 1]
    end

    BoundaryScalingFunctions(Y, b, phi, p, R, side(b))
end


function support_union(Phi::BoundaryScalingFunctions)
    p = vanishing_moments(Phi)
    support(Phi[p - 1])
end


function sorted_support_union(Phi::BoundaryScalingFunctions{LeftScalingFunction})
    support_union(Phi) |> reverse
end


function sorted_support_union(Phi::BoundaryScalingFunctions{RightScalingFunction})
    support_union(Phi)
end


function supports(Phi::BoundaryScalingFunctions)::Vector{Vector{DyadicRational}}
    p = vanishing_moments(Phi)
    [support(Phi[k]) for k = 0:p - 1]
end


function boundary_scaling_functions(b::BoundaryFilters, phi::InteriorScalingFunction)
    Phi = initialize_boundary_scaling_functions(b, phi, 0)

    p = vanishing_moments(Phi)

    support_values = sorted_support_union(Phi)

    for x in support_values
        for k in p-1:-1:0
            if x ∉ support(Phi[k])
                continue
            end

            phi_val = 0.0

            # Boundary contribution
            for l in 0:p - 1
                phi_val += filters(Phi)[k][l] * Phi[l](2*x)
            end
            
            # Interior contribution
            for m in p:p + 2k
                phi_val += filters(Phi)[k][m] * phi(2*x - interior_translation(m, side(b)))
            end

            Phi[k][x] = sqrt2 * phi_val
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


@inline function interior_translation(m::Integer, side::Sides)
    if side == LEFT
        return m
    elseif side == RIGHT
        return -m - 1
    end
end


function boundary_scaling_functions(side::Sides, p::Integer, R::Integer)
    h = interior_filter(p)
    phi = interior_scaling_function(h, R)

    b = boundary_filters(p, side)

    boundary_scaling_functions(b, phi, R)
end


function boundary_scaling_functions(b::BoundaryFilters, phi::InteriorScalingFunction, R::Integer)
    if R < 0
        throw(DomainError(R, "Resolution must be non-negative"))
    end

    Phi = boundary_scaling_functions(b::BoundaryFilters, phi::InteriorScalingFunction)

    for _ = 1:R
        Phi = increase_resolution(Phi)
    end

    compute_edge_value!(Phi)

    return Phi
end


function increase_resolution(Phi::BoundaryScalingFunctions)
    phi = interior(Phi)
    Phi2 = initialize_boundary_scaling_functions(filters(Phi), phi, resolution(Phi) + 1)

    p = vanishing_moments(Phi)

    support_values = support_union(Phi2)

    for (index, x) in enumerate(support_values)
        for k in p-1:-1:0
            if x ∉ support(Phi2[k])
                continue
            end

            if isodd(index)
                Phi2[k][x] = Phi[k][x]
                continue
            end

            phi_val = 0.0

            # Boundary contribution
            for l in 0:p - 1
                phi_val += filters(Phi)[k][l] * Phi[l](2*x)
            end
            
            # Interior contribution
            for m in p:p + 2k
                phi_val += filters(Phi)[k][m] * phi(2*x - interior_translation(m, side(Phi)))
            end

            Phi2[k][x] = sqrt2 * phi_val
        end
    end

    return Phi2
end


function compute_edge_value!(Phi)
    p = vanishing_moments(Phi)
    R = resolution(Phi)
    if R <= ceil(log2(p))
        return Phi
    end

    x = sign(side(Phi)) * DyadicRational.(1:p, R)
    representation_matrix = [Phi[j](x[i]) for i = 1:p, j = 0:p - 1]
    basis_coeff = representation_matrix \ ones(p)

    filter_matrix = [sqrt2 * filters(Phi)[i][j] for i = 0:p - 1, j = 0:p - 1]
    scaled_function_values = eigval1(filter_matrix)

    scaling_factor = LinearAlgebra.dot(basis_coeff, scaled_function_values)

    function_values = scaled_function_values ./ scaling_factor

    for (index, phi) in enumerate(Phi.functions)
        Phi.functions[index][DyadicRational(0,0)] = function_values[index]
    end

    return Phi
end


@recipe function f(phi::BoundaryScalingFunctions)
    for y in phi
        @series begin
            collect(y)
        end
    end
end

