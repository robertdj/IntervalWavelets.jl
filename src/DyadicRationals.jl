struct DyadicRational <: AbstractFloat
    numerator::Int64
    R::Int64

    function DyadicRational(numerator, R)
        if R < 0
            throw(DomainError(R, "Scale must be non-negative"))
        end

        largest_power_of_two_dividing_numerator = trailing_zeros(numerator)
        common_power_of_two = min(R, largest_power_of_two_dividing_numerator)

        new(numerator >> common_power_of_two, R - common_power_of_two)
    end
end

numerator(dr::DyadicRational) = dr.numerator
resolution(dr::DyadicRational) = dr.R

Base.:+(dr::DyadicRational, k::Integer) = DyadicRational(numerator(dr) + (k << resolution(dr)), resolution(dr))
Base.:-(dr::DyadicRational, k::Integer) = DyadicRational(numerator(dr) - (k << resolution(dr)), resolution(dr))
Base.:*(a::Integer, dr::DyadicRational) = DyadicRational(a * numerator(dr), resolution(dr))


function Base.show(io::IO, dr::DyadicRational)
    if resolution(dr) == 0
        print(io, numerator(dr))
    elseif resolution(dr) == 1
        print(io, numerator(dr), "/2")
    else
        print(io, numerator(dr), "/2^", resolution(dr))
    end
end


#= Base.convert(::Type{DyadicRational}, x::Integer) = Integer(x) =#
#= Base.promote_rule(::Type{DyadicRational}, ::Type{Integer}) = DyadicRational =#

#= Base.:<(x::DyadicRational, y::DyadicRational) = numerator(x) << resolution(y) < numerator(y) << resolution(x) =#
#= Base.:<=(x::DyadicRational, y::DyadicRational) = numerator(x) << resolution(y) <= numerator(y) << resolution(x) =#

function Base.Integer(dr::DyadicRational)
    if resolution(dr) != 0
        throw(InexactError(:Integer, Int64, dr))
    end

    numerator(dr)
end


function all_dyadic_rationals(left::Integer, right::Integer, R::Integer)
    if R < 0
        throw(DomainError(R, "Scale must be non-negative"))
    end

    DyadicRational.(left*2^R:right*2^R, R)
end


function AbstractFloat(dr::DyadicRational)
    numerator(dr) / 2.0^resolution(dr)
end

