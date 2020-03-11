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
scale(dr::DyadicRational) = dr.R

Base.:+(dr::DyadicRational, k::Integer) = DyadicRational(numerator(dr) + (k << scale(dr)), scale(dr))
Base.:-(dr::DyadicRational, k::Integer) = DyadicRational(numerator(dr) - (k << scale(dr)), scale(dr))
Base.:*(a::Integer, dr::DyadicRational) = DyadicRational(a * numerator(dr), scale(dr))


function Base.show(io::IO, dr::DyadicRational)
    if scale(dr) == 0
        print(io, numerator(dr))
    elseif scale(dr) == 1
        print(io, numerator(dr), "/2")
    else
        print(io, numerator(dr), "/2^", scale(dr))
    end
end


function Base.Integer(dr::DyadicRational)
    if scale(dr) != 0
        throw(InexactError(:Integer, Int64, dr))
    end

    numerator(dr)
end


function all_dyadic_rationals(left::Integer, right::Integer, R::Integer)
    if R < 0
        throw(DomainError(R, "Scale must be positive"))
    end

    DyadicRational.(left*2^R:right*2^R, R)
end


function AbstractFloat(dr::DyadicRational)
    numerator(dr) / 2.0^scale(dr)
end

