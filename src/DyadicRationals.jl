struct DyadicRational <: AbstractFloat
    numerator::Int64
    R::Int64

    function DyadicRational(numerator, R)
        if R < 0
            throw(DomainError(R, "Scale must be positive"))
        end

        new(numerator, R)
    end
end

numerator(dr::DyadicRational) = dr.numerator
scale(dr::DyadicRational) = dr.R

Base.:+(dr::DyadicRational, k::Integer) = DyadicRational(numerator(dr) + 2^scale(dr) * k, scale(dr))
Base.:-(dr::DyadicRational, k::Integer) = DyadicRational(numerator(dr) - 2^scale(dr) * k, scale(dr))
Base.:*(a::Integer, dr::DyadicRational) = DyadicRational(a * numerator(dr), scale(dr))


function Base.show(io::IO, dr::DyadicRational)
    print(io, numerator(dr), "/2^", scale(dr))
end


function Base.Integer(dr::DyadicRational)
    if scale(dr) != 0
        throw(InexactError(:Integer, Int64, dr))
    end

    numerator(dr)
end


"""
Reduce the fraction in a `DyadicRational`.
"""
function reduce(dr::DyadicRational)
    largest_power_of_two_dividing_numerator = trailing_zeros(numerator(dr))

    common_power_of_two = min(scale(dr), largest_power_of_two_dividing_numerator)

    DyadicRational(numerator(dr) >> common_power_of_two, scale(dr) - common_power_of_two)
end


function all_dyadic_rationals(left::Integer, right::Integer, R::Integer)
    if R < 0
        throw(DomainError(R, "Scale must be positive"))
    end

    DyadicRational.(left*2^R:right*2^R, R)
end


function AbstractFloat(dr::DyadicRational)
    # TODO: Can this be done with >> ?
    numerator(dr) / 2.0^scale(dr)
end


struct DyadicRationalVector
    numerator::AbstractVector{Int64}
    R::Int64

    function DyadicRationalVector(numerator, R)
        if R < 0
            throw(DomainError(R, "Scale must be positive"))
        end

        new(numerator, R)
    end
end

Base.getindex(dr::DyadicRationalVector, key) = DyadicRational(numerator(dr)[key], scale(dr))
Base.lastindex(dr::DyadicRationalVector) = lastindex(numerator(dr))

numerator(dr::DyadicRationalVector) = dr.numerator
scale(dr::DyadicRationalVector) = dr.R

function Base.show(io::IO, dr::DyadicRationalVector)
    print(io, dr[1], ", ..., ", dr[end])
end

