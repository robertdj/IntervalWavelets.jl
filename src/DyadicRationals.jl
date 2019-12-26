struct DyadicRational
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


function reduce(dr::DyadicRational)
    common_factor = gcd(numerator(dr), 2^scale(dr))

    if ispow2(common_factor)
        # TODO: prevpow instead of log2?
        return DyadicRational(numerator(dr) / common_factor, scale(dr) - Int(log2(common_factor)))
    end

    return dr
end


function all_dyadic_rationals(left::Integer, right::Integer, R::Integer)
    if R < 0
        throw(DomainError(R, "Scale must be positive"))
    end

    DyadicRational.(left*2^R:right*2^R, R)
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

