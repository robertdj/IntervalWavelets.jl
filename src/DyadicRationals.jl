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

function Base.show(io::IO, dr::DyadicRational)
    print(io, numerator(dr), "/2^", scale(dr))
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

