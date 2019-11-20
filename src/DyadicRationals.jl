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
    print(io, numerator(dr), "/2^", scale(dr), "\n")
end

