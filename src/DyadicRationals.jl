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

numerator(x::DyadicRational) = x.numerator
resolution(x::DyadicRational) = x.R


function Base.show(io::IO, x::DyadicRational)
    if resolution(x) == 0
        print(io, numerator(x))
    elseif resolution(x) == 1
        print(io, numerator(x), "/2")
    else
        print(io, numerator(x), "/2^", resolution(x))
    end
end


Base.convert(::Type{DyadicRational}, x::Integer) = DyadicRational(x, 0)
Base.promote_rule(::Type{DyadicRational}, ::Type{T}) where {T <: Integer} = DyadicRational

Base.:*(x::DyadicRational, y::DyadicRational) = DyadicRational(numerator(x) * numerator(y), resolution(x) + resolution(y))
Base.:*(x::DyadicRational, k::Integer) = DyadicRational(k * numerator(x), resolution(x))
Base.:*(k::Integer, x::DyadicRational) = x * k

for op in [:+, :-]
    @eval begin
        Base.$op(k::Integer, x::DyadicRational) = $op(x, k)

        function Base.$op(x::DyadicRational, k::Integer)
            new_numerator = $op(numerator(x), k << resolution(x))
            DyadicRational(new_numerator, resolution(x))
        end

        function Base.$op(a::DyadicRational, b::DyadicRational)
            common_numerator = $op(numerator(a) << resolution(b), numerator(b) << resolution(a))
            common_resolution = resolution(a) + resolution(b)

            DyadicRational(common_numerator, common_resolution)
        end
    end
end

for op in [:<, :<=]
    @eval begin
        Base.$op(x::DyadicRational, y::DyadicRational) = $op(numerator(x) << resolution(y), numerator(y) << resolution(x))
    end
end


function Base.Integer(x::DyadicRational)
    if resolution(x) != 0
        throw(InexactError(:Integer, Int64, x))
    end

    numerator(x)
end


function all_dyadic_rationals(left, right, R)
    all_dyadic_rationals(convert(DyadicRational, left), convert(DyadicRational, right), R)
end


function all_dyadic_rationals(left::DyadicRational, right::DyadicRational, R::Integer)
    if R < resolution(left) || R < resolution(right)
        throw(DomainError(R, "Resolution must be at least that of the endpoints"))
    end

    common_resolution = max(resolution(left), resolution(right), R)
    left_numerator = numerator(left) << (common_resolution - resolution(left))
    right_numerator = numerator(right) << (common_resolution - resolution(right))

    DyadicRational.(left_numerator:right_numerator, R)
end


function AbstractFloat(x::DyadicRational)
    numerator(x) / 2.0^resolution(x)
end

