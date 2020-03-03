# TODO: How to use this enum in a good way? In its own module?
@enum Sides begin
    LEFT
    RIGHT
end


sign(side::Sides) = side == LEFT ? 1 : -1

