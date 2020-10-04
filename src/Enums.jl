@enum Sides begin
    LEFT
    RIGHT
end


sign(side::Sides) = side == LEFT ? 1 : -1

