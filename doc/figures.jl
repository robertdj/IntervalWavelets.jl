using IntervalWavelets

# Interior
x, y = DaubScaling(4,10);
C = ifilter(4,true)
z = DaubScaling(C,10);
plot(x, y, x, z, "r")
savefig("interior.pdf")

# Left boundary
x, Y = DaubScaling(4, 'L', 10);
Y[1,:] = Y[2,:]
plot(x, Y)
legend( ["k=0"; "k=1"; "k=2"; "k=3"], (0.8,0.9) )
savefig("left.pdf")

# Right boundary
x, Y = DaubScaling(4, 'R', 10);
Y[end,:] = Y[end-1,:]
plot(x, Y)
legend( ["k=0"; "k=1"; "k=2"; "k=3"] )
savefig("right.pdf")

