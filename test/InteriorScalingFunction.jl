using IntervalWavelets
using Test


@testset "Interior scaling functions" begin
    @testset "Dyadic dilation matrix has the expected form" begin
        h1 = interior_filter(3, :min)
        H1 = IntervalWavelets.dyadic_dilation_matrix(h1)

        @test size(H1) == (length(h1) - 2, length(h1) - 2)

        @test H1 == sqrt(2) * 
        [ h1[1] h1[0]    0     0  ;
          h1[3] h1[2] h1[1] h1[0] ;
          h1[5] h1[4] h1[3] h1[2] ;
             0     0  h1[5] h1[4] ]


        h2 = interior_filter(3, :symmlet)
        H2 = IntervalWavelets.dyadic_dilation_matrix(h2)

        @test H2 == sqrt(2) * 
        [ h2[-1] h2[-2]      0      0  ;
          h2[1]  h2[0]   h2[-1] h2[-2] ;
          h2[3]  h2[2]   h2[1]  h2[0] ;
             0      0    h2[3]  h2[2] ]
    end

end

