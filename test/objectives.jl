@testset "Objectives" begin
    
@testset "LinearNonnegative" begin
    @test_throws ArgumentError LinearNonnegative(-ones(2))
    int_obj = LinearNonnegative([1, 1])
    @test int_obj.c isa Vector{Float64}
    obj = LinearNonnegative(ones(2))
    
    @test iszero(CR.f(obj, 2*ones(2)))
    @test isinf(CR.f(obj, 0.5*ones(2)))

    x = ones(2)
    CR.grad!(x, obj, 2*ones(2))
    @test iszero(x) 
    CR.grad!(x, obj, 0.5*ones(2))
    @test all(isinf.(x))
end

@testset "BasketLiquidation" begin
    @test_throws ArgumentError BasketLiquidation(0, [0.0, 1.0])
    obj = BasketLiquidation(1, [0, 1])
    @test obj.Δin isa Vector{Float64}
    
    @test CR.f(obj, [2, 3]) == 3
    @test isinf(CR.f(obj, 0.5*ones(2)))

    x = ones(2)
    CR.grad!(x, obj, 2*ones(2))
    @test iszero(x - [0, 1])
    CR.grad!(x, obj, 0.5*ones(2))
    @test all(isinf.(x))

end

end