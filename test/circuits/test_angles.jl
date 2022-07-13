@testset "Angles" begin
	
	@test Angle(1//2) == Angle(2//4)
	@test Angle(1//4) < Angle(1//2)
	@test Angle(1//2) == (Angle(1//4) + Angle(1//4))
	@test Angle(1//2) == (Angle(3//4) - Angle(1//4))
	@test Angle(1//1) == Angle(1//2)*2.0
	@test HiQuER._to_latex_raw(Angle(-3//2)) == "\\frac{-3\\pi}{2}"

	@test FloatAngle(0.5).value ≈ Angle(1//2).value

end