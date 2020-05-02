using HartreeFock: gaussian_product_center, norm2, boys, Vec3

@test gaussian_product_center(1.0, Vec3([0.0, 0.0, 0.0]), 1.5, Vec3([1.0, 1.0, 1.0])) == [0.6, 0.6, 0.6]
@test norm2(Vec3([1.0, 2.0, 3.0])) == 14.0
@test boys(0, 1.0) ≈ 0.7468241328
