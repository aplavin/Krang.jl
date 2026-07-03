@testset "Ray tracer Functions" begin

    @testset "Emission radius" begin
        a = 0.99
        met = Krang.Kerr(a)
        α = 10.0
        β = 1.0
        θo = π / 2
        @testset "$pixtype" for (pixtype, pix) in [
            ("Intensity Pixel", Krang.IntensityPixel(met, α, β, θo)),
            (
                "Cached Slow Light Intensity Pixel",
                Krang.SlowLightIntensityPixel(met, α, β, θo),
            ),
        ]

            @test !emission_radius(pix, π / 2, true, 2)[5]

            @testset "Case 2" begin
                α = 10.0
                β = 10.0
                θo = π / 4
                @testset "$pixtype" for (pixtype, pix) in [
                    ("Intensity Pixel", Krang.IntensityPixel(met, α, β, θo)),
                    (
                        "Cached Slow Light Intensity Pixel",
                        Krang.SlowLightIntensityPixel(met, α, β, θo),
                    ),
                ]

                    ηcase1 = η(met, α, β, θo)
                    λcase1 = λ(met, α, θo)
                    roots = get_radial_roots(met, ηcase1, λcase1)
                    _, _, _, root = roots
                    @test sum(Krang._isreal2.(roots)) == 4
                    rs = 1.1 * real(root)
                    τ1 = Ir(pix, true, rs)[1]
                    @test Krang.emission_radius(pix, τ1)[1] / rs ≈ 1 atol = 1e-5
                end
            end

            @testset "Case 3" begin
                α = 1.0
                β = 1.0
                θo = π / 4

                @testset "$pixtype" for (pixtype, pix) in [
                    ("Intensity Pixel", Krang.IntensityPixel(met, α, β, θo)),
                    (
                        "Cached Slow Light Intensity Pixel",
                        Krang.SlowLightIntensityPixel(met, α, β, θo),
                    ),
                ]

                    ηcase3 = η(met, α, β, θo)
                    λcase3 = λ(met, α, θo)
                    roots = get_radial_roots(met, ηcase3, λcase3)
                    @test sum(Krang._isreal2.(roots)) == 2
                    rs = 1.1horizon(met)
                    τ3 = Ir(pix, true, rs)[1]
                    @test Krang.emission_radius(pix, τ3)[1] / rs ≈ 1 atol = 1e-5
                end
            end
            @testset "Case 4" begin
                α = 0.1
                β = 0.1
                θo = π / 4
                @testset "$pixtype" for (pixtype, pix) in [
                    ("Intensity Pixel", Krang.IntensityPixel(met, α, β, θo)),
                    (
                        "Cached Slow Light Intensity Pixel",
                        Krang.SlowLightIntensityPixel(met, α, β, θo),
                    ),
                ]

                    ηcase4 = η(met, α, β, θo)
                    λcase4 = λ(met, α, θo)
                    roots = get_radial_roots(met, ηcase4, λcase4)
                    @test sum(Krang._isreal2.(roots)) == 0
                    rs = 1.1horizon(met)
                    τ4 = Ir(pix, true, rs)[1]
                    @test Krang.emission_radius(pix, τ4)[1] / rs ≈ 1 atol = 1e-5
                end
            end
        end
    end
    @testset "Emission inclination" begin
        # Make sure Emission inclination is the inverse of Gθ
        function coordinate_point(
            pix::Krang.AbstractPixel,
            geometry::Krang.ConeGeometry{T,A},
        ) where {T,A}
            n, rmin, rmax = geometry.attributes
            θs = geometry.opening_angle

            coords = ntuple(j -> zero(T), Val(4))

            isindir = false
            for _ = 1:2
                isindir ⊻= true
                ts, rs, ϕs = emission_coordinates(pix, geometry.opening_angle, isindir, n)
                if rs ≤ rmin || rs ≥ rmax
                    continue
                end
                coords = isnan(rs) ? observation : (ts, rs, θs, ϕs)
            end
            return coords
        end

        @testset "a=$a" for a in [0.99, 0.5, 0.01, -0.01, -0.5, -0.99]
            met = Krang.Kerr(a)
            @testset "$pixtype" for (pixtype, fcam) in [
                ("Intensity Pixel", Krang.IntensityCamera),
                ("Slow Light Intensity Pixel", Krang.SlowLightIntensityCamera),
            ]
                @testset "θs:$θs" for θs in
                                      [π / 5, π / 4, π / 3, π / 2, 2π / 3, 3π / 4, 4π / 5]
                    θo = π/4.1;
                    ρmax = 7.0;
                    sze = 400;
                    rmax = 1000.0; # maximum radius to be ray traced

                    pixels = Vector{Krang.AbstractPixel}(undef, 3)

                    cam = fcam(met, θo, -ρmax, ρmax, -ρmax, ρmax, sze);
                    geometries =
                        (Krang.ConeGeometry(θs, (i, Krang.horizon(met), rmax)) for i = 0:2)

                    for (i, geometry) in enumerate(geometries)
                        temp = false
                        pix = nothing
                        while !temp
                            pix = rand(cam.screen.pixels)
                            temp = coordinate_point(pix, geometry)[2] != 0
                        end
                        pixels[i] = pix
                    end

                    for n = 0:2
                        pix = pixels[n+1]
                        isindir = pix.screen_coordinate[2] > 0

                        ηcase1 = η(pix)
                        λcase1 = λ(pix)
                        roots = get_radial_roots(met, ηcase1, λcase1)
                        _, _, _, root = roots
                        τ1, _, _, _, _, issuccess = Krang.Gθ(pix, θs, isindir, n)
                        testθs, _, _, _, testn, _ = Krang.emission_inclination(pix, τ1)
                        @test testθs / θs ≈ 1 atol = 1e-5
                        @test testn == n
                    end
                end
            end
        end
    end
    @testset "Emission Coordinates" begin
        @testset "α:$α, β:$β, isindir:$isindir" for (α, β, isindir) in
                                                    [(10.0, 10.0, true), (1.0, -1.0, false)]
            θo = π / 4
            θs = π / 3
            a = 0.99
            met = Krang.Kerr(a)
            ηtemp = η(met, α, β, θo)
            λtemp = λ(met, α, θo)
            a2 = met.spin^2
            Δθ = (1.0 - (ηtemp + λtemp^2) / a2) / 2
            Δθ2 = Δθ^2
            desc = √(Δθ2 + ηtemp / a2)
            up = Δθ + desc
            θturning = acos(√up) * (1 + 1e-10)
            @testset "$pixtype" for (pixtype, pix) in [
                ("Intensity Pixel", Krang.IntensityPixel(met, α, β, θo)),
                (
                    "Cached Slow Light Intensity Pixel",
                    Krang.SlowLightIntensityPixel(met, α, β, θo),
                ),
            ]

                τ = Krang.Gθ(pix, θs, isindir, 0)[1]
                ts, testrs, ϕs, νr, νθ = Krang.emission_coordinates(pix, θs, isindir, 0)
                testrs2, ϕs2, νr2, νθ2 =
                    Krang.emission_coordinates_fast_light(pix, θs, isindir, 0)
                ts3, testrs3, testθs, ϕs3, νr3, νθ3 = Krang.emission_coordinates(pix, τ)
                testrs4, θs4, ϕs4, νr4, νθ4 = Krang.emission_coordinates_fast_light(pix, τ)

                testτ = Krang.Ir(pix, isindir, testrs)[1]
                @testset "Consistency between raytracing methods" begin
                    @test testrs / testrs2 ≈ 1.0 atol = 1e-5
                    @test ϕs2 / ϕs ≈ 1.0 atol = 1e-5
                    @test νr == νr2
                    @test νθ == νθ2
                    @test testrs / testrs3 ≈ 1.0 atol = 1e-5
                    @test ϕs / ϕs3 ≈ 1.0 atol = 1e-5
                    @test νr == νr3
                    @test νθ == νθ3
                    @test testrs / testrs4 ≈ 1.0 atol = 1e-5
                    @test θs4 / θs ≈ 1.0 atol = 1e-5
                    @test ϕs4 / ϕs ≈ 1.0 atol = 1e-5
                    @test νr == νr4
                    @test νθ == νθ4


                    @test testθs / θs ≈ 1.0 atol = 1e-5
                    @test testτ / τ ≈ 1.0 atol = 1e-5
                end

                fϕ(r, p) =
                    a *
                    (2r - a * λtemp) *
                    inv((r^2 - 2r + a^2) * √(r_potential(met, ηtemp, λtemp, r)))
                probϕ = IntegralProblem(fϕ, (testrs, Inf))
                solϕ = solve(probϕ, HCubatureJL(); reltol = 1e-8, abstol = 1e-8)
                Iϕ = solϕ.u

                gϕ(θ, p) = csc(θ)^2 * inv(√(Krang.θ_potential(met, ηtemp, λtemp, θ)))
                probϕs = IntegralProblem(gϕ, (θs, π / 2))
                solθs = solve(probϕs, HCubatureJL(); reltol = 1e-12, abstol = 1e-12)
                probϕo = IntegralProblem(gϕ, (θo, π / 2))
                solθo = solve(probϕo, HCubatureJL(); reltol = 1e-12, abstol = 1e-12)
                probϕ = IntegralProblem(gϕ, (θturning, π / 2))
                solθt = solve(probϕ, HCubatureJL(); reltol = 1e-12, abstol = 1e-12)

                Gϕ = 0
                if isindir
                    if sign(cos(θs) * cos(θo)) > 0
                        Gϕ = abs(2solθt.u - (solθo.u + solθs.u))
                    else
                        Gϕ = abs(2solθt.u + solθs.u - solθo.u)
                    end
                else
                    if sign(cos(θs) * cos(θo)) > 0
                        Gϕ = abs(solθo.u - solθs.u)
                    else
                        Gϕ = abs(solθo.u + solθs.u)
                    end
                end

                @test (Iϕ + λtemp * Gϕ) / ϕs ≈ 1.0 atol = 1e-3

                ft(r, p) = (
                    (r^2 * (r^2 - 2r + a^2) + 2r * (r^2 + a^2 - a * λtemp)) *
                    inv((r^2 - 2r + a^2) * √(r_potential(met, ηtemp, λtemp, r)))
                )
                probt = IntegralProblem(ft, (testrs, 1e6))
                solt = solve(probt, HCubatureJL(); reltol = 1e-8, abstol = 1e-8)
                It = solt.u

                gt(θ, p) = cos(θ)^2 * inv(√(Krang.θ_potential(met, ηtemp, λtemp, θ)))
                probts = IntegralProblem(gt, (θs, π / 2))
                solθs = solve(probts, HCubatureJL(); reltol = 1e-12, abstol = 1e-12)
                probto = IntegralProblem(gt, (θo, π / 2))
                solθo = solve(probto, HCubatureJL(); reltol = 1e-12, abstol = 1e-12)
                probt = IntegralProblem(gt, (θturning, π / 2))
                solθt = solve(probt, HCubatureJL(); reltol = 1e-12, abstol = 1e-12)

                Gt = 0
                if isindir
                    if sign(cos(θs) * cos(θo)) > 0
                        Gt = abs(2solθt.u - (solθo.u + solθs.u))
                    else
                        Gt = abs(2solθt.u + solθs.u - solθo.u)
                    end
                else
                    if sign(cos(θs) * cos(θo)) > 0
                        Gt = abs(solθo.u - solθs.u)
                    else
                        Gt = abs(solθo.u + solθs.u)
                    end
                end

                @test ((It - 1e6 - 2log(1e6)) + a^2 * Gt) / ts ≈ 1.0 atol = 1e-3
            end
        end
    end

    @testset "on-axis (θ_o→0,π) geodesics" begin
        m = Kerr(0.5)
        pix(α, β, θd) = Krang.SlowLightIntensityPixel(m, α, β, deg2rad(θd))
        coords(p, f) = Krang.emission_coordinates(p, f * p.total_mino_time)  # → (t,r,θ,φ,νr,νθ,ok)
        # r,θ continuous across the pole (vortical central pixel, ρ<a)
        p0, plim = pix(0.106, 0.106, 0.0), pix(0.106, 0.106, 1e-3)
        _, r0, θ0, _, _, _, ok0 = coords(p0, 0.45)
        _, rl, θl, _, _, _, _ = coords(plim, 0.45)
        @test ok0
        @test isapprox(r0, rl; rtol = 1e-4)
        @test isapprox(θ0, θl; rtol = 1e-3)
        # angular half-orbit nonzero on-axis (rejection would return 0)
        @test Krang._absGθo_Gθhat(m, deg2rad(1e-6), Krang.η(m, 0.106, 0.106, deg2rad(1e-6)), Krang.λ(m, 0.106, deg2rad(1e-6)))[2] > 0
        # φ continuous across θo=0 (τ samples span a pole crossing)
        for f in (0.1, 0.45, 0.85)
            @test isapprox(coords(pix(5.0, 2.0, 0.0), f)[4], coords(pix(5.0, 2.0, 1e-4), f)[4]; atol = 2e-3)
        end
        # screen azimuth not collapsed
        @test !isapprox(coords(pix(5.0, 2.0, 0.0), 0.45)[4], coords(pix(2.0, 5.0, 0.0), 0.45)[4]; atol = 1e-6)
        # south-pole azimuth matches its limit
        @test isapprox(coords(pix(5.0, 2.0, 180.0), 0.45)[4], coords(pix(5.0, 2.0, 180.0 - 1e-4), 0.45)[4]; atol = 2e-3)
    end

    @testset "α≈0 pixel (off-axis observer): pole-crossing azimuth" begin
        # λ = -α·sinθo → 0 geodesics cross the spin axis; φ must jump by π per crossing.
        # Degenerate at λ=0 (λ·Gϕ = 0·∞) and, pre-fix, numerically lost well before that
        # (1-up and the RJ p-argument cancel catastrophically).
        wrapφ(x) = mod(x + π, 2π) - π
        coords(p, f) = Krang.emission_coordinates(p, f * p.total_mino_time)
        # φ(α) continuous through α=0, F32 ≡ F64: compare against a trusted-α reference
        for T in (Float64, Float32), θod in (20.0, 60.0), β in (-5.0, 5.0)
            pref = Krang.SlowLightIntensityPixel(Kerr(0.9), 0.05, β, deg2rad(θod))
            for α in (0.0, 1e-6, 1e-4, 1e-3), f in (0.1, 0.35, 0.55, 0.8, 0.95)
                p = Krang.SlowLightIntensityPixel(Kerr(T(0.9)), T(α), T(β), T(deg2rad(θod)))
                c, cref = coords(p, T(f)), coords(pref, f)
                (c[7] && cref[7]) || continue
                min(c[3], π - c[3]) < 0.15 && continue  # skip the reference's own near-pole swing
                @test abs(wrapφ(Float64(c[4]) - cref[4])) < 0.1
            end
        end
        # φ jumps by ≈π across the first axis crossing (τ/τ_total: 0.305 → 0.318 spans it)
        for T in (Float64, Float32), α in (0.0, 1e-8, 1e-6)
            p = Krang.SlowLightIntensityPixel(Kerr(T(0.9)), T(α), T(-5.0), T(deg2rad(60)))
            Δφ = abs(wrapφ(Float64(coords(p, T(0.318))[4] - coords(p, T(0.305))[4])))
            @test 2.9 < Δφ < 3.4
        end
        # stable Π ≡ Legendre-API Π at moderate arguments (complete and incomplete)
        let up = 0.97, um = -3.2, m = up / um
            @test Krang._Pi_stable(up, m, 1.0, 0.0, 1 - up) ≈ Krang.JacobiElliptic.Pi(up, m) rtol = 1e-12
            @test Krang._Pi_stable(up, m, 0.6, 1 - 0.6^2, 1 - up * 0.6^2) ≈
                  Krang.JacobiElliptic.Pi(up, asin(0.6), m) rtol = 1e-12
        end
    end
end
