/************************************************************************************
*                      POLARIS: POLArized RadIation Simulator                       *
*                         Copyright (C) 2018 Stefan Reissl                          *
************************************************************************************/

#include "DustComponent.hpp"
#include "catch2/catch.hpp"

TEST_CASE("CDustComponent::henyeygreen", "[Dust][CDustComponent]")
{
    CDustComponent dust {};
    dust.setNrOfDustSpecies(1);
    dust.setNrOfWavelength(1);
    dust.initDustProperties();

    auto g = GENERATE(-0.99, -0.5, 0.0, 0.5, 0.99);
    dust.setHGg(0, 0, g);

    CRandomGenerator randgen {};
    randgen.init(123);

    photon_package pp{};
    pp.setWavelength(1.0, 0);
    pp.setPosition(Vector3D(0.0, 0.0, 0.0));
    Vector3D ez {Vector3D(0.0, 0.0, 1.0)};

    // mu = cos(theta)
    double mu {};
    double mu_mean {0.0};
    double mu_min {1.0};
    double mu_max {-1.0};
    constexpr std::size_t n_samples {1024 * 1024};

    for (std::size_t i = 0; i < n_samples; ++i) {
        pp.setDirection(ez);
        pp.initCoordSystem();
        dust.henyeygreen(&pp, 0, &randgen);
        mu = ez * pp.getDirection();
        mu_mean += mu;
        mu_min = std::min(mu, mu_min);
        mu_max = std::max(mu, mu_max);
    }

    CAPTURE(g);

    mu_mean /= static_cast<double>(n_samples);
    // mu should be between -1 and 1
    REQUIRE(mu_min >= -1.0);
    REQUIRE(mu_max <= 1.0);
    REQUIRE(mu_min < mu_max);
    // mean value should be approximately g
    REQUIRE(std::abs(mu_mean - g) < 0.001);
}

TEST_CASE("CDustComponent::drainehenyeygreen", "[Dust][CDustComponent]")
{
    CDustComponent dust {};
    dust.setNrOfDustSpecies(1);
    dust.setNrOfWavelength(1);
    dust.initDustProperties();

    auto g = GENERATE(-0.99, -0.5, 0.0, 0.5, 0.99);
    auto alpha = GENERATE(0.0, 0.33, 0.67, 1.0);
    dust.setHGg(0, 0, g);
    dust.setHGg2(0, 0, alpha);

    CRandomGenerator randgen {};
    randgen.init(123);

    photon_package pp {};
    pp.setWavelength(1.0, 0);
    pp.setPosition(Vector3D(0.0, 0.0, 0.0));
    Vector3D ez {Vector3D(0.0, 0.0, 1.0)};

    // mu = cos(theta)
    double mu {};
    double mu_mean {0.0};
    double mu_min {1.0};
    double mu_max {-1.0};
    constexpr std::size_t n_samples {1024 * 1024};

    for (std::size_t i = 0; i < n_samples; ++i) {
        pp.setDirection(ez);
        pp.initCoordSystem();
        dust.drainehenyeygreen(&pp, 0, &randgen);
        mu = ez * pp.getDirection();
        mu_mean += mu;
        mu_min = std::min(mu, mu_min);
        mu_max = std::max(mu, mu_max);
    }

    CAPTURE(g);
    CAPTURE(alpha);

    mu_mean /= static_cast<double>(n_samples);
    // mu should be between -1 and 1
    REQUIRE(mu_min >= -1.0);
    REQUIRE(mu_max <= 1.0);
    REQUIRE(mu_min < mu_max);
    // mean value should be approximately Eq. (A2) in Draine 2003, ApJ 598, 1017
    double first_moment {g * (1.0 + alpha * (3.0 + 2.0 * g * g) / 5.0) / (1.0 + alpha * (1.0 + 2.0 * g * g) / 3.0)};
    REQUIRE(std::abs(mu_mean - first_moment) < 0.01);
}

TEST_CASE("CDustComponent::threeparamhenyeygreen", "[Dust][CDustComponent]")
{
    CDustComponent dust {};
    dust.setNrOfDustSpecies(1);
    dust.setNrOfWavelength(1);
    dust.initDustProperties();

    auto g1 = GENERATE(0.0, 0.5, 0.99);
    auto g2 = GENERATE(-0.99, -0.5, 0.0);
    auto weight = GENERATE(0.0, 0.33, 0.67, 1.0);
    dust.setHGg(0, 0, g1);
    dust.setHGg2(0, 0, g2);
    dust.setHGg3(0, 0, weight);

    CRandomGenerator randgen {};
    randgen.init(123);

    photon_package pp {};
    pp.setWavelength(1.0, 0);
    pp.setPosition(Vector3D(0.0, 0.0, 0.0));
    Vector3D ez {Vector3D(0.0, 0.0, 1.0)};

    // mu = cos(theta)
    double mu {};
    double mu_mean {0.0};
    double mu_min {1.0};
    double mu_max {-1.0};
    constexpr std::size_t n_samples {1024 * 1024};

    for (std::size_t i = 0; i < n_samples; ++i) {
        pp.setDirection(ez);
        pp.initCoordSystem();
        dust.threeparamhenyeygreen(&pp, 0, &randgen);
        mu = ez * pp.getDirection();
        mu_mean += mu;
        mu_min = std::min(mu, mu_min);
        mu_max = std::max(mu, mu_max);
    }

    CAPTURE(g1);
    CAPTURE(g2);
    CAPTURE(weight);

    mu_mean /= static_cast<double>(n_samples);
    // mu should be between -1 and 1
    REQUIRE(mu_min >= -1.0);
    REQUIRE(mu_max <= 1.0);
    REQUIRE(mu_min < mu_max);
    // mean value should be approximately w*g1 + (1-w)*g2
    double first_moment {weight * g1 + (1.0 - weight) * g2};
    REQUIRE(std::abs(mu_mean - first_moment) < 0.01);
}
