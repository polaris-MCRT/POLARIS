/************************************************************************************
*                      POLARIS: POLArized RadIation Simulator                       *
*                         Copyright (C) 2018 Stefan Reissl                          *
************************************************************************************/

#include "RandomGenerator.hpp"
#include "catch2/catch.hpp"
#include <algorithm>

TEST_CASE("CRandomGenerator::getRND", "[MathFunctions][CRandomGenerator]") {
    CRandomGenerator rng{};
    rng.init(123);
    double r_mean{0};
    double r_min{2};
    double r_max{-1};
    constexpr std::size_t n_samples{1024 * 1024};
    for (std::size_t i = 0; i < n_samples; ++i) {
        double r{rng.getRND()};
        r_mean += r;
        r_min = std::min(r, r_min);
        r_max = std::max(r, r_max);
    }
    r_mean /= static_cast<double>(n_samples);
    // getRND returns a double between 0 and 1 with uniform distribution
    REQUIRE(r_min >= 0);
    REQUIRE(r_max <= 1);
    REQUIRE(r_min < r_max);
    // mean value should be approximately 0.5
    REQUIRE(std::abs(r_mean - 0.5) < 0.01);
}

TEST_CASE("CRandomGenerator::getRNDnormal", "[MathFunctions][CRandomGenerator]") {
    CRandomGenerator rng{};
    rng.init(123);
    double r_mean{0};
    double r_min{2};
    constexpr std::size_t n_samples{1024*1024};
    // getRNDnormal returns a normally distributed double
    // unless it would be negative, in which case it generates a new one
    // choose mu, sigma here such that this case will rarely occur
    auto mu = GENERATE(1.0, 2.3);
    auto sigma = GENERATE(0.02, 0.08);
    std::size_t n_within_1sd{0};
    for (std::size_t i = 0; i < n_samples; ++i) {
        double r{rng.getRNDnormal(mu, sigma)};
        r_mean += r;
        r_min = std::min(r, r_min);
        if(std::abs(r-mu) < sigma){
            ++n_within_1sd;
        }
    }
    r_mean /= static_cast<double>(n_samples);
    double r_frac_within_1sd{static_cast<double>(n_within_1sd)/n_samples};
    CAPTURE(mu);
    CAPTURE(sigma);
    // all values should be non-negative
    REQUIRE(r_min >= 0);
    // mean value should be approximately mu
    REQUIRE(std::abs(r_mean - mu) < 0.001);
    // ~68% of values should lie within mu +/- sigma
    constexpr double frac_within_1sd{0.6827};
    REQUIRE(std::abs(r_frac_within_1sd - frac_within_1sd) < 0.001);
}
