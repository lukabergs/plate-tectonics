#include "gtest/gtest.h"
#include "topography_codec.hpp"

TEST(TopographyCodec, MetricRoundTripAroundSeaLevel)
{
    const uint16_t sea_level_m = 32768;
    const uint16_t samples[] = {
        0,
        1,
        static_cast<uint16_t>(sea_level_m - 1),
        sea_level_m,
        static_cast<uint16_t>(sea_level_m + 1),
        TopographyCodec::kMaxHeightMeters,
    };

    for (uint16_t sample : samples) {
        const float internal = TopographyCodec::meters_to_internal(sample, sea_level_m);
        EXPECT_EQ(sample, TopographyCodec::internal_to_meters(internal, sea_level_m));
    }
}

TEST(TopographyCodec, OceanicClassificationMatchesCodecBoundary)
{
    const uint16_t sea_level_m = 32000;
    EXPECT_TRUE(TopographyCodec::is_oceanic_internal(
        TopographyCodec::meters_to_internal(static_cast<uint16_t>(sea_level_m - 1), sea_level_m)));
    EXPECT_FALSE(TopographyCodec::is_oceanic_internal(
        TopographyCodec::meters_to_internal(sea_level_m, sea_level_m)));
}
