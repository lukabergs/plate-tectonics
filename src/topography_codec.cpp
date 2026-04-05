#include "topography_codec.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace TopographyCodec {

namespace {

constexpr float kInternalLandReliefRange = 1.0f;

float oceanic_internal_ceiling()
{
    return std::nextafter(kContinentalBase, kOceanicBase);
}

float oceanic_internal_range()
{
    return oceanic_internal_ceiling() - kOceanicBase;
}

template <typename T>
T clamp_value(T value, T min_value, T max_value)
{
    return std::min(max_value, std::max(min_value, value));
}

} // namespace

float clamp_normalized(float normalized)
{
    return clamp_value(normalized, 0.0f, 1.0f);
}

float meters_to_internal(uint16_t meters, uint16_t sea_level_m)
{
    if (meters < sea_level_m) {
        if (sea_level_m <= 1) {
            return kOceanicBase;
        }

        const float ratio = static_cast<float>(meters) / static_cast<float>(sea_level_m - 1);
        return kOceanicBase + ratio * oceanic_internal_range();
    }

    if (sea_level_m >= kMaxHeightMeters) {
        if (meters >= kMaxHeightMeters) {
            return kContinentalBase + kInternalLandReliefRange;
        }
        return kOceanicBase +
               (static_cast<float>(meters) / static_cast<float>(kMaxHeightMeters - 1)) *
                   oceanic_internal_range();
    }

    const float ratio = static_cast<float>(meters - sea_level_m) /
                        static_cast<float>(kMaxHeightMeters - sea_level_m);
    return kContinentalBase + ratio * kInternalLandReliefRange;
}

uint16_t internal_to_meters(float internal, uint16_t sea_level_m)
{
    if (internal <= kOceanicBase) {
        return 0;
    }

    if (is_oceanic_internal(internal)) {
        if (sea_level_m <= 1) {
            return 0;
        }

        const float ratio =
            (clamp_value(internal, kOceanicBase, oceanic_internal_ceiling()) - kOceanicBase) /
            oceanic_internal_range();
        return static_cast<uint16_t>(
            std::lround(ratio * static_cast<float>(sea_level_m - 1)));
    }

    if (sea_level_m >= kMaxHeightMeters) {
        return kMaxHeightMeters;
    }

    const float ratio =
        (clamp_value(internal, kContinentalBase, kContinentalBase + kInternalLandReliefRange) -
         kContinentalBase) /
        kInternalLandReliefRange;
    const uint32_t meters = static_cast<uint32_t>(std::lround(
        static_cast<float>(sea_level_m) +
        ratio * static_cast<float>(kMaxHeightMeters - sea_level_m)));
    return static_cast<uint16_t>(clamp_value<uint32_t>(meters, sea_level_m, kMaxHeightMeters));
}

float normalized_to_internal(float normalized, float legacy_sea_threshold)
{
    const float clamped = clamp_normalized(normalized);
    return clamped > legacy_sea_threshold ? clamped + kContinentalBase : kOceanicBase;
}

bool is_oceanic_internal(float value)
{
    return value < kContinentalBase;
}

float legacy_raw_sea_ratio()
{
    return (kContinentalBase - kOceanicBase) /
           ((kContinentalBase + kInternalLandReliefRange) - kOceanicBase);
}

uint16_t legacy_raw_sea_level_m()
{
    return static_cast<uint16_t>(
        std::lround(legacy_raw_sea_ratio() * static_cast<float>(kMaxHeightMeters)));
}

float infer_normalized_sea_threshold(const float* values, size_t count, float ocean_coverage)
{
    if (values == nullptr || count == 0) {
        return 0.5f;
    }

    float threshold = 0.5f;
    float step = 0.5f;
    const float target = clamp_normalized(ocean_coverage);
    while (step > 0.00001f) {
        size_t below_count = 0;
        for (size_t i = 0; i < count; ++i) {
            below_count += clamp_normalized(values[i]) < threshold;
        }

        step *= 0.5f;
        if (static_cast<float>(below_count) / static_cast<float>(count) < target) {
            threshold += step;
        } else {
            threshold -= step;
        }
    }

    return clamp_normalized(threshold);
}

uint16_t infer_metric_sea_level(const uint16_t* values, size_t count, float ocean_coverage)
{
    if (values == nullptr || count == 0) {
        return static_cast<uint16_t>(std::lround(
            clamp_normalized(ocean_coverage) * static_cast<float>(kMaxHeightMeters)));
    }

    std::vector<uint16_t> sorted(values, values + count);
    std::sort(sorted.begin(), sorted.end());

    const float target = clamp_normalized(ocean_coverage);
    const size_t index = static_cast<size_t>(
        std::floor(target * static_cast<float>(sorted.size() - 1)));
    return sorted[index];
}

bool host_is_little_endian()
{
    const uint16_t value = 0x0100;
    const unsigned char* bytes = reinterpret_cast<const unsigned char*>(&value);
    return bytes[1] == 0x01;
}

} // namespace TopographyCodec
