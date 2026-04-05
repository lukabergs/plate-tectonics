#ifndef TOPOGRAPHY_CODEC_HPP
#define TOPOGRAPHY_CODEC_HPP

#include <cstddef>
#include <string>

#include "utils.hpp"

namespace TopographyCodec {

constexpr uint16_t kMaxHeightMeters = 65535;
constexpr int32_t kNoSeaLevelOverride = -1;
constexpr uint32_t kMetadataVersion = 1;
constexpr float kOceanicBase = 0.1f;
constexpr float kContinentalBase = 1.0f;
constexpr const char* kLittleEndian = "little";
constexpr const char* kBigEndian = "big";
constexpr const char* kRowMajorLayout = "row-major";
constexpr const char* kMetricFormatR16 = "r16";
constexpr const char* kMetricFormatPng16 = "png16";
constexpr const char* kMetricFormatFloat32 = "float32";

struct Metadata {
    uint32_t width = 0;
    uint32_t height = 0;
    uint16_t sea_level_m = kMaxHeightMeters / 2;
    uint16_t max_height_m = kMaxHeightMeters;
    uint32_t version = kMetadataVersion;
    std::string format;
    std::string endianness = kLittleEndian;
    std::string layout = kRowMajorLayout;
};

float meters_to_internal(uint16_t meters, uint16_t sea_level_m);
uint16_t internal_to_meters(float internal, uint16_t sea_level_m);
float normalized_to_internal(float normalized, float legacy_sea_threshold);
bool is_oceanic_internal(float value);

float clamp_normalized(float normalized);
float legacy_raw_sea_ratio();
uint16_t legacy_raw_sea_level_m();
float infer_normalized_sea_threshold(const float* values, size_t count, float ocean_coverage);
uint16_t infer_metric_sea_level(const uint16_t* values, size_t count, float ocean_coverage);
bool host_is_little_endian();

} // namespace TopographyCodec

#endif
