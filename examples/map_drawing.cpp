#ifndef NOMINMAX
#define NOMINMAX
#endif

#include "map_drawing.hpp"
#include "utils.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <vector>

using namespace std;

namespace {

using Byte = png_byte;

constexpr int kChannelsPerPixel = 3;
constexpr int kReadChannelsPerPixel = 4;
constexpr float kFaultNormalThreshold = 0.2f;
constexpr float kFaultTangentialThreshold = 0.35f;
constexpr float kInvSqrt2 = 0.70710677f;

inline float rgb_to_luminance(Byte r, Byte g, Byte b)
{
    return (0.299f * static_cast<float>(r) + 0.587f * static_cast<float>(g) +
            0.114f * static_cast<float>(b)) /
           255.0f;
}

const Byte* fault_color(FaultKind kind)
{
    static constexpr Byte convergent[3] = {255, 84, 84};
    static constexpr Byte divergent[3] = {80, 240, 255};
    static constexpr Byte transform[3] = {255, 80, 220};
    static constexpr Byte none[3] = {0, 0, 0};

    switch (kind) {
    case FaultKind::Convergent:
        return convergent;
    case FaultKind::Divergent:
        return divergent;
    case FaultKind::Transform:
        return transform;
    case FaultKind::None:
    default:
        return none;
    }
}

inline void setGray(Byte* ptr, int val)
{
    ptr[0] = static_cast<Byte>(val);
    ptr[1] = static_cast<Byte>(val);
    ptr[2] = static_cast<Byte>(val);
}

inline void setColor(Byte* ptr, Byte r, Byte g, Byte b)
{
    ptr[0] = r;
    ptr[1] = g;
    ptr[2] = b;
}

float sample_height(const float* heightmap, int width, int height, int x, int y)
{
    const int clamped_x = std::clamp(x, 0, width - 1);
    const int clamped_y = std::clamp(y, 0, height - 1);
    return heightmap[clamped_y * width + clamped_x];
}

float relief_shade(const float* heightmap, int width, int height, int x, int y)
{
    const float dx = sample_height(heightmap, width, height, x - 1, y) -
                     sample_height(heightmap, width, height, x + 1, y);
    const float dy = sample_height(heightmap, width, height, x, y - 1) -
                     sample_height(heightmap, width, height, x, y + 1);
    return std::clamp(dx * 0.35f + dy * 0.55f, -0.18f, 0.18f);
}

void apply_shade(Byte* ptr, float shade)
{
    const float factor = std::clamp(1.0f + shade, 0.7f, 1.3f);
    ptr[0] = static_cast<Byte>(std::clamp(std::round(static_cast<float>(ptr[0]) * factor), 0.0f, 255.0f));
    ptr[1] = static_cast<Byte>(std::clamp(std::round(static_cast<float>(ptr[1]) * factor), 0.0f, 255.0f));
    ptr[2] = static_cast<Byte>(std::clamp(std::round(static_cast<float>(ptr[2]) * factor), 0.0f, 255.0f));
}

int write_png(const char* filename, int width, int height, const std::vector<Byte>& pixels)
{
    if (filename == nullptr || width <= 0 || height <= 0) {
        return 1;
    }

    png_image image;
    memset(&image, 0, sizeof(image));
    image.version = PNG_IMAGE_VERSION;
    image.width = static_cast<png_uint_32>(width);
    image.height = static_cast<png_uint_32>(height);
    image.format = PNG_FORMAT_RGB;

    if (!png_image_write_to_file(&image, filename, 0, pixels.data(), 0, nullptr)) {
        fprintf(stderr, "Could not write PNG file %s: %s\n", filename,
                image.message != nullptr ? image.message : "unknown libpng error");
        png_image_free(&image);
        return 1;
    }

    png_image_free(&image);
    return 0;
}

float find_value_for_quantile(float quantile, const float* array, std::uint32_t size)
{
    float value = 0.5f;
    float th_step = 0.5f;

    while (th_step > 0.00001f) {
        uint32_t count = 0;
        for (uint32_t i = 0; i < size; ++i) {
            count += (array[i] < value);
        }

        th_step *= 0.5f;
        if (count / static_cast<float>(size) < quantile) {
            value += th_step;
        } else {
            value -= th_step;
        }
    }

    return value;
}

void gradient(Byte* ptr, Byte ra, Byte ga, Byte ba, Byte rb, Byte gb, Byte bb, float h, float ha,
              float hb)
{
    if (ha > hb) {
        throw runtime_error("gradient bounds are inverted");
    }

    if (hb <= ha) {
        setColor(ptr, rb, gb, bb);
        return;
    }

    const float clamped_h = std::clamp(h, ha, hb);
    const float simil_b = (clamped_h - ha) / (hb - ha);
    const float simil_a = 1.0f - simil_b;
    setColor(ptr,
             static_cast<Byte>(simil_a * static_cast<float>(ra) + simil_b * static_cast<float>(rb)),
             static_cast<Byte>(simil_a * static_cast<float>(ga) + simil_b * static_cast<float>(gb)),
             static_cast<Byte>(simil_a * static_cast<float>(ba) + simil_b * static_cast<float>(bb)));
}

void drawGrayImage(std::vector<Byte>& pixels, int width, int height, float* heightmap)
{
    const auto size = static_cast<std::uint32_t>(width * height);
    const float q01 = find_value_for_quantile(0.01f, heightmap, size);
    const float q995 = find_value_for_quantile(0.995f, heightmap, size);
    const float inv_range = (q995 > q01) ? 1.0f / (q995 - q01) : 1.0f;

    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            const std::size_t offset =
                static_cast<std::size_t>((y * width + x) * kChannelsPerPixel);
            const float h = heightmap[y * width + x];
            const float tone = std::clamp((h - q01) * inv_range, 0.0f, 1.0f);
            const float shaded =
                std::clamp(std::pow(tone, 0.82f) + relief_shade(heightmap, width, height, x, y), 0.0f, 1.0f);
            setGray(&pixels[offset], static_cast<int>(shaded * 255.0f));
        }
    }
}

void drawColorsImage(std::vector<Byte>& pixels, int width, int height, float* heightmap)
{
    const auto size = static_cast<std::uint32_t>(width * height);
    const float q15 = find_value_for_quantile(0.15f, heightmap, size);
    const float q70 = find_value_for_quantile(0.70f, heightmap, size);
    const float q75 = find_value_for_quantile(0.75f, heightmap, size);
    const float q90 = find_value_for_quantile(0.90f, heightmap, size);
    const float q95 = find_value_for_quantile(0.95f, heightmap, size);
    const float q99 = find_value_for_quantile(0.99f, heightmap, size);

    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            const std::size_t offset =
                static_cast<std::size_t>((y * width + x) * kChannelsPerPixel);
            const float h = heightmap[y * width + x];

            if (h < q15) {
                gradient(&pixels[offset], 4, 18, 84, 12, 64, 168, h, 0.0f, q15);
            } else if (h < q70) {
                gradient(&pixels[offset], 12, 64, 168, 72, 148, 222, h, q15, q70);
            } else if (h < q75) {
                gradient(&pixels[offset], 72, 148, 222, 188, 225, 216, h, q70, q75);
            } else if (h < q90) {
                gradient(&pixels[offset], 71, 132, 61, 161, 181, 84, h, q75, q90);
            } else if (h < q95) {
                gradient(&pixels[offset], 161, 181, 84, 214, 198, 112, h, q90, q95);
            } else if (h < q99) {
                gradient(&pixels[offset], 214, 198, 112, 145, 98, 62, h, q95, q99);
            } else {
                gradient(&pixels[offset], 145, 98, 62, 92, 62, 46, h, q99, 1.0f);
            }

            apply_shade(&pixels[offset], relief_shade(heightmap, width, height, x, y) * 0.8f);
        }
    }
}

void overlay_fault_map(std::vector<Byte>& pixels, const FaultKind* fault_map)
{
    if (fault_map == nullptr) {
        return;
    }

    const std::size_t pixel_count = pixels.size() / kChannelsPerPixel;
    for (std::size_t i = 0; i < pixel_count; ++i) {
        const FaultKind kind = fault_map[i];
        if (kind == FaultKind::None) {
            continue;
        }

        const Byte* color = fault_color(kind);
        const std::size_t offset = i * kChannelsPerPixel;
        pixels[offset + 0] = static_cast<Byte>(
            (static_cast<unsigned int>(pixels[offset + 0]) + 4U * color[0]) / 5U);
        pixels[offset + 1] = static_cast<Byte>(
            (static_cast<unsigned int>(pixels[offset + 1]) + 4U * color[1]) / 5U);
        pixels[offset + 2] = static_cast<Byte>(
            (static_cast<unsigned int>(pixels[offset + 2]) + 4U * color[2]) / 5U);
    }
}

FaultKind choose_fault_kind(FaultKind current, FaultKind candidate)
{
    return static_cast<int>(candidate) > static_cast<int>(current) ? candidate : current;
}

FaultKind classify_fault(uint32_t first_plate, uint32_t second_plate, float normal_x,
                         float normal_y, const float* plate_vx, const float* plate_vy,
                         uint32_t plate_count)
{
    if (first_plate == second_plate || first_plate >= plate_count || second_plate >= plate_count ||
        plate_vx == nullptr || plate_vy == nullptr) {
        return FaultKind::None;
    }

    const float rel_x = plate_vx[second_plate] - plate_vx[first_plate];
    const float rel_y = plate_vy[second_plate] - plate_vy[first_plate];
    const float normal_motion = rel_x * normal_x + rel_y * normal_y;
    const float tangential_motion = std::abs(rel_x * -normal_y + rel_y * normal_x);

    if (normal_motion <= -kFaultNormalThreshold) {
        return FaultKind::Convergent;
    }
    if (normal_motion >= kFaultNormalThreshold) {
        return FaultKind::Divergent;
    }
    if (tangential_motion >= kFaultTangentialThreshold) {
        return FaultKind::Transform;
    }

    return FaultKind::None;
}

void mark_fault_pixel(std::vector<FaultKind>& fault_map, int width, int x, int y, FaultKind kind)
{
    if (kind == FaultKind::None) {
        return;
    }

    const std::size_t index = static_cast<std::size_t>(y) * static_cast<std::size_t>(width) +
                              static_cast<std::size_t>(x);
    fault_map[index] = choose_fault_kind(fault_map[index], kind);
}

void drawImage(std::vector<Byte>& pixels, int width, int height, float* heightmap,
               void(drawFunction)(std::vector<Byte>&, int, int, float*))
{
    pixels.resize(static_cast<std::size_t>(width) * static_cast<std::size_t>(height) *
                  kChannelsPerPixel);
    drawFunction(pixels, width, height, heightmap);
}

int writeImageImpl(const char* filename, int width, int height, float* heightmap,
                   void(drawFunction)(std::vector<Byte>&, int, int, float*),
                   const FaultKind* fault_map)
{
    if (heightmap == nullptr || width <= 0 || height <= 0) {
        return 1;
    }

    std::vector<Byte> pixels;
    drawImage(pixels, width, height, heightmap, drawFunction);
    overlay_fault_map(pixels, fault_map);
    return write_png(filename, width, height, pixels);
}

} // namespace

int writeImageGray(const char* filename, int width, int height, float* heightmap,
                   [[maybe_unused]] const char* title)
{
    return writeImageImpl(filename, width, height, heightmap, drawGrayImage, nullptr);
}

int writeImageColors(const char* filename, int width, int height, float* heightmap,
                     [[maybe_unused]] const char* title)
{
    return writeImageImpl(filename, width, height, heightmap, drawColorsImage, nullptr);
}

int buildFaultMap(int width, int height, const uint32_t* platesmap, const float* plate_vx,
                  const float* plate_vy, uint32_t plate_count, std::vector<FaultKind>* fault_map)
{
    if (fault_map == nullptr || width <= 0 || height <= 0) {
        return 1;
    }

    fault_map->assign(static_cast<std::size_t>(width) * static_cast<std::size_t>(height),
                      FaultKind::None);

    if (platesmap == nullptr || plate_vx == nullptr || plate_vy == nullptr || plate_count == 0) {
        return 0;
    }

    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            const std::size_t index = static_cast<std::size_t>(y) * static_cast<std::size_t>(width) +
                                      static_cast<std::size_t>(x);
            const uint32_t current = platesmap[index];

            if (x + 1 < width) {
                mark_fault_pixel(*fault_map, width, x, y,
                                 classify_fault(current, platesmap[index + 1], 1.0f, 0.0f,
                                                plate_vx, plate_vy, plate_count));
            }
            if (y + 1 < height) {
                mark_fault_pixel(*fault_map, width, x, y,
                                 classify_fault(current,
                                                platesmap[index + static_cast<std::size_t>(width)],
                                                0.0f, 1.0f, plate_vx, plate_vy, plate_count));
            }
            if (x + 1 < width && y + 1 < height) {
                mark_fault_pixel(*fault_map, width, x, y,
                                 classify_fault(
                                     current, platesmap[index + static_cast<std::size_t>(width) + 1],
                                     kInvSqrt2, kInvSqrt2, plate_vx, plate_vy, plate_count));
            }
            if (x > 0 && y + 1 < height) {
                mark_fault_pixel(*fault_map, width, x, y,
                                 classify_fault(
                                     current, platesmap[index + static_cast<std::size_t>(width) - 1],
                                     -kInvSqrt2, kInvSqrt2, plate_vx, plate_vy, plate_count));
            }
        }
    }

    return 0;
}

int writeImageGrayWithFaultMap(const char* filename, int width, int height, float* heightmap,
                               const FaultKind* fault_map, [[maybe_unused]] const char* title)
{
    return writeImageImpl(filename, width, height, heightmap, drawGrayImage, fault_map);
}

int writeImageColorsWithFaultMap(const char* filename, int width, int height, float* heightmap,
                                 const FaultKind* fault_map, [[maybe_unused]] const char* title)
{
    return writeImageImpl(filename, width, height, heightmap, drawColorsImage, fault_map);
}

int writeImageGrayWithFaultLines(const char* filename, int width, int height, float* heightmap,
                                 const uint32_t* platesmap, const float* plate_vx,
                                 const float* plate_vy, uint32_t plate_count,
                                 [[maybe_unused]] const char* title)
{
    std::vector<FaultKind> fault_map;
    if (buildFaultMap(width, height, platesmap, plate_vx, plate_vy, plate_count, &fault_map) != 0) {
        return 1;
    }

    return writeImageGrayWithFaultMap(filename, width, height, heightmap, fault_map.data(), title);
}

int writeImageColorsWithFaultLines(const char* filename, int width, int height, float* heightmap,
                                   const uint32_t* platesmap, const float* plate_vx,
                                   const float* plate_vy, uint32_t plate_count,
                                   [[maybe_unused]] const char* title)
{
    std::vector<FaultKind> fault_map;
    if (buildFaultMap(width, height, platesmap, plate_vx, plate_vy, plate_count, &fault_map) != 0) {
        return 1;
    }

    return writeImageColorsWithFaultMap(filename, width, height, heightmap, fault_map.data(),
                                        title);
}

int loadHeightmapFromPng(const char* filename, uint32_t* width, uint32_t* height,
                         std::vector<float>* heightmap)
{
    if (filename == nullptr || width == nullptr || height == nullptr || heightmap == nullptr) {
        return 1;
    }

    png_image image;
    memset(&image, 0, sizeof(image));
    image.version = PNG_IMAGE_VERSION;

    if (!png_image_begin_read_from_file(&image, filename)) {
        fprintf(stderr, "Could not read PNG file %s: %s\n", filename,
                image.message != nullptr ? image.message : "unknown libpng error");
        png_image_free(&image);
        return 1;
    }

    image.format = PNG_FORMAT_RGBA;
    std::vector<Byte> pixels(PNG_IMAGE_SIZE(image));
    if (!png_image_finish_read(&image, nullptr, pixels.data(), 0, nullptr)) {
        fprintf(stderr, "Could not decode PNG file %s: %s\n", filename,
                image.message != nullptr ? image.message : "unknown libpng error");
        png_image_free(&image);
        return 1;
    }

    const std::size_t pixel_count =
        static_cast<std::size_t>(image.width) * static_cast<std::size_t>(image.height);
    heightmap->resize(pixel_count);
    for (std::size_t i = 0; i < pixel_count; ++i) {
        const std::size_t offset = i * kReadChannelsPerPixel;
        (*heightmap)[i] =
            rgb_to_luminance(pixels[offset + 0], pixels[offset + 1], pixels[offset + 2]);
    }

    *width = image.width;
    *height = image.height;
    png_image_free(&image);
    return 0;
}

int loadHeightmapFromImage(const char* filename, uint32_t* width, uint32_t* height,
                           std::vector<float>* heightmap)
{
    return loadHeightmapFromPng(filename, width, height, heightmap);
}
