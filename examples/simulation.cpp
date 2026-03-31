#ifndef NOMINMAX
#define NOMINMAX
#endif

#include <stdio.h>
#include "platecapi.hpp"
#include "sqrdmd.hpp"
#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>
#include "map_drawing.hpp"
#include <iostream>
#include <utils.hpp>

namespace {

struct FaultSnapshot {
    std::vector<uint32_t> platesmap;
    std::vector<float> plate_vx;
    std::vector<float> plate_vy;
    uint32_t plate_count = 0;
    bool valid = false;
};

struct FaultAccumulator {
    std::vector<uint32_t> convergent;
    std::vector<uint32_t> divergent;
    std::vector<uint32_t> transform;
    uint32_t samples = 0;
    bool valid = false;
};

void prune_fault_clusters(std::vector<FaultKind>& fault_map, int width, int height)
{
    const std::vector<FaultKind> base = fault_map;
    auto index_of = [width](int x, int y) -> size_t {
        return static_cast<size_t>(y) * static_cast<size_t>(width) + static_cast<size_t>(x);
    };

    for (int y = 1; y + 1 < height; ++y) {
        for (int x = 1; x + 1 < width; ++x) {
            const size_t index = index_of(x, y);
            const FaultKind kind = base[index];
            if (kind == FaultKind::None) {
                continue;
            }

            int same_neighbours = 0;
            for (int dy = -1; dy <= 1; ++dy) {
                for (int dx = -1; dx <= 1; ++dx) {
                    if (dx == 0 && dy == 0) {
                        continue;
                    }
                    same_neighbours += static_cast<int>(base[index_of(x + dx, y + dy)] == kind);
                }
            }

            if (same_neighbours < 5) {
                continue;
            }

            const bool horizontal = base[index_of(x - 1, y)] == kind &&
                                    base[index_of(x + 1, y)] == kind;
            const bool vertical = base[index_of(x, y - 1)] == kind &&
                                  base[index_of(x, y + 1)] == kind;
            const bool diagonal_a = base[index_of(x - 1, y - 1)] == kind &&
                                    base[index_of(x + 1, y + 1)] == kind;
            const bool diagonal_b = base[index_of(x - 1, y + 1)] == kind &&
                                    base[index_of(x + 1, y - 1)] == kind;
            const int axis_count = static_cast<int>(horizontal) + static_cast<int>(vertical) +
                                   static_cast<int>(diagonal_a) + static_cast<int>(diagonal_b);
            if (axis_count < 2) {
                continue;
            }

            const uint32_t hash = static_cast<uint32_t>(
                (x * 73856093) ^ (y * 19349663) ^ (static_cast<int>(kind) * 83492791));
            if ((hash & 1U) == 0U) {
                fault_map[index] = FaultKind::None;
            }
        }
    }
}

void roughen_fault_map(std::vector<FaultKind>& fault_map, int width, int height)
{
    const std::vector<FaultKind> base = fault_map;
    std::vector<FaultKind> warped(static_cast<size_t>(width) * static_cast<size_t>(height),
                                  FaultKind::None);
    auto index_of = [width](int x, int y) -> size_t {
        return static_cast<size_t>(y) * static_cast<size_t>(width) + static_cast<size_t>(x);
    };
    auto place_fault = [&](int x, int y, FaultKind kind) {
        if (x < 0 || x >= width || y < 0 || y >= height) {
            return;
        }
        const size_t target = index_of(x, y);
        if (warped[target] == FaultKind::None || warped[target] == kind) {
            warped[target] = kind;
        }
    };

    for (int x = 0; x < width; ++x) {
        place_fault(x, 0, base[index_of(x, 0)]);
        place_fault(x, height - 1, base[index_of(x, height - 1)]);
    }
    for (int y = 1; y + 1 < height; ++y) {
        place_fault(0, y, base[index_of(0, y)]);
        place_fault(width - 1, y, base[index_of(width - 1, y)]);
    }

    for (int y = 1; y + 1 < height; ++y) {
        for (int x = 1; x + 1 < width; ++x) {
            const size_t index = index_of(x, y);
            const FaultKind kind = base[index];
            if (kind == FaultKind::None) {
                continue;
            }

            const int horizontal = static_cast<int>(base[index_of(x - 1, y)] == kind) +
                                   static_cast<int>(base[index_of(x + 1, y)] == kind);
            const int vertical = static_cast<int>(base[index_of(x, y - 1)] == kind) +
                                 static_cast<int>(base[index_of(x, y + 1)] == kind);
            const int diagonal_a =
                static_cast<int>(base[index_of(x - 1, y - 1)] == kind) +
                static_cast<int>(base[index_of(x + 1, y + 1)] == kind);
            const int diagonal_b =
                static_cast<int>(base[index_of(x + 1, y - 1)] == kind) +
                static_cast<int>(base[index_of(x - 1, y + 1)] == kind);

            const int strongest =
                std::max(std::max(horizontal, vertical), std::max(diagonal_a, diagonal_b));
            if (strongest < 2) {
                place_fault(x, y, kind);
                continue;
            }

            const uint32_t hash = static_cast<uint32_t>(
                (x * 73856093) ^ (y * 19349663) ^ (static_cast<int>(kind) * 83492791));
            int offset = 0;
            if ((hash & 3U) != 0U) {
                offset = ((hash >> 2U) & 1U) == 0U ? -1 : 1;
            }
            const int drift = static_cast<int>((hash >> 4U) % 3U) - 1;

            int perp_dx = 0;
            int perp_dy = 0;
            int along_dx = 0;
            int along_dy = 0;
            if (strongest == horizontal) {
                perp_dy = 1;
                along_dx = 1;
            } else if (strongest == vertical) {
                perp_dx = 1;
                along_dy = 1;
            } else if (strongest == diagonal_a) {
                perp_dx = 1;
                perp_dy = -1;
                along_dx = 1;
                along_dy = 1;
            } else {
                perp_dx = 1;
                perp_dy = 1;
                along_dx = 1;
                along_dy = -1;
            }

            place_fault(x + perp_dx * offset + along_dx * drift,
                        y + perp_dy * offset + along_dy * drift, kind);
            if ((hash & 7U) < 5U) {
                place_fault(x, y, kind);
            }
        }
    }

    for (int y = 1; y + 1 < height; ++y) {
        for (int x = 1; x + 1 < width; ++x) {
            const size_t index = index_of(x, y);
            if (warped[index] != FaultKind::None) {
                continue;
            }

            const FaultKind horizontal = warped[index_of(x - 1, y)];
            if (horizontal != FaultKind::None && horizontal == warped[index_of(x + 1, y)]) {
                warped[index] = horizontal;
                continue;
            }

            const FaultKind vertical = warped[index_of(x, y - 1)];
            if (vertical != FaultKind::None && vertical == warped[index_of(x, y + 1)]) {
                warped[index] = vertical;
                continue;
            }

            const FaultKind diagonal_a = warped[index_of(x - 1, y - 1)];
            if (diagonal_a != FaultKind::None &&
                diagonal_a == warped[index_of(x + 1, y + 1)]) {
                warped[index] = diagonal_a;
                continue;
            }

            const FaultKind diagonal_b = warped[index_of(x - 1, y + 1)];
            if (diagonal_b != FaultKind::None &&
                diagonal_b == warped[index_of(x + 1, y - 1)]) {
                warped[index] = diagonal_b;
            }
        }
    }

    fault_map.swap(warped);
}

void produce_image_gray(float* heightmap, int width, int height, const char* filename)
{
    writeImageGray(filename, width, height, heightmap, "FOO");
}

void produce_image_colors(float* heightmap, int width, int height, const char* filename)
{
    writeImageColors(filename, width, height, heightmap, "FOO");
}

void produce_image_gray_with_faults(float* heightmap, int width, int height, const char* filename,
                                    const FaultKind* fault_map)
{
    writeImageGrayWithFaultMap(filename, width, height, heightmap, fault_map, "FOO");
}

void produce_image_colors_with_faults(float* heightmap, int width, int height, const char* filename,
                                      const FaultKind* fault_map)
{
    writeImageColorsWithFaultMap(filename, width, height, heightmap, fault_map, "FOO");
}

void build_output_filename(const char* filename, const char* suffix, char* destination,
                           size_t destination_size)
{
    const char* extension = strrchr(filename, '.');
    if (extension != nullptr) {
        snprintf(destination, destination_size, "%.*s%s%s",
                 static_cast<int>(extension - filename), filename, suffix, extension);
    } else {
        snprintf(destination, destination_size, "%s%s.png", filename, suffix);
    }
}

void ensure_parent_directory(const char* filename)
{
    if (filename == nullptr) {
        return;
    }

    const std::filesystem::path output_path(filename);
    const std::filesystem::path parent = output_path.parent_path();
    if (!parent.empty()) {
        std::filesystem::create_directories(parent);
    }
}

std::string replace_extension(const char* filename, const char* extension)
{
    std::filesystem::path output_path(filename);
    if (output_path.has_extension()) {
        output_path.replace_extension(extension);
    } else {
        output_path += extension;
    }
    return output_path.string();
}

std::string escape_powershell_single_quoted(const std::string& value)
{
    std::string escaped;
    escaped.reserve(value.size());
    for (const char ch : value) {
        escaped.push_back(ch);
        if (ch == '\'') {
            escaped.push_back('\'');
        }
    }
    return escaped;
}

int write_gif_with_powershell(const std::vector<std::string>& frame_files,
                              const std::string& output_file, uint32_t delay_cs)
{
#ifndef _WIN32
    (void)frame_files;
    (void)output_file;
    (void)delay_cs;
    fprintf(stderr, "GIF export is only supported on Windows in this build.\n");
    return 1;
#else
    if (frame_files.empty()) {
        return 0;
    }

    ensure_parent_directory(output_file.c_str());

    const std::filesystem::path output_path(output_file);
    const std::filesystem::path script_path = output_path.string() + ".encode.ps1";
    std::ofstream script(script_path, std::ios::binary | std::ios::trunc);
    if (!script) {
        fprintf(stderr, "Could not create GIF encoder script %s\n",
                script_path.string().c_str());
        return 1;
    }

    script << "$ErrorActionPreference = 'Stop'\n";
    script << "Add-Type -AssemblyName PresentationCore\n";
    script << "$files = @(\n";
    for (const std::string& frame : frame_files) {
        script << "  '" << escape_powershell_single_quoted(frame) << "'\n";
    }
    script << ")\n";
    script << "$encoder = New-Object System.Windows.Media.Imaging.GifBitmapEncoder\n";
    script << "$delay = [UInt16]" << std::max<uint32_t>(1U, delay_cs) << "\n";
    script << "foreach($file in $files){\n";
    script << "  $stream = [System.IO.File]::OpenRead($file)\n";
    script << "  try {\n";
    script << "    $decoder = New-Object System.Windows.Media.Imaging.PngBitmapDecoder($stream,[System.Windows.Media.Imaging.BitmapCreateOptions]::PreservePixelFormat,[System.Windows.Media.Imaging.BitmapCacheOption]::OnLoad)\n";
    script << "    $frame = $decoder.Frames[0]\n";
    script << "    $meta = New-Object System.Windows.Media.Imaging.BitmapMetadata 'gif'\n";
    script << "    $meta.SetQuery('/grctlext/Delay', $delay)\n";
    script << "    $meta.SetQuery('/grctlext/Disposal', [Byte]2)\n";
    script << "    $newFrame = [System.Windows.Media.Imaging.BitmapFrame]::Create($frame, $frame.Thumbnail, $meta, $frame.ColorContexts)\n";
    script << "    $encoder.Frames.Add($newFrame)\n";
    script << "  } finally {\n";
    script << "    $stream.Dispose()\n";
    script << "  }\n";
    script << "}\n";
    script << "$fs = [System.IO.File]::Create('" << escape_powershell_single_quoted(output_file)
           << "')\n";
    script << "try { $encoder.Save($fs) } finally { $fs.Dispose() }\n";
    script.close();

    const std::string command = "powershell -NoProfile -ExecutionPolicy Bypass -File \"" +
                                script_path.string() + "\"";
    const int result = std::system(command.c_str());
    std::error_code ec;
    std::filesystem::remove(script_path, ec);
    if (result != 0) {
        fprintf(stderr, "GIF export failed for %s\n", output_file.c_str());
        return 1;
    }

    return 0;
#endif
}

void append_saved_frames(const char* filename, bool show_fault_lines, const FaultKind* fault_map,
                         std::vector<std::string>& plain_frames,
                         std::vector<std::string>& fault_frames)
{
    plain_frames.emplace_back(filename);
    if (!show_fault_lines || fault_map == nullptr) {
        return;
    }

    char faults_filename[260];
    build_output_filename(filename, "_faults", faults_filename, sizeof(faults_filename));
    fault_frames.emplace_back(faults_filename);
}

void capture_fault_snapshot(void* p, int width, int height, FaultSnapshot& snapshot)
{
    const uint32_t plate_count = platec_api_get_plate_count(p);
    if (plate_count == 0) {
        snapshot.valid = false;
        snapshot.platesmap.clear();
        snapshot.plate_vx.clear();
        snapshot.plate_vy.clear();
        snapshot.plate_count = 0;
        return;
    }

    const size_t pixel_count = static_cast<size_t>(width) * static_cast<size_t>(height);
    snapshot.platesmap.resize(pixel_count);
    memcpy(snapshot.platesmap.data(), platec_api_get_platesmap(p), pixel_count * sizeof(uint32_t));

    snapshot.plate_vx.resize(plate_count);
    snapshot.plate_vy.resize(plate_count);
    for (uint32_t i = 0; i < plate_count; ++i) {
        snapshot.plate_vx[i] = platec_api_velocity_unity_vector_x(p, i);
        snapshot.plate_vy[i] = platec_api_velocity_unity_vector_y(p, i);
    }

    snapshot.plate_count = plate_count;
    snapshot.valid = true;
}

bool build_fault_map_from_snapshot(const FaultSnapshot& snapshot, int width, int height,
                                   std::vector<FaultKind>& fault_map)
{
    if (!snapshot.valid) {
        return false;
    }

    return buildFaultMap(width, height, snapshot.platesmap.data(), snapshot.plate_vx.data(),
                         snapshot.plate_vy.data(), snapshot.plate_count, &fault_map) == 0;
}

bool build_plate_boundary_map_from_snapshot(const FaultSnapshot& snapshot, int width, int height,
                                            std::vector<uint8_t>& boundary_map)
{
    if (!snapshot.valid) {
        return false;
    }

    const size_t pixel_count = static_cast<size_t>(width) * static_cast<size_t>(height);
    auto wrapped_index_of = [width, height](int x, int y) -> size_t {
        const int wrapped_x = (x % width + width) % width;
        const int wrapped_y = (y % height + height) % height;
        return static_cast<size_t>(wrapped_y) * static_cast<size_t>(width) +
               static_cast<size_t>(wrapped_x);
    };
    boundary_map.assign(pixel_count, 0);
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            const size_t index = wrapped_index_of(x, y);
            const uint32_t plate = snapshot.platesmap[index];
            if (plate >= snapshot.plate_count) {
                continue;
            }

            const uint32_t left = snapshot.platesmap[wrapped_index_of(x - 1, y)];
            const uint32_t right = snapshot.platesmap[wrapped_index_of(x + 1, y)];
            const uint32_t top = snapshot.platesmap[wrapped_index_of(x, y - 1)];
            const uint32_t bottom = snapshot.platesmap[wrapped_index_of(x, y + 1)];

            if ((left < snapshot.plate_count && left != plate) ||
                (right < snapshot.plate_count && right != plate) ||
                (top < snapshot.plate_count && top != plate) ||
                (bottom < snapshot.plate_count && bottom != plate)) {
                boundary_map[index] = 1;
            }
        }
    }

    return true;
}

void accumulate_faults(const FaultSnapshot& snapshot, int width, int height, FaultAccumulator& acc)
{
    std::vector<FaultKind> fault_map;
    if (!build_fault_map_from_snapshot(snapshot, width, height, fault_map)) {
        return;
    }

    const size_t pixel_count = static_cast<size_t>(width) * static_cast<size_t>(height);
    if (!acc.valid) {
        acc.convergent.assign(pixel_count, 0);
        acc.divergent.assign(pixel_count, 0);
        acc.transform.assign(pixel_count, 0);
        acc.valid = true;
    }

    for (size_t i = 0; i < pixel_count; ++i) {
        switch (fault_map[i]) {
        case FaultKind::Convergent:
            ++acc.convergent[i];
            break;
        case FaultKind::Divergent:
            ++acc.divergent[i];
            break;
        case FaultKind::Transform:
            ++acc.transform[i];
            break;
        case FaultKind::None:
        default:
            break;
        }
    }

    ++acc.samples;
}

bool build_render_fault_map(const FaultAccumulator& acc, int width, int height,
                            std::vector<FaultKind>& fault_map)
{
    if (!acc.valid || acc.samples == 0) {
        return false;
    }

    const size_t pixel_count = static_cast<size_t>(width) * static_cast<size_t>(height);
    std::vector<FaultKind> dominant(pixel_count, FaultKind::None);
    std::vector<uint32_t> dominant_count(pixel_count, 0);

    uint32_t max_count = 0;
    for (size_t i = 0; i < pixel_count; ++i) {
        const uint32_t convergent = acc.convergent[i];
        const uint32_t divergent = acc.divergent[i];
        const uint32_t transform = acc.transform[i];
        const uint32_t total = convergent + divergent + transform;

        uint32_t best = convergent;
        FaultKind kind = FaultKind::Convergent;
        if (divergent > best) {
            best = divergent;
            kind = FaultKind::Divergent;
        }
        if (transform > best) {
            best = transform;
            kind = FaultKind::Transform;
        }

        if (best == 0 || (total > 0 && best * 3U < total * 2U)) {
            continue;
        }

        dominant[i] = kind;
        dominant_count[i] = best;
        max_count = std::max(max_count, best);
    }

    if (max_count == 0) {
        return false;
    }

    const uint32_t activity_threshold = std::max(1U, max_count / 10U);
    fault_map.assign(pixel_count, FaultKind::None);

    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            const size_t index = static_cast<size_t>(y) * static_cast<size_t>(width) +
                                 static_cast<size_t>(x);
            if (dominant[index] == FaultKind::None || dominant_count[index] < activity_threshold) {
                continue;
            }

            uint32_t local_max = dominant_count[index];
            for (int ny = std::max(0, y - 1); ny <= std::min(height - 1, y + 1); ++ny) {
                for (int nx = std::max(0, x - 1); nx <= std::min(width - 1, x + 1); ++nx) {
                    const size_t neighbour =
                        static_cast<size_t>(ny) * static_cast<size_t>(width) +
                        static_cast<size_t>(nx);
                    if (dominant[neighbour] == dominant[index]) {
                        local_max = std::max(local_max, dominant_count[neighbour]);
                    }
                }
            }

            const uint32_t hash = static_cast<uint32_t>(
                (x * 73856093) ^ (y * 19349663) ^
                (static_cast<int>(dominant[index]) * 83492791));
            if (dominant_count[index] == local_max ||
                (dominant_count[index] + 1U == local_max && (hash & 1U) == 0U)) {
                fault_map[index] = dominant[index];
            }
        }
    }

    prune_fault_clusters(fault_map, width, height);
    roughen_fault_map(fault_map, width, height);
    return true;
}

uint32_t noise_hash(int x, int y, uint32_t seed)
{
    uint32_t h = static_cast<uint32_t>(x) * 374761393U + static_cast<uint32_t>(y) * 668265263U +
                 seed * 2246822519U;
    h = (h ^ (h >> 13U)) * 1274126177U;
    return h ^ (h >> 16U);
}

float signed_noise(int x, int y, uint32_t seed)
{
    const uint32_t h = noise_hash(x, y, seed);
    return (static_cast<float>(h & 1023U) / 511.5f) - 1.0f;
}

size_t wrapped_index(int x, int y, int width, int height)
{
    const int wrapped_x = (x % width + width) % width;
    const int wrapped_y = (y % height + height) % height;
    return static_cast<size_t>(wrapped_y) * static_cast<size_t>(width) +
           static_cast<size_t>(wrapped_x);
}

float sample_wrapped(const std::vector<float>& heightmap, int x, int y, int width, int height)
{
    return heightmap[wrapped_index(x, y, width, height)];
}

float sample_wrapped_bilinear(const std::vector<float>& heightmap, float x, float y, int width,
                              int height)
{
    const float floor_x = std::floor(x);
    const float floor_y = std::floor(y);
    const int x0 = static_cast<int>(floor_x);
    const int y0 = static_cast<int>(floor_y);
    const int x1 = x0 + 1;
    const int y1 = y0 + 1;
    const float tx = x - floor_x;
    const float ty = y - floor_y;

    const float s00 = sample_wrapped(heightmap, x0, y0, width, height);
    const float s10 = sample_wrapped(heightmap, x1, y0, width, height);
    const float s01 = sample_wrapped(heightmap, x0, y1, width, height);
    const float s11 = sample_wrapped(heightmap, x1, y1, width, height);

    const float top = s00 * (1.0f - tx) + s10 * tx;
    const float bottom = s01 * (1.0f - tx) + s11 * tx;
    return top * (1.0f - ty) + bottom * ty;
}

float smoothstep01(float value)
{
    const float clamped = std::clamp(value, 0.0f, 1.0f);
    return clamped * clamped * (3.0f - 2.0f * clamped);
}

std::vector<uint8_t> build_fault_distance_map(const FaultKind* fault_map, const uint8_t* seam_map,
                                              int width, int height, int max_radius)
{
    const size_t pixel_count = static_cast<size_t>(width) * static_cast<size_t>(height);
    std::vector<uint8_t> distances(pixel_count, static_cast<uint8_t>(max_radius + 1));
    if (fault_map == nullptr && seam_map == nullptr) {
        return distances;
    }

    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            const size_t index = wrapped_index(x, y, width, height);
            if ((fault_map != nullptr && fault_map[index] != FaultKind::None) ||
                (seam_map != nullptr && seam_map[index] != 0)) {
                distances[index] = 0;
                continue;
            }

            int best = max_radius + 1;
            for (int dy = -max_radius; dy <= max_radius && best > 0; ++dy) {
                for (int dx = -max_radius; dx <= max_radius; ++dx) {
                    const size_t neighbour = wrapped_index(x + dx, y + dy, width, height);
                    const bool has_fault =
                        fault_map != nullptr && fault_map[neighbour] != FaultKind::None;
                    const bool has_seam = seam_map != nullptr && seam_map[neighbour] != 0;
                    if (!has_fault && !has_seam) {
                        continue;
                    }

                    const int distance = std::max(std::abs(dx), std::abs(dy));
                    if (distance < best) {
                        best = distance;
                    }
                }
            }

            distances[index] = static_cast<uint8_t>(best);
        }
    }

    return distances;
}

void postprocess_imported_heightmap(float* heightmap, int width, int height,
                                    const FaultKind* fault_map, const uint8_t* seam_map,
                                    const float* original_heightmap)
{
    const size_t pixel_count = static_cast<size_t>(width) * static_cast<size_t>(height);
    std::vector<float> current(heightmap, heightmap + pixel_count);
    std::vector<float> next(pixel_count, 0.0f);
    const int fault_radius = 6;
    const int seam_radius = 4;
    const std::vector<uint8_t> fault_distance_map =
        build_fault_distance_map(fault_map, nullptr, width, height, fault_radius);
    const std::vector<uint8_t> seam_distance_map =
        build_fault_distance_map(nullptr, seam_map, width, height, seam_radius);
    std::vector<float> original;
    std::vector<float> original_blur;

    for (int iteration = 0; iteration < 3; ++iteration) {
        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                const size_t index = wrapped_index(x, y, width, height);
                float weighted_sum = 0.0f;
                float total_weight = 0.0f;
                for (int dy = -1; dy <= 1; ++dy) {
                    for (int dx = -1; dx <= 1; ++dx) {
                        const float weight = (dx == 0 && dy == 0) ? 4.0f :
                                             ((dx == 0 || dy == 0) ? 2.0f : 1.0f);
                        weighted_sum += sample_wrapped(current, x + dx, y + dy, width, height) * weight;
                        total_weight += weight;
                    }
                }

                const float average = weighted_sum / total_weight;
                float strength = (current[index] < 0.62f) ? 0.14f : 0.04f;
                const uint8_t fault_distance = fault_distance_map[index];
                const uint8_t seam_distance = seam_distance_map[index];
                const uint8_t distance = std::min(fault_distance, seam_distance);
                if (fault_distance <= 1) {
                    strength = 0.26f;
                } else if (seam_distance <= 1) {
                    strength = 0.18f;
                } else if (distance <= 3) {
                    strength = 0.18f;
                }

                next[index] = current[index] * (1.0f - strength) + average * strength;
            }
        }

        current.swap(next);
    }

    if (original_heightmap != nullptr) {
        original.assign(original_heightmap, original_heightmap + pixel_count);
        original_blur.assign(pixel_count, 0.0f);

        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                float weighted_sum = 0.0f;
                float total_weight = 0.0f;
                for (int dy = -2; dy <= 2; ++dy) {
                    for (int dx = -2; dx <= 2; ++dx) {
                        const float weight =
                            (dx == 0 && dy == 0) ? 6.0f :
                            ((dx == 0 || dy == 0) ? 3.0f : 1.0f);
                        weighted_sum += sample_wrapped(original, x + dx, y + dy, width, height) *
                                        weight;
                        total_weight += weight;
                    }
                }

                original_blur[wrapped_index(x, y, width, height)] = weighted_sum / total_weight;
            }
        }

        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                const size_t index = wrapped_index(x, y, width, height);
                const float detail = original[index] - original_blur[index];
                const float land_factor = std::clamp((current[index] - 0.42f) / 0.30f, 0.0f, 1.0f);
                const uint8_t distance =
                    std::min(fault_distance_map[index], seam_distance_map[index]);
                const float fault_factor =
                    (distance > fault_radius) ? 1.0f :
                    smoothstep01((static_cast<float>(distance) - 1.0f) /
                                 static_cast<float>(fault_radius - 1));
                current[index] =
                    std::clamp(current[index] + detail * 0.85f * land_factor * fault_factor, 0.0f,
                               1.0f);
            }
        }
    }

    if (fault_map != nullptr || seam_map != nullptr) {
        std::vector<float> warped = current;
        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                const size_t index = wrapped_index(x, y, width, height);
                const uint8_t fault_distance = fault_distance_map[index];
                const uint8_t seam_distance = seam_distance_map[index];
                const float fault_influence =
                    (fault_distance > fault_radius)
                        ? 0.0f
                        : smoothstep01((static_cast<float>(fault_radius + 1) -
                                        static_cast<float>(fault_distance)) /
                                       static_cast<float>(fault_radius + 1));
                const float seam_influence =
                    (seam_distance > seam_radius)
                        ? 0.0f
                        : smoothstep01((static_cast<float>(seam_radius + 1) -
                                        static_cast<float>(seam_distance)) /
                                       static_cast<float>(seam_radius + 1));
                if (fault_influence <= 0.0f && seam_influence <= 0.0f) {
                    continue;
                }

                const float jitter_x = 3.6f * signed_noise(x / 19, y / 19, 11U) +
                                       1.8f * signed_noise(x / 7, y / 7, 29U) +
                                       0.8f * signed_noise(x / 3, y / 3, 53U);
                const float jitter_y = 3.6f * signed_noise(x / 19, y / 19, 17U) +
                                       1.8f * signed_noise(x / 7, y / 7, 37U) +
                                       0.8f * signed_noise(x / 3, y / 3, 61U);
                const float warp_scale = 1.1f * seam_influence + 2.6f * fault_influence;
                const float sampled = sample_wrapped_bilinear(
                    current, static_cast<float>(x) + jitter_x * warp_scale,
                    static_cast<float>(y) + jitter_y * warp_scale, width, height);
                float relief_boost = 0.0f;
                if (!original.empty()) {
                    relief_boost = (original[index] - original_blur[index]) *
                                   (0.10f * seam_influence + 0.24f * fault_influence);
                }
                const float ridge_noise =
                    (0.016f * signed_noise(x / 4, y / 4, 71U) +
                     0.010f * signed_noise(x / 9, y / 9, 89U)) * fault_influence;
                const float blend = 0.18f * seam_influence + 0.34f * fault_influence;
                warped[index] = std::clamp(current[index] * (1.0f - blend) +
                                           sampled * blend + relief_boost +
                                           ridge_noise,
                                           0.0f, 1.0f);
            }
        }
        current.swap(warped);
    }

    if (fault_map != nullptr) {
        std::vector<float> softened = current;
        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                const size_t index = wrapped_index(x, y, width, height);
                if (fault_distance_map[index] > 1) {
                    continue;
                }

                float weighted_sum = 0.0f;
                float total_weight = 0.0f;
                for (int dy = -1; dy <= 1; ++dy) {
                    for (int dx = -1; dx <= 1; ++dx) {
                        const float weight = (dx == 0 && dy == 0) ? 4.0f :
                                             ((dx == 0 || dy == 0) ? 2.0f : 1.0f);
                        weighted_sum += sample_wrapped(current, x + dx, y + dy, width, height) *
                                        weight;
                        total_weight += weight;
                    }
                }

                const float average = weighted_sum / total_weight;
                softened[index] = current[index] * 0.82f + average * 0.18f;
            }
        }
        current.swap(softened);
    }

    memcpy(heightmap, current.data(), pixel_count * sizeof(float));
}

void save_image(void* p, const char* filename, const int width, const int height, bool colors,
                bool show_fault_lines, const FaultKind* fault_map, const uint8_t* seam_map,
                bool smooth_imported_output,
                const float* original_heightmap = nullptr)
{
    ensure_parent_directory(filename);
    const float* heightmap = platec_api_get_heightmap(p);
    float* copy = new float[width * height];
    memcpy(copy, heightmap, sizeof(float) * width * height);
    normalize(copy, width * height);
    if (smooth_imported_output) {
        postprocess_imported_heightmap(copy, width, height, fault_map, seam_map,
                                       original_heightmap);
        normalize(copy, width * height);
    }

    if (colors) {
        produce_image_colors(copy, width, height, filename);
    } else {
        produce_image_gray(copy, width, height, filename);
    }

    if (show_fault_lines && fault_map != nullptr) {
        char faults_filename[260];
        build_output_filename(filename, "_faults", faults_filename, sizeof(faults_filename));
        ensure_parent_directory(faults_filename);
        if (colors) {
            produce_image_colors_with_faults(copy, width, height, faults_filename, fault_map);
        } else {
            produce_image_gray_with_faults(copy, width, height, faults_filename, fault_map);
        }
    }

    delete[] copy;
}

} // namespace

typedef struct {
    uint32_t seed;
    uint32_t width;
    uint32_t height;
    bool colors;
    bool seed_set;
    bool dimensions_set;
    bool show_fault_lines;
    bool cycle_count_set;
    bool erosion_period_set;
    bool sea_level_set;
    bool folding_ratio_set;
    bool aggr_overlap_abs_set;
    bool aggr_overlap_rel_set;
    bool gif;
    bool num_plates_set;
    char* filename;
    char* input_image;
    float sea_level;
    float folding_ratio;
    float aggr_overlap_rel;
    uint32_t erosion_period;
    uint32_t aggr_overlap_abs;
    uint32_t cycle_count;
    uint32_t gif_delay_cs;
    uint32_t num_plates;
    uint32_t progress_interval;
    uint32_t step;
} Params;

char DEFAULT_FILENAME[] = "simulation";

void fill_params(Params& params, int argc, char* argv[])
{
    params.seed = 1;
    params.width = 600;
    params.height = 400;
    params.colors = true;
    params.seed_set = false;
    params.dimensions_set = false;
    params.show_fault_lines = false;
    params.cycle_count_set = false;
    params.erosion_period_set = false;
    params.sea_level_set = false;
    params.folding_ratio_set = false;
    params.aggr_overlap_abs_set = false;
    params.aggr_overlap_rel_set = false;
    params.gif = false;
    params.num_plates_set = false;
    params.filename = DEFAULT_FILENAME;
    params.input_image = nullptr;
    params.sea_level = 0.65f;
    params.folding_ratio = 0.02f;
    params.aggr_overlap_rel = 0.33f;
    params.erosion_period = 60;
    params.aggr_overlap_abs = 1000000;
    params.cycle_count = 2;
    params.gif_delay_cs = 10;
    params.num_plates = 10;
    params.progress_interval = 0;
    params.step = 0;

    int p = 1;
    while (p < argc) {
        if (0 == strcmp(argv[p], "--help") || 0 == strcmp(argv[p], "-h")) {
            printf(" -h --help           : show this message\n");
            printf(" -s SEED             : use the given SEED\n");
            printf(" --seed SEED         : same as -s\n");
            printf(" --dim WIDTH HEIGHT  : use the given width and height\n");
            printf(" --num-plates N      : use N tectonic plates (default: 10, image mode: 4)\n");
            printf(" --cycle-count N     : use N tectonic cycles (default: 2, image mode: 1)\n");
            printf(" --erosion-period N  : erode every N iterations (0 disables erosion)\n");
            printf(" --sea-level F       : target ocean coverage between 0 and 1\n");
            printf(" --folding-ratio F   : collision uplift ratio\n");
            printf(" --aggr-overlap-abs N: absolute aggregation threshold\n");
            printf(" --aggr-overlap-rel F: relative aggregation threshold\n");
            printf(" --colors            : generate a colors map\n");
            printf(" --grayscale         : generate a grayscale map\n");
            printf(" --filename FILENAME : generated map are named with the given filename (the extension is appended)\n");
            printf(" --input-image FILE  : use a PNG image as the initial terrain (RGB converted to luminance)\n");
            printf(" --show-fault-lines  : also write *_faults.png with convergent/divergent/transform overlays\n");
            printf(" --gif               : write an animated GIF from the saved simulation frames (defaults to --step 20)\n");
            printf(" --gif-delay N       : GIF frame delay in centiseconds (default: 10)\n");
            printf(" --progress-interval N: print progress every N simulation steps\n");
            printf(" --step X            : generate intermediate maps any given steps\n");
            exit(0);
        } else if (0 == strcmp(argv[p], "-s") || 0 == strcmp(argv[p], "--seed")) {
            if (p + 1 >= argc) {
                printf("error: a parameter should follow %s\n", argv[p]);
                exit(1);
            }
            long seed = atol(argv[p+1]);
            if (seed==0) {
                printf("error: not a number\n");
                exit(1);
            }
            params.seed = seed;
            params.seed_set = true;
            p += 2;
        } else if (0 == strcmp(argv[p], "--dim")) {
            if (p + 2 >= argc) {
                printf("error: two parameters should follow --dim\n");
                exit(1);
            }
            int width = atoi(argv[p+1]);
            int height = atoi(argv[p+2]);
            if (width==0 || height==0) {
                printf("error: not a number\n");
                exit(1);
            }
            if (width<5 || height<5) {
                printf("error: dimensions have to be positive and >= 5\n");
                exit(1);
            }
            params.width = width;
            params.height = height;
            params.dimensions_set = true;
            p += 3;
        } else if (0 == strcmp(argv[p], "--num-plates")) {
            if (p + 1 >= argc) {
                printf("error: a parameter should follow --num-plates\n");
                exit(1);
            }
            int num_plates = atoi(argv[p+1]);
            if (num_plates < 2) {
                printf("error: number of plates must be >= 2\n");
                exit(1);
            }
            params.num_plates = static_cast<uint32_t>(num_plates);
            params.num_plates_set = true;
            p += 2;
        } else if (0 == strcmp(argv[p], "--cycle-count")) {
            if (p + 1 >= argc) {
                printf("error: a parameter should follow --cycle-count\n");
                exit(1);
            }
            int cycle_count = atoi(argv[p+1]);
            if (cycle_count < 1) {
                printf("error: cycle count must be >= 1\n");
                exit(1);
            }
            params.cycle_count = static_cast<uint32_t>(cycle_count);
            params.cycle_count_set = true;
            p += 2;
        } else if (0 == strcmp(argv[p], "--erosion-period")) {
            if (p + 1 >= argc) {
                printf("error: a parameter should follow --erosion-period\n");
                exit(1);
            }
            int erosion_period = atoi(argv[p+1]);
            if (erosion_period < 0) {
                printf("error: erosion period must be >= 0\n");
                exit(1);
            }
            params.erosion_period = static_cast<uint32_t>(erosion_period);
            params.erosion_period_set = true;
            p += 2;
        } else if (0 == strcmp(argv[p], "--sea-level")) {
            if (p + 1 >= argc) {
                printf("error: a parameter should follow --sea-level\n");
                exit(1);
            }
            float sea_level = static_cast<float>(atof(argv[p+1]));
            if (sea_level < 0.0f || sea_level > 1.0f) {
                printf("error: sea level must be between 0 and 1\n");
                exit(1);
            }
            params.sea_level = sea_level;
            params.sea_level_set = true;
            p += 2;
        } else if (0 == strcmp(argv[p], "--folding-ratio")) {
            if (p + 1 >= argc) {
                printf("error: a parameter should follow --folding-ratio\n");
                exit(1);
            }
            float folding_ratio = static_cast<float>(atof(argv[p+1]));
            if (folding_ratio < 0.0f || folding_ratio > 1.0f) {
                printf("error: folding ratio must be between 0 and 1\n");
                exit(1);
            }
            params.folding_ratio = folding_ratio;
            params.folding_ratio_set = true;
            p += 2;
        } else if (0 == strcmp(argv[p], "--aggr-overlap-abs")) {
            if (p + 1 >= argc) {
                printf("error: a parameter should follow --aggr-overlap-abs\n");
                exit(1);
            }
            int aggr_overlap_abs = atoi(argv[p+1]);
            if (aggr_overlap_abs < 0) {
                printf("error: aggr overlap abs must be >= 0\n");
                exit(1);
            }
            params.aggr_overlap_abs = static_cast<uint32_t>(aggr_overlap_abs);
            params.aggr_overlap_abs_set = true;
            p += 2;
        } else if (0 == strcmp(argv[p], "--aggr-overlap-rel")) {
            if (p + 1 >= argc) {
                printf("error: a parameter should follow --aggr-overlap-rel\n");
                exit(1);
            }
            float aggr_overlap_rel = static_cast<float>(atof(argv[p+1]));
            if (aggr_overlap_rel < 0.0f || aggr_overlap_rel > 1.0f) {
                printf("error: aggr overlap rel must be between 0 and 1\n");
                exit(1);
            }
            params.aggr_overlap_rel = aggr_overlap_rel;
            params.aggr_overlap_rel_set = true;
            p += 2;
        } else if (0 == strcmp(argv[p], "--colors")) {
            params.colors = true;
            p += 1;
        } else if (0 == strcmp(argv[p], "--grayscale")) {
            params.colors = false;
            p += 1;
        } else if (0 == strcmp(argv[p], "--show-fault-lines")) {
            params.show_fault_lines = true;
            p += 1;
        } else if (0 == strcmp(argv[p], "--gif")) {
            params.gif = true;
            p += 1;
        } else if (0 == strcmp(argv[p], "--gif-delay")) {
            if (p + 1 >= argc) {
                printf("error: a parameter should follow --gif-delay\n");
                exit(1);
            }
            int gif_delay_cs = atoi(argv[p+1]);
            if (gif_delay_cs < 1) {
                printf("error: gif delay must be >= 1 centisecond\n");
                exit(1);
            }
            params.gif_delay_cs = static_cast<uint32_t>(gif_delay_cs);
            p += 2;
        } else if (0 == strcmp(argv[p], "--filename")) {
            if (p + 1 >= argc) {
                printf("error: a parameter should follow --filename\n");
                exit(1);
            }
            params.filename = argv[p+1];
            p += 2;
        } else if (0 == strcmp(argv[p], "--input-image")) {
            if (p + 1 >= argc) {
                printf("error: a parameter should follow --input-image\n");
                exit(1);
            }
            params.input_image = argv[p+1];
            p += 2;
        } else if (0 == strcmp(argv[p], "--step")) {
            if (p + 1 >= argc) {
                printf("error: a parameter should follow --step\n");
                exit(1);
            }
            int step = atoi(argv[p+1]);
            if (step==0) {
                printf("error: not a number\n");
                exit(1);
            }
            if (step<0) {
                printf("error: step have to be positive\n");
                exit(1);
            }
            params.step = step;
            p += 2;
        } else if (0 == strcmp(argv[p], "--progress-interval")) {
            if (p + 1 >= argc) {
                printf("error: a parameter should follow --progress-interval\n");
                exit(1);
            }
            int progress_interval = atoi(argv[p+1]);
            if (progress_interval < 1) {
                printf("error: progress interval must be >= 1\n");
                exit(1);
            }
            params.progress_interval = static_cast<uint32_t>(progress_interval);
            p += 2;
        } else {
            printf("Unexpected param '%s' use -h to display a list of params\n", argv[p]);
            exit(1);
        }
    }
}

/// Should take several parameters:
/// - colors/grayscale
/// - width and height
/// - if to generate the intermediate imaged and how frequently
/// - the output filename
/// - the seed
int main(int argc, char* argv[])
{
    Params params;
    fill_params(params, argc, argv);
    if (params.gif && params.step == 0) {
        params.step = 20;
    }
    const bool track_faults = params.show_fault_lines || params.input_image != nullptr;

    std::vector<float> input_heightmap;
    if (params.input_image != nullptr) {
        uint32_t image_width = 0;
        uint32_t image_height = 0;
        if (loadHeightmapFromImage(params.input_image, &image_width, &image_height,
                                 &input_heightmap) != 0) {
            return 1;
        }

        if (params.dimensions_set) {
            if (params.width != image_width || params.height != image_height) {
                printf("error: input image dimensions (%u x %u) do not match --dim (%u x %u)\n",
                       image_width, image_height, params.width, params.height);
                return 1;
            }
        } else {
            params.width = image_width;
            params.height = image_height;
        }

        if (!params.num_plates_set) {
            params.num_plates = 4;
        }
        if (!params.cycle_count_set) {
            params.cycle_count = 1;
        }
        if (!params.erosion_period_set) {
            params.erosion_period = 0;
        }
    }

    printf("Plate-tectonics simulation example\n");
    printf(" seed     : %d\n", params.seed);
    printf(" width    : %d\n", params.width);
    printf(" height   : %d\n", params.height);
    printf(" map      : %s\n", params.colors ? "colors" : "grayscale");
    printf(" filename : %s\n", params.filename);
    printf(" plates   : %u\n", params.num_plates);
    printf(" cycles   : %u\n", params.cycle_count);
    printf(" erosion  : %u\n", params.erosion_period);
    printf(" sea      : %.2f\n", params.sea_level);
    printf(" folding  : %.3f\n", params.folding_ratio);
    if (params.input_image != nullptr)
        printf(" input    : %s\n", params.input_image);
    if (params.show_fault_lines)
        printf(" faults   : convergent=red, divergent=cyan, transform=magenta\n");
    if (params.progress_interval != 0)
        printf(" progress : every %u steps\n", params.progress_interval);
    if (params.gif)
        printf(" gif      : yes (%u cs)\n", params.gif_delay_cs);
    if (params.step == 0)
        printf(" step     : no\n");
    else
        printf(" step     : %i\n", params.step);

    printf("\n");

    void* p = nullptr;
    if (params.input_image != nullptr) {
        p = platec_api_create_from_heightmap(params.seed, params.width, params.height,
                                             input_heightmap.data(), params.sea_level,
                                             params.erosion_period, params.folding_ratio,
                                             params.aggr_overlap_abs, params.aggr_overlap_rel,
                                             params.cycle_count, params.num_plates);
    } else {
        p = platec_api_create(params.seed, params.width, params.height, params.sea_level,
                              params.erosion_period, params.folding_ratio,
                              params.aggr_overlap_abs, params.aggr_overlap_rel,
                              params.cycle_count, params.num_plates);
    }

    FaultSnapshot initial_faults;
    FaultSnapshot latest_faults;
    FaultAccumulator fault_accumulator;
    std::vector<FaultKind> initial_fault_map;
    std::vector<FaultKind> rendered_fault_map;
    std::vector<uint8_t> initial_seam_map;
    std::vector<uint8_t> rendered_seam_map;
    std::vector<std::string> plain_frames;
    std::vector<std::string> fault_frames;
    if (track_faults) {
        capture_fault_snapshot(p, static_cast<int>(params.width), static_cast<int>(params.height),
                               initial_faults);
        latest_faults = initial_faults;
        build_fault_map_from_snapshot(initial_faults, static_cast<int>(params.width),
                                      static_cast<int>(params.height), initial_fault_map);
        build_plate_boundary_map_from_snapshot(initial_faults, static_cast<int>(params.width),
                                               static_cast<int>(params.height), initial_seam_map);
        accumulate_faults(initial_faults, static_cast<int>(params.width),
                          static_cast<int>(params.height), fault_accumulator);
    }

    char filenamei[250];
    build_output_filename(params.filename, "_initial", filenamei, sizeof(filenamei));
    save_image(p, filenamei, static_cast<int>(params.width), static_cast<int>(params.height),
               params.colors, params.show_fault_lines,
               initial_fault_map.empty() ? nullptr : initial_fault_map.data(),
               initial_seam_map.empty() ? nullptr : initial_seam_map.data(), false,
               params.input_image != nullptr ? input_heightmap.data() : nullptr);
    append_saved_frames(filenamei, params.show_fault_lines,
                        initial_fault_map.empty() ? nullptr : initial_fault_map.data(),
                        plain_frames, fault_frames);
    printf(" * initial map created\n");

    int step = 0;
    while (platec_api_is_finished(p) == 0) {
        step++;
        platec_api_step(p);

        if (params.progress_interval != 0 &&
            (step % static_cast<int>(params.progress_interval) == 0)) {
            printf(" * progress step=%i iter=%u cycle=%u plates=%u\n", step,
                   platec_api_get_iteration_count(p), platec_api_get_cycle_count(p),
                   platec_api_get_plate_count(p));
        }

        if (track_faults) {
            FaultSnapshot current_faults;
            capture_fault_snapshot(p, static_cast<int>(params.width), static_cast<int>(params.height),
                                   current_faults);
            if (current_faults.valid) {
                latest_faults = std::move(current_faults);
                accumulate_faults(latest_faults, static_cast<int>(params.width),
                                  static_cast<int>(params.height), fault_accumulator);
            }
        }

        if (params.step != 0 && (step % params.step == 0) ) {
            char filename[250];
            char suffix[64];
            snprintf(suffix, sizeof(suffix), "_%05i", step);
            build_output_filename(params.filename, suffix, filename, sizeof(filename));
            const FaultKind* fault_map = nullptr;
            const uint8_t* seam_map = nullptr;
            rendered_fault_map.clear();
            rendered_seam_map.clear();
            if (track_faults) {
                if (build_render_fault_map(fault_accumulator, static_cast<int>(params.width),
                                           static_cast<int>(params.height), rendered_fault_map) &&
                    !rendered_fault_map.empty()) {
                    fault_map = rendered_fault_map.data();
                }
                if (build_plate_boundary_map_from_snapshot(
                        latest_faults, static_cast<int>(params.width),
                        static_cast<int>(params.height), rendered_seam_map) &&
                    !rendered_seam_map.empty()) {
                    seam_map = rendered_seam_map.data();
                }
            }
            printf(" * step %i (filename %s)\n", step, filename);
            save_image(p, filename, static_cast<int>(params.width), static_cast<int>(params.height),
                       params.colors, params.show_fault_lines, fault_map, seam_map,
                       params.input_image != nullptr,
                       params.input_image != nullptr ? input_heightmap.data() : nullptr);
            append_saved_frames(filename, params.show_fault_lines, fault_map, plain_frames,
                                fault_frames);
        }
    }

    char filename[250];
    build_output_filename(params.filename, "", filename, sizeof(filename));
    const FaultKind* final_fault_map = nullptr;
    const uint8_t* final_seam_map = nullptr;
    rendered_fault_map.clear();
    rendered_seam_map.clear();
    if (track_faults &&
        build_render_fault_map(fault_accumulator, static_cast<int>(params.width),
                               static_cast<int>(params.height), rendered_fault_map) &&
        !rendered_fault_map.empty()) {
        final_fault_map = rendered_fault_map.data();
    }
    if (track_faults &&
        build_plate_boundary_map_from_snapshot(latest_faults, static_cast<int>(params.width),
                                               static_cast<int>(params.height),
                                               rendered_seam_map) &&
        !rendered_seam_map.empty()) {
        final_seam_map = rendered_seam_map.data();
    }
    save_image(p, filename, static_cast<int>(params.width), static_cast<int>(params.height),
               params.colors, params.show_fault_lines, final_fault_map, final_seam_map,
               params.input_image != nullptr,
               params.input_image != nullptr ? input_heightmap.data() : nullptr);
    append_saved_frames(filename, params.show_fault_lines, final_fault_map, plain_frames,
                        fault_frames);
    if (params.gif) {
        const std::string gif_filename = replace_extension(params.filename, ".gif");
        if (write_gif_with_powershell(plain_frames, gif_filename, params.gif_delay_cs) != 0) {
            platec_api_destroy(p);
            return 1;
        }
        if (params.show_fault_lines && !fault_frames.empty()) {
            char faults_gif[260];
            build_output_filename(gif_filename.c_str(), "_faults", faults_gif, sizeof(faults_gif));
            if (write_gif_with_powershell(fault_frames, faults_gif, params.gif_delay_cs) != 0) {
                platec_api_destroy(p);
                return 1;
            }
        }
    }
    printf(" * simulation completed (filename %s)\n", filename);
    platec_api_destroy(p);
    return 0;
}
