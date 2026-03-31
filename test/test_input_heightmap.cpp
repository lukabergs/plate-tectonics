/******************************************************************************
 *  plate-tectonics, a plate tectonics simulation library
 *  Copyright (C) 2012-2013 Lauri Viitanen
 *  Copyright (C) 2014-2015 Federico Tomassetti, Bret Curtis
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, see http://www.gnu.org/licenses/
 *****************************************************************************/

#include "gtest/gtest.h"
#include "map_drawing.hpp"
#include "platecapi.hpp"

#include <cmath>
#include <cstring>
#include <cstdio>
#include <vector>

namespace {

float quantized_gray(float value) {
    return static_cast<float>(static_cast<int>(value * 255.0f)) / 255.0f;
}

bool write_linear_grayscale_png(const char* filename, uint32_t width, uint32_t height,
                                const std::vector<float>& source) {
    std::vector<png_byte> pixels(static_cast<std::size_t>(width) * static_cast<std::size_t>(height) * 3U);
    for (std::size_t i = 0; i < source.size(); ++i) {
        const png_byte value = static_cast<png_byte>(quantized_gray(source[i]) * 255.0f);
        pixels[i * 3U + 0U] = value;
        pixels[i * 3U + 1U] = value;
        pixels[i * 3U + 2U] = value;
    }

    png_image image;
    std::memset(&image, 0, sizeof(image));
    image.version = PNG_IMAGE_VERSION;
    image.width = width;
    image.height = height;
    image.format = PNG_FORMAT_RGB;
    const bool ok = png_image_write_to_file(&image, filename, 0, pixels.data(), 0, nullptr) != 0;
    png_image_free(&image);
    return ok;
}

} // namespace

TEST(InputHeightmap, CreateFromHeightmapUsesProvidedTerrain)
{
    const uint32_t width = 8;
    const uint32_t height = 8;
    std::vector<float> input(width * height, 0.0f);
    for (uint32_t y = 0; y < height; ++y) {
        for (uint32_t x = width / 2; x < width; ++x) {
            input[y * width + x] = 1.0f;
        }
    }

    void* p = platec_api_create_from_heightmap(7, width, height, input.data(), 0.5f, 60, 0.02f,
                                               1000000, 0.33f, 2, 4);
    ASSERT_NE(p, nullptr);

    const float* heightmap = platec_api_get_heightmap(p);
    ASSERT_NE(heightmap, nullptr);
    EXPECT_NEAR(0.02f, heightmap[0], 0.001f);
    EXPECT_NEAR(0.02f, heightmap[width / 2 - 1], 0.001f);
    EXPECT_FLOAT_EQ(2.0f, heightmap[width / 2]);
    EXPECT_FLOAT_EQ(2.0f, heightmap[width - 1]);

    platec_api_destroy(p);
}

TEST(InputHeightmap, LoadImageReadsBackGrayscaleLuminance)
{
    const uint32_t width = 5;
    const uint32_t height = 5;
    const char* filename = "test_input_heightmap.png";
    std::vector<float> source = {
        0.0f, 0.25f, 0.5f, 0.75f, 1.0f,
        1.0f, 0.75f, 0.5f, 0.25f, 0.0f,
        0.0f, 0.5f, 1.0f, 0.5f, 0.0f,
        0.1f, 0.2f, 0.3f, 0.4f, 0.5f,
        0.9f, 0.8f, 0.7f, 0.6f, 0.5f,
    };

    ASSERT_TRUE(write_linear_grayscale_png(filename, width, height, source));

    uint32_t loaded_width = 0;
    uint32_t loaded_height = 0;
    std::vector<float> loaded;
    ASSERT_EQ(0, loadHeightmapFromImage(filename, &loaded_width, &loaded_height, &loaded));

    EXPECT_EQ(width, loaded_width);
    EXPECT_EQ(height, loaded_height);
    ASSERT_EQ(source.size(), loaded.size());

    for (std::size_t i = 0; i < source.size(); ++i) {
        EXPECT_NEAR(quantized_gray(source[i]), loaded[i], 1.0f / 255.0f);
    }

    std::remove(filename);
}
