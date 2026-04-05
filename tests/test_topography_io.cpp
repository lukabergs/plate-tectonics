#include "gtest/gtest.h"
#include "heightmap_io.hpp"

#include <filesystem>
#include <vector>

namespace fs = std::filesystem;

namespace {

struct ScopedDelete {
    explicit ScopedDelete(fs::path _path) : path(std::move(_path)) {}
    ~ScopedDelete()
    {
        std::error_code ec;
        fs::remove(path, ec);
    }

    fs::path path;
};

} // namespace

TEST(TopographyIo, RawR32RoundTrip)
{
    const std::vector<float> samples = {0.0f, 1.0f, -1.5f, 1024.25f, 32768.5f, 65535.0f};
    const fs::path path = fs::current_path() / "topography_roundtrip.r32";
    ScopedDelete cleanup(path);

    ASSERT_EQ(0, writeRawR32(path.string().c_str(), samples.data(), samples.size()));

    std::vector<float> loaded;
    ASSERT_EQ(0, readRawR32(path.string().c_str(), loaded, samples.size()));
    EXPECT_EQ(samples, loaded);
}

TEST(TopographyIo, Png32RoundTrip)
{
    const int width = 3;
    const int height = 2;
    const std::vector<float> samples = {0.0f, 1.0f, -0.25f, 1024.5f, 32768.25f, 65535.0f};
    const fs::path path = fs::current_path() / "topography_roundtrip.png";
    ScopedDelete cleanup(path);

    ASSERT_EQ(0, writeImageRgba32(path.string().c_str(), width, height, samples.data(), "roundtrip"));

    std::vector<float> loaded;
    int loaded_width = 0;
    int loaded_height = 0;
    ASSERT_EQ(0, readImageRgba32(path.string().c_str(), loaded, loaded_width, loaded_height));
    EXPECT_EQ(width, loaded_width);
    EXPECT_EQ(height, loaded_height);
    EXPECT_EQ(samples, loaded);
}

TEST(TopographyIo, MetadataRoundTrip)
{
    const fs::path path = fs::current_path() / "topography_roundtrip.r16.json";
    ScopedDelete cleanup(path);

    TopographyCodec::Metadata metadata;
    metadata.width = 64;
    metadata.height = 32;
    metadata.sea_level_m = 31000;
    metadata.max_height_m = TopographyCodec::kMaxHeightMeters;
    metadata.version = TopographyCodec::kMetadataVersion;
    metadata.format = TopographyCodec::kMetricFormatR16;
    metadata.endianness = TopographyCodec::kLittleEndian;
    metadata.layout = TopographyCodec::kRowMajorLayout;

    ASSERT_EQ(0, writeTopographyMetadataJson(path.string().c_str(), metadata));

    TopographyCodec::Metadata loaded;
    ASSERT_EQ(0, readTopographyMetadataJson(path.string().c_str(), loaded));
    EXPECT_EQ(metadata.width, loaded.width);
    EXPECT_EQ(metadata.height, loaded.height);
    EXPECT_EQ(metadata.sea_level_m, loaded.sea_level_m);
    EXPECT_EQ(metadata.max_height_m, loaded.max_height_m);
    EXPECT_EQ(metadata.version, loaded.version);
    EXPECT_EQ(metadata.format, loaded.format);
    EXPECT_EQ(metadata.endianness, loaded.endianness);
    EXPECT_EQ(metadata.layout, loaded.layout);
}
