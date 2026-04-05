#include "heightmap_io.hpp"
#include "utils.hpp"
#include <cctype>
#include <fstream>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using namespace std;

namespace {

bool extract_json_number(const std::string& content, const char* key, uint32_t& value)
{
    const std::string token = std::string("\"") + key + "\"";
    size_t pos = content.find(token);
    if (pos == std::string::npos) {
        return false;
    }

    pos = content.find(':', pos + token.size());
    if (pos == std::string::npos) {
        return false;
    }
    ++pos;
    while (pos < content.size() && std::isspace(static_cast<unsigned char>(content[pos])) != 0) {
        ++pos;
    }

    size_t end = pos;
    while (end < content.size() && std::isdigit(static_cast<unsigned char>(content[end])) != 0) {
        ++end;
    }

    if (end == pos) {
        return false;
    }

    char* parse_end = nullptr;
    const unsigned long parsed = std::strtoul(content.substr(pos, end - pos).c_str(), &parse_end, 10);
    if (parse_end == nullptr || *parse_end != '\0' ||
        parsed > static_cast<unsigned long>(std::numeric_limits<uint32_t>::max())) {
        return false;
    }

    value = static_cast<uint32_t>(parsed);
    return true;
}

bool extract_json_string(const std::string& content, const char* key, std::string& value)
{
    const std::string token = std::string("\"") + key + "\"";
    size_t pos = content.find(token);
    if (pos == std::string::npos) {
        return false;
    }

    pos = content.find(':', pos + token.size());
    if (pos == std::string::npos) {
        return false;
    }
    ++pos;
    while (pos < content.size() && std::isspace(static_cast<unsigned char>(content[pos])) != 0) {
        ++pos;
    }
    if (pos >= content.size() || content[pos] != '"') {
        return false;
    }
    ++pos;

    const size_t end = content.find('"', pos);
    if (end == std::string::npos) {
        return false;
    }

    value = content.substr(pos, end - pos);
    return true;
}

} // namespace

inline void setGray(png_byte *ptr, int val)
{
    ptr[0] = static_cast<png_byte>(val);
    ptr[1] = static_cast<png_byte>(val);
    ptr[2] = static_cast<png_byte>(val);
}

inline void setColor(png_byte *ptr, png_byte r, png_byte g, png_byte b)
{
    ptr[0] = r;
    ptr[1] = g;
    ptr[2] = b;
}

void gradient(png_byte *ptr, png_byte ra, png_byte ga, png_byte ba, png_byte rb, png_byte gb,
              png_byte bb, float h, float ha, float hb);

struct ColorStop {
    float position;
    png_byte r;
    png_byte g;
    png_byte b;
};

void spectralColor(png_byte* ptr, float h)
{
    static const ColorStop kStops[] = {
        {0.00f,  0,   8,  92},
        {0.12f,  0,  84, 196},
        {0.24f,  0, 232, 255},
        {0.38f,  0, 168,  74},
        {0.52f, 244, 235,  60},
        {0.66f, 255, 164,  32},
        {0.78f, 228,  42,  24},
        {0.88f, 118,  72,  42},
        {0.95f, 132,  86, 178},
        {1.00f, 255, 255, 255},
    };

    if (h <= kStops[0].position) {
        setColor(ptr, kStops[0].r, kStops[0].g, kStops[0].b);
        return;
    }

    for (size_t i = 1; i < sizeof(kStops) / sizeof(kStops[0]); ++i) {
        if (h <= kStops[i].position) {
            gradient(ptr, kStops[i - 1].r, kStops[i - 1].g, kStops[i - 1].b, kStops[i].r,
                     kStops[i].g, kStops[i].b, h, kStops[i - 1].position, kStops[i].position);
            return;
        }
    }

    setColor(ptr, 255, 255, 255);
}

int writeImageRgb(const char* filename, int width, int height, const png_byte* rgb, const char* title)
{
    volatile int code = 0;
    FILE * volatile fp = nullptr;
    png_structp volatile png_ptr = nullptr;
    png_infop volatile info_ptr = nullptr;

#ifdef _WIN32
    FILE* fp_temp = nullptr;
    errno_t err = fopen_s(&fp_temp, filename, "wb");
    fp = fp_temp;
    if (err != 0 || fp == nullptr) {
#else
    fp = fopen(filename, "wb");
    if (fp == nullptr) {
#endif
        fprintf(stderr, "Could not open file %s for writing\n", filename);
        code = 1;
        goto finalise;
    }

    png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
    if (png_ptr == nullptr) {
        fprintf(stderr, "Could not allocate write struct\n");
        code = 1;
        goto finalise;
    }

    info_ptr = png_create_info_struct(png_ptr);
    if (info_ptr == nullptr) {
        fprintf(stderr, "Could not allocate info struct\n");
        code = 1;
        goto finalise;
    }

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable: 4611)
#endif
    if (setjmp(png_jmpbuf(png_ptr))) {
        fprintf(stderr, "Error during png creation\n");
        code = 1;
        goto finalise;
    }
#ifdef _MSC_VER
#pragma warning(pop)
#endif

    png_init_io(png_ptr, fp);
    png_set_IHDR(png_ptr, info_ptr, width, height,
                 8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
                 PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

    if (title != nullptr) {
        png_text title_text;
        title_text.compression = PNG_TEXT_COMPRESSION_NONE;
        title_text.key = const_cast<char*>("Title");
        title_text.text = const_cast<char*>(title);
        png_set_text(png_ptr, info_ptr, &title_text, 1);
    }

    png_write_info(png_ptr, info_ptr);

    for (int y = 0; y < height; ++y) {
        png_write_row(png_ptr, const_cast<png_bytep>(rgb + static_cast<size_t>(y) * static_cast<size_t>(width) * 3U));
    }
    png_write_end(png_ptr, nullptr);

finalise:
    if (fp != nullptr) fclose(fp);
    if (png_ptr != nullptr) {
        if (info_ptr != nullptr) {
            png_free_data(png_ptr, info_ptr, PNG_FREE_ALL, -1);
        }
        {
            png_structp png_ptr_nv = png_ptr;
            png_destroy_write_struct(&png_ptr_nv, (png_infopp)nullptr);
        }
    }

    return code;
}

int writeImage(const char* filename, int width, int height, float *heightmap, const char* title,
               void (drawFunction)(png_structp&, png_bytep&, int, int, float*))
{
    volatile int code = 0;
    FILE * volatile fp = nullptr;
    png_structp volatile png_ptr = nullptr;
    png_infop volatile info_ptr = nullptr;
    png_bytep volatile row = nullptr;
    size_t row_bytes = 0;  // Declare early to avoid goto issues

    // Open file for writing (binary mode)
#ifdef _WIN32
    // fopen_s doesn't accept volatile pointer, so use a non-volatile temporary
    FILE* fp_temp = nullptr;
    errno_t err = fopen_s(&fp_temp, filename, "wb");
    fp = fp_temp;
    if (err != 0 || fp == nullptr) {
#else
    fp = fopen(filename, "wb");
    if (fp == nullptr) {
#endif
        fprintf(stderr, "Could not open file %s for writing\n", filename);
        code = 1;
        goto finalise;
    }

    // Initialize write structure
    png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
    if (png_ptr == nullptr) {
        fprintf(stderr, "Could not allocate write struct\n");
        code = 1;
        goto finalise;
    }

    // Initialize info structure
    info_ptr = png_create_info_struct(png_ptr);
    if (info_ptr == nullptr) {
        fprintf(stderr, "Could not allocate info struct\n");
        code = 1;
        goto finalise;
    }

    // Setup Exception handling
#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable: 4611)  // Interaction between setjmp and C++ object destruction is acceptable here
#endif
    if (setjmp(png_jmpbuf(png_ptr))) {
        fprintf(stderr, "Error during png creation\n");
        code = 1;
        goto finalise;
    }
#ifdef _MSC_VER
#pragma warning(pop)
#endif

    png_init_io(png_ptr, fp);

    // Write header (8-bit colour depth)
    png_set_IHDR(png_ptr, info_ptr, width, height,
                 8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
                 PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

    // Set title
    if (title != nullptr) {
        png_text title_text;
        title_text.compression = PNG_TEXT_COMPRESSION_NONE;
        title_text.key = const_cast<char*>("Title");
        title_text.text = (char*)title;
        png_set_text(png_ptr, info_ptr, &title_text, 1);
    }

    png_write_info(png_ptr, info_ptr);

    // Allocate memory for one row (3 bytes per pixel - RGB)
    row_bytes = 3 * width * sizeof(png_byte);
    row = (png_bytep) malloc(row_bytes);

    if (row == nullptr) {
        fprintf(stderr, "Could not allocate memory for one row\n");
        code = 1;
        goto finalise;
    }

    // Write image data
    // Need to create non-volatile references for function calls
    {
        png_structp png_ptr_nv = png_ptr;
        png_bytep row_nv = row;
        drawFunction(png_ptr_nv, row_nv, width, height, heightmap);
    }

    // End write
    png_write_end(png_ptr, nullptr);

finalise:
    if (fp != nullptr) fclose(fp);
    if (row != nullptr) free(row);
    if (png_ptr != nullptr) {
        if (info_ptr != nullptr) {
            png_free_data(png_ptr, info_ptr, PNG_FREE_ALL, -1);
        }
        {
            png_structp png_ptr_nv = png_ptr;
            png_destroy_write_struct(&png_ptr_nv, (png_infopp)nullptr);
        }
    }

    return code;
}

float find_value_for_quantile(const float quantile, const float* array, const uint32_t size)
{
    float value = 0.5;
    float th_step = 0.5;

    while (th_step > 0.00001)
    {
        uint32_t count = 0;
        for (uint32_t i = 0; i < size; ++i)
            count += (array[i] < value);

        th_step *= 0.5;
        if (count / (float)size < quantile)
            value += th_step;
        else
            value -= th_step;
    }
    return value;
}

void gradient(png_byte *ptr, png_byte ra, png_byte ga, png_byte ba, png_byte rb, png_byte gb, png_byte bb, float h, float ha, float hb)
{
    if (ha>hb) {
        printf("BAD1\n");
        throw runtime_error("BAD1\n");
    }
    if (hb<h) {
        printf("BAD2\n");
        throw runtime_error("BAD2\n");
    }
    if (ha>h) {
        printf("BAD3\n");
        throw runtime_error("BAD3\n");
    }
    float h_delta = hb - ha;
    float simil_b = (h - ha)/h_delta;
    float simil_a = (1.0f - simil_b);
    setColor(ptr,
             static_cast<png_byte>((float)simil_a * ra + (float)simil_b * rb),
             static_cast<png_byte>((float)simil_a * ga + (float)simil_b * gb),
             static_cast<png_byte>((float)simil_a * ba + (float)simil_b * bb));
}

void drawGrayImage(png_structp& png_ptr, png_bytep& row, int width, int height, float *heightmap)
{
    int x, y;
    for (y=0 ; y<height ; y++) {
        for (x=0 ; x<width ; x++) {

            float h = heightmap[(y*width + x)];
            float res = 0.0f;
            if (h <= 0.0f) {
                res = 0;
            } else if (h >= 1.0f) {
                res = 255;
            } else {
                res = (h * 255.0f);
            }

            setGray(&(row[x*3]), static_cast<int>(res));
        }
        png_write_row(png_ptr, row);
    }
}

void renderImageGrayRgb(int width, int height, float* heightmap, std::vector<png_byte>& rgb)
{
    rgb.resize(static_cast<size_t>(width) * static_cast<size_t>(height) * 3U);
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            const float h = heightmap[y * width + x];
            float res = 0.0f;
            if (h <= 0.0f) {
                res = 0.0f;
            } else if (h >= 1.0f) {
                res = 255.0f;
            } else {
                res = h * 255.0f;
            }

            png_byte* ptr = &rgb[(static_cast<size_t>(y) * static_cast<size_t>(width) + static_cast<size_t>(x)) * 3U];
            setGray(ptr, static_cast<int>(res));
        }
    }
}

void drawColorsImage(png_structp& png_ptr, png_bytep& row, int width, int height, float *heightmap)
{
    int x, y;
    for (y=0 ; y<height ; y++) {
        for (x=0 ; x<width ; x++) {
            float h = heightmap[(y*width + x)];
            spectralColor(&(row[x*3]), h);
        }
        png_write_row(png_ptr, row);
    }
}

void renderImageColorsRgb(int width, int height, float* heightmap, std::vector<png_byte>& rgb)
{
    rgb.resize(static_cast<size_t>(width) * static_cast<size_t>(height) * 3U);

    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            const float h = heightmap[y * width + x];
            png_byte* ptr = &rgb[(static_cast<size_t>(y) * static_cast<size_t>(width) + static_cast<size_t>(x)) * 3U];
            spectralColor(ptr, h);
        }
    }
}

int writeImageGray(const char* filename, int width, int height, float *heightmap, const char* title)
{
    return writeImage(filename, width, height, heightmap, title, drawGrayImage);
}

int writeImageColors(const char* filename, int width, int height, float *heightmap, const char* title)
{
    return writeImage(filename, width, height, heightmap, title, drawColorsImage);
}

int writeImageGray16(const char* filename, int width, int height, const uint16_t* heightmap,
                     const char* title)
{
    volatile int code = 0;
    FILE * volatile fp = nullptr;
    png_structp volatile png_ptr = nullptr;
    png_infop volatile info_ptr = nullptr;
    std::vector<png_byte> row;

#ifdef _WIN32
    FILE* fp_temp = nullptr;
    errno_t err = fopen_s(&fp_temp, filename, "wb");
    fp = fp_temp;
    if (err != 0 || fp == nullptr) {
#else
    fp = fopen(filename, "wb");
    if (fp == nullptr) {
#endif
        fprintf(stderr, "Could not open file %s for writing\n", filename);
        code = 1;
        goto finalise;
    }

    png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
    if (png_ptr == nullptr) {
        fprintf(stderr, "Could not allocate write struct\n");
        code = 1;
        goto finalise;
    }

    info_ptr = png_create_info_struct(png_ptr);
    if (info_ptr == nullptr) {
        fprintf(stderr, "Could not allocate info struct\n");
        code = 1;
        goto finalise;
    }

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable: 4611)
#endif
    if (setjmp(png_jmpbuf(png_ptr))) {
        fprintf(stderr, "Error during png creation\n");
        code = 1;
        goto finalise;
    }
#ifdef _MSC_VER
#pragma warning(pop)
#endif

    png_init_io(png_ptr, fp);
    png_set_IHDR(png_ptr, info_ptr, width, height, 16, PNG_COLOR_TYPE_GRAY,
                 PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

    if (title != nullptr) {
        png_text title_text;
        title_text.compression = PNG_TEXT_COMPRESSION_NONE;
        title_text.key = const_cast<char*>("Title");
        title_text.text = const_cast<char*>(title);
        png_set_text(png_ptr, info_ptr, &title_text, 1);
    }

    png_write_info(png_ptr, info_ptr);
    row.resize(static_cast<size_t>(width) * 2U);
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            const uint16_t sample =
                heightmap[static_cast<size_t>(y) * static_cast<size_t>(width) + static_cast<size_t>(x)];
            row[static_cast<size_t>(x) * 2U] = static_cast<png_byte>((sample >> 8) & 0xFFU);
            row[static_cast<size_t>(x) * 2U + 1U] = static_cast<png_byte>(sample & 0xFFU);
        }
        png_write_row(png_ptr, row.data());
    }
    png_write_end(png_ptr, nullptr);

finalise:
    if (fp != nullptr) fclose(fp);
    if (png_ptr != nullptr) {
        if (info_ptr != nullptr) {
            png_free_data(png_ptr, info_ptr, PNG_FREE_ALL, -1);
        }
        {
            png_structp png_ptr_nv = png_ptr;
            png_destroy_write_struct(&png_ptr_nv, (png_infopp)nullptr);
        }
    }

    return code;
}

int readImageGray16(const char* filename, std::vector<uint16_t>& heightmap, int& width, int& height)
{
    volatile int code = 1;
    FILE * volatile fp = nullptr;
    png_structp volatile png_ptr = nullptr;
    png_infop volatile info_ptr = nullptr;
    png_bytep volatile image_data = nullptr;
    png_bytep* volatile rows = nullptr;
    png_uint_32 png_width = 0;
    png_uint_32 png_height = 0;
    int bit_depth = 0;
    int color_type = 0;
    png_size_t row_bytes = 0;

#ifdef _WIN32
    FILE* fp_temp = nullptr;
    errno_t err = fopen_s(&fp_temp, filename, "rb");
    fp = fp_temp;
    if (err != 0 || fp == nullptr) {
#else
    fp = fopen(filename, "rb");
    if (fp == nullptr) {
#endif
        fprintf(stderr, "Could not open file %s for reading\n", filename);
        goto finalise;
    }

    png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
    if (png_ptr == nullptr) {
        fprintf(stderr, "Could not allocate read struct\n");
        goto finalise;
    }

    info_ptr = png_create_info_struct(png_ptr);
    if (info_ptr == nullptr) {
        fprintf(stderr, "Could not allocate info struct\n");
        goto finalise;
    }

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable: 4611)
#endif
    if (setjmp(png_jmpbuf(png_ptr))) {
        fprintf(stderr, "Error during png read\n");
        goto finalise;
    }
#ifdef _MSC_VER
#pragma warning(pop)
#endif

    png_init_io(png_ptr, fp);
    png_read_info(png_ptr, info_ptr);

    png_width = png_get_image_width(png_ptr, info_ptr);
    png_height = png_get_image_height(png_ptr, info_ptr);
    bit_depth = png_get_bit_depth(png_ptr, info_ptr);
    color_type = png_get_color_type(png_ptr, info_ptr);

    if (color_type == PNG_COLOR_TYPE_GRAY_ALPHA) {
        png_set_strip_alpha(png_ptr);
        color_type = PNG_COLOR_TYPE_GRAY;
    }

    if (bit_depth != 16 || color_type != PNG_COLOR_TYPE_GRAY) {
        fprintf(stderr, "Expected 16-bit grayscale PNG in %s\n", filename);
        goto finalise;
    }

    png_read_update_info(png_ptr, info_ptr);

    width = static_cast<int>(png_width);
    height = static_cast<int>(png_height);
    heightmap.resize(static_cast<size_t>(width) * static_cast<size_t>(height));

    row_bytes = png_get_rowbytes(png_ptr, info_ptr);
    image_data = static_cast<png_bytep>(malloc(row_bytes * png_height));
    rows = static_cast<png_bytep*>(malloc(sizeof(png_bytep) * png_height));
    if (image_data == nullptr || rows == nullptr) {
        fprintf(stderr, "Could not allocate PNG16 buffers\n");
        goto finalise;
    }

    for (png_uint_32 y = 0; y < png_height; ++y) {
        rows[y] = image_data + y * row_bytes;
    }

    {
        png_structp png_ptr_nv = png_ptr;
        png_infop info_ptr_nv = info_ptr;
        png_bytep* rows_nv = rows;
        png_read_image(png_ptr_nv, rows_nv);
        png_read_end(png_ptr_nv, info_ptr_nv);
    }

    for (int y = 0; y < height; ++y) {
        const png_bytep row = rows[y];
        for (int x = 0; x < width; ++x) {
            const size_t index =
                static_cast<size_t>(y) * static_cast<size_t>(width) + static_cast<size_t>(x);
            heightmap[index] =
                static_cast<uint16_t>((static_cast<uint16_t>(row[x * 2]) << 8) |
                                      static_cast<uint16_t>(row[x * 2 + 1]));
        }
    }

    code = 0;

finalise:
    if (fp != nullptr) fclose(fp);
    if (image_data != nullptr) free(const_cast<png_bytep>(image_data));
    if (rows != nullptr) free(const_cast<png_bytep*>(rows));
    if (png_ptr != nullptr) {
        if (info_ptr != nullptr) {
            png_free_data(png_ptr, info_ptr, PNG_FREE_ALL, -1);
        }
        {
            png_structp png_ptr_nv = png_ptr;
            png_destroy_read_struct(&png_ptr_nv, (png_infopp)nullptr, (png_infopp)nullptr);
        }
    }
    if (code != 0) {
        heightmap.clear();
        width = 0;
        height = 0;
    }
    return code;
}

int writeRawR16(const char* filename, const uint16_t* heightmap, size_t sample_count)
{
    std::ofstream output(filename, std::ios::binary);
    if (!output) {
        return 1;
    }

    for (size_t i = 0; i < sample_count; ++i) {
        const char bytes[2] = {
            static_cast<char>(heightmap[i] & 0xFFU),
            static_cast<char>((heightmap[i] >> 8) & 0xFFU),
        };
        output.write(bytes, sizeof(bytes));
        if (!output) {
            return 1;
        }
    }

    return 0;
}

int readRawR16(const char* filename, std::vector<uint16_t>& heightmap, size_t expected_sample_count)
{
    std::ifstream input(filename, std::ios::binary);
    if (!input) {
        return 1;
    }

    input.seekg(0, std::ios::end);
    const std::streamoff size = input.tellg();
    input.seekg(0, std::ios::beg);
    if (size != static_cast<std::streamoff>(expected_sample_count * sizeof(uint16_t))) {
        return 1;
    }

    heightmap.resize(expected_sample_count);
    for (size_t i = 0; i < expected_sample_count; ++i) {
        unsigned char bytes[2] = {0, 0};
        input.read(reinterpret_cast<char*>(bytes), sizeof(bytes));
        if (!input) {
            heightmap.clear();
            return 1;
        }
        heightmap[i] = static_cast<uint16_t>(bytes[0] | (static_cast<uint16_t>(bytes[1]) << 8));
    }

    return 0;
}

int writeTopographyMetadataJson(const char* filename, const TopographyCodec::Metadata& metadata)
{
    std::ofstream output(filename, std::ios::binary);
    if (!output) {
        return 1;
    }

    output << "{\n";
    output << "  \"width\": " << metadata.width << ",\n";
    output << "  \"height\": " << metadata.height << ",\n";
    output << "  \"sea_level_m\": " << metadata.sea_level_m << ",\n";
    output << "  \"max_height_m\": " << metadata.max_height_m << ",\n";
    output << "  \"format\": \"" << metadata.format << "\",\n";
    output << "  \"endianness\": \"" << metadata.endianness << "\",\n";
    output << "  \"layout\": \"" << metadata.layout << "\",\n";
    output << "  \"version\": " << metadata.version << "\n";
    output << "}\n";

    return output ? 0 : 1;
}

int readTopographyMetadataJson(const char* filename, TopographyCodec::Metadata& metadata)
{
    std::ifstream input(filename, std::ios::binary);
    if (!input) {
        return 1;
    }

    const std::string content((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
    uint32_t width = 0;
    uint32_t height = 0;
    uint32_t sea_level_m = 0;
    uint32_t max_height_m = 0;
    uint32_t version = 0;
    std::string format;
    std::string endianness;
    std::string layout;

    if (!extract_json_number(content, "width", width) ||
        !extract_json_number(content, "height", height) ||
        !extract_json_number(content, "sea_level_m", sea_level_m) ||
        !extract_json_number(content, "max_height_m", max_height_m) ||
        !extract_json_number(content, "version", version) ||
        !extract_json_string(content, "format", format) ||
        !extract_json_string(content, "endianness", endianness) ||
        !extract_json_string(content, "layout", layout) ||
        sea_level_m > TopographyCodec::kMaxHeightMeters ||
        max_height_m > TopographyCodec::kMaxHeightMeters) {
        return 1;
    }

    metadata.width = width;
    metadata.height = height;
    metadata.sea_level_m = static_cast<uint16_t>(sea_level_m);
    metadata.max_height_m = static_cast<uint16_t>(max_height_m);
    metadata.version = version;
    metadata.format = format;
    metadata.endianness = endianness;
    metadata.layout = layout;
    return 0;
}

int readImageNormalized(const char* filename, std::vector<float>& heightmap, int& width, int& height)
{
    volatile int code = 1;
    FILE * volatile fp = nullptr;
    png_structp volatile png_ptr = nullptr;
    png_infop volatile info_ptr = nullptr;
    png_bytep volatile image_data = nullptr;
    png_bytep* volatile rows = nullptr;
    png_uint_32 png_width = 0;
    png_uint_32 png_height = 0;
    int bit_depth = 0;
    int color_type = 0;
    png_size_t row_bytes = 0;

#ifdef _WIN32
    FILE* fp_temp = nullptr;
    errno_t err = fopen_s(&fp_temp, filename, "rb");
    fp = fp_temp;
    if (err != 0 || fp == nullptr) {
#else
    fp = fopen(filename, "rb");
    if (fp == nullptr) {
#endif
        fprintf(stderr, "Could not open file %s for reading\n", filename);
        goto finalise;
    }

    png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
    if (png_ptr == nullptr) {
        fprintf(stderr, "Could not allocate read struct\n");
        goto finalise;
    }

    info_ptr = png_create_info_struct(png_ptr);
    if (info_ptr == nullptr) {
        fprintf(stderr, "Could not allocate info struct\n");
        goto finalise;
    }

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable: 4611)
#endif
    if (setjmp(png_jmpbuf(png_ptr))) {
        fprintf(stderr, "Error during png read\n");
        goto finalise;
    }
#ifdef _MSC_VER
#pragma warning(pop)
#endif

    png_init_io(png_ptr, fp);
    png_read_info(png_ptr, info_ptr);

    png_width = png_get_image_width(png_ptr, info_ptr);
    png_height = png_get_image_height(png_ptr, info_ptr);
    bit_depth = png_get_bit_depth(png_ptr, info_ptr);
    color_type = png_get_color_type(png_ptr, info_ptr);

    if (bit_depth == 16) {
        png_set_strip_16(png_ptr);
    }
    if (color_type == PNG_COLOR_TYPE_PALETTE) {
        png_set_palette_to_rgb(png_ptr);
    }
    if (color_type == PNG_COLOR_TYPE_GRAY && bit_depth < 8) {
        png_set_expand_gray_1_2_4_to_8(png_ptr);
    }
    if (png_get_valid(png_ptr, info_ptr, PNG_INFO_tRNS)) {
        png_set_tRNS_to_alpha(png_ptr);
    }
    if (color_type == PNG_COLOR_TYPE_GRAY || color_type == PNG_COLOR_TYPE_GRAY_ALPHA) {
        png_set_gray_to_rgb(png_ptr);
    }
    if ((color_type & PNG_COLOR_MASK_ALPHA) != 0 || png_get_valid(png_ptr, info_ptr, PNG_INFO_tRNS)) {
        png_set_strip_alpha(png_ptr);
    }

    png_read_update_info(png_ptr, info_ptr);

    row_bytes = png_get_rowbytes(png_ptr, info_ptr);
    image_data = static_cast<png_bytep>(malloc(row_bytes * png_height));
    rows = static_cast<png_bytep*>(malloc(sizeof(png_bytep) * png_height));

    if (image_data == nullptr || rows == nullptr) {
        fprintf(stderr, "Could not allocate memory for image data\n");
        goto finalise;
    }

    for (png_uint_32 y = 0; y < png_height; ++y) {
        rows[y] = image_data + y * row_bytes;
    }

    {
        png_structp png_ptr_nv = png_ptr;
        png_infop info_ptr_nv = info_ptr;
        png_bytep* rows_nv = rows;
        png_read_image(png_ptr_nv, rows_nv);
        png_read_end(png_ptr_nv, info_ptr_nv);
    }

    width = static_cast<int>(png_width);
    height = static_cast<int>(png_height);
    heightmap.resize(static_cast<size_t>(width) * static_cast<size_t>(height));

    for (int y = 0; y < height; ++y) {
        const png_bytep row = rows[y];
        for (int x = 0; x < width; ++x) {
            const png_bytep ptr = &(row[x * 3]);
            const float luminance =
                (0.2126f * ptr[0] + 0.7152f * ptr[1] + 0.0722f * ptr[2]) / 255.0f;
            heightmap[static_cast<size_t>(y) * static_cast<size_t>(width) + static_cast<size_t>(x)] =
                luminance;
        }
    }

    code = 0;

finalise:
    if (fp != nullptr) fclose(fp);
    if (rows != nullptr) free(const_cast<png_bytep*>(rows));
    if (image_data != nullptr) free(const_cast<png_bytep>(image_data));
    if (png_ptr != nullptr) {
        if (info_ptr != nullptr) {
            png_free_data(png_ptr, info_ptr, PNG_FREE_ALL, -1);
        }
        {
            png_structp png_ptr_nv = png_ptr;
            png_destroy_read_struct(&png_ptr_nv, (png_infopp)nullptr, (png_infopp)nullptr);
        }
    }
    if (code != 0) {
        heightmap.clear();
        width = 0;
        height = 0;
    }
    return code;
}
