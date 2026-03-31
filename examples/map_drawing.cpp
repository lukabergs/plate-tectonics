#include "map_drawing.hpp"
#include "utils.hpp"
#include <stdexcept>
#include <vector>

using namespace std;

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

void drawColorsImage(png_structp& png_ptr, png_bytep& row, int width, int height, float *heightmap)
{
    float q15 = find_value_for_quantile(0.15f, heightmap, width * height);
    float q70 = find_value_for_quantile(0.70f, heightmap, width * height);
    float q75 = find_value_for_quantile(0.75f, heightmap, width * height);
    float q90 = find_value_for_quantile(0.90f, heightmap, width * height);
    float q95 = find_value_for_quantile(0.95f, heightmap, width * height);
    float q99 = find_value_for_quantile(0.99f, heightmap, width * height);

    int x, y;
    for (y=0 ; y<height ; y++) {
        for (x=0 ; x<width ; x++) {
            float h = heightmap[(y*width + x)];

            if (h < q15) {
                gradient(&(row[x*3]), 0, 0, 255, 0, 20, 200, h, 0.0f, q15);
                continue;
            }

            if (h < q70) {
                gradient(&(row[x*3]), 0, 20, 200, 50, 80, 225, h, q15, q70);
                continue;
            }

            if (h < q75) {
                gradient(&(row[x*3]), 50, 80, 225, 135, 237, 235, h, q70, q75);
                continue;
            }

            if (h < q90) {
                gradient(&(row[x*3]), 88, 173, 49, 218, 226, 58, h, q75, q90);
                continue;
            }

            if (h < q95) {
                gradient(&(row[x*3]), 218, 226, 58, 251, 252, 42, h, q90, q95);
                continue;
            }

            if (h < q99) {
                gradient(&(row[x*3]), 251, 252, 42, 91, 28, 13, h, q95, q99);
                continue;
            }

            gradient(&(row[x*3]), 91, 28, 13, 51, 0, 4, h, q99, 1.0f);
        }
        png_write_row(png_ptr, row);
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
