#ifndef MAP_DRAWING
#define MAP_DRAWING

#include <stdio.h>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <png.h>

// This function actually writes out the PNG image file. The string 'title' is
// also written into the image file
int writeImageGray(const char* filename, int width, int height, float *heightmap, const char* title);

int writeImageColors(const char* filename, int width, int height, float *heightmap, const char* title);

int writeImageRgb(const char* filename, int width, int height, const png_byte* rgb, const char* title);

void renderImageGrayRgb(int width, int height, float* heightmap, std::vector<png_byte>& rgb);

void renderImageColorsRgb(int width, int height, float* heightmap, std::vector<png_byte>& rgb);

int readImageNormalized(const char* filename, std::vector<float>& heightmap, int& width, int& height);

#endif
