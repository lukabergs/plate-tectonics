#ifndef MAP_DRAWING
#define MAP_DRAWING

#include <cstdint>
#include <png.h>
#include <vector>

enum class FaultKind : std::uint8_t {
    None = 0,
    Convergent = 1,
    Divergent = 2,
    Transform = 3,
};

int writeImageGray(const char* filename, int width, int height, float* heightmap,
                   const char* title);

int writeImageColors(const char* filename, int width, int height, float* heightmap,
                     const char* title);

int writeImageGrayWithFaultLines(const char* filename, int width, int height, float* heightmap,
                                 const uint32_t* platesmap, const float* plate_vx,
                                 const float* plate_vy, uint32_t plate_count,
                                 const char* title);

int writeImageColorsWithFaultLines(const char* filename, int width, int height, float* heightmap,
                                   const uint32_t* platesmap, const float* plate_vx,
                                   const float* plate_vy, uint32_t plate_count,
                                   const char* title);

int writeImageGrayWithFaultMap(const char* filename, int width, int height, float* heightmap,
                               const FaultKind* fault_map, const char* title);

int writeImageColorsWithFaultMap(const char* filename, int width, int height, float* heightmap,
                                 const FaultKind* fault_map, const char* title);

int buildFaultMap(int width, int height, const uint32_t* platesmap, const float* plate_vx,
                  const float* plate_vy, uint32_t plate_count,
                  std::vector<FaultKind>* fault_map);

int loadHeightmapFromPng(const char* filename, uint32_t* width, uint32_t* height,
                         std::vector<float>* heightmap);

int loadHeightmapFromImage(const char* filename, uint32_t* width, uint32_t* height,
                           std::vector<float>* heightmap);

#endif
