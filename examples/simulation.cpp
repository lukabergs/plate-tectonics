#ifndef NOMINMAX
#define NOMINMAX
#endif

#include "map_drawing.hpp"
#include "platecapi.hpp"
#include "sqrdmd.hpp"
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <filesystem>
#include <limits>
#include <string>
#include <vector>

namespace fs = std::filesystem;

#ifndef PLATE_TECTONICS_ROOT
#define PLATE_TECTONICS_ROOT "."
#endif

namespace {

constexpr float kSeaLevel = 0.65f;
constexpr uint32_t kErosionPeriod = 60;
constexpr float kFoldingRatio = 0.02f;
constexpr uint32_t kAggrOverlapAbs = 1000000;
constexpr float kAggrOverlapRel = 0.33f;
constexpr uint32_t kCycleCount = 2;
constexpr uint32_t kPlateCount = 10;

const fs::path kRepoRoot(PLATE_TECTONICS_ROOT);
const fs::path kDefaultInputDir = kRepoRoot / "img" / "in";
const fs::path kDefaultOutputDir = kRepoRoot / "img" / "out";

struct Params {
    uint32_t seed;
    uint32_t width;
    uint32_t height;
    bool colors;
    uint32_t step;
    bool custom_dimensions;
    std::string filename;
    std::string input_path;
};

[[noreturn]] void fail(const std::string& message)
{
    fprintf(stderr, "error: %s\n", message.c_str());
    exit(1);
}

uint32_t parse_u32(const char* flag, const char* value)
{
    char* end = nullptr;
    const unsigned long parsed = strtoul(value, &end, 10);
    if (end == value || *end != '\0' || parsed == 0 ||
        parsed > static_cast<unsigned long>(std::numeric_limits<uint32_t>::max())) {
        fail(std::string("invalid value for ") + flag + ": " + value);
    }
    return static_cast<uint32_t>(parsed);
}

void print_help()
{
    printf(" -h --help           : show this message\n");
    printf(" -s SEED             : use the given SEED\n");
    printf(" -i --input FILE     : load FILE from the given path or from %s\n",
           kDefaultInputDir.string().c_str());
    printf(" --dim WIDTH HEIGHT  : use the given width and height when no input image is used\n");
    printf(" --colors            : generate a colors map\n");
    printf(" --grayscale         : generate a grayscale map\n");
    printf(" --filename NAME     : output basename when no input image is used\n");
    printf(" --step X            : save intermediate maps every X steps\n");
}

Params fill_params(int argc, char* argv[])
{
    srand(static_cast<unsigned int>(time(nullptr)));

    Params params = {
        static_cast<uint32_t>(rand()),
        600,
        400,
        true,
        0,
        false,
        "simulation",
        ""
    };

    int p = 1;
    while (p < argc) {
        if (strcmp(argv[p], "--help") == 0 || strcmp(argv[p], "-h") == 0) {
            print_help();
            exit(0);
        } else if (strcmp(argv[p], "-s") == 0) {
            if (p + 1 >= argc) {
                fail("a parameter should follow -s");
            }
            params.seed = parse_u32("-s", argv[p + 1]);
            p += 2;
        } else if (strcmp(argv[p], "--dim") == 0) {
            if (p + 2 >= argc) {
                fail("two parameters should follow --dim");
            }
            params.width = parse_u32("--dim", argv[p + 1]);
            params.height = parse_u32("--dim", argv[p + 2]);
            if (params.width < 5 || params.height < 5) {
                fail("dimensions have to be >= 5");
            }
            params.custom_dimensions = true;
            p += 3;
        } else if (strcmp(argv[p], "--colors") == 0) {
            params.colors = true;
            p += 1;
        } else if (strcmp(argv[p], "--grayscale") == 0) {
            params.colors = false;
            p += 1;
        } else if (strcmp(argv[p], "--filename") == 0) {
            if (p + 1 >= argc) {
                fail("a parameter should follow --filename");
            }
            params.filename = argv[p + 1];
            p += 2;
        } else if (strcmp(argv[p], "--step") == 0) {
            if (p + 1 >= argc) {
                fail("a parameter should follow --step");
            }
            params.step = parse_u32("--step", argv[p + 1]);
            p += 2;
        } else if (strcmp(argv[p], "--input") == 0 || strcmp(argv[p], "-i") == 0) {
            if (p + 1 >= argc) {
                fail("a parameter should follow --input");
            }
            params.input_path = argv[p + 1];
            p += 2;
        } else {
            fail(std::string("unexpected param '") + argv[p] + "', use -h to display a list of params");
        }
    }

    if (!params.input_path.empty() && params.custom_dimensions) {
        fail("--dim cannot be used together with --input");
    }

    return params;
}

fs::path resolve_input_path(const std::string& input_path)
{
    const fs::path provided(input_path);
    std::vector<fs::path> candidates;

    if (provided.is_absolute()) {
        candidates.push_back(provided);
    } else {
        candidates.push_back(fs::current_path() / provided);
        candidates.push_back(kDefaultInputDir / provided);
    }

    if (provided.extension().empty()) {
        const fs::path with_png = provided.string() + ".png";
        if (with_png.is_absolute()) {
            candidates.push_back(with_png);
        } else {
            candidates.push_back(fs::current_path() / with_png);
            candidates.push_back(kDefaultInputDir / with_png);
        }
    }

    for (const fs::path& candidate : candidates) {
        if (fs::exists(candidate) && fs::is_regular_file(candidate)) {
            return fs::weakly_canonical(candidate);
        }
    }

    fail("input image not found: " + input_path);
}

uint32_t next_run_index(const fs::path& output_dir, const std::string& stem)
{
    uint32_t next_index = 0;
    const std::string prefix = stem + "_";

    if (!fs::exists(output_dir)) {
        return next_index;
    }

    for (const fs::directory_entry& entry : fs::directory_iterator(output_dir)) {
        if (!entry.is_regular_file() || entry.path().extension() != ".png") {
            continue;
        }

        const std::string name = entry.path().stem().string();
        if (name.rfind(prefix, 0) != 0) {
            continue;
        }

        const std::string suffix = name.substr(prefix.size());
        if (suffix.empty()) {
            continue;
        }

        char* end = nullptr;
        const unsigned long parsed = strtoul(suffix.c_str(), &end, 10);
        if (*end == '\0' &&
            parsed <= static_cast<unsigned long>(std::numeric_limits<uint32_t>::max()) &&
            parsed >= next_index) {
            next_index = static_cast<uint32_t>(parsed) + 1;
        }
    }

    return next_index;
}

void save_image(const float* heightmap, const fs::path& filename, int width, int height, bool colors)
{
    std::vector<float> copy(heightmap, heightmap + static_cast<size_t>(width) * static_cast<size_t>(height));
    normalize(copy.data(), width * height);

    const int result = colors
        ? writeImageColors(filename.string().c_str(), width, height, copy.data(), "Plate Tectonics")
        : writeImageGray(filename.string().c_str(), width, height, copy.data(), "Plate Tectonics");

    if (result != 0) {
        fail("failed to write image: " + filename.string());
    }
}

void save_image(void* simulation, const fs::path& filename, int width, int height, bool colors)
{
    save_image(platec_api_get_heightmap(simulation), filename, width, height, colors);
}

} // namespace

int main(int argc, char* argv[])
{
    Params params = fill_params(argc, argv);
    std::vector<float> input_heightmap;
    fs::path resolved_input_path;
    fs::path final_output_path;
    std::string output_label;
    uint32_t run_index = 0;

    if (!params.input_path.empty()) {
        resolved_input_path = resolve_input_path(params.input_path);

        int input_width = 0;
        int input_height = 0;
        if (readImageNormalized(resolved_input_path.string().c_str(), input_heightmap, input_width, input_height) != 0) {
            fail("failed to read input image: " + resolved_input_path.string());
        }

        params.width = static_cast<uint32_t>(input_width);
        params.height = static_cast<uint32_t>(input_height);

        fs::create_directories(kDefaultOutputDir);
        run_index = next_run_index(kDefaultOutputDir, resolved_input_path.stem().string());
        final_output_path =
            kDefaultOutputDir / (resolved_input_path.stem().string() + "_" + std::to_string(run_index) + ".png");
        output_label = final_output_path.string();
    } else {
        final_output_path = params.filename + ".png";
        output_label = final_output_path.string();
    }

    printf("Plate-tectonics simulation example\n");
    printf(" seed     : %u\n", params.seed);
    printf(" width    : %u\n", params.width);
    printf(" height   : %u\n", params.height);
    printf(" map      : %s\n", params.colors ? "colors" : "grayscale");
    if (!resolved_input_path.empty()) {
        printf(" input    : %s\n", resolved_input_path.string().c_str());
        printf(" output   : %s\n", output_label.c_str());
    } else {
        printf(" output   : %s\n", output_label.c_str());
    }
    if (params.step == 0) {
        printf(" step     : no\n");
    } else {
        printf(" step     : %u\n", params.step);
    }
    printf("\n");

    void* simulation = platec_api_create(params.seed, params.width, params.height, kSeaLevel, kErosionPeriod,
                                         kFoldingRatio, kAggrOverlapAbs, kAggrOverlapRel,
                                         kCycleCount, kPlateCount);

    if (!input_heightmap.empty()) {
        platec_api_load_heightmap(simulation, input_heightmap.data(), kSeaLevel);
    }

    uint32_t step = 0;
    while (platec_api_is_finished(simulation) == 0) {
        ++step;
        platec_api_step(simulation);

        if (params.step != 0 && step % params.step == 0) {
            fs::path step_output_path;
            if (!resolved_input_path.empty()) {
                step_output_path = kDefaultOutputDir /
                    (resolved_input_path.stem().string() + "_" + std::to_string(run_index) +
                     "_step_" + std::to_string(step) + ".png");
            } else {
                step_output_path = params.filename + "_" + std::to_string(step) + ".png";
            }
            save_image(simulation, step_output_path, static_cast<int>(params.width),
                       static_cast<int>(params.height), params.colors);
            printf(" * step %u (filename %s)\n", step, step_output_path.string().c_str());
        }
    }

    save_image(simulation, final_output_path, static_cast<int>(params.width),
               static_cast<int>(params.height), params.colors);
    printf(" * simulation completed (filename %s)\n", final_output_path.string().c_str());

    platec_api_destroy(simulation);
    return 0;
}
