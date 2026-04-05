#ifndef NOMINMAX
#define NOMINMAX
#endif

#include "heightmap_io.hpp"
#include "platecapi.hpp"
#include "sqrdmd.hpp"

#include <cmath>
#include <cctype>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <limits>
#include <string>
#include <system_error>
#include <vector>

namespace fs = std::filesystem;

#ifndef PLATE_TECTONICS_ROOT
#define PLATE_TECTONICS_ROOT "."
#endif

namespace {

constexpr float kSeaLevel = 0.65f;
constexpr uint32_t kErosionPeriod = 60;
constexpr float kDefaultErosionStrength = 1.0f;
constexpr float kFoldingRatio = 0.02f;
constexpr uint32_t kAggrOverlapAbs = 1000000;
constexpr float kAggrOverlapRel = 0.33f;
constexpr uint32_t kDefaultCycleCount = 2;
constexpr uint32_t kDefaultPlateCount = 10;
constexpr float kDefaultRotationStrength = 1.0f;
constexpr float kDefaultLandmassRotation = 0.0f;
constexpr float kDefaultSubductionStrength = 1.0f;
constexpr uint16_t kDefaultMinInitialHeight = TopographyCodec::kDefaultInitialMinHeightMeters;
constexpr uint16_t kDefaultMaxInitialHeight = TopographyCodec::kDefaultInitialMaxHeightMeters;
constexpr uint16_t kGifDelayCs = 8;

const fs::path kRepoRoot(PLATE_TECTONICS_ROOT);
const fs::path kDefaultInputDir = kRepoRoot / "img" / "in";
const fs::path kDefaultOutputDir = kRepoRoot / "img" / "out";
const fs::path kDefaultGifDir = kRepoRoot / "img" / "gif";
const fs::path kDefaultFramesDir = kRepoRoot / "img" / "frames";

enum class ToneMapMode {
    Linear,
    Log,
    Asinh
};

enum class InputMode {
    Procedural,
    LegacyNormalizedPng,
    MetricPng16,
    MetricR16
};

struct Params {
    uint32_t seed;
    uint32_t width;
    uint32_t height;
    bool colors;
    bool normalize_output;
    bool export_heightmap_f32;
    bool make_gif;
    bool delete_gif_frames;
    bool show_boundaries;
    uint32_t cycles;
    uint32_t plates;
    uint32_t aggregation_overlap_abs;
    float aggregation_overlap_rel;
    uint32_t erosion_period;
    float folding_ratio;
    float erosion_strength;
    float rotation_strength;
    float landmass_rotation;
    float subduction_strength;
    uint16_t min_initial_height_m;
    uint16_t max_initial_height_m;
    float display_min;
    float display_max;
    ToneMapMode tone_map;
    uint32_t step;
    bool custom_dimensions;
    bool has_sea_level_m;
    uint16_t sea_level_m;
    bool legacy_filename_override;
    InputMode input_mode;
    std::string filename;
    std::string input_path;
};

struct BoundaryOverlayData {
    std::vector<uint32_t> platesmap;
    std::vector<float> linear_vel_x;
    std::vector<float> linear_vel_y;
    std::vector<float> angular_velocity;
    std::vector<float> mass_center_x;
    std::vector<float> mass_center_y;
    uint32_t plate_count = 0;

    bool valid() const
    {
        return plate_count > 0 && !platesmap.empty() &&
               linear_vel_x.size() == plate_count &&
               linear_vel_y.size() == plate_count &&
               angular_velocity.size() == plate_count &&
               mass_center_x.size() == plate_count &&
               mass_center_y.size() == plate_count;
    }
};

struct BoundaryVelocity {
    float x;
    float y;
};

struct NormalizationRange {
    bool enabled = false;
    float min = 0.0f;
    float max = 1.0f;
    float output_min = 0.0f;
    float output_max = 1.0f;
};

[[noreturn]] void fail(const std::string& message)
{
    fprintf(stderr, "error: %s\n", message.c_str());
    exit(1);
}

std::string display_path(fs::path path)
{
    path.make_preferred();
    return path.string();
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

uint32_t parse_nonnegative_u32(const char* flag, const char* value)
{
    char* end = nullptr;
    const unsigned long parsed = strtoul(value, &end, 10);
    if (end == value || *end != '\0' ||
        parsed > static_cast<unsigned long>(std::numeric_limits<uint32_t>::max())) {
        fail(std::string("invalid value for ") + flag + ": " + value);
    }
    return static_cast<uint32_t>(parsed);
}

uint16_t parse_u16(const char* flag, const char* value)
{
    const uint32_t parsed = parse_nonnegative_u32(flag, value);
    if (parsed > TopographyCodec::kMaxHeightMeters) {
        fail(std::string("invalid value for ") + flag + ": " + value);
    }
    return static_cast<uint16_t>(parsed);
}

float parse_nonnegative_float(const char* flag, const char* value)
{
    char* end = nullptr;
    const float parsed = strtof(value, &end);
    if (end == value || *end != '\0' || !std::isfinite(parsed) || parsed < 0.0f) {
        fail(std::string("invalid value for ") + flag + ": " + value);
    }
    return parsed;
}

float parse_float(const char* flag, const char* value)
{
    char* end = nullptr;
    const float parsed = strtof(value, &end);
    if (end == value || *end != '\0' || !std::isfinite(parsed)) {
        fail(std::string("invalid value for ") + flag + ": " + value);
    }
    return parsed;
}

float parse_unit_float(const char* flag, const char* value)
{
    const float parsed = parse_nonnegative_float(flag, value);
    if (parsed > 1.0f) {
        fail(std::string("invalid value for ") + flag + ": " + value);
    }
    return parsed;
}

ToneMapMode parse_tone_map(const char* value)
{
    if (strcmp(value, "linear") == 0) {
        return ToneMapMode::Linear;
    }
    if (strcmp(value, "log") == 0) {
        return ToneMapMode::Log;
    }
    if (strcmp(value, "asinh") == 0) {
        return ToneMapMode::Asinh;
    }

    fail(std::string("invalid value for --tone-map: ") + value);
}

const char* tone_map_name(ToneMapMode mode)
{
    switch (mode) {
    case ToneMapMode::Linear:
        return "linear";
    case ToneMapMode::Log:
        return "log";
    case ToneMapMode::Asinh:
        return "asinh";
    default:
        return "linear";
    }
}

std::string escape_powershell_single_quotes(const std::string& value)
{
    std::string escaped;
    escaped.reserve(value.size());
    for (char ch : value) {
        escaped += ch;
        if (ch == '\'') {
            escaped += '\'';
        }
    }
    return escaped;
}

void print_help()
{
    printf(" -h --help           : show this message\n");
    printf(" -s SEED             : use the given SEED\n");
    printf(" -i --input FILE     : load legacy normalized PNG FILE from the given path or from %s\n",
           display_path(kDefaultInputDir).c_str());
    printf(" --input-png16 FILE  : load 16-bit grayscale metric PNG\n");
    printf(" --input-r16 FILE    : load little-endian metric .r16 raw data\n");
    printf(" --dim WIDTH HEIGHT  : use the given width and height when no input image is used\n");
    printf(" --sea-level-m N     : metric coastline threshold in [0, 65535]\n");
    printf(" --colors            : colorize preview frames/GIF only\n");
    printf(" --grayscale         : grayscale preview frames/GIF only\n");
    printf(" --no-output-normalization: disable preview-frame normalization\n");
    printf(" --min-initial-height N: minimum seeded metric height in meters (default %u)\n",
           static_cast<unsigned>(kDefaultMinInitialHeight));
    printf(" --max-initial-height N: maximum seeded metric height in meters (default %u)\n",
           static_cast<unsigned>(kDefaultMaxInitialHeight));
    printf(" --display-min X     : lower bound of the preview display window\n");
    printf(" --display-max X     : upper bound of the preview display window\n");
    printf(" --tone-map MODE     : preview tone map for frames/GIF: linear, log, asinh\n");
    printf(" --export-heightmap-f32: write the final raw float32 debug heightmap plus metadata\n");
    printf(" --filename NAME     : deprecated; saved filenames now use timestamp + seed/input name\n");
    printf(" --cycles N          : number of simulation cycles to run\n");
    printf(" --plates N          : number of tectonic plates to initialize\n");
    printf(" --aggregation-overlap-abs N: overlap pixels needed to aggregate continents\n");
    printf(" --aggregation-overlap-rel X: overlap ratio needed to aggregate continents\n");
    printf(" --folding-ratio X   : fraction of overlapping continental crust turned into uplift\n");
    printf(" --erosion-period N  : number of updates between erosion passes\n");
    printf(" --erosion-strength X: scale erosion amount (0 disables)\n");
    printf(" --rotation-strength X: scale angular plate motion\n");
    printf(" --landmass-rotation X: scale visible crust rotation (0 disables)\n");
    printf(" --subduction-strength X: scale oceanic crust removal during subduction\n");
    printf(" --step X            : save intermediate maps every X steps\n");
    printf(" --gif               : export an animated GIF using --step sampling, or every update if --step is omitted\n");
    printf(" --no-steps          : delete GIF frame PNGs after the GIF is created\n");
    printf(" --show-boundaries   : overlay convergent/divergent/transform boundaries\n");
    printf(" output files        : img/in/<YYMMDDHHMMSS>_<seed-or-input>.r16/.r16.json and PNG16/.png.json\n");
    printf("                       img/out/<YYMMDDHHMMSS>_<seed-or-input>.r16/.r16.json and PNG16/.png.json\n");
    printf("                       img/gif/<YYMMDDHHMMSS>_<seed-or-input>.gif when --gif is used\n");
    printf("                       img/frames/<YYMMDDHHMMSS>_<seed-or-input>_<frame>.png for preview frames\n");
}

Params fill_params(int argc, char* argv[])
{
    srand(static_cast<unsigned int>(time(nullptr)));

    Params params = {
        static_cast<uint32_t>(rand()),
        600,
        400,
        true,
        true,
        false,
        false,
        false,
        false,
        kDefaultCycleCount,
        kDefaultPlateCount,
        kAggrOverlapAbs,
        kAggrOverlapRel,
        kErosionPeriod,
        kFoldingRatio,
        kDefaultErosionStrength,
        kDefaultRotationStrength,
        kDefaultLandmassRotation,
        kDefaultSubductionStrength,
        kDefaultMinInitialHeight,
        kDefaultMaxInitialHeight,
        0.0f,
        1.0f,
        ToneMapMode::Linear,
        0,
        false,
        false,
        0,
        false,
        InputMode::Procedural,
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
        } else if (strcmp(argv[p], "--sea-level-m") == 0) {
            if (p + 1 >= argc) {
                fail("a parameter should follow --sea-level-m");
            }
            params.has_sea_level_m = true;
            params.sea_level_m = parse_u16("--sea-level-m", argv[p + 1]);
            p += 2;
        } else if (strcmp(argv[p], "--colors") == 0) {
            params.colors = true;
            p += 1;
        } else if (strcmp(argv[p], "--grayscale") == 0) {
            params.colors = false;
            p += 1;
        } else if (strcmp(argv[p], "--no-output-normalization") == 0) {
            params.normalize_output = false;
            p += 1;
        } else if (strcmp(argv[p], "--min-initial-height") == 0) {
            if (p + 1 >= argc) {
                fail("a parameter should follow --min-initial-height");
            }
            params.min_initial_height_m = parse_u16("--min-initial-height", argv[p + 1]);
            p += 2;
        } else if (strcmp(argv[p], "--max-initial-height") == 0) {
            if (p + 1 >= argc) {
                fail("a parameter should follow --max-initial-height");
            }
            params.max_initial_height_m = parse_u16("--max-initial-height", argv[p + 1]);
            p += 2;
        } else if (strcmp(argv[p], "--display-min") == 0) {
            if (p + 1 >= argc) {
                fail("a parameter should follow --display-min");
            }
            params.display_min = parse_float("--display-min", argv[p + 1]);
            p += 2;
        } else if (strcmp(argv[p], "--display-max") == 0) {
            if (p + 1 >= argc) {
                fail("a parameter should follow --display-max");
            }
            params.display_max = parse_float("--display-max", argv[p + 1]);
            p += 2;
        } else if (strcmp(argv[p], "--tone-map") == 0) {
            if (p + 1 >= argc) {
                fail("a parameter should follow --tone-map");
            }
            params.tone_map = parse_tone_map(argv[p + 1]);
            p += 2;
        } else if (strcmp(argv[p], "--export-heightmap-f32") == 0) {
            params.export_heightmap_f32 = true;
            p += 1;
        } else if (strcmp(argv[p], "--filename") == 0) {
            if (p + 1 >= argc) {
                fail("a parameter should follow --filename");
            }
            params.legacy_filename_override = true;
            params.filename = argv[p + 1];
            p += 2;
        } else if (strcmp(argv[p], "--cycles") == 0) {
            if (p + 1 >= argc) {
                fail("a parameter should follow --cycles");
            }
            params.cycles = parse_u32("--cycles", argv[p + 1]);
            p += 2;
        } else if (strcmp(argv[p], "--plates") == 0) {
            if (p + 1 >= argc) {
                fail("a parameter should follow --plates");
            }
            params.plates = parse_u32("--plates", argv[p + 1]);
            p += 2;
        } else if (strcmp(argv[p], "--aggregation-overlap-abs") == 0) {
            if (p + 1 >= argc) {
                fail("a parameter should follow --aggregation-overlap-abs");
            }
            params.aggregation_overlap_abs =
                parse_nonnegative_u32("--aggregation-overlap-abs", argv[p + 1]);
            p += 2;
        } else if (strcmp(argv[p], "--aggregation-overlap-rel") == 0) {
            if (p + 1 >= argc) {
                fail("a parameter should follow --aggregation-overlap-rel");
            }
            params.aggregation_overlap_rel =
                parse_unit_float("--aggregation-overlap-rel", argv[p + 1]);
            p += 2;
        } else if (strcmp(argv[p], "--folding-ratio") == 0) {
            if (p + 1 >= argc) {
                fail("a parameter should follow --folding-ratio");
            }
            params.folding_ratio = parse_nonnegative_float("--folding-ratio", argv[p + 1]);
            p += 2;
        } else if (strcmp(argv[p], "--erosion-period") == 0) {
            if (p + 1 >= argc) {
                fail("a parameter should follow --erosion-period");
            }
            params.erosion_period = parse_u32("--erosion-period", argv[p + 1]);
            p += 2;
        } else if (strcmp(argv[p], "--erosion-strength") == 0) {
            if (p + 1 >= argc) {
                fail("a parameter should follow --erosion-strength");
            }
            params.erosion_strength = parse_nonnegative_float("--erosion-strength", argv[p + 1]);
            p += 2;
        } else if (strcmp(argv[p], "--rotation-strength") == 0) {
            if (p + 1 >= argc) {
                fail("a parameter should follow --rotation-strength");
            }
            params.rotation_strength = parse_nonnegative_float("--rotation-strength", argv[p + 1]);
            p += 2;
        } else if (strcmp(argv[p], "--landmass-rotation") == 0) {
            if (p + 1 >= argc) {
                fail("a parameter should follow --landmass-rotation");
            }
            params.landmass_rotation =
                parse_nonnegative_float("--landmass-rotation", argv[p + 1]);
            p += 2;
        } else if (strcmp(argv[p], "--subduction-strength") == 0) {
            if (p + 1 >= argc) {
                fail("a parameter should follow --subduction-strength");
            }
            params.subduction_strength =
                parse_nonnegative_float("--subduction-strength", argv[p + 1]);
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
            if (params.input_mode != InputMode::Procedural) {
                fail("only one input mode can be used");
            }
            params.input_mode = InputMode::LegacyNormalizedPng;
            params.input_path = argv[p + 1];
            p += 2;
        } else if (strcmp(argv[p], "--input-png16") == 0) {
            if (p + 1 >= argc) {
                fail("a parameter should follow --input-png16");
            }
            if (params.input_mode != InputMode::Procedural) {
                fail("only one input mode can be used");
            }
            params.input_mode = InputMode::MetricPng16;
            params.input_path = argv[p + 1];
            p += 2;
        } else if (strcmp(argv[p], "--input-r16") == 0) {
            if (p + 1 >= argc) {
                fail("a parameter should follow --input-r16");
            }
            if (params.input_mode != InputMode::Procedural) {
                fail("only one input mode can be used");
            }
            params.input_mode = InputMode::MetricR16;
            params.input_path = argv[p + 1];
            p += 2;
        } else if (strcmp(argv[p], "--gif") == 0) {
            params.make_gif = true;
            p += 1;
        } else if (strcmp(argv[p], "--no-steps") == 0) {
            params.delete_gif_frames = true;
            p += 1;
        } else if (strcmp(argv[p], "--show-boundaries") == 0) {
            params.show_boundaries = true;
            p += 1;
        } else {
            fail(std::string("unexpected param '") + argv[p] + "', use -h to display a list of params");
        }
    }

    if (params.input_mode == InputMode::LegacyNormalizedPng && params.custom_dimensions) {
        fail("--dim cannot be used together with --input");
    }
    if (params.input_mode == InputMode::MetricPng16 && params.custom_dimensions) {
        fail("--dim cannot be used together with --input-png16");
    }
    if (params.input_path.empty() && params.input_mode != InputMode::Procedural) {
        fail("missing input path");
    }
    if (!params.input_path.empty() && params.input_mode == InputMode::Procedural) {
        fail("internal input parsing error");
    }
    if (params.input_mode == InputMode::LegacyNormalizedPng && params.has_sea_level_m) {
        fail("--sea-level-m is not supported with legacy --input");
    }
    if (params.delete_gif_frames && !params.make_gif) {
        fail("--no-steps requires --gif");
    }
    if (params.min_initial_height_m >= params.max_initial_height_m) {
        fail("--min-initial-height must be lower than --max-initial-height");
    }
    if (static_cast<uint32_t>(params.max_initial_height_m) -
            static_cast<uint32_t>(params.min_initial_height_m) <
        2U) {
        fail("--max-initial-height must be at least 2 meters above --min-initial-height");
    }
    if (params.input_mode == InputMode::Procedural && params.has_sea_level_m &&
        (params.sea_level_m <= params.min_initial_height_m ||
         params.sea_level_m >= params.max_initial_height_m)) {
        fail("--sea-level-m must lie strictly between --min-initial-height and --max-initial-height");
    }
    if (!(params.display_min < params.display_max)) {
        fail("--display-min must be smaller than --display-max");
    }

    return params;
}

fs::path resolve_input_path(const std::string& input_path, const char* default_extension = nullptr)
{
    const fs::path provided(input_path);
    std::vector<fs::path> candidates;

    if (provided.is_absolute()) {
        candidates.push_back(provided);
    } else {
        candidates.push_back(fs::current_path() / provided);
        candidates.push_back(kDefaultInputDir / provided);
    }

    if (provided.extension().empty() && default_extension != nullptr) {
        const fs::path with_extension = provided.string() + default_extension;
        if (with_extension.is_absolute()) {
            candidates.push_back(with_extension);
        } else {
            candidates.push_back(fs::current_path() / with_extension);
            candidates.push_back(kDefaultInputDir / with_extension);
        }
    }

    for (const fs::path& candidate : candidates) {
        if (fs::exists(candidate) && fs::is_regular_file(candidate)) {
            return fs::weakly_canonical(candidate);
        }
    }

    fail("input image not found: " + input_path);
}

void ensure_output_directories()
{
    fs::create_directories(kDefaultInputDir);
    fs::create_directories(kDefaultOutputDir);
    fs::create_directories(kDefaultGifDir);
    fs::create_directories(kDefaultFramesDir);
}

std::string sanitize_file_component(const std::string& value)
{
    std::string sanitized;
    sanitized.reserve(value.size());

    for (char ch : value) {
        const unsigned char byte = static_cast<unsigned char>(ch);
        if (std::isalnum(byte) != 0 || ch == '-' || ch == '_') {
            sanitized += ch;
        } else if (sanitized.empty() || sanitized.back() != '_') {
            sanitized += '_';
        }
    }

    while (!sanitized.empty() && sanitized.back() == '_') {
        sanitized.pop_back();
    }

    return sanitized.empty() ? "run" : sanitized;
}

std::string make_run_timestamp(std::time_t now)
{
    std::tm local_tm = {};
#ifdef _WIN32
    localtime_s(&local_tm, &now);
#else
    localtime_r(&now, &local_tm);
#endif

    char buffer[32];
    if (std::strftime(buffer, sizeof(buffer), "%y%m%d%H%M%S", &local_tm) == 0) {
        fail("failed to format run timestamp");
    }
    return buffer;
}

std::string make_run_id(const Params& params, const fs::path& resolved_input_path,
                        std::time_t started_at)
{
    const std::string label = resolved_input_path.empty()
        ? std::to_string(params.seed)
        : resolved_input_path.stem().string();
    return make_run_timestamp(started_at) + "_" + sanitize_file_component(label);
}

std::string frame_file_name(const std::string& run_id, uint32_t frame_index)
{
    char buffer[128];
    snprintf(buffer, sizeof(buffer), "%s_%05u.png", run_id.c_str(), frame_index);
    return buffer;
}

NormalizationRange capture_normalization_range(const float* heightmap, int size, bool normalize_output,
                                               float output_min, float output_max)
{
    NormalizationRange range;
    range.enabled = normalize_output;
    range.output_min = output_min;
    range.output_max = output_max;
    if (!normalize_output || size <= 0) {
        return range;
    }

    float min_value = heightmap[0];
    float max_value = heightmap[0];
    for (int i = 1; i < size; ++i) {
        min_value = min_value < heightmap[i] ? min_value : heightmap[i];
        max_value = max_value > heightmap[i] ? max_value : heightmap[i];
    }

    range.min = min_value;
    range.max = max_value;
    return range;
}

void apply_normalization_range(std::vector<float>& values, const NormalizationRange& range)
{
    if (!range.enabled || values.empty()) {
        return;
    }

    const float diff = range.max - range.min;
    if (diff <= 0.0f) {
        const float fill = 0.5f * (range.output_min + range.output_max);
        for (float& value : values) {
            value = fill;
        }
        return;
    }

    const float output_diff = range.output_max - range.output_min;
    for (float& value : values) {
        value = range.output_min + ((value - range.min) / diff) * output_diff;
    }
}

float apply_tone_curve(float value, ToneMapMode tone_map)
{
    constexpr float kLogStrength = 9.0f;
    constexpr float kAsinhStrength = 8.0f;

    if (value <= 0.0f) {
        return 0.0f;
    }
    if (value >= 1.0f) {
        return 1.0f;
    }

    switch (tone_map) {
    case ToneMapMode::Linear:
        return value;
    case ToneMapMode::Log:
        return static_cast<float>(std::log1p(kLogStrength * value) / std::log1p(kLogStrength));
    case ToneMapMode::Asinh:
        return static_cast<float>(std::asinh(kAsinhStrength * value) / std::asinh(kAsinhStrength));
    default:
        return value;
    }
}

void apply_display_mapping(std::vector<float>& values, float display_min, float display_max,
                           ToneMapMode tone_map)
{
    if (values.empty()) {
        return;
    }

    const float display_range = display_max - display_min;
    for (float& value : values) {
        float mapped = (value - display_min) / display_range;
        if (mapped <= 0.0f) {
            value = 0.0f;
            continue;
        }
        if (mapped >= 1.0f) {
            value = 1.0f;
            continue;
        }
        value = apply_tone_curve(mapped, tone_map);
    }
}

void export_heightmap_f32(const float* heightmap, const fs::path& filename, int width, int height)
{
    const size_t sample_count = static_cast<size_t>(width) * static_cast<size_t>(height);
    std::ofstream output(filename, std::ios::binary);
    if (!output) {
        fail("failed to write heightmap export: " + display_path(filename));
    }
    output.write(reinterpret_cast<const char*>(heightmap), static_cast<std::streamsize>(sample_count * sizeof(float)));
    if (!output) {
        fail("failed to write heightmap export: " + display_path(filename));
    }

    const fs::path metadata_path = fs::path(filename.string() + ".json");
    std::ofstream metadata(metadata_path, std::ios::binary);
    if (!metadata) {
        fail("failed to write heightmap metadata: " + display_path(metadata_path));
    }
    metadata << "{\n";
    metadata << "  \"width\": " << width << ",\n";
    metadata << "  \"height\": " << height << ",\n";
    metadata << "  \"format\": \"float32\",\n";
    metadata << "  \"endianness\": \"little\",\n";
    metadata << "  \"layout\": \"row-major\",\n";
    metadata << "  \"origin\": \"top-left\"\n";
    metadata << "}\n";
    if (!metadata) {
        fail("failed to write heightmap metadata: " + display_path(metadata_path));
    }
}

fs::path metadata_sidecar_path(const fs::path& filename)
{
    return fs::path(filename.string() + ".json");
}

bool try_read_sidecar_metadata(const fs::path& filename, TopographyCodec::Metadata& metadata)
{
    const fs::path sidecar = metadata_sidecar_path(filename);
    return fs::exists(sidecar) &&
           readTopographyMetadataJson(sidecar.string().c_str(), metadata) == 0;
}

std::vector<uint16_t> encode_metric_heightmap(const float* heightmap, int width, int height,
                                              uint16_t sea_level_m)
{
    const size_t sample_count = static_cast<size_t>(width) * static_cast<size_t>(height);
    std::vector<uint16_t> metric_heightmap(sample_count);
    for (size_t i = 0; i < sample_count; ++i) {
        metric_heightmap[i] = TopographyCodec::internal_to_meters(heightmap[i], sea_level_m);
    }
    return metric_heightmap;
}

void export_metric_heightmap(const float* heightmap, int width, int height, uint16_t sea_level_m,
                             const fs::path& raw_path, const fs::path& png_path)
{
    const std::vector<uint16_t> metric_heightmap =
        encode_metric_heightmap(heightmap, width, height, sea_level_m);

    if (writeRawR16(raw_path.string().c_str(), metric_heightmap.data(), metric_heightmap.size()) != 0) {
        fail("failed to write metric heightmap: " + display_path(raw_path));
    }

    TopographyCodec::Metadata raw_metadata;
    raw_metadata.width = static_cast<uint32_t>(width);
    raw_metadata.height = static_cast<uint32_t>(height);
    raw_metadata.sea_level_m = sea_level_m;
    raw_metadata.format = TopographyCodec::kMetricFormatR16;
    raw_metadata.endianness = TopographyCodec::kLittleEndian;
    if (writeTopographyMetadataJson(metadata_sidecar_path(raw_path).string().c_str(), raw_metadata) != 0) {
        fail("failed to write metric metadata: " + display_path(metadata_sidecar_path(raw_path)));
    }

    if (writeImageGray16(png_path.string().c_str(), width, height, metric_heightmap.data(),
                         "Plate Tectonics Metric Heightmap") != 0) {
        fail("failed to write metric PNG16: " + display_path(png_path));
    }

    TopographyCodec::Metadata png_metadata = raw_metadata;
    png_metadata.format = TopographyCodec::kMetricFormatPng16;
    png_metadata.endianness = TopographyCodec::kBigEndian;
    if (writeTopographyMetadataJson(metadata_sidecar_path(png_path).string().c_str(), png_metadata) != 0) {
        fail("failed to write PNG16 metadata: " + display_path(metadata_sidecar_path(png_path)));
    }
}

uint16_t infer_metric_sea_level_with_warning(const std::vector<uint16_t>& heightmap_m,
                                             const fs::path& source_path)
{
    fprintf(stderr,
            "warning: no sea_level_m metadata for %s; falling back to %.2f ocean coverage quantile\n",
            display_path(source_path).c_str(), kSeaLevel);
    return TopographyCodec::infer_metric_sea_level(heightmap_m.data(), heightmap_m.size(), kSeaLevel);
}

void save_image(const float* heightmap, const fs::path& filename, int width, int height, bool colors,
                const NormalizationRange& normalization_range, float display_min,
                float display_max, ToneMapMode tone_map)
{
    std::vector<float> copy(heightmap,
                            heightmap + static_cast<size_t>(width) * static_cast<size_t>(height));
    apply_normalization_range(copy, normalization_range);
    apply_display_mapping(copy, display_min, display_max, tone_map);

    const int result = colors
        ? writeImageColors(filename.string().c_str(), width, height, copy.data(), "Plate Tectonics")
        : writeImageGray(filename.string().c_str(), width, height, copy.data(), "Plate Tectonics");

    if (result != 0) {
        fail("failed to write image: " + display_path(filename));
    }
}

enum class BoundaryType {
    None,
    Convergent,
    Divergent,
    Transform
};

BoundaryType classify_boundary(float vel_ax, float vel_ay, float vel_bx, float vel_by, float nx, float ny)
{
    const float rel_x = vel_bx - vel_ax;
    const float rel_y = vel_by - vel_ay;
    const float normal_component = rel_x * nx + rel_y * ny;
    const float tangent_component = std::fabs(rel_x * (-ny) + rel_y * nx);
    constexpr float kNormalThreshold = 0.15f;
    constexpr float kTangentialThreshold = 0.15f;

    if (normal_component < -kNormalThreshold) {
        return BoundaryType::Convergent;
    }
    if (normal_component > kNormalThreshold) {
        return BoundaryType::Divergent;
    }
    if (tangent_component > kTangentialThreshold) {
        return BoundaryType::Transform;
    }
    return BoundaryType::Transform;
}

BoundaryType strongest_boundary(BoundaryType current, BoundaryType candidate)
{
    auto score = [](BoundaryType type) {
        switch (type) {
        case BoundaryType::Convergent:
            return 3;
        case BoundaryType::Divergent:
            return 2;
        case BoundaryType::Transform:
            return 1;
        case BoundaryType::None:
        default:
            return 0;
        }
    };

    return score(candidate) > score(current) ? candidate : current;
}

void apply_boundary_color(std::vector<png_byte>& rgb, size_t pixel_index, BoundaryType boundary_type)
{
    png_byte* ptr = &rgb[pixel_index * 3U];
    switch (boundary_type) {
    case BoundaryType::Convergent:
        ptr[0] = 255;
        ptr[1] = 64;
        ptr[2] = 64;
        break;
    case BoundaryType::Divergent:
        ptr[0] = 64;
        ptr[1] = 192;
        ptr[2] = 255;
        break;
    case BoundaryType::Transform:
        ptr[0] = 255;
        ptr[1] = 215;
        ptr[2] = 0;
        break;
    case BoundaryType::None:
    default:
        break;
    }
}

bool capture_boundary_overlay(void* simulation, int width, int height, BoundaryOverlayData& data)
{
    const uint32_t plate_count = platec_api_get_plate_count(simulation);
    if (plate_count == 0) {
        return false;
    }

    const uint32_t* platesmap = platec_api_get_platesmap(simulation);
    if (platesmap == nullptr) {
        return false;
    }

    data.plate_count = plate_count;
    data.platesmap.assign(platesmap, platesmap + static_cast<size_t>(width) * static_cast<size_t>(height));
    data.linear_vel_x.resize(plate_count);
    data.linear_vel_y.resize(plate_count);
    data.angular_velocity.resize(plate_count);
    data.mass_center_x.resize(plate_count);
    data.mass_center_y.resize(plate_count);

    for (uint32_t plate_index = 0; plate_index < plate_count; ++plate_index) {
        data.linear_vel_x[plate_index] = platec_api_velocity_vector_x(simulation, plate_index);
        data.linear_vel_y[plate_index] = platec_api_velocity_vector_y(simulation, plate_index);
        data.angular_velocity[plate_index] = platec_api_angular_velocity(simulation, plate_index);
        data.mass_center_x[plate_index] = platec_api_mass_center_x(simulation, plate_index);
        data.mass_center_y[plate_index] = platec_api_mass_center_y(simulation, plate_index);
    }
    return true;
}

float wrapped_delta(float point, float center, float period)
{
    float delta = point - center;
    const float half_period = period * 0.5f;
    if (delta > half_period) {
        delta -= period;
    } else if (delta < -half_period) {
        delta += period;
    }
    return delta;
}

BoundaryVelocity surface_velocity_at(const BoundaryOverlayData& data, uint32_t plate_index, int x,
                                     int y, int width, int height)
{
    const float dx = wrapped_delta(static_cast<float>(x), data.mass_center_x[plate_index],
                                   static_cast<float>(width));
    const float dy = wrapped_delta(static_cast<float>(y), data.mass_center_y[plate_index],
                                   static_cast<float>(height));
    const float omega = data.angular_velocity[plate_index];

    return BoundaryVelocity {
        data.linear_vel_x[plate_index] - dy * omega,
        data.linear_vel_y[plate_index] + dx * omega
    };
}

void overlay_boundaries(const BoundaryOverlayData& data, int width, int height, std::vector<png_byte>& rgb)
{
    const int offsets[4][2] = {
        {1, 0},
        {-1, 0},
        {0, 1},
        {0, -1}
    };

    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            const size_t pixel_index = static_cast<size_t>(y) * static_cast<size_t>(width) +
                                       static_cast<size_t>(x);
            const uint32_t plate_a = data.platesmap[pixel_index];
            if (plate_a >= data.plate_count) {
                continue;
            }

            BoundaryType boundary_type = BoundaryType::None;
            for (const auto& offset : offsets) {
                const int nx = (x + offset[0] + width) % width;
                const int ny = (y + offset[1] + height) % height;
                const size_t neighbor_index = static_cast<size_t>(ny) * static_cast<size_t>(width) +
                                              static_cast<size_t>(nx);
                const uint32_t plate_b = data.platesmap[neighbor_index];
                if (plate_b >= data.plate_count || plate_a == plate_b) {
                    continue;
                }

                const BoundaryVelocity velocity_a =
                    surface_velocity_at(data, plate_a, x, y, width, height);
                const BoundaryVelocity velocity_b =
                    surface_velocity_at(data, plate_b, nx, ny, width, height);
                const BoundaryType candidate =
                    classify_boundary(velocity_a.x, velocity_a.y,
                                      velocity_b.x, velocity_b.y,
                                      static_cast<float>(offset[0]), static_cast<float>(offset[1]));
                boundary_type = strongest_boundary(boundary_type, candidate);
            }

            if (boundary_type != BoundaryType::None) {
                apply_boundary_color(rgb, pixel_index, boundary_type);
            }
        }
    }
}

void save_image(void* simulation, const fs::path& filename, int width, int height, bool colors,
                const NormalizationRange& normalization_range, float display_min,
                float display_max, ToneMapMode tone_map, bool show_boundaries,
                const BoundaryOverlayData* fallback_boundaries = nullptr)
{
    if (!show_boundaries) {
        save_image(platec_api_get_heightmap(simulation), filename, width, height, colors,
                   normalization_range, display_min, display_max, tone_map);
        return;
    }

    std::vector<float> copy(platec_api_get_heightmap(simulation),
                            platec_api_get_heightmap(simulation) +
                                static_cast<size_t>(width) * static_cast<size_t>(height));
    apply_normalization_range(copy, normalization_range);
    apply_display_mapping(copy, display_min, display_max, tone_map);

    std::vector<png_byte> rgb;
    if (colors) {
        renderImageColorsRgb(width, height, copy.data(), rgb);
    } else {
        renderImageGrayRgb(width, height, copy.data(), rgb);
    }

    BoundaryOverlayData current_boundaries;
    if (capture_boundary_overlay(simulation, width, height, current_boundaries)) {
        overlay_boundaries(current_boundaries, width, height, rgb);
    } else if (fallback_boundaries != nullptr && fallback_boundaries->valid()) {
        overlay_boundaries(*fallback_boundaries, width, height, rgb);
    }

    if (writeImageRgb(filename.string().c_str(), width, height, rgb.data(), "Plate Tectonics") != 0) {
        fail("failed to write image: " + display_path(filename));
    }
}

void create_gif(const fs::path& gif_path, const std::vector<fs::path>& frame_paths)
{
    if (frame_paths.empty()) {
        fail("cannot create GIF without frames");
    }

    fs::path script_path = fs::temp_directory_path() /
                           ("plate-tectonics-gif-" + std::to_string(time(nullptr)) + ".ps1");
    std::ofstream script(script_path, std::ios::binary);
    if (!script) {
        fail("failed to create temporary PowerShell script for GIF export");
    }

    script << "Add-Type -AssemblyName WindowsBase\n";
    script << "Add-Type -AssemblyName PresentationCore\n";
    script << "$encoder = [System.Windows.Media.Imaging.GifBitmapEncoder]::new()\n";
    script << "$delay = [UInt16]" << kGifDelayCs << "\n";
    script << "$frames = @(\n";
    for (size_t i = 0; i < frame_paths.size(); ++i) {
        script << "    '" << escape_powershell_single_quotes(display_path(frame_paths[i])) << "'";
        script << (i + 1 < frame_paths.size() ? ",\n" : "\n");
    }
    script << ")\n";
    script << "foreach ($framePath in $frames) {\n";
    script << "    $bitmap = [System.Windows.Media.Imaging.BitmapImage]::new()\n";
    script << "    $bitmap.BeginInit()\n";
    script << "    $bitmap.CacheOption = [System.Windows.Media.Imaging.BitmapCacheOption]::OnLoad\n";
    script << "    $bitmap.UriSource = [System.Uri]::new($framePath)\n";
    script << "    $bitmap.EndInit()\n";
    script << "    $metadata = [System.Windows.Media.Imaging.BitmapMetadata]::new('gif')\n";
    script << "    $metadata.SetQuery('/grctlext/Delay', $delay)\n";
    script << "    $metadata.SetQuery('/grctlext/Disposal', [byte]2)\n";
    script << "    $frame = [System.Windows.Media.Imaging.BitmapFrame]::Create($bitmap, $bitmap.Thumbnail, $metadata, $bitmap.ColorContexts)\n";
    script << "    $encoder.Frames.Add($frame)\n";
    script << "}\n";
    script << "$stream = [System.IO.File]::Open('"
           << escape_powershell_single_quotes(display_path(gif_path))
           << "', [System.IO.FileMode]::Create, [System.IO.FileAccess]::Write)\n";
    script << "try { $encoder.Save($stream) } finally { $stream.Dispose() }\n";
    script.close();

    const std::string command =
        "powershell -NoProfile -NonInteractive -ExecutionPolicy Bypass -File \"" +
        script_path.string() + "\"";

    const int rc = system(command.c_str());
    std::error_code ec;
    fs::remove(script_path, ec);

    if (rc != 0) {
        fail("failed to create GIF: " + display_path(gif_path));
    }
}

void delete_files(const std::vector<fs::path>& paths)
{
    for (const fs::path& path : paths) {
        std::error_code ec;
        const bool removed = fs::remove(path, ec);
        if (!removed || ec) {
            fail("failed to delete frame: " + display_path(path));
        }
    }
}

} // namespace

int main(int argc, char* argv[])
{
    Params params = fill_params(argc, argv);
    std::vector<float> input_heightmap;
    std::vector<uint16_t> input_metric_heightmap;
    std::vector<fs::path> gif_frames;
    std::vector<fs::path> gif_temp_files;
    std::vector<fs::path> step_outputs;
    BoundaryOverlayData last_boundary_state;
    NormalizationRange normalization_range;
    TopographyCodec::Metadata input_metadata;
    bool has_input_metadata = false;
    fs::path resolved_input_path;
    fs::path initial_png_output_path;
    fs::path final_png_output_path;
    fs::path gif_output_path;
    fs::path debug_output_path;
    fs::path initial_r16_output_path;
    fs::path final_r16_output_path;
    std::string initial_png_label;
    std::string final_png_label;
    std::string gif_label;
    std::string debug_label;
    std::string initial_r16_label;
    std::string final_r16_label;
    std::string run_id;

    if (params.input_mode != InputMode::Procedural) {
        const char* default_extension =
            params.input_mode == InputMode::MetricR16 ? ".r16" : ".png";
        resolved_input_path = resolve_input_path(params.input_path, default_extension);

        if (params.input_mode == InputMode::LegacyNormalizedPng) {
            int input_width = 0;
            int input_height = 0;
            if (readImageNormalized(resolved_input_path.string().c_str(), input_heightmap, input_width,
                                    input_height) != 0) {
                fail("failed to read input image: " + display_path(resolved_input_path));
            }
            params.width = static_cast<uint32_t>(input_width);
            params.height = static_cast<uint32_t>(input_height);
        } else if (params.input_mode == InputMode::MetricPng16) {
            int input_width = 0;
            int input_height = 0;
            if (readImageGray16(resolved_input_path.string().c_str(), input_metric_heightmap,
                                input_width, input_height) != 0) {
                fail("failed to read metric PNG16: " + display_path(resolved_input_path));
            }
            params.width = static_cast<uint32_t>(input_width);
            params.height = static_cast<uint32_t>(input_height);
            has_input_metadata = try_read_sidecar_metadata(resolved_input_path, input_metadata);
            if (has_input_metadata &&
                (input_metadata.width != params.width || input_metadata.height != params.height)) {
                fail("metric PNG16 sidecar dimensions do not match the image");
            }
        } else if (params.input_mode == InputMode::MetricR16) {
            has_input_metadata = try_read_sidecar_metadata(resolved_input_path, input_metadata);
            if (has_input_metadata) {
                if (params.custom_dimensions &&
                    (input_metadata.width != params.width || input_metadata.height != params.height)) {
                    fail("--dim does not match the .r16 sidecar dimensions");
                }
                params.width = input_metadata.width;
                params.height = input_metadata.height;
            } else if (!params.custom_dimensions) {
                fail("--input-r16 requires .r16.json metadata or --dim WIDTH HEIGHT");
            }

            if (readRawR16(resolved_input_path.string().c_str(), input_metric_heightmap,
                           static_cast<size_t>(params.width) * static_cast<size_t>(params.height)) != 0) {
                fail("failed to read metric .r16: " + display_path(resolved_input_path));
            }
        }

    }

    ensure_output_directories();
    const std::time_t started_at = std::time(nullptr);
    run_id = make_run_id(params, resolved_input_path, started_at);
    initial_png_output_path = kDefaultInputDir / (run_id + ".png");
    final_png_output_path = kDefaultOutputDir / (run_id + ".png");
    gif_output_path = kDefaultGifDir / (run_id + ".gif");
    debug_output_path = kDefaultOutputDir / (run_id + ".f32");
    initial_r16_output_path = kDefaultInputDir / (run_id + ".r16");
    final_r16_output_path = kDefaultOutputDir / (run_id + ".r16");
    initial_png_label = display_path(initial_png_output_path);
    final_png_label = display_path(final_png_output_path);
    gif_label = display_path(gif_output_path);
    debug_label = display_path(debug_output_path);
    initial_r16_label = display_path(initial_r16_output_path);
    final_r16_label = display_path(final_r16_output_path);

    if (params.legacy_filename_override) {
        fprintf(stderr,
                "warning: --filename is ignored for output naming; using %s instead\n",
                run_id.c_str());
    }

    const int32_t sea_level_override_m =
        (params.input_mode == InputMode::Procedural && params.has_sea_level_m)
            ? static_cast<int32_t>(params.sea_level_m)
            : TopographyCodec::kNoSeaLevelOverride;

    void* simulation = platec_api_create(params.seed, params.width, params.height, kSeaLevel,
                                         params.erosion_period, params.folding_ratio,
                                         params.aggregation_overlap_abs,
                                         params.aggregation_overlap_rel, params.cycles,
                                         params.plates,
                                         params.erosion_strength, params.landmass_rotation,
                                         params.rotation_strength,
                                         params.subduction_strength, sea_level_override_m,
                                         params.min_initial_height_m, params.max_initial_height_m);

    if (!input_heightmap.empty()) {
        platec_api_load_heightmap_raw(simulation, input_heightmap.data());
    } else if (!input_metric_heightmap.empty()) {
        uint16_t metric_sea_level_m = 0;
        if (has_input_metadata) {
            metric_sea_level_m = input_metadata.sea_level_m;
            if (params.has_sea_level_m && metric_sea_level_m != params.sea_level_m) {
                fprintf(stderr,
                        "warning: %s metadata sea_level_m=%u overrides CLI value %u\n",
                        display_path(resolved_input_path).c_str(), metric_sea_level_m,
                        params.sea_level_m);
            }
        } else if (params.has_sea_level_m) {
            metric_sea_level_m = params.sea_level_m;
        } else {
            metric_sea_level_m =
                infer_metric_sea_level_with_warning(input_metric_heightmap, resolved_input_path);
        }

        platec_api_load_heightmap_u16(simulation, input_metric_heightmap.data(), metric_sea_level_m);
    }

    const uint16_t active_sea_level_m = platec_api_get_sea_level_m(simulation);

    printf("Plate-tectonics simulation example\n");
    printf(" seed     : %u\n", params.seed);
    printf(" width    : %u\n", params.width);
    printf(" height   : %u\n", params.height);
    printf(" preview  : %s\n", params.colors ? "colors" : "grayscale");
    printf(" cycles   : %u\n", params.cycles);
    printf(" plates   : %u\n", params.plates);
    printf(" aggregate: abs %u, rel %.3f\n", params.aggregation_overlap_abs,
           params.aggregation_overlap_rel);
    printf(" folding  : %.3f\n", params.folding_ratio);
    printf(" erosion  : every %u updates, strength %.2f\n", params.erosion_period,
           params.erosion_strength);
    printf(" rotation : motion %.2f, landmass %.2f\n", params.rotation_strength,
           params.landmass_rotation);
    printf(" subduct  : %.2f\n", params.subduction_strength);
    printf(" sea lvl  : %u m\n", active_sea_level_m);
    if (params.input_mode == InputMode::Procedural) {
        printf(" init rng : %u..%u m\n", static_cast<unsigned>(params.min_initial_height_m),
               static_cast<unsigned>(params.max_initial_height_m));
    }
    printf(" normalize: %s\n", params.normalize_output ? "preview frames only" : "no");
    printf(" display  : [%.2f, %.2f], tone=%s\n", params.display_min, params.display_max,
           tone_map_name(params.tone_map));
    if (!resolved_input_path.empty()) {
        printf(" input    : %s\n", display_path(resolved_input_path).c_str());
    }
    printf(" init r16 : %s\n", initial_r16_label.c_str());
    printf(" init png16: %s\n", initial_png_label.c_str());
    printf(" final r16: %s\n", final_r16_label.c_str());
    printf(" final png16: %s\n", final_png_label.c_str());
    if (params.export_heightmap_f32) {
        printf(" debug f32: %s\n", debug_label.c_str());
    }
    if (params.make_gif) {
        printf(" gif      : %s\n", gif_label.c_str());
        printf(" keep png : %s\n", params.delete_gif_frames ? "no" : "yes");
    }
    if (params.show_boundaries) {
        printf(" bounds   : convergent=red divergent=blue transform=yellow\n");
    }
    if (params.step == 0) {
        printf(" step     : no\n");
    } else {
        printf(" step     : %u\n", params.step);
    }
    if (params.step != 0 || params.make_gif) {
        printf(" frames   : %s\n", display_path(kDefaultFramesDir).c_str());
    }
    printf("\n");

    normalization_range =
        capture_normalization_range(platec_api_get_heightmap(simulation),
                                    static_cast<int>(params.width * params.height),
                                    params.normalize_output, 0.0f, 1.0f);

    if (params.show_boundaries) {
        capture_boundary_overlay(simulation, static_cast<int>(params.width),
                                 static_cast<int>(params.height), last_boundary_state);
    }

    export_metric_heightmap(platec_api_get_heightmap(simulation), static_cast<int>(params.width),
                            static_cast<int>(params.height), active_sea_level_m,
                            initial_r16_output_path, initial_png_output_path);

    uint32_t saved_frame_count = 0;
    const auto save_frame = [&](const BoundaryOverlayData* fallback_boundaries) {
        const fs::path frame_path = kDefaultFramesDir / frame_file_name(run_id, saved_frame_count);
        save_image(simulation, frame_path, static_cast<int>(params.width),
                   static_cast<int>(params.height), params.colors, normalization_range,
                   params.display_min, params.display_max, params.tone_map,
                   params.show_boundaries, fallback_boundaries);
        ++saved_frame_count;
        return frame_path;
    };

    if (params.make_gif && params.step == 0) {
        const fs::path frame_path = save_frame(&last_boundary_state);
        gif_frames.push_back(frame_path);
        gif_temp_files.push_back(frame_path);
    }

    uint32_t step = 0;
    while (platec_api_is_finished(simulation) == 0) {
        ++step;
        const BoundaryOverlayData boundary_before_step = last_boundary_state;
        platec_api_step(simulation);

        if (params.make_gif && params.step == 0) {
            const fs::path frame_path = save_frame(&boundary_before_step);
            gif_frames.push_back(frame_path);
            gif_temp_files.push_back(frame_path);
        }

        if (params.step != 0 && step % params.step == 0) {
            const fs::path step_output_path = save_frame(&boundary_before_step);
            step_outputs.push_back(step_output_path);
            if (params.make_gif) {
                gif_frames.push_back(step_output_path);
            }
            printf(" * frame %u (step %u, filename %s)\n",
                   saved_frame_count - 1, step, display_path(step_output_path).c_str());
        }

        if (params.show_boundaries) {
            capture_boundary_overlay(simulation, static_cast<int>(params.width),
                                     static_cast<int>(params.height), last_boundary_state);
        }
    }

    export_metric_heightmap(platec_api_get_heightmap(simulation), static_cast<int>(params.width),
                            static_cast<int>(params.height), active_sea_level_m,
                            final_r16_output_path, final_png_output_path);

    if (params.export_heightmap_f32) {
        export_heightmap_f32(platec_api_get_heightmap(simulation), debug_output_path,
                             static_cast<int>(params.width), static_cast<int>(params.height));
    }

    if (params.make_gif) {
        if (params.step != 0 && (gif_frames.empty() || step % params.step != 0)) {
            const fs::path frame_path = save_frame(&last_boundary_state);
            gif_frames.push_back(frame_path);
            if (params.delete_gif_frames) {
                step_outputs.push_back(frame_path);
            }
        }
        create_gif(gif_output_path, gif_frames);
        if (params.delete_gif_frames) {
            if (params.step == 0) {
                delete_files(gif_temp_files);
            } else {
                delete_files(step_outputs);
            }
        }
    }

    printf(" * initial metric heightmap exported (filename %s)\n", initial_r16_label.c_str());
    printf(" * initial metric PNG16 exported (filename %s)\n", initial_png_label.c_str());
    printf(" * final metric heightmap exported (filename %s)\n", final_r16_label.c_str());
    printf(" * final metric PNG16 exported (filename %s)\n", final_png_label.c_str());
    if (params.export_heightmap_f32) {
        printf(" * debug float32 exported (filename %s)\n", debug_label.c_str());
    }
    if (params.make_gif) {
        printf(" * gif created (filename %s)\n", gif_label.c_str());
    }

    platec_api_destroy(simulation);
    return 0;
}
