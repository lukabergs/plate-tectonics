#ifndef NOMINMAX
#define NOMINMAX
#endif

#include "map_drawing.hpp"
#include "platecapi.hpp"
#include "sqrdmd.hpp"

#include <cmath>
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
constexpr uint16_t kGifDelayCs = 8;

const fs::path kRepoRoot(PLATE_TECTONICS_ROOT);
const fs::path kDefaultInputDir = kRepoRoot / "img" / "in";
const fs::path kDefaultOutputDir = kRepoRoot / "img" / "out";

struct Params {
    uint32_t seed;
    uint32_t width;
    uint32_t height;
    bool colors;
    bool make_gif;
    bool delete_gif_frames;
    bool show_boundaries;
    uint32_t cycles;
    uint32_t plates;
    uint32_t erosion_period;
    float erosion_strength;
    float rotation_strength;
    float landmass_rotation;
    uint32_t step;
    bool custom_dimensions;
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

float parse_nonnegative_float(const char* flag, const char* value)
{
    char* end = nullptr;
    const float parsed = strtof(value, &end);
    if (end == value || *end != '\0' || !std::isfinite(parsed) || parsed < 0.0f) {
        fail(std::string("invalid value for ") + flag + ": " + value);
    }
    return parsed;
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
    printf(" -i --input FILE     : load FILE from the given path or from %s\n",
           display_path(kDefaultInputDir).c_str());
    printf(" --dim WIDTH HEIGHT  : use the given width and height when no input image is used\n");
    printf(" --colors            : generate a colors map\n");
    printf(" --grayscale         : generate a grayscale map\n");
    printf(" --filename NAME     : output basename when no input image is used\n");
    printf(" --cycles N          : number of simulation cycles to run\n");
    printf(" --plates N          : number of tectonic plates to initialize\n");
    printf(" --erosion-period N  : number of updates between erosion passes\n");
    printf(" --erosion-strength X: scale erosion amount (0 disables)\n");
    printf(" --rotation-strength X: scale angular plate motion\n");
    printf(" --landmass-rotation X: scale visible crust rotation (0 disables)\n");
    printf(" --step X            : save intermediate maps every X steps\n");
    printf(" --gif               : export an animated GIF using --step sampling, or every update if --step is omitted\n");
    printf(" --no-steps          : delete GIF frame PNGs after the GIF is created\n");
    printf(" --show-boundaries   : overlay convergent/divergent/transform boundaries\n");
}

Params fill_params(int argc, char* argv[])
{
    srand(static_cast<unsigned int>(time(nullptr)));

    Params params = {
        static_cast<uint32_t>(rand()),
        600,
        400,
        true,
        false,
        false,
        false,
        kDefaultCycleCount,
        kDefaultPlateCount,
        kErosionPeriod,
        kDefaultErosionStrength,
        kDefaultRotationStrength,
        kDefaultLandmassRotation,
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

    if (!params.input_path.empty() && params.custom_dimensions) {
        fail("--dim cannot be used together with --input");
    }
    if (params.delete_gif_frames && !params.make_gif) {
        fail("--no-steps requires --gif");
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

std::string step_file_name(const std::string& stem, uint32_t run_index, uint32_t step)
{
    char buffer[128];
    snprintf(buffer, sizeof(buffer), "%s_%u_step_%u.png", stem.c_str(), run_index, step);
    return buffer;
}

std::string frame_file_name(const std::string& stem, uint32_t run_index, uint32_t frame_index)
{
    char buffer[128];
    snprintf(buffer, sizeof(buffer), "%s_%u_frame_%05u.png", stem.c_str(), run_index, frame_index);
    return buffer;
}

void save_image(const float* heightmap, const fs::path& filename, int width, int height, bool colors)
{
    std::vector<float> copy(heightmap,
                            heightmap + static_cast<size_t>(width) * static_cast<size_t>(height));
    normalize(copy.data(), width * height);

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
                bool show_boundaries, const BoundaryOverlayData* fallback_boundaries = nullptr)
{
    if (!show_boundaries) {
        save_image(platec_api_get_heightmap(simulation), filename, width, height, colors);
        return;
    }

    std::vector<float> copy(platec_api_get_heightmap(simulation),
                            platec_api_get_heightmap(simulation) +
                                static_cast<size_t>(width) * static_cast<size_t>(height));
    normalize(copy.data(), width * height);

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
    std::vector<fs::path> gif_frames;
    std::vector<fs::path> gif_temp_files;
    std::vector<fs::path> step_outputs;
    BoundaryOverlayData last_boundary_state;
    fs::path resolved_input_path;
    fs::path final_output_path;
    fs::path gif_output_path;
    std::string output_label;
    std::string gif_label;
    std::string run_stem;
    uint32_t run_index = 0;

    if (!params.input_path.empty()) {
        resolved_input_path = resolve_input_path(params.input_path);

        int input_width = 0;
        int input_height = 0;
        if (readImageNormalized(resolved_input_path.string().c_str(), input_heightmap, input_width,
                                input_height) != 0) {
            fail("failed to read input image: " + display_path(resolved_input_path));
        }

        params.width = static_cast<uint32_t>(input_width);
        params.height = static_cast<uint32_t>(input_height);

        fs::create_directories(kDefaultOutputDir);
        run_stem = resolved_input_path.stem().string();
        run_index = next_run_index(kDefaultOutputDir, run_stem);
        final_output_path = kDefaultOutputDir / (run_stem + "_" + std::to_string(run_index) + ".png");
        gif_output_path = kDefaultOutputDir / (run_stem + "_" + std::to_string(run_index) + ".gif");
        output_label = display_path(final_output_path);
        gif_label = display_path(gif_output_path);
    } else {
        run_stem = params.filename;
        final_output_path = params.filename + ".png";
        gif_output_path = params.filename + ".gif";
        output_label = display_path(final_output_path);
        gif_label = display_path(gif_output_path);
    }

    printf("Plate-tectonics simulation example\n");
    printf(" seed     : %u\n", params.seed);
    printf(" width    : %u\n", params.width);
    printf(" height   : %u\n", params.height);
    printf(" map      : %s\n", params.colors ? "colors" : "grayscale");
    printf(" cycles   : %u\n", params.cycles);
    printf(" plates   : %u\n", params.plates);
    printf(" erosion  : every %u updates, strength %.2f\n", params.erosion_period,
           params.erosion_strength);
    printf(" rotation : motion %.2f, landmass %.2f\n", params.rotation_strength,
           params.landmass_rotation);
    if (!resolved_input_path.empty()) {
        printf(" input    : %s\n", display_path(resolved_input_path).c_str());
    }
    printf(" output   : %s\n", output_label.c_str());
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
    printf("\n");

    void* simulation = platec_api_create(params.seed, params.width, params.height, kSeaLevel,
                                         params.erosion_period, kFoldingRatio, kAggrOverlapAbs,
                                         kAggrOverlapRel, params.cycles, params.plates,
                                         params.erosion_strength, params.landmass_rotation,
                                         params.rotation_strength);

    if (!input_heightmap.empty()) {
        platec_api_load_heightmap(simulation, input_heightmap.data(), kSeaLevel);
    }

    if (params.show_boundaries) {
        capture_boundary_overlay(simulation, static_cast<int>(params.width),
                                 static_cast<int>(params.height), last_boundary_state);
    }

    if (params.make_gif && params.step == 0) {
        const fs::path frame_path =
            (!resolved_input_path.empty() ? kDefaultOutputDir : fs::path(".")) /
            frame_file_name(run_stem, run_index, 0);
        save_image(simulation, frame_path, static_cast<int>(params.width),
                   static_cast<int>(params.height), params.colors, params.show_boundaries,
                   &last_boundary_state);
        gif_frames.push_back(frame_path);
        gif_temp_files.push_back(frame_path);
    }

    uint32_t step = 0;
    while (platec_api_is_finished(simulation) == 0) {
        ++step;
        const BoundaryOverlayData boundary_before_step = last_boundary_state;
        platec_api_step(simulation);

        if (params.make_gif && params.step == 0) {
            const fs::path frame_path =
                (!resolved_input_path.empty() ? kDefaultOutputDir : fs::path(".")) /
                frame_file_name(run_stem, run_index, step);
            save_image(simulation, frame_path, static_cast<int>(params.width),
                       static_cast<int>(params.height), params.colors, params.show_boundaries,
                       &boundary_before_step);
            gif_frames.push_back(frame_path);
            gif_temp_files.push_back(frame_path);
        }

        if (params.step != 0 && step % params.step == 0) {
            const fs::path step_output_path =
                (!resolved_input_path.empty() ? kDefaultOutputDir : fs::path(".")) /
                step_file_name(run_stem, run_index, step);
            save_image(simulation, step_output_path, static_cast<int>(params.width),
                       static_cast<int>(params.height), params.colors, params.show_boundaries,
                       &boundary_before_step);
            step_outputs.push_back(step_output_path);
            if (params.make_gif) {
                gif_frames.push_back(step_output_path);
            }
            printf(" * step %u (filename %s)\n", step, display_path(step_output_path).c_str());
        }

        if (params.show_boundaries) {
            capture_boundary_overlay(simulation, static_cast<int>(params.width),
                                     static_cast<int>(params.height), last_boundary_state);
        }
    }

    save_image(simulation, final_output_path, static_cast<int>(params.width),
               static_cast<int>(params.height), params.colors, params.show_boundaries,
               &last_boundary_state);

    if (params.make_gif) {
        if (params.step != 0 && (gif_frames.empty() || step % params.step != 0)) {
            gif_frames.push_back(final_output_path);
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

    printf(" * simulation completed (filename %s)\n", output_label.c_str());
    if (params.make_gif) {
        printf(" * gif created (filename %s)\n", gif_label.c_str());
    }

    platec_api_destroy(simulation);
    return 0;
}
