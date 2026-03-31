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

#include "lithosphere.hpp"
#include "noise.hpp"
#include "plate.hpp"
#include "simplexnoise.hpp"
#include "sqrdmd.hpp"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <vector>

#define BOOL_REGENERATE_CRUST 1

using namespace std;

static const float SUBDUCT_RATIO = 0.5f;

static const float BUOYANCY_BONUS_X = 3;
static const uint32_t MAX_BUOYANCY_AGE = 20;
static const float MULINV_MAX_BUOYANCY_AGE = 1.0f / (float)MAX_BUOYANCY_AGE;

static const float RESTART_ENERGY_RATIO = 0.15f;
static const float RESTART_SPEED_LIMIT = 2.0f;
static const uint32_t RESTART_ITERATIONS = 600;
static const uint32_t NO_COLLISION_TIME_LIMIT = 10;

namespace {

uint32_t seam_hash(uint32_t x, uint32_t y, uint32_t iteration, uint32_t salt)
{
    uint32_t h = x * 73856093U;
    h ^= y * 19349663U;
    h ^= iteration * 83492791U;
    h ^= salt * 2654435761U;
    h = (h ^ (h >> 13U)) * 1274126177U;
    return h ^ (h >> 16U);
}

float seam_noise(uint32_t x, uint32_t y, uint32_t iteration, uint32_t salt)
{
    const uint32_t h = seam_hash(x, y, iteration, salt);
    return (static_cast<float>(h & 1023U) / 511.5f) - 1.0f;
}

void choose_regenerated_gap_fill(const WorldDimension& world_dimension, const IndexMap& current_imap,
                                 const IndexMap& previous_imap, const HeightMap& current_hmap,
                                 uint32_t x, uint32_t y, uint32_t num_plates,
                                 uint32_t iteration_count, std::vector<float>& owner_scores,
                                 std::vector<float>& owner_height_sums,
                                 std::vector<uint32_t>& owner_height_counts,
                                 uint32_t* owner_out, float* height_out)
{
    ASSERT(owner_out != nullptr && height_out != nullptr, "gap fill outputs must not be null");
    ASSERT(owner_scores.size() >= num_plates, "owner score buffer too small");
    ASSERT(owner_height_sums.size() >= num_plates, "owner height buffer too small");
    ASSERT(owner_height_counts.size() >= num_plates, "owner count buffer too small");

    std::fill(owner_scores.begin(), owner_scores.end(), 0.0f);
    std::fill(owner_height_sums.begin(), owner_height_sums.end(), 0.0f);
    std::fill(owner_height_counts.begin(), owner_height_counts.end(), 0U);

    const uint32_t world_index = world_dimension.indexOf(x, y);
    const uint32_t previous_owner = previous_imap[world_index];

    auto add_candidate = [&](uint32_t plate, float weight, float local_height, bool collect_height) {
        if (plate >= num_plates || weight <= 0.0f) {
            return;
        }

        owner_scores[plate] += weight;
        if (collect_height && local_height > 0.0f) {
            owner_height_sums[plate] += local_height;
            owner_height_counts[plate] += 1U;
        }
    };

    for (int dy = -2; dy <= 2; ++dy) {
        for (int dx = -2; dx <= 2; ++dx) {
            if (dx == 0 && dy == 0) {
                continue;
            }

            const int abs_dx = std::abs(dx);
            const int abs_dy = std::abs(dy);
            const int distance = abs_dx > abs_dy ? abs_dx : abs_dy;
            const bool cardinal = dx == 0 || dy == 0;
            const float current_weight =
                (distance == 1) ? (cardinal ? 7.0f : 5.0f) : (cardinal ? 3.5f : 2.2f);
            const uint32_t neighbour_index =
                world_dimension.normalizedIndexOf(static_cast<uint32_t>(x + dx),
                                                  static_cast<uint32_t>(y + dy));
            add_candidate(current_imap[neighbour_index], current_weight,
                          current_hmap[neighbour_index], true);
            add_candidate(previous_imap[neighbour_index], current_weight * 0.35f, 0.0f, false);
        }
    }

    if (previous_owner < num_plates) {
        owner_scores[previous_owner] += 1.4f;
    }

    uint32_t best_owner = previous_owner < num_plates ? previous_owner : 0xFFFFFFFFU;
    float best_score = -FLT_MAX;
    for (uint32_t plate = 0; plate < num_plates; ++plate) {
        if (owner_scores[plate] <= 0.0f) {
            continue;
        }

        float score = owner_scores[plate];
        if (previous_owner == plate) {
            score += 0.7f;
        }
        score += 0.45f * seam_noise(x, y, iteration_count, plate + 1U);

        if (owner_height_counts[plate] > 0U) {
            const float local_mean = owner_height_sums[plate] /
                                     static_cast<float>(owner_height_counts[plate]);
            const float capped_mean =
                local_mean < CONTINENTAL_BASE ? local_mean : CONTINENTAL_BASE;
            const float ocean_bias =
                std::clamp((CONTINENTAL_BASE - capped_mean) / CONTINENTAL_BASE, 0.0f, 1.0f);
            score += 0.8f * ocean_bias;
        }

        if (score > best_score) {
            best_score = score;
            best_owner = plate;
        }
    }

    if (best_owner >= num_plates) {
        best_owner = previous_owner;
    }
    if (best_owner >= num_plates) {
        for (int dy = -1; dy <= 1 && best_owner >= num_plates; ++dy) {
            for (int dx = -1; dx <= 1; ++dx) {
                if (dx == 0 && dy == 0) {
                    continue;
                }
                const uint32_t neighbour_index =
                    world_dimension.normalizedIndexOf(static_cast<uint32_t>(x + dx),
                                                      static_cast<uint32_t>(y + dy));
                if (current_imap[neighbour_index] < num_plates) {
                    best_owner = current_imap[neighbour_index];
                    break;
                }
                if (previous_imap[neighbour_index] < num_plates) {
                    best_owner = previous_imap[neighbour_index];
                }
            }
        }
    }

    if (best_owner >= num_plates) {
        *owner_out = previous_owner;
        *height_out = OCEANIC_BASE * BUOYANCY_BONUS_X;
        return;
    }

    float regenerated_height = OCEANIC_BASE * BUOYANCY_BONUS_X;
    float local_height_sum = 0.0f;
    uint32_t local_height_count = 0U;
    for (int dy = -1; dy <= 1; ++dy) {
        for (int dx = -1; dx <= 1; ++dx) {
            if (dx == 0 && dy == 0) {
                continue;
            }

            const uint32_t neighbour_index =
                world_dimension.normalizedIndexOf(static_cast<uint32_t>(x + dx),
                                                  static_cast<uint32_t>(y + dy));
            if (current_imap[neighbour_index] != best_owner) {
                continue;
            }

            local_height_sum += std::clamp(current_hmap[neighbour_index], OCEANIC_BASE,
                                           OCEANIC_BASE * BUOYANCY_BONUS_X);
            local_height_count += 1U;
        }
    }

    if (local_height_count > 0U) {
        regenerated_height = regenerated_height * 0.55f +
                             (local_height_sum / static_cast<float>(local_height_count)) * 0.45f;
    }

    regenerated_height *=
        0.88f + 0.16f * seam_noise(x, y, iteration_count, best_owner + 17U);
    regenerated_height = std::clamp(regenerated_height, OCEANIC_BASE * 1.6f,
                                    OCEANIC_BASE * 2.6f);

    *owner_out = best_owner;
    *height_out = regenerated_height;
}

} // namespace

uint32_t findBound(const uint32_t* map, uint32_t length, uint32_t x0, uint32_t y0, int dx, int dy);
uint32_t findPlate(plate** plates, float x, float y, uint32_t num_plates);

WorldPoint lithosphere::randomPosition() {
    return WorldPoint(_randsource.next() % _worldDimension.getWidth(),
                      _randsource.next() % _worldDimension.getHeight(), _worldDimension);
}

void lithosphere::initializeTopography(const float* height_samples,
                                       const WorldDimension& source_dimension, float sea_level,
                                       bool preserve_relief) {
    if (height_samples == nullptr) {
        throw runtime_error("Heightmap input cannot be null");
    }
    if (source_dimension.getWidth() < _worldDimension.getWidth() ||
        source_dimension.getHeight() < _worldDimension.getHeight()) {
        throw runtime_error("Source heightmap dimensions must cover the world dimensions");
    }

    const uint32_t area = source_dimension.getArea();
    std::vector<float> normalized(height_samples, height_samples + area);

    const auto [lowest_it, highest_it] = std::minmax_element(normalized.begin(), normalized.end());
    const float lowest = *lowest_it;
    const float highest = *highest_it;
    if (highest > lowest) {
        const float inv_range = 1.0f / (highest - lowest);
        for (uint32_t i = 0; i < area; ++i) {
            normalized[i] = (normalized[i] - lowest) * inv_range;
        }
    } else {
        std::fill(normalized.begin(), normalized.end(), 0.0f);
    }

    float sea_threshold = 0.5f;
    float th_step = 0.5f;
    while (th_step > 0.01f) {
        uint32_t count = 0;
        for (uint32_t i = 0; i < area; ++i) {
            count += (normalized[i] < sea_threshold);
        }

        th_step *= 0.5f;
        if (count / static_cast<float>(area) < sea_level) {
            sea_threshold += th_step;
        } else {
            sea_threshold -= th_step;
        }
    }

    if (preserve_relief) {
        const float sea_range = (sea_threshold > 0.0001f) ? sea_threshold : 0.0001f;
        const float land_range =
            (1.0f - sea_threshold > 0.0001f) ? (1.0f - sea_threshold) : 0.0001f;

        for (uint32_t i = 0; i < area; ++i) {
            if (normalized[i] <= sea_threshold) {
                const float ocean_relative = std::pow(normalized[i] / sea_range, 1.15f);
                normalized[i] = 0.02f + ocean_relative * 0.33f;
            } else {
                const float land_relative =
                    std::pow((normalized[i] - sea_threshold) / land_range, 0.85f);
                normalized[i] = CONTINENTAL_BASE + land_relative;
            }
        }
    } else {
        for (uint32_t i = 0; i < area; ++i) {
            normalized[i] =
                (normalized[i] > sea_threshold) ? normalized[i] + CONTINENTAL_BASE : OCEANIC_BASE;
        }
    }

    for (uint32_t y = 0; y < _worldDimension.getHeight(); ++y) {
        memcpy(&hmap[_worldDimension.lineIndex(y)], &normalized[source_dimension.lineIndex(y)],
               _worldDimension.getWidth() * sizeof(float));
    }
}

void lithosphere::initializePlates() {
    collisions.resize(max_plates);
    subductions.resize(max_plates);

    plates = new plate*[max_plates];
    for (uint32_t i = 0; i < max_plates; ++i) {
        plate_areas[i].border.reserve(8);
    }
    createPlates();
}

void lithosphere::createNoise(float* tmp, const WorldDimension& tmpDim, bool useSimplex) {
    ::createNoise(tmp, tmpDim, _randsource, useSimplex);
}

void lithosphere::createSlowNoise(float* tmp, const WorldDimension& tmpDim) {
    ::createSlowNoise(tmp, tmpDim, _randsource);
}

lithosphere::lithosphere(long seed, uint32_t width, uint32_t height, float sea_level,
                         uint32_t _erosion_period, float _folding_ratio, uint32_t aggr_ratio_abs,
                         float aggr_ratio_rel, uint32_t num_cycles,
                         uint32_t _max_plates) noexcept(false)
    : hmap(width, height), imap(width, height), prev_imap(width, height), amap(width, height),
      plates(nullptr), plate_areas(_max_plates), plate_indices_found(_max_plates),
      aggr_overlap_abs(aggr_ratio_abs), aggr_overlap_rel(aggr_ratio_rel), cycle_count(0),
      erosion_period(_erosion_period), folding_ratio(_folding_ratio), iter_count(0),
      max_cycles(num_cycles), max_plates(_max_plates), num_plates(0),
      _worldDimension(width, height), _randsource(seed), imported_heightmap_mode(false),
      _steps(0) {
    if (width < 5 || height < 5) {
        throw runtime_error("Width and height should be >=5");
    }

    WorldDimension tmpDim = WorldDimension(width + 1, height + 1);
    const uint32_t A = tmpDim.getArea();
    float* tmp = new float[A];

    createSlowNoise(tmp, tmpDim);
    initializeTopography(tmp, tmpDim, sea_level, false);
    delete[] tmp;
    initializePlates();
}

lithosphere::lithosphere(long seed, uint32_t width, uint32_t height, const float* heightmap,
                         float sea_level, uint32_t _erosion_period, float _folding_ratio,
                         uint32_t aggr_ratio_abs, float aggr_ratio_rel, uint32_t num_cycles,
                         uint32_t _max_plates) noexcept(false)
    : hmap(width, height), imap(width, height), prev_imap(width, height), amap(width, height),
      plates(nullptr), plate_areas(_max_plates), plate_indices_found(_max_plates),
      aggr_overlap_abs(aggr_ratio_abs), aggr_overlap_rel(aggr_ratio_rel), cycle_count(0),
      erosion_period(_erosion_period), folding_ratio(_folding_ratio), iter_count(0),
      max_cycles(num_cycles), max_plates(_max_plates), num_plates(0),
      _worldDimension(width, height), _randsource(seed), imported_heightmap_mode(true),
      _steps(0) {
    if (width < 5 || height < 5) {
        throw runtime_error("Width and height should be >=5");
    }

    initializeTopography(heightmap, _worldDimension, sea_level, true);
    initializePlates();
}

lithosphere::~lithosphere() throw() {
    clearPlates();
    delete[] plates;
    plates = 0;
}

void lithosphere::clearPlates() {
    for (uint32_t i = 0; i < num_plates; i++) {
        delete plates[i];
    }
    num_plates = 0;
}

void lithosphere::growPlates() {
    // "Grow" plates from their origins until surface is fully populated.
    uint32_t max_border = 1;
    uint32_t i;
    while (max_border) {
        for (max_border = i = 0; i < num_plates; ++i) {
            plateArea& area = plate_areas[i];
            const uint32_t N = (uint32_t)area.border.size();
            max_border = max_border > N ? max_border : N;

            if (N == 0) {
                continue;
            }
            const uint32_t j = _randsource.next() % N;
            const uint32_t p = area.border[j];
            const uint32_t cy = _worldDimension.yFromIndex(p);
            const uint32_t cx = _worldDimension.xFromIndex(p);

            const uint32_t lft = cx > 0 ? cx - 1 : _worldDimension.getWidth() - 1;
            const uint32_t rgt = cx < _worldDimension.getWidth() - 1 ? cx + 1 : 0;
            const uint32_t top = cy > 0 ? cy - 1 : _worldDimension.getHeight() - 1;
            const uint32_t btm = cy < _worldDimension.getHeight() - 1 ? cy + 1 : 0;

            const uint32_t n = top * _worldDimension.getWidth() + cx; // North.
            const uint32_t s = btm * _worldDimension.getWidth() + cx; // South.
            const uint32_t w = cy * _worldDimension.getWidth() + lft; // West.
            const uint32_t e = cy * _worldDimension.getWidth() + rgt; // East.

            if (imap[n] >= num_plates) {
                imap[n] = i;
                area.border.push_back(n);

                if (area.top == _worldDimension.yMod(top + 1)) {
                    area.top = top;
                    area.hgt++;
                }
            }

            if (imap[s] >= num_plates) {
                imap[s] = i;
                area.border.push_back(s);

                if (btm == _worldDimension.yMod(area.btm + 1)) {
                    area.btm = btm;
                    area.hgt++;
                }
            }

            if (imap[w] >= num_plates) {
                imap[w] = i;
                area.border.push_back(w);

                if (area.lft == _worldDimension.xMod(lft + 1)) {
                    area.lft = lft;
                    area.wdt++;
                }
            }

            if (imap[e] >= num_plates) {
                imap[e] = i;
                area.border.push_back(e);

                if (rgt == _worldDimension.xMod(area.rgt + 1)) {
                    area.rgt = rgt;
                    area.wdt++;
                }
            }

            // Overwrite processed point with unprocessed one.
            area.border[j] = area.border.back();
            area.border.pop_back();
        }
    }
}

void lithosphere::createPlates() {
    try {
        const uint32_t map_area = _worldDimension.getArea();
        num_plates = max_plates;

        // Initialize "Free plate center position" lookup table.
        // This way two plate centers will never be identical.
        for (uint32_t i = 0; i < map_area; ++i)
            imap[i] = i;

        // Select N plate centers from the global map.

        for (uint32_t i = 0; i < num_plates; ++i) {
            plateArea& area = plate_areas[i];

            // Randomly select an unused plate origin.
            const uint32_t p = imap[(uint32_t)_randsource.next() % (map_area - i)];
            const uint32_t y = _worldDimension.yFromIndex(p);
            const uint32_t x = _worldDimension.xFromIndex(p);

            area.lft = area.rgt = x; // Save origin...
            area.top = area.btm = y;
            area.wdt = area.hgt = 1;

            area.border.clear();
            area.border.push_back(p); // ...and mark it as border.

            // Overwrite used entry with last unused entry in array.
            imap[p] = imap[map_area - i - 1];
        }

        imap.set_all(0xFFFFFFFF);

        growPlates();

        // check all the points of the map are owned
        for (uint32_t i = 0; i < map_area; i++) {
            ASSERT(imap[i] < num_plates, "A point was not assigned to any plate");
        }

        // Extract and create plates from initial terrain.
        for (uint32_t i = 0; i < num_plates; ++i) {
            plateArea& area = plate_areas[i];

            area.wdt = _worldDimension.xCap(area.wdt);
            area.hgt = _worldDimension.yCap(area.hgt);

            const uint32_t x0 = area.lft;
            const uint32_t x1 = 1 + x0 + area.wdt;
            const uint32_t y0 = area.top;
            const uint32_t y1 = 1 + y0 + area.hgt;
            const uint32_t width = x1 - x0;
            const uint32_t height = y1 - y0;
            float* pmap = new float[width * height];

            // Copy plate's height data from global map into local map.
            for (uint32_t y = y0, j = 0; y < y1; ++y) {
                for (uint32_t x = x0; x < x1; ++x, ++j) {
                    uint32_t k = _worldDimension.normalizedIndexOf(x, y);
                    pmap[j] = hmap[k] * (imap[k] == i);
                }
            }
            // Create plate.
            // MK: The pmap array becomes owned by map, do not delete it
            plates[i] =
                new plate(_randsource.next(), pmap, width, height, x0, y0, i, _worldDimension);
        }

        iter_count = num_plates + MAX_BUOYANCY_AGE;
        peak_Ek = 0;
        last_coll_count = 0;

    } catch (const exception& e) {
        string msg = "Problem during createPlates: ";
        msg = msg + e.what();
        throw runtime_error(msg.c_str());
    }
}

uint32_t lithosphere::getPlateCount() const throw() {
    return num_plates;
}

const uint32_t* lithosphere::getAgeMap() const throw() {
    return amap.raw_data();
}

float* lithosphere::getTopography() const throw() {
    return hmap.raw_data();
}

bool lithosphere::isFinished() const {
    return getPlateCount() == 0;
}

// At least two plates are at same location.
// Move some crust from the SMALLER plate onto LARGER one.
void lithosphere::resolveJuxtapositions(const uint32_t& i, const uint32_t& j, const uint32_t& k,
                                        const uint32_t& x_mod, const uint32_t& y_mod,
                                        const float*& this_map, const uint32_t*& this_age,
                                        uint32_t& continental_collisions) {
    ASSERT(i < num_plates, "Given invalid plate index");

    // Record collisions to both plates. This also creates
    // continent segment at the collided location to plates.
    uint32_t this_area = plates[i]->addCollision(x_mod, y_mod);
    uint32_t prev_area = plates[imap[k]]->addCollision(x_mod, y_mod);

    if (this_area < prev_area) {
        plateCollision coll(imap[k], x_mod, y_mod, this_map[j] * folding_ratio);

        // Give some...
        hmap[k] += coll.crust;
        plates[imap[k]]->setCrust(x_mod, y_mod, hmap[k], this_age[j]);

        // And take some.
        plates[i]->setCrust(x_mod, y_mod, this_map[j] * (1.0f - folding_ratio), this_age[j]);

        // Add collision to the earlier plate's list.
        collisions[i].push_back(coll);
        ++continental_collisions;
    } else {
        plateCollision coll(i, x_mod, y_mod, hmap[k] * folding_ratio);

        plates[i]->setCrust(x_mod, y_mod, this_map[j] + coll.crust, amap[k]);

        plates[imap[k]]->setCrust(x_mod, y_mod, hmap[k] * (1.0f - folding_ratio), amap[k]);

        collisions[imap[k]].push_back(coll);
        ++continental_collisions;

        // Give the location to the larger plate.
        hmap[k] = this_map[j];
        imap[k] = i;
        amap[k] = this_age[j];
    }
}

// Update height and plate index maps.
// Doing it plate by plate is much faster than doing it index wise:
// Each plate's map's memory area is accessed sequentially and only
// once as opposed to calculating "num_plates" indices within plate
// maps in order to find out which plate(s) own current location.
void lithosphere::updateHeightAndPlateIndexMaps(uint32_t& oceanic_collisions,
                                                uint32_t& continental_collisions) {
    uint32_t world_width = _worldDimension.getWidth();
    uint32_t world_height = _worldDimension.getHeight();
    hmap.set_all(0);
    imap.set_all(0xFFFFFFFF);
    for (uint32_t i = 0; i < num_plates; ++i) {
        const uint32_t x0 = plates[i]->getLeftAsUint();
        const uint32_t y0 = plates[i]->getTopAsUint();
        const uint32_t x1 = x0 + plates[i]->getWidth();
        const uint32_t y1 = y0 + plates[i]->getHeight();

        const float* this_map;
        const uint32_t* this_age;
        plates[i]->getMap(&this_map, &this_age);

        uint32_t x_mod_start = (x0 + world_width) % world_width;
        uint32_t y_mod = (y0 + world_height) % world_height;

        // Copy first part of plate onto world map.
        // MK: These loops are ugly, but using modulus in here is a hog
        for (uint32_t y = y0, j = 0; y < y1;
             ++y, y_mod = ++y_mod >= world_height ? y_mod - world_height : y_mod) {
            const uint32_t y_width = y_mod * world_width;
            uint32_t x_mod = x_mod_start;

            for (uint32_t x = x0; x < x1;
                 ++x, ++j, x_mod = ++x_mod >= world_width ? x_mod - world_width : x_mod) {
                const uint32_t k = x_mod + y_width;

                if (this_map[j] < 2 * FLT_EPSILON) // No crust here...
                    continue;

                if (imap[k] >= num_plates) // No one here yet?
                {
                    // This plate becomes the "owner" of current location
                    // if it is the first plate to have crust on it.
                    hmap[k] = this_map[j];
                    imap[k] = i;
                    amap[k] = this_age[j];

                    continue;
                }

                // DO NOT ACCEPT HEIGHT EQUALITY! Equality leads to subduction
                // of shore that 's barely above sea level. It's a lot less
                // serious problem to treat very shallow waters as continent...
                const bool prev_is_oceanic = hmap[k] < CONTINENTAL_BASE;
                const bool this_is_oceanic = this_map[j] < CONTINENTAL_BASE;

                const uint32_t prev_timestamp = plates[imap[k]]->getCrustTimestamp(x_mod, y_mod);
                const uint32_t this_timestamp = this_age[j];
                const bool prev_is_buoyant =
                    (hmap[k] > this_map[j]) || ((hmap[k] + 2 * FLT_EPSILON > this_map[j]) &&
                                                (hmap[k] < 2 * FLT_EPSILON + this_map[j]) &&
                                                (prev_timestamp >= this_timestamp));

                // Handle subduction of oceanic crust as special case.
                if (this_is_oceanic && prev_is_buoyant) {
                    // This plate will be the subducting one.
                    // The level of effect that subduction has
                    // is directly related to the amount of water
                    // on top of the subducting plate.
                    const float sediment = SUBDUCT_RATIO * OCEANIC_BASE *
                                           (CONTINENTAL_BASE - this_map[j]) / CONTINENTAL_BASE;

                    // Save collision to the receiving plate's list.
                    plateCollision coll(i, x_mod, y_mod, sediment);
                    subductions[imap[k]].push_back(coll);
                    ++oceanic_collisions;

                    // Remove subducted oceanic lithosphere from plate.
                    // This is crucial for
                    // a) having correct amount of colliding crust (below)
                    // b) protecting subducted locations from receiving
                    //    crust from other subductions/collisions.
                    plates[i]->setCrust(x_mod, y_mod, this_map[j] - OCEANIC_BASE, this_timestamp);

                    if (this_map[j] <= 0)
                        continue; // Nothing more to collide.
                } else if (prev_is_oceanic) {
                    const float sediment = SUBDUCT_RATIO * OCEANIC_BASE *
                                           (CONTINENTAL_BASE - hmap[k]) / CONTINENTAL_BASE;

                    plateCollision coll(imap[k], x_mod, y_mod, sediment);
                    subductions[i].push_back(coll);
                    ++oceanic_collisions;

                    plates[imap[k]]->setCrust(x_mod, y_mod, hmap[k] - OCEANIC_BASE, prev_timestamp);
                    hmap[k] -= OCEANIC_BASE;

                    if (hmap[k] <= 0) {
                        imap[k] = i;
                        hmap[k] = this_map[j];
                        amap[k] = this_age[j];

                        continue;
                    }
                }

                resolveJuxtapositions(i, j, k, x_mod, y_mod, this_map, this_age,
                                      continental_collisions);
            }
        }
    }
}

void lithosphere::updateCollisions() {
    for (uint32_t i = 0; i < num_plates; ++i) {
        for (uint32_t j = 0; j < collisions[i].size(); ++j) {
            const plateCollision& coll = collisions[i][j];
            uint32_t coll_count, coll_count_i, coll_count_j;
            float coll_ratio, coll_ratio_i, coll_ratio_j;

            ASSERT(i != coll.index, "when colliding: SRC == DEST!");

            // Collision causes friction. Apply it to both plates.
            plates[i]->applyFriction(coll.crust);
            plates[coll.index]->applyFriction(coll.crust);

            plates[i]->getCollisionInfo(coll.wx, coll.wy, &coll_count_i, &coll_ratio_i);
            plates[coll.index]->getCollisionInfo(coll.wx, coll.wy, &coll_count_j, &coll_ratio_j);

            // Find the minimum count of collisions between two
            // continents on different plates.
            // It's minimum because large plate will get collisions
            // from all over whereas smaller plate will get just
            // a few. It's those few that matter between these two
            // plates, not what the big plate has with all the
            // other plates around it.
            coll_count = coll_count_i;
            coll_count -= (coll_count - coll_count_j) & -(coll_count > coll_count_j);

            // Find maximum amount of collided surface area between
            // two continents on different plates.
            // Like earlier, it's the "experience" of the smaller
            // plate that matters here.
            coll_ratio = coll_ratio_i;
            coll_ratio += (coll_ratio_j - coll_ratio) * (coll_ratio_j > coll_ratio);

            if ((coll_count > aggr_overlap_abs) | (coll_ratio > aggr_overlap_rel)) {
                float amount = plates[i]->aggregateCrust(plates[coll.index], coll.wx, coll.wy);

                // Calculate new direction and speed for the
                // merged plate system, that is, for the
                // receiving plate!
                plates[coll.index]->collide(*plates[i], amount);
            }
        }

        collisions[i].clear();
    }
}

// Remove empty plates from the system.
void lithosphere::removeEmptyPlates() {
    for (uint32_t i = 0; i < num_plates; ++i) {
        if (num_plates == 1)
            puts("ONLY ONE PLATE LEFT!");
        else if (plate_indices_found[i] == 0) {
            delete plates[i];
            plates[i] = plates[num_plates - 1];
            plate_indices_found[i] = plate_indices_found[num_plates - 1];

            // Life is seldom as simple as seems at first.
            // Replace the moved plate's index in the index map
            // to match its current position in the array!
            for (uint32_t j = 0; j < _worldDimension.getArea(); ++j)
                if (imap[j] == num_plates - 1)
                    imap[j] = i;

            --num_plates;
            --i;
        }
    }
}

void lithosphere::update() {
    try {
        _steps++;
        float totalVelocity = 0;
        float systemKineticEnergy = 0;

        for (uint32_t i = 0; i < num_plates; ++i) {
            totalVelocity += plates[i]->getVelocity();
            systemKineticEnergy += plates[i]->getMomentum();
        }

        if (systemKineticEnergy > peak_Ek) {
            peak_Ek = systemKineticEnergy;
        }

        // If there's no continental collisions during past iterations,
        // then interesting activity has ceased and we should restart.
        // Also if the simulation has been going on for too long already,
        // restart, because interesting stuff has most likely ended.
        if (totalVelocity < RESTART_SPEED_LIMIT ||
            systemKineticEnergy / peak_Ek < RESTART_ENERGY_RATIO ||
            last_coll_count > NO_COLLISION_TIME_LIMIT || iter_count > RESTART_ITERATIONS) {
            restart();
            return;
        }

        const uint32_t map_area = _worldDimension.getArea();
        // Keep a copy of the previous index map
        prev_imap.copy(imap);

        // Realize accumulated external forces to each plate.
        for (uint32_t i = 0; i < num_plates; ++i) {
            plates[i]->resetSegments();

            if (erosion_period > 0 && iter_count % erosion_period == 0)
                plates[i]->erode(CONTINENTAL_BASE);

            plates[i]->move();
        }

        uint32_t oceanic_collisions = 0;
        uint32_t continental_collisions = 0;

        updateHeightAndPlateIndexMaps(oceanic_collisions, continental_collisions);

        // Update the counter of iterations since last continental collision.
        last_coll_count = (last_coll_count + 1) & -(continental_collisions == 0);

        for (uint32_t i = 0; i < num_plates; ++i) {
            for (uint32_t j = 0; j < subductions[i].size(); ++j) {
                const plateCollision& coll = subductions[i][j];

                ASSERT(i != coll.index, "when subducting: SRC == DEST!");

                // Do not apply friction to oceanic plates.
                // This is a very cheap way to emulate slab pull.
                // Just perform subduction and on our way we go!
                plates[i]->addCrustBySubduction(coll.wx, coll.wy, coll.crust, iter_count,
                                                plates[coll.index]->getVelX(),
                                                plates[coll.index]->getVelY());
            }

            subductions[i].clear();
        }

        updateCollisions();

        fill(plate_indices_found.begin(), plate_indices_found.end(), 0);
        std::vector<uint32_t> regenerated_owner;
        std::vector<float> regenerated_height;
        std::vector<float> owner_scores;
        std::vector<float> owner_height_sums;
        std::vector<uint32_t> owner_height_counts;

        if (imported_heightmap_mode) {
            regenerated_owner.assign(map_area, 0xFFFFFFFFU);
            regenerated_height.assign(map_area, 0.0f);
            owner_scores.resize(num_plates);
            owner_height_sums.resize(num_plates);
            owner_height_counts.resize(num_plates);

            for (uint32_t y = 0, i = 0; y < _worldDimension.getHeight(); ++y) {
                for (uint32_t x = 0; x < _worldDimension.getWidth(); ++x, ++i) {
                    if (imap[i] < num_plates) {
                        continue;
                    }

                    choose_regenerated_gap_fill(
                        _worldDimension, imap, prev_imap, hmap, x, y, num_plates, iter_count,
                        owner_scores, owner_height_sums, owner_height_counts, &regenerated_owner[i],
                        &regenerated_height[i]);
                }
            }
        }

        // Fill divergent boundaries with new crustal material, molten magma.
        for (uint32_t y = 0, i = 0; y < BOOL_REGENERATE_CRUST * _worldDimension.getHeight(); ++y) {
            for (uint32_t x = 0; x < _worldDimension.getWidth(); ++x, ++i) {
                if (imap[i] >= num_plates) {
                    // The owner of this new crust is that neighbour plate
                    // who was located at this point before plates moved.
                    if (imported_heightmap_mode && i < regenerated_owner.size() &&
                        regenerated_owner[i] < num_plates) {
                        imap[i] = regenerated_owner[i];
                    } else {
                        imap[i] = prev_imap[i];
                    }

                    // If this is oceanic crust then add buoyancy to it.
                    // Magma that has just crystallized into oceanic crust
                    // is more buoyant than that which has had a lot of
                    // time to cool down and become more dense.
                    if (imported_heightmap_mode && imap[i] < num_plates) {
                        const float age_noise =
                            seam_noise(x, y, iter_count, static_cast<uint32_t>(imap[i]) + 31U);
                        const uint32_t age_offset =
                            4U + static_cast<uint32_t>((age_noise + 1.0f) * 3.5f);
                        amap[i] = iter_count > age_offset ? iter_count - age_offset : 0U;
                    } else {
                        amap[i] = iter_count;
                    }
                    if (imported_heightmap_mode && i < regenerated_height.size() &&
                        regenerated_owner[i] < num_plates) {
                        hmap[i] = regenerated_height[i];
                    } else {
                        hmap[i] = OCEANIC_BASE * BUOYANCY_BONUS_X;
                    }

                    // This should probably not happen
                    if (imap[i] < num_plates) {
                        plates[imap[i]]->setCrust(x, y, OCEANIC_BASE, iter_count);
                    }

                } else if (++plate_indices_found[imap[i]] && hmap[i] <= 0) {
                    puts("Occupied point has no land mass!");
                    exit(1);
                }
            }
        }

        removeEmptyPlates();

        // delete[] indexFound;

        // Add some "virginity buoyancy" to all pixels for a visual boost! :)
        for (uint32_t i = 0; i < (BUOYANCY_BONUS_X > 0) * map_area; ++i) {
            // Calculate the inverted age of this piece of crust.
            // Force result to be minimum between inv. age and
            // max buoyancy bonus age.
            uint32_t crust_age = iter_count - amap[i];
            crust_age = MAX_BUOYANCY_AGE - crust_age;
            crust_age &= -(crust_age <= MAX_BUOYANCY_AGE);

            hmap[i] += (hmap[i] < CONTINENTAL_BASE) * BUOYANCY_BONUS_X * OCEANIC_BASE * crust_age *
                       MULINV_MAX_BUOYANCY_AGE;
        }

        ++iter_count;
    } catch (const exception& e) {
        string msg = "Problem during update: ";
        msg = msg + e.what();
        cerr << msg << endl;
        throw runtime_error(msg.c_str());
    }
}

void lithosphere::restart() {
    try {
        const uint32_t map_area = _worldDimension.getArea();

        cycle_count += max_cycles > 0; // No increment if running forever.
        if (cycle_count > max_cycles)
            return;

        // Update height map to include all recent changes.
        hmap.set_all(0);
        for (uint32_t i = 0; i < num_plates; ++i) {
            const uint32_t x0 = plates[i]->getLeftAsUint();
            const uint32_t y0 = plates[i]->getTopAsUint();
            const uint32_t x1 = x0 + plates[i]->getWidth();
            const uint32_t y1 = y0 + plates[i]->getHeight();

            const float* this_map;
            const uint32_t* this_age;
            plates[i]->getMap(&this_map, &this_age);

            // Copy first part of plate onto world map.
            for (uint32_t y = y0, j = 0; y < y1; ++y) {
                for (uint32_t x = x0; x < x1; ++x, ++j) {
                    const uint32_t x_mod = _worldDimension.xMod(x);
                    const uint32_t y_mod = _worldDimension.yMod(y);
                    const float h0 = hmap[_worldDimension.indexOf(x_mod, y_mod)];
                    const float h1 = this_map[j];
                    const uint32_t a0 = amap[_worldDimension.indexOf(x_mod, y_mod)];
                    const uint32_t a1 = this_age[j];

                    const float h_sum = h0 + h1;
                    // Avoid division by zero: if both heights are zero, use the new age
                    amap[_worldDimension.indexOf(x_mod, y_mod)] =
                        (h_sum > 0.0f) ? static_cast<uint32_t>((h0 * a0 + h1 * a1) / h_sum) : a1;
                    hmap[_worldDimension.indexOf(x_mod, y_mod)] += this_map[j];
                }
            }
        }
        // Clear plate array
        clearPlates();

        // create new plates IFF there are cycles left to run!
        // However, if max cycle count is "ETERNITY", then 0 < 0 + 1 always.
        if (cycle_count < max_cycles + !max_cycles) {
            createPlates();

            // Restore the ages of plates' points of crust!
            for (uint32_t i = 0; i < num_plates; ++i) {
                const uint32_t x0 = plates[i]->getLeftAsUint();
                const uint32_t y0 = plates[i]->getTopAsUint();
                const uint32_t x1 = x0 + plates[i]->getWidth();
                const uint32_t y1 = y0 + plates[i]->getHeight();

                const float* this_map;
                const uint32_t* this_age_const;
                uint32_t* this_age;

                plates[i]->getMap(&this_map, &this_age_const);
                this_age = const_cast<uint32_t*>(this_age_const);

                for (uint32_t y = y0, j = 0; y < y1; ++y) {
                    for (uint32_t x = x0; x < x1; ++x, ++j) {
                        const uint32_t x_mod = _worldDimension.xMod(x);
                        const uint32_t y_mod = _worldDimension.yMod(y);

                        this_age[j] = amap[_worldDimension.indexOf(x_mod, y_mod)];
                    }
                }
            }

            return;
        }

        // Add some "virginity buoyancy" to all pixels for a visual boost.
        for (uint32_t i = 0; i < (BUOYANCY_BONUS_X > 0) * map_area; ++i) {
            uint32_t crust_age = iter_count - amap[i];
            crust_age = MAX_BUOYANCY_AGE - crust_age;
            crust_age &= -(crust_age <= MAX_BUOYANCY_AGE);

            hmap[i] += (hmap[i] < CONTINENTAL_BASE) * BUOYANCY_BONUS_X * OCEANIC_BASE * crust_age *
                       MULINV_MAX_BUOYANCY_AGE;
        }
    } catch (const exception& e) {
        std::string msg = "Problem during restart: ";
        msg = msg + e.what();
        throw runtime_error(msg.c_str());
    }
}

uint32_t lithosphere::getWidth() const {
    return _worldDimension.getWidth();
}

uint32_t lithosphere::getHeight() const {
    return _worldDimension.getHeight();
}

uint32_t* lithosphere::getPlatesMap() const throw() {
    return imap.raw_data();
}

const plate* lithosphere::getPlate(uint32_t index) const {
    ASSERT(index < num_plates, "invalid plate index");
    return plates[index];
}
