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

#include "platecapi.hpp"
#include "lithosphere.hpp"
#include "plate.hpp"
#include <stdio.h>
#include <stdlib.h>

#include <vector>

class platec_api_list_elem {
  public:
    platec_api_list_elem(uint32_t _id, lithosphere* _data) : data(_data), id(_id) {}

    lithosphere* data;
    uint32_t id;
};

extern lithosphere* platec_api_get_lithosphere(uint32_t);

static std::vector<platec_api_list_elem> lithospheres;
static uint32_t last_id = 1;

void* platec_api_create(long seed, uint32_t width, uint32_t height, float sea_level,
                        uint32_t erosion_period, float folding_ratio, uint32_t aggr_overlap_abs,
                        float aggr_overlap_rel, uint32_t cycle_count, uint32_t num_plates,
                        float erosion_strength, float crust_rotation_strength,
                        float rotation_strength) {
    /* Miten nykyisen opengl-mainin koodit refaktoroidaan tänne?
     *    parametrien tarkistus, kommentit eli dokumentointi, muuta? */

    lithosphere* litho =
        new lithosphere(seed, width, height, sea_level, erosion_period, folding_ratio,
                        aggr_overlap_abs, aggr_overlap_rel, cycle_count, num_plates,
                        erosion_strength, crust_rotation_strength, rotation_strength);

    platec_api_list_elem elem(++last_id, litho);
    lithospheres.push_back(elem);

    return litho;
}

void platec_api_destroy(void* litho) {
    for (uint32_t i = 0; i < lithospheres.size(); ++i)
        if (lithospheres[i].data == litho) {
            lithospheres.erase(lithospheres.begin() + i);
            break;
        }
}

const uint32_t* platec_api_get_agemap(uint32_t id) {
    lithosphere* litho = platec_api_get_lithosphere(id);
    if (!litho)
        return nullptr;

    return litho->getAgeMap();
}

float* platec_api_get_heightmap(void* pointer) {
    lithosphere* litho = static_cast<lithosphere*>(pointer);
    float* res = litho->getTopography();
    return res;
}

uint32_t* platec_api_get_platesmap(void* pointer) {
    lithosphere* litho = static_cast<lithosphere*>(pointer);
    uint32_t* res = litho->getPlatesMap();
    return res;
}

lithosphere* platec_api_get_lithosphere(uint32_t id) {
    for (uint32_t i = 0; i < lithospheres.size(); ++i)
        if (lithospheres[i].id == id)
            return lithospheres[i].data;

    return nullptr;
}

uint32_t platec_api_is_finished(void* pointer) {
    lithosphere* litho = static_cast<lithosphere*>(pointer);
    if (litho->isFinished()) {
        return 1;
    } else {
        return 0;
    }
}

void platec_api_step(void* pointer) {
    lithosphere* litho = static_cast<lithosphere*>(pointer);
    litho->update();
}

void platec_api_load_heightmap(void* pointer, const float* normalized_heightmap, float sea_level) {
    lithosphere* litho = static_cast<lithosphere*>(pointer);
    litho->importNormalizedHeightMap(normalized_heightmap, sea_level);
}

uint32_t platec_api_get_plate_count(void* pointer) {
    lithosphere* litho = static_cast<lithosphere*>(pointer);
    return litho->getPlateCount();
}

uint32_t lithosphere_getMapWidth(void* object) {
    return static_cast<lithosphere*>(object)->getWidth();
}

uint32_t lithosphere_getMapHeight(void* object) {
    return static_cast<lithosphere*>(object)->getHeight();
}

float platec_api_velocity_unity_vector_x(void* pointer, uint32_t plate_index) {
    lithosphere* litho = static_cast<lithosphere*>(pointer);
    return litho->getPlate(plate_index)->velocityUnitVector().x();
}

float platec_api_velocity_unity_vector_y(void* pointer, uint32_t plate_index) {
    lithosphere* litho = static_cast<lithosphere*>(pointer);
    return litho->getPlate(plate_index)->velocityUnitVector().y();
}

float platec_api_velocity_vector_x(void* pointer, uint32_t plate_index) {
    lithosphere* litho = static_cast<lithosphere*>(pointer);
    return litho->getPlate(plate_index)->linearVelocityVector().x();
}

float platec_api_velocity_vector_y(void* pointer, uint32_t plate_index) {
    lithosphere* litho = static_cast<lithosphere*>(pointer);
    return litho->getPlate(plate_index)->linearVelocityVector().y();
}

float platec_api_angular_velocity(void* pointer, uint32_t plate_index) {
    lithosphere* litho = static_cast<lithosphere*>(pointer);
    return litho->getPlate(plate_index)->getAngularVelocity();
}

float platec_api_mass_center_x(void* pointer, uint32_t plate_index) {
    lithosphere* litho = static_cast<lithosphere*>(pointer);
    return litho->getPlate(plate_index)->worldMassCenter().getX();
}

float platec_api_mass_center_y(void* pointer, uint32_t plate_index) {
    lithosphere* litho = static_cast<lithosphere*>(pointer);
    return litho->getPlate(plate_index)->worldMassCenter().getY();
}
