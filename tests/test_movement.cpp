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

#include "movement.hpp"
#include "mass.hpp"
#include "plate.hpp"
#include "gtest/gtest.h"

using namespace Platec;

TEST(Movement, Constructor)
{
    SimpleRandom sr(123);
    WorldDimension wd(5, 4);
    Movement mov(sr, wd);
    EXPECT_FLOAT_EQ(0.99992257f, mov.velX());
    EXPECT_FLOAT_EQ(0.01244594f, mov.velY());
    EXPECT_FLOAT_EQ(1.0f, mov.getVelocity());
}

TEST(Movement, ApplyFriction)
{
    SimpleRandom sr(456);
    WorldDimension wd(50, 40);
    Movement mov(sr, wd);

    EXPECT_FLOAT_EQ(0.9989379f, mov.velX());
    EXPECT_FLOAT_EQ(0.046077024f, mov.velY());

    mov.applyFriction(2.2f, 10.5f);
    EXPECT_FLOAT_EQ(0.9989379f, mov.velX());
    EXPECT_FLOAT_EQ(0.046077024f, mov.velY());
    EXPECT_FLOAT_EQ(0.58095241f, mov.getVelocity());

    mov.applyFriction(7.2f, 0.0f);
    EXPECT_FLOAT_EQ(0.0f, mov.getVelocity());
}

TEST(Movement, ApplyFrictionWithNullMass)
{
    SimpleRandom sr(456);
    WorldDimension wd(50, 40);
    Movement mov(sr, wd);

    mov.applyFriction(7.2f, 0.0f);
    EXPECT_FLOAT_EQ(0.0f, mov.getVelocity());
}

TEST(Movement, Move)
{
    SimpleRandom sr(789890);
    WorldDimension wd(500, 400);
    Movement mov(sr, wd);

    EXPECT_FLOAT_EQ(-0.29389676f, mov.velX());
    EXPECT_FLOAT_EQ(-0.95583719f, mov.velY());
    EXPECT_FLOAT_EQ(1.0f, mov.getVelocity());

    mov.move();
    EXPECT_FLOAT_EQ(-0.28745356f, mov.velX());
    EXPECT_FLOAT_EQ(-0.95779467f, mov.velY());
    EXPECT_FLOAT_EQ(1.0f, mov.getVelocity());
}

TEST(Movement, VelocityOnXNoParams)
{
    SimpleRandom sr(789890);
    WorldDimension wd(500, 400);
    Movement mov(sr, wd);

    EXPECT_FLOAT_EQ(-0.29389676f, mov.velX());
    EXPECT_FLOAT_EQ(1.0f, mov.getVelocity());
    EXPECT_FLOAT_EQ(-0.29389676f, mov.velocityOnX());
}

TEST(Movement, VelocityOnYNoParams)
{
    SimpleRandom sr(789890);
    WorldDimension wd(500, 400);
    Movement mov(sr, wd);

    EXPECT_FLOAT_EQ(-0.95583719f, mov.velY());
    EXPECT_FLOAT_EQ(1.0f, mov.getVelocity());
    EXPECT_FLOAT_EQ(-0.95583719f, mov.velocityOnY());
}

TEST(Movement, VelocityOnXOneParam)
{
    SimpleRandom sr(789890);
    WorldDimension wd(500, 400);
    Movement mov(sr, wd);

    EXPECT_FLOAT_EQ(-0.29389676f, mov.velX());
    EXPECT_FLOAT_EQ(1.0f, mov.getVelocity());
    EXPECT_FLOAT_EQ(-2.9389676f, mov.velocityOnX(10.0));
}

TEST(Movement, VelocityOnYOneParam)
{
    SimpleRandom sr(789890);
    WorldDimension wd(500, 400);
    Movement mov(sr, wd);

    EXPECT_FLOAT_EQ(-0.95583719f, mov.velY());
    EXPECT_FLOAT_EQ(1.0f, mov.getVelocity());
    EXPECT_FLOAT_EQ(-9.5583719f, mov.velocityOnY(10.0));
}

TEST(Movement, Dot)
{
    SimpleRandom sr(789890);
    WorldDimension wd(500, 400);
    Movement mov(sr, wd);

    EXPECT_FLOAT_EQ(-0.29389676f, mov.velX());
    EXPECT_FLOAT_EQ(-0.95583719f, mov.velY());
    EXPECT_FLOAT_EQ(-3.45530509f, mov.dot(2.0, 3.0));
}

class MockPlate : public IPlate {
public:
    MockPlate(const FloatVector& velocityUnitVector, float mass, const FloatPoint& massCenter)
        : MockPlate(velocityUnitVector, mass, massCenter, massCenter)
    { }

    MockPlate(const FloatVector& velocityUnitVector, float mass, const FloatPoint& massCenter,
              const FloatPoint& worldMassCenter)
        : _velocityUnitVector(velocityUnitVector),
          _decImpulseDelta(nullptr),
          _mass(mass),
          _massCenter(massCenter),
          _worldMassCenter(worldMassCenter)
    { }

    ~MockPlate() {
        if (_decImpulseDelta) delete _decImpulseDelta;
    }

    FloatVector velocityUnitVector() const override {
        return _velocityUnitVector;
    }

    void decImpulse(const FloatVector& delta) override {
        _decImpulseDelta = new FloatVector(delta);
    }

    FloatVector decImpulseDelta() {
        if (_decImpulseDelta == nullptr) throw runtime_error("(MockPlate::decImpulseDelta) Data not ready");
        return *_decImpulseDelta;
    }

    float getMass() const override {
        return _mass;
    }

    FloatPoint massCenter() const override {
        return _massCenter;
    }

    FloatPoint worldMassCenter() const override {
        return _worldMassCenter;
    }

private:
    FloatVector _velocityUnitVector;
    FloatVector* _decImpulseDelta;
    float _mass;
    FloatPoint _massCenter;
    FloatPoint _worldMassCenter;
};

TEST(Movement, Collide)
{
    SimpleRandom sr(789890);
    WorldDimension wd(500, 400);
    Movement mov(sr, wd);

    EXPECT_FLOAT_EQ(-0.29389676f, mov.velX());
    EXPECT_FLOAT_EQ(-0.95583719f, mov.velY());
    EXPECT_FLOAT_EQ(1.0f, mov.getVelocity());

    MockPlate thisPlate(FloatVector(mov.velX(), mov.velY()), 100.0f, FloatPoint(70.0f, 90.0f));
    FloatVector otherPlateVelocityUnitVector(0.0f, -1.0f);
    float otherPlateMass = 10000.0f;
    FloatPoint otherPlateMassCenter(20.0f, 120.0f);
    MockPlate otherPlate(otherPlateVelocityUnitVector, otherPlateMass, otherPlateMassCenter);
    mov.collide(thisPlate, otherPlate, 456.2f);

    EXPECT_FLOAT_EQ(0.0098300017f, otherPlate.decImpulseDelta().x());
    EXPECT_FLOAT_EQ(-0.0058980011f, otherPlate.decImpulseDelta().y());
}

TEST(Movement, CollideWrapsWorldMassCenter)
{
    SimpleRandom sr(123);
    WorldDimension wd(100, 60);
    Movement mov(sr, wd);

    MockPlate thisPlate(FloatVector(mov.velX(), mov.velY()), 100.0f, FloatPoint(8.0f, 12.0f),
                        FloatPoint(98.0f, 20.0f));
    MockPlate otherPlate(FloatVector(-1.0f, 0.0f), 1000.0f, FloatPoint(4.0f, 9.0f),
                         FloatPoint(2.0f, 20.0f));

    mov.collide(thisPlate, otherPlate, 250.0f);

    EXPECT_LT(otherPlate.decImpulseDelta().x(), -0.001f);
    EXPECT_LT(std::fabs(otherPlate.decImpulseDelta().y()), std::fabs(otherPlate.decImpulseDelta().x()));
}
