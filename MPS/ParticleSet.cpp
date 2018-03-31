#include "constants.h"
#include "ParticleSet.h"

ParticleSet::ParticleSet() :
    m_num(0)
{
}


ParticleSet::~ParticleSet()
{
}

void ParticleSet::initSet(int type, double x, double y, double v_x, double v_y)
{
    ++m_num;
    types.emplace_back(type);
    boundary_conditions.emplace_back(0);
    positions.emplace_back(x, y, 0);
    velocities.emplace_back(v_x, v_y, 0);
    accelerations.emplace_back(0, 0, 0);
    pressure_gradients.emplace_back(0, 0, 0);
    pressures.emplace_back(0);
    min_pressures.emplace_back(0);
    number_densities.emplace_back(0);
    flag_for_checking_boundary_conditions.emplace_back(0);
}

void ParticleSet::initSet(int type, double x, double y, double z, double v_x, double v_y, double v_z)
{
    ++m_num;
    types.emplace_back(type);

    if (type == GHOST || type == DUMMY_WALL)
    {
        boundary_conditions.emplace_back(GHOST_OR_DUMMY);
    }
    positions.emplace_back(x, y, z);
    velocities.emplace_back(v_x, v_y, v_z);
    accelerations.emplace_back(0, 0, 0);
    pressure_gradients.emplace_back(0, 0, 0);
    pressures.emplace_back(0);
    min_pressures.emplace_back(0);
    number_densities.emplace_back(0);
    flag_for_checking_boundary_conditions.emplace_back(0);
}
