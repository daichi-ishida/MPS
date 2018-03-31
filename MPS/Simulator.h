#pragma once
#include "ParticleSet.h"

class Simulator
{
public:
    Simulator(ParticleSet *particle_set);
    ~Simulator();

    void update();
    const int getNumberOfParticles() { return m_particle_set->getNumberOfParticles(); };

private:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    void  setNZeroAndLambda();
    double weight(double distance, double re);

    void moveToTentativePos();
    void collision();
    void calNumberDensity();
    void setBoundaryCondition();
    void calPressure();
    void calPressureGradient();
    void moveParticleUsingPressureGradient();

    ParticleSet* m_particle_set;

    const double m_re_forNumberDensity;
    const double m_re_forGradient;
    const double m_re_forLaplacian;

    const double m_fluid_density;
    const double m_collision_distance;

    double m_n0_forNumberDensity;
    double m_n0_forGradient;

    double m_n0_lambda;

};

