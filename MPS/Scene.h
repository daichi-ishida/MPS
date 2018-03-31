#pragma once
#include "ParticleSet.h"

class Scene
{
public:
    Scene(const double &time);
    ~Scene();

    void makeParticleSet();
    void writeData();

    ParticleSet* getParticleSetPtr() { return &m_particle_set; };

    ParticleSet  m_particle_set;


private:
    void initializeParticles_for2dim();
    void initializeParticles_for3dim();
    void writeData_inProfFormat();
    void writeData_inVtuFormat();

    int          m_file_num;
    const double &m_time;
};

