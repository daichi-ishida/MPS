#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip> 
#include <string>
#include <direct.h>
#include <Eigen/Core>
#include "constants.h"
#include "Scene.h"


Scene::Scene(const double &time):
    m_time(time)
{
}


Scene::~Scene()
{
}


void Scene::makeParticleSet()
{

    switch (DIMENSION)
    {
    case 2:
        initializeParticles_for2dim();
        break;
    case 3:
        initializeParticles_for3dim();
        break;
    default:
        std::cout << "invalid value : DIMENSION" << std::endl;
        exit(EXIT_FAILURE);
    }
}


void Scene::writeData()
{
    _mkdir("output");
    writeData_inVtuFormat();
    //writeData_inProfFormat();
    ++m_file_num;
}

/* private */
void Scene::writeData_inProfFormat()
{
    std::ostringstream sout;
    sout << std::setfill('0') << std::setw(3) << std::right << m_file_num;

    std::string file_name = "output/output_" + sout.str() + ".prof";
    std::ofstream ofs;
    ofs.open(file_name);
    if (!ofs)
    {
        std::cout << "ERROR : file open error at writing data in .prof format\n" << file_name << " cannot open" << std::endl;
        exit(EXIT_FAILURE);
    }
    ofs << m_time << std::endl;
    ofs << m_particle_set.getNumberOfParticles() << std::endl;
    for (unsigned int i = 0, n = m_particle_set.getNumberOfParticles(); i < n; ++i)
    {
        ofs << m_particle_set.types[i] << " "
            << m_particle_set.positions[i].x() << " "
            << m_particle_set.positions[i].y() << " "
            << m_particle_set.positions[i].z() << " "
            << m_particle_set.velocities[i].x() << " "
            << m_particle_set.velocities[i].y() << " "
            << m_particle_set.velocities[i].z() << " "
            << m_particle_set.pressures[i] << " "
            << m_particle_set.number_densities[i] << std::endl;
    }

    ofs.close();
}

void Scene::writeData_inVtuFormat()
{
    std::ostringstream sout;
    sout << std::setfill('0') << std::setw(3) << std::right << m_file_num;

    std::string file_name = "output/particle_" + sout.str() + ".vtu";
    std::ofstream ofs;
    ofs.open(file_name);
    if (!ofs)
    {
        std::cout << "ERROR : file open error at writing data in .vtu format\n" << file_name << " cannot open" << std::endl;
        exit(EXIT_FAILURE);
    }

    /* header */
    ofs << "<?xml version='1.0' encoding='UTF-8'?>" << std::endl;
    ofs << "<VTKFile xmlns='VTK' byte_order='LittleEndian' version='0.1' type='UnstructuredGrid'>" << std::endl;
    ofs << "<UnstructuredGrid>" << std::endl;

    ofs << "<Piece NumberOfCells='" << m_particle_set.getNumberOfParticles() << "' NumberOfPoints='" << m_particle_set.getNumberOfParticles() << "'>" << std::endl;

    /* point position */
    ofs << "<Points>" << std::endl;
    ofs << "<DataArray NumberOfComponents='3' type='Float32' Name='Position' format='ascii'>" << std::endl;
    for (unsigned int i = 0, n = m_particle_set.getNumberOfParticles(); i < n; ++i)
    {
        ofs << m_particle_set.positions[i].x() << " "
            << m_particle_set.positions[i].y() << " "
            << m_particle_set.positions[i].z() << std::endl;
    }
    ofs << "</DataArray>" << std::endl;
    ofs << "</Points>" << std::endl;

    /* point data */
    ofs << "<PointData>" << std::endl;

    // type
    ofs << "<DataArray NumberOfComponents='1' type='Int32' Name='ParticleType' format='ascii'>" << std::endl;
    for (unsigned int i = 0, n = m_particle_set.getNumberOfParticles(); i < n; ++i)
    {
        ofs << m_particle_set.types[i] << std::endl;
    }
    ofs << "</DataArray>" << std::endl;

    // velocity ( speed )
    ofs << "<DataArray NumberOfComponents='1' type='Float32' Name='Velocity' format='ascii'>" << std::endl;
    for (unsigned int i = 0, n = m_particle_set.getNumberOfParticles(); i < n; ++i)
    {
        double speed = m_particle_set.velocities[i].norm();
        ofs << speed << std::endl;
    }
    ofs << "</DataArray>" << std::endl;

    // pressures[i]
    ofs << "<DataArray NumberOfComponents='1' type='Float32' Name='Pressure' format='ascii'>" << std::endl;
    for (unsigned int i = 0, n = m_particle_set.getNumberOfParticles(); i < n; ++i)
    {
        ofs << m_particle_set.pressures[i] << std::endl;
    }
    ofs << "</DataArray>" << std::endl;

    // pressures[i] gradient
    ofs << "<DataArray NumberOfComponents='1' type='Float32' Name='Pressure Gradient' format='ascii'>" << std::endl;
    for (unsigned int i = 0, n = m_particle_set.getNumberOfParticles(); i < n; ++i)
    {
        double grad_magnitude = m_particle_set.pressure_gradients[i].norm();
        ofs << grad_magnitude << std::endl;
    }
    ofs << "</DataArray>" << std::endl;

    ofs << "</PointData>" << std::endl;

    ofs << "<Cells>" << std::endl;
    ofs << "<DataArray type='Int32' Name='connectivity' format='ascii'>" << std::endl;

    for (unsigned int i = 0, n = m_particle_set.getNumberOfParticles(); i < n; ++i)
    {
        ofs << i << std::endl;
    }
    ofs << "</DataArray>" << std::endl;
    ofs << "<DataArray type='Int32' Name='offsets' format='ascii'>" << std::endl;

    for (unsigned int i = 0, n = m_particle_set.getNumberOfParticles(); i < n; ++i)
    {
        ofs << i + 1 << std::endl;
    }
    ofs << "</DataArray>" << std::endl;
    ofs << "<DataArray type='UInt8' Name='types' format='ascii'>" << std::endl;

    for (unsigned int i = 0, n = m_particle_set.getNumberOfParticles(); i < n; ++i)
    {
        ofs << 1 << std::endl;
    }
    ofs << "</DataArray>" << std::endl;
    ofs << "</Cells>" << std::endl;

    ofs << "</Piece>" << std::endl;

    ofs << "</UnstructuredGrid>" << std::endl;
    ofs << "</VTKFile>" << std::endl;

    ofs.close();
}

void Scene::initializeParticles_for2dim()
{
    int	   type;
    double x, y;
    bool   canGenerate;

    /*
    /* x_max,y_max : total particle in each direction
    /* x : y = 1.0 : 0.6
    /* +5 and -4 is for wall and dummy particles
    */

    int x_max = (int)(1.0 / PARTICLE_DISTANCE) + 5;
    int y_max = (int)(0.6 / PARTICLE_DISTANCE) + 5;
    for (int i = -4; i < x_max; ++i)
    {
        for (int j = -4; j < y_max; ++j)
        {
            x = PARTICLE_DISTANCE * i;
            y = PARTICLE_DISTANCE * j;
            canGenerate = false;

            /* dummy wall region */
            if (((x > -4.0 * PARTICLE_DISTANCE + EPS) &&
                (x <= 1.00 + 4.0 * PARTICLE_DISTANCE + EPS)) &&
                ((y > 0.0 - 4.0 * PARTICLE_DISTANCE + EPS) &&
                (y <= 0.6 + EPS)))
            {
                type = DUMMY_WALL;
                canGenerate = true;
            }

            /* wall region */
            if (((x > -2.0 * PARTICLE_DISTANCE + EPS) &&
                (x <= 1.00 + 2.0 * PARTICLE_DISTANCE + EPS)) &&
                ((y > 0.0 - 2.0 * PARTICLE_DISTANCE + EPS) &&
                (y <= 0.6 + EPS)))
            {
                type = WALL;
                canGenerate = true;
            }

            /* wall region */
            if (((x > -4.0 * PARTICLE_DISTANCE + EPS) &&
                (x <= 1.00 + 4.0 * PARTICLE_DISTANCE + EPS)) &&
                ((y > 0.6 - 2.0 * PARTICLE_DISTANCE + EPS) &&
                (y <= 0.6 + EPS)))
            {
                type = WALL;
                canGenerate = true;
            }

            /* empty region */
            if (((x > 0.0 + EPS) &&
                (x <= 1.00 + EPS)) &&
                (y > 0.0 + EPS))
            {
                canGenerate = false;
            }

            /* fluid region */
            if (((x > 0.0 + EPS) &&
                (x <= 0.25 + EPS)) &&
                ((y > 0.0 + EPS) &&
                (y <= 0.50 + EPS)))
            {
                type = FLUID;
                canGenerate = true;
            }

            if (canGenerate)
            {
                m_particle_set.initSet(type, x, y, 0.0, 0.0);
            }
        }
    }
}

void Scene::initializeParticles_for3dim()
{
    int    type;
    double x, y, z;
    bool   canGenerate;

    /*
    /* x_max,y_max,z_max : total particle in each direction
    /* x : y : z = 1.0 : 0.6 : 0.3
    /* +5 and -4 is for wall and dummy particles
    */

    int x_max = (int)(1.0 / PARTICLE_DISTANCE) + 5;
    int y_max = (int)(0.6 / PARTICLE_DISTANCE) + 5;
    int z_max = (int)(0.3 / PARTICLE_DISTANCE) + 5;

    for (int i = -4; i < x_max; ++i)
    {
        for (int j = -4; j < y_max; ++j)
        {
            for (int k = -4; k < z_max; ++k)
            {
                x = PARTICLE_DISTANCE * i;
                y = PARTICLE_DISTANCE * j;
                z = PARTICLE_DISTANCE * k;
                canGenerate = false;

                /* dummy wall region */
                if ((((x > -4.0 * PARTICLE_DISTANCE + EPS) &&
                    (x <= 1.00 + 4.0 * PARTICLE_DISTANCE + EPS)) &&
                    ((y > 0.0 - 4.0 * PARTICLE_DISTANCE + EPS) &&
                    (y <= 0.6 + EPS))) &&
                        ((z > 0.0 - 4.0 * PARTICLE_DISTANCE + EPS) &&
                    (z <= 0.3 + 4.0 * PARTICLE_DISTANCE + EPS)))
                {
                    type = DUMMY_WALL;
                    canGenerate = true;
                }

                /* wall region */
                if ((((x > -2.0 * PARTICLE_DISTANCE + EPS) &&
                    (x <= 1.00 + 2.0 * PARTICLE_DISTANCE + EPS)) &&
                    ((y > 0.0 - 2.0 * PARTICLE_DISTANCE + EPS) &&
                    (y <= 0.6 + EPS))) &&
                        ((z > 0.0 - 2.0 * PARTICLE_DISTANCE + EPS) &&
                    (z <= 0.3 + 2.0 * PARTICLE_DISTANCE + EPS)))
                {
                    type = WALL;
                    canGenerate = true;
                }

                /* wall region */
                if ((((x > -4.0 * PARTICLE_DISTANCE + EPS) &&
                    (x <= 1.00 + 4.0 * PARTICLE_DISTANCE + EPS)) &&
                    ((y > 0.6 - 2.0 * PARTICLE_DISTANCE + EPS) &&
                    (y <= 0.6 + EPS))) &&
                        ((z > 0.0 - 4.0 * PARTICLE_DISTANCE + EPS) &&
                    (z <= 0.3 + 4.0 * PARTICLE_DISTANCE + EPS)))
                {
                    type = WALL;
                    canGenerate = true;
                }

                /* empty region */
                if ((((x > 0.0 + EPS) &&
                    (x <= 1.00 + EPS)) &&
                    (y > 0.0 + EPS)) &&
                    ((z > 0.0 + EPS) &&
                    (z <= 0.3 + EPS)))
                {
                    canGenerate = false;
                }

                /* fluid region */
                if ((((x > 0.0 + EPS) &&
                    (x <= 0.25 + EPS)) &&
                    ((y > 0.0 + EPS) &&
                    (y < 0.5 + EPS))) &&
                        ((z > 0.0 + EPS) &&
                    (z <= 0.3 + EPS)))
                {
                    type = FLUID;
                    canGenerate = true;
                }

                if (canGenerate)
                {
                    m_particle_set.initSet(type, x, y, z, 0.0, 0.0, 0.0);
                }
            }
        }
    }
}
