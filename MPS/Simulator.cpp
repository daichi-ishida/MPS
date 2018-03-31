#include <iostream>
#include <omp.h>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include "Simulator.h"
#include "constants.h"


Simulator::Simulator(ParticleSet *particle_set) :
    m_particle_set(particle_set),
    m_re_forNumberDensity(RADIUS_FOR_NUMBER_DENSITY),
    m_re_forGradient(RADIUS_FOR_GRADIENT),
    m_re_forLaplacian(RADIUS_FOR_LAPLACIAN),
    m_fluid_density(FLUID_DENSITY),
    m_collision_distance(COLLISION_DISTANCE)
{
    setNZeroAndLambda();
}

Simulator::~Simulator()
{
    delete m_particle_set;
}

void Simulator::update()
{
    for (int i = 0, n = m_particle_set->getNumberOfParticles(); i < n; ++i)
    {
        m_particle_set->accelerations[i].setZero();
        m_particle_set->pressures[i] = 0;
        m_particle_set->pressure_gradients[i].setZero();
    }

    moveToTentativePos();
    collision();
    calNumberDensity();
    setBoundaryCondition();
    calPressure();
    calPressureGradient();
    moveParticleUsingPressureGradient();
}


/* private */
void Simulator::setNZeroAndLambda()
{
    Eigen::Vector3d ri(0, 0, 0), rj(0, 0, 0);
    Eigen::Vector3i r_min;
    Eigen::Vector3i r_max;
    double distance;

    r_min.setConstant(-4);
    r_max.setConstant(5);

    if (DIMENSION == 2) {
        r_min.z() = 0;
        r_max.z() = 1;
    }

    for (int x = r_min.x(); x < r_max.x(); ++x)
    {
        for (int y = r_min.y(); y < r_max.y(); ++y)
        {
            for (int z = r_min.z(); z < r_max.z(); ++z)
            {
                if (((x == 0) && (y == 0)) && (z == 0))continue;
                rj = Eigen::Vector3d(x, y, z) * PARTICLE_DISTANCE;
                distance = (rj - ri).norm();

                m_n0_forNumberDensity += weight(distance, RADIUS_FOR_NUMBER_DENSITY);
                m_n0_forGradient += weight(distance, RADIUS_FOR_GRADIENT);
                m_n0_lambda += distance * distance * weight(distance, RADIUS_FOR_LAPLACIAN);
            }
        }
    }
}

double Simulator::weight(double distance, double re)
{
    //‰e‹¿”ÍˆÍ“à‚Ì‚Ýd‚Ý‚Ã‚¯
    return (distance < re) ? re / distance - 1.0 : 0.0;
}

void Simulator::moveToTentativePos()
{
    Eigen::Vector3d viscosity_term(0, 0, 0);
    Eigen::Vector3d r_ij(0, 0, 0), v_ij(0, 0, 0);
    double a = KINEMATIC_VISCOSITY * 2.0 * DIMENSION / m_n0_lambda;

    for (int i = 0, n = m_particle_set->getNumberOfParticles(); i < n; ++i)
    {
        if (m_particle_set->types[i] != FLUID) continue;

        m_particle_set->accelerations[i] << GRAVITY_X, GRAVITY_Y, GRAVITY_Z;

        for (int j = 0; j < n; ++j)
        {
            if (j == i || m_particle_set->types[j] == GHOST) continue;
            r_ij = m_particle_set->positions[j] - m_particle_set->positions[i];
            v_ij = m_particle_set->velocities[j] - m_particle_set->velocities[i];
            double distance = r_ij.norm();
            viscosity_term += v_ij * weight(distance, m_re_forLaplacian);
        }

        viscosity_term *= a;
        m_particle_set->accelerations[i] += viscosity_term;
        m_particle_set->velocities[i] += m_particle_set->accelerations[i] * DELTA_TIME;
        m_particle_set->positions[i] += m_particle_set->velocities[i] * DELTA_TIME;
        m_particle_set->accelerations[i].setZero();
        viscosity_term.setZero();
    }
}

void Simulator::collision()
{
    Eigen::Vector3d v_ij_after_collision(0, 0, 0);
    Eigen::Vector3d n_ij(0, 0, 0);
    Eigen::Vector3d r_ij(0, 0, 0), v_ij(0, 0, 0);
    double e = COEFFICIENT_OF_RESTITUTION;

    for (int i = 0, n = m_particle_set->getNumberOfParticles(); i < n; ++i)
    {
        if (m_particle_set->types[i] != FLUID) continue;

        for (int j = 0; j < n; ++j)
        {
            if (j == i || m_particle_set->types[j] == GHOST) continue;
            r_ij = m_particle_set->positions[j] - m_particle_set->positions[i];
            v_ij = m_particle_set->velocities[j] - m_particle_set->velocities[i];
            double distance = r_ij.norm();
            if (distance < m_collision_distance)
            {
                n_ij = r_ij.normalized();
                double speed_along_nij = v_ij.dot(n_ij);

                // if speed along_nij is negative,
                // i and j get closer than now (it causes unstable weight calculation)
                if (speed_along_nij < 0)
                {
                    v_ij_after_collision = (1 + e) * speed_along_nij * n_ij / 2;
                    m_particle_set->velocities[i] += v_ij_after_collision;
                    m_particle_set->positions[i] += v_ij_after_collision * DELTA_TIME;
                }
            }
        }
    }
}

void Simulator::calNumberDensity()
{
    Eigen::Vector3d r_ij(0, 0, 0);

    for (int i = 0, n = m_particle_set->getNumberOfParticles(); i < n; ++i)
    {
        m_particle_set->number_densities[i] = 0;
        if (m_particle_set->types[i] == GHOST) continue;

        for (int j = 0; j < n; ++j)
        {
            if (j == i || m_particle_set->types[j] == GHOST) continue;
            r_ij = m_particle_set->positions[j] - m_particle_set->positions[i];
            double distance = r_ij.norm();

            m_particle_set->number_densities[i] += weight(distance, m_re_forNumberDensity);
        }
    }
}

void Simulator::setBoundaryCondition()
{
    double beta = THRESHOLD_RATIO_OF_NUMBER_DENSITY;

    for (int i = 0, n = m_particle_set->getNumberOfParticles(); i < n; ++i)
    {
        if (m_particle_set->types[i] == GHOST || m_particle_set->types[i] == DUMMY_WALL) 
        {
            m_particle_set->boundary_conditions[i] = GHOST_OR_DUMMY;
        }

        else if (m_particle_set->number_densities[i] < beta * m_n0_forNumberDensity)
        {
            m_particle_set->boundary_conditions[i] = SURFACE_PARTICLE;
        }
        else
        {
            m_particle_set->boundary_conditions[i] = INNER_PARTICLE;
        }
    }
}

void Simulator::calPressure()
{
    double gamma = RELAXATION_COEFFICIENT_FOR_PRESSURE;
    //double gamma = 1;
    double a = 2.0 * DIMENSION / m_n0_lambda;
    Eigen::Vector3d r_ij(0, 0, 0);

    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
    Eigen::SparseMatrix<double, Eigen::RowMajor> SparseA(m_particle_set->getNumberOfParticles(), m_particle_set->getNumberOfParticles());

    Eigen::MatrixXd A(m_particle_set->getNumberOfParticles(), m_particle_set->getNumberOfParticles());
    Eigen::VectorXd x(m_particle_set->getNumberOfParticles());
    Eigen::VectorXd b(m_particle_set->getNumberOfParticles());
    A.setZero();
    x.setZero();
    b.setZero();

    for (int i = 0, n = m_particle_set->getNumberOfParticles(); i < n; ++i)
    {
        if (m_particle_set->boundary_conditions[i] == INNER_PARTICLE)
        {
            /*
            /* Set source term (right side of Ax=b)
            */

            b(i) = gamma * (1.0 / (DELTA_TIME * DELTA_TIME))*((m_particle_set->number_densities[i] - m_n0_forNumberDensity) / m_n0_forNumberDensity);

            /*
            /* Set matrix A (left side of Ax=b)
            */

            for (int j = 0; j < n; ++j)
            {
                if (m_particle_set->boundary_conditions[j] == GHOST_OR_DUMMY || j == i) continue;
                r_ij = m_particle_set->positions[j] - m_particle_set->positions[i];
                double distance = r_ij.norm();
                if (distance < m_re_forLaplacian)
                {
                    double coefficient_ij = a * weight(distance, m_re_forLaplacian) / m_fluid_density;
                    A.coeffRef(i, j) = -coefficient_ij;
                    A.coeffRef(i, i) += coefficient_ij;
                }
            }
            A.coeffRef(i, i) += COMPRESSIBILITY / (DELTA_TIME * DELTA_TIME);
        }

        /*
        /* Exceptional processing for boundary condition
        /*  - If tere is no Dirichlet boundary condition on the fluid,
        /*	  increase the diagonal terms of the matrix for an exception. This allows us to solve the matrix without Dirichlet boundary conditions.
        */

        switch (m_particle_set->boundary_conditions[i])
        {
            case GHOST_OR_DUMMY:
                m_particle_set->flag_for_checking_boundary_conditions[i] = GHOST_OR_DUMMY;
                break;

            case SURFACE_PARTICLE:
                m_particle_set->flag_for_checking_boundary_conditions[i] = DIRICHLET_BOUNDARY_IS_CONNECTED;
                break;

            default:
                m_particle_set->flag_for_checking_boundary_conditions[i] = DIRICHLET_BOUNDARY_IS_NOT_CONNECTED;
                break;
        }
    }

    /*
    /* Check boundary condition
    /*  - This procedure is repeated until the all fluid or wall particles (which have Dirhchlet boundary condition in the particle group)
    /*    are in the state of "DIRICHLET_BOUNDARY_IS_CHECKED".
    */

    int count;
    do
    {
        count = 0;
        for (int i = 0, n = m_particle_set->getNumberOfParticles(); i < n; ++i)
        {
            if (m_particle_set->flag_for_checking_boundary_conditions[i] == DIRICHLET_BOUNDARY_IS_CONNECTED)
            {
                for (int j = 0; j < n; ++j)
                {
                    if (m_particle_set->boundary_conditions[j] == GHOST_OR_DUMMY || j == i) continue;
                    if (m_particle_set->flag_for_checking_boundary_conditions[j] == DIRICHLET_BOUNDARY_IS_NOT_CONNECTED)
                    {
                        r_ij = m_particle_set->positions[j] - m_particle_set->positions[i];
                        double distance = r_ij.norm();
                        if (distance < m_re_forLaplacian)
                        {
                            m_particle_set->flag_for_checking_boundary_conditions[j] = DIRICHLET_BOUNDARY_IS_CONNECTED;
                        }
                    }
                }
                m_particle_set->flag_for_checking_boundary_conditions[i] = DIRICHLET_BOUNDARY_IS_CHECKED;
                ++count;
            }
        }
    } while (count != 0);

    // Increase diagonal term
    for (int i = 0, n = m_particle_set->getNumberOfParticles(); i < n; ++i)
    {
        if (m_particle_set->flag_for_checking_boundary_conditions[i] == DIRICHLET_BOUNDARY_IS_NOT_CONNECTED)
        {
            //std::cout << "WARNING: There is no dirichlet boundary condition for " << i << "-th particle." << std::endl;
            A.coeffRef(i, i) *= 2.0;
        }
    }

    /*
    /* Solve simultanious equations by bi conjugate gradient
    */

    SparseA = A.sparseView();
    solver.compute(SparseA);
    x = solver.solve(b);

    for (int i = 0, n = m_particle_set->getNumberOfParticles(); i < n; ++i)
    {
        if (x(i) < 0.0) x(i) = 0.0;
        m_particle_set->pressures[i] = x(i);
    }

}

void Simulator::calPressureGradient()
{
    Eigen::Vector3d r_ij(0, 0, 0);

    for (int i = 0, n = m_particle_set->getNumberOfParticles(); i < n; ++i)
    {
        if (m_particle_set->types[i] == FLUID || m_particle_set->types[i] == WALL)
        {
            m_particle_set->min_pressures[i] = m_particle_set->pressures[i];
            for (int j = 0; j < n; ++j)
            {
                if (m_particle_set->boundary_conditions[j] == GHOST_OR_DUMMY || j == i) continue;
                r_ij = m_particle_set->positions[j] - m_particle_set->positions[i];
                double distance = r_ij.norm();
                if (distance < m_re_forGradient)
                {
                    if (m_particle_set->min_pressures[i] > m_particle_set->pressures[j])
                    {
                        m_particle_set->min_pressures[i] = m_particle_set->pressures[j];
                    }
                    if (m_particle_set->types[i] == FLUID)
                    {
                        m_particle_set->pressure_gradients[i] += r_ij * (m_particle_set->pressures[j] - m_particle_set->min_pressures[i]) * weight(distance, m_re_forGradient) / (distance * distance);
                    }
                }
            }
            m_particle_set->pressure_gradients[i] *= DIMENSION / m_n0_forGradient;
            m_particle_set->accelerations[i] = -m_particle_set->pressure_gradients[i] / m_fluid_density;
        }
    }
}

void Simulator::moveParticleUsingPressureGradient()
{
    for (int i = 0, n = m_particle_set->getNumberOfParticles(); i < n; ++i)
    {
        if (m_particle_set->types[i] == FLUID)
        {
            m_particle_set->velocities[i] += m_particle_set->accelerations[i] * DELTA_TIME;
            m_particle_set->positions[i] += m_particle_set->accelerations[i] * DELTA_TIME * DELTA_TIME;
        }
    }
}
