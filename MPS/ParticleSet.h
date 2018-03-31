#pragma once
#include <vector>
#include <Eigen/Core>

class ParticleSet
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    ParticleSet();
    ~ParticleSet();

    void initSet(int type, double x, double y, double v_x, double v_y);
    void initSet(int type, double x, double y, double z, double v_x, double v_y, double v_z);

    const int getNumberOfParticles() { return m_num; };

    std::vector<int>             types;
    std::vector<int>             boundary_conditions;
    std::vector<Eigen::Vector3d> positions;
    std::vector<Eigen::Vector3d> velocities;
    std::vector<Eigen::Vector3d> accelerations;
    std::vector<Eigen::Vector3d> pressure_gradients;
    std::vector<double>          pressures;
    std::vector<double>          min_pressures;
    std::vector<double>          number_densities;
    std::vector<int>             flag_for_checking_boundary_conditions;

private:
    int m_num;
};

