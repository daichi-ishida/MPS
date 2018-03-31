#pragma once

/* for two-dimensional simulation */

constexpr unsigned int DIMENSION                     = 2;
constexpr double PARTICLE_DISTANCE                   = 0.025;
constexpr double DELTA_TIME                          = 0.001;
constexpr unsigned int OUTPUT_INTERVAL               = 20;


/* for three-dimensional simulation */
/*
constexpr unsigned int DIMENSION                     = 3;
constexpr double PARTICLE_DISTANCE                   = 0.075; //èâä˙ó±éqä‘ãóó£
constexpr double DELTA_TIME                          = 0.003;
constexpr unsigned int OUTPUT_INTERVAL               = 2;
*/

constexpr double FINISH_TIME                         = 1.0;
constexpr double KINEMATIC_VISCOSITY                 = 1.0E-6; // ìÆîSê´åWêî
constexpr double FLUID_DENSITY                       = 1000.0;
constexpr double GRAVITY_X                           = 0.0;
constexpr double GRAVITY_Y                           = -9.8;
constexpr double GRAVITY_Z                           = 0.0;
constexpr double RADIUS_FOR_NUMBER_DENSITY           = 2.1 * PARTICLE_DISTANCE;
constexpr double RADIUS_FOR_GRADIENT                 = 2.1 * PARTICLE_DISTANCE;
constexpr double RADIUS_FOR_LAPLACIAN                = 3.1 * PARTICLE_DISTANCE;
constexpr double COLLISION_DISTANCE                  = 0.5 * PARTICLE_DISTANCE;
constexpr double THRESHOLD_RATIO_OF_NUMBER_DENSITY   = 0.97;
constexpr double COEFFICIENT_OF_RESTITUTION          = 0.2; // îΩî≠åWêî
constexpr double COMPRESSIBILITY                     = 0.45E-9;
constexpr double EPS                                 = 0.01 * PARTICLE_DISTANCE;
constexpr double RELAXATION_COEFFICIENT_FOR_PRESSURE = 0.2;
constexpr int GHOST                                  = -1;
constexpr int FLUID                                  = 0;
constexpr int WALL                                   = 2;
constexpr int DUMMY_WALL                             = 3;
constexpr int GHOST_OR_DUMMY                         = -1;
constexpr int SURFACE_PARTICLE                       = 1;
constexpr int INNER_PARTICLE                         = 0;
constexpr int DIRICHLET_BOUNDARY_IS_NOT_CONNECTED    = 0;
constexpr int DIRICHLET_BOUNDARY_IS_CONNECTED        = 1;
constexpr int DIRICHLET_BOUNDARY_IS_CHECKED          = 2;