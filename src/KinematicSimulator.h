/* 
 * File:   KinematicSimulator.h
 * Author: tbabb
 *
 * Created on June 26, 2014, 1:17 AM
 */

/*
 
 on finite impulses:

 limitations
  - does not conserve angular momentum
  - friction not proportional to F_normal
  - Baumgaute term is tweaky
  - rotational integration is linear; not product.
  - constraints are pairwise
  - big mass ratios don't work well.

(can we "triangulate" an impulse for each V that does not change its 
magnitude? X no, it must be a V that does not add energy. which is different)

Can we explicitly add an energy constraint?
  is this perhaps added to each kind of constraint explicitly?
 
 
 */

#ifndef KINEMATICSIMULATOR_H
#define	KINEMATICSIMULATOR_H

#include <cmath>

#include "KinematicSolver.h"
#include "ParticleFilter.h"

void macc(KinematicState *o, double k, const KinematicState *d_ds, const KinematicState *s, size_t n);
void delta(KinematicState *d_dt, const KinematicState *s0, double t, size_t n, void *data);


inline real_t cycle(real_t t, real_t hz) {
    return std::sin(t * M_PI * 2 * hz);
}


class SensorSimulator {
    /*
     * This class keeps a ground-truth state and advances it according to physics
     * plus some random buffeting. The ground truth is used to generate simulated
     * sensor readings, which can then be fed to a KinematicSolver for state estimation.
     */
public:
    rng_t *rng;
    ParticleFilter<KinematicState> filter;
    KinematicState truth;
    
    SensorMagnetometer     s_mag;
    SensorAccelerometer    s_acc;
    SensorAccelerometer    s_acc2;
    SensorRateGyro         s_gyr;
    ConstraintAcceleration c_acc;
    real_t t;
    
    typename SensorMagnetometer::observation_t     mag_obs;
    typename SensorAccelerometer::observation_t    acc_obs;
    typename SensorRateGyro::observation_t         gyr_obs;
    typename ConstraintAcceleration::observation_t acc_cns;
    
    ObservationBase<KinematicState>* latest_obs;
    
    SensorSimulator(index_t n_particles);
    SensorSimulator(index_t n_particles, rng_t *rng);
    
    void init_state();
    void init_population();
    void init_sensors();
    void advance_dt(real_t dt);
    vec3 jerk(real_t t);
    vec3 alpha(real_t t);
    
};

#endif	/* KINEMATICSIMULATOR_H */

