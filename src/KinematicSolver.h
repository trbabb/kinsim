/* 
 * File:   KinematicSolver.h
 * Author: tbabb
 *
 * Created on June 24, 2014, 1:05 AM
 */

#ifndef KINEMATICSOLVER_H
#define	KINEMATICSOLVER_H

#include <geomc/linalg/AffineTransform.h>
#include "ParticleFilter.h"
#include "KinematicState.h"


// sample the value of a normal distribution pdf at x
real_t normal(real_t x, real_t mean, real_t variance);
real_t vector_likelihood(const vec3 &actual, const vec3 &measurement, const vec3 &variance);
vec3   draw_random_vector(const vec3 &ctr, const vec3 &variance, rng_t *rng);
real_t draw_vector_with_pdf(vec3 *out, const vec3 &ctr, const vec3 &variance, rng_t *rng);
quat   draw_quat(rng_t *rng);
real_t draw_quat_with_pdf(quat *q, rng_t *rng);


real_t clipped_vector_likelihood(const vec3  &actual, 
                                 const vec3  &measurement, 
                                 const vec3  &variance, 
                                 const rect3 &measurement_bounds);



/*******************************************
 * Sensor classes                          *
 *******************************************/


class Sensor3Axis : public Sensor<vec3, vec3> {
public:
    
    using observation_t = Sensor<vec3, vec3>::observation_t;
    
    Sensor3Axis();
    
    transform_t state2reading; // to initially obtain imperically
    rect3 clip_bounds;
    bool does_clip;
    vec3 variance; // could be a matrix?
                   // i.e. could linearly depend on reading
                   // ...and there could be a covariance component to the noise too.
    
    virtual real_t likelihood(const vec3 &state, const vec3 &observed);
    observation_t simulate(const vec3 &v, rng_t *rng);
    
};


class SensorSpatial : public Sensor<KinematicState, vec3> {
public:
    
    SensorSpatial(std::string s) : Sensor<KinematicState,vec3>(s) {}
    SensorSpatial() {}
    
    using observation_t = Sensor<KinematicState, vec3>::observation_t;
    Sensor3Axis body_space_sensor;
};



class SensorMagnetometer : public SensorSpatial {
public:
    
    SensorMagnetometer();
    
    
    transform_t ref2ecef;
    
    vec3 compute_inertial_field(vec3 ref);
    
    virtual real_t likelihood(const KinematicState &s, const vec3 &m);
    observation_t simulate(const KinematicState &s, rng_t *rng);
    
};


class SensorAccelerometer : public SensorSpatial {
public:
    
    SensorAccelerometer();
    vec3 c_e; // center of the earth in reference space.
    
    vec3 gravity_inertial(vec3 p);
    
    virtual real_t likelihood(const KinematicState &s, const vec3 &m);
    observation_t simulate(const KinematicState &s, rng_t *rng);
    
};


class SensorRateGyro : public SensorSpatial {
public:
    
    SensorRateGyro();
    
    virtual real_t likelihood(const KinematicState &s, const vec3 &m);
    observation_t simulate(const KinematicState &s, rng_t *rng);
    
};


// todo: measurement type will need to include the precision estimate.
// todo: oh yeah. gps fixes correspond to readings 1 second in the past.
//       that makes everything all fucked. we are going to have to increase
//       the likelihood of the current sample based on where it was a whole
//       second ago!
//       note that the history could be simplified by ramer-douglas-peucker.
class SensorGPSPositionXYZ : public Sensor<KinematicState, vec3> {
public:
    
    SensorGPSPositionXYZ();
    
    Sensor3Axis ECEF_space_sensor; // this will probably go away
    transform_t ref2ecef;
    
    virtual real_t likelihood(const KinematicState &s, const vec3 &m);
    
};


class SensorGPSVelocityXYZ : public Sensor<KinematicState, vec3> {
public:
    
    SensorGPSVelocityXYZ();
    
    Sensor3Axis ECEF_space_sensor; // this will probably go away
    transform_t ref2ecef;
    
    virtual real_t likelihood(const KinematicState &s, const vec3 &m);
    
};


// could be guided by control input
// models an observation that acceleration "should be" near a certain value
class ConstraintAcceleration : public Sensor<KinematicState, vec3> {
public:
    
    using observation_t = Sensor<KinematicState, vec3>::observation_t;
    
    ConstraintAcceleration();
    
    vec3 variance;
    
    virtual real_t likelihood(const KinematicState &s, const vec3 &m);
    
};


/*******************************************
 * Predictor class                         *
 *******************************************/


// todo: I think the pdf calc is still busted. it should not depend
//       on the number of samples. the weight is your dX, the pdf is 
//       your function value. sum(dX * pdf) for all dX should equal 1.
//       so the 1/N factor is accounted for.
class KinematicPredictor : public Predictor<KinematicState> {
public:
    // todo: expose jerk / alpha assumptions
    
    // this predictor could internally keep track of interesting things,
    // like the distribution of jerks/alphas that worked well.
    
    // for now, we shall assume the population size is the same as the 
    // new population size.
    void drawNextStates(std::vector< Particle<KinematicState> > &out,
                  const std::vector< Particle<KinematicState> > &population,
                  real_t dt,
                  rng_t *rng);
    
    // suppose an acceleration and a rotation rate; tell us where we'll be.
    void predict(KinematicState *s1, const KinematicState &s0, vec3 accel, vec3 omega, real_t dt);
    
};


#endif	/* KINEMATICSOLVER_H */

