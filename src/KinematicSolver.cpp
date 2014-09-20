/*
 * KinematicSolver.cpp
 * 
 * A particle filter implementation for solving kinematic state (i.e. position,
 * velocity, orientation, angular velocity).
 */

#include <vector>
#include <geomc/linalg/AffineTransform.h>
#include <geomc/shape/Rect.h>

#include "ParticleFilter.h"
#include "KinematicSolver.h"


// todo: figure out meaning/computation of "true" pdf value,
//       and its relationship to sample area.
// todo: ^^ in the above, the particle count is included in the pdf value;
//       that's wrong. wt is like `dx`; pdf is function value.
//       sum(dx * pdf) should equate to 1. beware of things like choices on 
//       the unit ball. do geometrical transforms hurt it?
// todo: better magnetometer model.
// todo: better gravity model.
// todo: add a pressure altimiter and model for it.
// todo: improve predictor to draw better samples
// todo: predictor may incorporate control inputs


KinematicState::KinematicState():orient(0,0,0,1) {}

/*******************************************
 * Helper functions                        *
 *******************************************/


// sample the value of a normal distribution pdf at x
real_t normal(real_t x, real_t mean, real_t variance) {
    variance = std::abs(variance);
    real_t p = x - mean;
    real_t q = (p * p) / (2 * variance); 
    return std::exp(-q) / std::sqrt(variance * 2 * M_PI);
}


vec3 draw_random_vector(const vec3 &ctr, const vec3 &variance, rng_t *rng) {
    vec3 o;
    for (index_t i = 0; i < 3; i++) {
        d_normal_t N = d_normal_t(ctr[i], std::sqrt(std::abs(variance[i])));
        o[i] = N(*rng);
    }
    return o;
}


// this motherfucker is busteeeeed:
real_t draw_vector_with_pdf(vec3 *out, const vec3 &ctr, const vec3 &variance, rng_t *rng) {
    *out = draw_random_vector(ctr, variance, rng);
    real_t pdf = 1;
    for (index_t i = 0; i < 3; i++) {
        pdf *= normal((*out)[i], ctr[i], variance[i]);
    }
    return pdf;
}

quat draw_quat(rng_t *rng) {
    quat q;
    for (index_t i = 0; i < 4; i++) {
        q[i] = d_normal_t(0, 1)(*rng);
    }
    return q.unit();
}


real_t draw_quat_with_pdf(quat *q, rng_t *rng) {
    *q = draw_quat(rng);
    return 1 / (2 * M_PI * M_PI);
}


real_t vector_likelihood(const vec3 &actual, const vec3 &measurement, const vec3 &variance) {
    real_t p = 1;
    for (index_t i = 0; i < 3; i++) {
        p *= normal(actual[i], measurement[i], variance[i]);
    }
    return p;
}


real_t clipped_vector_likelihood(const vec3  &actual, 
                                 const vec3  &measurement, 
                                 const vec3  &variance, 
                                 const rect3 &measurement_bounds) {
    // this pdf is somewhat degenerate because it has no indefinite integral.
    // however because we always renormalize our samples, and we always intersect
    // this pdf with others which *are* well defined, there is no issue.
    vec3 ctr = measurement_bounds.clamp(measurement);
    vec3 bnds[2] = { measurement_bounds.min(), measurement_bounds.max() };
    real_t p = 1;
    for (index_t i = 0; i < 3; i++) {
        // p(actual value is beyond measurement range) ==
        // area of pdf beyond measurement range.
        real_t edge;
        if (actual[i] < bnds[0][i]) {
            edge = bnds[0][i];
        } else if (actual[i] > bnds[1][i]) {
            edge = bnds[1][i];
        } else {
            p *= normal(actual[i], ctr[i], variance[i]);
            continue;
        }
        real_t dist_to_clip = std::abs(measurement[i] - edge);
        real_t likelihood_clipped = 1 - std::erf(dist_to_clip);
        p *= likelihood_clipped;
    }
    return p;
}


/*******************************************
 * Sensor classes                          *
 *******************************************/


///////////////// Sensor3Axis /////////////////


Sensor3Axis::Sensor3Axis():Observer<vec3, vec3>("3-Axis Sensor"), does_clip(false), variance(1) {};

real_t Sensor3Axis::likelihood(const vec3 &state, const vec3 &observed) {
    vec3 sensor_space = state2reading * state;
    if (does_clip) {
        return clipped_vector_likelihood(sensor_space, observed, variance * 2, clip_bounds);
    } else {
        return vector_likelihood(sensor_space, observed, variance * 2);
    }
}

Sensor3Axis::observation_t Sensor3Axis::simulate(const vec3 &v, rng_t *rng) {
    vec3 reading = state2reading * draw_random_vector(v, variance, rng);
    if (does_clip) {
        reading = clip_bounds.clamp(reading);
    }
    return observation_t(reading, this);
}


///////////////// SensorMagnetometer /////////////////


SensorMagnetometer::SensorMagnetometer():Observer<KinematicState,vec3>("Magnetometer") {}

vec3 SensorMagnetometer::compute_inertial_field(vec3 ref) {
    // todo: this
    // could use DoD World Magnetic Model
    // http://www.ngdc.noaa.gov/geomag/geomag.shtml
    return vec3(1,0,0);
}

real_t SensorMagnetometer::likelihood(const KinematicState &s, const vec3 &m) {
    // we don't care about field strength.
    vec3 field_in = compute_inertial_field(s.x);
    return body_space_sensor.likelihood(field_in.unit() * s.orient, m.unit());
}

SensorMagnetometer::observation_t SensorMagnetometer::simulate(const KinematicState &s, rng_t *rng) {
    // todo: really check this.
    vec3 field_in = compute_inertial_field(s.x);
    return observation_t(body_space_sensor.simulate(field_in * s.orient, rng).measurement, this);
}


///////////////// SensorAccelerometer /////////////////

// xxx: without some constraint on the likely inertial acceleration,
// the accelerometer tells you literally nothing about the orientation.
// ex: acceleration is likely 0, or some other value due to control input,
// or not larger than X. 
// the reason why is: picture the distribution of possible accelerations
// when you have no a priori information-- it's a uniform field over the reals.
// what happens when you translate the field by `g`? Nothing. Nothing at all.
// The moment you have a blob in there-- anywhere-- then your accel reading
// tells you something about your orientation, because some part of that blob
// will be closer to the origin, and certain orientations will change `g` and
// thus move it closer or farther away.

SensorAccelerometer::SensorAccelerometer():Observer<KinematicState,vec3>("Accelerometer") {}

vec3 SensorAccelerometer::gravity_inertial(vec3 p) {
    // todo: model 1/r^2
    // todo: gravity model
    // could later depend on ECEF coordinates, i.e. WGS84 or EGM2008
    // http://earth-info.nga.mil/GandG/wgs84/gravitymod/
    return (c_e - p).unit() * 9.80;
}

real_t SensorAccelerometer::likelihood(const KinematicState &s, const vec3 &m) {
    return body_space_sensor.likelihood((gravity_inertial(s.x) - s.a) * s.orient, m);
    //return body_space_sensor.likelihood(s.a * s.orient, m);
}

SensorAccelerometer::observation_t SensorAccelerometer::simulate(const KinematicState &s, rng_t *rng) {
    vec3 a_body = (gravity_inertial(s.x) - s.a) * s.orient;
    return observation_t(body_space_sensor.simulate(a_body, rng).measurement, this);
    //return observation_t(body_space_sensor.simulate(s.a, rng).measurement, this);
}


///////////////// SensorRateGyro /////////////////
    

SensorRateGyro::SensorRateGyro():Observer<KinematicState,vec3>("Rate Gyro") {}

real_t SensorRateGyro::likelihood(const KinematicState &s, const vec3 &m) {
    return body_space_sensor.likelihood(s.omega * s.orient, m);
}

SensorRateGyro::observation_t SensorRateGyro::simulate(const KinematicState &s, rng_t *rng) {
    return observation_t(body_space_sensor.simulate(s.omega * s.orient, rng).measurement, this);
}


///////////////// SensorGPSPositionXYZ /////////////////


SensorGPSPositionXYZ::SensorGPSPositionXYZ():Observer<KinematicState,vec3>("GPS Position") {}

real_t SensorGPSPositionXYZ::likelihood(const KinematicState &s, const vec3 &m) {
    // todo: your uncertainty might be different in alt, which would
    // require a special pdf computation. 
    vec3 ecef_position = ref2ecef * s.x;
    return ECEF_space_sensor.likelihood(ecef_position, m);
}


///////////////// SensorGPSVelocityXYZ /////////////////


SensorGPSVelocityXYZ::SensorGPSVelocityXYZ():Observer<KinematicState,vec3>("GPS Velocity") {}

real_t SensorGPSVelocityXYZ::likelihood(const KinematicState &s, const vec3 &m) {
    // todo: uncertainty calc is all fuxed.
    vec3 ecef_velocity = ref2ecef.applyVector(s.v);
    return ECEF_space_sensor.likelihood(ecef_velocity, m);
}


///////////////// ConstraintAcceleration /////////////////


ConstraintAcceleration::ConstraintAcceleration():Observer<KinematicState,vec3>("Acceleration constraint"), variance(1) {}

real_t ConstraintAcceleration::likelihood(const KinematicState& s, const vec3& m) {
    return vector_likelihood(s.a, ctr, variance);
}


/*******************************************
 * Kinematic Predictor                     *
 *******************************************/


void KinematicPredictor::drawNextStates(std::vector< Particle<KinematicState> > &out,
              const std::vector< Particle<KinematicState> > &population,
              real_t dt,
              rng_t *rng) {
    // the stupidest thing we could do, for now:
    // assume random jerk / angular acceleration about 0, with 
    // hard coded-variance.
    vec3 jerk_variance  = vec3(2);
    vec3 alpha_variance = vec3(2);
    for (index_t i = 0; i < out.size(); i++) {
        auto &p_i = population[i];
        auto &p_o = out[i];
        KinematicState tmp = p_i.state;
        // stupidly draw some process noise: choose a jerk, angular accel.
        vec3 a, omega;
        p_o = p_i;
        p_o.pdf *= draw_vector_with_pdf(&a,     vec3(0.), jerk_variance,  rng); 
        p_o.pdf *= draw_vector_with_pdf(&omega, vec3(0.), alpha_variance, rng);
        if (p_o.pdf == 0) std::cout << "pdf is 0\n";
        a     = p_i.state.a     + dt * a;
        omega = p_i.state.omega + dt * omega;
        
        // compute dynamics based on chosen noise.
        predict(&(out[i].state), p_i.state, a, omega, dt);
    }
}


// suppose an acceleration and a rotation rate; tell us where we'll be.
void KinematicPredictor::predict(KinematicState *s1, const KinematicState &s0, vec3 accel, vec3 omega, real_t dt) {
    // "perfect" dynamics for linearly changing acceleration:
    real_t dt2 = dt * dt;
    vec3 jerk = (accel - s0.a) / dt;
    s1->x = jerk * dt2 * dt / 6 + s0.a * dt2 / 2 + s0.v * dt + s0.x;
    s1->v = jerk * dt2      / 2 + s0.a * dt      + s0.v;
    s1->a = accel;

    // todo: better dynamics here.
    // i.e. compute angular acceleration for new omega, and get analytic result.
    // you could probably solve rk4 and get a closed form 4th order solution...
    // right now we use verlet integration. kinda qrappy.
    s1->orient = std::exp(quat(dt * s0.omega, 0)) * s0.orient;
    s1->omega  = omega;

    // renormalize if we've drifted appreciably.
    // empirically, it would take quite awhile to drift this far.
    if (std::abs(s1->orient.mag2() - 1) > 1e-7) {
        s1->orient = s1->orient.unit();
    }
}


/*******************************************
 * Notes                                   *
 *******************************************/

///////////// state prediction /////////////

// you want to pick a new, random state to advance by which is very likely
// given the particle we are advancing from. one way to do this is to take its 
// state and jitter it, but the question is "by how much and in what direction
// do we jitter?" so we might choose another nearby particle which is likely 
// given s0, and then jitter the current state according to the discrepancy between
// the two. for example, if we have two clusters around a particular, say, acceleration,
// a sample from cluster2 is not likely given s0 in cluster1; however, the spread/
// discrepancy within cluster1 gives us some idea about the uncertainty of s0's accel.

// we cannot simply select states from other particles, because as the values are traded
// around, the diversity of the population will decrease and trend to zero.
// you *must* add continuous noise to keep the population healthy and the sample
// space well-explored.

// we might also calculate a weighted average of state discrepancy along each 
// axis, based on likelihood, and use that discrepancy to generate jitter noise.
// this would be an n^2 calculation, though, which is likely too slow.

// you will probably need the pdf of the noise. i.e. how likely is the new state
// given the old state? i.e. a large accel is unlikely given a small delta v.
// and we want to pick noise with probability proportional to that likelihood.

// to jitter about existing states, you will need the precise pdf value of that
// state. This will tell you the approximate particles/unit area, which you can
// invert to get the area of the particle. This can then drive your "jitter".

// consider MIS-- samples based on the spread of our estimates (to ensure we
// finely sample the the region of high likelihood) and those based on sensor/time
// uncertainty-- to ensure we properly sample change.
//    i.e. take our supposition about jerk and propagate it through
//    to the other variables. we might pay attention to the average weight
//    of samples from either strategy and use the info to balance the two.
//    note that we should not every balance entirely to one side due to the
//    "stationary state" problem.

// consider ways to make your sampling stratified.

// ----

// note that we only have to explore 6 dimensions: 3 acceleration and 3 angular velocity.
// all other state derives precisely from those.

// we might draw acceleration jitter by choosing reasonable limits on jerk.
// me might get a better idea of reasonable jerk by histogramming delta a, and measured a.
// same with omega. this would require keeping P_prev with each particle.

// might also steal an accelerometer reading and use it ahead of time to guess `a`.
// diff between new and old `a` (and old reading and new reading) tell us 
// something about amount to jitter. similar logic applies for 'omega'.

// but be generous with jitter; imagine how we can "get out" of an estimated
// stationary state-- we must hypothesize a nonzero acceleration.
// however wide jitter creates inefficient sample space coverage, requiring
// more particles.

// finally, we may want to increase noise based on dt, since "more has changed"


///////////// further work /////////////


/*
 *  * other things you can take advantage of:
 *   - response of controls can inform airspeed/angle of attack
 *   - a microphone might even inform airspeed
 *   - meaning of pitot may depend on pressure
 *   - optical gyro could tell us about orientation
 *   - stars could tell us about orientation
 *   - a second low-gain accel/gyro could tell us about quick movements
 *     and refine detailed accel/gyro measurements
 * 
 *   - may want to map/learn space of control input/output so we can
 *     invert it.
 */