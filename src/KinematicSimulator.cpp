/* 
 * KinematicSimulator.cpp
 * 
 * Simulate the dynamics of a body, and attempt to estimate them with a 
 * particle filter. For verification/assessment of filter robustness.
 * 
 */

#include "KinematicSimulator.h"
#include "Solver.h"

#define ALPHA_SCALE 1

void macc(KinematicState *o, double k, const KinematicState *d_ds, const KinematicState *s, size_t n) {
    KinematicState *end = o + n;
    for (; o != end; o++, d_ds++, s++) {
        o->x = s->x + k * d_ds->x;
        o->v = s->v + k * d_ds->v;
        o->a = s->a + k * d_ds->a;
        
        o->orient = std::exp(k * d_ds->orient) * s->orient;
        o->omega  = s->omega + k * d_ds->omega;
    }
}


void delta(KinematicState *d_dt, const KinematicState *s0, double t, size_t n, void *data) {
    SensorSimulator *sim = (SensorSimulator*) data;
    for (index_t i = 0; i < n; i++) {
        d_dt[i].x = s0[i].v;
        d_dt[i].v = s0[i].a;
        d_dt[i].a = sim->jerk(t);
        
        d_dt[i].orient = quat(s0[i].omega, 0);
        d_dt[i].omega  = sim->alpha(t);
    }
}

// things to model:
// - displacement of sensor from center
// - delta over sample period

    
SensorSimulator::SensorSimulator(index_t n_particles):
        rng(new rng_t(11937294775LL)),
        filter(n_particles, new KinematicPredictor(), rng), 
        t(0),
        latest_obs(NULL) {
    init_sensors();
    init_state();
    init_population();
}


SensorSimulator::SensorSimulator(index_t n_particles, rng_t *rng):
        rng(rng),
        filter(n_particles, new KinematicPredictor(), rng), 
        t(0) {
    init_sensors();
    init_state();
    init_population();
}


void SensorSimulator::init_state() {
    truth.x = 0.;
    truth.v = 0.;
    truth.a = 0.;
    truth.orient = quat(0,0,0,1);
    truth.omega  = draw_random_vector(0., 0.25, rng);
}


void SensorSimulator::init_population() {
    for (index_t i = 0; i < filter.particles.size(); i++) {
        Particle<KinematicState> &p = filter.particles[i];
        p.state = KinematicState();
        vec3 dO;
        p.pdf *= draw_vector_with_pdf(&dO,            0., 0.01, rng);
        p.pdf *= draw_vector_with_pdf(&p.state.a,     0., 0.01,  rng);
        p.pdf *= draw_vector_with_pdf(&p.state.omega, 0., 1.5,   rng);
        
        p.state.orient = std::exp(quat(dO,0)) * p.state.orient;
    }
}


void SensorSimulator::init_sensors() {
    //s_acc.body_space_sensor.state2reading = rotation(draw_quat(rng)) * translation(vec3(-1)) * scale(vec3(0.1));
    //s_gyr.body_space_sensor.state2reading = rotation(draw_quat(rng)) * scale(vec3(0.5/M_PI));
    //s_mag.body_space_sensor.state2reading = rotation(draw_quat(rng));

    s_acc.c_e = s_acc2.c_e = vec3(0,0,-6378100);
    s_acc.body_space_sensor.clip_bounds = rect3(vec3(-6), vec3(6));
    //s_acc.body_space_sensor.does_clip   = true;
    s_acc.body_space_sensor.variance    = vec3(0.125 * 0.125); // convert from stddev
    s_acc2.body_space_sensor.variance    = vec3(0.0125);

    s_gyr.body_space_sensor.clip_bounds = rect3(vec3(-6), vec3(6));
    //s_gyr.body_space_sensor.does_clip   = true;
    s_gyr.body_space_sensor.variance    = vec3(0.025 * 0.025);

    s_mag.body_space_sensor.variance    = vec3(0.025 * 0.025);
    
    c_acc.ctr      = vec3(0.);
    c_acc.variance = vec3(0.1);
}


void SensorSimulator::advance_dt(real_t dt) {
    int n_subsims = 3;
    KinematicState s0 = truth;
    KinematicState s1;
    KinematicState b0, b1;

    // advance the ground truth position/orientation
    real_t dt_i = dt / n_subsims;
    for (int i = 0; i < n_subsims; i++) {
        real_t t_i = t + i / (real_t)n_subsims;
        rk4_advance(&s1, 1, t_i, dt_i, &s0, &delta, this, &b0, &b1);
        s0 = s1;
    }
    truth = s0;

    std::vector<ObservationBase<KinematicState>*>  obs(1);

    switch (std::uniform_int_distribution<int>(0,3)(*rng)) {
        case 0:
            mag_obs = s_mag.simulate(s0, rng);
            obs[0] = &mag_obs;
            break;
        case 1:
            // xxx: I don't know why this isn't tracking the true value.
            acc_obs = s_acc.simulate(s0, rng);
            obs[0] = &acc_obs;
            break;
        case 2:
            gyr_obs = s_gyr.simulate(s0, rng);
            obs[0] = &gyr_obs;
            break;
        case 3:
            acc_obs = s_acc2.simulate(s0, rng);
            obs[0] = &acc_obs;
            break;
    }
    latest_obs = obs[0];
    // xxx debug: arbitrarily constrain accel info.
    //acc_cns = ConstraintAcceleration::observation_t(0., &c_acc);
    //obs[1]  = &acc_cns;
 
    // update the filter's estimate, given some jittered sensor readings.
    filter.advance(obs, dt);

    t += dt;
}


vec3 SensorSimulator::jerk(real_t t) {
    //return vec3(cycle(t, 1), cycle(t + 0.77, 0.48776), cycle(t - 0.1693, 1.2284));
    return 0.;
}


vec3 SensorSimulator::alpha(real_t t) {
    return vec3(cycle(t - 2.14159, 0.89216), cycle(t + 0.14, 2.15216), cycle(t - 0.5411, 1.6623)) * ALPHA_SCALE;
}
