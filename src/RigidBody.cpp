/*
 * RigidBody.cpp
 *
 *  Created on: Dec 28, 2012
 *      Author: tbabb
 */

#include "RigidBody.h"

// multiply/accum a BodyMotion, doing "intermediate" quaternion math more sensibly.
void macc(BodyState *o, double k, const BodyState *d_ds, const BodyState *s, size_t n) {
    BodyState *end = o + n;
    for (; o != end; o++, d_ds++, s++) {
        o->x = s->x + k * d_ds->x;
        o->v = s->v + k * d_ds->v;
        o->L = s->L + k * d_ds->L;
        o->m = s->m + k * d_ds->m; // probably not changing, but why not?
        o->J = s->J + k * d_ds->J;
        
        // now do the quaternion magic.
        // (this is to approximate a product integral rather than a summation integral)
        o->o = std::exp(k * d_ds->o) * s->o;
        //o->o = s->o + k * d_ds->o; // david baraff version
    }
}

// would be better if it were clearer which variables were integrated
// and which were intermediate state / directly computed
void body_delta(BodyState *d_dt, const BodyState &state, double t, ForceSystem* fs) {
    Vec3d f, tau;
    
    fs->sumForces(state, t, &f, &tau);
    
    d_dt->x = state.v;
    d_dt->v = f / state.m;
    
    // J(t)                 // inertial space angular momentum
    // dR / dt = omega      // angular velocity.
    // domega/dt = alpha    // angular acceleration.
    // dL/dt = tau          // torque.
    // L = J(t) * omega
    
    // q * v := rotate v by q.
    // v * q := rotate v by q`.
    
    d_dt->o   = Quatd(state.get_omega(), 0); //0.5 * Quatd(state.get_omega(),0) * state.o; <- deb version
    d_dt->L   = tau;
    
    d_dt->J   = 0.0;
    d_dt->m   = 0.0;
    
    // J(t)  = R`(t) * J_0  * R(t)
    // J`(t) = R(t)  * J_0` * R`(t)
}

void body_delta_for_rk4(BodyState *d_dt, const BodyState *s0, double t, size_t n, void *data) {
    for (size_t i = 0; i < n; i++) {
        body_delta(d_dt + i, s0[i], t, (ForceSystem*)data);
    }
}