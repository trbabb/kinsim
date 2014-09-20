/*
 * RigidBody.h
 *
 *  Created on: Apr 11, 2012
 *      Author: tbabb
 */

#ifndef RIGIDBODY_H_
#define RIGIDBODY_H_


#include <geomc/linalg/Vec.h>
#include <geomc/linalg/Quaternion.h>

using namespace geom;

// integrate with runge-kutta
// beware that "perfect" integration is not energy-conserving for non-constant forces.
//   > I think the article's definition of "perfect" uses equations for constant
//     forces while allowing non-constant force inputs.
//   - I wonder how well your "lerped acceleration with kinematic equations" method
//     will compare to RK4.
// beware of body space vs. world space forces/measurements!
// you will have to abstract your integrator so you can play with cleverer methods
// you should allow:
//   world space constant/interpolated forces
//   world space analytic forces (computed dx/dv)
//   body space constant forces (rocket pinned to a body)
//   body space analytic forces (??? super hard)

// when computing forces, allow force applicators to generate dx and dt directly, rather than v and a, 
// so they can solve implicitly over the time interval if they need to.

// dR/dt = omega (angular velocity)
// L = I * omega (angular momentum)
// torque = dL / dt = I * alpha

// constructing off-center MI tensors: 
//   see parallel axis theorem
//   I_offs = I_cm - M * D^2
//   where d is a skew symmetric matrix constructed from the displacement vector.
//   (d x ...  matrix)


// todo: Can we represent moment of inertia in a quaternion-friendly way?

// todo: in the future, these will be long vectors of vectors representing
//       multiple bodies.
//       hopefully we could make addition/multiplication/powers/etc friendly

struct BodyState {
    static const size_t n_linear_values;
public:
    
    Vec3d x;   // position
    Vec3d v;   // velocity
    
    Quatd o;   // orientation
    Vec3d L;   // angular momentum, inertial space
    
    // unchanging state:
    
    double m;  // mass
    Vec3d  J;  // body-space moment of inertia
    
    inline Vec3d apply_mi(const Vec3d &v) const {
         // HAZARD: Multiplication not associative here!
        return ((o * v) / J) * o;
    }
    
    inline Vec3d apply_inv_mi(const Vec3d &v) const {
        return o * ((v * o) * J);
    }
    
    // get the instantaneous angular velocity in inertial space
    inline Vec3d get_omega() const {
        return apply_inv_mi(L);
    }
    
    // position relative to CM / body origin
    // directions in inertial space
    inline void accumulate_force_inertial(Vec3d *f_out, Vec3d *t_out, Vec3d p_in, Vec3d f_in) const {
        *f_out += f_in;
        *t_out += p_in ^ f_in;
        //TODO: this doesn't seem complete. see the top case.
    }

    // position relative to CM / body origin
    // directions in body space
    inline void accumulate_force_body(Vec3d *f_out, Vec3d *t_out, Vec3d p_body, Vec3d f_body) const {
        Vec3d r_in  = o * p_body;
        Vec3d f_in  = o * f_body;
        *f_out += f_in;
        *t_out += r_in ^ f_in;
    }
    
    inline void accumulate_torque_inertial(Vec3d *t_out, Vec3d t_in) const {
        *t_out += t_in;
    }
    
    inline void accumulate_torque_body(Vec3d *t_out, Vec3d t_in) const {
        *t_out += o * t_in;
    }
    
};

struct ForceSystem {
    virtual void sumForces(const BodyState &state, double t, Vec3d *a, Vec3d *tau) = 0;
};

void macc(BodyState *out, double k, const BodyState *x, const BodyState *b, size_t n);

void body_delta(BodyState *out_dBodyDt, const BodyState &state, double t);

void body_delta_for_rk4(BodyState *d_dt, const BodyState *s0, double t, size_t n, void *data);

#endif /* RIGIDBODY_H_ */
