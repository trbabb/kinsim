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

// constraints:
//   impose constraints after each integration step.
//   allow rigid body to adjust internal state given constraint adjustment
//     position constraints, e.g.
//     position- adjusting constraints may also choose to correct velocity to match.
//   see also: Sequential Impulse method

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
    
    Vec3d x;   // position, world space
    Vec3d v;   // velocity, world space
    
    Quatd o;   // orientation (body to world)
    Vec3d L;   // angular momentum, inertial space
    
    // unchanging state:
    
    double m;  // mass
    Vec3d  J;  // body-space (really, prinicple axis-space) moment of inertia
    
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
    inline void accumulate_force_inertial(Vec3d *f_out, Vec3d *t_out, Vec3d p_body, Vec3d f_in) const {
        Vec3d p_in = o * p_body;
        *f_out += f_in;
        *t_out += p_in ^ f_in;
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
    
    // return velocity of point, incl. rotation; inertial space
    inline Vec3d get_velocity_point_body(Vec3d p_body) const {
        return (get_omega() ^ (o * p_body)) + v;
    }
    
    // return world-oriented centripetal acceleration of p_body about CM.
    inline Vec3d get_centripetal_accel_point_body(Vec3d p_body) const {
        // find instantaneous v
        Vec3d p_in     = o * p_body;
        Vec3d omega_in = get_omega();
        Vec3d v_in     = omega_in ^ p_in;
        double o_mag2  = omega_in.mag2();
        
        // r perpendicular to axis
        double d       = p_in.dot(omega_in);
        Vec3d r_in     = p_in - (d / o_mag2) * omega_in;
        double r2      = r_in.mag2();
        
        // a = r_hat * v^2 / r
        return (r2 == 0 or o_mag2 == 0) ? ZERO_VEC3d : (-r_in * v_in.mag2() / r2);
    }
    
};

struct ForceSystem {
    virtual void sumForces(const BodyState &state, double t, Vec3d *a, Vec3d *tau) = 0;
    virtual void applyConstraints(BodyState *state, index_t n, double t, double dt) = 0;
};

void body_combine_add_mul(BodyState *out, double k, const BodyState *x, const BodyState *b, size_t n);
void body_combine_mul_exp(BodyState *out, double k, const BodyState *x, const BodyState *b, size_t n);
// void          (*accum)(State*,         T,        const State*,       const State*,       size_t)
void inline          macc(BodyState *out, double k, const BodyState *x, const BodyState *b, size_t n) { body_combine_mul_exp(out, k, x, b, n); }
void body_delta(BodyState *out_dBodyDt, const BodyState &state, double t);
void body_delta_for_rk4(BodyState *d_dt, const BodyState *s0, double t, size_t n, void *data);

#endif /* RIGIDBODY_H_ */
