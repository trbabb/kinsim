#include <OpenGL/gl.h>

#include "RigidBody.h"
#include "Solver.h"

#include "visible/VisBox3d.h"
#include "GLWindow.h"
#include "AnimTimer.h"
#include "Manipulator.h"

#define RESOLUTION 1280

/*
 * Purpose: 
 * 1. To assess whether the exponentiation method is valid for rotational dynamics.
 * 2. To compare the relative importance of using rk4 vs. euler, and product vs. 
 *    riemann integrals.
 * 
 * Results:
 * - Exponentiation (after a bug fix) matches summation integrals.
 * - Using RK4 is far more important than product integrals.
 *   The two RK4 integrators stick closely together.
 * - Exponentiation performs better in the long run; summation RK4 will drift to
 *   axial rotation. Exponentiation does this too, but far more slowly.
 * - There may be a phase difference between the two RK4 integrators.
 *   It might be related to starting conditions.
 * - There is a factor of 2 on the angle of quaternion magnitudes that must not 
 *   be ignored; this is due to the fact that quaternions are a double cover of 
 *   rotation space. E.g. q and -q are the same rotation.
 * - The coeffs in RK4 seem to be (at least, experimentally) valid for product
 *   integral-based solving.
 * 
 * Follow up:
 * - May want to numerically check advantage of product integrals within RK4.
 *   Earlier test suggests strong numerical stability. Let's check that again and 
 *   make sure it agrees with our new results.
 */


template <typename T>
Vec<T,3> moment_of_inertia(const Rect<T,3> &box) {
    Vec<T,3> dims = box.getDimensions();
    return Vec<T,3>(
            dims.y * dims.y + dims.z * dims.z,
            dims.x * dims.x + dims.z * dims.z,
            dims.x * dims.x + dims.y * dims.y
        )  / 12.0;
}


struct Top : virtual public Animated,
             virtual public Drawable,
             virtual public ForceSystem {
    
    // cm is at body space origin.
    
    BodyState s;
    Rect3d box;
    Vec3d color;
    int substeps;
    
    Top(Vec3d color=ONE_VEC3d, int substeps=1):
            box(Rect3d::fromCenter(ZERO_VEC3d, Vec3d(1, 1.33, 1/3.0))),
            color(color),
            substeps(substeps) {
        const Vec3d axis = Vec3d(0.25,0,1).unit();
        s.o = Quatd::rotDirectionAlign(Z_AXIS3d, axis);
        s.m = 1;
        s.J = s.m * moment_of_inertia(box);
        s.x = Vec3d(0, 0, 1/3.);
        s.L = (Z_AXIS3d * 30);
    }
    
    void sumForces(const BodyState &state, double t, Vec3d *f, Vec3d *tau) {
        // todo: so far, only torque-free rotation.
        *f   = 0.0;
        *tau = 0.0;
    }
    
    void applyConstraints(BodyState *state, index_t n, double t, double dt) {
        // do nothing.
    }
    
    virtual void do_advance(BodyState *s1, double t, double dt) {
        BodyState b0, b1;
        rk4_advance(s1, 1,                // one output state
                    t, dt,                // times
                    &s,                   // current state
                    body_delta_euler_for_rk4,   // delta function
                    body_combine_add_mul, // accumulation function
                    (ForceSystem*)this,   // extra data (convert the ptr now, while you still have type info)
                    &b0, &b1);            // working buffers
    }
    
    void init(double t) {}
    
    void update(double t, double dt) {
        const int n_steps = substeps;
        double dtt = dt / n_steps;
        
        for (int n = 0; n < n_steps; n++) {
            BodyState s1;

            // advance the simulation
            do_advance(&s1, t + n * dtt, dtt);
            
            s = s1;
        }
    }
    
    void draw() {
        glPushMatrix();
        glColor(color);
            glTranslate(s.x);
            glRotate(s.o);
            VisBox3d(box, true).draw();
            glBegin(GL_LINES);
                glVertex3d(0, 0, -1/3.0);
                glVertex3d(0, 0,  1/3.0);
                
                //glVertex3d(0,0,0);
                //glVertex((s.get_omega() * s.o) / 20);
            glEnd();
            
        glPopMatrix();
    }
    
};


struct TopExp : public Top {
    
    TopExp(Vec3d color=ONE_VEC3d, int substeps=1):Top(color,substeps) {}
    
    void do_advance(BodyState *s1, double t, double dt) {
        BodyState b0, b1;
        rk4_advance(s1, 1,                // one output state
                    t, dt,                // times
                    &s,                   // current state
                    body_delta_for_rk4,   // delta function
                    body_combine_mul_exp, // accumulation function
                    (ForceSystem*)this,   // extra data (convert the ptr now, while you still have type info)
                    &b0, &b1);            // working buffers
    }
    
};

struct TopEuler : public Top {
    
    TopEuler(Vec3d color=ONE_VEC3d, int substeps=1):Top(color,substeps) {}
    
    void do_advance(BodyState *s1, double t, double dt) {
        BodyState b0;
        euler_advance(s1, 1,              // one output state
                    t, dt,                // times
                    &s,                   // current state
                    body_delta_euler_for_rk4,   // delta function
                    body_combine_add_mul, // accumulation function
                    (ForceSystem*)this,   // extra data (convert the ptr now, while you still have type info)
                    &b0);                 // working buffer
    }
    
};


struct TopEulerExp : public Top {
    
    TopEulerExp(Vec3d color=ONE_VEC3d, int substeps=1):Top(color,substeps) {}
    
    void do_advance(BodyState *s1, double t, double dt) {
        BodyState b0;
        euler_advance(s1, 1,              // one output state
                    t, dt,                // times
                    &s,                   // current state
                    body_delta_for_rk4,   // delta function
                    body_combine_mul_exp, // accumulation function
                    (ForceSystem*)this,   // extra data (convert the ptr now, while you still have type info)
                    &b0);                 // working buffer
    }
    
};


int main(int argc, char *argv[]) {
    GLWindow win(&argc, argv, "rotsim", RESOLUTION, RESOLUTION);
    AnimTimer timer(&win);
    
    Camera &cam = win.cam;
    cam.setFov(0.4 * M_PI);
    cam.setFar(5000);
    cam.setNear(0.1);
    cam.setPosition(-X_AXIS3d * 5);
    cam.setUp(Z_AXIS3d);
    cam.setCenterOfInterest(ZERO_VEC3d);
    
    Manipulator manip(&win, &win.cam);
    
    Top         t_rk_add(Vec3d(1, 0, 1));
    TopExp      t_rk_exp(Vec3d(0, 0, 1));
    TopEuler    t_eu_add(Vec3d(1, 1, 0));
    TopEulerExp t_eu_exp(Vec3d(0, 1, 0));
    
    win.scene.push_back(&t_rk_add);
    win.scene.push_back(&t_rk_exp);
    //win.scene.push_back(&t_eu_add);
    //win.scene.push_back(&t_eu_exp);
    
    timer.anims.push_back(&t_rk_add);
    timer.anims.push_back(&t_rk_exp);
    timer.anims.push_back(&t_eu_add);
    timer.anims.push_back(&t_eu_exp);
    
    win.setClearColor(Vec3d(0.01));
    
    timer.begin();
    win.showAll();
    return 0;
    
}