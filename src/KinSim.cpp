#include <iostream>
#include <OpenGL/gl.h>
#include <algorithm>
#include <cmath>

#include "GLWindow.h"
#include "GUIListener.h"
#include "AnimTimer.h"
#include "Manipulator.h"
#include "glHelpers.h"
#include "visible/VisAxis.h"
#include "visible/StateUtils.h"
#include "visible/VisBall.h"
#include "visible/VisCallback.h"
#include "visible/VisBox3d.h"

#include "KinematicSimulator.h"


#define N_PARTICLES 256

/*
 * Algorithm input parameters:
 *   - jerk process variance
 *   - alpha process variance
 *     (the above should roughly match your jerk() and alpha() force functions)
 *   - variances of your initial population
 *     (true omega is chosen randomly, variance of omega should encapsulate that)
 *     (other variables should have enough variation to accommodate the initial
 *      timestep)
 *   - sensor variances
 *     (precise sensors might actually make it hard to estimate state, because
 *      the areas of nonzero pdf are too small to hit).
 *     (maybe MIS could be useful here)
 * 
 * Problems:
 *   - the gyro reading looks biased.
 *   - sometimes the initial state will get a wild, fast rotation that doesn't match
 *     at all. The estimates will then never converge on the true solution
 *     because the sampled pdf is essentially zero; so the estimator doesn't know
 *     where to look. (In this case, we set all pdfs to 1/n and allow the space to 
 *     sort itself out)
 *   - If we have a reasonable estimate of orientation (which we can get), then
 *     acceleration should be easy to get because we can subtract our estimated gravity.
 *     Without gravity at all, acceleration should track the truth relatively well.
 *     This is not borne out. 
 *   - It seems that adding the accelerometer measurement to the simulation actually 
 *     hurts the accuracy. (dimensionality problem?)
 *   - Imagine the case where the orientation estimate is dead on, but the position/accel
 *     estimate is off. Especially if our cdf is zero, we should prefer these solutions.
 *     currently, all will be treated equally.
 *   - our gyro, e.g. can tell us something directly about possible states to explore.
 *     > i.e. importance sampling.
 *   - would be nice if we sampled nearby state space directly and not just advancements from
 *     current configuration.
 * 
 * other thoughts:
 *     maybe more than sampling we need root-finding?
 *     if we know the likelihood function we have some idea of where to look.
 */


void draw_kstate(const KinematicState &ks, vec3 color) {
    glPushMatrix();
        //glTranslate(ks.x);
        glRotate(ks.orient);
        glColor(color);
        VisBox3d(Rect3d::fromCenter(ZERO_VEC3d, Vec3d(1,3,9))).draw_wireframe();
        //VisAxis().draw();
    glPopMatrix();
    
    glBegin(GL_POINTS);
        glColor3d(1,0,0);
        glVertex(ks.a);
        glColor3d(0,0,1);
        glVertex(ks.omega);
        glColor3d(0,1,1);
        glVertex(ks.v);
    glEnd();
}


class ProgramState : virtual public Animated, 
                     virtual public Drawable,
                     virtual public GUIListener {
public:
    
    SensorSimulator ssim;
    real_t time;
    
    ProgramState():ssim(N_PARTICLES),time(0) {}
    
    
    void draw() {
        KinematicState guess = ssim.filter.estimate();
        KinematicState truth = ssim.truth;
        glColor3d(0.125,0.77,1.0);
        draw_kstate(truth, vec3(1.));
        //glColor3d(1.0, 0, 0);
        draw_kstate(guess, vec3(1,0,0));
        //glColor3d(0.5,0.5,0.5);
        for (int i = 0; i < ssim.filter.particles.size(); i++) {
            draw_kstate(ssim.filter.particles[i].state, vec3(0.25));
        }
        
        draw_kstate(ssim.filter.estimate(), vec3(1,1,0));
        
        glBegin(GL_LINES);
            glColor3d(1,0,0);
            glVertex(vec3(0.));
            glVertex(truth.a);
            
            glColor3d(0,0,1);
            glVertex(vec3(0.));
            glVertex(truth.omega);
            
            glColor3d(0,1,1);
            glVertex(vec3(0.));
            glVertex(truth.v);
        glEnd();
        
        if (ssim.latest_obs) {
            Observation<KinematicState,vec3>* obs = (Observation<KinematicState,vec3>*)ssim.latest_obs;
            vec3 c;
            if      (obs->name() == "Magnetometer")  c = vec3(0,1,0);
            else if (obs->name() == "Accelerometer") c = vec3(1,0,0);
            else if (obs->name() == "Rate Gyro")     c = vec3(0,0,1);
            vec3 o = obs->measurement;
            //SensorSpatial *sp;
            //if (sp = dynamic_cast<SensorSpatial*>(obs->observer)) { o = o / sp->body_space_sensor.state2reading; }
            glColor(c);
            VisBall(o, 0.08).draw();
        }
    }
    
    void init(double t) {
        ssim.t = time = t;
    }
    
    void update(double t, double dt) {
        const int n = 4;
        for (int i = 0; i < n; i++) {
            ssim.advance_dt(dt / n);
        }
    }
    
    void reset() {
        ssim   = SensorSimulator(N_PARTICLES, ssim.rng);
        ssim.t = time;
    }
    
    bool keyEvent(GLWindow* window, unsigned char key, int keycode, bool down, bool special, int x, int y) {
        if (key == ' ' and down) {
            reset();
            return true;
        }
        return false;
    }

    
};


int main(int argc, char **argv) {
    GLWindow win(&argc, argv, "kinematic solver", 1280, 1280);
    AnimTimer timer(&win);
    
    Camera &cam = win.cam;
    cam.setPosition(-X_AXIS3d * 15);
    cam.setCenterOfInterest(ZERO_VEC3d);
    cam.setUp(Z_AXIS3d);
    cam.setFar(250);
    cam.setNear(0.1);
    
    Manipulator manip(&win, &win.cam);
    
    ProgramState s;
    
    win.scene.push_back(&s);
    timer.anims.push_back(&s);
    win.guiListeners.push_back(&s);
    
    timer.begin();
    win.showAll();
}