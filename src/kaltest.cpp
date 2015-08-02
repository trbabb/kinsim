
#include <iostream>
#include <OpenGL/gl.h>
#include <algorithm>
#include <cmath>
#include <random>
#include <time.h>

#include <geomc/function/Dual.h>

#define SID_ACC  1
#define SID_VEL  2
#define SID_GYR  3
#define SID_MAG  4

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

#include "KalmanFilter.h"
#include "Solver.h"

#define OMEGA_PROCESS_VARIANCE   0.25
#define ACCEL_PROCESS_VARIANCE   0.025

#define MAX_READINGS 12


class SensorSimulator;

/*
 * Observation:
 * Having *any* information, however crappy, about velocity, vastly
 * improves the estimate of everything.
 */


/********** simulation stuff **********/

vec3 cmp_accel(real_t t, rng_t *rng) {
    // white noise shall buffet our object
    return draw_random_vector(0.0, ACCEL_PROCESS_VARIANCE, rng);
}

vec3 cmp_omega(real_t t, vec3 base_omega, rng_t *rng) {
    // buffet with white noise
    return draw_random_vector(base_omega, OMEGA_PROCESS_VARIANCE, rng);
}

/********** simulation stuff **********/

struct delta_data {
    rng_t *rng;
    real_t dt;
    vec3 base_omega;
};

void macc(KinematicState<real_t> *o, double k, const KinematicState<real_t> *d_ds, const KinematicState<real_t> *s, size_t n) {
    KinematicState<real_t> *end = o + n;
    for (; o != end; o++, d_ds++, s++) {
        o->x = s->x + k * d_ds->x;
        o->v = s->v + k * d_ds->v;
        o->a = s->a + k * d_ds->a;
        
        o->orient = std::exp(k * d_ds->orient / 2) * s->orient;
        o->omega  = s->omega + k * d_ds->omega;
    }
}


void delta(KinematicState<real_t> *d_dt, const KinematicState<real_t> *s0, double t, size_t n, void *data) {
    delta_data *dat = (delta_data*) data;
    for (size_t i = 0; i < n; i++) {
        d_dt[i].x = s0[i].v;
        d_dt[i].v = s0[i].a;
        d_dt[i].a = (cmp_accel(t, dat->rng) - s0[i].a) / dat->dt;
        
        d_dt[i].orient = quat(s0[i].omega, 0);
        d_dt[i].omega  = (cmp_omega(t, dat->base_omega, dat->rng) - s0[i].omega) / dat->dt;
    }
}

/********** draw helpers **********/

void draw_kstate(const KinematicState<real_t> &ks, vec3 color) {
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


/********** simulator **********/


class SensorSimulator {
    /*
     * This class keeps a ground-truth state and advances it according to physics
     * plus some random buffeting. The ground truth is used to generate simulated
     * sensor readings, which can then be fed to a KinematicSolver for state estimation.
     */
public:
    rng_t *rng;
    KalmanFilter<real_t> *filter;
    KinematicState<real_t> truth;
    
    SensorMagnetometer<real_t>   s_mag;
    SensorAccelerometer<real_t>  s_acc;
    SensorGPSVelocity<real_t>    s_vel;
    SensorRateGyro<real_t>       s_gyr;
    real_t t;
    vec3 base_omega;
    
    real_t measure_buf[3 * 5];
    Measurement<real_t> latest_obs;
    bool measure_available;
    
    SensorSimulator():
            rng(new rng_t(11937294775LL)),
            filter(new KalmanFilter<real_t>(KINSTATE_SIZE, 
                                            MAX_READINGS,
                                            new KinematicPredictor<real_t>(
                                                    ACCEL_PROCESS_VARIANCE, 
                                                    OMEGA_PROCESS_VARIANCE))), 
            t(0) {
        init_sensors();
        init_state();
    }
    
    SensorSimulator(rng_t *rng):
            rng(rng),
            filter(new KalmanFilter<real_t>(KINSTATE_SIZE, 
                                            MAX_READINGS,
                                            new KinematicPredictor<real_t>(
                                                    ACCEL_PROCESS_VARIANCE, 
                                                    OMEGA_PROCESS_VARIANCE))), 
            t(0) {
        init_sensors();
        init_state();
    }
    
    void init_state() {
        measure_available = false;
        truth.x = 0.;
        truth.v = 0.;
        truth.a = 0.;
        truth.orient = quat(0,0,0,1);
        truth.omega  = base_omega = draw_random_vector(0., 0.25, rng);
        std::copy(truth.begin(), truth.end(), filter->x);
    }
    
    void init_sensors() {
        
        real_t sens_var = pow(0.025, 2);
        
        s_acc._id = SID_ACC;
        s_vel._id = SID_VEL;
        s_gyr._id = SID_GYR;
        s_mag._id = SID_MAG;

        s_acc.c_e = vec3(0,0,-6378100);
        s_acc.set_variance(sens_var);
        s_vel.set_variance(1);
        s_gyr.set_variance(sens_var);
        s_mag.set_variance(sens_var);
    }
    
    double advance_dt(real_t dt) {
        int n_subsims = 8;
        KinematicState<real_t> s0 = truth;
        KinematicState<real_t> s1;
        KinematicState<real_t> b0, b1;
        
        delta_data dat = {rng, dt, base_omega};

        // advance the ground truth position/orientation
        real_t dt_i = dt / n_subsims;
        for (int i = 0; i < n_subsims; i++) {
            real_t t_i = t + i / (real_t)n_subsims;
            rk4_advance(&s1, 1, t_i, dt_i, &s0, &delta, &macc, &dat, &b0, &b1);
            s0 = s1;
        }
        truth = s0;

        index_t n_samples_per_frame = 3;
        if (std::uniform_int_distribution<int>(0,120)(*rng) == 1) n_samples_per_frame++;
        Sensor<real_t> *sens;
        std::vector< Measurement<real_t> > obs_list(n_samples_per_frame);

        for (int i = 0; i < n_samples_per_frame; i++) {
            real_t *b = measure_buf + (i * 3);
            switch (i) { //(std::uniform_int_distribution<int>(0,3)(*rng)) {
                case 0:
                    sens = &s_mag;
                    break;
                case 1:
                    sens = &s_acc;
                    break;
                case 2:
                    sens = &s_gyr;
                    break;
                case 3:
                    sens = &s_vel;
                    break;
            }
            Measurement<real_t> obs;
            obs.data = b;
            obs.sensor = sens;

            sens->simulate(b, truth.begin(), rng);

            latest_obs = obs_list[i] = obs;
            measure_available = true;
        }
        
        // update the filter's estimate, given some jittered sensor readings.
        clock_t start = clock();
        filter->advance(obs_list.data(), obs_list.size(), t, dt);
        clock_t end = clock();
        
        // xxx debug
        KinematicState<real_t> &fs = *(KinematicState<real_t>*)filter->x;
        fs.orient = fs.orient.unit();

        t += dt;
        
        return (end-start) / (double)CLOCKS_PER_SEC;
    }
    
};


/********** program state **********/


KinematicState<real_t> get_estimate(real_t *x) {
    KinematicState<real_t> s;
    std::copy(x, x + KINSTATE_SIZE, s.begin());
    return s;
}


class ProgramState : virtual public Animated, 
                     virtual public Drawable,
                     virtual public GUIListener {
public:
    
    SensorSimulator ssim;
    real_t time;
    
    ProgramState():time(0) {}
    
    
    void draw() {
        KinematicState<real_t> &guess = *((KinematicState<real_t>*)ssim.filter->x);
        KinematicState<real_t> truth = ssim.truth;
        
        draw_kstate(truth, vec3(1.));
        draw_kstate(guess, vec3(1,0,0));
        
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
        
        glColor3d(1,0,0);
        VisBall(guess.a, 0.08).draw();
        glColor3d(0,0,1);
        VisBall(guess.omega, 0.08).draw();
        glColor3d(0,1,1);
        VisBall(guess.v, 0.08).draw();
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
        ssim   = SensorSimulator(ssim.rng);
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


void profile_kalman() {
    SensorSimulator ssim;
    double t_acc = 0;
    
    clock_t start = clock();
    for (index_t i = 0; i < 10000; i++) {
        t_acc += ssim.advance_dt(1/60.);
    }
    clock_t end = clock();
    double t_delt = (end-start) / (double)CLOCKS_PER_SEC;
    std::cout << "time: " << t_delt << " acc: " << t_acc << std::endl;
}


int not_main(int argc, char** argv) {
    profile_kalman();
    return 0;
}


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
    
    timer.fps = 60;
    timer.begin();
    win.showAll();
    return 0;
}

