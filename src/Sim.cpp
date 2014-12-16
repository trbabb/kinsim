#include <iostream>
#include <OpenGL/gl.h>
#include <algorithm>
#include <cmath>

#include <geomc/random/RandomTools.h>
#include <geomc/random/MTRand.h>
#include <geomc/function/Utils.h>

#include "GLWindow.h"
#include "Shader.h"
#include "AnimTimer.h"
#include "Manipulator.h"
#include "visible/StateUtils.h"
#include "visible/VisBall.h"
#include "visible/VisCallback.h"
#include "visible/VisBox3d.h"

// for Sky
#include "VertexBuffer.h"

#include "RigidBody.h"
#include "Solver.h"
#include "KeyIntegrator.h"

#define GAME_CAMERA 1

#define RESOLUTION 1920

using namespace geom;
using namespace std;

Shader surface;

template <typename T>
Vec<T,3> moment_of_inertia(const Rect<T,3> &box) {
    Vec<T,3> dims = box.getDimensions();
    return Vec<T,3>(
            dims.y * dims.y + dims.z * dims.z,
            dims.x * dims.x + dims.z * dims.z,
            dims.x * dims.x + dims.y * dims.y
        )  / 12.0;
}


// a cone with its tip at P
void cone(Vec3d p, Vec3d axis, double r, int segs) {
    Vec3d perp = axis.cross(Z_AXIS3d);
    double p2 = perp.mag2();
    perp = r * perp / sqrt(p2);
    if (p2 == 0) {
        perp = r * X_AXIS3d;
    }
    
    glBegin(GL_TRIANGLE_FAN);
    glVertex(p);
    for (int i = 0; i < segs + 1; i++) {
        double angle = 2 * M_PI * i / (double)segs;
        Vec3d seg = perp.rotate(axis, angle);
        glVertex(p + axis + seg);
    }
    glEnd();
    
    glBegin(GL_TRIANGLE_FAN);
    glVertex(p + axis);
    
    // copy paste. hooray
    for (int i = 0; i < segs + 1; i++) {
        double angle = 2 * M_PI * i / (double)segs;
        Vec3d seg = perp.rotate(axis, angle);
        glVertex(p + axis + seg);
    }
    glEnd();
}


struct Thruster {
    
    AffineTransform3d xf;
    double size;
    
    void draw(Vec3d torque, double normal_f) {
        Vec3d t_local = xf.applyInverseVector(torque);
        Vec3d r = ZERO_VEC3d / xf;
        const int spans = 10;
        
        glPushMatrix();
        glMultMatrix(xf);
            glScaled(size, size, size);
            glTranslated(0, 0, 0.25);

            glColor3d(1, 0.75, 0.1);
            for (int i = 0; i < 4; i++) {
                int axis =  i & 1;
                int sign = (i & 2) ? -1 : 1;
                Vec3d thrust_dir;
                thrust_dir[axis] = sign;
                double dot = (r ^ thrust_dir).dot(t_local);
                if (dot > 0) {
                    cone(thrust_dir * 0.5, thrust_dir * 2, 0.5, spans);
                }
            }
            
            glColor3d(1,1,1);
            VisBox3d(Rect3d::fromCenter(ZERO_VEC3d, Vec3d(1, 1, 0.5))).draw();
            
            if (normal_f < 0)
                cone(Z_AXIS3d * 0.5, 2* Z_AXIS3d, 1, spans);
        
        glPopMatrix();
    }
    
    inline Vec3d getNormal() {
        return xf.applyVector(Z_AXIS3d);
    }
};


Vec3f shitty_random_temperature(Sampler<float> &smp) {
    Vec3f pts[4] = { Vec3f(1,   0.2, 0),   // red
                     Vec3f(1,  0.85, 0),   // yellow 
                     Vec3f(1,     1, 1),   // white
                     Vec3f(0.5, 0.8, 1) }; // blue
    Vec2f p = smp.box<2>();
    // n-dimensional linear interp
    return interp_linear(p.begin(), pts, 2);
}

void useShader(int dummy) {
    surface.use();
}

void disableShader(int dummy) {
    surface.disable();
}

struct Star {
    Vec3f P;
    Vec3f color;
};


struct Sky : public Drawable, 
             public Animated,
             public GUIListener {
    /*
     * Class which draws a starfield. Stars are drawn out by rotational camera
     * blur. 
     */
    
    // TODO: draw (round?) points to prevent lines from disappearing in stillness.
    
    boost::shared_array<Star> stars;
    int nstars;
    VertexBuffer<Star> vbuf;
    BodyState *state;
    bool inited;
    double dt;
    Shader geo_shader;
    double shutter_angle;
    int w, h;
    
    Sky(int nstars, BodyState *state):
            stars(new Star[nstars]),
            nstars(nstars),
            state(state),
            inited(false),
            dt(1),
            shutter_angle(1),
            w(RESOLUTION),
            h(RESOLUTION) {
        MTRand rng;
        Sampler<float> smp(&rng);
        for (int i = 0; i < nstars; i++) {
            stars[i].P = smp.unit<3>() * 4500;
            stars[i].color = 0.75 * Vec3f(rng.rand<float>(0.1,1)) + 
                             0.25 * shitty_random_temperature(smp);
        }
        geo_shader = shaderFromFile("starfield", "resource/starfield", GL_POINTS, GL_LINE_STRIP, 2);
    }
    
    Quatd getCamOrient() {
        //TODO: obtain camera orientation in non-hard-coded way
        return state->o * Quatd::rotFromAxisAngle(X_AXIS3d, M_PI/2);
    }
    
    void sendUniforms() {
        Quatd q0   = getCamOrient();
        Quatd q1   = std::exp(Quatd(shutter_angle * dt * state->get_omega(), 0)) * q0;
        SimpleMatrix4d m0, m1;
        rotmat(&m0, q0.conj()); // inverse rotation because we're moving "world".
        rotmat(&m1, q1.conj());
        float buf[32];
        std::copy(m0.begin(), m0.end(), buf);
        std::copy(m1.begin(), m1.end(), buf+16);
        glUniformMatrix4fv(geo_shader.getUniformID("rotmat0"), 1, true, buf);
        glUniformMatrix4fv(geo_shader.getUniformID("rotmat1"), 1, true, buf+16);
        glUniform2f(geo_shader.getUniformID("res"), w, h);
    }
    
    void initGL() {
        vbuf.setVertexData(stars, nstars);
        if (nstars < std::numeric_limits<GLushort>::max()) {
            boost::shared_array<GLushort> idx(new GLushort[nstars]);
            for (int i = 0; i < nstars; i++) {
                idx[i] = i;
            }
            vbuf.setIndexData(idx, nstars);
        } else {
            boost::shared_array<GLuint> idx(new GLuint[nstars]);
            for (int i = 0; i < nstars; i++) {
                idx[i] = i;
            }
            vbuf.setIndexData(idx, nstars);
        }
        vbuf.defineAttribute("vertex", offsetof(Star, P));
        vbuf.defineAttribute("color",  offsetof(Star, color), GL_FLOAT, 3);
        
        inited = true;
    }
    
    void update(double t, double dt) {
        this->dt = dt;
    }
    
    void draw() {
        if (!inited) initGL();
        
        vbuf.use();
        geo_shader.use();
        sendUniforms();
        vbuf.bindShader(geo_shader);
        vbuf.drawElements(GL_POINTS, nstars, 0);
        geo_shader.disable();
        glPushMatrix();
        glLoadMatrix(getCamOrient().conj()); // no translation kthx
        vbuf.drawElements(GL_POINTS, nstars, 0);
        glPopMatrix();
        vbuf.disable();
    }
    
    bool windowReshaped(GLWindow *win, int w_new, int h_new) {
        w = w_new; h = h_new;
        return false;
    }
    
};


struct Top : virtual public Animated,
             virtual public Drawable,
             virtual public ForceSystem {
    
    // cm is at body space origin.
    
    BodyState s;
    Rect3d box;
    
    Top():
            box(Rect3d::fromCenter(ZERO_VEC3d, Vec3d(1, 1, 1/3.0))) {
        const Vec3d axis = Vec3d(0.15,0,1).unit();
        s.o = Quatd::rotDirectionAlign(Z_AXIS3d, axis);
        s.m = 1;
        s.J = s.m * moment_of_inertia(box);
        s.x = Vec3d(0, 0, 1/3.);
        s.L = s.o * (Z_AXIS3d * 0);
    }
    
    // sadly, this ain't work
    void sumForces(const BodyState &state, double t, Vec3d *f, Vec3d *tau) {
        *f   = 0.0;
        *tau = 0.0;
        
        const Vec3d tip_body(0, 0, -1/3.);
        const Vec3d g_in(0, 0, -9.80);
        Vec3d centr_accel = state.get_centripetal_accel_point_body(tip_body);
        Vec3d tip_accel = g_in + centr_accel; // needs also angular acceleration.
        //Vec3d tip_v = state.get_velocity_point_body(tip_body);

        // mg downward on cm.
        state.accumulate_force_inertial(f, tau, 
                                        state.x, 
                                        state.m * g_in);
        // oppose tip acceleration.
        state.accumulate_force_inertial(f, tau, 
                                        tip_body, 
                                        state.m * -tip_accel);
    }
    
    void applyConstraints(BodyState *state, index_t n, double t, double dt) {
        // do nothing.
    }
    
    void init(double t) {}
    
    void update(double t, double dt) {
        const int n_steps = 60;
        double dtt = dt / n_steps;
        
        for (int n = 0; n < n_steps; n++) {
            BodyState b0, b1;
            BodyState s1;

            // advance the simulation
            rk4_advance(&s1, 1,               // one output state
                        t + n * dtt, dtt,     // times
                        &s,                   // current state
                        body_delta_for_rk4,   // delta function
                        body_combine_mul_exp, // accumulation function
                        (ForceSystem*)this,   // extra data (convert the ptr now, while you still have type info)
                        &b0, &b1);            // working buffers
            s = s1;
        }
    }
    
    void draw() {
        glPushMatrix();
        
            glTranslate(s.x);
            glRotate(s.o);
            VisBox3d(box, true).draw();
            glBegin(GL_LINES);
                glVertex3d(0, 0, -1/3.0);
                glVertex3d(0, 0,  1/3.0);
            glEnd();
            
        glPopMatrix();
    }
    
};


class ProgramState : virtual public Animated, 
                     virtual public Drawable, 
                     virtual public GUIListener,
                     virtual public ForceSystem {
public:
    
    ProgramState(Camera *cam, Top *top) : 
        cam(cam),
        box(ZERO_VEC3d, Vec3d(3,4,2).unit() * 2), // {}
        top(top) {}
    
    Camera *cam;
    double t0;
    BodyState s;
    Rect3d box;
    Thruster thrusters[6];
    
    Top *top;
    
    // thruster state: [xyz][+/-]
    KeyIntegrator thrust_attitude_state[3][2];
    KeyIntegrator thrust_velocity_state[3][2];
    
    Vec3d F; // force
    Vec3d T; // torque
    
    void init(double t) {
        t0 = t;
        
        box.setCenter(ZERO_VEC3d);
        
        s.o = Quatd(ZERO_VEC3d, 1);
        s.m = 1;
        s.J = s.m * moment_of_inertia(box);
        
        Vec3d dims = box.getDimensions();
        int i = 0;
        for (int axis = 0; axis < 3; axis++){
            for (int sign = -1; sign <= 1; sign += 2, i++) {
                Vec3d N;
                N[axis] = sign;
                thrusters[i].xf = translation(N * dims[axis] * 0.5) * direction_align(Z_AXIS3d, N);
                thrusters[i].size = 0.125;
            }
        }
        
        surface = shaderFromFile("diffuse", "resource/diffuse");
        
        /*
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 2; j++) {
                thrust_attitude_state[i][j] = false;
                thrust_velocity_state[i][j] = false;
            }
        }
         */
    }
    
    void update(double t, double dt) {
        BodyState b0, b1;
        BodyState s1;
        
        updateKeys();
        
        // advance the simulation
        rk4_advance(&s1, 1,             // one output state
                    t, dt,              // times
                    &s,                 // current state
                    body_delta_for_rk4, // delta function
                    body_combine_mul_exp, // accum function
                    (ForceSystem*)this, // extra data (convert the ptr now, while you still have type info)
                    &b0, &b1);          // working buffers
        
        s = s1;
        
#if GAME_CAMERA
        setCamera();
#endif
    }
    
    void applyConstraints(BodyState *state, index_t n, double t, double dt) {
        // do nothing.
    }
    
    void setCamera() {
        Vec3d cam_pos = Vec3d(0,-1.975,0.5) * 3.5;
        //Vec3d v_body = s.o * s.v.unit();
        //Vec3d a = cam_pos.unit() ^ -v_body;
        //Quatd v_o = Quatd::rotFromAxisAngle(a.unit(), cam_pos.angleTo(-v_body) / 3);
        Quatd v_o = Quatd(0,0,0,1);
        cam->setPosition(s.x + s.o * v_o * cam_pos);
        cam->setDirection(s.o * v_o * Y_AXIS3d);
        cam->setUp(s.o * v_o * Z_AXIS3d);
    }
    
    void updateKeys() {
        // update thrust
        for (int i = 0; i < 3; i++) {
            double thr_pos = thrust_velocity_state[i][0].getFractionalTime();
            double thr_neg = thrust_velocity_state[i][1].getFractionalTime();
            F[i] = thr_pos - thr_neg;
        }
        // update torque
        for (int i = 0; i < 3; i++) {
            double thr_pos = thrust_attitude_state[i][0].getFractionalTime();
            double thr_neg = thrust_attitude_state[i][1].getFractionalTime();
            T[i] = thr_pos - thr_neg;
        }
    }
    
    inline Vec3d getThrust() {
        return F;
    }
    
    inline Vec3d getTorque() {
        return T;
    }
    
    void sumForces(const BodyState &state, double t, Vec3d *f, Vec3d *tau) {
        *f   = 0.0;
        *tau = 0.0;
        
        const double amt = 3.0;
        
        state.accumulate_torque_body(tau, amt * getTorque());
        state.accumulate_force_body(f, tau, ZERO_VEC3d, 4 * amt * getThrust());
    }
    
    void draw() {
        glPushMatrix();
            // box
            glTranslate(s.x);
            glPushMatrix();
                glRotate(s.o);
                glColor(ONE_VEC3d);
                VisBox3d(box, true).draw();
                
                for (int i = 0; i < 6; i++) {
                    thrusters[i].draw(getTorque(), thrusters[i].getNormal().dot(getThrust()));
                }
            glPopMatrix();
        glPopMatrix();
    }
    
    bool keyEvent(GLWindow* window, unsigned char key, int keycode, bool down, bool special, int x, int y) {
        //todo: event timing is very poor; probably due to glut implementation. perhaps SDL is better?
        int t_idx = -1;
        int t_neg = 0;
        
        if (key == 'a' or key == 'd') {
            t_idx = 2;
            t_neg = key == 'd';
        } else if (key == 'w' or key == 's') {
            t_idx = 0;
            t_neg = key == 's';
        } else if (key == 'q' or key == 'e') {
            t_idx = 1;
            t_neg = key == 'q';
        }
        
        int x_idx = -1;
        int x_neg;
        if (key == 'j' or key == 'l') {
            x_idx = 0;
            x_neg = key == 'j';
        } else if (key == 'i' or key == 'k') {
            x_idx = 1;
            x_neg = key == 'k';
        } else if (key == 'u' or key == 'h') {
            x_idx = 2;
            x_neg = key == 'h';
        } else if (key == ' ') {
            *top = Top();
        }
        
        if (t_idx != -1) {
            thrust_attitude_state[t_idx][t_neg].setDown(down);
        }
        
        if (x_idx != -1) {
            thrust_velocity_state[x_idx][x_neg].setDown(down);
        }
        
        return (x_idx != -1 or t_idx != -1);
    } 
    
};



void shaderParams(ProgramState *s) {
    glUniform1f(surface.getUniformID("Kd"), 1.0);
    glUniform1f(surface.getUniformID("Ks"), 0.0);
    glUniform1f(surface.getUniformID("gloss"), 100);
    
    glUniform3f(surface.getUniformID("Lcolor1"), 1,1,1);
    glUniform3f(surface.getUniformID("Lpos1"), 1,1,1);
    
    glUniform3f(surface.getUniformID("Lcolor2"), 0.18, 0.45, 0.80);
    glUniform3f(surface.getUniformID("Lpos2"), -1,-1,-1);
    
    float xf_buf[16];
    AffineTransform3d xf = s->cam->getMatrix();
    std::copy(xf.mat.begin(), xf.mat.end(), xf_buf);
    glUniformMatrix4fv(surface.getUniformID("w2cam"),   1, true,  xf_buf);
}


int main(int argc, char *argv[]) {
    GLWindow win(&argc, argv, "noodle", RESOLUTION, RESOLUTION);
    AnimTimer timer(&win);
    
    Camera &cam = win.cam;
    cam.setFov(0.4 * M_PI);
    cam.setFar(5000);
    cam.setNear(0.1);
    cam.setPosition(-X_AXIS3d * 5);
    cam.setUp(Z_AXIS3d);
    cam.setCenterOfInterest(ZERO_VEC3d);
    
#if !GAME_CAMERA
    Manipulator manip(&win, &win.cam);
#endif
    
    Top top;
    ProgramState s(&win.cam, &top);
    Sky sky(5000, &s.s);
    
    win.scene.push_back(&s);
    win.scene.push_back(&top);
    timer.anims.push_back(&s);
    timer.anims.push_back(&sky);
    timer.anims.push_back(&top);
    win.guiListeners.push_back(&s);
    win.scene.push_back(&sky);
    
    win.scene.push_back(new VisCallback<Vec4d>(&glDoSetColor4d, 2* Vec4d(0.125,1)));
    win.scene.push_back(new VisCallback<int>(&useShader, 1));
    win.scene.push_back(new VisCallback<ProgramState*>(&shaderParams, &s));
    Sampler<double> rv;
    for (int i = 0; i < 250; i++){
    	Vec3d c = rv.solidball<3>(300);
        VisBall *b = new VisBall(c, rv.rng->rand<double>(0.1, 18), 24, 12);
        win.scene.push_back(b);
    }
    win.scene.push_back(new VisCallback<int>(&disableShader, 1));
    
    win.setClearColor(Vec3d(0.01));
    
    timer.begin();
    win.showAll();
    return 0;
    
}