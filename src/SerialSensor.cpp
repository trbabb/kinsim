#include <OpenGL/gl.h>
#include <algorithm>
#include <geomc/function/Dual.h>
#include <string>
#include <list>

// io
#include <fcntl.h>
#include <unistd.h>
#include <sys/signal.h>
#include <sys/ioctl.h>
#include <termios.h>
#include <vector>

#include "Sensor.h"
#include "KalmanFilter.h"
#include "KinematicState.h"

#include "GLWindow.h"
#include "GUIListener.h"
#include "AnimTimer.h"
#include "Manipulator.h"
#include "glHelpers.h"
#include "visible/VisBox3d.h"

#define SID_ACC   0
#define SID_GYR   1
#define SID_MAG   2
#define SID_VEL   3
#define N_SENSORS 4


#define OMEGA_PROCESS_VARIANCE   0.25
#define ACCEL_PROCESS_VARIANCE   0.025

#define MEASUREMENT_COUNT 12

using namespace geom;
using namespace std;

// todo: allow for readings to report their own variance.
// todo: sensor dynamic range.
// todo: matrix pre-allocation.
// todo: binary message passing / sensor abstraction.
// todo: accelerometer filtering. 
//       > not available on your board.

void draw_kstate(const KinematicState<real_t> &ks, vec3 color) {
    glPushMatrix();
        //glTranslate(ks.x);
        glRotate(ks.orient);
        glColor(color);
        VisBox3d(Rect3d::fromCenter(ZERO_VEC3d, Vec3d(3,4,1))).draw_wireframe();
        //VisAxis().draw();
    glPopMatrix();
    
    glBegin(GL_LINES);
        glColor3d(1,0,0);
        glVertex3d(0,0,0);
        glVertex(ks.a);
        
        glColor3d(0,0,1);
        glVertex3d(0,0,0);
        glVertex(ks.omega);
        
        glColor3d(0,1,1);
        glVertex3d(0,0,0);
        glVertex(ks.v);
    glEnd();
}


class SerialSensor : public Drawable, public Animated {
public:
    KalmanFilter<real_t>   filter;
    KinematicState<real_t> state;
    
    // time:
    real_t last_update_t;
    
    // readings:
    vec3 readings[N_SENSORS];
    bool available[N_SENSORS];
    
    // sensor models:
    SensorAccelerometer<real_t>  s_acc;
    SensorRateGyro<real_t>       s_gyr;
    SensorMagnetometer<real_t>   s_mag;
    SensorGPSLocation<real_t>    s_loc;
    
    // incoming serial data:
    list<string> tokens;
    int serial_fd;
    bool serial_ok;
    string tok_buf;
    
    SerialSensor():
            filter(KINSTATE_SIZE,
                   MEASUREMENT_COUNT,
                   new KinematicPredictor<real_t>(
                        ACCEL_PROCESS_VARIANCE,
                        OMEGA_PROCESS_VARIANCE)),
            last_update_t(-1) {
        std::fill(available, available + 4, false);
        serial_ok = openSerial("/dev/tty.usbmodem1411", B115200);
        init_sensors();
        std::copy(state.begin(), state.end(), filter.x);
    }
    
    void init_sensors() {
        // problem: all this calibration data is meaningless when
        // we change sensor scale. 
        // solution: Actually, all you need to do is change the state2reading.
        // if the noise is in digit space and not state space, you're fine
        // ...but your bias is probably in world units. :(
        
        s_mag.set_variance(vec3(2.35506, 2.36629, 2.20033));
        s_acc.set_variance(vec3(15814.7, 13691.8, 26806.5));
        s_gyr.set_variance(vec3(665.779,  270.91, 468.631)); 
        s_loc.set_variance(vec3(25)); //  stddev of 5m.
        
        const real_t gravity = 9.80;
        const double gyr_deg_per_digit = 8.75 / 1000;
        const real_t mag_f_strength = 0.486; // gauss
        
        AffineTransform<real_t,3> xf;
        xf = scale(vec3(0.59411 * mag_f_strength / 300)) * translation(vec3( 94.2309, 143.86, 7.6891));
        s_mag.state2reading = xf.inverse();
        xf = scale(vec3(gravity / 16400.))               * translation(vec3(-0.396003, 0.154401, 0.108385));
        s_acc.state2reading = xf.inverse();
        xf = scale(vec3(M_PI * gyr_deg_per_digit / 180)) * translation(vec3(65.5869, 38.9592, 34.1288));
        s_gyr.state2reading = xf.inverse();
    }
    
    void draw() {
        draw_kstate(state, vec3(1));
    }
    
    void update(double t, double dt) {
        consumeSerial();
        while (tokens.size() > 0) {
            // we gots data, my friends.
            string tok = tokens.front();
            tokens.pop_front();
            vec3 v;
            char c;
            int redf = sscanf(tok.c_str(), "%c: %lf %lf %lf", &c, &v.x, &v.y, &v.z);
            int idx = -1;
            
            if (redf == 4) {
                switch (c) {
                    case 'a':
                        idx = SID_ACC; break;
                    case 'g':
                        idx = SID_GYR; break;
                    case 'm':
                        idx = SID_MAG; break;
                }
            } else {
                std::cout << "couldn't read " << tok << "\n";
            }
            
            if (idx == SID_MAG and v == readings[SID_MAG]) continue; // stale reading. 
            
            if (available[idx]) {
                // collision. Time to flush the sensor readings.
                flushSensors(t);
            }
            
            available[idx] = true;
            readings[idx]  = v;
        }
    }
    
    void flushSensors(double t) {
        Sensor<real_t> *sensors[N_SENSORS] = {&s_acc, &s_gyr, &s_mag, &s_loc};
        double dt = t - last_update_t;
        if (last_update_t < 0) {
            // first run through.
            dt = 0.1;
            last_update_t = t - dt;
        }
        std::vector< Measurement<real_t> > measurements;
        for (int i = 0; i < N_SENSORS; i++) {
            if (available[i]) {
                Measurement<real_t> m;
                m.data   = (readings + i)->begin();
                m.sensor = sensors[i];
                measurements.push_back(m);
            }
        }
        vec3 loc(0.0);
        Measurement<real_t> loc_cnstr;
        loc_cnstr.data   = loc.begin();
        loc_cnstr.sensor = &s_loc;
        measurements.push_back(loc_cnstr);
        
        filter.advance(measurements.data(), measurements.size(), last_update_t, dt);
        state = *((KinematicState<real_t>*)filter.x);
        last_update_t = t;
    }
    
    bool openSerial(const char *port, speed_t baud) {
        serial_fd = open(port, O_RDWR | O_NOCTTY | O_NDELAY);
    
        if (serial_fd < 0) {
            return false;
        }

        //block on read
        fcntl(serial_fd, F_SETFL, O_NDELAY);

        struct termios settings;
        tcgetattr(serial_fd, &settings);
        cfsetispeed(&settings, baud);
        cfsetospeed(&settings, baud);
        settings.c_cflag &= ~CSIZE;
        settings.c_cflag |=  CS8;    // 8bit chars
        settings.c_cflag &= ~PARENB; // no parity
        settings.c_cflag &= ~CSTOPB; // one stop bit

        tcsetattr(serial_fd, TCSANOW, &settings);
        return true;
    }
    
    // it's breakfast time, bitches
    bool consumeSerial() {
        if (serial_fd < 0) {
            return false;
        }
        
        // take whatever chunk is in the serial buffer.
        // break off a token if there's a complete one built up.
        char buf[128];
        int bytes_avail;
        ioctl(serial_fd, FIONREAD, &bytes_avail);
        bool found = false;
        while (bytes_avail > 0) {
            size_t red = read(serial_fd, &buf, 127);
            if (red == 0) return found;

            int start = 0;
            for (size_t i = 0; i < red; i++) {
                if (buf[i] == '\n') {
                    tok_buf.append(string(&buf[start], i - start + 1));
                    tokens.push_back(tok_buf);

                    found = true;
                    tok_buf = "";
                    start = i+1;
                }
            }
            string remaining = string(&buf[start], red-start);
            tok_buf += remaining;

            ioctl(serial_fd, FIONREAD, &bytes_avail);
        }
        return found;
    }
    
    
};


int main(int argc, char **argv) {
    GLWindow win(&argc, argv, "kinematic sensor", 1280, 1024);
    AnimTimer timer(&win);
    
    Camera &cam = win.cam;
    cam.setPosition(-X_AXIS3d * 15);
    cam.setCenterOfInterest(ZERO_VEC3d);
    cam.setUp(Z_AXIS3d);
    cam.setFar(250);
    cam.setNear(0.1);
    
    Manipulator manip(&win, &win.cam);
    
    SerialSensor s;
    
    win.scene.push_back(&s);
    timer.anims.push_back(&s);
    //win.guiListeners.push_back(&s);
    
    timer.fps = 60;
    timer.begin();
    win.showAll();
}