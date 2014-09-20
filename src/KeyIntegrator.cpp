/* 
 * File:   KeyIntegrator.cpp
 * Author: tbabb
 * 
 * Created on May 9, 2014, 1:00 AM
 */

#include "KeyIntegrator.h"
#include "Timing.h"

KeyIntegrator::KeyIntegrator() : last_check_t(-1), is_down(false), seconds_total(0)  {}
KeyIntegrator::~KeyIntegrator() {}

void KeyIntegrator::keyDown() {
    if (!is_down) {
        last_down_t = now();
        is_down = true;
    }
    if (last_check_t < 0) last_check_t = last_down_t;
}

void KeyIntegrator::keyUp() {
    if (is_down) {
        seconds_total += now() - last_down_t;
        is_down = false;
    }
}

void KeyIntegrator::setDown(bool down) {
    if (down) keyDown(); else keyUp();
}

double KeyIntegrator::dumpTotal() {
    double tot = seconds_total;
    double t = now();
    seconds_total = 0;
    if (is_down) {
        tot += t - last_down_t;
        last_down_t = t;
    }
    last_check_t = t;
    return tot;
}

double KeyIntegrator::getFractionalTime() {
    double t0 = getTimeOfLastCheck();
    double time_down = dumpTotal();
    double t1 = getTimeOfLastCheck();
    double dt = t1 - t0;
    if (dt == 0) return 0;
    else return time_down / (t1 - t0);
}

double KeyIntegrator::getTimeOfLastCheck() const {
    return last_check_t;
}

