/* 
 * File:   KeyIntegrator.h
 * Author: tbabb
 *
 * Created on May 9, 2014, 1:00 AM
 */

#ifndef KEYINTEGRATOR_H
#define	KEYINTEGRATOR_H

struct KeyIntegrator {
public:
    KeyIntegrator();
    virtual ~KeyIntegrator();
    
    void setDown(bool down);
    void keyDown();
    void keyUp();
    double dumpTotal();
    double getTimeOfLastCheck() const;
    double getFractionalTime();
    
    double last_down_t;
    double last_check_t;
    bool is_down;
    
private:
    
    double seconds_total;
};

#endif	/* KEYINTEGRATOR_H */

