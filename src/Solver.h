/* 
 * File:   Solver.h
 * Author: tbabb
 *
 * Created on October 29, 2013, 8:44 PM
 */

#ifndef SOLVER_H
#define	SOLVER_H

#include <cstddef>

/**
 * Multiply / Accumulate.
 * Function for computing `m * x + b` over arrays `x` and `b`.
 * @param out Output array of `T`.
 * @param m Gain. 
 * @param x Array of `T` to be scaled.
 * @param b Bias array of `T`.
 * @param n Number of items.
 */
template <typename T>
void macc(T *out, T m, const T *x, const T *b, size_t n) {
    T *end = out + n;
    for (; out != end; out++, x++, b++) {
        *out = m * (*x) + (*b);
    }
}

// note that d_dt(state) is not necessarily the same shape as state.
// also, some of state may not be integrated. it would be better if dSdt had
// its own type, and we allowed for direct computation of intermediate and
// state variables.


/**
 * Advance a state with 4th-order Runge-Kutta.
 * 
 * @param s_final Output state.
 * @param n Number of parallel states.
 * @param t0 Time from which to advance.
 * @param dt Length of time to advance.
 * @param s0 State at time `t0`.
 * @param delta Function accepting a state `s`, a time `t`, an item count, and 
 *        a pointer to arbitrary external data, generating a rate of change for 
 *        each state variable at `t`.
 * @param buf0 Optional pre-allocated buffer for intermediate work, of length `n`.
 * @param buf1 Optional pre-allocated buffer, distinct from `buf0`, of length `n`.
 */
template <typename T, typename State>
void rk4_advance(State *s_final, size_t n, 
                 T t0, T dt, 
                 const State *s0, 
                 void (*delta)(State*, const State*, T, size_t, void *), 
                 void *external_data=NULL, 
                 State *buf0=NULL, 
                 State *buf1=NULL) {
    bool alloc_buf = (buf0 == NULL or buf1 == NULL);
    if (alloc_buf) {
        buf0 = new State[n*2];
        buf1 = buf0 + n;
    }
    
    T t1 = t0 + dt / 2.0;
    T t2 = t0 + dt;
    
    // k1 <- dSdt(t0, s0):
    delta(buf0, s0, t0, n, external_data);   // buf0 = k1
    
    // k2 <- dSdt(t0 + dt/2, s0 + k1 * dt/2):
    macc(s_final, dt/6, buf0, s0, n);        // s_final <- dt/6 * buf0 + s0
    macc(buf0, dt/2, buf0, s0, n);           // buf0    <- dt/2 * buf0 + s0
    delta(buf1, buf0, t1, n, external_data); // buf1 = k2
    
    // k3 <- dSdt(t0 + dt/2, s0 + k2 * dt / 2):
    macc(s_final, dt/3, buf1, s_final, n);   // s_final <- dt/3 * buf1 + s_final
    macc(buf1, dt/2, buf1, s0, n);           // buf1    <- dt/2 * buf1 + s0
    delta(buf0, buf1, t1, n, external_data); // buf0 = k3
    
    // k4 <- dSdt(t0 + dt, s0 + k3 * dt):
    macc(s_final, dt/3, buf0, s_final, n);   // s_final <- dt/3 * buf0 + s_final
    macc(buf0, dt, buf0, s0, n);             // buf0    <- dt   * buf0 + s0
    delta(buf1, buf0, t2, n, external_data); // buf1 = k4
    
    macc(s_final, dt/6, buf1, s_final, n);   // s_final <- dt/6 * buf1 + s_final
    
    if (alloc_buf) delete [] buf0;
}


#endif	/* SOLVER_H */

