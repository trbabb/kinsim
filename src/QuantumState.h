/* 
 * File:   QuantumState.h
 * Author: tbabb
 *
 * Created on November 22, 2014, 1:46 PM
 */

#ifndef QUANTUMSTATE_H
#define	QUANTUMSTATE_H

#include "Grid.h"

#define HBAR 1.0

using namespace std;
using namespace geom;

// todo: how do we handle V() discontinuities?
// todo: edge handling is not terribly robust.
// todo: smarter convolutions
// todo: where can we use duals?
// todo: can we do something better with product calculus?
// todo: can we evolve in the momentum representation?
//       is the operator different?
// todo: what if we evolve it in a lagrangian coordsys (particles)
//       instead of eulerian (fixed grid)? 
//         > interesting, because those particles are very much like actual
//           physical particles.
//         > how do we importance sample? i.e. distribute particles in the most
//           interesting part of the wavefunction. is there a reason why
//           "most interesting" might be mag^2(psi)?
// future: we should not alloc on operations. use a temporary operation
//         object which evals on assignment (or after a certain stack depth)
// future: solve the dirac equation instead.
// question: is complex<dual<>> at play?
//           we know complex<complex<>> is at play (spinor)
//           is complex<complex<dual>> a dirac matrix, by any chance?
// todo: understand the relationship between complex numbers and derivatives better.
//       hint: think of the wavefunction as a field of oscillators.
// todo: understand "continuous derivatives". 
// future: come up with a "taylor series" for product derivatives.
// future: come up with a series method of approximating product_integral(f(x)) 
//         when analytic things are known about f(x), like its nth product derivative
//         or nth difference derivative.
//       > think about "tangent exponentials", i.e., curves of infinitesimal ratio.


template <typename T, index_t N>
struct Wavefunction {
    
    typedef complex<T> B;
    
    Grid<B,N> state;
    Grid<B,N> potential;
    T reduced_mass;
    
    Wavefunction():
            state(Vec<index_t,N>(1)),
            potential(Vec<index_t,N>(1)),
            reduced_mass(1) {}
    
    Wavefunction(Vec<index_t,N> dims):
            state(dims),
            potential(dims),
            reduced_mass(1) {}
    
    Grid<B,N> derivative() const {
        // really, you want the "product derivative" laplacian, since we're going
        // to reverse it with a multiplication, not addition.
        return complex<T>(HBAR,0) * state.laplace() / complex<T>(0, reduced_mass * 2) + 
               potential * state / (complex<T>(0, HBAR));
    }
    
};

template <typename T, index_t N> 
void wavefunction_delta(Wavefunction<T,N> *d_dt, const Wavefunction<T,N> *s0, double t, size_t n, void *data) {
    for (index_t i = 0; i < n; i++) {
        d_dt[i] = s0[i];
        d_dt[i].state = s0[i].derivative();
    }
}


template <typename T, index_t N>
void macc(Wavefunction<T,N> *out, double k, const Wavefunction<T,N> *d_ds, const Wavefunction<T,N> *s, size_t n) {
    for (index_t i = 0; i < n; i++) {
        out[i].state = s[i].state + complex<T>(k,0) * d_ds[i].state;
    }
}

#endif	/* QUANTUMSTATE_H */

