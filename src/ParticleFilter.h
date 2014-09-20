/* 
 * File:   ParticleFilter.h
 * Author: tbabb
 * 
 * A template class for estimating state given imperfect measurements.
 * 
 * Requires implementations of Observer objects which describe the likelihood
 * of a truth value given a measurement, and an implementation of a Predictor
 * object which can take a potential state and propagate it to the next likely
 * state. Predictor must be able to describe the likelihood of its prediction
 * given the position it's advancing from.
 *
 * Created on June 23, 2014, 9:20 PM
 */

#ifndef PARTICLEFILTER_H
#define	PARTICLEFILTER_H

#include <cmath>
#include <random>
#include <vector>

#include "lintypes.h"

typedef std::mt19937_64                        rng_t;
typedef std::normal_distribution<real_t>       d_normal_t; 
typedef std::uniform_real_distribution<real_t> d_uniform_t;

// fwd decl.
template <typename State, typename Measurement>
class Observation;

/**
 * @brief Models a sensor with uncertain relationship to a ground truth state.
 */
template <typename State, typename Measurement>
class Observer {
public:
    std::string name;
    typedef Observation<State, Measurement> observation_t;
    
    Observer() {}
    Observer(std::string name) : name(name) {}
    
    
    virtual real_t likelihood(const State &s, const Measurement &m) = 0;
};


/**
 * Base class for all observations.
 */
template <typename State>
class ObservationBase {
public:
    virtual real_t likelihood(const State &s) = 0;
    virtual std::string name() = 0;
};


/**
 * @brief Models a single reading with a particular type of sensor.
 */
template <typename State, typename Measurement> 
class Observation : public ObservationBase<State> {
public:
    Measurement measurement;
    Observer<State,Measurement> *observer;
    
    Observation(Measurement m, Observer<State,Measurement> *o):
            measurement(m),
            observer(o) {}
    Observation() : observer(NULL) {}
    
    virtual real_t likelihood(const State &s) {
        return observer->likelihood(s, measurement);
    }
    
    virtual std::string name() {
        return observer->name;
    }
};


/**
 * @brief Class for encapsulating a state and local pdf.
 */
template <typename State>
struct Particle {
    State  state;
    real_t wt;
    real_t pdf;
};


/**
 * @brief Class for drawing samples from the next time step of the simulation.
 */
template <typename State>
class Predictor {
    public:
        
    // this function should update the particle pdf with the process pdf
    // the integrator will do the renormalization for you.
    virtual void drawNextStates(std::vector< Particle<State> > &out,
                          const std::vector< Particle<State> > &population,
                          real_t dt,
                          rng_t *rng) = 0;
};


/**
 * @brief Stochastically estimates the most likely state given uncertain observations.
 */
template <typename State> 
class ParticleFilter {
public:
    
    typedef Particle<State>        particle_t;
    typedef ObservationBase<State> observation_t;
    
    std::vector<particle_t> particles;
    Predictor<State>        *predictor;
    rng_t                   *rng;
    
protected:
    std::vector<particle_t> buffer;
    
public:
    
    
    ParticleFilter(index_t n_particles, Predictor<State> *predictor, rng_t *rng):
            predictor(predictor),
            rng(rng) {
        particles.resize(n_particles);
        buffer.resize(n_particles);
        for (index_t i = 0; i < n_particles; i++) {
            particles[i].wt = particles[i].pdf = 1. / n_particles;
        }
    }
    
    
    void setParticleCount(index_t n) {
        buffer.resize(n);
        resample(buffer.data(), particles.data(), particles.size(), n);
        particles.resize(n);
        std::copy(buffer.begin(), buffer.end(), particles.begin());
    }
    
    
    static real_t likelihood(const State &s, const std::vector<observation_t*> &observations) {
        real_t p = 1;
        for (auto i = observations.begin(); i != observations.end(); i++) {
            observation_t *obs = *i;
            real_t k = obs->likelihood(s);
            p *= k;
        }
        return p;
    }
    
    
    void advance(const std::vector<observation_t*> &observations, real_t dt) {
        real_t  wt_tot = 0;
        real_t pdf_tot = 0;
        
        predictor->drawNextStates(buffer, particles, dt, rng);
        
        index_t      n = particles.size();
        particle_t *p0 = particles.data();
        particle_t *p1 = buffer.data();
        for (index_t i = 0; i < n; i++) {
            real_t transition_pdf = likelihood(p1[i].state, observations);
            p1[i].wt  *= transition_pdf;
            p1[i].pdf *= transition_pdf;
            wt_tot    += p1[i].wt;
            pdf_tot   += p1[i].pdf;
        }

        for (index_t i = 0; i < n; i++) {
            p1[i].wt  /= wt_tot;
            p1[i].pdf /= pdf_tot;
        }
        resample(p0, p1, n, n);
    }
    
    
    const State& estimate() {
        // todo: save ptr to this when you update.
        real_t  best_pdf = 0;
        index_t best_idx = 0; 
        for (index_t i = 0; i < particles.size(); i++) {
            if (particles[i].pdf > best_pdf) {
                best_pdf = particles[i].pdf;
                best_idx = i;
            }
        }
        return particles[best_idx].state;
    }
    
    
protected:
    
    void resample(particle_t *dests, particle_t *p0, index_t n_in, index_t n_out) {
        d_uniform_t rand01 = d_uniform_t(0, 1);
        real_t cdf = p0->wt;
        real_t pop_scale = n_in / n_out;
        particle_t *sampled_particle = p0;
        particle_t *end_sample = p0 + n_in;
        
        // use a statified sample: make n cells of uniform size, jitter within the cell
        for (index_t i = 0; i < n_out; i++) {
            real_t smp = (i + rand01(*rng)) / (real_t)n_out;
            while (smp > cdf and sampled_particle < end_sample - 1) {
                sampled_particle++;
                cdf += sampled_particle->wt;
            }
            dests[i]      = *sampled_particle;
            dests[i].wt   = 1.0 / n_out;
            dests[i].pdf *= pop_scale;
        }
        
        if (cdf != cdf or cdf == 0)  {
            std::cout << "WARNING: Zero cdf\n"; // xxx debug
            for (index_t i = 0; i < n_in; i++) {
                p0[i].wt = 1.0 / n_in;
            }
            // try again. all particles are equally likely.
            resample(dests, p0, n_in, n_out);
        }
    }
    
}; // end ParticleFilter class.


#endif	/* PARTICLEFILTER_H */

