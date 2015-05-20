/* 
 * File:   KalmanFilter.h
 * Author: tbabb
 *
 * Created on February 3, 2015, 6:50 PM
 */

#ifndef KALMANFILTER_H
#define	KALMANFILTER_H

#include <vector>
#include <geomc/function/Dual.h>
#include <geomc/linalg/Matrix.h>

#include "Sensor.h"
#include "Predictor.h"

// todo: can we get the process noise by looking at y~?

// todo: process noise can be a vector of different size.
//       it is an input to f().

/**
 * Estimates the true state of a possibly-multivariate system based on
 * gaussian-uncertain measurements of that state, using an extended Kalman 
 * filter (EKF).
 * 
 * One advantage of this class is that dynamically-varying numbers and types
 * of sensors are tolerated. The state can be updated as data becomes available,
 * and sensors with no new data may be omitted. Similarly, new sensors may be 
 * added to the estimation process as they come online.
 * 
 * This class works with Predictor and Sensor classes to automatically linearize
 * the state about the current estimate, in order to get accurate covariance
 * matrices.
 * 
 * If multiple sensors have noise with mutual covariance, they should be modeled 
 * as a single Sensor object.
 */
template <typename T>
class KalmanFilter {
  public:
    index_t n;               // number of state fields
    T *x;                    // most recent state estimate
    SimpleMatrix<T,0,0> P;   // most recent estimate covariance
    Predictor<T> *predictor;
        
    KalmanFilter(index_t n=1):
            n(n),
            x(new T[n]),
            P(n,n),
            predictor(NULL) {
        std::fill(x, x+n, 0);
    }
    
    KalmanFilter(index_t n, Predictor<T> *predictor):
            n(n),
            x(new T[n]),
            P(n,n),
            predictor(predictor) {
        std::fill(x, x+n, 0);
    }
    
    ~KalmanFilter() { delete [] x; }
    
    // todo: pass around the size for safety/generality.
    // todo: handle control input to prediction
    // todo: factor out intermediate variables into buffers to avoid realloc.
    // todo: split into predict(predictor, dt) and update(measurements)
    // todo: consider a factoring where all sensors are known ahead of time
    //       and are either on or off. 
    // todo: use x directly instead of copying to it (possible after mul fix)
    // todo: IEKF re-linearizes the process (and measurement) about the final 
    //       prediction, and repeats until convergence.
    // todo: some formulations factor the P matrix to maintain positive definite 
    //       symmetry. Can we just manually condition it? max(abs(elem, elem^T))?
    
    // options for pre-alloc: 
    //       * implement mul() with bare arrays
    //       ~ allow SimpleMatrices to be resized if dynamic
    //       - allow SimpleMatrices to accept existing storage/use an allocator
    
    
    // xxx debug {
    template <typename U>
    void print_arr(U *arr, index_t n) {
        for (index_t i = 0; i < n; i++) 
            std::cout << std::setw(12) << std::setprecision(5) << std::right << arr[i] << " ";
        std::cout << "\n";
    }
    // }
    
    
    void predict(T t, T dt) {
        Dual<T> *dual_state = new Dual<T>[n];
        Dual<T> *prediction = new Dual<T>[n];
        SimpleMatrix<T,0,0>      F(n,n); // Jacobian of predictor fn
        SimpleMatrix<T,0,0> tmp_nn(n,n); // temp
        SimpleMatrix<T,0,1>  x_hat(n,1); // predicted state
        
        std::copy(x, x + n, dual_state);
        // x_hat  <- f(x_{k-1}, u_{k-1})  (predict based on prev state and ctrl input)
        // F      <- J[f] dx              (compute jacobian of prediction with respect to state)
        
        // compute prediction and its covariance.
        for (index_t i = 0; i < n; i++) {
            // compute derivative in the ith direction
            dual_state[i].dx = 1;
            
            std::copy(dual_state, dual_state + n, prediction);
            predictor->predict(prediction, dual_state, t, dt);
            // copy derivatives to jacobian
            for (index_t j = 0; j < n; j++) {
                F[j][i] = prediction[j].dx;
            }
            
            dual_state[i].dx = 0;
        }
        
#ifdef DEBUG_KALMAN
        // xxx debug {
            std::cout << "predicted state: \n";
            T *fuckass = new T[n];
            for (index_t i = 0; i < n; i++) { fuckass[i]= prediction[i].x; }
            print_arr(fuckass,n);
            for (index_t i = 0; i < n; i++) { fuckass[i]= dual_state[i].x; }
            std::cout << "previous state: \n";
            print_arr(fuckass,n);
            delete [] fuckass;
        // }
#endif
            
        // xhat no longer needs to be dual
        for (index_t i = 0; i < n; i++) x_hat.begin()[i] = prediction[i].x;
        
        // P_k <- F * P_{k-1} * F^T + Q_k  (compute prediction covariance)
        
        mul(&tmp_nn, F, P);  // tmp <- F * P_{k-1}
        F.transpose();
        mul(&P, tmp_nn, F);  // P_k <- tmp * F^T
        
        // add covariance. We use F's memory as a temp.
        std::fill(F.begin(), F.end(), 0);
        predictor->covariance(F.begin(), x, t, dt);
        
        for (index_t i = 0; i < n * n; i++) P.begin()[i] += F.begin()[i];
        
        std::copy(x_hat.begin(), x_hat.end(), x);
        delete [] dual_state;
        delete [] prediction;
    }
    
    void update(const std::vector< Measurement<T> > &observations) {
        SimpleMatrix<T,0,0> tmp_nn(n,n);
        SimpleMatrix<T,0,1> x_hat(n,1); // predicted state
        Dual<T> *dual_state = new Dual<T>[n];
        std::copy(x, x+n, dual_state);
        std::copy(x, x+n, x_hat.begin());
        
        // how many readings are we dealing with?
        index_t m = 0;
        index_t m_max = 0;
        for (Measurement<T> z_i : observations) {
            index_t m_i = z_i.sensor->reading_size();
            m += m_i;
            m_max = std::max(m_max, m_i);
        }
        
        if (m > 0) {
            // compute predicted sensor readings and covariance
            Dual<T> *dual_z = new Dual<T>[m]; // for holding derivatives
            SimpleMatrix<T,0,0>  H_k(m,n);    // sensor jacobian
            SimpleMatrix<T,0,0> H_kT(n,m);    // transpose of H_k
            SimpleMatrix<T,0,1>    y(m,1);    // measurement residual
            SimpleMatrix<T,0,0>  S_k(m,m);    // residual covariance
            SimpleMatrix<T,0,0>  K_k(n,m);    // Kalman gain
            
            for (index_t j = 0; j < n; j++) {
                // compute the sensor reading vector, and its derivative
                // with respect to the jth state variable
                index_t i = 0;
                dual_state[j].dx = 1;
                for (Measurement<T> z_i : observations) {
                    index_t m_i = z_i.sensor->reading_size();
                    z_i.sensor->measure(dual_z + i, dual_state);
                    i += m_i;
                }
                // copy derivatives to jacobian
                for (index_t a = 0; a < m; a++) { H_k[a][j] = dual_z[a].dx; }
                dual_state[j].dx = 0;
            }

            // y <- z_k - h(x_hat)  (compute measurement residual)
            //                      (h was already invoked above and h(x_hat) remains in dual_z)

            index_t a = 0;
            for (Measurement<T> z_i : observations) {
                // copy measurments into y.
                index_t m_i = z_i.sensor->reading_size();
                std::copy(z_i.data, z_i.data + m_i , y.begin() + a);
                a += m_i;
            }
            
#ifdef DEBUG_KALMAN
            // xxx debug {
            std::cout << "readings: \n";
            print_arr(y.begin(),m);
            T *fucknuts = new T[m];
            for (index_t i = 0; i < m; i++) fucknuts[i] = dual_z[i].x;
            std::cout << "predicted readings: \n";
            print_arr(fucknuts,m);
            delete [] fucknuts;
            std::cout << "H:\n" << H_k;
            // }
#endif
            
            // subtract h(z).
            for (index_t i = 0; i < m; i++) {
                y[i][0] -= dual_z[i].x;
            }
            
            // S_k <- H_k * P_k * H_k^T + R_k  (compute residual covariance)
            
            transpose(&H_kT, H_k);
            mul(&K_k, P, H_kT);    // tmp <- P_k * H_k^T
            mul(&S_k, H_k, K_k);   // S_k <- H_k * tmp
            
            // add sensor covariance.
            // if sensors don't interfere with each other (an assumption),
            // the covariance matrix is block-diagonal.
            T *mtx_tmp = new T[m_max * m_max]; // (alloc)
            a = 0;
            for (Measurement<T> z_i : observations) {
                index_t m_i = z_i.sensor->reading_size();
                z_i.sensor->covariance(mtx_tmp, z_i.data);
                
                index_t j = 0;
                auto    i = S_k.region_begin(MatrixRegion(a, a + m_i));
                auto  end = i.end();
                for (; i != end; i++) {
                    *i += mtx_tmp[j];
                    j++;
                }
                a += m_i;
            }
            delete [] mtx_tmp;

            // K_k <- P_k * H_k^T * S_k^-1
            
#ifdef DEBUG_KALMAN
            // xxx debug {
            std::cout << "R_k:\n" << S_k;
            // }
#endif
            
            
            inv(&S_k, S_k);
            mul(&K_k, H_kT, S_k); // tmp <- H_k^T * S_k^-1 
            K_k = P * K_k;        // K_k <- P_k * tmp (alloc)
            
#ifdef DEBUG_KALMAN
            // xxx debug {
            std::cout << "S_k^-1:\n" << S_k;
            std::cout << "K_k * y:\n" << K_k * y;
            // }
#endif

            // x` <- x_hat + K_k * y  (compute new state estimate)

            SimpleMatrix<T,0,1> x_prime = x_hat + K_k * y; // (alloc x2)
            std::copy(x_prime.begin(), x_prime.end(), x);

            // P <- (I - K_k * H_k) * P_k  (compute new state covariance)

            mul(&tmp_nn, K_k, H_k);
            for (index_t i = 0; i < n*n; i++) tmp_nn.begin()[i] *= -1;
            for (index_t i = 0; i < n; i++) {
                tmp_nn[i][i] += 1;
            }
            mul(&P, tmp_nn, P); // (alloc)
        } else {
            // no sensor update; just use the predicted state.
            // P is already updated.
            std::copy(x_hat.begin(), x_hat.end(), x);
        }

#ifdef DEBUG_KALMAN
        // xxx debug {
        std::cout << "===================\n";
        // }
#endif
    }
    
    void advance(const std::vector< Measurement<T> > &observations, T t, T dt) {
        
        if (predictor) predict(t, dt);
        
        if (observations.size() > 0) {
            T *xx = new T[n];
            for (int i = 0; i < 10; i++) {
                // tah-dah, we are an iterated extended kalman filter.
                std::copy(x, x + n, xx);
                update(observations);
                T sum = 0;
                for (index_t i = 0; i < n; i++) {
                    sum += std::pow(x[i] - xx[i], 2);
                }
                if (std::sqrt(sum) < 0.001) break; // xxx: threshold is arbitrary
            }
            delete [] xx;
        }
    }
    
};


#endif	/* KALMANFILTER_H */

