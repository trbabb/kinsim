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
#include "KalmanBuffer.h"

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
    index_t  n;              // number of state fields
    T*       x;              // most recent state estimate
    SimpleMatrix<T,0,0> P;   // most recent estimate covariance
    Predictor<T> *predictor; // estimator of "next" state
    KalmanBuffer<T> pool;    // buffer allocator
        
    KalmanFilter(index_t n=1, index_t m=1):
            n(n),
            x(new T[n]),
            P(n,n),
            predictor(NULL),
            pool(n,m) {
        std::fill(x, x+n, 0);
    }
    
    KalmanFilter(index_t n, index_t m, Predictor<T> *predictor):
            n(n),
            x(new T[n]),
            P(n,n),
            predictor(predictor),
            pool(n,m) {
        std::fill(x, x+n, 0);
    }
    
    ~KalmanFilter() { delete [] x; }
    
    // todo: pass around the size for safety/generality.
    // todo: handle control input to prediction
    //    x: factor out intermediate variables into buffers to avoid realloc.
    //    x: split into predict(predictor, dt) and update(measurements)
    // todo: consider a factoring where all sensors are known ahead of time
    //       and are either on or off. 
    // todo: use x directly instead of copying to it (possible after mul fix)
    //    x: IEKF re-linearizes the process (and measurement) about the final 
    //       prediction, and repeats until convergence.
    // todo: some formulations factor the P matrix to maintain positive definite 
    //       symmetry. Can we just manually condition it? abs(max(elem, elem^T))?
    
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
        Dual<T>* x0 = pool.getX0(); // dual version of x
        Dual<T>* x1 = pool.getX1(); // predicted x
        WrapperMatrix<T,0,0>      F = pool.getPredictionJacobian();     // Jacobian of predictor fn
        WrapperMatrix<T,0,0> tmp_nn = pool.getPredictionMatrixBuffer(); // a temp
        
        std::copy(x, x + n, x0);
        // x_hat  <- f(x_{k-1}, u_{k-1})  (predict based on prev state and ctrl input)
        // F      <- J[f] dx              (compute jacobian of prediction with respect to state)
        
        // compute prediction and its covariance.
        for (index_t i = 0; i < n; i++) {
            // compute derivative in the ith direction
            x0[i].dx = 1;
            
            std::copy(x0, x0 + n, x1); // predictor might only update some variables.
            predictor->predict(x1, x0, t, dt);
            // copy derivatives to jacobian
            for (index_t j = 0; j < n; j++) {
                F[j][i] = x1[j].dx;
            }
            
            x0[i].dx = 0;
        }
        
#ifdef DEBUG_KALMAN
        // xxx debug {
            std::cout << "predicted state: \n";
            T *fuckass = new T[n];
            for (index_t i = 0; i < n; i++) { fuckass[i]= x1[i].x; }
            print_arr(fuckass,n);
            for (index_t i = 0; i < n; i++) { fuckass[i]= x0[i].x; }
            std::cout << "previous state: \n";
            print_arr(fuckass,n);
            delete [] fuckass;
        // }
#endif
         
        // P_k <- F * P_{k-1} * F^T + Q_k  (compute prediction covariance)
        
        mul(&tmp_nn, F, P);  // tmp <- F * P_{k-1}
        F.transpose();
        mul(&P, tmp_nn, F);  // P_k <- tmp * F^T
        
        // add covariance. We use F's memory as a temp.
        std::fill(F.begin(), F.end(), 0);
        predictor->covariance(F.begin(), x, t, dt);
        
        for (index_t i = 0; i < n * n; i++) P.begin()[i] += F.begin()[i];
        
        // xhat no longer needs to be dual
        for (index_t i = 0; i < n; i++) x[i] = x1[i].x;
    }
    
    
    void update(const std::vector< Measurement<T> > &observations) {
        // how many readings are we dealing with?
        index_t m = 0;
        index_t m_max = 0;
        for (Measurement<T> z_i : observations) {
            index_t m_i = z_i.sensor->reading_size();
            m += m_i;
            m_max = std::max(m_max, m_i);
        }
        
        if (m == 0) return; // no readings to use.
        
        // compute predicted sensor readings and covariance
        // todo: check initialization. who needs to be identity?
        Dual<T>* z = pool.getZ(); // for holding derivatives                 (length m)
        WrapperMatrix<T,0,0>    H_k = pool.getHk(m);     // sensor jacobian      (m x n)
        WrapperMatrix<T,0,0>   H_kT = pool.getHkT(m);    // transpose of H_k     (n x m)
        WrapperMatrix<T,0,1>      y = pool.getY(m);      // measurement residual (m x 1)
        WrapperMatrix<T,0,0>    S_k = pool.getSk(m);     // residual covariance  (m x m)
        WrapperMatrix<T,0,0>    K_k = pool.getKk(m);     // Kalman gain          (n x m)
        WrapperMatrix<T,0,0> tmp_nn = pool.getTempNn();  // Temporary mtx        (n x n)
        WrapperMatrix<T,0,1>  x_hat(x, n, 1);            // predicted state      (n x 1) (wraps x)
        
        // compute the sensor reading vector, and its derivative
        // with respect to the jth state variable. here we are linearizing
        // the readings about the measurement.
        Dual<T>* x0 = pool.getX0(); // state as a dual vector
        std::copy(x, x+n, x0);
        for (index_t j = 0; j < n; j++) {
            index_t i = 0;
            x0[j].dx  = 1; // query derivative in jth direction
            for (Measurement<T> z_i : observations) {
                index_t m_i = z_i.sensor->reading_size();
                z_i.sensor->measure(z + i, x0);
                i += m_i;
            }
            // copy derivatives to jacobian
            for (index_t a = 0; a < m; a++) { H_k[a][j] = z[a].dx; }
            x0[j].dx = 0;
        }

        // y <- z_k - h(x_hat)  (compute measurement residual)
        //                      (h was already invoked above and h(x_hat) remains in z)
        
        // copy measurments into y.
        index_t a = 0;
        for (Measurement<T> z_i : observations) {
            index_t m_i = z_i.sensor->reading_size();
            std::copy(z_i.data, z_i.data + m_i , y.begin() + a);
            a += m_i;
        }
        
#ifdef DEBUG_KALMAN
        // xxx debug {
        std::cout << "readings: \n";
        print_arr(y.begin(),m);
        T *fucknuts = new T[m];
        for (index_t i = 0; i < m; i++) fucknuts[i] = z[i].x;
        std::cout << "predicted readings: \n";
        print_arr(fucknuts,m);
        delete [] fucknuts;
        std::cout << "H:\n" << H_k;
        // }
#endif
        
        // subtract h(z).
        for (index_t i = 0; i < m; i++) {
            y[i][0] -= z[i].x;
        }
        
        // S_k <- H_k * P_k * H_k^T + R_k  (compute residual covariance)
        
        WrapperMatrix<T,0,0> tmp_nm = pool.getTempPHkT(m);
        transpose(&H_kT, H_k);
        mul(&tmp_nm, P, H_kT);    // tmp <- P_k * H_k^T
        mul(&S_k, H_k, tmp_nm);   // S_k <- H_k * tmp
        
        // xxx: todo: this. needs tweaking of sensors also.
        // add sensor noise covariance.
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
        
        WrapperMatrix<T,0,0> HkTSk_tmp = pool.getTempHkTSk(m); // (n x m)
        inv(&S_k, S_k);
        mul(&HkTSk_tmp, H_kT, S_k); // tmp <- H_k^T * S_k^-1 
        mul(&K_k, P, HkTSk_tmp);    // K_k <- P_k * tmp
        
#ifdef DEBUG_KALMAN
        // xxx debug {
        std::cout << "S_k^-1:\n" << S_k;
        std::cout << "K_k * y:\n" << K_k * y;
        // }
#endif

        // x` <- x_hat + K_k * y  (compute new state estimate)
        
        WrapperMatrix<T,0,1> x_prime = pool.getTempX();
        mul(&x_prime, K_k, y);
        add(&x_hat, x_hat, x_prime);

        // P <- (I - K_k * H_k) * P_k  (compute new state covariance)

        WrapperMatrix<T,0,0> tmp_kk = pool.getTempKk();
        mul(&tmp_nn, K_k, H_k);
        for (index_t i = 0; i < n*n; i++) tmp_nn.begin()[i] *= -1;
        for (index_t i = 0; i < n; i++) {
            tmp_nn[i][i] += 1;
        }
        mul(&tmp_kk, tmp_nn, P);
        mtxcopy(&P, tmp_kk);

#ifdef DEBUG_KALMAN
        // xxx debug {
        std::cout << "===================\n";
        // }
#endif
        
        // inventory: 
        // z of size m
        // y of size m
        // K_k, H_k and H_kT of size m*n
        // tmp_nn of size n*n
        // S_k which is m*m
        // one buffer of size (largest sensor output)^2 to receive sensor noise covariance
        // a temporary allocation to compute P * K_k (n x m)   // line 258
        // a temporary to compute tmp_nn * P. (n x n)          // line 279
        // x_hat can alias x.
        // temporaries to compute x + K_k * y; destination = x again. (n x 1) // line 269
    }
    
    void advance(const std::vector< Measurement<T> > &observations, T t, T dt, index_t iters=5) {
        if (predictor) predict(t, dt);
        
        if (observations.size() > 0) {
            T *xx = new T[n];
            for (int i = 0; i < iters; i++) {
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

