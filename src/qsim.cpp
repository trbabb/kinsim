#include "QuantumState.h"
#include "Solver.h"
#include "Image.h"

#include <geomc/function/Raster.h>

#define N 2

using namespace geom;

typedef Vec<index_t,N> ivec;
typedef Vec<float,N>   fvec;
typedef complex<float> cplex;

template <typename T>
Raster<T,T,2,3> to_image(const Grid<complex<T>,2> &g) {
    Vec<index_t,2> dims(g.shape);
    const index_t n = dims[0] * dims[1];
    complex<T> *src = g.data.get();
    T *data = new T[3*n];
    
    T max_val = 0;
    for (index_t i = 0; i < n; i++) {
        T a = src[i].real();
        T b = src[i].imag();
        //data[3*i] = data[3*i+1] = data[3*i+2] = a*a + b*b;
        data[3*i]     = (a + 1) / 2;
        data[3*i + 1] = (b + 1) / 2;
        data[3*i + 2] = 0.5;
        max_val = 1.0;//std::max(max_val, data[3*i]); // 1.0;
    }
    
    if (max_val > 0) {
        for (index_t i = 0; i < n * 3; i++) {
            data[i] /= max_val;
        }
    }
    
    return Raster<T,T,2,3>(dims, data);
}

template<typename T>
Raster<T,T,2,1> to_image(const Grid<T,2> &g) {
    ivec dims(g.shape);
    return Raster<T,T,2,1>(dims, g.data.get());
}

Wavefunction<float,N> make() {
    const index_t n = 512;
    ivec dim(512,512);
    Wavefunction<float,2> psi(dim);
    
    const float freq = 200;
    
    GridIterator<index_t,N> i(ivec(),dim);
    for (; i != i.end(); i++) {
        ivec p  = *i;
        fvec offs = fvec(-0.25,0);
        fvec uv = (2*(fvec)p) / ((fvec)dim) - fvec(1);
        float blob_val = exp(-(uv - offs).mag() * 20) * 5;
        psi.potential[p] = uv.mag();
        psi.state[p] = cplex(blob_val,0) * cplex(std::sin(uv[0]*freq), std::cos(uv[0]*freq));
    }
    return psi;
}

int main(int argc, char **argv) {
    Wavefunction<float,N> wf0 = make();
    Wavefunction<float,N> wf1(ivec(wf0.state.shape));
    Wavefunction<float,N>  b0(ivec(wf0.state.shape));
    Wavefunction<float,N>  b1(ivec(wf0.state.shape));
    
    char name[512];
    
    Raster<float,float,2,3> i0 = to_image(wf0.state);
    save_png(i0, "/Users/tbabb/test/qsim/electron000.png");
    
    const index_t n = 512;
    
    /*
    ivec dim(n,n);
    Grid<float,2> g(dim);
    for (index_t i = 0; i < n; i++) {
        for (index_t j = 0; j < n; j++) {
            ivec p  = ivec(j,i);
            fvec uv = fvec(j/(float)n,i/(float)n) - fvec(0.5,0.5);
            float d = uv.mag();
            float z[4] = {0,0,1,1};
            g[p] = clamp(100* pow(d * 2 - 0.5, 3), 0., 1.);//clamp(interp_cubic(d,z), 0.f, 1.f);
        }
    }
    
    //Grid<Vec<float,2>,2> grad = g.gradient();
    //save_png(Raster<float,float,2,2>(ivec(grad.shape), (const float*)grad.data.get()), "/Users/tbabb/test/qsim/laplace.png");
    save_png(to_image(100.f * g.laplace()), "/Users/tbabb/test/qsim/laplace.png");
     * */
    
    index_t every = 5;
    for (index_t i = 1; i < 1500; i++) {
        rk4_advance(&wf1, 1,   // out state, count
                    0.0, 1.,   // t0, dt
                    &wf0,      // start state
                    &wavefunction_delta<float,N>, // d/dt
                    NULL,      // extra data
                    &b0, &b1); // buffers
        cout << ".";
        cout.flush();
        if (i % every == 0) {
            snprintf(name, 512, "/Users/tbabb/test/qsim/electron%0.3ld.png", i/every);
            Raster<float,float,2,3> img = to_image(wf1.state);
            save_png(img, name);
        }
        std::swap(wf0,wf1);
    }
    cout << endl;
    
}
