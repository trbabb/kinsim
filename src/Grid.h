/* 
 * File:   Grid.h
 * Author: tbabb
 *
 * Created on November 22, 2014, 6:29 PM
 */

#ifndef GRID_H
#define	GRID_H

#include <complex>
#include <algorithm>
#include <boost/shared_array.hpp>

#include <geomc/linalg/Vec.h>

#define GEOMC_NONGRID_TYPE(U,V) \
typename boost::enable_if< \
       boost::is_base_of< Grid<typename U::elem_t, U::dim>, U >, \
       V >::type

// what if instead of a grid, we had particles? i.e. lagrangian frame
// and did something like sph for derivatives?

// todo: edge handling is not terribly robust.
// todo: smarter convolutions
// todo: where can we use duals?
// todo: can we do something better with product calculus?
// future: we should not alloc on operations. use a temporary operation
//         object which evals on assignment (or after a certain stack depth)

// todo: make this a proper class
// todo: give this an iterator

namespace geom {
    
namespace detail {
    
    // C++ is stupid and will not let you access the native operators as
    // function pointers. so we have to do this bullshit.
    
    template <typename T>
    inline T basicop_mul(const T &a, const T &b) {
        return a * b;
    }
    
    template <typename T>
    inline T basicop_div(const T &a, const T &b) {
        return a / b;
    }
    
    template <typename T>
    inline T basicop_add(const T &a, const T &b) {
        return a + b;
    }
    
    template <typename T>
    inline T basicop_sub(const T &a, const T &b) {
        return a - b;
    }
}

template <typename T, index_t N>
struct Grid {
    
    static const index_t DIM = N;
    typedef T elem_t;
    
    boost::shared_array<T> data;
    index_t n;
    index_t shape[N];  // this is stupid, make it a Vec<index_t,N>
    
protected:
    
    void init(index_t *extents) {
        n = 1;
        std::copy(extents, extents + N, shape);
        for (index_t i = 0; i < N; i++) {
            n *= shape[i];
        }
        data = boost::shared_array<T>(new T[n]);
    }
    
public:
    
    Grid(Vec<index_t,N> extents) { init(extents.begin()); }
    Grid(index_t *extents)       { init(extents); }
    
    inline T& operator[](const Vec<index_t,N> &p) {
        return data[index(p)];
    }
    
    inline T operator[](const Vec<index_t,N> &p) const {
        return data[index(p)];
    }
    
    index_t index(const Vec<index_t,N> &p) const {
        index_t i = 0;
        index_t blksz = 1;
        for (index_t j = 0; j < N; j++) {
            index_t c = std::max(std::min(p[j], shape[j] - 1), (index_t)0);
            i += blksz * c;
            blksz *= shape[j];
        }
        return i;
    }
    
    Grid<T,N> binop(const Grid<T,N> &g, T (*fn)(const T&, const T&)) const {
        // xxx todo: handle dimension mismatches
        Grid<T,N> out(shape);
        T *out_data = out.data.get();
        T *g0_data =   data.get();
        T *g1_data = g.data.get();
        for (index_t i = 0; i < n; i++) {
            out_data[i] = fn(g0_data[i], g1_data[i]);
        }
        return out;
    }
    
    inline Grid<T,N> operator+(const Grid<T,N> &g) const {
        return binop(g, detail::basicop_add<T>);
    }
    
    inline Grid<T,N> operator-(const Grid<T,N> &g) const {
        return binop(g, detail::basicop_sub<T>);
    }
    
    inline Grid<T,N> operator*(const Grid<T,N> &g) const {
        return binop(g, detail::basicop_mul<T>);
    }
    
    inline Grid<T,N> operator/(const Grid<T,N> &g) const {
        return binop(g, &detail::basicop_div<T>);
    }
    
    Vec<T,N> gradient(const Vec<index_t,N> &p) const {
        Vec<T,N> g;
        for (index_t i = 0; i < N; i++) {
            Vec<index_t,N> dx;
            dx[i] = 1;
            g[i] = (this->operator[](p + dx) - this->operator[](p - dx)) / (T)2;
        }
        return g;
    }
    
    inline Grid<T,N> dupe() const {
        Grid<T,N> out(shape);
        std::copy(data.get(), data.get() + n, out.data.get());
        return out;
    }
    
    Grid< Vec<T,N>, N > gradient() const {
        Grid< Vec<T,N>, N > out(shape);
        Vec<index_t,N> hi(shape);
        GridIterator<index_t,N> i(Vec<index_t,N>(), hi);
        for (; i != i.end(); i++) {
            Vec<index_t,N> p = *i;
            out[p] = gradient(p);
        }
        return out;
    }
    
    Grid<T,N> laplace() const {
        // less direct method, but more likely to be right. I think.
        Grid<T,N> out(shape);
        Grid<Vec<T,N>,N> grad = gradient();
        GridIterator<index_t,N> i = GridIterator<index_t,N>(Vec<index_t,N>(), Vec<index_t,N>(shape));
        GridIterator<index_t,N> end = i.end();
        for (; i != end; i++) {
            Vec<index_t,N> p = *i;
            T lapl = 0;
            for (index_t d = 0; d < N; d++) {
                Vec<index_t,N> dx;
                dx[d] = 1;
                lapl += (grad[p + dx][d] - grad[p - dx][d]) / (T)2;
            }
            out[p] = lapl;
        }
        return out;
    }
    
    /*
    Grid<T,N> laplace() const {
        // xxx: I guess this is wrong :(
        Grid<T,N> p0 = dupe();
        Grid<T,N> p1(shape);
        Vec<index_t,N> hi(shape);
        for (index_t axis = 0; axis < N; axis++) {
            GridIterator<index_t,N> i(Vec<index_t,N>(), hi);
            Vec<index_t,N> dx; 
            dx[axis] = 1;
            for (; i != i.end(); i++) {
                Vec<index_t,N> p = *i;
                p1[p] = p0[p + dx] - ((T)2) * p0[p] + p0[p - dx];
            }
            std::swap(p0, p1);
        }
        return p0;
    }*/
};

template <typename T, index_t N, typename U> 
Grid<T,N> grid_op(const Grid<T,N> &g, T v, T (*fn)(const T&, const U&)) {
    Grid<T,N> out(g.shape);
    T *out_data = out.data.get();
    T *g0_data  = g.data.get();
    for (index_t i = 0; i < g.n; i++) {
        out_data[i] = fn(g0_data[i], v);
    }
    return out;
}

template <typename T, index_t N, typename U> 
Grid<T,N> grid_op(U v, const Grid<T,N> &g, T (*fn)(const U&, const T&)) {
    Grid<T,N> out(g.shape);
    T *out_data = out.data.get();
    T *g0_data  =   g.data.get();
    for (index_t i = 0; i < g.n; i++) {
        out_data[i] = fn(v, g0_data[i]);
    }
    return out;
}

template <typename U, typename V>
struct NonGridType {
    typedef V type;
};

template <typename T, index_t N, typename V>
struct NonGridType<Grid<T,N>,V> {};


// we must specify a non-grid type so the free operators are not ambiguous with 
// the member operators on Grid<>. 

template <typename T, index_t N>
inline typename NonGridType< T,Grid<T,N> >::type operator*(const Grid<T,N> &g, T v) {
    return grid_op(g,v, detail::basicop_mul<T>);
}

template <typename T, index_t N>
inline typename NonGridType< T,Grid<T,N> >::type operator*(T v, const Grid<T,N> &g) {
    return grid_op(v,g, detail::basicop_mul<T>);
}

template <typename T, index_t N>
inline typename NonGridType< T,Grid<T,N> >::type operator/(const Grid<T,N> &g, T v) {
    return grid_op(g,v, detail::basicop_div<T>);
}

template <typename T, index_t N>
inline typename NonGridType< T,Grid<T,N> >::type operator/(T v, const Grid<T,N> &g) {
    return grid_op(v,g, detail::basicop_div<T>);
}

template <typename T, index_t N>
inline typename NonGridType< T,Grid<T,N> >::type operator+(const Grid<T,N> &g, T v) {
    return grid_op(g,v, detail::basicop_add<T>);
}

template <typename T, index_t N>
inline typename NonGridType< T,Grid<T,N> >::type operator+(T v, const Grid<T,N> &g) {
    return grid_op(v,g, detail::basicop_add<T>);
}

template <typename T, index_t N>
inline typename NonGridType< T,Grid<T,N> >::type operator-(const Grid<T,N> &g, T v) {
    return grid_op(g,v, detail::basicop_sub<T>);
}

template <typename T, index_t N>
inline typename NonGridType< T,Grid<T,N> >::type operator-(T v, const Grid<T,N> &g) {
    return grid_op(v,g, detail::basicop_sub<T>);
}

} //namespace geom

#endif	/* GRID_H */

