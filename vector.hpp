#pragma once

#include <stddef.h>                         // defines: size_t
#include <cmath>
#include <exception>

// valid arithmetic
// (x) point - point = vector
// (x) point + vector = point
// (x) vector + point = point
// (x) point - vector = point
// (x) vector - point = point
// (x) vector + vector = vector
// (x) vector - vector = vector
// (x) scalar * vector = vector
// (x) vector * scalar = vector
// (TODO: not possible with type system?) -vector = vector



typedef float scalar_t;
 
class OutOfBoundsException : public std::exception {};


//
template<typename T, size_t N, typename E>
struct HasElementBase {
        T operator[](size_t i) const { return static_cast<E const&>(*this)[i]; }
        size_t size() { return N; }
};

// Leave empty for specialization
template<typename T, size_t N, typename E>
struct PointExpression : public HasElementBase<T, N, E> {};

//
template<typename T, size_t N, typename E>
class VectorExpressionBase : public HasElementBase<T, N, E> {
        // Dot product expands at compile time to : u[N]*v[N] + u[N-1]*v[N-1] + ... + u[0]*v[0]
        template <typename E2, size_t W> 
        struct DotProduct  { 
                static inline T eval(const VectorExpressionBase<T, N, E>& u, const VectorExpressionBase<T, N, E2>& v) {
                        return  u[W]*v[W] + DotProduct<E2, W-1>::eval(u, v);  
                } 
        }; 

        template<typename E2> 
        struct DotProduct<E2, 0> { 
                static inline T eval(const VectorExpressionBase<T, N, E>& u, const VectorExpressionBase<T, N, E2>& v) { 
                        return u[0]*v[0];
                } 
        };
public:
        template<typename E2>
        inline T dot(VectorExpressionBase<T, N, E2>& b) const { 
                return DotProduct<E2, N-1>::eval(*this, b);
        }
        
    T length() {
                return sqrt( dot(*this) );
        }
};

// Leave empty for specialization
template<typename T, size_t N, typename E>
struct VectorExpression : public VectorExpressionBase<T, N, E> {};


template<typename T, size_t N1, size_t N2, typename E>
struct MatrixExpressionBase {
        VectorExpression<T, N2, E> operator[](size_t i) const { return static_cast<E const&>(*this)[i]; }
        size_t nRows() { return N1; }
        size_t nCols() { return N2; }
};

// Leave empty for specialization
template<typename T, size_t N1, size_t N2, typename E>
struct MatrixExpression : public MatrixExpressionBase<T, N1, N2, E> {};

//



// Define the possible vector operations
template <typename T, size_t N, typename E>
struct VectorScaled : public VectorExpression<T, N, VectorScaled<T, N, E> > {
        VectorScaled(T alpha, VectorExpression<T, N, E> const& v) : _alpha(alpha), _v(v) {}
        T operator[](size_t i) const { 
                if(i >= N || i < 0) throw OutOfBoundsException();
                return _alpha * _v[i]; 
        }

private:
        T _alpha; 
        VectorExpression<T, N, E> const& _v; 
};

template <typename T, size_t N, typename E1, typename E2>
struct VectorSum : public VectorExpression<T, N, VectorSum<T, N, E1, E2> > {
        VectorSum(VectorExpression<T, N, E1> const& u, VectorExpression<T, N, E2> const& v) : _u(u), _v(v) {};
        T operator[](size_t i) const { 
                if(i >= N || i < 0) throw OutOfBoundsException();
                return _u[i] + _v[i]; 
        }

private:
        VectorExpression<T, N, E1>  const& _u; 
        VectorExpression<T, N, E2>  const& _v; 
};

template <typename T, size_t N, typename E1, typename E2>
struct PointDiff : public VectorExpression<T, N, PointDiff<T, N, E1, E2> > {
        PointDiff(PointExpression<T, N, E1> const& u, PointExpression<T, N, E2> const& v) : _u(u), _v(v) {};
        T operator[](size_t i) const { 
                if(i >= N || i < 0) throw OutOfBoundsException();
                return _u[i] - _v[i]; 
        }

private:
        PointExpression<T, N, E1>  const& _u; 
        PointExpression<T, N, E2>  const& _v; 
};

template <typename T, size_t N, typename E1, typename E2>
struct PointVectorSum : public PointExpression<T, N, PointVectorSum<T, N, E1, E2> > {
        PointVectorSum(PointExpression<T, N, E1> const& u, VectorExpression<T, N, E2> const& v) : _u(u), _v(v) {};
        T operator[](size_t i) const { 
                if(i >= N || i < 0) throw OutOfBoundsException();
                return _u[i] + _v[i]; 
        }

private:
        PointExpression<T, N, E1>  const& _u; 
        VectorExpression<T, N, E2>  const& _v; 
};

template <typename T, size_t N, typename E1, typename E2>
struct PointVectorDiff : public PointExpression<T, N, PointVectorDiff<T, N, E1, E2> > {
        PointVectorDiff(PointExpression<T, N, E1> const& u, VectorExpression<T, N, E2> const& v) : _u(u), _v(v) {};
        T operator[](size_t i) const { 
                if(i >= N || i < 0) throw OutOfBoundsException();
                return _u[i] - _v[i]; 
        }

private:
        PointExpression<T, N, E1>  const& _u; 
        VectorExpression<T, N, E2>  const& _v; 
};


// NOTE: There are limitations to expression templates, notice how I'm only
//       using scalar T and not a more general typename T
//       if you wish to use something else you can simply redefine scalar_t
//       or copy all expressions with a new type
// TODO: There probably is a way to embed these definitions inside VectorBase
//       so that they are generic against all types T
template <typename E, size_t N>
VectorScaled<scalar_t, N, E> const
operator*(scalar_t alpha, VectorExpression<scalar_t, N, E> const& v) {
        return VectorScaled<scalar_t, N, E>(alpha,v);
}

template <typename E, size_t N>
VectorScaled<scalar_t, N, E> const
operator*(VectorExpression<scalar_t, N, E> const& v, scalar_t alpha) {
        return VectorScaled<scalar_t, N, E>(alpha,v);
}

template <typename E1, typename E2, size_t N>
VectorSum<scalar_t, N, E1, E2> const
operator+(VectorExpression<scalar_t, N, E1> const& u, VectorExpression<scalar_t, N, E2> const& v) {
        return VectorSum<scalar_t, N, E1, E2>(u,v);
}

template <typename E1, typename E2, size_t N>
PointDiff<scalar_t, N, E1, E2> const
operator-(PointExpression<scalar_t, N, E1> const& u, PointExpression<scalar_t, N, E2> const& v) {
        return PointDiff<scalar_t, N, E1, E2>(u,v);
}

template <typename E1, typename E2, size_t N>
PointVectorSum<scalar_t, N, E1, E2>
operator+(PointExpression<scalar_t, N, E1> const& u, VectorExpression<scalar_t, N, E2> const& v) {
        return PointVectorSum<scalar_t, N, E1, E2>(u,v);
}

template <typename E1, typename E2, size_t N>
PointVectorSum<scalar_t, N, E1, E2>
operator+(VectorExpression<scalar_t, N, E2> const& v, PointExpression<scalar_t, N, E1> const& u) {
        return PointVectorSum<scalar_t, N, E1, E2>(u,v);
}

template <typename E1, typename E2, size_t N>
PointVectorDiff<scalar_t, N, E1, E2>
operator-(PointExpression<scalar_t, N, E1> const& u, VectorExpression<scalar_t, N, E2> const& v) {
        return PointVectorDiff<scalar_t, N, E1, E2>(u,v);
}

template <typename E1, typename E2, size_t N>
PointVectorDiff<scalar_t, N, E1, E2>
operator-(VectorExpression<scalar_t, N, E2> const& v, PointExpression<scalar_t, N, E1> const& u) {
        return PointVectorDiff<scalar_t, N, E1, E2>(u,v);
}

// template specialization for VectorExpression with dim = 3
//       so that this: (e1 + e2).cross() expands properly
template<typename T,  typename E1>
class VectorExpression<T, 3, E1> : public VectorExpressionBase<T, 3, E1> {

        // cross product is only valid for 3rd dimension
        template <typename E2>
        struct VectorCross : public VectorExpression<T, 3, VectorCross<E2> > {
                VectorCross(E1 const& u, E2 const& v) : _u(u), _v(v) {};
                T operator[](size_t i) const { 
                        if(i == 0) return _u[1]*_v[2] - _u[2]*_v[1];
                        if(i == 1) return _u[2]*_v[0] - _u[0]*_v[2];
                        if(i == 2) return _u[0]*_v[1] - _u[1]*_v[0];
                        throw OutOfBoundsException();
                }

        private:
                E1 const& _u; 
                E2 const& _v; 
        };

public:         
        template<typename E2>
        inline VectorCross<VectorExpression<T, 3, E2> > const
        cross(VectorExpression<T, 3, E2> const& b) const {
                return VectorCross<VectorExpression<T, 3, E2> > ( static_cast<E1 const&>(*this), b );
        }
};

// Finally we can define the Vector class
// Define Vectors (they have storage)
template<typename T, size_t N>
struct VectorBase : public VectorExpression< T, N, VectorBase<T, N> > {
        VectorBase() {} 
        T operator[](size_t i) const {
                if(i >= N || i < 0) throw OutOfBoundsException();
                return _data[i];
        }
        
        T& operator[](size_t i) {
                if(i >= N || i < 0) throw OutOfBoundsException();
                return _data[i];
        }
        
        template<typename E>
        VectorBase& operator=(VectorExpressionBase< T, N, E> const& q) {            
            for(size_t i = 0; i < N; i++) {
                    _data[i] = q[i];
            }
            return *this;
        }
        
protected:
        T _data[N];
};

// template specialization to add nice constructors for useful classes
template<typename T, size_t N>
struct Vector : public VectorBase<T, N> {
        Vector() {}
        
        // Construct from any VectorExpression:
        template<typename E>
        Vector(const VectorExpression<T, N, E>& v) {
                for (size_t i = 0; i < N; i++)
                        this->_data[i] = v[i];
        }       
};

template<typename T>
struct Vector<T, 3> : public VectorBase<T, 3> {
        Vector() {}
        
        Vector(T x, T y, T z) {
                this->_data[0] = x;
                this->_data[1] = y;
                this->_data[2] = z;
        }
        // Construct from any VectorExpression:
        template<typename E>
        Vector(const VectorExpression<T, 3, E>& v) {
                for (size_t i = 0; i < 3; i++)
                        this->_data[i] = v[i];
        }       
};

template<typename T>
struct Vector<T, 4> : public VectorBase<T, 4> {
        Vector() {}
        Vector(T x, T y, T z, T w) {
                this->_data[0] = x;
                this->_data[1] = y;
                this->_data[2] = z;
                this->_data[3] = w;
        }
        // Construct from any VectorExpression:
        template<typename E>
        Vector(const VectorExpression<T, 4, E>& v) {
                for (size_t i = 0; i < 4; i++)
                        this->_data[i] = v[i];
        }       
};



// Start defining the point operations
template<typename T, size_t N>
struct PointBase : public  PointExpression<T, N, PointBase<T, N> >{
        PointBase() {}
        
        T operator[](size_t i) const {
                if(i >= N || i < 0) throw OutOfBoundsException();
                return _data[i];
        }
        
        T& operator[](size_t i) {
                if(i >= N || i < 0) throw OutOfBoundsException();
                return _data[i];
        }
        
        template<typename E>
        PointBase& operator=(PointExpression< T, N, E> const& q) {
                for(size_t i = 0; i < N; i++) {
                        _data[i] = q[i];
                }
                return *this;
        }
        

        
protected:
        T _data[N];
};

// template specialization to add nice constructors for useful classes
template<typename T, size_t N>
struct Point : public PointBase<T, N> {
        Point() {}
        // Construct from any PointExpression:
        template<typename E>
        Point(const PointExpression<T, N, E>& v) {
                for (size_t i = 0; i < N; i++)
                        this->_data[i] = v[i];
        }               
};

template<typename T>
struct Point<T, 3> : public PointBase<T, 3> {
        Point() {}
        Point(T x, T y, T z) {
                this->_data[0] = x;
                this->_data[1] = y;
                this->_data[2] = z;
        }
        // Construct from any PointExpression:
        template<typename E>
        Point(const PointExpression<T, 3, E>& v) {
                for (size_t i = 0; i < 3; i++)
                        this->_data[i] = v[i];
        }               
};

template<typename T>
struct Point<T, 4> : public PointBase<T, 4> {
        Point() {}
        Point(T x, T y, T z, T w) {
                this->_data[0] = x;
                this->_data[1] = y;
                this->_data[2] = z;
                this->_data[3] = w;
        }
        // Construct from any PointExpression:
        template<typename E>
        Point(const PointExpression<T, 4, E>& v) {
                for (size_t i = 0; i < 4; i++)
                        this->_data[i] = v[i];
        }               
};

template<typename T, size_t N1, size_t N2>
struct MatrixBase : public MatrixExpression<T, N1, N2, MatrixBase<T, N1, N2> > {
        MatrixBase() {}
        
        Vector<T, N2>& operator[](size_t i) const {
                if(i >= N1 || i < 0) throw OutOfBoundsException();
                return _data[i];
        }

protected:
        Vector<T, N2> _data[N1];
};

template<typename T, size_t N1, size_t N2>
struct Matrix : public MatrixBase<T, N1, N2> {
        Matrix() {}
        // Construct from any PointExpression:
        template<typename E>
        Matrix(const MatrixExpression<T, N1, N2, E>& v) {
                for (size_t i = 0; i < N1; i++)
                        this->_data[i] = v[i];
        }               
};

template<typename T>
struct Matrix<T, 4, 4> : public MatrixBase<T, 4, 4> {
        Matrix() {}
        // Construct from any PointExpression:
        template<typename E>
        Matrix(const MatrixExpression<T, 4, 4, E>& v) {
                for (size_t i = 0; i < 4; i++)
                        this->_data[i] = v[i];
        }               
        
        template<typename E>
        Matrix(VectorExpression<T, 4, E> &v1, VectorExpression<T, 4, E> &v2, VectorExpression<T, 4, E> &v3, VectorExpression<T, 4, E> &v4)
        {
                this->_data[0] = v1;
                this->_data[1] = v2;
                this->_data[2] = v3;
                this->_data[3] = v4;
        }
};


// Helper subclasses
typedef Vector<scalar_t, 3> Vector3d;
typedef Point<scalar_t, 3> Point3d;
typedef Vector<scalar_t, 4> Vector4d;
typedef Point<scalar_t, 4> Point4d;