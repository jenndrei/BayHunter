/******************************************************************
 *                                                                *
 *   Copyright by Joachim Saul <saul@gfz-potsdam.de>              *
 *                                                                *
 ******************************************************************/

# ifndef CMAT2_H
# define CMAT2_H
# include "Complex.h"

class Cmat2
{
    // defines complex 2 x 2 matrices and operations
  public:
    // constructors
    Cmat2 () {}
    Cmat2 (Complex const &f);
    Cmat2 (Complex const &f11, Complex const &f12,
	   Complex const &f21, Complex const &f22);
    // destructor
//    ~Cmat2 () {}

    void get (Complex *f11, Complex *f12,
	      Complex *f21, Complex *f22);

    // member operators
    void operator += (Cmat2   const &);
    void operator -= (Cmat2   const &);
    void operator *= (Cmat2   const &);
    void operator *= (Complex const &);
    Cmat2 &operator =  (Complex const &);

    // friend functions 
    friend Complex comp11      (Cmat2   const &m);

    friend Complex det         (Cmat2   const &m);
    friend Cmat2   inv         (Cmat2   const &m);

    // unary operators
    friend Cmat2   operator +  (Cmat2   const &m);
    friend Cmat2   operator -  (Cmat2   const &m);

    // binary operators
    friend Cmat2   operator +  (Cmat2   const &x, Cmat2   const &y);
    friend Cmat2   operator -  (Cmat2   const &x, Cmat2   const &y);
    friend Cmat2   operator *  (Cmat2   const &x, Cmat2   const &y);
    friend Cmat2   operator *  (Complex const &r, Cmat2   const &x);
    friend Cmat2   operator *  (Cmat2   const &x, Complex const &r);
    friend Cmat2   operator *  (double  const &r, Cmat2   const &x);
    friend Cmat2   operator *  (Cmat2   const &x, double  const &r);

//  private:
    Complex c11, c12, c21, c22;
};

inline
Cmat2::Cmat2 (Complex const &f)
{
    c11 = c22 = f;
    c12 = c21 = 0.;
}

inline 
Cmat2::Cmat2 (Complex const &f11, Complex const &f12,
	      Complex const &f21, Complex const &f22)
{
    c11=f11; c12=f12;
    c21=f21; c22=f22;
}

inline void 
Cmat2::get (Complex *f11, Complex *f12,
	    Complex *f21, Complex *f22)
{
    *f11=c11; *f12=c12;
    *f21=c21; *f22=c22;
}

inline void 
Cmat2::operator+= (Cmat2 const &other)
{
    c11 += other.c11; c12 += other.c12;
    c21 += other.c21; c22 += other.c22;
}


inline void 
Cmat2::operator-= (Cmat2 const &other)
{
    c11 -= other.c11;
    c12 -= other.c12;
    c21 -= other.c21;
    c22 -= other.c22;
}

inline Cmat2& 
Cmat2::operator= (Complex const &c)
{
    c11 = c22 = c;
    c12 = c21 = 0.;

    return *this;
}

inline Complex
comp11 (Cmat2 const &x)
{
   return x.c11;
}

inline Complex
det (Cmat2 const &x)
{
    return x.c11*x.c22 - x.c12*x.c21;
}

inline Cmat2
operator+ (Cmat2 const &x, Cmat2 const &y)
{
    return Cmat2 (x.c11 + y.c11, x.c12 + y.c12, 
		  x.c21 + y.c21, x.c22 + y.c22);
}

inline Cmat2
operator+ (Cmat2 const &x)
{
    return x;
}

inline Cmat2
operator- (Cmat2 const &x, Cmat2 const &y)
{
    return Cmat2 (x.c11 - y.c11, x.c12 - y.c12, 
		  x.c21 - y.c21, x.c22 - y.c22);
}

inline Cmat2
operator- (Cmat2 const &x)
{
    return Cmat2 (-x.c11, -x.c12, -x.c21, -x.c22);
}

inline Cmat2
inv (Cmat2 const &x)
{
    Complex q;

    q = 1. / (x.c11*x.c22 - x.c12*x.c21);

    return Cmat2 ( q*x.c22, -q*x.c12,
		  -q*x.c21,  q*x.c11);
}

inline void 
Cmat2::operator*= (Cmat2 const &other)
{
    // temporary matrix
    Complex x11, x12, x21, x22;  
  
    x11 = c11 * other.c11 + c12 * other.c21;
    x12 = c11 * other.c12 + c12 * other.c22;
    x21 = c21 * other.c11 + c22 * other.c21;
    x22 = c21 * other.c12 + c22 * other.c22;
 
    c11 = x11;
    c12 = x12;
    c21 = x21;
    c22 = x22;
}

inline Cmat2
operator* (Cmat2 const &x, Cmat2 const &y)
{
    return Cmat2 (x.c11 * y.c11 + x.c12 * y.c21,
		  x.c11 * y.c12 + x.c12 * y.c22,
		  x.c21 * y.c11 + x.c22 * y.c21,
		  x.c21 * y.c12 + x.c22 * y.c22);
}

inline Cmat2
operator* (Complex const &r, Cmat2 const &x)
{
    return Cmat2 (r*x.c11, r*x.c12,
		  r*x.c21, r*x.c22);
}

inline Cmat2
operator* (Cmat2 const &x, Complex const &r)
{
    return Cmat2 (r*x.c11, r*x.c12,
		  r*x.c21, r*x.c22);
}

inline Cmat2
operator* (double const &r, Cmat2 const &x)
{
    return Cmat2 (r*x.c11, r*x.c12,
		  r*x.c21, r*x.c22);
}

inline Cmat2
operator* (Cmat2 const &x, double const &r)
{
    return Cmat2 (r*x.c11, r*x.c12,
		  r*x.c21, r*x.c22);
}

inline void 
Cmat2::operator*= (Complex const &r)
{
    c11 *= r;
    c12 *= r;
    c21 *= r;
    c22 *= r;
}

# endif
