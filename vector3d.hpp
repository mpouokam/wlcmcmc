

#ifndef __vector3d_hpp__
#define __vector3d_hpp__

#include <iostream>
#include <sstream>

using std::ostream;
using std::cout;

#include <math.h>


class vector3d {

	public:
	// We need to make everything public and then we don't have to go crazy with friend declarations later.
	
	double x;
	double y;
	double z;


	vector3d()										{ x=0; y=0; z=0; }
	vector3d(double xin, double yin, double zin)	{ x=xin; y=yin; z=zin; }

	// vevtor addition
	inline vector3d & operator+=		(const vector3d & other)				{ x+=other.x; y+=other.y; z+=other.z; return *this; }
	inline vector3d operator+			(const vector3d & other)	const		{ return vector3d(*this)+=other; }
	inline vector3d & operator-=		(const vector3d & other)				{ x-=other.x; y-=other.y; z-=other.z; return *this; }
	inline vector3d operator-			(const vector3d & other)	const		{ return vector3d(*this)-=other; }
	//																	^^^^^ const because you can't for instance say (a-b)=3;

	inline vector3d operator-			()	const								{ return vector3d(*this)*=-1; }

	// scalar multiplication / division
	inline vector3d & operator*=		(const double k)						{ x*=k; y*=k; z*=k; return *this; }
	inline vector3d & operator/=		(const double k)						{ x/=k; y/=k; z/=k; return *this; }
	inline friend vector3d operator*	(const double k, const vector3d & vec)	{ return vector3d(vec)*=k; }
	inline friend vector3d operator*	(const vector3d & vec, const double k)	{ return vector3d(vec)*=k; }
	inline friend vector3d operator/	(const vector3d & vec, const double k)	{ return vector3d(vec)/=k; }

	inline friend bool operator ==		(const vector3d &a, const vector3d &b)	{ return (a.x==b.x && a.y==b.y && a.z==b.z);}
	inline friend bool operator !=		(const vector3d &a, const vector3d &b)	{ return !(a==b);}

	friend ostream & operator<<			(ostream & os, const vector3d & vec)	{ return os << vec.x << " " << vec.y << " " << vec.z; }
	// The following function is for python to use
	const char * __str__()				{std::ostringstream oos; oos << (*this); return oos.str().c_str();}

	inline double _dotprod				(const vector3d vec)			const	{ return x*vec.x+y*vec.y+z*vec.z; }
	inline vector3d _crossprod			(const vector3d v)				const	{ return vector3d(y*v.z-z*v.y,z*v.x-x*v.z,x*v.y-y*v.x); }


	inline double _norm					()								const	{ return sqrt(this->_dotprod(*this));}
	inline double _norm1					()							const	{ return this->x+this->y+this->z;}
	//																	^^^^^ needed if we are access a member function of a const reference
	// See this: http://stackoverflow.com/questions/98705/what-are-the-semantics-of-a-const-member-function

	inline vector3d & _normalize			()										{ return *this /= this->_norm(); }
//	double & operator[]	(const int index)					{	switch (index%3) {case 0: return x; case 1: return y; case 2: return z; } return 0.0;}

};

// I'm not sure why you can't make v1 as const here
inline const double dotprod(const vector3d & v1, const vector3d & v2 )
{
	return v1._dotprod(v2);
}

inline const vector3d crossprod(const vector3d & v1, const vector3d & v2)
{
	return v1._crossprod(v2);
}

inline const double det3x3(const vector3d & r0, const vector3d & r1, const vector3d & r2)
{
	return r0.x*(r1.y*r2.z-r1.z*r2.y)-r0.y*(r1.x*r2.z-r1.z*r2.x)+r0.z*(r1.x*r2.y-r1.y*r2.x);
}

inline double norm(const vector3d & vector)
{
	return vector._norm();
}

inline double norm1(const vector3d & vector)
{
	return vector._norm1();
}

inline const vector3d normalize(const vector3d &v)
{
	return vector3d(v)._normalize();
}

#endif