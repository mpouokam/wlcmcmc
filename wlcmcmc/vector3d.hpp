#ifndef __vector3d_hpp__
#define __vector3d_hpp__

#include <iostream>
#include <sstream>

using std::ostream;
using std::cout;

#include <math.h>


class vector3d {
	public:

	double x;
	double y;
	double z;


	vector3d()										{ x=0; y=0; z=0; }
	vector3d(double xin, double yin, double zin)	{ x=xin; y=yin; z=zin; }

	// vector addition
	vector3d & operator+=		(const vector3d & other)				{ x+=other.x; y+=other.y; z+=other.z; return *this; }
	vector3d operator+			(const vector3d & other)	const		{ return vector3d(*this)+=other; }
	vector3d & operator-=		(const vector3d & other)				{ x-=other.x; y-=other.y; z-=other.z; return *this; }
	vector3d operator-			(const vector3d & other)	const		{ return vector3d(*this)-=other; }

	vector3d operator-			()	const								{ return vector3d(*this)*=-1; }

	// scalar multiplication / division
	vector3d & operator*=		(const double k)						{ x*=k; y*=k; z*=k; return *this; }
	vector3d & operator/=		(const double k)						{ x/=k; y/=k; z/=k; return *this; }
	friend vector3d operator*	(const double k, const vector3d & vec)	{ return vector3d(vec)*=k; }
	friend vector3d operator*	(const vector3d & vec, const double k)	{ return vector3d(vec)*=k; }
	friend vector3d operator/	(const vector3d & vec, const double k)	{ return vector3d(vec)/=k; }

	friend bool operator ==		(const vector3d &a, const vector3d &b)	{ return (a.x==b.x && a.y==b.y && a.z==b.z);}
	friend bool operator !=		(const vector3d &a, const vector3d &b)	{ return !(a==b);}

	friend ostream & operator<<	(ostream & os, const vector3d & vec)	{ return os << vec.x << " " << vec.y << " " << vec.z; }

	// The following function is for python to use
	const char * __str__()				{std::ostringstream oos; oos << (*this); return oos.str().c_str();}

	double _dotprod				(const vector3d vec)			const	{ return x*vec.x+y*vec.y+z*vec.z; }
	vector3d _crossprod			(const vector3d v)				const	{ return vector3d(y*v.z-z*v.y,z*v.x-x*v.z,x*v.y-y*v.x); }


	double _norm					()								const	{ return sqrt(this->_dotprod(*this));}
	double _norm1					()							const	{ return this->x+this->y+this->z;}

	vector3d & _normalize			()										{ return *this /= this->_norm(); }

};


const double dotprod(const vector3d & v1, const vector3d & v2 )
{
	return v1._dotprod(v2);
}

const vector3d crossprod(const vector3d & v1, const vector3d & v2)
{
	return v1._crossprod(v2);
}

const double det3x3(const vector3d & r0, const vector3d & r1, const vector3d & r2)
{
	return r0.x*(r1.y*r2.z-r1.z*r2.y)-r0.y*(r1.x*r2.z-r1.z*r2.x)+r0.z*(r1.x*r2.y-r1.y*r2.x);
}

double norm(const vector3d & vector)
{
	return vector._norm();
}

double norm1(const vector3d & vector)
{
	return vector._norm1();
}

 const vector3d normalize(const vector3d &v)
{
	return vector3d(v)._normalize();
}

#endif