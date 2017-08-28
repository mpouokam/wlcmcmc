#include <iostream>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

#include "vector3d.hpp"
#include "mtrand/mtrand.h"

using namespace std;

MTRand_int32 mtrandint;
MTRand mtrand;


int rand_integer(int a, int b)
{
    return (mtrandint() % (b - a)) + a;
}

double rand_real(double a, double b)
{
    return (mtrand() * (b - a)) + a;
}


inline void _do_crank_shaft(vector3d* vertices, int numberofvertices, int m, int n, double theta)
{
	// if n (the number of segments to rotate) is congruent to 0, +1, -1 then do nothing
	switch ((n+1)%numberofvertices) { case 0: case 1: case 2: return; }
	
	vector3d axisofrotation = normalize(vertices[m] - vertices[(m+n)%numberofvertices]);

	for (unsigned i=1; i<n; i++) // we going to do this for n different vertices
	{
		vector3d tmp;
		int thisone = (i+m)%numberofvertices; // the index of the vertex we will rotate
		tmp = (vertices[thisone] - vertices[m]); // shift so that a point on the axis is now on the origin

		// rodrigues's formula
		vertices[thisone] =	tmp*cos(theta) + crossprod(axisofrotation,tmp) * sin(theta)
							+ axisofrotation * dotprod(axisofrotation,tmp) * (1-cos(theta))
							+ vertices[m];	//shift back away from origin
	}
}

void do_crank_shaft(double* coords, int numberofvertices, int ibe3, int m, int n, double theta)
{
	vector3d *vertices = reinterpret_cast<vector3d*>(coords);
	_do_crank_shaft(vertices, numberofvertices, m, n, theta);
}
inline void _do_random_crank_shaft(vector3d* vertices, int numberofvertices, double thetamax)
{
	// get random values
	int m = rand_integer(0,numberofvertices-1);
	int n = rand_integer(2,numberofvertices-2);
	double theta = rand_real(-thetamax,thetamax);

	_do_crank_shaft(vertices,numberofvertices,m,n,theta);
}

void do_random_crank_shaft(double* coords, int numberofvertices, int ibe3, double thetamax)
{
	vector3d *vertices = reinterpret_cast<vector3d*>(coords);
	_do_random_crank_shaft(vertices, numberofvertices, thetamax);
}


inline double _get_wlc_energy(vector3d* vertices, int n, double bendingrigidityconstant) // the unit being kBT
{
	if (bendingrigidityconstant==0.0)
		return 0;

	vector3d edges[n];


	edges[n-1]=normalize(vertices[0]-vertices[n-1]);
	for (unsigned i=0 ; i < n-1; i++)
		edges[i]=normalize(vertices[i+1]-vertices[i]);

	
	double tmp = acos(dotprod(edges[n-1],edges[0]));
	double output = tmp*tmp;
	for (unsigned i=0 ; i < n-1; i++)
	{
		tmp =  acos(dotprod(edges[i],edges[i+1]));
		output += tmp*tmp;
	}

	output *= 0.5*bendingrigidityconstant;
	return output;

}

double get_wlc_energy(double* coords, int numberofvertices, int ibe3, double bendingrigidityconstant)
{
	vector3d *vertices = reinterpret_cast<vector3d*>(coords);
	return _get_wlc_energy(vertices, numberofvertices, bendingrigidityconstant);
}

inline bool is_monotonic(const double& a,const double& b,const double& c)
{
	return ((a<=b && b<=c) || (c<=b && b<=a));
}

inline int sgn(double val) {
   return (double(0) < val) - (val < double(0));
}


inline double _distance_between_edges(vector3d * vertices, int numberofvertices, int _a1, int _a2, int _b1, int _b2)
{
	// maybe i could get these on the stack. this was faster when I had it in c++
	vector3d a1 = vertices[_a1];
	vector3d a2 = vertices[_a2];
	vector3d b1 = vertices[_b1];
	vector3d b2 = vertices[_b2];

	vector3d a = (a1+a2)/2; // translate a1 and a2 to midpoint a	
	if (a1!=a2)
	{
		vector3d m = a2-a1;

		//translate b1
		if (is_monotonic(dotprod(a1,m),dotprod(b1,m),dotprod(a2,m))) // ie, is in orthogonalspace in between
			b1 += m*dotprod(a-b1,m)/dotprod(m,m); // project onto orthogonal bisector
		else
			b1 += m/2*sgn(dotprod(a-b1,m)); // translate toward the projection

		// translate b2
		if (is_monotonic(dotprod(a1,m),dotprod(b2,m),dotprod(a2,m))) // ie, is in orthogonalspace in between
			b2 += m*dotprod(a-b2,m)/dotprod(m,m);
		else
			b2 += m/2*sgn(dotprod(a-b2,m));
	}

	vector3d b=(b1+b2)/2;
	if (b1!=b2)
	{
		vector3d n = b2-b1;

		// translate a
		if (is_monotonic(dotprod(b1,n),dotprod(a,n),dotprod(b2,n))) // ie, is in orthogonalspace in between
			a += n*dotprod(b-a,n)/dotprod(n,n);
		else
			a += n/2*sgn(dotprod(b-a,n));;
	}
	return norm(a-b);
}
double distance_between_edges(double* coords, int numberofvertices, int ibe3, int _a1, int _a2, int _b1, int _b2)
{
	return _distance_between_edges(reinterpret_cast<vector3d*>(coords),numberofvertices, _a1, _a2, _b1, _b2);
}

bool _is_collision(vector3d *vertices, int N, double diameter, double edgelength)
{
	if (diameter==0) return false;
	for (unsigned i=0; i<N-2; i++)
	{
		for (unsigned j=i+2; j<N-1; j++)
		{
			double gap = _distance_between_edges(vertices,N,i,i+1,j,j+1)-diameter;
			if (gap < 0)
				return true;
			else
				// j+=1;
				j+=int((gap-diameter)/(edgelength));

		}
		if (i>0 && ( 0 > (_distance_between_edges(vertices,N,i,i+1,N-1,0)-diameter)))
			return true;
	}

	return false;
}
bool is_collision(double* coords, int numberofvertices, int ibe3, double diameter, double edgelength)
{
	return _is_collision(reinterpret_cast<vector3d*>(coords),numberofvertices, diameter, edgelength);
}

inline int _do_monte_carlo_steps(vector3d* vertices, int n, int numberofsteps, double bendingrigidityconstant, double thetamax, double diameter, double edgelength)
{

	int acceptances = 0;
	vector3d candidatevertices[n];

	for (unsigned i=0; i<numberofsteps; i++)
	{
		//This is the code for a single MC step
		memcpy(candidatevertices,vertices,n*sizeof(vector3d));

		_do_random_crank_shaft(candidatevertices,n,thetamax);
		double deltaenergy = _get_wlc_energy(candidatevertices,n,bendingrigidityconstant) -
		                     _get_wlc_energy(vertices,n,bendingrigidityconstant);

		if ( (deltaenergy<0 || exp(-deltaenergy)>rand_real(0, 1)) && !_is_collision(candidatevertices,n,diameter,edgelength) )
		{
			// ACCEPT
			memcpy(vertices,candidatevertices,n*sizeof(vector3d));
			acceptances+=1;
		}
		else // REJECT
		{
		}

	}

	return acceptances;

}

int do_monte_carlo_steps(double* coords, int numberofvertices, int ibe3, int numberofsteps, double bendingrigidityconstant, double thetamax, double diameter, double edgelength)
{
	vector3d *vertices = reinterpret_cast<vector3d*>(coords);
	return _do_monte_carlo_steps(vertices, numberofvertices, numberofsteps, bendingrigidityconstant,thetamax,diameter,edgelength);
}

void set_random_seed()
{
    int seed = ::time(NULL) + ::clock() + ::getpid();
	cout << "Setting random seed: " << seed << endl;

	MTRand mt(seed);
}

int main(int argc, char const *argv[])
{
    MTRand(38);
}
