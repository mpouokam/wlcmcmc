#include <iostream>
#include "vector3d.hpp"
#include "random.h"
#include <string.h>

// TODO: Replace random.h with Python's MT routine
//#include "mtrand.h"

using namespace std;

inline void _docrankshaft(vector3d* vertices, int numberofvertices, int m, int n, double theta)
{
	// if n (the number of segments to rotate) is congruent to 0, +1, -1 then do nothing
	switch ((n+1)%numberofvertices) { case 0: case 1: case 2: return; }
	
	vector3d axisofrotation = normalize(vertices[m] - vertices[(m+n)%numberofvertices]);

	for (unsigned i=1; i<n; i++) // we going to do this for n different vertices
	{
		vector3d tmp;
		int thisone = (i+m)%numberofvertices; // the index of the vertex we will rotate
		tmp = (vertices[thisone] - vertices[m]); 	// shift so that a point on the axis is now on the origin
		// rodrigues's formula
		vertices[thisone] =		tmp*cos(theta) + crossprod(axisofrotation,tmp) * sin(theta)
								+ axisofrotation * dotprod(axisofrotation,tmp) * (1-cos(theta))
								+ vertices[m];	//shift back away from origin
	}
}

// TODO: Put this wrapper function in the SWIG interface file
void docrankshaft(double* coords, int numberofvertices, int ibe3, int m, int n, double theta)
{
	vector3d *vertices = reinterpret_cast<vector3d*>(coords);
	_docrankshaft(vertices, numberofvertices, m, n, theta);
}
inline void _dorandomcrankshaft(vector3d* vertices, int numberofvertices, double thetamax)
{
	// Call to random.cpp
	// TODO: Move this out of the innerloop

	// get random values
	int m = rand_integer(0,numberofvertices-1);
	int n = rand_integer(2,numberofvertices-2);
	double theta = rand_real(-thetamax,thetamax);

	_docrankshaft(vertices,numberofvertices,m,n,theta);
}

// TODO: Put this wrapper function in the SWIG interface file
void dorandomcrankshaft(double* coords, int numberofvertices, int ibe3, double thetamax)
{
	vector3d *vertices = reinterpret_cast<vector3d*>(coords);
	_dorandomcrankshaft(vertices, numberofvertices, thetamax);
}


inline double _getwlcenergy(vector3d* vertices, int n, double bendingrigidityconstant)
{
	if (bendingrigidityconstant==0.0) // Freely Jointed Chain
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

	output *= 0.5*bendingrigidityconstant; // The thing is that the K_BT will be cancelled out anyway in the Metropolis rejection step
	return output;

}

// TODO: Put this wrapper function in the SWIG interface file
double getwlcenergy(double* coords, int numberofvertices, int ibe3, double bendingrigidityconstant)
{
	vector3d *vertices = reinterpret_cast<vector3d*>(coords);
	return _getwlcenergy(vertices, numberofvertices, bendingrigidityconstant);
}

inline bool ismonotonic(const double& a,const double& b,const double& c)
{
	return ((a<=b && b<=c) || (c<=b && b<=a));
}

inline int sgn(double val) {
   return (double(0) < val) - (val < double(0));
}


inline double _distancebetweenedges(vector3d * vertices, int numberofvertices, int _a1, int _a2, int _b1, int _b2)
// b1 and b2 will be copied by value to be able to modify them
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
		if (ismonotonic(dotprod(a1,m),dotprod(b1,m),dotprod(a2,m))) // ie, is in orthogonalspace in between
			b1 += m*dotprod(a-b1,m)/dotprod(m,m); // project onto orthogonal bisector
		else
			b1 += m/2*sgn(dotprod(a-b1,m)); // translate toward the projection

		// translate b2
		if (ismonotonic(dotprod(a1,m),dotprod(b2,m),dotprod(a2,m))) // ie, is in orthogonalspace in between
			b2 += m*dotprod(a-b2,m)/dotprod(m,m);
		else
			b2 += m/2*sgn(dotprod(a-b2,m));
	}

	vector3d b=(b1+b2)/2;
	if (b1!=b2)
	{
		vector3d n = b2-b1;

		// translate a
		if (ismonotonic(dotprod(b1,n),dotprod(a,n),dotprod(b2,n))) // ie, is in orthogonalspace in between
			a += n*dotprod(b-a,n)/dotprod(n,n);
		else
			a += n/2*sgn(dotprod(b-a,n));;
	}
	return norm(a-b);
}
double distancebetweenedges(double* coords, int numberofvertices, int ibe3, int _a1, int _a2, int _b1, int _b2)
// b1 and b2 will be copied by value to be able to modify them
{
	return _distancebetweenedges(reinterpret_cast<vector3d*>(coords),numberofvertices, _a1, _a2, _b1, _b2);
}

bool _iscollision(vector3d *vertices, int N, double diameter, double edgelength)
{
	if (diameter==0) return false;
	for (unsigned i=0; i<N-2; i++)
	{
		for (unsigned j=i+2; j<N-1; j++)
		{
			double gap = _distancebetweenedges(vertices,N,i,i+1,j,j+1)-diameter;
			if (gap < 0)
				return true;
			else
				// j+=1;
				j+=int((gap-diameter)/(edgelength));

		}
		if (i>0 && ( 0 > (_distancebetweenedges(vertices,N,i,i+1,N-1,0)-diameter)))
			return true;
	}

	return false;
}
bool iscollision(double* coords, int numberofvertices, int ibe3, double diameter, double edgelength)
{
	return _iscollision(reinterpret_cast<vector3d*>(coords),numberofvertices, diameter, edgelength);
}



inline int _domontecarlosteps(vector3d* vertices, int n, int numberofsteps, double bendingrigidityconstant, double thetamax, double diameter, double edgelength)
{

	int acceptances = 0;
	vector3d candidatevertices[n];

	for (unsigned i=0; i<numberofsteps; i++)
	{
		//This is the code for a single MC step
		memcpy(candidatevertices,vertices,n*sizeof(vector3d));

		_dorandomcrankshaft(candidatevertices,n,thetamax);
		double deltaenergy = 	_getwlcenergy(candidatevertices,n,bendingrigidityconstant) -
								_getwlcenergy(vertices,n,bendingrigidityconstant);

		if ( (deltaenergy<0 || exp(-deltaenergy)>rand_uniform()) && !_iscollision(candidatevertices,n,diameter,edgelength) ) 
		{
			// ACCEPT
			memcpy(vertices,candidatevertices,n*sizeof(vector3d));
			acceptances+=1;
		}
		else // REJECT
		{
		}

		// if (	(deltaenergy<0 || exp(-deltaenergy)>rand_uniform()) &&
		// 		!(_iscollision(vertices,n,diameter,edgelength))
		// 	) // ACCEPT
		// {
		// 	memcpy(vertices,candidatevertices,n*sizeof(vector3d));
		// 	acceptances+=1;
		// }
		// else // REJECT
		// {
		// }

	}

	return acceptances;

}

// TODO: Put this wrapper function in the SWIG interface file
int domontecarlosteps(double* coords, int numberofvertices, int ibe3, int numberofsteps, double bendingrigidityconstant, double thetamax, double diameter, double edgelength)
{
	vector3d *vertices = reinterpret_cast<vector3d*>(coords);
	return _domontecarlosteps(vertices, numberofvertices, numberofsteps, bendingrigidityconstant,thetamax,diameter,edgelength);
}


		// double EnergyChange = candidate->getenergy()-this->getenergy(); // the slowest part. Perhaps you could optimize within the class

		// // you can put debugging messages here.
		// if (verbose) std::cout << m << " " << n << " " << theta << "\t" << candidate->getenergy() - this->getenergy() << "\t";

void setrandomseedtoclocktime()
{
	set_sRand_seed_to_clocktime ();
}

int main(int argc, char const *argv[])
{

	double array[] = {0,0,1,0,0,-1,1,1,0,-1,1,0};


	cout << distancebetweenedges(array,4,3,0,1,0,1);
	// vector3d *vertices = reinterpret_cast<vector3d*>(array);

	// cout << vertices[0] << vertices[1] << vertices[2] << vertices[3] << endl;
	// domontecarlosteps(array,4,3,0,0,0);
	// cout << vertices[0] << vertices[1] << vertices[2] << vertices[3] << endl;


	// MTRand(37); // seed with 37
	// MTRand_int32 randint; // 32-bit int generator
	// MTRand rand01; // double in [0, 1) generator, already init'd

	// cout << randint() << " " << randint() << " " << randint() << " " << randint() << endl; 

	// string str;
	// cin >> str;

	// cout << *(str.c_str());
	// cstring cstr(str.c_str());

	// vector<double> vertices (&str,&str+sizeof(str) / sizeof(str[0]));
	
	// vector3d* vertices = reinterpret_cast<vector3d*>(&str)

	// vertices[0]

	// wormlikechain *wlc = new wormlikechain();
	// cout << wlc->vertices[0];

	// delete wlc;

 //    return 0;
}
