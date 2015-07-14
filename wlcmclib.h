#ifndef __WLCMC_H__
#define __WLCMC_H__

void docrankshaft(double* coords, int numberofvertices, int ibe3, int vert1, int numtonextvert, double theta);
void dorandomcrankshaft(double* coords, int numberofvertices, int ibe3, double thetamax);
double getwlcenergy(double* coords, int numberofvertices, int ibe3, double bendingrigidityconstant);
int domontecarlosteps(double* coords, int numberofvertices, int ibe3, int numberofsteps, double bendingrigidityconstant, double thetamax, double diameter, double edgelength);
double distancebetweenedges(double* coords, int numberofvertices, int ibe3, int _a1, int _a2, int _b1, int _b2);
bool iscollision(double* coords, int numberofvertices, int ibe3, double diameter, double edgelength);
void setrandomseedtoclocktime();

#endif 