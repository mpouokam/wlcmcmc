#ifndef __WLCMCMC_H__
#define __WLCMCMC_H__

void do_crank_shaft(double* coords, int numberofvertices, int ibe3, int vert1, int numtonextvert, double theta);
void do_random_crank_shaft(double* coords, int numberofvertices, int ibe3, double thetamax);
double get_wlc_energy(double* coords, int numberofvertices, int ibe3, double bendingrigidityconstant);
int do_monte_carlo_steps(double* coords, int numberofvertices, int ibe3, int numberofsteps, double bendingrigidityconstant, double thetamax, double diameter, double edgelength);
double distance_between_edges(double* coords, int numberofvertices, int ibe3, int _a1, int _a2, int _b1, int _b2);
bool is_collision(double* coords, int numberofvertices, int ibe3, double diameter, double edgelength);
void set_random_seed();

#endif 