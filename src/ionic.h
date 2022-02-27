#ifndef __IONIC__
#define __IONIC__

#include <armadillo>
#include "potential.h"

class ionic;

class ionic
{
public:
	int istate;
	int nhops;
	double dt;
	int hop_period;
	double mass, ek, etot, x, v;
	double x_t1, x_t2;
	double xendl, xendr;
	//
	void init(double Mass, double X, double V, int state, double Dt, int Period, double Xendl, double Xendr);
	void move(potential& HH);
	void try_hop(potential& HH, arma::cx_mat& rho, arma::mat& hop_bath);
	int check_stop();
};

#endif
