#ifndef __ELECTRONIC_
#define __ELECTRONIC_

#include <armadillo>
#include "potential.h"

class electronic;

class electronic
{
private:
	arma::mat H_t_t1, H_t_t2, H_t;
	int not_evolve;
	double x_2, x_1, x_0;
public:
	int sz_s, sz_t, sz_f;
	// N_t in diabats, N_s in adiabats, rho_fock in adiabats. N_t and N_s are really Tr[rho d_i\dagger d_j]
	arma::cx_mat N_t, N_s, rho_fock, rho_fock_old; 
	arma::cx_mat drho, drho_2fit;
	// TODO: dt here is twice the step of dt in ionic
	// need to run two ionic steps then run electronic step
	double beta, damping, dt;
	arma::vec N_b;
	arma::mat hop_bath;
	//
	void init_rho(arma::mat N0_s, potential& HH, double Beta, double Damping, double Dt);
	void evolve(potential& HH, double x_last2, double x_last, double x_now);
	void fit_drho(potential& HH, int method);
	//void try_decoherence(ionic& AA);
private:
	void construct_rho_fock();
	arma::cx_mat N_dot(arma::mat H0, arma::cx_mat N0);
	//void fit_drho_v1(potential& HH); // impose detailed balance, fit only diagonal
	//void fit_drho_v2(potential& HH); // impose detailed balance, fit full matrix
	void fit_drho_v3(); // impose single orbital the same, fit only diagonal
	void fit_drho_v3_1imp(); // impose single orbital the same, fit only diagonal
};
#endif
