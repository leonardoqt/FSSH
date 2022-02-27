#ifndef __POTENTIAL__
#define __POTENTIAL__

#include <armadillo>

class potential; // time independent part

class potential
{
private:
	const int nbath = 100;
	const double dep_bath = 0.3;
	const double dx = 1e-3;
public:
	const int sz_s = 1;
	const int sz_f = 1<<sz_s; // this is 2^sz_s
	const int sz_t = sz_s + nbath;
	arma::vec Eb;
	//
	double omega, g0, ed, gamma;
	void init_H(double Omega, double G0, double Ed, double Gamma);
	arma::mat Hs(double x);
	arma::mat Hf(double x);
	arma::mat Ht(double x);
	void ionic(double x, double& Eion, double& Fion);
	void adiab_H_f(double x, arma::vec& E_f, arma::vec& F_f, arma::mat& U_f);
	void adiab_H_s(double x, arma::vec& E_s, arma::mat& U_s);
	arma::cx_mat ddt_f(double x1, double x2);
	//
};

#endif
