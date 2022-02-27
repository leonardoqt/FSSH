#include "potential.h"

using namespace arma;

void potential::init_H(double Omega, double G0, double Ed, double Gamma)
{
	omega = Omega;
	g0 = G0;
	ed = Ed;
	gamma = Gamma;
	Eb = linspace(-dep_bath,dep_bath,nbath);
}

mat potential::Hs(double x)
{
	mat HH(sz_s,sz_s,fill::zeros);
	HH(0,0) = sqrt(2)*g0*x + ed;
	return HH;
}

mat potential::Hf(double x)
{
	mat HH(sz_f,sz_f,fill::zeros);
	HH(span(1,sz_s),span(1,sz_s)) = Hs(x);
	return HH;
}

mat potential::Ht(double x)
{
	double vsb = sqrt( gamma/(datum::pi*nbath/dep_bath) );
	mat HH(sz_t,sz_t,fill::zeros);
	HH.diag() = join_vert(zeros<vec>(sz_s),Eb);
	HH(span(0,sz_s-1),span(0,sz_s-1)) = Hs(x);
	for(int t1 = sz_s; t1<sz_t; t1++)
	{
		HH(0,t1) = vsb;
		HH(t1,0) = vsb;
	}
	return HH;
}

void potential::ionic(double x, double& Eion, double& Fion)
{
	Eion = 0.5*omega*x*x;
	Fion = -omega*x;
}

void potential::adiab_H_f(double x, vec& E_f, vec& F_f, mat& U_f)
{
	vec e1,e2;
	eig_sym(e1,U_f,Hf(x-dx));
	eig_sym(e2,U_f,Hf(x+dx));
	F_f = (e1-e2)/2/dx;
	eig_sym(E_f,U_f,Hf(x));
}

void potential::adiab_H_s(double x, vec& E_s, mat& U_s)
{
	eig_sym(E_s,U_s,Hs(x));
}

cx_mat potential::ddt_f(double x1, double x2)
{
	mat U1, U2;
	vec E,F;
	adiab_H_f(x1,E,F,U1);
	adiab_H_f(x2,E,F,U2);
	for(int t1=0; t1<sz_f; t1++)
		if (dot(U1.col(t1),U2.col(t1)) < 0)
			U2.col(t1) *= -1;
	return logmat(U1.t()*U2);
}
