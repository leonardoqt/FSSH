#include "ionic.h"

using namespace arma;

void ionic::init(double Mass, double X, double V, int state, double Dt, int Period, double Xendl, double Xendr)
{
	dt = Dt;
	hop_period = Period;
	mass = Mass;
	x_t2 = x_t1 = x = X;
	v = V;
	ek = mass*v*v/2;
	istate = state;
	xendl = Xendl;
	xendr = Xendr;
	nhops = 0;
}

void ionic::move(potential& HH)
{
	// Velocity Verlet
	vec E_f, F_f;
	mat U;
	double Eion, Fion, a1, a2;
	//
	x_t2 = x_t1;
	x_t1 = x;
	//
	HH.adiab_H_f(x,E_f,F_f,U);
	HH.ionic(x,Eion,Fion);
	a1 = (F_f(istate)+Fion) / mass;
	x += v*dt + a1*dt*dt/2;
	//
	HH.adiab_H_f(x,E_f,F_f,U);
	HH.ionic(x,Eion,Fion);
	a2 = (F_f(istate)+Fion) / mass;
	v += (a1+a2)/2*dt;
	ek = mass*v*v/2;
	//cout<<x<<'\t'<<v<<'\t'<<a1<<'\t'<<a2<<endl;
}

void ionic::try_hop(potential &HH, arma::cx_mat &rho, arma::mat &hop_bath)
{
	cx_mat T = HH.ddt_f(x_t2, x);
	vec rate_s(HH.sz_f), rate_b(HH.sz_f);
	//
	// TODO: rethink about this part
	//*************
	/*
	// rho_ii_dot can be splitted into two terms, the ration decides whether it undergoes the derivative coupling procedue or relaxation procedure
	cx_mat rho_dot1 = (T*rho - rho*T)/dt;
	double rho_ii_dot1 = real(rho_dot1(istate,istate));
	double rho_ii_dot2 = sum(hop_bath.row(istate))*real(rho(istate,istate));
	*/
	//*************
	//
	double Eion,dtmp;
	vec Ef,Ftmp;
	mat Utmp;
	HH.adiab_H_f(x,Ef,Ftmp,Utmp);
	HH.ionic(x,Eion,dtmp);
	for (int t1=0; t1<HH.sz_f; t1++)
	{
		if (t1==istate)
		{
			rate_s(t1) = 0;
			rate_b(t1) = 0;
		}
		else
		{
			rate_s(t1) = real( T(istate,t1)*rho(t1,istate) ) * 2 / real(rho(istate,istate));
			rate_b(t1) = ( hop_bath(istate,t1)*real(rho(istate,istate))-hop_bath(t1,istate)*real(rho(t1,t1)) )/real(rho(istate,istate)) * dt*hop_period;
		}
		if (rate_s(t1) < 0)
			rate_s(t1) = 0;
		if (rate_b(t1) < 0)
			rate_b(t1) = 0;
		if (ek + Ef(istate) < Ef(t1))
			rate_s(t1) = 0;
	}
	//rate_s.t().print();
	//rate_b.t().print();
	//cout<<endl;
	//
	vec rate = join_vert(rate_s,rate_b);
	for (int t1=1; t1<HH.sz_f*2; t1++)
		rate(t1) += rate(t1-1);
	//
	vec tmp(1,fill::randu);
	int new_state = istate;
	int from_bath = 0;
	for (int t1=0; t1<HH.sz_f*2; t1++)
		if( tmp(0) < rate(t1) )
		{
			new_state = t1 % HH.sz_f;
			from_bath = t1 / HH.sz_f;
			break;
		}
	//
	// adjust velocity
	if (from_bath == 0)
	{
		double ek_new = ek + Ef(istate) - Ef(new_state);
		v = v * sqrt(ek_new / ek);
		ek = ek_new;
	}
	if (istate != new_state) nhops++;
	istate = new_state;
	etot = ek + Ef(istate) + Eion;
}

int ionic::check_stop()
{
	if (x < xendl)
		return -1;
	else if (x > xendr)
		return 1;
	else
		return 0;
}
