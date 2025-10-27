#include <armadillo>
#include "statespace.hpp"
#include "pdf.hpp"
using namespace arma;

/** Generating model for example 1 **/
extern "C" colvec SIRB(const colvec & parameters , const colvec & state, const double time, const double dt ) {
	const double gamma= 1.0/14.0;
	const double q = 1.0;
	
	double rho = parameters[0];

	double sk = state[0];
	double ik = state[1];
	double rk = state[2];
	double betak = state[3];

	colvec temp = state;

	temp[0] = sk + dt*(-betak*sk*ik);	 
	temp[1] = ik + dt*(betak*sk*ik - gamma*ik);
	temp[2] = rk + dt*(gamma*ik);
	temp[3] = betak;                                                                                                              	//
	return temp;
}

static double staticTime;

/*************************************************************
**************** SIRS Model Fully Connected  ***************
*************************************************************/

colvec modelSIRS_Full( const colvec & state, const double dt, double time, const colvec & parameters ) {
	const double alpha = parameters[0];
	const double gamma = parameters[1];
	const double delta = parameters[2];
	const double q = parameters[3];
	const double beta0 = parameters[4];
		
	double sk = state[0];
	double ik = state[1];
	double rk = state[2];
	double betak = state[3];
	
	if (time == 0 ) {
		betak = beta0;
	}
	
	colvec temp = state;
		
	temp[0] = sk + dt*(-betak*sk*ik + alpha*ik + delta*rk);	 
	temp[1] = ik + dt*(betak*sk*ik - gamma*ik - alpha*ik);
	temp[2] = rk + dt*(gamma*ik - delta*rk);
	temp[3] = betak;                                                                                                              	//                                                                              	//
	return temp;
}

colvec hSIRS_Full = zeros<colvec>(1);
colvec _hSIRS_Full( const colvec & state, const mat & covariance, const colvec & parameters ) {
	double ik = state[1];

	hSIRS_Full[0] = ik;
	return hSIRS_Full;
}

// Declare all the jacobians for EKF
mat dfdxSIRS_Full (const colvec & state, const double dt, const colvec & parameters ) {
	const double alpha = parameters[0];
	const double gamma = parameters[1];
	const double delta = parameters[2];
	const double q = parameters[3];
	const double beta0 = parameters[4];
		
	double sk = state[0];
	double ik = state[1];
	double rk = state[2];
	double betak = state[3];
	
	if (time == 0 ) {
		betak = beta0;
	}
	
	colvec temp = state;      

	return {{1.0 + dt*(-betak*ik), dt*(-betak*sk + alpha), dt*(delta), dt*(-sk*ik)},
		{dt*(betak*ik), 1.0 + dt*(betak*sk - gamma - alpha), 0.0, dt*(sk*ik)},
		{0.0, dt*gamma, 1.0 - dt*delta, 0.0},
		{0.0, 0.0, 0.0, 1.0}};
}

mat dfdeSIRS_Full (const colvec & state, const double dt, const colvec & parameters ) {
	double q = parameters[3];
	
	return {{0.0, 0.0, 0.0, 0.0},
		{0.0, 0.0, 0.0, 0.0},
		{0.0, 0.0, 0.0, 0.0},
		{q*std::sqrt(dt), 0.0, 0.0, 0.0}};
}	

mat dhdxSIRS_Full (const colvec& , const colvec & parameters ) {
	return {0.0,1.0,0.0,0.0};

}

mat dhdeSIRS_Full (const colvec & state, const colvec & parameters ) {
	double ik = state[1];
	
	return {ik + 1.0e-10};
}


/*************************************************************
******************** Nested SIRS Model *********************
*************************************************************/

colvec modelSIRS_SIRS( const colvec & state, const double dt, double time, const colvec & parameters ) {
	const double gamma = parameters[0];
	const double delta = parameters[1];
	const double q = parameters[2];
	const double beta0 = parameters[3];
		
	double sk = state[0];
	double ik = state[1];
	double rk = state[2];
	double betak = state[3];
	
	if (time == 0 ) {
		betak = beta0;
	}
	
	colvec temp = state;
		
	temp[0] = sk + dt*(-betak*sk*ik + delta*rk);	 
	temp[1] = ik + dt*(betak*sk*ik - gamma*ik);
	temp[2] = rk + dt*(gamma*ik - delta*rk);
	temp[3] = betak;                                                                                                              	                                                                         	//
	return temp;
}

colvec hSIRS_SIRS = zeros<colvec>(1);
colvec _hSIRS_SIRS( const colvec & state, const mat & covariance, const colvec & parameters ) {
	double ik = state[1];

	hSIRS_SIRS[0] = ik;
	return hSIRS_SIRS;
}

mat dfdxSIRS_SIRS (const colvec & state, const double dt, const colvec & parameters ) {
	const double gamma = parameters[0];
	const double delta = parameters[1];
	const double q = parameters[2];
	const double beta0 = parameters[3];
		
	double sk = state[0];
	double ik = state[1];
	double rk = state[2];
	double betak = state[3];
	
	if (time == 0 ) {
		betak = beta0;
	}
	
	colvec temp = state;      

	return {{1.0 + dt*(-betak*ik), dt*(-betak*sk), dt*(delta), dt*(-sk*ik)},
		{dt*(betak*ik), 1.0 + dt*(betak*sk - gamma), 0.0, dt*(sk*ik)},
		{0.0, dt*gamma, 1.0 - dt*delta, 0.0},
		{0.0, 0.0, 0.0, 1.0}};
}

mat dfdeSIRS_SIRS (const colvec & state, const double dt, const colvec & parameters ) {
	double q = parameters[2];
	
	return {{0.0, 0.0, 0.0, 0.0},
		{0.0, 0.0, 0.0, 0.0},
		{0.0, 0.0, 0.0, 0.0},
		{q*std::sqrt(dt), 0.0, 0.0, 0.0}};
}	

mat dhdxSIRS_SIRS (const colvec& , const colvec & parameters ) {
	return {0.0,1.0,0.0,0.0};

}

mat dhdeSIRS_SIRS (const colvec & state, const colvec & parameters ) {
	double ik = state[1];
	
	return {ik + 1.0e-10};
}


/*************************************************************
******************** Nested SIR Model **********************
*************************************************************/

colvec modelSIRS_SIR( const colvec & state, const double dt, double time, const colvec & parameters ) {
	const double gamma = parameters[0];
	const double q = parameters[1];
	const double beta0 = parameters[2];
		
	double sk = state[0];
	double ik = state[1];
	double rk = state[2];
	double betak = state[3];
	
	if (time == 0 ) {
		betak = beta0;
	}
	
	colvec temp = state;
		
	temp[0] = sk + dt*(-betak*sk*ik);	 
	temp[1] = ik + dt*(betak*sk*ik - gamma*ik);
	temp[2] = rk + dt*(gamma*ik);
	temp[3] = betak;                                                                                                              	                                                                            	//
	return temp;
}

colvec hSIRS_SIR = zeros<colvec>(1);
colvec _hSIRS_SIR( const colvec & state, const mat & covariance, const colvec & parameters ) {
	double ik = state[1];

	hSIRS_SIR[0] = ik;
	return hSIRS_SIR;
}

mat dfdxSIRS_SIR (const colvec & state, const double dt, const colvec & parameters ) {
	const double gamma = parameters[0];
	const double q = parameters[1];
	const double beta0 = parameters[2];
		
	double sk = state[0];
	double ik = state[1];
	double rk = state[2];
	double betak = state[3];
	
	if (time == 0 ) {
		betak = beta0;
	}
	
	colvec temp = state;      

	return {{1.0 + dt*(-betak*ik), dt*(-betak*sk), 0.0, dt*(-sk*ik)},
		{dt*(betak*ik), 1.0 + dt*(betak*sk - gamma), 0.0, dt*(sk*ik)},
		{0.0, dt*gamma, 1.0, 0.0},
		{0.0, 0.0, 0.0, 1.0}};
}

mat dfdeSIRS_SIR (const colvec & state, const double dt, const colvec & parameters ) {
	double q = parameters[1];
	
	return {{0.0, 0.0, 0.0, 0.0},
		{0.0, 0.0, 0.0, 0.0},
		{0.0, 0.0, 0.0, 0.0},
		{q*std::sqrt(dt), 0.0, 0.0, 0.0}};
}	

mat dhdxSIRS_SIR (const colvec& , const colvec & parameters ) {
	return {0.0,1.0,0.0,0.0};

}

mat dhdeSIRS_SIR (const colvec & state, const colvec & parameters ) {
	double ik = state[1];
	
	return {ik + 1.0e-10};
}

/*************************************************************
******************** Nested SIS Model ****I*****************
*************************************************************/

colvec modelSIRS_SIS( const colvec & state, const double dt, double time, const colvec & parameters ) {
	const double alpha = parameters[0];
	const double q = parameters[1];
	const double beta0 = parameters[2];
		
	double sk = state[0];
	double ik = state[1];
	double betak = state[2];
	
	if (time == 0 ) {
		betak = beta0;
	}
	
	colvec temp = state;
		
	temp[0] = sk + dt*(-betak*sk*ik + alpha*ik);	 
	temp[1] = ik + dt*(betak*sk*ik - alpha*ik);
	temp[2] = betak;                                                                                                              	                                                                         	//
	return temp;
}

colvec hSIRS_SIS = zeros<colvec>(1);
colvec _hSIRS_SIS( const colvec & state, const mat & covariance, const colvec & parameters ) {
	double ik = state[1];

	hSIRS_SIS[0] = ik;
	return hSIRS_SIS;
}

mat dfdxSIRS_SIS (const colvec & state, const double dt, const colvec & parameters ) {
	const double alpha = parameters[0];
	const double q = parameters[1];
	const double beta0 = parameters[2];
		
	double sk = state[0];
	double ik = state[1];
	double betak = state[2];
	
	if (time == 0 ) {
		betak = beta0;
	}
	
	colvec temp = state;      

	return {{1.0 + dt*(-betak*ik), dt*(-betak*sk + alpha), dt*(-sk*ik)},
		{dt*(betak*ik), 1.0 + dt*(betak*sk - alpha), dt*(sk*ik)},
		{0.0, 0.0, 1.0}};
}

mat dfdeSIRS_SIS (const colvec & state, const double dt, const colvec & parameters ) {
	double q = parameters[1];
	
	return {{0.0, 0.0, 0.0},
		{0.0, 0.0, 0.0},
		{q*std::sqrt(dt), 0.0, 0.0}};
}	

mat dhdxSIRS_SIS (const colvec& , const colvec & parameters ) {
	return {0.0,1.0,0.0};

}

mat dhdeSIRS_SIS (const colvec & state, const colvec & parameters ) {
	double ik = state[1];
	
	return {ik + 1.0e-10};
}

/*************************************************************
******************* Build all the statespaces*****************
*************************************************************/

extern "C" statespace mSIRS_SIR = statespace(modelSIRS_SIR, dfdxSIRS_SIR, dfdeSIRS_SIR, _hSIRS_SIR, dhdxSIRS_SIR, dhdeSIRS_SIR, 4, 1);
extern "C" statespace mSIRS_SIRS = statespace(modelSIRS_SIRS, dfdxSIRS_SIRS, dfdeSIRS_SIRS, _hSIRS_SIRS, dhdxSIRS_SIRS, dhdeSIRS_SIRS, 4, 1);
extern "C" statespace mSIRS_SIS = statespace(modelSIRS_SIS, dfdxSIRS_SIS, dfdeSIRS_SIS, _hSIRS_SIS, dhdxSIRS_SIS, dhdeSIRS_SIS, 3, 1);
extern "C" statespace mSIRS_Full = statespace(modelSIRS_Full, dfdxSIRS_Full, dfdeSIRS_Full, _hSIRS_Full, dhdxSIRS_Full, dhdeSIRS_Full, 4, 1);
