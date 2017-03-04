/*
 * MinloKeithDeclarations.h
 *
 *  Created on: 2 Feb 2016
 *      Author: daniel
 */

#ifndef MINLOKEITHDECLARATIONS_H_
#define MINLOKEITHDECLARATIONS_H_

// depending on the compiler and flags there is a trailing underscore
// or not.


#define I_PWHG_ALPHAS FC_UNDERSORE(i_pwhg_alphas)
#define LO_SUDAKOV FC_UNDERSORE(lo_sudakov)
#define EXPSUDAKOV FC_UNDERSORE(expsudakov)
#define SUDAKOV_EXPONENT FC_UNDERSORE(sudakov_exponent)
#define NLL_SUDAKOV FC_UNDERSORE(nll_sudakov)
#define I_SETLOCALSCALES FC_UNDERSORE(i_setlocalscales)
#define PDFFROMLHAPDF FC_UNDERSORE(pdffromlhapdf)
#define FC_UNDERSORE(fname) fname ## _

// can't be const so they can be passed to fortran
const double      KEITH_rad_charmthr2_CST=2.25;   //! Charm  (mass) threshold squared: (1.5 GeV)^2
const double      KEITH_rad_bottomthr2_CST=25.0;  //! Bottom (mass) threshold squared: (5.0 GeV)^2
const int         KEITH_st_nlight_CST=5;       //! Nominal number of light flavs in computation (5)
const double      KEITH_st_lambda5MSB_CST=0.226 /*same as in sudakov.f*/;   //! Lambda_5 in the MS bar scheme


extern "C" {

	double  LO_SUDAKOV(const  double &q20,const double &q2h,const double &q2l,const int &flav);
	double EXPSUDAKOV(const double &q20,const double &q2h,const double &q2l,const int &flav);

	double  SUDAKOV_EXPONENT(const  double &q20,const double &q2h,const double &q2l,const int &flav);
	double  NLL_SUDAKOV(const  double &q20,const double &q2h,const double &q2l,const int &flav);

	double I_SETLOCALSCALES(
//			C - New extra inputs for independence from includes
		     int*      flg_bornonly,    //! Are we feeding through only Born stuff, no NLO?
		     int&      imode,           //! imode=1 for Born, imode=2 for all NLO contribs
		     int&      isReal,           //! Set isReal=1 for real kinematics, 0 otherwise.
			     int&      nlegborn,      //  ! Number of final state particle in Born plus 2.
			     double*      kn_cmpborn,      //! Born particles momenta in C.O.M frame
			     int*      flav,            //! Born particle flavours (1,2) is incoming
			     int&      raisingscales,   //! Freeze clustering scale if it attempts to
			                            //! go down - no more sudakovs/subtractions, only
			                            //! coupling reweighting (for all remaining nodes).
			     double&      st_alpha,        //! Fixed aS used in MEs before they got here
			     double&      st_muren2,       //! The mu_R^2 value that was used to evaluate virt
			     int&      st_bornorder,    //! Number of powers of aS in the Born (ZJ=1, ZJJ=2)
			     double&      rad_charmthr2,   //! Charm  (mass) threshold squared: (1.5 GeV)^2
			     double&      rad_bottomthr2,  //! Bottom (mass) threshold squared: (5.0 GeV)^2
			     int&      st_nlight,       //! Nominal number of light flavs in computation (5)
			     double&      st_lambda5MSB,   //! Lambda_5 in the MS bar scheme
			     double&      st_renfact,      //! Scaling factor for mu_R (not mu_R^2): 0.5,1,2.0
			     double&      st_facfact,      //! Scaling factor for mu_F (not mu_F^2): 0.5,1,2.0
			     double&      rescfac,         //! Comb. of suds, subtractions & aS weights to
			                            //! rescale the weight of the contrib assoc. to this
			                            //! mom. and flav. configuration.
			     double&      st_mufact2,      //! The mu_F^2 value that the PDFs should be
			                             //! re-evaluated at now
//			C - Old inputs and output
			     double&      basicfac,        //! Comb. of suds, subtractions & aS weights to
			     double&      bornfac,         //! Comb. of suds, subtractions & aS weights to
			     double&      nlofac,          //! Comb. of suds, subtractions & aS weights to
//			C - Extra factorisation scale output as it normally lives in pwhg common block
//////////////// not used any more
				 //			C - Momenta tha tare only used for debugging
//				 			     double*      kn_pborn,        //! Born particles momenta in LAB frame
int& alltheway
	);


  double FC_UNDERSORE(i_pwhg_alphas)(double &Q2,double &st_lambda5MSB, int &st_nlight,double &rad_charmthr2,double &rad_bottomthr2);

}




#endif /* MINLOKEITHDECLARATIONS_H_ */
