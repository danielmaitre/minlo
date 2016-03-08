// File generated automatically by ./GenIModeIncludes.pl 
// from: ../../common/alg_declarations.f90.
// Provides definitions of IDs for all jet algorithms
#ifndef __JETIMODES_H__
#define __JETIMODES_H__
const int itype_ee = 1000 ; // e+e- type kt-algo;
const int itype_pp = 4000 ; // pp type kt-algo;
const int iangle_ang    = 100 ; // min(E_i,E_j)^2 (1-cos theta_ij);
const int iangle_DeltaR = 200 ; // min(pt_i,pt_j)^2 (Delta R_ij)^2;
const int iangle_DeltaR_flav1  = 300;
const int iangle_DeltaR_flav2  = 400;
const int iangle_DeltaR_flavs1 = 500;
const int iangle_DeltaR_flavs2 = 600;
const int iangle_DeltaRs       = 700;
const int irecom_P   = 1;
const int irecom_P0  = 2;
const int irecom_E   = 3;
const int irecom_E0  = 4;
const int irecom_Pt  = 5;
const int irecom_Pt2 = 6;
const int irecom_Et  = 7;
const int irecom_Et2 = 8;
const int ialg_pp_ktLI_E=itype_pp+iangle_DeltaR+irecom_E;
const int ialg_pp_ktf1_E=itype_pp+iangle_DeltaR_flav1+ irecom_E;
const int ialg_pp_ktf2_E=itype_pp+iangle_DeltaR_flav2+ irecom_E;
const int ialg_pp_ktfs1_E=itype_pp+iangle_DeltaR_flavs1+ irecom_E;
const int ialg_pp_ktfs2_E=itype_pp+iangle_DeltaR_flavs2+ irecom_E;
const int ialg_ee_kt_E  =itype_ee+iangle_ang+irecom_E;
const int ialg_ee_kt_P  =itype_ee+iangle_ang+irecom_P;
const int ialg_ee_kt_E0  =itype_ee+iangle_ang+irecom_E0;
#endif   // __JETIMODES_H__
