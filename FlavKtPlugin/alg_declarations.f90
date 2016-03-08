module alg_declarations
  implicit none 
  private 
  
  ! different types of jet-algorithms... 
  integer, public, parameter :: Durham    = 1
  integer, public, parameter :: Jade      = 2
  integer, public, parameter :: Geneva    = 3
  integer, public, parameter :: Durham_flav   = 4 ! flavour with E suprsn
  integer, public, parameter :: Durham_flav2  = 5 ! flavour with E^2 suprsn
  integer, public, parameter :: Durham_flavs  = 6 ! flav selective, E sup
  integer, public, parameter :: Durham_flavs2 = 7 ! flav selective, E^2 sup
  integer, public, parameter :: Durham_s      = 8 ! flav selective plain Durham
  integer, public, parameter :: AODurh_flav   = 9 ! flavour with E suprsn
  integer, public, parameter :: AODurh_flav2  = 10! flavour with E^2 suprsn
  integer, public, parameter :: AODurh_flavs  = 11! flav selective, E sup
  integer, public, parameter :: AODurh_flavs2 = 12! flav selective, E^2 sup
  integer, public, parameter :: AODurh_s      = 13! AO with flav.selection

  ! different types of recombinations 
  integer, public, parameter :: Pscheme  = 1
  integer, public, parameter :: P0scheme = 2
  integer, public, parameter :: Escheme  = 3
  integer, public, parameter :: E0scheme = 4

  ! different types of flavour recombination schemes
  integer, public, parameter :: FlavScheme_expb  = 1 
  integer, public, parameter :: FlavScheme_net   = 2
  integer, public, parameter :: FlavScheme_mod2  = 3


  !---- HH ---------------------------------------
  integer, public, parameter :: itype_ee = 1000 !! e+e- type kt-algo
  integer, public, parameter :: itype_pp = 4000 !! pp type kt-algo
  integer, public, parameter :: iangle_ang    = 100 !! min(E_i,E_j)^2 (1-cos theta_ij)
  integer, public, parameter :: iangle_DeltaR = 200 !! min(pt_i,pt_j)^2 (Delta R_ij)^2
  ! flavour variants of DeltaR
  integer, public, parameter :: iangle_DeltaR_flav1  = 300
  integer, public, parameter :: iangle_DeltaR_flav2  = 400
  integer, public, parameter :: iangle_DeltaR_flavs1 = 500
  integer, public, parameter :: iangle_DeltaR_flavs2 = 600
  integer, public, parameter :: iangle_DeltaRs       = 700
  ! different types of recombinations 
  integer, public, parameter :: irecom_P   = 1
  integer, public, parameter :: irecom_P0  = 2
  integer, public, parameter :: irecom_E   = 3
  integer, public, parameter :: irecom_E0  = 4
  ! new recombinations for the kt clustering algorithm
  integer, public, parameter :: irecom_Pt  = 5
  integer, public, parameter :: irecom_Pt2 = 6
  integer, public, parameter :: irecom_Et  = 7
  integer, public, parameter :: irecom_Et2 = 8

  ! some shortcuts...
  ! NB: all declarations must stay on a single line to help f90->C++ 
  !     conversions.
  integer,public,parameter :: ialg_pp_ktLI_E=itype_pp+iangle_DeltaR+irecom_E
  integer,public,parameter :: ialg_pp_ktf1_E=itype_pp+iangle_DeltaR_flav1+ irecom_E
  integer,public,parameter :: ialg_pp_ktf2_E=itype_pp+iangle_DeltaR_flav2+ irecom_E
  integer,public,parameter :: ialg_pp_ktfs1_E=itype_pp+iangle_DeltaR_flavs1+ irecom_E
  integer,public,parameter :: ialg_pp_ktfs2_E=itype_pp+iangle_DeltaR_flavs2+ irecom_E
  integer,public,parameter :: ialg_ee_kt_E  =itype_ee+iangle_ang+irecom_E
  integer,public,parameter :: ialg_ee_kt_P  =itype_ee+iangle_ang+irecom_P
  integer,public,parameter :: ialg_ee_kt_E0  =itype_ee+iangle_ang+irecom_E0


  ! for variants of the flavs2 algs:
  integer, parameter, public :: modelA = 10000
  integer, parameter, public :: modelB = 10001
  integer, parameter, public :: modelC = 10002
  integer, parameter, public :: modelD = 10003
  integer, parameter, public :: modelE = 10004
  integer, parameter, public :: modelF = 10005
  integer, parameter, public :: modelG = 10006
  integer, parameter, public :: modelH = 10007
  integer, parameter, public :: modelI = 10008
  integer, parameter, public :: modelJ = 10009
  integer, public :: flavs2_model = modelA
  
  public :: flav_add

contains
  
  !-----------------------------------------------------------
  !! Given flavour arrays a,b return the sum of the combination according
  !! to flav_scheme
  elemental function flav_add(a,b,flav_scheme) result(res)
    integer, intent(in) :: a, b, flav_scheme
    integer             :: res
    !-----------------------------------------
    select case(flav_scheme)
    case(FlavScheme_expb)
       res = min(1,abs(a)+abs(b))
    case(FlavScheme_net)
       res = a + b
    case(FlavScheme_mod2)
       res = abs(mod(a+b,2))
    case default
       res = 0
    end select
  end function flav_add
  

end module alg_declarations
