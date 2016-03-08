
!!
!! Module generated automatically by NameSelect.pl
!! Provides conversion from string codes to integer codes
!!
module NameSelect
  use alg_declarations
  
  implicit none
  private

  public :: CodeOfName
  public :: NameOfCode
  public :: code_val_opt

  integer, parameter :: maxlen_name = 40
  integer, parameter :: maxlen_longname = 120
  integer, parameter, public :: NameSelect_maxlen_name = maxlen_name
  integer, parameter, public :: NameSelect_maxlen_longname = maxlen_longname

contains

  !!
  !! Returns the code associated with 'name'; if status dummy arg
  !! is present a non-zero result indicates failure to find the name.
  !! Otherwise failure results in a hard error.
  !!
  function CodeOfName(name, status) result(code)
    character(len=*), intent(in)   :: name
    integer                        :: code
    integer, optional, intent(out) :: status

    if (present(status)) status = 0

    ! now follows code generated specifically for above used modules
    select case(trim(name))
    case('Durham')
      code = Durham
    case('Jade')
      code = Jade
    case('Geneva')
      code = Geneva
    case('Durham_flav')
      code = Durham_flav
    case('Durham_flav2')
      code = Durham_flav2
    case('Durham_flavs')
      code = Durham_flavs
    case('Durham_flavs2')
      code = Durham_flavs2
    case('Durham_s')
      code = Durham_s
    case('AODurh_flav')
      code = AODurh_flav
    case('AODurh_flav2')
      code = AODurh_flav2
    case('AODurh_flavs')
      code = AODurh_flavs
    case('AODurh_flavs2')
      code = AODurh_flavs2
    case('AODurh_s')
      code = AODurh_s
    case('Pscheme')
      code = Pscheme
    case('P0scheme')
      code = P0scheme
    case('Escheme')
      code = Escheme
    case('E0scheme')
      code = E0scheme
    case('FlavScheme_expb')
      code = FlavScheme_expb
    case('FlavScheme_net')
      code = FlavScheme_net
    case('FlavScheme_mod2')
      code = FlavScheme_mod2
    case('itype_ee')
      code = itype_ee
    case('itype_pp')
      code = itype_pp
    case('iangle_ang')
      code = iangle_ang
    case('iangle_DeltaR')
      code = iangle_DeltaR
    case('iangle_DeltaR_flav1')
      code = iangle_DeltaR_flav1
    case('iangle_DeltaR_flav2')
      code = iangle_DeltaR_flav2
    case('iangle_DeltaR_flavs1')
      code = iangle_DeltaR_flavs1
    case('iangle_DeltaR_flavs2')
      code = iangle_DeltaR_flavs2
    case('iangle_DeltaRs')
      code = iangle_DeltaRs
    case('irecom_P')
      code = irecom_P
    case('irecom_P0')
      code = irecom_P0
    case('irecom_E')
      code = irecom_E
    case('irecom_E0')
      code = irecom_E0
    case('irecom_Pt')
      code = irecom_Pt
    case('irecom_Pt2')
      code = irecom_Pt2
    case('irecom_Et')
      code = irecom_Et
    case('irecom_Et2')
      code = irecom_Et2
    case('ialg_pp_ktLI_E')
      code = ialg_pp_ktLI_E
    case('ialg_pp_ktf1_E')
      code = ialg_pp_ktf1_E
    case('ialg_pp_ktf2_E')
      code = ialg_pp_ktf2_E
    case('ialg_pp_ktfs1_E')
      code = ialg_pp_ktfs1_E
    case('ialg_pp_ktfs2_E')
      code = ialg_pp_ktfs2_E
    case('ialg_ee_kt_E')
      code = ialg_ee_kt_E
    case('ialg_ee_kt_P')
      code = ialg_ee_kt_P
    case('ialg_ee_kt_E0')
      code = ialg_ee_kt_E0
    case('modelA')
      code = modelA
    case('modelB')
      code = modelB
    case('modelC')
      code = modelC
    case('modelD')
      code = modelD
    case('modelE')
      code = modelE
    case('modelF')
      code = modelF
    case('modelG')
      code = modelG
    case('modelH')
      code = modelH
    case('modelI')
      code = modelI
    case('modelJ')
      code = modelJ
    case default
      if (present(status)) then
        status = -1
      else
        write(0,*) 'ERROR. CodeOfName: unrecognized name "'//name//'"'
        stop
      endif
    end select
  end function CodeOfName


  ! standard code (independent of used modules)
  !!
  !! looks for the command-line argument 'option' and if it is
  !! present its value (optionally prefixed with 'prefix') is
  !! fed to CodeOfName
  !!
  function code_val_opt(option, default, prefix) result(code)
    use sub_defs_io
    character(len=*),           intent(in) :: option
    integer, optional,          intent(in) :: default
    character(len=*), optional, intent(in) :: prefix
    integer                                :: code
    !---------------------------------------
    character(len=maxlen_name) :: opt_val

    if (log_val_opt(option)) then
       opt_val = string_val_opt(option)
       if (present(prefix)) then
         code = CodeOfName(prefix//trim(opt_val))
       else
         code = CodeOfName(trim(opt_val))
       end if
    else
       if (present(default)) then
          code = default
       else
          write(0,*) 'Error in code_val_opt: command-line option '&
               &//option//' absent and no default provided'
       end if
    end if
  end function code_val_opt

  !!
  !! Returns the name of the given integer code which has
  !! the (optional) specified prefix
  !!
  function NameOfCode(code,prefix,longname) result(name)
    integer,          intent(in)            :: code
    character(len=*), intent(in),  optional :: prefix
    character(len=*), intent(out), optional :: longname
    character(len=maxlen_name)              :: name
    !----------------------------------------------
    integer :: nocc

    nocc = 0
    if (present(longname)) longname = ''

    ! code specific to used modules starts here
    if (PrefixMatches('Durham',prefix) .and. code == Durham) then
       nocc = nocc + 1; name = 'Durham'
    end if

    if (PrefixMatches('Jade',prefix) .and. code == Jade) then
       nocc = nocc + 1; name = 'Jade'
    end if

    if (PrefixMatches('Geneva',prefix) .and. code == Geneva) then
       nocc = nocc + 1; name = 'Geneva'
    end if

    if (PrefixMatches('Durham_flav',prefix) .and. code == Durham_flav) then
       nocc = nocc + 1; name = 'Durham_flav'
    end if

    if (PrefixMatches('Durham_flav2',prefix) .and. code == Durham_flav2) then
       nocc = nocc + 1; name = 'Durham_flav2'
    end if

    if (PrefixMatches('Durham_flavs',prefix) .and. code == Durham_flavs) then
       nocc = nocc + 1; name = 'Durham_flavs'
    end if

    if (PrefixMatches('Durham_flavs2',prefix) .and. code == Durham_flavs2) then
       nocc = nocc + 1; name = 'Durham_flavs2'
    end if

    if (PrefixMatches('Durham_s',prefix) .and. code == Durham_s) then
       nocc = nocc + 1; name = 'Durham_s'
    end if

    if (PrefixMatches('AODurh_flav',prefix) .and. code == AODurh_flav) then
       nocc = nocc + 1; name = 'AODurh_flav'
    end if

    if (PrefixMatches('AODurh_flav2',prefix) .and. code == AODurh_flav2) then
       nocc = nocc + 1; name = 'AODurh_flav2'
    end if

    if (PrefixMatches('AODurh_flavs',prefix) .and. code == AODurh_flavs) then
       nocc = nocc + 1; name = 'AODurh_flavs'
    end if

    if (PrefixMatches('AODurh_flavs2',prefix) .and. code == AODurh_flavs2) then
       nocc = nocc + 1; name = 'AODurh_flavs2'
    end if

    if (PrefixMatches('AODurh_s',prefix) .and. code == AODurh_s) then
       nocc = nocc + 1; name = 'AODurh_s'
    end if

    if (PrefixMatches('Pscheme',prefix) .and. code == Pscheme) then
       nocc = nocc + 1; name = 'Pscheme'
    end if

    if (PrefixMatches('P0scheme',prefix) .and. code == P0scheme) then
       nocc = nocc + 1; name = 'P0scheme'
    end if

    if (PrefixMatches('Escheme',prefix) .and. code == Escheme) then
       nocc = nocc + 1; name = 'Escheme'
    end if

    if (PrefixMatches('E0scheme',prefix) .and. code == E0scheme) then
       nocc = nocc + 1; name = 'E0scheme'
    end if

    if (PrefixMatches('FlavScheme_expb',prefix) .and. code == FlavScheme_expb) then
       nocc = nocc + 1; name = 'FlavScheme_expb'
    end if

    if (PrefixMatches('FlavScheme_net',prefix) .and. code == FlavScheme_net) then
       nocc = nocc + 1; name = 'FlavScheme_net'
    end if

    if (PrefixMatches('FlavScheme_mod2',prefix) .and. code == FlavScheme_mod2) then
       nocc = nocc + 1; name = 'FlavScheme_mod2'
    end if

    if (PrefixMatches('itype_ee',prefix) .and. code == itype_ee) then
       nocc = nocc + 1; name = 'itype_ee'
       if (present(longname)) longname = 'e+e- type kt-algo'
    end if

    if (PrefixMatches('itype_pp',prefix) .and. code == itype_pp) then
       nocc = nocc + 1; name = 'itype_pp'
       if (present(longname)) longname = 'pp type kt-algo'
    end if

    if (PrefixMatches('iangle_ang',prefix) .and. code == iangle_ang) then
       nocc = nocc + 1; name = 'iangle_ang'
       if (present(longname)) longname = 'min(E_i,E_j)^2 (1-cos theta_ij)'
    end if

    if (PrefixMatches('iangle_DeltaR',prefix) .and. code == iangle_DeltaR) then
       nocc = nocc + 1; name = 'iangle_DeltaR'
       if (present(longname)) longname = 'min(pt_i,pt_j)^2 (Delta R_ij)^2'
    end if

    if (PrefixMatches('iangle_DeltaR_flav1',prefix) .and. code == iangle_DeltaR_flav1) then
       nocc = nocc + 1; name = 'iangle_DeltaR_flav1'
    end if

    if (PrefixMatches('iangle_DeltaR_flav2',prefix) .and. code == iangle_DeltaR_flav2) then
       nocc = nocc + 1; name = 'iangle_DeltaR_flav2'
    end if

    if (PrefixMatches('iangle_DeltaR_flavs1',prefix) .and. code == iangle_DeltaR_flavs1) then
       nocc = nocc + 1; name = 'iangle_DeltaR_flavs1'
    end if

    if (PrefixMatches('iangle_DeltaR_flavs2',prefix) .and. code == iangle_DeltaR_flavs2) then
       nocc = nocc + 1; name = 'iangle_DeltaR_flavs2'
    end if

    if (PrefixMatches('iangle_DeltaRs',prefix) .and. code == iangle_DeltaRs) then
       nocc = nocc + 1; name = 'iangle_DeltaRs'
    end if

    if (PrefixMatches('irecom_P',prefix) .and. code == irecom_P) then
       nocc = nocc + 1; name = 'irecom_P'
    end if

    if (PrefixMatches('irecom_P0',prefix) .and. code == irecom_P0) then
       nocc = nocc + 1; name = 'irecom_P0'
    end if

    if (PrefixMatches('irecom_E',prefix) .and. code == irecom_E) then
       nocc = nocc + 1; name = 'irecom_E'
    end if

    if (PrefixMatches('irecom_E0',prefix) .and. code == irecom_E0) then
       nocc = nocc + 1; name = 'irecom_E0'
    end if

    if (PrefixMatches('irecom_Pt',prefix) .and. code == irecom_Pt) then
       nocc = nocc + 1; name = 'irecom_Pt'
    end if

    if (PrefixMatches('irecom_Pt2',prefix) .and. code == irecom_Pt2) then
       nocc = nocc + 1; name = 'irecom_Pt2'
    end if

    if (PrefixMatches('irecom_Et',prefix) .and. code == irecom_Et) then
       nocc = nocc + 1; name = 'irecom_Et'
    end if

    if (PrefixMatches('irecom_Et2',prefix) .and. code == irecom_Et2) then
       nocc = nocc + 1; name = 'irecom_Et2'
    end if

    if (PrefixMatches('ialg_pp_ktLI_E',prefix) .and. code == ialg_pp_ktLI_E) then
       nocc = nocc + 1; name = 'ialg_pp_ktLI_E'
    end if

    if (PrefixMatches('ialg_pp_ktf1_E',prefix) .and. code == ialg_pp_ktf1_E) then
       nocc = nocc + 1; name = 'ialg_pp_ktf1_E'
    end if

    if (PrefixMatches('ialg_pp_ktf2_E',prefix) .and. code == ialg_pp_ktf2_E) then
       nocc = nocc + 1; name = 'ialg_pp_ktf2_E'
    end if

    if (PrefixMatches('ialg_pp_ktfs1_E',prefix) .and. code == ialg_pp_ktfs1_E) then
       nocc = nocc + 1; name = 'ialg_pp_ktfs1_E'
    end if

    if (PrefixMatches('ialg_pp_ktfs2_E',prefix) .and. code == ialg_pp_ktfs2_E) then
       nocc = nocc + 1; name = 'ialg_pp_ktfs2_E'
    end if

    if (PrefixMatches('ialg_ee_kt_E',prefix) .and. code == ialg_ee_kt_E) then
       nocc = nocc + 1; name = 'ialg_ee_kt_E'
    end if

    if (PrefixMatches('ialg_ee_kt_P',prefix) .and. code == ialg_ee_kt_P) then
       nocc = nocc + 1; name = 'ialg_ee_kt_P'
    end if

    if (PrefixMatches('ialg_ee_kt_E0',prefix) .and. code == ialg_ee_kt_E0) then
       nocc = nocc + 1; name = 'ialg_ee_kt_E0'
    end if

    if (PrefixMatches('modelA',prefix) .and. code == modelA) then
       nocc = nocc + 1; name = 'modelA'
    end if

    if (PrefixMatches('modelB',prefix) .and. code == modelB) then
       nocc = nocc + 1; name = 'modelB'
    end if

    if (PrefixMatches('modelC',prefix) .and. code == modelC) then
       nocc = nocc + 1; name = 'modelC'
    end if

    if (PrefixMatches('modelD',prefix) .and. code == modelD) then
       nocc = nocc + 1; name = 'modelD'
    end if

    if (PrefixMatches('modelE',prefix) .and. code == modelE) then
       nocc = nocc + 1; name = 'modelE'
    end if

    if (PrefixMatches('modelF',prefix) .and. code == modelF) then
       nocc = nocc + 1; name = 'modelF'
    end if

    if (PrefixMatches('modelG',prefix) .and. code == modelG) then
       nocc = nocc + 1; name = 'modelG'
    end if

    if (PrefixMatches('modelH',prefix) .and. code == modelH) then
       nocc = nocc + 1; name = 'modelH'
    end if

    if (PrefixMatches('modelI',prefix) .and. code == modelI) then
       nocc = nocc + 1; name = 'modelI'
    end if

    if (PrefixMatches('modelJ',prefix) .and. code == modelJ) then
       nocc = nocc + 1; name = 'modelJ'
    end if

    ! common-code resumes
    if (nocc == 0) then
       if (present(prefix)) then
          write(0,*) 'Error in NameOfCode: could not find code',&
               &code,' with prefix "'//prefix//'"'
       else
          write(0,*) 'Error in NameOfCode: could not find code',&
               &code,' (without prefix)'
       end if
       stop
    end if

    if (nocc > 1) then
       if (present(prefix)) then
          write(0,*) 'Error in NameOfCode: several meanings for code',&
               &code,' with prefix "'//prefix//'"'
       else
          write(0,*) 'Error in NameOfCode: several meanings for code',&
               &code,' (without prefix)'
       end if
       stop
    end if

  end function NameOfCode



  !!
  !! for establishing whether a prefix matches a given string
  !!
  logical function PrefixMatches(string,prefix)
    character(len=*), intent(in)           :: string
    character(len=*), intent(in), optional :: prefix

    if (present(prefix)) then
       PrefixMatches = (index(string,prefix) == 1)
    else
       PrefixMatches = .true.
    end if
  end function PrefixMatches


end module NameSelect
