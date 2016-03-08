!! File generated automatically by autogen.pl from 
!! general precision template ../Common/genPREC/PRECktclus_gen.f90.
!! subsequently modified by hand...


! attempt to move towards hh version, with analysis of how to keep
! things quick, while not passing around vast amounts of information
! from one subroutine to the next.
!
! In practice it seems that the real "killer" in terms of speed
! (slowdown of 1.5 with intel compiler) is when one goes from having p
! passed along as a dummy argument, to when it is an element either 
! of a "event" derived type, or when it is a global variable; the only
! explanation I can find is that when it is passed as a dummy argument, 
! some memory prefetching has gone on which makes things more efficient.
! 
! Now updated to deal also with hh -- on vect-proc.dat it works perfectly, 
! while with the herwig runs it does seem to have some difference relative
! to the official version (though not in actual y23 results) -- this
! may be because of rounding errors (in ktclusalt?)?
module dpktclus_gen
  use types; use consts_dp
  use warnings_and_errors; use assertions
  use alg_declarations 
  use dpsorter
  implicit none
  private
  
  type iy
     integer  :: i
     real(dp) :: y
  end type 
  
  real(dp), public , parameter :: ten_pow200 = 1e200_dp
  real(dp) :: Evis, Evis2
  ! -- for normal e+e- case (iangle_ang)
  integer, parameter  :: i_modp3 = 5
  ! -- for hh case (iangle_DeltaR)
  integer, parameter  :: i_ptsq  = 5
  integer, parameter  :: i_eta   = 6
  integer, parameter  :: i_phi   = 7
  integer, parameter  :: maxdim1p = 7
  integer :: nyres

  !------------ internal variables ------------
  type(iy), allocatable :: iy_min_atj(:,:)
  integer               :: sub_size
  integer               :: iangle_glb
  ! the distance between particles and beam jets
  real(dp), allocatable :: ybeam(:,:)
  real(dp), allocatable :: ybeam_calc(:,:)
  integer, allocatable  :: indx(:)
  integer, allocatable  :: flavinfo(:,:)
  integer, allocatable  :: flavbeam(:,:)
  real(dp), parameter :: y_impossible = 1e100_dp   
  real(dp), parameter :: ptsq_infinitesimal = 1e-200_dp
  logical               :: use_ptbeam_glb
  
  type node
     integer    :: prnt1, prnt2
     real(dp) :: y
  end type node
  type nodelist
     integer  :: itype, iangle, irecom, nmin, n
     logical  :: etacut, monotonic, use_local_ptbeam
     real(dp) :: ptbeam_rescale, R, R2, invR2
     ! momenta run from 1:2n, with n+1:2n being the original n 
     ! particles and 1:n being the momenta at the recombination 
     ! nodes; expires has same numbering
     real(dp),   pointer :: p(:,:) => null()
     integer,    pointer :: expires(:) => null()
     ! nodes run from nmin:n
     type(node), pointer :: nodes(:) => null()
  end type nodelist

  type(nodelist), pointer :: ndlist_ptr

  public :: node, nodelist
  
  public :: ktclus_genR
  public :: ktclus_gen, ktclus_gen_nonodes, delete_nodelist, get_m_jets
  public :: init_nodelist, recombine_nodes, get_ym, get_m, get_all_jets
  
  interface get_phi
     module procedure get_phi_0d, get_phi_1d
  end interface
  interface get_pt2
     module procedure get_pt2_0d, get_pt2_1d
  end interface
  interface get_Et2
     module procedure get_Et2_0d, get_Et2_1d
  end interface
  interface get_pt
     module procedure get_pt_0d, get_pt_1d
  end interface
  interface get_Et
     module procedure get_Et_0d, get_Et_1d
  end interface
  interface get_rap
     module procedure get_rap_0d, get_rap_1d
  end interface
  public :: get_phi, get_pt2, get_Et2, get_pt, get_Et, get_rap


contains

  !======================================================================
  !! This is a legacy interface to the jet clustering. See below
  !! (ktclus_genR) for an explanation of all arguments
  !!
  subroutine ktclus_gen(pp,ndlist,imode,flavour,beamflavour,&
       &ptbeam_rescale, use_ptbeam, use_local_ptbeam, etacut)
    use alt_flavclus
    !use sort
    real(dp),     intent(in)  :: pp(:,:)
    type(nodelist), intent(inout), target :: ndlist
    integer,        intent(in)  :: imode
    integer,  optional, intent(in)  :: flavour(:,:)
    integer,  optional, intent(in)  :: beamflavour(:,:)
    real(dp), optional, intent(in) :: ptbeam_rescale
    logical,  optional, intent(in)  :: use_ptbeam, use_local_ptbeam
    logical,  optional, intent(in)  :: etacut
    call ktclus_genR(pp, ndlist, imode, R=one,&
        &            flavour          = flavour,&
        &            beamflavour      = beamflavour,&
        &            ptbeam_rescale   = ptbeam_rescale, &
        &            use_ptbeam       = use_ptbeam, &
        &            use_local_ptbeam = use_local_ptbeam, &
        &            etacut           = etacut)
  end subroutine ktclus_gen
  

  !======================================================================
  !! 
  !! Routine that performs KT clustering and flavour variants. 
  !!
  !! The momentta are to be supplied in the array pp and the output is
  !! the ndlist structure which allows to one to efficiently extract
  !! information about the clustering structure.
  !!
  !! The parameter R is the usual inclusive kt-algorithm R-parameter,
  !! i.e. dij_R = dij_(R=1) / R^2 and diB is independent of R.
  !!
  !! You MUST call delete_nodelist(ndlist) when you have finished
  !! using it, so as to free up dynamically allocated memory. [Some
  !! attempt will be made at automatic cleanup of an already allocated
  !! ndlist -- but this will only work if calling subroutine gives the
  !! "save" attribute to ndlist.]
  !!
  !! For flavour algorithms the flavour array must be supplied. It has
  !! dimensions (1:nflav,1:n) where nflav is the number of active
  !! flavours being considered and n is the number of particles. For a
  !! particle i, flavour(:,i) indicates the net amount of each kind of
  !! flavour in the particle. E.g.
  !!    flavour(:,i) = (/-1,1,0,0,0,0/)
  !! would be a dbar-u meson-like particle.
  !! 
  !! For DIS (not implemented) and hh algorithms beamflavour is the
  !! flavour of the beams, and is an array of size (nflav,-nbeam:-1). For
  !! hh, we assume that the beams are along the z-axis and that beam -2
  !! has negative pz, while beam -1 has positive pz.
  !!
  !! ptbeam_rescale and use_ptbeam modify behaviour of hh flavour
  !! algorithms (these algorithms use a special "ptbeam" rather than
  !! pt^2(i) for the d_iB of certain beam recombinations:
  !! ptbeam_rescale renormalises ptbeam and use_ptbeam=.false. turns
  !! off this part of the algorithm).
  !!
  !! etacut=.false. kills all recombination with beams (it is designed
  !! for hh case where only central particles are supplied)
  !!
  !---------------------------------------------------------------------
  !!
  !! Notes on algorithms:
  !! - the underlying kt-algorithm should be at worst n^{5/2} and
  !!   hopefully in practice will be n^2
  !!                                            k   ,    j
  !!   Idea is to store info about yij in an sqrt(N) by N array, where
  !!   in the short direction (k) we store the minimum of the yij, with
  !!   i running between jbase+1:jbase+sub_size, where jbase = 
  !!   (k-1)*sub_size; idea is, for a recombination ij, to reduce to
  !!   O(sqrt(N)) operations each of the updates of a column that is
  !!   neither i nor j (which tend to dominate the run-time).
  !!
  !! - flavour variants of the algorithms are probably a lot slower:
  !!   still N^2 for e+e- (not yet implemented!) but N^2 ln N for hh.
  !!   The ln N factor can be eliminated, but this has not yet been done.
  !!
  !! - note that methods for N ln N implementation of the plain kt
  !!   algo can probably not be immediately applied to flavour 
  !!   variants.
  !!
  subroutine ktclus_genR(pp,ndlist,imode,R,flavour,beamflavour,&
       &ptbeam_rescale, use_ptbeam, use_local_ptbeam, etacut)
    use alt_flavclus
    !use sort
    real(dp),     intent(in)  :: pp(:,:)
    type(nodelist), intent(inout), target :: ndlist
    integer,        intent(in)  :: imode
    real(dp), optional, intent(in)  :: R
    integer,  optional, intent(in)  :: flavour(:,:)
    integer,  optional, intent(in)  :: beamflavour(:,:)
    real(dp), optional, intent(in)  :: ptbeam_rescale
    logical,  optional, intent(in)  :: use_ptbeam, use_local_ptbeam
    logical,  optional, intent(in)  :: etacut
    !---------------------------------
    real(dp) :: p(maxdim1p, size(pp,dim=2)), y
    real(dp) :: ranarray(size(pp,dim=2))
!    integer  :: indx(size(pp,dim=2))
    integer  :: inodes(size(pp,dim=2))
    logical  :: reorder
    !type(iy), allocatable :: iy_min_atj(2:size(pp,dim=2))
    
    integer  :: n, i, j, ii, nsub, jbase, k, kmx
    integer    :: jj     ! FOR IFORT WORKAROUND
    real(dp) :: minvalue ! FOR IFORT WORKAROUND
    integer  :: ibeam    
    logical  :: hh
    !integer, parameter :: imode_debug = 4203
    !integer, parameter :: imode_debug = 4303
    !integer, parameter :: imode_debug = 4603
    !integer, parameter :: imode_debug = 4503
    !integer, parameter :: imode_debug = 4403
    !integer, parameter :: imode_debug = 4703
    integer, parameter :: imode_debug = -1
    logical :: useflav

    call init_nodelist(ndlist,pp,imode,R,etacut,ptbeam_rescale,use_local_ptbeam)
    ndlist_ptr => ndlist
    ! globals to be set here
    iangle_glb = ndlist%iangle
    use_ptbeam_glb = default_or_opt(.true., use_ptbeam)

    ! consistency checks on current status of implementation
    if ((ndlist%etacut .and. ndlist%itype /= itype_pp) .or. &
         &(ndlist%itype == itype_pp .and. ndlist%iangle > iangle_DeltaRs)) then
       call wae_error('ktclus_gen','Unimplemented imode:',intval=imode)
    end if

    

    ! tracking
    nyres = 0

    n        = ndlist%n
    sub_size = floor(sqrt(one*n))
    !sub_size = floor((one*n)**0.6_dp) ! tried playing with 0.3-0.7
    !                                  ! but 0.5 seems best...
    nsub     = kofn(n-1)
    allocate (iy_min_atj(0:nsub,2:n))
    ! initialise info for use of nodelist
    do i = 1, n
       inodes(i) = n+i ! node associated with pseudoparticle i in local list
    end do

    
    
    !-- reorder the p to see if it helps: gives a marginal (3 percent)
    !   improvement, so not clear whether it's worth the effort...
    ! reorder = .false.
    ! if (reorder) then
    !    write(0,*) 'reordering'
    !    call random_number(ranarray)
    !    call indexx(ranarray,indx)
    !    p(:4,:) = pp(:4,indx)
    ! else
    !    p(:4,:) = pp(:4,:)
    ! end if

    p(:4,:) = pp(:4,:)
    do i = 1, size(pp,2)
       call setup_p_info(p(:,i))
    end do

    hh = .not.ndlist%etacut .and. ndlist%itype == itype_pp

    ! behave differently for flavour and non-flavour algos

    if (present(flavour)) then 
       allocate (flavinfo(size(flavour,dim=1),1:n))
       flavinfo = flavour
       ! NB: think what to do about DIS?
       if (hh) then 
          allocate(flavbeam(size(flavour,dim=1),-size(beamflavour,dim=2):-1))
          flavbeam = beamflavour
       end if
    end if

    if (ndlist%iangle == iangle_ang .and. ndlist%R /= one) then
       call wae_error("ktclus_genR",&
            & "for iangle=iangle_ang, R must be 1, instead it was",&
            & dbleval = R)
    end if
    

    select case(ndlist%iangle)
    case(iangle_DeltaR_flav1 ,&
       & iangle_DeltaR_flav2 ,&
       & iangle_DeltaR_flavs1,&
       & iangle_DeltaR_flavs2,&
       & iangle_DeltaRs      )

       useflav = .true.

       if (.not. present(flavour)) call wae_error('ktclus_gen',&
            &'Requested flavour algorithm, but no flavour info provided')
       
       !allocate (flavinfo(size(flavour,dim=1),1:n))
       !flavinfo = flavour

       if (hh) then
          if (.not. present(beamflavour)) call wae_error('ktclus_gen',&
          &'Requested hh flavour algorithm, but no beam flavour info provided')

          if (size(flavour,dim=1) /= size(beamflavour,dim=1)) call wae_error(&
               & 'ktclus_gen','size of dim-1 of flavour and beamflavour&
               & differ')

          !! NB: think what to do about DIS?
          !allocate(flavbeam(size(flavour,dim=1),-size(beamflavour,dim=2):-1))
          !flavbeam = beamflavour

       end if
       
    case default
       ! DO NOT REQUIRE ABSENCE OF FLAVOUR INFO: there are natural
       ! cases (considering many algos) when it could be supplied 
       ! but not used.
       useflav = .false.
    end select

    !---- for debugging -----------------
    if (imode == imode_debug .and. present(flavour)) then
       call do_alt_flavclus(pp,flavinfo,flavbeam,imode)
    end if


    ! set up the basic array...
    do j = 2, size(pp,2)
       call set_iy_min_j(p,j)
    end do

    ! set up the distances between particles and the beam jets
    ! -2 proton (infinite negative rapidity)
    ! -1 antiproton (infinite positive rapidity)
    ! in ybeam(:,-2:-1) particles are ordered in rapidity
    ! indx(i) is the id of the i-th particle in ybeam(1:n,-2:-1)
    if (hh) then
       allocate(indx(1:n),ybeam(1:n,-2:-1))
       call set_ybeam(p)
       !! for testing and display purposes
       !call writeout_ptbeam_curves(p)
    end if


    do n = size(pp,2),ndlist%nmin,-1

       !-- find minimum of the yij values
       if (n >= 2) then
          ! add 1 because iy_min_atj starts from 2!
          ! GIVES COMPILER ERROR with ifort 8.1.21
          j = 1 + sum(minloc(iy_min_atj(0,2:n)%y))
          ! WORKAROUND for intel 8.1.21 internal compiler error
          !minvalue = iy_min_atj(0,2)%y; j = 2
          !do jj = 3, n
          !   if (iy_min_atj(0,jj)%y < minvalue) then
          !      j = jj; minvalue = iy_min_atj(0,jj)%y
          !   end if
          !end do
          ! for testing ...
          !write(0,*) j, 1 + sum(minloc(iy_min_atj(0,2:n)%y))
          y = iy_min_atj(0,j)%y 
          i = iy_min_atj(0,j)%i
       end if

       if (hh) then
          ! new version
          ii = sum(minloc(ybeam(1:n,-1))); ibeam = -1
          jj = sum(minloc(ybeam(1:n,-2)))
          if (ybeam(jj,-2) <= ybeam(ii,-1)) then
             ii = jj; ibeam = -2
          end if
          if (ybeam(ii,ibeam) <= y .or. n==1) then
             y = ybeam(ii,ibeam)
             j = indx(ii)
             i = ibeam
          end if

       end if

       call recombine(p,i,j,n,ndlist%irecom)
       call recombine_nodes(p,i,j,n,y,ndlist,inodes)

       if (useflav) then
          call recombine_flavour(i,j,n,hh)
       end if

       if (n > 1) then
          call renew_iy_min_atj(p,i,j,n)
          if (.not.ndlist%etacut .and. ndlist%itype == itype_pp) then
             call renew_ybeam(p,i,j,n)
          end if
       end if
       if (imode == imode_debug) then
          write(21,'(3i9,es20.7)') n,i,j,y
       end if
    end do

    select case(ndlist%itype)
    case(itype_ee)
       Evis2 = sum(pp(4,:))**2
       ndlist%nodes%y = two*ndlist%nodes%y / Evis2
    case(itype_pp)
       ! do nothing
    case default
       call wae_error('ktclus_gen: Unsupported itype', intval=ndlist%itype)
    end select


    deallocate(iy_min_atj)
    deallocate(indx,ybeam)

!    if (useflav) then
    if (present(flavour)) then
       deallocate(flavinfo)
       if (hh) deallocate(flavbeam)
    end if

  end subroutine ktclus_genR
  



  ! attempt at more efficient ktclus algorithm -- this time
  ! using an algorithm that should be at worst n^{5/2} and hopefully
  ! in practice will be n^2 with a somewhat smaller coefficient than
  ! v1.
  !                                         k   ,    j
  ! Idea is to store info about yij in an sqrt(N) by N array, where
  ! in the short direction (k) we store the minimum of the yij, with
  ! i running between jbase+1:jbase+sub_size, where jbase = 
  ! (k-1)*sub_size; idea is, for a recombination ij, to reduce to
  ! O(sqrt(N)) operations each of the updates of a column that is
  ! neither i nor j (which tend to dominate the run-time).
  subroutine ktclus_gen_nonodes(pp,yvals,imode,etacut)
    !use sort
    real(dp), intent(in)  :: pp(:,:)
    real(dp), intent(out) :: yvals(:)
    integer, intent(in)   :: imode
    logical,  intent(in)  :: etacut
    !---------------------------------
    real(dp) :: p(maxdim1p, size(pp,dim=2))
    real(dp) :: ranarray(size(pp,dim=2))
    integer    :: irecom, itype,nmin
!    integer  :: indx(size(pp,dim=2))
    logical  :: reorder
    !type(iy), allocatable :: iy_min_atj(2:size(pp,dim=2))
    
    integer  :: n, i, j, ii, nsub, jbase, k, kmx
    integer    :: jj     ! FOR IFORT WORKAROUND
    real(dp) :: minval ! FOR IFORT WORKAROUND
    
    !-- first sort out the kind of algo
    itype  = (imode/1000)*1000
    iangle_glb = (mod(imode,1000)/100)*100
    irecom = mod(imode,10)

    ! consistency checks on current status of implementation
    if ((etacut .and. itype /= itype_pp) .or. &
         &(itype == itype_pp .and. iangle_glb /= iangle_DeltaR)) then
       call wae_error('ktclus_gen','Unimplemented imode:',intval=imode)
    end if
    

    ! direct consequences...
    if (etacut .or. itype /= itype_pp) then
       nmin = 2
    else
       nmin = 1
    end if

    ! tracking
    nyres = 0

    n        = size(pp,dim=2)
    sub_size = floor(sqrt(one*n))
    !sub_size = floor((one*n)**0.6_dp) ! tried playing with 0.3-0.7
    !                                  ! but 0.5 seems best...
    nsub     = kofn(n-1)
    allocate (iy_min_atj(0:nsub,2:n))
    
    
    
    !-- reorder the p to see if it helps: gives a marginal (3 percent)
    !   improvement, so not clear whether it's worth the effort...
    ! reorder = .false.
    ! if (reorder) then
    !    write(0,*) 'reordering'
    !    call random_number(ranarray)
    !    call indexx(ranarray,indx)
    !    p(:4,:) = pp(:4,indx)
    ! else
    !    p(:4,:) = pp(:4,:)
    ! end if

    p(:4,:) = pp(:4,:)
    do i = 1, size(pp,2)
       call setup_p_info(p(:,i))
    end do

    Evis2 = sum(sqrt(p(4,:)))**2

    ! set up the basic array...
    do j = 2, size(pp,2)
       call set_iy_min_j(p,j)
    end do
    
    do n = size(pp,2),nmin,-1
       !-- find minimum of the yij values
       if (n >= 2) then
          ! add 1 because iy_min_atj starts from 2!
          ! GIVES COMPILER ERROR with ifort 8.1.21
          !j = 1 + sum(minloc(iy_min_atj(0,2:n)%y))
          ! WORKAROUND for intel 8.1.21 internal compiler error
          minval = iy_min_atj(0,2)%y; j = 2
          do jj = 3, n
             if (iy_min_atj(0,jj)%y < minval) then
                j = jj; minval = iy_min_atj(0,jj)%y
             end if
          end do
          ! for testing ...
          !write(0,*) j, 1 + sum(minloc(iy_min_atj(0,2:n)%y))
          yvals(n) = iy_min_atj(0,j)%y 
          i = iy_min_atj(0,j)%i
       end if
       
       ! deal in a not very general manner with the minimum of the pt**2
       if (.not.etacut .and. itype == itype_pp) then
          ii = sum(minloc(p(i_ptsq,:n)))
          if (p(i_ptsq,ii) <= yvals(n) .or. n==1) then
             yvals(n) = p(i_ptsq,ii)
             j = ii
             i = -1
          end if

       end if
       
       call recombine(p,i,j,n,irecom)
       if (n > 1) then
          call renew_iy_min_atj(p,i,j,n)
       end if
    end do

    select case(itype)
    case(itype_ee)
       yvals = two*yvals / Evis2
    case(itype_pp)
       ! do nothing
    case default
       call wae_error('ktclus_gen: Unsupported itype', intval=itype)
    end select

    deallocate (iy_min_atj)
  end subroutine ktclus_gen_nonodes



  !----------------------------------------------
  ! recombines particles i and j into position i
  ! and moves n to position j
  subroutine recombine(p,i,j,n,irecom)
    real(dp), intent(inout) :: p(:,:)
    integer,  intent(in)    :: i, j, n
    integer,  intent(in)    :: irecom    
    !--------------------------------------------
    real(dp) :: ptsum, pti,ptj, pmod3
    real(dp) :: dphi
    real(dp) :: pcopy(size(p,1), size(p,2))

    if (i > 0) then

       select case(irecom)
       case(irecom_E)
          p(:4,i) = p(:4,i) + p(:4,j)
          call setup_p_info(p(:,i))
       case(irecom_P)
          p(:3,i) = p(:3,i) + p(:3,j)
          p(4,i)  = sqrt(p(1,i)**2+p(2,i)**2+p(3,i)**2)
          call setup_p_info(p(:,i))
       case(irecom_E0)
          p(:4,i) = p(:4,i) + p(:4,j)
          pmod3   = sqrt(p(1,i)**2+p(2,i)**2+p(3,i)**2)
          p(:3,i) = p(:3,i)*p(4,i)/pmod3
          call setup_p_info(p(:,i))
          ! in these two schemes we do not update the four vectors 
       case(irecom_Pt)
          pcopy(:4,i) = p(:4,i) + p(:4,j)
          call setup_p_info(pcopy(:,i))

          pti = sqrt(p(i_ptsq,i))
          ptj = sqrt(p(i_ptsq,j))
          ptsum = pti+ptj
          p(i_ptsq,i) = ptsum**2
!          write(6,*) '------------------------------------------'
!          write(6,*) 'recombining', i,j, dble(pti),dble(ptj)
! 
!
!          write(6,*) 'eta', dble(p(i_eta,i)), dble(p(i_eta,j)),&
!               &dble(p(i_eta,i)-p(i_eta,j))
!          write(6,*) 'dcut', log(dble(ptj**2*(dphi**2+(p(i_eta,i)-p(i_eta,j))**2)/(two*pti)))

          p(i_eta,i)  = (pti*p(i_eta,i)+ptj*p(i_eta,j))/ptsum

          if (abs(p(i_phi,i)-p(i_phi,j)) > pi) then
             p(i_phi,j) =  p(i_phi,j) - twopi*sign(one,p(i_phi,j))
          end if

!         write(6,*) 'phi', real(p(i_phi,i)), dble(p(i_phi,j)),&
!               &dble(abs(p(i_phi,i)-p(i_phi,j)))

          p(i_phi,i)  = (pti*p(i_phi,i)+ptj*p(i_phi,j))/ptsum

          if (abs(p(i_phi,i)) > pi) then
             p(i_phi,i)=p(i_phi,i)-twopi*sign(one,p(i_phi,i))
          end if

!          write(6,*) 'phi_new, eta_new', dble(ptj/ptsum*twopi), dble(p(i_phi,i)+pi),&
!               & dble(p(i_eta,i))

       case(irecom_Pt2)
          ptsum = p(i_ptsq,i)+p(i_ptsq,j)

          p(i_eta,i) = (p(i_ptsq,i)*p(i_eta,i)+p(i_ptsq,j)*p(i_eta,j))/ptsum

          if (abs(p(i_phi,i)-p(i_phi,j)) > pi) then
             p(i_phi,j) =  p(i_phi,j) - twopi*sign(one,p(i_phi,j))
          end if

          p(i_phi,i) = (p(i_ptsq,i)*p(i_phi,i)+p(i_ptsq,j)*p(i_phi,j))/ptsum

          if (abs(p(i_phi,i)) > pi) then
             p(i_phi,i)=p(i_phi,i)-twopi*sign(one,p(i_phi,i))
          end if

          p(i_ptsq,i) = (sqrt(p(i_ptsq,i))+sqrt(p(i_ptsq,j)))**2
       case default
          call wae_error("recombine", "Unrecognised recombination scheme")
       end select
    end if

    if (j < n) then
       p(:,j)  = p(:,n)
    end if
    
  end subroutine recombine

  !---------------------------------------------------
  ! carry out all initialisation and a number of consistency
  ! checks and extraction of information (e.g. imode -> itype,...)
  ! for the nodelist
  subroutine init_nodelist(ndlist,pp,imode,R,etacut,ptbeam_rescale,use_local_ptbeam)
    type(nodelist), intent(inout)  :: ndlist
    real(dp),     intent(in)       :: pp(:,:)
    integer,        intent(in)     :: imode
    real(dp), optional, intent(in) :: R
    logical,  optional, intent(in) :: etacut,use_local_ptbeam
    real(dp), optional, intent(in) :: ptbeam_rescale
    !-----------------------------------------

    ! avoid memory leaks...
    if (associated(ndlist%nodes)) call delete_nodelist(ndlist)

    !-- first sort out the kind of algo
    ndlist%itype  = (imode/1000)*1000
    ndlist%iangle = (mod(imode,1000)/100)*100
    ndlist%irecom = mod(imode,10)
    ndlist%R      = default_or_opt(one,R)
    ndlist%R2     = R**2
    ndlist%invR2  = one/ndlist%R2
    ndlist%etacut = default_or_opt(.false.,etacut)
    ndlist%monotonic = .true. ! if it is not this will be realised
                              ! during run
    ndlist%ptbeam_rescale = default_or_opt(one,ptbeam_rescale)
    ndlist%use_local_ptbeam = default_or_opt(.false.,use_local_ptbeam)

    ! direct consequences...
    if (ndlist%etacut .or. ndlist%itype /= itype_pp) then
       ndlist%nmin = 2
    else
       ndlist%nmin = 1
    end if

    ndlist%n = size(pp,dim=2)
    allocate(ndlist%p(4,2*ndlist%n), ndlist%expires(2*ndlist%n), &
         &   ndlist%nodes(ndlist%n))
    ndlist%p(:4,ndlist%n+1:2*ndlist%n) = pp(:4,1:ndlist%n)
    ndlist%expires = 0
    
    ! AB initialise the first node for e+e- annihilation
    if (ndlist%itype == itype_ee) then
       ndlist%nodes(1:ndlist%n)%prnt1 = 0
       ndlist%nodes(1:ndlist%n)%prnt2 = 0
    end if
       
  end subroutine init_nodelist
  
  !-----------------------------------------------------------------
  ! Carry out the part of the recombination procedure associated 
  ! with the node structure
  subroutine recombine_nodes(p,i,j,n,y,ndlist,inodes)
    real(dp),     intent(in)    :: p(:,:)
    integer,        intent(in)    :: i, j, n
    real(dp),     intent(in)    :: y
    type(nodelist), intent(inout) :: ndlist
    integer,        intent(inout) :: inodes(:)
    !-----------------------------------------

    !ndlist%nodes(n)%y = y
    ndlist%nodes(n)%y = y * ndlist%invR2 ! GPS-R: include R scale

    if (j <= 0) call wae_error('recombine_nodes',&
         &'j had illegal value:',intval=j)
    if (i < 0) then
       ndlist%nodes(n)%prnt1 = i
       ndlist%p(:,:n)        = zero
       !implement the fact that when a (pseudo)particle is recombined with 
       ! the beam the corresponding jet is eliminated from the jet list 
       ndlist%expires(n)     = n
    else
       ndlist%nodes(n)%prnt1     = inodes(i)
       ndlist%expires(inodes(i)) = n
       ndlist%p(:4,n)            = p(:4,i)
       inodes(i)                 = n
    end if
    ndlist%nodes(n)%prnt2 = inodes(j)
    ndlist%expires(inodes(j)) = n
    
    if (j < n) inodes(j) = inodes(n)
    if (n < ndlist%n) then
       ndlist%monotonic = ndlist%monotonic .and. &
            & (ndlist%nodes(n)%y < ndlist%nodes(n+1)%y)
    end if
  end subroutine recombine_nodes

  !----------------------------------------------
  ! Deallocate things so as to avoid memory leaks
  subroutine delete_nodelist(ndlist)
    type(nodelist), intent(inout) :: ndlist
    deallocate(ndlist%p,ndlist%nodes,ndlist%expires)
  end subroutine delete_nodelist
  
  !-------------------------------------------------------
  ! Returns momenta of the state when clustered to m jets
  ! Should function in a time O(m)
  subroutine get_m_jets(ndlist, m, p)
    type(nodelist), intent(in)  :: ndlist
    integer,        intent(in)  :: m
    real(dp),     intent(out) :: p(:,:)
    !--------------------------------
    integer :: i, njets
    if (size(p,1) < 4 .or. size(p,2) < m) call wae_error('get_m_jets',&
         &'array for output jet momenta is too small')

    if (m > ndlist%n .or. m<ndlist%nmin) call wae_error('get_m_jets', &
         &'m > n  or m<nmin (this is not allowed); m was:',intval=m)

    njets = 0
    do i = ndlist%nmin, m
       if (ndlist%nodes(i)%prnt1 > m) then
          njets = njets + 1
          p(:4,njets) = ndlist%p(:4,ndlist%nodes(i)%prnt1)
       end if
       if (ndlist%nodes(i)%prnt2 > m) then
          njets = njets + 1
          p(:4,njets) = ndlist%p(:4,ndlist%nodes(i)%prnt2)
       end if
    end do

    if (njets /= m) call wae_error('get_m_jets',&
         &'internal inconsistency, njets /= m; njets was:', intval=njets)
  end subroutine get_m_jets
  
  !-----------------------------------------------------------
  ! Returns momenta of all jets in the "inclusive" sense of a
  ! longitudinally-invariant kt algorithm, with m indicating the
  ! number of jets.
  !
  ! If one of sortPt (sortEt) is present and .true. the jets are
  ! sorted into decreasing order in transverse Pt (Et). Otherwise
  ! they are returned in the same order as which they appear in
  ! nodes (inverse order of jet clustering).
  !
  subroutine get_all_jets(ndlist, p, m, sortPt, sortEt)
    use dpsorter
    type(nodelist), intent(in)  :: ndlist
    real(dp),     intent(out) :: p(:,:)
    integer,        intent(out) :: m
    logical, optional, intent(in):: sortPt
    logical, optional, intent(in):: sortEt
    !------------------------------------
    integer :: i, jets(1:ndlist%n)
    integer,    allocatable :: idx(:)

    if (ndlist%iangle < iangle_DeltaR) call wae_error('get_all_jets',&
         &'This routine should only be called if iangle>=iangle_DeltaR. &
         &Instead it was', intval=ndlist%iangle)
    
    m = 0
    do i = ndlist%nmin, ndlist%n
       if (ndlist%nodes(i)%prnt1 < 0) then
          m = m + 1
          jets(m) = ndlist%nodes(i)%prnt2
       end if
    end do

    if (size(p,dim=2) < m) call wae_error('get_all_jets',&
         &'dim2 of p is smaller than m, which is',intval=m)
    
    if (default_or_opt(.false.,sortPt)) then
       allocate(idx(1:m))
       call indexx(get_pt2(ndlist%p(:4,jets(1:m))),idx)
       p(:4,1:m) = ndlist%p(:4,jets(idx(m:1:-1)))
       deallocate(idx)
    else if (default_or_opt(.false.,sortEt)) then
       allocate(idx(1:m))
       call indexx(get_Et2(ndlist%p(:4,jets(1:m))),idx)
       p(:4,1:m) = ndlist%p(:4,jets(idx(m:1:-1)))
       deallocate(idx)
    else
       p(:4,1:m) = ndlist%p(:4,jets(1:m))
    end if
    
  end subroutine get_all_jets
  
  
  !-----------------------------------------------------------------
  ! Returns maximum y such that event has m jets (or potentially 
  ! more if non-monotonic) -- i.e. what is usually called y_{m-1,m}
  function get_ym(ndlist,m) result(ym)
    type(nodelist), intent(in) :: ndlist
    integer,        intent(in) :: m
    real(dp)                 :: ym

    if (m > ndlist%n) then
       ym = zero
       return
    end if
    if (m < ndlist%nmin) call wae_error('get_ym', &
         &'m<nmin (this is not allowed); m was:',intval=m)

    if (ndlist%monotonic) then
       ym = ndlist%nodes(m)%y
    else
       ym = maxval(ndlist%nodes(m:)%y)
    end if
  end function get_ym
  
  
  !-----------------------------------------------------------------
  ! Returns number of jets (m) present for a given value of y
  function get_m(ndlist, y) result(m)
    type(nodelist), intent(in) :: ndlist
    real(dp),     intent(in) :: y
    integer                    :: m
    !-------------------------------
    ! algo should be same even for non-monotonic situation...
    do m = ndlist%n, ndlist%nmin,-1
       if (ndlist%nodes(m)%y > y) exit
    end do
  end function get_m
  
  !----------------------------------------------------
  ! some overloaded "service" functions: not necessarily always
  ! the fastest approach if "extra" information is available
  function get_phi_0d(p) result(phi)
    real(dp), intent(in) :: p(:)
    real(dp)             :: phi
    phi = atan2(p(2),p(1))
  end function get_phi_0d
  function get_phi_1d(p) result(phi)
    real(dp), intent(in) :: p(:,:)
    real(dp)             :: phi(size(p,dim=2))
    phi(:) = atan2(p(2,:),p(1,:))
  end function get_phi_1d

  function get_pt2_0d(p) result(pt2)
    real(dp), intent(in) :: p(:)
    real(dp)             :: pt2
    pt2 = (p(1)**2 + p(2)**2)
  end function get_pt2_0d
  !
  function get_pt2_1d(p) result(pt2)
    real(dp), intent(in) :: p(:,:)
    real(dp)             :: pt2(size(p,dim=2))
    pt2(:) = (p(1,:)**2 + p(2,:)**2)
  end function get_pt2_1d
  !
  function get_pt_0d(p) result(pt)
    real(dp), intent(in) :: p(:)
    real(dp)             :: pt
    pt = sqrt(p(1)**2 + p(2)**2)
  end function get_pt_0d
  !
  function get_pt_1d(p) result(pt)
    real(dp), intent(in) :: p(:,:)
    real(dp)             :: pt(size(p,dim=2))
    pt(:) = sqrt(p(1,:)**2 + p(2,:)**2)
  end function get_pt_1d
  !
  function get_Et2_0d(p) result(Et2)
    real(dp), intent(in) :: p(:)
    real(dp)             :: Et2
    Et2 = p(4)**2 * ((p(1)**2 + p(2)**2) / (p(1)**2+p(2)**2+p(3)**2))
  end function get_Et2_0d
  !
  function get_Et2_1d(p) result(Et2)
    real(dp), intent(in) :: p(:,:)
    real(dp)             :: Et2(size(p,dim=2))
    Et2(:) = (p(1,:)**2 + p(2,:)**2) &
         &* p(4,:)**2/(p(1,:)**2+p(2,:)**2+p(3,:)**2)
  end function get_Et2_1d
  !
  function get_Et_0d(p) result(Et)
    real(dp), intent(in) :: p(:)
    real(dp)             :: Et
    Et = p(4) * sqrt((p(1)**2 + p(2)**2) / (p(1)**2+p(2)**2+p(3)**2))
  end function get_Et_0d
  !
  function get_Et_1d(p) result(Et)
    real(dp), intent(in) :: p(:,:)
    real(dp)             :: Et(size(p,dim=2))
    Et(:) = p(4,:)*sqrt((p(1,:)**2 + p(2,:)**2) &
         &              / (p(1,:)**2+p(2,:)**2+p(3,:)**2))
  end function get_Et_1d
  !
  function get_rap_0d(p) result(eta)
    real(dp), intent(in) :: p(:)
    real(dp)             :: eta
    if (abs(p(3)) /= p(4)) then
       eta = half*log((p(4)+p(3))/(p(4)-p(3)))
    else
       eta = sign(ten_pow200,p(3))
    end if
  end function get_rap_0d
  !
  function get_rap_1d(p) result(eta)
    real(dp), intent(in) :: p(:,:)
    real(dp)             :: eta(size(p,dim=2))
    where(.not.(abs(p(3,:)) == p(4,:)))
       eta(:) = half*log((p(4,:)+p(3,:))/(p(4,:)-p(3,:)))
    elsewhere
       ! had problems 
       eta(:) = sign(ten_pow200,p(3,:))
    end where
  end function get_rap_1d
  
  !-----------------------------------------------------------------
  ! does the bookkeeping and new calculations needed
  ! to keep iy_min_atj uptodate after the recombination of i,j 
  ! and the move of n to position j.
  subroutine renew_iy_min_atj(p,i,j,n)
    real(dp), intent(in)    :: p(:,:)
    integer,  intent(in)    :: i, j, n
    !----------------------------------
    integer :: jj, k, kmx

    if (j > 1) then
       ! as much as possible just copy info from the n column to j
       kmx = kofn(j-1)
       iy_min_atj(1:kmx-1,j) = iy_min_atj(1:kmx-1,n)
       if (iy_min_atj(kmx,n)%i < j) then
          iy_min_atj(kmx,j) = iy_min_atj(kmx,n)
       else
          iy_min_atj(kmx,j) = iy_min_kj(p,kmx,j)
       end if
       ! reestablish minimum (here too could look if minimum was already well-
       ! positioned?)
       k = sum(minloc(iy_min_atj(1:kmx,j)%y))
       iy_min_atj(0,j) = iy_min_atj(k,j)
    end if
    !-- now worry about all points above j
    call renew_iy_min_abovel(p,j,n-1)
    
    !write(0,*) i,j
    if (i>0) then
       ! first redo column i
       if (i>1) call set_iy_min_j(p,i)
       ! then worry about points above i
       call renew_iy_min_abovel(p,i,n-1)
    end if
  end subroutine renew_iy_min_atj


  !-----------------------------------------------------------------
  ! returns the vector of y-values corresponding to tests between
  ! p(:,j) and p(:,imin:imax)
  !
  ! NB: a leading factor (here two/Evis**2) is left out for speed,
  !     and is to be put in later (saving is about 7% on total 
  !     program speed for large numbers of partons).
  function yres(p,j,imin,imax)
    real(dp), intent(in) :: p(:,:)
    integer , intent(in) :: j, imin, imax
    real(dp)             :: yres(imin:imax)
    !--------------------------------------
    integer  :: i
    real(dp) :: costh, dphi, dist, bigX2
    integer  :: jmin, jmax

    !!! IF PUTTING IN DIFFERENT YIJ OPTIONS (E.G. DIFFERENT ALGOS)
    !!! DO SELECT/CASE OUTSIDE LOOP(?) IN ORDER TO PRESERVE 
    !!! SPEED?

    nyres = nyres + imax-imin+1

    select case(iangle_glb)
    case(iangle_ang)
       if (p(i_modp3,j) == zero) then
          yres = zero
          return
       end if
       
       do i = imin, imax
          if (p(i_modp3,i) == zero) then
             yres(i) = zero
          else
             !--- tried some equivalent alternatives...
             !costh = sum(p(:3,i)*p(:3,j))/(p(i_modp3,i)*p(i_modp3,j))
             !costh = (p(1,i)*p(1,j) + p(2,i)*p(2,j) + p(3,i)*p(3,j))&
             !     &/(p(i_modp3,i)*p(i_modp3,j))
             !---------------------------------------------------------
             ! this is the fastest version (with ifort on pentium IV).
             costh = dot_product(p(:3,i),p(:3,j))/(p(i_modp3,i)*p(i_modp3,j))
             yres(i) = min(p(4,i),p(4,j))**2*(one-costh)
          end if
       end do
    case(iangle_DeltaR,iangle_DeltaRs)
       if (p(i_ptsq,j) == zero) then
          yres = zero
          return
       end if
       
       do i = imin, imax
          if (p(i_ptsq,i) == zero) then
             yres(i) = zero
          else
             dphi = abs(p(i_phi,j) - p(i_phi,i))
             if (dphi > pi) dphi = twopi - dphi 
             yres(i) = min(p(i_ptsq,i),p(i_ptsq,j)) * (&
                  &dphi**2 + (p(i_eta,j) - p(i_eta,i))**2)
          end if

          if (iangle_glb == iangle_DeltaRs) then
             if(sum(abs(flavinfo(:,i)+flavinfo(:,j))) > 1) &
                  &yres(i) = y_impossible
          end if
       end do

       
    case(iangle_DeltaR_flav1, iangle_DeltaR_flav2, iangle_DeltaR_flavs1, &
         &iangle_DeltaR_flavs2)

       !if (p(i_ptsq,j) == zero) then
       !   yres = zero
       !   return
       !end if
       
       do i = imin, imax

          ! Compute first the distance DeltaR
          dphi = abs(p(i_phi,j) - p(i_phi,i))
          if (dphi > pi) dphi = twopi - dphi 
          dist = (dphi**2 + (p(i_eta,j) - p(i_eta,i))**2)
          if (p(i_ptsq,i) < p(i_ptsq,j)) then
             jmin = i; jmax = j
          else
             jmin = j; jmax = i
          end if

          !if (max(i,j) == 3 .and. min(i,j) == 1) write(0,*)&
          !     &'Hey there:', dist

          ! put in the flavour correction...
          ! AB ????
          !if (any(flavinfo(:,imin) /= 0)) then
          if (sum(abs(flavinfo(:,jmin))) /= 0) then

             !NB: if dist is zero then leave it as zero so as to
             !    avoid division by zero
             !    (BUT: this is not right for flav2)
             select case(iangle_glb)
             case(iangle_DeltaR_flav1)
                yres(i) = dist * sqrt(p(i_ptsq,jmax)*p(i_ptsq,jmin))
             case(iangle_DeltaR_flav2)
                yres(i) = dist * p(i_ptsq,jmax)
             case(iangle_DeltaR_flavs1)
                yres(i) = dist * sqrt(p(i_ptsq,jmax)*p(i_ptsq,jmin))
                if (sum(abs(flavinfo(:,jmin)+flavinfo(:,jmax))) > 1) &
                     & yres(i) = y_impossible


             case(iangle_DeltaR_flavs2)
                yres(i) = dist * p(i_ptsq,jmax)
                ! START G&G HACKING
                if (sum(abs(flavinfo(:,jmax))) == 1) then
                   !-MODELS B-J
                   if (flavs2_model /= modelA) &
                        &yres(i) = dist * sqrt(p(i_ptsq,jmax)*p(i_ptsq,jmin))
                else
                   bigX2 = min(dist,one)
                   select case(flavs2_model)
                   case(modelD)
                      yres(i) = dist * &
                        &((one-dist)*sqrt(p(i_ptsq,jmax)*p(i_ptsq,jmin)) + &
                        &  dist*p(i_ptsq,jmax))
                   case(modelE)
                      yres(i) = dist * &
                        &((one-bigX2)*sqrt(p(i_ptsq,jmax)*p(i_ptsq,jmin)) + &
                        &  bigX2*p(i_ptsq,jmax))
                   case(modelH)
                      yres(i) = dist * p(i_ptsq,jmax)**(half*(one+bigX2)) * &
                        &           p(i_ptsq,jmin)**(half*(one-bigX2))
                   case(modelF)
                      yres(i) = dist * &
                        &((one-min(sqrt(dist),one))&
                        &*sqrt(p(i_ptsq,jmax)*p(i_ptsq,jmin)) + &
                        &  min(sqrt(dist),one)*p(i_ptsq,jmax))
                   case(modelI)
                      yres(i) = dist * sqrt(p(i_ptsq,jmax)*p(i_ptsq,jmin)) * &
                        &           (one + dist)
                   case(modelJ)
                      yres(i) = dist * sqrt(p(i_ptsq,jmax)*p(i_ptsq,jmin)) * &
                        &           (one + dist*(p(i_ptsq,jmax)-p(i_ptsq,jmin))&
                        &                       /(p(i_ptsq,jmax)+p(i_ptsq,jmin)))
                   end select
                end if
                ! END G&G HACKING
                if (sum(abs(flavinfo(:,jmin)+flavinfo(:,jmax))) > 1) then
                   yres(i) = y_impossible
                end if
             case default
                call wae_error('yres',&
                     &'unrecognized algorithm',intval=iangle_glb)
             end select
          else
             yres(i) = dist * p(i_ptsq,jmin)
          end if
          !if (max(i,j) == 3 .and. min(i,j) == 1) write(0,*)&
          !     &'Hey again:', yres(i)
          
       end do
       

    case default
       write(0,*) 'yres: unsupported iangle = ', iangle_glb
       stop
    end select
    
  end function yres
  
  
!  !-------------------------------------------------------
!  !! return .true. if both particles are beam-particles (i.e.
!  !! have ptsq = ptsq_infinitesimal) and additionally are in the 
!  !! same beam
!  logical function same_beam(i,j)
!    integer, intent(in) :: i, j
!    same_beam = p(i_ptsq,j) == ptsq_infinitesimal .and. &
!         &      p(i_ptsq,i) == ptsq_infinitesimal .and.&
!         &      sign(one,p(3,i)) == sign(one,p(3,j))
!  end function same_beam
  

  !-----------------------------------------------------------------
  ! return the iy combination correspond to the minimum of the vector
  function miniy(vect)
    real(dp), intent(in) :: vect(:)
    type(iy)             :: miniy
    integer :: i(1)
    i = minloc(vect)
    !if (i(1) == 0) write(0,*) 'whoah', vect
    miniy%i = i(1)
    miniy%y = vect(miniy%i)
  end function miniy

  !-------------------------------------------------------------------
  ! return the k value (index of sub divided yij ) for a given
  ! n
  function kofn(n)
    integer, intent(in) :: n
    integer             :: kofn
    kofn = (n-1)/sub_size + 1
  end function kofn

  !---------------------------------------------------------------------
  ! returns the iy_min corresonding to the k,j stretch
  function iy_min_kj(p,k,j)
    real(dp), intent(in) :: p(:,:)
    integer,  intent(in) :: k,j
    type(iy)             :: iy_min_kj
    !-----------------
    integer :: jbase
    jbase = (k-1)*sub_size
    if (-(jbase+1- min(jbase + sub_size, j-1))<0) write(0,*) 'asdfhjlk',jbase,j
    iy_min_kj = miniy(yres(p,j,jbase+1, min(jbase + sub_size, j-1)))
    iy_min_kj%i = iy_min_kj%i + jbase
  end function iy_min_kj
  

  !---------------------------------------------------------------------
  ! returns the iy_min(0:? corresonding to the j column
  subroutine set_iy_min_j(p,j)
    real(dp), intent(in)  :: p(:,:)
    integer,  intent(in)  :: j
    
    !-----------------
    integer :: kmx, k
    kmx = kofn(j-1)
    do k = 1, kmx
       iy_min_atj(k,j) = iy_min_kj(p,k,j)
    end do
    k = sum(minloc(iy_min_atj(1:kmx,j)%y))
    iy_min_atj(0,j) = iy_min_atj(k,j)
  end subroutine set_iy_min_j
  

  !---------------------------------------------------------------------
  ! given a change in the state of parton l, update the parts of 
  ! iy_min_atj(:,l+1:n)
  subroutine renew_iy_min_abovel(p,l,n)
    real(dp), intent(in)    :: p(:,:)
    integer,  intent(in)    :: l,n
    !--------------------------------------------
    integer  :: k, jj, kofl
    real(dp) :: yl(l+1:n)

    yl = yres(p,l,l+1,n)

    kofl = kofn(l)
    do jj = l+1, n
       !-- if l was a local "minimum" previously, then need
       !   to redetermine the local minimum and the global
       !   minimum as well
       if (iy_min_atj(kofl,jj)%i == l) then
          iy_min_atj(kofl,jj) = iy_min_kj(p,kofl,jj)
          k = sum(minloc(iy_min_atj(1:kofn(jj-1),jj)%y))
          iy_min_atj(0,jj) = iy_min_atj(k,jj)
       end if
       !-- check whether new i might become a local minimum
       if (yl(jj) < iy_min_atj(kofl,jj)%y) then
          iy_min_atj(kofl,jj)%y = yl(jj)
          iy_min_atj(kofl,jj)%i = l
          if (iy_min_atj(kofl,jj)%y < iy_min_atj(0,jj)%y)&
               & iy_min_atj(0,jj) = iy_min_atj(kofl,jj)
       end if
    end do

  end subroutine renew_iy_min_abovel
  
  !-- according to the situation being dealt with, set up the "extra"
  !   information in a different manner...
  subroutine setup_p_info(p_col)
    real(dp), intent(inout) :: p_col(:)
    real(dp) :: m2
    real(dp) :: etamax

    etamax = 1e90_dp
    select case(iangle_glb)
    case(iangle_ang)
       p_col(i_modp3) = sqrt(p_col(1)**2 + p_col(2)**2 + p_col(3)**2)
    case(iangle_DeltaR, iangle_DeltaRs, iangle_DeltaR_flav1, &
         &iangle_DeltaR_flav2, iangle_DeltaR_flavs1, iangle_DeltaR_flavs2)
       p_col(i_ptsq) = p_col(1)**2 + p_col(2)**2
       if (p_col(i_ptsq) /= zero) then
          p_col(i_phi) = atan2(p_col(1),p_col(2))
          ! avoid rounding errors this way?
          m2 =  max(zero,p_col(4)**2 - p_col(3)**2 - p_col(i_ptsq))
          p_col(i_eta) = sign(half*log((m2 + p_col(i_ptsq))&
               &/(p_col(4)+abs(p_col(3)))**2), p_col(3))
          !p_col(i_eta) = sign(half*log((p_col(4)**2 - p_col(3)**2)&
          !     &/(p_col(4)+abs(p_col(3)))**2), p_col(3))
          !p_col(i_eta) = half*log((p_col(4)-p_col(3))/(p_col(4)+p_col(3)))
       else
          p_col(i_phi) = zero
          p_col(i_eta) = sign(etamax,p_col(3))
          ! GZ & GPS (23/5/05) -- never allow ptsq to be
          ! truly zero because it leads to confusion in other
          ! contexts(?)
          p_col(i_ptsq) = ptsq_infinitesimal
       end if
    case default
       write(0,*) 'Unrecognized iangle', iangle_glb
    end select
  end subroutine setup_p_info

!OBS  !=================================================================
!OBS  !!
!OBS  !! set up the flavour of all particles, including the beam jets
!OBS  subroutine setup_flav(flav,hh)
!OBS    integer, intent(in) :: flav(:)
!OBS    logical, intent(in) :: hh
!OBS    integer :: i
!OBS    integer :: idqq(2)
!OBS
!OBS    flavinfo = 0
!OBS    if (hh) then
!OBS       ! assume ppbar collision
!OBS       ! -2  = proton
!OBS       flavbeam(:,-2) = (/ 1,2,0,0,0,0 /)
!OBS       ! -1 = antiproton
!OBS       flavbeam(:,-1) = (/-1,-2,0,0,0,0/)
!OBS    end if
!OBS
!OBS    do i = 1, size(flav)
!OBS       ! flavour of quarks
!OBS       if (abs(flav(i)) > 0 .and. abs(flav(i)) <= 6) then
!OBS          flavinfo(abs(flav(i)),i) = sign(1,flav(i))
!OBS      ! flavour of diquarks
!OBS       else if (abs(flav(i)/100) <= 66 .and. &
!OBS            &abs(flav(i)/100) >= 11) then
!OBS           idqq(1) = flav(i)/1000
!OBS           flavinfo(abs(idqq(1)),i) = sign(1,idqq(1))
!OBS           
!OBS           idqq(2) = (flav(i)-idqq(1)*1000)/100
!OBS           if (abs(idqq(1)) == abs(idqq(2))) then
!OBS              flavinfo(abs(idqq(2)),i) = &
!OBS                   & flavinfo(abs(idqq(1)),i)+ sign(1,idqq(2))
!OBS           else
!OBS              flavinfo(abs(idqq(2)),i) = sign(1,idqq(2))
!OBS           end if
!OBS        end if
!OBS    end do
!OBS
!OBS  end subroutine setup_flav
!OBS


  ! redetermine the flavour after each recombination
  subroutine recombine_flavour(i,j,n,hh)
    integer, intent(in)    :: i,j,n
    logical, intent(in)    :: hh

    if (i > 0) then
       flavinfo(:,i) = flavinfo(:,i) + flavinfo(:,j) 
    else
       ! flavour of the beam jets
       if (hh) then
          flavbeam(:,i) = flavbeam(:,i) - flavinfo(:,j) 
       end if
    end if

    if (j < n) then
       flavinfo(:,j)  = flavinfo(:,n)
    end if

  end subroutine recombine_flavour

  
  !======================================================================
  !! order particles in rapidity and returns two auxiliary arrays, to avoid 
  !! computing sqrt's and exp's many times. Note that it also sets the global 
  !! array indx(:) which is such that p(indx(i),i_eta) < p(indx(i+1),i_eta).
  subroutine order_pjet(p,pt,expdeta)
    real(dp), intent(in)  :: p(:,:)
    real(dp), intent(out) :: pt(:),expdeta(:)
    integer :: i,n

    n = size(p,dim=2)

    ! order particles in rapidity
    call indexx(p(i_eta,:), indx(1:n))
    
    pt(1) = sqrt(p(i_ptsq,indx(1)))
    expdeta(1) = one

    do i=2,n
       pt(i) = sqrt(p(i_ptsq,indx(i)))
       expdeta(i) = exp(p(i_eta,indx(i-1))-p(i_eta,indx(i)))
    end do

  end subroutine order_pjet

  ! determines the pt of the beam jets according to eq.(8) and (9) in 
  ! jetflav.tex
  ! ptf(i,:) contribution to ptbeam from (forward) particles with eta > eta(i) 
  ! ptb(i,:) contribution to ptbeam from (backward)particles with eta > eta(i) 
  function get_ptbeam(pt,expdeta) result(ptbeam)
    real(dp), intent(in) :: pt(:), expdeta(:)
    real(dp) :: ptbeam(1:size(pt),-2:-1)
    real(dp) :: ptf(1:size(pt),-2:-1)
    real(dp) :: ptb(1:size(pt),-2:-1)
    integer :: i,n, ii    

    n = size(pt)

    ! set up the pt of the beam
    ptf = zero; ptb = zero

    ptb(1,-2) = pt(1); ptf(n,-1) = pt(n)
    do i=2,n
       ii = i
       ptb(ii,-1) = (ptb(ii-1,-1)+ pt(ii-1))*expdeta(ii)
       ptb(ii,-2) = ptb(ii-1,-2)+ pt(ii)

       ii = n-i+1
       ptf(ii,-1) = ptf(ii+1,-1) + pt(ii)
       ptf(ii,-2) =(ptf(ii+1,-2)+pt(ii+1))*expdeta(ii+1)
    end do

    ptbeam = ptf+ptb

    ! GPS & GZ: it is useful to examine sensitivity to the arbitrary
    ! rescaling of ptbeam...
    ptbeam = ndlist_ptr%ptbeam_rescale * ptbeam
  end function get_ptbeam


  !======================================================================
  !! Determines the pt of the beam jets according to equation in CCN26-37.
  !! The momenta are taken to be ordered in rapidity, with pt(i) the 
  !! transverse momentum of particle i and expdeta = exp(eta_{i-1} - eta_i), 
  !! where i corresponds elsewhere to indx(i).
  !!
  !! Our approach here is as follows:
  !! Define PTb(i) = sum_j pt_j exp(eta_i - eta_j) Theta(eta_j - eta_i)
  !! then ptbeam_b(i) = max(maxval(PTb(1:i)), 
  !!                               maxval(PTb(i:n)*exp(eta_i-eta(i:n))) )
  !!                         ^^^                  ^^^
  !!                          |                    |
  !!                    ptb_max_left          ptb_max_right
  !!
  !! To ensure that we have O(N) operations, we do not actually
  !! repeat the maxval for each i, but rather work it out as we go
  !! along.
  function get_ptbeam_ccn26_37(p, pt,expdeta) result(ptbeam)
    real(dp), intent(in)  :: p(:,:)
    real(dp), intent(in), target :: pt(:), expdeta(:)
    real(dp), target :: ptbeam(1:size(pt),-2:-1)
    real(dp) :: ptb(1:size(pt))
    real(dp) :: ptb_max_left(1:size(pt)), ptb_max_right(1:size(pt))
    !real(dp) :: ptb(1:size(pt),-2:-1)
    integer :: i,n, ii    
    real(dp), pointer :: ptbeam_p(:), pt_p(:), expdeta_p(:)
    ! testing
    !T real(dp) :: eta(1:size(pt)), ptbeam_here, eta_here, ptf(size(pt)), ptfeam_here

    n = size(pt)

    do ii = 1, 2
       if (ii == 1) then
          ! backwards
          ptbeam_p  => ptbeam(1:n,-2)
          expdeta_p => expdeta(2:n)  ! only 2:n contains anything useful
          pt_p      => pt(1:n)
       else
          ! forwards (reverse it all!)
          ptbeam_p  => ptbeam(n:1:-1,-1)
          expdeta_p => expdeta(n:2:-1) ! only 2:n contains anything useful
          pt_p      => pt(n:1:-1)
       end if
       
       ! calculate quantities as if for beam towards negative eta
       ! (above pointers will sort out other case too)
       ptb(n) = pt_p(n)
       ptb_max_right(n) = ptb(n)
       do i = n-1, 1, -1
          ptb(i) = ptb(i+1) * expdeta_p(i) + pt_p(i)
          ptb_max_right(i) = max(ptb_max_right(i+1)*expdeta_p(i), ptb(i))
       end do
       ptb_max_left(1) = ptb(1)
       do i = 2, n
          ptb_max_left(i) = max(ptb_max_left(i-1), ptb(i))
       end do
       ptbeam_p = max(ptb_max_left, ptb_max_right)
    end do

    ! GPS & GZ: it is useful to examine sensitivity to the arbitrary
    ! rescaling of ptbeam...
    ptbeam = ndlist_ptr%ptbeam_rescale * ptbeam
    
    !T ! testing -- doing things the SLOW way!
    !T eta(:) = p(i_eta,indx(:))
    !T forall(i=1:n)
    !T    ptb(i) = sum(pt(i:n)*exp(eta(i) - eta(i:n)))
    !T    ptf(i) = sum(pt(1:i)*exp(eta(1:i) - eta(i)))
    !T end forall
    !T 
    !T do i = 1, n
    !T    ptbeam_here = max(maxval(ptb(1:i)),&
    !T         &      maxval(ptb(i:n)*exp(eta(i)-eta(i:n))))
    !T    !write(0,'(i4,4f10.4)') i, -sum(log(expdeta(2:i))), pt(i), ptbeam(i,:)
    !T    write(0,'(i4,6es12.4)') i, expdeta(i), p(i_eta,indx(i)), pt(i), ptbeam(i,:), ptbeam_here
    !T end do
    !T 
    !T do i = 0, 200
    !T    eta_here = -6.0_dp + i*(12.0_dp/200)
    !T    ptbeam_here = max(maxval(ptb(:),mask=eta<=eta_here),&
    !T         &      maxval(ptb(:)*exp(eta_here-eta(:)),mask=eta>=eta_here))
    !T    ptfeam_here = max(maxval(ptf(:),mask=eta>=eta_here),&
    !T         &      maxval(ptf(:)*exp(eta(:)-eta_here),mask=eta<=eta_here))
    !T    write(0,'(5es15.7)') eta_here, ptbeam_here, ptfeam_here,&
    !T         &sum(pt,eta<eta_here)+sum(pt*exp(eta_here-eta),eta>=eta_here), &
    !T         &sum(pt,eta>eta_here)+sum(pt*exp(eta-eta_here),eta<=eta_here)
    !T end do
    !T 
    !T stop
  end function get_ptbeam_ccn26_37


  subroutine set_ybeam(p)
    real(dp), intent(in)  :: p(:,:)
    real(dp):: pt(size(p,dim=2)), expdeta(size(p,dim=2))
    real(dp):: ptbeam(1:size(p,dim=2),-2:-1)
    integer :: i, n

    n = size(p,dim=2)

    call order_pjet(p, pt, expdeta)

    if (ndlist_ptr%use_local_ptbeam) then
       ptbeam = get_ptbeam_ccn26_37(p, pt, expdeta)
    else
       ptbeam = get_ptbeam(pt, expdeta)
    end if

    ! initialise ybeam in such a way that particles with negative (positive) 
    ! rapidities are preferentially recombined with the proton (antiproton)
    !forall(i=1:n)
    do i=1, n
       !ybeam(i,-2) = p(i_ptsq,indx(i))*&
       !     &(one+y_impossible*(one+sign(one,p(i_eta,indx(i)))))
       !ybeam(i,-1) = p(i_ptsq,indx(i))*&
       !     &(one+y_impossible*(one-sign(one,p(i_eta,indx(i)))))
       ybeam(i,:) = p(i_ptsq,indx(i))
       !if (p(i_eta,indx(i)) >= 0) then
       if (ptbeam(i,-2) > ptbeam(i,-1)) then
          ybeam(i,-2) = y_impossible
       else
          ybeam(i,-1) = y_impossible
       end if
    end do
    !end forall
    

    ! correct ybeam according to the algorithm
    if (use_ptbeam_glb) then
       call ybres(p,ptbeam,pt,1,n)
    end if

    ! GPS-R: include R scale (include it in beam, and then 
    ! later, when recording things in ndlist, divide it out from
    ! yiB and yij, so that one is left with 1/R2 just in the yij,
    ! as per the original definition by Ellis & Soper).
    ybeam(:n,:) = ybeam(:n,:) * ndlist_ptr%R2
  end subroutine set_ybeam

  ! Update ybeam after a recombination. For the moment it just recomputes 
  ! ybeam using the new momenta after a recombination
  subroutine renew_ybeam(p,i,j,n)
    real(dp), intent(in) :: p(:,:)
    integer, intent(in)  :: i,j,n
    
    call set_ybeam(p(:,:n-1))

  end subroutine renew_ybeam

  subroutine ybres(p, ptbeam, pt,imin,imax)
    real(dp), intent(in) :: p(:,:)
    real(dp), intent(in) :: ptbeam(:,-2:), pt(:)
    integer, intent(in)  :: imin, imax
    integer              :: i, ibeam
    integer :: flavsum(-2:-1)
    logical :: isbeam(imin:imax)

    ! only do this if the particle has a decent pt -- otherwise
    ! we consider it to be "predestined" to its beam.
    isbeam(:) = p(i_ptsq,indx(imin:imax))== ptsq_infinitesimal

    select case(iangle_glb)
    case(iangle_DeltaR)
    case(iangle_DeltaRs)
       do i=imin,imax
          if (isbeam(i)) cycle
          do ibeam = -2, -1
             flavsum(ibeam) = &
                  &sum(abs(flavbeam(:,ibeam)-flavinfo(:,indx(i))),dim=1)
          end do
          ! If any of the beams is impossible because of flavour issues
          ! we need to reconsider our purely kinematic choice of beam (recall
          ! the more distant beam was already labelled impossible). This
          ! should hopefully avoid both beam recombinations beig labelled
          ! impossible (unless they truly are for flavour reasons).
          !
!          if (any(flavsum > 1) .and. p(i_ptsq,indx(i))/=ptsq_infinitesimal) then
          if (any(flavsum > 1)) then
             do ibeam = -2, -1
                !if (flavsum(ibeam) > 1 .and. &
                !     &ndlist_ptr%p(i_ptsq,indx(i)) /= ptsq_infinitesimal) then
                if (flavsum(ibeam) > 1) then
                   ybeam(i,ibeam) = y_impossible
                else
                   ybeam(i,ibeam) = pt(i)**2
                end if
             end do
          end if
       end do
    case(iangle_DeltaR_flav1)
       do i=imin,imax
          if (sum(abs(flavinfo(:,indx(i))),dim=1) /= 0) then
             ybeam(i,-2:) = pt(i)*ptbeam(i,-2:)
          end if
       end do
    case(iangle_DeltaR_flav2)
       do i=imin,imax
          if (sum(abs(flavinfo(:,indx(i))),dim=1) /= 0) then
             ybeam(i,-2:) = ptbeam(i,-2:)**2
          end if
       end do
    case(iangle_DeltaR_flavs1)
       do i=imin,imax
          if (isbeam(i)) cycle
          if (sum(abs(flavinfo(:,indx(i))),dim=1) /= 0) then
             do ibeam = -2,-1
                if (sum(abs(flavbeam(:,ibeam)-flavinfo(:,indx(i))),dim=1) > 1)&
                     & then 
                   ybeam(i,ibeam) = y_impossible
                else
                   ybeam(i,ibeam) = pt(i)*ptbeam(i,ibeam)
                end if
             end do
          end if
       end do
    case(iangle_DeltaR_flavs2)
       do i=imin,imax
          if (isbeam(i)) cycle
          ! START G&G HACKING (UNLABELLED MODEL!)
          !ybeam(i,:) = pt(i)*ptbeam(i,:)
          ! END G&G HACKING
          if (sum(abs(flavinfo(:,indx(i))),dim=1) /= 0) then
             do ibeam = -2,-1
                if (sum(abs(flavbeam(:,ibeam)-flavinfo(:,indx(i))),dim=1) > 1)&
                     & then 
                   ybeam(i,ibeam) = y_impossible
                else
                   ybeam(i,ibeam) = ptbeam(i,ibeam)**2
                   select case(flavs2_model) 
                   case(modelC:)
                      ! START G&G HACKING (MODELS C-I)
                      ybeam(i,ibeam) = pt(i)*ptbeam(i,ibeam)
                      ! END G&G HACKING
                   end select
                end if
             end do
          end if
       end do
    case default
       call wae_error('dpktclus_gen','Algorithm not yet implemented')
    end select
  end subroutine ybres


  !======================================================================
  !! Calculates and writes out ptbeam; since this calls the set_ptbeam 
  !! routine, it shares some side-effects with that routine (e.g. setting
  !! up indx array).
  subroutine writeout_ptbeam_points(p)
    real(dp), intent(in)  :: p(:,:)
    integer :: i, n
    real(dp):: pt(size(p,dim=2)), expdeta(size(p,dim=2))
    real(dp):: ptbeam(1:size(p,dim=2),-2:-1)

    n = size(p,dim=2)

    call order_pjet(p, pt, expdeta)

    if (ndlist_ptr%use_local_ptbeam) then
       ptbeam = get_ptbeam_ccn26_37(p, pt, expdeta)
    else
       ptbeam = get_ptbeam(pt, expdeta)
    end if

    write(6,*) '# Results obtained by calling writeout_ptbeam_points &
         &(probably hacked in by hand)'
    write(6,*) '# eta, pt, ptbeam_1, ptbeam_2'
    do i = 1, n
       write(6,'(2i5,6es15.6E3)') i, indx(i), p(i_eta, indx(i)), &
            &       sqrt(p(i_ptsq, indx(i))), &
            &       ptbeam(i,-1), ptbeam(i,-2)
    end do
    write(6,*)
    write(6,*)
    
  end subroutine writeout_ptbeam_points


  !======================================================================
  !! Calculates and writes out ptbeam as a continuous function of eta; 
  !! since this calls the order_pjet
  !! routine, it shares some side-effects with that routine (e.g. setting
  !! up indx array).
  !!
  !! Note that it is EXTREMELY inefficient since it recalculates the
  !! componentss of ptbeam (as well as kt histograms) with an N^2 algorithm!
  subroutine writeout_ptbeam_curves(p)
    use sub_defs_io
    real(dp), intent(in)  :: p(:,:)
    integer :: i, n
    real(dp):: pt(size(p,dim=2)), expdeta(size(p,dim=2))
    real(dp):: ptbeam(1:size(p,dim=2),-2:-1)
    !--------------------------------------
    integer,  parameter :: neta = 120
    real(dp), parameter :: eta_min=-6.0_dp, eta_max = 6.0_dp
    real(dp), parameter :: deta = (eta_max - eta_min)/neta
    real(dp), parameter :: first_last_eta = -1000.0_dp
    real(dp) :: eta, ptbeam_curve(0:neta,-2:-1)
    real(dp) :: kthist(0:neta)
    integer  :: ieta
    real(dp), save :: last_eta = first_last_eta
    integer,  save :: iev = -1

    if (last_eta == first_last_eta) then
       rewind(6)
       write(6,'(a)') '# '//trim(command_line())
    end if
    
    if (p(i_eta,2) == last_eta) then
       return
    else
       last_eta = p(i_eta,2)
       iev = iev + 1
    end if
     
    
    n = size(p,dim=2)

    call order_pjet(p, pt, expdeta)

    write(6,'(a)') '# Results obtained calling writeout_ptbeam_curves &
         &(probably hacked in by hand)'
    write(6,'(a)') '# eta, ptbeam_1, ptbeam_2, pt (in bin eta-deta..eta)'
    write(6,'(a,i5)') '# iev = ', iev

    do ieta = 0, neta
       eta = ieta*deta + eta_min
       ptbeam_curve(ieta,-1) = sum(pt,p(i_eta,indx(:)) > eta) + &
            &   sum(pt*exp(p(i_eta,indx(:))-eta),p(i_eta,indx(:)) <= eta)
       ptbeam_curve(ieta,-2) = sum(pt,p(i_eta,indx(:)) < eta) + &
            &   sum(pt*exp(eta-p(i_eta,indx(:))),p(i_eta,indx(:)) >= eta)
       if (ieta > 0) then
          kthist(ieta) = sum(pt,(p(i_eta,indx(:)) < eta .and. &
                     &           p(i_eta,indx(:)) >= eta-deta))
       else
          kthist(ieta) = 0
       end if
       
       write(6,'(6es15.6E3)') eta, &
            &                 ptbeam_curve(ieta,-1), ptbeam_curve(ieta,-2), &
            &                 kthist(ieta)
    end do
    write(6,*)
    write(6,*)
    
  end subroutine writeout_ptbeam_curves
  
end module dpktclus_gen
