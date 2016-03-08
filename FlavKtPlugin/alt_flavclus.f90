!======================================================================
! Module for an alternative implementation of flavour clustering in 
! hadron-hadron, intended as a cross-check on current routines.
!======================================================================
module alt_flavclus
  use types; use consts_dp
  implicit none
  private

  
  real(dp), allocatable :: p(:,:), phi(:), eta(:), ktsq(:)
  real(dp), allocatable :: ptbeam(:,:)
  integer , allocatable :: flav(:,:)

  real(dp), parameter :: y_huge = 1e200_dp
  real(dp), parameter :: etamax = 1e50_dp
  real(dp), parameter :: ktsq_infinitesimal = 1e-200_dp

  integer, parameter :: algo_kt = 1
  integer, parameter :: algo_flav1 = 4
  integer, parameter :: algo_flav2 = 5
  integer, parameter :: algo_flavs1 = 6
  integer, parameter :: algo_flavs2 = 7
  integer, parameter :: algo_kts = 8
  integer, save      :: algo = algo_kt

  public :: do_alt_flavclus
contains

  !-------------------------------------------------
  ! bbeamflav is flavour of incoming beam. 
  ! -2 is beam whose remnant is at negative infinity
  ! -1 is beam whose remnant is at positive infinity
  !
  ! Note that internally we will use -bbeamflav (since in general 
  ! out internal flavour is always an outgoing flavour).
  !
  subroutine do_alt_flavclus(pp,fflav,bbeamflav,imode)
    real(dp), intent(in) :: pp(:,:)
    integer,  intent(in) :: fflav(:,:), bbeamflav(:,-2:), imode
    !-----------------------------------
    integer  :: i, j, n, ibeam
    integer  :: ii, jj
    real(dp) :: ymin, y

    select case(imode)
    case(4203)
       algo = algo_kt
    case(4303)
       algo = algo_flav1
    case(4403)
       algo = algo_flav2
    case(4503)
       algo = algo_flavs1
    case(4603)
       algo = algo_flavs2
    case(4703)
       algo = algo_kts
    case default
       stop "Unsupported imode"
    end select
    
    allocate(p(size(pp,dim=1),size(pp,dim=2)))
    allocate(flav(size(fflav,dim=1),-2:size(fflav,dim=2)))
    allocate(phi(size(pp,dim=2)),eta(size(pp,dim=2)))
    allocate(ktsq(size(pp,dim=2)))
    allocate(ptbeam(-2:-1,size(pp,dim=2)))

    p    = pp
    flav(:,1:) = fflav
    flav(:,-2:-1) = -bbeamflav(:,-2:-1)

    do n = size(p,dim=2), 1, -1
       ktsq(:n) = p(1,:n)**2 + p(2,:n)**2
       where (ktsq/= zero)
          phi(:n) = atan2(p(1,:n),p(2,:n))
          eta(:n) = half*log((p(4,:n)+p(3,:n))/(p(4,:n)-p(3,:n)))
       elsewhere
          phi(:n) = zero
          eta(:n) = sign(etamax,p(3,:n))
          ktsq(:n) = ktsq_infinitesimal
       end where
       call calculate_ptbeam(n)
       !if (n == 2) then
       !   write(0,'(6i1," ",6i1," ",6i1," ",6i1," ")') &
       !        &flav(:,-2), flav(:,-1),flav(:,1), flav(:,2)
       !   write(0,*) ktsq(1:2)
       !   write(0,*) eta(1:2)
       !   write(0,'(4es20.10)') ptbeam(-2,1),ptbeam(-2,2),&
       !        &                ptbeam(-1,1),ptbeam(-1,2)
       !   write(0,'(4es20.10)') yibeam(1,-2),yibeam(2,-2),&
       !        &               yibeam(1,-1),yibeam(2,-1)
       !   write(0,*) yibeam(1,-2) - yibeam(2,-1)
       !end if
       
       ! find smallest yij
       ymin = y_huge
       do i = 1, n-1
          do j = i+1,n
             y = yij(i,j)
             if (y <= ymin) then
                ii = i ; jj = j; ymin = y
             end if
          end do
       end do
       
       ! see if there is a smaller yibeam
       do ibeam = -1,-2,-1
          do j = 1, n
             y = yibeam(j,ibeam)
             if (y <= ymin) then
                ii = ibeam; jj = j; ymin = y
             end if
          end do
       end do

       !-- do a dumb recombination...
       if (ii > 0) then
          p(:,ii) = p(:,ii) + p(:,jj)
       end if
       ! NB: flavour recomb holds even for beam
       flav(:,ii) = flav(:,ii) + flav(:,jj)
       p(:,jj)    = p(:,n)
       flav(:,jj) = flav(:,n)

    end do

    
    deallocate(p)
    deallocate(flav)
    deallocate(phi,eta)
    deallocate(ktsq,ptbeam)
  end subroutine do_alt_flavclus
  

  !-------------------------------------------------
  real(dp) function yij(i,j) 
    integer,  intent(in) :: i, j
    !-----------
    integer :: jmin, jmax
    real(dp) :: dphi,deta, dist
    dphi = abs(phi(i) - phi(j))
    dphi = min(dphi, twopi - dphi)
    deta = (eta(i) - eta(j))
    dist = deta**2 + dphi**2

    if (ktsq(i) < ktsq(j)) then
       jmin = i; jmax = j
    else
       jmin = j; jmax = i
    end if
    
    if (algo /= algo_kt .and. any(flav(:,jmin) /= 0)) then
       select case(algo)
       case(algo_flav1, algo_flavs1)
          yij = sqrt(ktsq(jmin)*ktsq(jmax))*dist
       case(algo_flav2, algo_flavs2)
          yij = ktsq(jmax)*dist
       case(algo_kts)
          yij = ktsq(jmin)*dist
       end select

       select case(algo)
       case (algo_flavs1, algo_flavs2, algo_kts)
          if (sum(abs(flav(:,i)+flav(:,j))) > 1) yij = y_huge
       end select

    else
       yij = ktsq(jmin)*dist
    end if

    
  end function yij
  

  !-------------------------------------------------
  real(dp) function yibeam(i,ibeam)
    integer,  intent(in) :: i, ibeam
    real(dp) :: flavsum(-2:-1)
    integer  :: ib

    if (algo /= algo_kt.and. any(flav(:,i) /= 0)) then
       select case(algo)
       case(algo_flav1, algo_flavs1)
          yibeam = sqrt(ktsq(i))*ptbeam(ibeam,i)
       case(algo_flav2, algo_flavs2)
          yibeam = ptbeam(ibeam,i)**2
       case(algo_kts)
          do ib = -2, -1
             flavsum(ib) = sum(abs(flav(:,i)+flav(:,ib)))
          end do
          if (any(flavsum > 1)) then
             if (flavsum(ibeam) > 1) then
                yibeam = y_huge
             else
                yibeam = ktsq(i)
             end if
          else
             yibeam = yibeam_plainkt(i,ibeam)
          end if
       end select

       if (algo == algo_flavs1 .or. algo == algo_flavs2) then
          if (sum(abs(flav(:,i)+flav(:,ibeam))) > 1) yibeam = y_huge
       end if
    else
       yibeam = yibeam_plainkt(i,ibeam)
    end if
  end function yibeam

  !------------------------------------------------
  real(dp) function yibeam_plainkt(i,ibeam)
    integer,  intent(in) :: i, ibeam
    logical :: thisbeam
    thisbeam = (ibeam==-1).eqv.(ptbeam(-1,i) < ptbeam(-2,i))
    if (thisbeam) then
       yibeam_plainkt = ktsq(i)
    else
       yibeam_plainkt = y_huge
       end if
  end function yibeam_plainkt


  !-------------------------------------------------
  subroutine calculate_ptbeam(n)
    integer, intent(in) :: n
    integer :: i, j
    ptbeam = zero
    do i = 1, n
       do j = 1, n
          if (eta(j) < eta(i)) then
             ptbeam(-2,i) = ptbeam(-2,i) + sqrt(ktsq(j))
             ptbeam(-1,i) = ptbeam(-1,i) + sqrt(ktsq(j))*exp(eta(j)-eta(i))
          else
             ptbeam(-1,i) = ptbeam(-1,i) + sqrt(ktsq(j))
             ptbeam(-2,i) = ptbeam(-2,i) + sqrt(ktsq(j))*exp(eta(i)-eta(j))
          end if
       end do
    end do
    
  end subroutine calculate_ptbeam
  
end module alt_flavclus
