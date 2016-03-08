

!!
!! Cluster an event according to the "imode" algorithm.
!!
!! The event is taken to have np outgoing partices with momenta
!! p(1:4,1:np), and up to n different flavours for particles, stored
!! in flav(1:nflav,1:np). 
!!
!! The beam flavours are given in beamflav(:nflav,-2:-1), with the
!! assumption that beam(:,-2) has negative pz, while beam(:,-1) has
!! positive pz.
!!
!! If njet > 0, the event is clustered to exactly njet outgoing jets
!! (this may change in the future) and njet is unchanged at the end;
!! dcut is set to the (unnormalised) 3-jet resolution 
!!
!! If njet == 0 the event is clustered up to the scale dcut and on output
!! njet is set equal to the number of outgoing jets.
!!
!! if njet < 0, an inclusive clustering is carried out, and on output
!! njet is set to the number of outgoing jets, while dcut is left 
!! unchanged
!!
!!
!! 
!! the jet momenta are stored on output in jetpp(1:4,1:njet), the jet
!! flavours in jetflav(1:nflav,1:njet) the incoming flavours in
!! incomingflav(1:nflav,-1:-2).
!!
!! d3 is set to the 3-jet resolution parameter (not normalised)
!!
subroutine ktflavf77(imode, rparam,np, pp, nflav, flav, beamflav, &
     &                           njet, jetpp, jetflav, incomingflav, dcut)
  use types; use consts_dp
  use warnings_and_errors
  use dpktclus_gen
  use ktflav
  implicit none
  !-- input
  integer,  intent(in)  :: imode, np
  real(dp), intent(in)  :: rparam, pp(4,np)
  integer,  intent(in)  :: nflav, flav(nflav,np), beamflav(nflav,-2:-1)
  !-- output
  integer,  intent(inout) :: njet
  real(dp), intent(out) :: jetpp(4,np)
  integer,  intent(out) :: jetflav(nflav,np), incomingflav(nflav,-2:-1)
  real(dp), intent(inout) :: dcut
  !----------------------------------------
  type(nodelist) :: ndlist
  integer        :: labels(np)
  integer        :: ip, ijet

  call ktclus_genR(pp,ndlist,imode,r=rparam,flavour=flav,beamflavour=beamflav)
  
  ! get two exclusive jets
  if (njet > 0) then
     call labelMjets(njet, ndlist%nodes, labels)
     ! get the n-jet resolution parameter (unnormalized)
     dcut = get_ym(ndlist, njet+1)
  else if (njet == 0) then
     njet = get_m(ndlist, dcut)
     call labelMjets(njet, ndlist%nodes, labels)
  else ! njet < 0, inclusive case
     call labelInclusiveJets(ndlist%nodes,labels,njet)
  end if
  

  jetpp(:,:njet)   = zero
  jetflav(:,:njet) = zero
  incomingflav     = beamflav

  do ip = 1, np
     ijet = labels(ip)
     if (ijet > 0) then
        jetpp(:,ijet)   = jetpp(:,ijet)   + pp(:,ip)
        jetflav(:,ijet) = jetflav(:,ijet) + flav(:,ip)
     else if (ijet < 0) then
        incomingflav(:,ijet) = incomingflav(:,ijet) - flav(:,ip)
     else 
        call wae_error('ktflavf77', &
             &'ijet=0 is internally inconsistent')
     end if
     
  end do

  
  call delete_nodelist(ndlist)
end subroutine ktflavf77

!=================================================================
!!
!! Filters the flavour array so that all flavours except id are set to
!! zero.
!!
!! It assumes that flavour array (dim 1) starts from index 1.
subroutine filterflavour(np,nflav,flav,id)
  use warnings_and_errors
  integer, intent(in)    :: np, nflav
  integer, intent(inout) :: flav(nflav,np)
  integer, intent(in)    :: id
  integer :: i
  if (id <= 0 .or. id > size(flav,dim=1)) &
       &call wae_error('filterflavour_',&
       &'id value is out of range of size of flavour; id = ',intval=id)
  do i = 1, nflav
     if (i == id) cycle
     flav(i,:) = 0
  end do
end subroutine filterflavour


!======================================================================
module ktf77_storage
  use types; use consts_dp
  use warnings_and_errors
  use dpktclus_gen
  use ktflav
  implicit none
  type(nodelist), save :: ndlist
end module ktf77_storage


!======================================================================
subroutine ktf77MakeNodeList(imode, rparam, np, pp)
  use ktf77_storage
  implicit none
  integer,  intent(in)  :: imode, np
  real(dp), intent(in)  :: rparam, pp(4,np)
  !integer,  intent(in)  ::nflav, flav(nflav,np), beamflav(nflav,-2:-1)
  call ktclus_genR(pp,ndlist,imode,r=rparam)
end subroutine ktf77MakeNodeList

!======================================================================
subroutine ktf77MakeNodeListFlav(imode, rparam, np, pp, nflav, flav, beamflav)
  use ktf77_storage
  implicit none
  integer,  intent(in)  :: imode, np
  real(dp), intent(in)  :: rparam, pp(4,np)
  integer,  intent(in)  :: nflav
  integer,  intent(in)  :: flav(nflav,np), beamflav(nflav,-2:-1)
  call ktclus_genR(pp,ndlist,imode,r=rparam, flavour=flav, beamflavour=beamflav)
end subroutine ktf77MakeNodeListFlav

!----------------------------------------------------------------------
subroutine ktf77DeleteNodeList()
  use ktf77_storage
  implicit none
  call delete_nodelist(ndlist)
end subroutine ktf77DeleteNodeList

!----------------------------------------------------------------------
subroutine ktf77ReadNodeList(i,p,parent1,parent2,y)
  use ktf77_storage
  implicit none
  integer,  intent(in)  :: i
  real(dp), intent(out) :: p(4)
  integer,  intent(out) :: parent1, parent2
  real(dp), intent(out) :: y
  p = ndlist%p(:,i)
  parent1 = ndlist%nodes(i)%prnt1
  parent2 = ndlist%nodes(i)%prnt2
  y     = ndlist%nodes(i)%y
end subroutine ktf77ReadNodeList
