module ktflav
  use types; use consts_dp; use warnings_and_errors; use dpktclus_gen
  use alg_declarations 
  implicit none 
  private

  public :: labelMjets, labelInclusiveJets

  !NB: in various places we assume that a label of 0 has no other
  !    meaning -- this is potentially dangerous?

contains

  !--------------------------------------------------
  ! Resolve the event into m (non-beam) jets and label each 
  ! particle according to the jet it belongs to.
  !
  ! For nodes where parent 1 is < 0, the particles associated
  ! with parent2 are all associated to a jet whose number=parent1
  ! (i.e. if prnt1 = -1, then all particles that go into
  ! that node get jet number = -1) -- this helps us reconstruct
  ! the beam jets, i.e. it gives us the exclusive kt-algorithm
  subroutine labelMjets(m, nodes,labels)
    integer,    intent(in)  :: m
    type(node), intent(in)  :: nodes(:)
    integer,    intent(out) :: labels(:)
    integer :: ijet, nn!, modm
    
    if (size(labels) /= size(nodes)) call wae_error('labelMjets',&
         &'size(labels) /= size(nodes)')
    if (m < 2) call wae_error('labelMjets',&
         &'m must be >=2. it was',intval=m)

    ! will help sanity checks...
    labels = 0

    ijet = 0

    do nn = 1, m
       if (nodes(nn)%prnt2 > m) then
          ijet = ijet + 1
          call RecLabel(ijet,nodes(nn)%prnt2,nodes,labels)
       end if
       if (nodes(nn)%prnt1 > m) then
          ijet = ijet + 1
          call RecLabel(ijet,nodes(nn)%prnt1,nodes,labels)
       end if
    end do

    ! now label remaining partons that make up beam jets...
    do nn = m+1, size(nodes)
       if (nodes(nn)%prnt1 < 0) then
          !if (nn <= modm) modm = modm + 1
          call RecLabel(nodes(nn)%prnt1,nodes(nn)%prnt2,nodes,labels)
       end if
    end do
    
    ! sanity checks (maybe too expensive?)
    if (ijet /= m) call wae_error('labelMjets',&
         &'Found inconsistent number of jets.')
    if (any(labels == 0)) call wae_error('labelMjets',&
         &'Some partons not labelled')
  end subroutine labelMjets


  !-------------------------------------------------------
  ! Resolve the event into inclusive jets and label each 
  ! particle according to the jet it belongs to.
  subroutine labelInclusiveJets(nodes,labels,nincl)
    type(node), intent(in)  :: nodes(:)
    integer,    intent(out) :: labels(:)
    integer,    optional, intent(out) :: nincl
    integer :: ijet, nn

    if (size(labels) /= size(nodes)) call wae_error('labelInclusiveJets',&
         &'size(labels) /= size(nodes)')

    ! will help sanity checks...
    labels = 0

    ijet = 0
    do nn = 1, size(labels)
       if (nodes(nn)%prnt1 < 0) then
          ijet = ijet + 1
          call RecLabel(ijet,nodes(nn)%prnt2,nodes,labels)
       end if
    end do
    if (present(nincl)) nincl = ijet
  end subroutine labelInclusiveJets
  


  !--------------------------------------------------------------
  ! Follow node inode and label all subsiduary particles as being
  ! part of jet ijet
  recursive subroutine RecLabel(ijet,inode,nodes,labels)
    integer,    intent(in)    :: ijet, inode
    type(node), intent(in)    :: nodes(:)
    integer,    intent(inout) :: labels(:)
    !--------------
    integer :: ipart, nn
    
    nn = size(nodes)
    if (inode > nn) then
       ipart = inode - nn 
       if (labels(ipart) /= 0) call wae_error('RecLabel',&
            &'labels(ipart/=0) for ipart',intval=ipart)
       labels(ipart) = ijet
    else
       call RecLabel(ijet,nodes(inode)%prnt1,nodes,labels)
       call RecLabel(ijet,nodes(inode)%prnt2,nodes,labels)
    end if
  end subroutine RecLabel

end module ktflav
