      subroutine i_setlocalscales(
C - Inputs (first the serious ones)
     $      flg_bornonly,    ! Are we doing an LO (flg_bornonly=1) or NLO
     $                       ! computation (flg_bornonly=0). LO means only Born.
     $      imode,           ! Set imode =1 for Born, imode=2 for all NLO contribs
     $      isReal,          ! Set isReal=1 for real kinematics, 0 otherwise.
     $      nlegborn,        ! Number of f.s. particles in the Born PLUS 2 (ZJJ=6) 
     $      cmp,             ! nlegborn/nlegborn+1 momenta in C.O.M frame
     $                       ! Format: cmp(0-3,i), for i'th particle 0-3=E,px,py,pz
     $                       ! Particle order shouldnt matter but we adopt convention
     $                       ! incoming +/- (1/2), leptons (3,4), partons (5,6,...)) 
     $      flav,            ! nlegborn/nlegborn+1 flavours; same order as momenta.
     $      raisingscales,   ! Freeze clustering scale if it attempts to
     $                       ! go down - no more sudakovs/subtractions, only
     $                       ! coupling reweighting (for all remaining nodes).
     $                       ! (def: raisingscales=1 i.e. do this)
     $      st_alpha,        ! Value of fixed aS used in MEs before getting here
     $      st_muren2,       ! The mu_R^2 value that was used to evaluate virt
     $      st_bornorder,    ! Number of powers of aS in the Born (ZJ=1,ZJJ=2)
C - Input used for computing aS in this code
     $      rad_charmthr2,   ! Charm  (mass) threshold squared. def=(1.5 GeV)^2
     $      rad_bottomthr2,  ! Bottom (mass) threshold squared. def=(5.0 GeV)^2
     $      st_nlight,       ! Nominal number of light flavs in computation (def=5)
     $      st_lambda5MSB,   ! Lambda_5 in the MS bar scheme 
C - Input factors used to recale mu_R and mu_F if performing scale variations
     $      st_renfact,      ! Scaling factor for mu_R (not mu_R^2): 0.5,1 (def),2.0
     $      st_facfact,      ! Scaling factor for mu_F (not mu_F^2): 0.5,1 (def),2.0
C - Outputs
     $      rescfac,         ! Comb. of sudakovs, subtractions & aS weights to
     $                       ! rescale the weight of the contrib assoc. to this
     $                       ! mom. and flav. configuration.
     $      st_mufact2,      ! The MiNLO mu_F^2 value that the PDFs should be 
     $                       ! (re-)evaluated at now
C - Additional output information on make-up of rescfac for debugging
     $      basicfac,        ! Comb. of suds & aS ratios multiplying all contribs
     $      bornfac,         ! Comb. of subtractions to multiply Born term only by 
     $      nlofac)          ! Multiplies only NLO terms (adjusts NLO term's aS).
     $                       ! So, for Born terms rescfac=basicfac*bornfac
     $                       !     for NLO  terms rescfac=basicfac*nlofac
     $                       ! If flg_bornonly=1 we are considering only LO
     $                       ! so then for Born terms rescfac=basicfac only.
c$$$     $      indie_ub,
c$$$     $      pout)
      implicit none

C - maxprocborn and maxprocreal must be greater or equal to the number of
C - independent flavour structures for the born and real process.
      integer  maxprocborn
      parameter(maxprocborn=999)
      integer  nlegborn
      integer  isReal
      real * 8 st_alpha
      integer  st_bornorder
      real * 8 st_muren2
      integer  flg_bornonly
      real * 8 rad_charmthr2
      real * 8 rad_bottomthr2
      integer  st_nlight
      real * 8 st_lambda5MSB
      real * 8 st_renfact
      real * 8 st_facfact
      integer  raisingscales
      real * 8 cmp(0:3,20)

      integer  iuborn,imode
      real * 8 rescfac
      real * 8 st_mufact2

      integer  flav(20),indie_ub(20)
      real * 8 pout(0:3,20)
      real * 8 basicfac,bornfac,nlofac

c$$$CXXXX DEBUGGING XXXXC
      if(.false.) then
         write(*,*) ''
         write(*,*) ''
         write(*,*) ''
         write(*,*) ''
         write(*,*) 'New born event! Top level! (indy setlocalscales)'
         write(*,*) 'nlegborn ',nlegborn
         write(*,*) 'st_alpha ',st_alpha
         write(*,*) 'st_bornorder ',st_bornorder
         write(*,*) 'st_muren2 ',st_muren2
         write(*,*) 'rad_charmthr2,rad_bottomthr2 ',
     $               rad_charmthr2,rad_bottomthr2
         write(*,*) 'st_nlight,st_lambda5MSB ',st_nlight,st_lambda5MSB 
         write(*,*) 'st_renfact,st_facfact,raisingscales ',
     $               st_renfact,st_facfact,raisingscales
         write(*,*) 'flg_bornonly = ',flg_bornonly
         write(*,*) ''
      endif
c$$$CXXXX DEBUGGING XXXXC      

      call i_setlocalscales0(
     $     flg_bornonly,isReal,nlegborn,
     $     cmp,flav,raisingscales,st_alpha,st_muren2,st_bornorder,
     $     rad_charmthr2,rad_bottomthr2,st_nlight,st_lambda5MSB,
     $     st_renfact,st_facfact,
     $     st_mufact2,
     $     basicfac,bornfac,nlofac)
c$$$     $     indie_ub,pout)

      if(imode.eq.1) then
         if(flg_bornonly.eq.1) then
            rescfac=basicfac
         else
C - Include also Sudakov subtractions (for terms already included in the
C - NLO correction) and explicit scale logs included in the virtual.
            rescfac=basicfac*bornfac
         endif
      elseif(imode.eq.2) then
         rescfac=basicfac*nlofac
      endif

      end

C ---------------------------------------------------------- C
C - Inputs:                                                - C
C - *******                                                - C
C - nlegborn  - No. external legs of Born proc             - C
C - st_nlight - The default number of light flavours under - C
C -             consideration i.e. max no. light flavours  - C
C -             (should be 5, unless you know what you're  - C
C -              doing, really).                           - C
C -                                                        - C
C ---------------------------------------------------------- C
      subroutine i_setlocalscales0
     1    (flg_bornonly,isReal,nlegborn,
     $     pin,flav,raisingscales,st_alpha,st_muren2,st_bornorder,
     $     rad_charmthr2,rad_bottomthr2,st_nlight,st_lambda5MSB,
     $     st_renfact,st_facfact,
     $     st_mufact2,
     $     basicfac,bornfac,nlofac)
c$$$     $     indie_ub,pout)

      implicit none
      integer  nlegborn
      integer  nleg
      integer  ixx
      integer  flav(20),indie_ub(20)
      real * 8 pout(0:3,20)
      real * 8 pin(0:3,20),basicfac,bornfac,nlofac
      integer  onem
      parameter (onem=1000000)
      real * 8 scales(20),p(0:3,20),ptot(0:3),
     1         lscalej,lscalek
      integer  i,j,k,lflav(20),jmerge,kmerge,inlofac
      integer  mergedfl
      real * 8 q2merge,q2merge0,renfac2,facfact2,alphas,mu2,muf2
      real * 8 i_sudakov,i_expsudakov,i_pwhg_alphas,b0,powheginput
      external i_sudakov,i_expsudakov,i_pwhg_alphas,powheginput
      real * 8 q2mergeMAX
      integer  raisingscales
      logical  ini
      save ini
      data ini/.true./

C - maxprocborn and maxprocreal must be greater or equal to the number of
C - independent flavour structures for the born and real process.
      integer maxprocborn
      parameter(maxprocborn=999)
      logical debuggingInfo
      real * 8 st_alpha
      integer  st_bornorder
      real * 8 st_lambda5MSB
      integer  st_nlight
      real * 8 rad_charmthr2,rad_bottomthr2
      integer  flg_bornonly
      real * 8 pi
      parameter (pi=3.141592653589793d0)
      real * 8 st_muren2
      real * 8 st_mufact2
      real * 8 st_renfact
      real * 8 st_facfact

      real * 8 cmpin(0:3,20)
      real * 8 cmpout(0:3,20)
      integer  isReal

CXXXX - Following variables and includes only needed if debuggingInfo = .true.
      integer  flst_born(20,maxprocborn)
      real * 8 tmp_Sudakov,tmp_alphaSratio
      real * 8 tmp_bornfac
CXXXX - Variables and includes above only needed if debuggingInfo = .true.

      if(isReal.eq.1) then
         nleg=nlegborn+1
      else
         nleg=nlegborn
      endif
      do ixx=1+nleg,20
         flav(ixx)    =1000000
         pin(:,ixx)   =0d0
      enddo
      do ixx=1,20
         indie_ub(ixx)=1000000
         lflav(ixx)   =1000000
         pout(:,ixx)  =0d0
         p(:,ixx)     =0d0
         cmpout(:,ixx)=0d0
         cmpin(:,ixx) =0d0
         scales(ixx)  =0d0
      enddo
      ptot(:)=0d0

C ----------------------------- C
C - To debug or not to debug: - C 
      debuggingInfo=.false.
C ----------------------------- C

c$$$CXXXX DEBUGGING XXXXC
      if(debuggingInfo) then
         write(*,*) ''
         write(*,*) ''
         write(*,*) ''
         write(*,*) ''
         write(*,*) 'New born event (indy setlocalscales0)'
         write(*,*) 'nleg ',nleg
         write(*,*) 'st_alpha ',st_alpha
         write(*,*) 'st_bornorder ',st_bornorder
         write(*,*) 'st_muren2 ',st_muren2
         write(*,*) 'rad_charmthr2,rad_bottomthr2 ',
     $               rad_charmthr2,rad_bottomthr2
         write(*,*) 'st_nlight,st_lambda5MSB ',st_nlight,st_lambda5MSB 
         write(*,*) 'st_renfact,st_facfact,raisingscales ',
     $               st_renfact,st_facfact,raisingscales
         write(*,*) 'flg_bornonly = ',flg_bornonly
         write(*,*) ''
      endif
c$$$CXXXX DEBUGGING XXXXC      

      if(debuggingInfo) then
         tmp_bornfac     =-999d0
         tmp_Sudakov     =-999d0
         tmp_alphaSratio =-999d0
      endif
      renfac2=st_renfact**2
      facfact2=st_facfact**2
      lflav=flav
      p=pin
      scales=0
      q2mergeMAX=-1d10
      if(isReal.eq.1) then
         call i_findNearestNeighbours(p,lflav,nleg,st_nlight,
     $                                jmerge,kmerge,mergedfl,q2merge)
c$$$         if(q2merge.lt.1d10) then
c$$$            write(*,*) ''
c$$$            write(*,*) 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
c$$$            write(*,*) 'findNearest neighbours merged: ',
c$$$     $           lflav(jmerge),'(',jmerge,')',
c$$$     $           lflav(kmerge),'(',kmerge,')','--->',mergedfl
c$$$         endif
         if(q2merge.lt.1d10) then
            if(jmerge.gt.2) then
               cmpin=p
               call i_findubfsr(nleg,jmerge,kmerge,cmpin,cmpout)
               p=cmpout
            else
               cmpin=p
               call i_findubisr(nleg,kmerge,cmpin,cmpout)
               p=cmpout
            endif
            lflav(kmerge)=onem
            lflav(jmerge)=mergedfl
         endif
         indie_ub=lflav
         pout=p
         if(q2merge.ge.1d10) then
            goto 99
         endif
      endif
      do ixx=1,nleg
         call i_findNearestNeighbours(p,lflav,nleg,st_nlight,
     $                                jmerge,kmerge,mergedfl,q2merge)
         if(q2merge.lt.1d10) then
c     perform the merging
            if(q2merge.gt.q2mergeMAX) q2mergeMAX=q2merge
            lscalej=scales(jmerge)
            lscalek=scales(kmerge)
            scales(jmerge)=q2merge
            if(lscalej.eq.0) then
c     This is the first merge; it sets the low scale for
c     all partons; no Sudakov factor or reweighting is introduced
               do j=1,nleg
                  scales(j)=q2merge
               enddo
c save this scale; it is the Q_0 scale that appears in all Sudakovs
               q2merge0=q2merge
               bornfac=0
c     Provide alpha_S reweighting for the first merge
               alphas=i_pwhg_alphas(max(q2merge*renfac2,1d0),
     1                              st_lambda5MSB,st_nlight,
     $                              rad_charmthr2,rad_bottomthr2)
               basicfac=alphas/st_alpha
               if(debuggingInfo) then
                  tmp_alphaSratio=alphas/st_alpha
               endif
               nlofac=basicfac
               mu2=max(q2merge*renfac2,1d0)
               muf2=max(q2merge*facfact2,1d0)
               inlofac=1
            else
c provide Sudakov
               basicfac=basicfac*
     1              i_sudakov(q2merge0,q2merge,lscalej,lflav(jmerge),
     $                        st_lambda5MSB,st_nlight,
     $                        rad_charmthr2,rad_bottomthr2)
               basicfac=basicfac*
     1              i_sudakov(q2merge0,q2merge,lscalek,lflav(kmerge),
     $                        st_lambda5MSB,st_nlight,
     $                        rad_charmthr2,rad_bottomthr2)
               bornfac=bornfac+
     1              i_expsudakov(q2merge0,q2merge,lscalej,lflav(jmerge),
     $                           st_lambda5MSB,st_nlight,
     $                           flg_bornonly)
               bornfac=bornfac+
     1              i_expsudakov(q2merge0,q2merge,lscalek,lflav(kmerge),
     $                           st_lambda5MSB,st_nlight,
     $                           flg_bornonly)
c provide alpha_S reweighting
               alphas=i_pwhg_alphas(max(q2merge*renfac2,1d0),
     1                              st_lambda5MSB,st_nlight,
     $                              rad_charmthr2,rad_bottomthr2)
               basicfac=basicfac*alphas/st_alpha
               mu2=mu2*max(q2merge*renfac2,1d0)
               nlofac=nlofac+alphas/st_alpha
               inlofac=inlofac+1
            endif
            if(jmerge.gt.2) then
               p(:,jmerge)=p(:,jmerge)+p(:,kmerge)
            else
               p(3,jmerge)=p(3,jmerge)-p(3,kmerge)
               p(0,jmerge)=abs(p(3,jmerge))
            endif
            lflav(kmerge)=onem
            lflav(jmerge)=mergedfl
         else
            goto 99
         endif
      enddo
 99   continue
c     No more merging is possible. Define the initial scale as
c     the invariant mass of the remaining system
      ptot=0
      do j=3,nleg
         if(lflav(j).ne.onem) then
            ptot=ptot+p(:,j)
         endif
      enddo
      if(raisingscales.eq.1) then
        q2merge=max(q2mergeMAX,
     $              ptot(0)**2-ptot(1)**2-ptot(2)**2-ptot(3)**2)
      else
        q2merge=ptot(0)**2-ptot(1)**2-ptot(2)**2-ptot(3)**2
      endif
      if(scales(1).gt.0) then
         do j=1,nleg
            if(abs(lflav(j)).le.st_nlight) then
               basicfac=basicfac*
     1              i_sudakov(q2merge0,q2merge,scales(j),lflav(j),
     $                        st_lambda5MSB,st_nlight,
     $                        rad_charmthr2,rad_bottomthr2)
               bornfac=bornfac+
     1              i_expsudakov(q2merge0,q2merge,scales(j),lflav(j),
     $                           st_lambda5MSB,st_nlight,
     $                           flg_bornonly)
            endif
         enddo
      else
c If scales(1)=0 no merge has taken place: no sudakovs.
         mu2=1
         muf2=max(q2merge*facfact2,1d0)
         inlofac=0
         bornfac=0
         basicfac=1
         nlofac=0
      endif
      if(debuggingInfo) then
         if(scales(1).gt.0) then
            tmp_Sudakov=1d0
            do j=1,nleg
               if(abs(lflav(j)).le.st_nlight) then
                  tmp_Sudakov=tmp_Sudakov*
     1                 i_sudakov(q2merge0,q2merge,scales(j),lflav(j),
     $                           st_lambda5MSB,st_nlight,
     $                           rad_charmthr2,rad_bottomthr2)
               endif
            enddo
         else
            tmp_Sudakov=1d0
         endif
         tmp_bornfac=bornfac
      endif
      if(debuggingInfo) then
         write(*,*) ''
         write(*,*) ''
         write(*,*) 'st_bornorder = ',st_bornorder
         write(*,*) 'inlofac      = ',inlofac
         write(*,*) ''
         write(*,*) ''
      endif
      if(st_bornorder.gt.inlofac) then
         alphas=i_pwhg_alphas(max(q2merge*renfac2,1d0),
     1                        st_lambda5MSB,st_nlight,
     $                        rad_charmthr2,rad_bottomthr2)

c         do j=inlofac+1,st_bornorder
c            mu2=mu2*max(q2merge*renfac2,1d0)
c            nlofac=nlofac+alphas/st_alpha
c            basicfac=basicfac*alphas/st_alpha
c         enddo
         mu2=mu2*max(q2merge*renfac2,1d0)**(st_bornorder-inlofac)
         nlofac=nlofac+alphas/st_alpha*(st_bornorder-inlofac)
         basicfac=basicfac*(alphas/st_alpha)**(st_bornorder-inlofac)
         inlofac=st_bornorder
      endif 
      nlofac=nlofac/inlofac
      mu2=mu2**(1d0/inlofac)
      b0=(33-2*st_nlight)/(12*pi)
      bornfac=1+st_alpha*nlofac*
     1     (bornfac+st_bornorder*b0*log(mu2/st_muren2))
      st_mufact2=muf2

c$$$CXXXX DEBUGGING XXXXC
      if(debuggingInfo) then
         write(*,*) ''
         write(*,*) ''
         write(*,*) 'Back in setlocalscales ...'
         write(*,*) ''
         write(*,*) 'Flavours: ',(flav(i),i=1,10)
         write(*,*) ''
         write(*,*) 'Born kinematics'
         write(*,*) 'Particle 1      : ',(pin(i,1),i=0,3)
         write(*,*) 'Particle 2      : ',(pin(i,2),i=0,3)
         write(*,*) 'Particle 3      : ',(pin(i,3),i=0,3)
         write(*,*) 'Particle 4      : ',(pin(i,4),i=0,3)
         write(*,*) 'Particles 3+4   : ',
     $        (pin(i,3)+pin(i,4),i=0,3)
         write(*,*) 'Particle 5      : ',(pin(i,5),i=0,3)
         write(*,*) 'Particle 6      : ',(pin(i,6),i=0,3)
         if(isReal.eq.1) write(*,*)'Particle 7      : ',(pin(i,7),i=0,3)
         write(*,*) 'Particle 3+4 pT : ',
     $        sqrt((pin(1,3)+pin(1,4))**2
     $        +(pin(2,3)+pin(2,4))**2)
         write(*,*) 'Particle 5 pT   : ',
     $        sqrt(pin(1,5)**2+pin(2,5)**2)
         write(*,*) 'Particle 6 pT   : ',
     $        sqrt(pin(1,6)**2+pin(2,6)**2)
         write(*,*) ''
         write(*,*) 'Post jet-finder kinematics'
         write(*,*) 'Particle 1      : ',(p(i,1),i=0,3)
         write(*,*) 'Particle 2      : ',(p(i,2),i=0,3)
         write(*,*) 'Particle 3      : ',(p(i,3),i=0,3)
         write(*,*) 'Particle 4      : ',(p(i,4),i=0,3)
         write(*,*) 'Particles 3+4   : ',(p(i,3)+p(i,4),i=0,3)
         write(*,*) 'Particle 5      : ',(p(i,5),i=0,3)
         write(*,*) 'Particle 6      : ',(p(i,6),i=0,3)
         if(isReal.eq.1) write(*,*)'Particle 7      : ',(p(i,7),i=0,3)
         write(*,*) 'Particle 3+4 pT : ',
     $        sqrt((p(1,3)+p(1,4))**2+(p(2,3)+p(2,4))**2)
         write(*,*) 'Particle 5 pT   : ',
     $        sqrt(p(1,5)**2+p(2,5)**2)
         write(*,*) 'Particle 6 pT   : ',
     $        sqrt(p(1,6)**2+p(2,6)**2)
         write(*,*) ''
         write(*,*) 'st_alpha        = ',st_alpha
         write(*,*) 'sqrt(st_muren2) = ',sqrt(st_muren2)
         write(*,*) 'sqrt(st_mufact2)= ',sqrt(st_mufact2)
         write(*,*) 'st_bornorder    = ',st_bornorder
         write(*,*) 'sqrt(q2merge)   = ',sqrt(q2merge)
         write(*,*) 'sqrt(scales(1)) = ',sqrt(scales(1))
         write(*,*) 'sqrt(mu2)       = ',sqrt(mu2)
         write(*,*) 'basicfac        = ',basicfac
         write(*,*) 'Sudakov         = ',tmp_Sudakov
         write(*,*) 'alphaS new      = ',tmp_alphaSratio*st_alpha
         write(*,*) '-999d0*st_alpha = ',-999d0*st_alpha
         write(*,*) 'Sudakov subtr w/o aS = ',tmp_bornfac
         write(*,*) 'full bornfac    = ',bornfac
         write(*,*) 'nlofac          = ',nlofac
         write(*,*) 'inlofac         = ',inlofac
         write(*,*) 'st_lambda5MSB   = ',st_lambda5MSB
         write(*,*) 'st_nlight       = ',st_nlight
         write(*,*) 'rad_charmthr2   = ',rad_charmthr2
         write(*,*) 'rad_bottomthr2  = ',rad_bottomthr2
         write(*,*) 'flg_bornonly    = ',flg_bornonly
      endif
c$$$CXXXX DEBUGGING XXXXC      
      end

C ------------------------------------------------ C
C - Inputs:                                      - C
C - *******                                      - C
C -                                              - C
C - p        - Underlying born momenta           - C
C -                                              - C
C - lflav    - Flavour list derived from         - C
C -            flst_born by subjecting it to     - C
C -            repeated QCD clusterings.         - C
C -                                              - C
C - nleg     - No. external legs                 - C
C -                                              - C
C - st_nlight - The default number of light      - C 
C -             flavours under consideration     - C
C -              i.e. max no. light flavours     - C
C -             (should be 5, unless you know    - C
C -              what you're doing, really).     - C
C -                                              - C
C -                                              - C
C - Outputs:                                     - C
C - ********                                     - C
C -                                              - C
C - jmerge   - Index in lflav of one of the two  - C
C -            closest partons.                  - C
C -                                              - C
C - kmerge   - Index in lflav of the             - C
C -            corresponding parton.             - C
C -                                              - C
C - mergedfl - Flavour of parton resulting from  - C 
C -            combination.                      - C
C -                                              - C
C - q2merge  - pT^2 scale associated to the      - C
C -            merging of jmerge and kmerge.     - C
C -                                              - C
C - checked 24/03/12                             - C
C ------------------------------------------------ C
      subroutine i_findNearestNeighbours(p,lflav,nleg,st_nlight,
     $                                   jmerge,kmerge,mergedfl,q2merge)
      implicit none
C - Input / output:
      real * 8 p(0:3,20)
      integer  lflav(20)
      integer  jmerge,kmerge,mergedfl
      real * 8 q2merge
C - Local variables:
      real * 8 ycm
      integer  onem
      parameter (onem=1000000)
      integer  i,j,k
      integer  fl1,fl2,fl
      real * 8 yj,phij,q2j
      real * 8 yk,phik,q2k
      real * 8 dphi
      real * 8 q2

      integer nleg
      integer st_nlight
      real * 8 pi
      parameter (pi=3.141592653589793d0)

      logical debuggingInfo
C ----------------------------- C
C - To debug or not to debug: - C 
      debuggingInfo=.false.
C ----------------------------- C

c$$$CXXXX DEBUGGING XXXXC
      if(debuggingInfo) then
         write(*,*) 'jet-finder inputs: '
         write(*,*) 'jet-finder flavours: ',(lflav(i),i=1,10)
         write(*,*) 'initial jet-finder kinematics'
         write(*,*) 'Particle 1      : ',(p(i,1),i=0,3)
         write(*,*) 'Particle 2      : ',(p(i,2),i=0,3)
         write(*,*) 'Particle 3      : ',(p(i,3),i=0,3)
         write(*,*) 'Particle 4      : ',(p(i,4),i=0,3)
         write(*,*) 'Particles 3+4   : ',(p(i,3)+p(i,4),i=0,3)
         write(*,*) 'Particle 5      : ',(p(i,5),i=0,3)
         write(*,*) 'Particle 6      : ',(p(i,6),i=0,3)
         write(*,*) 'Particle 3+4 pT : ',
     $        sqrt((p(1,3)+p(1,4))**2+(p(2,3)+p(2,4))**2)
         write(*,*) 'Particle 5 pT   : ',
     $        sqrt(p(1,5)**2+p(2,5)**2)
         write(*,*) 'Particle 6 pT   : ',
     $        sqrt(p(1,6)**2+p(2,6)**2)
         write(*,*) ''
         write(*,*) 'nleg   : ',nleg
         write(*,*) 'st_nlight  : ',st_nlight
         write(*,*) ''
      endif
c$$$CXXXX DEBUGGING XXXXC      

      q2merge=1d10
      ycm=log(p(0,1)/p(0,2))/2
      mergedfl=onem
      do j=3,nleg
         if(abs(lflav(j)).gt.st_nlight) goto 11
         yj=0.5d0*log((p(0,j)+p(3,j))/(p(0,j)-p(3,j)))
         if(yj.gt.ycm) then
            call i_validmergeisr(nleg,st_nlight,lflav,1,j,fl1)
            if(fl1.ne.onem) then
               q2j = p(1,j)**2+p(2,j)**2
               if(q2j.lt.q2merge) then
                  q2merge=q2j
                  jmerge=1
                  kmerge=j
                  mergedfl=fl1
               endif
            endif
         else
            call i_validmergeisr(nleg,st_nlight,lflav,2,j,fl2)
            if(fl2.ne.onem) then
               q2j = p(1,j)**2+p(2,j)**2
               if(q2j.lt.q2merge) then
                  q2merge=q2j
                  jmerge=2
                  kmerge=j
                  mergedfl=fl2
               endif
            endif
         endif
         do k=j+1,nleg
            if(abs(lflav(k)).gt.st_nlight) goto 12
            call i_validmergefsr(nleg,st_nlight,lflav,j,k,fl)
            if(fl.ne.onem) then
               yk=0.5d0*log((p(0,k)+p(3,k))/(p(0,k)-p(3,k)))
               call i_phipt2(p(:,k),phik,q2k)
               call i_phipt2(p(:,j),phij,q2j)
               dphi=abs(phik-phij)
               if(dphi.gt.2*pi) dphi=dphi-2*pi
               if(dphi.gt.pi) dphi=2*pi-dphi
               q2=((yk-yj)**2+dphi**2)*min(q2k,q2j)
               if(q2.lt.q2merge) then
                  q2merge=q2
                  jmerge=j
                  kmerge=k
                  mergedfl=fl
               endif
            endif
 12         continue
         enddo
 11      continue
      enddo

c$$$CXXXX DEBUGGING XXXXC
      if(debuggingInfo) then
         write(*,*) 'jet-finder output: '
         write(*,*) 'jmerge        : ',jmerge
         write(*,*) 'kmerge        : ',kmerge
         write(*,*) 'mergedfl      : ',mergedfl
         write(*,*) 'sqrt(q2merge) : ',sqrt(q2merge)
      endif
c$$$CXXXX DEBUGGING XXXXC
      end

C ------------------------------------------------ C
C - Inputs:                                      - C
C - *******                                      - C
C -                                              - C
C - q20  - q0 scale                              - C
C -                                              - C
C - q2h  - Upper node scale / bound on Sudakov   - C
C -                                              - C
C - q2l  - Lower node scale / bound on Sudakov   - C
C -                                              - C
C - flav - flavour index for the evolving parton - C
C -                                              - C
C - st_lambda5MSB - Lambda 5 in the MS bar       - C
C -                 scheme (by default we are    - C
C -                 getting this directly        - C
C -                 from LHAPDF.                 - C
C -                                              - C
C - st_nlight     - The default number of light  - C 
C -                 flavours under consideration - C
C -                 i.e. max no. light flavours  - C
C -                (should be 5, unless you know - C
C -                 what you're doing, really).  - C
C -                                              - C
C - rad_charmthr2 - The mass of the charm quark  - C 
C -                 (threshold). Should be 1.5**2- C
C -                                              - C
C - rad_bottomthr2 - The mass of the bottom quark- C 
C -                 (threshold). Should be 5.0**2- C
C -                                              - C
C - Outputs:                                     - C
C - ********                                     - C
C - sudakov - The Sudakov form factor.           - C
C -                                              - C
C ------------------------------------------------ C
      function i_sudakov(q20,q2h,q2l,flav,
     $                   st_lambda5MSB,st_nlight,
     $                   rad_charmthr2,rad_bottomthr2)
      implicit none
      real * 8 i_sudakov,q2h,q2l,q20
      integer  flav
      real * 8 lam2
      logical  isQuark
      real * 8 theExponentN,theExponentD
      logical  ini
      data     ini/.true./
      save     ini

      real * 8 st_lambda5MSB
      integer  st_nlight
      real * 8 rad_charmthr2,rad_bottomthr2

      if(ini) then
         call i_sudakov_plotter(st_lambda5MSB,st_nlight,
     $                          rad_charmthr2,rad_bottomthr2)
         ini=.false.
      endif
      lam2=st_lambda5MSB**2
      if(q20.le.lam2.or.q2l.lt.lam2.or.q2h.lt.lam2) then
         i_sudakov=0
         goto 999
      endif
      if(q2l.ge.q2h.or.q2h.le.q20) then
         i_sudakov=1
         goto 999
      endif
      if(flav.eq.0) then
         isQuark=.false.
      else
         isQuark=.true.
      endif
      if(q2l.le.q20) then
        call i_sudakov_exponent(q20,q2h,q2h,theExponentN,
     $                          isQuark,2,.true.,
     $                          st_lambda5MSB,st_nlight,
     $                          rad_charmthr2,rad_bottomthr2)
        i_sudakov=exp(theExponentN)
      else
         call i_sudakov_exponent(q20,q2h,q2h,theExponentN,
     $                           isQuark,2,.true.,
     $                           st_lambda5MSB,st_nlight,
     $                           rad_charmthr2,rad_bottomthr2)
         call i_sudakov_exponent(q20,q2l,q2l,theExponentD,
     $                           isQuark,2,.true.,
     $                           st_lambda5MSB,st_nlight,
     $                           rad_charmthr2,rad_bottomthr2)
         i_sudakov=exp(theExponentN-theExponentD)
      endif
 999  continue
      end


C ------------------------------------------------ C
C - Inputs:                                      - C
C - *******                                      - C
C -                                              - C
C - q20  - q0 scale                              - C
C -                                              - C
C - q2h  - Upper node scale / bound on Sudakov   - C
C -                                              - C
C - q2l  - Lower node scale / bound on Sudakov   - C
C -                                              - C
C - flav - flavour index for the evolving parton - C
C -                                              - C
C - st_lambda5MSB - Lambda 5 in the MS bar       - C
C -                 scheme (by default we are    - C
C -                 getting this directly        - C
C -                 from LHAPDF.                 - C
C -                                              - C
C - st_nlight     - The default number of light  - C 
C -                 flavours under consideration - C
C -                 i.e. max no. light flavours  - C
C -                (should be 5, unless you know - C
C -                 what you're doing, really).  - C
C -                                              - C
C - Outputs:                                     - C
C - ********                                     - C
C -                                              - C
C - LO_sudakov    - The LO Sudakov form factor.  - C
C -                 Uses 1-loop alpha, no CMW    - C
C -                 scheme etc.                  - C
C -                                              - C
C ------------------------------------------------ C
      function i_LO_sudakov(q20,q2h,q2l,flav,
     $                      st_lambda5MSB,st_nlight)
      implicit none
      real * 8 i_LO_sudakov,q2h,q2l,q20
      integer  flav
      real * 8 b0,c,b,lam2

      real * 8 st_lambda5MSB
      integer  st_nlight

      lam2=st_lambda5MSB**2
      if(q20.le.lam2.or.q2l.lt.lam2.or.q2h.lt.lam2) then
         i_LO_sudakov=0
         goto 999
      endif
      if(q2l.ge.q2h.or.q2h.le.q20) then
         i_LO_sudakov=1
         goto 999
      endif
      b0=(33-2*st_nlight)/12d0
      if(flav.eq.0) then
         c=3
         b=b0/3
      else
         c=4d0/3
         b=3d0/4
      endif
      if(q2l.le.q20) then
         i_LO_sudakov= exp(
     1        -c/b0*(  log(log(q2h/lam2)/log(q20/lam2))
     2        *(0.5d0*log(q2h/lam2)-b)
     3         -0.5d0*log(q2h/q20)
     4         ))
      else
         i_LO_sudakov= exp(
     1        -c/b0*(  log(log(q2h/lam2)/log(q20/lam2))
     2        *(0.5d0*log(q2h/lam2)-b)
     3         -0.5d0*log(q2h/q20)
     4         )
     5        +c/b0*(  log(log(q2l/lam2)/log(q20/lam2))
     6        *(0.5d0*log(q2l/lam2)-b)
     7         -0.5d0*log(q2l/q20)
     8         ))
      endif
 999  continue
      end

C ------------------------------------------------ C
C - Inputs:                                      - C
C - *******                                      - C
C -                                              - C
C - q20  - q0 scale                              - C
C -                                              - C
C - q2h  - Upper node scale / bound on Sudakov   - C
C -                                              - C
C - q2l  - Lower node scale / bound on Sudakov   - C
C -                                              - C
C - flav - flavour index for the evolving parton - C
C -                                              - C
C - st_lambda5MSB - Lambda 5 in the MS bar       - C
C -                 scheme (by default we are    - C
C -                 getting this directly        - C
C -                 from LHAPDF.                 - C
C -                                              - C
C - st_nlight     - The default number of light  - C 
C -                 flavours under consideration - C
C -                 i.e. max no. light flavours  - C
C -                (should be 5, unless you know - C
C -                 what you're doing, really).  - C
C -                                              - C
C - flg_bornonly  - 1 if we are just consider    - C
C -                 -ing an LO computation.      - C
C -                                              - C
C - Outputs:                                     - C
C - ********                                     - C
C -                                              - C
C - expsudakov - The Sudakov form factor's expon - C
C -              -ent MODULO a factor of minus   - C
C -              alphaS, integrated with alphaS  - C
C -              fixed. Summed over with the     - C
C -              relevant alphaS factors this is - C
C -              used in compensating the NLO    - C
C -              correction induced when the     - C
C -              Sudakov multiplies the Born.    - C
C -                                              - C
C ------------------------------------------------ C
      function i_expsudakov(q20,q2h,q2l,flav,
     $                      st_lambda5MSB,st_nlight,flg_bornonly)
      implicit none
      real * 8 i_expsudakov,q2h,q2l,q20
      integer  flav
      real * 8 b0,c,b,lam2

      real * 8 st_lambda5MSB
      integer  st_nlight
      integer  flg_bornonly
      real * 8 pi
      parameter (pi=3.141592653589793d0)

C ----------------------------- C
C - To debug or not to debug: - C 
      logical debuggingInfo
      debuggingInfo=.false.
C ----------------------------- C

c$$$CXXXX DEBUGGING XXXXC
      if(debuggingInfo) then
         write(*,*) ''
         write(*,*) 'i_expsudakov debugging'
         write(*,*) 'sqrt(q20)  = ',sqrt(q20)
         write(*,*) 'sqrt(lam2) = ',sqrt(lam2)
         write(*,*) 'sqrt(q2l)  = ',sqrt(q2l)
         write(*,*) 'sqrt(q2h)  = ',sqrt(q2h)
         write(*,*) 'flg_bornonly = ',flg_bornonly
         write(*,*) ''
      endif
c$$$CXXXX DEBUGGING XXXXC      

      lam2=st_lambda5MSB**2
      if(q20.le.lam2.or.q2l.lt.lam2.or.q2h.lt.lam2) then
c in this case everything is zero, irrelevant
         i_expsudakov=0
         return
      endif
      if(q2l.ge.q2h.or.q2h.le.q20.or.(flg_bornonly.eq.1)) then
         i_expsudakov=0
         return
      endif
      b0=(33-2*st_nlight)/12d0
      if(flav.eq.0) then
         c=3
         b=b0/3
      else
         c=4d0/3
         b=3d0/4
      endif
      if(q2l.le.q20) then
         i_expsudakov=
     1        c/pi*(0.25d0*log(q2h/q20)**2 - b*log(q2h/q20))
      else
         i_expsudakov=
     1        c/pi*(0.25d0*log(q2h/q20)**2 - b*log(q2h/q20))
     2      - c/pi*(0.25d0*log(q2l/q20)**2 - b*log(q2l/q20))
      endif
      end

C ---------------------------------------------- C
C - Inputs:                                    - C
C - *******                                    - C
C - nlegborn  - No. external legs of Born proc - C
C -                                            - C
C - st_nlight - The default number of light    - C 
C -             flavours under consideration   - C
C -              i.e. max no. light flavours   - C
C -             (should be 5, unless you know  - C
C -              what you're doing, really).   - C
C -                                            - C
C - flav - flavour list derived from flst_born - C
C -        by subjecting it to repeated QCD    - C
C -        compatible clusterings.             - C
C -                                            - C
C -  i   - index of i-th final-state partICLE  - C
C -        in flav: hence i = 1 or 2 only.     - C
C -                                            - C
C -  j   - index of j-th final-state partICLE  - C
C -        particle in flav.                   - C
C -                                            - C
C - Outputs:                                   - C
C - ********                                   - C
C - fl   - Would-be PDG code of spacelike      - C
C -        "mother" parton obtained by merging - C
C -        (~on-shell) incoming parton i with  - C
C -        outgoing particle j:                - C
C -          i -> fl + j                       - C
C -        Note gluons have id=0 in Powheg-Box - C
C -        instead of 21. If the splitting is  - C
C -        not possible in QCD, fl=1000000 ;   - C
C -        this setting signals to the rest of - C
C -        the algorithm that this is not a    - C
C -        candidate pair for combination.     - C
C -                                            - C
C ---------------------------------------------- C
      subroutine i_validmergeisr(nlegborn,st_nlight,flav,i,j,fl)
      implicit none
      integer  onem
      parameter (onem=1000000)
      integer  flav(20),i,j,fl
      integer  lflav(20)
      logical  i_validflav
      external i_validflav

      integer nlegborn
      integer st_nlight

      if(i.gt.2.or.j.le.2) then  ! Remove when development is finished.
         write(*,*) 'validmergeisr: fatal error'
         write(*,*) 'Routine demands an i.s. and f.s. particle'
         write(*,*) 'index for the 2nd and 3rd input values   '
         write(*,*) 'respectively. Quitting.'
         call exit(-1)
      endif
      if(abs(flav(i)).gt.st_nlight.or.abs(flav(j)).gt.st_nlight) then
         fl=onem
         return
      endif
      if(flav(i).eq.flav(j)) then
c g -> g g or q -> g q
         fl=0
         goto 999
      endif
      if(flav(j).eq.0) then
c q -> q g
         fl=flav(i)
         goto 999
      endif
      if(flav(i).eq.0) then
c g -> qbar q
         fl=-flav(j)
         goto 999
      endif
      fl=onem
      return
 999  continue
C - Check that the flavour list that results from the merging
C - is acceptable e.g. check that for HJJ you don't get back to
C - qqbar->H; if you do then set fl to 1000000, as if the 
C - branching were not possible in QCD s.t. it will be neglected
C - as a candidate for clustering.
      lflav=flav
      lflav(j)=onem
      lflav(i)=fl
      if(.not.i_validflav(lflav)) then
         fl=onem
      endif
      end


C ---------------------------------------------- C
C - Inputs:                                    - C
C - *******                                    - C
C -                                            - C
C - nlegborn  - No. external legs of Born proc - C
C -                                            - C
C - st_nlight - The default number of light    - C 
C -             flavours under consideration   - C
C -              i.e. max no. light flavours   - C
C -             (should be 5, unless you know  - C
C -              what you're doing, really).   - C
C -                                            - C
C - flav - flavour list derived from flst_born - C
C -        by subjecting it to repeated QCD    - C
C -        compatible clusterings.             - C
C -                                            - C
C -  i   - index of i-th final-state partICLE  - C
C -        in flav: hence i = 1 or 2 only.     - C
C -                                            - C
C -  j   - index of j-th final-state partICLE  - C
C -        particle in flav.                   - C
C -                                            - C
C - Outputs:                                   - C
C - ********                                   - C
C - fl   - Would-be PDG code of timelike       - C
C -        "mother" parton obtained by merging - C
C -        outgoing particles i and j:         - C
C -          fl -> i + j                       - C
C -        Note gluons have id=0 in Powheg-Box - C
C -        instead of 21. If the splitting is  - C
C -        not possible in QCD, fl=1000000 ;   - C
C -        this setting signals to the rest of - C
C -        the algorithm that this is not a    - C
C -        candidate pair for combination.     - C
C -                                            - C
C ---------------------------------------------- C
      subroutine i_validmergefsr(nlegborn,st_nlight,flav,i,j,fl)
      implicit none
      integer  onem
      parameter (onem=1000000)
      integer  flav(20),i,j,fl
      integer  lflav(20)
      logical  i_validflav
      external i_validflav

      integer nlegborn
      integer st_nlight

      if(i.le.2.or.j.le.2) then  ! Remove when development is finished.
         write(*,*) 'validmergefsr: fatal error'
         write(*,*) 'Routine demands an f.s. and f.s. particle'
         write(*,*) 'index for the 2nd and 3rd input values   '
         write(*,*) 'respectively. Quitting.'
         call exit(-1)
      endif
      if(abs(flav(i)).gt.st_nlight.or.abs(flav(j)).gt.st_nlight) then
         fl=onem
         return
      endif
      if(flav(i).eq.-flav(j)) then
c g -> g g or g -> q qbar
         fl=0
         goto 999
      endif
      if(flav(j).eq.0) then
c q -> q g
         fl=flav(i)
         goto 999
      endif
      if(flav(i).eq.0) then
c q -> g q
         fl=flav(j)
         goto 999
      endif
      fl=onem
      return
 999  continue
C - Check that the flavour list that results from the merging
C - is acceptable e.g. check that for HJJ you don't get back to
C - qqbar->H; if you do then set fl to 1000000, as if the 
C - branching were not possible in QCD s.t. it will be neglected
C - as a candidate for clustering.
      lflav=flav
      lflav(j)=onem
      lflav(i)=fl
      if(.not.i_validflav(lflav)) then
         fl=onem
      endif
      end


C ---------------------------------------------- C
C - Inputs:                                    - C
C - *******                                    - C
C - p    - p(0) = Energy, p(3) = p_Z           - C
C -                                            - C
C - Outputs:                                   - C
C - ********                                   - C
C - y    - Rapidity                            - C
C - phi  - phi                                 - C
C - q2   - pT^2 w.r.t the beam                 - C
C -                                            - C
C ---------------------------------------------- C
      subroutine i_phipt2(p,phi,q2)
      implicit none
      real * 8 p(0:3),phi,q2
      q2=p(1)**2+p(2)**2
      phi=atan2(p(2),p(1))
      end


C ---------------------------------------------- C
C - validflav:                                 - C
C - ==========                                 - C
C - Forbid us clustering back to an impossible - C
C - Born configuration e.g. gg->Z+nothing.     - C
C - Although the validflavisr and validflavfsr - C
C - routines check that clusterings conserve   - C
C - QCD quantum numbers, they do not check if  - C
C - weak quantum numbers are conserved overall - C
C - This routine is intended to stop us bwds   - C
C - clustering to nonsense like gg->Z+nothing. - C
C - Probably it would be better if validflav   - C
C - was made more general and checked that the - C
C - weak quantum numbers in the initial state  - C
C - equalled those in the final state, to      - C
C - handle the more complicated (higher        - C
C - multiplicity processes).                   - C
C ---------------------------------------------- C
      function i_validflav(lflav)
      implicit none
      logical  i_validflav
      integer  lflav(20)
      integer  ixx
      integer  el_charge_in,el_charge_out

      el_charge_in  = 0
      el_charge_out = 0

      do ixx=1,20
         if(abs(lflav(ixx)).ge.1.and.abs(lflav(ixx)).le.6) then
            if(ixx.le.2) then
               if(mod(abs(lflav(ixx)),2).eq.1) then
                  el_charge_in=el_charge_in-sign(1,lflav(ixx))
               else
                  el_charge_in=el_charge_in+sign(2,lflav(ixx))
               endif
            else
               if(mod(abs(lflav(ixx)),2).eq.1) then
                  el_charge_out=el_charge_out-sign(1,lflav(ixx))
               else
                  el_charge_out=el_charge_out+sign(2,lflav(ixx))
               endif
            endif
         elseif(abs(lflav(ixx)).ge.11.and.abs(lflav(ixx)).le.16) then
            if(ixx.le.2) then
               if(mod(abs(lflav(ixx)),2).eq.1) then
                  el_charge_in=el_charge_in-sign(3,lflav(ixx))
               endif
            else
               if(mod(abs(lflav(ixx)),2).eq.1) then
                  el_charge_out=el_charge_out-sign(3,lflav(ixx))
               endif
            endif
         endif
      enddo

      if(el_charge_in.ne.el_charge_out) then
         i_validflav = .false.
      else
         i_validflav = .true.
      endif

      if(lflav(1).eq.0.and.lflav(2).eq.0.and.
     $   lflav( 5).eq.1000000.and.lflav( 6).eq.1000000.and.
     $   lflav( 7).eq.1000000.and.lflav( 8).eq.1000000.and.
     $   lflav( 9).eq.1000000.and.lflav(10).eq.1000000.and.
     $   lflav(11).eq.1000000.and.lflav(12).eq.1000000) then
         i_validflav = .false.
      else
         i_validflav = .true.
      endif

      end


C ********* DDT / Ellis-Veseli / Nason-Ridolfi Sudakov ************ C
C -                                                               - C
C - Output:                                                       - C
C - ========                                                      - C
C - The value of the Sudakov exponent defined as the integral,    - C
C - from  Log [ ql^2/Lambda^2 ]  up to Log [  qh^2/Lambda^2 ], of - C
C -                                                               - C
C -    d Log[ q^2/Lambda^2 ]                                      - C
C -  - {                                                          - C
C -      aSBar*A1*Log[m^2/q^2] + aSBar^2*A2*Log[m^2/q^2]          - C
C -    + aSBar*B1              + aSBar^2*B2                       - C
C -    }                                                          - C
C -                                                               - C
C - where aSBar = aS/2/Pi.                                        - C
C -                                                               - C
C - For m2=qh2, except for an overall factor of two this is the   - C
C - Sudakov form factor of eq. 32 in the Ellis-Veseli paper - in  - C
C - that paper they have two quark lines to consider while here   - C
C - we only want to consider one line at a time. The factor of    - C
C - two is manifest in the code below as our A1, A2, B1, B2       - C
C - coefficients defined to be HALF of the Ellis-Veseli ones.     - C
C -                                                               - C
C - The m2 dependence is a relic of the Nason & Ridolfi form of   - C
C - the Sudakov form factor, which has the numerator in the large - C
C - log equal to mZZ but the upper bound on the Sudakov integral  - C
C - is Q^2. There doesn't seem to be any problem arising when     - C
C - you just call the routine with m^2=Q^2, but maybe if this     - C
C - gets resolved I can re-do the mathematica integral.           - C
C -                                                               - C
C - The analytic integral was done in Mathematica assuming no     - C
C - flavour thresholds. When q^2 is below the b or c quark        - C 
C - flavour thresholds a numerical integration is done instead    - C
C - using h_dgauss. The numerical integration and analytic results- C
C - agree very well above these thresholds - try resetting        - C
C - debuggingEpsilon below.                                       - C
C -                                                               - C 
C - For the default values of A1, A2, B1, B2 in the code below    - C
C - the Sudakov should correspond to that of Nason and Ridolfi,   - C
C - which has an effective B2 term by virtue of the fact that     - C
C - the CMW alpha_S is used to multiply the leading & subleading  - C
C - term. At least with the calculation done in the way it is we  - C
C - can easily play around with the coefficients.                 - C 
C -                                                               - C 
C - To Use:                                                       - C
C - =======                                                       - C
C - q2l = The scale of the lower  clustered node.                 - C
C - q2h = The scale of the higher clustered node.                 - C
C - m2  = The boson mass squared (argument of the                 - C
C -       log in the exponent of N.R. eq 4.8).                    - C
C - theExponent                                                   - C 
C -     = The value of the curly brackets in N.R. eq 4.8.         - C
C - isQuark = .true. for a quark propagator                       - C 
C - theAccuracy = 0 for 1-loop alphaS and A2=B2=0,                - C 
C -             = 1 for 2-loop alphaS and Powheg A & B coeffs     - C 
C -             = 2 for 2-loop alphaS and NLL A & B coeffs        - C 
C -                                                               - C 
C - Notes:                                                        - C
C - ======                                                        - C
C - Details for the integration in the sudakov exponent can be    - C
C - found in the Mathematica notebook: menlops/DDT_exponent.nb .  - C 
C -                                                               - C 
C - The Mathematica notebook shows plots in which the 5-flavour   - C
C - and 4-flavour alphaS differ by <2% at pT=3 GeV and 4% at      - C
C - pT=2 GeV. Using the C.M.W. alphaS (aS -> aS*(1+aS*K/2*pi))    - C
C - increases these differences but they remain small: 2% at pT=4 - C
C - GeV and 6% at pT=2 GeV. Note well that since the program      - C
C - matches alphaS at flavour thresholds, not                     - C
C - alphaS*(1+alphaS*K/2*pi), since K too actually depends on     - C
C - the number of flavours, this means the 3,4 and 5 flavour      - C
C - C.M.W. alphaS*(1+alphaS*K/2*pi) DO NOT match at the flavour   - C
C - thresholds in pT! Whereas alphaS 4 and 5 flavour couplings    - C
C - match at 5 GeV, the nf dependence of K means that the 4 and 5 - C
C - flavour alphaS*(1+alphaS*K/2*pi) actually meet at about 9 GeV - C
C - instead.                                                      - C
C -                                                               - C
C - Inputs:                                                       - C
C - =======                                                       - C
C - st_lambda5MSB    - Lambda 5 in the MS bar scheme (by default  - C
C -                    we are getting this directly from LHAPDF.  - C
C - st_nlight        - The default number of light flavours under - C
C -                    consideration i.e. max no. light flavours  - C
C -                    (should be 5, unless you know what you're  - C
C -                    doing, really).                            - C
C -                    we are getting this directly from LHAPDF.  - C
C - rad_charmthr2    - The mass of the charm quark (threshold)    - C 
C -                    Should be 1.5**2.                          - C 
C - rad_bottomthr2   - The mass of the bottom quark (threshold)   - C 
C -                    Should be 5.0**2.                          - C 
C -                                                               - C
C ***************************************************************** C
      subroutine i_sudakov_exponent(q2l,q2h,m2,theExponent,isQuark,
     $                              theAccuracy,fixed_nf,
     $                              st_lambda5MSB,st_nlight,
     $                              rad_charmthr2,rad_bottomthr2)
      implicit none

C - The inputs:
      real * 8 q2l,q2h,m2
      logical  isQuark
      integer  theAccuracy
      logical  fixed_nf ! Assume the number of active flavours in 
                        ! is always st_nlight when evaluating the
                        ! Sudakov.
C - The output:
      real * 8 theExponent

C - Auxiliary variables:
      real * 8 st_lambda5MSB
      integer  st_nlight
      real * 8 rad_charmthr2,rad_bottomthr2
      real * 8 CF,CA
      save     CF,CA
      logical  ini
      data     ini/.true./
      save     ini
      real * 8 pi
      parameter (pi=3.141592653589793d0)

C - More auxiliary variables:
      integer  nf
      real * 8 bnf,bpnf,K
      real * 8 A1,B1,A2,B2,zeta3
      real * 8 Lq2l,Lq2h,Lm2
      real * 8 i_pwhg_alphas
      external i_pwhg_alphas
      real * 8 aSbar
      real * 8 theA1coeff,theA2coeff,theB1coeff,theB2coeff
      real * 8 eps
      logical  i_isQuarkLine
      integer  i_accuracy
      real * 8 i_m2_common
      common/i_sudakov_integral/i_isQuarkLine,i_accuracy,i_m2_common
      real * 8 h_dgauss,i_sudakov_exponent_integrand
      external h_dgauss,i_sudakov_exponent_integrand
      logical  i_pwhg_isfinite
      external i_pwhg_isfinite
      real * 8 debuggingEpsilon,tmp

      if(ini) then
         CF=4d0/3d0
         CA=3d0
         ini=.false.
      endif

C - Fractional difference (%) between analytic and numerical Sudakov
C - exponent integration which leads to an error on the screen if
C - it is exceeded. Making it negative (recommended) deactivates
C - this debugging.
      debuggingEpsilon = -999d0 !  1d-9

      i_isQuarkLine  = isQuark

      i_accuracy = theAccuracy

      if(fixed_nf) then
         nf=st_nlight
      else
         if(q2l.lt.rad_charmthr2) then
            nf=3
         elseif(q2l.lt.rad_bottomthr2) then
            nf=4
         else
            nf=5
         endif
      endif

      bnf  = (11d0*CA-2d0*nf)/12/Pi
      bpnf = (153 - 19d0*nf) / Pi / 2 / (33 - 2*nf)

      K  =  (67d0/18-Pi**2/6)*CA-5d0/9*nf
      if(i_isQuarkLine) then
         A1 =   Cf
         A2 =   Cf*K
         B1 =  -3d0/2*Cf
C - Powheg spurious B2:
         B2 =  -3d0/2*Cf*K
      else
         A1 =  CA
         A2 =  CA*K
         B1 = -2*Pi*bnf
C - Powheg spurious B2:
         B2 = -2*Pi*bnf*K
      endif

      Lq2l = Log(q2l/st_lambda5MSB/st_lambda5MSB)
      Lq2h = Log(q2h/st_lambda5MSB/st_lambda5MSB)
      Lm2  = Log(m2/st_lambda5MSB/st_lambda5MSB)

      aSbar = i_pwhg_alphas(q2l,st_lambda5MSB,-1,
     $                      rad_charmthr2,rad_bottomthr2)/2/Pi

      if(i_accuracy.eq.0) then
         A2   = 0d0 ! NLL coefficient
         B2   = 0d0 ! NNLL coefficient
         bpnf = 0d0 ! Reduce 2-loop to 1-loop alpha in the calculation
      else if(i_accuracy.eq.2) then
         B2   = 0d0 ! NNLL coefficient
      else if(i_accuracy.eq.3) then
         B1   = 0d0 ! NLL coefficient
         A2   = 0d0 ! NLL coefficient
         B2   = 0d0 ! NNLL coefficient
         bpnf = 0d0 ! Reduce 2-loop to 1-loop alpha in the calculation
      endif

      if(q2l.ge.rad_bottomthr2.or.fixed_nf) then

         theA1coeff =
     $        ( (Lq2h - Lq2l)
     $           - Lm2*Log(Lq2h/Lq2l)
     $        )/(2*bnf*Pi)
     $      + bpnf*( 2*Lm2*(Lq2h - Lq2l)
     $             + 2*Lm2*Lq2h*Log(Lq2l)
     $             + Lq2l*Lq2h*Log(Lq2l)**2 
     $             - 2*Lm2*Lq2l*Log(Lq2h)
     $             - Lq2l*Lq2h*Log(Lq2h)**2
     $             )/(4*bnf**2*Lq2l*Lq2h*Pi)

         theA2coeff = 
     $         ( Lm2*(Lq2l - Lq2h)
     $         - Lq2l*Lq2h*Log(Lq2l)
     $         + Lq2l*Lq2h*Log(Lq2h)
     $         )/(4d0*bnf**2*Lq2l*Lq2h*Pi**2)
     $       + bpnf*( 0.5*(Lq2h-Lq2l)*( Lm2*(Lq2l + Lq2h)
     $                                - 4*Lq2l*Lq2h
     $                                )
     $              + (Lm2 - 2*Lq2l)*Lq2h**2*Log(Lq2l)
     $              - Lq2l**2*(Lm2 - 2*Lq2h)*Log(Lq2h)
     $              )/(4d0*bnf**3*Lq2l**2*Lq2h**2*Pi**2)
     $       + bpnf**2*(  (Lq2l - Lq2h)*( 8*Lm2*( Lq2l**2
     $                                          + Lq2l*Lq2h
     $                                          + Lq2h**2
     $                                          )
     $                                  - 27*Lq2l*Lq2h*(Lq2l + Lq2h)
     $                                  )
     $                 - 6*(4*Lm2 - 9*Lq2l)*Lq2h**3*Log(Lq2l)
     $                 + 18*(-2*Lm2 + 3*Lq2l)*Lq2h**3*Log(Lq2l)**2
     $                 + 6*Lq2l**3*(4*Lm2 - 9*Lq2h)*Log(Lq2h)
     $                 + 18*Lq2l**3*(2*Lm2 - 3*Lq2h)*Log(Lq2h)**2
     $                 )/(432.*bnf**4*Lq2l**3*Lq2h**3*Pi**2)

         theB1coeff = 
     $       - Log(Lq2h/Lq2l)/(2.*bnf*Pi)
     $       - bpnf*( Lq2l - Lq2h - Lq2h*Log(Lq2l)
     $              + Lq2l*Log(Lq2h)
     $              )/(2.*bnf**2*Lq2l*Lq2h*Pi) 

         theB2coeff =
     $       - (1/Lq2l - 1/Lq2h)/(4.*bnf**2*Pi**2)
     $       - bpnf*(   Lq2l**2 - Lq2h**2 - 2*Lq2h**2*Log(Lq2l)
     $              + 2*Lq2l**2*Log(Lq2h)
     $              )/(8.*bnf**3*Lq2l**2*Lq2h**2*Pi**2)
     $       - bpnf**2*( 2*Lq2h**3 - 2*Lq2l**3
     $                 + 3*Lq2h**3*Log(Lq2l)*(2 + 3*Log(Lq2l)) 
     $                 - 3*Lq2l**3*Log(Lq2h)*(2 + 3*Log(Lq2h))
     $                 )/(108.*bnf**4*Lq2l**3*Lq2h**3*Pi**2)

         theExponent = A1*theA1coeff + A2*theA2coeff
     $               + B1*theB1coeff + B2*theB2coeff

         if(debuggingEpsilon.gt.0d0) then
            tmp=h_dgauss(i_sudakov_exponent_integrand,Lq2l,Lq2h,eps,
     $                   st_lambda5MSB,rad_charmthr2,rad_bottomthr2)
            tmp=100d0*abs((tmp-theExponent)/(tmp+theExponent))
            if(tmp.gt.debuggingEpsilon) then
               write(*,*) ''
               write(*,*) 'sudakov_exponent: debug mode'
               write(*,*) '============================'
               write(*,*) 'sqrt(q2l) (GeV) = ',sqrt(q2l)
               write(*,*) 'sqrt(q2h) (GeV) = ',sqrt(q2h)
               write(*,*) 'sqrt(m2) (GeV)  = ',sqrt(m2)
               write(*,*) 'analytic - numerical exponent (%) = ',tmp
               write(*,*) 'analytic exponent = ',theExponent
               write(*,*) 'A1 term = ',A1*theA1coeff
               write(*,*) 'A2 term = ',A1*theA2coeff
               write(*,*) 'B1 term = ',A1*theB1coeff
               write(*,*) 'B2 term = ',A1*theB2coeff
            endif
         endif
      else
         eps = 1d-6
         i_m2_common=m2
         theExponent
     $       = h_dgauss(i_sudakov_exponent_integrand,Lq2l,Lq2h,eps,
     $                  st_lambda5MSB,rad_charmthr2,rad_bottomthr2)
      endif

      if(.not.i_pwhg_isfinite(theExponent)) then
         write(6,*) ' '
         write(6,*) 'Warning: sudakov_exponent is weird.'
         write(6,*) 'theExponent = ',theExponent
         write(6,*) 'exp(theExponent) = ',exp(theExponent)
         write(6,*) 'q_low   = ',sqrt(q2l)
         write(6,*) 'q_hi    = ',sqrt(q2h)
         write(6,*) 'm       = ',sqrt(m2)
      endif

      end


C ***************************************************************** C
C - The integrand in Sudakov exponent times q^2 i.e. we effectiv  - C
C - ely integrate in log(q2/lambda^2) (Lq2) to better sample      - C
C - the Sudakov peak region. This can be used to replace the      - C
C - above analytic result in the region below the b- and c-quark  - C
C - thresholds when a variable flavour number should come into    - C
C - effect (not taken account of in the analytic result!). We see - C
C - excellent agreement above the b-quark threshold between       - C
C - numerical and analytic results (basically exact agreement),   - C
C - while below pT=5 GeV we see differences of at most 5%, rising - C
C - to 10% and 50% as pT goes below the charm threshold (~2 GeV). - C
C - I suggest using this routine AS AN ENHANCEMENT of the one     - C
C - above, i.e. to call the one here with a low value of EPS for  - C
C - events in which pT < b-quark threshold.                       - C
C -                                                               - C
C - Inputs:                                                       - C
C - =======                                                       - C
C - Lq2              - log(q2/lambda^2) with q2 being the "usual" - C
C -                    integration variable in the Sudakov.       - C
C - st_lambda5MSB    - Lambda 5 in the MS bar scheme (by default  - C
C -                    we are getting this directly from LHAPDF.  - C
C - rad_charmthr2    - The mass of the charm quark (threshold)    - C 
C -                    Should be 1.5**2.                          - C 
C - rad_bottomthr2   - The mass of the bottom quark (threshold)   - C 
C -                    Should be 5.0**2.                          - C 
C ***************************************************************** C
      function i_sudakov_exponent_integrand(Lq2,
     $                                   st_lambda5MSB,
     $                                   rad_charmthr2,rad_bottomthr2)
      implicit none
      real * 8 i_sudakov_exponent_integrand
      real * 8 Lq2
      logical  i_isQuarkLine
      integer  i_accuracy
      real * 8 i_m2_common
      common/i_sudakov_integral/i_isQuarkLine,i_accuracy,i_m2_common
      real * 8 q2,aSbar
      integer  nf
      real * 8 bnf,bpnf,K
      real * 8 A1,B1,A2,B2
      real * 8 i_pwhg_alphas,i_pwhg_alphas0
      external i_pwhg_alphas,i_pwhg_alphas0

      real * 8 st_lambda5MSB
      real * 8 rad_charmthr2,rad_bottomthr2
      real * 8 CF,CA
      save     CF,CA
      logical  ini
      data     ini/.true./
      save     ini
      real * 8 pi
      parameter (pi=3.141592653589793d0)
      if(ini) then
         CF=4d0/3d0
         CA=3d0
         ini=.false.
      endif

      q2=st_lambda5MSB**2*exp(Lq2)

      if(q2.lt.rad_charmthr2) then
         nf=3
      elseif(q2.lt.rad_bottomthr2) then
         nf=4
      else
         nf=5
      endif

      if(i_accuracy.eq.0) then    ! Do 1-loop running AND keep nf=5
         nf = 5
         aSbar = i_pwhg_alphas0(q2,st_lambda5MSB,nf)/2/Pi
      elseif(i_accuracy.eq.3) then! Do 1-loop running AND keep nf=5
         nf = 5
         aSbar = i_pwhg_alphas0(q2,st_lambda5MSB,nf)/2/Pi
      else                     ! Do 2-loop running
         aSbar = i_pwhg_alphas(q2,st_lambda5MSB,-1,
     $                         rad_charmthr2,rad_bottomthr2)/2/Pi
      endif

      bnf  = (11d0*CA-2d0*nf)/12/Pi
      bpnf = (153 - 19d0*nf) / Pi / 2 / (33 - 2*nf)

      K  =  (67d0/18-Pi**2/6)*CA-5d0/9*nf
      if(i_isQuarkLine) then
         A1 =   Cf
         A2 =   Cf*K
         B1 =  -3d0/2*Cf
C - Powheg spurious B2:
         B2 =  -3d0/2*Cf*K
      else
         A1 =  CA
         A2 =  CA*K
         B1 = -2*Pi*bnf
C - Powheg spurious B2:
         B2 = -2*Pi*bnf*K
      endif

      if(i_accuracy.eq.0) then
         A2  = 0d0 ! NLL coefficient
         B2  = 0d0 ! NNLL coefficient
      else if(i_accuracy.eq.2) then
         B2  = 0d0 ! NNLL coefficient
      else if(i_accuracy.eq.3) then
         B1  = 0d0 ! NLL coefficient
         A2  = 0d0 ! NLL coefficient
         B2  = 0d0 ! NNLL coefficient
      endif

      i_sudakov_exponent_integrand=
     $     - ( (A1*aSbar+A2*aSbar*aSbar)*log(i_m2_common/q2)
     $       + (B1*aSbar+B2*aSbar*aSbar)
     $       )

      end


C ***************************************************************** C
C - H_DGAUSS integration routine from cernlib (not sure precisely - C
C - where I got it anymore).                                      - C
C ***************************************************************** C
*
* $Id: gauss64.F,v 1.1.1.1 1996/04/01 15:02:13 mclareni Exp $
*
* $Log: gauss64.F,v $
* Revision 1.1.1.1  1996/04/01 15:02:13  mclareni
* Mathlib gen
*
*
      FUNCTION H_DGAUSS(F,A,B,EPS,
     $                  st_lambda5MSB,rad_charmthr2,rad_bottomthr2)

*
* $Id: imp64.inc,v 1.1.1.1 1996/04/01 15:02:59 mclareni Exp $
*
* $Log: imp64.inc,v $
* Revision 1.1.1.1  1996/04/01 15:02:59  mclareni
* Mathlib gen
*
*
* imp64.inc
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER NAME*(*)
      PARAMETER (NAME = 'H_DGAUSS')
      real * 8 st_lambda5MSB
      real * 8 rad_charmthr2,rad_bottomthr2
*
* $Id: gausscod.inc,v 1.1.1.1 1996/04/01 15:02:13 mclareni Exp $
*
* $Log: gausscod.inc,v $
* Revision 1.1.1.1  1996/04/01 15:02:13  mclareni
* Mathlib gen
*
*
*
* gausscod.inc
*
      DIMENSION W(12),X(12)
      PARAMETER (Z1 = 1, HF = Z1/2, CST = 5*Z1/1000)
      DATA X( 1) /9.6028985649753623D-1/, W( 1) /1.0122853629037626D-1/
      DATA X( 2) /7.9666647741362674D-1/, W( 2) /2.2238103445337447D-1/
      DATA X( 3) /5.2553240991632899D-1/, W( 3) /3.1370664587788729D-1/
      DATA X( 4) /1.8343464249564980D-1/, W( 4) /3.6268378337836198D-1/
      DATA X( 5) /9.8940093499164993D-1/, W( 5) /2.7152459411754095D-2/
      DATA X( 6) /9.4457502307323258D-1/, W( 6) /6.2253523938647893D-2/
      DATA X( 7) /8.6563120238783174D-1/, W( 7) /9.5158511682492785D-2/
      DATA X( 8) /7.5540440835500303D-1/, W( 8) /1.2462897125553387D-1/
      DATA X( 9) /6.1787624440264375D-1/, W( 9) /1.4959598881657673D-1/
      DATA X(10) /4.5801677765722739D-1/, W(10) /1.6915651939500254D-1/
      DATA X(11) /2.8160355077925891D-1/, W(11) /1.8260341504492359D-1/
      DATA X(12) /9.5012509837637440D-2/, W(12) /1.8945061045506850D-1/
      H=0
      IF(B .EQ. A) GO TO 99
      CONST=CST/ABS(B-A)
      BB=A
    1 AA=BB
      BB=B
    2 C1=HF*(BB+AA)
      C2=HF*(BB-AA)
      S8=0
      DO 3 I = 1,4
      U=C2*X(I)
    3 S8=S8+W(I)*(F(C1+U,st_lambda5MSB,rad_charmthr2,rad_bottomthr2)
     $           +F(C1-U,st_lambda5MSB,rad_charmthr2,rad_bottomthr2))
      S16=0
      DO 4 I = 5,12
      U=C2*X(I)
    4 S16=S16+W(I)*(F(C1+U,st_lambda5MSB,rad_charmthr2,rad_bottomthr2)
     $             +F(C1-U,st_lambda5MSB,rad_charmthr2,rad_bottomthr2))
      S16=C2*S16
      IF(ABS(S16-C2*S8) .LE. EPS*(1+ABS(S16))) THEN
       H=H+S16
       IF(BB .NE. B) GO TO 1
      ELSE
       BB=C1
       IF(1+CONST*ABS(C2) .NE. 1) GO TO 2
       H=0
       CALL i_MTLPRT(NAME,'D103.1','TOO HIGH ACCURACY REQUIRED')
       GO TO 99
      END IF
   99 H_DGAUSS=H
      RETURN
      END

C ***************************************************************** C
C - This routine plots the sudakov according to what was in the   - C
C - sudakov exponent routine above.                               - C
C -                                                               - C
C - Inputs:                                                       - C
C - =======                                                       - C
C - st_lambda5MSB    - Lambda 5 in the MS bar scheme (by default  - C
C -                    we are getting this directly from LHAPDF.  - C
C - st_nlight        - The default number of light flavours under - C
C -                    consideration i.e. max no. light flavours  - C
C -                    (should be 5, unless you know what you're  - C
C -                    doing, really).                            - C
C -                    we are getting this directly from LHAPDF.  - C
C - rad_charmthr2    - The mass of the charm quark (threshold)    - C 
C -                    Should be 1.5**2.                          - C 
C - rad_bottomthr2   - The mass of the bottom quark (threshold)   - C 
C -                    Should be 5.0**2.                          - C 
C -                                                               - C
C ***************************************************************** C

      subroutine i_sudakov_plotter(st_lambda5MSB,st_nlight,
     $                             rad_charmthr2,rad_bottomthr2)

      implicit none
      real * 8 the_q,the_Lq2,the_qh
      real * 8 the_Sudakov,one_minus_alpha
      real * 8 the_step,the_start,the_end
      integer  no_steps,ixx,iun
      logical  isQuark

      real * 8 st_lambda5MSB
      integer  st_nlight
      real * 8 rad_charmthr2,rad_bottomthr2

      no_steps=1000

      the_qh=120d0

      call i_newunit(iun)
      open(unit=iun,file='Sudakov.top')

      write(iun,*) 'SET FONT DUPLEX'
      write(iun,*) 'SET WINDOW X 1 12 Y 2.5 9'
      write(iun,*) 'SET INTENSITY 5'
      write(iun,*) 'SET LABELS SIZE 2.5'
      write(iun,*) 'SET TITLE TOP SCALE 1.0'
      write(iun,*) 'SET TITLE SIZE 2.0'
      write(iun,*) 'TITLE TEXT 7.0 6.0 "quark A1,B1"'
      write(iun,*) 'CASE               "           "'
      write(iun,*) '5.5 6.0'
      write(iun,*) '6.5 6.0'
      write(iun,*) 'JOIN TEXT RED'
      write(iun,*) 'TITLE TEXT 7.0 5.5 "quark A1,B1,A2"'
      write(iun,*) 'CASE               "              "'
      write(iun,*) '5.5 5.5'
      write(iun,*) '6.5 5.5'
      write(iun,*) 'JOIN TEXT GREEN'
      write(iun,*) 'TITLE TEXT 7.0 5.0 "quark Powheg"'
      write(iun,*) 'CASE               "            "'
      write(iun,*) '5.5 5.0'
      write(iun,*) '6.5 5.0'
      write(iun,*) 'JOIN TEXT BLUE'
      write(iun,*) 'TITLE TEXT 7.0 4.5 "gluon A1,B1"'
      write(iun,*) 'CASE               "           "'
      write(iun,*) '5.5 4.5'
      write(iun,*) '6.5 4.5'
      write(iun,*) 'JOIN TEXT RED DASHES'
      write(iun,*) 'TITLE TEXT 7.0 4.0 "gluon A1,B1,A2"'
      write(iun,*) 'CASE               "              "'
      write(iun,*) '5.5 4.0'
      write(iun,*) '6.5 4.0'
      write(iun,*) 'JOIN TEXT GREEN DASHES'
      write(iun,*) 'TITLE TEXT 7.0 3.5 "gluon Powheg"'
      write(iun,*) 'CASE               "            "'
      write(iun,*) '5.5 3.5'
      write(iun,*) '6.5 3.5'
      write(iun,*) 'JOIN TEXT GREEN DASHES'
      write(iun,*) 'TITLE BOTTOM "p0T1 (GeV)"'
      write(iun,*) 'CASE         " X X S   S"'
      write(iun,*) 'SET ORDER X Y DY'

C - Plot Sudakov vs kT
      the_start= 0.5d0
      the_end  = 
     $     st_lambda5MSB
     $   * exp(0.5d0*(Log(the_qh**2/st_lambda5MSB**2)-1.5d0))
      the_step = (the_end-the_start)/dble(no_steps)

      isQuark = .true.
      do ixx=1,no_steps
         the_q = the_start + ixx*the_step
         call i_sudakov_exponent(the_q**2,the_qh**2,the_qh**2,
     $                           the_Sudakov,isQuark,0,.true.,
     $                           st_lambda5MSB,st_nlight,
     $                           rad_charmthr2,rad_bottomthr2)

         the_Sudakov=exp(the_Sudakov)
         write(iun,'(F8.5,A,F8.6)') the_q,'  ',the_Sudakov
      enddo
      write(iun,*) 'JOIN RED'
      isQuark = .true.
      do ixx=1,no_steps
         the_q = the_start + ixx*the_step
         call i_sudakov_exponent(the_q**2,the_qh**2,the_qh**2,
     $                           the_Sudakov,isQuark,2,.true.,
     $                           st_lambda5MSB,st_nlight,
     $                           rad_charmthr2,rad_bottomthr2)

         the_Sudakov=exp(the_Sudakov)
         write(iun,'(F8.5,A,F8.6)') the_q,'  ',the_Sudakov
      enddo
      write(iun,*) 'JOIN GREEN'
      isQuark = .true.
      do ixx=1,no_steps
         the_q = the_start + ixx*the_step
         call i_sudakov_exponent(the_q**2,the_qh**2,the_qh**2,
     $                           the_Sudakov,isQuark,1,.true.,
     $                           st_lambda5MSB,st_nlight,
     $                           rad_charmthr2,rad_bottomthr2)

         the_Sudakov=exp(the_Sudakov)
         write(iun,'(F8.5,A,F8.6)') the_q,'  ',the_Sudakov
      enddo
      write(iun,*) 'JOIN BLUE'
      isQuark = .false.
      do ixx=1,no_steps
         the_q = the_start + ixx*the_step
         call i_sudakov_exponent(the_q**2,the_qh**2,the_qh**2,
     $                           the_Sudakov,isQuark,0,.true.,
     $                           st_lambda5MSB,st_nlight,
     $                           rad_charmthr2,rad_bottomthr2)
         the_Sudakov=exp(the_Sudakov)
         write(iun,'(F8.5,A,F8.6)') the_q,'  ',the_Sudakov
      enddo
      write(iun,*) 'JOIN RED DASHES'
      isQuark = .false.
      do ixx=1,no_steps
         the_q = the_start + ixx*the_step
         call i_sudakov_exponent(the_q**2,the_qh**2,the_qh**2,
     $                           the_Sudakov,isQuark,2,.true.,
     $                           st_lambda5MSB,st_nlight,
     $                           rad_charmthr2,rad_bottomthr2)
         the_Sudakov=exp(the_Sudakov)
         write(iun,'(F8.5,A,F8.6)') the_q,'  ',the_Sudakov
      enddo
      write(iun,*) 'JOIN GREEN DASHES'
      isQuark = .false.
      do ixx=1,no_steps
         the_q = the_start + ixx*the_step
         call i_sudakov_exponent(the_q**2,the_qh**2,the_qh**2,
     $                           the_Sudakov,isQuark,1,.true.,
     $                           st_lambda5MSB,st_nlight,
     $                           rad_charmthr2,rad_bottomthr2)
         the_Sudakov=exp(the_Sudakov)
         write(iun,'(F8.5,A,F8.6)') the_q,'  ',the_Sudakov
      enddo
      write(iun,*) 'JOIN BLUE DASHES'

      write(iun,*) ''
      write(iun,*) 'NEW FRAME'
      write(iun,*) 'SET FONT DUPLEX'
      write(iun,*) 'SET WINDOW X 1 12 Y 2.5 9'
      write(iun,*) 'SET INTENSITY 5'
      write(iun,*) 'SET LABELS SIZE 2.5'
      write(iun,*) 'SET TITLE TOP SCALE 1.0'
      write(iun,*) 'SET TITLE SIZE 2.0'
      write(iun,*) 'TITLE TEXT 3.0 8.0 "quark A1,B1"'
      write(iun,*) 'CASE               "           "'
      write(iun,*) '1.5 8.0'
      write(iun,*) '2.5 8.0'
      write(iun,*) 'JOIN TEXT RED'
      write(iun,*) 'TITLE TEXT 3.0 7.5 "quark A1,B1,A2"'
      write(iun,*) 'CASE               "              "'
      write(iun,*) '1.5 7.5'
      write(iun,*) '2.5 7.5'
      write(iun,*) 'JOIN TEXT GREEN'
      write(iun,*) 'TITLE TEXT 3.0 7.0 "quark Powheg"'
      write(iun,*) 'CASE               "            "'
      write(iun,*) '1.5 7.0'
      write(iun,*) '2.5 7.0'
      write(iun,*) 'JOIN TEXT BLUE'
      write(iun,*) 'TITLE TEXT 3.0 6.5 "gluon A1,B1"'
      write(iun,*) 'CASE               "           "'
      write(iun,*) '1.5 6.5'
      write(iun,*) '2.5 6.5'
      write(iun,*) 'JOIN TEXT RED DASHES'
      write(iun,*) 'TITLE TEXT 3.0 6.0 "gluon A1,B1,A2"'
      write(iun,*) 'CASE               "              "'
      write(iun,*) '1.5 6.0'
      write(iun,*) '2.5 6.0'
      write(iun,*) 'JOIN TEXT GREEN DASHES'
      write(iun,*) 'TITLE TEXT 3.0 5.5 "gluon Powheg"'
      write(iun,*) 'CASE               "            "'
      write(iun,*) '1.5 5.5'
      write(iun,*) '2.5 5.5'
      write(iun,*) 'JOIN TEXT BLUE DASHES'
      write(iun,*) 'TITLE BOTTOM "Log(p0T1223/L051223)"'
      write(iun,*) 'CASE         "   S X XX X FX XX XS"'
      write(iun,*) 'SET ORDER X Y DY'

C - Plot Sudakov vs Log(kT/Lambda)
      the_start= Log(0.5d0**2/st_lambda5MSB**2)
      the_end  = Log(the_qh**2/st_lambda5MSB**2)-1.5d0
      the_step = (the_end-the_start)/dble(no_steps)

      isQuark = .true.
      do ixx=1,no_steps
         the_Lq2 = the_start + ixx*the_step
         call i_sudakov_exponent(exp(the_Lq2)*st_lambda5MSB**2,
     $                           the_qh**2,the_qh**2,
     $                           the_Sudakov,isQuark,0,.true.,
     $                           st_lambda5MSB,st_nlight,
     $                           rad_charmthr2,rad_bottomthr2)

         the_Sudakov=exp(the_Sudakov)
         write(iun,'(F8.5,A,F8.6)') the_Lq2,'  ',the_Sudakov
      enddo
      write(iun,*) 'JOIN RED'
      isQuark = .true.
      do ixx=1,no_steps
         the_Lq2 = the_start + ixx*the_step
         call i_sudakov_exponent(exp(the_Lq2)*st_lambda5MSB**2,
     $                           the_qh**2,the_qh**2,
     $                           the_Sudakov,isQuark,2,.true.,
     $                           st_lambda5MSB,st_nlight,
     $                           rad_charmthr2,rad_bottomthr2)
         the_Sudakov=exp(the_Sudakov)
         write(iun,'(F8.5,A,F8.6)') the_Lq2,'  ',the_Sudakov
      enddo
      write(iun,*) 'JOIN GREEN'
      isQuark = .true.
      do ixx=1,no_steps
         the_Lq2 = the_start + ixx*the_step
         call i_sudakov_exponent(exp(the_Lq2)*st_lambda5MSB**2,
     $                           the_qh**2,the_qh**2,
     $                           the_Sudakov,isQuark,1,.true.,
     $                           st_lambda5MSB,st_nlight,
     $                           rad_charmthr2,rad_bottomthr2)
         the_Sudakov=exp(the_Sudakov)
         write(iun,'(F8.5,A,F8.6)') the_Lq2,'  ',the_Sudakov
      enddo
      write(iun,*) 'JOIN BLUE'
      isQuark = .false.
      do ixx=1,no_steps
         the_Lq2 = the_start + ixx*the_step
         call i_sudakov_exponent(exp(the_Lq2)*st_lambda5MSB**2,
     $                           the_qh**2,the_qh**2,
     $                           the_Sudakov,isQuark,0,.true.,
     $                           st_lambda5MSB,st_nlight,
     $                           rad_charmthr2,rad_bottomthr2)

         the_Sudakov=exp(the_Sudakov)
         write(iun,'(F8.5,A,F8.6)') the_Lq2,'  ',the_Sudakov
      enddo
      write(iun,*) 'JOIN RED DASHES'
      isQuark = .false.
      do ixx=1,no_steps
         the_Lq2 = the_start + ixx*the_step
         call i_sudakov_exponent(exp(the_Lq2)*st_lambda5MSB**2,
     $                           the_qh**2,the_qh**2,
     $                           the_Sudakov,isQuark,2,.true.,
     $                           st_lambda5MSB,st_nlight,
     $                           rad_charmthr2,rad_bottomthr2)
         the_Sudakov=exp(the_Sudakov)
         write(iun,'(F8.5,A,F8.6)') the_Lq2,'  ',the_Sudakov
      enddo
      write(iun,*) 'JOIN GREEN DASHES'
      isQuark = .false.
      do ixx=1,no_steps
         the_Lq2 = the_start + ixx*the_step
         call i_sudakov_exponent(exp(the_Lq2)*st_lambda5MSB**2,
     $                           the_qh**2,the_qh**2,
     $                           the_Sudakov,isQuark,1,.true.,
     $                           st_lambda5MSB,st_nlight,
     $                           rad_charmthr2,rad_bottomthr2)
         the_Sudakov=exp(the_Sudakov)
         write(iun,'(F8.5,A,F8.6)') the_Lq2,'  ',the_Sudakov
      enddo
      write(iun,*) 'JOIN BLUE DASHES'

      close(iun)

      end 


C ***************************************************************** C
C -                                                               - C
C ---------------- EXTERNAL NON-MINLO ROUTINES -------------------- C
C ---------------- EXTERNAL NON-MINLO ROUTINES -------------------- C
C ---------------- EXTERNAL NON-MINLO ROUTINES -------------------- C
C -                                                               - C
C ***************************************************************** C


C ***************************************************************** C
C - ALPHA_S QCD at LO                                             - C
C - Program to calculate alfa strong with inf flavours,           - C
C - at leading order.                                             - C
C ***************************************************************** C
      function i_pwhg_alphas0(q2,xlam,inf)
      implicit none
      real * 8 i_pwhg_alphas0,q2,xlam
      integer inf
      real * 8 pi
      parameter (pi=3.141592653589793d0)
      real * 8 b0
      b0=(33-2*inf)/(12*pi)
      i_pwhg_alphas0=1/(b0*log(q2/xlam**2))
      end

C ***************************************************************** C
C - ALPHA_S QCD                                                   - C
C - Program to calculate alfa strong with nf flavours,            - C
C - as a function of lambda with 5 flavors.                       - C
C - The value of alfa is matched at the thresholds q = mq.        - C
C - When invoked with nf < 0 it chooses nf as the number of       - C
C - flavors with mass less then q.                                - C
C -                                                               - C
C - rad_charmthr2    - The mass of the charm quark (threshold)    - C 
C -                    Should be 1.5**2.                          - C 
C - rad_bottomthr2   - The mass of the bottom quark (threshold)   - C 
C -                    Should be 5.0**2.                          - C 
C ***************************************************************** C
      function i_pwhg_alphas(q2,xlam,inf,
     $                       rad_charmthr2,rad_bottomthr2)
      implicit none
      real * 8 i_pwhg_alphas,q2,xlam
      integer inf
      real * 8 pi
      parameter (pi=3.141592653589793d0)
      real * 8 olam,b5,bp5,b4,bp4,b3,bp3,xlc,xlb,xllc,xllb,c45,c35,
     # xmc,xmb
      real * 8 q,xlq,xllq
      integer nf
      data olam/0.d0/
      save olam,b5,bp5,b4,bp4,b3,bp3,xlc,xlb,xllc,xllb,c45,c35,xmc,xmb

      real * 8 rad_charmthr2,rad_bottomthr2

      if(xlam.ne.olam) then
        olam = xlam
        xmc=sqrt(rad_charmthr2)
        xmb=sqrt(rad_bottomthr2)
        b5  = (33-2*5)/pi/12
        bp5 = (153 - 19*5) / pi / 2 / (33 - 2*5)
        b4  = (33-2*4)/pi/12
        bp4 = (153 - 19*4) / pi / 2 / (33 - 2*4)
        b3  = (33-2*3)/pi/12
        bp3 = (153 - 19*3) / pi / 2 / (33 - 2*3)
        xlc = 2 * log(xmc/xlam)
        xlb = 2 * log(xmb/xlam)
        xllc = log(xlc)
        xllb = log(xlb)
        c45  =  1/( 1/(b5 * xlb) - xllb*bp5/(b5 * xlb)**2 )
     #        - 1/( 1/(b4 * xlb) - xllb*bp4/(b4 * xlb)**2 )
        c35  =  1/( 1/(b4 * xlc) - xllc*bp4/(b4 * xlc)**2 )
     #        - 1/( 1/(b3 * xlc) - xllc*bp3/(b3 * xlc)**2 ) + c45
      endif
      q   = sqrt(q2)
      xlq = 2 * log( q/xlam )
      xllq = log( xlq )
      nf = inf
      if( nf .lt. 0) then
        if( q .gt. xmb ) then
          nf = 5
        elseif( q .gt. xmc ) then
          nf = 4
        else
          nf = 3
        endif
      endif
      if    ( nf .eq. 5 ) then
        i_pwhg_alphas = 1/(b5 * xlq) -  bp5/(b5 * xlq)**2 * xllq
      elseif( nf .eq. 4 ) then
        i_pwhg_alphas =
     #    1/( 1/(1/(b4 * xlq) - bp4/(b4 * xlq)**2 * xllq) + c45 )
      elseif( nf .eq. 3 ) then
        i_pwhg_alphas =
     #    1/( 1/(1/(b3 * xlq) - bp3/(b3 * xlq)**2 * xllq) + c35 )
      else
        print *,'error in alfa: unimplemented # of light flavours',nf
        call exit(1)
      endif
      return
      end


C ***************************************************************** C
C -                                                               - C
C - NaN and inf detection.                                        - C
C -                                                               - C
C ***************************************************************** C
      logical function i_pwhg_isfinite(x)
      implicit none
      real * 8 x
c According to ieee standards, a NaN is the only real not
c satisfying x.eq.x.
      if (.not.(x.eq.x)) then
         i_pwhg_isfinite = .false.
         call i_increasecnt('NaN exception')
         return
      endif
c Put constraint to avoid denormals
      if(x.gt.1.or.x.lt.-1) then
         if (1/x.eq.0) then
            i_pwhg_isfinite = .false.
            call i_increasecnt('Inf exception')
            return
         endif
      endif
      i_pwhg_isfinite = .true.
      end


C ***************************************************************** C
C -                                                               - C
C - 3 CERNROUTINES and one ENTRY needed for CERNROUTINE DGAUSS    - C
C -                                                               - C
C ***************************************************************** C
c# 10 "mtlset.F" 2
      SUBROUTINE i_MTLSET(ERC,NLG,MXM,MXR)

      PARAMETER (KTE = 132)
      CHARACTER*6 ERC,CODE(KTE)
      LOGICAL LMF,LRF
      DIMENSION KNTM(KTE),KNTR(KTE)

      DATA ILG /0/

C     renumber the data statements after putting new codes in Unix with:
C     awk -F'[()]' '{ printf"%s(%s)%s(%s)%s(%s)%s\n",$1,NR,$3,NR,$5,NR,$7 }'
C     and modify KTE to the number of lines below

      DATA CODE(1),KNTM(1),KNTR(1) / 'B100.1', 255, 255 /
      DATA CODE(2),KNTM(2),KNTR(2) / 'B300.1', 255, 255 /
      DATA CODE(3),KNTM(3),KNTR(3) / 'B300.2', 255, 255 /
      DATA CODE(4),KNTM(4),KNTR(4) / 'C200.0', 255, 255 /
      DATA CODE(5),KNTM(5),KNTR(5) / 'C200.1', 255, 255 /
      DATA CODE(6),KNTM(6),KNTR(6) / 'C200.2', 255, 255 /
      DATA CODE(7),KNTM(7),KNTR(7) / 'C200.3', 255, 255 /
      DATA CODE(8),KNTM(8),KNTR(8) / 'C201.0', 255, 255 /
      DATA CODE(9),KNTM(9),KNTR(9) / 'C202.0', 255, 255 /
      DATA CODE(10),KNTM(10),KNTR(10) / 'C202.1', 255, 255 /
      DATA CODE(11),KNTM(11),KNTR(11) / 'C202.2', 255, 255 /
      DATA CODE(12),KNTM(12),KNTR(12) / 'C205.1', 255, 255 /
      DATA CODE(13),KNTM(13),KNTR(13) / 'C205.2', 255, 255 /
      DATA CODE(14),KNTM(14),KNTR(14) / 'C207.0', 255, 255 /
      DATA CODE(15),KNTM(15),KNTR(15) / 'C208.0', 255, 255 /
      DATA CODE(16),KNTM(16),KNTR(16) / 'C209.0', 255, 255 /
      DATA CODE(17),KNTM(17),KNTR(17) / 'C209.1', 255, 255 /
      DATA CODE(18),KNTM(18),KNTR(18) / 'C209.2', 255, 255 /
      DATA CODE(19),KNTM(19),KNTR(19) / 'C209.3', 255, 255 /
      DATA CODE(20),KNTM(20),KNTR(20) / 'C210.1', 255, 255 /
      DATA CODE(21),KNTM(21),KNTR(21) / 'C302.1', 255, 255 /
      DATA CODE(22),KNTM(22),KNTR(22) / 'C303.1', 255, 255 /
      DATA CODE(23),KNTM(23),KNTR(23) / 'C304.1', 255, 255 /
      DATA CODE(24),KNTM(24),KNTR(24) / 'C305.1', 255, 255 /
      DATA CODE(25),KNTM(25),KNTR(25) / 'C306.1', 255, 255 /
      DATA CODE(26),KNTM(26),KNTR(26) / 'C307.1', 255, 255 /
      DATA CODE(27),KNTM(27),KNTR(27) / 'C312.1', 255, 255 /
      DATA CODE(28),KNTM(28),KNTR(28) / 'C313.1', 255, 255 /
      DATA CODE(29),KNTM(29),KNTR(29) / 'C315.1', 255, 255 /
      DATA CODE(30),KNTM(30),KNTR(30) / 'C316.1', 255, 255 /
      DATA CODE(31),KNTM(31),KNTR(31) / 'C316.2', 255, 255 /
      DATA CODE(32),KNTM(32),KNTR(32) / 'C320.1', 255, 255 /
      DATA CODE(33),KNTM(33),KNTR(33) / 'C321.1', 255, 255 /
      DATA CODE(34),KNTM(34),KNTR(34) / 'C323.1', 255, 255 /
      DATA CODE(35),KNTM(35),KNTR(35) / 'C327.1', 255, 255 /
      DATA CODE(36),KNTM(36),KNTR(36) / 'C328.1', 255, 255 /
      DATA CODE(37),KNTM(37),KNTR(37) / 'C328.2', 255, 255 /
      DATA CODE(38),KNTM(38),KNTR(38) / 'C328.3', 255, 255 /
      DATA CODE(39),KNTM(39),KNTR(39) / 'C330.1', 255, 255 /
      DATA CODE(40),KNTM(40),KNTR(40) / 'C330.2', 255, 255 /
      DATA CODE(41),KNTM(41),KNTR(41) / 'C330.3', 255, 255 /
      DATA CODE(42),KNTM(42),KNTR(42) / 'C331.1', 255, 255 /
      DATA CODE(43),KNTM(43),KNTR(43) / 'C331.2', 255, 255 /
      DATA CODE(44),KNTM(44),KNTR(44) / 'C334.1', 255, 255 /
      DATA CODE(45),KNTM(45),KNTR(45) / 'C334.2', 255, 255 /
      DATA CODE(46),KNTM(46),KNTR(46) / 'C334.3', 255, 255 /
      DATA CODE(47),KNTM(47),KNTR(47) / 'C334.4', 255, 255 /
      DATA CODE(48),KNTM(48),KNTR(48) / 'C334.5', 255, 255 /
      DATA CODE(49),KNTM(49),KNTR(49) / 'C334.6', 255, 255 /
      DATA CODE(50),KNTM(50),KNTR(50) / 'C336.1', 255, 255 /
      DATA CODE(51),KNTM(51),KNTR(51) / 'C337.1', 255, 255 /
      DATA CODE(52),KNTM(52),KNTR(52) / 'C338.1', 255, 255 /
      DATA CODE(53),KNTM(53),KNTR(53) / 'C340.1', 255, 255 /
      DATA CODE(54),KNTM(54),KNTR(54) / 'C343.1', 255, 255 /
      DATA CODE(55),KNTM(55),KNTR(55) / 'C343.2', 255, 255 /
      DATA CODE(56),KNTM(56),KNTR(56) / 'C343.3', 255, 255 /
      DATA CODE(57),KNTM(57),KNTR(57) / 'C343.4', 255, 255 /
      DATA CODE(58),KNTM(58),KNTR(58) / 'C344.1', 255, 255 /
      DATA CODE(59),KNTM(59),KNTR(59) / 'C344.2', 255, 255 /
      DATA CODE(60),KNTM(60),KNTR(60) / 'C344.3', 255, 255 /
      DATA CODE(61),KNTM(61),KNTR(61) / 'C344.4', 255, 255 /
      DATA CODE(62),KNTM(62),KNTR(62) / 'C345.1', 255, 255 /
      DATA CODE(63),KNTM(63),KNTR(63) / 'C346.1', 255, 255 /
      DATA CODE(64),KNTM(64),KNTR(64) / 'C346.2', 255, 255 /
      DATA CODE(65),KNTM(65),KNTR(65) / 'C346.3', 255, 255 /
      DATA CODE(66),KNTM(66),KNTR(66) / 'C347.1', 255, 255 /
      DATA CODE(67),KNTM(67),KNTR(67) / 'C347.2', 255, 255 /
      DATA CODE(68),KNTM(68),KNTR(68) / 'C347.3', 255, 255 /
      DATA CODE(69),KNTM(69),KNTR(69) / 'C347.4', 255, 255 /
      DATA CODE(70),KNTM(70),KNTR(70) / 'C347.5', 255, 255 /
      DATA CODE(71),KNTM(71),KNTR(71) / 'C347.6', 255, 255 /
      DATA CODE(72),KNTM(72),KNTR(72) / 'C348.1', 255, 255 /
      DATA CODE(73),KNTM(73),KNTR(73) / 'C349.1', 255, 255 /
      DATA CODE(74),KNTM(74),KNTR(74) / 'C349.2', 255, 255 /
      DATA CODE(75),KNTM(75),KNTR(75) / 'C349.3', 255, 255 /
      DATA CODE(76),KNTM(76),KNTR(76) / 'D101.1', 255, 255 /
      DATA CODE(77),KNTM(77),KNTR(77) / 'D103.1', 255, 255 /
      DATA CODE(78),KNTM(78),KNTR(78) / 'D104.1', 255, 255 /
      DATA CODE(79),KNTM(79),KNTR(79) / 'D104.2', 255, 255 /
      DATA CODE(80),KNTM(80),KNTR(80) / 'D105.1', 255, 255 /
      DATA CODE(81),KNTM(81),KNTR(81) / 'D105.2', 255, 255 /
      DATA CODE(82),KNTM(82),KNTR(82) / 'D107.1', 255, 255 /
      DATA CODE(83),KNTM(83),KNTR(83) / 'D110.0', 255, 255 /
      DATA CODE(84),KNTM(84),KNTR(84) / 'D110.1', 255, 255 /
      DATA CODE(85),KNTM(85),KNTR(85) / 'D110.2', 255, 255 /
      DATA CODE(86),KNTM(86),KNTR(86) / 'D110.3', 255, 255 /
      DATA CODE(87),KNTM(87),KNTR(87) / 'D110.4', 255, 255 /
      DATA CODE(88),KNTM(88),KNTR(88) / 'D110.5', 255, 255 /
      DATA CODE(89),KNTM(89),KNTR(89) / 'D110.6', 255, 255 /
      DATA CODE(90),KNTM(90),KNTR(90) / 'D113.1', 255, 255 /
      DATA CODE(91),KNTM(91),KNTR(91) / 'D201.1', 255, 255 /
      DATA CODE(92),KNTM(92),KNTR(92) / 'D202.1', 255, 255 /
      DATA CODE(93),KNTM(93),KNTR(93) / 'D401.1', 255, 255 /
      DATA CODE(94),KNTM(94),KNTR(94) / 'D601.1', 255, 255 /
      DATA CODE(95),KNTM(95),KNTR(95) / 'E210.1', 255, 255 /
      DATA CODE(96),KNTM(96),KNTR(96) / 'E210.2', 255, 255 /
      DATA CODE(97),KNTM(97),KNTR(97) / 'E210.3', 255, 255 /
      DATA CODE(98),KNTM(98),KNTR(98) / 'E210.4', 255, 255 /
      DATA CODE(99),KNTM(99),KNTR(99) / 'E210.5', 255, 255 /
      DATA CODE(100),KNTM(100),KNTR(100) / 'E210.6', 255, 255 /
      DATA CODE(101),KNTM(101),KNTR(101) / 'E210.7', 255, 255 /
      DATA CODE(102),KNTM(102),KNTR(102) / 'E211.0', 255, 255 /
      DATA CODE(103),KNTM(103),KNTR(103) / 'E211.1', 255, 255 /
      DATA CODE(104),KNTM(104),KNTR(104) / 'E211.2', 255, 255 /
      DATA CODE(105),KNTM(105),KNTR(105) / 'E211.3', 255, 255 /
      DATA CODE(106),KNTM(106),KNTR(106) / 'E211.4', 255, 255 /
      DATA CODE(107),KNTM(107),KNTR(107) / 'E406.0', 255, 255 /
      DATA CODE(108),KNTM(108),KNTR(108) / 'E406.1', 255, 255 /
      DATA CODE(109),KNTM(109),KNTR(109) / 'E407.0', 255, 255 /
      DATA CODE(110),KNTM(110),KNTR(110) / 'E408.0', 255, 255 /
      DATA CODE(111),KNTM(111),KNTR(111) / 'E408.1', 255, 255 /
      DATA CODE(112),KNTM(112),KNTR(112) / 'F500.0', 255, 255 /
      DATA CODE(113),KNTM(113),KNTR(113) / 'F500.1', 255, 255 /
      DATA CODE(114),KNTM(114),KNTR(114) / 'F500.2', 255, 255 /
      DATA CODE(115),KNTM(115),KNTR(115) / 'F500.3', 255, 255 /
      DATA CODE(116),KNTM(116),KNTR(116) / 'G100.1', 255, 255 /
      DATA CODE(117),KNTM(117),KNTR(117) / 'G100.2', 255, 255 /
      DATA CODE(118),KNTM(118),KNTR(118) / 'G101.1', 255, 255 /
      DATA CODE(119),KNTM(119),KNTR(119) / 'G101.2', 255, 255 /
      DATA CODE(120),KNTM(120),KNTR(120) / 'G105.1', 255, 255 /
      DATA CODE(121),KNTM(121),KNTR(121) / 'G106.1', 255, 255 /
      DATA CODE(122),KNTM(122),KNTR(122) / 'G106.2', 255, 255 /
      DATA CODE(123),KNTM(123),KNTR(123) / 'G116.1', 255, 255 /
      DATA CODE(124),KNTM(124),KNTR(124) / 'G116.2', 255, 255 /
      DATA CODE(125),KNTM(125),KNTR(125) / 'H101.0', 255, 255 /
      DATA CODE(126),KNTM(126),KNTR(126) / 'H101.1', 255, 255 /
      DATA CODE(127),KNTM(127),KNTR(127) / 'H101.2', 255, 255 /
      DATA CODE(128),KNTM(128),KNTR(128) / 'H301.1', 255, 255 /
      DATA CODE(129),KNTM(129),KNTR(129) / 'U501.1', 255, 255 /
      DATA CODE(130),KNTM(130),KNTR(130) / 'V202.1', 255, 255 /
      DATA CODE(131),KNTM(131),KNTR(131) / 'V202.2', 255, 255 /
      DATA CODE(132),KNTM(132),KNTR(132) / 'V202.3', 255, 255 /

c# 175 "mtlset.F"

      ILG=NLG
      L=0
      IF(ERC .NE. ' ') THEN
       DO 10 L = 1,6
       IF(ERC(1:L) .EQ. ERC) GOTO 12
   10  CONTINUE
   12  CONTINUE
      ENDIF
      DO 14 I = 1,KTE
      IF(L .EQ. 0 .OR. CODE(I)(1:L) .EQ. ERC(1:L)) THEN
       IF(MXM .GE. 0) KNTM(I)=MXM
       IF(MXR .GE. 0) KNTR(I)=MXR
      ENDIF
   14 CONTINUE
      RETURN

      ENTRY i_MTLMTR(ERC,MLG,LMF,LRF)

      MLG=ILG
      DO 20 I = 1,KTE
      IF(ERC .EQ. CODE(I))  GOTO 21
   20 CONTINUE
      WRITE(*,100) ERC
      CALL i_ABEND
      RETURN

   21 LMF=KNTM(I) .GE. 1
      LRF=KNTR(I) .GE. 1
      IF(LMF .AND. KNTM(I) .LT. 255)  KNTM(I)=KNTM(I)-1
      IF(LRF .AND. KNTR(I) .LT. 255)  KNTR(I)=KNTR(I)-1
      IF(.NOT.LRF) THEN
       IF(ILG .LT. 1) WRITE(  *,101) CODE(I)
       IF(ILG .GE. 1) WRITE(ILG,101) CODE(I)
      ENDIF
      RETURN
  100 FORMAT(7X,'***** CERN N002 MTLSET ... ERROR N002: ',
     1'ERROR CODE ',A6,' NOT RECOGNIZED BY ERROR MONITOR. RUN ABORTED.')
  101 FORMAT(7X,'***** CERN N002 MTLSET ... ERROR NOO2.1: ',
     1'RUN TERMINATED BY LIBRARY ERROR CONDITION ',A6)
      END

c# 10 "mtlprt.F" 2
      SUBROUTINE i_MTLPRT(NAME,ERC,TEXT)
      CHARACTER*(*) NAME,ERC,TEXT
      LOGICAL LMF,LRF

      IF(ERC(5:6).NE.'.0') THEN
        CALL i_MTLMTR(ERC,MLG,LMF,LRF)
      ELSE
        LMF=.TRUE.
        LRF=.FALSE.
      ENDIF
      IF(LMF) THEN
        LT=i_LENOCC(TEXT)
        IF(MLG .LT. 1) WRITE(  *,100) ERC(1:4),NAME,ERC,TEXT(1:LT)
        IF(MLG .GE. 1) WRITE(MLG,100) ERC(1:4),NAME,ERC,TEXT(1:LT)
      ENDIF
      IF(.NOT.LRF) CALL i_ABEND
      RETURN
100   FORMAT(7X,'***** CERN ',A,1X,A,' ERROR ',A,': ',A)
      END

      SUBROUTINE i_ABEND
C
C CERN PROGLIB# Z035    ABEND           .VERSION KERNFOR  4.31  911111
C ORIG.  8/02/88  JZ
C

      STOP  7
      END

c# 1 "lenocc.F"
c# 1 "<built-in>"
c# 1 "<command line>"
c# 1 "lenocc.F"
*
* $Id: lenocc.F,v 1.1.1.1 1996/02/15 17:49:49 mclareni Exp $
*
* $Log: lenocc.F,v $
* Revision 1.1.1.1  1996/02/15 17:49:49  mclareni
* Kernlib
*
*
c# 1 "/usr/local/home/video/cernlib/2005/include/kerngen/pilot.h" 1
c# 94 "/usr/local/home/video/cernlib/2005/include/kerngen/pilot.h"





c# 10 "lenocc.F" 2
      FUNCTION i_LENOCC (CHV)
C
C CERN PROGLIB# M507    LENOCC          .VERSION KERNFOR  4.21  890323
C ORIG. March 85, A.Petrilli, re-write 21/02/89, JZ
C
C-    Find last non-blank character in CHV

      CHARACTER    CHV*(*)

      N = LEN(CHV)

      DO 17  JJ= N,1,-1
      IF (CHV(JJ:JJ).NE.' ') GO TO 99
   17 CONTINUE
      JJ = 0

   99 i_LENOCC = JJ
      RETURN
      END

C ***************************************************************** C
C -                                                               - C
C - Finds a new unit number to open a file on.                    - C
C -                                                               - C
C ***************************************************************** C
      subroutine i_newunit(iun)
      implicit none
      integer iun
      logical ok
      integer j
      do j=10,99
         inquire(unit=j,opened=ok)
         if(.not.ok) then
            iun=j
            return
         endif
      enddo
      end

      subroutine i_increasecnt(string)
      implicit none
      character *(*) string
      integer maxnum
      parameter (maxnum=100)
      character * 100 i_keywords(maxnum)
      real * 8 i_counters(maxnum)
      integer i_ncounters
      common/i_ccounters/i_keywords,i_counters,i_ncounters
      integer ini,j
      call i_initcnt
      do j=1,i_ncounters
         if(string.eq.i_keywords(j)) then
            i_counters(j)=i_counters(j)+1
            return
         endif
      enddo
c not found
      if(i_ncounters.eq.maxnum) then
         write(*,*) 'ERROR: increasecnt too many counters requested'
         stop
      endif
      i_ncounters=i_ncounters+1
      i_keywords(i_ncounters)=string
      i_counters(i_ncounters)=1
      end

      subroutine i_initcnt
      implicit none
      integer maxnum
      parameter (maxnum=100)
      logical ini
      character * 100 i_keywords(maxnum)
      real * 8 i_counters(maxnum)
      integer i_ncounters
      common/i_ccounters/i_keywords,i_counters,i_ncounters
      data ini/.true./
      save ini
      if(ini) then
         i_ncounters=0
         ini=.false.
      endif
      end


      subroutine i_findubfsr(nleg,i,j,cmpin,cmpout)
      implicit none
      integer nleg
      integer i,j
      real * 8 cmpin(0:3,20)
      real * 8 cmpout(0:3,20)
      real * 8 krecv(3),q0,q2,krec,beta,
     1 k0rec,k,vec(3)
      cmpout(0:3,1)=cmpin(0:3,1)
      cmpout(0:3,2)=cmpin(0:3,2)
      q0=2*cmpout(0,1)
      q2=q0**2
c recoil system momentum 
      k0rec=q0-cmpin(0,i)-cmpin(0,j)
      krecv=-cmpin(1:3,i)-cmpin(1:3,j)
      krec=sqrt(krecv(1)**2+krecv(2)**2+krecv(3)**2)
      beta=(q2-(k0rec+krec)**2)/(q2+(k0rec+krec)**2)
      vec=krecv/krec
      call i_mboost(20-2,vec,beta,
     1     cmpin(:,3:20),cmpout(:,3:20))
      k=(q2-(k0rec**2-krec**2))/(2*q0)
      cmpout(0,i)=k
      cmpout(1:3,i)=-vec*k
      cmpout(:,j)=0
      end

      subroutine i_findubisr(nleg,j,cmpin,cmpout)
      implicit none
      integer nleg
      integer ixx,j
      real * 8 cmpin(0:3,20)
      real * 8 cmpout(0:3,20)
      real * 8 krecv(3),q0,q2,krec,k0rec,
     1     krecperp,mrec2,beta,vec(3)
      real * 8 En,pZ
      cmpout(0:3,1)=cmpin(0:3,1)
      cmpout(0:3,2)=cmpin(0:3,2)
      q0=2*cmpout(0,1)
      q2=q0**2
c recoil system momentum 
      k0rec=q0-cmpin(0,j)
      krecv=-cmpin(1:3,j)
      krec=sqrt(krecv(1)**2+krecv(2)**2+krecv(3)**2)
      mrec2=(k0rec**2-krec**2)
      beta=-krecv(3)/k0rec
      vec(1)=0
      vec(2)=0
      vec(3)=1
      call i_mboost(20-2,vec,beta,
     1     cmpin(:,3:20),cmpout(:,3:20))
c Now the transverse boost
      krecperp=sqrt(krecv(1)**2+krecv(2)**2)
      vec(3)=0
      vec(1:2)=krecv(1:2)/krecperp
      beta=-krecperp/sqrt(mrec2+krecperp**2)
      call i_mboost(20-2,vec,beta,
     1     cmpout(:,3:20),cmpout(:,3:20))
      cmpout(:,j)=0
C - I.S. momenta
      En = 0d0
      pz = 0d0
      do ixx=3,nleg
         En = En + cmpout(0,ixx)
         pz = pz + cmpout(3,ixx)
      enddo
      cmpout(0,1) =  0.5d0*(En+pz)
      cmpout(1,1) =  0d0
      cmpout(2,1) =  0d0
      cmpout(3,1) =  0.5d0*(En+pz)
      cmpout(0,2) =  0.5d0*(En-pz)
      cmpout(1,2) =  0d0
      cmpout(2,2) =  0d0
      cmpout(3,2) = -0.5d0*(En-pz)
      end

      subroutine i_mboost(m,vec,beta,vin,vout)
c     boosts the m vectors vin(0:3,m) into the vectors vout(0:3,m) (that can
c     be the same) in the direction of vec(3) (|vec|=1) with velocity
c     beta.  Lorents convention: (t,x,y,z).
      implicit none
      integer m
      real * 8 vec(3),beta,vin(0:3,m),vout(0:3,m)
      real * 8 betav,gamma
      real * 8 vdotb
      integer ipart,idim
      gamma=1/sqrt(1-beta**2)
      do ipart=1,m
         vdotb=vin(1,ipart)*vec(1)
     #         +vin(2,ipart)*vec(2)+vin(3,ipart)*vec(3)
         do idim=1,3
            vout(idim,ipart)=vin(idim,ipart)
     #           +vec(idim)*((gamma-1)*vdotb
     #           +gamma*beta*vin(0,ipart))
         enddo
         vout(0,ipart)=gamma*(vin(0,ipart)+vdotb*beta)
      enddo
      end

C ---------------------------------------------------------- C 
C - TEST CODE SETUP ---------------------------------------- C
C - TEST CODE SETUP ---------------------------------------- C
C - TEST CODE SETUP ---------------------------------------- C
C ---------------------------------------------------------- C 

C - How to hook up and test code in the Powheg-Box:
C - Replace subrountine setlocalscales with the version 
C - below, which hooks in the independent code here after
C - it calls ist own native version (full of dependencies,
C - includes and the like). You can then compare what comes
C - out for rescfac, basicfac, nlofac and bornfac from the
C - original version to what comes out in the independent
C - version (i_rescfac, i_basicfac, i_nlofac and i_bornfac).
C - I found them all to be identical, using Z2jets, which
C - I take to mean the independent version is a genuine
C - replica.

c$$$      subroutine setlocalscales(iuborn,imode,rescfac)
c$$$c returns the rescaling factor including sudakov form factors and
c$$$c coupling rescaling, for born (imode=1) and NLO corrections (imode=2)
c$$$      implicit none
c$$$      integer iuborn,imode
c$$$      real * 8 rescfac
c$$$      include 'nlegborn.h'
c$$$      include 'pwhg_flst.h'
c$$$      include 'pwhg_kn.h'
c$$$      include 'pwhg_flg.h'
c$$$      include 'pwhg_par.h'
c$$$      include 'pwhg_st.h'
c$$$      integer j,k,mu
c$$$      logical samekin,sameflv
c$$$      integer flav(nlegborn)
c$$$      real * 8 op(0:3,nlegborn),basicfac,bornfac,nlofac
c$$$      real * 8 savbasicfac(maxprocborn),savbornfac(maxprocborn),
c$$$     1     savnlofac(maxprocborn),savmuf2(maxprocborn)
c$$$      logical valid(maxprocborn)
c$$$      save savbasicfac,savbornfac,savnlofac,savmuf2,valid
c$$$      save op
c$$$      data valid/maxprocborn*.false./
c$$$C - Below is just for indy code:
c$$$      include 'pwhg_rad.h'
c$$$      integer  i_flav(20)
c$$$      real * 8 i_kn_pborn(0:3,20),i_kn_cmpborn(0:3,20)
c$$$      real * 8 i_rescfac,i_basicfac,i_bornfac,i_nlofac,i_st_mufact2
c$$$C - Above is just for indy code.
c$$$      do j=1,nlegborn
c$$$         do mu=0,3
c$$$            if(op(mu,j).ne.kn_cmpborn(mu,j)) then
c$$$               do k=1,flst_nborn
c$$$                  valid(k)=.false.
c$$$               enddo
c$$$               op=kn_cmpborn
c$$$            endif
c$$$         enddo
c$$$      enddo
c$$$ 33   continue
c$$$      if(valid(iuborn)) then
c$$$         basicfac=savbasicfac(iuborn)
c$$$         nlofac=savnlofac(iuborn)
c$$$         bornfac=savbornfac(iuborn)
c$$$         st_mufact2=savmuf2(iuborn)
c$$$      else
c$$$         flav=flst_born(:,iuborn)
c$$$         call setlocalscales0(flav,kn_cmpborn,
c$$$     1        basicfac,bornfac,nlofac)
c$$$         savbasicfac(iuborn)=basicfac
c$$$         savnlofac(iuborn)=nlofac
c$$$         savbornfac(iuborn)=bornfac
c$$$         savmuf2(iuborn)=st_mufact2
c$$$         valid(iuborn)=.true.
c$$$      endif
c$$$      if(imode.eq.1) then
c$$$         rescfac=basicfac*bornfac
c$$$      elseif(imode.eq.2) then
c$$$         rescfac=basicfac*nlofac
c$$$      endif
c$$$
c$$$C - Now hook on independent code version and see what happens:
c$$$      do j=1,20
c$$$         i_flav(j)=-999
c$$$         if(j.le.nlegborn) i_flav(j)=flst_born(j,iuborn)
c$$$         if(j.le.nlegborn) then
c$$$            do mu=0,3
c$$$               i_kn_pborn(mu,j)=kn_pborn(mu,j)
c$$$               i_kn_cmpborn(mu,j)=kn_cmpborn(mu,j)
c$$$            enddo
c$$$         else
c$$$            do mu=0,3
c$$$               i_kn_pborn(mu,j)=0d0
c$$$               i_kn_cmpborn(mu,j)=0d0
c$$$            enddo
c$$$         endif
c$$$      enddo
c$$$      call i_setlocalscales(
c$$$C - New extra inputs for independence from includes
c$$$     $      nlegborn,        ! Number of final state particle in Born plus 2.
c$$$     $      st_alpha,        ! Fixed aS used in MEs before they got here
c$$$     $      st_bornorder,    ! Number of powers of aS in the Born (ZJ=1, ZJJ=2)
c$$$     $      st_muren2,       ! The mu_R^2 value that was used to evaluate virt
c$$$     $      flg_bornonly,    ! Are we feeding through only Born stuff, no NLO?
c$$$     $      rad_charmthr2,   ! Charm  (mass) threshold squared: (1.5 GeV)^2
c$$$     $      rad_bottomthr2,  ! Bottom (mass) threshold squared: (5.0 GeV)^2
c$$$     $      st_nlight,       ! Nominal number of light flavs in computation (5)
c$$$     $      st_lambda5MSB,   ! Lambda_5 in the MS bar scheme
c$$$     $      st_renfact,      ! Scaling factor for mu_R (not mu_R^2): 0.5,1,2.0
c$$$     $      st_facfact,      ! Scaling factor for mu_F (not mu_F^2): 0.5,1,2.0
c$$$     $      .true.,          ! Freeze clustering scale if it attempts to
c$$$     $                       ! go down - no more sudakovs/subtractions, only
c$$$     $                       ! coupling reweighting (for all remaining nodes).
c$$$C - Old inputs and output
c$$$     $      i_kn_cmpborn,    ! Born particles momenta in C.O.M frame
c$$$     $      i_flav,          ! Born particle flavours (1,2) is incoming
c$$$     $      imode,           ! imode=1 for Born, imode=2 for all NLO contribs
c$$$     $      i_rescfac,       ! Comb. of suds, subtractions & aS weights to
c$$$     $                       ! rescale the weight of the contrib assoc. to this
c$$$     $                       ! mom. and flav. configuration.
c$$$     $      i_basicfac,      ! Comb. of suds, subtractions & aS weights to
c$$$     $      i_bornfac,       ! Comb. of suds, subtractions & aS weights to
c$$$     $      i_nlofac,        ! Comb. of suds, subtractions & aS weights to
c$$$C - Extra factorisation scale output as it normally lives in pwhg common block
c$$$     $      i_st_mufact2)    ! The mu_F^2 value that the PDFs should be 
c$$$                             ! re-evaluated at now
c$$$
c$$$c$$$      if(abs((rescfac   -i_rescfac )  /rescfac )  .gt.1d-127.or.
c$$$c$$$     $   abs((basicfac  -i_basicfac)  /basicfac)  .gt.1d-127.or.
c$$$c$$$     $   abs((nlofac    -i_nlofac)    /nlofac  )  .gt.1d-127.or.
c$$$c$$$     $   abs((bornfac   -i_bornfac)   /bornfac )  .gt.1d-127.or.
c$$$c$$$     $   abs((st_mufact2-i_st_mufact2)/st_mufact2).gt.1d-127) then
c$$$         write(*,*) '--------------------------------------'
c$$$         write(*,*) 'result of subsequent call to indy setlocalscales:'
c$$$         write(*,*) 'rescfac    = ',rescfac,  i_rescfac
c$$$         write(*,*) 'basicfac   = ',basicfac, i_basicfac
c$$$         write(*,*) 'bornfac    = ',bornfac,  i_bornfac
c$$$         write(*,*) 'nlofac     = ',nlofac,   i_nlofac
c$$$         write(*,*) 'st_mufact2 = ',st_mufact2,i_st_mufact2
c$$$c$$$      endif
c$$$
c$$$      end
