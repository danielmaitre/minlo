
C ------------------------------------------------ C
C - Inputs:                                      - C
C - *******                                      - C
C - q2h  - Upper node scale / bound on Sudakov   - C
C - q2l  - Lower node scale / bound on Sudakov   - C
C - flav - flavour index for the evolving parton - C
C -                                              - C
C - Outputs:                                     - C
C - ********                                     - C
C - sudakov - The Sudakov form factor.           - C
C -           Uses 1-loop alpha, no CMW scheme   - C
C -           etc.                               - C
C -                                              - C
C ------------------------------------------------ C
      function LO_sudakov(q20,q2h,q2l,flav)
      implicit none
      real * 8 LO_sudakov,q2h,q2l,q20
      integer flav
C      include 'pwhg_st.h'
      integer st_nlight
      real * 8 st_lambda5MSB
C      include 'pwhg_math.h'
      real * 8 pi
      real * 8 b0,c,b,lam2
      st_nlight = 5
      pi = 3.1415926535
      st_lambda5MSB = 0.2
      lam2=st_lambda5MSB**2
      if(q20.le.lam2.or.q2l.lt.lam2.or.q2h.lt.lam2) then
         LO_sudakov=0
         goto 999
      endif
      if(q2l.ge.q2h.or.q2h.le.q20) then
         LO_sudakov=1
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
         LO_sudakov= exp(
     1        -c/b0*(  log(log(q2h/lam2)/log(q20/lam2))
     2        *(0.5d0*log(q2h/lam2)-b)
     3        -0.5d0*log(q2h/q20)
     4        ))
      else
         LO_sudakov= exp(
     1        -c/b0*(  log(log(q2h/lam2)/log(q20/lam2))
     2        *(0.5d0*log(q2h/lam2)-b)
     3        -0.5d0*log(q2h/q20)
     4        )
     5        +c/b0*(  log(log(q2l/lam2)/log(q20/lam2))
     6        *(0.5d0*log(q2l/lam2)-b)
     7        -0.5d0*log(q2l/q20)
     8        ))
      endif
 999  continue
      end

C ------------------------------------------------ C
C - Inputs:                                      - C
C - *******                                      - C
C - q2h  - Upper node scale / bound on Sudakov   - C
C - q2l  - Lower node scale / bound on Sudakov   - C
C - flav - flavour index for the evolving parton - C
C -                                              - C
C - Outputs:                                     - C
C - ********                                     - C
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
      function expsudakov(q20,q2h,q2l,flav)
      implicit none
      real * 8 expsudakov,q2h,q2l,q20
      integer flav
C      include 'pwhg_st.h'
      integer st_nlight
      real * 8 st_lambda5MSB
C     include 'pwhg_math.h'
      real * 8 pi
C      include 'pwhg_flg.h'
      logical flg_bornonly
      real * 8 b0,c,b,lam2
      flg_bornonly = .false.
      st_lambda5MSB = 0.2
      pi = 3.1415926535
      st_nlight = 5
      lam2=st_lambda5MSB**2
      if(q20.le.lam2.or.q2l.lt.lam2.or.q2h.lt.lam2) then
c in this case everything is zero, irrelevant
         expsudakov=0
         return
      endif
      if(q2l.ge.q2h.or.q2h.le.q20.or.flg_bornonly) then
         expsudakov=0
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
         expsudakov=
     1        c/pi*(0.25d0*log(q2h/q20)**2 - b*log(q2h/q20))
      else
         expsudakov=
     1        c/pi*(0.25d0*log(q2h/q20)**2 - b*log(q2h/q20))
     2      - c/pi*(0.25d0*log(q2l/q20)**2 - b*log(q2l/q20))
      endif
      end


