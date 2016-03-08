

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
C - using dgauss. The numerical integration and analytic results  - C
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
C ***************************************************************** C
      subroutine sudakov_exponent(q2l,q2h,m2,theExponent,isQuark)
      implicit none
      real * 8 q2l,q2h,m2,theExponent
      logical  isQuark
      integer  nf
      real * 8 bnf,bpnf,K
      real * 8 A1,B1,A2,B2
      real * 8 Lq2l,Lq2h,Lm2
      real * 8 aSbar
      real * 8 theA1coeff,theA2coeff,theB1coeff,theB2coeff
      real * 8 st_lambda5MSB
      real * 8 pi
      parameter (pi=3.141592653589793238462643383279502884197D0)
      real * 8 CF,CA
      parameter (CF=4d0/3, CA=3d0)

C - This function is the only external dependency, all it is is the 2-loop
C - alpha_S:
      real * 8 pdfFromLHAPDF
      external pdfFromLHAPDF

C - This is what comes out when I run the code with CTEQ66 using LHAPDF.
C - Apparently it corresponds to alpha_s(Mz) = 0.11798. I put this here to
C - avoid using the Powheg-Box include which uses it. Powheg-Box will
C - get this from the PDF in general i.e. in the code I'm using it's not
C - hardwired.
      st_lambda5MSB=0.226000000000003

C - In practice the Powheg-Box Minlo Sudakov is always using 5 flavours.
      nf = 5

      bnf  = (11d0*CA-2d0*nf)/12/Pi
      bpnf = (153 - 19d0*nf) / Pi / 2 / (33 - 2*nf)

      K  =  (67d0/18-Pi**2/6)*CA-5d0/9*nf

      if(isQuark) then
         A1 =  Cf
         A2 =  Cf*K
         B1 = -3d0/2*Cf
         B2 =  0d0
      else
         A1 =  CA
         A2 =  CA*K
         B1 = -2*Pi*bnf
         B2 =  0d0
      endif

      Lq2l = Log(q2l/st_lambda5MSB/st_lambda5MSB)
      Lq2h = Log(q2h/st_lambda5MSB/st_lambda5MSB)
      Lm2  = Log(m2/st_lambda5MSB/st_lambda5MSB)

      aSbar = pdfFromLHAPDF(q2l)/2/Pi

      theA1coeff =
     $     ( (Lq2h - Lq2l)
     $       - Lm2*Log(Lq2h/Lq2l)
     $     )/(2*bnf*Pi)
     $     + bpnf*( 2*Lm2*(Lq2h - Lq2l)
     $            + 2*Lm2*Lq2h*Log(Lq2l)
     $            + Lq2l*Lq2h*Log(Lq2l)**2 
     $            - 2*Lm2*Lq2l*Log(Lq2h)
     $            - Lq2l*Lq2h*Log(Lq2h)**2
     $            )/(4*bnf**2*Lq2l*Lq2h*Pi)

      theA2coeff = 
     $     ( Lm2*(Lq2l - Lq2h)
     $     - Lq2l*Lq2h*Log(Lq2l)
     $     + Lq2l*Lq2h*Log(Lq2h)
     $     )/(4d0*bnf**2*Lq2l*Lq2h*Pi**2)
     $     + bpnf*( 0.5*(Lq2h-Lq2l)*( Lm2*(Lq2l + Lq2h)
     $                              - 4*Lq2l*Lq2h
     $                              )
     $            + (Lm2 - 2*Lq2l)*Lq2h**2*Log(Lq2l)
     $            - Lq2l**2*(Lm2 - 2*Lq2h)*Log(Lq2h)
     $            )/(4d0*bnf**3*Lq2l**2*Lq2h**2*Pi**2)
     $     + bpnf**2*(  (Lq2l - Lq2h)*( 8*Lm2*( Lq2l**2
     $                                        + Lq2l*Lq2h
     $                                        + Lq2h**2
     $                                        )
     $                                - 27*Lq2l*Lq2h*(Lq2l + Lq2h)
     $                                )
     $               - 6*(4*Lm2 - 9*Lq2l)*Lq2h**3*Log(Lq2l)
     $               + 18*(-2*Lm2 + 3*Lq2l)*Lq2h**3*Log(Lq2l)**2
     $               + 6*Lq2l**3*(4*Lm2 - 9*Lq2h)*Log(Lq2h)
     $               + 18*Lq2l**3*(2*Lm2 - 3*Lq2h)*Log(Lq2h)**2
     $               )/(432.*bnf**4*Lq2l**3*Lq2h**3*Pi**2)

      theB1coeff = 
     $     - Log(Lq2h/Lq2l)/(2.*bnf*Pi)
     $     - bpnf*( Lq2l - Lq2h - Lq2h*Log(Lq2l)
     $            + Lq2l*Log(Lq2h)
     $            )/(2.*bnf**2*Lq2l*Lq2h*Pi) 

      theB2coeff =
     $     - (1/Lq2l - 1/Lq2h)/(4.*bnf**2*Pi**2)
     $     - bpnf*(   Lq2l**2 - Lq2h**2 - 2*Lq2h**2*Log(Lq2l)
     $            + 2*Lq2l**2*Log(Lq2h)
     $            )/(8.*bnf**3*Lq2l**2*Lq2h**2*Pi**2)
     $     - bpnf**2*( 2*Lq2h**3 - 2*Lq2l**3
     $               + 3*Lq2h**3*Log(Lq2l)*(2 + 3*Log(Lq2l)) 
     $               - 3*Lq2l**3*Log(Lq2h)*(2 + 3*Log(Lq2h))
     $               )/(108.*bnf**4*Lq2l**3*Lq2h**3*Pi**2)

      theExponent = A1*theA1coeff + A2*theA2coeff
     $            + B1*theB1coeff + B2*theB2coeff

c$$$      write(*,*) 'theExponent = ',theExponent
c$$$      write(*,*) 'theA1coeff  = ',theA1coeff
c$$$      write(*,*) 'theA2coeff  = ',theA2coeff
c$$$      write(*,*) 'theB1coeff  = ',theB1coeff
c$$$      write(*,*) 'theB2coeff  = ',theB2coeff
c$$$      write(*,*) 'A1,B1,A2,B2 = ',A1,B1,A2,B2
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
C - sudakov - The Sudakov form factor.           - C
C -                                              - C
C ------------------------------------------------ C
      function nll_sudakov(q20,q2h,q2l,flav)
      implicit none
      real * 8 nll_sudakov,q2h,q2l,q20
      integer flav
      real * 8 lam2
      logical isQuark
      real * 8 theExponentN,theExponentD
      logical ini
      data ini/.true./
      save ini
      real * 8 st_lambda5MSB
      st_lambda5MSB=0.226000000000003

      lam2=st_lambda5MSB**2
      if(q20.le.lam2.or.q2l.lt.lam2.or.q2h.lt.lam2) then
         nll_sudakov=0
         goto 999
      endif
      if(q2l.ge.q2h.or.q2h.le.q20) then
         nll_sudakov=1
         goto 999
      endif
      if(flav.eq.0) then
         isQuark=.false.
      else
         isQuark=.true.
      endif
      if(q2l.le.q20) then
        call sudakov_exponent(q20,q2h,q2h,theExponentN,
     $                         isQuark)
        nll_sudakov=exp(theExponentN)
      else
         call sudakov_exponent(q20,q2h,q2h,theExponentN,
     $                         isQuark)
         call sudakov_exponent(q20,q2l,q2l,theExponentD,
     $                         isQuark)
         nll_sudakov=exp(theExponentN-theExponentD)
      endif
 999  continue
      end
