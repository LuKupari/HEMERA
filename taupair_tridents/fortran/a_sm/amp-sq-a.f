      subroutine projections(v2,t,capY,sx,mx2,dp,dm,
     x ppdp,pmdp,
     x b11,b12,b21,b22,b11bsm,b12bsm,b21bsm,b22bsm)
c      implicit real*8 (a-z)
      implicit none
      real*8 m,mmu,mp,mhadmin,dm2,mmu2,M2,v2,t,capY,sx,mx2,dp,dm,
     x         F2,ppdp,pmdp,b11,b12,b21,b22,capT,ppdq,pmdq,
     x	       b11bsm,b12bsm,b21bsm,b22bsm,ppdbigQ,pmdbigQ
      common/masses/m,mmu,mp,mhadmin,dm2,mmu2,M2
      common/bsm/F2
      
      capT = M2-t-mx2
      
      if (dp.eq.0.d0.or.dm.eq.0.d0)then
         b11 = 0.
         b12 = 0.
         b21 = 0.
         b22 = 0.
         b11bsm = 0.
         b12bsm = 0.
         b21bsm = 0.
         b22bsm = 0.
         write(6,*)'dp,dm',dp,dm
      else

c              B11        Terms in output =         39
c              B12        Terms in output =         22
c              B21        Terms in output =         37
c              B22        Terms in output =         35
c              B11BSM       Terms in output =       93
c              B12BSM        Terms in output =      53
c              B21BSM        Terms in output =      119
c              B22BSM        Terms in output =      98

         ppdq = (-t-dp)/2.
         pmdq = (-t-dm)/2.
         ppdbigQ = (-capY-dm)/2.
         pmdbigQ = (-capY-dp)/2.

      B11 =
     &  + t*Dp**(-2) * (  - m**2*M2 - 1./2.*M2*v2 + 2*ppdp*pmdp )
      B11 = B11 + t*Dp**(-1)*Dm**(-1) * (  - 2*m**2*M2 + M2*v2 )
      B11 = B11 + t*Dm**(-2) * (  - m**2*M2 - 1./2.*M2*v2 + 2*ppdp*pmdp
     &     )
      B11 = B11 + Dp**(-2)*capT * ( 2*m**2*pmdp + 2*pmdp*ppdq )
      B11 = B11 + Dp**(-2) * ( m**2*M2*v2 - 2*m**2*M2*ppdq - 2*m**2*M2*
     &    pmdq - 4*m**2*ppdp*pmdp - 2*M2*ppdq*pmdq )
      B11 = B11 + Dp**(-1)*Dm**(-1)*capT * ( 2*m**2*ppdp + 2*m**2*pmdp
     &     - v2*ppdp - v2*pmdp )
      B11 = B11 + Dp**(-1)*Dm**(-1)*capT**2 * (  - m**2 )
      B11 = B11 + Dp**(-1)*Dm**(-1) * ( 2*m**2*M2*v2 - 4*m**2*M2*ppdq
     &     - 4*m**2*M2*pmdq - 8*m**2*ppdp*pmdp + 2*M2*v2*ppdq + 2*M2*v2
     &    *pmdq - M2*v2**2 + 4*v2*ppdp*pmdp - 4*ppdp*pmdp*ppdq - 4*ppdp
     &    *pmdp*pmdq + 4*ppdp**2*pmdq + 4*pmdp**2*ppdq )
      B11 = B11 + Dm**(-2)*capT * ( 2*m**2*ppdp + 2*ppdp*pmdq )
      B11 = B11 + Dm**(-2) * ( m**2*M2*v2 - 2*m**2*M2*ppdq - 2*m**2*M2*
     &    pmdq - 4*m**2*ppdp*pmdp - 2*M2*ppdq*pmdq )

      B12 =
     &  + t*Dp**(-2) * (  - 6*m**2 - v2 )
      B12 = B12 + t*Dp**(-1)*Dm**(-1) * (  - 4*m**2 + 4*v2 )
      B12 = B12 + t*Dm**(-2) * (  - 6*m**2 - v2 )
      B12 = B12 + Dp**(-2) * ( 2*m**2*v2 - 8*m**2*ppdq - 4*m**2*pmdq + 
     &    4*m**4 - 4*ppdq*pmdq )
      B12 = B12 + Dp**(-1)*Dm**(-1) * (  - 4*m**2*ppdq - 4*m**2*pmdq + 
     &    8*m**4 + 4*v2*ppdq + 4*v2*pmdq - 2*v2**2 )
      B12 = B12 + Dm**(-2) * ( 2*m**2*v2 - 4*m**2*ppdq - 8*m**2*pmdq + 
     &    4*m**4 - 4*ppdq*pmdq )

      B21 =
     &  + t*Dp**(-2) * ( 1./4.*M2**2*v2 )
      B21 = B21 + t*Dp**(-1)*Dm**(-1) * (  - 1./2.*M2**2*v2 )
      B21 = B21 + t*Dm**(-2) * ( 1./4.*M2**2*v2 )
      B21 = B21 + Dp**(-2)*capT * ( 1./2.*M2*v2*ppdp - M2*ppdp*pmdq - 
     &    M2*pmdp*ppdq - 4*ppdp**2*pmdp )
      B21 = B21 + Dp**(-2)*capT**2 * ( ppdp*pmdp )
      B21 = B21 + Dp**(-2) * (  - M2*v2*ppdp**2 + 2*M2*ppdp*pmdp*ppdq
     &     + 2*M2*ppdp**2*pmdq + M2**2*ppdq*pmdq + 4*ppdp**3*pmdp )
      B21 = B21 + Dp**(-1)*Dm**(-1)*capT * (  - 1./2.*M2*v2*ppdp - 1./2.
     &    *M2*v2*pmdp + 2*M2*ppdp*pmdq + 2*M2*pmdp*ppdq + 4*ppdp*
     &    pmdp**2 + 4*ppdp**2*pmdp )
      B21 = B21 + Dp**(-1)*Dm**(-1)*capT**2 * (  - 2*ppdp*pmdp )
      B21 = B21 + Dp**(-1)*Dm**(-1) * ( 2*M2*v2*ppdp*pmdp - 2*M2*ppdp*
     &    pmdp*ppdq - 2*M2*ppdp*pmdp*pmdq - 2*M2*ppdp**2*pmdq - 2*M2*
     &    pmdp**2*ppdq - 2*M2**2*ppdq*pmdq - 8*ppdp**2*pmdp**2 )
      B21 = B21 + Dm**(-2)*capT * ( 1./2.*M2*v2*pmdp - M2*ppdp*pmdq - 
     &    M2*pmdp*ppdq - 4*ppdp*pmdp**2 )
      B21 = B21 + Dm**(-2)*capT**2 * ( ppdp*pmdp )
      B21 = B21 + Dm**(-2) * (  - M2*v2*pmdp**2 + 2*M2*ppdp*pmdp*pmdq
     &     + 2*M2*pmdp**2*ppdq + M2**2*ppdq*pmdq + 4*ppdp*pmdp**3 )

      B22 =
     &  + t*Dp**(-2) * ( m**2*M2 + 1./2.*M2*v2 - 2*ppdp*pmdp )
      B22 = B22 + t*Dp**(-1)*Dm**(-1) * ( 2*m**2*M2 - M2*v2 )
      B22 = B22 + t*Dm**(-2) * ( m**2*M2 + 1./2.*M2*v2 - 2*ppdp*pmdp )
      B22 = B22 + Dp**(-2)*capT * ( 2*m**2*ppdp + v2*ppdp - 2*ppdp*pmdq
     &     )
      B22 = B22 + Dp**(-2) * (  - 4*m**2*ppdp**2 + 2*M2*ppdq*pmdq - 2*
     &    v2*ppdp**2 - 4*ppdp*pmdp*ppdq + 4*ppdp**2*pmdq )
      B22 = B22 + Dp**(-1)*Dm**(-1)*capT * (  - 2*m**2*ppdp - 2*m**2*
     &    pmdp - v2*ppdp - v2*pmdp )
      B22 = B22 + Dp**(-1)*Dm**(-1)*capT**2 * ( m**2 )
      B22 = B22 + Dp**(-1)*Dm**(-1) * ( 8*m**2*ppdp*pmdp + 4*v2*ppdp*
     &    pmdp - 4*ppdp*pmdp*ppdq - 4*ppdp*pmdp*pmdq + 4*ppdp**2*pmdq
     &     + 4*pmdp**2*ppdq )
      B22 = B22 + Dm**(-2)*capT * ( 2*m**2*pmdp + v2*pmdp - 2*pmdp*ppdq
     &     )
      B22 = B22 + Dm**(-2) * (  - 4*m**2*pmdp**2 + 2*M2*ppdq*pmdq - 2*
     &    v2*pmdp**2 - 4*ppdp*pmdp*pmdq + 4*pmdp**2*ppdq )
     
      B11BSM =
     &  + Dp**(-2)*F2 * ( Sx*pmdp*ppdq - 2*M2*ppdq*pmdbigQ - 3*M2*
     &    ppdq**2 - 3*M2*pmdq*ppdq - M2*capY*ppdq - M2*v2*ppdq - m**2*
     &    Sx*pmdp - m**2*Sx*ppdp + 2*m**2*M2*pmdbigQ + 2*m**2*M2*
     &    ppdbigQ - m**2*M2*capY - m**2*M2*v2 )
      B11BSM = B11BSM + Dp**(-2)*F2*capT * ( 3*pmdp*ppdq + 1./2.*Sx*
     &    ppdq + 1./2.*m**2*Sx )
      B11BSM = B11BSM + Dp**(-1)*Dm**(-1)*F2 * (  - Sx*pmdp*pmdq - Sx*
     &    ppdp*ppdq + v2*Sx*pmdp + v2*Sx*ppdp + 2*M2*ppdq*ppdbigQ + M2*
     &    ppdq**2 + 2*M2*pmdq*pmdbigQ + 2*M2*pmdq*ppdq + M2*pmdq**2 - 2
     &    *M2*v2*pmdbigQ - 2*M2*v2*ppdbigQ + M2*v2*capY + M2*v2**2 - 2*
     &    m**2*Sx*pmdp - 2*m**2*Sx*ppdp + 4*m**2*M2*pmdbigQ + 4*m**2*M2
     &    *ppdbigQ - 2*m**2*M2*capY - 2*m**2*M2*v2 )
      B11BSM = B11BSM + Dp**(-1)*Dm**(-1)*F2*capT * ( 2*pmdp*ppdq - 
     &    pmdp*pmdq - ppdp*ppdq + 2*ppdp*pmdq + 1./2.*capY*pmdp + 1./2.
     &    *capY*ppdp + 1./2.*v2*pmdp + 1./2.*v2*ppdp - 1./2.*v2*Sx + 
     &    m**2*Sx )
      B11BSM = B11BSM + Dp**(-1)*Dm**(-1)*F2*capT**2 * (  - 1./2.*
     &    pmdbigQ - 1./2.*ppdbigQ - 1./2.*v2 )
      B11BSM = B11BSM + Dm**(-2)*F2 * ( Sx*ppdp*pmdq - 2*M2*pmdq*
     &    ppdbigQ - 3*M2*pmdq*ppdq - 3*M2*pmdq**2 - M2*capY*pmdq - M2*
     &    v2*pmdq - m**2*Sx*pmdp - m**2*Sx*ppdp + 2*m**2*M2*pmdbigQ + 2
     &    *m**2*M2*ppdbigQ - m**2*M2*capY - m**2*M2*v2 )
      B11BSM = B11BSM + Dm**(-2)*F2*capT * ( 3*ppdp*pmdq + 1./2.*Sx*
     &    pmdq + 1./2.*m**2*Sx )
      B11BSM = B11BSM + t*Dp**(-2)*F2 * ( 6*ppdp*pmdp + Sx*pmdp + 1./2.
     &    *Sx*ppdp - 2*M2*pmdbigQ - M2*ppdbigQ - 5./2.*M2*ppdq + 3./2.*
     &    M2*pmdq - 3./2.*M2*v2 - m**2*M2 )
      B11BSM = B11BSM + t*Dp**(-2)*F2*capT * (  - 3./2.*pmdp )
      B11BSM = B11BSM + t*Dp**(-1)*Dm**(-1)*F2 * (  - 4*pmdp**2 + 4*
     &    ppdp*pmdp - 4*ppdp**2 - 1./2.*Sx*pmdp - 1./2.*Sx*ppdp + M2*
     &    pmdbigQ + M2*ppdbigQ + 2*M2*v2 - 2*m**2*M2 )
      B11BSM = B11BSM + t*Dp**(-1)*Dm**(-1)*F2*capT * ( 2*pmdp + 2*ppdp
     &     )
      B11BSM = B11BSM + t*Dm**(-2)*F2 * ( 6*ppdp*pmdp + 1./2.*Sx*pmdp
     &     + Sx*ppdp - M2*pmdbigQ - 2*M2*ppdbigQ + 3./2.*M2*ppdq - 5./2.
     &    *M2*pmdq - 3./2.*M2*v2 - m**2*M2 )
      B11BSM = B11BSM + t*Dm**(-2)*F2*capT * (  - 3./2.*ppdp )
      
      B12BSM =
     &  + Dp**(-2)*F2 * (  - 6*ppdq*pmdbigQ - 12*ppdq**2 - 6*pmdq*ppdq
     &     - 3*capY*ppdq - 3*v2*ppdq + 6*m**2*pmdbigQ + 6*m**2*ppdbigQ
     &     - 3*m**2*capY - 3*m**2*v2 )
      B12BSM = B12BSM + Dp**(-1)*Dm**(-1)*F2 * ( 6*ppdq*ppdbigQ + 2*
     &    ppdq**2 + 6*pmdq*pmdbigQ + 16*pmdq*ppdq + 2*pmdq**2 + capY*
     &    ppdq + capY*pmdq - 6*v2*pmdbigQ - 6*v2*ppdbigQ + v2*ppdq + v2
     &    *pmdq + 3*v2*capY + 3*v2**2 + 12*m**2*pmdbigQ + 12*m**2*
     &    ppdbigQ - 6*m**2*capY - 6*m**2*v2 )
      B12BSM = B12BSM + Dm**(-2)*F2 * (  - 6*pmdq*ppdbigQ - 6*pmdq*ppdq
     &     - 12*pmdq**2 - 3*capY*pmdq - 3*v2*pmdq + 6*m**2*pmdbigQ + 6*
     &    m**2*ppdbigQ - 3*m**2*capY - 3*m**2*v2 )
      B12BSM = B12BSM + t*Dp**(-2)*F2 * (  - 6*pmdbigQ - 3*ppdbigQ - 9*
     &    ppdq + 3*pmdq - 3*v2 - 9*m**2 )
      B12BSM = B12BSM + t*Dp**(-1)*Dm**(-1)*F2 * ( 5*pmdbigQ + 5*
     &    ppdbigQ + 4*ppdq + 4*pmdq + 11*v2 - 18*m**2 )
      B12BSM = B12BSM + t*Dm**(-2)*F2 * (  - 3*pmdbigQ - 6*ppdbigQ + 3*
     &    ppdq - 9*pmdq - 3*v2 - 9*m**2 )


      B21BSM =
     &  + Dp**(-2)*F2 * ( Sx*ppdp**2*pmdp + Sx*ppdp**3 - 2*M2*ppdp**2*
     &    pmdbigQ - 2*M2*ppdp**2*ppdbigQ + M2*capY*ppdp**2 + M2*v2*
     &    ppdp**2 + M2**2*ppdq**2 + M2**2*pmdq*ppdq + 1./2.*M2**2*capY*
     &    ppdq + 1./2.*M2**2*v2*ppdq )
      B21BSM = B21BSM + Dp**(-2)*F2*capT * (  - 1./2.*Sx*ppdp*pmdp - Sx
     &    *ppdp**2 - 1./2.*M2*pmdp*ppdq + M2*ppdp*pmdbigQ + M2*ppdp*
     &    ppdbigQ - 1./2.*M2*ppdp*ppdq - 1./2.*M2*capY*ppdp - 1./4.*M2*
     &    Sx*ppdq - 1./2.*M2*v2*ppdp )
      B21BSM = B21BSM + Dp**(-2)*F2*capT**2 * (  - 1./2.*ppdp*pmdp + 1./
     &    4.*Sx*ppdp - 1./4.*M2*pmdq + 1./8.*M2*v2 )
      B21BSM = B21BSM + Dp**(-2)*F2*capT**3 * ( 1./4.*pmdp )
      B21BSM = B21BSM + Dp**(-1)*Dm**(-1)*F2 * (  - 2*Sx*ppdp*pmdp**2
     &     - 2*Sx*ppdp**2*pmdp + 4*M2*ppdp*pmdp*pmdbigQ + 4*M2*ppdp*
     &    pmdp*ppdbigQ - 2*M2*capY*ppdp*pmdp - 2*M2*v2*ppdp*pmdp - 
     &    M2**2*ppdq**2 - 2*M2**2*pmdq*ppdq - M2**2*pmdq**2 - 1./2.*
     &    M2**2*capY*ppdq - 1./2.*M2**2*capY*pmdq - 1./2.*M2**2*v2*ppdq
     &     - 1./2.*M2**2*v2*pmdq )
      B21BSM = B21BSM + Dp**(-1)*Dm**(-1)*F2*capT * ( 1./2.*Sx*pmdp**2
     &     + 2*Sx*ppdp*pmdp + 1./2.*Sx*ppdp**2 - M2*pmdp*pmdbigQ - M2*
     &    pmdp*ppdbigQ + 1./2.*M2*pmdp*ppdq + 1./2.*M2*pmdp*pmdq - M2*
     &    ppdp*pmdbigQ - M2*ppdp*ppdbigQ + 1./2.*M2*ppdp*ppdq + 1./2.*
     &    M2*ppdp*pmdq + 1./2.*M2*capY*pmdp + 1./2.*M2*capY*ppdp + 1./4.
     &    *M2*Sx*ppdq + 1./4.*M2*Sx*pmdq + 1./2.*M2*v2*pmdp + 1./2.*M2*
     &    v2*ppdp )
      B21BSM = B21BSM + Dp**(-1)*Dm**(-1)*F2*capT**2 * ( 1./2.*pmdp**2
     &     + 1./2.*ppdp**2 - 1./4.*Sx*pmdp - 1./4.*Sx*ppdp + 1./4.*M2*
     &    ppdq + 1./4.*M2*pmdq - 1./4.*M2*v2 )
      B21BSM = B21BSM + Dp**(-1)*Dm**(-1)*F2*capT**3 * (  - 1./4.*pmdp
     &     - 1./4.*ppdp )
      B21BSM = B21BSM + Dm**(-2)*F2 * ( Sx*pmdp**3 + Sx*ppdp*pmdp**2 - 
     &    2*M2*pmdp**2*pmdbigQ - 2*M2*pmdp**2*ppdbigQ + M2*capY*pmdp**2
     &     + M2*v2*pmdp**2 + M2**2*pmdq*ppdq + M2**2*pmdq**2 + 1./2.*
     &    M2**2*capY*pmdq + 1./2.*M2**2*v2*pmdq )
      B21BSM = B21BSM + Dm**(-2)*F2*capT * (  - Sx*pmdp**2 - 1./2.*Sx*
     &    ppdp*pmdp + M2*pmdp*pmdbigQ + M2*pmdp*ppdbigQ - 1./2.*M2*pmdp
     &    *pmdq - 1./2.*M2*ppdp*pmdq - 1./2.*M2*capY*pmdp - 1./4.*M2*Sx
     &    *pmdq - 1./2.*M2*v2*pmdp )
      B21BSM = B21BSM + Dm**(-2)*F2*capT**2 * (  - 1./2.*ppdp*pmdp + 1./
     &    4.*Sx*pmdp - 1./4.*M2*ppdq + 1./8.*M2*v2 )
      B21BSM = B21BSM + Dm**(-2)*F2*capT**3 * ( 1./4.*ppdp )
      B21BSM = B21BSM + t*Dp**(-2)*F2 * (  - 2*M2*ppdp*pmdp + M2*
     &    ppdp**2 - 1./4.*M2*Sx*pmdp - 1./4.*M2*Sx*ppdp + 1./2.*M2**2*
     &    pmdbigQ + 1./2.*M2**2*ppdbigQ + M2**2*ppdq - 1./2.*M2**2*pmdq
     &     + 1./2.*M2**2*v2 )
      B21BSM = B21BSM + t*Dp**(-2)*F2*capT * ( 3./4.*M2*pmdp - 3./4.*M2
     &    *ppdp )
      B21BSM = B21BSM + t*Dp**(-1)*Dm**(-1)*F2 * ( 2*M2*pmdp**2 - 2*M2*
     &    ppdp*pmdp + 2*M2*ppdp**2 + 1./2.*M2*Sx*pmdp + 1./2.*M2*Sx*
     &    ppdp - M2**2*pmdbigQ - M2**2*ppdbigQ - 1./2.*M2**2*ppdq - 1./
     &    2.*M2**2*pmdq - M2**2*v2 )
      B21BSM = B21BSM + t*Dm**(-2)*F2 * ( M2*pmdp**2 - 2*M2*ppdp*pmdp
     &     - 1./4.*M2*Sx*pmdp - 1./4.*M2*Sx*ppdp + 1./2.*M2**2*pmdbigQ
     &     + 1./2.*M2**2*ppdbigQ - 1./2.*M2**2*ppdq + M2**2*pmdq + 1./2.
     &    *M2**2*v2 )
      B21BSM = B21BSM + t*Dm**(-2)*F2*capT * (  - 3./4.*M2*pmdp + 3./4.
     &    *M2*ppdp )

      B22BSM =
     &  + Dp**(-2)*F2 * (  - 6*ppdp**2*pmdQ - 6*ppdp**2*ppdQ + 3*capY*
     &    ppdp**2 - 3*Sx*ppdp*ppdq + 3*v2*ppdp**2 + 4*M2*ppdq**2 + 2*M2
     &    *pmdq*ppdq + 3./2.*M2*capY*ppdq + 3./2.*M2*v2*ppdq )
      B22BSM = B22BSM + Dp**(-2)*F2*capT * (  - pmdp*ppdq + 3*ppdp*pmdQ
     &     + 3*ppdp*ppdQ - 2*ppdp*ppdq - 3./2.*capY*ppdp - 3./2.*v2*
     &    ppdp )
      B22BSM = B22BSM + Dp**(-2)*F2*capT**2 * (  - 1./2.*pmdq + 1./4.*
     &    v2 + 1./2.*m**2 )
      B22BSM = B22BSM + Dp**(-1)*Dm**(-1)*F2 * ( 12*ppdp*pmdp*pmdQ + 12
     &    *ppdp*pmdp*ppdQ - 6*capY*ppdp*pmdp + 3*Sx*pmdp*ppdq + 3*Sx*
     &    ppdp*pmdq - 6*v2*ppdp*pmdp - 4*M2*pmdq*ppdq - 1./2.*M2*capY*
     &    ppdq - 1./2.*M2*capY*pmdq - 1./2.*M2*v2*ppdq - 1./2.*M2*v2*
     &    pmdq )
      B22BSM = B22BSM + Dp**(-1)*Dm**(-1)*F2*capT * (  - 3*pmdp*pmdQ - 
     &    3*pmdp*ppdQ - pmdp*pmdq - 3*ppdp*pmdQ - 3*ppdp*ppdQ - ppdp*
     &    ppdq + capY*pmdp + capY*ppdp - 1./2.*Sx*ppdq - 1./2.*Sx*pmdq
     &     + v2*pmdp + v2*ppdp )
      B22BSM = B22BSM + Dp**(-1)*Dm**(-1)*F2*capT**2 * ( 1./2.*pmdQ + 1.
     &    /2.*ppdQ + 1./2.*ppdq + 1./2.*pmdq - 1./2.*v2 + m**2 )
      B22BSM = B22BSM + Dm**(-2)*F2 * (  - 6*pmdp**2*pmdQ - 6*pmdp**2*
     &    ppdQ + 3*capY*pmdp**2 - 3*Sx*pmdp*pmdq + 3*v2*pmdp**2 + 2*M2*
     &    pmdq*ppdq + 4*M2*pmdq**2 + 3./2.*M2*capY*pmdq + 3./2.*M2*v2*
     &    pmdq )
      B22BSM = B22BSM + Dm**(-2)*F2*capT * ( 3*pmdp*pmdQ + 3*pmdp*ppdQ
     &     - 2*pmdp*pmdq - ppdp*pmdq - 3./2.*capY*pmdp - 3./2.*v2*pmdp
     &     )
      B22BSM = B22BSM + Dm**(-2)*F2*capT**2 * (  - 1./2.*ppdq + 1./4.*
     &    v2 + 1./2.*m**2 )
      B22BSM = B22BSM + t*Dp**(-2)*F2 * ( 3*ppdp**2 - 3./2.*Sx*ppdp + 3.
     &    /2.*M2*pmdQ + 3./2.*M2*ppdQ + 7./2.*M2*ppdq - M2*pmdq + M2*v2
     &     + 2*m**2*M2 )
      B22BSM = B22BSM + t*Dp**(-2)*F2*capT * (  - 1./2.*pmdp - 5./2.*
     &    ppdp )
      B22BSM = B22BSM + t*Dp**(-1)*Dm**(-1)*F2 * (  - 2*pmdp**2 - 10*
     &    ppdp*pmdp - 2*ppdp**2 + 1./2.*Sx*pmdp + 1./2.*Sx*ppdp - M2*
     &    pmdQ - M2*ppdQ - 3./2.*M2*ppdq - 3./2.*M2*pmdq - 2*M2*v2 + 4*
     &    m**2*M2 )
      B22BSM = B22BSM + t*Dp**(-1)*Dm**(-1)*F2*capT * ( 5./2.*pmdp + 5./
     &    2.*ppdp )
      B22BSM = B22BSM + t*Dm**(-2)*F2 * ( 3*pmdp**2 - 3./2.*Sx*pmdp + 3.
     &    /2.*M2*pmdQ + 3./2.*M2*ppdQ - M2*ppdq + 7./2.*M2*pmdq + M2*v2
     &     + 2*m**2*M2 )
      B22BSM = B22BSM + t*Dm**(-2)*F2*capT * (  - 5./2.*pmdp - 1./2.*
     &    ppdp )


 
      endif
      return
      end
