c*     june 2025
c* mu(k)+p(p)-> mu(k1)+tau-(p-)+tau+(p+)+p(p2)
c* q = p-p2, k+q -> kap + k1
c* Dkq = (k+q)^2-mmu^2
c* Dk1q = (k1-q)^2-mmu^2
c* t=-q^2
c* u = p/M
c     * s = (k+p)^2
c     m2 is pair mass, Mp2 is target mass (both squared)
c             M11         Terms in output =         43
c             M21         Terms in output =         20
c             M12         Terms in output =        145
c             M22         Terms in output =         68


      subroutine projections(m2,Mp2,ein,v2,t,Y,k1dpp,kdpp,
     x       kdq,k1dq,
     x       ppdp,X,B11,B12,B21,B22)
      implicit none
      real*8 M11,M12,M21,M22,B11,B12,B21,B22,m2,mmu2,Mp2,v2,
     x     t,Y,k1dpp,kdpp,kdq,k1dq,Dkq,Dk1q,W1,W2,s0,ein,X,
     x     ppdp

      W1 = 1.
      W2 = 1.
      mmu2 = (0.10566)**2
C      Mp2 = 0.94**2 this is input in projections 
      
      s0 = 2.*dsqrt(Mp2)*ein +mmu2+Mp2
      
      Dkq = 2.*kdq - t
      Dk1q = -2.*k1dq - t

      M11 = 0.
      M12 = 0.
      M21 = 0.
      M22 = 0.

      goto 100

 100  continue
     
       M11 =
     &  + W1*Dkq**(-2) * (  - 4*t*k1dpp**2 - 2*v2*t*k1dpp + 8*mmu2*
     &    k1dpp**2 + 4*mmu2*v2*k1dpp + 2*m2*t*k1dq - 4*m2*mmu2*k1dq + 
     &    m2*Y*t - 2*m2*Y*mmu2 )
      M11 = M11 + W1*Dkq**(-1)*Dk1q**(-1) * (  - 4*t*k1dpp**2 - 8*t*
     &    kdpp*k1dpp - 4*t*kdpp**2 + 16*mmu2*kdpp*k1dpp - 4*mmu2*v2*
     &    k1dpp + 4*mmu2*v2*kdpp - 2*mmu2*v2**2 + 4*Y*k1dpp**2 + 4*Y*
     &    kdpp**2 + 2*Y*v2*k1dpp - 2*Y*v2*kdpp + 4*m2*mmu2*t + 2*m2*Y*t
     &     - 4*m2*Y*mmu2 - 2*m2*Y**2 )
      M11 = M11 + W1*Dkq**(-1) * (  - 4*k1dpp**2 - 4*kdpp**2 - 2*v2*
     &    k1dpp + 2*m2*k1dq + 6*m2*mmu2 + 2*m2*Y )
      M11 = M11 + W1*Dk1q**(-2) * (  - 4*t*kdpp**2 + 2*v2*t*kdpp + 8*
     &    mmu2*kdpp**2 - 4*mmu2*v2*kdpp - 2*m2*t*kdq + 4*m2*mmu2*kdq + 
     &    m2*Y*t - 2*m2*Y*mmu2 )
      M11 = M11 + W1*Dk1q**(-1) * (  - 4*k1dpp**2 - 4*kdpp**2 + 2*v2*
     &    kdpp - 2*m2*kdq + 6*m2*mmu2 + 2*m2*Y )



 200  continue


      M21 =
     &  + W1*Dkq**(-2) * ( v2*t*k1dq - 2*mmu2*v2*k1dq - mmu2*v2*t + 2*
     &    mmu2**2*v2 + 1./2.*Y*v2*t - Y*mmu2*v2 )
      M21 = M21 + W1*Dkq**(-1)*Dk1q**(-1) * ( 4*mmu2**2*v2 - Y**2*v2 )
      M21 = M21 + W1*Dkq**(-1) * ( v2*k1dq + 3*mmu2*v2 + Y*v2 )
      M21 = M21 + W1*Dk1q**(-2) * (  - v2*t*kdq + 2*mmu2*v2*kdq - mmu2*
     &    v2*t + 2*mmu2**2*v2 + 1./2.*Y*v2*t - Y*mmu2*v2 )
      M21 = M21 + W1*Dk1q**(-1) * (  - v2*kdq + 3*mmu2*v2 + Y*v2 )
     



 300  continue
C     problem in M12


      M12 =
     &  + Mp2**(-1)*W2*Dkq**(-2) * (  - 2*k1dpp**2*s0**2 + 2*t*k1dpp**2
     &    *s0 - v2*k1dpp*s0**2 + v2*t*k1dpp*s0 + 4*mmu2*k1dpp**2*s0 - 2
     &    *mmu2*t*k1dpp**2 + 2*mmu2*v2*k1dpp*s0 - mmu2*v2*t*k1dpp - 2*
     &    mmu2**2*k1dpp**2 - mmu2**2*v2*k1dpp + m2*k1dq*s0**2 - m2*t*
     &    k1dq*s0 - 2*m2*mmu2*k1dq*s0 + m2*mmu2*t*k1dq + m2*mmu2**2*
     &    k1dq + 1./2.*m2*Y*s0**2 - 1./2.*m2*Y*t*s0 - m2*Y*mmu2*s0 + 1./
     &    2.*m2*Y*mmu2*t + 1./2.*m2*Y*mmu2**2 )
      M12 = M12 + Mp2**(-1)*W2*Dkq**(-1)*Dk1q**(-1) * (  - 2*t*kdpp*
     &    k1dpp*s0 - 2*t**2*ppdp**2 + 2*t**2*k1dpp*ppdp - 2*t**2*kdpp*
     &    ppdp + 2*t**2*kdpp*k1dpp + v2*t*ppdp*s0 + 2*mmu2*t*kdpp*k1dpp
     &     - mmu2*v2*t*ppdp - 4*X*kdpp*k1dpp*s0 + 2*X*t*kdpp*k1dpp + X*
     &    v2*k1dpp*s0 - X*v2*kdpp*s0 - X*v2*t*ppdp + 1./2.*X*v2**2*s0
     &     + 4*X*mmu2*kdpp*k1dpp - X*mmu2*v2*k1dpp + X*mmu2*v2*kdpp - 1.
     &    /2.*X*mmu2*v2**2 + 2*Y*t*ppdp**2 - 2*Y*t*k1dpp*ppdp + 2*Y*t*
     &    kdpp*ppdp - Y*v2*t*ppdp - 1./2.*m2*t*s0**2 + 1./2.*m2*t**2*s0
     &     + m2*mmu2*t*s0 - 1./2.*m2*mmu2*t**2 - 1./2.*m2*mmu2**2*t - 1.
     &    /2.*m2*X*t**2 - 1./2.*m2*X**2*t + 1./2.*m2*Y*t*s0 - 1./2.*m2*
     &    Y*t**2 - 1./2.*m2*Y*mmu2*t + m2*Y*X*s0 - 1./2.*m2*Y*X*t - m2*
     &    Y*X*mmu2 )
      M12 = M12 + Mp2**(-1)*W2*Dkq**(-1) * ( 4*k1dpp*ppdp*s0 - 2*t*
     &    ppdp**2 - 2*t*kdpp*ppdp + v2*ppdp*s0 - 4*mmu2*k1dpp*ppdp - 
     &    mmu2*v2*ppdp - 1./2.*m2*s0**2 + 1./2.*m2*t*s0 + m2*mmu2*s0 - 
     &    1./2.*m2*mmu2*t - 1./2.*m2*mmu2**2 - m2*X*s0 + m2*X*mmu2 )
      M12 = M12 + Mp2**(-1)*W2*Dk1q**(-2) * (  - 2*X*t*kdpp**2 + X*v2*t
     &    *kdpp - 2*X**2*kdpp**2 + X**2*v2*kdpp - m2*X*t*kdq - m2*X**2*
     &    kdq + 1./2.*m2*Y*X*t + 1./2.*m2*Y*X**2 )
      M12 = M12 + Mp2**(-1)*W2*Dk1q**(-1) * (  - 2*t*ppdp**2 + 2*t*
     &    k1dpp*ppdp + 4*X*kdpp*ppdp - X*v2*ppdp - m2*X*s0 - 1./2.*m2*X
     &    *t + m2*X*mmu2 - 1./2.*m2*X**2 )
      M12 = M12 + Mp2**(-1)*W2 * (  - 2*ppdp**2 )
      M12 = M12 + W2*Dkq**(-2) * ( 4*k1dpp**2*s0 + 2*v2*k1dpp*s0 - 4*
     &    mmu2*k1dpp**2 - 2*mmu2*v2*k1dpp - 2*m2*k1dq*s0 + 2*m2*mmu2*
     &    k1dq - m2*Y*s0 + m2*Y*mmu2 )
      M12 = M12 + W2*Dkq**(-1)*Dk1q**(-1) * ( 2*t*k1dpp**2 + 2*t*kdpp*
     &    k1dpp + 2*t*kdpp**2 - v2*t*ppdp + v2*t*k1dpp - v2*t*kdpp + 4*
     &    X*kdpp*k1dpp - X*v2*k1dpp + X*v2*kdpp - 1./2.*X*v2**2 - 2*Y*
     &    k1dpp**2 + 4*Y*kdpp*k1dpp - 2*Y*kdpp**2 - 2*Y*v2*k1dpp + 2*Y*
     &    v2*kdpp - 1./2.*Y*v2**2 + m2*t*s0 + 1./2.*m2*t**2 - m2*mmu2*t
     &     - 3./2.*m2*Y*t - m2*Y*X )
      M12 = M12 + W2*Dkq**(-1) * (  - 4*k1dpp*ppdp + 2*k1dpp**2 - 4*
     &    kdpp*k1dpp + 2*kdpp**2 - v2*ppdp + v2*k1dpp - v2*kdpp + m2*s0
     &     - m2*k1dq + 1./2.*m2*t - m2*mmu2 + m2*X )
      M12 = M12 + W2*Dk1q**(-2) * ( 2*t*kdpp**2 - v2*t*kdpp + m2*t*kdq
     &     - 1./2.*m2*Y*t )
      M12 = M12 + W2*Dk1q**(-1) * ( 2*k1dpp**2 - 4*kdpp*k1dpp + 2*
     &    kdpp**2 + v2*k1dpp - v2*kdpp + m2*kdq + m2*t + m2*X )
      M12 = M12 + W2 * ( m2 )
      M12 = M12 + Mp2*W2*Dkq**(-2) * (  - 2*k1dpp**2 - v2*k1dpp + m2*
     &    k1dq + 1./2.*m2*Y )
      M12 = M12 + Mp2*W2*Dkq**(-1)*Dk1q**(-1) * (  - 1./2.*m2*t )
      M12 = M12 + Mp2*W2*Dkq**(-1) * (  - 1./2.*m2 )

      

 400  continue
C     problem in M22
      
      
      M22 =
     &  + Mp2**(-1)*W2*Dkq**(-2) * ( 1./2.*v2*k1dq*s0**2 - 1./2.*v2*t*
     &    k1dq*s0 - 1./2.*mmu2*v2*s0**2 - mmu2*v2*k1dq*s0 + 1./2.*mmu2*
     &    v2*t*s0 + 1./2.*mmu2*v2*t*k1dq + mmu2**2*v2*s0 + 1./2.*
     &    mmu2**2*v2*k1dq - 1./2.*mmu2**2*v2*t - 1./2.*mmu2**3*v2 + 1./
     &    4.*Y*v2*s0**2 - 1./4.*Y*v2*t*s0 - 1./2.*Y*mmu2*v2*s0 + 1./4.*
     &    Y*mmu2*v2*t + 1./4.*Y*mmu2**2*v2 )
      M22 = M22 + Mp2**(-1)*W2*Dkq**(-1)*Dk1q**(-1) * (  - 1./4.*v2*t*
     &    s0**2 + 1./2.*mmu2*v2*t**2 + 1./4.*mmu2**2*v2*t - 1./2.*X*v2*
     &    t*s0 - X*mmu2*v2*s0 + X*mmu2*v2*t + X*mmu2**2*v2 - 1./4.*X**2
     &    *v2*t + 1./4.*Y*v2*t*s0 - 1./4.*Y*mmu2*v2*t + 1./2.*Y*X*v2*s0
     &     - 1./4.*Y*X*v2*t - 1./2.*Y*X*mmu2*v2 )
      M22 = M22 + Mp2**(-1)*W2*Dkq**(-1) * (  - 1./4.*v2*s0**2 + 1./2.*
     &    mmu2*v2*s0 - 1./4.*mmu2**2*v2 - 1./2.*X*v2*s0 + 1./2.*X*mmu2*
     &    v2 )
      M22 = M22 + Mp2**(-1)*W2*Dk1q**(-2) * (  - 1./2.*X*v2*t*kdq - 1./
     &    2.*X*mmu2*v2*t - 1./2.*X**2*v2*kdq - 1./2.*X**2*mmu2*v2 + 1./
     &    4.*Y*X*v2*t + 1./4.*Y*X**2*v2 )
      M22 = M22 + Mp2**(-1)*W2*Dk1q**(-1) * (  - 1./2.*X*v2*s0 + 1./2.*
     &    X*mmu2*v2 - 1./4.*X**2*v2 )
      M22 = M22 + W2*Dkq**(-2) * (  - v2*k1dq*s0 + mmu2*v2*s0 + mmu2*v2
     &    *k1dq - mmu2**2*v2 - 1./2.*Y*v2*s0 + 1./2.*Y*mmu2*v2 )
      M22 = M22 + W2*Dkq**(-1)*Dk1q**(-1) * ( 1./2.*v2*t*s0 + mmu2*v2*t
     &     + 1./2.*X*v2*t + X*mmu2*v2 + 1./4.*Y*v2*t - 1./2.*Y*X*v2 )
      M22 = M22 + W2*Dkq**(-1) * ( 1./2.*v2*s0 - 1./2.*v2*k1dq - 1./2.*
     &    mmu2*v2 + 1./2.*X*v2 )
      M22 = M22 + W2*Dk1q**(-2) * ( 1./2.*v2*t*kdq + 1./2.*mmu2*v2*t - 
     &    1./4.*Y*v2*t )
      M22 = M22 + W2*Dk1q**(-1) * ( 1./2.*v2*kdq + 1./2.*X*v2 )
      M22 = M22 + Mp2*W2*Dkq**(-2) * ( 1./2.*v2*k1dq - 1./2.*mmu2*v2 + 
     &    1./4.*Y*v2 )
      M22 = M22 + Mp2*W2*Dkq**(-1)*Dk1q**(-1) * (  - 1./4.*v2*t )
      M22 = M22 + Mp2*W2*Dkq**(-1) * (  - 1./4.*v2 )


 500  continue
c     KK = L/4, M = B/4 so multiply everything times 16 for total contractions
C form factor definitions will include the spin average of the proton      
      B11 = M11*4.*4.
      B12 = M12*4.*4.
      B21 = M21*4.*4.
      B22 = M22*4.*4.

      return
      end 
