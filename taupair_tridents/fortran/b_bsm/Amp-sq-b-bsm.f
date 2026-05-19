c*     june 2025
c* mu(k)+p(p)-> mu(k1)+tau-(p-)+tau+(p+)+p(p2)
c* q = p-p2, k+q -> kap + k1
c* Dkq = (k+q)^2-mmu^2
c* Dk1q = (k1-q)^2-mmu^2
c* t=-q^2
c* u = p/M
c     * s = (k+p)^2
c     m2 is pair mass, Mp2 is target mass (both squared)
c             M11         Terms in output =         54
c             M22         Terms in output =         154


      subroutine projections(m2,Mp2,ein,v2,t,Y,k1dpp,kdpp,
     x       kdq,k1dq,
     x       ppdp,X,B1,B2)
      implicit none
      real*8 M11,M12,M21,M22,B11,B12,B21,B22,m2,mmu2,Mp2,v2,
     x     B1,B2,
     x     t,Y,k1dpp,kdpp,kdq,k1dq,Dkq,Dk1q,W1,W2,s0,ein,X,
     x     ppdp

      W1 = 1.
      W2 = 1.
      mmu2 = (0.106d0)**2
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
     &  + W1*Dkq**(-2) * (  - 8*m2*t*k1dq + 16*m2*mmu2*k1dq - 16*m2*
     &    mmu2*t + 32*m2*mmu2**2 - 4*m2*Y*t + 8*m2*Y*mmu2 )
      M11 = M11 + W1*Dkq**(-2)*v2 * ( 2*t*k1dq - 4*mmu2*k1dq + 4*mmu2*t
     &     - 8*mmu2**2 + Y*t - 2*Y*mmu2 )
      M11 = M11 + W1*Dkq**(-1)*Dk1q**(-1) * ( 8*m2*t**2 - 80*m2*mmu2*t
     &     + 64*m2*mmu2**2 - 16*m2*Y*t + 48*m2*Y*mmu2 + 8*m2*Y**2 )
      M11 = M11 + W1*Dkq**(-1)*Dk1q**(-1)*v2 * (  - 2*t**2 + 20*mmu2*t
     &     - 16*mmu2**2 + 4*Y*t - 12*Y*mmu2 - 2*Y**2 )
      M11 = M11 + W1*Dkq**(-1) * (  - 8*m2*k1dq + 8*m2*t - 24*m2*mmu2
     &     - 8*m2*Y )
      M11 = M11 + W1*Dkq**(-1)*v2 * ( 2*k1dq - 2*t + 6*mmu2 + 2*Y )
      M11 = M11 + W1*Dk1q**(-2) * ( 8*m2*t*kdq - 16*m2*mmu2*kdq - 16*m2
     &    *mmu2*t + 32*m2*mmu2**2 - 4*m2*Y*t + 8*m2*Y*mmu2 )
      M11 = M11 + W1*Dk1q**(-2)*v2 * (  - 2*t*kdq + 4*mmu2*kdq + 4*mmu2
     &    *t - 8*mmu2**2 + Y*t - 2*Y*mmu2 )
      M11 = M11 + W1*Dk1q**(-1) * ( 8*m2*kdq + 8*m2*t - 24*m2*mmu2 - 8*
     &    m2*Y )
      M11 = M11 + W1*Dk1q**(-1)*v2 * (  - 2*kdq - 2*t + 6*mmu2 + 2*Y )
      M11 = M11 + W1 * ( 8*m2 )
      M11 = M11 + W1*v2 * (  - 2 )
      
c 200  continue


       M22 =
     &  + Mp2**(-1)*W2*Dkq**(-2) * (  - 4*m2*k1dq*s0**2 + 4*m2*t*k1dq*
     &    s0 - 8*m2*mmu2*s0**2 + 8*m2*mmu2*k1dq*s0 + 8*m2*mmu2*t*s0 - 4
     &    *m2*mmu2*t*k1dq + 16*m2*mmu2**2*s0 - 4*m2*mmu2**2*k1dq - 8*m2
     &    *mmu2**2*t - 8*m2*mmu2**3 - 2*m2*Y*s0**2 + 2*m2*Y*t*s0 + 4*m2
     &    *Y*mmu2*s0 - 2*m2*Y*mmu2*t - 2*m2*Y*mmu2**2 )
      M22 = M22 + Mp2**(-1)*W2*Dkq**(-2)*v2 * ( k1dq*s0**2 - t*k1dq*s0
     &     + 2*mmu2*s0**2 - 2*mmu2*k1dq*s0 - 2*mmu2*t*s0 + mmu2*t*k1dq
     &     - 4*mmu2**2*s0 + mmu2**2*k1dq + 2*mmu2**2*t + 2*mmu2**3 + 1./
     &    2.*Y*s0**2 - 1./2.*Y*t*s0 - Y*mmu2*s0 + 1./2.*Y*mmu2*t + 1./2.
     &    *Y*mmu2**2 )
      M22 = M22 + Mp2**(-1)*W2*Dkq**(-1)*Dk1q**(-1) * ( 2*m2*t*s0**2 - 
     &    2*m2*t**2*s0 - 12*m2*mmu2*t*s0 + 10*m2*mmu2*t**2 + 10*m2*
     &    mmu2**2*t + 2*m2*X*t**2 - 16*m2*X*mmu2*s0 + 8*m2*X*mmu2*t + 
     &    16*m2*X*mmu2**2 + 2*m2*X**2*t - 2*m2*Y*t*s0 + 2*m2*Y*t**2 + 2
     &    *m2*Y*mmu2*t - 4*m2*Y*X*s0 + 2*m2*Y*X*t + 4*m2*Y*X*mmu2 )
      M22 = M22 + Mp2**(-1)*W2*Dkq**(-1)*Dk1q**(-1)*v2 * (  - 1./2.*t*
     &    s0**2 + 1./2.*t**2*s0 + 3*mmu2*t*s0 - 5./2.*mmu2*t**2 - 5./2.
     &    *mmu2**2*t - 1./2.*X*t**2 + 4*X*mmu2*s0 - 2*X*mmu2*t - 4*X*
     &    mmu2**2 - 1./2.*X**2*t + 1./2.*Y*t*s0 - 1./2.*Y*t**2 - 1./2.*
     &    Y*mmu2*t + Y*X*s0 - 1./2.*Y*X*t - Y*X*mmu2 )
      M22 = M22 + Mp2**(-1)*W2*Dkq**(-1) * ( 2*m2*s0**2 - 2*m2*t*s0 - 4
     &    *m2*mmu2*s0 + 2*m2*mmu2*t + 2*m2*mmu2**2 + 4*m2*X*s0 - 4*m2*X
     &    *mmu2 )
      M22 = M22 + Mp2**(-1)*W2*Dkq**(-1)*v2 * (  - 1./2.*s0**2 + 1./2.*
     &    t*s0 + mmu2*s0 - 1./2.*mmu2*t - 1./2.*mmu2**2 - X*s0 + X*mmu2
     &     )
      M22 = M22 + Mp2**(-1)*W2*Dk1q**(-2) * ( 4*m2*X*t*kdq - 8*m2*X*
     &    mmu2*t + 4*m2*X**2*kdq - 8*m2*X**2*mmu2 - 2*m2*Y*X*t - 2*m2*Y
     &    *X**2 )
      M22 = M22 + Mp2**(-1)*W2*Dk1q**(-2)*v2 * (  - X*t*kdq + 2*X*mmu2*
     &    t - X**2*kdq + 2*X**2*mmu2 + 1./2.*Y*X*t + 1./2.*Y*X**2 )
      M22 = M22 + Mp2**(-1)*W2*Dk1q**(-1) * ( 4*m2*X*s0 + 2*m2*X*t - 4*
     &    m2*X*mmu2 + 2*m2*X**2 )
      M22 = M22 + Mp2**(-1)*W2*Dk1q**(-1)*v2 * (  - X*s0 - 1./2.*X*t + 
     &    X*mmu2 - 1./2.*X**2 )
      M22 = M22 + W2*Dkq**(-2) * ( 8*m2*k1dq*s0 + 16*m2*mmu2*s0 - 8*m2*
     &    mmu2*k1dq - 16*m2*mmu2**2 + 4*m2*Y*s0 - 4*m2*Y*mmu2 )
      M22 = M22 + W2*Dkq**(-2)*v2 * (  - 2*k1dq*s0 - 4*mmu2*s0 + 2*mmu2
     &    *k1dq + 4*mmu2**2 - Y*s0 + Y*mmu2 )
      M22 = M22 + W2*Dkq**(-1)*Dk1q**(-1) * (  - 4*m2*t*s0 - 2*m2*t**2
     &     + 28*m2*mmu2*t + 16*m2*X*mmu2 + 6*m2*Y*t + 4*m2*Y*X )
      M22 = M22 + W2*Dkq**(-1)*Dk1q**(-1)*v2 * ( t*s0 + 1./2.*t**2 - 7*
     &    mmu2*t - 4*X*mmu2 - 3./2.*Y*t - Y*X )
      M22 = M22 + W2*Dkq**(-1) * (  - 4*m2*s0 + 4*m2*k1dq - 2*m2*t + 4*
     &    m2*mmu2 - 4*m2*X )
      M22 = M22 + W2*Dkq**(-1)*v2 * ( s0 - k1dq + 1./2.*t - mmu2 + X )
      M22 = M22 + W2*Dk1q**(-2) * (  - 4*m2*t*kdq + 8*m2*mmu2*t + 2*m2*
     &    Y*t )
      M22 = M22 + W2*Dk1q**(-2)*v2 * ( t*kdq - 2*mmu2*t - 1./2.*Y*t )
      M22 = M22 + W2*Dk1q**(-1) * (  - 4*m2*kdq - 4*m2*t - 4*m2*X )
      M22 = M22 + W2*Dk1q**(-1)*v2 * ( kdq + t + X )
      M22 = M22 + W2 * (  - 4*m2 )
      M22 = M22 + W2*v2 * ( 1 )
      M22 = M22 + Mp2*W2*Dkq**(-2) * (  - 4*m2*k1dq - 8*m2*mmu2 - 2*m2*
     &    Y )
      M22 = M22 + Mp2*W2*Dkq**(-2)*v2 * ( k1dq + 2*mmu2 + 1./2.*Y )
      M22 = M22 + Mp2*W2*Dkq**(-1)*Dk1q**(-1) * ( 2*m2*t )
      M22 = M22 + Mp2*W2*Dkq**(-1)*Dk1q**(-1)*v2 * (  - 1./2.*t )
      M22 = M22 + Mp2*W2*Dkq**(-1) * ( 2*m2 )
      M22 = M22 + Mp2*W2*Dkq**(-1)*v2 * (  - 1./2. )

 300  continue
 500  continue
c     KK = L/4, M = B/4 so multiply everything times 16 for total contractions
C form factor definitions will include the spin average of the proton      
      B1 = M11*4.*1.
      B2 = M22*4.*1.
      
      return
      end 
