c this is in the v2 rest frame
      subroutine v2rest(y,sx,v2,t,cthe,phie,phiq,
     x k1dpp,kdpp,kdq,k1dq,ppdp)
      implicit real*8 (a-z)

      dimension kappavec(4),pplab(4),pmlab(4)
c      common/inputs/emu
      common/inputs/emu,Gamma
      common/masses/mtau,mmu,mp,dm2,mmu2,mp2,mx
c      common/scalar/mx2,gams,gmu2,gtau2


c      print *, y,x,v2,mx2,cthe,phie,phiq

c     rest frame momenta      

      v = dsqrt(v2)

      

c        New ---------- c

      s = mp*emu

     
     
      
      ep = v/2.
      ppvec2 = ep**2-mtau**2
      if (ppvec2.gt.0.d0)then
         ppvec = dsqrt(ep**2-mtau**2)
      else
         ppvec = 0.d0
      endif


      sthe = dsqrt(1.d0-cthe**2)
      ppx = ppvec*sthe*dcos(phie)
      ppy = ppvec*sthe*dsin(phie)
      ppz = ppvec*cthe

       

      pmx = -ppx
      pmy = -ppy
      pmz = -ppz
      em = ep

c    Lab Frame Momenta --------------------------------------

c      iemu = 3

      absklab = dsqrt(emu**2-mmu2)
c      klab = (/0.d0,0.d0,k3lab,emu/)

      klab0 = emu
      klab1 = 0.d0
      klab2 = 0.d0
      klab3 = absklab


      E1 = (2.d0*s-sx)/2.d0/mp

      

      absk1lab = dsqrt(E1**2-mmu2)

      costh1 = (2.d0*emu*E1-y-2.d0*mmu2)/(2.d0*absklab*absk1lab)
      sinth1 = dsqrt(1.d0-costh1**2)

      

c      k1lab = (/k13lab*dsin(th1),0.d0,k13lab*dcos(th1),E1/)

      k1lab0 = E1
      k1lab1 = absk1lab*sinth1
      k1lab2 = 0.d0
      k1lab3 = absk1lab*costh1

      


      plab0 = mp
      plab1 = 0.d0
      plab2 = 0.d0
      plab3 = 0.d0


c      Z' frame

      rho = dsqrt(sx**2/4.d0+y*mp2)
      costhq = (v2+y+t+t*sx/2.d0/mp2)/(rho*dsqrt(t**2+
     x 4.d0*mp2*t)/mp2)
      

c      if(costhq.gt.1.d0) then
c       costhq = 1.d0
c       print *, 'costhq >1', costhq
c      endif

      sinthq = dsqrt(1.d0-costhq**2)


      

      
      qf0 = t/(2.d0*mp) 
      absqf = dsqrt(t**2+4.d0*mp2*t)/(2.d0*mp)

      if( Isnan(absqf)) then
       write(6,*) 'qvec prob',absqf, t
      endif 


      qfn0 = t/(2.d0*mp)            
      qfn1 = absqf*sinthq*dcos(phiq)
      qfn2 = absqf*sinthq*dsin(phiq)
      qfn3 = absqf*costhq

     

c     back to Z frame

c      abskmk1 = dsqrt((absk1lab*sinth1)**2+(absklab-absk1lab*costh1)**2)    !|k-k1|
      abskmk1 = rho/mp
      cosalpha =  (absklab-absk1lab*costh1)/(abskmk1)
      sinalpha = dsqrt(1.d0-cosalpha**2)

      

      qflab0 = qf0
      qflab1 = qfn1*cosalpha-qfn3*sinalpha
      qflab2 = qfn2
      qflab3 = qfn1*sinalpha+qfn3*cosalpha

      q0 = -qflab0
      q1 = -qflab1
      q2 = -qflab2
      q3 = -qflab3


      kappa0 = klab0-k1lab0+q0   !kappa = k-k1+qF Check!!
      kappa1 = klab1-k1lab1+q1 
      kappa2 = klab2-k1lab2+q2
      kappa3 = klab3-k1lab3+q3 

c      write(6,*) kappa0**2-kappa1**2-kappa2**2-kappa3**2,v2

      kappavec(1) = kappa1
      kappavec(2) = kappa2
      kappavec(3) = kappa3
      kappavec(4) = kappa0

      masskappa = dsqrt(v2)

      pplab(1) = ppx
      pplab(2) = ppy
      pplab(3) = ppz
      pplab(4) = ep

      pmlab(1) = pmx
      pmlab(2) = pmy
      pmlab(3) = pmz
      pmlab(4) = em

c      write(6,*) kappavec(1),kappavec(2),kappavec(3),kappavec(4)

      call BOOST(kappavec,masskappa,pplab)
      call BOOST(kappavec,masskappa,pmlab)

      ppxlab = pplab(1)
      ppylab = pplab(2)
      ppzlab = pplab(3)
      eplab = pplab(4)


      pmxlab = pmlab(1)
      pmylab = pmlab(2)
      pmzlab = pmlab(3)
      emlab = pmlab(4)

c      write(6,*) ppxlab,ppylab,ppzlab

      

      
      kdq = klab0*q0-klab1*q1-
     x klab2*q2-klab3*q3
      k1dq = k1lab0*q0-k1lab1*q1-
     x k1lab2*q2-k1lab3*q3
      kdpp =klab0*eplab-klab1*ppxlab-
     x klab2*ppylab-klab3*ppzlab 
      kdpm =klab0*emlab-klab1*pmxlab-
     x klab2*pmylab-klab3*pmzlab
      k1dpp =k1lab0*eplab-k1lab1*ppxlab-
     x k1lab2*ppylab-k1lab3*ppzlab  
      k1dpm =k1lab0*emlab-k1lab1*pmxlab-
     x k1lab2*pmylab-k1lab3*pmzlab 

      
      

      ppdp = eplab*plab0;

      if(Isnan(costh1)) then
      write(6,*) 'q1 dot product issue'
      endif


c      write(6,*) kdq,k1dq,kdpp,k1dpp,ppdp



      v2lab0 = epslab+emlab
      v2lab1 = ppxlab+pmxlab
      v2lab2 = ppylab+pmylab
      v2lab3 = ppzlab+pmzlab

      absv2lab = (v2lab1**2+v2lab2**2+v2lab3**2)
      v2lab =  v2lab0**2-absv2lab  

c      write(6,*) ppdp,kdq
c      write(6,*) kdq,k1dq,ppdp,kdpp,k1dpp



      return
      end

      subroutine normal_dist(x)
            implicit real*8 (a-z)
            real,parameter :: pi=3.14159265
            call random_stduniform(u1)
            call random_stduniform(u2)
            x = sqrt(-2*log(u1))*cos(2*pi*u2)
      end subroutine normal_dist

      subroutine random_stduniform(u)
            implicit real*8 (a-z)
            call random_number(r)
            u = 1 - r
      end subroutine random_stduniform

C--------------------------------------------------------------
      SUBROUTINE BOOST(P,MP,Q)
C     ---------------------------------------------------
C THIS SUBROUTINE CALCULATES BOOSTED 4-MOMENTA
C Q TO FRAME OF PARTICLE P
C EQ=Q(4)
C PQ=Q(I) ; I=1,3

C MP = MASS OF PARTICLE P -H.BAER-
C--------------------------------------------------

      REAL*8 P(4),Q(4),G,VQ,HQ,MP
      G=P(4)/MP
      VQ=(Q(1)*P(1)+Q(2)*P(2)+Q(3)*P(3))/P(4)
      HQ=VQ*G*G/(1.+G)/P(4)+Q(4)/MP
      Q(4)=G*(Q(4)+VQ)
      DO I=1,3
      Q(I)=Q(I)+HQ*P(I)
      END DO
      RETURN
      END
C---------------------------------------------------------------



