      INCLUDE 'vegas-wgt-ncall.f'
      INCLUDE 'vec-lab-a.f' !vecs and limits
      INCLUDE 'amp-sq-a.f'

      
C     f2py wrapper for the a-diagram contribution.

      subroutine compute_a_sm(z_in, a_in, mp_in, mmu_in, me_in,
     x   emu_in, nbins, e1l_arr, e1u_arr, ncall_in, itmx_in,
     x   nprn_in, results, avgi_arr, sd_arr)

Cf2py intent(in)  :: z_in, a_in, mp_in, mmu_in, me_in, emu_in
Cf2py intent(in)  :: nbins, e1l_arr, e1u_arr, ncall_in, itmx_in
Cf2py intent(in)  :: nprn_in
Cf2py intent(out) :: results, avgi_arr, sd_arr

      implicit real*8 (a-z)
      integer ndim,nprn,igraph,ib,nbins,ncall_in,itmx_in,nprn_in
      real*8 z_in,a_in,mp_in,mmu_in,me_in,emu_in
      real*8 e1l_arr(nbins),e1u_arr(nbins)
      real*8 results(nbins),avgi_arr(nbins),sd_arr(nbins)
      dimension x(8)
      common/masses/me,mmu,mp,mhadmin,dm2,mmu2,mp2
      common/inputs/emu
      common/intvar/capy,capx,mx2,t,v2,dp,dm,ppdp,pmdp
      common/checkvars/denom,xjacp,dpair,lams,dsig
      common/result/s1,s2,s3,s4
      common/tmx/tmin,tmax
      common/vegascalls/calls
      common/chrg/z,a
      common/bsm/F2
      common/E1lim/e1l,e1u
      external sxcint0
      data cm2_scale /1.d-27/

      F2  = 0.d0
      z   = z_in
      a   = a_in
      mp  = mp_in
      mmu = mmu_in
      me  = me_in
      emu = emu_in
      
      if (a.ge.1.d0)then
         mp = mp*a
      endif

      mmu2 = mmu**2
      mp2 = mp**2

      mhadmin = mp 
      dm2 = (mhadmin+2.d0*me)**2-mp**2

      acc = 0.d0
      ndim = 8
      nprn = nprn_in
      igraph = 0

      do ib = 1,nbins

         e1l = e1l_arr(ib)
         e1u = e1u_arr(ib)

         call vegas(sxcint0,acc,ndim,ncall_in,itmx_in,nprn,igraph)

         results(ib)  = s1*cm2_scale
         avgi_arr(ib) = s1*cm2_scale
         sd_arr(ib)   = s2*cm2_scale
      end do

      return
      end

      



      function sxcint0(x) 
      implicit none
      integer ix
      real*8 x(8),me,mmu,mp,mhadmin,dm2,mmu2,mp2,emu,capy,capx,mx2,t,
     x   v2,dp,dm,ppdp,pmdp,lams,lamq,lamY,pi,xjac,skp,ymin,ymax,
     x   dly,dlt,tmin,tmax,mx2min,mx2max,v2min,v2max,sxcint0,cthe,phie,
     x   cthq,phiq,denom,vdsigdv0,dpair,capT,lamt,xmin,xmax,xjacp,dsig,
     x   xlams,capyinv,tinv,alpha,sx,sxmin,sxmax,s,tx,sum,xinv,
     x   xmid,wgt,yx,e1,Ex,e1l,e1u
      common/masses/me,mmu,mp,mhadmin,dm2,mmu2,mp2
      common/inputs/emu
      common/vegaswgt/wgt
      common/intvar/capy,capx,mx2,t,v2,dp,dm,ppdp,pmdp
      common/checkvars/denom,xjacp,dpair,xlams,dsig
      common/lambdas/lamY,lams,lamq
      common/tmx/tmin,tmax
      Common/bins/xmid(1000),xinv(1000),sum
      common/E1lim/e1l,e1u
      
      data pi/3.14159d0/,alpha/7.2993d-3/

      
      xjac = 2.d0*pi !for final phi integral
      skp = 2.*mp*emu + mp2 + mmu2 
      s = 2.*mp*emu

      call ssxminmax(sxmin,sxmax)
      sx = sxmin+(sxmax-sxmin)*x(1)
      capx = s-sx

      e1 = capx/2.d0/mp



      xjac = xjac*(sxmax-sxmin)
      call syminmax(sx,ymin,ymax)
      if (ymin.ge.ymax)then
         write(6,*)'yminmax',ymin,ymax
         sxcint0=0.d0
         return
      endif


      dly = dlog(ymin) + (dlog(ymax)-dlog(ymin))*x(2)
      capy = dexp(dly)
      xjacp = xjac * (dlog(ymax)-dlog(ymin))/capy 


      mx2 = mp2


      call tminmax(skp,capy,capx,mp2,tmin,tmax)
      if (tmin.ge.tmax)then
         sxcint0 = 0.
         return
         write(6,*)'tminmax=',tmin,tmax
      endif
      if (tmin.le.0.d0)write(6,*)'tmin problem',tmin

c     takes care of one power of t in denominator
      dlt = dlog(tmin) + (dlog(tmax)-dlog(tmin))*x(4)
      t = dexp(dlt)
      xjacp = xjacp*(dlog(tmax)-dlog(tmin))/t


      call v2minmax(skp,capy,capx,mp2,t,v2min,v2max)
      v2 = v2min + (v2max-v2min)*x(5)
      if (v2min.ge.v2max)then
         capt = mp2-t-mx2
         capt = -t
         lamt = capt**2+4.*mp2*t 
         sxcint0 = 0.d0
         return
      endif
      xjacp = xjacp*(v2max-v2min)

      phiq = 2.*pi*x(6)
      xjacp = xjacp*2.*pi

      phie = 2.*pi*x(7)
      cthe = -1.+2.*x(8)
      xjacp = xjacp*4.*pi
      dpair = 1./8.*dsqrt(1.d0-4.d0*me**2/v2)

      call v2rest(capy,capx,v2,t,mp2,cthe,phie,ppdp,pmdp,Dp,Dm)

      lamY = (skp - mp2 - mmu2 - capx)**2 + 4.*mp2*capy
      lams = (skp - mp2 - mmu2)**2 - 4.d0*mmu2*mp2
      xlams = lams
      lamq =(v2**2+t**2+capy**2-2.d0*(t*capY-capY*v2-t*v2)) 
      denom = 16.d0*dsqrt(lamY*lams)


      sxcint0 = xjacp*dpair/denom/2./dsqrt(lams) 
      dsig = vdsigdv0(capy,capx,v2,t,mp2,dp,dm,ppdp,pmdp)*s/sx
      sxcint0 = sxcint0*dsig * 8.* alpha**4/pi**4/6.022e+23  
      
       

 10   continue

c      If(e1.gt.e1out.or.e1.lt.e1in) sxcint0 = 0.d0

     

c     -------------------------------------------------------------------------------
      Ex = e1/emu
      ix = Ex*100+1   
      if (ix.gt.100)ix = 100 

      if (e1.le.e1l .or. e1.ge.e1u) then
         sxcint0 = 0.d0
      endif
C     -------------------------------------------------------------------------------


      return
      end

C this function below if for the elastic portion of the cross section
      function vdsigdv0(y,x,v2,t,mx2,dp,dm,ppdp,pmdp)
      implicit none
      real*8 vdsigdv0,y,x,v2,t,capy,mx2,phiq,phie,cthe,me,mmu,
     x   mp,mhadmin,dm2,mmu2,mp2,emu,ep,ppvec,q0,qz,p0,px,pz,
     x   lams,lamq,lamY,pi,alpha,xbj,f2,f1,blockf2,
     x   pro1,pro2,fa,fb,af,dsig,capT,s,sx,ans11,ans12,
     x   ans21,ans22,conv,dp,dm,ppdp,pmdp,bb11,bb12,
     x   bb21,bb22,bb11bsm,bb12bsm,bb21bsm,bb22bsm,tau,ge,gm,z,
     x   a,agev,aegev,fe
      common/masses/me,mmu,mp,mhadmin,dm2,mmu2,mp2
      common/inputs/emu
      common/chrg/z,a
c these variables are all in the eplus-eminus rest frame

      common/lambdas/lamY,lams,lamq !these are defined in xint

      capT = -t
      S = 2.*mp*emu
      sx = 2.*mp*emu - x
      p0 = (capT+sx)/2./dsqrt(v2)

      pz = (v2*(capT-sx)+(t-y)*(capT+sx))/2./dsqrt(v2)/dsqrt(lamq)
      px = dsqrt(p0**2-mp2-pz**2)
      qz = dsqrt(lamq)/2./dsqrt(v2)
      q0 = (v2-t+Y)/2./dsqrt(v2)


      vdsigdv0 = 0.
      xbj = 1.d0
      if (t.le.0.d0)then
         write(6,*)'t problems',t
         return
      endif

      tau = t/4./mp2
      if (a.lt.1.d0)then
         ge = 1./(t/0.71+1.)**2
         gm = ge*2.79           
         fe = 0.d0
      endif
      if (a.eq.1.d0)then
      ge = 1./(t/0.71+1.)**2
      gm = ge*2.79 !check mu
      aegev = 137./0.511d-3
      fe = 1./(1.+aegev**2*t/4.)**2
      else if (a.gt.1.d0)then
         agev =( 0.58+0.82*a**(1.d0/3.d0))*5.07
         ge = z/(agev**2*t/12.+1)**2
         gm =  0.d0
         aegev = 111.7/z**(1.d0/3.d0)/0.511d-3
         fe = z/(1.+aegev**2*t)
      endif
      f1 = gm**2/2.*t !t factor to convert delta(x-1) to
      f2 = (ge**2 + tau*gm**2)/(1.+tau)*t   ! delta(Mx2-Mp2)
      f1 = gm**2/2.*t
      f2 = ((ge-fe)**2+tau*gm**2)/(1.+tau)*t
      

      if (t.le.0.d0)write(6,*)'t.le.0!!'
      if (lamy.le.0.d0)write(6,*)'lamy.le.0!!'
      if (y.le.0.d0)write(6,*)'y.le.0!!'

      call projections(v2,t,y,sx,mp2,dp,dm,ppdp,pmdp,bb11,bb12,bb21,
     x  bb22,bb11bsm,bb12bsm,bb21bsm,bb22bsm)

      pro1 = f1*(bb12+bb12bsm) + 2.*xbj/t*f2*(bb22+bb22bsm)
      pro2 = f1*(bb11+bb11bsm) + 2.*xbj/t*f2*(bb21+bb21bsm)

      fa = 0.5*(pro1-4.*y/lamy*pro2)   
      fb = 2.*y/lamy*(-pro1+12.*y/lamy*pro2)  
      af = (4.*mmu**2-2.*y)*fa +(s*x-mp2*y)*fb
      dsig = 2.    

      conv = 0.3894*6.022d+23
      vdsigdv0 = dsig * sx/s *af * conv

      return
      end

      function dot4(a,b)
      implicit none
      real*8 a(0:3),b(0:3),dot4
      dot4 = a(0)*b(0)-a(1)*b(1)-a(2)*b(2)-a(3)*b(3)
      return
      end

