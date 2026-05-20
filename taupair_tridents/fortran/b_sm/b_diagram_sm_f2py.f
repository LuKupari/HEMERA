      INCLUDE 'vegas-wgt-ncall.f'
      INCLUDE 'Amp-sq-b.f'
      INCLUDE 'vec-lab-b-sm.f'
      
C     f2py wrapper for the SM b-diagram contribution.

      subroutine compute_b_sm(z_in, a_in, mp_in, mmu_in, me_in,
     x   emu_in, nbins, e1l_arr, e1u_arr, ncall_in, itmx_in,
     x   nprn_in, distribution_in, results, avgi_arr, sd_arr)

Cf2py intent(in)  :: z_in, a_in, mp_in, mmu_in, me_in, emu_in
Cf2py intent(in)  :: nbins, e1l_arr, e1u_arr, ncall_in, itmx_in
Cf2py intent(in)  :: nprn_in, distribution_in
Cf2py intent(out) :: results, avgi_arr, sd_arr

      implicit real*8 (a-z)
      integer ndim,nprn,igraph,ib,nbins,ncall_in,itmx_in,nprn_in,i
      integer distribution_in,nbins_common
      logical :: emuout
      real*8 z_in,a_in,mp_in,mmu_in,me_in,emu_in
      real*8 e1l_arr(nbins),e1u_arr(nbins)
      real*8 results(nbins),avgi_arr(nbins),sd_arr(nbins)
      dimension x(8)
      common/masses/me,mmu,mp,dm2,mmu2,mp2
      common/inputs/emu
      common/intvar/capsx,capy,t,v2,kdq,k1dq,k1dpp,kdpp,ppdp
      common/checkvars/denom,lams,dsig
      common/result/s1,s2,s3,s4
      common/chrg/z,a
      common/vegascalls/calls
      external sxcint0
      common/bins/xmid(1000),xinv(1000),xinv2(1000),sum,
     x   nbins_common
      common/E1lim/e1l,e1u
      common/InvMasslim/mttl,mttu
      common/arg/emuout
      data cm2_scale /1.d-27/


      emu = emu_in
      mp  = mp_in
      mmu = mmu_in
      me  = me_in
      z   = z_in
      a   = a_in
      
      if (a.ge.1.d0)then
         mp = mp*a
      endif

      dm2 = 4.d0*me**2+4.d0*me*mp
      mmu2 = mmu**2
      mp2 = mp**2
      me2 = me**2

      acc = 0.d0
      ndim = 8
      nprn = nprn_in
      igraph = 0

      do i = 1,100
         xmid(i) = dfloat(i)
      end do

      emuout = distribution_in .eq. 0
      nbins_common = nbins

      if(emuout) then

         do ib = 1,nbins

            e1l = e1l_arr(ib)
            e1u = e1u_arr(ib)

            do i = 1,100
               xinv(i) = 0.d0
            end do

            call vegas(sxcint0,acc,ndim,ncall_in,itmx_in,nprn,igraph)

            xnorm = itmx_in*calls
            sum = 0.d0
            do i = 1,100
               sum = sum + xinv(i)/xnorm
            end do
            if (Isnan(sum)) sum = 0.d0

            if (Isnan(s1)) s1 = 0.d0
            if (Isnan(s2)) s2 = 0.d0
            results(ib)  = sum*cm2_scale
            avgi_arr(ib) = s1*cm2_scale
            sd_arr(ib)   = s2*cm2_scale
         end do


      else

!        New --------------------------------------------
         do ib = 1,nbins

            mttl = e1l_arr(ib)
            mttu = e1u_arr(ib)
            do i = 1,100
               xinv2(i) = 0.d0
            end do

            call vegas(sxcint0,acc,ndim,ncall_in,itmx_in,nprn,igraph)

            xnorm = itmx_in*calls
            sum = 0.d0
            do i = 1,100
                  sum = sum + xinv2(i)/xnorm
            end do
            if (Isnan(sum)) sum = 0.d0

            if (Isnan(s1)) s1 = 0.d0
            if (Isnan(s2)) s2 = 0.d0
            results(ib)  = sum*cm2_scale
!         avgi_arr(ib) = s1*cm2_scale
!         sd_arr(ib)   = s2*cm2_scale
         end do

      endif
!------------------------------------------------------

      return
      end

      
      function sxcint0(x) !this has a delta function for xbj
      implicit none
      integer ix,iy
      integer nbins_common
      real*8 x(8),me,mmu,mp,mhadmin,dm2,mmu2,mp2,emu,capy,
     x   capsx,mx2,t,
     x   v2,dp,dm,k1dpp,kdpp,lams,lamq,lamY,pi,xjac,s0,
     x   ymin,ymax,ymax1,
     x   dly,dlt,tmin,tmax,mx2min,mx2max,v2min,v2max,
     x   sxcint0,cthe,phie,
     x   cthq,phiq,denom,capT,dsig,sx,sxmin,sxmax,s,
     x   rho,kdq,k1dq,m2,me2,sum,dlsx,
     x   v2smr,ppdp,dsig0,kvec,domegakap,dphi1,
     x   argsx,bigX,arg,t1,dt,xmid,xinv,wgt,
     x   t2,u,xjacp,alpha,xjacp1,v2lab,vx,e1,lamqarg,
     x   Ex,e1l,e1u,mttl,mttu,xinv2
      logical emuout
      common/masses/me,mmu,mp,dm2,mmu2,mp2
      common/inputs/emu
      common/vegaswgt/wgt
      common/intvar/sx,capy,t,v2,kdq,k1dq,k1dpp,kdpp,ppdp
      common/checkvars/denom,lams,dsig
      common/bins/xmid(1000),xinv(1000),xinv2(1000),sum,
     x   nbins_common
      common/E1lim/e1l,e1u
      common/InvMasslim/mttl,mttu
      common/arg/emuout

      data pi/3.14159d0/,alpha/7.2993d-3/

      
      me2 = me**2
      m2 = 4.d0*me**2
      s0 = 2.*mp*emu + mp2 + mmu2 !sprime
      s = mp*emu
      s0 = 2.*s + mmu2 + mp2
      lams = 4.*s**2 - 4.*mmu2*mp2
     
      argsx = lams*(lams - 4.d0*s*dm2 - 4.d0*mmu2*dm2+dm2**2)
      if (argsx.lt.0.d0) then
         sxcint0 = 0.d0
         return
      endif

      sxmin = (lams +2.d0*dm2*(s+mp2) - dsqrt(argsx))/2.d0/s0
      sxmax = 2.*(s - 2.d0*mmu*mp)

      if (sxmin.lt.4.*me*(me+mp)) sxmin = 4.*me*(me+mp)

      if (sxmin.ge.sxmax)then
         write(6,*)'sx problem',sxmin,sxmax
      endif

      capsx = sxmin+(sxmax-sxmin)*x(1)
      xjacp =(sxmax-sxmin)

      if(sxmax.lt.sxmin) then
         sxcint0 = 0.d0
         return
      endif

      if(Isnan(capsx)) then
         sxcint0 = 0.d0
         return
      endif

      bigX = 2.d0*s-capsx
      t1 = (lams-2.d0*s*capsx)/(2.d0*mp2)
      arg = (lams-2.d0*s*capsx)**2-4.d0*mmu2*mp2*capsx**2

      if (arg.lt.0.d0)then
         sxcint0 = 0.d0
         return
      endif

      ymin = t1-(1/2.d0/mp2)*dsqrt(arg)
      ymax = capsx-4.d0*me*(me+mp)

      ymax1 = s*(2.d0*s-capsx)/mp2-2.d0*mmu2+
     x dsqrt(lams*(bigX**2/4.d0-mmu2*mp2))/mp2 !check 

      if(ymax.gt.ymax1) ymax = ymax1

C this takes care of one power of capy in the denominator
      capy = ymin+(ymax-ymin)*x(2)
      xjacp = xjacp*(ymax-ymin)

      dphi1 = 2.d0*pi*x(3)
      xjacp = xjacp*2.d0*pi

      e1 = bigX/2./mp
      
      t1 = capsx*(capsx-capy)+2.d0*mp2*capy-m2*(capsx+2.d0*mp2)
      t2 = 2.d0*(capsx-capy+mp2)
      u = (capsx**2+4.d0*capy*mp2)*((capsx-capy)**2-
     x 2.d0*m2*(capsx-capy)+4.d0*m2*(me2-mp2))

      if (u.lt.0.d0)then
         sxcint0 = 0.d0
         return
      endif

      tmin = (t1 - dsqrt(u))/t2
c      tmax = (t1 + dsqrt(u))/t2
      tmax = 0.1

      if (tmin.ge.tmax)then
         sxcint0 = 0.
         return
         write(6,*)'tminmax=',tmin,tmax
      endif

      if (tmin.le.0.d0) then
      tmin = 1.e-10
      endif

      dlt = dlog(tmin) + (dlog(tmax)-dlog(tmin))*x(4)
      t = dexp(dlt)
      dt = t*(dlog(tmax)-dlog(tmin))

      xjacp = xjacp*dt

c      write(6,*)  tmax,tmin

      phiq = 2.*pi*x(5)
      xjacp= xjacp*2.*pi

      v2min = 4.*me2

      rho = dsqrt(capsx**2/4.+capy*mp2)
      v2max = rho*dsqrt(t**2+4.d0*mp2*t)/mp2-(capy+t+t*capsx/(2.*mp2))

      if(v2max.lt.v2min) v2max = v2min

      v2 = v2min + (v2max-v2min)*x(6)
      xjacp = xjacp*(v2max-v2min)

      if(v2.le.0.d0) then
         sxcint0 = 0.d0
         return
      endif

      phie = 2.*pi*x(7)
      xjacp = xjacp*2.*pi

      cthe = -1.+2.*x(8)
      xjacp = xjacp*2.d0

      domegakap = 2.d0*pi*2.d0

      xjacp = xjacp
      
      lamqarg = v2**2+me2**2+me2**2
     x   -2.d0*(me2*v2+me2*v2+me2*me2)
      if(lamqarg.le.0.d0) then
         sxcint0 = 0.d0
         return
      endif
      lamq = dsqrt(lamqarg)

      kvec = dsqrt(emu**2-me**2)
      rho = dsqrt(capsx**2/4.d0+capy*mp2)
      denom = 8.d0*v2*8.d0*rho*8.d0*mp*kvec

      if(denom.le.0.d0) then
         sxcint0 = 0.d0
         return
      endif

      sxcint0 = xjacp*lamq/denom 

      if(Isnan(lamq)) then
         sxcint0 = 0.d0
         return
      endif

       
      call v2rest(capy,capsx,v2,t,cthe,phie,phiq,k1dpp,kdpp,kdq,    
     x k1dq,ppdp)



      dsig = dsig0(m2,v2,t,capsx,capy,k1dpp,kdpp,kdq,
     x k1dq,ppdp)

      if(Isnan(dsig)) then
         sxcint0 = 0.d0
         return
      endif

      sxcint0 = sxcint0*dsig
      if(Isnan(sxcint0)) sxcint0 = 0.d0

      Ex = e1/emu
      ix = Ex*100+1
      if (ix.gt.100) ix = 100
      if (ix.lt.1) ix = 1
      if (e1.gt.e1l .and. e1.lt.e1u) then
         xinv(ix) = xinv(ix) + sxcint0*wgt
      endif

!-------------------------New-------------
      vx = dsqrt(v2)
      iy = vx-0.5d0
      if (iy.gt.nbins_common) iy = nbins_common
      if (iy.lt.1) iy = 1
      If(vx.gt.mttl.and.vx.lt.mttu) then
         xinv2(iy) = xinv2(iy) + sxcint0*wgt
      endif
!-------------------------------------------      


      return
      end

      function dsig0(m2,v2,t,Sx,Y,k1dpp,kdpp,kdq,k1dq,ppdp)
      implicit none
      real*8 m2,mmu2,Mp2,v2,t,Sx,Y,k1dpp,kdpp,kdq,k1dq,B11,B12,emu,
     x     B21,B22,z,a,tau,ge,gm,fe,aegev,f1,f2,lams,Ein,pi,alpha,
     x     agev,dsig0,ppdp,dm2,me,mmu,mp,X,s,conv
      common/inputs/emu
      common/chrg/z,a
      common/masses/me,mmu,mp,dm2,mmu2,mp2
      data pi/3.14159d0/,alpha/7.2993d-3/      
      
      
      s = mp*emu
      mmu2 = (0.10566)**2  
      conv = 0.38938    
C need elastic form factors instead, with delta function jacobian
      tau = t/4./mp2
      lams = (2.*emu*dsqrt(mp2))**2 -4.*mmu2*Mp2

      if (a.lt.1.d0)then            !!from hallsie
         ge = 1./(t/0.71+1.)**2
         gm = ge*2.79           !check mu
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

      f1 = 4.d0*mp2*tau*gm**2
      f2 = 4.d0*mp2*(((ge-fe)**2+tau*gm**2)/(1.+tau))

      X = 2.d0*s-Sx

      call projections(m2/4.d0,Mp2,emu,v2,t,Y,k1dpp,kdpp,                        
     x       kdq,k1dq,
     x       ppdp,X,B11,B12,B21,B22)


      dsig0 = ((B11 + B21)*f1 !factor 4 from leptonic part
     x     + (B12 + B22)*f2) 

      dsig0 =dsig0*(4.d0*pi*alpha)**4/2./dsqrt(lams)/(2.d0*pi)**8                    
      dsig0 = dsig0/v2**2/t**2*conv

      return
      end

      function dot4(a,b)
      implicit none
      real*8 a(0:3),b(0:3),dot4
      dot4 = a(0)*b(0)-a(1)*b(1)-a(2)*b(2)-a(3)*b(3)
      return
      end



 
