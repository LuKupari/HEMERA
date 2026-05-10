c this is in the v2 rest frame
      subroutine v2rest(y,x,v2,t,mx2,cthe,phie,ppdp,pmdp,Dp,Dm)
      implicit real*8 (a-z)
      common/inputs/en
      common/masses/me,mmu,mp,mhadmin,dm2,mmu2,mp2
      v = dsqrt(v2)
      capT = mp2-t-mx2
      s = 2.*mp*en
      sx = s-x
      q0 = (v2-t+y)/2./v
      lamq = v2**2+t**2+y**2-2.*(t*y-t*v2-y*v2)
      if (lamq.le.0.d0)lamq = 0.d0
      qz = dsqrt(lamq)/2./v
      p0 = (capT+sx)/2./v
      pz = (v2*(capT-sx)+(t-y)*(capT+sx))/2./v/dsqrt(lamq)
      pz = -pz  
      px2 = p0**2-pz**2-mp2
      if (px2.le.0.d0)then

         px = 0.
      else
         px = dsqrt(px2)
      endif
      ep = v/2.
      ppvec2 = ep**2-me**2
      if (ppvec2.gt.0.d0)then
         ppvec = dsqrt(ep**2-me**2)
      else
         ppvec = 0.d0
      endif
      the = dacos(cthe)
      sthe = dsin(the)
      ppx = ppvec*sthe*dcos(phie)
      ppy = ppvec*sthe*dsin(phie)
      ppz = ppvec*cthe
      pmx = - ppx
      pmy = - ppy
      pmz = - ppz
      em = ep
      ppdp = ep*p0-px*ppx-pz*ppz
      pmdp = em*p0-px*pmx-pz*pmz
      qdpp = q0*ep - qz*ppz
      qdpm = q0*em - qz*pmz
      Dp = -t - 2.*qdpp
      Dm = -t - 2.*qdpm
      return
      end


      subroutine ssxminmax(sxmin,sxmax)
      implicit none
      real*8 sxmin,sxmax,me,mmu,mp,mhadmin,dm2,mmu2,mp2,en,s,
     x  xlams,arg,skp,sxmax0,t1,t2
      common/masses/me,mmu,mp,mhadmin,dm2,mmu2,mp2
      common/inputs/en

      s = 2.*mp*en
      skp = s + mmu2 + mp2
      xlams = s**2 - 4.*mmu2*mp2
      arg = xlams*(xlams - 2.d0*s*dm2 - 4.d0*mmu2*dm2+dm2**2)
c      write(6,*)'arg =',arg
      sxmin = 0.d0 
      sxmax = s - 2.d0*mmu*mp
      if (arg.ge.0.d0)then
         t1 = xlams + dm2*(s+2.d0*mp2)
         t2 = dsqrt(arg)
         sxmin = (xlams +dm2*(s+2.d0*mp2) - dsqrt(arg))/2.d0/skp
      endif

      if (sxmin.ge.sxmax)then
         write(6,*)'sx problem',sxmin,sxmax,sxmax0
      endif
      return
      end

      subroutine syminmax(sx,ymin,ymax)
      implicit none
      real*8 ymin,ymax,me,mmu,mp,mhadmin,dm2,mmu2,mp2,en,s,
     x  xlams,arg,skp,t1,t2,sx,ymax0,t1fac,targfac,yminck,arg2
      common/masses/me,mmu,mp,mhadmin,dm2,mmu2,mp2
      common/inputs/en
      s = 2.*mp*en
      skp = s + mmu2 + mp2
      xlams = s**2 - 4.*mmu2*mp2

      arg = (xlams - s*sx)**2 - 4.*mp2*mmu2*sx**2
      arg2 = xlams *(xlams+sx**2-2.*s*sx)
      t1 = (xlams - s*sx)/2./mp2
      t2 = dsqrt(arg) /2./mp2 
      ymin = (xlams - s*sx)/2./mp2 - dsqrt(arg)/2./mp2 
      yminck = -4.*mmu2*mp2 + 2.*mmu2*mp2*(s/(s-sx)+(s-sx)/s)
      yminck = yminck/2.d0/mp2

      yminck = sx**2*mmu2/s**2  
      yminck = (2.*mmu2*mp2*sx**2/(xlams - s*sx)*
     x      (1.+mmu2*mp2*sx**2/(xlams-s*sx)**2))/2./mp2
      if (ymin.lt.1.d-4)ymin = yminck
      ymax = sx - dm2
      ymax0 = t1 + t2 
c      write(6,*)'ymax,ymax0',ymax,ymax0
      if (ymax.gt.ymax0)ymax = ymax0
      if (ymin.ge.ymax)then
         write(86,*)'yminmax',ymin,ymax,arg,sx
         ymin = ymax
      endif
      return
      end



      
      subroutine tminmax(skp,y,x,mhad2,tmin,tmax)!here mx2 = mp2
      implicit none
      real*8 me,mmu,mp,mhadmin,dm2,mmu2,mp2
      real*8 skp,y,x,mhad2,tmin,tmax,me2,sp,s,w2,sx,mx2,t1,
     x      lamy,lamp,t2,diff
      common/masses/me,mmu,mp,mhadmin,dm2,mmu2,mp2

      me2 = me**2
      sp = skp
      s = sp - mp**2 - mmu**2
      w2 = mp2 + s - x - y
      sx = s - x
      mx2 = mhad2
      t1 = sx*(w2 - mx2)+2.*mx2*y - 4.*me2*(sx+2.*mp2)
      lamy = sx**2 + 4.*mp2*y
      lamp = w2**2+(4.*me2)**2+mx2**2 
     x     - 2.*(w2*(4.*me2+mx2)+4.*me2*mx2)
      if (lamy*lamp.lt.0.d0)then
         t2 = 0.d0
      else
         t2 = dsqrt(lamy*lamp)
      endif
      tmin = (t1 - t2)/2./w2
      tmax = (t1 + t2)/2./w2

c      if (tmin.gt.tmax)then
c         write(6,*)'tmin.gt.tmax'
c         tmin = tmax
c      endif

      if (tmax.le.0.d0)then
         write(6,*)'tmax.le.0',tmax
         tmax = tmin
      endif
      if (tmin.le.0.d0)then

c         write(6,*)'tmin neg=',tmin
c         write(6,*)'t1,t2,tmax vals=',t1,t2,tmax
c         write(6,*)'w2,mx2,y,sx',w2,mx2,y,sx
         tmin = 2.*y*(mx2-mp2*y/sx*(w2-mx2))/2./w2
         tmin = 2.*(y+4.*me2)*(mx2-mp2)/2./w2 
         tmin = 2.*mp2*(y-4.*me2)+sx/w2*8.*me2*mp2 
     x        + 2.*y*mp2*(mp2+4.*me2)/sx -2.*y*mp2*w2/sx
     x      + sx/w2**2*(2.*me2*(mp2**2+4.*me2*mp2-16.*me2**2)
     x      - mp2**3/2.)

         tmin = 2.d0*mp2*(y + 4.d0*me2)**2/sx 
         tmin = tmin + 2.d0*mp2*(mp2+4.d0*me2)*(y+4.d0*me2)**2/sx**2 
         tmin = tmin/2.d0/w2
c         write(6,*)sx,y,me2 
      endif
      if (tmin.ge.tmax)then
         diff = 2.*t2/2./w2
         write(6,*)'diff = ',diff,lamy,lamp,w2,sx,y
         tmax = tmin+diff
      endif
      if (t2.eq.0.d0)then
      tmax = tmin
      endif
      return
      end

      subroutine V2minmax(skp,y,x,mhad2,t,v2min,v2max)
      implicit none
      real*8 me,mmu,mp,mhadmin,dm2,mmu2,mp2
      real*8 skp,y,x,mhad2,t,v2min,v2max,me2,sp,s,w2,sx,mx2,
     x   lamt,lamy,bigt,dlamt,dlamy,expand1,expand2
      common/masses/me,mmu,mp,mhadmin,dm2,mmu2,mp2

      me2 = me**2
      v2min = 4.*me2
      sp = skp
      s = sp - mp**2 - mmu**2
      w2 = mp2 + s - x - y
      sx = s - x
      mx2 = mhad2
      bigt = mp2 - t - mx2
      lamt = bigt**2 + 4.*mp2*t
      dlamt = dsqrt(lamt)
      lamy = sx**2 + 4.*mp2*y
      dlamy = dsqrt(lamy)
      expand1 = lamy/sx**2-1.d0
      expand2 = lamt/bigt**2-1.d0
      if (expand2.lt.1.d-5)then  
         v2max = (-sx/bigt - 1.d0)*t +(-bigt/sx-1.d0)*y
      else
         v2max = (sx*bigt + dlamt*dlamy)/2./mp2 - t - y
      endif
      if (v2max.le.v2min)then
       v2max=v2min  
      endif
      return
      end


