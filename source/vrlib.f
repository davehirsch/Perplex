c-----------------------------------------------------------------------

c RLIB - a library of subprograms called by VERTEX, FRENDLY, ACTCOR
c BUILD, ADIABAT, and ISOCHOR.

c Copyright (c) 1990 by James A. D. Connolly, Institute for Mineralogy
c & Petrography, Swiss Federal Insitute of Technology, CH-8092 Zurich,
c SWITZERLAND. All rights reserved.
 
c Please do not distribute this source.

c-----------------------------------------------------------------------

      character*10 function gname (id)
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer id
 
      character names*8, fname*10
 
      common/ cst8  /names(k1)/ csta7 /fname(h9)

      if (id.lt.0) then
         gname = names(-id)
      else if (id.gt.0) then 
         gname = fname(id)
      end if 

      end 

      subroutine gexces (id,dg)
c-----------------------------------------------------------------------
c gexces evaluates the contribution to the gibbs energy of a pseudocompound
c arising from configurational entropy, excess properties, and dqf corrections
c for special cases (internal EoS's or Van Laar gexces does not return
c the excess energy, see routines fexces, gvlaar, and toop).
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer id

      double precision p,t,xco2,u1,u2,tr,pr,r,ps,dg
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      integer ixp,ifp
      double precision sxs,exces
      common/ cst304 /sxs(k13),exces(m3,k1),ixp(k1),ifp(k1)
c-----------------------------------------------------------------------
      
      dg = exces(1,id) + t * exces(2,id) + p * exces(3,id)

      end 

      double precision function gfluid (y)
c-----------------------------------------------------------------------
c gfluid returns the fugacities computed from a binary phase described 
c by an internal EoS referenced by cfluid, the composition of the phase
c is y (the amount of component 2).
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision fo2, fs2, y

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      double precision f
      common/ cst11 /f(2)
c-----------------------------------------------------------------------
      xco2 = y

      call cfluid (fo2,fs2)

      gfluid = r*t*((1d0-y)*f(1) + y*f(2))

      end 


      subroutine fexces (id,dg)
c-----------------------------------------------------------------------
c gexces evaluates the contribution to the gibbs energy of a pseudocompound
c arising from dqf corrections, and excess properties computed from an
c internal fluid EoS
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer id

      double precision dg,cp,f,fo2,fs2

      common/ cst11 /f(2)/ cst12 /cp(k5,k1)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      integer ixp,ifp
      double precision sxs,exces
      common/ cst304 /sxs(k13),exces(m3,k1),ixp(k1),ifp(k1)
c-----------------------------------------------------------------------

      dg = exces(1,id) + t * exces(2,id) + p * exces(3,id)

      xco2 = cp(2,id)

      call cfluid (fo2,fs2)

      dg = dg + r*t*(cp(1,id)*f(1) + cp(2,id)*f(2))

      end 

      double precision function gzero (id)
c----------------------------------------------------------------------
c gzero computes the reference pressure free energy of a compound 
c with no transitions identified by the argument 'id'. see gcpd. 
c----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer id,j,jfct,jmct,jprct

      double precision thermo,uf,us
 
      common/ cst1 /thermo(k4,k10),uf(2),us(h5)
     *      / cst307 /jfct,jmct,jprct

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      double precision mu
      common/ cst39 /mu(i6)

      double precision vnumu
      common/ cst44 /vnumu(i6,k10)
c----------------------------------------------------------------------
 
      gzero = thermo(1,id)
c                                 -sdt
     *      + t * (thermo(2,id) - thermo(4,id) * dlog(t)
     *      - t * (thermo(5,id) + thermo(7,id) * t))
     *      - (thermo(6,id) + thermo(10,id) / t) / t
     *      + thermo(8,id) * dsqrt(t) + thermo(9,id)*dlog(t)
 
      do j = 1, jmct      
c                                 -ndu   
         gzero = gzero - vnumu(j,id) * mu(j)
      end do

      end 

      subroutine gcpd (id,gval)
c-----------------------------------------------------------------------
c gcpd computes the gibbs free energy of a compound identified by
c the arguement 'id' from the thermochemical parameters stored
c in the array 'thermo' which is located in common block 'cst1'.
c the parameters are: g(pr,tr) [thermo(1,id)], s(pr,tr) [thermo
c (2,id)], v(pr,tr) [thermo(3,id)], the coefficients of the extended
c maier-kelley heat capacity equation:
 
c           cp(pr,t)=a+b*t+c/(t*t)+d*t*t+e/t**(1/2)+f/t+g/t**3
 
c [thermo(4-10,id)], and the coefficients of the volumetric 
c equations given in the program documentation (Eqs 2.1-2.3):

c        b8 = 0 =>
 
c           v(p,t) = v(pr,tr) + b2*(t-tr) + b4*(p-pr)
c                           + b6*(p-pr)**2 + b7*(t-tr)**2

c        -3d0 < b8 < 0 =>

c           v(p,t) = v(pr,tr) * exp[ b3*(t-tr) + b8*(p-pr) ]

c        b8 < -3d0 =>

c           p = 3 * k * f * (2f + 1)^(5/2) * (1 - 3/2 * f * (4 + b8))
c           f = 1/2 * ((v0/v)^(2/3) - 1)
c           v0 = v(pr,tr) + b1 * (t - tr)
c           k = -v0 / (b2 + b3 * (t - tr))

c        b8 > 0 =>

c           v(p,t) = v(pr,t) * [1 + b8*p/(b8*p+Kt)]^(1/b8)
c               Kt = b6 + b7 * t
c            alpha = b1 + b2*t + b3/t - b4/t^2 + b5/t^(1/2)      

c [thermo(11-18,id)].
c---------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer id,j,jfct,jmct,jprct,jd,idis,lmda,ltyp

      double precision thermo,uf,us,
     *       gmurn, gsixtr, gstxgi,therdi,therlm,f,fs2,fo2,
     *       gss, del, dg, rkt, gval, gmake
 
      common/ cst1 /thermo(k4,k10),uf(2),us(h5)/ cst11 /f(2)
     *      / cst203 /therdi(m8,m9),therlm(m7,m6,k9)
     *      / cst204 /ltyp(k1),lmda(k1),idis(k1)
     *      / cst307 /jfct,jmct,jprct

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      double precision mu
      common/ cst39 /mu(i6)

      integer make
      common / cst335 /make(k10)

      integer isfp
      common/ cst303 /isfp(k10)

      double precision vnumu
      common/ cst44 /vnumu(i6,k10)
c--------------------------------------------
      save rkt

      data rkt/0d0/
c---------------------------------------------------------------------

      if (make(id).ne.0) then 
c                                 the phase is a made phase, compute 
c                                 and sum the component g's.
         gval = gmake (id)

         goto 999

      else if (thermo(3,id).lt.0d0) then
c                                 identify sixtrude's EoS by negative
c                                 entropy term and 
c                                 identify the version of stxrude's
c                                 EoS by sign on number of atoms
         if (thermo(2,id).gt.0) then 
c                                 use this call to evalute the 
c                                 debye T as in stixrude JGR '05
            gval = gsixtr (id)
         else 
c                                 and this call to evalute the 
c                                 debye T as in stixrude JGI '05
            gval = gstxgi (id) 
         end if 

         goto 999

      end if 
c                                 all other EoS's with Cp function
      gval = thermo(1,id)
c                                 -sdt
     *       + t * (thermo(2,id) - thermo(4,id) * dlog(t)
     *       - t * (thermo(5,id) + thermo(7,id) * t))
     *       - (thermo(6,id) + thermo(10,id) / t) / t
     *       + thermo(8,id) * dsqrt(t) + thermo(9,id)*dlog(t)
c                                 vdp-ndu term:
      if (thermo(18,id).eq.0d0) then 
c                                 normal polynomial:
         del =  p * (thermo(3,id) 
     *               + (thermo(17,id) * t + thermo(12,id)) * t
     *               + (thermo(16,id) * p + thermo(14,id)) * p)

      else if (thermo(18,id).gt.0d0) then
c                                 HP-murnaghan EoS:
c                                 vdp = V(1,T)*KT**(1/K')/(K'-1)
c                                 * [(KT+K'*p)**(1-1/K') -
c                                 (KT+K'*pr)**(1-1/K')]
         rkt = thermo(16,id) + thermo(17,id)*t
c                                 a ****wit has entered a ridiculous
c                                 temperature
         if (rkt.lt.0d0) call error (46,t,id,'GCPD') 

         del = (thermo(3,id) + (thermo(11,id) + thermo(12,id)*t)*t +
     *          thermo(13,id)*dlog(t) + thermo(14,id)/t + 
     *          thermo(15,id)*dsqrt(t)) * rkt**thermo(21,id)
     *          /thermo(22,id) 
     *          *((rkt+thermo(18,id)*p)**thermo(19,id)
     *           -(rkt+thermo(20,id))**thermo(19,id))

      else if (thermo(18,id).lt.-3d0) then
c                                 Ghiorso-murnaghan

         del = gmurn (id)

      else 
c                                 gottschalk.
         del = thermo(11,id)*dexp(thermo(13,id)*t)*
     *               (1d0 - dexp(thermo(18,id)*(p - pr)))

      end if
c                                 -ndu term
      do j = 1, jmct      
         del = del - vnumu(j,id) * mu(j)
      end do

      gval = gval + del 
c                                 check for lambda transitions:
      if (ltyp(id).ne.0) then
 
         if (ltyp(id).lt.4) then
c                                 ubc-type transitions
            call lamubc (p,t,dg,lmda(id),ltyp(id))
            gval = gval + dg
 
         else if (ltyp(id).lt.7) then
c                                 standard transitions
            call lamhel (p,t,gval,del,lmda(id),ltyp(id))
 
         else if (ltyp(id).lt.10) then
c                                 supcrt q/coe lambda transition not correct 
c                                 as of 12/07 cause del = vdp - ndu, should 
c                                 be just ndu. JADC
            call lamqtz (p,t,gval,del,lmda(id),id)
 
         else if (ltyp(id).gt.9 .and. ltyp(id).le.12 ) then
c                                 holland and powell (landau model)
            call lamlan (dg,rkt,lmda(id),ltyp(id))
            gval = gval + dg
 
         end if
      end if
c                                 check for temperature dependent
c                                 order/disorder:
      if (idis(id).ne.0) then
         jd = idis(id)
         call disord (t,therdi(1,jd),therdi(2,jd),therdi(3,jd),
     *                  therdi(4,jd),therdi(5,jd),therdi(6,jd),
     *                  therdi(7,jd),therdi(8,jd),therdi(9,jd),gss)
         gval = gval + gss
      end if
c                                 fluids in the saturated
c                                 or thermodynamic composition space
      if (isfp(id).ne.0) then 

         if (isfp(id).eq.1) then 
            xco2 = 0d0 
            call cfluid (fo2,fs2)
            gval = gval + r*t*f(1)
         else
            xco2 = 1d0
            call cfluid (fo2,fs2)
            gval = gval + r*t*f(2)  
         end if          

      end if
 
999   end


      subroutine zeroi (iarray,index,ivalue,n)

      implicit none

      integer n, iarray(n), index, ivalue, i

      do i = 1, index
         iarray(i) = ivalue
      end do 

      end

      subroutine xchk (xmin, xmax, xinc, tname) 

      implicit none

      double precision xmin, xmax, xinc

      character tname*10

      if (xmax.gt.1d0) then 
         call warn (21,xmax,1,tname)
         xmax = 1d0
      end if 

      if (xmin.lt.0d0) then 
         call warn (22,xmin,1,tname) 
         xmin = 0d0
      end if 

      if (xmax.lt.xmin) then
         call warn (23,xmax,1,tname)
         xmax = 1d0
         xmin = 0d0
      end if 

      if (xinc.le.0d0) then 
         call warn (23,xinc,1,tname)
         xinc = 1d0       
      end if

      end 

      subroutine zchk (z, tname)

      implicit none

      double precision z

      character tname*10

      if (z.gt.0d0.and.z.le.1d-4) then
         z = 0d0
      else if (z.ge.0.9999d0.and.z.lt.1d0) then
         z = 1d0
      else if (z.ge.0d0.and.z.le.1d0) then 
         return
      else if (z.gt.1d0.and.z.lt.1.0001d0) then
         z = 1d0
      else if (z.gt.-1d-3) then
         z = 0d0
      else 
         call error (125,z,1,tname)
      end if

      end 

      subroutine assptx
c---------------------------------------------------------------------
c assptx assigns the current value of v1 and v2 to the array ptx and
c increments the ordinate counter ipt2.
c---------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer ipt2

      double precision ptx,v,tr,pr,r,ps
 
      common/ cst32 /ptx(l5),ipt2
     *      / cst5 /v(l2),tr,pr,r,ps

      integer ipot,jv,iv1,iv2,iv3,iv4,iv5
      common / cst24 /ipot,jv(l2),iv1,iv2,iv3,iv4,iv5

  
      ipt2 = ipt2 + 2

      if (ipt2.gt.l5) ipt2 = l5

      ptx(ipt2-1) = v(iv1)
      ptx(ipt2)   = v(iv2)
 
      end
 
      subroutine concrt
c---------------------------------------------------------------------
c concrt determines convergence criteria and limits for reasonable
c solutions for the routine univeq. the array delt is also used in
c the routine slope.
 
c input:  dv - array of default search increments.
c         vmax,vmin - arrays containing the maximum and minimum values
c           of the independent and dependent intensive variables.
c         iv - array indexes variables in vmax, vmin, and output arrays
c output: ulim, blim - arrays with upper and lower limits for reasonable
c           solutions (vmax,vmin+/-4*dv)
c         delt - array containing the finite difference increments
c           (vmax-vmin/10**5), delt also serves as the test value
c           for convergence.
c---------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i

      double precision ddv
 
      double precision blim, ulim, dgr
      common/ cxt62 /blim(l2),ulim(l2),dgr

      double precision vmax,vmin,dv
      common/ cst9 /vmax(l2),vmin(l2),dv(l2)
c---------------------------------------------------------------------
      do i = 1, l2
         if (dv(i).lt.0d0) call error (34,dv(i),i,'CONCRT')
         if (i.eq.3) then 
            ulim(i) = vmax(i)
            blim(i) = vmin(i)
         else if (i.eq.1.or.i.eq.2) then
            ulim(i) = vmax(i)+dv(i)
            blim(i) = vmin(i)-dv(i)
            if (blim(i).lt.0d0) blim(i) = 1d0 
         else 
            ulim(i) = vmax(i)+dv(i)
            blim(i) = vmin(i)-dv(i)
         end if 
         ddv = vmax(i)-vmin(i)
         if (ddv.lt.0d0) call error (35,ddv,i,'CONCRT')
      end do 
 
      end
 
      subroutine conver (g,s,v,a,b,c,d,e,f,gg,b1,b2,b3,b4,b5,b6,b7,b8,
     *                   b9,b10,b11,b12,tr,pr,r)
c---------------------------------------------------------------------
c convert thermodynamic equation of state from a pr tr reference state
c to get polynomial coefficients
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision g,s,v,a,b,c,d,e,f,gg,b1,b2,b3,b4,b5,b6,b7,b8,
     *                 b9,b10,b11,b12,tr,pr,n,v0,k00,k0p,gamma0,q0,
     *                 etas0,g0,g0p,r

      double precision emodu
      common/ cst318 /emodu(k15)
c----------------------------------------------------------------------
c                              this routine uses special characteristics
c                              of the data to identify the different types
c                              of EoS's in use within perple_x
      if (v.lt.0d0) then 
c                              negative "volume" signals one of the 
c                              stixrude EoS's, for these EoS's "s" is
c                              +/-n - number of atoms pfu.
        if (s.gt.0) then
c                              stixrude & bukowinski JGR '93
            n = s

         else
c                              stixrude & lithgow-bertelloni GJI '05
            n = -s

         end if 

         v0     = -v 
         k00    = a
         k0p    = b
         gamma0 = d
         q0     = e
         etas0  = f
         g0     =  emodu(1)
         g0p    =  emodu(2)            
c                                 nr9
         b1 = 9d0*n*r
c                                 c1
         b2 = 9d0*k00*v0
c                                 c2
         b3 = k0p/2d0-2d0
c                                 c3
         b4 = 3d0*b2*b3
c                                 aii
         b5 = 6d0*gamma0
c                                 aiikk
         b6 = -12d0*gamma0 + 36d0*gamma0**2 - 18d0*q0*gamma0
c                                 as 
         b7 = -(gamma0 + etas0)
c                                 aiikk2 
         b8 = b6/2d0
c                                 aii2   
         b9 = b5/2d0
c                                 nr9t0  
         b10 = b1*tr

         b11 = (3d0*k00*g0p-5d0*g0)

         b12 = ((6d0*g0p-24d0+4.5d0*k0p)*k00-14d0*g0)

      else if (b8.eq.0d0) then
c                              normal polynomial vdp term:
 
c                              first compute G*
         g  = g
c                              add the integral of sdt at tr,pr:
     *       + s * tr - a * tr - b * tr * tr / 2d0 + c / tr
     *       - d * tr**3 / 3d0 - 2d0 * e * dsqrt(tr)
     *       - f * dlog(tr) + gg / tr / tr / 2d0 + f
c                              add the constant components of the
c                              integral -vdp at (t,pr):
     *       - v * pr + b2 * tr * pr + b4 * pr * pr / 2d0
c                              sign on b7 corrected June 16, 2004. PJ Gorman. 
     *       - b6 * pr**3/3d0 - b7 * tr * tr * pr
c                              S* (T)
         s  = a - b2*pr - s + a*dlog(tr) + b*tr - c/tr/tr/2d0
     *        + d * tr * tr / 2d0 - 2d0 * e / dsqrt(tr) - f/tr
     *        - gg / tr**3 / 3.d0
c                              b7 term added June 16, 2004. PJ Gorman.
     *        + b7 * 2d0 * pr * tr
c                              V* (P)
         v  = v - b2 * tr - b4 * pr  
c                              b6 term added June 16, 2004. PJ Gorman.
     *          + b6 * pr * pr
c                              sign on b7 corrected June 16, 2004. PJ Gorman.
     *          + b7 * tr * tr
c                              a*  (-T log T)
c        a  = a
c                              b*  (-T*T)
         b  = b7 * pr + b / 2d0
c                              c*  (-1/T)
         c  = c / 2d0
c                              d*  (-T**3)
         d  = d / 6d0
c                              e*  (sqrt T)
         e  = 4d0 * e
c                              f*  (log T)
c        f  = f
c                              gg* (-1/T**2)
         gg = gg/6d0
c                              b2* (PT)
         b2 = b2
c                              b7 term added June 16, 2004. PJ Gorman.
     *          - 2d0 * b7 * tr 
c                              b3* ((P-Pr) exp (T / c1))
c        b3 = b3
c                              b4* (P**2)
         b4 = b4/2d0
c                              b6 term added June 16, 2004. PJ Gorman.
     *          - b6 * pr 
c                              b5* (exp (P/ c2))
c        b5 = c2 * b5
c                              b6* (P**3)
         b6 = b6/3d0
c                              b7* (P*T**2)
c        b7 = b7

      else
c                              thermal terms:
         g  = g
c                              add the integral of sdt at tr,pr:
     *       + s * tr - a * tr - b * tr * tr / 2d0 + c / tr
     *       - d * tr**3 / 3d0 - 2d0 * e * dsqrt(tr)
     *       - f * dlog(tr) + gg / tr / tr / 2d0 + f
c                              S* (T)
         s  = a - s + a * dlog(tr) + b * tr - c / tr / tr /2d0
     *        + d * tr * tr / 2d0 - 2d0 * e / dsqrt(tr) - f / tr
     *        - gg / tr**3 / 3d0
c                              these terms are as above
         b  = b / 2d0
         c  = c / 2d0
         d  = d / 6.d0
         e  = 4d0 * e
         gg = gg/6.d0
c                              volume terms:
c                              b8 is K' in the HP-murnaghan EoS
c                              b8 is -K' in the Ghiorso-murnaghan EoS
c                              or -beta in exponential volumetric
c                              functions.
         if (b8.gt.0d0) then
c                              HP-murnaghan vdp
            b1 = v*b1 
            b2 = v*b2/2d0 
            b3 = v*b3 
            b4 = -v*b4 
            b5 = 2d0*v*b5 
            v  = v - (b1-b2*tr)*tr - b3*dlog(tr) - b4/tr - b5*dsqrt(tr) 
c                                 these are unnecessary, but may cut time?
            b9 = (1d0-1d0/b8) 
            b10 = b8*pr 
            b11 = 1d0/b8 
            b12 = b8 - 1d0

         else if (b8.lt.-3d0) then
c                              Ghiorso Murnaghan



         else 
c                              gottschalk, b3 = alpha, -b8 = beta: 
            b1 = -v / b8 / dexp (b3*tr)
      
         end if 

      end if 
 
      end

      subroutine disord (tk,d0,d1,d2,d3,d4,d5,d6,tr,tmax,g)
c----------------------------------------------------------------------
c compute t-dependent disorder contribution, g is the
c gibbs energy of disordering, tr is the t of onset of
c disordering, tmax is the t of complete disorder.
c---------------------------------------------------------------------
      implicit none

      double precision tk,d0,d1,d2,d3,d4,d5,d6,tr,tmax,g,dh,tt,ds,gss

      double precision p,v2,v3,v4,v5,trr,pr,r,ps
      common/ cst5 /p,v2,v3,v4,v5,trr,pr,r,ps
 
      g = 0d0
      if (tk .lt. tr) return
      tt  = tk
      if (tk.gt.tmax) tt = tmax
 
      dh = d0 * (tt - tr)
     *     + 2d0 * d1 * (dsqrt(tt) - dsqrt(tr))
     *     - d2 * (1d0 / tt - 1d0 / tr)
     *     + d4 * dlog(tt/tr)
     *     + d5 * (tt*tt - tr*tr) / 2d0
     *     + d6 * (tt**3 - tr**3) / 3d0
 
      ds = d0 * dlog(tt/tr)
     *     - 2d0 * d1 * (tt**(-0.5d0) - tr**(-0.5d0))
     *     - d2 * (1d0/tt/tt - 1d0/tr/tr) / 2d0
     *     - d4 * (1d0/tt - 1d0 / tr)
     *     + d5 * (tt - tr)
     *     + d6 * (tt*tt - tr*tr) / 2d0
 
      gss = dh - (tk * ds)
      if (d3.ne.0d0) gss = gss + dh/d3 * (p - pr)
      end
 
      subroutine chkphi (ichk,name,good)
c-----------------------------------------------------------------------
c ichk = 0 and 2 -> test for saturated entities
c ichk = 1 and 3 -> test for non-saturated entities
c ichk > 1 -> do not compare against excluded list (for make definitions). 
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'
 
      character*8 name
 
      integer ichk,i,j

      logical good
 
      integer icomp,istct,iphct,icp
      common/ cst6 /icomp,istct,iphct,icp

      integer ic
      common/ cst42 /ic(k5)

      integer jfct,jmct,jprct
      common/ cst307 /jfct,jmct,jprct

      integer ids,isct,icp1,isat,io2
      common/ cst40 /ids(h5,h6),isct(h5),icp1,isat,io2

      integer ifct,idfl
      common/ cst208 /ifct,idfl

      integer ixct,iexyn,ifact
      common/ cst37 /ixct,iexyn,ifact 

      character*8 exname,afname
      common/ cst36 /exname(h8),afname(2)

      integer idh2o,idco2,ikind,icmpn,icout
      double precision therm,comp,tot
      common/ cst43 /therm(k4),comp(k0),tot,icout(k0),idh2o,idco2,
     *               ikind,icmpn
c-----------------------------------------------------------------------
      good = .true.
c                               reject the data if excluded in input:
      if (ichk.lt.2) then 
         do i = 1, ixct
            if (name.eq.exname(i)) goto 90
         end do
      end if 
c                               reject phases which have a component
c                               not included in the input:
      do i= 1, icmpn
         if (icout(i).eq.0.and.comp(i).ne.0d0) goto 90
      end do
c                               reject phases with negative/zero compositions
      tot = 0d0

      do j = 1, icmpn
         if (comp(j).lt.0d0.and.comp(j).gt.-1d-14) then
            comp(j) = 0d0
         else if (comp(j).lt.0d0) then 
c                               use ichk to avoid multiple messages
            if (ichk.eq.1) call warn (13,tot,j,name)
            goto 90
         else 
            tot = tot + comp(j)
         end if
      end do 

      if (tot.eq.0d0) goto 90 
c                               do a check to make sure that the phase does
c                               not consist of just mobile components
      if (jmct.gt.0) then

         tot = 0d0

         do j = 1, icp + isat + ifct
            tot = tot + comp(ic(j))
         end do  

         if (tot.eq.0d0) goto 90

      end if 
c                               the following is not executed by build:
c                               if ichk ne 0 reject phases that consist entirely
c                               of saturated components: 

      if (ichk.eq.0.or.ichk.eq.2) return
c                               reject phases with null composition, in case
c                               a user puts one in by accident
      tot = 0d0
      do j = 1, icp
         tot = tot + comp(ic(j))
      end do 

      if (tot.eq.0d0) goto 90

      goto 99

90    good = .false.
 
99    end

      logical function findph (igood)
c-----------------------------------------------------------------------
c findph checks if a phase loaded by getphi consists entirely of 
c component igood, if it does it returns true. 
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'
 
      integer igood, i

      integer idh2o,idco2,ikind,icmpn,icout
      double precision therm,comp,tot
      common/ cst43 /therm(k4),comp(k0),tot,icout(k0),idh2o,idco2,
     *               ikind,icmpn
c-----------------------------------------------------------------------

      if (comp(igood).ne.0d0) then 
c                                 the phase has the component,
c                                 does it have others?
         do i= 1, icmpn
            if (i.ne.igood.and.comp(i).ne.0d0) then
               findph = .false.
               return
            end if  
         end do

         findph = .true.

      else 

         findph = .false.

      end if 

      end
 
      subroutine lamlan (dg,rkt,ld,ltype)
c---------------------------------------------------------------------
c     calculate the extra energy of a lamdba transition using  the
c     Landau model (Holland and Powell).
 
c     input variables
 
c                        t = temperature in k
c                        p = pressure in bar
c                        ld = pointer to phase in therlm
c                        rkt = bulk modulus
 
c     returned - dg - difference in free energy due to the transition
c---------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer lct,ltype,i,ld

      double precision therdi,therlm,dg,tc,smax,q2,rkt
 
      common/ cst203 /therdi(m8,m9),therlm(m7,m6,k9)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      lct = ltype - 9
c                              don't do this stuff until TJBH
c                              has a consistent policy:

c                              create a constant which contains
c                              the 0 to Tr contribution to the
c                              free energy:
c               q = (1.-(tr/tc))**0.25
c               if (tr.gt.tc) q = 0.
c                              maximum enthalpy of disorder:
c               therlm (3,j,lamin) = 2./3.*smax*tc
c                              enthalpy of disorder at tr:
c               therlm (4,j,lamin) = 2.*smax*tc * 
c     *                              (q**6/6. - q*q/2. + 1./3.)
c                              entropy of disorder at tr:
c               therlm (5,j,lamin) = 

 
      dg = 0d0

      do i = 1, lct
 
         tc = therlm(1,i,ld) + therlm(3,i,ld)*(p-pr)
 
         smax = therlm (2,i,ld)

         if (t.gt.tc) then
c                                 fully disordered:
            dg = dg + therlm(4,i,ld) - t*therlm(5,i,ld)
 
         else
c                                 partially disordered:
            q2 = dsqrt(1d0-t/tc)
c                                 bulk modulus:
            dg = dg   
c                                 hr-t*sr
     *           + therlm(4,i,ld) - t*therlm(5,i,ld)   
c                                 G(landau)
     *           + smax * (q2*(t-tc) + tc*q2**3/3d0)
 
         end if
c                                 if pressure dependent
c                                 add the vdp term:
c                                 pressure dependent 
c                                 landau transitions are
c                                 only in the murnaghan 
c                                 versions of HP data.
         if (therlm(3,i,ld).ne.0d0) then
            dg = dg 
c                                 vt*kt/3
     *           + (therlm(6,i,ld) + therlm(7,i,ld)*t
     *              + therlm(8,i,ld)*dsqrt(t))*rkt/3d0
c                                 *[(1+4p/kt)^(3/4) - 1]
c                                 this approximation is only
c                                 good for pr near 0.
     *             *((1d0 + 4d0*p/rkt)**(0.75d0) - 1d0)
         end if 

 
      end do
 
      end

      double precision function gtrans (id,j)
c-----------------------------------------------------------------------
c gtrans computes the reference pressure free energy of a compound 
c aboves its jth transition. 
c---------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer id,j
 
      double precision therdi,therlm
      common/ cst203 /therdi(m8,m9),therlm(m7,m6,k9)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps
 
      gtrans = therlm(12,j,id)
c                                 -sdt
     *       + t * (therlm(3,j,id) - therlm(5,j,id) * dlog(t)
     *       - t * (therlm(6,j,id) + therlm(8,j,id) * t))
     *       - (therlm(7,j,id) + therlm(11,j,id) / t) / t
     *       + therlm(9,j,id) * dsqrt(t) + therlm(10,j,id)*dlog(t)

      end 
 
      subroutine lamhel (p,t,g,vdp,ld,ltype)
c---------------------------------------------------------------------
c     calculate the extra energy of a lamdba transition using model
c     of helgeson et al 1978 (AJS).
 
c     there is something seriously wrong in this routine!!!!
 
c     input variables
c                        p = pressure in bars
c                        t = temperature in k
c                        ld = pointer to phase in therlm
c                        g = initial free energy
 
c     returned - g - modified phase free energy
c---------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer lct,ltype,ld,i,jtran
      double precision t,g,gtrans,vdp,trtp,p,dt,pstar
 
      double precision therlm,therdi
      common/ cst203 /therdi(m8,m9),therlm(m7,m6,k9)

      double precision v,pr,u
      common/ cst5 /v(6),pr,u(2)
 
      lct = ltype - 3
c                                 T<T lowest transition, ignore
c                                 possibility of < clapeyron slope
c                                 and exit:
      if (t.lt.therlm(1,1,ld)) goto 99

      do i = 1, lct 
         if (t.lt.therlm(1,i,ld)) then 

            if (i.eq.1) then
c                                 T<T lowest transition, ignore
c                                 possibility of < clapeyron slope
c                                 and exit:
                goto 99
             else
c                                 at the i-1 th transition:
                jtran = i - 1
                exit
             end if 
          else 
             jtran = i 
          end if 
      end do 

      g = gtrans(ld,jtran) + vdp 
c                                 add vtrans*dp terms, this is
c                                 only set up for one transition
c                                 with a non-zero clapeyron.
c                                 should write warning. 
      if (therlm(2,1,ld).eq.0d0) goto 99

      trtp = therlm(1,1,ld) + (p-pr)/therlm(2,1,ld)
      dt = t - therlm(1,1,ld)

      if (t .gt. trtp) then
c                                 1 bar polymorph is stable at p,
c                                 p*thermo(3,id) is 0 bar standard
c                                 state pdv polymorph contribution &
c                                 therlm(4,j,ld) is Delta V of trans
c                                 contribution.
         pstar = dt*therlm(2,1,ld) + pr
 
         g = g + (p-pstar)*therlm(4,1,ld)
 
      else
c                                 the 1 bar polymorph isn't stable.
         g = g + dt * therlm(2,1,ld) * therlm(4,1,ld)
 
      end if
 
99    end
 
      subroutine lamqtz (p,t,g,du,ld,id)
c---------------------------------------------------------------------
c     calculate the extra energy of a lamdba transition using model
c     of helgeson et al 1978 (AJS). eq. (114) (corrected for the
c     misplaced parenthesis and the sign on c sub alpha), this is
C     probably bogus for coesite note that only one transition
c     is allowed.
 
c     input variables
c                        p = pressure in bars
c                        t = temperature in k
c                        ld = pointer to phase in therlm
c                        g = initial free energy
 
c     constants
c                Tr = 298.15 K
c                Pr = 1. bar
c                S(lope) = 38.5 bar/K
c                ba = 0.066 j/bar
c                aa = 549.82 K
c                ca = -4.973d-6 j/bar/bar
c                trt = 848. K
 
c     returned - g - modified phase free energy
c---------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer id,ld

      double precision p,t,ps,g,du,pdv
 
      double precision thermo,uf,us
      common/ cst1 /thermo(k4,k10),uf(2),us(h5)

      double precision therdi,therlm
      common/ cst203 /therdi(m8,m9),therlm(m7,m6,k9)

      double precision trt, tr, pr, s, aa, ba, ca , vb
 
      save trt, tr, pr, s, aa, ba, ca , vb
 
      data trt, tr, pr, s, aa, ba, ca, vb / 848., 298.15, 1d0,
     *                  38.5, 549.82, 0.066, -4.973d-6, 2.372 /
 
c      trtp = trt + (p-pr) / s
 
      if (t.gt.trt) then
         ps = (t-trt) * therlm(2,1,ld) + pr
      else
         ps = pr
      end if
c                                 if above the ref P transition T
c                                 do the cp integration:
      if (t.gt.trt) g = therlm(8,1,ld) + (p-ps) * thermo(3,id)
     *                - therlm(3,1,ld) * (t-trt) + du
     *     + therlm(5,1,ld) * (t - trt - t * dlog(t/trt))
     *     - (therlm(7,1,ld) + therlm(6,1,ld) * t * trt * trt)
     *     * (t - trt)**2 / 2d0 / t / trt / trt
 
c                                 now add in pdv terms to
c                                 the free energy, note that
c                                 pdv term for the ref polymorph
c                                 is already in:
 
      pdv  =  vb * (ps-pr)
     *      - ca * ( (2d0 * pr * (p-ps) - (p * p - ps * ps) ) / 2d0
     *              + s * (t-tr) * (p-ps) )
     *      + s * (ba + aa*ca*s) * (t-tr)
     *                           * dlog ((aa+p/s) / (aa + ps/s))
 
      g = g + pdv
 
      end
 
      subroutine lamubc (p,t,gspk,k,ltype)
c---------------------------------------------------------------------
c     calculate the extra energy of a lamdba transition using model
c     of berman and brown (1985, contribs. miner. petro.)
 
c     input variables
c                        p = pressure in bars
c                        t = temperature in k
c                        k = pointer to phase in therlm
 
c     returned - gspk - energy because of lambda spike
 
c     modified from perkins et al 1986 (computers and geosciences).
c---------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer j,ltype,k

      double precision gspk,ctrans,aspk2,bspk2,ct2,ct3,a1,b1,c1,t92,t93,
     *                 tr92,tr93,dhspk,dsspk,t9,tr,teq,tq1bar,p,tr9,
     *                 dvdtr,dvdp,dstr,abspk,t
 
      double precision therdi, therlm
      common/ cst203 /therdi(m8,m9),therlm(m7,m6,k9)
c---------------------------------------------------------------------
      gspk=0d0
 
      do j = 1, ltype
 
         tq1bar = therlm(3,j,k)
 
         if (tq1bar.eq.0d0) cycle
 
            tr    = therlm(7,j,k)
            teq   = therlm(4,j,k) * (p-1d0) + tq1bar
            ctrans = tq1bar - teq
            tr9    = tr - ctrans
 
            if (t.lt.tr9) cycle 
 
            aspk2 = therlm(1,j,k)
            bspk2 = therlm(2,j,k)
            dvdtr = therlm(5,j,k)
            dvdp  = therlm(6,j,k)
            dstr  = therlm(8,j,k) / tq1bar
            abspk = therlm(9,j,k)
 
            if (t .gt. teq) then
               t9 = teq
            else
               t9 = t
            end if
 
            ct2 = ctrans * ctrans
            ct3 = ct2 * ctrans
 
            a1 = aspk2 * ctrans + 2d0 * abspk * ct2 + bspk2 * ct3
            b1 = aspk2 + 4d0 * abspk * ctrans + 3d0 * bspk2 * ct2
            c1 = 2d0 * abspk + 3d0 * ctrans * bspk2
 
            t92 = t9 * t9
            t93 = t9 *t92
            tr92 = tr9 * tr9
            tr93 = tr92 * tr9
 
            dhspk = a1 * (t9 - tr9)
     *            + b1 * (t92 - tr92) / 2d0
     *            + c1 * (t93 - tr93) / 3d0
     *            + bspk2 * (t9*t93 - tr93*tr9) / 4d0
 
            dsspk = a1 * (dlog(t9) - dlog(tr9))
     *            + b1 * (t9 - tr9)
     *            + c1 * (t92 - tr92) / 2d0
     *            + bspk2 * (t93 - tr93) / 3d0
 
            gspk = gspk - (t9 * dsspk) + dhspk
 
            if (t.gt.teq) gspk = gspk - (dstr + dsspk ) * (t-teq)
 
            gspk = gspk + dvdtr * (p-1d0) * (t9-tr)
     *                  + dvdp * ((p*p-1d0)/2d0 - (p-1d0))
 
      end do 

      end
 
      subroutine loadit (id)
c---------------------------------------------------------------------
c loadit loads descriptive data for phases and species (name,comp,
c and therm) into the appropriate arrays (names,comps,thermo,vf,
c and vs).  the arguement 'id' indexes the phase in the arrays.
c note that loadit also computes reference state constants which
c are dependent on the state function being used and on its
c analytical expression.
c---------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'
 
      character*8 names, name

      integer ipoint,imyn,iff,
     *        idss,ifug,ifyn,isyn,icomp,istct,iphct,icp,
     *        ltyp,lmda,idis,id,i,j,k,lct

      double precision z(13),thermo,uf,us,cp,therdi,therlm,
     *                 smax,t0,a0,qr2,qr4,vmax,gph
 
      common/ cst60  /ipoint,imyn
     *      / cst1   /thermo(k4,k10),uf(2),us(h5)
     *      / cst10  /iff(2),idss(h5),ifug,ifyn,isyn
     *      / cst8   /names(k1)
     *      / cst12 /cp(k5,k1)
     *      / cst6   /icomp,istct,iphct,icp/ csta6 /name
     *      / cst203 /therdi(m8,m9),therlm(m7,m6,k9)
     *      / cst204 /ltyp(k1),lmda(k1),idis(k1)

      double precision emodu
      common/ cst318 /emodu(k15)

      integer ic
      common/ cst42 /ic(k5)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      integer ilam,idiso,lamin,idsin
      double precision tm,td
      common/ cst202 /tm(m7,m6),td(m8),ilam,idiso,lamin,idsin

      integer idh2o,idco2,ikind,icmpn,icout
      double precision therm,comp,tot
      common/ cst43 /therm(k4),comp(k0),tot,icout(k0),idh2o,idco2,
     *               ikind,icmpn

      integer ikp
      common/ cst61 /ikp(k1)

      integer iemod,kmod
      logical smod
      double precision emod
      common/ cst319 /emod(k15,k10),smod(h9),iemod(k10),kmod

      integer ixp,ifp
      double precision sxs,exces
      common/ cst304 /sxs(k13),exces(m3,k1),ixp(k1),ifp(k1)

      logical gflu,aflu,fluid,shear,lflu,volume
      common/ cxt20 /gflu,aflu,fluid(k5),shear,lflu,volume

      integer isfp
      common/ cst303 /isfp(k10)

      character dname*40,cmpnt*5
      common/ csta5 /dname,cmpnt(k0)

      integer ifct,idfl
      common/ cst208 /ifct,idfl

      integer junk
      double precision del, rand
      common/ cst321 /del(11,k10),rand(12,k10),junk

      double precision delta
      common/ cst325 /delta(11)

      double precision vnumu
      common/ cst44 /vnumu(i6,k10)

      integer jfct,jmct,jprct
      common/ cst307 /jfct,jmct,jprct
c---------------------------------------------------------------------

      if (id.gt.k10) call error (1,0d0,id,'k10')

      ipoint = iphct
      
c                               load name and phase flag
      names(id) = name
c                               shear modulus flag
      if (therm(3).lt.0) then 
         iemod(id) = 2
      else 
c                               this is gonna be bad news
c                               if anyone ever take the 
c                               seimic velocity tutorial too
c                               earnestly.
         iemod(id) = ikind
      end if 

      ikp(id) = 0
      ifp(id) = 0
      isfp(id) = 0 

      if (name.eq.cmpnt(idh2o)) then 
c                                 set fluid flag for gcpd, this flag
c                                 is used to indicate the endmember is
c                                 a fluid and the id of the fluid component.
         if (ifyn.ne.0) then 
c                                 no saturated phase
            isfp(id) = 1

         else if (idfl.ne.2) then 
c                                 saturated phase, and it's not CO2, ergo
c                                 will be computed by ufluid.
            isfp(id) = 0 

         end if

c                                 set fluid flag, this flag is 
c                                 used only to match fluid endmembers
c                                 with fluid pseudocompounds
         ifp(id) = 1  
c                                 gflu used to indicate whether a fluid is 
c                                 in the calculation.  
         if (ifyn.eq.1) gflu = .true. 

      else if (name.eq.cmpnt(idco2)) then 

         if (ifyn.ne.0) then 
c                                 no saturated phase
            isfp(id) = 2

         else if (idfl.ne.1) then 
c                                 saturated phase, and it's not H2O, ergo
c                                 will be computed by ufluid.
            isfp(id) = 0 

         end if

         ifp(id) = 1    
         if (ifyn.eq.1) gflu = .true. 

      end if 
c                               load stoichiometry of components.
      do i = 1, icomp
         cp(i,id) = comp(ic(i))
      end do 
c                               and just mobile components
      do i = 1, jmct 
         vnumu(i,id) = comp(ic(i+jprct))
      end do 
c                               load elastic props if present
      if (iemod(id).ne.0) then 
c                               kmod initialized 0 in main.
         kmod = 1        
         do i = 1, k15
            emod(i,id) = emodu(i)
         end do

      end if 
c                               compute reference state constants,
c                               etc..
      do i = 1, k4
         thermo(i,id) = therm(i)
      end do 
c                               load errors for MC calculations
      do i = 1, 11
         del(i,id) = delta(i)
      end do 
 
      call conver (thermo(1,id),thermo(2,id),thermo(3,id),thermo(4,id),
     *             thermo(5,id),thermo(6,id),thermo(7,id),thermo(8,id),
     *             thermo(9,id),thermo(10,id),thermo(11,id),
     *             thermo(12,id),thermo(13,id),thermo(14,id),
     *             thermo(15,id),thermo(16,id),thermo(17,id),
     *             thermo(18,id),thermo(19,id),thermo(20,id),
     *             thermo(21,id),thermo(22,id),tr,pr,r)
 
      if (tr.eq.0d0) then
         thermo(1,id) = therm(1)
         thermo(2,id) = therm(2)
      end if
c                              lmda transitions:
      idis(id) = 0
      lmda(id) = 0
      if (ilam.ne.0) then

         lamin=lamin + 1
         if (lamin.gt.k9) call error (1,0d0,lamin,'k9')
 
         if (ilam.ge.10.and.ilam.le.12) then
c                              holland and powell, landau model:
            lct = ilam - 9
            do j = 1, lct       

               smax = tm(2,j)
               t0 = tm(1,j)
               vmax = tm(3,j)
               a0 = therm(11)

               therlm(1,j,lamin) = t0
               therlm(2,j,lamin) = smax
c                              this makes therlm(3) dt/dp
               therlm(3,j,lamin) = tm(3,j)/tm(2,j)
c                              qr2,qr4
               
               qr2 = dsqrt (1d0 - tr/t0)
               qr4 = (1d0 - tr/t0)**(1.5d0)/3.d0
c                              href:
               therlm(4,j,lamin) = smax*t0*(qr2-qr4)
c                              sref:
               therlm(5,j,lamin) = smax*qr2
c                              vT coefficients:  
c                              constant:
               therlm(6,j,lamin) = vmax*qr2*
     *                             (1d0-a0*tr+20d0*a0*dsqrt(tr))
c                              *T coefficient:
               therlm(7,j,lamin) = vmax*qr2*a0
c                              *sqrt(T) coeff:
               therlm(8,j,lamin) = -20d0*vmax*qr2*a0
c                              TJBH said that he is not yet
c                              consistent in correcting the 
c                              reference state props for the
c                              landau effect, accordingly I
c                              assume here that the effect is
c                              not incorporated. 

c                              this stuff is no longer possible with
c                              hollands variable tc formulation:

c                              create a constant which contains
c                              the 0 to Tr contribution to the
c                              free energy:
c               q = (1.-(tr/tc))**0.25
c               if (tr.gt.tc) q = 0.
c                              maximum enthalpy of disorder:
c               therlm (3,j,lamin) = 2./3.*smax*tc
c                              enthalpy of disorder at tr:
c               therlm (4,j,lamin) = 2.*smax*tc * 
c     *                              (q**6/6. - q*q/2. + 1./3.)
c                              entropy of disorder at tr:
c               therlm (5,j,lamin) = smax * (1. - q*q)
            end do 

 
 
         else if (ilam.le.3) then
c                              ubc:
            do j = 1, ilam
 
               therlm(1,j,lamin)=tm(1,j)*tm(1,j)
               therlm(2,j,lamin)=tm(2,j)*tm(2,j)
               therlm(9,j,lamin)=tm(1,j)*tm(2,j)
 
               do k = 3, 8
                  therlm(k,j,lamin)=tm(k,j)
               end do 
            end do 
 
         else if (ilam.gt.3.and.ilam.lt.8) then
c                              helgeson:
            p = pr
c                              now convert paramters:
            do k = 1, ilam - 3
c                              load into therlm:
               therlm(1,k,lamin) = tm(1,k)
               therlm(2,k,lamin) = tm(2,k)

               do j = 4, 10
                  therlm(j+1,k,lamin) = tm(j,k)
               end do 

               t = tm(1,k)
               if (k.eq.1) then
c                              set transition type to null
c                              for call to gphase
                  lmda(id) = 0
               else 
                  lmda(id) = lamin
                  ltyp(id) = 2 + k
               end if
c                              g at trt:
               call gphase (id,gph)
               therlm(12,k,lamin) = gph
c                              delta v trans:
               therlm(4,k,lamin) = 0d0
               if (tm(2,k).ne.0d0) therlm(4,k,lamin) = tm(3,k)/tm(2,k)
c                              s + st at trt:
               t = t + 1d-3
               call gphase (id,gph)
               therlm(3,k,lamin) = tm(3,k) -
     *                             (gph - therlm(12,k,lamin))/1d-3
c                              streamline the eos:
               do j = 1, 13
                  z(j) = 0d0 
               end do 

               call conver (therlm(12,k,lamin),therlm(3,k,lamin),
     *          z(1),therlm(5,k,lamin),therlm(6,k,lamin),
     *          therlm(7,k,lamin),therlm(8,k,lamin),therlm(9,k,lamin),
     *          therlm(10,k,lamin),therlm(11,k,lamin),
     *          z(2),z(3),z(4),z(5),z(6),z(7),z(8),
     *          z(9),z(10),z(11),z(12),z(13),tm(1,k),pr,r)
            end do 

         end if
 
         lmda(id)=lamin
         ltyp(id)=ilam
 
      end if
c                              t dependent order: load berman and brown
c                              parameters, (this should be cleaned up)
      if (idiso.ne.0) then
         idsin = idsin + 1
         idis(id) = idsin
         do j = 1, m8
            therdi(j,idsin) = td(j)
         end do 
      end if
 
      end

      subroutine univeq (i,ier)
c---------------------------------------------------------------------
c univeq computes the equilibrium condition of a univariant reaction.
c using a finite difference newton-raphson approximation. the method
c will fail when the derivative of the state function, with respect to
c the state variable v(i), goes to infinity.
c---------------------------------------------------------------------
      implicit none

      double precision v,tr,pr,r,ps,delt,dtol,utol,ptol,
     *                 vi,delv,del,u,gr,ddv,vmax,vmin,b

      integer i,j,ier
 
      include 'perplex_parameters.h'

      double precision blim, ulim, dgr
      common/ cxt62 /blim(l2),ulim(l2),dgr
 
      common/ cst5  /v(l2),tr,pr,r,ps
     *      / cst87 /delt(l2),dtol,utol,ptol
     *      / cst9 /vmax(l2),vmin(l2),ddv(l2)
 
      ier = 0
 
      vi = v(i)
      del = delt(i)
      b = blim(i)
      u = ulim(i)
c                                 phase composition special case:
      if (i.eq.3) then 
         if (vi.lt.1d1*del) then 
            del = dabs(vi)/1d1
         else if (1d0-vi.lt.1d1*del) then 
            del = dabs(1d0-vi)/1d1
         end if 
      end if 
 
      if (vi+dabs(del).gt.u.or.vi-dabs(del).lt.b) goto 30
 
      do j = 1, 100

         call grxn (gr)

         v(i) = vi + del
         call incdep (i)

         call grxn (dgr) 
         dgr = dgr - gr

         if (dgr.eq.0d0) exit
 
         delv = gr*del/dgr

         if (delv/ddv(i).gt.1d0) then 
            v(i) = vi - del
            del = del*ddv(i)/delv/2d0
            if (del/delt(i).lt.1d-6) goto 30 
            cycle
         end if 

         vi = vi - delv

         if (vi+dabs(del).gt.u.or.vi-dabs(del).lt.b) goto 30

         v(i) = vi

         call incdep (i)
c                                 used to be on the basis of utol,
c                                 but this allows uniform relative error.
         if (dabs(delv).lt.del) goto 99

      end do 
c                                 error code 1: iv and dv must be
c                                 switched
      ier = 1
      return
c                                 error code 2: value too far out of
c                                 range, refine iv increment.
30    ier = 2
99    end

      subroutine unver (g,s,v,a,b,c,d,e,f,gg,
     *                  b1,b2,b3,b4,b5,b6,b7,b8,tr,pr)
c---------------------------------------------------------------------
c convert thermodynamic equation of state from a 0-0 reference state
c to a pr-tr reference state.

c corrections made corresponding to PJ's corrections in conver, 
c June 16, 2004.
c---------------------------------------------------------------------
 
      implicit none

      double precision v,gg,a,b,c,d,e,f,g,b1,b2,b3,b4,b5,b6,b7,b8,pr,
     *                 tr,s

c double precision (a-g,o-z),integer (h-n)
c                               stixrudes EoS, exit without
c                               doing anything
      if (v.lt.0d0) goto 99

      gg = 6d0 * gg
      e  = e / 4d0
      d  = 6d0 * d
      c  = 2d0 * c

      if (b8.eq.0d0) then 
c                                normal vdp term:
         b6 = 3d0 * b6
         b4 = 2d0 * (b4 + b6 * pr)
         b2 = b2 + 2d0 * b7 * tr
         b  = 2d0 * (b - b7*pr)
         v  = v + b2 * tr + b4 * pr - b6 * pr * pr - b7 * tr * tr
         s  = -1d0 * ( s - ( a - b2 * pr  + a * dlog(tr) + b * tr
     *        - c / tr / tr /2d0
     *        + d * tr * tr / 2d0 - 2d0 * e / dsqrt(tr) - f / tr
     *        - gg / tr**3 / 3d0 + b7 * 2d0 * pr * tr ) )
         g  = g - ( s * tr - a * tr - b * tr * tr / 2d0 + c / tr
     *        - d * tr**3 / 3d0 - 2d0 * e * dsqrt(tr)
     *        - f * dlog(tr) + gg / tr / tr / 2d0 + f
     *        - v * pr + b2 * tr * pr + b4 * pr * pr / 2d0 
     *        - b6 * pr**3/3d0 - b7 * tr * tr * pr )
      else 

         b  = 2d0 * b

         s  = -1d0 * ( s - ( a + a * dlog(tr) + b * tr
     *        - c / tr / tr /2d0
     *        + d * tr * tr / 2d0 - 2d0 * e / dsqrt(tr) - f / tr
     *        - gg / tr**3 / 3d0 ) )
         g  = g - ( s * tr - a * tr - b * tr * tr / 2d0 + c / tr
     *        - d * tr**3 / 3d0 - 2d0 * e * dsqrt(tr)
     *        - f * dlog(tr) + gg / tr / tr / 2d0 + f)

         if (b8.gt.0d0) then 
c                                 b8 is K' 
c                                 HP-murnaghan vdp:
            v  = v + (b1+b2*tr)*tr + b3*dlog(tr) + b4/tr + b5*dsqrt(tr)
            b1 = b1 / v
            b2 = 2d0 * b2/v
            b3 = b3 / v
            b4 = -b4 / v
            b5 = b5 / v / 2d0  

         else if (b8.lt.-3d0) then
c                                 b8 is -K' 
c                                 Ghirso-murnaghan

         else 
c                                 v = f(exponential(beta))
c                                 zero b1, so users don't get confused.
            b1 = 0d0
         end if  

      end if 
 
99    end
 
c-----------------------------------------------------------------------
 
c PLIB -  subroutines common to FRENDLY and VERTEX.
 
c-----------------------------------------------------------------------
 
      subroutine slope (iv1,iv2,s)
c---------------------------------------------------------------------
c slope computes the slope (div(1)/div(2)) of a univariant
c equilibria by approximating the derivatives dg/div(1) and
c dg/div(2) by finite differences.
c---------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'
 
      integer iv(2),iv1,iv2,i

      double precision dg(2),v,tr,pr,r,ps,delt,dtol,utol,ptol,gr,s,gval
 
      common/ cst5 /v(l2),tr,pr,r,ps/ cst87 /delt(l2),dtol,utol,ptol
 
      iv(1) = iv1
      iv(2) = iv2
      call grxn (gr)

      do i = 1, 2
         v(iv(i)) = v(iv(i)) + delt(iv(i))
         call incdep (iv(i))
         call grxn (gval)
         dg(i) = (gval-gr)/delt(iv(i))
c                                 note the possibility exists here
c                                 for the value of an intensive
c                                 parameter to excede allowed
c                                 limits.  although this is very
c                                 unlikely a limit test like the
c                                 one done in 'univeq' could be
c                                 used.
         v(iv(i)) = v(iv(i)) - delt(iv(i))
         call incdep (iv(i))
      end do 

      s = -dg(2) / dg(1)
 
      end
 
      subroutine switch (div,ivi,ivd,jer)
c---------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'
 
      integer ivi,ivd,jer,iovd

      double precision div,vo,s,odiv

      double precision v,tr,pr,r,ps,vmax,vmin,dv
      common/ cst5 /v(l2),tr,pr,r,ps/ cst9 /vmax(l2),vmin(l2),dv(l2)
c---------------------------------------------------------------------
c                                 reset intensive variables
      call reptx 
      vo = v(ivi)
c                                 determine the sign for the
c                                 ivd increment
      call slope(ivd,ivi,s)
      jer = 0
c
      if (s.eq.0d0) then
c                                 a zero slope shouldn't occur
c                                 someplace is a bug
         jer = 1
         return
      end if
c
      odiv = div
      div = dv(ivd)
      if (odiv.gt.0d0) goto 10
      if (s.gt.0d0) div = -dv(ivd)
      goto 20
10    if (s.lt.0d0) div = -dv(ivd)
c                                 estimate a new value for v(ivi)
20    v(ivi) = v(ivi)+div/s
      call incdep (ivi) 
c                                 switch variables
      if ((v(ivi).gt.vmin(ivi)).and.(v(ivi).lt.vmax(ivi))) goto 30
      div = div/5d0

      if (dabs(div).lt.dv(ivd)/1d6) then
         jer = 1
         return
      end if 

      v(ivi) = vo
      call incdep (ivi) 
      goto 20

30    iovd = ivd
      ivd = ivi
      ivi = iovd

      end

      subroutine unblnk (text)
c------------------------------------------------------------------- 
c unblnk - subroutine to remove blanks from text
 
c     text - character string 
c     jchar - length of unblanked character string, 0 on input 
c             if unknown.
c-------------------------------------------------------------------
      implicit none

      integer i,ict,nchar

      character text*(*), bitsy(400)*1 

      nchar = len(text)
 
      read (text,1000) (bitsy(i), i = 1, nchar)
c                                 scan for blanks:
      ict = 0

      do i = 1, nchar
         if (bitsy(i).eq.' ') cycle 
         ict = ict + 1
         bitsy(ict) = bitsy(i)
      end do 

      write (text,1000) (bitsy(i), i = 1, ict)
 
1000  format (400a1)
      end

      subroutine enblnk (text)
c---------------------------------------------------------------------- 
c enblnk - scan text to first blank and cut remaining text.
 
c     text - character string 
c     jchar - length of unblanked character string, 0 on input 
c             if unknown.
c----------------------------------------------------------------------
      implicit none

      integer i,ict,nchar
 
      character text*(*), bitsy(400)*1 

      nchar = len(text) 
      read (text,1000) (bitsy(i), i = 1, nchar)
c                                 scan for blanks:
      ict = 0

      do i = 1, nchar
         if (bitsy(i).eq.' ') exit
         ict = ict + 1
      end do 

      write (text,1000) (bitsy(i), i = 1, ict),(' ',i = ict+1, nchar)
 
1000  format (400a1)
      end

      subroutine mertxt (text,text1,text2)
c------------------------------------------------------------------- 
c mertxt - merge text strings cutting trailing blanks, if text1 eq 
c          blank then pad by an initial 40 blanks.
 
c     text - character string 
c-------------------------------------------------------------------
      implicit none

      integer i, nchar1, nchar2, nchar3
 
      character text*(*), text1*(*), text2*(*), bitsy(400)*1 

      if (text1.eq.' ') then
         write (text,1010) text2
      else 
         nchar1 = len(text1) 
         read (text1,1000) (bitsy(i), i = 1, nchar1)  
         do i = 2, nchar1
            if (bitsy(i).eq.' '.and.bitsy(i-1).eq.' ') then 
               nchar1 = i - 1
               exit 
            end if 
         end do 
c                                 put 5 blanks between 1st and 
c                                 2nd strings
         do i = nchar1+1, nchar1 + 6
            bitsy(i) = ' '
         end do 

         nchar1 = nchar1 + 5 
         nchar2 = nchar1 + len(text2) 
         nchar3 = len(text)

         if (nchar2.gt.nchar3) then 
            nchar2 = nchar3 - nchar1
         end if 

         if (nchar2.gt.400) then 
            nchar2 = 400 - nchar1
         end if 

         read (text2,1000) (bitsy(i), i = nchar1, nchar2)  
         write (text,1000) (bitsy(i), i = 1, nchar2),
     *                     (' ',i = nchar2+1, nchar3)
      end if 

1010  format (40x,a)
1000  format (400a1)

      end

      subroutine reblnk (text)
c------------------------------------------------------------------- 
c reblnk - subroutine to replace blanks followed by a character
c          with an underscore.
 
c     text - character string 
c     jchar - length of unblanked character string, 0 on input 
c             if unknown.
c----------------------------------------------------------------------
      implicit none 

      integer i,ict
 
      character*8 text, bitsy(8)*1 
 
      read (text,1000) bitsy
c                                 scan for blanks:
      ict = 0

      do i = 1, 7

         if (i.eq.1.and.bitsy(i).eq.' ') cycle

         if (bitsy(i).eq.' '.and.bitsy(i+1).ne.' ') then
            ict = ict + 1
            bitsy(ict) = '_'
         else if (bitsy(i).eq.' ') then 
            cycle
         else
            ict = ict + 1
            bitsy(ict) = bitsy(i)
         end if

      end do 
 
      bitsy(ict+1) = bitsy(8)

      write (text,1000) (bitsy(i), i = 1, ict + 1)
 
1000  format (8a1)
      end


      subroutine reptx
c-----------------------------------------------------------------------
c reset - resets indendent potentials to the last known stable condition 
c along a univariant curve 
c----------------------------------------------------------------------- 
      implicit none

      double precision v,tr,pr,r,ps,ptx

      integer ipt2

      include 'perplex_parameters.h'
    
      common/ cst5  /v(l2),tr,pr,r,ps/ cst32 /ptx(l5),ipt2

      integer ipot,jv,iv
      common/ cst24 /ipot,jv(l2),iv(l2)

      v(iv(1)) = ptx(ipt2-1)
      v(iv(2)) = ptx(ipt2)
      call incdp0
   
      end

      subroutine toop (id,gex)
c-----------------------------------------------------------------------
c giulio, here is a dummy routine to serve as a template, be careful
c of the type conventions that i use (i.e., all reals are implicit
c double) and also not to change any of the variables in the common
c blocks. you have to keep the three common blocks cts6, cst5 and cst12
c for input, the variables that will concern you are 

c on input:
c           r - the gas constant, J/mol-K
c           t - temperature (K)
c           p - pressure (bar)
c           tr - reference t, as in the thermodynamic data file header.
c           pt - reference p, dido.
c           icp - the number of components.
c           id - a pointer to the phase composition
c           cp(i,id) - the mol fraction of the ith component in composition id.
c on output:
c           gex - the excess free energy in J (i.e., in excess of a mechanical
c                 mixture of the endmembers) per mol of components.

c N.B. the components appear in the same order that they are entered in 
c the input file and this must also be identical to the order in the 
c solution model data file. i think build will always make SIO2 the 
c first component, overriding your choice. so to be safe anywhere you
c enter components expect SIO2 to be the 1st composition, then enter 
c the remaining components in order of the cation atomic weight, e.g., 
c cp(1,id) will always be X(SIO2) in the id'th composition, cp(2,id) will
c be X(CA0) if there are no other components.
c----------------------------------------------------------------------- 
      implicit none

      include 'perplex_parameters.h'

      integer i,id
     
      double precision gex

      double precision cp
      common/ cst12 /cp(k5,k1)

      integer icomp,istct,iphct,icp
      common/ cst6 /icomp,istct,iphct,icp

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      gex = 0d0

      do i = 1, icp
         if (cp(i,id).ne.0d0) gex = gex + r*t*cp(i,id)*log(cp(i,id))
      end do 

      end 
         
      subroutine unlam (tm,id,lct)
c---------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer ilam,id,jd,lct,i,j,k
  
      double precision tm(m7,m6),z(9),g1,gph,g0,s0

      double precision therdi, therlm
      common/ cst203 /therdi(m8,m9),therlm(m7,m6,k9)

      integer ltyp,lmda,idis
      common/ cst204 /ltyp(k1),lmda(k1),idis(k1)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps
c-----------------------------------------------------------------------
                                 
      ilam = ltyp(id)
      jd = lmda(id)
      lct = 0 

      do i = 1, m7
         do j = 1, m6
            tm(i,j) = 0d0
         end do 
      end do 

      if (ilam.ne.0.and.(ilam.lt.4.or.ilam.gt.7)) then
   
         ilam = ltyp(id)
         jd = lmda(id)
         if (ilam.le.3) then
            lct = ilam
            do j = 1, ilam
               tm(1,j) = dsqrt (therlm(1,j,jd))
               tm(2,j) = dsqrt (therlm(2,j,jd)) 
            end do         
         else if (ilam.gt.9.and.ilam.lt.13) then
            lct = ilam - 9
            do j = 1, ilam - 9
               tm(1,j) = therlm(1,j,jd)
               tm(2,j) = therlm(2,j,jd)
               tm(3,j) = therlm(3,j,jd) *therlm(2,j,jd) 
            end do 
         end if
      end if

      if (ilam.gt.3.and.ilam.lt.7) then

         lct = ilam - 3
         p = pr
         do i = lct, 1 , -1
c                              get s transition:
c                              load into therlm:
            tm(1,i) = therlm(1,i,jd) 
            tm(2,i) = therlm(2,i,jd) 

            do j = 5, 11
               tm(j-1,i) = therlm(j,i,jd)  
            end do 
            t = tm(1,i)

            if (i.eq.1) then
c                              set transition type to null
c                              for call to gphase 
               ltyp(id) = 0
            else 
               ltyp(id) = 2 + i
            end if
c    c                         -s at trt:

            call gphase (id,g1)
            t = t + 1d-3

            call gphase (id,gph)
            tm(3,i) =  (gph - g1)/1d-3

            g0 = therlm(12,i,jd)
            s0 = therlm(3,i,jd)
            do k = 1, 9
               z(k) = 0d0
            end do 

            call unver (g0,s0,z(1),tm(4,i),tm(5,i),tm(6,i),tm(7,i),
     *                  tm(8,i),tm(9,i),tm(10,i),z(2),z(3),z(4),z(5),
     *                  z(6),z(7),z(8),z(9),tm(1,i),pr) 

            tm(3,i) = s0 + tm(3,i)
         end do 

         ltyp(id) = ilam
         lmda(id) = jd

      end if
      
      end 

      subroutine outlam (tm,jd,lct)
c---------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'
      
      integer jd,lct,i,j
      double precision tm(m7,m6)

      double precision therdi,therlm
      common/ cst203 /therdi(m8,m9),therlm(m7,m6,k9)
c-----------------------------------------------------------------------
      if (lct.ne.0) write (n2,1000) ((tm(i,j),i=1,10),j= 1,lct)
      if (jd.ne.0) write (n2,1000) (therdi(j,jd), j = 1, m8)
1000  format (5(g13.7,1x))
      end 
 
      subroutine incdep (ind) 
c-----------------------------------------------------------------------
c either indep or incdp0 are called whenever any primary potential
c variables are changed to reevaluate secondary variables. this 
c cumbersome structure is necessitated by the fact that computational
c variables were mapped directly to the thermodynamic variables in 
c cst5. a more rational strategy, which is used in werami/pssect, is
c to separate the computational and thermodynamic variables.

c incdep:

c  1) if ind = iind, computes the dependent variable if P(T) or T(P)

c  ind  - is the index of the current variable. 
c  iind - is the index of the independent variable that v(idep) is a 
c         function of.
c  idep - is the index of the dependent variable.

c  2) if jmct > 0, assigns the independent chemical potentials
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ind
       
      double precision var

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      integer iind, idep
      double precision c0,c1,c2,c3,c4,c5
      common/ cst316 /c0,c1,c2,c3,c4,c5,iind,idep

      integer jfct,jmct,jprct
      common/ cst307 /jfct,jmct,jprct
c----------------------------------------------------------------------
      if (ind.eq.iind) then 
         var = v(iind)
         v(idep) = c0 + var*(c1 + var*(c2 
     *                          + var*(c3 + var*c4))) 
c    *                           var*(c3 + var*(c4 + c5*var)))) 
      end if 

      if (jmct.gt.0) call subinc 

      end 

      subroutine subinc 
c-----------------------------------------------------------------------
c assigns the independent chemical potentials, called by incdep and 
c incdp0
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,id
       
      double precision gref, xp

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      integer jfct,jmct,jprct
      common/ cst307 /jfct,jmct,jprct

      integer imaf,idaf
      common/ cst33 /imaf(i6),idaf(i6)

      double precision mu
      common/ cst39 /mu(i6)
c----------------------------------------------------------------------
      do i = 1, jmct

            if (imaf(i).eq.1) then 
c                                 the primary variable is a chemical
c                                 potential.
               mu(i) = v(3+i)

            else 
c                                 the primary variable is a fugacity or
c                                 an activity.
               if (imaf(i).eq.2) then 
c                                 fugacity
                  xp = v(1)
                  v(1) = pr
                  call gcpd (idaf(i),gref)
                  v(1) = xp

               else 
c                                 activity
                  call gcpd (idaf(i),gref)

               end if 

               mu(i) = gref + r*v(2)*v(3+i)*2.302585093d0

             end if 

      end do 
c                                 this line is here to prevent an optimization 
c                                 bug with compaq visual fortran. 
      xp = xp * xp

      end 

      subroutine incdp0 
c----------------------------------------------------------------------
c incdep 1) conditionally computes the dependent variable if one exists 
c (idep ne 0), 2) assigns independent chemical potentials.

c  idep - is the index of the dependent variable 
c  iind - is the index of the independent variable that v(idep) is a 
c         function of.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'
       
      double precision var 

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      integer iind, idep
      double precision c0,c1,c2,c3,c4,c5
      common/ cst316 /c0,c1,c2,c3,c4,c5,iind,idep

      integer jfct,jmct,jprct
      common/ cst307 /jfct,jmct,jprct
c----------------------------------------------------------------------
      if (idep.ne.0) then 
         var = v(iind)
         v(idep) = c0 + var*(c1 + var*(c2 
     *                          + var*(c3 + var*c4))) 
c    *                            var*(c3 + var*(c4 + c5*var)))) 
      end if 

      if (jmct.gt.0) call subinc

      end 

      double precision function depvar (var) 
c--------------------------------------------------------------------
c depvar computes the dependent variable from the independent variable
c var
c--------------------------------------------------------------------
      implicit none

      double precision c0,c1,c2,c3,c4,c5,var

      integer iind,idep

      common/ cst316 /c0,c1,c2,c3,c4,c5,iind,idep

      depvar = c0 + var*(c1 + var*(c2 + var*(c3 + var*c4)))

c      depvar = c0 + var*(c1 + var*(c2 + var*(c3 + var*(c4 + c5*var))))

      end 

      integer function match (idim,ier,name)
c----------------------------------------------------------------------
c find the index of endmember identified by name
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer idim, ier 

      character*8 mname, name

      common / cst18a /mname(m4)

      ier = 0

      do match = 1, idim
         if (name.eq.mname(match)) exit
      end do  

      if (match.gt.idim) ier = 1
      
      end 

      subroutine readn (idim,tname)
c----------------------------------------------------------------------
c readn - read idim endmember names expected format

c  name_1 name_2 ... | comment
c  name_n+1....      | comment

c end_of_data is either a "|" or the end of the record.
c----------------------------------------------------------------------
 
      implicit none

      include 'perplex_parameters.h'

      integer ibeg, iend, len, ier, idim, ict, i

      character chars(240)*1, mname*8, name*8, tname*10 

      common / cst18a /mname(m4)

      ier = 0 

      call readcd (n9,len,ier,chars)
      if (ier.ne.0) goto 90

      ibeg = 1
      iend = len
      ict = 0  

      do while (ict.lt.idim) 
c                                 find the name
         call readnm (ibeg,iend,len,ier,name,chars)
         if (ier.ne.0) goto 90
         ict = ict + 1
         mname(ict) = name

         if (ibeg.ge.len.and.ict.lt.idim) then
            call readcd (n9,len,ier,chars)
            ibeg = 1
            if (ier.ne.0) goto 90
         end if 

      end do

      goto 99      

90    write (*,1000) tname,(chars(i),i=1,len),name
      stop
      
1000  format ('**error ver200** READN bad data, currently ',
     *        'reading solution model: ',a,' data was:',/,240a,/,
     *        'last name read was: ',a,/)

99    end 


      subroutine readda (rnums,idim,tname)
c----------------------------------------------------------------------
c readda - reads idim numbers, discarding all data in records 
c          beyond the comment marker "|" or the end of the record.
c----------------------------------------------------------------------

      implicit none

      include 'perplex_parameters.h'

      integer len, idim, kdim, jdim, i, isnum, ier

      character chars(240)*1, tname*10, nums*240

      double precision rnums(100)
c----------------------------------------------------------------------
c                                 read card scans for non blank data
c                                 card:
      kdim = 1
      jdim = 0
      isnum = 0 

      do while (jdim.lt.idim)

         call readcd (n9,len,ier,chars)
         if (ier.ne.0) goto 91
c                                 got data, count how many "numbers"
         do i = 1, len
            if (chars(i).ne.' '.and.isnum.eq.0) then
               jdim = jdim + 1
               isnum = 1
            else if (chars(i).eq.' ') then
               isnum = 0
            end if 
         end do
c                                 if there is no comment
c                                 marker there can be more 
c                                 data than anticipated:
         if (jdim.gt.idim) jdim = idim
c                                 write numbers to string
         write (nums,*) (chars(i), i = 1, len),' '
         
         read (nums,*,err=90) (rnums(i), i = kdim, jdim)

         kdim = jdim + 1

      end do 

      goto 99 

90    write (*,1000) tname, (chars(i),i = 1, len)

      goto 99

91    write (*,1010) tname

      goto 99 

1000  format (' **error ver209** READDA bad data, currently',
     *        ' reading solution model: ',a,' data was:',/,240a,/,
     *        ' READDA was expecting numeric data but found',
     *        ' character data.',/)

1010  format (' **error ver210** READDA read to end of file',
     *        ' reading solution model: ',a,/,
     *        ' READDA was expecting numeric data.',/)

99    end 

      subroutine readx (idim,tname)
c----------------------------------------------------------------------
c readx - read excess function for a solution model, assumes
c data on one line of less than 240 characters, the expected format

c        W(name-name-...) number number number

c end_of_data is either a "|" or the end of the record.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ibeg, jend, len, ier, iscan, lord, imax, match, idim

      character chars(240)*1, name*8, begin*5, eod*3, tname*10

      integer i,j

      integer iend,isub,imd,insp,ist,isp,isite,iterm,iord,istot,jstot,
     *        kstot
      double precision wg,xmn,xmx,xnc
      common/ cst108 /wg(m1,m3),xmn(mst,msp),xmx(mst,msp),xnc(mst,msp),
     *      iend(m4),isub(m1,m2,2),imd(ms1,mst),insp(m4),ist(mst),
     *      isp(mst),isite,iterm,iord,istot,jstot,kstot
c----------------------------------------------------------------------

      iterm = 0 
      iord = 0 

      call readcd (n9,len,ier,chars)

      write (begin,'(5a)') (chars(i),i=1, 5)

      if (begin.eq.'ideal') then
         goto 99
      else if (begin.ne.'begin') then 
         goto 90
      end if 

      do i = 1, m1
         do j = 1, m2
            isub(i,j,1) = 0 
         end do 
      end do 
 
      eod = ' '

      do while (eod.ne.'end')  

         call readcd (n9,len,ier,chars)
         if (ier.ne.0) goto 90

         write (eod,'(3a)') (chars(i),i=1,3) 
c                                 find brackets
         ibeg = iscan (1,len,chars,'(') + 1
         imax = iscan (1,len,chars,')') - 1

         if (ibeg.gt.len.or.imax.gt.len) cycle
c                                 data found
         iterm = iterm + 1
         if (iterm.gt.m1) call error (48,wg(1,1),m1,tname)

         lord = 0

         do while (ibeg.lt.imax)

            jend = imax

            call readnm (ibeg,jend,imax,ier,name,chars)
            if (ier.ne.0) goto 90

            lord = lord + 1
            if (lord.gt.m2) call error (49,wg(1,1),m2,tname)

            isub(iterm,lord,1) = match (idim,ier,name)  

            if (ier.ne.0) goto 90         

         end do 

         if (lord.gt.iord) iord = lord

         ibeg = imax + 2
         jend = len

         do i = 1, 3

            call readfr (wg(iterm,i),ibeg,jend,len,ier,chars)
     
            if (ier.ne.0) goto 90 
    
         end do 
      end do

      goto 99

90    write (*,1000) tname,(chars(i),i=1,len),name
      stop
      
1000  format ('**error ver200** READX bad data, currently ',
     *        'reading solution model: ',a,' data was:',/,240a,/,
     *        'last name read was: ',a,/,
     *        'usually this error is due to a mispelled ',
     *        'endmember name.',/)

99    end 

      subroutine readop (idim,jlaar,kstot,maxit,tname)
c----------------------------------------------------------------------
c readop - tail of solution model to find optional dqf and
c          van laar size parameters

c readop - reads data until it finds an     "end_of_model" record

c          van laar data is identified by a "begin_van_la" record
c          dqf data is identified by a      "begin_dqf_co" record
c          or the iteration limit is        "max_iteratio" record

c readop returns:

c          jlaar = 1 if van laar data found (else 0).
c          idqf  > 0 if dqf data found (else 0)
c          maxit >= 0 if max_iteration is specified (else -1)
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer len, ier, idim, jlaar, maxit

      character chars(240)*1, begin*12, tname*10

      integer kstot,i,indq,idqf
      double precision dqf
      common/ cst222 /dqf(m3,m4),indq(m4),idqf
c----------------------------------------------------------------------

      jlaar = 0 
      idqf = 0 
      maxit = 99

      do 

         call readcd (n9,len,ier,chars)

         write (begin,'(12a)') (chars(i),i=1,12)

         if (begin.eq.     'end_of_model') then 
            
            exit

         else if (begin.eq.'begin_model ') then
c                              found new model, current
c                              model does not of end_of_model keyword
            write (*,1000) tname,(chars(i),i=1,len)
            stop

         else if (begin.eq.'begin_van_la') then 
c                              read van laar data: 
            jlaar = 1
            call readvl (idim,kstot,tname)
            cycle 

         else if (begin.eq.'begin_dqf_co') then 
c                              read dqf data:
            call readdq (idim,tname)
            cycle

         else if (begin.eq.'max_iteratio') then
c                              read the max # of iterations a solutions
c                              composition will be refined.
            call readmx (maxit,chars,tname)
            cycle 

         else 

            write (*,1010) tname,(chars(i),i=1,len)
            write (*,1020)
            stop

         end if

      end do  

1000  format (/,' **error ver200** READOP missing "end_of_model"',
     *          ' keyword at end',' of solution model:',a,/)
      
1010  format (/,' **error ver210** READOP bad data, currently',
     *          ' reading solution model: ',a,' data was:',/,240a)

1020  format (/,' This error is most probably due to an out-of-date',
     *          ' solution model file.',//,
     *          ' Copy the current version from:',//,
     *          ' www.perplex.ethz.ch/perplex/datafiles/',
     *          'solut_08.dat',//)

99    end 

      subroutine readmx (maxit,chars,tname)
c----------------------------------------------------------------------
c readvl - read van laar volumes for a solution models endmembers, assumes
c data on one line of less than 240 characters, the expected format

c        alpha(name) number number number

c end_of_data is either a "|" or the end of the record.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ibeg,len,iscan,i,iscnlt,maxit,ier

      character chars(240)*1, tname*10

      len = 240
c                                 locate first blank      
      ibeg = iscan (13,len,chars,' ') + 1
c                                 find the number (0-9)
      i = iscnlt (ibeg,len,chars,' ')
c                                 read the number 
      read (chars(i),'(i1)',iostat=ier) maxit 

      if (ier.ne.0) then 

        write (*,1000) tname,(chars(i),i=1,len)
        stop

      end if 

1000  format (' **error ver200** READMX bad data, currently',
     *        ' reading solution model: ',a,' data was:',/,240a)

      end 

      subroutine readvl (idim,kstot,tname)
c----------------------------------------------------------------------
c readvl - read van laar volumes for a solution models endmembers, assumes
c data on one line of less than 240 characters, the expected format

c        alpha(name) number number number

c end_of_data is either a "|" or the end of the record.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ibeg, len, ier, iscan, imax, match, idim, index

      character chars(240)*1, name*8, eod*3, tname*10

      integer kstot,jend,i,ict,jsmod
      double precision vlaar
      common/ cst221 /vlaar(m3,m4),jsmod
c----------------------------------------------------------------------

      ict = 0 
 
      eod = ' '

      do while (eod.ne.'end')  

         call readcd (n9,len,ier,chars)
         if (ier.ne.0) goto 90

         write (eod,'(3a)') (chars(i),i=1,3) 
c                                 find brackets
         ibeg = iscan (1,len,chars,'(') + 1
         imax = iscan (1,len,chars,')') - 1

         if (ibeg.gt.len.or.imax.gt.len) cycle
c                                 data found
         ict = ict + 1
         if (ict.gt.m4) goto 91
 
         jend = imax

         call readnm (ibeg,jend,imax,ier,name,chars)
         if (ier.ne.0) goto 90

         index = match (idim,ier,name)  

         if (ier.ne.0) goto 90         

         ibeg = imax + 2
         jend = len

         do i = 1, m3

            call readfr (vlaar(i,index),ibeg,jend,len,ier,chars)
     
            if (ier.ne.0) goto 90 
    
         end do 

      end do

      if (ict.lt.kstot) goto 91 

      goto 99

90    write (*,1000) tname,(chars(i),i=1,len),vlaar(i,index)
      write (*,1001)
      stop

91    write (*,1010) tname
      stop
      
1000  format ('**error ver200** READVL bad data, currently ',
     *        'reading solution model: ',a,' data was:',/,240a,/,
     *        'last number (or real equivalent) was: ',g12.6,/)

1001  format (/,'usually this error is caused by a mispelled ',
     *          'endmember name.',/)

1010  format (' **error ver201** READVL bad data, currently',
     *        ' reading solution model: ',a,/,
     *        ' this model requires 1 size parameter for',
     *        ' each independent endmember, READVL found ',i2,
     *        ' parameters.',/)

99    end 

      subroutine readdq (idim,tname)
c----------------------------------------------------------------------
c readvl - read dqf corrections for solution models endmembers, assumes
c data on one line of less than 240 characters, the expected format

c        dqf(name) number number number

c        output:

c          idqf           - the number of endmembers with dqf corrections
c          indq(idqf)     - pointer to corrected endmember in the phase name array
c          dqf(1..3,idqf) - the dqf parameters for the correction

c end_of_data is either a "|" or the end of the record.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ibeg, len, ier, iscan, imax, match, idim

      character chars(240)*1, name*8, eod*3, tname*10

      integer jend,i,idqf,indq
      double precision dqf
      common/ cst222 /dqf(m3,m4),indq(m4),idqf
c----------------------------------------------------------------------

      eod = ' '

      do while (eod.ne.'end')  

         call readcd (n9,len,ier,chars)
         if (ier.ne.0) goto 90

         write (eod,'(3a)') (chars(i),i=1,3) 
c                                 find brackets
         ibeg = iscan (1,len,chars,'(') + 1
         imax = iscan (1,len,chars,')') - 1

         if (ibeg.gt.len.or.imax.gt.len) cycle
c                                 data found
         idqf = idqf + 1
 
         jend = imax

         call readnm (ibeg,jend,imax,ier,name,chars)
         if (ier.ne.0) goto 90

         indq(idqf) = match (idim,ier,name)  

         if (ier.ne.0) goto 90         

         ibeg = imax + 2
         jend = len

         do i = 1, m3

            call readfr (dqf(i,idqf),ibeg,jend,len,ier,chars)
     
            if (ier.ne.0) goto 90 
    
         end do 

      end do

      goto 99

90    write (*,1000) tname,(chars(i),i=1,len),dqf(i,idqf)
      write (*,1001)
      stop
      
1000  format ('**error ver200** READDQ bad data, currently',
     *        'reading solution model: ',a,' data was:',/,240a,/,
     *        'last number (or real equivalent) was: ',g12.6,/)
1001  format (/,'usually this error is caused by a mispelled ',
     *          'endmember name.',/)

99    end 

      subroutine readr (coeffs,enth,inds,idim,nreact,tname)
c----------------------------------------------------------------------
c readr - read stoichiometric reaction data for a solution model, assumes
c data on one line of less than 240 characters, the expected format

c        name "=" (acoef(i), mame(i), i= 2..nreact) Enthalpy_value

c nreact >= 4 + nord for reciprocal reactions
c nreact = 3 for ordered species
c enthalpy_value is only read if on input nreact = 3

c end_of_data is either a "|" or the end of the record.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ibeg, iend, len, ier, iscan, i, nreact, inds(k7), 
     *        idim, match

      double precision coeffs(k7), enth, rnum

      character chars(240)*1, mname*8, name*8, tname*10 

      common / cst18a /mname(m4)
c----------------------------------------------------------------------
      ier = 0 

      call readcd (n9,len,ier,chars)
      if (ier.ne.0) goto 90

      ibeg = 1
      iend = len
c                                 first name
      call readnm (ibeg,iend,len,ier,name,chars)

      if (ier.ne.0) goto 90

      if (nreact.eq.3) then 
c                                 if nreact = 3, new name
         idim = idim + 1
         mname(idim) = name
         inds(1) = idim

      else 
c                                 reciprocal sol, get index
         inds(1) = match(idim,ier,name)
         if (ier.ne.0) goto 90
              
      end if 
c                                 find marker '='
      ibeg = iscan (1,len,chars,'=') + 1

      iend = len

      i = 2

      do 
c                                 find a stoichiometric coeff
         call readfr (rnum,ibeg,iend,len,ier,chars)

         coeffs(i) = rnum

         if (ier.ne.0) exit 
c                                 find the name
         call readnm (ibeg,iend,len,ier,name,chars)

         if (ier.ne.0) goto 90

         if (i.gt.k7) call error (1,0d0,i,'k7')

         inds(i) = match(idim,ier,name)

         if (ier.ne.0) goto 90

         if (nreact.gt.0.and.i.eq.nreact) exit 
c                                 increment counter and find 
c                                 the next coeff+name
         i = i + 1
           
      end do

      if (nreact.eq.3.and.i.eq.3) then      
c                                 ordered compound, read
c                                 enthalpy, find marker '='
         ibeg = iscan (ibeg,len,chars,'=') + 1

         call readfr (rnum,ibeg,iend,len,ier,chars)

         enth = rnum

         if (ier.ne.0) goto 90

      else if (i.lt.3) then 
c                                 a reaction with only 2 species 
c                                 is unexpected, write error
         goto 90 

      else 

         nreact = i - 1

      end if 

      goto 99      

90    write (*,1000) tname,(chars(i),i=1,len),name,rnum
      stop
      
1000  format ('**error ver200** READR bad data, currently ',
     *        'reading solution model: ',a,' data was:',/,240a,
     *        'last name read was: ',a,/,
     *        'last number (or real equivalent) was: ',g12.6,/)

99    end 

      subroutine readz (coeffs,inds,ict,idim,tname,tag)
c----------------------------------------------------------------------
c readz - read site fraction data for a solution model, assumes
c data on one line of less than 240 characters, the expected format

c        comments "=" a0 (acoef(i), mame(i), i= 1...end_of_data)

c end_of_data is either a "|" or the end of the record.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ibeg, iend, len, ier, iscan, ict, inds(k7),
     *        match, idim, i, iscnlt

      double precision rnum, coeffs(k7)

      character chars(240)*1, name*8, tname*10, tag*3
c----------------------------------------------------------------------
      ict = 0 

      call readcd (n9,len,ier,chars)
      if (ier.ne.0) goto 90 
c                                 this first segment is only done for
c                                 the readlm routine:
      ibeg = 1
c                                 read the fist word
      call readnm (ibeg,iend,len,ier,name,chars)
      read (name,'(a)') tag

      if (tag.eq.'end') then 
c                                 if called from readlm tag may be the
c                                 end of data marker
         goto 99

      else 

         i = match(idim,ier,name)

         if (ier.eq.0) then 
c                                 the "comment" is a name for
c                                 readlm
            ict = ict + 1
            inds(ict) = i

         end if
 
      end if 
c                                 find start of data marker '='
      ibeg = iscan (iend,len,chars,'=') + 1
      ict = ibeg
c                                 find a number
      call readfr (rnum,ibeg,iend,len,ier,chars)
      if (ier.ne.0) goto 90 
c                                 find the next non-blank chracter
c                                 if it's text, then the expression
c                                 has no constant, else save the 
c                                 constant.
      if (chars(iscnlt(iend+1,240,chars,'/')).lt.'A') then
c                                 assuming ascii collating sequence,
c                                 the next character is numeric, so 
c                                 the number read previously is the a0 constant
         coeffs(1) = rnum

      else 
c                                 no constant, reset character pointer.
         coeffs(1) = 0d0 
         ibeg = ict

      end if 
c                                 the rest of the data should
c                                 consist of coefficients followed
c                                 by names
      ict = 1 

      do while (ibeg.lt.len) 
c                                 find the number
         call readfr (rnum,ibeg,iend,len,ier,chars)

         if (ier.ne.0) then 
c                                 if called from readlm may be the 
c                                 legitimate text string "delta"
            call readnm (ibeg,iend,len,ier,name,chars)

            if (name.eq.'delta') then 

               ibeg = iscan (iend,len,chars,'=') + 1
               call readfr (rnum,ibeg,iend,len,ier,chars)
               if (ier.ne.0) goto 90
               coeffs(ict+1) = rnum
               exit 

            else 
c                                 invalid data 
               goto 90

            end if 

         end if

         if (ier.ne.0) goto 90
c                                 find the name
         call readnm (ibeg,iend,len,ier,name,chars)
         if (ier.ne.0) goto 90

         ict = ict + 1
         coeffs(ict) = rnum 
         inds(ict) = match(idim,ier,name)

         if (ier.ne.0) then 
            write (*,1010) name,tname,(chars(i),i=1,len)
            stop
         end if 
           
      end do

      goto 99

90    write (*,1000) tname,(chars(i),i=1,len),name,rnum
      stop
      
1010  format (/,'**error ver201** invalid name: ',a,' in an expression',
     *        ' for solution model: ',a,/,' data was:',/,240a)
1000  format (/,'**error ver200** READZ bad data, currently',
     *        ' reading solution model: ',a,' data was:',/,240a,/,
     *        'last name read was: ',a,/,
     *        'last number (or real equivalent) was: ',g12.6,/)

99    end 
 
      subroutine xrmode (tname)
c---------------------------------------------------------------------
c xrmode - reads old style solution models from LUN n9.
c---------------------------------------------------------------------
      implicit none
  
      include 'perplex_parameters.h'
 
      character*10 mname*8,tname

      integer i,j,k,l,m,inew,ijk(mst),h,jmsol,kdsol

      common/ cst18a /mname(m4)
     *      / cst142 /jmsol(m4,mst),kdsol(m4)/ cst140 /inew

      integer nsub1,nsub2,nttyp,nterm,nspm1,nsite
      double precision acoef,smult,a0
      common/ cst107 /a0(m10,m11),acoef(m10,m11,m0),smult(m10),
     *      nsite,nspm1(m10),nterm(m10,m11),nsub1(m10,m11,m0,m12),
     *      nsub2(m10,m11,m0,m12),nttyp(m10,m11,m0)

      integer iend,isub,imd,insp,ist,isp,isite,iterm,iord,istot,jstot,
     *        kstot
      double precision wg,xmn,xmx,xnc
      common/ cst108 /wg(m1,m3),xmn(mst,msp),xmx(mst,msp),xnc(mst,msp),
     *      iend(m4),isub(m1,m2,2),imd(ms1,mst),insp(m4),ist(mst),
     *      isp(mst),isite,iterm,iord,istot,jstot,kstot

      integer mdep,idep,jdep,ndph
      double precision nu,y2p
      common/ cst146 /nu(m15,j4),y2p(m4,m15),mdep,jdep(m15),
     *                idep(m15,j4),ndph(m15)
c----------------------------------------------------------------------

      mdep = 0 
      istot = 0
c                               read number of independent sites:
      read (n9,*,err=90) isite

      do i = 1, mst
         isp(i) = 1
         ist(i) = 1
      end do 
c                               read number of species on each site:
      read (n9,*,err=90) (isp(i),ist(i),i=1,isite)
c                               total number of endmembers:
      istot = 1
      do i = 1, isite
         istot = istot*isp(i)
      end do 
c                               read endmember names:
      read (n9,1000) (mname(i), i = 1, istot)
c                               read endmember flags:
      read (n9,*,err=90) (iend(i), i= 1, istot)
c                               read composition limits:
      do i = 1, isite

         if (inew.eq.0) then 
c                               old style file
            read (n9,*,err=90) (xmn(i,j), xmx(i,j), xnc(i,j), 
     *                          j = 1, isp(i) - 1), imd(1,i)
            do j = 2, isp(i) - 1
               imd(j,i) = imd(1,i)
            end do 
         else
            read (n9,*,err=90) (xmn(i,j), xmx(i,j), xnc(i,j), 
     *                          imd(j,i), j = 1, isp(i) - 1)
         end if 
      end do 
c                                create bragg-williams indexes
      do i = 2, isite
         ijk(i) = 1
      end do 
 
      ijk(1) = 0  

      do 20 l = 1, istot

         do m = 1, isite

            if (ijk(m).lt.isp(m)) then 

               ijk(m) = ijk(m) + 1
c                                increment only one index per endmember
               do i = 1, isite
                  jmsol(l,i) = ijk(i)
               end do 

               if (ijk(m).eq.isp(m).and.l.lt.istot-isp(1)+1) then
c                                saturated counter, increment first
c                                unsaturated counter and reset all
c                                lower site counters to first index,
                  do j = 2, isite
                     if (ijk(j).lt.isp(j)) then
                        ijk(j) = ijk(j) + 1
                        do k = 2, j-1
                           ijk(k) = 1
                        end do 
                        ijk(1) = 0
                        goto 20 
                     end if
                  end do  
               end if 

               goto 20
            end if
         end do 
20    continue 
c                               read total number of terms, and order
c                               of highest order term, in mixing model
c                               ideal models: iterm = iord = 0.
      read (n9,*,err=90) iterm, iord
      if (iterm.ne.0) then 
         if (iterm.gt.m1) call error (48,wg(1,1),m1,tname)
         if (iord.gt.m2) call error (49,wg(1,1),m2,tname)
c                               read parameters for margules
c                               expression of the excess free
c                               energy of the solution.
         do h = 1, iterm
c                               weigthing by bragg-williams fractions
            read (n9,*,err=90) ((isub(h,i,j), j = 1, 2),i = 1, iord)
            read (n9,*,err=90) (wg(h,k), k = 1, m3)
         end do 
      end if 
c                              expansion for S(configurational)
      read (n9,*,err=90) nsite

      if (nsite.ne.0) call smodel (tname) 
c                              read p-t dependent size parameter for 
c                              each endmember
c     if (inew.ne.0) read (n9,*,err=90) idummy

      goto 99

90    call error (20,wg(1,1),i,tname)

1000  format (3(a,1x))

99    end 

      subroutine smodel (tname)
c-----------------------------------------------------------------
c smodel reads models for configurational entropy
c of solutions with dependent site mixing or disordered
c endmembers, see documentation section 1.3.1, equation (8b)
c-----------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      character*10 tname

      integer i,j,knsp,k,l

      integer jsmod
      double precision vlaar
      common/ cst221 /vlaar(m3,m4),jsmod

      integer nsub1,nsub2,nttyp,nterm,nspm1,nsite
      double precision acoef,smult,a0
      common/ cst107 /a0(m10,m11),acoef(m10,m11,m0),smult(m10),
     *      nsite,nspm1(m10),nterm(m10,m11),nsub1(m10,m11,m0,m12),
     *      nsub2(m10,m11,m0,m12),nttyp(m10,m11,m0)

      integer iend,isub,imd,insp,ist,isp,isite,iterm,iord,istot,jstot,
     *        kstot
      double precision wg,xmn,xmx,xnc
      common/ cst108 /wg(m1,m3),xmn(mst,msp),xmx(mst,msp),xnc(mst,msp),
     *      iend(m4),isub(m1,m2,2),imd(ms1,mst),insp(m4),ist(mst),
     *      isp(mst),isite,iterm,iord,istot,jstot,kstot
c-----------------------------------------------------------------

      if (nsite.gt.m10) call error (31,a0(1,1),isite,tname)
c
      if (nsite.lt.isite) call error (30,a0(1,1),isite,tname)
c                                 for each site
      do i = 1, nsite 
c                                 read # of species, and site 
c                                 multiplicty.
         read (n9,*,err=90) knsp,smult(i)
         nspm1(i) = knsp - 1
c                                 for each species, read
c                                 function to define the
c                                 site fraction of the species:
         do j = 1, nspm1(i)
c                                 read # of terms in the 
c                                 site fraction function and a0.
            read (n9,*,err=90) nterm(i,j), a0(i,j)
            if (nterm(i,j).gt.m0) call error (33,a0(1,1),m0,tname)
c                                 for each term:
            do k = 1, nterm(i,j)
c                                 read term type:
               read (n9,*,err=90) nttyp(i,j,k)
               if (nttyp(i,j,k).lt.m12) then
                  if (jsmod.eq.1) then 
                     read (n9,*,err=90) acoef(i,j,k),(nsub1(i,j,k,l),
     *                            nsub2(i,j,k,l),l = 1, nttyp(i,j,k))
                  else
                     read (n9,*,err=90) acoef(i,j,k),(nsub1(i,j,k,l),
     *                                  l = 1, nttyp(i,j,k))                     
                  end if 
               else
                  call error (29,a0(1,1),nttyp(i,j,k),tname)
               end if
            end do 

         end do 
      end do 

      return

90    call error (20,a0(1,1),i,tname)

      end 

      double precision function gsixtr (id)
c-----------------------------------------------------------------------
c gsixtr computes G from the EoS formulated by Sixtrude & Bukowski '90, 
c JGR 95:19311-19325.

c Sixtrude parameters:

c    F0    n  -v(J/bar) K0(bar) K0'   theta0 gamma0 q etaS0 Smag

c Perple_X parameters (1st index in array thermo):

c    1     2     3        4     5      6      7    8    9    10

c and for shear modulii

c    G0(bar)   G0'

c are in the Perple_X array emod(i,id) with i

c      1        2
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer id, itic, izap, izap1

      double precision turd2, a0, v0, v, df, f, dv, root, nr9t0, nr9t,
     *                 gamma0, k00, plg, c1, c2, c3, f1, gvq, 
     *                 n, q, vq, v23, theta0, tol, k0p, a, ethv

      double precision nr9, qm1, d2f, tht, tht0, etht, etht0, ltht, df1,
     *                 ltht0, dtht, dtht0, d2tht, d2tht0, g, g0, dg, 
     *                 dg0, d2g, d2g0, dfc, d2fc, dft, d2ft, dft0, d2ft0
 
      double precision thermo, uf, us
      common/ cst1 /thermo(k4,k10),uf(2),us(h5)

      double precision smu
      common/ cst323 /smu

      character*8 names 
      common/ cst8   /names(k1)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      integer iemod,kmod
      logical smod
      double precision emod
      common/ cst319 /emod(k15,k10),smod(h9),iemod(k10),kmod

      save turd2, izap, izap1 

      data turd2, izap, izap1 /0.6666666666666666667d0, 0, 0/
c----------------------------------------------------------------------
c                                 assign local variables:
      a0     =  thermo(1,id)
      n      =  thermo(2,id)
      v0     = -thermo(3,id)
      k00    =  thermo(4,id)
      k0p    =  thermo(5,id)
      theta0 =  thermo(6,id)
      gamma0 =  thermo(7,id)
      q      =  thermo(8,id)

      nr9    = thermo(11,id)
      c1     = thermo(12,id)
      c2     = thermo(13,id)
      c3     = thermo(14,id)
      qm1    = q - 1d0
      nr9t   = nr9*t
      nr9t0  = thermo(20,id)
c                                 initial guess for volume:
c                                 taylor(diff(FTH,v),v=v0,1)
c                                 JADC Feb 26, 2008. 
c                                 the dfth0 could be loaded as a 
c                                 constant. 
      tht    = theta0/t
      tht0   = theta0/tr
      k00    = thermo(4,id)
      k0p    = thermo(5,id)

      dft   = nr9t*gamma0/v0*(3d0*plg(tht)/tht**3 
     *        - dlog(1d0-exp(-tht)))
      dft0  = nr9t0*gamma0/v0*(3d0*plg(tht0)/tht0**3 
     *        - dlog(1d0-exp(-tht0)))
c                                 taylor(diff(FC,v),v=v0,2)
c     v       = (k00-dft+dft0-p)/k00*v0
c                                 taylor(diff(FC,v),v=v0,3)
      root = k00*((2d0+2d0*k0p)*(p+dft-dft0) + k00)

      if (root.gt.0d0) then 
         v = ((2d0+k0p)-dsqrt(root)/k00)*v0/(1d0+k0p)
         if (v.lt.v0/2d0.or.v.gt.v0*2d0) v = v0
      else 
         v = v0
      end if 

      f1 = 1d9
      itic = 0 
c                                 change to use relative tolerance
c                                 JADC March 1, 2005. formerly 1d-1 bar.
      tol = 1d-6*p

      do while (dabs(f1).gt.tol)

         itic = itic + 1

         vq = (v/v0)**q
         gvq = gamma0*vq
         v23 = (v0/v)**turd2
         f = 0.5d0*v23 - 0.5d0
         df = -v23/v/3d0
         d2f = 5d0/9d0*v23/v**2

         tht  =  theta0*dexp(-gamma0*((v/v0)**q-1d0)/q)/t
         if (tht.lt.1d-10) goto 90 
         tht0 =  tht*t/tr
         etht  = dexp(-tht )
         etht0 = dexp(-tht0)
         ltht  = dlog(1d0 - etht )
         ltht0 = dlog(1d0 - etht0)
c                                 diff(theta/T,v)
         dtht  = -gvq/v*tht
         dtht0 = -gvq/v*tht0
c                                 diff(theta/T,v,v)
         d2tht  = gvq*tht /v**2*(gvq - qm1)
         d2tht0 = gvq*tht0/v**2*(gvq - qm1)

         g   = plg(tht )
         g0  = plg(tht0)
         dg  = tht **2*ltht *dtht
         dg0 = tht0**2*ltht0*dtht0
         d2g  = ((2d0*ltht  + tht *etht /(1d0-etht ))*dtht **2 + 
     *          tht *ltht *d2tht )*tht
         d2g0 = ((2d0*ltht0 + tht0*etht0/(1d0-etht0))*dtht0**2 + 
     *          tht0*ltht0*d2tht0)*tht0

         dfc = (c3*f+c1)*f*df
         d2fc = (2d0*c3*f+c1)*df**2+(c3*f+c1)*f*d2f

         dft  = nr9t /tht **3*(dg  -3d0/tht *g *dtht )
         dft0 = nr9t0/tht0**3*(dg0 -3d0/tht0*g0*dtht0)

         d2ft =  nr9t /tht **3*(3d0/tht *(dtht *(4d0/tht *g *dtht  
     *                        - 2d0*dg ) - g *d2tht ) + d2g )
         d2ft0 = nr9t0/tht0**3*(3d0/tht0*(dtht0*(4d0/tht0*g0*dtht0 
     *                        - 2d0*dg0) - g0*d2tht0) + d2g0)

         f1  = -dfc - dft + dft0 - p

         df1 = -d2fc - d2ft + d2ft0 

         dv = f1/df1
         v = v - dv

         if (v.le.0d0.or.v/v0.gt.2d1.or.
     *       itic.gt.100.or.dabs(f1).gt.1d40) goto 90
 
      end do
c                                 if everything is ok, now get 
c                                 helmoltz energy:
      goto 10 
c                                 if we get here, failed to converge
90    if (izap.lt.10) then
         write (*,1000) t,p,names(id)
      else if (izap.eq.100) then 
         write (*,1010)
      end if 

      izap = izap + 1
c                                 destabilize the phase:
      gsixtr = 1d10
      goto 99 
         
10    vq = (v/v0)**q
      f = 0.5d0*(v0/v)**turd2 - 0.5d0
      tht  =  theta0*dexp(-gamma0*(vq-1d0)/q)/t
      tht0 =  tht*t/tr 

      a = a0 + c1*f**2*(0.5d0 + c2*f) 
     *  + nr9*(t/tht**3*plg(tht ) -tr/tht0**3*plg(tht0))

      gsixtr = a + p*v - t*thermo(10,id)
c                                 thermal energy/v
      ethv = (dft0-dft)/gamma0/vq
c                                 etas0 = thermo(9,id)
c                                 g0 = emod(1,id)
c                                 g0p = emod(2,id)
c                                 etas = thermo(9,id)*v/v0
c                                 adiabatic shear modulus
      smu = (1d0 + 2d0*f)**(2.5d0)*
     *      (emod(1,id)*(1d0 - 5d0*f) + f*emod(2,id)*3d0*k00)
     *    -  thermo(9,id)*v/v0*ethv

      if (tht.gt.5d0.or.tht0.gt.5d0.and.izap1.lt.100) then
         write (*,*) 'danger will robinson, danger, danger!!',
     *               id,tht0,tht
         izap1 = izap1 + 1
      end if 
                   
1000  format (/,'**warning ver369** failed to converge at T= ',f8.2,' K'
     *       ,' P=',f9.1,' bar',/,'Using Sixtrude EoS.',
     *        ' Phase ',a,' will be destabilized.',/)
1010  format (/,'**warning ver369** will not be repeated for future ',
     *          'instances of this problem',/)

99    end 

      double precision function gstxgi (id)
c-----------------------------------------------------------------------
c gstxgi computes G from the EoS formulated by Sixtrude & Lithgow-B

c Stxrude parameters:

c    F0    -n  -v(J/bar) K0(bar) K0'   theta0 gamma0 q etaS0 Smag

c Perple_X parameters (1st index in array thermo):

c    1     2     3        4     5      6      7    8    9    10

c and for shear modulii

c    G0(bar)   G0'

c are in the Perple_X array emod(i,id) with i

c      1        2
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer id, itic, izap, izap1

      logical bad

      double precision v0, v, df, f, dv, gamma0, k00, k0p,
     *           plg, c1, c2, c3, f1, aiikk, aiikk2, nr9t,
     *           root, aii, etas, a, ethv, gamma, da, nr9t0,
     *           fpoly, fpoly0, letht, letht0, z, aii2, f23,
     *           v23, tol, t1, t2, a2f, f59

      double precision nr9, d2f, tht, tht0, etht, etht0, df1,
     *                 dtht, dtht0, d2tht, d2tht0, 
     *                 dfc, d2fc, dfth, d2fth, dfth0, d2fth0
 
      double precision thermo, uf, us
      common/ cst1 /thermo(k4,k10),uf(2),us(h5)

      double precision smu
      common/ cst323 /smu

      character*8 names 
      common/ cst8   /names(k1)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      integer iemod,kmod
      logical smod
      double precision emod
      common/ cst319 /emod(k15,k10),smod(h9),iemod(k10),kmod

      save f23,f59, izap, izap1 

      data f23, f59, izap, izap1 /0.66666666666666666667d0,
     *                            0.55555555555555555556d0,0, 0/
c----------------------------------------------------------------------
c                                 assign local variables:
      v0     = -thermo(3,id)
      nr9    = thermo(11,id)
      c1     = thermo(12,id)
      c2     = thermo(13,id)
      c3     = thermo(14,id)
      aii    = thermo(15,id)
      aiikk  = thermo(16,id)
      aiikk2 = thermo(18,id)
      aii2   = thermo(19,id)
      nr9t0  = thermo(20,id) 

      t1     = thermo(6,id)/t
      t2     = t/tr
      nr9t   = nr9*t
c                                 initial guess for volume:
c                                 taylor(diff(FTH,v),v=v0,1)
c                                 JADC Feb 26, 2008. 
c                                 the dfth0 could be loaded as a 
c                                 constant. 
      tht    = t1
      tht0   = tht*t2
      gamma0 = thermo(7,id)
      k00    = thermo(4,id)
      k0p    = thermo(5,id)

      dfth   = nr9t*gamma0/v0*(3d0*plg(tht)/tht**3 
     *         - dlog(1d0-exp(-tht)))
      dfth0  = nr9t0*gamma0/v0*(3d0*plg(tht0)/tht0**3 
     *         - dlog(1d0-exp(-tht0)))
c                                 taylor(diff(FC,v),v=v0,2)
c     v       = (k00-dfth+dfth0-p)/k00*v0
c                                 taylor(diff(FC,v),v=v0,3)
      root = k00*((2d0+2d0*k0p)*(p+dfth-dfth0) + k00)

      if (root.gt.0d0) then 
         v = ((2d0+k0p)-dsqrt(root)/k00)*v0/(1d0+k0p)
         if (v.lt.v0/2d0.or.v.gt.v0*2d0) v = v0
      else 
         v = v0
      end if 

      itic = 0 
c                                 change to use relative tolerance
c                                 JADC March 1, 2005. formerly 1d-1 bar.
      tol = 1d-6*p

      do 

         itic = itic + 1
c                                 f, and derivatives
         v23 = (v0/v)**f23
         f = 0.5d0*v23 - 0.5d0
         df = -v23/v/3d0
         d2f = f59*v23/v**2
c                                 cold part derivatives
         dfc = (c3*f+c1)*f*df
         d2fc = (2d0*c3*f+c1)*df**2+(c3*f+c1)*f*d2f
c                                 debye T/T (tht)
         z  = 1d0+(aii+aiikk2*f)*f

         if (z.lt.0d0) then
            bad = .true.
            exit 
         end if 

         root = dsqrt(z)

         tht   = t1*root
         tht0  =  tht*t/tr
c                                 tht derivatives
         a2f   = aii2+aiikk2*f
         da    = a2f/root
         dtht  = t1*da*df
         d2tht = t1*((aiikk2/root-a2f**2/z**1.5d0)*df**2
     *               + da*d2f) 

         dtht0 = dtht*t2
         d2tht0 = d2tht*t2
c                                 polylog functions:
         fpoly   = 3d0*plg(tht )/tht**3
         fpoly0  = 3d0*plg(tht0)/tht0**3
c                                 thermal part derivatives:
         etht  = dexp(-tht )

         if (1d0-etht.lt.0d0) then
            bad = .true.
            exit 
         end if 

         letht = dlog(1d0-etht)

         dfth = (letht-fpoly)*nr9t*dtht/tht
         d2fth = ((4d0*dtht**2/tht-d2tht)*(fpoly-letht)
     *         + dtht**2*etht/(1d0-etht))*nr9t/tht

         etht0 = dexp(-tht0)

         if (1d0-etht0.lt.0d0) then
            bad = .true.
            exit 
         end if 
 
         letht0 = dlog(1d0-etht0)

         dfth0 = (letht0-fpoly0)*nr9t0*dtht0/tht0
         d2fth0 = ((4d0*dtht0**2/tht0-d2tht0)*(fpoly0-letht0)
     *          + dtht0**2*etht0/(1d0-etht0))*nr9t0/tht0

         f1  = -dfc - dfth + dfth0 - p

         df1 = -d2fc - d2fth + d2fth0 

         dv = f1/df1

         if (dabs(dv).gt.1d-2) dv = 1d-2*dv/dabs(dv)

         v = v - dv
c                                 the limit on f1 should be 
c                                 machine dependent
         if (v.le.0d0.or.itic.gt.100.or.dabs(f1).gt.1d40) then 
            bad = .true.
            exit 
         else if (dabs(f1).lt.tol) then
            bad = .false.
            exit 
         end if 
 
      end do

      if (bad) then  
c                                 if we get here, failed to converge
        if (izap.lt.10) then
            write (*,1000) t,p,names(id)
            if (itic.gt.100) write (*,*) 'failed cause of it limit'
         else if (izap.eq.100) then 
            write (*,1010)
         end if 

         izap = izap + 1
c                                 destabilize the phase.
         gstxgi  = 1d10

      else 

c                                 everything is ok, now get 
c                                 helmoltz energy:
         f = 0.5d0*(v0/v)**f23 - 0.5d0 
         z = 1d0+(aii+aiikk2*f)*f
         root = dsqrt(z)
c                                 final estimate for tht
         tht   = t1*root
         tht0  = tht*t2
c                                 helmholtz enery
         a = thermo(1,id) + c1*f**2*(0.5d0 + c2*f) 
     *     + nr9*(t/tht**3*plg(tht ) -tr/tht0**3*plg(tht0))

         gstxgi = a + p*v - t*thermo(10,id)
c                                 z = (theta/theta0)^2
         gamma = (2d0*f+1d0)*(aii+aiikk*f)/6d0/z
         etas = - gamma - thermo(17,id)/z*(2d0*f + 1d0)**2
c                                 thermal energy/V, based on
c                                 previous v estimate
         ethv = (dfth0-dfth)/gamma
c                                 adiabatic shear modulus
         smu = (1d0 + 2d0*f)**(2.5d0)*(
     *         emod(1,id) + f*(thermo(21,id) + thermo(22,id)*f))
     *       - etas*ethv    
     
         if (tht.gt.5d0.or.tht0.gt.5d0.and.izap1.lt.100) then
            write (*,*) 'danger will robinson, danger, danger!!',
     *                  id,tht0,tht
            izap1 = izap1 + 1
         end if
      end if  
                   
1000  format (/,'**warning ver369** failed to converge at T= ',f8.2,' K'
     *       ,' P=',f9.1,' bar',/,'Using Sixtrude EoS.',
     *        ' Phase ',a,' will be destabilized.',/)
1010  format (/,'**warning ver369** will not be repeated for future ',
     *          'instances of this problem',/)
99    end 

      double precision function plg (t)
c-----------------------------------------------------------------------
c 29th order series expansion of polylog terms in sixtrude's EoS about
c t = 0, good to 0.0001 accuracy to t = 5.
c f := int((ln(1-exp(-t))*t^2),t=0..TD);
c gg := convert(series(f,TD=0,29),polynom);
c-----------------------------------------------------------------------
      implicit none

      double precision c2, t

      save c2 

      data c2 /0.1111111111111111111d0/

c                              13th order expansion good to t ~ 3
c       c1 = 1/3
c     plg  = t**3*((dlog(t)*c1-c2) - t/8d0 + t**2/120d0     
c    *     - t**4/20160d0 + t**6/1632960d0 - t**8/106444800d0)

      plg  = (t**3*((dlog(t)/3d0-c2) - t/8d0 + t**2/120d0     
     *     - t**4/20160d0 + t**6/1632960d0 - t**8/106444800d0
     *     + t**10/6227020800d0 - 0.2935661189d-11*t**12 
     *     + 0.562291451d-13*t**14 - 0.1115026413d-14*t**16
     *     + 0.2271444989d-16*t**18 - 0.4727975432d-18*t**20
     *     + 0.1001636878d-19*t**22 - 0.2153466772d-21*t**24))

      end 

      double precision function gmurn (id)
c-----------------------------------------------------------------------
c gmurn computes the vdp integral of a compound identified by id
c that is described by ghiorso's pmelt murnaghan EoS.

c Ghiorso parameters:

c v(J/bar) dv/dt(J/K/bar) dv/dp(J/bar^2) dv/dp/dt(J/bar^2/K) -K'

c Perple_X parameters (1st index in array thermo):

c    3            11           12             13              18

c N.B. for some reason Ghiorso chooses a reference T = 1673 for tbe
c mechanical EoS parameters, this is hardwired here as variable trv.

c modified may 11, 2006 to use holland and powell's alpha and kt expresions
c if thermo(id,16) ne 0.

c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer id, itic, jerk

      double precision thermo, uf, us, turd2, trv,
     *                 v0, k, turd, rat, rat2, c0, c1, 
     *                 c2, c3, c4, c5, a0, a1, v, df, f, dv, kprime
 
      common/ cst1 /thermo(k4,k10),uf(2),us(h5)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      save turd, turd2, trv, jerk 

      data turd, turd2, trv, jerk /0.3333333333333333333d0,
     *                      0.6666666666666666667d0, 1673.15d0, 0/
c----------------------------------------------------------------------
      if (thermo(16,id).eq.0) then 
c                                 assume ghiorso's expressions for v0 and k (KT)
c                                 v0 = Volume at 1 bar, T
         v0 = thermo(3,id) + thermo(11,id) * (t - trv)

         k = -v0 / (thermo(12,id) + thermo(13,id) * (t - trv) )
      
      else
c                                 use hopo's expressions
         k = thermo(16,id) + thermo(17,id)*t
c                                 a ****wit has entered a ridiculous
c                                 temperature
         if (k.lt.0d0) call error (46,t,id,'GMURN') 

         v0 = (thermo(3,id) + (thermo(11,id) + thermo(12,id)*t)*t +
     *          thermo(13,id)*dlog(t) + thermo(14,id)/t + 
     *          thermo(15,id)*dsqrt(t)) 

      end if 
c                                 NOTE kprime = -K'
      kprime = thermo(18,id)
c                                 constants:
      a0 = 0.375d0 * v0 * k
      a1 = -0.125d0 * v0**2 * k
      c0 = (-28d0 -6d0 * kprime) * v0 * a0
      c1 = (12d0 + 3d0 * kprime) * v0**2 * a0
      c2 = (16d0 + 3d0 * kprime) * a0
      c3 = a1 * v0 * (-196d0 - 42d0 * kprime)
      c4 = a1 * (80d0 + 15d0 * kprime)
      c5 = a1 * v0 * (108d0 + 27d0 * kprime)
c                                 initial guess for volume:
      v = k * v0 / (p + k)
      dv = 1d0
      itic = 0 

      do while (dabs(dv).gt.1d-5)

         itic = itic + 1
         rat = (v0/v)**turd
         rat2 = rat**2
         f = p  + ((c0*v*rat+c1+c2*v**2*rat2)/v**3)
         df = (c3/rat2+c4*v/rat+c5)/v**4
         dv = f/df
         v = v - dv
         if (v.le.0d0.or.v.gt.1d6.or.itic.gt.20) then

            if (jerk.lt.10) then 
               jerk = jerk + 1
               write (*,1000) t,p
               if (jerk.eq.50) write (*,1010)
            end if 

            if (k*(k-2d0*p*(1d0+kprime)).gt.0d0) then 

               v = 0.5d0/((3d0+kprime)*k+2d0*p)*((2d0+kprime)*2d0*k
     *             +2d0*dsqrt(k*(k-2d0*p*(1d0+kprime))))*v0 

               if (v.lt.0d0) v = 0.5d0/((3d0+kprime)*k+2d0*p)*
     *                       ((2d0+kprime)*2d0*k
     *                       -2d0*dsqrt(k*(k-2d0*p*(1d0+kprime))))*v0 
            else 

               v = k * v0 / (p + k)

            end if 

            exit

         end if 

1000  format (/,'**warning ver369** failed to converge at T= ',f8.2,' K'
     *       ,' P=',f9.1,' bar',/,'Using Murnaghan ',
     *        'EoS, probably for Ghiorso et al. MELTS/PMELTS endmember',
     *        ' data.',/,
     *        'Volume estimated using 3rd order taylor series.',/)
1010  format (/,'**warning ver369** will not be repeated for future ',
     *          'instances of this problem',/)

      end do
c                                 and the vdp integral is:
      f = 0.5d0*((v0/v)**turd2-1d0)

      gmurn = p*v - v0*(pr-4.5d0*k*f**2*(1d0-f*(4d0+kprime))) 

      end 

      subroutine cartes (ksite,sname,ids)
c---------------------------------------------------------------------
c subroutine to cartesian or transform subdivision on site ksite of 
c solution ids (or sname). called by subdiv.
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer mres
 
      parameter (mres=3000)

      integer mode, i, ind(ms1), iy(ms1), jsp, ksite, indx, iexit, 
     *        ieyit, j, ids

      double precision tol, y(ms1,mres), ycum, ymax, dy, ync, res, ylmn,
     *                 ylmx, yloc, x, unstch, strtch, yreal

      logical odd

      character sname*10

      integer ntot,npairs
      double precision yy,xy
      common/ cst86 /xy(mdim,k21),yy(k21,mst,ms1),ntot,npairs

      integer iend,isub,imd,insp,ist,isp,isite,iterm,iord,istot,jstot,
     *        kstot
      double precision wg,xmn,xmx,xnc
      common/ cst108 /wg(m1,m3),xmn(mst,msp),xmx(mst,msp),xnc(mst,msp),
     *      iend(m4),isub(m1,m2,2),imd(ms1,mst),insp(m4),ist(mst),
     *      isp(mst),isite,iterm,iord,istot,jstot,kstot
c                                 interval limits conformal transformation
      integer intv
      double precision yint, yfrc
      common/ cst47 /yint(5,ms1,mst,h9),yfrc(4,ms1,mst,h9),intv(4)
c                                 x coordinate description
      integer istg, ispg, imlt, imdg
      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)

      character fname*10
      common/ csta7 /fname(h9)

      logical good(h9)

      save tol, good

      data tol, good/1d-5,h9*.true./
c----------------------------------------------------------------------
      ycum = 0d0
      jsp = isp(ksite) - 1

      if (jsp.eq.0) then 
c                                 a relict site with only one species
c                                 left over from reform, i have know idea
c                                 if this works
         y(1,1) = xmn(ksite,1)
         npairs = 1
         goto 99 
      end if 

      do i = 1, jsp
c                                 generate coordinates for i'th component
         iy(i) = 1
         y(i,1) = xmn(ksite,i)
         ync = xnc(ksite,i)

         if (ync.eq.0d0) cycle

         mode = imdg(i,ksite,ids)
c                                 avoid impossible compositions 'cause a min > 0
         if (i.gt.1) then 

            ycum = ycum + xmn(ksite,i-1)
c                                 1-ycum is the smallest fraction possible
            if (1d0-ycum.lt.0) then 
c                                 inconsistent limits
               call error (999,ycum,jsp,'cartes')

            else
c                                 the smallest fraction possible is lt
c                                 than xmax
               ymax = xmx(ksite,i)

            end if 
         else 
            ymax = xmx(ksite,i)
         end if 
c                                 two means of extracting y-range, cartesian
c                                 imod = 0 and transformation imod = 1
         if (mode.eq.0) then 
c                                 check limit
            if (ync.lt.tol.and.good(ids)) then
               good(ids) = .false.
               write (*,1000) fname(ids)
            end if
               
c                                 cartesian
            do while (y(i,iy(i)).lt.ymax)
               iy(i) = iy(i) + 1
               if (iy(i).gt.mres) call error (999,ync,mres,sname)
               y(i,iy(i)) = y(i,iy(i)-1) + ync
               if (dabs(y(i,iy(i))-ymax).lt.tol) then
                  y(i,iy(i)) = ymax
               else if (y(i,iy(i)).gt.ymax) then
                  y(i,iy(i)) = ymax
               end if 
            end do

         else 
c                                 conformal x is the cartesian coordinate
c                                 y is the real coordinate.
            if (mode.lt.4) then 
               odd = .false.
            else
               odd = .true.
            end if 

            res = 0d0
c                                 there are as many as intv(mode)
c                                 intervals to cycle through
            do j = 1, intv(mode)
c                                 odd or even interval?
               odd = .not.odd
c                                 interval limits              
               ylmn = yint(j,i,ksite,ids)
               ylmx = yint(j+1,i,ksite,ids)
c                                 which interval are we starting from?
               if (y(i,iy(i)).gt.ylmx-tol) cycle
c
               dy = ylmx - ylmn
c                                 pathological case y = ylmn
               if (dabs(y(i,iy(i))-ylmn).lt.tol) y(i,iy(i)) = ylmn + tol

               if (res.eq.0d0) then 
c                                 the current value is in interval j
c                                 convert to raw y (varies from 0 ->1 
c                                 over the local interval)
                  yloc = (y(i,iy(i))-ylmn) / dy
c                                 convert to conformal x
                  if (odd) then 
                     x = unstch(yloc)
                  else 
                     x = 1d0 - unstch(1d0-yloc)
                  end if 

               else
c                                 have jumped from an earlier interval
                  x = res - ync / yfrc(j-1,i,ksite,ids)
c                 if (x.lt.0d0) x = 0d0

               end if                 
c                                 now generate all compositions in
c                                 local interval
               do while (x.le.1d0) 
c                                 increment conformal x
                  x = x + ync / yfrc(j,i,ksite,ids)
c                                 compute yreal
                  if (x.le.1d0) then 
                     if (odd) then 
                        yreal = ylmn + strtch(x) * dy
                     else
                        yreal = ylmx - strtch(1d0-x) * dy
                     end if 
c                                 cycle if below tolerance
                     if (yreal-y(i,iy(i)).lt.tol.and.good(ids)) then
                        good(ids) = .false.
                        write (*,1010) fname(ids)
                     end if 
 
                     iy(i) = iy(i) + 1
                     if (iy(i).gt.mres) call error (999,ync,mres,sname)
c                                 check if in bounds
                     if (dabs(yreal-ymax).lt.tol.or.yreal.gt.ymax) then
                        res = 0d0
                        y(i,iy(i)) = ymax
                        exit
                     else 
                        y(i,iy(i)) = yreal
                     end if 

                  else if (x.gt.1d0.and.j.eq.intv(mode)) then
c                                 at the last interval
                     iy(i) = iy(i) + 1
                     y(i,iy(i)) = ymax                
                     exit

                  else 
c                                the point is in the next interval
                     res = x - 1d0 
                     exit 

                  end if 
c                                 coordinate generating end do 
               end do 
c                                 if y is at ymax exit
               if (y(i,iy(i)).ge.ymax) exit
c                                 interval loop end do 
            end do 

         end if 
c                                 add last point if necessary, can it be? 
c                                 certainly not for conformal. 
         if (y(i,iy(i)).lt.ymax) then
            iy(i) = iy(i) + 1
            if (iy(i).gt.mres) call error (999,ync,mres,sname)
            y(i,iy(i)) = ymax
         end if          
 
      end do
c                                  
      do i = 1, jsp
         ind(i) = 1
      end do 
c                                 assign the first point
      npairs = 1

      do i = 1, jsp
         xy(i,1) = y(i,ind(i))
      end do
c                                 now make the array index run over all
c                                 values increasing the last index fastest
      iexit = 0 
      ieyit = 0 
      dy = 0d0

      do while (iexit.eq.0)
c                                 figure out which index to increment
         do i = jsp, 1, -1
            if (ind(i).lt.iy(i).and.ieyit.eq.0) then
c                                 this is the one to increment
               ind(i) = ind(i) + 1
               indx = i 
               exit 
            else if (i.gt.1) then 
c                                 saturated the index
               ind(i) = 1
               ieyit = 0 
               
            else
c                                 saturated first index, done.
                goto 99 

            end if 
         end do 
c                                 ok now we have the indices, check
c                                 the composition
         ycum = 0d0

         do i = 1, jsp 
            ycum = ycum + y(i,ind(i))  
         end do 

         if (ycum.gt.1d0) then

            ieyit = 1
c                                 here is where it gets messy:
            if (indx.eq.1) then
c                                 we're at the first point, and already
c                                 over the top
               iexit = 1
               cycle

            else if ( y(indx,ind(indx)) - y(indx,ind(indx)-1)
     *               - ycum + 1d0    .gt. tol ) then          
c                                 reset the current variable
c                                 and max loop index
               dy =  1d0 - ycum 

            else
c                                 must have just hit on the last increment
               cycle 

            end if 
         end if 

         npairs = npairs + 1

         if (npairs.gt.k21) call error (180,ycum,k21,'CARTES')

         do i = 1, jsp
            xy(i,npairs) = y(i,ind(i))
         end do

         xy(indx,npairs) = xy(indx,npairs) + dy

         dy = 0d0 

      end do 

1000  format (/,'Warning: a composition of solution ',a,' has been ',
     *          'refined to a level that',/,'may destabilize ',
     *          'calculations. If excessive failure occurs, use less ',
     *        /,'extreme iteration parameters.',/)
1010  format (/,'Warning: a composition of solution ',a,' has been ',
     *          'refined to a level that',/,'may destabilize ',
     *          'calculations. If excessive failure occurs use less ',
     *        /,'extreme iteration parameters or ',
     *          'increase the stretching factor.',/)
99    end


      subroutine blanko (text,chars,nchar,ilen)
c------------------------------------------------------------------- 
c unblnk - find the last non-blank character in text, 
c     text - character string 
c-------------------------------------------------------------------
      implicit none

      integer nchar, i, ilen

      character text*(*), chars(*)*1
 
      read (text,1000) (chars(i),i=1,ilen)
c                                 scan for blanks:
      do nchar = ilen, 1, -1
         if (chars(nchar).gt.' ') exit
      end do 

1000  format (400a1)

      end

      subroutine makecp (inames,mnames,vertex)
c----------------------------------------------------------------------
c makecp reads the thermodynamic to compute the composition vector of 
c made phases, called by vertex. programs without composition checking
c use smakcp.

c output to console if vertex = .true.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer inames, i, j, k,ict, id, incomp(k0), jct 

      logical inph(k16*k17), inmk(k16), eof, good, vertex

      double precision mcp(k16*k17,k0)
      character name*8, mnames(k16*k17)*8

      character dname*40,cmpnt*5
      common/ csta5 /dname,cmpnt(k0)

      double precision mcomp
      character mknam*8
      integer nmak
      logical mksat
      common / cst333 /mcomp(k16,k0),nmak,mksat(k16),mknam(k16,k17)

      double precision mkcoef, mdqf
      integer mknum, mkind
      common / cst334 /mkcoef(k16,k17),mdqf(k16,k17),mkind(k16,k17),
     *                 mknum(k16)

      integer idh2o,idco2,ikind,icmpn,icout
      double precision therm,comp,tot
      common/ cst43 /therm(k4),comp(k0),tot,icout(k0),idh2o,idco2,
     *               ikind,icmpn
c----------------------------------------------------------------------
c                                 make a list of all definition names:
      inames = 0 

      do i = 1, nmak

         inmk(i) = .true.

         do 10 j = 1, mknum(i)

            do k = 1, inames
               if (mnames(k).eq.mknam(i,j)) goto 10
            end do 

            inames = inames + 1
            mnames(inames) = mknam(i,j) 

10       end do 

      end do 
c                                 array to find valid makes:
      do i = 1, inames
         inph(i) = .false.
      end do 
c                                 now get the composition vectors for 
c                                 the mnames phases:
      eof = .false.

      do 

         call getphi (name,eof)

         if (eof) exit

         do i = 1, inames
            if (name.eq.mnames(i)) then 

               do j = 1, icmpn
                  mcp(i,j) = comp(j)
               end do 

               inph(i) = .true.

               exit
            end if 
         end do 

      end do 
c                                 find valid makes:
      do 40 i = 1, nmak
         do j = 1, mknum(i)
            do k = 1, inames
               if (mnames(k).eq.mknam(i,j).and.(.not.inph(k))) then 
                  inmk(i) = .false.
                  call warn (51,tot,icmpn,mknam(i,mknum(i)+1))
                  goto 40 
               else if (mnames(k).eq.mknam(i,j)) then 
                  mkind(i,j) = k
               end if 
            end do
         end do 
40    end do 
c                                 compute the composition for each 
c                                 made entitity and check if it's 
c                                 valid
      ict = 0 

      do i = 1, nmak

         if (inmk(i)) then

            name = mknam(i,mknum(i)+1)

            do j = 1, icmpn
               mcomp(i,j) = 0d0
            end do 
 
            do j = 1, mknum(i)
               id = mkind(i,j)
               do k = 1, icmpn
                  mcomp(i,k) = mcomp(i,k) + mkcoef(i,j)*mcp(id,k)
               end do 
            end do 
c                                 test the composition vector
c                                 is it a normal phase (i.e.,
c                                 non-zero thermodynamic components)
            do k = 1, icmpn
               comp(k) = mcomp(i,k) 
            end do 

            call chkphi (1,name,good)

            if (good) then 
               
               mksat(i) = .false.

            else 
c                                 no thermo componnents, sat comps?
               call chkphi (0,name,good)

               mksat(i) = .true.

               if (.not.good) then

                  inmk(i) = .false.

c                 call warn (52,tot,icmpn,mknam(i,mknum(i)+1))

                  cycle

               end if 

            end if 

         end if 
      end do 
c                                 clean up arrays:
      ict = 0

      do i = 1, icmpn
         incomp(i) = 0
      end do 
   

      do i = 1, nmak

         if (inmk(i)) then

            ict = ict + 1

            mksat(ict) = mksat(i)

            mknum(ict) = mknum(i)

            do j = 1, mknum(ict)+1
               mknam(ict,j) = mknam(i,j)
            end do 

            do j = 1, mknum(ict)
               mkcoef(ict,j) = mkcoef(i,j)
            end do 

            do j = 1, k17
               mdqf(ict,j) = mdqf(i,j)
            end do 

            do j = 1, icmpn
               mcomp(ict,j) = mcomp(i,j)
c                                get list of used components
               if (mcomp(ict,j).ne.0d0.and.incomp(j).eq.0) incomp(j) = j

            end do 

         end if

      end do  

      jct = 0 
      do j = 1, icmpn
         if (incomp(j).ne.0) then 
            jct = jct + 1
            incomp(jct) = incomp(j)
         end if
      end do 

      nmak = ict

      if (nmak.gt.0.and.vertex) 
     *               write (*,1010) (cmpnt(incomp(j)),j=1,jct)
c                                remake list of phases required for 
c                                makes:
      inames = 0 

      do i = 1, nmak

         do 20 j = 1, mknum(i)

            do k = 1, inames
               if (mnames(k).eq.mknam(i,j)) goto 20
            end do 

            inames = inames + 1
            mnames(inames) = mknam(i,j) 
   
20       end do
c                                write list of valid makes:
         if (vertex) write (*,1000) mknam(i,mknum(i)+1),
     *                             (mcomp(i,incomp(j)),j=1,jct)

      end do 

      if (nmak.gt.0.and.vertex) write (*,1020)

1000  format (a,1x,15(f5.2,1x))
1010  format (/,'Summary of valid make definitions:',//,10x,15(a,1x)) 
1020  format (/)

      end

      subroutine sattst (ifer)
c----------------------------------------------------------------------
c sorts phases into the appropriate saturated phase list called by 
c input2.             
c----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'
 
      character name*8

      integer iff,idss,ifug,ifyn,isyn,icomp,
     *        istct,iphct,icp,j,ifer,idc

      common/ cst10 /iff(2),idss(h5),ifug,ifyn,isyn/ csta6 /name
     *      / cst6 /icomp,istct,iphct,icp

      integer ic
      common/ cst42 /ic(k5)

      character cmpnt*5, dname*40
      common/ csta5 /dname,cmpnt(k0)

      integer idh2o,idco2,ikind,icmpn,icout
      double precision therm,comp,tot
      common/ cst43 /therm(k4),comp(k0),tot,icout(k0),idh2o,idco2,
     *               ikind,icmpn

      integer ids,isct,icp1,isat,io2
      common/ cst40 /ids(h5,h6),isct(h5),icp1,isat,io2
c-----------------------------------------------------------------------

         if (ifyn.eq.0) then
c                               check for fluid species data
            if (name.eq.cmpnt(idco2)) then
               j = 2
               ifer = ifer + 1
            else if (name.eq.cmpnt(idh2o)) then
               j = 1
               ifer = ifer + 1
            else
               goto 70
            end if
 
            call loadit (j)

         end if
 
70       if (isyn.eq.0) then 
c                               check for saturated phases:
c                               reject the phase if it contains
c                               a thermodynamic component:
            do j = 1, icp
               if (comp(ic(j)).ne.0d0) goto 99
            end do
c                               now load the phase if it has
c                               the saturated component idc:
            do j = isat, 1, -1
               idc = icp + j
               if (comp(ic(idc)).ne.0d0) then
                  isct(j) = isct(j) + 1 
                  if (isct(j).gt.h6) call error (17,1d0,h6,'SATTST')
                  iphct = iphct + 1
                  if (iphct.gt.k1) call error (180,1d0,k1,'SATTST')
                  ids(j,isct(j)) = iphct
                  call loadit (iphct)
                  goto 99
               end if
            end do 
         end if 

99    end 

      double precision function gmake (id)
c-----------------------------------------------------------------------
c gmake computes and sums the component g's for a make definition.
c the component g's may be calculated redundantly because gmake is
c called by gcpd, which in turn may be called by routines that call
c for a single g (e.g., gphase). 
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i, id, jd

      double precision g, dg

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      double precision mkcoef, mdqf

      integer mknum, mkind

      common / cst334 /mkcoef(k16,k17),mdqf(k16,k17),mkind(k16,k17),
     *                 mknum(k16)
c--------------------------------------------
      integer make

      common / cst335 /make(k10)
c-----------------------------------------------------------------------

      jd = make(id)

      g = 0d0
c                                compute the sum of the component g's
      do i = 1, mknum(jd)

         call gcpd (mkind(jd,i),dg)
         g = g + mkcoef(jd,i)*dg 

      end do 
c                                add the dqf correction
      gmake = g + mdqf(jd,1) + t*mdqf(jd,2) + p*mdqf(jd,3)

      end 

      subroutine smakcp (inames,mnames)
c----------------------------------------------------------------------
c smakcp reads the thermodynamic data file to compute the composition 
c vector of made phases without composition checking. vertex uses
c makecp
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer inames, i, j, k, ict, id, incomp(k5), jct 

      logical inph(k16*k17), inmk(k16), eof

      double precision mcp(k16*k17,k0)

      character name*8, mnames(k16*k17)*8

      character dname*40,cmpnt*5
      common/ csta5 /dname,cmpnt(k0)

      double precision mcomp
      character mknam*8
      integer nmak
      logical mksat
      common / cst333 /mcomp(k16,k0),nmak,mksat(k16),mknam(k16,k17)

      double precision mkcoef, mdqf
      integer mknum, mkind
      common / cst334 /mkcoef(k16,k17),mdqf(k16,k17),mkind(k16,k17),
     *                 mknum(k16)

      integer idh2o,idco2,ikind,icmpn,icout
      double precision therm,comp,tot
      common/ cst43 /therm(k4),comp(k0),tot,icout(k0),idh2o,idco2,
     *               ikind,icmpn
c----------------------------------------------------------------------
c                                 make a list of all definition names:
      inames = 0 

      do i = 1, nmak

         inmk(i) = .true.

         do 10 j = 1, mknum(i)

            do k = 1, inames
               if (mnames(k).eq.mknam(i,j)) goto 10
            end do 

            inames = inames + 1
            mnames(inames) = mknam(i,j) 

10       end do 

      end do 
c                                 array to find valid makes:
      do i = 1, inames
         inph(i) = .false.
      end do 
c                                 now get the composition vectors for 
c                                 the mnames phases:
      eof = .false.

      do 

         call getphi (name,eof)

         if (eof) exit

         do i = 1, inames
            if (name.eq.mnames(i)) then 

               do j = 1, icmpn
                  mcp(i,j) = comp(j)
               end do 

               inph(i) = .true.
 
               exit
            end if 
         end do 

      end do 
c                                 find valid makes:
      do 40 i = 1, nmak
         do j = 1, mknum(i)
            do k = 1, inames
               if (mnames(k).eq.mknam(i,j).and.(.not.inph(k))) then 
                  inmk(i) = .false.
                  call warn (51,tot,icmpn,mknam(i,mknum(i)+1))
                  goto 40 
               else if (mnames(k).eq.mknam(i,j)) then 
                  mkind(i,j) = k
               end if 
            end do
         end do 
40    end do 
c                                 compute the composition for each 
c                                 made entitity and check if it's 
c                                 valid
      ict = 0 

      do i = 1, nmak

         if (inmk(i)) then

            do j = 1, icmpn
               mcomp(i,j) = 0d0
            end do 
 
            do j = 1, mknum(i)
               id = mkind(i,j)
               do k = 1, icmpn
                  mcomp(i,k) = mcomp(i,k) + mkcoef(i,j)*mcp(id,k)
               end do 
            end do 
c                                 test the composition vector
c                                 is it a normal phase (i.e.,
c                                 non-zero thermodynamic components)
            do k = 1, icmpn
               comp(k) = mcomp(i,k) 
            end do 

         end if 
      end do 
c                                 clean up arrays:
      ict = 0

      do i = 1, icmpn
         incomp(i) = 0
      end do 
   

      do i = 1, nmak

         if (inmk(i)) then

            ict = ict + 1

            mknum(ict) = mknum(i)

            do j = 1, mknum(ict)+1
               mknam(ict,j) = mknam(i,j)
            end do 

            do j = 1, mknum(ict)
               mkcoef(ict,j) = mkcoef(i,j)
            end do 

            do j = 1, k17
               mdqf(ict,j) = mdqf(i,j)
            end do 

            do j = 1, icmpn
               mcomp(ict,j) = mcomp(i,j)
c                                get list of used components
               if (mcomp(ict,j).ne.0d0.and.incomp(j).eq.0) incomp(j) = j

            end do 

         end if

      end do  

      jct = 0 
      do j = 1, icmpn
         if (incomp(j).ne.0) then 
            jct = jct + 1
            incomp(jct) = incomp(j)
         end if
      end do 

      nmak = ict

      if (nmak.gt.0) write (*,1010) (cmpnt(incomp(j)),j=1,jct)
c                                remake list of phases required for 
c                                makes:
      inames = 0 

      do i = 1, nmak

         do 20 j = 1, mknum(i)

            do k = 1, inames
               if (mnames(k).eq.mknam(i,j)) goto 20
            end do 

            inames = inames + 1
            mnames(inames) = mknam(i,j) 
   
20       end do
c                                write list of valid makes:
         write (*,1000) mknam(i,mknum(i)+1),(mcomp(i,incomp(j)),j=1,jct)

      end do 

      write (*,1020)

1000  format (a,1x,15(f5.2,1x))
1010  format (/,'Summary of valid make definitions:',//,10x,15(a,1x)) 
1020  format (/)

      end

c       the x-files library for adaptive minimization

c local variables:
 
c insp(jstot/jstot+1) - pointer to original index of first kstot independent
c         endmembers, followed by jstot-kstot dependent endmembers, followed
c         by the ordered species, if present (jstot+1).
c istot - total number of endmembers (excluding ordered species) used
c         to formulate the solution model in the input file
c jmsol(jstot,msp) - chemical species present on m'th site of jstot'th
c          endmember, original indexing.
c jstot - total number of endmembers (excluding ordered species) used
c         in the computation (i.e., excluding missing endmembers), but
c         including dependent endmembers.
c kdsol(istot/istot+1) - endmember index, 0 if missing, -2 if dependent,
c          -1 if ordered species (istot+1/jstot+1), original indexing.
c kstot - total number of endmembers (excluding ordered species) 
c         used in the computation (i.e., excluding missing endmembers)

c global variables

c indx(i,j,k,l) - for solution i, pointer to the l'th original endmember
c                 index with species k on site j. 
c mstot(i) - istot globally
c jgsol(i,j,k) - k species indices of endmember j in solution i (jmsol globally) 

      subroutine reform (sname,im,vertex)
c---------------------------------------------------------------------
c reform - counts the number of species that can be respresented for a 
c solution given the present endmembers.
c---------------------------------------------------------------------
      implicit none
  
      include 'perplex_parameters.h'

      character*10 sname

      logical vertex

      integer kill,ikill,jkill,kill1,i,j,kosp(mst,msp),kill2,
     *        jmsol,k,kdsol,im,idsp,ksp(mst)

      common/ cst142 /jmsol(m4,mst),kdsol(m4)

      integer iend,isub,imd,insp,ist,isp,isite,iterm,iord,istot,jstot,
     *        kstot
      double precision wg,xmn,xmx,xnc
      common/ cst108 /wg(m1,m3),xmn(mst,msp),xmx(mst,msp),xnc(mst,msp),
     *      iend(m4),isub(m1,m2,2),imd(ms1,mst),insp(m4),ist(mst),
     *      isp(mst),isite,iterm,iord,istot,jstot,kstot
c----------------------------------------------------------------------
      kill = 1

      if (vertex.and.isite.gt.1) call warn (50,wg(1,1),isite,sname)

      do while (kill.lt.99) 

         do i = 1, isite
            ksp(i) = 1
            do j = 1, isp(i)
               kosp(i,j) = 0
            end do 
            do j = 1, isite
c                                 the number of endmembers that
c                                 should be present to represent
c                                 site i if all endmembers are present. 
               if (i.ne.j) ksp(i) = ksp(i)*isp(j)
            end do 
         end do
c                                 count the number of species
c                                 missing on each site
         do i = 1, istot
            if (kdsol(i).eq.0) then
               do j = 1, isite 
                  kosp(j,jmsol(i,j)) = kosp(j,jmsol(i,j)) + 1
               end do 
            end if 
         end do 
c                                 find the species that is missing
c                                 the most from the model
         kill = 99 
         kill1 = 0 

         do i = 1, isite
            if (isp(i).gt.1) then 
               do j = 1, isp(i)
c                                 idsp is the the number of species
c                                 possible - the number of missing
c                                 endmembers
                  idsp = ksp(i) -kosp(i,j)
                 
                  if (idsp.lt.kill) then 
                     ikill = i
                     jkill = j 
                     kill = idsp
                  else if (idsp.eq.kill.and.kosp(i,j).gt.0) then
c                                 kill the species that will kill the
c                                 most dependent endmembers
                     kill2 = 0 

                     do k = 1, istot
c                                 count the number that will be killed
                        if (jmsol(k,ikill).eq.jkill.and.kdsol(k).eq.-2) 
     *                     kill2 = kill2 + 1
                     end do

                     if (kill2.gt.kill1) then 
c                                 this is more than before (kill1)
                        kill1 = kill2 
                        ikill = i
                        jkill = j 
                        kill = idsp
                     end if 
                  end if 
               end do 
            end if 
         end do      
c                                 kill the species jkill on site ikill
c                                 and reformulate the model (this is 
c                                 inefficient, but we don't care here). 
         call killsp (ikill,jkill)
c                                 check exit conditions
         if (istot.lt.2) then
c                                 failed, rejected too many endmembers
            kill = 99 

            im = im - 1
            if (vertex) call warn (25,wg(1,1),jstot,sname)
            jstot = 0

         else if (istot.eq.jstot) then  
c                                 succeeded
            kill = 99 
c                                 reorder the insp array so that the 
c                                 the first kstot endmembers are 
c                                 independent, followed by the 
c                                 dependent endmember, followed by
c                                 the ordered species. 
         end if 

      end do 

      end 


      subroutine killsp (ikill,jkill)
c---------------------------------------------------------------------
c killsp - eliminates species jkill from site ikill in a solution model
c and reformulates the model accordingly
c---------------------------------------------------------------------
      implicit none
  
      include 'perplex_parameters.h'

      logical skip, bad

      integer jsp,jtic,morder,
     *        i,j,ikill,jkill,kill,kdep,jdqf,ktic,jold,
     *        iwas(m4),jmsol,i2ni(m4),kwas(m4),
     *        iddeps,jsmod,k,l,itic,kdsol,ijkill(m4),norder,
     *        j2oj(msp),j2nj(msp),i2oi(m4),maxord,mord

      double precision depvnu,denth,vlaar
c                                 dqf variables
      integer indq,idqf
      double precision dqf
      common/ cst222 /dqf(m3,m4),indq(m4),idqf
c
      logical dead,depend,laar,order,fluid,macro,specil
      common/ cst160 /depend,laar,order,fluid,macro,specil
c                                 local input variables
      common/ cst221 /vlaar(m3,m4),jsmod
     *      / cst141 /depvnu(2,j3),denth(j3),iddeps(2,j3),norder
     *      / cst142 /jmsol(m4,mst),kdsol(m4)

      integer nsub1,nsub2,nttyp,nterm,nspm1,nsite
      double precision acoef,smult,a0
      common/ cst107 /a0(m10,m11),acoef(m10,m11,m0),smult(m10),
     *      nsite,nspm1(m10),nterm(m10,m11),nsub1(m10,m11,m0,m12),
     *      nsub2(m10,m11,m0,m12),nttyp(m10,m11,m0)

      integer iend,isub,imd,insp,ist,isp,isite,iterm,iord,istot,jstot,
     *        kstot
      double precision wg,xmn,xmx,xnc
      common/ cst108 /wg(m1,m3),xmn(mst,msp),xmx(mst,msp),xnc(mst,msp),
     *      iend(m4),isub(m1,m2,2),imd(ms1,mst),insp(m4),ist(mst),
     *      isp(mst),isite,iterm,iord,istot,jstot,kstot

      integer iorig,jnsp,iy2p
      common / cst159 /iorig(m4),jnsp(m4),iy2p(m4)

      integer mdep,idep,jdep,ndph
      double precision nu,y2p
      common/ cst146 /nu(m15,j4),y2p(m4,m15),mdep,jdep(m15),
     *                idep(m15,j4),ndph(m15)

      double precision yin
      common/ cst50 /yin(ms1,mst)

      integer limn,limt,limid,jimid,jimt
      double precision limc,jimc
      common/ cxt30 /limc(j6+2,j5,j3),limid(m0,j5,j3),jimid(j3,j5,j3),
     *               limn(j3),limt(j5,j3),jimc(j3,j5,j3),jimt(j5,j3)
c----------------------------------------------------------------------
      do i = 1, isite

         if (i.ne.ikill) then
c                                 nothing happens  
            do j = 1, isp(i)  
              j2oj(j) = j 
            end do
         else 
c                                 on a site where we kill a species
            jsp = isp(i) - 1
c                                 should also check and store subdivsion
c                                 model here (i.e., some ternary model
c                                 choices are not valid as binary):

c                                 make a pointer from the new species
c                                 index to the old index
            jtic = 0 

            do j = 1, isp(i)
               if (j.ne.jkill) then 
                  jtic = jtic + 1 
c                              pointer from new j to old j
                  j2oj(jtic) = j
c                              pointer from old j to new j
                  j2nj(j) = jtic
               end if 
            end do
c                              now reload
            isp(i) = jsp

            if (jsp.gt.1) then   
c                              now shift subdivision ranges
               do j = 1, jsp - 1
                  xmn(i,j) = xmn(i,j2oj(j))
                  xmx(i,j) = xmx(i,j2oj(j))
                  xnc(i,j) = xnc(i,j2oj(j))
                  imd(j,i) = imd(j2oj(j),i)

                  if (imd(j,i).gt.0) yin(j,i) = yin(j2oj(j),i)

               end do
            else 
               xmn(i,j) = 1d0
               xmx(i,j) = 1d0
               xnc(i,j) = 1d0                 
            end if
         end if 
      end do 

      kdep = 0 

      do i = 1, istot

         if (depend.and.kdsol(i).eq.-2) then 
c                                create an array which gives the
c                                original locations of the dependent
c                                endmembers, need this to be able to
c                                reorder the y2p array:
            kdep = kdep + 1
            iwas(kdep) = i 
         end if 
c                                 kill endmembers with the species
c                                 to be deleted:
         if (jmsol(i,ikill).eq.jkill) kdsol(i) = -3

      end do 
c                                 check for dependent endmembers
      call redep (-3)
c                                 at this point all ordered endmembers to 
c                                 be killed are flagged by kdsol = -3.

c                                 now check the ordered species
      morder = 0 

      if (order) then 
c                                 first check if the ordered endmember
c                                 may be stable 
         do k = 1, norder      
c                                 check if a missing constituent 
            bad = .false.

            do j = 1, 2
               if (kdsol(iddeps(j,k)).eq.-3) then 
                  bad = .true.
                  exit 
               end if 
            end do

            if (bad) then 
c                                 add species to the kill list
               kdsol(istot+k) = -3
            
            else 

               morder = morder + 1
               kdsol(istot+k) = -1
               kwas(morder) = k

            end if 

         end do  

      end if 
c                                figure out which dependent endmembers have
c                                been killed:
      if (depend) then 

         if (kdep.ne.mdep) call error (54,dqf(1,1),mdep,'KILLSP')

         mdep = 0 

         do i = 1, kdep
            if (kdsol(iwas(i)).ne.-3) then 
               mdep = mdep + 1
c                                 iwas is now the original index of the 
c                                 dependent endmember, and mdep is the reset
c                                 value of the dependent endmember counter
               iwas(mdep) = i                
            end if
         end do 

      end if 
 
      itic = 0 
      jtic = 0
      ktic = 0
      kill = 0 
      
      do i = 1, istot + norder 

         if (kdsol(i).ge.-2) then 
c                                 replacement for istot (itic)
            itic = itic + 1
c                                 pointers from new to old endmember index (i2oi)
            i2oi(itic) = i
c                                 pointers from new to old endmember index (i2ni)
            i2ni(i) = itic                                
c                                 pointer to original species index
            iorig(itic) = iorig(i)
c                                 number of missing endmembers (jtic)
            if (kdsol(i).eq.0) jtic = jtic + 1
c                                reset the kdsol array
            kdsol(itic) = kdsol(i)

            if (i.gt.istot) cycle
c                                 reset the species pointers (jmsol) 
            do j = 1, isite
               if (j.eq.ikill) then
                  jmsol(itic,j) = j2nj(jmsol(i,j))
               else
                  jmsol(itic,j) = jmsol(i,j)
               end if 
            end do
         else 
c                                 kill records the killed endmembers
            kill = kill + 1
            ijkill(kill) = i

         end if 
      end do  
c                                reset total and present counters
      istot = itic - morder 
c            
      jstot = itic - jtic - morder
c                                --------------------------------------
c                                excess terms:
      itic = 0 
      maxord = 0 
    
      do i = 1, iterm 
c                                check for forbidden terms (i.e., terms
c                                with a missing endmember
         skip = .false.

         if (jsmod.eq.1) then 
c                                microscopic formulation
            do j = 1, iord
               if (isub(i,j,1).eq.ikill.and.isub(i,j,2).eq.jkill) then
                  skip = .true.
                  exit
               end if 
            end do 

            if (skip) cycle
c                                 the term is acceptable
            itic = itic + 1

            mord = 0 

            do j = 1, iord

               isub(itic,j,1) = isub(i,j,1)

               if (isub(i,j,1).eq.ikill) then 
                  isub(itic,j,2) = j2nj(isub(i,j,2))
               else
                  isub(itic,j,2) = isub(i,j,2)
               end if 

               if (isub(i,j,1).eq.0) exit

               mord = mord + 1

            end do 
         else 
c                                 macroscopic formulation
            do j = 1, kill
c                                 check if subscript points to a killed 
c                                 endmember
               do k = 1, iord
                  if (isub(i,k,1).eq.ijkill(j)) then
                     skip = .true.
                     exit 
                  end if 
               end do 

               if (skip) exit

            end do 

            if (skip) cycle 
c                               the term is acceptable
            itic = itic + 1

            mord = iord

            do j = 1, iord
               isub(itic,j,1) = i2ni(isub(i,j,1))
            end do 

         end if 
c                                save the coefficient
         do j = 1, m3
            wg(itic,j) = wg(i,j)
         end do 
c                                find highest order term
         if (mord.gt.maxord) maxord = mord

      end do     
c                                reset counters, iord is not reset
      iterm = itic
      iord = maxord                         
c                                --------------------------------------
c                                van laar volume functions
      if (laar) then
         do i = 1, istot + morder 
            do j = 1, m3
               vlaar(j,i) = vlaar(j,i2oi(i))
            end do 
         end do  
      end if 
c                                 --------------------------------------
c                                 dqf corrections, this is sloppy since
c                                 uses istot instead of kstot
      if (idqf.gt.0) then 

         jdqf = 0 
c                                 check if a retained species has a dqf
c                                 correction
         do j = 1, idqf
c                                 the itoi index must be in the inner loop
c                                 in case the values of indq are not sequential
            do i = 1, istot
               if (indq(j).eq.i2oi(i)) then 
c                                 found a dqf'd endmember
                  jdqf = jdqf + 1
                  indq(jdqf) = i
                  do k = 1, m3
                     dqf(k,jdqf) = dqf(k,j)
                  end do 
                  exit 
               end if
            end do 

            if (jdqf.eq.idqf) exit
 
         end do 

         idqf = jdqf 

      end if 
c                                 --------------------------------------
c                                 configurational entropy model

c                                 site fractions as a function of bulk
c                                 y's and dependent species y:
      do i = 1, nsite 
c                                 for each species, read function to define 
c                                 the site fraction of the species and eliminate
c                                 killed species

c                                 species counter is incremented in advance
c                                 and must be decremented before saving the 
c                                 final value:
         jtic = 1 

         do j = 1, nspm1(i)

            ktic = 0 
c                                 for each term:
            do k = 1, nterm(i,j)

               if (jsmod.eq.1) then 
c                                 microscopic formulation
                  dead = .false.

                  do l = 1, nttyp(i,j,k)
                     if (nsub1(i,j,k,l).eq.ikill) then 
c                                 on the site where we kill species jkill
                        if (nsub2(i,j,k,l).eq.jkill) then 
c                                 the term has the killed species
                           dead = .true.
                           exit 
                        end if 
                     end if 
                  end do
c                                 if not dead, count the term
                  if (.not.dead) then 
                     ktic = ktic + 1
c                                 nttyp cannot change
                     nttyp(i,jtic,ktic) = nttyp(i,j,k)
                     acoef(i,jtic,ktic) = acoef(i,j,k)
                     do l = 1, nttyp(i,j,ktic)
                        nsub1(i,jtic,ktic,l) = nsub1(i,j,k,l)
                        if (nsub1(i,j,ktic,l).eq.ikill) then 
                           nsub2(i,jtic,ktic,l) = j2nj(nsub2(i,j,k,l))
                        else 
                           nsub2(i,jtic,ktic,l) = nsub2(i,j,k,l)
                        end if 
                     end do 
                  end if 

               else 
c                                 macroscopic formulation:
c                                 note: 4th index (nttyp) is only used
c                                 for bragg-williams models.
                  dead = .false.
                  do l = 1, kill 
                     if (nsub1(i,j,k,1).eq.ijkill(l)) then 
                        dead = .true.
                        exit 
                     end if 
                  end do 

                  if (.not.dead) then
c                                 the term has survived (and therefore 
c                                 also the species):
c                                 don't save nttyp since this is always
c                                 1 for the macroscopic formulation
                     ktic = ktic + 1
c                                 but my dear peanut brained friend, do 
c                                 not forget to move the pointer:
                     nsub1(i,jtic,ktic,1) = i2ni(nsub1(i,j,k,1)) 
                     acoef(i,jtic,ktic) = acoef(i,j,k)
                  end if 
               end if 
            end do 
c                                 ktic is the number of terms representing
c                                 the jth species, we won't count species 
c                                 with no terms because the endmember configurational
c                                 entropy is assumed to be implicit. 
            if (ktic.gt.0) then 
c                                 increment the species counter
               nterm(i,jtic) = ktic 
               a0(i,jtic) = a0(i,j)
               jtic = jtic + 1
            end if 

         end do

         nspm1(i) = jtic - 1

      end do 
c                                 ---------------------------------------
c                                 ordered species:
      if (order) then 

         norder = morder

         if (morder.eq.0) then 
c                                 there are no ordered species left
            order = .false.

            if (depend) then
               jsmod = 7
            else 
               jsmod = 2
            end if 

         else 
c                                 shift the ordered species pointers 
c                                 and data to eliminate kill ordered
c                                 species.    
            do j = 1, morder 

               jold = kwas(j)
               denth(j) = denth(jold)

               do i = 1, 2
                  iddeps(i,j) = i2ni(iddeps(i,jold))
                  depvnu(i,j) = depvnu(i,jold)
               end do

               itic = 1 

               do i = 1, limn(jold)
c                                 eliminate absent species from 
c                                 stoichiometric p0 limits
                  ktic = 0 

                  do k = 1, limt(i,jold)

                     skip = .false.

                     do l = 1, kill
c                                 check if limid points to a killed 
c                                 endmember
                        if (limid(k,i,jold).eq.ijkill(l)) then
                           skip = .true.
                           exit 
                        end if
  
                     end do  

                     if (skip) cycle

                     ktic = ktic + 1

                     limid(ktic,itic,j) = i2ni(limid(k,i,jold))
                     limc(ktic,itic,j) = limc(k,i,jold)

                  end do 

                  if (ktic.eq.0) cycle

                  limc(ktic+1,itic,j) = limc(k,i,jold)
                  limc(ktic+2,itic,j) = limc(k+1,i,jold)
                  limt(itic,j) = ktic
c                                 now check the p terms, this assumes
c                                 there are no p terms if there are no
c                                 p0 terms, which maynot be true?
                  ktic = 0 

                  do k = 1, jimt(i,jold)

                     skip = .false.

                     do l = 1, kill
c                                 check if limid points to a killed 
c                                 endmember
                        if (jimid(k,i,jold).eq.ijkill(l)) then
                           skip = .true.
                           exit 
                        end if
  
                     end do  

                     if (skip) cycle

                     ktic = ktic + 1

                     jimid(ktic,itic,j) = i2ni(jimid(k,i,jold))
                     jimc(ktic,itic,j) = jimc(k,i,jold)

                  end do 

                  jimt(itic,j) = ktic

                  itic = itic + 1

               end do 

               limn(j) = itic - 1

            end do 

         end if  

      end if 
c                                 --------------------------------------
c                                 dependent endmember properties, the
      if (depend) then 
c                                 dependent endmembers have been reordered
c                                 in redep, but are still expressed in
c                                 terms of the old indices, so reset the 
c                                 indices:
         do i = 1, mdep
            jdep(i) = i2ni(jdep(i))
            do j = 1, ndph(i)
               idep(i,j) = i2ni(idep(i,j))
            end do 
         end do 

      end if 

      end

      subroutine cmodel (im,idsol,tname,ibuild,x1,x2,first)
c---------------------------------------------------------------------
c cmodel - checks to see if solution models contain valid endmembers.
c modified to allow saturated phase/component endmembers, 10/25/05.
c---------------------------------------------------------------------
      implicit none
  
      include 'perplex_parameters.h'

      logical first
 
      character*10 mname*8,tname,fname,names*8,missin(m4)*8,x1*8,x2*8

      integer imiss,im,idsol,ibuild,i,j,imyn,
     *        ipoint,h,isoct,jmsol,jsmod,kdsol

      double precision vlaar
 
      common/ cst8 /names(k1)/ cst60 /ipoint,imyn/ cst79 /isoct
     *      / csta7 /fname(h9)   
     *      / cst221 /vlaar(3,m4),jsmod
     *      / cst142 /jmsol(m4,mst),kdsol(m4)/ cst18a /mname(m4)


      integer iend,isub,imd,insp,ist,isp,isite,iterm,iord,istot,jstot,
     *        kstot
      double precision wg,xmn,xmx,xnc
      common/ cst108 /wg(m1,m3),xmn(mst,msp),xmx(mst,msp),xnc(mst,msp),
     *      iend(m4),isub(m1,m2,2),imd(ms1,mst),insp(m4),ist(mst),
     *      isp(mst),isite,iterm,iord,istot,jstot,kstot

      integer mdep,idep,jdep,ndph
      double precision nu,y2p
      common/ cst146 /nu(m15,j4),y2p(m4,m15),mdep,jdep(m15),
     *                idep(m15,j4),ndph(m15)

      integer ikp
      common/ cst61 /ikp(k1)
c----------------------------------------------------------------------
      jstot = 0
c                              if called by build (ibuild = 1) skip the
c                              name check:
      if (ibuild.eq.1) goto 60
c                              check that solution name matches a 
c                              name requested in input from n1
      do i = 1, isoct

         if (tname.eq.fname(i)) then 
c                              got a match, goto 60:
            idsol = i 
            im = im + 1
            goto 60 
         end if 
c
      end do 
c                              didn't find a match, read a new name:
      goto 99 
c                              check that the endmembers match with data
c                              from n2:   
60    do 70 i = 1, istot

         kdsol(i) = 0

         if (jsmod.eq.7.or.jsmod.eq.8) then
c                              reciprocal solution, check if endmember i
c                              is a dependent endmember, if it is flag
c                              it by setting kdsol(i) = -2
            do j = 1, mdep
               if (jdep(j).eq.i) then
                  kdsol(i) = -2
                  goto 70
               end if 
            end do 
         end if 

         do h = 1, ipoint
c                                 in build, check for forbidden choices:
            if (ibuild.eq.1) then
               if (mname(i).eq.x1.or.mname(i).eq.x2) then
                  write (*,1010) tname, x1, x2
                  jstot = 0
                  goto 99
               end if
            end if

            if (names(h).eq.mname(i)) then
c                                 got a valid endmember, count and                     
               jstot = jstot + 1 
               kstot = jstot
c                                 create arrays of convenience, where j = 1, jstot
               kdsol(i) = h                                 
c                                 set kill flag: 
c                                 don't reset ikp if it has been set for 
c                                 another model.
               if (iend(i).ne.0.and.ikp(h).eq.0) ikp(h) = -1
c                                 found all endmembers:
               if (jstot.eq.istot) goto 99

               goto 70
            end if 
         end do
70    continue

      call redep (0)
c                                done if nothing is missing:
      if (jstot.eq.istot) goto 99
c                                missing endmember warnings:
      if (jstot.lt.2) then
         im = im - 1
         if (first) call warn (25,wg(1,1),jstot,tname)
         jstot = 0 
      else 
         imiss = 0 
         do i = 1, istot
            if (kdsol(i).eq.0) then
               imiss = imiss + 1
               missin(imiss) = mname(i)
            end if 
         end do 
         if (first) write (*,1000) tname,(missin(i), i = 1, imiss)
      end if 
 

1000  format (/,'**warning ver114** the following endmembers',
     *          ' are missing for ',a,/,4(8(2x,a),/),/)
1010  format (/,'**warning ver113** ',a,' is not a valid model',
     *        ' because component ',a,' or ',a,/,
     *        ' is constrained or missing.',/)
99    end

      subroutine redep (jkill)
c----------------------------------------------------------------------
c redep does reordering of dependent endmembers

c jkill is the value of kdsol that indicates invalid (missing) endmembers,
c this is 0 from cmodel and -3 from reform.
c----------------------------------------------------------------------
      implicit none
  
      include 'perplex_parameters.h'

      integer i,j,l,jmsol,ndep,jsmod,kdsol,k,jkill

      double precision vlaar

      logical depend,laar,order,fluid,macro,specil

      common/ cst160 /depend,laar,order,fluid,macro,specil
 
      common/ cst221 /vlaar(3,m4),jsmod
     *      / cst142 /jmsol(m4,mst),kdsol(m4)

      integer iend,isub,imd,insp,ist,isp,isite,iterm,iord,istot,jstot,
     *        kstot
      double precision wg,xmn,xmx,xnc
      common/ cst108 /wg(m1,m3),xmn(mst,msp),xmx(mst,msp),xnc(mst,msp),
     *      iend(m4),isub(m1,m2,2),imd(ms1,mst),insp(m4),ist(mst),
     *      isp(mst),isite,iterm,iord,istot,jstot,kstot

      integer iddeps,norder 
      double precision depvnu,denth
      common/ cst141 /depvnu(2,j3),denth(j3),iddeps(2,j3),norder

      integer mdep,idep,jdep,ndph
      double precision nu,y2p
      common/ cst146 /nu(m15,j4),y2p(m4,m15),mdep,jdep(m15),
     *                idep(m15,j4),ndph(m15)
c----------------------------------------------------------------------

c                                check for dependent endmembers, necessary?
      if (depend) then

         ndep = 0  

         do 100 i = 1, mdep

            do j = 1, ndph(i)

               if (idep(i,j).le.istot) then
c                                looking for a dependent endmember component: 
                  if (kdsol(idep(i,j)).eq.jkill.and.
     *                kdsol(jdep(i)).ne.-3) then
c                                a component is missing in cmodel, reset kdsol to 0
                     kdsol(jdep(i)) = 0
                     goto 100 
                  end if 

               else
c                                looking for an ordered species
                  do l = 1, norder 
                     do k = 1, 2
                        if (kdsol(iddeps(k,l)).eq.jkill) then
                           kdsol(jdep(i)) = 0 
                           goto 100
                        end if 
                     end do 
                  end do 

               end if 
            end do 
c                                dependent endmember is ok
            ndep = ndep + 1
            jdep(ndep) = jdep(i)
            ndph(ndep) = ndph(i)

            do j = 1, ndph(i)
               idep(ndep,j) = idep(i,j)
               nu(ndep,j) = nu(i,j)
            end do 

            jstot = jstot + 1

100      continue 

         mdep = ndep

         if (mdep.eq.0) then

            depend = .false.

            if (jsmod.eq.7) then
               jsmod = 2
               if (laar) jsmod = 3
            else 
               jsmod = 6
            end if 
         end if  

      end if 

      end 

      subroutine rmodel (tname,bad,vertex)
c---------------------------------------------------------------------
c rmodel - reads solution models from LUN n9.
c---------------------------------------------------------------------
      implicit none
  
      include 'perplex_parameters.h'
 
      character*10 tname, tag*3

      integer nreact,i,j,k,l,m,jlaar,jmsol,jsmod,kdsol,idim

      logical bad,vertex

      double precision vlaar,coeffs(k7),rnums(100),enth

      integer ijk(mst),inds(k7),ict

      common/ cst221 /vlaar(m3,m4),jsmod
     *      / cst142 /jmsol(m4,mst),kdsol(m4)

      integer nsub1,nsub2,nttyp,nterm,nspm1,nsite
      double precision acoef,smult,a0
      common/ cst107 /a0(m10,m11),acoef(m10,m11,m0),smult(m10),
     *      nsite,nspm1(m10),nterm(m10,m11),nsub1(m10,m11,m0,m12),
     *      nsub2(m10,m11,m0,m12),nttyp(m10,m11,m0)

      integer iend,isub,imd,insp,ist,isp,isite,iterm,iord,istot,jstot,
     *        kstot
      double precision wg,xmn,xmx,xnc
      common/ cst108 /wg(m1,m3),xmn(mst,msp),xmx(mst,msp),xnc(mst,msp),
     *      iend(m4),isub(m1,m2,2),imd(ms1,mst),insp(m4),ist(mst),
     *      isp(mst),isite,iterm,iord,istot,jstot,kstot

      logical depend, laar,order,fluid,macro,specil
      common/ cst160 /depend,laar,order,fluid,macro,specil

      integer iorig,jnsp,iy2p
      common / cst159 /iorig(m4),jnsp(m4),iy2p(m4)

      integer iddeps,norder 
      double precision depvnu,denth
      common/ cst141 /depvnu(2,j3),denth(j3),iddeps(2,j3),norder

      integer mdep,idep,jdep,ndph
      double precision nu,y2p
      common/ cst146 /nu(m15,j4),y2p(m4,m15),mdep,jdep(m15),
     *                idep(m15,j4),ndph(m15)

      integer maxit 
      common/ cst149 /maxit

      integer inew
      common/ cst140 /inew

      double precision yin
      common/ cst50 /yin(ms1,mst)
c----------------------------------------------------------------------

      mdep = 0 
      istot = 0
      ist(1) = 0 
c                               set logical variables
      specil = .false.
      macro = .false.
      order = .false.
      laar = .false.
      depend = .false.
      fluid = .false.
      bad = .false.

      do 
          read (n9, 1000, end = 99) tname
          if (tname.ne.' '.and.tname.ne.'begin_mode') exit
      end do 

      call readda (rnums,1,tname)
c                               model type flag
      jsmod = idint(rnums(1))   
c                               correct jsmod for old versions    
      if (jsmod.eq.3) jsmod = 2  


      if (jsmod.eq.0) fluid = .true.
      if (jsmod.gt.1) macro = .true.
      if (jsmod.eq.6.or.jsmod.eq.8) order = .true.
      if (jsmod.eq.7.or.jsmod.eq.8) depend = .true. 
      if (jsmod.gt.20) specil = .true.  

      if (jsmod.eq.1) then 
c                               old style solution models
c                               to put old style solution models
c                               in current version delete jfix line at end
c                               and add model type (1) after name. 
         call xrmode (tname)
c                               look for the end of model marker
c                               added in perplex04  
         call readop (idim,jlaar,istot,maxit,tname)

      else 
c                               new style:

c                               read number of independent sites:
      if (inew.eq.0.or.jsmod.eq.7.or.jsmod.eq.8) then 
         call readda (rnums,1,tname)    
         isite = idint(rnums(1)) 
      else
         isite = 1
      end if 
c                               read number of species on each site:
      call readda (rnums,isite,tname)  
      do i = 1, isite
         isp(i) = idint(rnums(i))
      end do 
c                               total number of endmembers:
      istot = 1
      do i = 1, isite
         istot = istot*isp(i)
      end do 
c                               counter for character read routines
      idim = istot
c                               read endmember names:
      call readn (idim,tname)
c                               compound formation models
      norder = 0 

      if (order) then 
c                               get the number of ordered species
         if (inew.eq.0) then
            norder = 1
         else 
            call readda (rnums,1,tname)    
            norder = idint(rnums(1)) 
         end if 

         if (norder.gt.j3) call error (5,rnums(1),norder,tname)
c                               nreact can only be 3 for ordered
c                               species:
         nreact = 3
c                               get ordering reaction and name
c                               of ordered species:   
         do i = 1, norder   

            call readr (coeffs,enth,inds,idim,nreact,tname)

            denth(i) = enth
  
            do j = 1, 2
               depvnu(j,i) = coeffs(j+1)
               iddeps(j,i) = inds(j+1)
            end do

         end do  
c                               read the limit equations for the 
c                               amount of the ordered endmembers
         call readlm (tname,bad,vertex)

      end if 
c                               read dependent endmembers
      if (depend) then
c                               number of dependent endmembers
         call readda (rnums,1,tname)    
         mdep = idint(rnums(1))
         if (mdep.gt.m15) call error (1,xmn(1,1),mdep,'m15')

         do i = 1, mdep
c                               nreact is returned by readr
            nreact = 0 

            call readr (coeffs,enth,inds,idim,nreact,tname)

            jdep(i) = inds(1)
            ndph(i) = nreact - 1

            do j = 1, ndph(i)
               nu(i,j) = coeffs(j+1)
               idep(i,j) = inds(j+1)
            end do 

         end do 
      end if 
c                               read endmember flags:
      call readda (rnums,istot,tname)  

      do i = 1, istot
         iend(i) = idint(rnums(i))
      end do 
c                               read composition limits, subdivision type:
      m = 0

      do i = 1, isite
c                               number of ranges to be read
         m = m + isp(i) - 1
      end do 
c                               get the numbers
      if (inew.eq.0) then 
c                               old format
         call readda (rnums,3*m+isite,tname) 
c                               assign the numbers:
         l = 1

         do i = 1, isite
         
            do j = 1, isp(i) - 1
               xmn(i,j) = rnums(l)
               xmx(i,j) = rnums(l+1)
               xnc(i,j) = rnums(l+2)
               l = l + 3
            end do 
c                                 fill the imd array
            do j = 1, isp(i) - j
               imd(j,i) = idint(rnums(l))
            end do 

            l = l + 1

         end do

      else 
c                                 new format
         l = 1

         do i = 1, isite
         
            do j = 1, isp(i) - 1

               call readda (rnums,4,tname)

               xmn(i,j) = rnums(1)
               xmx(i,j) = rnums(2)
               xnc(i,j) = rnums(3)
               imd(j,i) = idint(rnums(4))

               if (imd(j,i).eq.3) then
c                                 read extra parm
                  call readda (rnums,1,tname)
                  yin(j,i) = rnums(1)
               else if (imd(j,i).gt.4) then
                  call error (169,rnums(1),imd(j,i),tname)
               end if 

            end do 

         end do

      end if 
c                                create bragg-williams indexes
      do i = 2, isite
         ijk(i) = 1
      end do 
 
      ijk(1) = 0  

      do 20 l = 1, istot

         do m = 1, isite

            if (ijk(m).lt.isp(m)) then 

               ijk(m) = ijk(m) + 1
c                                increment only one index per endmember
               do i = 1, isite
                  jmsol(l,i) = ijk(i)
               end do 

               if (ijk(m).eq.isp(m).and.l.lt.istot-isp(1)+1) then
c                                saturated counter, increment first
c                                unsaturated counter and reset all
c                                lower site counters to first index,
                  do j = 2, isite
                     if (ijk(j).lt.isp(j)) then
                        ijk(j) = ijk(j) + 1
                        do k = 2, j-1
                           ijk(k) = 1
                        end do 
                        ijk(1) = 0
                        goto 20 
                     end if
                  end do  
               end if 

               goto 20
            end if
         end do 
20    continue 
c                              read excess function
      call readx (idim,tname)
c                              expansion for S(configurational)
      call readda (rnums,1,tname)    

      nsite = idint(rnums(1))

      if (nsite.gt.m10) call error (31,a0(1,1),isite,tname)
c                                 for each site
      do i = 1, nsite 
c                                 read # of species, and site 
c                                 multiplicty.
         call readda (rnums,2,tname)    
         nspm1(i) = idint(rnums(1)) - 1
         smult(i) = rnums(2)
c                                 for each species, read
c                                 function to define the
c                                 site fraction of the species:
         do j = 1, nspm1(i)
c                                 read expression for site
c                                 fraction of species j on 
c                                 site i.
            call readz (coeffs,inds,ict,idim,tname,tag)

            a0(i,j) = coeffs(1)
            nterm(i,j) = ict - 1
            if (nterm(i,j).gt.m0) call error (33,a0(1,1),m0,tname)
c                                 for each term:
            do k = 2, ict
c                                 all terms 1 order type, this
c                                 saved for compatability with 
c                                 old versions:
               nttyp(i,j,k-1)   = 1
               acoef(i,j,k-1)   = coeffs(k)
               nsub1(i,j,k-1,1) = inds(k)
            end do 
         end do 
      end do 
c                              look for van laar and/or dqf parameters
c                              or the end of model marker

      call readop (idim,jlaar,istot-mdep,maxit,tname)

      if (jlaar.ne.0) then

         laar = .true.
c                                 high order terms not allowed for
c                                 van laar.
         if (iord.gt.2.and.laar) call error (999,coeffs(1),800,'RMODEL')
c                                 re-set jsmod flag only for jsmod 2
         if (jsmod.eq.2) jsmod = 3

      end if 
c                                 save original indices, need this for 
c                                 melt models etc that have species specific
c                                 equations of state.
      end if 

      do i = 1, istot + norder
         iorig(i) = i 
      end do

1000  format (3(a,1x))

99    end 

      subroutine nmodel 
c---------------------------------------------------------------------
c nmodel - creates some counters and conversion arrays and stores
c local solution model parameters in global arrays. 
c---------------------------------------------------------------------
      implicit none
  
      include 'perplex_parameters.h'
 
      integer itic,i,j,k

      logical depend,laar,order,fluid,macro,specil

      common/ cst160 /depend,laar,order,fluid,macro,specil

      integer iend,isub,imd,insp,ist,isp,isite,iterm,iord,istot,jstot,
     *        kstot
      double precision wg,xmn,xmx,xnc
      common/ cst108 /wg(m1,m3),xmn(mst,msp),xmx(mst,msp),xnc(mst,msp),
     *      iend(m4),isub(m1,m2,2),imd(ms1,mst),insp(m4),ist(mst),
     *      isp(mst),isite,iterm,iord,istot,jstot,kstot

      integer jmsol,kdsol
      common/ cst142 /jmsol(m4,mst),kdsol(m4)

      integer iorig,jnsp,iy2p
      common / cst159 /iorig(m4),jnsp(m4),iy2p(m4)

      integer mdep,idep,jdep,ndph
      double precision nu,y2p
      common/ cst146 /nu(m15,j4),y2p(m4,m15),mdep,jdep(m15),
     *                idep(m15,j4),ndph(m15)

      integer iddeps,norder 
      double precision depvnu,denth
      common/ cst141 /depvnu(2,j3),denth(j3),iddeps(2,j3),norder
c----------------------------------------------------------------------
c                                 make the insp arrays, 

c                                 insp orders endmembers; first independent
c                                 then dependent, last ordered species.
 
c                                 jnsp points only to independent+ordered
c                                 species (used for reciprocal solutions 
c                                 with ordering).

c                                 iy2p points from an endmember in the y
c                                 array (1..istot+norder) to the endmember
c                                 in the p array (1..kstot+norder)

      if (depend) then       

         itic = 0 

         do i = 1, istot
            if (kdsol(i).gt.0) then 
               itic = itic + 1
               insp(itic) = i 
               jnsp(itic) = i 
               iy2p(i) = itic
            end if 
         end do 

         do i = 1, mdep 
            insp(itic+i) = jdep(i)
         end do 

         do i = 1, norder 
            insp(istot+i) = istot + i
            jnsp(itic+i) = istot + i 
            iy2p(istot+i) = itic + i 
         end do 

      else 

         do i = 1, istot + norder
            insp(i) = i 
            jnsp(i) = i 
            iy2p(i) = i 
         end do 

      end if 

      if (depend) then 
c                                 make y2p array       
         do i = 1, kstot + norder 

            do j = 1, mdep

               y2p(i,j) = 0d0

               do k = 1, ndph(j)
                  if (idep(j,k).eq.jnsp(i)) y2p(i,j) = nu(j,k)
               end do 

            end do
         end do

      end if

      end 

      subroutine gphase (id,gph)
c-----------------------------------------------------------------------
c gphase computes the gibbs free energy of a compound identified by index id.
c gphase does not assume that the properties of a pseudocompound endmembers
c are known (unlike gall) it is thus less efficient than gall.
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer ipoint,imyn,k,id,ids

      double precision gph,gzero,dg

      common/ cst60  /ipoint,imyn

      integer ikp
      common/ cst61 /ikp(k1)

      integer ixp,ifp
      double precision sxs,exces
      common/ cst304 /sxs(k13),exces(m3,k1),ixp(k1),ifp(k1)

      integer jend
      common/ cxt23 /jend(h9,k12)
c                                 working arrays
      double precision z, pa, p0a, x, w, y
      common/ cxt7 /x(m4),y(m4),pa(m4),p0a(m4),z(mst,msp),w(m1)
c                                 new global arrays, 10/25/05:
c                                 bookkeeping variables
      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)
c                                 model type
      logical lorder, lexces, llaar, lrecip
      common/ cxt27 /lorder(h9),lexces(h9),llaar(h9),lrecip(h9)
c----------------------------------------------------------------------
      ids = ikp(id)
 
      if (id.le.ipoint) then
c                                 phase is an endmember compound
         call gcpd (id,gph)

      else if (lorder(ids).and.lrecip(ids)) then 
c                                 reciprocal solution speciation model 
c                                 with nord order parameters
c                                 load x's from sxs array.                              
         do k = 1, nstot(ids) 
            p0a(k) = sxs(ixp(id)+k)
            pa(k) = p0a(k)
         end do 
c                                 compute margules coefficients
         call setw (ids)
c                                 get the speciation energy effect
         call specis (gph,ids)
c                                 get endmember dqf's
         call gexces (id,dg) 

         gph = gph + dg 
c                                 get gmech
         do k = 1, lstot(ids)
            call gcpd (jend(ids,2+k),dg)
            gph = gph + dg * p0a(k) 
         end do 

      else if (lorder(ids)) then
c                                 initialize gph and get any dqf corrections
         call gexces (id,gph)
c                                 get gmech
         do k = 1, lstot(ids)
            p0a(k) = sxs(ixp(id)+k)
            pa(k) = p0a(k)
            call gcpd (jend(ids,2+k),dg)
            gph = gph + dg*p0a(k)
         end do 
c                                 compute margules coefficients
         call setw (ids) 
c                                 get the speciation energy effect
         call specis (dg,ids)
c                                 add in entropy effect pseudocompound version
         gph = gph + dg

      else 
c                              a pseudocompound without speciation:
         if (ifp(id).ne.1) then
c                              get the excess and/or ideal mixing effect
c                              and/or dqf corrections:
            if (ifp(id).eq.23) then 
               call toop(id,gph)
            else if (ifp(id).eq.26) then 
               call hcneos (gph,sxs(ixp(id)+1),
     *                      sxs(ixp(id)+2),sxs(ixp(id)+3))
            else
               call gexces (id,gph)
            end if 
c                              compute mech mix G for 
c                              all models except fluid 
            do k = 1, nstot(ids) 
               call gcpd (jend(ids,2+k),dg)
               gph = gph + dg*sxs(ixp(id)+k)
            end do 

         else 
c                              get the excess and/or ideal mixing effect
c                              and/or dqf corrections:
            call fexces (id,gph)
c                              excess props don't include vdp:
            do k = 1, nstot(ids) 
               gph = gph + gzero(jend(ids,2+k))*sxs(ixp(id)+k)
            end do 
         end if 
c                              for van laar get fancier excess function         
         if (llaar(ids)) then
            call setw(ids) 
            call gvlaar (ikp(id),id,gph)
         end if 
      end if 

      end

      double precision function dgdp (id)
c----------------------------------------------------------------------
c subroutine to compute the G derivative of solution id with respect to 
c the concentration of the 1st ordered species. assumes the excess function 
c is not higher than second order.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,i1,i2,id

      double precision gex,dgex,dsconf,tphi,dtphi

      external dsconf

      double precision r,v,tr,pr,ps
      common/ cst5   /v(l2),tr,pr,r,ps
c                                 working arrays
      double precision z, pa, p0a, x, w, y
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1)

      double precision alpha,dt
      common/ cyt0  /alpha(m4),dt(j3)
c                                 excess energy variables
      integer jterm, jord, jsub
      common/ cxt2i /jterm(h9),jord(h9),jsub(m2,m1,h9)

      double precision dvnu,deph,dydy
      common/ cxt3r /dvnu(m4,j3,h9),deph(j3,h9),dydy(m4,j3,h9)
c                                 model type
      logical lorder, lexces, llaar, lrecip
      common/ cxt27 /lorder(h9),lexces(h9),llaar(h9),lrecip(h9)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)
c----------------------------------------------------------------------
c                                 get derivative of excess function
      dgex = 0d0

      if (lexces(id)) then 

         if (.not.llaar(id)) then 
c                                 regular excess function
            do i = 1, jterm(id)
               i1 = jsub(1,i,id)
               i2 = jsub(2,i,id)
               dgex = dgex + w(i)*( pa(i1)*dydy(i2,1,id) 
     *                            + pa(i2)*dydy(i1,1,id))
            end do
         else 
c                                 h&p van laar
            tphi = 0d0
            dtphi = 0d0

            do i = 1, nstot(id)
               tphi = tphi + alpha(i)* y(i) 
               dtphi = dtphi + alpha(i)*dydy(i,1,id)
            end do 

            gex = 0d0 

            do i = 1, jterm(id)
c                                 assume holland powell form, all terms regular
              i1 = jsub(1,i,id)
              i2 = jsub(2,i,id)

              gex = gex + w(i) * y(i1) * y(i2)   
              dgex = dgex + w(i) * (y(i1)*dydy(i2,1,id) 
     *                            + y(i2)*dydy(i1,1,id))

            end do 
c                                note the excess energy is gex/tphi
            dgex = (dgex - dtphi*gex/tphi)/tphi

         end if 

      end if 
c                                 now get dg/dy(jd)
      dgdp = deph(1,id) + dgex - v(2)*dsconf(id)

      end

      double precision function dsconf (id)
c----------------------------------------------------------------------
c subroutine to the derivative of the configurational entropy of a 
c solution with respect to the proportion of a dependent species.

c THIS DOES NOT INCLUDE ENDMEMBER CONFIGURATION ENTROPY DERIVATIVES!
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,j,k,id

      double precision dlnz,dscon,zt,q,dzdy,z,dd
c                                 working arrays
      double precision zz, pa, p0a, x, w, y
      common/ cxt7 /zz(m4),y(m4),pa(m4),p0a(m4),x(mst,msp),w(m1)
c                                 configurational entropy variables:
      integer msite, ksp, lterm, ksub
      common/ cxt1i /msite(h9),ksp(m10,h9),lterm(m11,m10,h9),
     *               ksub(m0,m11,m10,h9)

      double precision qmult, d0, dcoef, scoef      
      common/ cxt1r /qmult(m10,h9),d0(m11,m10,h9),dcoef(m0,m11,m10,h9),
     *               scoef(m4,h9)

      double precision dppp,d2gx,sdzdp
      common/ cxt28 /dppp(j3,j3,m1,h9),d2gx(j3,j3),sdzdp(j3,m11,m10,h9)
c----------------------------------------------------------------------
      dscon = 0d0
c                                 for each site
      do i = 1, msite(id)

         zt = 0d0
         dd = 0d0 
         q = qmult(i,id)
c                                 get site fractions
         do j = 1, ksp(i,id)

            z = d0(j,i,id)
c                                 for each term:
            do k = 1, lterm(j,i,id)
               z = z + dcoef(k,j,i,id) * pa(ksub(k,j,i,id))
            end do 

            zt = zt + z
c                                 sdzdp is (dz(i,j)/dp(k))
            dzdy = sdzdp(1,j,i,id)

            if (dzdy.eq.0d0) cycle

            if (z.gt.0d0) then 
               dlnz = 1d0 + dlog(z)
            else 
c                                 the term should goto -Inf
               dlnz = 1d0
            end if 

            dscon = dscon - q * dzdy * dlnz

            dd = dd + dzdy

         end do 

         z = 1d0 - zt 

         if (z.gt.0d0) then 
            dlnz = 1d0 + dlog(z)
         else 
            dlnz = 1d0
         end if 

         dscon = dscon - q * sdzdp(1,j,i,id) * dlnz

      end do 

      dsconf = dscon

      end

      subroutine gvlaar (jd,id,dg)
c-----------------------------------------------------------------------
c gvlaar evaluates the contribution to the gibbs energy of a pseudocompound
c arising from van laar excess properties. 
c
c now redundant with function gex except gvlaar uses sxs instead of y
c
c input:
c      jd - solution pointer
c      id - compound pointer
c      dg - ideal free energy
c output:
c      dg - ideal free energy + excess energy     
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,id,jd,i1,i2

      double precision phi(m4), tphi, dg

      integer ixp,ifp
      double precision sxs,exces
      common/ cst304 /sxs(k13),exces(m3,k1),ixp(k1),ifp(k1)
c                                 excess energy variables
      integer jterm, jord, jsub
      common/ cxt2i /jterm(h9),jord(h9),jsub(m2,m1,h9)

      double precision z, pa, p0a, x, w, y
      common/ cxt7 /y(m4),x(m4),pa(m4),p0a(m4),z(mst,msp),w(m1)

      double precision alpha,dt
      common/ cyt0  /alpha(m4),dt(j3)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)
c-----------------------------------------------------------------------
c                                 first compute "volumes"
      tphi = 0d0

      do i = 1, nstot(jd)


c                                 these phi's are the numerator of
c                                 holland & powells "phi's"
         phi(i) =  alpha(i)* sxs(ixp(id)+i)
c                                 tphi is the denominator of holland & powell's
c                                 phi's
         tphi = tphi + phi(i) 
      end do 
c                                 rearrange tphi
      tphi = 2d0/tphi
c                                 dg is initialized as gph in the calling 
c                                 program
      do i = 1, jterm(jd)
c                                 assume holland powell form, all terms regular
         i1 = jsub(1,i,jd)
         i2 = jsub(2,i,jd)
         dg = dg + w(i)  
     *           * tphi * phi(i1) * phi(i2) / (alpha(i1) + alpha(i2))   

      end do  

      end 


      double precision function hpmelt (im)
c----------------------------------------------------------------------
c subroutine to evaluate the configurational entropy of Holland & Powell's
c haplogranite melt model, dlnw is S.

c modified to use global arrays, 10/26/05.
c----------------------------------------------------------------------

      implicit none

      include 'perplex_parameters.h'

      integer i, im

      double precision dlnw, ytot, yh2o, yfo, yfa
c                                 global arrays:
      double precision t, p, xco2, u1, u2, tr, pr, r, ps

      common/ cst5   /p,t,xco2,u1,u2,tr,pr,r,ps
c                                 bookkeeping variables
      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)
c                                 special model endmember indexing
      integer ispec

      common/ cxt8 /ispec(h9,m4)
c                                 working arrays
      double precision z, pa, p0a, x, w, y
      common/ cxt7 /y(m4),x(m4),pa(m4),p0a(m4),z(mst,msp),w(m1)
c----------------------------------------------------------------------  
      dlnw = 0d0

      if (ispec(im,1).eq.1) then

         yh2o = y(1)
c                                 the quadratic water term
         if (yh2o.ne.0d0) dlnw = -2d0 * yh2o * dlog(yh2o)

      else 
 
         yh2o = 0d0 

      end if 
c                                 the basic entropy
      do i = ispec(im,4), mstot(im)

         if (y(i).le.0d0) cycle
         dlnw = dlnw - y(i) * dlog(y(i)*(1d0-yh2o))
 
      end do 
c                                 the fe-mg fudge factor
      yfo = 0d0
      yfa = 0d0
      if (ispec(im,2).ne.0) yfo = y(ispec(im,2))   
      if (ispec(im,3).ne.0) yfa = y(ispec(im,3))
      ytot = yfo + yfa

      if (ytot.ne.0d0) then 
         if (yfo.ne.0d0) dlnw = dlnw - 4d0 * yfo * dlog(yfo/ytot)
         if (yfa.ne.0d0) dlnw = dlnw - 4d0 * yfa * dlog(yfa/ytot)
      end if 

      hpmelt = dlnw*r

      end 

      double precision function gmelt (im)
c----------------------------------------------------------------------
c subroutine to evaluate the configurational entropy of Ghiorso's
c pMELTS/MELTS model, dlnw is S.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, im
      double precision dlnw, yh2o
c                                 global arrays:
      double precision t, p, xco2, u1, u2, tr, pr, r, ps
      common/ cst5   /p,t,xco2,u1,u2,tr,pr,r,ps
c                                 bookkeeping variables
      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)
c                                 special model endmember indexing
      integer ispec
      common/ cxt8 /ispec(h9,m4)
c                                 working arrays
      double precision z, pa, p0a, x, w, y
      common/ cxt7 /y(m4),x(m4),pa(m4),p0a(m4),z(mst,msp),w(m1)
c----------------------------------------------------------------------
      dlnw = 0d0

      if (ispec(im,1).eq.1) then

         yh2o = y(1)
c                                 the quadratic water term, cf hp model. 
         if (yh2o.ne.0d0) dlnw = -(yh2o * dlog(yh2o) 
     *                             + (1d0-yh2o)*dlog(1d0-yh2o))
      else 

         yh2o = 0d0  
       
      end if 
c                                 i wouldn't expect to count water here yet again, but that
c                                 is what eqs 1 & 2 of Ghiorso et al. (2002) imply

c                                 in any event equation 1 in Ghiorso et al. 2002
c                                 is almost certainly wrong, so the entropy here
c                                 is computed as in ghiorso & sack 1995 (melts).
c                                 this correction made Nov 4, 2005. This model seems
c                                 to be consistent with Nicholls 1980 CMP 74:211 

c                                 the basic entropy
      do i = 1, mstot(im)

         if (y(i).le.0d0) cycle
         dlnw = dlnw - y(i) * dlog(y(i)*(1d0-yh2o))
 
      end do 

      gmelt = r*dlnw

      end  

      subroutine ytox (ids)
c----------------------------------------------------------------------
c subroutine to convert endmember fractions (y) to geometric coordinates 
c (x) for solution model ids.
c----------------------------------------------------------------------
      implicit none 

      include 'perplex_parameters.h'

      integer ids, i, j
c                                 convert y -> x array
      integer indx

      common/ cxt5i /indx(h9,mst,msp)
c                                 working arrays
      double precision z, pa, p0a, x, w, y

      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1)

c                                 x coordinate description
      integer istg, ispg, imlt, imdg

      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)
c----------------------------------------------------------------------
      do i = 1, istg(ids)
c                                 this conditional is necessary
c                                 because reform can generate a 
c                                 one species site. 
         if (ispg(ids,i).gt.1) then 
c                                 initialize
            do j = 1, ispg(ids,i)
               x(i,j) = 0d0
            end do 
  
            do j = 1, ispg(ids,i)
               x(i,j) = x(i,j) + y(indx(ids,i,j))
            end do 
         else 
c                                 one species site. 
            x(i,1) = 1d0
         end if 
      end do
     
      end 

      subroutine getolx (ids,id)
c----------------------------------------------------------------------
c subroutine to recover geometric reciprocal solution compositions (x(i,j))
c from the xcoor (reciprocal) or sxs (single site) arrays loaded in soload
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, id, ids, jcoor
c                                 working arrays
      double precision z, pa, p0a, x, w, y

      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1)
c                                 x coordinate description
      integer istg, ispg, imlt, imdg

      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)
c                                 stored x coordinate
      double precision xcoor
      integer icoor
      common/ cxt10 /xcoor(k18),icoor(k1)
c                                 single site solution coordinates:
      integer ixp,ifp
      double precision sxs,exces
      common/ cst304 /sxs(k13),exces(m3,k1),ixp(k1),ifp(k1)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)
c----------------------------------------------------------------------
      if (istg(ids).gt.1) then 

         jcoor = icoor(id)

         do i = 1, istg(ids)
            do j = 1, ispg(ids,i)
               jcoor = jcoor + 1
               x(i,j) = xcoor(jcoor)
            end do 
         end do 

      else 

         do j = 1, nstot(ids)
            x(1,j) = sxs(ixp(id)+j) 
         end do 

      end if 

      end 

      subroutine xtoy (ids)
c----------------------------------------------------------------------
c subroutine to convert geometric reciprocal solution compositions (x(i,j))
c to geometric endmember fractions (y) for solution model ids.
c----------------------------------------------------------------------
      implicit none 

      include 'perplex_parameters.h'
c                                 -------------------------------------
c                                 local variables:
      integer ids, l, m
c                                 -------------------------------------
c                                 global variables:
c                                 bookkeeping variables
      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)
c                                 working arrays
      double precision z, pa, p0a, x, w, y
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1)
c                                 x coordinate description
      integer istg, ispg, imlt, imdg
      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)

      integer ksmod, ksite, kmsol, knsp
      common/ cxt0  /ksmod(h9),ksite(h9),kmsol(h9,m4,mst),knsp(m4,h9)
c----------------------------------------------------------------------

      do l = 1, mstot(ids)

         y(l) = 1d0

         do m = 1, istg(ids)
            y(l) = y(l)*x(m,kmsol(ids,l,m))
         end do

      end do   

      end 

      double precision function gsol (id)
c-----------------------------------------------------------------------
c gsolut computes the total (excess+ideal) free energy of solution 
c for a solution identified by index ids and composition y(m4) input
c from cxt7, the composition y is the independent endmember fractions
c for all model types except reciprocal solutions with ordering, in 
c which case it is the y's for the full reciprocal model.

c gsol assumes the endmember g's have not been calculated by gall,
c gsol1 is identical to gsol but can only been called after gall. 
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer k,id

      double precision omega, gproj, hpmelt, gmelt, gfluid, gzero, g, 
     *                 dg, gex

      external gex

      integer jend
      common/ cxt23 /jend(h9,k12)

      double precision r,tr,pr,ps,p,t,xco2,u1,u2
      common/ cst5   /p,t,xco2,u1,u2,tr,pr,r,ps
c                                 bookkeeping variables
      integer ksmod, ksite, kmsol, knsp
      common/ cxt0  /ksmod(h9),ksite(h9),kmsol(h9,m4,mst),knsp(m4,h9)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)
c                                 working arrays
      double precision z, pa, p0a, x, w, y
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1)
c                                 model type
      logical lorder, lexces, llaar, lrecip
      common/ cxt27 /lorder(h9),lexces(h9),llaar(h9),lrecip(h9)

      integer ispec
      common/ cxt8 /ispec(h9,m4)
c----------------------------------------------------------------------
      if (id.lt.0) then 

         call gphase (-id,g)
         gsol = g

      else 

         g = 0d0
c                                 evaluate dqf coefficients
         call setdqf (id)

         call setw (id) 

         if (ksmod(id).eq.1.or.ksmod(id).eq.2.or.ksmod(id).eq.3) then 
c                                 -------------------------------------
c                                 macroscopic formulation for normal solutions.
            call gdqf (id,g,y) 
c                                 add entropy and excess contributions
            g = g - t * omega(id,y) + gex(id,y)
c                                 get mechanical mixture contribution
            do k = 1, mstot(id) 
               g = g + y(k) * gproj (jend(id,2+k))
            end do 

         else if (lrecip(id).and.lorder(id)) then 
c                                 -------------------------------------
c                                 convert y coordinates to independent p coordinates
            call y2p0 (id)
c                                 get the speciation, excess and entropy effects.
            call specis (g,id)

            do k = 1, lstot(id) 
c                                 compute mechanical g from these z's, 
c                                 specip adds a correction for the ordered species.
               g = g + gproj(jend(id,2+k)) * p0a(k)
            end do 
c                                 get the dqf, this assumes the independent reactants
c                                 are not dqf'd. gex not neccessary as computed in specip
            call gdqf (id,g,p0a)

         else if (lorder(id)) then 
c                                 -------------------------------------
c                                 non-reciprocal speciation.
            do k = 1, lstot(id)  
               pa(k) = y(k)
               p0a(k) = y(k)
               g = g + y(k) * gproj (jend(id,2+k))
            end do 
c                                 get the speciation energy effect
            call specis (dg,id)

            g = g + dg 
c                                 get dqf corrections
            call gdqf (id,g,p0a) 
 
         else if (lrecip(id)) then 
c                                 -------------------------------------
c                                 macroscopic reciprocal solution 

c                                 initialize
c                                 convert y's to p's (p0a here), this loop 
c                                 could be eliminated by saving the p's. the reason
c                                 for naming the p's z here is so that they can be
c                                 subsequently treated the same way as is the case
c                                 for order-disorder reciprocal solution (ksmod=8)
            call y2p0 (id)

            do k = 1, lstot(id)
               g = g + gproj (jend(id,2+k)) * p0a(k) 
            end do 
c                                 get the dqf
            call gdqf (id,g,p0a)
c                                 and excess contributions
            g = g - t * omega(id,p0a) + gex(id,p0a)

         else if (ksmod(id).eq.23) then 

             write (*,*) 'toop samis model not coded'

         else if (ksmod(id).eq.24) then 
c                                 -------------------------------------
c                                 hp melt model         
            call gdqf (id,g,y) 
            g = g - t * hpmelt(id) + gex(id,y)
c                                 get mechanical mixture contribution
            do k = 1, mstot(id)  
               g = g + y(k) * gproj (jend(id,2+k))
            end do 

         else if (ksmod(id).eq.25) then 
c                                 -------------------------------------
c                                 ghiorso pmelt model  
            call gdqf (id,g,y) 

            g = g - t * gmelt(id) + gex(id,y)
c                                 get mechanical mixture contribution
            do k = 1, mstot(id)  
               g = g + y(k) * gproj (jend(id,2+k))
            end do 

         else if (ksmod(id).eq.26) then 
c                                 ------------------------------------
c                                 andreas salt model
            call hcneos (g,y(1),y(2),y(3))

            do k = 1, 3
               g = g + y(k) * gproj (jend(id,2+k))
            end do 

         else if (ksmod(id).eq.0) then 
c                                 ------------------------------------
c                                 internal fluid eos
            do k = 1, 2
               g = g + gzero(jend(id,2+k))*y(k)
            end do 

            g = g + gfluid(y(ispec(id,1)))

         else 

            write (*,*) 'what the **** am i doing here?'
            stop

         end if 

      gsol = g 

      end if 

      end

      double precision function gsol1 (id)
c-----------------------------------------------------------------------
c gsolut computes the total (excess+ideal) free energy of solution 
c for a solution identified by index ids and composition y(m4) input
c from cxt7, the composition y is the independent endmember fractions
c for all model types except reciprocal solutions with ordering, in 
c which case it is the y's for the full reciprocal model.

c gsol1 assumes the endmember g's have been calculated by gall. 
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer k,id

      double precision omega, hpmelt, gmelt, gfluid, gzero, gg, dg, gex

      external gex

      integer jend
      common/ cxt23 /jend(h9,k12)

      double precision g
      common/ cst2 /g(k1)

      double precision r,tr,pr,ps,p,t,xco2,u1,u2
      common/ cst5   /p,t,xco2,u1,u2,tr,pr,r,ps
c                                 bookkeeping variables
      integer ksmod, ksite, kmsol, knsp
      common/ cxt0  /ksmod(h9),ksite(h9),kmsol(h9,m4,mst),knsp(m4,h9)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)
c                                 working arrays
      double precision z, pa, p0a, x, w, y
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1)
c                                 model type
      logical lorder, lexces, llaar, lrecip
      common/ cxt27 /lorder(h9),lexces(h9),llaar(h9),lrecip(h9)

      integer ispec
      common/ cxt8 /ispec(h9,m4)
c----------------------------------------------------------------------
         gg = 0d0
c                                 evaluate margules coefficients
         call setw (id)
c                                 evaluate dqf coefficients
         call setdqf(id)

         if (ksmod(id).eq.1.or.ksmod(id).eq.2.or.ksmod(id).eq.3) then 
c                                 -------------------------------------
c                                 macroscopic formulation for normal solutions.
            call gdqf (id,gg,y) 

            gg = gg - t * omega(id,y) + gex(id,y)
c                                 get mechanical mixture contribution
            do k = 1, mstot(id) 
               gg = gg + y(k) * g(jend(id,2+k))
            end do 

         else if (lrecip(id).and.lorder(id)) then 
c                                 -------------------------------------
c                                 initialize p's
            call y2p0 (id)
c                                 get the speciation, excess and entropy effects.
            call specis (gg,id)
c                                 decompose the ordered species into 
c                                 the independent disordered species
c                                 i.e., the p0a array becomes the p's if the 
c                                 abundance of the ordered species is 0.
            do k = 1, lstot(id)
c                                 compute mechanical g from these z's, 
c                                 specip adds a correction for the ordered species.
               gg = gg + g(jend(id,2+k)) * p0a(k)
            end do 
c                                 get the dqf, this assumes the independent reactants
c                                 are not dqf'd. gex not neccessary as computed in specip
            call gdqf (id,gg,p0a)

         else if (lorder(id)) then 
c                                 -------------------------------------
c                                 macroscopic formulation for ordering solutions.

c                                 get mechanical mixture contribution
            do k = 1, lstot(id)  
               pa(k) = y(k)
               p0a(k) = y(k)
               gg = gg + y(k) * g(jend(id,2+k))
            end do 
c                                 get the speciation energy effect
            call specis (dg,id)

            gg = gg + dg 
c                                 get dqf corrections
            call gdqf (id,gg,p0a) 

         else if (lrecip(id)) then 
c                                 -------------------------------------
c                                 macroscopic reciprocal solution
c                                 initialize p's
            call y2p0 (id)

            do k = 1, lstot(id)
               gg = gg + g(jend(id,2+k)) * p0a(k) 
            end do 
c                                 get the dqf, this assumes the independent reactants
c                                 are not dqf'd
            call gdqf (id,gg,p0a)

            gg = gg - t * omega(id,p0a) + gex(id,p0a)

         else if (ksmod(id).eq.23) then 

             write (*,*) 'toop samis model not coded'

         else if (ksmod(id).eq.24) then 
c                                 -------------------------------------
c                                 hp melt model         
            call gdqf (id,gg,y) 

            gg = gg - t * hpmelt(id) + gex(id,y)
c                                 get mechanical mixture contribution
            do k = 1, mstot(id)  
               gg = gg + y(k) * g(jend(id,2+k))
            end do 

         else if (ksmod(id).eq.25) then 
c                                 -------------------------------------
c                                 ghiorso pmelt model 
            call gdqf (id,gg,y) 

            gg = gg - t * gmelt(id) + gex(id,y)
c                                 get mechanical mixture contribution
            do k = 1, mstot(id)  
               gg = gg + y(k) * g(jend(id,2+k))
            end do 

         else if (ksmod(id).eq.26) then 
c                                 ------------------------------------
c                                 andreas salt model
            call hcneos (gg,y(1),y(2),y(3))

            do k = 1, 3
               gg = gg + y(k) * g(jend(id,2+k))
            end do 

         else if (ksmod(id).eq.0) then 
c                                 ------------------------------------
c                                 internal fluid eos
            gg = gfluid(y(ispec(id,1)))
            
            do k = 1, 2
               gg = gg + gzero(jend(id,2+k))*y(k)
            end do 

         else 

            write (*,*) 'what the **** am i doing here?'
            stop

         end if 

      gsol1 = gg 

      end

      subroutine zchk0 (x,tname)
c----------------------------------------------------------------------
c subroutine to site fractions computed from equations entered by 
c user for configurational entropy (bragg williams form).
c----------------------------------------------------------------------

      implicit none

      include 'perplex_parameters.h'

      character*10 tname

      double precision x(mst,msp), z, zt, zad

      integer i,j,k,l

      integer nsub1,nsub2,nttyp,nterm,nspm1,nsite
      double precision acoef,smult,a0
      common/ cst107 /a0(m10,m11),acoef(m10,m11,m0),smult(m10),
     *      nsite,nspm1(m10),nterm(m10,m11),nsub1(m10,m11,m0,m12),
     *      nsub2(m10,m11,m0,m12),nttyp(m10,m11,m0)
c----------------------------------------------------------------------
c                                 for each site
      do i = 1, nsite
         zt = 0d0
c                                 get site fractions
         do j = 1, nspm1(i)
            z = a0(i,j)
            if (nterm(i,j).ne.0) then 
c                                 for each term:
               do k = 1, nterm(i,j)
                  zad = acoef(i,j,k)
                  do l = 1, nttyp(i,j,k)
                     zad = zad * x(nsub1(i,j,k,l),nsub2(i,j,k,l))
                  end do 
                  z = z + zad
               end do 
            end if 

            call zchk (z,tname)
            zt = zt + z
            
         end do 
         z = 1d0 - zt
         call zchk (z,tname)
      end do 

      end

      subroutine zchk1 (x,ids,tname)
c----------------------------------------------------------------------
c subroutine to site fractions computed from equations entered by 
c user for configurational entropy (macroscopic form).
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      character*10 tname

      double precision x(m4),z,zt

      integer i,j,k,ids
c                                 configurational entropy variables:
      integer msite, ksp, lterm, ksub
      common/ cxt1i /msite(h9),ksp(m10,h9),lterm(m11,m10,h9),
     *               ksub(m0,m11,m10,h9)

      double precision qmult, d0, dcoef, scoef      
      common/ cxt1r /qmult(m10,h9),d0(m11,m10,h9),dcoef(m0,m11,m10,h9),
     *               scoef(m4,h9)
c----------------------------------------------------------------------
c                                 for each site
      do i = 1, msite(ids)

         zt = 0d0
c                                 get site fractions
         do j = 1, ksp(i,ids)

            z = d0(j,i,ids)
c                                 for each term:
            do k = 1, lterm(j,i,ids)
               z = z + dcoef(k,j,i,ids) * x(ksub(k,j,i,ids))
            end do 

            call zchk (z,tname)
            zt = zt + z
            
         end do
 
         z = 1d0 - zt

         call zchk (z,tname)

      end do 

      end

      double precision function omega (id,y)
c----------------------------------------------------------------------
c subroutine to evaluate the configurational entropy of a solution
c with composition y, including the correction for endmember 
c configurational negentropies. reciprocal end-member composition version.

c note this version doesn't evaluate ordered endmember S, this is fine
c so long as the ordered and its stoichiometrically equivalent disordered 
c endmembers have the same S (as is currently the case for HP models, 
c 10/20/05).
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision z,zt,dlnw,dlnz,y(m4)

      integer i,j,k,id
c                                 configurational entropy variables:
      integer msite, ksp, lterm, ksub
      common/ cxt1i /msite(h9),ksp(m10,h9),lterm(m11,m10,h9),
     *               ksub(m0,m11,m10,h9)

      double precision qmult, d0, dcoef, scoef      
      common/ cxt1r /qmult(m10,h9),d0(m11,m10,h9),dcoef(m0,m11,m10,h9),
     *               scoef(m4,h9)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)

      double precision v,tr,pr,r,ps
      common / cst5 /v(l2),tr,pr,r,ps
c----------------------------------------------------------------------
      dlnw = 0d0
c                                 for each site
      do i = 1, msite(id)

         zt = 0d0
         dlnz = zt
c                                 get site fractions
         do j = 1, ksp(i,id)
            z = d0(j,i,id) 
c                                 for each term:
            do k = 1, lterm(j,i,id)
               z = z + dcoef(k,j,i,id) * y(ksub(k,j,i,id))
            end do 

            zt = zt + z
            if (z.gt.0d0) dlnz = dlnz - z * dlog (z)
            
         end do 

         z = 1d0 - zt
         if (z.gt.0d0) dlnz = dlnz - z * dlog(z)
         dlnw = dlnw + qmult(i,id)*dlnz 

      end do 
c                                 endmember corrections
      do i = 1, nstot(id)
         dlnw = dlnw - y(i)*scoef(i,id)
      end do 

      omega = dlnw 

      end

      subroutine snorm (id,tname)
c------------------------------------------------------------------------
c compute endmember configurational entropies
c site fractions expressed as a function of end-member proportions.
c see deleted snorm0 for bragg-williams site fraction version.
c------------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      character*10 tname

      integer h,id,j

      double precision omega

      integer iorig,jnsp,iy2p
      common / cst159 /iorig(m4),jnsp(m4),iy2p(m4)
c                                 working arrays
      double precision z, pa, p0a, x, w, y
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1)
c                                 bookkeeping variables
      integer ksmod, ksite, kmsol, knsp
      common/ cxt0  /ksmod(h9),ksite(h9),kmsol(h9,m4,mst),knsp(m4,h9)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)

      double precision qmult, d0, dcoef, scoef      
      common/ cxt1r /qmult(m10,h9),d0(m11,m10,h9),dcoef(m0,m11,m10,h9),
     *               scoef(m4,h9)
c----------------------------------------------------------------------
c                                 get normalization constants
c                                 for each endmember
      do h = 1, nstot(id) 
c                                 zero pa-array
         do j = 1, nstot(id) 
            pa(j) = 0d0 
         end do

         pa(h) = 1d0        
c                                 check pa's now set
         call zchk1 (pa,id,tname)
c                                 evaluate S
         scoef(h,id) = omega(id,pa)

      end do 

      end 

      subroutine setdqf (id)
c---------------------------------------------------------------------
c setdqf - evalutaes dqf coefficients at p and t
c---------------------------------------------------------------------
      implicit none
  
      include 'perplex_parameters.h'

      integer i, id

      double precision t, p, xco2, u1, u2, tr, pr, r, ps
      common/ cst5   /p,t,xco2,u1,u2,tr,pr,r,ps
c                                 dqf corrections
      integer jndq, jdqf, iq
      double precision dqfg, dq 
      common/ cxt9 /dqfg(m3,m4,h9),dq(m4),jndq(m4,h9),jdqf(h9),iq(m4)
c----------------------------------------------------------------------

      do i = 1, jdqf(id)
c                                 index points to the endmember in the full
c                                 model:
         iq(i) = jndq(i,id)

         dq(i) = dqfg(1,i,id) + t*dqfg(2,i,id) + p*dqfg(3,i,id)
         
      end do 

      end 

      subroutine gdqf (id,g,y)
c----------------------------------------------------------------------
c subroutine to evaluate the excess G of solution id with macroscopic
c composition y. assumes a previous call to setdqf. When setdqf is 
c not called use 
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, id
c                                 working arrays
      double precision g, y(m4)
c                                 dqf corrections
      integer jndq, jdqf, iq

      double precision dqfg, dq 

      common/ cxt9 /dqfg(m3,m4,h9),dq(m4),jndq(m4,h9),jdqf(h9),iq(m4)
c----------------------------------------------------------------------
      do i = 1, jdqf(id)

         g = g + y(iq(i))*dq(i)

      end do  

      end 

      double precision function gproj (id)
c-----------------------------------------------------------------------
c gproj computes the projected molar free energy of the
c phase with index id.
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer iff,idss,ifug,ifyn,isyn,icomp,istct,iphct,icp,j,id

      double precision thermo,uf,us,cp,gph

      common/ cst1 /thermo(k4,k10),uf(2),us(h5)/ cst12 /cp(k5,k1)
     *      / cst10 /iff(2),idss(h5),ifug,ifyn,isyn
     *      / cst6  /icomp,istct,iphct,icp

      integer ifct,idfl
      common/ cst208 /ifct,idfl

      integer ids,isct,icp1,isat,io2
      common/ cst40 /ids(h5,h6),isct(h5),icp1,isat,io2

c-----------------------------------------------------------------------
      call gphase (id,gph)
c                                 this is a screw up solution
c                                 necessary cause uf(1) and uf(2)
c                                 are allocated independent of ifct!
      if (ifct.gt.0) then 
         do j = 1, 2
            if (iff(j).ne.0) gph = gph - cp(iff(j),id)*uf(j)
          end do 
      end if 

      do j = 1, isat
         gph = gph - cp(icp+j,id) * us(j)
      end do 

      gproj = gph

      end

      subroutine setw (id)
c---------------------------------------------------------------------
c setw - evaluates margules coeffiecients and, for laar models, size
c parameters for solution model id at p and t
c---------------------------------------------------------------------
      implicit none
  
      include 'perplex_parameters.h'

      integer i, k, l, i1, i2, id

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps
c                                 excess energy variables
      integer jterm, jord, jsub
      common/ cxt2i /jterm(h9),jord(h9),jsub(m2,m1,h9)

      double precision wgl,vlar
      common/ cxt2r /wgl(m3,m1,h9),vlar(m3,m4,h9)
c                                 working arrays
      double precision z, pa, p0a, x, w, y
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1)
c                                 local alpha
      double precision alpha,dt
      common/ cyt0  /alpha(m4),dt(j3)

      double precision dvnu,deph,dydy
      common/ cxt3r /dvnu(m4,j3,h9),deph(j3,h9),dydy(m4,j3,h9)
c                                 model type
      logical lorder, lexces, llaar, lrecip
      common/ cxt27 /lorder(h9),lexces(h9),llaar(h9),lrecip(h9)
c                                 bookkeeping variables
      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)

      double precision dppp,d2gx,sdzdp
      common/ cxt28 /dppp(j3,j3,m1,h9),d2gx(j3,j3),sdzdp(j3,m11,m10,h9)
c----------------------------------------------------------------------
      do i = 1, jterm(id)
         w(i) = wgl(1,i,id) + t*wgl(2,i,id) + p*wgl(3,i,id)
      end do

      if (llaar(id)) then 

         do i = 1, nstot(id) 

            alpha(i) = vlar(1,i,id) 
     *               + t * vlar(2,i,id) + p * vlar(3,i,id)

         end do

         do i = 1, jterm(id)
            i1 = jsub(1,i,id)
            i2 = jsub(2,i,id)
            w(i) = 2d0 * w(i)
     *                 * alpha(i1)*alpha(i2) / (alpha(i1) + alpha(i2))
         end do

      end if 

      if (lorder(id)) then 
c                                 set higher order derivatives for 
c                                 speciation
         do k = 1, nord(id) 
            dt(k) = 0d0
            do l = 1, nord(id)
               d2gx(l,k) = 0d0
            end do  
         end do
c                                 for both laar and regular need
c                                 the d(gex)/dy(k)/dy(l)
         do i = 1, jterm(id)

            i1 = jsub(1,i,id)
            i2 = jsub(2,i,id)

            do k = 1, nord(id) 
               do l = 1, nord(id)
                  d2gx(l,k) = d2gx(l,k) + w(i) * dppp(l,k,i,id)
               end do  
            end do 
         end do


         if (llaar(id)) then 
c                                 for laar also need:
            do i = 1, nstot(id)
               do k = 1, nord(id)  
c                                 dt, derivative of sum(phi)
                  dt(k) = dt(k) + alpha(i)*dydy(i,k,id)
               end do 
            end do 

         end if 

      end if 

      end 

      double precision function gex (ids,y)
c-----------------------------------------------------------------------
c evaluate the excess function for solution model ids. assuming no prior
c call to set coefficients (as in function gexces). 
c input:
c      ids - solution pointer
c      y - composition array

c------------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,j,ids

      double precision y(m4), tphi, xpr

      double precision z, pa, p0a, x, w, yy
      common/ cxt7 /yy(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1)
c                                 bookkeeping variables
      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)
c                                 excess energy variables
      integer jterm, jord, jsub
      common/ cxt2i /jterm(h9),jord(h9),jsub(m2,m1,h9)
c                                 model type
      logical lorder, lexces, llaar, lrecip
      common/ cxt27 /lorder(h9),lexces(h9),llaar(h9),lrecip(h9)
c                                 local alpha
      double precision alpha,dt
      common/ cyt0  /alpha(m4),dt(j3)
c----------------------------------------------------------------------
      gex = 0d0 
 
      if (lexces(ids)) then 

         if (llaar(ids)) then 
c                                 holland & powells version of the van laar
c                                 first compute "volumes"
            tphi = 0d0

            do i = 1, nstot(ids)
c                                 tphi is the sum of holland & powell's
c                                 phi's
               tphi = tphi + alpha(i) * y(i)

             end do 
c                                 dg is initialized as gph in the calling 
c                                 program
            do i = 1, jterm(ids)
c                                 holland powell form, all terms regular
               gex = gex + w(i) * y(jsub(1,i,ids)) * y(jsub(2,i,ids))  

            end do 

            gex = gex/tphi 

         else 
c                                 macroscopic margules formulation by default
            do i = 1, jterm(ids)

               xpr = 1d0

               do j = 1, jord(ids)
                  if (jsub(j,i,ids).eq.0d0) exit
                  xpr = xpr * y(jsub(j,i,ids))
               end do 

               gex = gex + xpr * w(i) 

            end do  

         end if 
      end if 

      end 

      subroutine endcp (jd,id,ids)
c------------------------------------------------------------------------
c compute the composition of endmember id, for solution ids and load it
c into the jd'th position of the x3 array.
c------------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer jd, h, i, j, id, ids
c                                 bookkeeping variables
      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)
c                                 x coordinate description
      integer istg, ispg, imlt, imdg
      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)

      integer jend
      common/ cxt23 /jend(h9,k12)

      double precision x3
      common/ cxt16 /x3(k21,mst,msp)

      integer ksmod, ksite, kmsol, knsp
      common/ cxt0  /ksmod(h9),ksite(h9),kmsol(h9,m4,mst),knsp(m4,h9)
c----------------------------------------------------------------------
c                                 figure out which endmember we
c                                 are looking at:
      do h = 1, mstot(ids)
         if (id.eq.jend(ids,2+h)) exit
      end do          
c                                 zero x-array
      do i = 1, istg(ids)
         do j = 1, ispg(ids,i)
            x3(jd,i,j) = 0d0
         end do
c                                 now assign endmember fractions
         x3(jd,i,kmsol(ids,knsp(h,ids),i)) = 1d0 

      end do   

      end 

      subroutine gmodel (im,tname)
c---------------------------------------------------------------------
c qmodel - stores ALL solution model parameters in global arrays
c---------------------------------------------------------------------
      implicit none
  
      include 'perplex_parameters.h'

      character tname*10, sname*10

      logical add

      integer im,nloc,i,j,ind,jmsol,id,jsmod,k,l,itic,kdsol,ii,imatch
c                                 LOCAL SOLUTION PARAMETERS
      double precision dinc,xsym,dzt,vlaar
                             
      common/ cst221 /vlaar(m3,m4),jsmod
     *      / cst142 /jmsol(m4,mst),kdsol(m4)

      integer nsub1,nsub2,nttyp,nterm,nspm1,nsite
      double precision acoef,smult,a0
      common/ cst107 /a0(m10,m11),acoef(m10,m11,m0),smult(m10),
     *      nsite,nspm1(m10),nterm(m10,m11),nsub1(m10,m11,m0,m12),
     *      nsub2(m10,m11,m0,m12),nttyp(m10,m11,m0)

      integer iend,isub,imd,insp,ist,isp,isite,iterm,iord,istot,jstot,
     *        kstot
      double precision wg,xmn,xmx,xnc
      common/ cst108 /wg(m1,m3),xmn(mst,msp),xmx(mst,msp),xnc(mst,msp),
     *      iend(m4),isub(m1,m2,2),imd(ms1,mst),insp(m4),ist(mst),
     *      isp(mst),isite,iterm,iord,istot,jstot,kstot

      integer iddeps,norder 
      double precision depvnu,denth
      common/ cst141 /depvnu(2,j3),denth(j3),iddeps(2,j3),norder

      integer mdep,idep,jdep,ndph
      double precision nu,y2p
      common/ cst146 /nu(m15,j4),y2p(m4,m15),mdep,jdep(m15),
     *                idep(m15,j4),ndph(m15)

      double precision t, p, xco2, u1, u2, tr, pr, r, ps
      common/ cst5   /p,t,xco2,u1,u2,tr,pr,r,ps

      integer iorig,jnsp,iy2p
      common / cst159 /iorig(m4),jnsp(m4),iy2p(m4)

      logical depend,laar,order,fluid,macro,specil
      common/ cst160 /depend,laar,order,fluid,macro,specil
c                                 GLOBAL SOLUTION PARAMETERS:
      integer ksmod, ksite, kmsol, knsp
      common/ cxt0  /ksmod(h9),ksite(h9),kmsol(h9,m4,mst),knsp(m4,h9)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)
c                                 configurational entropy variables:
      integer msite, ksp, lterm, ksub
      common/ cxt1i /msite(h9),ksp(m10,h9),lterm(m11,m10,h9),
     *               ksub(m0,m11,m10,h9)

      double precision qmult, d0, dcoef, scoef      
      common/ cxt1r /qmult(m10,h9),d0(m11,m10,h9),dcoef(m0,m11,m10,h9),
     *               scoef(m4,h9)
c                                 excess energy variables
      integer jterm, jord, jsub
      common/ cxt2i /jterm(h9),jord(h9),jsub(m2,m1,h9)

      double precision wgl,vlar
      common/ cxt2r /wgl(m3,m1,h9),vlar(m3,m4,h9)

      double precision dppp,d2gx,sdzdp
      common/ cxt28 /dppp(j3,j3,m1,h9),d2gx(j3,j3),sdzdp(j3,m11,m10,h9)

      integer ideps,icase
      common/ cxt3i /ideps(2,j3,h9),icase(h9)

      double precision dvnu,deph,dydy
      common/ cxt3r /dvnu(m4,j3,h9),deph(j3,h9),dydy(m4,j3,h9)
      
      double precision y2pg
      common/ cxt4  /y2pg(m15,m4,h9)
c                                 convert y -> x array
      integer indx
      common/ cxt5i /indx(h9,mst,msp)
c                                 endmember pointers
      integer jend
      common/ cxt23 /jend(h9,k12)
c                                 x coordinate description
      integer istg, ispg, imlt, imdg
      double precision xmng, xmxg, xncg, xmno, xmxo
      common/ cxt6r /xmng(h9,mst,msp),xmxg(h9,mst,msp),xncg(h9,mst,msp),
     *               xmno(h9,mst,msp),xmxo(h9,mst,msp)
      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)
c                                 special model endmember indexing
      integer ispec
      common/ cxt8 /ispec(h9,m4)

      double precision cp
      common/ cst12 /cp(k5,k1)
c                                 dqf parameters
      integer idqf,indq
      double precision dqf
      common/ cst222 /dqf(m3,m4),indq(m4),idqf

      integer jndq, jdqf, iq
      double precision dqfg, dq 
      common/ cxt9 /dqfg(m3,m4,h9),dq(m4),jndq(m4,h9),jdqf(h9),iq(m4)

      integer iemod,kmod
      logical smod
      double precision emod
      common/ cst319 /emod(k15,k10),smod(h9),iemod(k10),kmod

      integer iopt
      double precision nopt
      common / opts /nopt(i10),iopt(i10)

      integer maxit
      common/ cst149 /maxit

      integer ncoor,maxitg
      common/ cxt24 /ncoor(h9),maxitg(h9)

      integer isec,icopt,ifull,imsg,io3p
      common/ cst103 /isec,icopt,ifull,imsg,io3p
c                                 temporary stretching coordinate
c                                 paramters
      double precision yin
      common/ cst50 /yin(ms1,mst)
c                                 parameters for autorefine
      logical stable,limit
      double precision xlo,xhi
      common/ cxt11 /xlo(m4,mst,h9),xhi(m4,mst,h9),stable(h9),limit(h9)

      logical refine
      common/ cxt26 /refine
c                                 interval limits conformal transformation
      integer intv
      double precision yint, yfrc
      common/ cst47 /yint(5,ms1,mst,h9),yfrc(4,ms1,mst,h9),intv(4)
c                                 model type
      logical lorder, lexces, llaar, lrecip
      common/ cxt27 /lorder(h9),lexces(h9),llaar(h9),lrecip(h9)

      integer killct,killid
      common/ cst51/killct,killid(20)

      integer limn,limt,limid,jimid,jimt
      double precision limc,jimc
      common/ cxt30 /limc(j6+2,j5,j3),limid(m0,j5,j3),jimid(j3,j5,j3),
     *               limn(j3),limt(j5,j3),jimc(j3,j5,j3),jimt(j5,j3)

      integer ln,lt,lid,jt,jid
      double precision lc, l0c, jc
      common/ cxt29 /lc(j6,j5,j3,h9),l0c(2,j5,j3,h9),lid(j6,j5,j3,h9),
     *               ln(j3,h9),lt(j5,j3,h9),jc(j3,j5,j3,h9),
     *               jid(j3,j5,j3,h9),jt(j5,j3,h9)
c----------------------------------------------------------------------
c                                 auto_refine changes
      if (refine) then
c                                 check for consistent auto-refine data
         read (n10,'(a)') sname
         if (tname.ne.sname) call error (63,r,i,'GMODEL')

      end if 
c                                 initialize autorefine arrays
      stable(im) = .false.
      limit(im) = .false.

      do i = 1, isite 
         do j = 1, isp(i) - 1

            xlo(j,i,im) = 1d0
            xhi(j,i,im) = 0d0

         end do 
      end do  
c                                 endmember counters:
      if (im.gt.h9) call error (52,dq(1),idqf,'GMODEL')
c                                 number of endmembers
      mstot(im) = istot
c                                 number of independent + ordered endmebers
      nstot(im) = kstot + norder 
c                                 number of independent disordered endmembers
      lstot(im) = kstot 
c                                 chemical mixing sites
      istg(im) = isite
c                                 number of ordered species
      nord(im) = norder 
c                                 max iteration counter
      if (maxit.lt.0) then
         maxitg(im) = iopt(10)
      else 
         maxitg(im) = maxit
      end if 
c                                 number of species and multiplicity and
c                                 site ranges
      ncoor(im) = 0 

      do i = 1, isite 

         ispg(im,i) = isp(i)
         imlt(im,i) = ist(i)
         ncoor(im) = ncoor(im) + isp(i)

         do j = 1, isp(i) - 1
c                                 subdivision override (iopt(13))
            if (iopt(13).eq.1) then
c                                 make linear
               imd(j,i) = 0

            else if (iopt(13).eq.2) then 
c                                 make all stretch (if not already)
               if (imd(j,i).eq.0) imd(j,i) = 2

            end if 
c                                 use default initial resolution, perturbed
c                                 by a per-mil scale increment to reduce 
c                                 compositional degeneracies. 
            if (nopt(13).gt.0d0) 
     *         xnc(i,j) = (1d0 + 5d-3*float(im-5)) * nopt(13)
c                                 save original values as hard limits for
c                                  
            xmno(im,i,j) = xmn(i,j)
            xmxo(im,i,j) = xmx(i,j)
c                                 ------------------------------------
c                                 auto_refine segment
            if (refine) then 

               read (n10,*) xmn(i,j),xmx(i,j)

               if (icopt.lt.4) then 
c                                 set slop to the initial spacing
                  dinc = xnc(i,j)
                
                  if (icopt.eq.1) then
c                                 Schreinemakers use refine factor III
                     xnc(i,j) = xnc(i,j)/dfloat(iopt(9))
                  else  
c                                 non-adaptive use refine factor II
                     xnc(i,j) = xnc(i,j)/dfloat(iopt(8))
                  end if 

               else 
c                                 use fractional slop, with a minimum 
c                                 corresponding to the initial compositional
c                                 resolution
                  dinc = nopt(10)/1d1 + (xmx(i,j) - xmn(i,j)) * nopt(3)
c                                 adaptive use refine factor I
                  xnc(i,j) = xnc(i,j)/dfloat(iopt(7)) 

               end if 

               xmn(i,j) = xmn(i,j) - dinc
               if (xmn(i,j).lt.0d0) xmn(i,j) = 0d0
               xmx(i,j) = xmx(i,j) + dinc
               if (xmn(i,j).lt.0d0) xmn(i,j) = 0d0
 
            end if 
c                                 -------------------------------------
c                                 stretching stuff
            if (imd(j,i).gt.0) then 
c                                 check for old subdivision schemes
               if (imd(j,i).gt.0) then
                  if (xnc(i,j).le.0d0.or.xnc(i,j).gt.1d0)
     *               call error (62,nopt(13),imd(j,i),tname)
               end if

               xsym = (xmxo(im,i,j)+xmno(im,i,j))/2d0

               dinc = xmx(i,j) - xmn(i,j)
c                                 rescale the increment if the limits 
c                                 are reduced
               xnc(i,j) = xnc(i,j)/(xmxo(im,i,j) - xmno(im,i,j))
c                                 get the intervals
               yint(1,j,i,im) = xmn(i,j)    
   
               if (imd(j,i).eq.1.or.imd(j,i).eq.4) then 
c                                 one interval
                  yfrc(1,j,i,im) = 1d0

                  yint(2,j,i,im) = xmx(i,j)

               else if (imd(j,i).eq.2.and.(.not.refine)) then 
c                                 two intervals, symmetry 
c                                 at (xmx(i,j) + xmn(i,j))/2d0.
                  yfrc(1,j,i,im) = (xsym-xmno(im,i,j))/dinc
                  yfrc(2,j,i,im) = 1d0 - yfrc(1,j,i,im)

                  yint(2,j,i,im) = xsym                 
                  yint(3,j,i,im) = xmx(i,j)

               else if (imd(j,i).eq.2) then 
c                                 in autorefine cycle, check if limits 
c                                 are on one-side of original symmetry 
c                                 axis:
                  if (xmx(i,j).le.xsym) then 
c                                 change to assymetric stretch toward xmax
                     imd(j,i) = 1
                     yfrc(1,j,i,im) = 1d0
                     yint(2,j,i,im) = xmx(i,j)

                  else if (xmn(i,j).ge.xsym) then 
c                                 change to assymetric stretch toward xmin
                     imd(j,i) = 4
                     yfrc(1,j,i,im) = 1d0
                     yint(2,j,i,im) = xmx(i,j)                                   

                  else 
c                                 compositions span axis, recompute 
c                                 fractional lengths of intervals:
                     yfrc(1,j,i,im) = (xsym-xmno(im,i,j))/dinc
                     yfrc(2,j,i,im) = 1d0 - yfrc(1,j,i,im)   
                     yint(2,j,i,im) = xsym          
                     yint(3,j,i,im) = xmx(i,j)                     

                  end if 

               else 
c                                 four intervals
                  yint(3,j,i,im) = yin(ms1,mst)
                  yint(5,j,i,im) = xmx(i,j)
                  yint(2,j,i,im) = (yint(1,j,i,im)+yint(3,j,i,im))/2d0
                  yint(4,j,i,im) = (yint(5,j,i,im)+yint(3,j,i,im))/2d0

                  yfrc(1,j,i,im) = (yint(3,j,i,im)-yint(1,j,i,im))/2d0
                  yfrc(2,j,i,im) = yfrc(1,j,i,im) 
                  yfrc(3,j,i,im) = (yint(5,j,i,im)-yint(3,j,i,im))/2d0
                  yfrc(4,j,i,im) = yfrc(3,j,i,im)

               end if 
            end if 

            imdg(j,i,im) = imd(j,i)
            xmng(im,i,j) = xmn(i,j)
            xmxg(im,i,j) = xmx(i,j)
            xncg(im,i,j) = xnc(i,j)

         end do 
      end do 
c                                 -------------------------------------
c                                 classify the model
      ksmod(im) = jsmod
c                                 -------------------------------------
      if (jsmod.eq.1) then 
         write (*,1000) tname
         stop
      end if 

1000  format (/,'Solution model: ',a,' is in microscopic format',/,
     *  /,'to use this model reformulate for macroscopic format',/)

c                                 -------------------------------------
c                                 save the excess terms.                              
      jterm(im) = iterm
      jord(im) = iord

      do i = 1, iterm 
         do j = 1, m3
            wgl(j,i,im) = wg(i,j)
         end do 

         do j = 1, iord
c                                 isub points to the position in the list
c                                 of endmembers potentially including dependent 
c                                 species. use iy2p to convert to independent
c                                 endmember pointers.
            if (isub(i,j,1).eq.0) then 
               jsub(j,i,im) = 0 
            else 
               jsub(j,i,im) = iy2p(isub(i,j,1))                              
            end if 

         end do 
      end do 
c                                 ----------------------------------------------
c                                 configurational entropy models

c                                 site fractions as a function of bulk
c                                 y's and dependent species y:
      nloc = 0

      do i = 1, nsite 
c                                 eliminate sites with 1 species
         if (nspm1(i).eq.0) cycle 

         nloc = nloc + 1
c                                 # of species, and site r*multiplicty.
         qmult(nloc,im) = r*smult(i)
         ksp(nloc,im) = nspm1(i) 
c                                 for each species, read
c                                 function to define the
c                                 site fraction of the species:
         do j = 1, nspm1(i)
c                                 # of terms in the 
c                                 site fraction function and a0.
            lterm(j,nloc,im) = nterm(i,j)
            d0(j,nloc,im) = a0(i,j)
c                                 for each term:
            do k = 1, nterm(i,j)
c                                 term coefficient amd species index:
               dcoef(k,j,nloc,im) = acoef(i,j,k) 
               ksub(k,j,nloc,im) = nsub1(i,j,k,1)

            end do 
         end do
      end do  
c                                 number of distinct identisites for entropy
      msite(im) = nloc

      do i = 1, mstot(im) + nord(im)
c                                 insp points to the original position 
c                                 of endmember i in the solution model input:
         knsp(i,im) = insp(i)
      end do 
c                                 -------------------------------------
c                                 kmsol points to the species on the j'th site
c                                 of the i'th endmember, used for the xtoy
c                                 conversion      
      do i = 1, mstot(im)
         do j = 1, isite
            kmsol(im,i,j) = jmsol(i,j)
         end do 
      end do 
c                                 --------------------------------------
c                                 van laar volumes, and pointers for "in" endmembers
      do i = 1, nstot(im)
c                                 if the solution is van laar save
c                                 the "volume" function.
         if (laar) then  
            do l = 1, m3
               vlar(l,i,im) = vlaar(l,jnsp(i))
            end do 
         end if 
c                                 initialize scoef's to zero for config 
c                                 entropy calculation (done by snorm).
         scoef(i,im) = 0d0

      end do   
c                                 -------------------------------------
      if (depend) then 
c                                 save y -> p array 
         ndep(im) = mdep
         do i = 1, nstot(im)
            do j = 1, mdep
               y2pg(j,i,im) = y2p(i,j)
            end do
         end do
c                                 for reasons of stupidity, convert the z(y) 
c                                 expressions to z(p). 
         call y2p4z (im)

      end if 

c                                 -------------------------------------
c                                 dqf parameters
      jdqf(im) = idqf

      do i = 1, idqf 
c                                 shift pointer from y array to p array
         jndq(i,im) = iy2p(indq(i))
         do j = 1, m3
            dqfg(j,i,im) = dqf(j,i)
         end do 
      end do 
c                                 -------------------------------------
c                                 if nsite ne 0 get "normalization" constants (endmember
c                                 configurational entropies) for entropy model:
      if (nsite.ne.0) call snorm (im,tname)
c                                 -------------------------------------
      if (order) then 
c                                 models with speciation: 
         do j = 1, norder 

            deph(j,im) = denth(j) 

            do i = 1, 2
               ideps(i,j,im) = iy2p(iddeps(i,j))
            end do 
c                                 stoichiometric limits on ordered species
            if (limn(j).gt.0) then

               ln(j,im) = limn(j)
 
               do k = 1, limn(j)
c                                 for each limit, number of p0 terms
                  lt(k,j,im) = limt(k,j)

                  do i = 1, limt(k,j)
c                                 for each p0 term
                     lc(i,k,j,im) = limc(i,k,j)
c                                 convert index from y to p pointer and save
                     lid(i,k,j,im) = iy2p(limid(i,k,j))

                  end do 
c                                 the constant and delta
                  l0c(1,k,j,im) = limc(i,k,j)
                  l0c(2,k,j,im) = limc(i+1,k,j)
c                                 and the number of p terms...
                  jt(k,j,im) = jimt(k,j)

                  do i = 1, jimt(k,j)
c                                 for each p0 term
                     jc(i,k,j,im) = jimc(i,k,j)
c                                 convert index from y to p pointer and save
                     jid(i,k,j,im) = iy2p(jimid(i,k,j))

                  end do 

               end do 

            end if 

         end do 
c                                 classify multiple species models according
c                                 to whether the disordered reactants are 
c                                 partially or completely correlated, assume
c                                 anti-correlation and no correlation are not 
c                                 possible cases.
         icase(im) = 0

         if (norder.gt.1) then 

            imatch = 0

            do j = 1, 2 

               id = ideps(j,1,im)
     
               do i = 1, 2          
                  if (id.eq.ideps(i,2,im)) then
                     imatch = imatch + 1
                     exit
                  end if 
               end do 

            end do 
 
            if (imatch.eq.1) then 
c                                 if match = 1 one species didn't match
c                                 assume partial correlation
               icase(im) = 2
            else if (imatch.eq.2) then 
               icase(im) = 1
            end if 

         end if  
c                                first create derivatives of endmember
c                                fractions with respect to the ordered
c                                species:
         do j = 1, norder 

            do i = 1, nstot(im)
               dydy(i,j,im) = 0d0
            end do
c                                derivative of the ordered species with
c                                respect to itself:
            dydy(kstot+j,j,im) = 1d0
c                                each ordered species decomposes to
c                                two disordered species iddeps(1-2,j)
c                                depvnu is the stoichiometric coefficient
c                                of the disordered species in the ordered
c                                species.

c                                derivatives of the consituent species 
c                                with respect to the ordered species
            do i = 1, 2
               dydy(ideps(i,j,im),j,im) = dydy(ideps(i,j,im),j,im) 
     *                                  - depvnu(i,j)
            end do

         end do 
c                                evaluate the second derivative of each
c                                pi*pj term in the excess function with
c                                respect to kth species
         do i = 1, iterm 
            do j = 1, norder
               do k = 1, norder

                  dppp(k,j,i,im) =  dydy(jsub(1,i,im),k,im)
     *                             *dydy(jsub(2,i,im),j,im)
     *                           + 
     *                              dydy(jsub(2,i,im),k,im)
     *                             *dydy(jsub(1,i,im),j,im)  
               end do                     
            end do 
         end do 

         if (depend) then 
c                                make an array to get the the fractions
c                                of the disordered species in the ordered
c                                species. this is essentially identical to 
c                                the dydy array formed above.                                 
            do i = 1, kstot

               do k = 1, norder

                  dvnu(i,k,im) = 0d0

                  do j = 1, 2
                     if (i.ne.ideps(j,k,im)) cycle
                     dvnu(i,k,im) = dvnu(i,k,im) + depvnu(j,k)
                  end do

               end do  
            end do  

         end if
c                                 site fractions as a function of bulk
c                                 y's and dependent species y:
         do i = 1, msite(im)
c                                 for each species, read
c                                 function to define the
c                                 site fraction of the species:
            do k = 1, norder
               sdzdp(k,ksp(i,im)+1,i,im) = 0d0 
            end do 

            do j = 1, ksp(i,im)
c                                 # of terms in the 
c                                 site fraction function and a0.
               do l = 1, norder
                  sdzdp(l,j,i,im) = 0d0
               end do 
c                                 for each term:
               do k = 1, lterm(j,i,im)
c                                 endmember indexes  
                  ind = ksub(k,j,i,im) 
c                                 get derivatives, of species fractions
c                                 with respect to ordered species
                  do l = 1, norder 
                     itic = 0 
                     do ii = 1, 2
                        if (ind.eq.ideps(ii,l,im)) then
                           sdzdp(l,j,i,im) = sdzdp(l,j,i,im) 
     *                     + dydy(ideps(ii,l,im),l,im)*dcoef(k,j,i,im)  
                           itic = itic + 1
c                                 high order terms not allowed
                           if (itic.gt.1) 
     *                        call error (999,r,801,'GMODEL')
                        end if                       
                     end do
c                                 the derivative of a term with the 
c                                 ordered species.
                     if (ind.eq.kstot+l) 
     *                  sdzdp(l,j,i,im) = sdzdp(l,j,i,im) 
     *                                  + dcoef(k,j,i,im) 

                  end do   
               end do 
            end do 
         end do
c                                 multiply each dzdp by qmult (R*site
c                                 multiplicity) to reduce operation 
c                                 count in evaluating derivatives. 
         do k = 1, norder 
            do i = 1, nsite 

               dzt = 0d0 

               do j = 1, ksp(i,im) 
                  if (dabs(sdzdp(k,j,i,im)).lt.1d-5) 
     *                     sdzdp(k,j,i,im) = 0d0
                  dzt = dzt + sdzdp(k,j,i,im)
               end do 

               if (dabs(dzt).lt.1d-5) dzt = 0d0
               sdzdp(k,j,i,im) = -dzt

            end do 
         end do 

      end if 
c                                 ----------------------------------------------
c                                 models with special endmember indexing:  
      if (jsmod.eq.24.or.jsmod.eq.25) then 
c                                 hp & ghiroso models:
         do i = 1, 3
            ispec(im,i) = 0 
         end do 
c                                 set start index assuming no water:
         ispec(im,4) = 1  

         if (iorig(1).eq.1) then 
c                                 h2o is present:
            ispec(im,1) = 1
c                                 set start index to avoid h2o:
            ispec(im,4) = 2
            if (iorig(2).eq.2) then
               ispec(im,2) = 2
               if (iorig(3).eq.3) ispec(im,3) = 3
            else if (iorig(2).eq.3) then 
               ispec(im,3) = 2
            end if 
         else if (iorig(1).eq.2) then
c                                 h2o absent, fo (in hp) is first endmember
            ispec(im,2) = 1
            if (iorig(2).eq.3) ispec(im,3) = 2
         else if (iorig(1).eq.3) then 
c                                 h2o and fo absent, fa (in hp) is first endmember
            ispec(im,3) = 1
         end if     
      else if (jsmod.eq.0) then
c                                 fluid eos, make pointer to co2
         do i = 1, 2
            id = kdsol(insp(i))
            if (cp(2,id).ne.0d0) then
               ispec(im,1) = i
               exit
            end if 
         end do 

      end if 
c                                 -------------------------------------
c                                 create a y -> x array, this array is 
c                                 to be used to convert endmember fractions (y's)
c                                 back to geometric coordinates (x's).
      do i = 1, isite
         do j = 1, isp(i)
c                                 now find all endmembers with
c                                 species j on site i, this method
c                                 is inefficient but idependent of
c                                 endmember order.   
            do k = 1, istot
               if (jmsol(k,i).eq.j) indx(im,i,j) = k
            end do 
         end do 
      end do 

      if (istot+2.gt.k12) call error (39,0d0,k12,'INPUT9')    

      smod(im) = .true.

      do i = 1, kstot
c                                 pointer to endmember
         id = kdsol(insp(i))

         jend(im,2+i) = id
c                                 set shear modulus flag
         if (iemod(id).eq.0) smod(im) = .false.
c                                 look for endmembers to be killed
         if (iend(insp(i)).ne.2) cycle

         add = .true. 

         do j = 1, killct

            if (killid(j).eq.id) then
               add = .false.
               exit
            end if 
         end do 

         if (add) then 

            killct = killct + 1
            if (killct.gt.10) call error (999,wg(1,1),killct,tname)
            killid(killct) = id

         end if 

      end do 

      if (laar.and.ksmod(im).ne.3) ksmod(im) = 7 
c                                 set type flags, presently no provision for 
c                                 bw summation
      llaar(im) = .false.
      lexces(im) = .false.
      lorder(im) = .false.
      lrecip(im) = .false.
      
      if (iterm.gt.0) then 
         lexces(im) = .true.
         if (laar) llaar(im) = .true.
      end if 

      if (order) lorder(im) = .true.
c                                 the ksmod(im) test is made because
c                                 reform may dump the dependent endmembers
c                                 setting depend = .false., while retaining
c                                 a dummy site with no mixing. reform should
c                                 be redone to truly reformulate multiple
c                                 models to single site models. 
      if (depend.or.ksmod(im).eq.7) lrecip(im) = .true. 

      if (iopt(15).eq.0) then 
c                                 hard limits are off, set limits to 0/1
         do i = 1, isite 
            do j = 1, isp(i) - 1

               xmxo(im,i,j) = 1d0
               xmno(im,i,j) = 0d0

            end do 
         end do
      end if   

      end 

      double precision function strtch (y)
c----------------------------------------------------------------------
c get the x-value from the unstretched coordinate y
c----------------------------------------------------------------------
      implicit none

      double precision y, t

      double precision bp1,bm1,bpm,lbpm
      common/ cst46 /bp1,bm1,bpm,lbpm

      t = bpm**(1d0-y)

      strtch = (bp1 - bm1*t)/(1d0 + t)

      end 

      double precision function unstch (x)
c----------------------------------------------------------------------
c get the y-value from the stretched coordinate x
c----------------------------------------------------------------------
      implicit none

      double precision x

      double precision bp1,bm1,bpm,lbpm
      common/ cst46 /bp1,bm1,bpm,lbpm

      unstch = 1d0 - dlog((bp1-x)/(bm1+x))/lbpm

      end 

      double precision function ydinc (y,xinc,mode,i,ksite,ids)
c----------------------------------------------------------------------
c get the new value of the stretched coordinate (y) from an
c increment in the cartesian coordinate (xinc)
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ids, i, j, mode, ksite

      double precision x, xinc, strtch, yreal, yloc, unstch, dy, ylmn,
     *                 ylmx, y
    
      logical odd
c                                 interval limits conformal transformation
      integer intv
      double precision yint, yfrc
      common/ cst47 /yint(5,ms1,mst,h9),yfrc(4,ms1,mst,h9),intv(4)
c----------------------------------------------------------------------
c                                 y is the real coordinate.
      if (mode.lt.4) then 
         odd = .false.
      else
         odd = .true. 
      end if 
c                                 there are as many as intv(mode)
c                                 intervals to cycle through
      do j = 1, intv(mode)
c                                 odd or even interval?
         odd = .not.odd
c                                 interval limits              
         ylmx = yint(j+1,i,ksite,ids)
c                                 which interval are we starting from?
         if (y.gt.ylmx) cycle

         ylmn = yint(j,i,ksite,ids)
         dy = ylmx - ylmn
c                                 the current value is in interval j
c                                 convert to raw y (varies from 0 ->1 
c                                 over the local interval)
         if (dy.ne.0d0) then 
            yloc = (y-ylmn) / dy
         else 
            yloc = y 
         end if 
c                                 convert to cartesian x and add increment

         if (odd) then 
            x = unstch(yloc) + xinc/yfrc(j,i,ksite,ids)
         else
            x = 1d0 - unstch(1d0-yloc) + xinc/yfrc(j,i,ksite,ids)
         end if 

         if (x.gt.1d0) then
c                                 jumped to next interval
            if (j.lt.intv(mode)) then
c                                 add the residual (opposite odd)
               x = (x - 1d0)*yfrc(j,i,ksite,ids)/yfrc(j+1,i,ksite,ids)
               ylmn = ylmx
               ylmx = yint(j+2,i,ksite,ids)
               dy = ylmx - ylmn

               if (odd) then 
                  yreal = ylmx - strtch(1d0-x) * dy
               else
                  yreal = ylmx + strtch(x) * dy 
               end if 
            else 
               yreal = ylmx
            end if 

         else if (x.ge.0d0) then
c                                 within interval j
            if (odd) then 
               yreal = ylmn + strtch(x) * dy 
            else
               yreal = ylmx - strtch(1d0-x) * dy 
            end if               

         else
c                                 jumped to an earlier interval
            if (j.gt.1) then
c                                 add the residual (opposite odd)
               x = (x + 1d0)*yfrc(j,i,ksite,ids)/yfrc(j-1,i,ksite,ids)
               ylmx = ylmn
               ylmn = yint(j-1,i,ksite,ids)
               dy = ylmx - ylmn

               if (odd) then 
                  yreal = ylmx - strtch(1d0-x) * dy
               else
                  yreal = ylmn + strtch(x) * dy 
               end if 
            else 
               yreal = ylmn
            end if 
         end if

         exit 
 
      end do 

      ydinc = yreal 
         
      end 

      subroutine factor (a,n,ipvt,ier)
c-----------------------------------------------------------------------
c factor is a subroutine which calculates the triangular
c decompositions of the matrix 'a'. factor is modified from
c the subroutine of the same name given by conte and de boor
c in 'elementary numerical analysis', mcgraw-hill, 1980.
c factor uses scaled partial pivoting.
c
c input     a- an n by n array containing the elements of matrix a.
c           n- the dimension of the matrix a.
c output    a- an n by n array containing the upper, u, and lower, l,
c              triangular decompositions of input matrix a.
c        ipvt- a vector indicating that row ipvt(k) was used to
c              eliminate the a(n,k).
c         ier- a flag, zero if a is of rank = n, and 1 if a is of
c              lower rank.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

     
      double precision a(k8,k8),d(k8),rmax,tmax,temp,ratio

      integer ipvt(k8),i,j,k,ier,ip1,n,istr,nm1
c-----------------------------------------------------------------------
      ier = 0
c                            initialize ipvt,d
      do i = 1, n
         ipvt(i) = i
         rmax = 0d0
         do j = 1, n
            rmax = dmax1(rmax,dabs(a(i,j)))
         end do 
c                            ax = b is singular if rmax = 0
         if (dabs(rmax).lt.1d-5) goto 9000
         d(i) = rmax
      end do 
c                            begin decomposition: 
      nm1 = n-1
c
      do i = 1, nm1
         ip1 = i+1
c                            determine pivot row (istr).
         rmax = dabs(a(i,i))/d(i)
         istr = i
         do j = ip1, n
            tmax = dabs(a(j,i))/d(j)
            if (tmax.gt.rmax) then
               rmax = tmax
               istr = j
            end if 
         end do

         if (dabs(rmax).lt.1d-5) goto 9000
c                            if istr gt i, make i the pivot row
c                            by interchanging it with row istr.
         if (istr.gt.i) then 
            j = ipvt(istr)
            ipvt(istr) = ipvt(i)
            ipvt(i) = j
            temp = d(istr)
            d(istr) = d(i)
            d(i) = temp

            do j = 1, n
               temp = a(istr,j)
               a(istr,j) = a(i,j)
               a(i,j) = temp
            end do 
         end if 
c                            eliminate x(k) from rows k+1,...,n.
         do j = ip1,n
         a(j,i) = a(j,i)/a(i,i)
         ratio = a(j,i)
            do k = ip1,n
               a(j,k) = a(j,k)-ratio*a(i,k)
            end do 
         end do
      end do 

      if (dabs(a(n,n)).lt.1d-5) ier = 1

      return
c                           algoritmic singularity.
9000  ier = 1
 
      end

      subroutine y2p0 (id)
c-----------------------------------------------------------------------
c y2p0 converts the y array of disordered dependent and independent 
c species abundance to the p0 array of the independent (ordered and 
c disordered) species. the p0 array gives the minimum possible 
c concentrations of the ordered species (the stable abundances are 
c determined by solving the speciation problem). 

c for non-reciprocal solutions the y and p0 arrays are identical.
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer id,k,l

      double precision y2pg
      common/ cxt4  /y2pg(m15,m4,h9)

      double precision z, pa, p0a, x, w, y
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)

      integer ksmod, ksite, kmsol, knsp
      common/ cxt0  /ksmod(h9),ksite(h9),kmsol(h9,m4,mst),knsp(m4,h9)
c-----------------------------------------------------------------------
c                                 convert y's to p's
c                                 initialize ordered species
      do k = 1, nord(id)
         p0a(lstot(id)+k) = 0d0
      end do        

      do k = 1, nstot(id)
c                                 initialize the independent species
c                                 other then the ordered species
         if (k.le.lstot(id)) p0a(k) =  y(knsp(k,id))
c                                 convert the dependent species to
c                                 idependent species
         do l = 1, ndep(id)
            p0a(k) = p0a(k) + y2pg(l,k,id) * y(knsp(lstot(id)+l,id))
         end do 

         pa(k) = p0a(k)

      end do

      end 

      subroutine p0dord (id)
c-----------------------------------------------------------------------
c decomposes p0a values that specify the abundances of independent 
c ordered endmembers to their stoichiometric equivalent disordered
c species. 
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer id,k,l,ind 

      double precision dvnu,deph,dydy
      common/ cxt3r /dvnu(m4,j3,h9),deph(j3,h9),dydy(m4,j3,h9)

      integer ideps,icase
      common/ cxt3i /ideps(2,j3,h9),icase(h9)

      double precision z, pa, p0a, x, w, y
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)
c-----------------------------------------------------------------------

      do k = 1, nord(id)
         do l = 1, 2
            ind = ideps(l,k,id)
            p0a(ind) = p0a(ind) + dvnu(ind,k,id) * p0a(lstot(id)+k)
         end do 
      end do 

      end 

      subroutine specip (g,id)
c----------------------------------------------------------------------
c subroutine to solve speciation of a solution with one ordering parameter
c and disordered composition p0a by halving. the speciated composition is returned 
c in array pa. 
c    id identifies the solution.
c    g is the change in G for the stable speciation relative to a mechanical
c      mixture of the endmembers. this version uses endmember proportions, see
c      specis for the equivalent problem in terms of endmember fractions.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,id,jd,ind

      double precision g,sign,dy,omega,gex,
     *                 dgdp,odgdp,yt,ymax,tol,sdgdp

      double precision r,tr,pr,ps,p,t,xco2,u1,u2
      common/ cst5   /p,t,xco2,u1,u2,tr,pr,r,ps
c                                 working arrays
      double precision z, pa, p0a, x, w, y
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1)
c                                 new global arrays, 10/25/05:
      integer ideps,icase
      common/ cxt3i /ideps(2,j3,h9),icase(h9)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)

      double precision dvnu,deph,dydy
      common/ cxt3r /dvnu(m4,j3,h9),deph(j3,h9),dydy(m4,j3,h9)
c                                 model type
      logical lorder, lexces, llaar, lrecip
      common/ cxt27 /lorder(h9),lexces(h9),llaar(h9),lrecip(h9)

      save tol

      data tol/1d-5/
c----------------------------------------------------------------------
c                                 two cases, ksmod = 6 non-reciprocal
c                                 and ksmod = 8 reciprocal
      jd = lstot(id) + 1 

      if (.not.lrecip(id)) then 
         do i = jd, nstot(id)
            p0a(i) = 0d0
         end do 
      end if 
c                                 find the maximum proportion of the 
c                                 ordered species
      ymax = 1d0
c                                 get stoichiometric limit on 
c                                 fraction of the ordered species
      do i = 1, 2
         yt = -p0a(ideps(i,1,id))/dydy(ideps(i,1,id),1,id)
         if (yt.lt.ymax) ymax = yt
      end do 
c                                 no ordering possible if yt = 0
      if (ymax-p0a(jd).lt.tol) then 
         pa(jd) = p0a(jd)
         goto 90
      end if 

      ymax = ymax - tol/1d3
c                                 solve by halving
      pa(jd) = p0a(jd) + ymax

      do i = 1, 2
         ind = ideps(i,1,id)
         pa(ind) = p0a(ind) + dydy(ind,1,id)*ymax
      end do 

      odgdp = dgdp(id)
c                                 if dgdp > 0 must be fully ordered
      if (odgdp.lt.0d0) goto 90
c                                 initialize at halfway point
      dy = ymax/2d0
      pa(jd) = pa(jd) - dy
c                                 iteration loop:
10    do i = 1, 2
         ind = ideps(i,1,id)
         pa(ind) = p0a(ind) + dydy(ind,1,id)*(pa(jd)-p0a(jd))
      end do 

      sdgdp = dgdp(id)

      sign = odgdp*sdgdp

      if (sign.gt.0d0) then 

         odgdp = sdgdp
         if ((pa(jd)-p0a(jd)).eq.tol) goto 90
         pa(jd) = pa(jd) - dy
         if ((pa(jd)-p0a(jd)).lt.tol) pa(jd) = p0a(jd) + tol
         goto 10 

      else if (dy.gt.tol) then 

         dy = dy / 2d0
         pa(jd) = pa(jd) + dy
         goto 10

      end if 
c                                 sum the enthalpic + entropic  + excess contributions
90    g = pa(jd)*deph(1,id) - t*omega(id,pa) + gex(id,pa)
c                                 convert the ordered endmember fractions to 
c                                 disordered fractions (stored in the p0a array).      
c     if (lrecip(id)) call p0dord (id)

      end 

      subroutine y2p4z (id)
c----------------------------------------------------------------------
c subroutine to convert site fraction expressions in terms of 
c y's to p's, this conversion consists simply of eliminating 
c the dependent endmembers and coverting the pointer array 
c jsub is changed from a pointer to the y array to a pointer
c to the p array.

c the reason for writing z(y) is it is a more intuitive input,
c i.e., users do not have to understand the distinction between
c dependent and independent endmembers.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision ncoef(m0)

      integer i,j,k,id,ind,jconv,mterm(m0)
c                                 -------------------------------------
      integer msite, ksp, lterm, ksub
      common/ cxt1i /msite(h9),ksp(m10,h9),lterm(m11,m10,h9),
     *               ksub(m0,m11,m10,h9)

      double precision qmult, d0, dcoef, scoef      
      common/ cxt1r /qmult(m10,h9),d0(m11,m10,h9),dcoef(m0,m11,m10,h9),
     *               scoef(m4,h9)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)

      integer jmsol,kdsol
      common/ cst142 /jmsol(m4,mst),kdsol(m4)

      integer iorig,jnsp,iy2p
      common / cst159 /iorig(m4),jnsp(m4),iy2p(m4)
c----------------------------------------------------------------------
c                                 for each site
      do i = 1, msite(id)
c                                 for each species
         do j = 1, ksp(i,id)
 
            jconv = 0 

            do k = 1, lterm(j,i,id)

               ind = ksub(k,j,i,id)

               if (kdsol(ind).ne.-2) then
c                                 index points to an independent endmember
                  jconv = jconv + 1
                  mterm(jconv) = ind
                  ncoef(jconv) = dcoef(k,j,i,id)

               end if 

            end do 
c                                 load the reformulated function
c                                 into lterm, ksub, and dcoef
            lterm(j,i,id) = jconv
 
            do k = 1, jconv

               dcoef(k,j,i,id) = ncoef(k)
               ind = mterm(k)

               if (ind.le.mstot(id)) then 
                  ksub(k,j,i,id) = iy2p(ind)
               else 
                  ksub(k,j,i,id) = mterm(k) - ndep(id)
               end if 

            end do                       

         end do 
 
      end do  

      end

      subroutine specis (g,id)
c----------------------------------------------------------------------
c subroutine to speciation of a solution with disordered composition p0a. 
c the speciated composition is returned in array pa. 
c    id identifies the solution.
c    g is the change in G for the stable speciation relative to a mechanical
c      mixture of the endmembers.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,id

      logical error

      double precision g, gdord, omega, gex

      external omega, gex

      double precision z, pa, p0a, x, w, y
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)

      logical lorder, lexces, llaar, lrecip
      common/ cxt27 /lorder(h9),lexces(h9),llaar(h9),lrecip(h9)

      double precision r,tr,pr,ps,p,t,xco2,u1,u2
      common/ cst5   /p,t,xco2,u1,u2,tr,pr,r,ps

      double precision dvnu,deph,dydy
      common/ cxt3r /dvnu(m4,j3,h9),deph(j3,h9),dydy(m4,j3,h9)

      integer iopt
      double precision nopt
      common / opts /nopt(i10),iopt(i10)
c----------------------------------------------------------------------
      g = 0d0 

      if (.not.lrecip(id)) then 
         do i = lstot(id)+1, nstot(id)
            p0a(i) = 0d0
            pa(i) = p0a(i)
         end do 
      end if 
c                                 as most models are single species and
c                                 there is so much overhead in computing
c                                 multiple speciation, use a special routine
c                                 for single species models:
      if (nord(id).eq.1) then
 
         call speci1 (g,id,1,error) 

      else 

         call speci2 (g,id,error)

      end if 

      if (error.or.iopt(17).ne.0) then 
c                                 if speciation returns error, or order_check is on,
c                                  i.e., iopt(17).ne.0, compute disordered g.                               
         gdord =  gex(id,p0a) - t*omega(id,p0a)
 
         if (lrecip(id)) then 
            do i = 1, nord(id)
               gdord = gdord + p0a(lstot(id)+i)*deph(i,id)
            end do 
         end if 

         if (error) then 

            g = gdord

         else

            if (gdord.lt.g) g = gdord

         end if 

      end if 
c                                 convert the ordered endmember fractions to 
c                                 disordered fractions (stored in the p0a array).      
      if (lrecip(id)) call p0dord (id)

      end 

      subroutine gderiv (id,g,dp,error)
c----------------------------------------------------------------------
c subroutine to compute the g of a solution (id) and it's 1st and 2nd 
c derivatives with respect to the oncentrations of nord(id) ordered 
c species. the formulation assumes atomic site fractions are linear 
c functions of the ordered species concentrations (p's) and that the 
c excess function is second order.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical error

      integer i,k,l,i1,i2,id,norder,ipvt(j3)

      double precision g,dp(j3),t,s,ds(j3),d2s(j3,j3)
c                                 working arrays
      double precision z, pa, p0a, x, w, y
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1)

      double precision alpha,dt
      common/ cyt0  /alpha(m4),dt(j3)
c                                 excess energy variables
      integer jterm, jord, jsub
      common/ cxt2i /jterm(h9),jord(h9),jsub(m2,m1,h9)

      double precision dvnu,deph,dydy
      common/ cxt3r /dvnu(m4,j3,h9),deph(j3,h9),dydy(m4,j3,h9)
c                                 model type
      logical lorder, lexces, llaar, lrecip
      common/ cxt27 /lorder(h9),lexces(h9),llaar(h9),lrecip(h9)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)

      double precision dppp,d2gx,sdzdp
      common/ cxt28 /dppp(j3,j3,m1,h9),d2gx(j3,j3),sdzdp(j3,m11,m10,h9)

      double precision p,tk,xc,u1,u2,tr,pr,r,ps
      common/ cst5 /p,tk,xc,u1,u2,tr,pr,r,ps

      double precision dg,d2g
      common/ cyt1 /dg(j3),d2g(j3,j3)

      logical pin
      common/ cyt2 /pin(j3)
c----------------------------------------------------------------------
c                                 initialize, d2gx has been set in setw
      g = 0d0

      norder = nord(id)

      do k = 1, norder
         dg(k) = 0d0
         do l = k, norder
            d2g(l,k) = d2gx(l,k)
         end do 
      end do 

      if (lexces(id)) then

         do i = 1, jterm(id)
c                                 assuming regular terms
           i1 = jsub(1,i,id)
           i2 = jsub(2,i,id)

           g = g + w(i) * pa(i1) * pa(i2)

           do k = 1, norder

              if (.not.pin(k)) cycle

              dg(k) = dg(k) + w(i) * (pa(i1)*dydy(i2,k,id) 
     *                              + pa(i2)*dydy(i1,k,id))
           end do 

         end do  
c                                 get derivative of excess function
         if (llaar(id)) then 

            t = 0d0 
c                                 h&p van laar
            do i = 1, nstot(id)
               t = t + alpha(i)* pa(i)
            end do 
c                                 coming out of this loop g, dg, and
c                                 d2g  are not the complete functions 
c                                 because of the "tphi" term in the 
c                                 van laar.

            do k = 1, norder

               if (.not.pin(k)) cycle
c                                 convert dg and d2g to the full derivative
               dg(k) = (dg(k) - g*dt(k)/t)/t
               do l = k, norder
                  d2g(l,k) = (d2g(l,k) - 2d0*dt(k)*dg(k))/t
               end do 
            end do 
c                                 and the full excess energy 
            g = g/t

         end if 

      end if 
c                                 get the configurational entropy derivatives
      call sderiv (id,s,ds,d2s)

      do k = 1, norder

         if (.not.pin(k)) cycle

         g = g + deph(k,id) * pa(lstot(id)+k)
c                                 dg is the negative of the differential of g 
c                                 with respect to the kth species.
         dg(k) = -(deph(k,id) + dg(k) - tk*ds(k))
         do l = k, norder 
            d2g(l,k) = d2g(l,k) - tk*d2s(l,k)
         end do 
      end do
c                                 determininats, to check for a saddle point
c      if (norder.eq.2) then 
c         detg = d2g(1,1)*d2g(2,2)-d2g(2,1)**2
c      else 
c         detg = d2g(1,1)*(d2g(2,2)*d2g(3,3)-d2g(3,2)**2) 
c     *        - d2g(2,2)*d2g(1,3)**2
c     *        + 2d0*d2g(2,1)*d2g(3,2)*d2g(3,1)-d2g(2,1)**2*d2g(3,3)
c      end if 

      g = g - tk*s 
c                                 copy dg and d2g into dp and d2s
      do k = 1, norder
         if (pin(k)) then 
            dp(k) = dg(k)
            d2s(k,k) = d2g(k,k)
            do l = k+1, norder
               if (pin(l)) then  
                  d2s(l,k) = d2g(l,k)
                  d2s(k,l) = d2g(l,k) 
               end if 
            end do 
         end if 
      end do 

      do k = 1, norder
         if (.not.pin(k)) then
            dp(k) = 1d0 
            d2s(k,k) = 1d0 
            do l = 1, norder
               if (l.eq.k) cycle 
               d2s(l,k) = 0d0
               d2s(k,l) = 0d0 
            end do 
         end if 
      end do     
c                                 get the newton-raphson increments:
c                                 this is a general factorization routine, should
c                                 exploit that d2g is symmetric.
      call factr2 (d2s,j3,norder,ipvt,error)
c                                 solve for the increments by back-substitution, 
c                                 this routine is also not efficient and should 
c                                 be re written. 
      if (.not.error) call subst2 (d2s,ipvt,j3,norder,dp,error) 
c                                 substitute replaces the values of dg with the 
c                                 newton-raphson increments for the ordered species
c                                 compositions. 
      end

      subroutine sderiv (id,s,dsy,dsyy)
c----------------------------------------------------------------------
c subroutine to the derivative of the configurational entropy of a 
c solution with respect to the proportion of a dependent species.

c THIS DOES NOT INCLUDE ENDMEMBER CONFIGURATION ENTROPY DERIVATIVES!
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,j,k,l,id

      logical inf

      double precision zt,dzdy,s,dsy(j3),dsyy(j3,j3),q,zl,
     *                 z(m11,m10),ds,d2s,s0,ztemp,zlnz,
     *                 dsinf(j3),d2sinf(j3,j3)
c                                 working arrays
      double precision zz, pa, p0a, x, w, y
      common/ cxt7 /zz(m4),y(m4),pa(m4),p0a(m4),x(mst,msp),w(m1)
c                                 configurational entropy variables:
      integer msite, ksp, lterm, ksub
      common/ cxt1i /msite(h9),ksp(m10,h9),lterm(m11,m10,h9),
     *               ksub(m0,m11,m10,h9)

      double precision qmult, d0, dcoef, scoef      
      common/ cxt1r /qmult(m10,h9),d0(m11,m10,h9),dcoef(m0,m11,m10,h9),
     *               scoef(m4,h9)

      double precision dppp,d2gx,sdzdp
      common/ cxt28 /dppp(j3,j3,m1,h9),d2gx(j3,j3),sdzdp(j3,m11,m10,h9)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)

      logical pin
      common/ cyt2 /pin(j3)
c----------------------------------------------------------------------
      s = 0d0
c                                 for each site
      do i = 1, msite(id)

         zt = 0d0
         s0 = zt 
c                                 get site fractions
         do j = 1, ksp(i,id)

            ztemp = d0(j,i,id)
c                                 for each term:
            do k = 1, lterm(j,i,id)
               ztemp = ztemp + dcoef(k,j,i,id) * pa(ksub(k,j,i,id))
            end do  

            zt = zt + ztemp
            if (ztemp.gt.0d0) s0 = s0 - ztemp * dlog (ztemp)
            z(j,i) = ztemp

         end do 

         ztemp = 1d0 - zt
         if (ztemp.gt.0d0) s0 = s0 - ztemp * dlog (ztemp)
         z(j,i) = ztemp
         s = s + qmult(i,id) * s0

      end do 
c                                 endmember corrections
      do i = 1, nstot(id)
         s = s - pa(i)*scoef(i,id)
      end do 

c                                 initialize derivatives:
      inf = .false.

      do k = 1, nord(id)
         dsy(k) = 0d0 
         dsinf(k) = 0d0 
         do l = k, nord(id) 
            dsyy(l,k) = 0d0 
            d2sinf(l,k) = 0d0 
         end do 
      end do  
c                                 evaluate derivatives:
      do i = 1, msite(id)

         ds = 0d0  
         d2s = 0d0
         q = qmult(i,id)
 
         do j = 1, ksp(i,id) + 1   

            zl = z(j,i) 

            if (zl.gt.0d0) then 
               zlnz = 1d0 + dlog(zl)
            else
               zlnz = 1d0 
            end if 

            do k = 1, nord(id)        
c                                 skip species not in the model
               if (.not.pin(k)) cycle
c                                 sdzdp is (dz(i,j)/dp(k))
               dzdy = sdzdp(k,j,i,id)

               if (dzdy.eq.0d0) cycle 

               if (zl.gt.0d0) then 
c                                 the first derivative is
                  dsy(k) = dsy(k) - q * dzdy * zlnz
c                                 and the jacobians are
                  do l = k, nord(id)

                     if (.not.pin(l)) cycle 
                     dsyy(l,k) = dsyy(l,k) 
     *                         - q * dzdy * sdzdp(l,j,i,id) / zl
                  end do 

               else if (zl.lt.-1d-6) then 

                   write (*,*) 'wacka boom',zl,j,i
                   write (*,*) (p0a(l),l=1,8)
                   write (*,*) (pa(l),l=1,8)
                   write (*,*) 

               else 

                  inf = .true.
c                                 a species with a non-zero
c                                 derivative is zero, the 
c                                 first will be sign(dzdy)*infinity
                  dsinf(k) = dsinf(k) + dsign(q,dzdy)

                  do l = k, nord(id)
c                                 the 2nd will be -sign of 
c                                 cross term * infinity
                     if (.not.pin(l)) cycle 
                     d2sinf(l,k) = dsinf(k) -  
     *                             dsign(q,dzdy * sdzdp(l,j,i,id))
                  end do 
               
               end if 
            end do  

         end do 

      end do 

      if (inf) then

         do k = 1, nord(id)

            if (.not.pin(k)) cycle 
            if (dabs(dsinf(k)).gt.1d-5) dsy(k) = 1d8*dsinf(k)

            do l = k, nord(id)
               if (.not.pin(l)) cycle 
               if (dabs(d2sinf(l,k)).gt.1d-5) 
     *                                  dsyy(l,k) = 1d10*d2sinf(l,k)
            end do  
 
         end do 

      end if 

      end

      subroutine factr2 (a,m,n,ipvt,error)
c-----------------------------------------------------------------------
c factr2 is a subroutine which calculates the triangular
c decompositions of the matrix 'a'. factor is modified from
c the subroutine of the same name given by conte and de boor
c in 'elementary numerical analysis', mcgraw-hill, 1980.
c factor uses scaled partial pivoting.
c
c input     a- an n by n array containing the elements of matrix a.
c           n- the actual dimension of matrix a
c           m- the phsical dimension of matrix a
c output    a- an n by n array containing the upper, u, and lower, l,
c              triangular decompositions of input matrix a.
c        ipvt- a vector indicating that row ipvt(k) was used to
c              eliminate the a(n,k).
c       error- false if a is of rank = n, and true if a is of
c              lower rank.
c-----------------------------------------------------------------------
      implicit none

      logical error 

      include 'perplex_parameters.h'

      integer m,ipvt(m),i,j,k,ier,ip1,n,istr
    
      double precision a(m,m),d(m),rmax,tmax,temp,ratio
c-----------------------------------------------------------------------
      error = .false.
c                            initialize ipvt,d
      do i = 1, n

         ipvt(i) = i
         rmax = 0d0

         do j = 1, n
            rmax = dmax1(rmax,dabs(a(i,j)))
         end do 
c                            ax = b is singular if rmax = 0
         if (dabs(rmax).lt.1d-5) then 
            error = .true.
            return 
         end if 

         d(i) = rmax

      end do 
c                            begin decomposition: 
      do i = 1, n-1
c                            determine pivot row (istr).
         ip1 = i+1
         rmax = dabs(a(i,i))/d(i)
         istr = i

         do j = ip1, n

            tmax = dabs(a(j,i))/d(j)

            if (tmax.gt.rmax) then
               rmax = tmax
               istr = j
            end if 

         end do

         if (dabs(rmax).lt.1d-5) then 
            ier = 1
            return
         end if 
c                            if istr gt i, make i the pivot row
c                            by interchanging it with row istr.
         if (istr.gt.i) then 

            j = ipvt(istr)
            ipvt(istr) = ipvt(i)
            ipvt(i) = j
            temp = d(istr)
            d(istr) = d(i)
            d(i) = temp

            do j = 1, n
               temp = a(istr,j)
               a(istr,j) = a(i,j)
               a(i,j) = temp
            end do 

         end if 
c                            eliminate x(k) from rows k+1,...,n.
         do j = ip1, n

            a(j,i) = a(j,i)/a(i,i)
            ratio = a(j,i)

            do k = ip1, n
               a(j,k) = a(j,k)-ratio*a(i,k)
            end do 

         end do

      end do 

      if (dabs(a(n,n)).lt.1d-5) ier = 1
 
      end

      subroutine subst2 (a,ipvt,m,n,b,error)
c-----------------------------------------------------------------------
c subst uses the lu decomposition of the matrix 'a' contained
c in the array 'a' to solve ax = b for x. subst is modified from the
c the subroutine of the same name listed by conte and de boor
c in 'elementary numerical analysis', mcgraw-hill, 1980.
c factor uses scaled partial pivoting.
 
c input     a- an n by n array containing the non-zero elements of
c              the u and l decompositions of a, as output by factor.
c           n- the actual dimension of the matrix a.
c           m- the physical dimension of a.
c        ipvt- a vector indicating that row ipvt(k) was used to
c              eliminate the coefficient a(n,k).
c           b- the vector b.
c output    b- the solution vector x.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical error 
 
      integer m,ipvt(m),ip,i,j,n,ii

      double precision a(m,m),b(m),x(m),sum
c----------------------------------------------------------------------
c                                 solve ly = b for y:
      ip = ipvt(1)
      x(1) = b(ip)

      do i = 2, n

         sum = 0d0

         do j = 1, i - 1
            sum = a(i,j)*x(j) + sum
         end do 

         ip = ipvt(i)
         x(i) = b(ip)-sum

      end do 
c                                 solve ux = y for x:
      if (a(n,n).eq.0d0) then
c                                 this check should be superfluous. should check
c                                 what's with factor. 
         error = .true.
         return
      end if 

      x(n) = x(n)/a(n,n)

      do ii = 1, n - 1

         i = n-ii

         sum = 0d0

         do j = i + 1, n
            sum = a(i,j)*x(j)+sum
         end do 

         if (a(i,i).eq.0d0) then
c                                 as above.
            error = .true.
            return
         end if 

         x(i) = (x(i)-sum)/a(i,i)
         b(i) = x(i)

      end do 

      b(n) = x(n)

      end

      subroutine speci1 (g,id,k,error)
c----------------------------------------------------------------------
c subroutine to speciation of a solution with a single ordering parameter
c and disordered composition p0a by newton raphson. the speciated 
c composition is returned 
c in array pa. 
c    id identifies the solution.
c    g is the change in G for the stable speciation relative to a mechanical
c      mixture of the endmembers.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i1,i2,id,jd,k,itic

      logical error 

      double precision g,pt,pmax,pmin,dy1,dy2,dp,dpmax,tol,
     *                 omega,gex,dg,d2g

      external omega, gex

      double precision z, pa, p0a, x, w, y
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1)

      double precision r,tr,pr,ps,p,t,xco2,u1,u2
      common/ cst5   /p,t,xco2,u1,u2,tr,pr,r,ps

      integer ideps,icase
      common/ cxt3i /ideps(2,j3,h9),icase(h9)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)

      double precision dvnu,deph,dydy
      common/ cxt3r /dvnu(m4,j3,h9),deph(j3,h9),dydy(m4,j3,h9)

      logical pin
      common/ cyt2 /pin(j3)

      save tol
      data tol/1d-5/
c----------------------------------------------------------------------
      i1 = ideps(1,k,id)
      i2 = ideps(2,k,id)
      dy1 = dydy(i1,k,id)
      dy2 = dydy(i2,k,id)
      jd = lstot(id) + k 
      error = .false.
c                                 find the maximum proportion of the 
c                                 ordered species cannot be > the amount 
c                                 of reactant initially present
      dpmax = dmin1(-pa(i1)/dy1,-pa(i2)/dy2)
c                                 to avoid singularity set the initial 
c                                 composition to the max - tol, at this
c                                 condition the first derivative < 0, 
c                                 and the second derivative > 0 (otherwise
c                                 the root must lie at p > pmax - tol.                
      if (dpmax.gt.0d0) then

         pin(k) = .true.
         dp = dpmax - tol
         pmax = pa(jd) + dp
         pmin = pa(jd) + tol
c                                 get starting end for the search
c                                 first try the maximum
         pa(jd) = pa(jd) + dp 
         pa(i1) = pa(i1) + dy1*dp
         pa(i2) = pa(i2) + dy2*dp 

         call gderi1 (k,id,dg,d2g)

         if (dg.gt.0d0.and.d2g.gt.0d0) then 
c                                 at the maximum concentration, the 
c                                 first derivative is positive, if 
c                                 the second is also > 0 then we're 
c                                 business
            dp = -dg/d2g

         else if (dg.lt.0d0) then
c                                 then saturated with the ordered 
c                                 species
            pa(jd) = p0a(jd) + dpmax
            pa(i1) = p0a(i1) + dy1*dpmax
            pa(i2) = p0a(i2) + dy2*dpmax    

            goto 90       

         else
c                                 try the min
            pa(jd) = pmin
            pa(i1) = p0a(i1) + dy1*tol
            pa(i2) = p0a(i2) + dy2*tol
 
            call gderi1 (k,id,dg,d2g)

            if (dg.lt.0d0.and.d2g.gt.0d0) then 
c                                 ok
               dp = -dg/d2g
            
            else                
c                                 full disordered
               error = .true.
               goto 99              

            end if 
         end if 

         pt = pa(jd) + dp
c                                 check bounds 
         if (pt.lt.pmin) then
                             
            pa(jd) = pa(jd) + (pmin - pa(jd))/2d0

         else if (pt.gt.pmax) then
 
            pa(jd) = pa(jd) + (pmax - pa(jd))/2d0
          
         else 

            pa(jd) = pt         

         end if
c                                 iteration counter to escape
c                                 infinite loops
         itic = 0 
c                                 newton raphson iteration
         do 

            call gderi1 (k,id,dg,d2g)

            dp = -dg/d2g 

            pt = pa(jd) + dp 

            if (pt.lt.pmin) then
c                                 increment would make p < pmin
c                                 switch the starting guess to pmax 
               dp = (pmin - pa(jd))/2d0
               pa(jd) = pa(jd) + dp

            else if (pt.gt.pmax) then
c                                 increment would make p > pmax
c                                 switch the starting guess to pmin 
               dp = (pmax - pa(jd))/2d0
               pa(jd) = pa(jd) + dp 

            else if (pt.eq.pmin.or.pt.eq.pmax) then 

               exit 

            else 

               pa(jd) = pt         

            end if

            pa(i1) = pa(i1) + dy1*dp
            pa(i2) = pa(i2) + dy2*dp 

            if (dabs(dp).lt.tol) then 
               exit
            else 
               itic = itic + 1
               if (itic.gt.20) then
                  exit
               end if 
            end if 

         end do

      end if  

90    g = pa(jd)*deph(k,id) - t*omega(id,pa) + gex(id,pa)

99    end

      subroutine speci2 (g,id,error)
c----------------------------------------------------------------------
c subroutine to multiple speciation of a solution with disordered composition 
c p0a. the speciated composition is returned in array pa. 
c    id identifies the solution.
c    g is the change in G for the stable speciation relative to a mechanical
c      mixture of the endmembers.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical error

      integer i,k,id,lord,itic,ifail,igood

      double precision g,dp(j3),tdp,tol,omega,gex,gold,xtdp

      external omega, gex

      double precision z, pa, p0a, x, w, y
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1)

      double precision r,tr,pr,ps,p,t,xco2,u1,u2
      common/ cst5   /p,t,xco2,u1,u2,tr,pr,r,ps

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)

      double precision dvnu,deph,dydy
      common/ cxt3r /dvnu(m4,j3,h9),deph(j3,h9),dydy(m4,j3,h9)

      logical lorder, lexces, llaar, lrecip
      common/ cxt27 /lorder(h9),lexces(h9),llaar(h9),lrecip(h9)

      logical pin
      common/ cyt2 /pin(j3)

      integer igood0,jtic
      common/junk/ifail,igood,igood0,jtic

      save tol
      data tol,ifail,igood,igood0/1d-5,0,0,0/
c----------------------------------------------------------------------
c                                 initialize limit expressions
      call p0limt (id)
c                                 get initial p values
      call pinc0 (id,lord) 
c                                 lord is the number of possible species
      if (lord.eq.1) then 

         do i = 1, nord(id)
            if (pin(i)) then 
               call speci1 (g,id,i,error)
               exit 
            end if 
         end do

      else if (lord.gt.1) then 

         itic = 0  
         gold = 0d0
         xtdp = 0d0

         do 

            call gderiv (id,g,dp,error)

            if (error) exit

            tdp = 0d0 

            do k = 1, nord(id)
               
               if (.not.pin(k)) cycle

               call pinc (dp(k),k,id)

               tdp = tdp + dabs(dp(k))

            end do 

            if (tdp.lt.tol.and.gold-g.lt.1d2) then 

               if (itic.lt.2) then 
                  igood0 = igood0 + 1
               end if 

               jtic = jtic + itic
               igood = igood + 1

               exit

            end if 

            if (g.gt.gold.and.gold.ne.0d0) then 
 
               if (tdp.eq.xtdp) then 
                  ifail = ifail + 1
                  exit 
               end if 

               xtdp = tdp

            end if 

            gold = g 

            itic = itic + 1

            if (itic.eq.16) then
               ifail = ifail + 1
               exit 
            end if 

         end do

      else 
c                                 no speciation possible, but still need
c                                 to calculate g (setting error will do this). 
         error = .true.

      end if 

      end                  

      subroutine pinc (dp,k,id)
c----------------------------------------------------------------------
c subroutine to increment the k'th species of solution id, if the increment
c violates a stoichiometric limit, it's set to half it's maximum value.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer k,i1,i2,id,jd

      double precision dp,pmx,pmn,tol
c                                 working arrays
      double precision z, pa, p0a, x, w, y
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1)

      integer ideps,icase
      common/ cxt3i /ideps(2,j3,h9),icase(h9)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)

      double precision dvnu,deph,dydy
      common/ cxt3r /dvnu(m4,j3,h9),deph(j3,h9),dydy(m4,j3,h9)

      save tol

      data tol/1d-5/
c----------------------------------------------------------------------
c                                 given dp check if it violates
c                                 stoichiometric constraints
      i1 = ideps(1,k,id)
      i2 = ideps(2,k,id)
      jd = lstot(id) + k 

      tol = 1d-7

      call plimit (pmn,pmx,k,id)       

      if (pa(jd)+dp.gt.pmx) then 
         dp = pmx - pa(jd) - tol
      else if (pa(jd)+dp.lt.pmn) then 
         dp = pmn - pa(jd) + tol
      end if  
c                                 adjust the composition by the increment
      pa(i1) = pa(i1) + dydy(i1,k,id)*dp
      pa(i2) = pa(i2) + dydy(i2,k,id)*dp 
      pa(jd) = pa(jd) + dp

      end 

      subroutine pinc0 (id,lord)
c----------------------------------------------------------------------
c subroutine set initial species concentrations to half their 
c stoichiometric limit.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,j,k,i1,i2,id,jd,lord,iout,ibad(m4)

      double precision dp,pmn,pmx,tol,dpp(j3),dinc,tinc
c                                 working arrays
      double precision z, pa, p0a, x, w, y
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1)

      integer ideps,icase
      common/ cxt3i /ideps(2,j3,h9),icase(h9)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)

      double precision dvnu,deph,dydy
      common/ cxt3r /dvnu(m4,j3,h9),deph(j3,h9),dydy(m4,j3,h9)

      logical pin
      common/ cyt2 /pin(j3)
c----------------------------------------------------------------------

      lord = 0
      tol = 1d-6

      if (icase(id).eq.1) then 
c                                 case 1: fully correlated
         dinc = 0.9d0/dfloat(nord(id))
         tinc = dinc

         do k = 1, nord(id)

            call plimit (pmn,pmx,k,id) 

            if (pmn.ge.pmx) then 
               pin(k) = .false.
               cycle 
            else 
               pin(k) = .true.
               lord = lord + 1
            end if 

            jd = lstot(id) + k 
            i1 = ideps(1,k,id)
            i2 = ideps(2,k,id)

            dp = pmn + (pmx - pmn) * tinc - pa(jd)
c                                 adjust the composition by the first increment
            pa(i1) = pa(i1) + dydy(i1,k,id)*dp
            pa(i2) = pa(i2) + dydy(i2,k,id)*dp 
            pa(jd) = pa(jd) + dp

            tinc = tinc + dinc

         end do

      else if (icase(id).eq.2) then 
c                                 case 2: positive partial correlation
         do i = 1, 5

            do k = 1, nord(id)

               call plimit (pmn,pmx,k,id) 

               if (i.eq.1) then 

                  if (pmn.ge.pmx) then 
                     pin(k) = .false.
                     cycle 
                  else 
                     pin(k) = .true.
                     lord = lord + 1
                 end if 

               end if 
c                                 adjust the composition by the first increment
               jd = lstot(id) + k 
               dp = pmx - pa(jd) 
               pa(jd) = pa(jd) + dp
               dpp(k) = pa(jd) - p0a(jd)

            end do 

            if (dp.lt.1d-2) exit 

         end do 
c                                 back off from maximum for final assignements
         do k = 1, nord(id)

            if (.not.pin(k)) cycle

            jd = lstot(id) + k 
            pa(jd) = p0a(jd)
            i1 = ideps(1,k,id)
            i2 = ideps(2,k,id)

            dp = dpp(k)*0.9
c                                 adjust the composition by the first increment
            pa(i1) = pa(i1) + dydy(i1,k,id)*dp
            pa(i2) = pa(i2) + dydy(i2,k,id)*dp 
            pa(jd) = pa(jd) + dp

         end do 

      else 
c                                 uncorrelated,anticorrelated? 
c       pause 

      end if
c                                 check for degenerate compositions
      if (lord.gt.0) then 

         iout = 0 

         do i = 1, lstot(id)
            if (p0a(i).eq.0d0) then
               iout = iout + 1
               ibad(iout) = i 
            end if
         end do 
c                                 the indices of the present components are igood(1..in)
         if (iout.gt.0) then 
            do k = 1, nord(id)
               if (pin(k)) then    
c                                 check that the ordered species are in the subcomposition
                  do j = 1, 2                                    
                     do i = 1, iout
                        if (ideps(j,k,id).eq.ibad(i)) then 
                           lord = 0 
                           goto 99 
                        end if 
                     end do 
                  end do 
               end if
            end do 
         end if 
      end if  

99    end 

      subroutine gderi1 (k,id,dg,d2g)
c----------------------------------------------------------------------
c subroutine to compute the 1st and 2nd derivatives of the g of 
c solution (id)  with respect to the concentrations of the kth ordered. 
c the formulation assumes atomic site fractions are linear 
c functions of the ordered species concentrations (p's) and that the 
c excess function is second order.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,k,i1,i2,id,norder

      double precision g,dg,d2g,t,ds,d2s,omega
c                                 working arrays
      double precision z, pa, p0a, x, w, y
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1)

      double precision alpha,dt
      common/ cyt0  /alpha(m4),dt(j3)
c                                 excess energy variables
      integer jterm, jord, jsub
      common/ cxt2i /jterm(h9),jord(h9),jsub(m2,m1,h9)

      double precision dvnu,deph,dydy
      common/ cxt3r /dvnu(m4,j3,h9),deph(j3,h9),dydy(m4,j3,h9)
c                                 model type
      logical lorder, lexces, llaar, lrecip
      common/ cxt27 /lorder(h9),lexces(h9),llaar(h9),lrecip(h9)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)

      double precision dppp,d2gx,sdzdp
      common/ cxt28 /dppp(j3,j3,m1,h9),d2gx(j3,j3),sdzdp(j3,m11,m10,h9)

      double precision v,tr,pr,r,ps
      common / cst5 /v(l2),tr,pr,r,ps
c----------------------------------------------------------------------
c                                 initialize, d2gx has been set in setw
      g = 0d0
      dg = 0d0
      d2g = d2gx(k,k)
      norder = nord(id)

      if (lexces(id)) then

         do i = 1, jterm(id)
c                                 assuming regular terms
           i1 = jsub(1,i,id)
           i2 = jsub(2,i,id)

           g = g + w(i) * pa(i1) * pa(i2)
           dg = dg + w(i) * (pa(i1)*dydy(i2,k,id) 
     *                     + pa(i2)*dydy(i1,k,id))

         end do  
c                                 get derivative of excess function
         if (llaar(id)) then 

            t = 0d0 
c                                 h&p van laar
            do i = 1, nstot(id)
               t = t + alpha(i)* pa(i)
            end do 
c                                 coming out of this loop g, dg, and
c                                 d2g  are not the complete functions 
c                                 because of the "tphi" term in the 
c                                 van laar.

c                                 convert dg and d2g to the full derivative
            dg = (dg - g*dt(k)/t)/t
            d2g = (d2g - 2d0*dt(k)*dg)/t

         end if 

      end if 
c                                 get the configurational entropy derivatives
      call sderi1 (k,id,ds,d2s)

      dg  = dg + deph(k,id)  - v(2)*ds
      d2g = d2g - v(2)*d2s

      t = omega(id,pa)

      end

      subroutine sderi1 (l,id,ds,d2s)
c----------------------------------------------------------------------
c subroutine to the derivative of the configurational entropy of a 
c solution with respect to the proportion of the lth ordered species.

c THIS DOES NOT INCLUDE ENDMEMBER CONFIGURATION ENTROPY DERIVATIVES!
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,j,k,l,id

      double precision zt,dzdy,dzy,dzyy,zl,ds,d2s,zlnz,dzdyt,dsinf
c                                 working arrays
      double precision zz, pa, p0a, x, w, y
      common/ cxt7 /zz(m4),y(m4),pa(m4),p0a(m4),x(mst,msp),w(m1)
c                                 configurational entropy variables:
      integer msite, ksp, lterm, ksub
      common/ cxt1i /msite(h9),ksp(m10,h9),lterm(m11,m10,h9),
     *               ksub(m0,m11,m10,h9)

      double precision qmult, d0, dcoef, scoef      
      common/ cxt1r /qmult(m10,h9),d0(m11,m10,h9),dcoef(m0,m11,m10,h9),
     *               scoef(m4,h9)

      double precision dppp,d2gx,sdzdp
      common/ cxt28 /dppp(j3,j3,m1,h9),d2gx(j3,j3),sdzdp(j3,m11,m10,h9)

      logical pin
      common/ cyt2 /pin(j3)
c----------------------------------------------------------------------

      ds = 0d0 
      d2s = 0d0 

      do i = 1, msite(id)

         dzy = 0d0  
         dzyy = 0d0
 
         zt = 0d0 
         dzdyt = 0d0 
         dsinf = 0d0 
  
         do j = 1, ksp(i,id)    

            zl = d0(j,i,id)
c                                 for each term:
            do k = 1, lterm(j,i,id)
               zl = zl + dcoef(k,j,i,id) * pa(ksub(k,j,i,id))
            end do 

            if (zl.gt.0d0) then 
               zt = zt + zl 
               zlnz = 1d0 + dlog(zl)
            else
               zlnz = 1d0 
            end if         
c                                 sdzdp is (dz(i,j)/dp(l))
            dzdy = sdzdp(l,j,i,id)

            if (dzdy.eq.0d0) cycle 

            if (zl.gt.0d0) then 
c                                 the first derivative is
               dzy = dzy - dzdy * zlnz
c                                 and the jacobians are

               dzyy = dzyy  - dzdy**2 / zl

            else 
c                                 a species with a non-zero
c                                 derivative is zero, the s
c                                 derivative may be +/-infinite
               dsinf = dsinf + dsign(1d0,dzdy)
                             
            end if 

         end do 
c                                 add the contibution from the ksp(i,id)+1th
c                                 species:
         zl = 1d0 - zt

         if (zl.gt.0d0) then 
            zlnz = 1d0 + dlog(zl)
         else
            zlnz = 1d0 
         end if 

         dzdy = sdzdp(l,j,i,id)

         if (dzdy.ne.0d0) then 
            if (zl.gt.0d0) then 
c                                 the first derivative is
               dzy = dzy - dzdy * zlnz
c                                 and the second is 
               dzyy = dzyy  - dzdy**2 / zl

            else 
c                                 
c                                 a species with a non-zero
c                                 derivative is zero, the s
c                                 derivative may be +/-infinite
               dsinf = dsinf + dsign(1d0,dzdy) 
               
            end if 

         end if 

         if (dabs(dsinf).lt.1d-5) then 
            ds = ds + qmult(i,id)*dzy
            d2s = d2s + qmult(i,id)*dzyy
         else 
            ds = ds + qmult(i,id)*dsinf*1d8
            d2s = d2s - qmult(i,id)*dabs(dsinf)*1d8
         end if 

      end do 

      end

      subroutine p0limt (id)
c----------------------------------------------------------------------
c subroutine to compute the sums of the p0 terms in ordered species 
c limit expressions for solution id.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,j,k,id
c                                 working arrays
      double precision z, pa, p0a, x, w, y
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)

      integer ln,lt,lid,jt,jid
      double precision lc, l0c, jc
      common/ cxt29 /lc(j6,j5,j3,h9),l0c(2,j5,j3,h9),lid(j6,j5,j3,h9),
     *               ln(j3,h9),lt(j5,j3,h9),jc(j3,j5,j3,h9),
     *               jid(j3,j5,j3,h9),jt(j5,j3,h9)

      double precision tsum
      common/ cxt31 /tsum(j5,j3)
c----------------------------------------------------------------------
      do k = 1, nord(id) 
c                                 for ordered species k
         do i = 1, ln(k,id)
c                                 for limit i
            tsum(i,k) = l0c(1,i,k,id)  

            do j = 1, lt(i,k,id)
c                                 for term j
               tsum(i,k) = tsum(i,k) + lc(j,i,k,id)*p0a(lid(j,i,k,id))

            end do
 
         end do 

      end do 

      end  

      subroutine plimit (pmn,pmx,k,id)
c----------------------------------------------------------------------
c subroutine to compute minimum and maximum concentration of ordered
c species k in solution id from site fraction constraints, assumes the
c p0 terms have been accumulated in tsum (routine p0limt)
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,j,k,id

      double precision pmn,pmx,mini
c                                 working arrays
      double precision z, pa, p0a, x, w, y
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)

      integer ln,lt,lid,jt,jid
      double precision lc, l0c, jc
      common/ cxt29 /lc(j6,j5,j3,h9),l0c(2,j5,j3,h9),lid(j6,j5,j3,h9),
     *               ln(j3,h9),lt(j5,j3,h9),jc(j3,j5,j3,h9),
     *               jid(j3,j5,j3,h9),jt(j5,j3,h9)

      double precision tsum
      common/ cxt31 /tsum(j5,j3)
c----------------------------------------------------------------------
      pmx = 1d99
      pmn = -1d99

      do i = 1, ln(k,id)

         mini =  tsum(i,k) 

         do j = 1, jt(i,k,id)

            mini = mini + jc(j,i,k,id)*pa(jid(j,i,k,id))

         end do 

         if (mini.gt.pmn) pmn = mini
         if (l0c(2,i,k,id)+mini.lt.pmx) pmx = mini + l0c(2,i,k,id)

      end do 

      end  

      subroutine readlm (tname,bad,vertex)
c---------------------------------------------------------------------
c readlm - reads stoichiometric limits on ordered species concentrations
c---------------------------------------------------------------------
      implicit none
  
      include 'perplex_parameters.h'

      integer j,k,l,jd,len,inds(k7),ier,ict 

      double precision coeffs(k7)

      logical bad,vertex

      character chars(240)*1, begin*5, tag*3, tname*10

      integer limn,limt,limid,jimid,jimt
      double precision limc,jimc
      common/ cxt30 /limc(j6+2,j5,j3),limid(m0,j5,j3),jimid(j3,j5,j3),
     *               limn(j3),limt(j5,j3),jimc(j3,j5,j3),jimt(j5,j3)

      double precision wg,xmn,xmx,xnc
      integer iend,isub,imd,insp,ist,isp,isite,iterm,iord,istot,jstot,
     *        kstot
      common/ cst108 /wg(m1,m3),xmn(mst,msp),xmx(mst,msp),xnc(mst,msp),
     *      iend(m4),isub(m1,m2,2),imd(ms1,mst),insp(m4),ist(mst),
     *      isp(mst),isite,iterm,iord,istot,jstot,kstot

      integer iddeps,norder 
      double precision depvnu,denth
      common/ cst141 /depvnu(2,j3),denth(j3),iddeps(2,j3),norder
c----------------------------------------------------------------------
c                               initialize limit counter
      do j = 1, norder
         limn(j) = 0
      end do 

      if (norder.eq.1) return

      call readcd (n9,len,ier,chars)

      write (begin,'(5a)') (chars(j),j=1, 5)

      if (begin.eq.'begin') then 

         do 
c                                 read the limit equations for the 
c                                 amount of the ordered endmembers
            call readz (coeffs,inds,ict,istot+norder,tname,tag)

            if (tag.eq.'end') then 
               exit 
            else if (ict.eq.1) then
               bad = .true.
               exit
            end if 
c                                 convert the endmember index to the 
c                                 ordered species index
            jd = inds(1) - istot
            limn(jd) = limn(jd) + 1
            k = limn(jd) 
            if (k.gt.j5) call error (999,coeffs(1),k,tname)
                   
            limt(k,jd) = 0
            jimt(k,jd) = 0 

            do l = 2, ict
c                                 four cases
               if (    (inds(l).le.istot) 
     *             .or.(inds(l) - istot.eq.jd)
     *             .or.(inds(l).eq.inds(l-1).and.
     *                  inds(l) - istot.ne.jd)) then
c                                 1) disordered species go in p0 array
c                                 2) the species is limit species, the 
c                                    species must be in the p0 array.
c                                 3) the species is an ordered species
c                                    and it's the second occurence
                  limt(k,jd) = limt(k,jd) + 1
                  j = limt(k,jd)
                  if (j.gt.j6) call error (33,coeffs(1),j,tname)
                  limid(j,k,jd) = inds(l)
                  limc(j,k,jd) = coeffs(l)

               else if (inds(l).gt.istot) then 
c                                 4) is an ordered species p-term
                  jimt(k,jd) = jimt(k,jd) + 1
                  j = jimt(k,jd)
                  if (j.gt.j3) call error (33,coeffs(1),j,tname)
                  jimid(j,k,jd) = inds(l)
                  jimc(j,k,jd) = coeffs(l)

               else 

                  bad = .true.

               end if 

            end do 
c                                 the constant and delta (max-min) are:
            j = limt(k,jd) + 1
            limc(j  ,k,jd) = coeffs(1)
            limc(j+1,k,jd) = coeffs(ict+1)

         end do 

      else 

         bad = .true.

      end if 

      if (bad) then 
         if (vertex) then 
            write (*,1000) tname,(chars(j),j=1,len)
            write (*,1010)
         end if  
         backspace (n9)
      end if 

1000  format ('**warning ver203** READLM missing or invalid format for '
     *       ,'stoichiometric limit of ordered species',/,'currently ',
     *        'reading (and rejecting) solution model: ',a,
     *      /,'last record was:',/,240a1)
1010  format (/,'This error may indicate you are using an out-of-date '
     *         ,'solution model file.',/
     *         ,'The current version is: '
     *         ,'www.perplex.ethz.ch/perplex/datafiles/solut_08.dat',/)

      end 

      subroutine input9 (vertex,first,output)
c-----------------------------------------------------------------------
c given a list of solution phase names (fname(h9)) input9 searches a
c data file (on unit n9) for the relevant data and subdivides the
c solutions into pseudo-compounds.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer istot,jstot,iend,isub,imd,jmsol,kdsol,icoct,imyn,
     *        ist,i,j,iord,iterm,isite,isp,insp,kstot,ier,
     *        h,im,icky,id,ipoint,icpct,jsmod,idsol,ixct

      logical output, first, bad, vertex
 
      character*10 tname, uname(2)*8, sname(h9), fname, names*8, new*3,
     *             mname*8

      double precision xnc,xmx,xmn,wg,zt,vlaar

      common/ cst108 /wg(m1,m3),xmn(mst,msp),xmx(mst,msp),xnc(mst,msp),
     *      iend(m4),isub(m1,m2,2),imd(ms1,mst),insp(m4),ist(mst),
     *      isp(mst),isite,iterm,iord,istot,jstot,kstot
     *      / cst60 /ipoint,imyn/ csta7 /fname(h9)
     *      / cst18a /mname(m4)
     *      / cst8  /names(k1)
     *      / cst221 /vlaar(m3,m4),jsmod
     *      / cst142 /jmsol(m4,mst),kdsol(m4)

      integer ntot,npairs
      double precision y,xy
      common/ cst86 /xy(mdim,k21),y(k21,mst,ms1),ntot,npairs

      integer iorig,jnsp,iy2p
      common / cst159 /iorig(m4),jnsp(m4),iy2p(m4)

      integer isoct
      common/ cst79 /isoct

      integer io3,io4,io9
      common / cst41 /io3,io4,io9

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  

      integer ixp,ifp
      double precision sxs,exces
      common/ cst304 /sxs(k13),exces(m3,k1),ixp(k1),ifp(k1)

      integer jend
      common/ cxt23 /jend(h9,k12)

      integer ikp
      common/ cst61 /ikp(k1)

      integer inew
      common/ cst140 /inew

      integer ksmod, ksite, kmsol, knsp
      common/ cxt0  /ksmod(h9),ksite(h9),kmsol(h9,m4,mst),knsp(m4,h9)

      double precision pa, p0a, xx, w, yy, z
      common/ cxt7 /yy(m4),xx(m4),pa(m4),p0a(m4),z(mst,msp),w(m1)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)

      save uname

      data uname/' ',' '/
c-----------------------------------------------------------------------
c                                 initialize counters
      ixct = 0 
c                                 gloabl coordinate counter for xcoor (cxt10)
      icoct = 0  
c                                 initialize model counter
      im = 0
c                                 no request for solutions
      if (io9.eq.1) then 
         isoct = 0 
         return 
      end if 
c                                 open pseudocompund list file
      open (99,file="freenergy3carb.dat")

      if (output) open (n8,file="pseudocompound_glossary.dat")
c                                 format test line
      read (n9,1030) new

      if (new.eq.'new'.or.new.eq.'NEW') then
         inew = 0 
      else if (new.eq.'007'.or.new.eq.'008') then
         inew = 1
      else 
         call error (3,xmx(1,1),inew,new)
      end if

      do 
c                                 -------------------------------------
c                                 read the solution name
         call rmodel (tname,bad,vertex)

         if (bad) cycle 
c                                 istot is zero, if eof: 
         if (istot.eq.0) then 
c                                 then at least one solution phase referenced
c                                 in the input is not present in the
c                                 solution phase data file, write warning:
            if (vertex) call warn (43,zt,isoct-im,'INPUT9')
            exit
         end if 
c                                 -------------------------------------
c                                 check the solution model:
         call cmodel (im,idsol,tname,0,uname(1),uname(2),first)

         if (jstot.lt.2) cycle 
c                                 -------------------------------------
c                                 reformulate the model so that it has 
c                                 no missing endmembers:             
         if (jstot.lt.istot) call reform (tname,im,first)

         if (istot.lt.2) cycle  
c                                 -------------------------------------
c                                 make various book keeping arrays (y2p,
c                                 jmsol, dydz, .....)
         call nmodel
c                                 save solution name
         sname(im) = tname
c           
c                                 save found solutions in global solution 
c                                 model arrays
         call gmodel (im,tname)
c                                 generate pseudocompound compositions.
c                                 subdiv returns the total
c                                 number of pseudocompounds (ipcps) and 
c                                 array y, of which element y(h,i,j) is
c                                 the site fraction of the jth species on
c                                 the ith site of the hth pseudocompound.
         if (vertex) then  
c                                 global pseudo-cpd counter for sxs
            icpct = 0           
      
            do 
c                                 load the composition into
c                                 a the site fraction array:
               read (99,*,iostat=ier) (z(1,j),j=1,3)

               if (ier.ne.0) exit 
c                               generate the pseudocompound:
               call soload (im,icoct,icpct,ixct,tname,icky,im)

            end do 

            if (icpct.gt.0) then 
 
               write (*,1100) icpct, tname

               if (output) then

                  if (jsmod.ne.6.and.jsmod.ne.8) then
                     write (n8,1060) tname, 
     *                         (names(jend(im,2+i)), i =1, lstot(im))
                  else if (jsmod.eq.6) then 
                     write (n8,1060) tname, 
     *                         (mname(iorig(knsp(i,im))), 
     *                         i = 1, lstot(im)),'* see footnote 1'
                  else if (jsmod.eq.8) then 
                     write (n8,1060) tname, 
     *                         (mname(iorig(knsp(i,im))), 
     *                         i = 1, lstot(im)),'* see footnotes 1 & 2'
                  end if

                  do i = iphct-icpct+1, iphct
                     write (n8,1070) names(i),
     *                            (sxs(ixp(i)+j), j = 1, lstot(im))
                  end do

                  if (jsmod.eq.6) then 
                     write (n8,1120)
                  else if (jsmod.eq.8) then
                     write (n8,1120)
                     write (n8,1130)
                  end if 
               end if 

            end if 

            jend(im,2) = icpct

         end if 

         if (im.eq.isoct) exit 
c                               read next solution
      end do 

c                               scan for "killed endmembers"
      do i = 1, ipoint
c                               reset ikp
         if (ikp(i).lt.0) ikp(i) = 0
      end do 

      if (vertex) write (*,1110) iphct - ipoint

      if (isoct.gt.0) then 
         if (io3.eq.0.and.output.and.vertex) then 
            write (n3,1020)
            write (n3,1010) (fname(i), i = 1, isoct)
            if (im.gt.0) then 
               write (n3,1000)
               write (n3,1010) (sname(i), i = 1, im)
            else if (output) then 
               write (n3,1040) 
            end if 
         end if 

         do i = 1, im
            fname(i) = sname(i)
         end do 
      end if 

      isoct = im
c                              close pseudocompound list
      if (output) close (n8)
c                              close solution model file
      close (n9)

      first = .false.

1000  format (/,'the following solution models will be considered:',/)
1010  format (7(2x,a10))
1020  format (/,'Of the requested solution models:',/)
1030  format (a)
1040  format (/,'no models will be considered.',/)
1060  format (/,'Solution: ',a,/,12x,'Endmember fractions:',
     *        /,12x,20(a,1x))
1070  format (a,2x,20(1x,f6.3,2x))
1100  format (i8,' pseudocompounds generated for: ',a)
1110  format (/,'Total number of pseudocompounds:',i8)
1120  format (/,'1 - Although the bulk composition of pseudocompounds'
     *        ,' for this solution is fixed,',/,' the proportions of'
     *        ,' its endmembers may vary due to respeciation.',/)
1130  format (/,'2 - Proportions output here may sum to <1 ', 
     *          'because the ordered species',/,'may have non-zero ',
     *          'initial proportions.',/)
      end 

      subroutine subdiv (tname,ids)
c---------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'
 
      character*10 tname

      integer last,i,j,np1,h,index,ids

      integer ntot,npairs
      double precision y,xy
      common/ cst86 /xy(mdim,k21),y(k21,mst,ms1),ntot,npairs

      double precision wg,xmn,xmx,xnc
      integer iend,isub,imd,insp,ist,isp,isite,iterm,iord,istot,jstot,
     *        kstot
      common/ cst108 /wg(m1,m3),xmn(mst,msp),xmx(mst,msp),xnc(mst,msp),
     *      iend(m4),isub(m1,m2,2),imd(ms1,mst),insp(m4),ist(mst),
     *      isp(mst),isite,iterm,iord,istot,jstot,kstot
c---------------------------------------------------------------------
c                                 do the first site:
      call cartes (1,tname,ids)

      last = isp(1) - 1
      if (isp(1).eq.1) last = 1

      do h = 1, npairs
         do i = 1, last
            y(h,1,i) = xy(i,h)
         end do 
      end do 

      np1 = npairs
      ntot = np1 

      if (isite.eq.1) goto 99
c                                 do the second site:
      call cartes (2,tname,ids)

      index = 2 
      last = isp(index) - 1
      if (isp(1).eq.1) last = 1
c                                 there will be a total of
c                                 (npairs-1)*np1 compositions,
c                                 copy the first site 2 distribution
c                                 into the first np1 compositions:
      do h = 1, np1 
         do i = 1, last
c                                 this could be an invalid compostion for
c                                 a 3 site model.
            y(h,2,i) = xy(i,1)
         end do 
      end do 

      do h = 2, npairs
c                                 for each site 2 composition,
c                                 duplicate the range of site 1
c                                 compositions.
         do i = 1, np1

            ntot = ntot + 1

            if (ntot.gt.k21) call error (41,xy(1,1),k21,'SUBDIV_2')

            do j = 1, isp(1) - 1
               y(ntot,1,j) = y(i,1,j)
            end do 
c                                 put in the new site 2 compositions:
            do j = 1, isp(2) - 1 
               y(ntot,2,j) = xy(j,h)
            end do 
         end do 
      end do 
c                                 do the third site:
c                                 this hardwires the array dimensions to "mst"
      if (isite.eq.2) goto 99

      np1 = (npairs-1) * np1

      call cartes (3,tname,ids)
c                                 the use of "index" is necessary in case mst<3.
      index = 3 
      last = isp(index) - 1
      if (isp(1).eq.1) last = 1
c                                 copy the first site 3 distribution
c                                 into the first np1*np2 compositions:
      do h = 1, ntot
         do i = 1, last
            y(h,mst,i) = xy(i,1)
         end do 
      end do 
c                                 for each site 3 composition,
c                                 duplicate the range of site 1 and 
c                                 site 2 compositions.
      do h = 2, npairs

         do i = 1, np1
            ntot = ntot + 1
            do j = 1, isp(1) - 1
               y(ntot,1,j) = y(i,1,j)
            end do 
            do j = 1, isp(2) - 1
               y(ntot,2,j) = y(i,1,j)
            end do 
            do j = 1, isp(mst) - 1
               y(ntot,mst,j) = xy(j,h)
            end do 
         end do 
      end do 
 
99    end

      subroutine soload (isoct,icoct,icpct,ixct,tname,icky,im)
c--------------------------------------------------------------------------
c soload - loads/requires solution properties: 

c   jend(h9,k12)  - h9 is the maximum number of solutions
c                   k12 is the maximum number of endmembers pers
c                   solution plus two.
c   jend(i,1)     - OBSOLETE! is the number of endmembers in solution i.
c   jend(i,2)     - is the number of pseudocompounds of solution i.
c   jend(i,3-3+j) - are the indices of the j endmembers in solution i.
c   sxs(k13)      - contains the mole fractions of the endmembers
c                   in the pseudocompounds.
c   ixp(i)        - a pointer that locates the first mole fraction (- 1)
c                   of the ith pseudocompound, the remaining mole fractions
c                   follow sequentially (as in jend(i,3-3+j)).
c   ikp(i)        - the index of the solution corresponding to pseudocompound i.
c   ifp(i)        - flag, 0  - margules/ideal
c                         1  - internal eos fluid                         
c                         3  - van laar
c                         23 - toop, internal
c   exces(j,i)    - the excess function of pseudocompound i, accounts for 
c                   excess properties and configurational entropy as a function
c                   of pressure and temperature:

c                       gexces(i) = exces(1) + exces(2)*T + exces(3)*P 
c--------------------------------------------------------------------------

      implicit none
  
      include 'perplex_parameters.h'

      character*10 tname, names*8
 
      double precision zpr,hpmelt,gmelt,wg,xmn,xmx,xnc,a0,acoef,smult,
     *                 ctot,vlaar,smix,esum,ctotal,omega

      integer iend,isub,imd,insp,ist,jtic,id,im,
     *        isp,isite,iterm,iord,istot,jstot,kstot,nsite,nspm1,nterm,
     *        nsub1,nsub2,nttyp,jsmod,jmsol,kdsol,
     *        h,i,j,l,m,icpct,isoct,ixct,icky,index,icoct

      common/ cst108 /wg(m1,m3),xmn(mst,msp),xmx(mst,msp),xnc(mst,msp),
     *     iend(m4),isub(m1,m2,2),imd(ms1,mst),insp(m4),ist(mst),
     *     isp(mst),isite,iterm,iord,istot,jstot,kstot
     *      / cst107 /a0(m10,m11),acoef(m10,m11,m0),smult(m10),
     *     nsite,nspm1(m10),nterm(m10,m11),
     *     nsub1(m10,m11,m0,m12),nsub2(m10,m11,m0,m12),nttyp(m10,m11,m0)
     *      / cst8 /names(k1)
     *      / cst3 /ctot(k1)
     *      / cst221 /vlaar(m3,m4),jsmod
     *      / cst142 /jmsol(m4,mst),kdsol(m4)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)

      integer iddeps,norder 
      double precision depvnu,denth
      common/ cst141 /depvnu(2,j3),denth(j3),iddeps(2,j3),norder

      double precision qmult, d0, dcoef, scoef      
      common/ cxt1r /qmult(m10,h9),d0(m11,m10,h9),dcoef(m0,m11,m10,h9),
     *               scoef(m4,h9)

      integer jndq, jdqf, iq
      double precision dqfg, dq 
      common/ cxt9 /dqfg(m3,m4,h9),dq(m4),jndq(m4,h9),jdqf(h9),iq(m4)

      double precision cp
      common/ cst12 /cp(k5,k1)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      integer iff,idss,ifug,ifyn,isyn
      common/ cst10  /iff(2),idss(h5),ifug,ifyn,isyn

      integer jterm, jord, jsub
      common/ cxt2i /jterm(h9),jord(h9),jsub(m2,m1,h9)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  

      integer ikp
      common/ cst61 /ikp(k1)

      integer ixp,ifp
      double precision sxs,exces
      common/ cst304 /sxs(k13),exces(m3,k1),ixp(k1),ifp(k1)

      logical depend,laar,order,fluid,macro,specil
      common/ cst160 /depend,laar,order,fluid,macro,specil

      double precision pa, p0a, zp, w, y, x, z
      common/ cxt7 /y(m4),zp(m4),pa(m4),p0a(m4),z(mst,msp),w(m1)

      double precision xcoor
      integer icoor
      common/ cxt10 /xcoor(k18),icoor(k1)

      double precision wgl,vlar
      common/ cxt2r /wgl(m3,m1,h9),vlar(m3,m4,h9)

      integer ksmod, ksite, kmsol, knsp
      common/ cxt0  /ksmod(h9),ksite(h9),kmsol(h9,m4,mst),knsp(m4,h9)

      integer ideps,icase
      common/ cxt3i /ideps(2,j3,h9),icase(h9)

      double precision y2pg
      common/ cxt4  /y2pg(m15,m4,h9)
c----------------------------------------------------------------------
c                              eliminate end-member compositions 
      do l = 1, istot
         p0a(l) = 0d0
         y(l) = 1d0
         do m = 1, isite
c                              check for invalid compositions
            call zchk (z(m,jmsol(l,m)),tname)
            y(l) = y(l)*z(m,jmsol(l,m))
         end do
      end do   
c                              move site fractions into array indexed 
c                              only by independent disordered endmembers:
      do i = 1, istot
         zp(i) = y(knsp(i,im))
      end do
c                              the composition is acceptable.
      iphct = iphct + 1
      icpct = icpct + 1 

      if (iphct.gt.k1) call error (41,z(1,1),k1,'SOLOAD')

      ikp(iphct) = isoct
      ixp(iphct) = ixct
      icky = 0 

      do i = 1, m3
         exces(i,iphct) = 0d0
      end do
c                               classify model:
      if (jsmod.eq.0) then
c                               fluid uses internal eos
         ifp(iphct) = 1
c                               don't allow gonzoids to treat fluid
c                               as a solution if the user has also 
c                               specified fluid saturation
         if (ifyn.eq.0) call error (43,r,i,tname)  
      else if ((jsmod.eq.1.or.jsmod.eq.2.or.jsmod.eq.7.or.jsmod.eq.24
     *          .or.jsmod.eq.25).and.(.not.laar)) then
c                               no special treatment in computational
c                               routines.
         ifp(iphct) = 0 
      else if (jsmod.eq.26) then 
c                               andreas's H2O-CO2-NaCl internal EoS
         ifp(iphct) = 26
c                               other special cases (internal solution
c                               models)
      else if (jsmod.ge.2) then
c                               Toop, internal EoS
c                               Van Laar ala HP.
         ifp(iphct) = jsmod    
      else 
c                               Ideal or Margules.
         ifp(iphct) = 0
      end if 
c                                encode a name
      if (isite.eq.2.and.istot.eq.4) then
c                                special case 1, bin-bin reciprocal solution
         write (names(iphct),1020) tname,
     *                             (idint(1d2*z(j,1)), j = 1, 2)
      else if (isite.eq.2.and.istot.eq.6.and.isp(1).eq.3) then
c                                special case 2, tern-bin reciprocal solution
         write (names(iphct),1060) tname,
     *                             (idint(1d2*z(1,j)), j = 1, 2),
     *                              idint(1d2*z(2,1))
      else if (isite.eq.2.and.istot.eq.6.and.isp(1).eq.2) then
c                                special case 3, bin-tern reciprocal solution
         write (names(iphct),1060) tname,
     *                              idint(1d2*z(1,1)),
     *                             (idint(1d2*z(2,j)), j = 1, 2)
      else if (isite.eq.2.and.istot.eq.9) then
c                                special case 4, tern-tern reciprocal solution
         write (names(iphct),1010) (idint(1d2*z(1,j)), j = 1, 2),
     *                             (idint(1d2*z(2,j)), j = 1, 2)

      else if (istot.eq.2) then 
c                                binary solutions
         if (zp(1).ge.0.04d0.and.zp(1).le.0.96d0) then
            write (names(iphct),1030) names(kdsol(knsp(1,im))),
     *                                idint(1d2*zp(1))
         else    
            write (names(iphct),1070) 
     *             names(kdsol(knsp(1,im))), 1d2*zp(1)             
         end if
      else if (istot.eq.3) then 
c                                ternary solutions
         write (names(iphct),1040) (names(kdsol(knsp(j,im))),
     *                                idint(1d2*zp(j)), j = 1, 2)
      else if (istot.eq.4) then 
c                                quaternary solutions
         icky = 1
         write (names(iphct),1060) tname,
     *                    (idint(1d2*zp(j)), j = 1, 3)
      else
c                                all the rest:
         icky = 1

         if (iphct.lt.1000000) then 
            write (names(iphct),1080) tname, iphct
         else if (iphct.lt.10000000) then
            write (names(iphct),1100) tname, iphct
         else
            write (names(iphct),1110) iphct
         end if 

      end if 
c                                get blanks out of name:
      if (istot.lt.4) then 
         call unblnk (names(iphct)) 
      else 
         call reblnk (names(iphct)) 
      end if   
c                                 initialize constants:
      smix = 0d0
      esum = 0d0

      do i = 1, icomp
         cp(i,iphct) = 0d0
      end do 
c                                 load constants:
      ctotal = 0d0
      jtic = 0 
c                                 load xcoors if reciprocal
      if (isite.gt.1) then 

         icoor(iphct) = icoct 

         do i = 1, isite
            do j = 1, isp(i)
               icoct = icoct + 1
               xcoor(icoct) = z(i,j)
            end do 
         end do 

      end if 

      do h = 1, lstot(im)
c                               do not count the mole
c                               fractions of absent endmembers
         id = kdsol(knsp(h,im)) 

         ixct = ixct + 1
         if (ixct.gt.k13) call error (40,y(1),k13,'SOLOAD')

         sxs(ixct) = zp(h)
     
         if (depend) then
c                               convert y's to p's
            do j = 1, ndep(im)
               sxs(ixct) = sxs(ixct) 
     *                   + y2pg(j,h,im) * y(knsp(lstot(im)+j,im))
            end do 

         end if
               
         p0a(h) = sxs(ixct) 

         if (sxs(ixct).ne.0d0) then
c                              composition vector
            do l = 1, icomp
               cp(l,iphct) = cp(l,iphct) + sxs(ixct) * cp(l,id)
               if (l.le.icp) ctotal = ctotal + sxs(ixct) * cp(l,id)
            end do 
c                              accumulate endmember configurational entropy
            esum = esum + sxs(ixct) * scoef(h,im)

            jtic = jtic + 1

         end if  

      end do  

      if (order.and.depend) then 

         do i = 1, norder 
c                              compute the fraction of the i'th ordered species
c                              required by the decomposition of the dependent 
c                              disordered species:
            h = lstot(im) + i
            p0a(h) = 0d0
c                               convert y's to p's
            do j = 1, ndep(im)
               p0a(h) = p0a(h) 
     *                + y2pg(j,h,im) * y(knsp(lstot(im)+j,im))
            end do 

            ixct = ixct + 1
            if (ixct.gt.k13) call error (40,y(1),k13,'SOLOAD')
            sxs(ixct) = p0a(h)
c                              split these fraction into the fractions of the
c                              consituent disordered species:
            do j = 1, 2

               x = depvnu(j,i)*p0a(h)
               id = kdsol(knsp(ideps(j,i,im),im))
c                              composition vector
               do l = 1, icomp
                  cp(l,iphct) = cp(l,iphct) + x * cp(l,id)
                  if (l.le.icp) ctotal = ctotal + x * cp(l,id)
               end do 

            end do 

         end do  

      end if 

      do l = 1, icomp
         if (cp(l,iphct).gt.-1d-8.and.cp(l,iphct).lt.0d0) then 
            cp(l,iphct) = 0d0 
         else if (cp(l,iphct).lt.0d0) then
            call error (228,cp(l,iphct),l,tname)
         end if 
      end do 
c                              check if the phase consists
c                              entirely of saturated components:
      if (ctotal.eq.0d0) then
         call satsrt 
c                                 to prevent nan's in the compositional coordinates:
         ctot(iphct) = 1d0
      else
         ctot(iphct) = ctotal
      end if 
c                              compute ideal configurational negentropy:
      if (jsmod.eq.1.and.nsite.eq.0.and.ist(1).ne.0) then  

         do i = 1, isite 
            do j = 1, isp(i)
               if (z(i,j).ne.0d0) smix = smix + 
     *                            dfloat(ist(i))*z(i,j)*dlog(z(i,j)) 
            end do 
         end do     

      else if (order) then
c                                 for cpd formation models, configurational entropy
c                                 is evaluated from speciation.
         smix = 0d0

      else if (depend) then 

         smix = -omega(im,p0a)

      else if (nsite.ne.0) then
 
         smix = -omega(im,zp)

      else if (jsmod.eq.24) then 
c                              hp melt model, use internal routine to get entropy
         smix = -hpmelt(im)

      else if (jsmod.eq.25) then 
c                              ghiorso melt model, use internal routine to get entropy
         smix = -gmelt(im)

      end if 
c                              save it:
      exces(2,iphct) = smix 
c                              load excess terms, if not Laar or ordered:
      if ((.not.laar).and.(.not.order)) then 
         do i = 1, iterm
               zpr = 1d0

               do j = 1, iord
                  if (jsmod.eq.1) then 
c                                 use bw site fraction expansion, redundant
c                                 index means only works for single site models
                     if (isub(i,j,1).ne.0) zpr = zpr 
     *                            * z(jsub(j,i,im),jsub(j,i,im))
                  else  
c                                 use independent endmember expansion
                     if (isub(i,j,1).ne.0) zpr = 
     *                                     zpr * p0a(jsub(j,i,im))
                  end if 
               end do  

               do j = 1, m3
                  exces(j,iphct) = exces(j,iphct) + zpr * wgl(j,i,im)
               end do 

         end do  
      end if 
c                              dqf corrections are also be saved in the
c                              exces array this implies that speciation
c                              does not effect the amount of the dqf'd
c                              endmembers.

c                              p0dord converts the p0 of any ordered species
c                              to it's disordered equivalents, as necessary
c                              for the dqf.
      if (depend.and.order) call p0dord (im)

      do i = 1, jdqf(im)
c                              index points to the endmember in the full
c                              model:
         index = jndq(i,im)
c                              user has made a dqf to an ordered species
c                              or a dependent endmember
         if (kdsol(knsp(index,im)).lt.0) 
     *                        call error (227,exces(1,1),index,tname)

         if (depend) then
            do j = 1, m3
               exces(j,iphct) = exces(j,iphct) + p0a(index)*dqfg(j,i,im)
            end do 
         else if (jsmod.ne.1) then
            do j = 1, m3
               exces(j,iphct) = exces(j,iphct) + y(index)*dqfg(j,i,im)
            end do 
         end if 
      end do

1010  format (i2,i2,i2,i2)
1020  format (a2,i2,'_',i2)
1030  format (a2,i2)
1040  format (a1,i2,a1,i2)
1060  format (a2,i2,i2,i2)
1070  format (a2,f4.1)
1080  format (a2,i6)
1100  format (a1,i7)
1110  format (i8)

999   end

      subroutine satsrt 
c---------------------------------------------------------------------
c routine to sort pseudocompounds consisting entirely of saturated
c components.
c---------------------------------------------------------------------
      implicit none
  
      include 'perplex_parameters.h'

      integer j,idc
 
      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  

      double precision cp
      common/ cst12 /cp(k5,k1)

      integer ids,isct,icp1,isat,io2
      common/ cst40 /ids(h5,h6),isct(h5),icp1,isat,io2
 
      do j = isat, 1, -1
         idc = icp + j
         if (cp(idc,iphct).ne.0d0) then
            isct(j) = isct(j) + 1
            if (isct(j).gt.h6) call error (17,cp(1,1),h6,'SATSRT')
            if (iphct.gt.k1) call error (180,cp(1,1),k1,'SATSRT')
            ids(j,isct(j)) = iphct
            exit
         end if
      end do 

      end