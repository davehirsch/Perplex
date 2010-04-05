      subroutine resub (id,ids,iref,iter)
c----------------------------------------------------------------------
c subroutine to generate new pseudocompound compositions around the
c pseudocompound id of solution ids in iteration iter. ifst is the 
c pointer to the first pseudocompound of the solution ids and ilst
c points to the last. 
c----------------------------------------------------------------------
      implicit none 

      include 'perplex_parameters.h'
c                                 -------------------------------------
c                                 local variables
      logical good

      double precision xxnc, ysum

      integer i, j, k, ids, id, iter, kcoct, iref
c                                 -------------------------------------
c                                 functions
      double precision gsol1, ydinc
c                                 -------------------------------------
c                                 global variables:
c                                 adaptive g and compositions
      integer jphct
      double precision g2, cp2
      common/ cxt12 /g2(k21),cp2(k5,k21),jphct
c                                 adaptive z coordinates
      integer jcoct, jcoor, jkp
      double precision zcoor
      common/ cxt13 /zcoor(k20),jcoor(k21),jkp(k21),jcoct
c                                 working arrays
      double precision z, pa, p0a, x, w, y
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1)
c                                 x coordinate description
      integer istg, ispg, imlt, imdg
      double precision xmng, xmxg, xncg, xmno, xmxo
      common/ cxt6r /xmng(h9,mst,msp),xmxg(h9,mst,msp),xncg(h9,mst,msp),
     *               xmno(h9,mst,msp),xmxo(h9,mst,msp)
      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)
c                                 temporary subdivision limits:
      double precision wg,xmn,xmx,xnc
      integer iend,isub,imd,insp,ist,isp,isite,iterm,iord,istot,
     *        jstot,kstot

      common/ cst108 /wg(m1,m3),xmn(mst,msp),xmx(mst,msp),xnc(mst,msp),
     *      iend(m4),isub(m1,m2,2),imd(ms1,mst),insp(m4),ist(mst),
     *      isp(mst),isite,iterm,iord,istot,jstot,kstot
c                                 coordinates output by subdiv
      double precision xy,yy
      integer ntot,npairs
      common/ cst86 /xy(mdim,k21),yy(k21,mst,ms1),ntot,npairs
c                                 max number of refinements for solution h9
      integer ncoor,maxitg
      common/ cxt24 /ncoor(h9),maxitg(h9)
c                                 option values
      integer iopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10)
c----------------------------------------------------------------------
      
      if (iter.eq.1) then
c                                first iteration id array points to 
c                                original compound arrays:
         call getolx (ids,id)

      else
c                                on subsequent iterations get the y's
c                                stored in the ycoor array by routine 
c                                saver, these are reindexed copies of the
c                                coordinates originally saved in zcoor
c                                below. 
         call getxy0 (ids,id)

      end if

      if (iter.le.maxitg(ids).and.iter.le.iopt(10)) then 
c                                load the subdivision limits into
c                                temporary limit arrays:
         isite = istg(ids)
      
         do i = 1, isite

            isp(i) = ispg(ids,i)

            do j = 1, isp(i) - 1

               imd(j,i) = imdg(j,i,ids)

               xxnc = nopt(14)*xncg(ids,i,j)/iopt(11)**(iter-1)

               if (imd(j,i).eq.0) then 
c                                 cartesian
                  xmn(i,j) = x(i,j) - xxnc
                  xmx(i,j) = x(i,j) + xxnc

               else
c                                 conformal
                  xmn(i,j) = ydinc (x(i,j),-xxnc,imd(j,i),j,i,ids)
                  xmx(i,j) = ydinc (x(i,j),xxnc,imd(j,i),j,i,ids)

               end if 

               xnc(i,j) = 2d0*nopt(14)*xncg(ids,i,j)/iopt(11)**iter
               if (xmn(i,j).lt.xmng(ids,i,j)) xmn(i,j) = xmng(ids,i,j)
               if (xmx(i,j).gt.xmxg(ids,i,j)) xmx(i,j) = xmxg(ids,i,j)

            end do 
         end do 
                            
         call subdiv ('characters',ids) 

         do i = 1, ntot 

            jphct = jphct + 1
            if (jphct.gt.k21) call error (58,x(1,1),k21,'resub')
c                                 convert to compositional corrdinates 
c                                 required by routine gsol, y coordinates
c                                 are placed in first array of cxt7,
c                                 store a copy of x coordinates in 
c                                 1-d array zcoor
            jkp(jphct) = ids
            jcoor(jphct) = jcoct - 1
            kcoct = jcoct + ncoor(ids)
c                                 counter for number of non 0 or 1 compositions
            good = .false.

            if (kcoct.gt.k20) call error (59,x(1,1),k20,'resub')

            do j = 1, isite
               ysum = 0d0
               do k = 1, isp(j) - 1
                  x(j,k) = yy(i,j,k)
                  ysum = ysum + x(j,k)
                  zcoor(jcoct) = x(j,k)
                  if (x(j,k).gt.0d0.and.x(j,k).lt.1d0) good = .true.
                  jcoct = jcoct + 1
               end do 
               x(j,isp(j)) = 1d0 - ysum
               zcoor(jcoct) = x(j,isp(j))
               jcoct = jcoct + 1
            end do 

            if (.not.good) then 
               jphct = jphct - 1
               jcoct = jcoct - ncoor(ids)
               cycle
            end if 

            call xtoy (ids)
c                                 call gsol to get g of the solution, gsol also
c                                 computes the p compositional coordinates
            g2(jphct) = gsol1(ids)
c                                 use the coordinates to compute the composition 
c                                 of the solution
            call csol (ids)

            iref = iref + 1
     
         end do 

      else
c                                  simply load the coordinate of the single composition
         jphct = jphct + 1

         jkp(jphct) = ids

         jcoor(jphct) = jcoct - 1

         kcoct = jcoct + ncoor(ids)

         do j = 1, istg(ids)
            do k = 1, ispg(ids,j)
               zcoor(jcoct) = x(j,k)
               jcoct = jcoct + 1
            end do 
         end do 

         call xtoy (ids)
c                                 call gsol to get g of the solution, gsol also
c                                 computes the p compositional coordinates
         g2(jphct) = gsol1(ids)
c                                 use the coordinates to compute the composition 
c                                 of the solution
         call csol (ids)

      end if 

      end 

      subroutine csol (id)
c-----------------------------------------------------------------------
c csol computes chemical composition of solution id from the macroscopic
c endmember fraction array y or p0a (cxt7), these arrays are prepared by a prior
c call to function gsol. the composition is loaded into the array cp2 at
c position jphct.
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i,j,id

      double precision ctot2
c                                 -------------------------------------
c                                 global variables:
      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      double precision cp
      common/ cst12 /cp(k5,k1)

      double precision ctot
      common/ cst3  /ctot(k1)
c                                 adaptive coordinates
      integer jphct
      double precision g2, cp2
      common/ cxt12 /g2(k21),cp2(k5,k21),jphct
c                                 bookkeeping variables
      integer ksmod, ksite, kmsol, knsp
      common/ cxt0  /ksmod(h9),ksite(h9),kmsol(h9,m4,mst),knsp(m4,h9)
c                                 working arrays
      double precision z, pa, p0a, x, w, y
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1)

      integer jend
      common/ cxt23 /jend(h9,k12)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)
c----------------------------------------------------------------------

      ctot2 = 0d0

      do i = 1, icp
         cp2(i,jphct) = 0d0
      end do  

      if (ksmod(id).lt.6.or.ksmod(id).gt.9) then 
c                                 general case (y coordinates)
         do i = 1, mstot(id)

            do j = 1, icp 
               cp2(j,jphct) = cp2(j,jphct) + y(i) * cp(j,jend(id,2+i))
            end do

            ctot2 = ctot2 + y(i)*ctot(jend(id,2+i)) 

         end do 

      else 
c                                 cpd-formation (ksmod = 6) &/or reciprocal 
c                                 solutions with dependent endmembers, p0a 
c                                 contains the p's. for ksmod=8 these are a 
c                                 reformulation of the p's to eliminate the ordered 
c                                 endmembers. p0a is constructed in function gsol.
         do i = 1, lstot(id) 
            do j = 1, icp 
               cp2(j,jphct) = cp2(j,jphct) + p0a(i) * cp(j,jend(id,2+i))
            end do 
            ctot2 = ctot2 + p0a(i)*ctot(jend(id,2+i))
         end do  
         
      end if 
c                                  normalize the composition and free energy
      g2(jphct) = g2(jphct)/ctot2

      do j = 1, icp 
         cp2(j,jphct) = cp2(j,jphct)/ctot2
      end do  

      end

      subroutine reopt (idead,jdv,npt,js)
c-----------------------------------------------------------------------
c reopt - given the results of an initial optimization for lpopt, reopt
c iteratively refines the solution by generating pseudocompounds in the
c neighborhood of the initial optimization.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer liw, lw, iter, iref, i, j, k, npt, 
     *        id, idead, ids, jstart, inc

      parameter (liw=2*k21+3,lw=2*(k5+1)**2+7*k21+5*k5)  

      double precision  ax(k5), x(k21), clamda(k21+k5), w(lw)

      integer is(k21+k5), iw(liw), jdv(k19), js(k19)

c                                 why save anything? augbug
c     save ax, x, is, iw
c                                 -------------------------------------
c                                 global variables
c                                 adaptive coordinates
      integer jphct
      double precision g2, cp2
      common/ cxt12 /g2(k21),cp2(k5,k21),jphct

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      double precision xa,b,xc
      common/ cst313 /xa(k5,k1),b(k5),xc(k1)

      integer ikp
      common/ cst61 /ikp(k1)

      double precision g
      common/ cst2  /g(k1)

      double precision cp
      common/ cst12 /cp(k5,k1)

      double precision ctot
      common/ cst3  /ctot(k1)
c                                 adaptive x(i,j) coordinates
      integer jcoct, jcoor, jkp
      double precision zcoor
      common/ cxt13 /zcoor(k20),jcoor(k21),jkp(k21),jcoct

      integer iopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10)

      integer lcoor,lkp
      double precision ycoor
      common/ cxt14 /ycoor(k22),lcoor(k19),lkp(k19)

      integer ipoint,imyn
      common/ cst60 /ipoint,imyn
c-----------------------------------------------------------------------
c                                 the pseudocompounds to be refined
c                                 are identified in jdv(1..npt)
      iter = 1
      jphct = 0
      iref = 0 
      jcoct = 1
      inc = istct - 1
c                                 --------------------------------------
c                                 first iteration
      do i = 1, npt
         id = jdv(i) + inc
         if (id.le.ipoint) then
c                                 the point is a true compound
            jphct = jphct + 1
            jkp(jphct) = -id
            g2(jphct) = g(id)/ctot(id)

            do j = 1, icp
               cp2(j,jphct) = cp(j,id)/ctot(id)
            end do 

         else 
c                                 the point is a pseudocompound 
c                                 to be refined
            call resub (id,ikp(id),iref,1)

         end if
c                                 reset jdv to point to the new
c                                 indexing, this is only for the
c                                 case of no iteration (i.e., resub
c                                 increments jphct by one.
         jdv(i) = jphct

      end do 

      if (iref.eq.0) goto 90
 
      do while (iref.gt.0.and.iter.le.iopt(10))
c                                 nothing has been refined

         iter = iter + 1
c                                 cold start
         jstart = 0 
c                                 set idead = 0 to prevent lpnag from
c                                 overwriting warm start parameters
         idead = 0 
c                                 do the optimization
         call lpnag (jphct,icp,cp2,k5,b,g2,is,x,ax,
     *               clamda,iw,liw,w,lw,idead,l6,jstart)
c                                 warn if severe error
         if (idead.gt.0) then 
            call lpwarn (idead,'REOPT ')
            goto 99
         end if 
c                                 analyze solution, get refinement points
         call yclos2 (clamda,x,is,jdv,js,npt,iter,idead)

         if (idead.gt.0) goto 99
c                                 save the id and compositions
c                                 of the refinement points, this
c                                 is necessary because resub rewrites
c                                 the xcoor array.
         call saver (jdv,npt)

         if (iter.gt.iopt(10)) exit 

         jphct = 0 
         iref = 0 
         jcoct = 1
c                                 generate new pseudocompounds
         do i = 1, npt

            ids = lkp(i)

            if (ids.lt.0) then 
c                                 the point is a true compound
               jphct = jphct + 1
               jkp(jphct) = ids
               ids = -ids
               g2(jphct) = g(ids)/ctot(ids)

               do j = 1, icp
                  cp2(j,jphct) = cp(j,ids)/ctot(ids)
               end do 

            else

c              write (*,*) 'iteration',iter
c                                 the point is a pseudocompound 
               call resub (i,ids,iref,iter)

            end if
c                                 reset jdv in case we're going to exit
            jdv(i) = i 

         end do 
      
      end do 

90    jphct = 0

      do k = 1, npt
         if (js(k).ne.1) then
             jphct = jphct + 1
             jdv(jphct) = jdv(k)
          end if 
      end do

      npt = jphct

      if (npt.gt.icp) idead = 1
c                                 get chemical potentials, this is
c                                 just for output.
      call getmus (jdv,2,idead)       

99    end 

      subroutine sortin (ind,k,n)
c-----------------------------------------------------------------------
c sort the first k values of ind
c-----------------------------------------------------------------------
      implicit none

      integer i, j, k, n, ind(n), imin

      do j = 1, k-1

         imin = ind(j)

         do i = j+1, k

            if (ind(i).lt.imin) then 
               imin = ind(i)
               ind(i) = ind(j)
               ind(j) = imin
            end if
 
         end do 

      end do 

      end 

      subroutine lpwarn (idead,char)
c----------------------------------------------------------------------
c write warning messages from lpnag as called by routine 'char',
c set flags ier and idead, the optimization is a total fail if
c idead set to 1.
c----------------------------------------------------------------------
      implicit none

      integer idead

      character*6 char     

      double precision c
c----------------------------------------------------------------------
c                                             look for errors                                            
      if (idead.eq.2.or.idead.gt.4) then 
c                                             unbounded solution, or
c                                             other programming error.
         call warn (91,c,idead,char) 

      else if (idead.eq.3) then 
c                                             no feasible solution
         call warn (42,c,idead,char)

      else if (idead.eq.4) then 
c                                             iteration count exceeded,
c                                             probable cause no feasible
c                                             solution.
         call warn (90,c,idead,char) 

      end if

      end 

      subroutine saver (jdv,npt)
c----------------------------------------------------------------------
c subroutine to save a copy of adaptive pseudocompound x(i,j) compositions
c in the temporary array ycoor (also lcoor) used by resub to generate
c the new zcoor array for the subsequent iteration.  
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'
c                                 -------------------------------------
c                                 local variables
      integer jdv(k19), npt, i, j, k, kcoct, id, ids, itic
c                                 -------------------------------------
c                                 global variables:
c                                 adaptive z coordinates
      integer jcoct, jcoor, jkp
      double precision zcoor
      common/ cxt13 /zcoor(k20),jcoor(k21),jkp(k21),jcoct
c                                 interim storage array
      integer lcoor,lkp
      double precision ycoor
      common/ cxt14 /ycoor(k22),lcoor(k19),lkp(k19)
c                                 x coordinate description
      integer istg, ispg, imlt, imdg
      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)
c----------------------------------------------------------------------
      kcoct = 0

      do i = 1, npt

         id = jdv(i)
         ids = jkp(id)
         lkp(i) = ids
c                                 cycle on a compound
         if (ids.lt.0) cycle
c                                 it's a solution:
         lcoor(i) = kcoct
         itic = 0

         do j = 1, istg(ids)
            do k = 1, ispg(ids,j)
               itic = itic + 1
               if (kcoct+itic.gt.k22) 
     *             call error (60,ycoor(1),k22,'saver')
               ycoor(lcoor(i)+itic) = zcoor(jcoor(id)+itic)
            end do 
         end do 

         kcoct = kcoct + itic

      end do 

      end 

      subroutine getxz (jd,id,ids)
c----------------------------------------------------------------------
c subroutine to recover geometric reciprocal solution compositions (x(i,j))
c from the zcoor array loaded in resub.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, jd, id, ids, icoor
c                                 working arrays
      double precision z, pa, p0a, x, w, y
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1)
c                                 x coordinate description
      integer istg, ispg, imlt, imdg
      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)
c                                 adaptive x(i,j) coordinates
      integer jcoct, jcoor, jkp
      double precision zcoor
      common/ cxt13 /zcoor(k20),jcoor(k21),jkp(k21),jcoct
c                                  xcoordinates for the final solution
      double precision x3
      common/ cxt16 /x3(k21,mst,msp)
c----------------------------------------------------------------------
      icoor = jcoor(id)

      do i = 1, istg(ids)
         do j = 1, ispg(ids,i)
            icoor = icoor + 1
            x(i,j) = zcoor(icoor)
            x3(jd,i,j) = x(i,j)
         end do 
      end do 

      end 

      subroutine getxy0 (ids,id)
c----------------------------------------------------------------------
c subroutine to recover geometric reciprocal solution compositions (x(i,j))
c from the ycoor array loaded in saver
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
c                                 interim storage array
      integer lcoor,lkp
      double precision ycoor
      common/ cxt14 /ycoor(k22),lcoor(k19),lkp(k19)
c----------------------------------------------------------------------
      jcoor = lcoor(id)

      do i = 1, istg(ids)
         do j = 1, ispg(ids,i)
            jcoor = jcoor + 1
            x(i,j) = ycoor(jcoor)
         end do 
      end do 

      end 

      subroutine rebulk (jdv,idead)
c----------------------------------------------------------------------
c rebulk computes the amounts of the stable compounds and eliminates 
c those with zero modes (if requested).
c----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer jdv(k19),idead,i,j,k,icomp,istct,iphct,icp,iff,idss,
     *        ifug,ifyn,isyn,ipvt,npt

      double precision a,b,cp
 
      common/ cst6   /icomp,istct,iphct,icp
     *      / cst10  /iff(2),idss(h5),ifug,ifyn,isyn
     *      / cst301 /a(k5,k5),b(k5),ipvt(k5)
     *      / cst12 /cp(k5,k1)

      integer jbulk
      double precision cblk
      common/ cst300 /cblk(k5),jbulk
c                                 adaptive x(i,j) coordinates
      integer jcoct, jcoor, jkp
      double precision zcoor
      common/ cxt13 /zcoor(k20),jcoor(k21),jkp(k21),jcoct
c                                 compositions of stable adaptive
c                                 coordinates (and solution ids).
      integer kkp,np,ncpd,ntot
      double precision cp3,ctot3
      common/ cxt15 /cp3(k5,k5),ctot3(k5),kkp(k5),np,ncpd,ntot
c                                 options from perplex_option.dat
      integer iopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10)

      integer ikp
      common/ cst61 /ikp(k1)
c----------------------------------------------------------------------
c                                 load the transpose of the
c                                 concentration matrix of the pseudo-
c                                 invariant assemblage.
      npt = icp 

      do i = 1, jbulk
         if (i.le.icp) then
c                                 getcmp assigns kkp, so jkp is uneccesary
c                                 from here on
            call getcmp (i,jdv(i),jkp(jdv(i)))
            do j = 1, jbulk
               a(j,i) = cp3(j,i)
            end do
         else
            npt = npt + 1
            k = idss(i-icp)
c                                  set case for solution in saturated component
c                                  space, the endmember composition is not set,
c                                  this is gonna cause problems, at least for 
c                                  meemum
            if (ikp(k).eq.0) then 
               kkp(i) = -k
            else 
               kkp(i) = ikp(k)
            end if 

            do j = 1, jbulk
               a(j,i) = cp(j,k)
            end do
         end if
      end do  
c                                 factor the matrix
      call factr1 (jbulk,idead)
      if (idead.eq.1) goto 99 
c                                 test for bounded compositions
      do i = 1, jbulk
c                                 load composition vector, the alpha
c                                 vector is returned in the same array
         b(i) = cblk(i)
      end do
c                                 solve for the alpha vector
      call subst1 (jbulk)
c                                 check for and eliminate zero mode
c                                 phases
      ntot = jbulk

      if (nopt(9).gt.0d0) call zmode 

99    end

      logical function solvs1 (id1,id2)
c-----------------------------------------------------------------------
c function to test if a solvus separates two pseudocompounds .
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      double precision dx

      integer icomp, iphct, icp, i, id1, id2, istct

      common/ cst6 /icomp,istct,iphct,icp

      integer iopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10)
c                                 composition and model flags
c                                 for final adaptive solution
      integer kkp,np,ncpd,ntot
      double precision cp3, ctot3
      common/ cxt15 /cp3(k5,k5),ctot3(k5),kkp(k5),np,ncpd,ntot
c-----------------------------------------------------------------------

      dx = 0d0

      do i = 1, icp
         dx = dx + (cp3(i,id1)/ctot3(id1) - cp3(i,id2)/ctot3(id2))**2
      end do 

      if (dsqrt(dx).gt.nopt(8)) then
         solvs1 = .true.
      else
         solvs1 = .false.
      end if 

      end 

      logical function solvs2 (id1,id2)
c-----------------------------------------------------------------------
c function to test if a solvus separates two pseudocompounds .
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      double precision dx

      integer icomp, iphct, icp, i, id1, id2, istct

      common/ cst6 /icomp,istct,iphct,icp

      integer iopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10)

      integer jphct
      double precision g2, cp2
      common/ cxt12 /g2(k21),cp2(k5,k21),jphct
c-----------------------------------------------------------------------

      dx = 0d0

      do i = 1, icp
         dx = dx + (cp2(i,id1) - cp2(i,id2))**2
      end do 

      if (dsqrt(dx).gt.nopt(8)) then
         solvs2 = .true.
      else
         solvs2 = .false.
      end if 

      end 

      subroutine avrger 
c----------------------------------------------------------------------
c avrger combines discretization points into a single solution
c composition. on output

c     np  - is the number of solutions, 
c     ncpd - is the number of true compounds
c     ntot - np+ncpd

c this routine is unecessarily complicated, because it assumes
c pseudocompounds are not ordered by solution (but they are
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical solvs1, check
c                                 -------------------------------------
c                                 local variables
      integer idsol(k5),kdsol(k5,k5),ids,isite,xidsol,xkdsol,irep,
     *        i,j,jdsol(k5,k5),jd,k,l,nkp(k5),xjdsol(k5)

      double precision bsol(k5,k5),cpnew(k5,k5),xx,ct,xb(k5), 
     *                 bnew(k5),xnew(k21,mst,msp)
c                                 -------------------------------------
c                                 global variables:
c                                 x coordinate description
      integer istg, ispg, imlt, imdg
      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      integer ipvt
      double precision a,b
      common/ cst301 /a(k5,k5),b(k5),ipvt(k5)
c                                 composition and model flags
c                                 for final adaptive solution
      integer kkp, np, ncpd, ntot
      double precision cp3, ctot3
      common/ cxt15 /cp3(k5,k5),ctot3(k5),kkp(k5),np,ncpd,ntot
c                                  x-coordinates for the final solution
      double precision x3
      common/ cxt16 /x3(k21,mst,msp)

      integer ikp
      common/ cst61 /ikp(k1)

      double precision cp
      common/ cst12 /cp(k5,k1)
c                                 working arrays
      double precision z, pa, p0a, x, w, y
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1)
c-----------------------------------------------------------------------
c                                first check if solution endmembers are
c                                among the stable compounds:
      do i = 1, ntot
         if (kkp(i).lt.0) then 
            if (ikp(-kkp(i)).ne.0) then 
c                                we have an endmember
               nkp(i) = ikp(-kkp(i))
            else
               nkp(i) = kkp(i)
            end if 
         else 
            nkp(i) = kkp(i)
         end if 
      end do 
c                                check if any solutions
      do i = 1, ntot
         if (nkp(i).gt.0) goto 10
      end do 

      np = 0
      ncpd = ntot
c                                the assemblage consists entirely 
c                                of true compounds, nothing to do, but load
c                                the cp3 array
      do i = 1, ncpd
         l = -kkp(i)
         do j = 1, icomp
            cp3(j,i) = cp(j,l)
         end do 
      end do 

      goto 99
c                                figure out how many solutions
c                                are present:
10    np = 0
      ncpd = 0

      do 30 i = 1, ntot
         if (nkp(i).lt.0) then
c                                 the pseudocompound is a true compound
            ncpd = ncpd + 1 
            idsol(ntot) = ncpd
            bsol(ntot,ncpd) = b(i)
            kdsol(ntot,ncpd) = nkp(i)       
            jdsol(ntot,ncpd) = i   
         else 
            do j = 1, np
c                                 compare the compound to the np solutions 
c                                 identfied so far:        
               if (kdsol(j,1).eq.nkp(i)) then 
c                                 if match check for a solvus
                  if (.not.(solvs1(i,jdsol(j,idsol(j))))) then
c                                 the pseudocompound matches a solution
c                                 found earlier.
                     idsol(j) = idsol(j) + 1
                     bsol(j,idsol(j)) = b(i)
                     jdsol(j,idsol(j)) = i  
                     goto 30 
                  end if 
               end if 
            end do
c                                 the pseudocompound is a new solution 
c                                 phase.
            np = np + 1
            idsol(np) = 1
            kdsol(np,1) = nkp(i)
            jdsol(np,1) = i 
            bsol(np,1) = b(i)

         end if    
30    continue  
c                                 check if a solution occurs more than once
c                                 but the occurences are not sequential (this
c                                 can only occur if an endmember is immiscible 
c                                 with a general composition
      if (np.gt.2) then
 
         do i = 1, np

            check = .false.
            irep = 0

            do j = i+1, np
               if (kdsol(j,1).ne.kdsol(i,1)) then
                  check = .true.
               else 
                  irep = irep + 1
               end if 
            end do 

            if (check.and.irep.gt.0) then 

               l = i + 1

               if (kdsol(l,1).ne.kdsol(i,1)) then 
c                                 not in sequence, find the next occurence
                  do j = i+2, np 
                     if (kdsol(i,1).eq.kdsol(j,1)) exit
                  end do 
c                                 swap phase at i+1 with the one at j
                  xidsol = idsol(l)
                  xkdsol = kdsol(l,1)
                  do k = 1, xidsol
                     xb(k) = bsol(l,k)
                     xjdsol(k) = jdsol(l,k)
                  end do 

                  idsol(l) = idsol(j)
                  kdsol(l,1) = kdsol(j,1)
                  do k = 1, idsol(j)
                     bsol(l,k) = bsol(j,k)
                     jdsol(l,k) = jdsol(j,k)
                  end do 

                  idsol(j) = xidsol
                  kdsol(j,1) = xkdsol
                  do k = 1, xidsol
                     bsol(j,k) = xb(k)
                     jdsol(j,k) = xjdsol(k)
                  end do 

               end if 
            end if 
         end do 
      end if 
c                                 if a solution is represented by
c                                 more than one pseudocompound get
c                                 the everage composition
      do i = 1, np 
c                                 initialize
         bnew(i) = 0d0

         do j = 1, icomp
            cpnew(j,i) = 0d0
         end do 

         ids = kdsol(i,1)
         isite = istg(ids)

         do j = 1, isite
            do k = 1, ispg(ids,j)
               xnew(i,j,k) = 0d0
            end do 
         end do 

         do j = 1, idsol(i)
            bnew(i) = bnew(i) + b(jdsol(i,j))
         end do 

         do j = 1, idsol(i)

            jd = jdsol(i,j)
c                                conditional in case zero mode
c                                is off:
            if (bnew(i).gt.0d0) then 

               xx =  b(jd)/bnew(i)
c                                save the new compositions
               do k = 1, icomp
                  cpnew(k,i) = cpnew(k,i) + xx*cp3(k,jd)
               end do 

               do k = 1, isite
                  do l = 1, ispg(ids,k)
                     xnew(i,k,l) = xnew(i,k,l) + xx*x3(jd,k,l)
                  end do 
               end do 
            
            else 
c                               
               do k = 1, icomp
                  cpnew(k,i) = cp3(k,jd)
               end do 

               do k = 1, isite
                  do l = 1, ispg(ids,k)
                     xnew(i,k,l) = x3(jd,k,l)
                  end do 
               end do 

            end if 

         end do 

      end do
c                                now reform the arrays kdv and b
      do i = 1, ncpd
         k = np + i
         l = kdsol(ntot,i)
         b(k) = bsol(ntot,i)
         kkp(k) = l
c                                for the sake of completeness load
c                                compound composition into cp3 array
         do j = 1, icomp
            cp3(j,k) = cp(j,-l)
         end do 
      end do 

      do i = 1, np

         b(i) = bnew(i)
         kkp(i) = kdsol(i,1)
         ids = kkp(i)

         ct = 0d0 

         do j = 1, icomp
            ct = ct + cpnew(j,i)
            cp3(j,i) = cpnew(j,i)
         end do

         ctot3(i) = ct

         do j = 1, istg(ids)
            do k = 1, ispg(ids,j)
c                                 set x's for sollim
               x(j,k) = xnew(i,j,k)
c                                 set x's for global storage
               x3(i,j,k) = x(j,k) 

            end do 
         end do 
c                                 check composition against solution model ranges
c                                 if auto_refine is on:
         call sollim (ids)

      end do

      ntot = np + ncpd

99    end 

      subroutine sollim (ids)
c----------------------------------------------------------------------
c subroutine to extract compositional range of endmembers in stable phases
c for auto_refine option.
c----------------------------------------------------------------------
      implicit none 

      include 'perplex_parameters.h'

      integer ids,i,j

      double precision z, pa, p0a, x, w, y
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1)
c                                 x coordinate description
      integer istg, ispg, imlt, imdg
      double precision xmng, xmxg, xncg, xmno, xmxo
      common/ cxt6r /xmng(h9,mst,msp),xmxg(h9,mst,msp),xncg(h9,mst,msp),
     *               xmno(h9,mst,msp),xmxo(h9,mst,msp)
      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)
c                                 solution limits and stability
      logical stable,limit
      double precision xlo,xhi
      common/ cxt11 /xlo(m4,mst,h9),xhi(m4,mst,h9),stable(h9),limit(h9)

      character fname*10
      common/ csta7 /fname(h9)
c                                 interval limits conformal transformation
      integer intv
      double precision yint, yfrc
      common/ cst47 /yint(5,ms1,mst,h9),yfrc(4,ms1,mst,h9),intv(4)

      integer iopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10)
c----------------------------------------------------------------------
c                                 set stable flag
      stable(ids) = .true.
c                                 check x-ranges
      do i = 1, istg(ids)
         do j = 1, ispg(ids,i) - 1
c                                 low limit:
            if (x(i,j).lt.xlo(j,i,ids)) then

               xlo(j,i,ids) = x(i,j)
c                                 check if solution is at an unnatural limit
               if (x(i,j).gt.xmno(ids,i,j).and.
     *             x(i,j).le.xmng(ids,i,j)) then
c                                 relax limits according to subdivsion model
                  if (imdg(j,i,ids).eq.0) then 
c                                 cartesian
                     xmng(ids,i,j) = xmng(ids,i,j) - nopt(10)
                     if (xmng(ids,i,j).lt.0d0) xmng(ids,i,j) = 0d0
                     write (*,1000) fname(ids),x(i,j),i,j,xmng(ids,i,j)

                  else if (imdg(j,i,ids).eq.1.or.imdg(j,i,ids).eq.4)then 
c                                 assymmetric stretching towards xmin
                     yint(1,j,i,ids) = yint(1,j,i,ids) - nopt(10)
                     if (yint(1,j,i,ids).lt.0d0) yint(1,j,i,ids) = 0d0
                     write (*,1000) fname(ids),x(i,j),i,j,
     *                              yint(1,j,i,ids)
                     xmng(ids,i,j) =  yint(1,j,i,ids)

                  else 
c                                 symmetric modes, don't reset
                     write (*,1010) fname(ids),x(i,j),i,j,imdg(j,i,ids)
c                                 set xmn to prevent future warnings
                     xmng(ids,i,j) = 0d0 
                  end if 

                  limit(ids) = .true.

               end if 
            end if 
c                                 high limit:
            if (x(i,j).gt.xhi(j,i,ids)) then
               xhi(j,i,ids) = x(i,j)
c                                 check if solution is at an unnatural limit
               if (x(i,j).lt.xmxo(ids,i,j).and.
     *             x(i,j).ge.xmxg(ids,i,j)) then
c                                 relax limits according to subdivsion model
                  if (imdg(j,i,ids).eq.0) then 
c                                 cartesian
                     xmxg(ids,i,j) = xmxg(ids,i,j) + nopt(10)
                     if (xmxg(ids,i,j).gt.1d0) xmxg(ids,i,j) = 1d0
                     write (*,1000) fname(ids),x(i,j),i,j,xmxg(ids,i,j)

                  else if (imdg(j,i,ids).eq.1.or.imdg(j,i,ids).eq.4)then 
c                                 assymmetric stretching
                     yint(2,j,i,ids) = yint(2,j,i,ids) + nopt(10)
                     if (yint(2,j,i,ids).gt.1d0) yint(2,j,i,ids) = 1d0
                     write (*,1000) fname(ids),x(i,j),i,j,
     *                              yint(2,j,i,ids)
                     xmxg(ids,i,j) = yint(2,j,i,ids)

                  else 
c                                 symmetric modes, don't reset
                     write (*,1010) fname(ids),x(i,j),i,j,imdg(j,i,ids)
c                                 set xmx to prevent future warnings
                     xmxg(ids,i,j) = 1d0 

                  end if 

                  limit(ids) = .true.

               end if 
            end if 
         end do 
      end do  

1000  format (/,'WARNING: composition of solution ',a,' has reached an',
     *          ' internal limit (',f5.3,')',/,
     *          'on site ',i1,' for species ',i2,'.',/,'If this warning'
     *         ,' occurs during auto-refinement, the problem can be ',
     *          'circumvented',/,'by increasing the auto_refine_slop ',
     *          'specified in perplex_option.dat.',//,
     *          'For the remainder of this calculation the internal ',
     *          'limit has been relaxed to: ',f5.3,/)

1010  format (/,'WARNING: composition of solution ',a,' has reached an',
     *          ' internal limit (',f5.3,')',/,
     *          'on site ',i1,' for species ',i2,'.',/,'If this warning'
     *         ,' occurs during auto-refinement, the problem can be ',
     *          'circumvented',/,'by increasing the auto_refine_slop ',
     *          'specified in perplex_option.dat.',//,
     *          'The limit for this subdivision scheme (imd=',i1,
     *          ') cannot be relaxed automatically.',/)
      end 

      subroutine outlim 
c----------------------------------------------------------------------
c subroutine to extract compositional range of endmembers in stable phases
c for auto_refine option.
c----------------------------------------------------------------------
      implicit none 

      include 'perplex_parameters.h'
c                                 -------------------------------------
c                                 local variables:
      integer i,j,k,ibad1,ibad2,igood

      logical bad1,bad2,good
c                                 -------------------------------------
c                                 global variables:
c                                 working arrays
c                                 x coordinate description
      integer istg, ispg, imlt, imdg
      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)
c                                 solution limits and stability
      logical stable,limit
      double precision xlo,xhi
      common/ cxt11 /xlo(m4,mst,h9),xhi(m4,mst,h9),stable(h9),limit(h9)
c                                 solution model counter
      integer isoct
      common/ cst79 /isoct
c                                 solution model names
      character*10 fname
      common/ csta7 /fname(h9)
c                                 endmember pointers
      integer jend
      common/ cxt23 /jend(h9,k12)
c                                 endmember names
      character names*8
      common/ cst8  /names(k1)

      logical refine
      common/ cxt26 /refine
c----------------------------------------------------------------------
      ibad1 = 0 
      ibad2 = 0 
      igood = 0 
      rewind (n10)
      rewind (n11)

      if (isoct.eq.0) goto 99

      bad1 = .false.
      bad2 = .false.
      good = .false.

      do i = 1, isoct

         if (.not.stable(i)) then
            bad1 = .true.
            ibad1 = ibad1 + 1
         else
            good = .true.
            igood = igood + 1
         end if
 
         if (limit(i)) bad2 = .true.

      end do 

      if (.not.refine) write (n10,*) ibad1,ibad2,igood
c                                 write solutions present that are 
c                                 not stable
      if (bad1) then 

         write (*,1000)
         write (n11,1000)
     
         do i = 1, isoct
            if (.not.stable(i)) then 
               write (*,'(5x,a)') fname(i) 
               if (.not.refine) write (n10,'(a)') fname(i)
               write (n11,'(5x,a)') fname(i) 
            end if 
         end do
      end if 

      if (.not.good) goto 99
c                                 write solutions that are on an internal
c                                 limit
      if (bad2) then 
         write (*,1010) 
         write (n11,1010) 
         do i = 1, isoct
            if (limit(i)) then
               write (*,'(5x,a)') fname(i) 
               write (n11,'(5x,a)') fname(i)
            end if  
         end do
      end if

      do i = 1, isoct
         if (.not.stable(i)) cycle

         if (.not.refine) then

            write (n10,'(a)') fname(i)

            do j = 1, istg(i)
               do k = 1, ispg(i,j)-1
                  write (n10,*) xlo(k,j,i),xhi(k,j,i)
               end do
            end do 

         end if 

         if (istg(i).eq.1) then 
c                                 single site solution
            write (*,1020) fname(i)
            write (n11,1020) fname(i)
            do j = 1, ispg(i,1) - 1
               write (*,1030) names(jend(i,2+j)),xlo(j,1,i),xhi(j,1,i)
               write (n11,1030) names(jend(i,2+j)),xlo(j,1,i),xhi(j,1,i)
            end do 
         else
c                                 reciprocal solution
            write (*,1040) fname(i)
            write (n11,1040) fname(i)

            do j = 1, istg(i)
               
               write (*,1050) j
               write (n11,1050) j

               if (ispg(i,j).eq.1) then 
                  write (*,1060)
                  write (n11,1060) 
               else
                  do k = 1, ispg(i,j) - 1
                     write (*,1070) k,xlo(k,j,i),xhi(k,j,i)
c    *                     ,names(jend(i,2+indx(i,j,k)))
                     write (n11,1070) k,xlo(k,j,i),xhi(k,j,i)
c    *                     ,names(jend(i,2+indx(i,j,k)))
                  end do 
               end if 
            end do
         end if 
      end do 

99    close (n10)
      close (n11)

1000  format (/,'WARNING: The following solutions were input, but are',
     *          ' not stable:',/)
1010  format (/,'WARNING: The following solutions have compositions on',
     *          ' an internal limit (i.e., 0<x<1)',/,'(see ranges ',
     *          'below to determine which limits should be relaxed or',
     *        /,'if executing in auto_refine mode inrease auto_refine',
     *          '_slop in perplex_option.dat):',/)
1020  format (/,'Endmember compositional ranges for model: ',a,//,5x,
     *        'Endmember   Minimum   Maximum')
1030  format (5x,a8,4x,f7.5,3x,f7.5)
1040  format (/,'Site fraction ranges for multisite model: ',a)
1050  format (/,'  Site ',i1,/,5x,'Species   Minimum   Maximum   ')
c     *          'Endmember with this species')
1070  format (8x,i1,6x,f7.5,3x,f7.5,3x,12(a8,1x))
1060  format (8x,'Dummy site generated by model reformulation',/)

      end 

      subroutine sorter (kdbulk,ico,jco,output)
c----------------------------------------------------------------------
c sorter compares assemblages to those already defined and reorders 
c the phases if the assemblage has been identified earlier
c----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i,j,k,l,m,kdbulk,ico,jco,ids,ioct,inct

      logical output 

      double precision cpt(k5,k5),xt(k5,mst,msp),bt(k5),ct3(k5)
c                                 x-coordinates for the final solution
      double precision x3
      common/ cxt16 /x3(k21,mst,msp)

      integer iap,ibulk
      common/ cst74  /iap(k2),ibulk

      integer idasls,iavar,iasct,ias
      common/ cst75  /idasls(k5,k3),iavar(3,k3),iasct,ias
c                                 composition and model flags
c                                 for final adaptive solution
      integer kkp, np, ncpd, ntot
      double precision cp3, ctot3
      common/ cxt15 /cp3(k5,k5),ctot3(k5),kkp(k5),np,ncpd,ntot

      integer ipvt
      double precision a,b
      common/ cst301 /a(k5,k5),b(k5),ipvt(k5)
c                                 x coordinate description
      integer istg, ispg, imlt, imdg
      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp
c----------------------------------------------------------------------
c                                 look for a match with known assemblages
      do 110 i = 1, iasct

         if (np.ne.iavar(1,i).or.ncpd.ne.iavar(2,i)) cycle 

         do 120 j = 1, ntot
            do k = 1, ntot 
               if (idasls(k,i).eq.kkp(j)) then 
c                                 check that the phase occurs the same 
c                                 number of times in each assemblage:
                   inct = 0 
                   ioct = 0 
                   do l = 1, ntot
                      if (kkp(l).eq.kkp(j)) inct = inct + 1
                      if (idasls(l,i).eq.kkp(j)) ioct = ioct + 1
                   end do 

                   if (ioct.ne.inct) goto 110 

                   goto 120

               end if 
            end do
c                                 no match with compound j:
c                                 do next assemblage
            goto 110 

120      continue  

         if (ibulk.gt.k2) call error (183,0d0,k2,'SORTER')
         ibulk = ibulk + 1
         iap(ibulk) = i
         kdbulk = ibulk
c                                 reorder the result arrays of the
c                                 current occurence to match initial 
c                                 occurence:
         do j = 1, ntot

            do k = 1, ntot

               if (kkp(k).eq.idasls(j,i)) then
c                                 load temporary array
                  bt(j) = b(k)

                  if (kkp(k).gt.0) then 

                     do l = 1, icomp
                        cpt(l,j) = cp3(l,k)
                     end do

                     ct3(j) = ctot3(k)

                     do l = 1, istg(kkp(k))
                        do m = 1, ispg(kkp(k),l)
                           xt(j,l,m) = x3(k,l,m) 
                        end do 
                     end do 
                  end if 
c                                 this eliminates immiscible phases
                  kkp(k) = 0

                  exit 
 
               end if 

            end do 

         end do
c                                 reload final arrays from temporary
         do j = 1, ntot

            b(j) = bt(j)
            ids = idasls(j,i)
            kkp(j) = ids

            if (ids.gt.0) then 

               do k = 1, icomp
                  cp3(k,j) = cpt(k,j)
               end do

               ctot3(j) = ct3(j)

               do k = 1, istg(ids)
                  do l = 1, ispg(ids,k)
                     x3(j,k,l) = xt(j,k,l) 
                  end do 
               end do 
            end if 
         end do 

         goto 98 

110   continue 
c                                 the assemblage is new:
      iasct = iasct + 1
      if (iasct.gt.k3) call error (184,0d0,k3,'BLKMAT')

      do i = 1, ntot
         idasls(i,iasct) = kkp(i)
      end do

      ibulk = ibulk + 1
      if (ibulk.gt.k2) call error (183,0d0,k2,'BLKMAT')
      kdbulk = ibulk 
      iap(ibulk) = iasct 

      iavar(1,iasct) = np
      iavar(2,iasct) = ncpd
      iavar(3,iasct) = np + ncpd
c                                
98    if (output) call outbl1 (ico,jco)

     
      end 

      subroutine outbl1 (ico,jco)
c----------------------------------------------------------------------
c output data for compositions and phases of assemblage ibulk
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ico,jco,i,j,k,ids
c                                 -------------------------------------
c                                 global variables
      integer iap,ibulk
      common/ cst74  /iap(k2),ibulk
c                                 composition and model flags
c                                 for final adaptive solution
      integer kkp, np, ncpd, ntot
      double precision cp3, ctot3
      common/ cxt15 /cp3(k5,k5),ctot3(k5),kkp(k5),np,ncpd,ntot
c                                 x-coordinates for the final solution
      double precision x3
      common/ cxt16 /x3(k21,mst,msp)
c                                 x coordinate description
      integer istg, ispg, imlt, imdg
      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)
c                                 molar amounts (b)
      integer ipvt
      double precision a,b
      common/ cst301 /a(k5,k5),b(k5),ipvt(k5)
c                                 i/o
      integer io3,io4,io9
      common / cst41 /io3,io4,io9

      double precision mu
      common/ cst330 /mu(k8)

      integer jtest,jpot
      common/ debug /jtest,jpot

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp
c----------------------------------------------------------------------
      if (io4.eq.1) return
c                                graphics output  
      write (n5,1000) ico,jco,iap(ibulk)
c                                phase molar amounts
      write (n5,1010) (b(i),i=1,np+ncpd)
c                                solution phase compositions
      do i = 1, np
         ids = kkp(i)
         write (n5,1010) ((x3(i,j,k),k=1,ispg(ids,j)),j=1,istg(ids))
      end do 
c                                dependent potentials
      if (jpot.ne.1) write (n5,1010) (mu(i),i=1,icp)

1000  format (10(i8,1x))
1010  format (6(g16.8,1x))

      end 

      subroutine yclos1 (clamda,x,is,jphct,jdv,js,npt,idead)
c----------------------------------------------------------------------
c subroutine to identify pseudocompounds close to the solution for 
c subsequent refinement for iteration 1.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer jphct, jdv(k19), js(k19), npt, i, j, is(k1+k5), imin,
     *        idsol(k5), kdv(k1+k5), nsol, mpt, iam, id, lpt,idead,
     *        left, right, inc, jdsol(k5,k5), kdsol(k5), max

      external ffirst

      logical smart, solvus

      double precision clamda(k1+k5), tol, cmin, dlamda(k1+k5),
     *                 cmax, dlam, x(k1)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      integer iopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10)

      integer ikp
      common/ cst61 /ikp(k1)
c----------------------------------------------------------------------

      npt = 0 
      nsol = 0 
      cmin = 1d99
      cmax = 0d0
      imin = 0 
      inc = istct - 1

      do i = 1, jphct

         if (is(i).eq.0.or.is(i).eq.2) then 
c                                 make a list of found phases:
            id = i + inc

            if (ikp(id).ne.0) then 
               do j = 1, nsol
                  if (ikp(id).eq.idsol(j)) then 
                     kdsol(j) = kdsol(j) + 1
                     jdsol(kdsol(j),j) = id
                     goto 10
                  end if 
               end do
c                                 new phase, add to list
               nsol = nsol + 1
               idsol(nsol) = ikp(id)
               jdsol(1,nsol) = id
               kdsol(nsol) = 1

            end if 
c                                 new point, add to list
10          npt = npt + 1
            jdv(npt) = i

         else 
c                                 find the nearest (and furthest)
c                                 phase
            dlam = dabs(clamda(i))

            if (dlam.lt.cmin) then  
               cmin = dlam
               imin = i
            else if (dlam.gt.cmax) then 
               cmax = dlam
            end if 

          end if 

      end do

      if (npt.gt.icp) then 
c                                 look through for zero phases
         do i = 1, npt
            if (x(jdv(i)).eq.0d0) idead = idead + 1
         end do

         if (idead.ge.npt-icp) then 

            idead = npt - icp 
            mpt = 0
c                                 eliminate zero phase
            do i = 1, npt 

               if (x(jdv(i)).eq.0d0.and.idead.gt.0) then 
                  idead = idead - 1
c                                 reset the active flag
c                 is(jdv(i)) = 1
                  cycle
               end if 

               mpt = mpt + 1
               jdv(mpt) = jdv(i)

            end do 

            npt = mpt 

         else 

            idead = 1
            goto 999

         end if 

      else if (npt.lt.icp) then 

         idead = 1 
         goto 999 

      end if 

      if (iopt(12).eq.1) then 
c                                 only saving the closest metastable
c                                 point
         if (imin.ne.0) then 
            npt = npt + 1
            jdv(npt) = imin
         end if 
         goto 999
      end if 


c     nsol = 0 
c                                 look for solvi and turn off any
c                                 phase with a solvus
      mpt = 0 
 
      do 15 i = 1, nsol
         if (kdsol(i).gt.1) then 
            do j = 2, kdsol(i)
               if (solvus(jdsol(1,i),jdsol(j,i))) goto 15 
            end do 
         end if 
c                                 no solvus
         mpt = mpt + 1
         idsol(mpt) = idsol(i)
15    continue 

      nsol = mpt
c                                 make a list of non-active
c                                 constraints that do not match
c                                 the phases of the active list
      mpt = 0 

      do 20 i = 1, jphct

         if (is(i).ne.0.and.is(i).ne.2) then 

            iam = ikp(i + inc)

            if (iam.ne.0) then                                  
               do j = 1, nsol
                  if (iam.eq.idsol(j)) goto 20 
               end do
            end if 

            mpt = mpt + 1
            kdv(mpt) = i
            dlamda(mpt) = dabs(clamda(i))

         end if 

20    continue 
c                                 smart sort:
      smart = .true.

      if (mpt.eq.0) then
c                                 no points other than active
         goto 999

      else if (smart) then 

         left = 1

         right = mpt

         max = iopt(12)

         if (mpt.lt.max) max = mpt

         call ffirst (dlamda,kdv,left,right,max,k1+k5,ffirst)
 
         do i = 1, max

            npt = npt + 1

            jdv(npt) = kdv(i)

         end do
c                                 sort the phases
         call sortin (jdv,npt,k19)

      else 
c                                 find the closest non-active 
c                                 phases:
         cmax = cmin + cmax * dfloat(iopt(12))/mpt

30       tol = cmin + (cmax-cmin)/2d0

         lpt = 0 
         cmax = 0d0

         do i = 1, mpt

            if (dlamda(i).lt.tol) then 

               lpt = lpt + 1
               jdv(npt+lpt) = kdv(i)
               if (dlamda(i).gt.cmax) cmax = dlamda(i)

               if (lpt.gt.iopt(12)) exit

            end if 

         end do 

         if (lpt.eq.0.and.imin.ne.0) then
c                                 didn't find anything, add
c                                 minimum point and go:
            npt = npt + 1
            jdv(npt) = imin

         else if (cmax-cmin.gt.1d-16.and.lpt.gt.iopt(12)) then 

            goto 30 

         else 

            npt = npt + lpt

         end if 
c                                 sort the phases
         call sortin (jdv,npt,k19)

      end if 

999   do i = 1, npt
         js(i) = is(jdv(i))
      end do 

      end 

      subroutine yclos2 (clamda,x,is,jdv,js,npt,iter,idead)
c----------------------------------------------------------------------
c subroutine to identify pseudocompounds close to the solution for 
c subsequent refinement for iteration 1.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer jdv(k19), js(k19), npt, i, j, is(k21+k5), iter,
     *        idsol(k5), kdv(k21+k5), nsol, mpt, iam, lpt, imin, left,
     *        right,max,idead, jdsol(k5,k5), kdsol(k5)

      logical smart, solvs2

      double precision clamda(k21+k5), tol, cmin, dlamda(k21+k5),
     *                 dlam, cmax, x(k21)

      external ffirst

      integer iopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10)
c                                 adaptive x(i,j) coordinates
      integer jcoct, jcoor, jkp
      double precision zcoor
      common/ cxt13 /zcoor(k20),jcoor(k21),jkp(k21),jcoct

      integer jphct
      double precision g2, cp2
      common/ cxt12 /g2(k21),cp2(k5,k21),jphct

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp
c----------------------------------------------------------------------

      npt = 0 
      nsol = 0 
      cmin = 1d99 
      cmax = 0d0 
      imin = 0 

      do i = 1, jphct

         if (is(i).eq.0.or.is(i).eq.2) then 

            if (jkp(i).gt.0) then 
c                                 make a list of found phases: 
               do j = 1, nsol
                  if (jkp(i).eq.idsol(j)) then 
                     kdsol(j) = kdsol(j) + 1
                     jdsol(kdsol(j),j) = i
                     goto 10
                  end if 
               end do
c                                 new phase, add to list
               nsol = nsol + 1
               idsol(nsol) = jkp(i)
               jdsol(1,nsol) = i
               kdsol(nsol) = 1

            end if 
c                                 new point, add to list
10          npt = npt + 1

            if (npt.gt.icp) then 
               idead = 1
               goto 99 
            end if 

            jdv(npt) = i

          else 

            dlam = dabs(clamda(i))

            if (dlam.lt.cmin) then  
               cmin = dlam
               imin = i
            else if (dlam.gt.cmax) then 
               cmax = dlam
            end if 

          end if 

      end do

      if (npt.gt.icp) then 
c                                 look through for zero phases
         do i = 1, npt
            if (x(jdv(i)).eq.0d0) idead = idead + 1
         end do

         if (idead.ge.npt-icp) then 

            idead = npt - icp 
            mpt = 0
c                                 eliminate zero phase
            do i = 1, npt 

               if (x(jdv(i)).eq.0d0.and.idead.gt.0) then 
                  idead = idead - 1
c                                 reset the active flag
c                 is(jdv(i)) = 1
                  cycle
               end if 

               mpt = mpt + 1
               jdv(mpt) = jdv(i)

            end do 

            npt = mpt 

         else 

            idead = 1
            goto 999

         end if 

      else if (npt.lt.icp) then 

         idead = 1 
         goto 999 

      end if 
 

      if (iter.gt.iopt(10)) goto 99

      if (iopt(12).eq.1) then 
         if (imin.ne.0) then
            npt = npt + 1
            jdv(npt) = imin
         end if 
         goto 99
      end if 
c                                 now look for solvi and turn off any
c                                 phase with a solvus
      mpt = 0 
 
      do 15 i = 1, nsol
         if (kdsol(i).gt.1) then 
            do j = 2, kdsol(i)
               if (solvs2(jdsol(1,i),jdsol(j,i))) goto 15
            end do 
         end if 
c                                 no solvus
         mpt = mpt + 1
         idsol(mpt) = idsol(i)

15    continue 

      nsol = mpt
c                                 make a list of non-active
c                                 constraints that do not match
c                                 the phases of the active list
      mpt = 0

      do 20 i = 1, jphct

         if (is(i).ne.0.and.is(i).ne.2) then 

            iam = jkp(i)

            if (iam.gt.0) then                                   
               do j = 1, nsol
                  if (iam.eq.idsol(j)) goto 20 
               end do
            end if

            mpt = mpt + 1
            kdv(mpt) = i
            dlamda(mpt) = dabs(clamda(i))

         end if 

20    continue 
c                                 smart sort:
      smart = .true.

      if (mpt.eq.0) then
c                                 no points other than active
         goto 999

      else if (iopt(12).eq.1.and.imin.ne.0) then 

         npt = npt + 1
         jdv(npt) = imin

      else if (smart) then 

         left = 1

         right = mpt

         max = iopt(12)

         if (mpt.lt.max) max = mpt

         call ffirst (dlamda,kdv,left,right,max,k1+k5,ffirst)
 
         do i = 1, max

            npt = npt + 1

            jdv(npt) = kdv(i)

         end do 
c                                 sort the phases
         call sortin (jdv,npt,k19)

      else 
c                                 find the closest non-active 
c                                 phases:
         cmax = cmin + cmax * dfloat(iopt(12))/mpt

30       tol = cmin + (cmax-cmin)/2d0

         lpt = 0 
         cmax = 0d0

         do i = 1, mpt

            if (dlamda(i).lt.tol) then 

               lpt = lpt + 1
               jdv(npt+lpt) = kdv(i)
               if (dlamda(i).gt.cmax) cmax = dlamda(i)

               if (lpt.gt.iopt(12)) exit

            end if 

         end do 

         if (lpt.eq.0.and.imin.ne.0) then
c                                 didn't find anything, add
c                                 minimum point and go:
            npt = npt + 1
            jdv(npt) = imin

         else if (cmax-cmin.gt.1d-16.and.lpt.gt.iopt(12)) then 

            goto 30 

         else 

            npt = npt + lpt

         end if 
c                                 sort the phases
         call sortin (jdv,npt,k19)

      end if 

99    do i = 1, npt
         js(i) = 0
      end do 

999   end

      subroutine ffirst (a, ind, left, right, k, n, dumsub)
c-----------------------------------------------------------------------
c find the k smallest values of array between indices left and right
c from http://en.wikipedia.org/wiki/Selection_algorithm
c-----------------------------------------------------------------------
      implicit none

      integer left, right, k, n, pivot, opivot, partit, ind(n)

      external dumsub

      double precision a(n)

      if (right.gt.left) then 

         opivot = left + (right-left)/2
         pivot = partit (a, ind, left, right, opivot, n)

         if (pivot.gt.k) then 
            call dumsub (a,ind,left,pivot-1,k,n,dumsub)
         else if (pivot.lt.k) then 
            call dumsub (a,ind,pivot+1,right,k-pivot,n,dumsub)
         end if 

      end if 

      end 

      integer function partit (a, ind, left, right, opivot, n)

      implicit none

      integer left, right, n, pivot, opivot, iold, ind(n), i

      double precision a(n), value, oldval

      value = a(opivot)
c                                 swap a(opivot) with a(right)
      iold = ind(opivot)
      a(opivot) = a(right)
      ind(opivot) = ind(right)
      a(right) = value
      ind(right) = iold

      pivot  = left

      do i = left, right-1

         if (a(i).le.value) then
 
            iold = ind(pivot)
            oldval = a(pivot)
            a(pivot) = a(i)
            ind(pivot) = ind(i)
            a(i) = oldval
            ind(i) = iold
            
            pivot = pivot + 1
         end if 
      end do 
c                                 swap a(right) with a(pivot)
      iold = ind(pivot)
      oldval = a(pivot)
      a(pivot) = a(right)
      ind(pivot) = ind(right)
      a(right) = oldval 
      ind(right) = iold
      
      partit = pivot

      end 

      subroutine subst (a,ipvt,n,b,ier)
c-----------------------------------------------------------------------
c subst uses the lu decomposition of the matrix 'a' contained
c in the array 'a' to solve ax = b for x. subst is modified from the
c the subroutine of the same name listed by conte and de boor
c in 'elementary numerical analysis', mcgraw-hill, 1980.
c factor uses scaled partial pivoting.
 
c input     a- an n by n array containing the non-zero elements of
c              the u and l decompositions of a, as output by factor.
c           n- the dimension of the matrix a.
c        ipvt- a vector indicating that row ipvt(k) was used to
c              eliminate the coefficient a(n,k).
c           b- the vector b.
c output    b- the solution vector x.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'
 
      double precision a(k8,k8),b(k8),x(k8),sum

      integer ipvt(k8),ip,i,j,n,ii,ier
c----------------------------------------------------------------------
c                                 solve ly = b for y:
      ip = ipvt(1)
      x(1) = b(ip)
      do i = 2, n

         sum = 0d0

         do j = 1, i - 1
            sum = a(i,j)*x(j)+sum
         end do 

         ip = ipvt(i)
         x(i) = b(ip)-sum

      end do 
c                                 solve ux = y for x:
      if (a(n,n).eq.0d0) then
c                                 this check should be superfluous,
c                                 but reopt requires it. should check
c                                 what's with factor. 
         ier = 1
         goto 99
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
            ier = 1
            goto 99
         end if 

         x(i) = (x(i)-sum)/a(i,i)
         b(i) = x(i)

      end do 
      b(n) = x(n)
 
99    end

      subroutine subst1 (n)
c-----------------------------------------------------------------------
c subst uses the lu decomposition of the matrix 'a' contained
c in the array 'a' to solve ax = b for x. subst is modified from the
c the subroutine of the same name listed by conte and de boor
c in 'elementary numerical analysis', mcgraw-hill, 1980.
c factor uses scaled partial pivoting.
 
c input     a- an n by n array containing the non-zero elements of
c              the u and l decompositions of a, as output by factor.
c           n- the dimension of the matrix a.
c        ipvt- a vector indicating that row ipvt(k) was used to
c              eliminate the coefficient a(n,k).
c           b- the vector b.
c output    b- the solution vector x.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision x(k5), sum

      integer n, i, j, im1, ip1, nm1, ii, ip

      integer ipvt
      double precision a,b
      common/ cst301 /a(k5,k5),b(k5),ipvt(k5) 

c                            solve ly = b for y:
      ip = ipvt(1)
      x(1) = b(ip)

      do i = 2, n
         sum = 0d0
         im1 = i-1
         do j = 1, im1
            sum = a(i,j)*x(j)+sum
         end do 
         ip = ipvt(i)
         x(i) = b(ip)-sum
      end do 
c                            solve ux = y for x:
      x(n) = x(n)/a(n,n)
      nm1 = n-1

      do ii = 1, nm1
         i = n-ii
         ip1 = i+1
         sum = 0d0
         do j = ip1, n
            sum = a(i,j)*x(j)+sum
         end do
         x(i) = (x(i)-sum)/a(i,i)
         b(i) = x(i)
      end do 

      b(n) = x(n)
 
      end

      subroutine factr1 (n,ier)
c-----------------------------------------------------------------------
c factr1 is a subroutine which calculates the triangular
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

      integer i,j,k,n,ip1,istr,ier

      double precision temp,ratio,tmax,rmax

      integer ipvt
      double precision a,d
      common/ cst301 /a(k5,k5),d(k5),ipvt(k5)
c-----------------------------------------------------------------------
      ier = 0
c                            initialize ipvt,d
      do i = 1, n
         ipvt(i) = i
         rmax = 0d0
         do j = 1,n
            rmax = dmax1(rmax,dabs(a(i,j)))
         end do 
c                            ax = b is singular if rmax = 0
         if (dabs(rmax).lt.1d-5) goto 9000
         d(i) = rmax
      end do 
c                            begin decomposition:
      do i = 1, n - 1
c                            determine pivot row (istr).
         rmax = dabs(a(i,i))/d(i)
         istr = i
         ip1 = i + 1

         do j = ip1, n
            tmax = dabs(a(j,i))/d(j)
            if (tmax.le.rmax) cycle
            rmax = tmax
            istr = j
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
            do k = ip1, n
               a(j,k) = a(j,k)-ratio*a(i,k)
            end do 
         end do 
 
      end do 
     
      if (dabs(a(n,n)).lt.1d-5) ier = 1

      return
c                           algoritmic singularity.
9000  ier = 1
 
      end

      subroutine ufluid (fo2)
c----------------------------------------------------------------------
c subroutine ufluid computes the potential of the components
c of a saturated fluid phase. if the mole fraction of a component is les
c less than 1.d-38 the chemical potential is set to -9.9d09.
c ufluid may call one of three molecular fluid equations of state, or
c alternatively users may supply their own routines, however,
c the routines currently in use return the log of a components fugacity
c which is then added to the reference state potential computed by the
c function gphase.
c----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i
 
      double precision xf(2),fo2,fs2,gph

      double precision thermo,uf,us 
      common/ cst1 /thermo(k4,k10),uf(2),us(h5)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      double precision f
      common/ cst11 /f(2)

      integer ifct,idfl
      common/ cst208 /ifct,idfl

      integer iff,idss,ifug,ifyn,isyn
      common/ cst10  /iff(2),idss(h5),ifug,ifyn,isyn
c-----------------------------------------------------------------------
c                           compute the chemical potentials of
c                           fluid components in fluid saturated
c                           systems.
      call cfluid (fo2,fs2)

      if (idfl.ne.0) then
         call gphase (idfl,gph)
         uf(idfl) = gph + r * t * f(idfl)
      else
         xf(1) = 1d0 - xco2
         xf(2) = xco2
 
         do i = 1, 2
            if (iff(i).ne.0) then 
               if (xf(i).lt.1d-38) then 
                  uf(i) = -1d10
               else 
                  call gphase (i,gph)
                  uf(i) = gph + r * t * f(i)
               end if
            end if 
         end do

      end if 
 
      end

      subroutine uproj
c----------------------------------------------------------------------
c subroutine uproj computes the potentials of saturated phase components
c and saturated components.
  
c the energies of saturated components are projected through
c saturated volatile components.
c----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'
 
      integer i,j,k,l,ict,ll,i1,id

      double precision uss(h6),fo2,gph,u

      double precision g
      common/ cst2 /g(k1)

      double precision cp
      common/ cst12 /cp(k5,k1)

      double precision thermo,uf,us 
      common/ cst1 /thermo(k4,k10),uf(2),us(h5)

      integer ifct,idfl
      common/ cst208 /ifct,idfl

      integer iff,idss,ifug,ifyn,isyn
      common/ cst10  /iff(2),idss(h5),ifug,ifyn,isyn

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  

      integer ids,isct,icp1,isat,io2
      common/ cst40 /ids(h5,h6),isct(h5),icp1,isat,io2

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps
c-----------------------------------------------------------------------
      fo2 = 0d0
c                                 compute the chemical potentials of
c                                 saturated phase components.
      if (ifyn.ne.1) call ufluid (fo2)

      do i = 1, isat
c                                 determine stable saturated composants
c                                 and the corresponding chemical potentials
         ict = isct(i)

         ll = icp+i

         do j = 1, ict

            k = ids(i,j)
            call gphase (k,gph)
            
            if (ifct.gt.0) then 
               do l = 1, 2
c                                 legendre transform for saturared phase
c                                 component potentials
                  if (iff(l).ne.0) gph = gph - cp(iff(l),k)*uf(l)
               end do 
            end if 

            uss(j) = gph 

            if (i.gt.1) then 
c                                 if multiple component saturation constraints
c                                 apply saturation hierarchy legendre transform:
               i1 = i-1
               do l = 1, i1
                  uss(j) = uss(j)-cp(icp+l,k)*us(l)
               end do
            end if 

            g(k) = uss(j)
            uss(j) = uss(j)/cp(ll,k)
         end do 
c                                 if O2, check if fo2 has been 
c                                 determined by a fluid phase routine,
c                                 if so, add the transform:
         if (io2.eq.i) then 
            do j = 1, ict 
               uss(j) = uss(j) + r*t*fo2
            end do 
         end if 
c                           now find stable "composant":

         u = uss(1)

         id = 1

         if (ict.ne.1) then 
            do j = 2, ict
               if (uss(j).gt.u) cycle  
               id = j
               u = uss(j)
            end do 
         end if 
c                               save the id of the stable composant.
         idss(i) = ids(i,id)
c                               and its chemical potential.
         us(i) = u
c                               in case a phase in the component
c                               saturation space is an endmember of
c                               a solution transform the endmember G's:
         do j = 1, ict
            k = ids(i,j)
            g(k) = g(k) - cp(icp+i,k)*u
         end do 

      end do 

      end

      subroutine initlp 
c--------------------------------------------------------------------
c initialize arrays and constants for lp minimizarion
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,j

      double precision ctotal

      integer hcp,idv
      common/ cst52  /hcp,idv(k7)  

      double precision ctot
      common/ cst3  /ctot(k1)

      double precision cp
      common/ cst12 /cp(k5,k1)

      double precision a,b,c
      common/ cst313 /a(k5,k1),b(k5),c(k1)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  

      integer jphct,istart
      common/ cst111 /jphct,istart

      integer jbulk
      double precision cblk
      common/ cst300 /cblk(k5),jbulk
c                                 stuff used in lpnag 
      double precision wmach(9)
      common /ax02za/wmach

      integer ldt,ldq
      common /be04nb/ldt,ldq

      double precision epspt8,epspt9
      common/ce04nb/epspt8,epspt9
c-----------------------------------------------------------------------
c                                 load arrays for lp solution
      hcp = icp
      ctotal = 0d0

      do i = 1, icp
         ctotal = ctotal + cblk(i)
      end do 
c                                 load arrays for lp solution
      jphct = iphct - istct + 1
c                                 composition constraint
      do i = 1, icp
         b(i) = cblk(i)/ctotal
      end do 

      do i = 1, jphct
         do j = 1, icp
            a(j,i) = cp(j,istct+i-1)/ctot(istct+i-1)
         end do
      end do
c                                 cold start istart = 0
      istart = 0
c                                 stuff for lpnag
      wmach(3) = 1.11022302462516d-16
      wmach(4) = dsqrt(wmach(3))
      wmach(5) = 2.22507385850721d-308
      wmach(7) = 1d0/wmach(5)
      wmach(8) = dsqrt(wmach(7))
      wmach(9) = max(1d0/wmach(4),1d2)
      wmach(2) = wmach(3)**0.8d0
      wmach(1) = wmach(3)**0.9d0

      ldt = icp + 1
      ldq = icp + 1

      end 

      subroutine gall 
c-----------------------------------------------------------------------
c subroutine gall computes molar free energies of all phases.
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i,istct,ipoint,j,icp,k,iphct,imyn,id,icomp

      double precision dg1,gval,dg,gzero,g0(k5),gsol,gsol1,gex

      common/ cst6  /icomp,istct,iphct,icp
     *      / cst60 /ipoint,imyn

      double precision g
      common/ cst2 /g(k1)

      double precision cp
      common/ cst12 /cp(k5,k1)

      double precision thermo,uf,us 
      common/ cst1 /thermo(k4,k10),uf(2),us(h5)

      integer ifct,idfl
      common/ cst208 /ifct,idfl

      integer iff,idss,ifug,ifyn,isyn
      common/ cst10  /iff(2),idss(h5),ifug,ifyn,isyn

      integer isoct
      common/ cst79 /isoct

      integer ids,isct,icp1,isat,io2
      common/ cst40 /ids(h5,h6),isct(h5),icp1,isat,io2

      integer ixp,ifp
      double precision sxs,exces
      common/ cst304 /sxs(k13),exces(m3,k1),ixp(k1),ifp(k1)

      integer jend
      common/ cxt23 /jend(h9,k12)
c                                 working arrays
      double precision z, pa, p0a, x, w, y
      common/ cxt7 /x(m4),y(m4),pa(m4),p0a(m4),z(mst,msp),w(m1)

      double precision dvnu,deph,dydy
      common/ cxt3r /dvnu(m4,j3,h9),deph(j3,h9),dydy(m4,j3,h9)
c                                 model type
      logical lorder, lexces, llaar, lrecip
      common/ cxt27 /lorder(h9),lexces(h9),llaar(h9),lrecip(h9)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      double precision numbs(21)

      integer index
c-----------------------------------------------------------------------
c                                 compute the chemical potential
c                                 of the projected components.
      call uproj
c                                 first do the endmembers:
      do id = istct, ipoint

         call gcpd(id,gval)

         g(id) = +1d2

      end do 
c                                 index to g-value based on victors data format index = T/100 - 2
      index = int(t/100) - 2
      rewind (99)
c                                 now do solutions:
      do i = 1, isoct

         do j = 1, jend(i,2) 

            read (99,*) numbs

            g(id) = numbs(index)*1d3

            id = id + 1

         end do 

      end do 

      end

      logical function solvus (id1,id2)
c-----------------------------------------------------------------------
c function to test if a solvus separates two pseudocompounds .
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      double precision ctot, dx

      integer icomp, iphct, icp, i, id1, id2, istct

      common/ cst6 /icomp,istct,iphct,icp
     *      / cst3 /ctot(k1)

      double precision cp
      common/ cst12 /cp(k5,k1)

      integer iopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10)
c-----------------------------------------------------------------------
      if (nopt(8).eq.1d0) then 
c                                 this conditional allows user to homogenize
c                                 immiscible phases be setting the tolerance to 1.  
         solvus = .false.
         return 
      end if

      dx = 0d0

      do i = 1, icp
         dx = dx + ( cp(i,id1)/ctot(id1) 
     *                     - cp(i,id2)/ctot(id2))**2
      end do 

      if (dx.ne.0) dx = dsqrt(dx)
c                                 to qualify as a solvus (as this is only
c                                 called by yclos1) require that the gap
c                                 must be greater than both the solvus tolerance
c                                 (nopt(8)) and the initial resolution (nopt(13)).
      if (dx.gt.nopt(8).and.dx.gt.nopt(13)) then
         solvus = .true.
      else
         solvus = .false.
      end if 

      end 

      subroutine ftopen (n2name,n3name,n4name,n9name,jbulk,icp,icopt,j)
c-----------------------------------------------------------------------
c open files for subroutine input1.
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      logical first

      integer ierr,icopt,jbulk,icp,j
 
      character*100 blank*1,n2name,n3name,n4name,n5name,n9name

      integer io3,io4,io9
      common / cst41 /io3,io4,io9

      save first, blank

      data first,blank/.true.,'  '/
c----------------------------------------------------------------------
c                             open thermodynamic data file
      open (n2, file = n2name, iostat = ierr, status = 'old')
      if (ierr.ne.0) call error (120,0d0,n2,n2name) 
      if (first) write (*,1170) n2name
c                                 open files if requested
      if (n3name.ne.blank) then 
         io3 = 0 
         open (n3, file = n3name)
         if (first) write (*,1180) n3name
      else
         io3 = 1
         if (first) write (*,1180) 'none requested'
      end if

      if (n4name.ne.blank) then
         io4 = 0
         open (n4, file = n4name)
         if (first) write (*,1190) n4name
      else
         io4 = 1
         if (first) write (*,1190) 'none requested'
      end if

      if (n9name.ne.blank) then
         io9 = 0 
c                             open solution model file
         open (n9,file = n9name,iostat = ierr,status = 'old')
         if (ierr.ne.0) call error (120,0d0,n9,n9name)
         if (first) then 
            write (*,1200)
            write (*,1210) n9name
         end if 
      else
         io9 = 1
         if (first) write (*,1210) 'none requested'
       end if

      if (jbulk.ge.icp.and.io4.ne.1) then
c                                 create special plot output file
c                                 name = b + plotfile name
         n5name = n4name
         call inblnk (n5name,'b')
         open (n5, file = n5name)
         if (first) write (*,1220) n5name

      else if (jbulk.ge.icp.and.io4.eq.1) then 
         if (first) write (*,1220) 'none requested'
      end if
c                              if jtest = 3, write a list of reactions for
c                              stefano
      if ((icopt.eq.1.or.icopt.eq.3).and.j.eq.3) then 
         open (n6,file='reaction_list.dat')
         if (first) write (*,1250)
      end if 

      first = .false.

1170  format (/,'Reading thermodynamic data from file: ',a)
1180  format ('Writing print output to file: ',a)
1190  format ('Writing plot output to file: ',a)
1200  format ('Writing pseudocompound glossary to file: ',
     *          'pseudocompound_glossary.dat')
1210  format ('Reading solution models from file: ',a)
1220  format ('Writing bulk composition plot output to file: ',a)
1250  format ('Writing complete reaction list to file: ',
     *          'reaction_list.dat')
      end 

      subroutine grxn (gval) 
c-----------------------------------------------------------------------
c grxn computes the free energy of univariant equilibria
c defined by the data in commonn block cst21 which is initialized
c in the subprogram balanc.  grxn is partially redundant with
c the function gphase but because of the frequency that these
c these routines are used a significant increase in efficiency is
c gained by maintaining separate functions.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer j

      double precision gval,gproj

      integer iffr,isr
      double precision vuf,vus
      common/ cst201 /vuf(2),vus(h5),iffr,isr

      integer iff,idss,ifug,ifyn,isyn
      common/ cst10  /iff(2),idss(h5),ifug,ifyn,isyn

      integer idr,ivct
      double precision vnu
      common/ cst25 /vnu(k7),idr(k7),ivct
c-----------------------------------------------------------------------
c                                 compute potentials of saturated phases
c                                 and components, note that in this
c                                 version of vertex the stoichiometry of
c                                 such components may vary.

c                                 no saturated phase components and no
c                                 saturated components:
      if (iffr.eq.1.and.isyn.eq.1) goto 10
c                                 note that this call to uproj makes a
c                                 subsequent call in gall redundant if
c                                 sfol1 is used to trace a univariant
c                                 curve.
      call uproj
c                                 compute free energy change of the rxn
10    gval = 0d0

      do j = 1, ivct
         gval = gval + vnu(j) * gproj(idr(j))
      end do 

      end

      subroutine lpopt0 (idead)
c-----------------------------------------------------------------------
c lpopt0 - calls lp minimization after a call to initlp. lpopt0
c does the minimization, writes error messages if necessary.

c this is an utterly stupid formulation of the lp problem because i modified
c the lp code to impose the implicit constraint that phase amounts were between
c 0 and 1, but did not impose the constraint the sum of the amounts must be
c be 1 (which would be unwise in any case). this requires that compositions and
c g's must be normalized to the total number of moles (and therefore lots of 
c extra bookkeeping). 
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer liw,lw,k,idead,npt,inc

      parameter (liw=2*k1+3,lw=2*(k5+1)**2+7*k1+5*k5)  

      double precision ax(k5),x(k1),clamda(k1+k5),w(lw),oldt

      integer is(k1+k5),iw(liw),jdv(k19),js(k19)
c                                 options from perplex_option.dat
      integer iopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10)

      double precision ctot
      common/ cst3  /ctot(k1)

      double precision g
      common/ cst2 /g(k1)

      integer jphct,istart
      common/ cst111 /jphct,istart

      double precision a,b,c
      common/ cst313 /a(k5,k1),b(k5),c(k1)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  
c                                 solution model counter
      integer isoct
      common/ cst79 /isoct

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      save ax, x, clamda, w, is, iw
c-----------------------------------------------------------------------
      inc = istct - 1

      oldt = t

      if (t.lt.nopt(12)) t = nopt(12)

      call gall

      do k = 1, jphct
         c(k) = g(k+inc)/ctot(k+inc)
      end do
c                                 idead = -1 tells lpnag to save parameters
c                                 for subsequent warm starts
      idead = -1
c                                 optimize by nag
      call lpnag (jphct,icp,a,k5,b,c,is,x,ax,
     *            clamda,iw,liw,w,lw,idead,l6,istart)

      if (idead.gt.0) then
c                                 look for severe errors                                            
         call lpwarn (idead,'LPOPT ')
c                                 on severe error do a cold start.
c                                 necessary?
         istart = 0
         return

      end if 

      if (icp.eq.1.or.iopt(10).eq.0.or.isoct.eq.0) then 
c                                 no refinement, final processing
c                                 with yclos0
         call yclos0 (x,is,jphct,jdv,idead) 
c                                 get chemical potentials, this is
c                                 just for output.
         call getmus (jdv,1,idead)

      else
c                                 find discretization points
c                                 for refinement
         call yclos1 (clamda,x,is,jphct,jdv,js,npt,idead)

         if (idead.gt.0) return 
c                                 reoptimize with refinement
         call reopt (idead,jdv,npt,js)
c                                 coming out of reopt the stable points
c                                 are indexed by jdv
         if (idead.gt.0) return 
c                                 get amounts of pseudocompounds
         call rebulk (jdv,idead)

      end if 
c                                 test for solvi and average
      if (idead.eq.0) call avrger       

      t = oldt

      end 

      subroutine yclos0 (x,is,jphct,jdv,idead)
c----------------------------------------------------------------------
c subroutine to save optimization results for non-iterative refinement
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer jdv(k19),inc,i,j,iff,idss,idead,jphct,mpt,
     *        ifug,ifyn,isyn,ipvt,npt,is(k1+k5),id

      double precision a,b,cp,x(k1)
 
      common/ cst10  /iff(2),idss(h5),ifug,ifyn,isyn
     *      / cst301 /a(k5,k5),b(k5),ipvt(k5)
     *      / cst12 /cp(k5,k1)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      double precision ctot
      common/ cst3  /ctot(k1)

      integer jbulk
      double precision cblk
      common/ cst300 /cblk(k5),jbulk
c                                 compositions of stable adaptive
c                                 coordinates (and solution ids).
      integer kkp,np,ncpd,ntot
      double precision cp3,ctot3
      common/ cxt15 /cp3(k5,k5),ctot3(k5),kkp(k5),np,ncpd,ntot
c                                 options from perplex_option.dat
      integer iopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10)

      integer ipoint,imyn
      common/ cst60 /ipoint,imyn

      integer ikp
      common/ cst61 /ikp(k1)
c----------------------------------------------------------------------

      npt = 0 
      inc = istct - 1

      do i = 1, jphct

         if (is(i).ne.0.and.is(i).ne.2) cycle  
c                                 new point, add to list, two acceptable
c                                 cases is = 0 inactive 
c                                       is = 2 active at upper bound (this
c                                            is really only possible for a 
c                                            one component system (though
c                                            conceivably might occur in 
c                                            fractionation calculations).
            npt = npt + 1

            jdv(npt) = i 
 
      end do

      if (icp.eq.1.and.npt.gt.1) then 

         npt = 0 

         do i = 1, jphct

            if (is(i).ne.2) cycle 

            npt = npt + 1
            jdv(npt) = i

         end do

      end if 

      if (npt.gt.icp) then 
c                                 look through for zero phases
         do i = 1, npt
            if (x(jdv(i)).eq.0d0) idead = idead + 1
         end do

         if (idead.ge.npt-icp) then 

            idead = npt - icp 
            mpt = 0
c                                 eliminate zero phase
            do i = 1, npt 

               if (x(jdv(i)).eq.0d0.and.idead.gt.0) then 
                  idead = idead - 1
c                                 reset the active flag
c                 is(jdv(i)) = 1
                  cycle
               end if 

               mpt = mpt + 1
               jdv(mpt) = jdv(i)

            end do 

            npt = mpt 

         else 

            idead = 1

         end if 

      else if (npt.lt.icp) then 

         idead = 1 

      end if 

      if (idead.eq.1) then 
        call warn (88,0d0,idead,'LPOPT') 
        goto 99
      end if 
c                                the assemblage is ok, get the modes and 
c                                save the data, this means creating the kkp
c                                array, the x3 array (for output), the cp3
c                                array (only for meemum)
      do i = 1, jbulk

         if (i.le.icp) then 
c                                normal phase
            id = jdv(i) + inc
         else
c                                saturated phase
            id = idss(i-icp)
         end if 
c                                set identifier flag
         if (i.le.ipoint) then
            kkp(i) = -id
         else  
            kkp(i) = ikp(id)
         end if 
c                                load the composition matrix
         do j = 1, jbulk
            a(j,i) = cp(j,id)
            cp3(j,i) = cp(j,id)
         end do
c                                total component amounts
         ctot3(i) = ctot(id)
c                                set the x3 array
         if (ikp(id).ne.0) call setx3 (i,id,ikp(id))

      end do 
c                                 factor the matrix
      call factr1 (jbulk,idead)

      if (idead.eq.1) goto 99 
c                                 get mode
      do i = 1, jbulk
c                                 load composition vector, the alpha
c                                 vector is returned in the same array
         b(i) = cblk(i)
      end do
c                                 solve for the alpha vector
      call subst1 (jbulk)

      ntot = jbulk
c                                 test for zero modes
      if (nopt(9).gt.0d0) call zmode

99    end 

      subroutine zmode
c----------------------------------------------------------------------
c rebulk computes the amounts of the stable compounds and eliminates 
c those with zero modes (if requested).
c----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i,j,k,l,ipvt,jzero(k5)

      double precision a,b
      common/ cst301 /a(k5,k5),b(k5),ipvt(k5) 

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      integer jbulk
      double precision cblk
      common/ cst300 /cblk(k5),jbulk
c                                 compositions of stable adaptive
c                                 coordinates (and solution ids).
      integer kkp,np,ncpd,ntot
      double precision cp3,ctot3
      common/ cxt15 /cp3(k5,k5),ctot3(k5),kkp(k5),np,ncpd,ntot
c                                 options from perplex_option.dat
      integer iopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10)
c                                  x-coordinates for the final solution
      double precision x3
      common/ cxt16 /x3(k21,mst,msp)
c                                 x coordinate description
      integer istg, ispg, imlt, imdg
      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)
c----------------------------------------------------------------------

      ntot = 0 

      do i = 1, jbulk
         if (b(i).lt.-nopt(9)) then
            write (*,*) ' a negative mode, should warn'
            write (*,*) (kkp(j),j=1,jbulk)
            write (*,*) (b(j),j=1,jbulk)
         else if (b(i).gt.nopt(9)) then 
            ntot = ntot + 1
            jzero(ntot) = i 
         end if 
      end do 

      if (jbulk.ne.ntot) then

         do i = 1, ntot

            j = jzero(i)
c                                 jkp and jdv are probably unnecessary.
            kkp(i) = kkp(j)
            b(i) = b(j)

            do k = 1, icomp
               cp3(k,i) = cp3(k,j)
            end do 

            ctot3(i) = ctot3(j)

            if (kkp(i).gt.0) then 

               do k = 1, istg(kkp(i))
                  do l = 1, ispg(kkp(i),k)
                     x3(i,k,l) = x3(j,k,l)
                  end do 
               end do 
            end if 
         end do 
      end if

      end

      subroutine setx3 (ind,id,ids)
c----------------------------------------------------------------------
c subroutine to recover geometric reciprocal solution compositions (x(i,j))
c from the xcoor (reciprocal) or sxs (single site) arrays loaded in soload
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, id, ids, jcoor, ind
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
c                                  x-coordinates for the final solution
      double precision x3
      common/ cxt16 /x3(k21,mst,msp)

      integer ipoint,imyn
      common/ cst60 /ipoint,imyn

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)
c----------------------------------------------------------------------
      if (id.gt.ipoint) then 
c                                 a normal solution
         if (istg(ids).gt.1) then 

            jcoor = icoor(id)

            do i = 1, istg(ids)
               do j = 1, ispg(ids,i)
                  jcoor = jcoor + 1
                  x3(ind,i,j) = xcoor(jcoor)
               end do 
            end do 

         else 

            do j = 1, nstot(ids)
               x3(ind,1,j) = sxs(ixp(id)+j) 
            end do 

         end if 
      else 
c                                 an endmember 
         call endcp (ind,id,ids)
      end if

      end 

      subroutine getmus (jdv,iam,idead)
c-----------------------------------------------------------------------
c getmus computes chemical potentials from a pseudocompound assemblage
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, idead,  ier, iam

      double precision comp(k8,k8)

      integer ipvt(k8), jerk

      integer jdv(k19)

      integer jphct
      double precision g2, cp2
      common/ cxt12 /g2(k21),cp2(k5,k21),jphct

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      double precision g
      common/ cst2  /g(k1)

      double precision mu
      common/ cst330 /mu(k8)

      integer iopt
      double precision nopt
      common / opts /nopt(i10),iopt(i10)

      integer jtest,jpot
      common/ debug /jtest,jpot

      double precision cp
      common/ cst12 /cp(k5,k1)

      save jerk 

      data jerk/0/ 
c-----------------------------------------------------------------------
      if (jpot.eq.1) return

      if (idead.eq.0) then 
c                                 compute chemical potentials if requested:
         if (iam.eq.1) then
c                                 non-adaptive solution
            do i = 1, icp
               mu(i) = g(jdv(i))
               do j = 1, icp
                  comp(i,j) = cp(j,jdv(i))
               end do
            end do

         else
c                                 adaptive solution
            do i = 1, icp
               mu(i) = g2(jdv(i))
               do j = 1, icp
                  comp(i,j) = cp2(j,jdv(i))
               end do
            end do

         end if 
c                                 compute chemical potentials if requested:
c                                 reload last solution (no test for uniqueness)
         call factor (comp,icp,ipvt,ier)

         if (ier.eq.1.and.jerk.lt.5) then 

            jerk = jerk + 1
            write (*,*) 'Deader than a doornail, but dont worry' 
     
         else 
 
            call subst (comp,ipvt,icp,mu,ier)

            if (ier.eq.1.and.jerk.lt.5) then 
               jerk = jerk + 1
               write (*,*) 'Deader than 2 doornails, but dont worry' 
            end if 

         end if 

      else
c                                 optimization failed.
         do i = 1, icp
            mu(i) = nopt(7)
         end do

      end if 

      end 

