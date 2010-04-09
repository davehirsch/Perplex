
c Please do not distribute any part of this source.
 
c Copyright (c) 1998 by James A. D. Connolly, Institute for Mineralogy
c & Petrography, Swiss Federal Insitute of Technology, CH-8092 Zurich,
c SWITZERLAND. All rights reserved.

      program vertx        
c----------------------------------------------------------------------
c                       ************************
c                       *                      *
c                       *  vertex.sep.10.1998  *
c                       *                      *
c                       ************************
c----------------------------------------------------------------------
 
c references:
 
c       connolly jad, kerrick dm (1987) bounds:an algorithm and program
c         for calculating composition diagrams, calphad 11, 1-57.
c       connolly jad (1990) multivariable phase diagrams: an algorithm
c         based on generalized thermodynamics, am j sci 290, p 666-718.
 
c-----------------------------------------------------------------------
c file structure:
 
c i/o logical unit numbers are specified in perplex_parameters.h 
c-----------------------------------------------------------------------
c this version of vertex has been modified to be compatible with cdc
c fortran 5 compilers.  memory requirements can be reduced on most
c other compilers by changing the implicit variable type integer to
c integer*2, on ibm compilers integer variables beginning with z must
c be integer*4, i.e., for ibm compilers the current implicit type
c statement 'implicit integer (h-n,z)' should be changed to 'implicit
c integer (z), integer*2 (h-n).  the necessity for double precision real
c variables has not been established, users interested in reducing 
c memory requirements may wish to experiment with single precision 
c variables.
c-----------------------------------------------------------------------
c functional(1) tolerance(2) constant(3) or zero (4) dependent
c     subprograms:                         (this list is not current)

c                 testit (4)      sreset (2)      asschk (4)
c itestc   (4)    loadit (1)      flipit (2)      wway   (2)
c                 solod  (1,3)    schk   (4)      concrt (2,3)
c miscib   (4)    newass (4)      pchk   (2,4)    univeq (2,3)
c balanc   (4)    uproj  (1,2,3)  grxn   (1)      gphase (1)
c mrk      (1,3)  hsmrk  (1,3)    newrap (1,3)    fug    (1,3)

c units, tolerances, and constants:

c function specific
c routines for fluid properties (mrk, trkmrk, hsmrk) require units of
c kelvins (v(2)) and bars (v(3)) and joules.
c-----------------------------------------------------------------------
      implicit none
c-----------------------------------------------------------------------

c parameters are assigned in "perplex_parameter.h"
 
c-----------------------------------------------------------------------
      include 'perplex_parameters.h'

      logical vertex, output, first,pots  
      
      integer ier99

      character n4name*100

      integer io3,io4,io9
      common / cst41 /io3,io4,io9

      integer isec,icopt,ifull,imsg,io3p
      common/ cst103 /isec,icopt,ifull,imsg,io3p

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      integer iemod,kmod
      logical smod
      double precision emod
      common/ cst319 /emod(k15,k10),smod(h9),iemod(k10),kmod

      integer jtest,jpot
      common/ debug /jtest,jpot

      logical  debug
      common / debugblk / debug

      save first,output,vertex,pots
      data vertex,output,first/.true.,.false.,.true./
c----------------------------------------------------------------------- 
c   Look for the "debug_yes" file to turn on debugging messages
      open (99,iostat=ier99,file='debug_yes',status='old')
      if (ier99.eq.0) then
          debug=.TRUE.
          close (99)
      else
          debug=.FALSE.
      end if

c                                 elastic modulii flag
      kmod = 0 
c                                 this do loop is a cheap strategy to automate
c                                 "auto_refine"
      do
c                                 -------------------------------------
c                                 open statements for units n1-n6 and n9
c                                 are in subroutine input1
 
c                                 read input from unit n1 (terminal/disk).
c                                 input1 also initializes: conditions,
c                                 equilibrium counters; units n2 n4 and n6;
c                                 and the limits for numerical results.
         call input1 (first,output,n4name)
c                                 read thermodynamic data on unit n2:
         call input2 (first)
c                                 read/set autorefine dependent parameters, 
c                                 it would be logical to output context specific 
c                                 parameter settings here instead of the generic 
c                                 blurb dumped by redop1
         call setau1 (vertex,output)
c                                 read data for solution phases on n9:
         call input9 (vertex,first,output)

         call setau2 (output)
c                                 initialize potentials
         call inipot 
c                                 -------------------------------------
c                                 at this point the problem is fully 
c                                 configured, 
c                                 -------------------------------------
         if (output) then 
c                                 header info for print and graphics files
            call topout (icopt)
c                                 inform user of 1st stage
            if (iopt(6).ne.0) write (*,1000) 'auto_refine'
c                                 turn of printing of potentials if no
c                                 print file request.
            if (.not.first.and.pots) jpot = 0
            if (icopt.lt.5.and.io3.eq.1) jpot = 1

         else 
c                                 inform user of 2nd stage
            if (iopt(6).ne.0) write (*,1000) 'exploratory'
c                                 suppress output to graphics and print files
c                                 (these flags are reset by input1). 
            io4 = 1
            io3 = 1

            if (jpot.ne.1) then 
               pots = .true.
            else 
               pots = .false.
            end if 

            jpot = 1

         end if 
            

         if (icopt.eq.0) then
c                                 calculate composition phase diagrams
c                                 calculations and remaining output
            call chmcal (output)

         else if (icopt.eq.1.or.icopt.eq.3) then                            
c                                 phase diagram projection or mixed variable
c                                 diagram 
            call newhld (output)

         else if (icopt.eq.4) then 
c                                 generate pseudo-compound file
            call swash 

            exit 
 
         else if (icopt.eq.5) then 
c                                 optimization on a 2-d grid.
            call wav2d1 (output)

         else if (icopt.eq.7.or.icopt.eq.10) then 
c                                 fractionation on a 1-d path.
            call frac1d (output)

         else if (icopt.eq.8) then 
c                                 a g-x data output format for manual optimization
            call gwash 

            exit 

         else if (icopt.eq.9) then 
c                                 fractionation on a 2-d path (space-time) 
            call frac2d (output)

         else      
c                                 disabled stability field calculation
            call error (32,0d0,k2,'MAIN')

         end if 
c                                 output compositions for autorefine
         call outlim 

         if (output) exit

         output = .true.   
         first = .false.

      end do 

1000  format (/,'** Starting ',a,' computational stage **',/)

      end

      subroutine topout (icopt)
c-----------------------------------------------------------------------
c call various routines to generate the header section of graphics and
c print files depending on computational model (icopt)
c-----------------------------------------------------------------------
      implicit none

      integer icopt 

      integer io3,io4,io9
      common / cst41 /io3,io4,io9
c-----------------------------------------------------------------------
c                              write header to graphics files (n4, n5):
      if (io4.ne.1) then
 
         if (icopt.eq.1) then 
c                              computing schreinemakers projection
            call header (icopt)

         else if (icopt.le.3) then 
c                              computing chemography or mixed-var.
            call outhed

         end if

      end if 
c                              title page for print file:
      if (io3.ne.1) call outtit

      end 

      subroutine frac1d (output)
c-----------------------------------------------------------------------
c a main program template to illustrate how to call the minimization 
c procedure. the example here is 1d fractionation.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical output

      integer i,j,idead,ier

      double precision ctotal

      character*8 cname*5, y*1, xname, vname
      common/ csta4 /cname(k5)/ csta2 /xname(k5),vname(l2)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      character*100 n1name,cfname
      common/ cst228 /n1name,cfname

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      integer io3,io4,io9
      common / cst41 /io3,io4,io9

      double precision dcomp
      common/ frct2 /dcomp(k5)

      double precision a,b,c
      common/ cst313 /a(k5,k1),b(k5),c(k1)

      integer ipot,jv,iv
      common/ cst24 /ipot,jv(l2),iv(l2)

      integer jbulk
      double precision cblk
      common/ cst300 /cblk(k5),jbulk

      integer jlow,jlev,loopx,loopy,jinc
      common/ cst312 /jlow,jlev,loopx,loopy,jinc

      integer idasls,iavar,iasct,ias
      common/ cst75  /idasls(k5,k3),iavar(3,k3),iasct,ias

      integer iap,ibulk
      common/ cst74 /iap(k2),ibulk

      integer ifrct,ifr
      logical gone
      common/ frct1 /ifrct,ifr(k23),gone(k5)

      integer isec,icopt,ifull,imsg,io3p
      common/ cst103 /isec,icopt,ifull,imsg,io3p
c-----------------------------------------------------------------------

      iasct = 0
      ibulk = 0  

      loopy = jlow
c                                 get phases to be fractionated
      call frname 
c                                 call initlp to initialize arrays 
c                                 for optimization.
      call initlp 
c                                 two cases icopt = 10, input from
c                                 file, else analytical path function
c                                 (icopt = 7)
      if (icopt.eq.10) then 

         open (n8,file=cfname,status='old',iostat=ier)

         if (ier.ne.0) call error (6,v(1),i,cfname)

         j = 0 

         do 

            read (n8,*,iostat=ier) (v(jv(i)), i = 1, ipot)
            if (ier.ne.0) exit

            j = j + 1

            ctotal = 0d0
c                                 get total moles to compute mole fractions             
            do i = 1, icp
               ctotal = ctotal + cblk(i)
            end do

            do i = 1, icp 
               b(i) = cblk(i)/ctotal
            end do
c                                 lpopt does the minimization and outputs
c                                 the results to the print file.
            call lpopt (1,j,idead,output)
c                                 fractionate the composition:
            if (idead.eq.0) call fractr (output)

         end do 

         close (n8)

         loopy = j 

      else
c                                 initialize sectioning variables
         call setvar 
c                                 loopx = 99 indicates interactive
c                                 mode, otherwise traverse on  
         if (loopx.eq.99) then 

            do 

               write (*,1070) (vname(jv(i)), i = 1, ipot)
               read (*,*) (v(jv(i)), i = 1, ipot)
               if (v(jv(1)).eq.0d0) exit
c                                 the systems molar composition is in 
c                                 the array cblk.
               ctotal = 0d0
c                                 get total moles to compute mole fractions             
               do i = 1, icp
                  ctotal = ctotal + cblk(i)
               end do

               do i = 1, icp 
                  b(i) = cblk(i)/ctotal
               end do
c                                 lpopt does the minimization and outputs
c                                 the results to the print file.
               call lpopt (1,j,idead,output)

               write (*,1000)
               do i = 1, jbulk
                  write (*,1010) cname(i), cblk(i)
               end do 
c                                 the bulk composition can be modified here, 
c                                 e.g., to add back previously fractionated
c                                 phases.
               write (*,1060) 
               read (*,1050) y
               if (y.eq.'y'.or.y.eq.'Y') then 
                  write (*,1020) 
                  read (*,*) (dcomp(i), i = 1, jbulk)
                  do i = 1, jbulk
                     cblk(i) = cblk(i) + dcomp(i)
                  end do 
               end if 

            end do

         else
c                                 automated mode:
            do j = 1, loopy

               ctotal = 0d0
c                                 get total moles to compute mole fractions             
               do i = 1, icp
                  ctotal = ctotal + cblk(i)
               end do

               do i = 1, icp 
                  b(i) = cblk(i)/ctotal
               end do

               call setvr0 (j,j)  
c                                 lpopt does the minimization and outputs
c                                 the results to the print file.
               call lpopt (1,j,idead,output)
c                                 fractionate the composition:
               if (idead.eq.0) call fractr (output)

            end do 

         end if 
      end if 
c                                 output 
      if (output) then 

        if (io4.eq.0) call outgrd (1,loopy,1) 
c                                 close fractionation data files
         do i = 1, ifrct
            close (n0+i)
         end do

         close (n0)

      end if 

1000  format (/,'Composition is now:',/)
1010  format (1x,a,1x,g14.6)
1020  format (/,'Enter molar amounts of the components to be added ',
     *        '(ordered as above:')
1050  format (a)
1060  format (/,'Modify composition (y/n)?')
1070  format (/,'Enter (zeroes to quit) ',7(a,1x))

      end 

      subroutine frac2d (output)
c-----------------------------------------------------------------------
c a subprogram template to illustrate how to do 2-d (space-time) 
c fractionation. 

c modified oct 15 for more rational input and variable gradients.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer maxbox,maxlay

      parameter (maxbox=1760,maxlay=6) 

      integer i,j,k,idead

      double precision ctotal,gblk(maxbox,k5),dz,p0

      logical output
      character*8 happy*3

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      integer io3,io4,io9
      common / cst41 /io3,io4,io9

      double precision a,b,c
      common/ cst313 /a(k5,k1),b(k5),c(k1)

      double precision dcomp
      common/ frct2 /dcomp(k5)

      integer idasls,iavar,iasct,ias
      common/ cst75  /idasls(k5,k3),iavar(3,k3),iasct,ias

      integer ipot,jv,iv
      common/ cst24 /ipot,jv(l2),iv(l2)

      double precision vmax,vmin,dv
      common/ cst9  /vmax(l2),vmin(l2),dv(l2)  

      integer jbulk
      double precision cblk
      common/ cst300 /cblk(k5),jbulk

      integer jlow,jlev,loopx,loopy,jinc
      common/ cst312 /jlow,jlev,loopx,loopy,jinc 

      integer iap,ibulk
      common/ cst74 /iap(k2),ibulk

      integer gloopy,ilay,irep
      double precision a0,a1,a2,a3,b0,b1,b2,b3,c0,c1,c2,c3,dv1dz,
     *               zbox,iblk
      common/ cst66 /a0,a1,a2,a3,b0,b1,b2,b3,c0,c1,c2,c3,dv1dz,
     *               zbox,iblk(maxlay,k5),gloopy,ilay,irep(maxlay)

      logical first

      save first,ctotal

      data first/.true./
c-----------------------------------------------------------------------
c                                 initialize
      iasct = 0 
      ibulk = 0 
c                                 jlow is set from 1dpath in perplex_option.dat
      loopy = jlow
c                                 set the number of independent variables
c                                 to 1 (the independent path variable must
c                                 be variable jv(1), and the dependent path
c                                 variable must be jv(2), the path variables
c                                 can only be pressure and temperature
      ipot = 1

      if (first) then 
c                                 set first to prevent re-reading of the 
c                                 input in auto_refine                                
         first = .false.
c                                 get the phase to be fractionated
         call frname 

      end if 
c                                 initialize compositional array gblk:
      loopx = 0 

      do i = 1, ilay

         do j = 1, irep(i)

            loopx = loopx + 1

            if (loopx.gt.maxbox) then
               write (*,*) 'increase maxbox in routine DUMMY1'
               stop
            end if 

            do k = 1, icp 
               gblk(loopx,k) = iblk(i,k)
            end do 

         end do

      end do

      if (loopx*loopy.gt.k2) then
         write (*,*) ' parameter k2 must be >= loopx*loopy'
         write (*,*) ' increase parameter k2 for routine DUMMY1'
         write (*,*) ' or increase box size (zbox) or decrease'
         write (*,*) ' number of path increments (loopy) or try'
         write (*,*) ' the large parameter version of VERTEX'
         write (*,*) ' k2 = ',k2 
         write (*,*) ' loopx * loopy = ',loopx*loopy
         stop
      end if 
c                                 call initlp to initialize arrays 
c                                 for optimization.
      call initlp 
c                                 initialize path variables
      call setvar 

      dv(jv(1)) = (vmax(jv(1)) - vmin(jv(1)))/(loopy-1)
c                                 loopy is the number of steps along the subduction
c                                 path:
      do j = 1, loopy

         happy = 'no'
c                                 j loop varies the pressure at the top of the 
c                                 column (p0), set p0:
         p0 = vmin(1) + (j-1)*dv(1)
c                                 happy is set to yes when we've finally finished 
c                                 equilibrating at point j on the path
         do while (happy.eq.'no') 

            happy = 'yes'
c                                 for each box in our column compute phase 
c                                 relations
            do k = 1, loopx
c                                 k-loop varies depth within the column (dz)
               dz = (loopx - k) * zbox  
c                                 set the p-t conditions
               call fr2dpt (p0,dz)
c                                 now put the bulk composition from the global
c                                 array into the local array and get the total
c                                 number of moles (ctotal)
               ctotal = 0d0
c                                 get total moles to compute mole fractions             
               do i = 1, icp
                  dcomp(i) = 0 
                  cblk(i) = gblk(k,i)
                  ctotal = ctotal + cblk(i)
               end do

               do i = 1, icp 
                  b(i) = cblk(i)/ctotal
               end do
c                                 lpopt does the minimization and outputs
c                                 the results to the print file.
               if (io3.eq.0) write (n3,*) 'step ',j,' on path, box ',k

               call lpopt (j,k,idead,output)

               call fractr (output)
c                                 at this point we've computed the stable
c                                 assemblage at each point in our column
c                                 and could do mass transfer, etc etc 
c                                 here we'll simply fractionate the fluid 
c                                 phase 
               do i = 1, icp 
c                                 subtract the fluid from the current composition
                  gblk(k,i) = gblk(k,i) - dcomp(i) 
                  if (gblk(k,i).lt.0d0) then
                     write (*,*) 'negative - bulk at k,i',k,i,gblk(k,i)
                      gblk(k,i) = 0d0
                  end if 
c                                 and add it to the overlying composition
                  if (k.lt.loopx) then
                      gblk(k+1,i) = gblk(k+1,i) + dcomp(i) 
                      if (gblk(k,i).lt.0d0) then
                         write (*,*) 'bulk < 0, at k,i',k,i,gblk(k,i)
                         gblk(k,i) = 0d0
                      end if 
                  end if 
               end do
c                                 reset initialize top layer if gloopy = 999
               if (gloopy.eq.999.and.k.gt.loopx-irep(ilay)) then 
                  do i = 1, icp
                     gblk(k,i) = iblk(ilay,i)
                  end do 
               end if 
c                                 end of the k index loop
            end do 
c                                 end of the happy loop. 
         end do 
c                                 end of j index loop
      end do 

      if (output.and.io4.eq.0) call outgrd (loopy,loopx,1) 

      end 

      subroutine chmcal (output)
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'
 
      integer jpoly, ier

      logical output

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      integer io3,io4,io9
      common / cst41 /io3,io4,io9
c-----------------------------------------------------------------------
      jpoly = 0 
 
      do 

         read (n1,*,iostat=ier) v

         call incdp1 
c                                 for old file formats:
         if (ier.ne.0.or.v(1).lt.0d0) exit

         jpoly = jpoly + 1
         write (*,1010) jpoly
c                                 compute phase diagram
         call gall
         call combin

         if (output) then 
c                                 graphics file
            if (io4.ne.1) call outgrf
c                                 print file:
            if (io3.ne.1) call outchm

         end if 

      end do 

      close (n1)

1010  format ('Computing the compositional phase relations at',
     *        ' condition ',i2)
 
      end

      subroutine outchm
c-------------------------------------------------------------------
c outchm writes new chemographies to the print file.
c-------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'
 
      character*8 names,vname,xname

      integer ipoint,imyn,jasmbl,iasmbl,i,j

      common/ cst8  /names(k1)
     *      / cst60 /ipoint,imyn/ cst27 /jasmbl(j9),iasmbl(j9)
     *      / csta2 /xname(k5),vname(l2)

      integer iff,idss,ifug,ifyn,isyn
      common/ cst10  /iff(2),idss(h5),ifug,ifyn,isyn

      integer idcf,icfct
      common/ cst96 /idcf(k5,j9),icfct

      integer ipot,jv,iv1,iv2,iv3,iv4,iv5
      common/ cst24 /ipot,jv(l2),iv1,iv2,iv3,iv4,iv5

      integer isoct
      common/ cst79 /isoct

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  

      integer ids,isct,icp1,isat,io2
      common/ cst40 /ids(h5,h6),isct(h5),icp1,isat,io2
c-----------------------------------------------------------------------
c                             output independent potential variable
c                             constraints:
      write (n3,1160)
      write (n3,1000)
      write (n3,1010) (vname(jv(i)), v(jv(i)), i = 1, ipot)
      write (n3,1015)
c                             next the stable assemblages:
      if (icp.gt.4) then 
c                             output high order chmographies
         do i = 1, icfct
            write (n3,7020) iasmbl(i), (names(idcf(j,i)), j = 1, icp)
         end do
      else if (icp.eq.1) then  
c                             output unary chemographies.
         write (n3,1020) names(idcf(1,1))
      else if (icp.eq.2) then 
c                             output binary chemographies.
         write (n3,2020) ((names(idcf(j,i)), j = 1, icp), iasmbl(i), 
     *                 i = 1, icfct)
      else if (icp.eq.3) then 
c                             output ternary chemographies:
         write (n3,3020) ((names(idcf(j,i)), j = 1, icp), iasmbl(i), 
     *                i = 1, icfct)
      else if (icp.eq.4) then 
c                             output quaternary chemographies:
         write (n3,4020) ((names(idcf(j,i)), j = 1, icp), iasmbl(i), 
     *                i = 1, icfct)
      end if 
c                             output phases consistent with component
c                             saturation constraints:
      if (isyn.eq.1) goto 9000

      write (n3,6000)
      write (n3,2030) (names(idss(i)), i = 1, isat)
c
9000  if (icp.gt.1.and.isoct.gt.0) then

         write (n3,'(/)')
c                             write miscibility blurb
         if (imyn.eq.1) then
            write (n3,6010) 
         else
            write (n3,6030) 
         end if
      end if 
 
1000  format ('the stable assemblages at:',/)
1010  format (25x,a,' = ',g12.6)
1015  format (/,'are (variance flag in parentheses):',/)
1020  format (25x,a)
1160  format (/,80('-'),/)
2020  format (3(a,'-',a,'(',i1,')',3x))
2030  format (6(1x,a))
3020  format (2(2(a,'-'),a,'(',i1,')',2x))
4020  format (2(3(a,'-'),a,'(',i1,')',2x))
6000  format (/,'these assemblages are compatible with the followi',
     *        'ng phases or species',/,'determined by component ',
     *        'saturation or buffering constraints:',/)
6010  format ('** no immiscibility occurs in the stable solution ',
     *        'phases **',/)
6030  format ('** immiscibility occurs in one or more of the ',
     *        'stable solution phases **',/)
7020  format ('(',i1,')',12(1x,a8))
 
      end
 
      subroutine outhed
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'
 
      character*8 names,vname,xname,fname*10

      integer ipoint,imyn,i,j

      double precision ctot

      character*162 title
      common/ csta8 /title(4)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  
 
      common/ cst3   /ctot(k1)/ cst8 /names(k1)
     *      / csta2  /xname(k5),vname(l2)
     *      / cst60  /ipoint,imyn/ csta7 /fname(h9)

      double precision cp
      common/ cst12 /cp(k5,k1)

      integer iff,idss,ifug,ifyn,isyn
      common/ cst10  /iff(2),idss(h5),ifug,ifyn,isyn

      integer isoct
      common/ cst79 /isoct

      integer ipot,jv,iv
      common/ cst24 /ipot,jv(l2),iv(l2)

      integer isec,icopt,ifull,imsg,io3p
      common/ cst103 /isec,icopt,ifull,imsg,io3p

      integer ikp
      common/ cst61 /ikp(k1)
c-----------------------------------------------------------------------
c                             number of components, phase counters,
c                             assemblage counter, fluid saturation flag,
c                             component saturation flag, ipot = number
c                             of independent potential variables:
      write (n4,1020) icopt
c
      write (n4,1020) icp,istct,iphct,ipoint,ifyn,isyn,ipot,isoct
c                             write graphics code variable names:
      write (n4,1000) (vname(jv(i)), i = 1, ipot)
c                             write a blank record as left caption
      write (n4,1000) title(1)
c                             write phase names
      write (n4,1015) (names(i), i = 1, iphct)
c                             write phase coordinates
      write (n4,1025) ((cp(j,i)/ctot(i), j = 2, icp), i = istct, iphct)
c                             write solution phase flags:
      write (n4,1090) (ikp(i), i = 1, iphct)
c                             solution names
      if (isoct.ne.0) write (n4,1080) (fname(i), i = 1, isoct)
c                             write graphics code names for x variables:
      write (n4,1015) (xname(i), i = 1, icp)

1000  format (a)
1015  format (10a8)
1020  format (20(i8,1x))
1025  format (9(f8.6,1x))
1080  format (8a10)
1090  format (20(i4,1x))
 
      end

      subroutine outgrf
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i,j
 
      integer iff,idss,ifug,ifyn,isyn
      common/ cst10  /iff(2),idss(h5),ifug,ifyn,isyn

      integer jasmbl,iasmbl
      common/ cst27  /jasmbl(j9),iasmbl(j9)

      integer idcf,icfct
      common/ cst96 /idcf(k5,j9),icfct

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      integer ids,isct,icp1,isat,io2
      common/ cst40 /ids(h5,h6),isct(h5),icp1,isat,io2

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  

      integer ipot,jv,iv
      common/ cst24 /ipot,jv(l2),iv(l2)
c-----------------------------------------------------------------------
c                             stable configurations, phases are
c                             labelled by the index 'i' in the
c                             list of phases.
c                             write values of the independent potentials
      write (n4,1040) (v(jv(i)), i = 1, ipot)

      if (icp.ne.2) then
         write (n4,*) icfct
      else
         write (n4,*) icfct+1
      end if 

      if (icp.eq.2) then 
c                             binary is a special case (1-d)
         write (n4,1020) (idcf(1,j), j = 1, icfct),idcf(2,icfct)

      else if (icp.eq.1) then 
         goto 99
      else 
c                             higher order:
         write (n4,1020) ((idcf(j,i), j = 1, icp), i = 1, icfct)
      end if 
c                             write assemblage flags
      if (icp.gt.2) write (n4,1020) (iasmbl(j), j = 1, icfct)

      if (isyn.eq.1) goto 99

      write (n4,1020) isat
      write (n4,1020) (idss(i), i = 1, isat)

1020  format (20(i8,1x))
1040  format (5(g12.6,1x))

99    end
 
      subroutine prtpot
c---------------------------------------------------------------------
c prtpot outputs the dependent potentials for phase assemblages
c determined by the combin subroutine. it is called if the flag
c jpot, in common debug, is not equal 1.
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      character*8 names
      
      integer i

      integer iar
      double precision a,u
      common/ cst23 /a(k8*k8),u(k8),iar(2*k8+4)

      common/ cst8 /names(k1)

      integer h,id
      common/ cst52 /h,id(k7)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  

      call abload (*99)

      if (icp.gt.5) then 
         write (n3,1060) (names(id(i)),i=1,icp)
         write (n3,1070) (u(i),i=1,icp)
      else if (icp.eq.1) then 
         write (n3,1010) names(id(1)),u(1)
      else if (icp.eq.2) then 
         write (n3,1020) (names(id(i)),i=1,icp),(u(i),i=1,icp)
      else if (icp.eq.3) then 
         write (n3,1030) (names(id(i)),i=1,icp),(u(i),i=1,icp)
      else if (icp.eq.4) then 
         write (n3,1040) (names(id(i)),i=1,icp),(u(i),i=1,icp)
      else if (icp.eq.5) then 
         write (n3,1050) (names(id(i)),i=1,icp),(u(i),i=1,icp)
      end if 

1010  format (1x,a8,1x,g14.7)          
1020  format (2(1x,a8),2(1x,g14.7))
1030  format (3(1x,a8),3(1x,g14.7))
1040  format (4(1x,a8),4(1x,g14.7))
1050  format (5(1x,a8),5(1x,g14.7))
1060  format (16(1x,a8))
1070  format (16(1x,g14.7))

99    end 

      subroutine newhld (output)
c-----------------------------------------------------------------------
c newhld replaces nohold. newhld is designed to do boundary traverses
c for individual equilibrium configurations and in certain cases,
c specifically for systems with fewer than 4 components, is 
c slower than nohold. in contrast to newhld, nohold simultaneously
c monitored all assemblages on the thermodynamic surface of the system.
c a copy of the nohold code is stored in check.fortran on tape xyw003.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      character*8 vname,xname*8,cname*5

      logical output

      integer iedge,ivfl,ipot,jv,iv1,iv2,iv3,iv4,iv5,irchk,
     *        jok,kok,index,i,j,irend,knct,ier,iste,jnct,ivi,ivd

      double precision vt,delt,dtol,utol,ptol,vti,div    

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  

      common/ csta4 /cname(k5)
     *      / csta2 /xname(k5),vname(l2)/ cst102 /ivfl
     *      / cst24 /ipot,jv(l2),iv1,iv2,iv3,iv4,iv5/ cst801 /irchk(k2)
     *      / cst65 /vt(j9),vti(j9),jok(j9),kok(j9),index 
     *      / cst87 /delt(l2),dtol,utol,ptol 

      integer ipvt,idv,iophi,idphi,iiphi,iflg1
      double precision a,b
      common/ cst23 /a(k8,k8),b(k8),ipvt(k8),idv(k8),
     *               iophi,idphi,iiphi,iflg1

      integer iff,idss,ifug,ifyn,isyn
      common/ cst10  /iff(2),idss(h5),ifug,ifyn,isyn

      integer irct,ird
      double precision vn
      common/ cst31 /vn(k2,k7),irct,ird

      integer idcf,icfct
      common/ cst96 /idcf(k5,j9),icfct

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      integer io3,io4,io9
      common / cst41 /io3,io4,io9

      integer ids,isct,icp1,isat,io2
      common/ cst40 /ids(h5,h6),isct(h5),icp1,isat,io2

      integer isec,icopt,ifull,imsg,io3p
      common/ cst103 /isec,icopt,ifull,imsg,io3p

      integer ismax,igot
      double precision value
      common/ cst49 /value,ismax,igot 

      integer ipid,ipct
      common/ cst29 /ipid(k2,k8),ipct

      double precision vmax,vmin,dv
      common/ cst9  /vmax(l2),vmin(l2),dv(l2)
c-----------------------------------------------------------------------
c                                 initialize invariant and univariant
c                                 field counters:
      ipct = 0
      irct = 0
 
      do i = 1, k2
         irchk(i) = 0
      end do 
c                                 set starting values for search
      v(iv1) = vmin(iv1)
      v(iv2) = vmin(iv2)
      call incdp0
c                                 determine the stable chemography.
      call gall 
      call combin

      if (icopt.eq.3) then 
         iedge = 1
         if (io4.eq.0) then 
c                                 write chemographies for mixed variable
c                                 diagrams:
            call outgrf
c                                 write min and max values for the
c                                 independent potentials:
            write (n4,*) vmin(iv1),vmax(iv1)
         end if 
      else 
         iedge = 4 
      end if

      if (io3.eq.0) then 
c                                 output chemography if requested.
c                                 this option is archaic with newhld.
         call outchm  

         write (n3,1160)
         write (n3,1170) vname(iv1),vname(iv2)
c                                 write potential variable sectioning
c                                 constraints:
         if (ipot.gt.2) then
            write (n3,1180)
            write (n3,1190) (vname(jv(i)),v(jv(i)), i = 3, ipot)
         end if
c                                 write saturated and buffered
c                                 component names:
         if (isyn.eq.0) write (n3,1200) (cname(i+icp), i= 1, isat)
c                                 write dependent extensities blurb
         write (n3,1210)
         write (n3,1230) vname(ivfl)

         write (n3,1160)

         if (io3p.eq.1) then
            write (n3,1240)
            write (n3,1160)
         end if

      end if 
c                                 initialize start parms:
      do i = 1, icfct
         vt(i)  = vmin(iv1)
         vti(i) = vmin(iv2)
         kok(i) = 1
         jok(i) = 1
      end do 
      
      write (*,1000) icfct
c                                 look for reactions involving each
c                                 assemblage.
      knct = 1
40    jnct = icfct
c
      do i = knct, jnct
c                                 save the current position for chkass
         index = i 

         if (imsg.eq.0) write (*,1220) i,icfct-i
c                                 load id's into array idv:
         do j = 1, icp
            idv(j) = idcf(j,i)
         end do 
c                                 get the lower/upper decomposition
c                                 of the transpose of the concentration
c                                 matrix:
         call pivots (ier)
         if (ier.eq.1) then
c                                 error, a degenerate assemblage in
c                                 pivots due to a programming error.
            call warn (29,v(1),ier,'NEWHLD')
            cycle
         end if 
c                                 set flag for maxend:
         ismax = jok(i) 
         igot = 0
         if (ismax.eq.1) then 
            value = vmin(iv1)
         else if (ismax.eq.2) then 
            value = vmin(iv2)
         else if (ismax.eq.3) then
            value = vmax(iv1) 
         else 
            value = vmax(iv2)
         end if 
c                                 test stability along an edge of the
c                                 diagrams coordinate frame:
30       kok(i) = jok(i)
         call search (vt(i),vti(i),jok(i),iste,ivi,ivd,div,iedge,ier)
c                                 ier =2 error from search.
c                                 ier =1 stable at all conditions.
         if (ier.eq.1.or.ier.eq.2) cycle
c                                 call coface to delineate all the
c                                 equilibria that are cofacial with
c                                 the simplex defined by the vertices
c                                 idv.
         call coface (ivd,ivi,div,iste,irend)
c                                 after tracing the net attempt
c                                 to resume testing the stability
c                                 of the assemblage:
         if (igot.eq.1) then
            if (ismax.eq.1) then
c                                 resume search of travserse 1:
               vt(i) = value + 1d1*delt(iv1)
               vti(i) = vmin(iv2)
            else if (ismax.eq.2) then
c                                 resume search of travserse 2:
               jok(i) = 2
               vt(i) = value + 1d1*delt(iv2)
               vti(i) = vmax(iv1)
            else if (ismax.eq.3) then
c                                 resume search of travserse 3:
               jok(i) = 3
               vt(i) = value - 1d1*delt(iv1)
               vti(i) = vmax(iv2)
            else
c                                 resume search of travserse 4:
               jok(i) = 4
               vt(i) = value - 1d1*delt(iv2)
               vti(i) = vmin(iv1)
            end if
c                                 test for stability before search? 
            do j = 1, icp
               idv(j) = idcf(j,i)
            end do 
            call pivots (ier)
            igot = igot + 1
            goto 30 
         end if 

      end do 

      if (icfct.eq.jnct) goto 90
      knct = jnct + 1
      goto 40
c                                 return for icopt=1
90    if (icfct.eq.j9) call warn (205,r,j9,'NEWHLD')

c                                 summarize print and graphics output
      if (output) then 
         if (icopt.eq.3) call onedim
         call outier
      end if 
      
1000  format ('Initial number of',
     *        ' divariant assemblages to be tested is:',i6)
1160  format (/,80('-'),/)
1170  format ('The ',a,'-',a,' loci of (pseudo-) univariant fields'
     *        ,' follow:')
1180  format (/,'the fields are subject to the constraint(s):'/)
1190  format (25x,a,'=',g12.6)
1200  format (/,'These fields are consistent with saturation or ',
     *            'buffering constraints',/,'on the component(s):',
     *            7(a,2x))
1210  format (/,'NOTE: For each field the values of the dependent ',
     *          ' extensities are output',/,'for the first',
     *          ' equilibrium condition, in general these properties',
     *          ' vary with',/,'the independent potentials.',/)
1220  format ('Testing divariant assemblage ',i6,', ',i6,
     *        ' assemblages remaining to be tested.')
1230  format (/,'Reaction equations are written such that the high ',
     *          a,/,'assemblage is on the right of the = sign',/)
1240  format (//,'You have chosen the short form print file. ',/,
     *           'If you want equilibrium coordinates output select',
     *           /,'the long form when you run BUILD.',//)
      end

      subroutine gwash 
c-----------------------------------------------------------------------
c swash outputs thermodynamic data for pseudocompounds. 
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      character*8 names, fname*10, cname*5

      integer i,j

      common/ cst8 /names(k1)
     *      / csta7 /fname(h9)/ csta4 /cname(k5)

      double precision g
      common/ cst2 /g(k1)

      double precision cp
      common/ cst12 /cp(k5,k1)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  

      double precision vmax,vmin,dv
      common/ cst9  /vmax(l2),vmin(l2),dv(l2)

      integer ikp
      common/ cst61 /ikp(k1)
c----------------------------------------------------------------------

      p = vmin(1)
      t = vmin(2)
      xco2 = vmin(3)

      call gall

      open (n2,file='components.dat')
      write (n2,'(a)') (cname(i),i=1,icp)
      close (n2)

      open (n2,file='names.dat')
      write (n2,'(a)') (names(i),i=1,iphct)
      close (n2)

      open (n2,file='g.dat')
      do i = 1, iphct
         write (n2,*) g(i)
      end do 
      close (n2)

      open (n2,file='comp.dat')
      do i = 1, iphct
         write (n2,'(15(g15.7,1x))') (cp(j,i),j=1,icp)
      end do 
      close (n2)

      open (n2,file='solution_name.dat')
      do i = 1, iphct
         if (ikp(i).eq.0) then 
            write (n2,*) names(i)
         else 
            write (n2,*) fname(ikp(i))
         end if 
      end do 
      close (n2)

      stop

      end

      subroutine swash 
c-----------------------------------------------------------------------
c swash outputs thermodynamic data for pseudocompounds. 
c-----------------------------------------------------------------------
      implicit none

      character*8 names

      include 'perplex_parameters.h'

      integer iwarn9,iwrn48,iwrn49,jd,id,jdis,i,j,klam,ld,k,kd,h,lct

      double precision tm(m7,m6),xkd

      double precision thermo,uf,us
      common/ cst1 /thermo(k4,k10),uf(2),us(h5)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      double precision therdi,therlm      
      common/ cst203 /therdi(m8,m9),therlm(m7,m6,k9)

      integer ltyp,lmda,idis
      common/ cst204 /ltyp(k1),lmda(k1),idis(k1)

      common/ cst8 /names(k1)

      integer ipoint,imyn
      common/ cst60 /ipoint,imyn

      integer ic
      common/ cst42 /ic(k5)

      double precision cp
      common/ cst12 /cp(k5,k1)

      integer idh2o,idco2,ikind,icmpn,icout
      double precision therm,comp,tot
      common/ cst43 /therm(k4),comp(k0),tot,icout(k0),idh2o,idco2,
     *               ikind,icmpn

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  

      integer ifp
      common/ cxt32 /ifp(k1)

      integer ixp
      double precision sxs,exces
      common/ cst304 /sxs(k13),exces(m3,k1),ixp(k1)

      integer ikp
      common/ cst61 /ikp(k1)

      integer jend
      common/ cxt23 /jend(h9,k12)
c----------------------------------------------------------------------

      open (n2,file='swash.dat')

      iwarn9 = 0
      iwrn49 = 0 
      iwrn48 = 0
      
      do 10 jd = ipoint + 1, iphct
c                              skip fluid cpds
         if (ifp(jd).eq.1) iwrn48 = 1
c                              create mock-entry iphct + 1 for output
c                              initialize paramaters
         id = k1
         idis(k1) = 0
         ltyp(k1) = 0
         jdis = 0
         klam = 0
         names(k1) = names(jd)  
 
         do i = 1, k4
            therm(i) = 0d0
         end do 

         do i = 1, k0
            comp(i) = 0d0  
         end do 
c                                 now combine endmember props

         do j = 1, jend(ikp(jd),1)
c                                 endmember index
            kd = jend(ikp(jd),2+j)
            xkd = sxs(ixp(jd)+j) 
c                                 test for mughnahan EoS
            if (thermo(18,kd).ne.0d0) then
               iwrn49 = 1 
               goto 10
            end if             

            do i = 1, 18
               therm(i) = therm(i) + xkd * thermo(i,kd)
            end do 

            do i = 1, icomp
               comp(ic(i)) = comp(ic(i)) + xkd * cp(i,kd)
            end do 

            if (idis(kd).ne.0) then

               idis(k1) = m9
               jdis = jdis + 1
               ld = idis(kd) 
               if (jdis.gt.1) goto 91
               do i = 1, 7      
                  therdi(i,m9) = xkd * therdi(i,ld)
               end do 
               therdi(8,m9) = therdi(8,ld)
               therdi(9,m9) = therdi(9,ld)

            end if

            if (ltyp(k1).ne.0) then

               ltyp(k1) = ltyp(kd)
               lmda(k1) = k9
               ld = lmda(kd)
               klam = klam + 1

               if (ltyp(kd).eq.10) then
                  if (klam.gt.3) goto 91
                  ltyp(k1) = 9 + klam
                  therlm(1,klam,k9) = therlm(1,1,ld)
                  therlm(2,klam,k9) = xkd * therlm(2,1,ld)
               else if (ltyp(kd).lt.4) then
                  if (klam.gt.1) goto 91
                  do i = 1, ltyp(kd)
                     therlm(1,i,k9) = xkd * therlm(1,i,ld)
                     therlm(2,i,k9) = xkd * therlm(2,i,ld)
                     therlm(5,i,k9) = xkd * therlm(5,i,ld)
                     therlm(6,i,k9) = xkd * therlm(6,i,ld)
                     therlm(3,i,k9) = therlm(3,i,ld)
                     therlm(4,i,k9) = therlm(4,i,ld)
                     therlm(7,i,k9) = therlm(7,i,ld)
                     therlm(8,i,k9) = therlm(8,i,ld)
                  end do
               else if (ltyp(kd).lt.8) then

                  if (klam.gt.1) goto 91

                  do k = 1, ltyp(kd) - 3

                     therlm(1,k,k9) = therlm(1,k,ld)
                     therlm(2,k,k9) = therlm(2,k,ld)

                     do h = 3, 12
                        therlm(h,k,k9) = xkd * therlm(h,k,ld) 
                     end do
                  end do 
               end if
            end if
         end do

         call unlam (tm,id,lct)

         call unver (therm(1),therm(2),therm(3),therm(4),
     *               therm(5),therm(6),therm(7),therm(8),
     *               therm(9),therm(10),therm(11),
     *               therm(12),therm(13),therm(14),
     *               therm(15),therm(16),therm(17),
     *               therm(18),tr,pr)
c                                 output the data 
         write (n2,1100) names(k1),1,ikp(jd),ltyp(k1),idis(k1)
         write (n2,1080) (comp(i),i=1,icmpn)
         write (n2,1090) (therm(i),i=1,18)
         call outlam (tm,idis(k1),lct)

         goto 10 

91       iwarn9 = 1

10    continue  

      if (iwrn49.eq.1) call warn (49,thermo(18,kd),iphct,names(kd))

      if (iwrn48.eq.1) call warn (48,thermo(18,kd),iphct,names(kd))

      if (iwarn9.eq.1) call warn (9,t,i,names(k1))

      close (n2)

      stop

1080  format (12(f7.4,1x))
1090  format (5(g13.7,1x))
1100  format (a8,i2,i2,i2,i2)

      end

      subroutine abload (*)
c---------------------------------------------------------------------
c abload assembles the matrix 'a' and vector 'b' for icp component
c sytems and then solves the equation ax = b, the vector x is
c returned in 'b'.
 
c   referenced by: combin
c   references to: subst,factor
c----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i,j,ier
 
      double precision a,b
      integer ipvt,idv,iophi,idphi,iiphi,iflg1
      common/ cst23  /a(k8,k8),b(k8),ipvt(k8),idv(k8),iophi,idphi,
     *                iiphi,iflg1

      double precision g
      common/ cst2  /g(k1)

      double precision cp
      common/ cst12 /cp(k5,k1)
 
      integer hcp,id
      common/ cst52 /hcp,id(k7)
c-----------------------------------------------------------------------
      do i = 1, hcp
         do j = 1, hcp
            a(i,j) = cp(j,id(i))
         end do
      end do

      call factor (a,hcp,ipvt,ier)
 
      if (ier.eq.1) goto 99
 
      do i = 1,  hcp
         b(i) = g(id(i))
      end do 
 
      call subst (a,ipvt,hcp,b,ier)
c                                 don't use the ier flag from 
c                                 subst 'cause it only became
c                                 necessary with "reopt". (not
c                                 great logic).
      return
99    return 1
      end

      subroutine asschk
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i,j,ier

      double precision delt,dtol,utol,ptol,gphi,dg

      common/ cst87 /delt(l2),dtol,utol,ptol

      double precision g
      common/ cst2 /g(k1)

      double precision cp
      common/ cst12 /cp(k5,k1)

      integer iflag
      common/ cst7 /iflag
 
      double precision a,b
      integer ipvt,idv,iophi,idphi,iiphi,iflg1
      common/ cst23  /a(k8,k8),b(k8),ipvt(k8),idv(k8),iophi,idphi,
     *                iiphi,iflg1

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  
c-----------------------------------------------------------------------
c                                 determine chemical potentials
      iflag = 0
      do i = 1, icp
         b(i) = g(idv(i))
      end do
c
      call subst (a,ipvt,icp,b,ier)
c                                 test phases against the
c                                 assemblage idv
      do 20 i = istct, iphct
         gphi = 0d0
         do j = 1, icp
            gphi = gphi + cp(j,i) * b(j)
         end do 
         dg = g(i) - gphi
         if (dg.gt.dtol) goto 20
c                                check that a phase is not metastable
c                                with respect to itself, this
c                                could be remedied by changing the 
c                                value of dtol.
         do j = 1, icp
            if (idv(j).eq.i) goto 20
         end do 
c                                assemblage is metastable with respect
c                                to phase idphi.
         iflag = iflag + 1
         idphi = i
c                                this would force a refinement
c                                of condtions by the calling
c                                routine. since utol << dtol
c        if (dabs(dg).gt.utol) iflag = iflag + 1
         if (iflag.gt.1) goto 99
20    continue

99    end

      subroutine assir (bad)
c----------------------------------------------------------------------
c assir is called only for phase diagram projection calculations. 

c irct - counter of univariant rxns.
c ivarrx(irct) - variance flag of the iirct(th) rxn (set in balanc as
c                 ivar in cst61.
c irchk(irct) - a flag indicating if the iirct(th) rxn occurs on the
c                edge of a diagram (initialized in newhld and set in
c                coface or sfol1).
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,j,k,irnms,isol,np

      logical solvs1, bad

      common/ csta1 /irnms(k2,k7)

      integer idr,ivct
      double precision vnu
      common/ cst25 /vnu(k7),idr(k7),ivct

      integer irct,ird
      double precision vn
      common/ cst31 /vn(k2,k7),irct,ird

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      integer irv
      common/ cst35 /irv(k2)

      integer ikp
      common/ cst61 /ikp(k1)

      integer ivarrx,ivarip,isudo,ivar
      common/ cst62 /ivarrx(k2),ivarip(k2),isudo,ivar
c-----------------------------------------------------------------------
      bad = .false.
c                                 check if the reaction is known
      do 70 i = 1, irct

         if (ivct.ne.irv(i)) cycle

         do 80 j = 1, ivct
            do k = 1, ivct
               if (irnms(i,k).eq.idr(j)) goto 80
            end do 
c                                 no match with irnms(i,k) try i + 1.
         goto 70

80       continue
c                                 match found with irnms(i,k)
         ird = i
         ivar = ivarrx(i)

         goto 99

70    continue
c                                 the reaction is new:
      irct = irct + 1
      isol = 0 
      ivar = 1
      ird = irct
      if (irct.gt.k2) call error (181,r,k2,'ASSIR')
      irv(irct) = ivct

      do i = 1, ivct
         if (ikp(idr(i)).gt.0) isol = isol + 1
      end do 

      if (isol.gt.0) then 
c                                 compositional limits
         call sollm0 (k7,ivct,idr)

         np = ivct 
c                                 classify the reaction
         if (isol.gt.1) call miscb0 (k7,ivct,np,solvs1,idr)
c                                 set ivar, just in case it's used someplace
         ivar = 1 + ivct - np

      end if 

      if (ivar.gt.isudo) then 
c                                 throw out reactions with variance > isudo
         bad = .true.
         irct = irct - 1
         goto 99

      end if 

      ivarrx(irct) = ivar

      do i = 1, ivct
         vn(irct,i) = vnu(i)
         irnms(irct,i) = idr(i)
      end do 

99    end

      subroutine assri (ier)
c----------------------------------------------------------------------
c assri is a special version of assir called only for mixed variable diagrams

c irct         - reaction counter
c irv(irct)    - number of phases in the irct(th) rxn
c                (set in balanc as ivct cst25).
c ivarrx(irct) - variance flag of the iirct(th) rxn (set in balanc as
c                 ivar in cst61.
c irchk(irct)  - a flag indicating if the iirct(th) rxn occurs on the
c                edge of a diagram (initialized in newhld and set in
c                coface or sfol1).
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,j,k,ier,np,isol

      logical solvs1

      integer irct,ird
      double precision vn
      common/ cst31 /vn(k2,k7),irct,ird

      double precision delt,dtol,utol,ptol
      common/ cst87 /delt(l2),dtol,utol,ptol

      integer idr,ivct
      double precision vnu
      common/ cst25 /vnu(k7),idr(k7),ivct

      integer irnms
      common/ csta1 /irnms(k2,k7)
 
      double precision vip
      common/ cst28 /vip(l2,k2)

      integer ipot,jv,iv1,iv2,iv3,iv4,iv5
      common/ cst24 /ipot,jv(l2),iv1,iv2,iv3,iv4,iv5

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      integer irv
      common/ cst35 /irv(k2)

      integer ivarrx,ivarip,isudo,ivar
      common/ cst62 /ivarrx(k2),ivarip(k2),isudo,ivar

      integer ikp
      common/ cst61 /ikp(k1)
c-----------------------------------------------------------------------
c                         check if the reaction is known
      do 70 i = 1, irct
         if (ivct.ne.irv(i)) goto 70
         do 80 j = 1, ivct
            do k = 1, ivct
               if (irnms(i,k).eq.idr(j)) goto 80
            end do 
c                         no match with irnms(i,k) try i + 1.
            goto 70
80       continue
c                         match found with irnms(i,k), now check
c                         that the conditions are different, if
c                         ier = 1, univeq failed, skip the check

         ivar = ivarrx(i)

         if (ier.eq.0) then 
            if (dabs(vip(iv1,i)-v(iv1)).le.delt(iv1)) then 
c                         since there mat be multiple occurences
c                         of the same reaction, go through the entire
               ird = i
               ier = 1
               goto 99 
            end if 
         end if 

70    continue
c                                 the reaction is new:
      irct = irct + 1
      ier = 0
      isol = 0 
      ivar = 0
      ird = irct
      if (irct.gt.k2) call error (181,r,k2,'ASSIR')
      irv(irct) = ivct

      do i = 1, ivct
         if (ikp(idr(i)).gt.0) isol = isol + 1
      end do 

      if (isol.gt.0) then 
c                                 compositional limits
         call sollm0 (k7,ivct,idr)

         np = ivct
c                                 classify the reaction
         if (isol.gt.1) call miscb0 (k7,ivct,np,solvs1,idr)
c                                 set ivar, just in case it's used someplace
         ivar = ivct - np

      end if 

      ivarrx(irct) = ivar

      do i = 1, ivct
         vn(irct,i) = vnu(i)
         irnms(irct,i) = idr(i)
      end do 

99    end

      subroutine balanc (b,idv,idphi,ier)
c-----------------------------------------------------------------------
c balanc balances reactions between the phases in idv and idphi
c if the system is saturated with any phases these are taken
c into account.  once the reaction is determined balanc determines
c the thermodynamic parameters needed for subroutine univeq and
c saves them in the array 'rxn'.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ipvt(k8),idv(k8),i,j,ip,idphi,ier

      double precision a(k8,k8),b(k8)

      integer iffr,isr
      double precision vuf,vus
      common/ cst201 /vuf(2),vus(h5),iffr,isr

      double precision cp
      common/ cst12 /cp(k5,k1)

      integer iff,idss,ifug,ifyn,isyn
      common/ cst10  /iff(2),idss(h5),ifug,ifyn,isyn

      integer idr,ivct
      double precision vnu
      common/ cst25 /vnu(k7),idr(k7),ivct

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  

      integer ids,isct,icp1,isat,io2
      common/ cst40 /ids(h5,h6),isct(h5),icp1,isat,io2

      integer ivarrx,ivarip,isudo,ivar
      common/ cst62 /ivarrx(k2),ivarip(k2),isudo,ivar
c-----------------------------------------------------------------------
c                                 assemble transpose of the composition
c                                 matrix for an assemblage of icp phases
      ivar = 1
                                                
      do i = 1, icp
         do j = 1, icp
            a(j,i) = cp(j,idv(i))
         end do 
         idr(i) = idv(i)
      end do 
c                                 assemble the composition vectore
c                                 for the icp+1 th phase tangent to
c                                 the free energy plane.
      ip = icp + 1
      do i = 1, icp
         b(i) = cp(i,idphi)
      end do 
      idr(ip) = idphi
      b(ip) = -1d0
c                                 balance the equilibrium
      call factor (a,icp,ipvt,ier)
      if (ier.eq.1) goto 999

      call subst (a,ipvt,icp,b,ier)
c                                 don't use the ier flag from 
c                                 subst 'cause it only became
c                                 necessary with "reopt". (not
c                                 great logic).

c                                 eliminate phases with vnu= 0
      ivct = 0

      do i = 1, ip

         if (dabs(b(i)).gt.1d-05) then 
            ivct = ivct + 1
            vnu(ivct) = b(i)
            idr(ivct) = idr(i)
         end if 

      end do 

      if (isyn.eq.1) goto 40
c                                determine stoichiometric coefficients
c                                of saturated components:
c                                (these aren't used for anything real)
      isr = 1
      do j = 1, isat
         vus(j) = 0d0
         do i = 1, ivct
            vus(j) = vus(j) + vnu(i) * cp(icp+j,idr(i))
         end do 
c                                flag isr is not used.
         if (vus(j).ne.0d0) isr = 0   
      end do 
c                                ok, now vus has the total deltas
   
c                                determine stoichiometric coefficients
c                                of saturated phase components:
      iffr = 1
40    do j = 1, 2
         vuf(j) = 0d0
         if (iff(j).ne.0d0) then 
            do i = 1, ivct
               vuf(j) = vuf(j) + vnu(i) * cp(iff(j),idr(i))
            end do 
c                                this is only the total of saturated
c                                phase components for the reaction
c                                (doesn't include saturated phases).
c                                flag iffr is used
            if (vuf(j).ne.0d0) iffr = 0
         end if
      end do 
c                                et finis
999   end

      subroutine combin
c----------------------------------------------------------------------------
c combinatorial algorithm for mfes. 

c for constrained bulk compositions (jbulk.ge.icp)  
c----------------------------------------------------------------------------
      implicit none

      integer ipoint,imyn,jcp,iclose,istart,jstart,i,j,k,l,m,igen,itic,
     *        idno,lastk,iend,icase

      double precision sum,dgphc,sign

      include 'perplex_parameters.h'

      integer hcp,id
      common/ cst52 /hcp,id(k7)/ cst60 /ipoint,imyn

      double precision cp
      common/ cst12 /cp(k5,k1)

      integer idcf,icfct
      common/ cst96 /idcf(k5,j9),icfct

      integer jtest,jpot
      common/ debug /jtest,jpot

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  

      integer jbulk
      double precision cblk
      common/ cst300 /cblk(k5),jbulk
c----------------------------------------------------------------------
c                                 initialization
      imyn = 1
      icfct = 0
      jcp = icp - 1

      iclose = 0
c                                 output thermo data 
      if (jtest.eq.2) call outdat
c                                 identify the seed stable assemblage
      call findas 
c                                 save the assemblage
      call assdc (icase)
c                                 finished if one component
      if (icp.eq.1) goto 99
c                                 now use the c-1 phase simplicial
c                                 facets generated from the seed
c                                 assemblage to identify the remaining
c                                 assemblages.
      istart = 1
      jstart = 0
      igen = 1

40    iend = icfct 
      write (*,*) 'cycle ',igen,istart,iend
      igen = igen + 1

      do 50 i = istart, iend
c                             each facet may generate icp-1 new psf's
         do 60 j = 1, jcp
c                             set vertices identities:
            itic = 0
c                             this loop generates the psf consisting of
c                             the icp-1 phases, saving the icp'th phase
            do k = 1, icp 
               if (k.eq.j) then 
                  idno = idcf(k,i) 
               else 
                  itic = itic + 1
                  id(itic) = idcf(k,i)
               end if 
            end do
c                              check that the assemblage is not degenerate
c                              a valid psf must span the icp dimensional 
c                              composition space, count the components:
            do k = 1, icp
               sum = 0d0
               do l = 1, jcp
c                             sum the amounts of component k in the psf 
                  sum = sum + cp(k,id(l))
               end do 
c                             if the component is 0, the psf is degenerate
               if (sum.lt.1d-5) goto 60
            end do
c                             next check if the psf matches any permutation
c                             of icp-1 phases in a facet tested earlier:
c                             this test is not essential, but it should
c                             save time. 
c           do 80 k = 1, i-1
            do 80 k = jstart+1, i-1
c                             unless there is some weird hyperdimensional
c                             effect that i don't understand, the index k
c                             should only need to run over the facets identified
c                             in the last cycle up to i-1. but tests seem to
c                             give more assemblages if k runs over all the
c                             assemblages, a problem to be investigated.
               do 70 l = 1, jcp
                  do m = 1, icp
                     if (idcf(m,k).eq.id(l)) goto 70
                  end do
c                             no match with assemblage k:
                  goto 80
70             continue
c                             matched assemblage k, reject the psf:
               goto 60
80          continue
c                             now identify the first potential vertice:

            do 10 k = istct, iphct
c                             don't allow phases of the original simplex.
               do l = 1, icp
                  if (k.eq.idcf(l,i)) goto 10
               end do 
c                             try to find a phase that is stable
c                             with respect to phase idno:
               id(icp) = k

               call abload (*10)

               if (dgphc(idno).gt.-1d-05) then
c                              the plane is stable
                  lastk = k + 1
                  goto 20   
               end if 
10          continue
c                               if no acceptable phase is found
c                               the assemblage must lie on an edge
c                               of the compositional space.
            if (iclose.eq.0) call warn (28,sign,1,'COMBIN')
            iclose = 1
            goto 60
c                               now do thermodynamic testing
20          if (lastk.gt.iphct) then
               call assdc (icase)
               goto 60
            end if          

            do k = lastk, iphct
c                                 checkd replaces the trial vertice
c                                 with phase j, if the assemblage is 
c                                 found to be metastable with respect 
c                                 to phase j. given the speed of this
c                                 test, i doubt geometric testing would
c                                 be worthwhile here.
c                                 don't allow return of the vertice
c                                 idpsf(icp,i)
               if (k.ne.idno) call checkd (k)
            end do 
c                                 save the new assemblage:
            call assdc (icase)

60       continue
50    continue 

      if (icfct.gt.iend) then
c                                 reset counters for next cycle
            jstart = istart
            istart = iend + 1
            goto 40

      end if
c                                 check that something was found.
      if (icfct.eq.0) then
         if (iclose.eq.1.and.jbulk.ne.0) call error (42,sign,2,'COMBIN')
         call error (999,sign,2,'COMBIN')
      end if 

99    end 

      subroutine assdc (icase)
c-----------------------------------------------------------------------
c assdc counts (icfct) and assigns (idcf) c-component facets. the
c vertices are identified by the array id(icp).

c   icase = 0 => no new assemblage counted.
c   icase = 1 => new assemblage.

c   referenced by: findas
c   references to: outlm0, miscib, error
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i,j,k,icase,ivar,isol,np

      logical solvs1
 
      integer jasmbl,iasmbl
      common/ cst27 /jasmbl(j9),iasmbl(j9)

      integer ht,id
      common/ cst52 /ht,id(k7)

      integer idcf,icfct
      common/ cst96 /idcf(k5,j9),icfct 

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  

      integer jtest,jpot
      common/ debug /jtest,jpot

      integer ikp
      common/ cst61 /ikp(k1)
c-----------------------------------------------------------------------
      icase = 0
c                             test for equivalence with earlier
c                             assemblage:
      do 10 i = 1, icfct
         do 20 j = 1, icp
            do k = 1, icp
              if (idcf(j,i).eq.id(k)) goto 20
            end do 
            goto 10
20       continue 
c                             matches earlier assemblage, return.
         return
10    continue
c                             unique (and bounding) assemblage:   
      icfct = icfct + 1
      if (icfct.gt.j9) call error (204,0d0,j9,'ASSDC')
      jasmbl(icfct) = 0
      icase = 1
c                             assign the assemblage:
      isol = 0  
      do i = 1, icp
         idcf(i,icfct) = id(i)
         if (ikp(id(i)).gt.0) isol = isol + 1
      end do 

      ivar = 0

      if (isol.gt.0) then 
c                             compositional limits
         call sollm0 (k5,icp,id)
c                             do miscibility check:
         if (isol.gt.1) then

            call  miscb0 (k5,icp,np,solvs1,id)
            ivar = icp - np

         end if 
      end if 

      iasmbl(icfct) = ivar
c                             print potentials:
      if (jpot.eq.0) call prtpot

      end

      subroutine sollm0 (idim,ntot,id)
c----------------------------------------------------------------------
c subroutine to extract compositional range of endmembers from stable 
c fixed pseudocompounds (icopt=0-3, cf. sollim) for auto_refine option.

c referenced by: assdc.
c references to:

c----------------------------------------------------------------------
      implicit none 

      include 'perplex_parameters.h'
c                                 -------------------------------------
c                                 local variables:
      integer jd,ids,i,j,k,idim,ntot,id(idim)
c                                 -------------------------------------
c                                 global variables:
c                                 working arrays
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

      integer ikp
      common/ cst61 /ikp(k1)
c----------------------------------------------------------------------
      do k = 1, ntot
         jd = id(k)
         ids = ikp(jd)
         if (ids.le.0) cycle
c                                 set stable flag
         stable(ids) = .true.
c                                 get composition
         call getolx (ids,jd)
c                                 check x-ranges
         do i = 1, istg(ids)
            do j = 1, ispg(ids,i)-1
c                                 low limit:
               if (x(i,j).lt.xlo(j,i,ids)) then
                  xlo(j,i,ids) = x(i,j)
c                                 check if solution is at an unnatural limit
                  if (x(i,j).gt.xmno(ids,i,j).and.
     *                x(i,j).le.xmng(ids,i,j)
     *               .and..not.limit(ids)) then
                     write (*,1000) fname(ids),x(i,j),i,j
                     limit(ids) = .true.
                  end if 
               end if 
c                                 high limit:
               if (x(i,j).gt.xhi(j,i,ids)) then
                  xhi(j,i,ids) = x(i,j)
c                                 check if solution is at an unnatural limit
                  if (x(i,j).lt.xmxg(ids,i,j).and.
     *                x(i,j).ge.xmxg(ids,i,j)
     *               .and..not.limit(ids)) then
                     write (*,1000) fname(ids),x(i,j),i,j
                     limit(ids) = .true.
                  end if 
               end if 
            end do 
         end do
      end do   

1000  format (/,'WARNING: composition of solution ',a,' has reached an',
     *          ' internal limit (',f5.3,')',/,
     *          'on site ',i1,' for species ',i2,'. If this warning'
     *         ,' occurs during auto-refinement the problem',/,'can ',
     *          'be circumvented by increasing auto_refine_slop ',
     *          'in perplex_option.dat or editing the ',/,
     *          'auto_refine_XXX file.',/)

      end 

      subroutine chkass (jdv,ivd,ivi,iste)
c-----------------------------------------------------------------------
c chkass determines if a divariant assemblage (jdv) generated by
c newass has already been identified. there are three possible cases:
c in the reaction as determined by routine balanc.
c (1) the assemblage is unique: the indices of the phases, the initial
c     conditions of stability, and the coordinate frame edge are recorde
c (2) the assemblage has already been defined and it's stability field
c     determined: quit
c (3) the assemblage has already been defined but it's stability field
c     has not been determined: reset the initial conditions of stability
c     to shorten traverses by the routine search, if possible. 

c input:  jdv(icp) - a vector containing the indices of phases in a
c                    divariant assemblage.
c         ipart    - 1 if testing original divariant assemblages, 2 if
c                    testing new divariant assemblages.
c         index    - index of the divariant assemblage in array idcf
c                    (ipart=1) or idns (ipart=2) currently being tested.
  
c referenced by: newass
c references to: no external references
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer jdv(k8), icrap, j, k, 
     *        i,jok,kok,index,iste,ivd,ivi

      double precision vt,vti

      common/ cst65 /vt(j9),vti(j9),jok(j9),kok(j9),index    

      integer idcf,icfct
      common/ cst96 /idcf(k5,j9),icfct

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  

      double precision vmax,vmin,dv
      common/ cst9  /vmax(l2),vmin(l2),dv(l2)

      save icrap

      data icrap /0/
c-----------------------------------------------------------------------

      do 40 i = 1, icfct
         do 50 j = 1, icp
            do k = 1, icp
               if (jdv(k).eq.idcf(j,i)) goto 50
            end do 
c                                 no match with assemblage i
            goto 40
c                                 
50       continue
c                                 first check if within known stability
c                                 field of the assemblage
         if (iste.lt.kok(i)) then 
            goto 99
         else if (iste.eq.kok(i)) then 
            if (iste.lt.3) then
               if (vt(i).ge.v(ivd)-1d-2*dv(ivd)) goto 99
            else 
               if (vt(i).le.v(ivd)+1d-2*dv(ivd)) goto 99
            end if 
c                                 it extends known stability field:
            if (i.gt.index) then
c                                 the assemblage hasn't been tested
c                                 then reset start conditions 
               call sreset (jok(i),iste,vt(i),v(ivd),vti(i),v(ivi))
               goto 99
            else 
c                                 the assemblage has been tested, so 
c                                 skip the check and look for another
c                                 occurence, if no other occurence the
c                                 assemblage will be saved as new. 
               goto 40

            end if 

         else if (i.gt.index) then 
            call sreset (jok(i),iste,vt(i),v(ivd),vti(i),v(ivi))
            goto 99
         end if
c                                 otherwise save as a new assemblage
         goto 100
  
40    continue
c                                 the assemblage is unique, store
c                                 parameters:
100   if (icrap.ne.1) then 

         icfct = icfct + 1

         if (icfct.gt.j9) then
            call warn (205,r,j9,'CHKASS')
            icfct = j9
            icrap = 1
         end if 

         do i = 1, icp
            idcf(i,icfct) = jdv(i)
         end do 

         vt(icfct) = v(ivd)
         vti(icfct) = v(ivi)
         jok(icfct) = iste
         kok(icfct) = iste

      end if 

99    end

      subroutine coface (iovd,iovi,odiv,iste,irend)
c----------------------------------------------------------------------
c once nohold has established that a stable equilibrium occurs on
c the boundary of a diagram, coface traces the equilibrium condiitons
c within the diagram.  if the equilibrium terminates at an
c invariant point coface traces the equilibria which emanate
c from the invariant point and any invariant points encountered
c subsequently.  coface does not detect indifferent points.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer iflag,irchk,iopct,icter,irend,iovd,
     *        ivi,ivd,iovi,jflg,ier,iste,i,inpct,jer,ikwk

      logical bad

      double precision odiv,div

      common/ cst7   /iflag/ cst801 /irchk(k2)

      integer ipid,ipct
      common/ cst29 /ipid(k2,k8),ipct
 
      double precision vip
      common/ cst28 /vip(l2,k2)

      integer ipvt,idv,iophi,idphi,iiphi,iflg1
      double precision a,b
      common/ cst23 /a(k8,k8),b(k8),ipvt(k8),idv(k8),
     *               iophi,idphi,iiphi,iflg1

      integer ipt2
      double precision ptx
      common/ cst32 /ptx(l5),ipt2

      integer ibug
      common/ cst105 /ibug(k2)

      integer irct,ird
      double precision vn
      common/ cst31 /vn(k2,k7),irct,ird

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      integer isec,icopt,ifull,imsg,io3p
      common/ cst103 /isec,icopt,ifull,imsg,io3p
c-----------------------------------------------------------------------
      ivi = iovi                                                  
      ivd = iovd
      jflg = 1

      call balanc (b,idv,idphi,ier)
c                                 singular concentration matrix
      if (ier.eq.1) call error (68,v(1),ier,'COFACE')
c                                 call newass to define new (pseudo-)
c                                 divariant assemblages generated by
c                                 the reaction defined by balanc:
      call newass (b,idv,idphi,ivd,ivi,iste)
c                                 call assir to determine if the
c                                 reaction has already been identified
c                                 and if so quit, there is a small
c                                 possibility that as a result of
c                                 quitting at this stage some equilibria
c                                 will not be identified, if this is a
c                                 significant concern the following
c                                 computed goto should be deleted.
      if (icopt.ne.1) then  
c                                 allow duplicated reactions in 
c                                 mixed variable diagrams:
         call univeq (ivd,ier)
c                                 assign reaction:
         call assri (ier)

         if (ier.eq.1) goto 9999

         call delrxn 

         vip(1,ird) = v(1)
         vip(2,ird) = v(2)
         vip(3,ird) = v(3)
         vip(4,ird) = v(4)
         vip(5,ird) = v(5)

         call outrxn (0,0)

         goto 9999

      end if 

      call assir (bad)
c                                 skip curves of variance > isudo
      if (bad) goto 9999
c                                 skip equilibria that have already
c                                 been defined on an edge, there is some
c                                 possibility that this may result in an
c                                 an incomplete diagram, but it saves
c                                 some time and redundancy in output.
c                                 check if the univariant rxn has alread
c                                 been found on the edge of a diagram.
c                                 set flags irchk and ier:
c      if (irchk(ird).eq.1) goto 9999
c                                 call svrend to see if the end point 
c                                 has already been found, if it has ier
c                                 is 1, and traverse can be skipped. 
      call svrend (ird,irend,ier)
      if (ier.eq.1) goto 9999

      irchk(ird) = 1
      call univeq (ivd,ier)
c                                 get deltas for dependent extensities:
      call delrxn
c                                 if univeq fails on a bounding edge
c                                 write error message and return:
      if (ier.eq.1.or.ier.eq.2) goto 9000
c                                 save the index of the c+1th phase
c                                 as iophi.
      iophi =idphi
c                                 factor the concentration matrix:
      call pivots (ier)
      if (ier.eq.1) call error (69,v(1),ier,'COFACE')
c                                 set the increment for the iv
      div = odiv
c                                 initialize counters
      ipt2 = 0
      icter = 0
c                                 assign the 1st point
      call assptx
c                                 follow the equilibrium
60    ikwk = 0 
      jflg = 0
      call sfol1 (ivd,ivi,ier,div,ikwk,irend)
c                                 sfol1 returns ier = 1
      if (ier.eq.1.or.ier.eq.2) goto 70

      ivi = iovi
      ivd = iovd

      if (iflag.eq.0) then 
c                                 iflag = 0 the univariant
c                                 curve was traced from one
c                                 edge to another, reset
c                                 the starting conditions
c                                 and set flag jer for newhld
         goto 9999
      else if (iflag.eq.1) then
c                                 an ip was found:
         goto 10
      end if 

70    call switch (div,ivi,ivd,jer)
      if (jer.eq.1) goto 75

      icter = icter + 1
      if (icter.lt.4) goto 60
      
75    call warn (10,v(1),ier,'COFACE')

      call outrxn (ipct,2)
      ibug(irct) = 1
      ivi =iovi
      ivd=iovd
c                                 return on error
      goto 9999                
c                                 a new invariant point has been
c                                 encountered, generate the new 
c                                 univariant fields.
c                           iopct=the total number of ip's encountered
c                                 before each call to sfol2.
c                           inpct=the total number after each call to 
c                                 sfol2.
10    iopct = ipct
      call sfol2 (iovi,iovd,iopct,irend)
      if (iopct.eq.ipct) goto 9999
      iopct = iopct + 1
      inpct = ipct
c                                 loop for new invariant points
30    do i = iopct, inpct
         call sfol2 (iovi,iovd,i,irend)
      end do 
      if (ipct.eq.inpct) goto 9999
      iopct = inpct + 1
      inpct = ipct
      goto 30
c                                 error in univeq:
9000  call warn (79,v(1),ird,'COFACE')
      ipt2 = 0
      call outrxn (ipct,2)
      ibug(irct) = 1

9999  iflg1= 0 

      if (jflg.eq.1.or.io3p.eq.1) goto 999

      write (n3,1020) 
      
1020  format ('Network traced, resuming boundary search.',/)

999   end

      subroutine delvar (dv,iflag,iflg1)
c-----------------------------------------------------------------------
c delvar determines if a variable increment 'dv' should
c be incremented or decremented on a traverse to locate
c to locate conditions where an equilibrium is metastable with respect
c to exactly one phase (iflag = 1). the flag iflg1 indicates if the
c equilibrium was stable (iflg1 = 0) or metastable (iflg1 = 1) at the
c last conditions tested.
c-----------------------------------------------------------------------
      implicit none
  
      integer iflag,iflg1

      double precision dv

c                            the equilibrium is stable (iflag = 0):
      if (iflag.eq.1.or.iflag.eq.2) then 
c                            the equilibrium is metastable (iflag = 2):
        if (iflg1.ne.1) then
c                            and was stable before (iflg1 = 0) then
c                            switch sign and halve dv:
            dv = -dv/2d0
            iflg1 = 1
         end if 

      else 
c                            and was stable before (iflg1 = 0) then
c                            continue looking in the same direction
         if (iflg1.eq.1) then
c                            and was metastable before (iflg1 = 1) then
c                            switch sign and halve dv:
            dv = -dv/2d0
            iflg1 = 0
         end if 
      end if 

      end

      subroutine flipit (ddv,vst,ivd,ist,inow,jer)
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      character*8 names

      integer iflag,ivd,ist,inow,jer,j

      double precision v,x,vst,ddv

      common/ cst5 /v(l2),x(4)/ cst7 /iflag/ cst8 /names(k1)

      integer ipvt,idv,iophi,idphi,iiphi,iflg1
      double precision a,b
      common/ cst23 /a(k8,k8),b(k8),ipvt(k8),idv(k8),
     *               iophi,idphi,iiphi,iflg1

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  

      double precision vmax,vmin,dv
      common/ cst9  /vmax(l2),vmin(l2),dv(l2)
c-----------------------------------------------------------------------
c                                 initialization
      call gall 
      call asschk
c                                 the goto 20 is a dummy to avoid
c                                 compiler errors, change 90 to 20 to
c                                 get backup search, ha ha
      if (iflag.ne.0.and.v(ivd).eq.vst.and.ist.eq.inow) goto 90

      goto (10,20,20),iflag
c                                 the assemblage is stable, return
      jer = 0
      return
c                                 the assemblage is metastable with
c                                 respect to only one phase:
10    jer = 1
      return
c                                 the assemblage is metastable with
c                                 respect to >1 phase, therefore
c                                 it must be stable at some condition
c                                 between the present condition and the
c                                 conditions of the univariant
c                                 equilibrium that generated the
c                                 assemblage. flipit does a reverse
c                                 traverse to find this condition.
c                                 initialize:
20    iflg1 = 1
      ddv = -ddv
c                                 begin search:
100   v(ivd) = v(ivd) + ddv
      call incdep (ivd) 
c                                 out of range?:
      if (ist.lt.3) then 
         if (v(ivd).lt.vmin(ivd)) v(ivd) = vmin(ivd)
         if (v(ivd).lt.vst) goto 130
         ddv = -dabs(ddv)/2d0
         v(ivd) = vst
         call incdep (ivd) 
         iflg1 = 0
         goto 100
      else 
         if (v(ivd).gt.vmax(ivd)) v(ivd) = vmax(ivd)
         if (v(ivd).gt.vst) goto 130
         ddv = dabs(ddv)/2d0
         v(ivd) = vst
         call incdep (ivd) 
         iflg1= 0
         goto 100
      end if 
c                                 calculate phase energies:
130   call gall 
c                                 test stability of the assemblage:
      call asschk
c                                 check if search is in range:
      if (ist.lt.3) then 
         if (v(ivd).le.vmin(ivd).and.iflag.gt.0) goto 60
      else
         if (v(ivd).ge.vmax(ivd).and.iflag.gt.0) goto 60
      end if 
c                                 iflag=1, found a stable equilibrium
c                                 in bounds:
      if (iflag.eq.1) write (n3,*) 'flipit worked please tell me!'
c
      if (iflag.eq.1) goto 10
c                                 refine increment if necessary:
      call delvar (ddv,iflag,iflg1)
c                                 check if increment has been refined
c                                 beyond the tolerance which is
c                                 arbitrarily set at 1.d-8.
      if (dabs(ddv).lt.1d-8) goto 9000
c                                 increment conditions again:
      goto 100
c                                 next traverse.
60    call warn (76,vst,ist,'FLIPIT')

      goto 90
c                                 write warning message:
9000  write (n3,1010) ivd,dv(ivd),(names(idv(j)),j = 1, icp)
90    jer = 2
c                                 done:
1010  format (/,'**warning ver045** FLIPIT: > 1 equilibrium',
     *          ' occurs within the',/,'minimum search increment for',
     *          ' variable: ',i1,', this often occurs as YCO2 => 1',
     *          ' or => 0, you may be able to correct this',/,
     *          'by reducing the default increment for this variable',
     *          ' (',g12.3,') in perplex_option.dat.',/,
     *          'Equilibria involving the following assemblage may',
     *          ' not be delineated:',/,7(1x,a8)) 
      end

      subroutine outdel
c-----------------------------------------------------------------------
c output deltas of dependent extensities for a reaction
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      character*8 cname*5,vname,xname,exten(2)*7,names

      integer i

      common/ csta2  /xname(k5),vname(l2)/ csta4 /cname(k5) 
     *      / cst8 /names(k1)

      integer iffr,isr
      double precision vuf,vus
      common/ cst201 /vuf(2),vus(h5),iffr,isr

      integer jfct,jmct,jprct
      common/ cst307 /jfct,jmct,jprct

      integer jds,ifr
      double precision du,dv
      common/ cst21 /du(2),dv(2),jds(h5),ifr

      integer iff,idss,ifug,ifyn,isyn
      common/ cst10  /iff(2),idss(h5),ifug,ifyn,isyn

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  

      integer ids,isct,icp1,isat,io2
      common/ cst40 /ids(h5,h6),isct(h5),icp1,isat,io2

      save exten 
c                                 this is a bullshit trick and
c                                 will cause errors of someone
c                                 uses a function other than G.
      data exten/'-V(j/b)','S(j/k)'/
c                                 composant stoichiometry:
 
      do i = 1, isat 
         write (n3,1000) cname(icp+i),vus(i),names(jds(i))
      end do 
c                                 fluid stoichiometry:
      if (ifyn.eq.0) then 
         do i = 1, 2
            if (iff(i).ne.0) write (n3,1010) names(i),vuf(i)
         end do 
      end if 
c                                 mobile components:
      do i = 1, jmct
         write (n3,1020) cname(jprct+i),du(i),vname(3+i)
      end do
c                                 dependent extensities 1 and 2, 
c                                 normally entropy and negative 
c                                 volume:
      do i = 1, 2
         write (n3,1030) exten(i), dv(i),vname(i)
      end do

1000  format (10x,'Delta(',a7,') =',g9.3,1x,
     *            '(saturated composant=',a8,')')
1010  format (10x,'Delta(',2x,a5,') =',g9.3,1x,
     *            '(saturated phase component)')
1020  format (10x,'Delta(',a7,') =',g9.3,1x,
     *            '(dependent conjugate of ',a8,')')
1030  format (10x,'Delta(',a7,') =',g9.3,1x,
     *            '(dependent conjugate of ',a8,')')
      end 

      subroutine delrxn
c-----------------------------------------------------------------------
c delrxn computes the change in the dependent extensities for a reaction
c these are returned in cst21, delrxn is always called after UNIVEQ, so
c the free energy at the present conditions is zero, if called 
c independently, then grxn must be called first.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,j,ivfl,invert

      double precision delt,dtol,utol,ptol,og,v,tr,pr,r,ps,gval

      common/ cst5   /v(l2),tr,pr,r,ps 
     *      / cst87  /delt(l2),dtol,utol,ptol/ cst102 /ivfl

      integer iffr,isr
      double precision vuf,vus
      common/ cst201 /vuf(2),vus(h5),iffr,isr

      integer jfct,jmct,jprct
      common/ cst307 /jfct,jmct,jprct

      double precision cp
      common/ cst12 /cp(k5,k1)

      integer ifct,idfl
      common/ cst208 /ifct,idfl

      integer jds,ifr
      double precision du,dv
      common/ cst21 /du(2),dv(2),jds(h5),ifr

      integer iff,idss,ifug,ifyn,isyn
      common/ cst10  /iff(2),idss(h5),ifug,ifyn,isyn

      integer idr,ivct
      double precision vnu
      common/ cst25 /vnu(k7),idr(k7),ivct

      integer irct,ird
      double precision vn
      common/ cst31 /vn(k2,k7),irct,ird

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  

      integer ids,isct,icp1,isat,io2
      common/ cst40 /ids(h5,h6),isct(h5),icp1,isat,io2

      integer ipot,jv,iv
      common/ cst24 /ipot,jv(l2),iv(l2)
c-----------------------------------------------------------------------
c                                 save dependent extensities

c                                 it is necessary to calculate
c                                 the initial G cause iptol may
c                                 be large relative to the 
c                                 change in G due to delt
      call grxn (og)

      do i = 1, 2
         v(i) = v(i) + delt(i)
         call grxn (gval)   
         dv(i) = (og - gval) / delt(i)
         v(i) = v(i) - delt(i) 
      end do 
c                                 get deltas on mobile comps
      if (jmct.gt.0) then
         do i = 1, jmct
            du(i) = 0d0
         end do

         do i = 1, ivct 
            do j = 1, jmct 
               du(j) = du(j) + vnu(i) * cp(jprct+j,idr(i))
            end do 
         end do
      end if 
c                                 do inverse projection to 
c                                 get stoichiometries of 
c                                 constrained components:
      if (isat.gt.0) then
         if (isat.gt.1) then  
            do j = 1, isat - 1
               do i = j + 1, isat
                  vus(j) = vus(j) - vus(i) * cp(icp+j,idss(i))
               end do 
            end do
         end if 
c                                 now do the effects of 
c                                 the inverse projection on
c                                 saturated phase components:
         if (ifct.ne.0) then 
            do i = 1, isat
               do j = 1, ifct
                  vuf(j) = vuf(j) - vus(i) * cp(jfct+j,idss(i))
               end do
            end do 
         end if 
c                                 now do the effects of 
c                                 the inverse projection on
c                                 mobile components:
         if (jmct.ne.0) then 
            do i = 1, isat 
               do j = 1, jmct
                  du(j) = du(j) - vus(i) * cp(jprct+j,idss(i))
               end do 
            end do 
         end if 
      end if 
c                                 see if coefficients should 
c                                 change sign:
      invert = 1
      if (ivfl.eq.1.and.dv(1).gt.0d0) then 
         invert = 0
      else if (ivfl.eq.2.and.dv(2).gt.0d0) then 
         invert = 0
      else if (ivfl.gt.2) then 
         if (du(ivfl-3).gt.0d0) then 
c                                 diagram must have a mu variable 
            invert = 0
         else if (du(ivfl-3).eq.0d0) then
            if (iv(2).eq.3) then
               if (vuf(2).lt.0d0) invert = 0
            else if (iv(1).eq.3) then
               if (vuf(2).lt.0d0) invert = 0
            else if (du(iv(1)-3).gt.0d0) then
               invert = 0
            end if
         end if
      end if

      if (invert.eq.1) then 
         do i = 1, ivct 
            vnu(i) = -vnu(i)
c                                 why different signs?
            vn(ird,i) = vnu(i)
         end do

         do i = 1, isat
            vus(i) = -vus(i)
         end do 

         do i = 1, ifct
            vuf(i) = -vuf(i)
         end do

         do i = 1, jmct 
            du(i) = -du(i)
         end do

         dv(1) = -dv(1)
         dv(2) = -dv(2)

      end if   
c                                 save id's of saturated comps
      do i = 1, isat
         jds(i) = idss(i)
      end do

      end     

      subroutine header (icopt)
c-----------------------------------------------------------------------
c header writes the graphics file header for schreinemakers diagrams
c to lun n4
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      character*8 names,vname,xname,fname*10

      integer iind,idep,icopt,i,ipres,itemp,jind

      double precision c0,c1,c2,c3,c4,c5

      character*162 title
      common/ csta8 /title(4)

      common/ cst8   /names(k1)
     *      / csta2  /xname(k5),vname(l2)
     *      / csta7 /fname(h9)
     *      / cst316 /c0,c1,c2,c3,c4,c5,iind,idep

      integer iff,idss,ifug,ifyn,isyn
      common/ cst10  /iff(2),idss(h5),ifug,ifyn,isyn

      integer ipot,jv,iv1,iv2,iv3,iv4,iv5
      common/ cst24 /ipot,jv(l2),iv1,iv2,iv3,iv4,iv5

      integer iso
      common/ cst79 /iso

      double precision vmax,vmin,dv
      common/ cst9  /vmax(l2),vmin(l2),dv(l2)

      integer ikp
      common/ cst61 /ikp(k1)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  

      double precision atwt
      common/ cst45 /atwt(k0)
c-----------------------------------------------------------------------
c                              value to be read as icopt
      write (n4,1020) icopt
c                              phase and solution count:
      write (n4,1020) iphct, iso
c                              component and volatile counters:
      if (ifyn.eq.0.or.isyn.eq.0) then 
         write (n4,1020) 1, icp
      else
         write (n4,1020) 0, icp
      end if 
c                              phase names:
      write (n4,1010) (names(i), i = 1, iphct)
c                              solution phase flags:
      write (n4,1060) (ikp(i), i = 1, iphct)  
c                              solution names:
      if (iso.ne.0) write (n4,1040) (fname(i), i = 1, iso)
c                               make and output title blurb
      call maktit

      write (n4,1000) title                                
c                               find out which variables are p and t:
      ipres = 0 
      itemp = 0 
      jind = 0 

      do i = 1, ipot
         if (jv(i).eq.1) then
            ipres = i
         else if (jv(i).eq.2) then
            itemp = i
         end if
      end do 

      if (idep.eq.1) then 
         jind = itemp
      else if (idep.eq.2) then
         jind = ipres
      end if 
         


      write (n4,1020) ipot, (jv(i), i = 1, ipot), ipres, itemp
      write (n4,1110) jind,idep,c0,c1,c2,c3,c4
      write (n4,1030) (vmax(jv(i)), vmin(jv(i)), i= 1, ipot)
      write (n4,1000) (vname(jv(i)), i= 1, ipot)

1000  format (a)
1010  format (10a8)
1020  format (20(i8,1x))
1030  format (6(d11.5,1x))
1040  format (8a10)
1060  format (26(i2,1x))
1110  format (2(i1,1x),5(g12.6,1x))

      end

      subroutine lchk (lphi,lchkl)
c-----------------------------------------------------------------------
c a function subprogram to determine if a phase, lphi, lies
c below a g-x plane defined by the vertices idv.
c lchk has a value of 0 if the plane is stable, and 1 if not.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,lphi,lchkl,j,ier

      double precision gproj, gphi

      double precision g
      common/ cst2 /g(k1)

      double precision cp
      common/ cst12 /cp(k5,k1)

      integer ipvt,idv,iophi,idphi,iiphi,iflg1
      double precision a,b
      common/ cst23 /a(k8,k8),b(k8),ipvt(k8),idv(k8),
     *               iophi,idphi,iiphi,iflg1

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  
c-----------------------------------------------------------------------
      call uproj

      do i = 1, icp
        b(i) = gproj(idv(i))
      end do 

      g(lphi) = gproj(lphi)
      lchkl = 0

      call subst (a,ipvt,icp,b,ier)

      gphi = 0d0
      do j = 1, icp
         gphi = gphi + cp(j,lphi) * b(j)
      end do

      if (gphi.lt.g(lphi)) goto 99

      lchkl = 1

99    end

      subroutine newass (b,idv,idphi,ivd,ivi,iste)
c-----------------------------------------------------------------------
c newass determines the (pseudo-) divariant assemblages generated by
c a reaction based on the stoichiometric coefficients of the phases
c in the reaction as determined by routine balanc.

c input:  idv(icp) - a vector containing the indices of phases in a
c           divariant assemblage.
c         idphi - the index of the phase which is involved with the
c           the divariant assemblage (idv) in a univariant equilibrium.
c         b(icp) - the stoichiometric coefficient of the phases inedxed
c           by idv in the reaction of these phases to form the phase
c           idphi. by convention the stoichiometric coeffeicient of
c           idphi is -1.0.

c output: jdv(icp) - a vector containing the indices of phases in a
c           divariant assemblage generated by the reaction between the
c           phases idphi and idv this vector is used by the routine
c           chkass. this assemblage is stable on the opposite 'side' of
c           the univariant equilibrium from the assemblage idv.

c referenced by: coface
c references to: chkass

c algorithm: phases with null or negative coefficients are always
c   compatible with eachother on the 'side' of a univariant reaction
c   which idv is metastable. the k phases with positive coefficients
c   define a chemographic simplex surrounding phases with null or
c   negative coefficients. consequently a univariant reaction will
c   generate k new divariant assemblages consisting of k-1 of the phases
c   in the original assemblage (idv) and all the phases with null or
c   negative coefficients. if only one phase has a positive coefficient
c   this phase is unstable and only one new assemblage will be generated
c   which consists of all the other phases (including idphi).
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer idv(k8),idpos(k5),jdv(k8),
     *        ipos,i,j,ineg,idphi,jtic,k,ivd,ivi,iste
     
      double precision b(k8)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  
c-----------------------------------------------------------------------
c                                 initialize counter:
      ipos = 0

      do i = 1, icp
c                                 is the b(i) coefficient zero/negative?
         if (b(i).ge.1d-5) then 
c                                 no, save in the array idpos and count.
            ipos = ipos + 1
            idpos(ipos) = idv(i)
         else 
c                                 yes, load into the array jdv
           jdv(i-ipos) = idv(i)
         end if 
      end do 
c                                 load idphi
      ineg = icp - ipos+1
      jdv(ineg) = idphi
c                                 generate the divariant assemblages
      do j = 1, ipos
c                                 make ipos permutations of ipos - 1 phases
         jtic = 0
         do k = 1, ipos
            if (k.ne.j) then
               jtic = jtic + 1
               jdv(ineg+jtic) = idpos(k)
            end if 
         end do            
c                                 send to the routine chkass to store th
c                                 the assemblage if unique.
         call chkass (jdv,ivd,ivi,iste)
      end do
c                                 done
      end

      subroutine onedim
c----------------------------------------------------------------------
c for phase diagrams with one extrinsic variable (iv1) onedim sorts
c equilibria, identified by newhld and located by coface, by veq(iv1).
c the sorted equilibria are then output to units n3 and n4 as requested.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ir,irct,ird,i,ist,ist1,itemp 

      double precision vn,vst

      common/ cst13 /ir(k2)
     *      / cst31 /vn(k2,k7),irct,ird   
 
      double precision vip
      common/ cst28 /vip(l2,k2)

      integer ipot,jv,iv1,iv2,iv3,iv4,iv5
      common/ cst24 /ipot,jv(l2),iv1,iv2,iv3,iv4,iv5
c-----------------------------------------------------------------------
c                             any rxns to be sorted?
      if (irct.le.1) goto 30
c                             initialize:
      ist = 2

      do i = 1, irct
         ir(i) = i
      end do 
c                             begin sort loop:
20    ist1 = ist - 1
      itemp = ir(ist1)
      vst = vip(iv1,itemp)

      do i = ist, irct
         if (vip(iv1,ir(i)).le.vst) then
            ir(ist1) = ir(i)
            ir(i) = itemp
            itemp = ir(ist1)
            vst = vip(iv1,itemp)
         end if 
      end do

      ist = ist + 1
      if (ist.gt.irct) goto 30
      goto 20
c                             end of sort loop
30    if (irct.eq.0) goto 99
c                             output invariant reactions and conditions:
      call outirn 

99    end
            
      subroutine outier
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      character vname*8, xname*8, text(kd2)*1

      integer icp2,i,j,iend

      common/ csta2 /xname(k5),vname(l2)/ cst81 /icp2

      integer ipid,ipct
      common/ cst29 /ipid(k2,k8),ipct
 
      double precision vip
      common/ cst28 /vip(l2,k2)

      integer io3,io4,io9
      common / cst41 /io3,io4,io9

      integer ipot,jv,iv
      common/ cst24 /ipot,jv(l2),iv(l2)

      integer isec,icopt,ifull,imsg,io3p
      common/ cst103 /isec,icopt,ifull,imsg,io3p

      integer ivarrx,ivarip,isudo,ivar
      common/ cst62 /ivarrx(k2),ivarip(k2),isudo,ivar
c-----------------------------------------------------------------------
c                                 write end of ptx data for univariant
c                                 equilibria marker to graphics file.
      if (io4.ne.1) then
         write (n4,1090)
         write (n4,1000) ipct,icp2
         do j = 1, ipct
c                                 n4:
            write (n4,1000) j, ivarip(j), (ipid(j,i),i = 1, icp2)
            write (n4,1060) (vip(jv(i),j), i = 1, ipot)
         end do     
      end if 
c                                 write ip names and locations to
c                                 units n6 and n3.
      if (ipct.gt.0.and.io3.ne.1) then 
         write (n3,1050)
         write (n3,1030)
         if (io3p.ne.1) then 
            do j = 1, ipct
         
               call iptext (text,iend,j)
               write (n3,1080) j, ivarip(j), (text(i), i= 1, iend)
               write (n3,1020)
               write (n3,1040) (vname(jv(i)), vip(jv(i),j), i= 1, ipot)
            end do 
         end if 
      end if
c                                 originally outier allowed output of multiple
c                                 sections, the lines below would have to be 
c                                 moved out to recover this functionality.
      if (io3.ne.1.and.icopt.eq.1) then

         write (n3,1160)
c                                 output cumulative equilibrium lists
         call outlst

      end if 

1000  format (14(i6,1x))
1020  format (15x,'occurs at:')
1030  format ('(pseudo-) invariant points are summarized below:')
1040  format (25x,a8,'=',g12.6)
1050  format (/,80('-'),/)
1060  format (5(1x,g12.6))
1080  format (/,' (',i6,'-',i1,') ',380a1)
1090  format (' 1 1 1 1 1 1 1 1 1 EOR',/,' 1.0 EOR')
1160  format (/,80('-'),/)

      end

      subroutine outirn
c-----------------------------------------------------------------------
c outirn writes the identity and equilibrium conditions of invariant
c equilibria to units n3 and n4 after they have been sorted by onedim.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      character*8 vname, xname*8, rxnstr*(kd2)

      integer i,j,k,l,irct,ird,ir,irnms

      double precision vn

      common/ cst31 /vn(k2,k7),irct,ird/ cst104 /rxnstr(k2)
     *      / cst13 /ir(k2)/ csta1 /irnms(k2,k7)
     *      / csta2 /xname(k5),vname(l2)
 
      double precision vip
      common/ cst28 /vip(l2,k2)

      integer ipot,jv,iv1,iv2,iv3,iv4,iv5
      common/ cst24 /ipot,jv(l2),iv1,iv2,iv3,iv4,iv5

      integer io3,io4,io9
      common / cst41 /io3,io4,io9

      integer irv
      common/ cst35 /irv(k2)

      integer ivarrx,ivarip,isudo,ivar
      common/ cst62 /ivarrx(k2),ivarip(k2),isudo,ivar
c-----------------------------------------------------------------------
c                             write invariant reactions:
      do i = 1, irct
         j = ir(i)
c                             the index on these two guys used to be j
c                             must be wrong? if so why did it work 
c                             before? changed to i. you know sometimes
c                             i'm so stupid, i just can't believe it.
         k = irv(j)
         l = ivarrx(j)
         if (io3.ne.1) then 
c                             output to print file:
            if (l.eq.1) then
               write (n3,1030) j, l, rxnstr(j)
            else
               write (n3,1040) j, l, rxnstr(j)
            end if 
            write (n3,1070) vname(iv1), vip(iv1,j), vname(iv2), 
     *                                  vip(iv2,j)
            write (n3,1160)
         end if 
c                             output to graphics file
         if (io4.ne.1) then 
            write (n4,1000) j, k, l, vip(iv1,j), (irnms(j,l), l= 1, k)
            write (n4,1050) (vn(j,l), l = 1, k)
         end if 
      end do 

1000  format (i6,1x,2(i6,1x),g12.5,1x,6(i6,1x))
1030  format ('The equilibrium of the invariant reaction:',//,
     *       ' (',i6,'-',i1,') ',a)
1040  format ('The equilibrium of the pseudoinvariant reaction:',//,
     *       ' (',i6,'-',i1,') ',a)
1050  format (12x,6(g9.2,1x))
1070  format (/,'occurs at ',a,'=',g12.6,' and ',a,'=',g12.6)
1160  format ('      ----------------------------------------')

      end

      subroutine outlst
c----------------------------------------------------------------------
c outlst writes cumulative lists of the univariant and invariant
c equilibria identified during a single execution of the program
c to the graphics and print files. 
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      character*8 names,rxnstr*(kd2)

      integer jbug,i,irct,ird,k

      double precision vn

      common/ cst104 /rxnstr(k2)
     *      / cst31 /vn(k2,k7),irct,ird
     *      / cst8  /names(k1)

      integer icp2
      common/ cst81 /icp2

      integer ipid,ipct
      common/ cst29 /ipid(k2,k8),ipct

      integer ibug
      common/ cst105 /ibug(k2)

      integer io3,io4,io9
      common / cst41 /io3,io4,io9

      integer ivarrx,ivarip,isudo,ivar
      common/ cst62 /ivarrx(k2),ivarip(k2),isudo,ivar
c------------------------------------------------------------------------
c                             printing to n4 is suppressed for the
c                             erlanger system 
      jbug = 0

      do i = 1, irct
         if (ibug(i).ne.0) jbug = 1
      end do    

      if (jbug.ne.0) then 
         write (*,1160)
         write (*,1190)
         do i = 1, irct
            if (ibug(i).ne.0) write (*,1000) i,ivarrx(i),rxnstr(i)
         end do 
         write (*,1160)
      end if 

      if (io3.eq.1) goto 99
c                             invariant point lists:
      if (ipct.ne.0) then 
c                             output print file heading:
c                             output ip lists:
         write (n3,1030)
         do i = 1, ipct
            write (n3,1150) i,ivarip(i),(names(ipid(i,k)),k=1,icp2)
         end do 
      end if 
c                              write end of invariant point list
c                              marker to graphics file.
      write (n3,1160)
c                              write list of stable univariant
c                              equilibria to units n3 and n4.
      write (n3,1070)

      do i = 1, irct
         write (n3,1000) i,ivarrx(i),rxnstr(i)
      end do    
c                              if itic= 0 no equilibria were identified
c                              goto 9000 to write no eq mess.
      if (irct.eq.0) goto 9000

      write (n3,1160)
c                              jbug ne 0, one or more equilibria
c                              may not have completely defined, write
c                              list of possible failures:
      if (jbug.ne.0) then 
         write (n3,1160)
         write (n3,1190)
         do i = 1, irct
            if (ibug(i).ne.0) write (n3,1000) i,ivarrx(i),rxnstr(i)
         end do 
         write (n3,1160)
      end if 

      goto 99
c                              write no equilibrium messages:
9000  write (n3,1180)

1000  format (' (',i6,'-',i1,') ',a)
1030  format ('(pseudo-) invariant points are summarized below:',/)
1070  format ('(pseudo-) univariant equilibria are summarized ',
     *       'below:',/)
1150  format (' (',i6,'-',i1,') ',12(a,1x))
1160  format (//,80('-'),//)
1180  format ('no univariant or invariant equilibria occur.')
1190  format ('WARNING!! The stability fields of the following',
     *        ' equilibria may',/,'have been entirely or',
     *        ' partially skipped in the calculation: ',/)
99    end

      subroutine iptext (text,iend,jp)
c-----------------------------------------------------------------------
c subprogram to write a text label for an invariant point 
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      character*8 names, text(kd2)*1, fname*10

      integer i,j,ist,id,iend,jp

      common/ cst8  /names(k1)
     *      / csta7 /fname(h9)

      integer icp2
      common/ cst81 /icp2

      integer ipid,ipct
      common/ cst29 /ipid(k2,k8),ipct

      integer ikp
      common/ cst61 /ikp(k1)
c----------------------------------------------------------------------
c                             dump the names into the array text
      ist = 1

      do 20 i = 1, icp2
         id = ipid(jp,i)
         if (ikp(id).eq.0) then
c                             simple compound:
            read (names(id),1000) (text(j), j = ist, ist + 7)
            ist = ist + 8
         else  
c                             solution phases:
            read (fname(ikp(id)),1000) (text(j), j = ist, ist + 9)
            text(ist + 10) = '('
            read (names(id),1000) (text(j), j = ist + 11, ist + 18)
            text(ist + 19) = ')'
            ist = ist + 20
         end if 
         text(ist) = ' '
20       ist = ist + 1 
c                             now filter out blanks
      iend = 1
      do 40 i = 2, ist - 1
         if (text(i).eq.' '.and.text(i + 1).eq.' ') then 
            goto 40 
         else if (text(i).eq.' '.and.text(i + 1).eq.')') then 
            goto 40
         else if (text(i).eq.' '.and.text(i + 1).eq.'(') then 
            goto 40 
         end if 
      iend = iend + 1
      text(iend) = text(i)
40    continue

1000  format (20a1)

      end 

      subroutine rxntxt (iend,jend,text,alpha)
c-----------------------------------------------------------------------
c subprogram to write a text label for a reaction 
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      character*8 names,text(kd2)*1,fname*10,alpha(130)*1,
     *            string*(kd2),rxnstr*(kd2)

      integer iplus(k5),iminus(k5),jplus(k5),jminus(k5),ip,im,i,j,ist,
     *        ione(k7),jone(k7),iend,jst,is,jend,id

      common/ csta7 /fname(h9)/ cst8 /names(k1)
     *      / cst104 /rxnstr(k2)

      integer irct,ird
      double precision vn
      common/ cst31 /vn(k2,k7),irct,ird

      integer idr,ivct
      double precision vnu
      common/ cst25 /vnu(k7),idr(k7),ivct

      integer ikp
      common/ cst61 /ikp(k1)
c----------------------------------------------------------------------

      ip = 0 
      im = 0
      jend = 1
c                             load id's of phases with + or - coeffs 
      do i = 1, ivct 
         if (vnu(i).gt.0d0) then 
            ip = ip + 1
            iplus(ip) = idr(i)
            jplus(ip) = i
         else 
            im = im + 1
            iminus(im) = idr(i)
            jminus(im) = i
         end if 
      end do 

      do i = 1, im 
         jone(i) = jminus(i)
         ione(i) = iminus(i)
      end do 

      do i = 1, ip
         j = i + im
         jone(j) = jplus(i)
         ione(j) = iplus(i)
      end do 
c                             dump coefficients into string:
      write (string,1010) (vnu(jone(i)),i = 1, ivct)
c                             read alpha from string:
      jst = 6 + ivct * 11
      read (string,1000) (alpha(i), i = 1, jst)
c                             filter blanks and add right )
      do 50 i = 1, jst 
         if (alpha(i).eq.' '.and.alpha(i + 1).eq.' ') then 
            goto 50
         else if (alpha(i).eq.' '.and.alpha(i + 1).eq.',') then
            goto 50
         end if 
         alpha(jend) = alpha(i)
         jend = jend + 1
50    continue

      if (alpha(jend-1).eq.' '.or.alpha(jend-1).eq.',') 
     *                                       jend = jend - 1  
      alpha(jend) = ')'
c                             now dump the names into the array text
      ist = 1
      is = 1
80    do i = is, im 
         id = ione(i)
         if (ikp(id).eq.0) then
c                             simple compound:
            read (names(id),1000) (text(j), j = ist, ist + 7)
            text(ist+8) = ' '
            ist = ist + 9
         else  
c                             solution phases:
            read (fname(ikp(id)),1000) (text(j), j = ist, ist + 9)
            text(ist + 10) = '('
            read (names(id),1000) (text(j), j = ist + 11, ist + 18)
            text(ist + 19) = ')'
            text(ist+20) = ' '
            ist = ist + 21
         end if 
      end do 

      if (is.eq.1) then 
         text(ist) = '='
         text(ist + 1) = ' '
         ist = ist + 2
         is = im + 1
         im = ivct 
         goto 80 
      else 
         text(ist) = ' ' 
      end if 
c                             now filter out blanks
      iend = 1
      do 40 i = 2, ist
         if (text(i).eq.' '.and.text(i + 1).eq.' ') then 
            goto 40 
         else if (text(i).eq.' '.and.text(i + 1).eq.')') then 
            goto 40
         else if (text(i).eq.' '.and.text(i + 1).eq.'(') then 
            goto 40 
         end if 
         iend = iend + 1
         text(iend) = text(i)
40    continue
c                             this is redundant if the 
c                             reaction was found before. 
      write (rxnstr(ird),1000) (text(i),i=1,iend)

1000  format (434a1)
1010  format ('Alpha(',13(g9.3,', '))

      end 

      subroutine outrxn (ip,ier)
c-----------------------------------------------------------------------
c outrxn writes the identity and v1-v2 loci of univariant equilibria
c to units n3 and n4.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      character text(kd2)*1, alpha(130)*1

      integer i,ier,jend,iend,ip

      integer ipt2
      double precision ptx
      common/ cst32 /ptx(l5),ipt2

      integer ibug
      common/ cst105 /ibug(k2)

      integer idr,ivct
      double precision vnu
      common/ cst25 /vnu(k7),idr(k7),ivct

      integer iflag
      common/ cst7 /iflag

      integer irct,ird
      double precision vn
      common/ cst31 /vn(k2,k7),irct,ird

      integer io3,io4,io9
      common / cst41 /io3,io4,io9

      integer jtest,jpot
      common/ debug /jtest,jpot

      integer isec,icopt,ifull,imsg,io3p
      common/ cst103 /isec,icopt,ifull,imsg,io3p

      integer ikp
      common/ cst61 /ikp(k1)

      integer ivarrx,ivarip,isudo,ivar
      common/ cst62 /ivarrx(k2),ivarip(k2),isudo,ivar
c----------------------------------------------------------------------

      do i = 1, ivct 
         if (ikp(idr(i)).lt.0) goto 99
      end do 
c                                set bug flag to indicate that
c                                the equilibrium was at least 
c                                partially located.
      if (ier.lt.2) ibug(irct) = 0
c                                write stefano's list to 
c                                reaction_list.dat
      if (ird.eq.irct.and.jtest.eq.3) call stetxt

      if (ifull.eq.0) then 
         call rxntxt (iend,jend,text,alpha)
      else  
         call fultxt (iend,text)
      end if 

      if (icopt.eq.3) goto 99
         
      if (imsg.eq.0) write (*,1110) ird,(text(i),i = 1, iend)

      if (io3p.eq.1) goto 10

      write (n3,1120) ird,ivarrx(ird),(text(i),i = 1, iend)

      if (ifull.ne.0) goto 80

      write (n3,1130) (alpha(i),i = 1,jend)

      if (ipt2.lt.3) then 
         write (n3,*)
         goto 99
      end if 

      call outdel 

80    write (n3,*)
      write (n3,1000) (ptx(i), i = 1, ipt2)
      write (n3,*)

      if (ier.ne.0) goto 10
      if (iflag.eq.1) goto 40 
      goto 10

40    write (n3,1070) ip        
      write (n3,1040)

10    if (io4.eq.1) goto 99

      write (n4,1080) ipt2,ird,ivar,ivct,(idr(i),i=1,ivct)
      write (n4,1090) (vnu(i),i = 1, ivct)
      write (n4,1000) (ptx(i),i = 1, ipt2)

1000  format (1x,g12.6,1x,g12.6,2x,g12.6,1x,g12.6,2x,g12.6,1x,g12.6)
1040  format (/)
1070  format ('the equilibrium extends to invariant point (',i6,')')
1080  format (16(i6,1x))
1090  format (6(g12.4,1x)) 
1110  format ('finished with equilibrium (',i6,') ',434a1)
1120  format (' (',i6,'-',i1,') ',434a1)
1130  format (/,10x,90a1)

99    end

      subroutine pchk (dg,igo)
c-----------------------------------------------------------------------
c a subprogram which returns the difference in free energy
c beteen a g-x plane, defined by the assemblage idv, and a phase lphi.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer igo,i,j,ier

      double precision delt,dtol,utol,ptol,gproj,dg

      common/ cst87 /delt(l2),dtol,utol,ptol

      double precision cp
      common/ cst12 /cp(k5,k1)

      integer ipvt,idv,iophi,idphi,iiphi,iflg1
      double precision a,b
      common/ cst23 /a(k8,k8),b(k8),ipvt(k8),idv(k8),
     *               iophi,idphi,iiphi,iflg1

      integer iflag
      common/ cst7 /iflag

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  
c-----------------------------------------------------------------------
c                                 compute energies:
      igo = 0

      call uproj

      do i = 1, icp
         b(i) = gproj(idv(i))
      end do 

      dg = gproj(idphi)
c                                 solve for chemical potentials:
      call subst (a,ipvt,icp,b,ier)
c                                 compute energy difference:
      do j = 1, icp
         dg = dg - cp(j,idphi)*b(j)
      end do 

      if (dabs(dg).lt.ptol) then
        igo = 1
        iflag = 1

        call ssaptx 
      else
        if (dg.gt.dtol) then
          iflag= 0

          call ssaptx
        else
          iflag= 1
        end if
      end if

      end

      subroutine pivots (ier)
c--------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer j,k,ier

      double precision cp
      common/ cst12 /cp(k5,k1)

      integer ipvt,idv,iophi,idphi,iiphi,iflg1
      double precision a,b
      common/ cst23 /a(k8,k8),b(k8),ipvt(k8),idv(k8),
     *               iophi,idphi,iiphi,iflg1

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  
c--------------------------------------------------------------------
c                              calculate pivots for the matrix 'a'
      do k = 1, icp
         do j = 1, icp
            a(k,j) = cp(j,idv(k))
         end do 
      end do  

      call factor (a,icp,ipvt,ier)

      end

      subroutine schk (lphi)
c-----------------------------------------------------------------------
c schk my ___
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,ipvt,iophi,lphi,j,idphi,iiphi,iflg1,idv,ier

      double precision u,a,dtol,dg,utol,delt,ptol

      common/ cst23 /a(k8,k8),u(k8),ipvt(k8),idv(k8),
     *               iophi,idphi,iiphi,iflg1
     *      / cst87 /delt(l2),dtol,utol,ptol

      double precision g
      common/ cst2 /g(k1)

      double precision cp
      common/ cst12 /cp(k5,k1)

      integer iflag
      common/ cst7 /iflag

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  
c-----------------------------------------------------------------------
c                                 determine chemical potentials
      iflag = 0
      do i = 1, icp
         u(i) = g(idv(i))
      end do 

      call subst (a,ipvt,icp,u,ier)
c                                 test phases, not=iophi, against the
c                                 assemblage idv(i),i = 1, icp
      do 20 i = istct, iphct
         if (i.eq.iophi.or.i.eq.lphi) goto 20
         dg = g(i)

         do j = 1, icp
            dg = dg - cp(j,i) * u(j)
         end do 

         if (dg.gt.dtol) goto 20
c                                check that a phase is not metastable
c                                with respect to itself, this
c                                could be remedied by changing the 
c                                value of dtol.
         do j = 1, icp
            if (idv(j).eq.i) goto 20
         end do 

         iflag = iflag + 1
         idphi = i
         goto (20,40,40), iflag

20    continue

40    end

      subroutine nschk 
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision a,u,delt,dtol,utol,ptol

      integer ipvt,idv,iophi,idphi,iiphi,iflg1,i,j,ier

      double precision blim, ulim, dg
      common/ cxt62 /blim(l2),ulim(l2),dg

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  

      common/ cst23 /a(k8,k8),u(k8),ipvt(k8),idv(k8),
     *               iophi,idphi,iiphi,iflg1
     *      / cst87 /delt(l2),dtol,utol,ptol

      double precision g
      common/ cst2 /g(k1)

      double precision cp
      common/ cst12 /cp(k5,k1)

      integer iflag
      common/ cst7 /iflag
c-----------------------------------------------------------------------
c                                 determine chemical potentials
      iflag = 0

      do i = 1, icp
         u(i) = g(idv(i))
      end do 

      call subst (a,ipvt,icp,u,ier)
c                                 test phases, not=iophi, against the
c                                 assemblage idv(i),i = 1, icp
      do i = istct, iphct

         if (i.eq.iophi) cycle

         dg = g(i)

         do j = 1, icp
            dg = dg - cp(j,i) * u(j)
         end do 

         if (dg.gt.dtol) cycle
c                                check that a phase is not metastable
c                                with respect to itself, this
c                                could be remedied by changing the 
c                                value of dtol.
         do j = 1, icp
            if (idv(j).eq.i) goto 90
         end do 

         iflag = iflag + 1
         idphi = i

         if (iflag.eq.1) then 
            cycle
         else 
            goto 99
         end if 
 
90    end do

99    end

      subroutine search (vst,vsti,ist,iste,ivi,ivd,div,iedge,ier)
c-----------------------------------------------------------------------
c search looks for an equilibrium around the periphery of the coordinate
c frame specified for a spd calculation. if a stable equilibrium is
c found ier = 0, if an error occurs ier =2, and if no equilibrium is
c found ier =1.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      character*8 names

      integer iedge,ier,iugh,i,ist,iste,ivi,ivd,jer,j,k

      double precision v,x,vst,ddv,div,vlast,vsti

      common/ cst5 /v(l2),x(4)/ cst8 /names(k1)

      integer ipvt,idv,iophi,idphi,iiphi,iflg1
      double precision a,b
      common/ cst23 /a(k8,k8),b(k8),ipvt(k8),idv(k8),
     *               iophi,idphi,iiphi,iflg1

      integer iflag
      common/ cst7 /iflag

      integer ipot,jv,iv1,iv2,iv3,iv4,iv5
      common/ cst24 /ipot,jv(l2),iv1,iv2,iv3,iv4,iv5

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  

      double precision vmax,vmin,dv
      common/ cst9  /vmax(l2),vmin(l2),dv(l2)
c-----------------------------------------------------------------------
c                                 initialization
      iflg1 = 0
      ier = 0
      iugh = 0

      v(iv1) = vst
      v(iv2) = vst

      do i = ist, iedge

         iste = i
c                                 set default dependent and independent
c                                 variables and increments for each edge
c                                 to be searched.
         if (i.eq.1) then 
c                                 traverse 1.
            ivi = iv2
            ivd = iv1
            ddv = dv(ivd)
            div = dv(ivi)
            v(ivi) = vmin(ivi)

         else if (i.eq.2) then 
c                                 traverse 2.
            ivi = iv1
            ivd = iv2
            ddv = dv(ivd)
            div = -dv(ivi)
            v(ivi) = vmax(ivi)

         else if (i.eq.3) then 
c                                 traverse 3.
            ivi = iv2
            ivd = iv1
            ddv = -dv(ivd)
            div = -dv(ivi)
            v(ivi) = vmax(ivi)

         else 
c                                 traverse 4.
            ivi = iv1
            ivd = iv2
            ddv = -dv(ivd)
            div = dv(ivi)
            v(ivi) = vmin(ivi)

         end if 

         call incdep (ivi) 

c                                 in some cases newass may generate
c                                 a metastable assemblage call flipit
c                                 to check. if the assemblage is
c                                 stable flipit will return jer = 0,
c                                 if it is metastable flipit tries a
c                                 reverse traverse and returns jer =1
c                                 if it is succesful, and jer =2 if not.
         if (i.eq.ist) then
            v(ivd) = vst
            v(ivi) = vsti
            call incdp0 
         end if 

         call flipit (ddv,vst,ivd,ist,i,jer)

         if (jer.eq.2) then 
            write (*,1000) (names(idv(k)), k= 1, icp)
            write (*,1010) v
            ier = 1
            goto 99
         end if 

         vlast = v(ivd)

         if (jer.eq.1) goto 99
c                                 begin search:
100      v(ivd) = v(ivd) + ddv
         call incdep (ivd) 
c                                 out of range?:
         if (i.lt.3) then 

            if (v(ivd).gt.vmax(ivd)) then
               v(ivd) = vmax(ivd)
               ddv = dabs(vmax(ivd) - vlast)
            else if (i.eq.ist) then
               if (v(ivd).gt.vst) goto 130
               ddv = dabs(ddv)/2d0
               v(ivd) = vst
               iflg1= 0
               goto 100
            end if

         else 

            if (v(ivd).lt.vmin(ivd)) then
               v(ivd) = vmin(ivd)
               ddv = -dabs(vlast - vmin(ivd))
            else if (i.eq.ist) then
               if (v(ivd).lt.vst) goto 130
               ddv = -dabs(ddv)/2d0
               v(ivd) = vst
               iflg1 = 0
               goto 100
            end if

         end if 

         call incdep (ivd)
c                                 calculate phase energies:
130      call gall 
c                                 test stability of the assemblage:
         call asschk

         if (iflag.eq.0) vlast = v(ivd)
c                                 iflag= 0 still stable
c                                 iflag=1 metastable with respect to one
c                                 phase indexed by idphi.
c                                 iflag=2 metastable with respect to
c                                 multiple phases refine the
c                                 search increment.
         if (iflag.eq.1) goto 99
c                                 check if search is in range:
         if (i.lt.3) then 
            if (v(ivd).ge.vmax(ivd).and.iflag.eq.0) cycle
         else 
            if (v(ivd).le.vmin(ivd).and.iflag.eq.0) cycle
         end if 
c                                 refine increment if necessary:
         call delvar (ddv,iflag,iflg1)
c                                 check if increment has been refined
c                                 beyond the tolerance which is
c                                 arbitrarily set at 1.d-8.
         if (dabs(ddv).lt.1d-8) then

            iugh = iugh + 1

            if (i.lt.3) then 
c                                 traverses 1,2.
               ddv = dv(ivd)
            else  
c                                 traverses 3,4.
               ddv = -dv(ivd)
            end if

            v(ivd) = vlast
            call incdep (ivd) 
            ddv = ddv/1d1
            iflg1 = 0
   
            if (iugh.gt.3) then
               write (n3,1020) ivd,dv(ivd),(names(idv(j)),j = 1, icp)
               write (n3,*) ' '
               ier = 2
               goto 99
            end if 
      
            goto 100 

         end if 
c                                 increment conditions again:
         goto 100
c                                 next traverse.
      end do 

      ier = 1

1000  format ('Metastable assemblage in FLIPIT, ',
     *        'the assemblage is:',/,4x,12(1x,a))
1010  format ('v =',5(g12.6,1x))
1020  format (/,'**warning ver047** SEARCH: > 1 equilibrium',
     *         ' occurs within the minimum search',/,'increment for',
     *         ' variable: ',i1,'. This often occurs as YCO2 =>',
     *         ' 1 or => 0.',/,'You may be able to correct this by',
     *         ' assigning slightly different limits',/,'for YCO2 or',
     *         ' by reducing the default increment for this variable',
     *         ' (',g12.3,')',/,'in the perplex_option.dat.',/,
     *         'Equilibria involving the following assemblage may',
     *         ' not be delineated:',/,7(1x,a8)) 
99    end

      subroutine sfol1 (ivd,ivi,ier,dv,ikwk,irend)
c----------------------------------------------------------------------
c eliminated computed goto's introduced do loop, 1/9/09.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ikwk,irend,ivd,ivi,iswtch,jswtch,ier,jer,ip

      integer irchk
      common/ cst801 /irchk(k2)

      double precision delt,dtol,utol,ptol,dv
      common/ cst87 /delt(l2),dtol,utol,ptol  

      integer iflag
      common/ cst7 /iflag   

      integer irct,ird
      double precision vn
      common/ cst31 /vn(k2,k7),irct,ird

      double precision v,tr,pr,r,ps
      common/ cst5 /v(l2),tr,pr,r,ps

      double precision vmax,vmin,ddv
      common/ cst9 /vmax(l2),vmin(l2),ddv(l2)
c----------------------------------------------------------------------

      if (ikwk.eq.1) then
         irchk(ird) = 1
         goto 900
      end if

      iswtch = 0
      jswtch = 0 
      iflag = 0 

      do
c                                 begin traverse:
          v(ivi) = v(ivi) + dv
c                                 is search in range?
         if (v(ivi).gt.vmax(ivi)) then
            v(ivi) = v(ivi) - dv
            dv = (vmax(ivi) - v(ivi))
            v(ivi) = vmax(ivi)
         else if (v(ivi).lt.vmin(ivi)) then
            v(ivi) = v(ivi) - dv
            dv = (vmin(ivi) - v(ivi))
            v(ivi) = vmin(ivi)
         end if

         call incdep (ivi) 
c                                 solve for the equilibrium conditions:
         call univeq (ivd,ier)
c                                 on error return:
c                                 calling routine will switch variables.
         if (ier.ne.0) goto 999

         if (jswtch.eq.0.and.
     *      (v(ivd).gt.vmax(ivd).or.v(ivd).lt.vmin(ivd))) then 
            call switch (dv,ivi,ivd,jer)
            jswtch = 1

            if (jer.eq.0) cycle

         end if 
c                                 test if the equilibrium is stable:
         call gall 
         call nschk 

         if (iflag.eq.0) then 
c                                 iflag= 0 stable:
c                                 dependent v in range? if > vmax
c                                 or < vmin reset conditions, and
c                                 switch independent/dependent variables
            if (v(ivd).gt.vmax(ivd)) then

               v(ivd) = vmax(ivd)
               call incdep (ivd) 

            else if (v(ivd).lt.vmin(ivd)) then

               v(ivd) = vmin(ivd)
               call incdep (ivd)  

            else

               call assptx
c                                 if independent variable has reached a
c                                 limit of the diagram output traverse
c                                 otherwise increment the variable and
c                                 continue the traverse
               if (v(ivi).eq.vmax(ivi).or.v(ivi).eq.vmin(ivi)) then
c                                 flag irchk indicates that the
c                                 traverse ended on an edge.
                  irchk(ird) = 1

                  call maxend

                  exit

               end if
c                                 try for next point
               cycle

            end if
c                                 solve for the equilibrium with
c                                 the switched variables:
            call univeq (ivi,ier)
c                                 flag irchk indicates that the
c                                 traverse ended on an edge.
            irchk(ird) = 1
c                                 on error output traverse anyway
            if (ier.ne.0) exit 
c                                 assign final value
            call assptx 

         else if (iflag.eq.1) then
c                                 iflag=1 metastable to 1 phase, i.e, 
c                                 an invariant equilibrium:

c                                 iterate independent variable to refine
c                                 the conditions of the invariant point:
            call findjp (ivi,ivd,dv,ip)


         else if (iflag.eq.2) then   
c                                 iflag=2 metastable to > 1 phase:
c                                 reset to last stable condition
c                                 and redefine increment dv
            call reptx

            dv = dv/2d0
            if (dabs(dv).gt.delt(ivi)) cycle
c                                 refined increment down to dtol, try
c                                 switching iv's, set iswitch flag
            if (iswtch.eq.0) then 
               call switch (dv,ivi,ivd,jer)
               iswtch = 1
               if (jer.eq.0) cycle
            end if 

            call warn (73,delt(ivi),ivi,'SFOL1 ')

            iflag = 0 

         end if 

         exit 

      end do                       
c                                 output the traversed equilibrium:
900   ier = 0

      call outrxn (ip,ier)    
c                                 save conditions of the endpoint
      call svrend (ird,irend,jer)

999   end                         

      subroutine sfol2 (iovi,iovd,ip,irend)
c----------------------------------------------------------------------
c sfol2 generates and computes the loci of stable equilibria
c emanating from an ip indexed by the array ipid and the
c arguement (i).  since one of these equilibria is calculated
c by sfol1 or by previous calls to sfol2 only c+1 (at most)
c curves are calibrated.
c
c the phases identified by idv define a plane in g-x space, idphi
c is the phase which is tangent to this plane, and lphi is the
c ,or one of the, phases which is present at the ip but absent
c along the univariant point.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      character text(kd2)*1, alpha(130)*1

      logical bad

      integer ikwk,irend,ivd,ivi,ier,jer,ip,jchk,icp3,
     *        iovi,iovd,iend,i,j,itic,jtic,lphi,idno,icter,jend

      double precision div

      integer ipid,ipct
      common/ cst29 /ipid(k2,k8),ipct
 
      integer icp2
      common/ cst81 /icp2

      double precision vip
      common/ cst28 /vip(l2,k2)

      integer ipvt,idv,iophi,idphi,iiphi,iflg1
      double precision a,b
      common/ cst23 /a(k8,k8),b(k8),ipvt(k8),idv(k8),
     *               iophi,idphi,iiphi,iflg1

      integer ipt2
      double precision ptx
      common/ cst32 /ptx(l5),ipt2

      integer ibug
      common/ cst105 /ibug(k2)

      integer irct,ird
      double precision vn
      common/ cst31 /vn(k2,k7),irct,ird

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      integer io3,io4,io9
      common / cst41 /io3,io4,io9

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  

      integer isec,icopt,ifull,imsg,io3p
      common/ cst103 /isec,icopt,ifull,imsg,io3p

      integer ids,isct,icp1,isat,io2
      common/ cst40 /ids(h5,h6),isct(h5),icp1,isat,io2
c-----------------------------------------------------------------------
      jchk = 0 
      icp3 = icp + 3
      ivi = iovi
      ivd = iovd

      goto (150), io3
      goto (150), io3p

      call iptext (text,iend,ip)

      write (n3,1040)
      write (n3,1020) ip, (text(i),i = 1, iend)
      write (n3,1050)

150   do j = 2, icp2
         jtic = 0
         idno  = icp3 - j
         lphi  = ipid(ip,idno)
         idphi = ipid(ip,icp2)

         do i = 1, icp1
            if (i.ne.idno) then 
               jtic = jtic + 1
               idv(jtic) = ipid(ip,i)
            end if 
         end do 
c                                 balance the equilibrium, ier = 1
c                                 if the plane idv is degenerate
         itic= 0
c                                 set invariant point conditions
         v(ivd) = vip(ivd,ip)
         v(ivi) = vip(ivi,ip)
         call incdp0   

110      call balanc (b,idv,idphi,ier)
         if (ier.eq.1) goto 80
c                                 call pivots for lchk.
         call pivots (ier)
         if (ier.eq.0) goto 90
c                                 in some cases singularity of the
c                                 transpose of the concentration matrix
c                                 may not be detected by balanc, but,
c                                 will be detected for the concentration
c                                 matrix by pivots, in which case the
c                                 columns are permuted and the stoichiometry
c                                 of the reaction is redetermined.
80       if (itic.eq.icp-1) cycle
         icter =idphi
c                                 if a degenerate facet has been
c                                 generated by permutation permute
c                                 idphi with the vertices idv until
c                                 balanc succeeds.  for a c component
c                                 system no more than c-1 permutations
c                                 will be necessary.
         idphi = idv(icp-itic)
         idv(icp-itic) = icter
         itic = itic + 1
         goto 110

90       call assir (bad)
c                                 skip curves of variance > isudo
         if (bad) cycle
c                                 TESTING 
c                                 test endpoint.
         call svrend (ird,irend,ier)
         if (ier.eq.1) cycle
c                                 calculate the loci of the equilibia,
         ipt2 = 0
         call delrxn 
         call assptx
c                                 determine direction of traverse and
c                                 dependent and independent variables
         call wway (div,lphi,ivi,ivd,ikwk,ier)
c                                 begin traverse
         if (ier.ne.0) then 
c                                 call rxntext to save the equilibrium
c                                 in the list in case it is not found 
c                                 again:
            call rxntxt (iend,jend,text,alpha)
c                                 set the bug flag:
            ibug(irct) = 1
            cycle 

         end if 

         icter = 0

60       call sfol1 (ivd,ivi,ier,div,ikwk,irend)
         jchk = jchk + 1

         if (ier.eq.0) then 
            ivd = iovd
            ivi = iovi
            cycle
         end if 

         call switch (div,ivi,ivd,jer)  
  
         if (jer.ne.1) then 
            icter = icter + 1
            if (icter.lt.4) goto 60
         end if 
      
         call warn (20,v(ivi),ivi,'SFOL2 ')

         call outrxn (ip,1)
         ivi =iovi
         ivd=iovd
      end do 

      if (jchk.eq.0) call warn (74,v(ivi),ivi,'SFOL2 ')
      if (io3p.eq.0) write (n3,1040)

1020  format ('equilibria about invariant point (',i6,'):',//,
     *        3x,200a1)
1040  format ('                       ------')
1050  format (/,'are listed below:',/)

      end

      subroutine sreset (jok,iste,vt,v,vti,vi)
c----------------------------------------------------------------------
      implicit none

      double precision vt,v,vti,vi

      integer jok, iste
c                                 if the assemblages conditions
c                                 of stability are already known (jok=5)
c                                 or are known to extend beyond the
c                                 current condition (jok>iste) quit:
      if (jok.gt.iste) goto 99
c                                 if the known conditions of stability
c                                 have been extened, i.e., either:
c                                 jok=iste  or v>vt
      if (jok.ne.iste) then
         jok=iste
      else if (iste.lt.3) then
         if (vt.gt.v) goto 99
      else
         if (vt.lt.v) goto 99
      end if

      vt = v
      vti = vi

99    end  

      subroutine maxend 
c----------------------------------------------------------------------
c save the most advanced condition for the stability of an 
c assemblage along the periphery of a diagram. 
c ismax - is the traverse
c value - is the dependent variable
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      integer ipot,jv,iv1,iv2,iv3,iv4,iv5
      common/ cst24 /ipot,jv(l2),iv1,iv2,iv3,iv4,iv5

      integer ismax,igot
      double precision value
      common/ cst49 /value,ismax,igot

      double precision vmax,vmin,ddv
      common/ cst9  /vmax(l2),vmin(l2),ddv(l2)
c-----------------------------------------------------------------------
      if (igot.eq.0) igot = 1

      if (v(iv2).eq.vmin(iv2).and.ismax.eq.1) then
c                                 traverse 1
c                                 now don't allow backing up:
         if (value.lt.v(iv1)) value = v(iv1)

      else if (v(iv1).eq.vmax(iv1).and.ismax.lt.3) then
c                                 traverse 2
         if (ismax.eq.1) then

            ismax = 2
            value = v(iv2)

         else if (value.lt.v(iv2)) then

            value = v(iv2)

         end if

      else if (v(iv2).eq.vmax(iv2).and.ismax.lt.4) then
c                                 traverse 3
         if (ismax.lt.3) then

            ismax = 3
            value = v(iv1)

         else if (value.gt.v(iv1)) then

            value = v(iv1)

         end if 

      else if (v(iv1).eq.vmin(iv1)) then
c                                 traverse 2
         if (ismax.lt.4) then
            ismax = 4
            value = v(iv2)
         else if (value.gt.v(iv1)) then
            ismax = 4
            value = v(iv2)
         end if 

      end if 

      end 

      subroutine svrend (ird,irend,ier)
c----------------------------------------------------------------------
c save end points of univariant equilibria and 
c check if the intersection has already been recorded (ier =1).

c two cases:

c 1) reaction occurs on an edge, if the intersection has already 
c    occurred then ignore future occurrences.
c 2) reaction occurs at an IP only ignore occurences being traced in 
c    the same direction if the reaction is degenerate.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ird, irend, ier, irends(k2), j, jct

      double precision rends(2,16,k2), dx, dy, x, y, xx, yy

      integer ids,isct,icp1,isat,io2
      common/ cst40 /ids(h5,h6),isct(h5),icp1,isat,io2

      integer irv
      common/ cst35 /irv(k2)

      integer ipot,jv,iv
      common/ cst24 /ipot,jv(l2),iv(l2)

      integer isec,icopt,ifull,imsg,io3p
      common/ cst103 /isec,icopt,ifull,imsg,io3p

      double precision vmax,vmin,ddv
      common/ cst9  /vmax(l2),vmin(l2),ddv(l2)

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      save rends, irends
c-----------------------------------------------------------------------
      ier = 0 
c                                 if isec = 1, don't do any checks
      if (isec.eq.1) goto 99

      dx =  ddv(iv(1))
      dy =  ddv(iv(2))

      x = v(iv(1))
      y = v(iv(2))

      if (x.ne.vmin(iv(1)).and.x.ne.vmax(iv(1)) .and.
     *    y.ne.vmin(iv(2)).and.y.ne.vmax(iv(2))) then
c                                 the endpoint is not on an edge:
          if (isec.eq.2) then 
c                                 if isec = 2, don't check unless
c                                 curve intersects at an edge.
             goto 99
          else if (isec.eq.3) then 
c                                 if isec = 3, don't do a check for 
c                                 degenerate interior endpoints:
             if (irv(ird).lt.icp1) goto 99
          end if 
      else 
c                                 the endpoint is on an edge:
c                                 if isec = 4, don't do a check for
c                                 degenerate equilibria.
          if (isec.eq.4.and.irv(ird).lt.icp1) goto 99
      end if 
c                                 isec = 5, do check on all ends
      if (ird.le.irend) then 
c                                 reaction has been found already
c                                 look to see if its at the same 
c                                 endpoint (tolerance +/- dv).
         jct = irends(ird)
         do j = 1, jct 
            xx = rends(1,j,ird)
            yy = rends(2,j,ird)
            if (xx.gt.x-dx.and.xx.lt.x+dx.and.
     *          yy.gt.y-dy.and.yy.lt.y+dy) then 
c                                 condition match has occurred
               ier = 1
               goto 99
            end if 
         end do 
c                                 no condition match has occurred
         jct = jct + 1
         if (jct.gt.16) then 
            call warn (999,x,ier,'SVREND')
            jct = 16 
         end if  

      else 
c                                 reaction has not been recorded
         jct = 1
         irend = irend + 1

      end if 

      if (ird.gt.k2) call error (206,x,k2,'SVREND')
      irends(ird) = jct 
      rends(1,jct,ird) = x
      rends(2,jct,ird) = y

99    end 

      subroutine wway (div,lphi,ivi,ivd,ikwk,ier)
c-----------------------------------------------------------------------
c wway determines from which side of an invariant point a specific
c univariant equilibrium emanates. the appropriate sign for the
c increment of the independent variable (dv) is also determined.
c where necessary wway also refines dv and/or may switch the
c independent and independent intensive variables (ivi,ivd). vo1
c and vo2 are the values of the intensive variables iv1 and iv2 at
c the invariant point. lphi is the id of the absent phase.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      character*8 names

      double precision vo(l2),vvi(2),vvd(2),ddv(2),
     *                 odiv,del,div,order

      integer l,ikwk,isign,lchkl,i,iv,jfail,ier,ivi,
     *        ivd,icter,lphi,jswit,ifail

      common/ cst8 /names(k1)

      integer ipvt,idv,iophi,idphi,iiphi,iflg1
      double precision a,b
      common/ cst23 /a(k8,k8),b(k8),ipvt(k8),idv(k8),
     *               iophi,idphi,iiphi,iflg1

      integer idr,ivct
      double precision vnu
      common/ cst25 /vnu(k7),idr(k7),ivct

      integer iflag
      common/ cst7 /iflag

      integer irct,ird
      double precision vn
      common/ cst31 /vn(k2,k7),irct,ird

      integer ipot,jv,iv1,iv2,iv3,iv4,iv5
      common/ cst24 /ipot,jv(l2),iv1,iv2,iv3,iv4,iv5

      double precision vmax,vmin,dv
      common/ cst9  /vmax(l2),vmin(l2),dv(l2)

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps
c-----------------------------------------------------------------------
      iophi =idphi
      ikwk = 0
      icter = 0
      ier = 0
      vo(iv1) = v(iv1)
      vo(iv2) = v(iv2)
c                                 first try iv1 as the iv
      ivi = iv2
      ivd = iv1

      ifail= 0

5     call univeq (ivd,ier)

      if (ier.eq.0) then
c                                 is the equilibrium condition
c                                 for this reaction consistent
c                                 with the invariant point?
         if ( dabs( (v(ivd)-vo(ivd))/dv(ivd) ).gt.0.1d0) then 
            vo(iv1) = v(iv1)
            vo(iv2) = v(iv2)
            call warn (18,v(ivd),ird,'WWAY')
         end if 

      else
         goto (9800),ifail
         ivi = iv1
         ivd = iv2
c                                 added april 2001 without check!!
         div = dv(ivi)

         v(iv1) = vo(iv1)
         v(iv2) = vo(iv2)
         call incdp0

         ifail = ifail + 1
         goto 5
      end if
c                                 test that the ip isn't on an edge
      if (v(iv1).eq.vmin(iv1).or.v(iv1).eq.vmax(iv1).or.
     *    v(iv2).eq.vmin(iv2).or.v(iv2).eq.vmax(iv2)) goto 9600

      jswit = 0

190   order = 1d0
c                                 set iv increment
60    div = dv(ivi) / order

      ifail = 0

      do i = 1, 2
c                                 because of numerical slop 
c                                 an equilibrium may appear stable
c                                 on its metastable extension close to
c                                 an ip, because of this wway always
c                                 tests for stability on both extensions
c                                 if wway cannot solve for equilibrium
c                                 conditions on both sides of an ip
c                                 the increment is decreased and if this
c                                 fails the independent and dependent
c                                 variables are switched.
         jfail = 0

15       v(ivi) = vo(ivi) + div
c                                 check if the increment is too large
         if (v(ivi).gt.vmax(ivi)) then
            v(ivi) = vmax(ivi)
         else if (v(ivi).lt.vmin(ivi)) then 
            v(ivi) = vmin(ivi)
         end if 
c                                 set the dependent variable to ip:
         v(ivd) = vo(ivd)

         call incdp0 
c                                 solve for the equilibrium
         call univeq (ivd,ier)

         if (jfail.le.5.and.ier.eq.0.and.(v(ivd).gt.vmax(ivd)
     *       .or.v(ivd).lt.vmin(ivd))) then
c                                 increment to big:
            div = div / 1d1
            jfail = jfail + 1
            goto 15
  
         else if (jswit.eq.0.and.ier.gt.0) then 
c                                 on error switch dependent and
c                                 independent variables
            iv  = ivi
            ivi = ivd
            ivd = iv
            jswit = 1
            jfail = 0

            goto 190

         else if (jswit.gt.0.and.ier.gt.0) then

            call warn (999,dv(1),ivi,'WWAY')
            ier = 1
            return

         end if 

         vvi(i) = v(ivi)
         vvd(i) = v(ivd)
         ddv(i) = div
c                                 determine if the equilibrium is
c                                 metastable with respect to lphi
c                                 (lchkl=1):
         call lchk (lphi,lchkl)
         if (lchkl.eq.0) then
            ifail = ifail + 1
            isign = i
         end if
c                                 end of looop
         div = -div
      end do 
c                                 at this stage there are three
c                                 possibilities indicated by ifail
      if (ifail.eq.1) goto 50
c                                 ifail= 0 or 2 the equilibrium was stable
c                                 or metastable on both sides, the most 
c                                 plausible scenario is the increment is
c                                 too big and overshoots the stability,
c                                 in this scenario ifail = 2; however the
c                                 ip could be located inside of, rather than
c                                 on the edge of the stability field, in 
c                                 which case there is a small chance the
c                                 increment is too small to cross the border.
c                                 
c                                 so first decrease the increment:
         if (order.ge.1d0) then 
            order = order * 10d0 
c                                 stop decreasing after four orders 
c                                 of magnitude, and try increasing
c                                 increment:
            if (order.gt.1d4) order = 1d-1
         else 
            order = order / 10d0
c                                 stop increasing after two orders 
c                                 of magnitude, cause there should 
c                                 never be a default increment much less 
c                                 than 1 % of the variable range.
            if (order.lt.1d-2) goto 9500
         end if 
c                                 why the **** was this here?
c         ivi = iv1
c         ivd = iv2

      goto 60
c                                 ifail=1, stable extension identified.
c                                 determine the sign of the increment
c                                 and its magnitude:
50    if (isign.eq.1) then
         div = dv(ivi)
         del= vmax(ivi)-vo(ivi)
         if (del.lt.div) div = del/3.d0
      else
         div = -dv(ivi)
         del= vo(ivi)-vmin(ivi)
         if (del.lt.-div) div = -del/3.d0
      end if
c                                 reset variable values and save
      v(ivi) = vvi(isign)
      v(ivd) = vvd(isign)
      call incdp0
  
      odiv = ddv(isign)
      icter = 0
c                                 now check the stability with respect
c                                 to all phases
      order = 1d0
70    call gall 
      call schk (lphi)
c                                 stable if iflag= 0 from schk
c                                 now check if dependent variable is in
c                                 range.
      if (iflag.eq.0) goto 9999
c                                 metastable try refining the increment
      order = order * 5.d0
      div = div / order

      v(ivi) = vo(ivi) + div
      v(ivd) = vo(ivd)
      call incdp0
 
100   call univeq (ivd,ier)

      if (ier.eq.0) goto 90
c                                 error from univeq, more than once quit
      goto (9400), icter
      div =odiv
      icter = icter + 1
      goto 100
c                                 no error from univeq
90    if (order.gt.10.d6) goto 9500
      goto 70
c                                 all systems go
9999  call assptx
      if (dabs(div).lt.dv(ivi)) then
         div = dv(ivi)*div/dabs(div)
      end if 

      if (v(ivi).eq.vmax(ivi).or.v(ivi).eq.vmin(ivi)) then
         ikwk= 1
         call maxend
      end if 

      return

9400  call warn (87,dv(1),ivi,'WWAY')
      ier = 1
      return  
       
9500  call warn (24,dv(ivi)/order,ivi,'WWAY')
      ier = 1
      return

9600  call warn  (63,dv(1),ivi,'WWAY')     
      ier = 1
      return
  
9800  call warn (58,dv(1),ier,'WWAY')
      if (ivct.gt.4) goto 9820
      write (n3,1050) ird,(vnu(l),names(idr(l)), l= 1, ivct)
      goto 9830
9820  write (n3,1050) ird,(vnu(l),names(idr(l)), l= 1, 4)
      write (n3,1060) (vnu(l),names(idr(l)), l= 5, ivct)
9830  write (n3,*) 
      ier = 1

1050  format (1x,'(',i6,')',4(1x,g9.3,1x,a))
1060  format (6x,4(1x,g9.3,1x,a),/,6x,4(1x,g9.3,1x,a))
      end

      subroutine fultxt (nchar,rtxt)
c-----------------------------------------------------------------------
c subprogram to write text for complete reaction reactions
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      character*8 names, text(kd2,2)*1, rtxt(kd2)*1,
     *            ppart(k8)*31, rxnstr*(kd2), char8, 
     *            mpart(k8)*31,cname*5, exten(2)*8, 
     *            ptext*(kd2), mtext*(kd2)

      integer jchar(2),ip,im,i,id,j,nchar

      common/ csta4 /cname(k5)/ cst104 /rxnstr(k2)/ cst8 /names(k1)

      integer iffr,isr
      double precision vuf,vus
      common/ cst201 /vuf(2),vus(h5),iffr,isr

      integer jfct,jmct,jprct
      common/ cst307 /jfct,jmct,jprct

      integer iff,idss,ifug,ifyn,isyn
      common/ cst10  /iff(2),idss(h5),ifug,ifyn,isyn

      integer idr,ivct
      double precision vnu
      common/ cst25 /vnu(k7),idr(k7),ivct

      integer jds,ifr
      double precision du,dv
      common/ cst21 /du(2),dv(2),jds(h5),ifr

      integer irct,ird
      double precision vn
      common/ cst31 /vn(k2,k7),irct,ird

      integer isec,icopt,ifull,imsg,io3p
      common/ cst103 /isec,icopt,ifull,imsg,io3p

      integer ikp
      common/ cst61 /ikp(k1)

      integer ids,isct,icp1,isat,io2
      common/ cst40 /ids(h5,h6),isct(h5),icp1,isat,io2
c----------------------------------------------------------------------

      save exten 
c                                 this is a bullshit trick and
c                                 will cause errors if someone
c                                 uses a function other than G.
      data exten/'-V(j/b) ','S(j/k) '/

      ip = 0 
      im = 0

      do i = 1, 8
         mpart(i) = '                              '
         ppart(i) = '                              '
      end do
c                             load id's of phases with + or - coeffs 
      do i = 1, ivct 
c           
         id = idr(i)   
         if (vnu(i).lt.0d0) then
            im = im + 1
            call wrpart (-vnu(i),ikp(id),names(id),mpart(im))
         else 
            ip = ip + 1
            call wrpart (vnu(i),ikp(id),names(id),ppart(ip))
         end if
      end do
c                             composant stoichiometry: 
      do i = 1, isat 
         if (vus(i).gt.0d0) then
            im = im + 1
            call wrpart (vus(i),0,names(jds(i)),mpart(im))
         else if (vus(i).lt.0) then
            ip = ip + 1
            call wrpart (-vus(i),0,names(jds(i)),ppart(ip))
         end if 
      end do 
c                                 fluid stoichiometry:
      if (ifyn.eq.0) then 
         do i = 1, 2
            if (iff(i).ne.0) then
               if (vuf(i).gt.0d0) then
                  im = im + 1
                  call wrpart (vuf(i),0,names(i),mpart(im))
               else if (vuf(i).lt.0d0) then
                  ip = ip + 1
                  call wrpart (-vuf(i),0,names(i),ppart(ip))
               end if 
            end if 
         end do
      end if 
c                                 mobile components: 
      do i = 1, jmct
         char8 = cname(jprct+i)
         if (du(i).gt.0d0) then
            im = im + 1
            call wrpart (du(i),0,char8,mpart(im))
         else if (du(i).lt.0d0) then
            ip = ip + 1
            call wrpart (-du(i),0,char8,ppart(ip))
         end if
      end do

      if (ifull.eq.1.or.ifull.eq.3) goto 55

      do i = 1, 2
         if (dv(i).gt.0d0) then
            im = im + 1
            call wrpart (dv(i),0,exten(i),mpart(im))
         else if (dv(i).lt.0d0) then
            ip = ip + 1
            call wrpart (-dv(i),0,exten(i),ppart(ip))
         end if
      end do 

55    if (ip.gt.k8.or.im.gt.k8) call error (997,du(1),ip,'FULTXT')

      write (ptext,1000) (ppart(i), i = 1, ip)
      write (mtext,1000) (mpart(i), i = 1, im)

      jchar(1) = im * 31
      read (mtext,1010) (text(i,1), i = 1, jchar(1))
      jchar(2) = ip * 31
      read (ptext,1010) (text(i,2), i = 1, jchar(2))
c                             filter out blanks:
      do i = 1, 2
         nchar = 1
         do j = 2, jchar(i)
            if (text(j-1,i).eq.' '.and.text(j,i).eq.' ') cycle
            nchar = nchar + 1
            text(nchar,i) = text(j,i)
         end do 
c                             filter out blanks before parentheses:
         jchar(i) = nchar
         nchar = 0
         do j = 1, jchar(i)-1
            if ((text(j,i).eq.' '.and.text(j+1,i).eq.')').or.
     *          (text(j,i).eq.' '.and.text(j+1,i).eq.'(')) cycle
            nchar = nchar + 1
            text(nchar,i) = text(j,i)
         end do  
         jchar(i) = nchar
      end do
c                             concatenate parts into rxnstr:
      if (jchar(1)+jchar(2).gt.(k8)*31-3) 
     *                    call error (208,du(1),ip,'FULTXT')
 
      nchar = 0 

      do i = 1, 2
         do j = 1, jchar(i)
            nchar = nchar + 1
            rtxt(nchar) = text(j,i)
         end do
         if (i.eq.2) cycle
         nchar = nchar + 3
         rtxt(nchar-2) = ' '
         rtxt(nchar-1) = '='
         rtxt(nchar)   = ' '
      end do 
c                             this is redundant if the reaction 
c                             was found before
      write (rxnstr(ird),1010) (rtxt(i),i=1,nchar)

1000  format (14(a,1x))
1010  format (434a1)

      end 


      subroutine stetxt 
c-----------------------------------------------------------------------
c subprogram to write complete reactions to a file called 
c "reaction_list.dat" for import into excel, etc. 
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision rcoef(k8)

      integer i,ict

      character*8 rname(k8), exten(2), names, cname*5

      common/ csta4 /cname(k5)/ cst8 /names(k1)

      integer iffr,isr
      double precision vuf,vus
      common/ cst201 /vuf(2),vus(h5),iffr,isr

      integer jfct,jmct,jprct
      common/ cst307 /jfct,jmct,jprct

      integer jds,ifr
      double precision du,dv
      common/ cst21 /du(2),dv(2),jds(h5),ifr

      integer iff,idss,ifug,ifyn,isyn
      common/ cst10  /iff(2),idss(h5),ifug,ifyn,isyn

      integer idr,ivct
      double precision vnu
      common/ cst25 /vnu(k7),idr(k7),ivct

      integer irct,ird
      double precision vn
      common/ cst31 /vn(k2,k7),irct,ird

      integer ivarrx,ivarip,isudo,ivar
      common/ cst62 /ivarrx(k2),ivarip(k2),isudo,ivar

      integer ids,isct,icp1,isat,io2
      common/ cst40 /ids(h5,h6),isct(h5),icp1,isat,io2

      save exten 
c                                 this is a bullshit trick and
c                                 will cause errors if someone
c                                 uses a function other than G.
      data exten/'-V(j/b) ','S(j/k) '/
c                                 ouput data as follows
c                                 V, S, mobile components, fluid components, 
c                                 saturated component, thermodynamic components
c                                 extensities:
      do i = 1, 2
         rcoef(i) = -dv(i)
         rname(i) = exten(i)
      end do 

      ict = 2
c                                 mobile components:
      do i = 3, 2 + jmct
         rcoef(i) = -du(i-2)
         rname(i) = cname(jprct+i-2)
      end do

      ict = ict + jmct
c                                 saturated phase components:
      if (ifyn.eq.0) then 
         do i = 1, 2
            if (iff(i).ne.0) then
               ict = ict + 1
               rcoef(ict) = -vuf(i)
               rname(ict) = names(i) 
            end if 
         end do
      end if 
c                                 saturated components: 
      do i = 1, isat 
         if (vus(i).ne.0d0) then 
            ict = ict + 1
            rcoef(ict) = -vus(i)
            rname(ict) = names(jds(i))
         end if 
      end do
c                                 thermodynamic components
      do i = 1, ivct         
            
         ict = ict + 1
         rcoef(ict) = vnu(i)
         rname(ict) = names(idr(i))
      end do
 
      write (n6,1000) irct,ivar,(rname(i), rcoef(i), i = 1, ict)

1000  format (i5,1x,i1,1x,20(a,1x,g14.8,1x))

      end 

      subroutine wrpart (vnu,ikp,name,part)
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ikp

      double precision vnu

      character*8 name, fname*10, part*30

      common/ csta7 /fname(h9)

      integer isec,icopt,ifull,imsg,io3p
      common/ cst103 /isec,icopt,ifull,imsg,io3p

      if (ikp.eq.0) then
         if (ifull.gt.2) then  
            write (part,1000) vnu,name
         else 
            write (part,2000) name
         end if 
      else 
         if (ifull.gt.2) then 
            write (part,1010) vnu,fname(ikp),name
         else 
            write (part,2010) fname(ikp),name
         end if 
      end if 

1000  format (g9.3,1x,a)
1010  format (g9.3,1x,a,'(',a,')')
2000  format (a)
2010  format (a,'(',a,')')

      end

      subroutine detest (ier)
c----------------------------------------------------------------------
c subroutine used by lpopt to test if an assemblage is non-degenerate
c ier = 0 if non-degenerate, otherwise 1.
c----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer ier,i,j
 
      integer ipvt
      double precision a,b
      common/ cst301 /a(k5,k5),b(k5),ipvt(k5) 

      integer hcp,idv
      common/ cst52  /hcp,idv(k7) 

      double precision cp
      common/ cst12 /cp(k5,k1)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  
c----------------------------------------------------------------------
      ier = 0
c                                load the transpose of the
c                                concentration matrix of the pseudo-
c                                invariant assemblage.
      do i = 1, icp
         do j = 1, icp
            a(j,i) = cp(j,idv(i))
         end do
      end do  
c                                factor the matrix
      call factr1 (icp,ier)

      end 



      subroutine money
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'
      
      integer imoney, ierr

      open (n8, file = 'money', iostat = ierr, status = 'old')
      if (ierr.ne.0) then
c                                 system could not find the file
         write (*,1000) 
         stop
1000     format (/,'There is no money file, make one and',
     *             ' deposit at least a dime.',/)
 
      end if

      read (n8,*) imoney

      rewind (n8)

      if (imoney.le.0) then

         write (*,1010) 
         stop
1010     format (/,'You are out of money, deposit at least a dime'/)

      else 

         imoney = imoney - 1
         write (n8,*) imoney

      end if

      close (n8)

      end

      subroutine testit (ier)
c-----------------------------------------------------------------------
c ier = 0 the plane is stable
c ier = 1 the plane is metastable
c ier = 2 the assemblage is degenerate
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i,j,ier,k

      double precision delt,dtol,utol,ptol,dg
 
      common/ cst87 /delt(l2),dtol,utol,ptol

      double precision g
      common/ cst2 /g(k1)

      double precision cp
      common/ cst12 /cp(k5,k1)
 
      double precision a,b
      integer ipvt,idv,iophi,idphi,iiphi,iflg1
      common/ cst23  /a(k8,k8),b(k8),ipvt(k8),idv(k8),iophi,idphi,
     *                iiphi,iflg1

      integer h,id
      common/ cst52 /h,id(k7)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  
c                              optional routine for testing
c                              for metastable g-x planes.
      ier = 0

      call abload (*991)
 
      do j = istct, iphct
 
         dg = g(j)
 
         do i = 1, icp
            dg = dg - cp(i,j) * b(i)
         end do
 
         if (dg.gt.dtol) cycle
c                                check that a phase is not metastable
c                                with respect to itself, this
c                                could be remedied by changing the 
c                                value of dtol.
         do k = 1, icp
            if (id(k).eq.j) goto 10
         end do
c                                the assemblage is metastable with
c                                respect to phase j
         ier = 1
 
10    end do 

      goto 99

991   ier = 2
 
99    end

      subroutine outtit
c-----------------------------------------------------------------------
c outtit writes title information and a brief description of the
c chemical system for each calculation requested.
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'
 
      character*8 cname*5,names,vname,xname

      integer i,j,jasmbl,iasmbl,itot,ic11

      double precision ctot

      character*162 title
      common/ csta8 /title(4)
 
      common/ cst3  /ctot(k1)/ cst8 /names(k1)
     *      / csta4 /cname(k5)
     *      / csta2 /xname(k5),vname(l2)
     *      / cst27 /jasmbl(j9),iasmbl(j9)

      character cmpnt*5, dname*40
      common/ csta5 /dname,cmpnt(k0)

      integer ixct,iexyn,ifact
      common/ cst37 /ixct,iexyn,ifact 

      character*8 exname,afname
      common/ cst36 /exname(h8),afname(2)

      double precision cp
      common/ cst12 /cp(k5,k1)

      integer ifct,idfl
      common/ cst208 /ifct,idfl

      integer iff,idss,ifug,ifyn,isyn
      common/ cst10  /iff(2),idss(h5),ifug,ifyn,isyn

      integer ipot,jv,iv1,iv2,iv3,iv4,iv5
      common/ cst24 /ipot,jv(l2),iv1,iv2,iv3,iv4,iv5

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  

      integer ids,isct,icp1,isat,io2
      common/ cst40 /ids(h5,h6),isct(h5),icp1,isat,io2

      logical gflu,aflu,fluid,shear,lflu,volume
      common/ cxt20 /gflu,aflu,fluid(k5),shear,lflu,volume
c----------------------------------------------------------------------
      write (n3,1000)
c                          title:
      write (n3,1190) title(1)
c                          data base
      write (n3,1210) dname
c                          fluid
      if (ifyn.eq.0.or.gflu) call rfluid (2,ifug) 
c                          independent potentials:
      write (n3,1070) (vname(jv(j)), j = 1, ipot)
c                          saturated phase components:
      if (ifyn.eq.0) then 
         j = icp + isat
         write (n3,1200) (cname(j+i), i = 1, ifct)
      end if 
c                          saturated components:
      if (isyn.eq.0) then 
         j = icp + isat
         write (n3,1180) (cname(i), i = icp1, j)
      end if 
c                          unconstrained components
      write (n3,1080) (cname(i), i = 1, icp)
c                          phases
      if (icp.gt.3) then 
         write (n3,1150) (cname(i), i = 1, icp)
         do i = istct, iphct
            write (n3,1160) names(i), (cp(j,i)/ctot(i), j = 1, icp)
         end do 
      else if (icp.eq.3) then 
         write (n3,1090) (cname(i), i = 2, 3)
         write (n3,1100) (names(i), cp(2,i)/ctot(i), cp(3,i)/ctot(i),
     *                                             i = istct, iphct)
      else if (icp.eq.2) then 
         write (n3,1040) cname(2)
         write (n3,1030) (names(i), cp(2,i)/ctot(i), i = istct, iphct)
      else if (icp.eq.1) then 
         write (n3,1130)
         write (n3,1110) (names(i), i = istct, iphct)
      end if 
c                          saturation composant phases,
c                          load the indexes into iasmbl:
      if (isat.gt.0) then 
         itot = 0
         do i = 1, isat
            ic11 = isct(i)
c                          do a check that there is at least
c                          one composant for each saturated component:
c                          why here? well, why not?
            if (ic11.eq.0) call error (279,cp(1,1),i,'OUTTIT')
 
            do j = 1, ic11
               iasmbl(itot+j) = ids(i,j)
            end do 
            itot = itot+ic11
         end do 
 
         write (n3,1170)
         write (n3,1110) (names(iasmbl(i)), i = 1, itot)
      end if 
c                          excluded phases
      if (iexyn.eq.0) then 
         write (n3,1140)
         write (n3,1110) (exname(i), i = 1, ixct)
      end if 
 
      write (n3,1000)
 
1000  format (/,80('-'),/)
1030  format (4(2x,a8,1x,f5.3))
1040  format (/,'Phases and (projected) mol fraction ',a,':',/)
1070  format (/,'Independently constrained potentials:',//,3x,8(a,1x))
1080  format (/,'Components with unconstrained potent'
     *       ,'ials:',//,3x,10(a5,3x))
1090  format (/,'Phases and (projected) composition with respect to '
     *         ,a5,' and ',a5,':',/)
1100  format (3(1x,a,2x,f5.3,2x,f5.3,5x))
1110  format (7(1x,a,1x))
1130  format (/,'Phases:',/)
1140  format (/,'Excluded phases:',/)
1150  format (/,'Phases and (projected) compositions:',//,
     *        11x,12(1x,a5,1x),/)
1160  format (3x,a8,12(1x,f5.3,1x))
1170  format (/,'Phases on saturation and buffering surfaces:',/)
1180  format (/,'Saturated or buffered components:',//,3x,7(a,3x))
1190  format (/,'Problem title: ',a,/)
1200  format (/,'Saturated phase components:',//,3x,5(a,3x))
1210  format ('Thermodynamic data base from: ',a)
      end

      subroutine findas
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      character*5 cname

      integer kmax(k1),i,j

      double precision r,dg
      
      integer hcp,id
      common/ cst52  /hcp,id(k7)

      double precision g
      common/ cst2 /g(k1)

      double precision a,u
      integer ipvt,idv,iophi,idphi,iiphi,iflg1
      common/ cst23  /a(k8,k8),u(k8),ipvt(k8),idv(k8),iophi,idphi,
     *                iiphi,iflg1

      double precision cp
      common/ cst12 /cp(k5,k1)
       
      common/ csta4  /cname(k5)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  
c----------------------------------------------------------------------
      do i = istct, iphct
c                                 sort phases into subcompositions
         do j = 1, icp
            if (cp(j,i).gt.1d-5) kmax(i) = j
         end do
      end do

      do hcp = 1, icp
c                             formerly ic(hcp) = hcp, 9/27/08.
         id(hcp) = hcp
c                             make the first guess for the
c                             stable hcp component phase:
         do i = istct, iphct
            if (kmax(i).eq.hcp) then
               id(hcp) = i
               goto 20
            end if
         end do
c                             missing composant error
         call error (15,r,i,cname(hcp))
c
20       call abload (*992)
c                             start test loop:
         do j = istct, iphct 
            if (kmax(j).gt.hcp) cycle
            dg = g(j)
            do i = 1, hcp
               if (id(i).eq.j) goto 10
               dg = dg - cp(i,j)*u(i)
            end do

            if (dg.lt.-1d-8) then
c                             plane is metastable
c                             replace the most recently
c                             added test phase:
               id(hcp) = j
               call abload (*992)
            end if
10       end do 

      end do
c                                 this is to reset the component
c                                 counter in ABLOAD via cst52. 
      hcp = icp

      goto 99

992   call error (999,dg,i,'FINDAS')

99    end 
         

      subroutine checkd (jd)
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer idold,jd,i

      double precision dgphc

      integer hcp,id
      common / cst52 /hcp,id(k7)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  
c----------------------------------------------------------------------
c                           istabl = 1 if the original
c                           assemblage was metastable, istable = 0
c                           otherwise.
      do i = 1, icp 
         if (id(i).eq.jd) goto 99
      end do

      if (dgphc(jd).gt.-1d-05) goto 99

      idold = id(icp)
      id(icp) = jd 

      call abload (*110)

      goto 99
c                           restore old vertex
110   id(icp) = idold

      call abload (*9000)
      
      goto 99

9000  call error (19,1d0,jd,'CHECKD')

99    end 

      double precision function dgphc (ld)
c-----------------------------------------------------------------------
c dgphc is the g difference between a phase
c identified by 'ld', and the plane defined by 'u'.
c   referenced by: 
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i,is,ld

      double precision dg,a,u
 
      common/ cst23 /a(k8*k8),u(k8),is(2*k8+4)

      double precision g
      common/ cst2 /g(k1)

      double precision cp
      common/ cst12 /cp(k5,k1)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  
 
      dg = g(ld)

      do i = 1, icp 
         dg = dg - cp(i,ld)*u(i)
      end do 

      dgphc = dg

      end

      subroutine miscb0 (idim,ntot,np,solvs1,id)
c----------------------------------------------------------------------
c miscb0 counts the number of distinct phases in an equilibrium of
c ntot pseudocompounds, the routine returns solvus = .true. if a solvus
c occurs in one or more of the phases.

c     np  - is the number of phases

c this routine is unecessarily complicated, because it assumes
c pseudocompounds are not ordered by solution (but they are
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical solvs1, solvus 
c                                 -------------------------------------
c                                 local variables
      integer idsol(k8,k8),jdsol(k8),ids,np,np1,np2,ncpd,idim,
     *        i,j,kdsol(k8),k,ntot,id(idim)

      integer ikp
      common/ cst61 /ikp(k1)
c-----------------------------------------------------------------------
      np = 0   
      ncpd = 0 
      solvs1 = .false.

      do 10 i = 1, ntot
         ids = ikp(id(i))
         if (ids.le.0) then 
            ncpd = ncpd + 1
            cycle
         else   
            do j = 1, np
               if (ikp(idsol(1,j)).eq.ids) then
                  jdsol(j) = jdsol(j) + 1
                  idsol(jdsol(j),j) = id(i)
                  goto 10 
               end if 
            end do                   
            np = np + 1
            jdsol(np) = 1
            idsol(1,np) = id(i)
         end if
10    continue  
c                                now check if the np solutions
c                                are homogeneous
      np1 = 0 

      do i = 1, np

         np2 = 0 

         if (jdsol(i).eq.1) then 
c                                 solution represented by a single
c                                 compound
            np1 = np1 + 1
            cycle 

         else 
c                                 check for immiscibility
            do 20 j = 1, jdsol(i)
c                                 how many phases are there of this
c                                 solution?
               do k = 1, np2
                  if (.not.solvus(kdsol(k),idsol(j,i))) then
                     goto 20 
                  else
                     solvs1 = .true.
                  end if 
               end do 
c                                 if here, it's a new phase
               np2 = np2 + 1
               kdsol(np2) = idsol(j,i)
               
20          continue
             
            np1 = np1 + np2
   
         end if
 
      end do 

      np = ncpd + np1

      end 

      subroutine findjp (ivi,ivd,odv,ip)
c----------------------------------------------------------------------
c findjp - called by sfol1 to iteratitvely refine the location of 
c an invariant point, the invariant point is located within an 
c energy tolerence specified by the variable ptol in routine pchk
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer iwarn,ier,ip,ivi,ivd,igo

      double precision delt,dtol,utol,ptol,odv,dg

      common/ cst87 /delt(l2),dtol,utol,ptol 

      integer iflag
      common/ cst7 /iflag

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      double precision vmax,vmin,dv
      common/ cst9  /vmax(l2),vmin(l2),dv(l2)
c----------------------------------------------------------------------    
c                                 iterate independent variable to refine
c                                 the conditions of the invariant point:
      iwarn = 0 
      ier = 0 
      ip = 0
c                                 step back to stable condition
20    call reptx
      odv = odv / 2d0
c                                 check if the point is already known
      call sameip (ip)

      if (ip.ne.0) goto 99

      if (dabs(odv).lt.delt(ivi)) goto 30
 
10    v(ivi) = v(ivi) + odv 
      call incdep (ivi)   

      call univeq (ivd,ier)
      if (ier.ne.0) goto 30 

      call pchk (dg,igo)

      if (igo.eq.1) then
         goto 40 
      else if (iflag.eq.1) then
         goto 20
      else  
         if (v(ivi).gt.vmax(ivi).or.v(ivi).lt.vmin(ivi)) goto 30
         goto 10
      end if 
c                                 univeq blew up, assign ip with 
c                                 warning:
30    call reptx 
      iwarn = 1 
c                                 assign the invariant assemblage:
40    call assip (ip)
      
      if (iwarn.eq.1) call warn (47,ptol,ip,'FINDJP')

99    end 

      subroutine ssaptx
c---------------------------------------------------------------------
c ssaptx assigns the current value of v1 and v2 to the array ptx and
c increments the ordinate counter ipt2 if the value of either variable
c has increased by more than 1% of the default increment. 
c---------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer ipt2
      double precision ptx
      common/ cst32 /ptx(l5),ipt2

      integer ipot,jv,iv1,iv2,iv3,iv4,iv5
      common/ cst24 /ipot,jv(l2),iv1,iv2,iv3,iv4,iv5

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      double precision vmax,vmin,dv
      common/ cst9  /vmax(l2),vmin(l2),dv(l2)

      if (dabs((v(iv1)-ptx(ipt2-1))/dv(iv1)).gt.1d-2.or.
     *    dabs((v(iv2)-ptx(ipt2))/dv(iv2)).gt.1d-2) then 
  
         ipt2 = ipt2 + 2

         if (ipt2.gt.l5) ipt2 = l5

         ptx(ipt2-1) = v(iv1)
         ptx(ipt2)   = v(iv2)

      else 

         ptx(ipt2-1) = v(iv1)
         ptx(ipt2)   = v(iv2)

      end if 
 
      end

      subroutine assip (ip)
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer j,ip,isol,np,icp2

      logical solvs1

      integer ipid,ipct
      common/ cst29 /ipid(k2,k8),ipct
 
      double precision vip
      common/ cst28 /vip(l2,k2)

      integer ipvt,idv,iophi,idphi,iiphi,iflg1
      double precision a,b
      common/ cst23 /a(k8,k8),b(k8),ipvt(k8),idv(k8),
     *               iophi,idphi,iiphi,iflg1

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  

      integer ivarrx,ivarip,isudo,ivar
      common/ cst62 /ivarrx(k2),ivarip(k2),isudo,ivar

      integer ikp
      common/ cst61 /ikp(k1)
c----------------------------------------------------------------------
c                                 assign the new elements of ipid
      ipct = ipct + 1
      if (ipct.gt.k2) call error (182,r,k2,'ASSIP')

      icp2 = icp + 2

      idv(icp+1) = iophi
      idv(icp2) = idphi

      isol = 0 

      do j = 1, icp2
         ipid(ipct,j) = idv(j)
         if (ikp(idv(j)).gt.0) isol = isol + 1
      end do 
c                                 get true variance
      if (isol.gt.1) then 
c                                 classify the reaction
         if (isol.gt.1) call miscb0 (k7,icp2,np,solvs1,idv)
c                                 set ivar, just in case it's used someplace
         ivarip(ipct) = icp2 - np

      else 

         ivarip(ipct) = 0 

      end if 

c                                 save conditions
      vip(1,ipct) = v(1)
      vip(2,ipct) = v(2)
      vip(3,ipct) = v(3)
      vip(4,ipct) = v(4)
      vip(5,ipct) = v(5)

      ip = ipct

      end

      subroutine sameip (ip)
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ip,icp2,jps,ips,i,j,k

      integer ipid,ipct
      common/ cst29 /ipid(k2,k8),ipct
 
      double precision vip
      common/ cst28 /vip(l2,k2)

      integer ipvt,idv,iophi,idphi,iiphi,iflg1
      double precision a,b
      common/ cst23 /a(k8,k8),b(k8),ipvt(k8),idv(k8),
     *               iophi,idphi,iiphi,iflg1

      integer iflag
      common/ cst7 /iflag

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  

      integer ipot,jv,iv
      common/ cst24 /ipot,jv(l2),iv(l2)

      double precision vmax,vmin,ddv
      common/ cst9  /vmax(l2),vmin(l2),ddv(l2)
c-----------------------------------------------------------------------
      ip = 0
c
      if (ipct.eq.0) goto 99
 
      icp2 = icp + 2
      idv(icp+1) = iophi
      idv(icp2) = idphi
c                                 determine if a new ip was found
      jps = 1
5     ips = jps 

      do 10 i = ips, ipct

         do 20 j = 1, icp2
            do 30 k= 1, icp2
30             if (idv(k).eq.ipid(i,j)) goto 20
c                                 check if invariant point is equivalent
c                                 to the ith ip already found
               goto 10
20       continue
c                      the new ip assemblage is equivalent to one found 
c                      earlier, check if the conditions match:
         do 25 j = 1, 5
            if (dabs(vip(j,i)-v(j)).gt.ddv(j).and.
     *                                 ddv(j).ne.0d0) then
               if (i.lt.ipct) then 
                  jps = i + 1
                  goto 5
               end if 
               goto 99 
            end if 
25       continue 
c                                 conditions and assemblage match
         v(iv(1)) = vip(iv(1),i)
         v(iv(2)) = vip(iv(2),i)
         call incdp0
         call assptx
c                                 reset iflag so sfol2 doesn't try 
c                                 anything stupid
         iflag = 0
         ip = i
         goto 99

10    continue

99    end

      subroutine lpopt (i,j,idead,output)
c-----------------------------------------------------------------------
c lpopt - is a shell to maintain backwards compatibility used for gridded
c minimization calculations. it calls lpopt0 to do a minimization, then 
c does some minor bookkeeping for node i,j. 
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical output

      integer i,j,k,idead

      integer igrd
      common/ cst311 /igrd(l7,l7)

      double precision a,b,c
      common/ cst313 /a(k5,k1),b(k5),c(k1)

      integer iap,ibulk
      common/ cst74  /iap(k2),ibulk

      integer hcp,idv
      common/ cst52  /hcp,idv(k7)  

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)
c-----------------------------------------------------------------------
c                                 check for positive bulk
      idead = 0 

      do k = 1, hcp
         if (b(k).gt.0d0) then 
            cycle
         else if (dabs(b(k)).lt.nopt(11)) then
            b(k) = 0d0
         else 
            idead = 2
            exit 
         end if 
      end do 

      if (idead.eq.0) call lpopt0 (idead)
c                                 if idead = 0 optimization was ok
      if (idead.eq.0) then 
c                                 at this point the compositions of
c                                 the np solutions are in cp3, ctot3, x3 indexed
c                                 by np and the original indices of the 
c                                 ncpd compounds are -kkp(np+1..ntot), and
c                                 the molar amounts of the phases are in b.
         call sorter (igrd(i,j),i,j,output)

      else 
 
         if (idead.ne.2) write (*,1000) i,j
         igrd(i,j) = k2
         iap(k2) = k3

      end if  

1000  format (/,'WARNING: minimization failed at node i=',i6,' j= ',i6,/
     *         'the null assemblage will be assigned to this node',/)

      end 

      subroutine outdat
c----------------------------------------------------------------------------
c output data for special uses
c----------------------------------------------------------------------------
      implicit none

      character*8 names

      include 'perplex_parameters.h'

      integer icomp,istct,iphct,icp,i,j

      double precision ctot

      common/ cst6   /icomp,istct,iphct,icp   
     *      / cst8  /names(k1)/ cst3  /ctot(k1)

      double precision g
      common/ cst2 /g(k1)

      double precision cp
      common/ cst12 /cp(k5,k1)

      integer jbulk
      double precision cblk
      common/ cst300 /cblk(k5),jbulk
          
1000  format (12(g16.8,1x))  
                   
      open (n8,file='xsystem.dat')

      do j = 1, icp
         write (n8,1000) cblk(j)
      end do 

      close (n8)
      
      open (n8,file='phase.dat')

      do i = 1, iphct
         write (n8,1000) names(i),g(i)/ctot(i),
     *                    (cp(j,i)/ctot(i),j=1,icp)
      end do

      close (n8)
      
      stop
      end 

      subroutine amihot (i,j,jhot,jinc)
c--------------------------------------------------------------------
c check if cell with lower left corner at node ij is homogeneous
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer jinc, i, j, jhot

      integer iap,ibulk
      common/ cst74 /iap(k2),ibulk

      integer igrd
      common/ cst311 /igrd(l7,l7)

      jhot = 1

      if (iap(igrd(i,j)).eq.iap(igrd(i,j+jinc)).and.
     *    iap(igrd(i,j)).eq.iap(igrd(i+jinc,j+jinc)).and.
     *    iap(igrd(i,j)).eq.iap(igrd(i+jinc,j))) jhot = 0
     
      end  

      subroutine filler (i,j,kinc)
c--------------------------------------------------------------------
c fill nodes between identical edges or diagonals, this could
c be modified to a less biased (symmetrical) strategy as in 
c routine aminot. 
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer kinc, i, j, ii, jj, k

      integer igrd
      common/ cst311 /igrd(l7,l7)

      if (kinc.eq.1) return 

      if (igrd(i,j).eq.igrd(i+kinc,j+kinc)) then
c                                right diagonal
         do k = 1, kinc-1
c                                the conditional prevents 
c                                overwriting of identities if 
c                                compression is off and searching
c                                for true boundaries.
            if (igrd(i+k,j+k).eq.0) igrd(i+k,j+k) = igrd(i,j)
         end do 
      else if (igrd(i+kinc,j).eq.igrd(i,j+kinc)) then
c                                left diagonal
         do k = 1, kinc-1
c                                the conditional prevents 
c                                overwriting of identities if 
c                                compression is off and searching
c                                for true boundaries.
            if (igrd(i+k,j+kinc-k).eq.0) 
     *          igrd(i+k,j+kinc-k) = igrd(i,j+kinc)
         end do 

c        write (*,*) 'left diagonal'

      end if
c                                bottom and top edge
      do jj = j, j+kinc, kinc
         if (igrd(i,jj).eq.igrd(i+kinc,jj)) then 
            do k = 1, kinc-1
c                                the conditional prevents 
c                                overwriting of identities if 
c                                compression is off and searching
c                                for true boundaries.
               if (igrd(i+k,jj).eq.0) igrd(i+k,jj) = igrd(i,jj)
            end do 
         end if
      end do 
c                                left and right edges
      do ii = i, i+kinc, kinc
         if (igrd(ii,j).eq.igrd(ii,j+kinc)) then 
            do k = 1, kinc-1
c                                the conditional prevents 
c                                overwriting of identities if 
c                                compression is off and searching
c                                for true boundaries.
               if (igrd(ii,j+k).eq.0)  igrd(ii,j+k) = igrd(ii,j)
            end do 
         end if
      end do 

      end

      subroutine wav2d1 (output)
c--------------------------------------------------------------------
c wav2d does constrained minimization on a 2 dimensional multilevel
c grid (ith column of a 2-d grid), lowest resolution is the default, 
c and phase boundaries are located at the highest level

c the phase assemblage at node(i,j) of the grid is identified by the 
c pointer to a phase assemblage.

c the compound grid data igrd is output to the "plot" file.
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical output 

      integer kinc, jinc(l8), iind(4), jind(4), iiind(4,2),htic,
     *        jjind(4,2), hotij(l7*l7,2), kotij(l7*l7,2), icind(4),
     *        jcind(4), lhot(4), ieind(5), jeind(5),i,ihot,
     *        iil,jjl,ll,je,ie,icell,
     *        ii,jj,kk,hh,hhot,jcent,icent,jjc,iic,h,jtic,khot,k,
     *        ktic,j,jhot,klow,kinc2,kinc21,idead

      integer ipot,jv,iv1,iv2,iv3,iv4,iv5
      common/ cst24 /ipot,jv(l2),iv1,iv2,iv3,iv4,iv5

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      integer igrd
      common/ cst311 /igrd(l7,l7)

      integer iap,ibulk
      common/ cst74 /iap(k2),ibulk

      integer jlow,jlev,loopx,loopy,jinc1
      common/ cst312 /jlow,jlev,loopx,loopy,jinc1

      integer idasls,iavar,iasct,ias
      common/ cst75  /idasls(k5,k3),iavar(3,k3),iasct,ias

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      save iind, jind, iiind, jjind, icind, jcind

      data iind, jind   /0,0,1,1, 0,1,1,0/
      data iiind, jjind /0,0,1,1,1,1,2,2, 1,1,2,0,0,2,1,1/
      data icind, jcind /0,1,2,1, 1,2,1,0/
      data ieind, jeind /0,0,2,2,0, 0,2,2,0,0/
c-----------------------------------------------------------------------
c                               initialize assemblage counter
      iasct = 0 
      ibulk = 0 
c                               call old routine for 1d grid
      if (loopx.eq.1) then
         call wavgrd (output)
         return
      end if 
c                               load arrays for lp solution
      call initlp 
c                               jlow is the number of nodes
c                               at the lowest level, the number of
c                               nodes is
      klow = jlow - 1
c                               first level:
      loopy = klow * 2**(jlev-1) + 1 

      loopx = (loopx-1) * 2**(jlev-1) + 1 

      write (*,1040)

      if (loopy.gt.l7) then
         call warn (92,v(iv1),loopy,'y_node')
         klow = (l7 - 1)/2**(jlev-1)
         loopy = klow * 2**(jlev-1) + 1 
      end if          

      if (loopx.gt.l7) then
         call warn (92,v(iv1),loopx,'x_node')
         loopx = l7
      end if  
c                               initialize igrd (this is critical 
c                               for auto_refine).
      do j = 1, loopy
         do i = 1, loopx
            igrd(i,j) = 0
         end do 
      end do 
c                               could check here if loopx*loopy, the
c                               theoretical max number of assemblages
c                               is > k2, but in practice the number of
c                               assemblages << k2, so only test when 
c                               actually set.
      do j = 1, k2
         iap(j) = 0 
      end do 
c                               increments at each level
      do j = 1, jlev
         jinc(j) = 2**(jlev-j)
      end do 

      kinc = jinc(1)
      jinc1 = kinc

      call setvar 
c                               do all points on lowest level
      do i = 1, loopx, kinc
         do j = 1, loopy, kinc
            call setvr0 (i,j)
            call lpopt (i,j,idead,output)
         end do
c                               progress info
         write (*,1030) dfloat(i)/dfloat(loopx)*1d2
      end do 
c                               get hot points
      ihot = 0 
      kinc2 = kinc/2
      kinc21 = kinc2 + 1

      do i = 1, loopx-1, kinc
         do j = 1, loopy-1, kinc
            call amihot (i,j,jhot,kinc)
            if (jhot.ne.0) then 
               ihot = ihot + 1
               hotij(ihot,1) = i
               hotij(ihot,2) = j 
c                               cell is heterogeneous
c                               fill in homogeneous diagonals
c                               and edges
               if (iopt(18).ne.0.and.kinc.gt.1) call filler (i,j,kinc)
            else 
c                               cell is homogeneous
c                               fill in entire cell
               if (kinc.gt.1) call aminot (i,j,kinc,kinc2,kinc21)
            end if 
         end do 
      end do

      if (ihot.eq.0) goto 10 

      ktic = (loopx/kinc+1)*(loopy/kinc+1)

      if (jlev.gt.1) write (*,1050) 

      htic = 0 
c                              now refine on all higher levels:
      do k = 2, jlev
c                              set new hot cell counter
         khot = 0 
         jtic = 0 
c                              now working on new level
         kinc = jinc(k)
         kinc2 = kinc/2
         kinc21 = kinc2 + 1
c
         write (*,1060) ihot,k
c                              compute assemblages at refinement
c                              points
         do h = 1, ihot
c                              cell corner
            iic = hotij(h,1) 
            jjc = hotij(h,2)
c                              first compute the central node of
c                              the hot cell
            icent = iic + kinc
            jcent = jjc + kinc
            if (igrd(icent,jcent).eq.0) then 
               call setvr0 (icent,jcent)
               call lpopt (icent,jcent,idead,output)                       
               jtic = jtic + 1
               htic = htic + 1
            end if 
c                              now determine which of the diagonals
c                              has a change
            hhot = 0 

            do hh = 1, 4
            
               i = iic + iind(hh)*2*kinc
               j = jjc + jind(hh)*2*kinc
               lhot(hh) = 0 

               if (iap(igrd(i,j)).ne.iap(igrd(icent,jcent))) then 
c                              cell is hot
                  khot = khot + 1
                  hhot = hhot + 1
                  kotij(khot,1) = iic + iind(hh)*kinc
                  kotij(khot,2) = jjc + jind(hh)*kinc
                  lhot(hh) = 1
c                              compute assemblages at new nodes
                  do kk = 1, 2
                     ii = iic + iiind(hh,kk)*kinc
                     jj = jjc + jjind(hh,kk)*kinc
                     if (igrd(ii,jj).eq.0) then 
                        call setvr0 (ii,jj)
                        call lpopt (ii,jj,idead,output)
                        jtic = jtic + 1
                        htic = htic + 1
                     end if 
                  end do 
               end if 
            end do 
c                              if less than 3 hot sub-cells check
c                              edges
            if (hhot.lt.4.and.hhot.gt.1) then 

               do hh = 1, 4
c                              index the edge node
                  ii = iic + icind(hh)*kinc
                  jj = jjc + jcind(hh)*kinc
                  if (igrd(ii,jj).ne.0) then 
c                              could have a second hot cell, check
c                              both corners
                     icell = hh

                     do kk = 1, 2
                        ie = iic + ieind(hh+kk-1)*kinc
                        je = jjc + jeind(hh+kk-1)*kinc
                        icell = icell + kk - 1
                        if (icell.gt.4) icell = 1
                        
                        if (iap(igrd(ii,jj)).ne.iap(igrd(ie,je)).and.
     *                     lhot(icell).eq.0) then 
c                               new cell
                           khot = khot + 1
                           hhot = hhot + 1
                           lhot(icell) = 1
c                                cell index is 
                           ii = iic + iind(icell)*kinc
                           jj = jjc + jind(icell)*kinc
                           kotij(khot,1) = ii
                           kotij(khot,2) = jj
c                                compute assemblage at cell nodes
                           do ll = 1, 4
                              iil = ii + iind(ll)*kinc
                              jjl = jj + jind(ll)*kinc
                              if (igrd(iil,jjl).eq.0) then              
                                 call setvr0 (iil,jjl)
                                 call lpopt (iil,jjl,idead,output)
                                 jtic = jtic + 1
                                 htic = htic + 1
                              end if 
                           end do  
                        end if
                     end do 
                  end if 
               end do 
            end if      

            do hh = 1, 4

               i = iic + iind(hh)*kinc
               j = jjc + jind(hh)*kinc
               if (i.lt.loopx.and.j.lt.loopy) then 
                  if (lhot(hh).eq.0) then 
c                                fill cold cells
                     call aminot1 (icent,jcent,i,j,kinc)
                  else 
c                                fill hot cells
                     if (iopt(18).ne.0) call filler (i,j,kinc)
                  end if 
               end if 
            end do 
         
            if (htic.gt.500) then 
               write (*,1090) jtic
               htic = 0 
            end if
 
         end do

         write (*,1070) k,jtic
         ktic = ktic + jtic
 
         write (*,1080) ktic,(loopx/kinc+1)*(loopy/kinc+1)

         if (khot.eq.0.or.k.eq.jlev) exit 
c                             now switch new and old hot list
         ihot = khot
         do i = 1, khot
            hotij(i,1) = kotij(i,1)
            hotij(i,2) = kotij(i,2)
         end do 

      end do 
c                                 ouput grid data
10    if (output) call outgrd (loopx,loopy,jinc(1))    

1030  format (f5.1,'% done with low level grid.')
1040  format (/)
1050  format (/,'Beginning grid refinement stage.',/)
1060  format (i6,' grid cells to be refined at grid level ',i1)
1070  format (7x,'refinement at level ',i1,' involved ',i6,
     *        ' minimizations')
1080  format (i6,' minimizations required of the ',
     *        'theoretical limit of ',i6)
1090  format (7x,'...working (',i6,' minimizations done)')

      end 

      subroutine outgrd (loopx,loopy,jinc)
c----------------------------------------------------------------------
c output grid data to the plot file
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer loopx,loopy,jinc,i,j,jst,kst,kd,ltic

      integer igrd
      common/ cst311 /igrd(l7,l7)

      integer idasls,iavar,iasct,ias
      common/ cst75  /idasls(k5,k3),iavar(3,k3),iasct,ias
c----------------------------------------------------------------------

      write (n4,*) loopx, loopy, jinc
c                                 fill in grid
      do i = 1, loopx

         if (i.ne.1.and.igrd(i,1).eq.0) igrd(i,1) = igrd(i-1,1)

         kst = 1

20       jst = kst
         if (i.ne.1.and.igrd(i,jst).eq.0) igrd(i,jst) = igrd(i-1,jst)
         kd = igrd(i,jst)
         ltic = -1

         do j = jst, loopy

            if (i.ne.1.and.igrd(i,j).eq.0) igrd(i,j) = igrd(i-1,j)

            if (igrd(i,j).eq.0.or.igrd(i,j).eq.kd) then
               ltic = ltic + 1
               if (j.eq.loopy) then 
                  write (n4,*) ltic,kd
                  kst = 1
               end if 
            else 
               write (n4,*) ltic,kd
               kst = j
               goto 20 
            end if 
         end do 
      end do         
c                                 write assemblage list
      write (n4,*) iasct

      do i = 1, iasct
         write (n4,*) iavar(1,i),iavar(2,i),iavar(3,i)
         write (n4,*) (idasls(j,i), j = 1, iavar(3,i))
      end do 

      end 

      subroutine wavgrd (output)
c----------------------------------------------------------------------
c wavgrd does constrained minimization on a 1 dimensional multilevel
c grid (ith column of a 2-d grid), lowest resolution is the default, 
c and phase boundaries are located at the highest level

c the phase assemblage at node(i,j) of the grid is identified by the 
c pointer igrd(i,j), igrd is a pointer to a phase assemblage.
c the compound grid data igrd is output to the "plot" file.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical output

      integer kinc,jinc(l8),i,j,k,klast,kmax,kd,idead,jtic,jmin,
     *        jmax,klev,klow

      integer ipot,jv,iv1,iv2,iv3,iv4,iv5
      common/ cst24 /ipot,jv(l2),iv1,iv2,iv3,iv4,iv5

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      integer igrd
      common/ cst311 /igrd(l7,l7)

      integer iap,ibulk
      common/ cst74 /iap(k2),ibulk

      integer jlow,jlev,loopx,loopy,jinc1
      common/ cst312 /jlow,jlev,loopx,loopy,jinc1
c-----------------------------------------------------------------------
c                               load arrays for lp solution
      call initlp
c                               jlow is the number of nodes
c                               at the lowest level, the number of
c                               nodes is
      klow = jlow - 1
c                               first level:
      loopy = klow * 2**(jlev-1) + 1 

      if (loopy.gt.l7) then
         call warn (92,v(iv1),loopy,'1dpath')
         klow = (l7 - 1)/2**(jlev-1)
         loopy = klow * 2**(jlev-1) + 1 
      end if        
c                               initialize igrd/iap (this is critical 
c                               for auto_refine).
      do j = 1, loopy
         igrd(1,j) = 0
         iap(j) = 0 
      end do        
c                                increments
      do j = 1, jlev
         jinc(j) = 2**(jlev-j)
      end do 

      jinc1 = jinc(1)
c                                initialize variables
      call setvar 
c                                initialize grid vars
      kinc = jinc(1)
      klev = 1

      i = 1

      j = 1
      jmax = 1 
      jmin = 0
      jtic = 0 

      do while (j.le.loopy) 

         if (igrd(i,j).eq.0) then
c                                call to incvr2 moved here (formerly
c                                at end of loop), april 5, 2006. JADC
            call setvr0 (j,1)
            call lpopt (i,j,idead,output)
            jtic = jtic + 1
         end if 

         kd = iap(igrd(i,j))
c                                jmax should always be the last known node
c                                at lowest level of resolution
         if (j.gt.jmax) then 
            jmax = j
            jmin = j - jinc(1) 
            kmax = kd
         end if 

         if (j.eq.1) then
 
            klast = kd
               
         else if (kd.ne.klast.and.klev.ne.jlev) then
c                                crossed a boundary and not
c                                at max resolution, increment resolution
            klev = klev + 1
            kinc = -jinc(klev)

         else if (kd.eq.klast.and.klev.ne.jlev) then 
c                                backed into field?
            if (kd.eq.kmax) then
c                                go to next low res node
               kinc = jinc(1)
               j = jmax
               klev = 1
               goto 10
            end if 
            
            klev = klev + 1
            kinc = jinc(klev)

         else if (klev.eq.jlev) then
c                                boundary has been located at max
c                                resolution
              
c                                jmin is the minimum value at which 
c                                subsequent searchs can go
            do k = j, jmax - 1
               if (igrd(i,k).eq.0) then 
                  jmin = k - 1
                  goto 20 
               end if 
            end do 
          
20          klast = iap(igrd(i,jmin))
c                                check if there is a grid point with
c                                this assemblage at a higher level
            do k = jmin, jmax - 1

               if (igrd(i,k).eq.0) cycle
               if (iap(igrd(i,k)).eq.klast) jmin = k

            end do 

            if (klast.eq.kmax.or.j.eq.jmax-1) then
c                                the grid can be filled in between j and jmax
c                                go on to the next low level grid point:
               klast = iap(igrd(i,jmax))
               kinc = jinc(1)
               klev = 1
               j = jmax
                  
            else 
c                                find the level to search for next boundary
c                                must be > level 1
               do k = 2, jlev
                  if (jmax-jmin.gt.jinc(k)) then
                     j = jmax
                     klev = k
                     kinc = -jinc(k)
                     goto 10 
                  end if
               end do 
            end if 
         end if 

10       j = j + kinc

         if (j.lt.jmin+1) then 
            j = j - kinc 
            klev = klev + 1
            kinc = -jinc(klev)
            goto 10 
         end if 

      end do 
c                                 output graphics data
      if (output) call outgrd (i,loopy,jinc(1))

      end 

      subroutine aminot (i,j,kinc,kinc2,kinc21)
c--------------------------------------------------------------------
c fill in nodes of a homogeneous cell, called only for lowest level
c grid, aminot1 is called for higher level cells. kinc is always
c even and > 1. modfied to fill from each corner if kinc > 2 in case
c compression is off. JADC sept '04. 

c kinc   - the current grid increment
c kinc2  - kinc/2 (integer division)
c kinc21 - kinc2 + 1
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, kinc, ii, jj, kinc2, kinc21

      integer igrd
      common/ cst311 /igrd(l7,l7)
c                                 the commented version fills
c                                 with just assemblage i,j this
c                                 would be fine if compression is
c                                 on, but is biased otherwise
c      do ii = i, i + kinc - 1
c         do jj = j, j + kinc - 1
c                                 the conditional prevents 
c                                 overwriting of identities if 
c                                 compression is off and searching
c                                 for true boundaries.
c            if (igrd(ii,jj).eq.0) igrd(ii,jj) = igrd(i,j)
c         end do
c      end do 
c                                 first fill from lower left (i,j)
      do ii = i, i + kinc2 
         do jj = j, j + kinc2
            if (igrd(ii,jj).eq.0) igrd(ii,jj) = igrd(i,j)
         end do
      end do
c                                 then from lower right (i+kinc,j)
      do ii = i + kinc21, i + kinc
         do jj = j, j + kinc2
            if (igrd(ii,jj).eq.0) igrd(ii,jj) = igrd(i+kinc,j)
         end do
      end do
c                                 then from upper left (i,j+kinc)
      do ii = i, i + kinc2
         do jj = j + kinc21, j + kinc
            if (igrd(ii,jj).eq.0) igrd(ii,jj) = igrd(i,j+kinc)
         end do
      end do
c                                 then from upper right (i+kinc,j+kinc)
      do ii = i + kinc21, i + kinc
         do jj = j + kinc2 + 1, j + kinc
            if (igrd(ii,jj).eq.0) igrd(ii,jj) = igrd(i+kinc,j+kinc)
         end do
      end do
     
      end  

      subroutine aminot1 (icent,jcent,i,j,kinc)
c--------------------------------------------------------------------
c fill in nodes of a homogeneous cell, called only for high level
c grids, aminot is called for lowest level cells. kinc is always
c even and > 1.  JADC sept '04. 

c kinc   - the current grid increment
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, ii, jcent, icent, jj, kinc

      integer igrd
      common/ cst311 /igrd(l7,l7)

      do ii = i, i + kinc
         do jj = j, j + kinc
c                                 the conditional prevents 
c                                 overwriting of identities if 
c                                 compression is off and searching
c                                 for true boundaries.
            if (igrd(ii,jj).eq.0) igrd(ii,jj) = igrd(icent,jcent)
         end do
      end do 
     
      end  

      subroutine fractr (output)  
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'
 
      character vname*8, xname*8, gname*10

      integer i,j,k

      logical there, warn, output
 
      common/ csta2 /xname(k5),vname(l2)

      double precision dcomp
      common/ frct2 /dcomp(k5)

      integer ifrct,ifr
      logical gone
      common/ frct1 /ifrct,ifr(k23),gone(k5)
 
      integer ipot,jv,iv
      common/ cst24 /ipot,jv(l2),iv(l2)

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      integer jbulk
      double precision cblk
      common/ cst300 /cblk(k5),jbulk

      character*5 cname
      common/ csta4 /cname(k5)

      integer kkp, np, ncpd, ntot
      double precision cp3, ctot3
      common/ cxt15 /cp3(k5,k5),ctot3(k5),kkp(k5),np,ncpd,ntot

      integer ipvt
      double precision a,b
      common/ cst301 /a(k5,k5),b(k5),ipvt(k5)

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)
c----------------------------------------------------------------------- 
c                                 fractionation effects:
      do i = 1, jbulk
         dcomp(i) = 0d0
      end do 

      do i = 1, ifrct

         there = .false.

         do j = 1, ntot

            if (kkp(j).eq.ifr(i).and.b(j).ge.0) then 

               there = .true.
c                                 the phase to be fractionated
c                                 is present
               do k = 1, jbulk 
                  dcomp(k) = dcomp(k) + b(j)*cp3(k,j)
               end do 
c                                 write to console
               write (*,1185) vname(iv(2)),v(iv(2)),
     *                        vname(iv(1)),v(iv(1))
               write (*,1190) b(j),gname(ifr(i)),
     *                        (b(j)*cp3(k,j),k=1,jbulk)
c                                 write to file
               if (output) write (n0+i,1200) v(iv(1)),b(j),
     *                     gname(ifr(i)),(b(j)*cp3(k,j),k=1,jbulk)
            end if
         end do 

         if (.not.there) then
            write (*,1185) vname(iv(2)),v(iv(2)),vname(iv(1)),v(iv(1))
            write (*,1210) gname(ifr(i))
         end if 
 
      end do   

      warn = .false.     

      do i = 1, jbulk
         cblk(i) = cblk(i) - dcomp(i)

         if (cblk(i).le.nopt(11)) then 

            if (.not.gone(i)) then
               warn = .true.
               gone(i) = .true. 
               write (*,1000) cname(i)
            end if 

            if (cblk(i).lt.nopt(11)) then 
               cblk(i) = 0d0
            end if 
         end if 

      end do

      if (output.and.ifrct.gt.0) 
     *   write (n0,1220) v(iv(1)),(cblk(k),k=1,jbulk)
 
      if (warn) then 
         do i = 1, jbulk 
            if (.not.gone(i)) write (*,1010) cname(i),cblk(i)
         end do
      end if  

1000  format (3(/,5('*** WARNING ***')),
     *     //,'Fractionation has completely'
     *       ,' depleted ',a,' from the bulk composition.',//,
     *        'This ',
     *        'depletion may destabilize the minimization algorithm ',
     *        'if excessive',/,'minimization errors occur, restart the',
     *        'fractionation calculation at the',/,'current point on ',
     *        'the fractionation path with the current molar bulk ', 
     *        'composition:',/)
1010  format (4x,a,2x,g12.6)       
1185  format (/,'At ',a,'=',g12.6,' and ',a,'=',g12.6)
1190  format ('fractionating ',g12.6,' moles of ',a,'; changes bulk by:'
     *        ,/,15(1x,g12.6))
1200  format (g12.6,1x,g12.6,1x,a,15(1x,g12.6))
1210  format (a,' is not stable.')
1220  format (17(g12.6,1x))
      end 

      subroutine frname 
c----------------------------------------------------------------------
c frname - prompt for names of phases to be fractionated

c          ifrct - number of phases to be fractionated
c          ifr(ifract) - 0 if cpd, else ikp (solution pointer)
c          jfr(ifract) - cpd pointer
c----------------------------------------------------------------------

      implicit none
 
      include 'perplex_parameters.h'

      logical first 

      integer iam, jfrct, i

      character*10 phase(k23), fname*14, y*1


      integer ifrct,ifr
      logical gone
      common/ frct1 /ifrct,ifr(k23),gone(k5)
 

      save first,phase

      data first/.true./
c----------------------------------------------------------------------

      if (first) then 

         first = .false.

         write (*,1030) 
         read (*,1020) y
         if (y.ne.'y'.and.y.ne.'Y') goto 99 

         ifrct = 1

         do 

            write (*,1040)
            read (*,1020) phase(ifrct)

            if (phase(ifrct).eq.' ') exit

            call matchj (phase(ifrct),ifr(ifrct))

            if (ifr(ifrct).eq.0) then
               write (*,1100) phase(ifrct)
               cycle 
            end if

            ifrct = ifrct + 1

            if (ifrct.gt.k23) call error (1,0d0,ifrct,'k23')

         end do 

         ifrct = ifrct - 1

      else 
c                                 must be in autorefine cycle, make
c                                 new phase list from old list:
         jfrct = ifrct
         ifrct = 0  
         
         do i = 1, jfrct 

            call matchj (phase(i),iam)

            if (iam.eq.0) cycle 
       
            ifrct = ifrct + 1
            ifr(ifrct) = iam 

         end do 

      end if 
c                                 initialize "gone" flags
      do i = 1, k5
         gone(i) = .false.
      end do 
c                                 open fractionation files
      open (n0,file='fractionated_bulk.dat',status='unknown')
      write (*,1050)

      do i = 1, ifrct 

         write (fname,1000) phase(i),'.dat'
         call unblnk (fname)

         write (*,1010) phase(i), fname

         open (n0+ifrct,file=fname,status='unknown')

      end do 

1000  format (a,a)
1010  format (/,'The fractionated amount and composition of ',a,/,
     *          'will be written to file: ',a,/)
1020  format (a)
1030  format ('Fractionate phases (y/n)?')
1040  format (/,'Enter the name of a phase to be fractionated',
     *        /,'(left justified, <cr> to finish): ')
1050  format (/,'The fractionated bulk composition will be ',
     *          'written to file: fractionated_bulk.dat',/)
1100  format (/,'No such entity as ',a,', try again: ')

99    end 