 
c Copyright (c) 1990 by James A. D. Connolly, Institute for Mineralogy and
c Petrography, Swiss Federal Insitute of Technology, CH-8092 Zurich,
c SWITZERLAND. All rights reserved.
 
c PSVDRAW - a program to generate postscript code for x-y plots and chemographic
c phase diagrams.  The input data format is consistent with the files generated
c by the programs VERTEX and FRENDLY.
 
c Please do not distribute any part of this source.
 
      PROGRAM PSPTS

      implicit none

      include 'perplex_parameters.h'

      integer ier, ier99
 
      character*100 fname

      integer  iop0 
      common / basic /iop0

      logical  debug
      common / debugblk / debug
c----------------------------------------------------------------------
c   Look for the "debug_yes" file to turn on debugging messages

      open (97,iostat=ier99,file='debug_yes',status='old')
      if (ier99.eq.0) then
          debug=.TRUE.
      else
          debug=.FALSE.
      end if
      close(97);

      if (debug) PRINT *,'IN PSPTS'

c                                 set iop0 to 1 to allow 
c                                 drafting ptompts
      iop0 = 1

      do 

         write (*,1020) 
         read (*,1000) fname
 
         call getfil (fname,n4,ier)

         if (ier.eq.0) exit 

      end do
c                                 read plot option file, set
c                                 default transformation
      call rdopt 
c                                 open output file 
      call psopen (fname)

      call psxypl 
 
      call psclos
 
      close (n4)
 
1000  format (a)
1020  format (/,'Enter the POINT plot file name: ',/)
 
      end

      subroutine psxypl
c--------------------------------------------------------------------- 
c psxypl - subroutine to output x-y plot.
 
      implicit none
 
      include 'perplex_parameters.h'

      integer iop1, jop0, ier, isym

      double precision x, y
 
      integer jvar
      double precision var,dvr,vmn,vmx
      common/ cxt18 /var(l3),dvr(2),vmn(l3),vmx(l3),jvar

      integer  iop0 
      common / basic /iop0

      logical  debug
      common / debugblk / debug

      character vnm*8
      common/ cxt18a /vnm(l3)   
c--------------------------------------------------------------------- 
      if (debug) PRINT *,'IN PSXYPL'

      jvar = 2
      vmn(1) = 1d30
      vmx(1) = -1d30
      vmn(2) = 1d30
      vmx(2) = -1d30
      vnm(1) = 'x axis'
      vnm(2) = 'y axis'
c                                 read the data to get the range  
      do                             
         read (n4,*,iostat=ier) isym, x, y
         if (ier.ne.0) exit 
         if (x.lt.vmn(1)) vmn(1) = x
         if (x.gt.vmx(1)) vmx(1) = x
         if (y.gt.vmx(2)) vmx(2) = y
         if (y.lt.vmn(2)) vmn(2) = y
      end do 
c                                 get some options and
c                                 set up transformations
      call psaxop (1,jop0,iop1)
 
      call psipts 
 
      call psaxes (jop0)

      end

      subroutine psipts 
c---------------------------------------------------------------
c psipts - subprogram to x-y points.

      implicit none

      include 'perplex_parameters.h'

      integer ier, isym

      double precision r,x,y

      double precision xmin,xmax,ymin,ymax,dcx,dcy,xlen,ylen
      common/ wsize /xmin,xmax,ymin,ymax,dcx,dcy,xlen,ylen 
c---------------------------------------------------------------  
      rewind (n4)

      do 

         read (n4,*,iostat=ier) isym,x,y
         if (ier.ne.0) exit

         r = 0.78d0

 
         if (isym.lt.4) then
            if (isym.eq.0) then

               call pselip (x,y,r*dcx,r*dcy,1d0,0d0,7)

            else if (isym.eq.1) then
 
               call psrect (x-r*dcx,x+r*dcx,y-r*dcy,y+r*dcy,1d0,0d0,7)

            else if (isym.eq.2) then

               call pselip (x,y,r*dcx,r*dcy,1d0,0d0,1)

            else if (isym.eq.3) then

               call psrect (x-r*dcx,x+r*dcx,y-r*dcy,y+r*dcy,1d0,0d0,1)

            end if 

         else 
 
            r = 0.38d0

            if (isym.eq.4) then

               call pselip (x,y,r*dcx,r*dcy,1d0,0d0,7)

            else if (isym.eq.5) then
 
               call psrect (x-r*dcx,x+r*dcx,y-r*dcy,y+r*dcy,1d0,0d0,7)

            else if (isym.eq.6) then

               call pselip (x,y,r*dcx,r*dcy,1d0,0d0,1)

            else if (isym.eq.7) then

               call psrect (x-r*dcx,x+r*dcx,y-r*dcy,y+r*dcy,1d0,0d0,1)

            end if 

         end if 

      end do 
 
      end
