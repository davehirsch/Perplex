      PROGRAM PSCONTOUR 

      implicit none

      include 'perplex_parameters.h'

      integer ier, ier99
 
      character*100 fname, yes*1

      integer  iop0 
      common / basic /iop0

      logical  debug
      common / debugblk / debug
c----------------------------------------------------------------------
c   Look for the "debug_yes" file to turn on debugging messages

      open (97,iostat=ier99,file='debug_yes',status='old')
      if (ier99.eq.0) then
          debug=.TRUE.
          close (99)
      else
          debug=.FALSE.
      end if
      close(97);

      if (debug) PRINT *,'IN PSCONTOR'

      do 
c                                 get plot input file 
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
c                                 allow drafting options prompt
      iop0 = 0
      write (*,1030) 
      read (*,1000) yes
      if (yes.eq.'y'.or.yes.eq.'Y') iop0 = 1

      call psxypl
 
      call psclos
 
      close (n4)
 
1000  format (a)
1020  format (/,'Enter the CONTOUR plot file name:')
1030  format (/,'Modify the default plot (y/n)?')

      end

c---------------------------------------------------------------------
      subroutine psxypl 
 
c psxypl - subroutine to output x-y plot.

      implicit none

      include 'perplex_parameters.h'

      character y*1, fname*162

      integer nx,ny,i,j,iox,ioy,jmn,imn,imx,jop0,ncon,jmx,iop1

      double precision dx,dy,xpmn,xpmx,cmin,cmax,dcon,ypmx,ypmn,
     *                 z0min,z0max
 
      parameter (nx=1000,ny=1000)

      integer jvar
      double precision var,dvr,vmn,vmx
      common/ cxt18 /var(l3),dvr(2),vmn(l3),vmx(l3),jvar

      character vnm*8
      common/ cxt18a /vnm(l3)  

      integer  iop0 
      common / basic /iop0

      logical  debug
      common / debugblk / debug

      integer ix,iy
      double precision z,zt 
      common/ dim   /z(nx,ny),zt(nx,ny),ix,iy

      double precision zmin,zmax
      common/ stuff /zmax,zmin

      double precision xmin,xmax,ymin,ymax,dcx,dcy,xlen,ylen
      common/ wsize /xmin,xmax,ymin,ymax,dcx,dcy,xlen,ylen
c----------------------------------------------------------------------
      if (debug) PRINT *,'IN PSXYPL'
      read (n4,'(a)') fname
      if (debug) PRINT *,'PSXYPL.  fname=', fname
      read (n4,*) ix,iy,xmin,ymin,dx,dy
      if (debug) PRINT *,'PSXYPL.  ix=', ix
      if (debug) PRINT *,'PSXYPL.  iy=', iy
      if (debug) PRINT *,'PSXYPL.  xmin=', xmin
      if (debug) PRINT *,'PSXYPL.  ymin=', ymin
      if (debug) PRINT *,'PSXYPL.  dx=', dx
      if (debug) PRINT *,'PSXYPL.  dy=', dy
      read (n4,1010) (vnm(i),i=1,2)
      if (debug) PRINT *,'PSXYPL.  vnm=', vnm(1), vnm(2)
      if (ix.gt.nx) call error (1,dx,nx,'NX, PSXYPL')
      if (iy.gt.ny) call error (1,dx,ny,'NY, PSXYPL')
      read (n4,*) ((z(i,j), i = 1, ix), j = 1, iy)

      if (iop0.eq.1) then 
         write (*,1050) 
         read (*,'(a)') y
         if (y.eq.'y') then 
            do j = 1, ny
               do i = 1, nx
                  if (z(i,j).ne.0d0) z(i,j) = dlog10(dabs(z(i,j)))
               end do 
            end do  
         end if 
      end if 

      ypmn = ymin
      xpmn = xmin
      ypmx = ymin + (iy-1)*dy
      xpmx = xmin + (ix-1)*dx    
      ymax = ypmx
      xmax = xpmx     

      write (*,1060) 
      read (*,'(a)') y
      if (y.eq.'y') then 
         write (*,1070) ypmx,ypmn,xpmx,xpmn
         read (*,*) ypmx,ypmn,xpmx,xpmn
      end if 

      iox = ix
      ioy = iy

      if (y.eq.'y') then 
         jmn = int(ypmn/dy) + 1
         jmx = int(ypmx/dy) + 1
         imn = int(xpmn/dx) + 1
         imx = int(xpmx/dx) + 1
         iy = (jmx-jmn+1)
         ix = (imx-imn+1)
         ymax = ypmn + (iy-1)*dy
         ymin = ypmn 
         xmin = xpmn
         xmax = xpmn + (ix-1)*dx
      end if 
c                                      reload mini matrix:
      if (y.eq.'y') then
         do i = 1, ix
            do j = 1, iy
               z(i,j) = z(i+imn-1,j+jmn-1)
            end do
         end do
      end if 

      vmn(1) = xmin
      vmn(2) = ymin
      vmx(1) = xmax
      vmx(2) = ymax
c                                 get some options and
c                                 set up transformations
      call psaxop (1,jop0,iop1)
        
      zmin = 1d9
      zmax = -1d9
      z0min = 1d30
      z0max = -1d30
c                                      set up contour intervals                                      
      do i = 1, ix
         do j = 1, iy 
            if (z(i,j).lt.zmin) zmin = z(i,j)
            if (z(i,j).gt.zmax) zmax = z(i,j)
            if (z(i,j).lt.z0min.and.z(i,j).ne.0d0) z0min = z(i,j)
            if (z(i,j).gt.z0max.and.z(i,j).ne.0d0) z0max = z(i,j)
         end do 
      end do 
c                                      set up contour intervals
      write (*,1020) zmin, zmax, z0min, z0max 
      read (*,'(a)') y

      if (y.eq.'y') then
         write (*,1030) 
         read (*,*) cmin, cmax, dcon
         ncon = int((cmax-cmin)/dcon) + 1
      else 
         dcon = (zmax-zmin)/11.
         cmax = zmax - 0.5d0 * dcon
         cmin = zmin + 0.5d0 * dcon
         ncon = 11
      end if 

      call pscontor (cmin,ncon,dcon)
 
      call psaxes (jop0)
 
1010  format (10a8)
1020  format ('Contoured variable range:',g14.6,'->',g14.6,/,
     *        'Range excluding zero values:',g14.6,'->',g14.6,/,
     *        'Modify default contour interval (y/n)?')
1030  format ('Enter min, max and interval for contours:')
1050  format ('Contour log10 of the z-value (y/n)?')
1060  format (/,'Reset plot limits (y/n)?')
1070  format (/,'Old values were: ',4(g12.4),/,'Enter new values:')

      end

