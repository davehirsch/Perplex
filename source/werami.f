
      PROGRAM WERAMI
c----------------------------------------------------------------------
c this version modified by M. Caddick 3/11/06 to permit switching 
c between modes.
c---------------------------------------------------------------------- 
      implicit real (a-g,o-z),integer (h-n)

      include 'perplex_parameters.h'

      real xy(2), xyp(2,2), prp, dxy(2), coef(0:10),
     *     tmin(2), tmax(2)

      real prmin, prmax, xx(5*l5), yy(5*l5)

      integer icx, jcx, lop
 
      character*100 n4name, fname, vname*8, yes*1, dname

      save icoors

      common/ vs    /vmin(l3),vmax(l3),iv/ vsa /vname(l3)
      common/ basic /n4,ifont,icopt,jwidth,iop0
     *      / cst84 /icx,jcx,lop
 
      n4 = 24
c                                 read options if present:
      call readop 
c                                 get plot file name
      write (*,1250) 
      read (*,1050) fname
c                                 read header info from 
c                                 plot file
      open (n4,err=91,file=fname,status='old')
      call plinp (tmin,tmax,pmin,pmax)
      close (n4)
c                                 read assemblage data
c                                 from b-plot file 
      call bplinp (fname,kpoly)
c                                 read polygons from p-plot file:
      if (icopt.eq.6) call pplinp (fname)

20    write (*,1020) 
      read (*,*,err=20) imode

      if (imode.eq.5) then 
         goto 99
         
      else if (imode.eq.1) then 

         n4name = fname 
         call inblnk (n4name,'r')
         write (*,1030) n4name
         open (n4,file=n4name) 

10       write (*,1000) (vname(i), i = 1, 2)
         read (*,*,) x, y

         if (x+y.eq.0e0) goto 20

         xy(1) = x
         xy(2) = y

         do i = 1, 2
            if (vmin(i).lt.vmax(i)) then
               if (xy(i).lt.vmin(i).or.xy(i).gt.vmax(i)) then  
                  write (*,1010) vname(i),vmin(i),vmax(i)
                  goto 10 
               end if 
            else 
               if (xy(i).lt.vmax(i).or.xy(i).gt.vmin(i)) then  
                  write (*,1010) vname(i),vmin(i),vmax(i)
                  goto 10 
               end if 
            end if 
         end do 
 
         call amiin (x,y,iam)
         if (iam.ne.0) then 
            call calpr1 (x,y,iam,6)
            call calpr1 (x,y,iam,n4)
         else 
            write (*,1070)
         end if 
         
         goto 10 

      else if (imode.eq.2) then 

         nplot = 0

25       nplot = nplot + 1

         call chsprp (lop,icx,jcx)
c                                         set up coordinates etc
         if (nplot.eq.1) then 
            write(*,1040)
            read (*,1050) yes 
            if (yes.eq.'y'.or.yes.eq.'Y') then 
               do i = 1, 2
30                write (*,1060) vname(i),vmin(i),vmax(i)
                  read (*,*,err=30) tmin(i),tmax(i)
               end do 
            else 
               tmin(1) = vmin(1)
               tmin(2) = vmin(2)
               tmax(1) = vmax(1)
               tmax(2) = vmax(2)
            end if 

            write (*,1080) 
            read (*,*) nx,ny
         
            do i = 1, 2
               tmin(i) = tmin(i) + (tmax(i)-tmin(i))*1e-6
               tmax(i) = tmax(i) - (tmax(i)-tmin(i))*1e-6
            end do 

            dx = (tmax(1)-tmin(1))/float(nx-1)
            dy = (tmax(2)-tmin(2))/float(ny-1)

            n4name = fname 
            call inblnk (n4name,'c')
         else 
            call inblnk (n4name,'c')
         end if  

         write (*,1110) n4name
         open (n4,file=n4name)

         write (n4,1050) fname
c        write (n4,*) nx,ny,vmin(1),vmin(2),dx,dy
c                                              write weird format
c                                              for matlab scripts based
c                                              on windows version:
         write (n4,1260) nx,ny,tmin(1),tmin(2),dx
         write (n4,1270) dy         

         prmin = 1e16
         prmax = -1e16

         write (n4,1100) (vname(i),i=1,2)
         do j = 1, ny 
            y = tmin(2)+dy*float(j-1)
            do i = 1, nx 
               call polprp (prp,tmin(1)+dx*float(i-1),y)
               if (prp.gt.0.and.prp.gt.prmax) prmax = prp
               if (prp.gt.0.and.prp.lt.prmin) prmin = prp
               write (n4,*) prp
            end do 
         end do 

         close (n4)

         write (*,1280) prmin, prmax

         write (*,1230)
         read (*,1050) yes

         if (yes.eq.'y'.or.yes.eq.'Y') goto 25

      else if (imode.eq.3) then 

         call chsprp (lop,icx,jcx)

         n4name = fname 
         call inblnk (n4name,'c')
         write (*,1130) n4name
         open (n4,file=n4name) 

65       write (*,1200) 
         read (*,1050) yes 

         icurve = 0 
         idxy = 0 
         ivi = 1
         ivd = 2

         if (yes.eq.'y'.or.yes.eq.'Y') then 
            icurve = 1
            write (*,1210) vname(ivd),vname(ivi)
            read (*,*) iord
            do i = 0, iord
               write (*,1220) i
               read (*,*) coef(i)
            end do
            write (*,1340) vname(ivd),(coef(i),vname(ivi),i,i=0,iord)

            write (*,1320)
            read (*,1050) yes 
            if (yes.eq.'y'.or.yes.eq.'Y') goto 65

            dxy(ivi) = vmax(ivi)-vmin(ivi)
            dxy(ivd) = vmax(ivd)-vmin(ivd)
            xyp(ivi,1) = vmin(ivi)
            xyp(ivd,1) = vmin(ivd)

         else 

45          do i = 1, 2
40             write (*,1140) i,vname(1),vname(2)
               read (*,*,err=40) xyp(1,i),xyp(2,i)
               do j = 1, 2
                  if (vmin(j).lt.vmax(j)) then 
                     if (xyp(j,i).lt.vmin(j).or.xyp(j,i).gt.vmax(j)) 
     *               then  
                        write (*,1010) vname(j),vmin(j),vmax(j)
                        goto 40
                     end if 
                  else
                     if (xyp(j,i).lt.vmax(j).or.xyp(j,i).gt.vmin(j)) 
     *               then  
                        write (*,1010) vname(j),vmin(j),vmax(j)
                        goto 40
                     end if 
                  end if 
               end do 
            end do 

            do j = 1, 2
               dxy(j) = xyp(j,2) - xyp(j,1)
               if (dxy(j).eq.0d0) idxy = j
            end do 

            if (dxy(1).eq.0d0.and.dxy(2).eq.0d0) then
               write (*,*) 
     *               'initial and final coordinates cannot be identical'
               goto 45
            end if   
         end if 

         write (*,1150) 
         read (*,*) ipts

         if (idxy.eq.0.and.icurve.eq.0) then 
            write (*,1160) (i,vname(i),i= 1, 2)
            read (*,*) ivi
            ivd = 2
            if (ivi.eq.2) ivd = 1 

            s = dxy(ivd)/dxy(ivi)
         else if (icurve.eq.0) then
            ivi = 1
            ivd = idxy 
            if (idxy.eq.1) ivi = 2
            s = 0d0 
         end if 

         d = dxy(ivi)/float(ipts - 1)

         write (*,1180) vname(ivi),vname(ivd)

         nprop = 0

80       jpts = 0 

         do i = 1, ipts
            xy(ivi) = xyp(ivi,1) + float(i-1)*d
            if (icurve.eq.0) then 
               xy(ivd) = xyp(ivd,1) + s*(xy(ivi)-xyp(ivi,1))
               jpts = 1
            else 
               xy(ivd) = 0d0
               do j = 0, iord
                  xy(ivd) = xy(ivd) + coef(j)*xy(ivi)**j
               end do
               if (xy(ivd).le.vmax(ivd).and.xy(ivd).ge.vmin(ivd)) then 
                  jpts = jpts + 1
               else 
                  cycle 
               end if 
            end if 

            call polprp(p,xy(1),xy(2))
            write (n4,1170) nprop,xy(ivi),p,xy(ivd)
            write (*,1190) xy(ivi),xy(ivd),p
         end do 

         if (jpts.eq.0) then 
            write (*,1330) 
            goto 65
         end if 

         write (*,1230)
         read (*,1050) yes
         nprop = nprop + 1
         if (yes.eq.'y'.or.yes.eq.'Y') then 
            call chsprp (lop,icx,jcx)
            goto 80
         end if 

      else if (imode.eq.4) then 

         nprop = 0

         call chsprp (lop,icx,jcx)

         n4name = fname 
         call inblnk (n4name,'c')
         write (*,1130) n4name
         open (n4,file=n4name) 
          
         write (*,1290) 
         read (*,*) jmode
          
         write (*,1300) 
         read (*,'(a)') dname

         if (jmode.eq.1) then 
c                                  points from a polynomial
            open (33,file=dname,err=120)
            read (33,*) pmin,pmax,ixy
110         read (33,*,end=99) tmin1,tmax1,dt1,a1,b1,c1,d1

            if (dt1.lt.0d0) then 
               x0 = tmax1
               x1 = tmin1
            else 
               x0 = tmin1
               x1 = tmax1
            end if 

90          x = x0

            do i = 1, 1000
               y = a1 + b1*x + c1*x**2 + d1*x**3
               if (y.le.pmax.and.y.ge.pmin) then 
                  if ((dt1.gt.0d0.and.x.le.x1).or.
     *                (dt1.lt.0d0.and.x.ge.x1)) then
                     if (ixy.eq.0) then 
                        call polprp(p,x,y)
                     else
                        call polprp(p,y,x)
                     end if 
                     write (n4,1170) nprop,x,p,y
                     write (*,1190) x,y,p
                  else 
                     goto 100
                  end if 
               end if 
               x = x + dt1
            end do  

100         write (*,1230)
            read (*,1050) yes
            nprop = nprop + 1
            if (yes.eq.'y'.or.yes.eq.'Y') then 
               call chsprp (lop,icx,jcx)
               goto 90
            end if 
c                                   more than one geotherm may
c                                   be entered in a file
            write (*,1320)
            read (*,1050) yes
            if (yes.eq.'y'.or.yes.eq.'Y') goto 110  
            goto 20
         else 
c                                   points from a data file:
            open (33,file=dname,err=120)
            icoors = 1
            do 
               read (33,*,end=290) xx(icoors),yy(icoors)
               icoors = icoors + 1
               if (icoors.gt.5*l5) then 
                  write (*,*) '**error** too many points, ',
     *                        'increase parameter l5.'
                  goto 99
               end if 
            end do 
               
290         icoors = icoors - 1
            write (*,1310) icoors
            read (*,*) inc
            
300         do i = 1, icoors, inc
               call polprp(p,xx(i),yy(i))
               write (n4,1170) nprop,yy(i),p,xx(i)
               write (*,1190) xx(i),yy(i),p
            end do  

            write (*,1230)
            read (*,1050) yes
            nprop = nprop + 1
            if (yes.eq.'y'.or.yes.eq.'Y') then 
               call chsprp (lop,icx,jcx)
               goto 300
            end if 
         end if 
      else 
  
         goto 20 

      end if 

      goto 20
120   write (*,*) 'No such data file as: ',dname

      goto 99

91    write (*,1240) fname

1000  format (/,'Enter ',a,' and ',a,' (zeroes to quit):')
1010  format (/,'The plot file range for ',a,' is ',g12.4,' - ',g12.4,
     *        /,'Try again:',/)
1020  format (/,'Select operational mode:',/,
     *        4x,'1 - compute properties at specified conditions',/,
     *        4x,'2 - create a property grid (plot with pscontor)',/,
     *        4x,'3 - compute properties along a curve or line',
     *        ' (plot with pspts)',/,
     *        4x,'4 - as in 3, but input from file',/,
     *        4x,'5 - EXIT',/)
1030  format (/,'Console output will be echoed in file: ',a,/)
1040  format (/,'Change default variable range (y/n)?')
1050  format (a)
1060  format (/,'Current limits on ',a,' are: ',g12.6,'->',g12.6,/,
     *        'Enter new values:')
1070  format (/,'No polygon bounds the specified conditions.',/,
     *        'Data missing, or you have specified conditions on an',
     *        ' edge of a polygon.',/)
1080  format (/,'Enter number of nodes in the x and y directions:')
1100  format (2a8)
1110  format (/,'Writing grid data to file: ',a,/)
1130  format (/,'Writing profile to file: ',a,/)
1140  format (/,'Enter endpoint ',i1,' (',a,'-',a,') coordinates:')
1150  format (/,'How many points along the profile?')
1160  format (/,'Select independent variable: ',2(/,1x,i1,' - ',a))
1170  format (1x,i1,3(1x,g12.6))
1180  format (/,3x,a8,5x,a8,2x,'  Property   ',/)
1190  format (3(1x,g12.6))
1200  format (/,'Construct a non-linear profile (y/n)?')
1210  format (/,'Profile must be described by the function',/,a,
     *        ' = Sum ( c(i) * ',a,' ^i, i = 0..n)',/,'Enter n (<10)')
1220  format (/,'Enter c(',i2,')')
1230  format (/,'Evaluate additional properties (y/n)?',/)
1240  format (/,'No such file as: ',a,/)
1250  format (/,'Enter the VERTEX plot file name:')
1260  format (i5,1x,i5,4(1x,g14.6))
1270  format (g14.6)
1280  format (/,'Range of >0 data is: ',g14.6,' -> ',g14.6,/)
1290  format (/,'Path will be described by:',/,
     *          '   1 - a file containing a polynomial function',/,
     *          '   2 - a file containing a list of x-y points',/,
     *          'Enter 1 or 2 (default=2):'/) 
1300  format (/,'Enter the file name:',/)
1310  format (/,'File contains ',i5,' points',/,
     *          'every nth plot will be plotted, enter n:',/)
1320  format (/,'Change the geotherm (Y/N)?')
1340  format (/,'Your polynomial is:',/,
     *        a,'=',5('+(',g12.6,')','*',a,'^',i1))
1330  format (/,'Your polynomial does not yield conditions within',
     *          'computational coordinate frame.',/,'Try again.',/)

99    end
