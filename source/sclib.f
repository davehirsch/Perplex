      subroutine setau1 (vertex,output)
c----------------------------------------------------------------------
c setau1 sets autorefine dependent parameters. vertex is true if vertex
c is the calling program. output is set to false if autorefine mode is 
c not auto (i.e., iopt(6) = 2) or it is auto and in the second cycle.
c----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      logical vertex, output
 
      character*8 y*1

      character*10 badnam(h9)

      integer ibad2,ibad1,igood,i,j,ierr

      character n10nam*100,n11nam*100,n12nam*100

      character*100 n1name,cfname
      common/ cst228 /n1name,cfname
c                                 solution model counter
      integer isoct
      common/ cst79 /isoct
c                                 solution model names
      character*10 fname
      common/ csta7 /fname(h9)

      integer isec,icopt,ifull,imsg,io3p
      common/ cst103 /isec,icopt,ifull,imsg,io3p

      integer iopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10)

      logical refine
      common/ cxt26 /refine
c-----------------------------------------------------------------------
      refine = .false.
c                                 only use autorefine if solutions
c                                 are present and it is requested.
      if (isoct.ne.0) then 

         call mertxt (n10nam,'auto_refine_',n1name)
         call unblnk (n10nam)
         open (n10, file = n10nam, iostat = ierr, status = 'old')

         call mertxt (n12nam,n10nam,'_true_or_false')
         call unblnk (n12nam)

         if (vertex) then

            open (n8, file = n12nam, status = 'unknown')
c                                 user friendly text version 
c                                 to unit n11
            call mertxt (n11nam,n10nam,'.txt')
            call unblnk (n11nam)
            open (n11, file = n11nam, status = 'unknown')

            ibad1 = 0 
            igood = 0 

            if (ierr.ne.0) then 
c                                 no auto_refine data
               write (*,1020) n10nam
               open (n10, file = n10nam, status = 'unknown')

            else 
                   
               read (n10,*,iostat=ierr) ibad1, ibad2, igood
               if (ibad1.gt.0) read (n10,'(a)') (badnam(i),i=1,ibad1)

               if (iopt(6).ne.2.or.output) write (*,1030) n10nam

               if (iopt(6).eq.1) then 
c                                 manual mode, allow reinitialization
c                                 or suppression.
                  write (*,1060) 
                  read (*,'(a)') y

                  if (y.eq.'y'.or.y.eq.'Y') then

                     iopt(6) = 0
                     igood = 0

                  else 

                     refine = .true.  

                  end if

                  output = .true.
 
               else if (output) then  
c                                 second cycle of automated mode
                  refine = .true.

               end if  

               write (n8,*) refine

            end if 

         else 
c                                 werami/pssect if refine, get the 
c                                 solution models to be rejected
            open (n8, file = n12nam, iostat=ierr, status = 'old')
        
            if (ierr.eq.0) then 
c                                 write a flag to indicate if auto-refine
c                                 has been used, this is necessary so that other
c                                 perplex programs know whether to reject the
c                                 badnam phases:
               read (n8,*,iostat=ierr) refine
c                                 read phases to be rejected if in auto-refine
               if (refine) then 
                  read (n10,*,iostat=ierr) ibad1, ibad2, igood
                  if (ibad1.gt.0) read (n10,'(a)') (badnam(i),i=1,ibad1)
               end if 

            end if 

         end if 

      end if 

      close (n8)
c                                 just to be sure
      if (iopt(6).eq.0) refine = .false.

      if (refine) then 
c                                 reject solution models that were 
c                                 not found to be stable and set parameters 
c                                 that depend on refinement
         ibad2 = 0 

         do 50 i = 1, isoct

            do j = 1, ibad1
               if (fname(i).eq.badnam(j)) then
                  if (vertex) write (*,1070) fname(i)
                  goto 50
               end if 
            end do 

            ibad2 = ibad2 + 1
            fname(ibad2) = fname(i)

50       continue 

         isoct = ibad2 

         write (*,'(/)')

      end if

      if (iopt(6).eq.2.and..not.refine) then 
         output = .false.
      else
         output = .true.
      end if 
c                                 check that solvus tolerance is not too small:
      if (nopt(8).gt.0d0.and.isoct.gt.0) then 

         if (icopt.eq.1.and.nopt(13)/dfloat(iopt(9)).ge.nopt(8)) then
c                                 non-adaptive minimization Schreinemakers diagram
            nopt(8) = 1.2d0 * nopt(13)/dfloat(iopt(9))
            write (*,1080) nopt(8)

         else if (icopt.le.3.and.nopt(13)/dfloat(iopt(8)).ge.nopt(8))
     *                                                           then 
c                                 other non-adaptive minimization calculations
            nopt(8) = 1.2d0 * nopt(13)/dfloat(iopt(8))
            write (*,1080) nopt(8) 

         else if (icopt.ge.5.and.nopt(10)/1d1.ge.nopt(8)) then 
c                                 adaptive minimization
            nopt(8) = 1.2d0 * nopt(10)/1d1
            write (*,1080) nopt(8) 

         end if 
      end if 

1020  format (/,'Writing data for auto-refinement to file: ',a,/)
1030  format (/,'Reading data for auto-refinement from file: ',a,/)
1060  format ('Suppress or reinitialize auto-refinement (y/n)?')
1070  format ('Eliminating solution model: ',a,' in auto-refinement.')
1080  format (/,'WARNING: the compositional resolution specified for ',
     *          'the auto-refine stage is < the solvus',/,'tolerance, ',
     *          'the solvus tolerance will be increased to: ',f5.3,/,
     *          'To maintain the original solvus tolerance adjust the ',
     *          'compositional resolution specified in',/,
     *          'perplex_option.dat',/)
      end 

      subroutine setau2 (output)
c----------------------------------------------------------------------
c setau2 sets/resets autorefine parameters after the solution models have
c been read. setau1 must be called first.

c output is set to true if autorefine mode is auto (i.e., iopt(6) = 2) 
c but no solutions are present (isoct = 0). 
c----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      logical output

      integer i,index
c                                 solution model counter
      integer isoct
      common/ cst79 /isoct

      integer isec,icopt,ifull,imsg,io3p
      common/ cst103 /isec,icopt,ifull,imsg,io3p

      integer ivarrx,ivarip,isudo,ivar
      common/ cst62 /ivarrx(k2),ivarip(k2),isudo,ivar

      integer jlow,jlev,loopx,loopy,jinc
      common/ cst312 /jlow,jlev,loopx,loopy,jinc

      integer grid
      double precision rid 
      common/ cst327 /grid(5,2),rid(2)

      integer ipot,jv,iv
      common/ cst24 /ipot,jv(l2),iv(l2)

      double precision vmax,vmin,dv
      common/ cst9  /vmax(l2),vmin(l2),dv(l2)

      logical oned
      common/ cst82 /oned

      logical refine
      common/ cxt26 /refine
c-----------------------------------------------------------------------
      if (isoct.eq.0) then 
     
         index = 2
         output = .true.

      else if (.not.output) then

         index = 1

      else 

          if (refine) then

             index = 2

          else 

             index = 1

          end if 

      end if 
c                                 set auto-refine dependent parameters
      if (icopt.eq.5) then 
c                                 gridded minimization
         if (oned) then 
            jlow = grid(4,index)
            loopx = 1
         else 
            jlow = grid(1,index)
            loopx = grid(2,index) 
         end if

         jlev = grid(3,index) 
          
      else if (icopt.gt.5) then 
c                                 1d/2d phase fractionation
         jlow = grid(4,index)

      else if (icopt.eq.1) then 
c                                 schreinemakers diagrams

c                                 max variance of curves to be traced
          isudo = grid(5,index)
c                                 default variable tracing increment
          do i = 1, 2
             dv(iv(i)) = (vmax(iv(i)) - vmin(iv(i)))*rid(index)
          end do 

      else if (icopt.eq.3) then 
c                                 mixed variable diagrams 

c                                 no variance restriction
          isudo = 99
c                                 default search increment
          dv(iv(1)) = (vmax(iv(1)) - vmin(iv(1)))*rid(index)

      end if 

      end 

      subroutine input1 (output,n4name)
c-----------------------------------------------------------------------
c input1 reads data from a file on unit n1, this data controls the
c computational options and is modified frequently.
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'
 
      logical output, eof

      character*8 fname*10,blank*1,string(3),
     *            cname*5,rname*5,name,sname*10,
     *            vname,xname,strg*80,
     *            n2name*100,n3name*100,n4name*100,
     *            n9name*100,tcname*5,y*1,zname*5

      integer icomp,istct,iphct,icp,isoct,icp2,idum,nstrg,
     *        iff,idss,ifug,ifyn,isyn,iwt,ictr,itrans,ivfl,
     *        ibuf,hu,hv,hw,hx,jfct,jmct,jprct,
     *        icont,i,j,ierr,icmpn,iind,idep,jcont,kct

      double precision buf,v,tr,pr,r,ps,vmax,vmin,dv,ctrans,
     *                 dlnfo2,elag,gz,gy,gx,dblk,cx,c0,c1,c2,c3,
     *                 c4,c5,dip

      integer ipot,jv,iv
      common/ cst24 /ipot,jv(l2),iv(l2)

      integer imaf,idaf
      common/ cst33 /imaf(i6),idaf(i6)

      character*162 title
      common/ csta8 /title(4)
 
      common/ cst6  /icomp,istct,iphct,icp/ cst79 /isoct
     *      / csta7 /fname(h9) 
     *      / cst81 /icp2/ cst112 /buf(5)
     *      / csta2 /xname(k5),vname(l2)
     *      / cst5  /v(l2),tr,pr,r,ps/ cst9 /vmax(l2),vmin(l2),dv(l2)    
     *      / cst10 /iff(2),idss(h5),ifug,ifyn,isyn
     *      / csta4 /cname(k5)
      common/ cst209 /iwt/ cst209a /zname
     *      / cst207 /ctrans(k0,k5),ictr(k5),itrans
     *      / csta9 /tcname(k5)/ cst102 /ivfl
     *      / cst100 /dlnfo2,elag,gz,gy,gx,ibuf,hu,hv,hw,hx     
     *      / cst307 /jfct,jmct,jprct
     *      / cst314 /dblk(3,k5),cx(2),icont
     *      / cst316 /c0,c1,c2,c3,c4,c5,iind,idep

      integer ifct,idfl
      common/ cst208 /ifct,idfl

      integer ixct,iexyn,ifact
      common/ cst37 /ixct,iexyn,ifact 

      character*100 n1name,cfname
      common/ cst228 /n1name,cfname

      character*8 exname,afname
      common/ cst36 /exname(h8),afname(2)

      integer ids,isct,icp1,isat,io2
      common/ cst40 /ids(h5,h6),isct(h5),icp1,isat,io2

      integer isec,icopt,ifull,imsg,io3p
      common/ cst103 /isec,icopt,ifull,imsg,io3p

      integer jbulk
      double precision cblk
      common/ cst300 /cblk(k5),jbulk

      integer ivarrx,ivarip,isudo,ivar
      common/ cst62 /ivarrx(k2),ivarip(k2),isudo,ivar

      integer jlow,jlev,loopx,loopy,jinc
      common/ cst312 /jlow,jlev,loopx,loopy,jinc

      integer jtest,jpot
      common/ debug /jtest,jpot

      logical oned
      common/ cst82 /oned

      integer io3,io4,io9
      common / cst41 /io3,io4,io9

      logical usv
      common/ cst54 /usv

      character*8 eoscmp
      common/ cst98 /eoscmp(2)

      save blank

      data blank/' '/
c-----------------------------------------------------------------------
c                             output = .false. then in 1st cycle of
c                             autorefine.
      if (.not.output) then 
c                             read computational option file name
         do 

            write (*,1050)
            read (*,1000) n1name
            call enblnk (n1name)

            open (n1, file = n1name, iostat = ierr, status = 'old')

            if (ierr.ne.0) then
c                             system could not find the file
               write (*,1040) n1name
               read (*,1000) y
               if (y.ne.'Y'.and.y.ne.'y')
     *                        call error (120,r,n1,n1name)
c                                 try again
            else 
               exit 
            end if

         end do 
      
      else 

         open (n1, file = n1name, iostat = ierr, status = 'old')
         if (ierr.ne.0) call error (120,r,n1,n1name)

      end if 
c                                 begin reading input:

c                                 read name of thermodynamic data file
      read (n1,1000) n2name
      call enblnk (n2name)
c                                 read print and graphic file names
      read (n1,1000) n3name
      call enblnk (n3name)
      read (n1,1000) n4name
      call enblnk (n4name)
      read (n1,1000) n9name
      call enblnk (n9name)
c
      do i = 1, 4
         title(i) = ' '
      end do 
c                                 read title for the calculation:
      read (n1,1000) title(1)
      read (n1,*,err=998) icopt 
c                                 if fractionation path from data 
c                                 file, get name:
      if (icopt.eq.10) then 
         read (n1,1000) cfname
         call enblnk (cfname)
      end if 
c                                 dummy variable place holders
      read (n1,*,err=998) idum
      read (n1,*,err=998) idum 
      read (n1,*,err=998) idum 
      read (n1,*,err=998) idum
      read (n1,*,err=998) idum 
      read (n1,*,err=998) idum
      read (n1,*,err=998) idum
      read (n1,*,err=998) idum
      read (n1,*,err=998) idum

      read (n1,*,err=998) itrans
      read (n1,*,err=998) icmpn
c                                 read new component definitions:
      do i = 1, itrans
         read (n1,1120) tcname(i), ictr(i)
         read (n1,*) (ctrans(j,i), j = 1, icmpn)
      end do

      read (n1,*,err=998) iwt
c                                 dummy variable place holders
      read (n1,*,err=998) idum
      read (n1,*,err=998) idum  
      read (n1,*,err=998) idum
c                                 read code for choice of fluid equation
c                                 of state from terminal. 
      read (n1,*,err=998) ifug
      if (ifug.ge.7.and.ifug.le.12.and.ifug.ne.9.or.ifug.eq.19.or.
     *    ifug.eq.14.or.ifug.eq.16.or.ifug.eq.17.or.ifug.eq.24.or.
     *    ifug.eq.20.or.ifug.eq.25) 
     *                  read (n1,*,err=998) ibuf,hu,dlnfo2,elag
      if (ibuf.eq.5) read (n1,*,err=998) buf     

      if (hu.eq.1) then 
c                                 hardwired fluid EoS endmember names
         eoscmp(1) = 'H2      '
         eoscmp(2) = 'O2      '

      else 

         eoscmp(1) = 'H2O     '
         eoscmp(2) = 'CO2     '

      end if 
c                                 no dependent variable
      iind = 0 
c                                 dummy variable
      read (n1,*,err=998) loopx
c                                 here loopx is just a 1d/2d flag for 
c                                 gridded minimization, for backwards 
c                                 compatibility set the to 2d if > 2 or < 1.
      if (loopx.eq.1) then 
         oned = .true.
      else
         oned = .false.
      end if 

      read (n1,*,err=998) idep
      read (n1,*,err=998) c0,c1,c2,c3,c4

      if (idep.eq.1) then 
         iind = 2
      else if (idep.eq.2) then 
         iind = 1
      end if 
c                                 decode thermodynamic components
c                                 read to the beginning of the component list
      do 
         read (n1,1000,end=998) rname
         if (rname.eq.'begin') exit
      end do 
c                                 count (icp) and save names (cname)
      icp = 0
      jbulk = 0 
  
      do 

         read (n1,1000) rname,strg

         if (rname.eq.'end t') then 
c                                 finished, could check for no components
            if (icp.eq.0) then
               write (*,*) 'No thermodynamic components'
               goto 998
            else if (icopt.eq.5.and.jbulk.lt.icp) then 
               write (*,*) 'All thermodynamic components must be ',
     *                     'constrained.'
               goto 998
            end if 
        
            exit 

         else if (rname.eq.blank) then 
 
            cycle 

         else if (rname.eq.'V'.or.rname.eq.'S') then

            usv = .true.

         else

            icp = icp + 1
            icp1 = icp + 1
            icp2 = icp + 2
            cname(icp) = rname
c                                 encode a graphics names for the
c                                 compositional variables, this is kind of
c                                 pointless, but it looks good.
            write (xname(icp),1010) 'x(',rname,')'
c                                 unblank the name
            call unblnk (xname(icp))
            if (icp.gt.k5) call error (197,r,k5,'INPUT1')

         end if 
c                                 check for compositional constraints
         read (strg,*,err=998) icont

         if (icont.ne.0) then 
            jbulk = jbulk + 1
            read (strg,*,err=998) j, (dblk(i,jbulk), i = 1, icont)    
         end if 

      end do           
c                                 decode saturated components    
c                                 isat is the saturated component counter,
c                                 isyn = 1 if isat = 0, else isyn = 0
      isyn = 1
      isat = 0
      io2  = 0 

      do 
         read (n1,1000,end=998) rname
         if (rname.eq.'begin') exit
      end do 

      do 

         read (n1,1000) rname,strg
         if (rname.eq.blank) cycle 

         if (rname.eq.'end s') then 

            if (isat.ne.0) isyn = 0
            icomp = icp + isat
            exit 

         else if (rname.eq.blank) then 

            cycle 

         end if 
c                                 check for compositional constraints
         read (strg,*,err=998) jcont

         if (jcont.ne.0) then 
            jbulk = jbulk + 1
            read (strg,*,err=998) j, (dblk(i,jbulk), i = 1, icont)    
c                                 override variance flag choice, why here?
            isudo = 0    
         end if  

         isat = isat + 1
         if (isat.gt.h5) call error (15,r,i,'BUILD')
         cname(icp+isat) = rname
         if (rname.eq.'O2') io2 = isat

      end do 
c                                 decode saturated phase components
      do 
         read (n1,1000,end=998) rname
         if (rname.eq.'begin') exit
      end do 
c                                 ifct is the fluid component counter,
c                                 ifyn = 1 if ifct = 0.
      ifyn = 1
      ifct = 0

      do 

         read (n1,1000) rname

         if (rname.eq.'end s') then 
            if (ifct.ne.0) ifyn = 0
            icomp = icomp + ifct
            exit 
         else if (rname.eq.blank) then 
            cycle 
         end if 
      
         ifct = ifct + 1
         if (ifct.gt.2) call error (44,r,i,'BUILD')
c                                 save the component if only one
c                                 for use in input2.
         if (ifct.eq.1) zname = rname
         cname(icomp+ifct) = rname
      end do 
c                                  decode mobile components
c                                  jmct - mobile component counter
      jmct = 0 
      ifact = 0 

      do 

         call rdstrg (n1,nstrg,string,eof)

         if (eof) then 

            goto 998

         else if (string(1).eq.'begin') then

            cycle 

         else if (string(1).eq.'end') then

            icomp = icomp + jmct
            exit 

         else 

            read (string(1),'(a5)') rname
            jmct = jmct + 1
            if (jmct.gt.2) call error (45,r,i,'BUILD')
            cname(icomp+jmct) = rname

            if (nstrg.eq.1) then 
c                                 old format, create variable name
               write (vname(3+jmct),1000) 'mu_',rname
               imaf(jmct) = 1

            else 
c                                 new format
               read (string(2),'(a1)') y
               vname(3+jmct) = string(2)
               afname(jmct) = string(3)

               if (y.eq.'m') then 
c                                 chemical potential
                  imaf(jmct) = 1

               else if (y.eq.'f') then 

                  imaf(jmct) = 2

               else if (y.eq.'a') then 

                  imaf(jmct) = 3

               end if 

               if (imaf(jmct).gt.1) ifact = ifact + 1 

            end if 
               
         end if 

      end do 
c                             the ifct flag can probably be set later if fluid
c                             is in the thermodynamic composition space.   
      jfct = icp + isat 
c                             jprct+1..icomp -> (jmct.ne.0) mobile components 
      jprct = icomp - jmct 
c                             excluded phases
      iexyn = 1
      ixct = 0
c                             decode excluded phases
      do 
         read (n1,1000,end=998) name
         if (name.eq.'begin ex') exit
      end do

      do 

        read (n1,1000) name

         if (name.eq.'end excl') then 
            if (ixct.ne.0) iexyn = 0 
            exit
         else if (name.eq.blank) then 
            cycle 
         end if 

         ixct = ixct + 1
         if (ixct.gt.h8) call error (13,r,i,'BUILD')
         exname(ixct) = name

      end do  
c                             solution phases:
      do 
         read (n1,1000,end=998) sname
         if (sname.eq.'begin solu') exit
      end do
c                             isoct - solution phase counter,
c                             io9 is a flag = 0 no solution file
      isoct = 0

      do 

         read (n1,1000) sname
 
         if (sname.eq.'end soluti') then 
            if (io9.eq.1) isoct = 0 
            exit 
         else if (name.eq.blank) then 
            cycle  
         end if 

         isoct = isoct + 1
         if (isoct.gt.h9) call error (25,r,i,'BUILD')
         fname(isoct) = sname

      end do  
c                             read the maximum pressure, temper-
c                             ature, xco2, u1, and u2; the minimum
c                             pressure temperature, xco2, u1, and u2;
c                             and the default pressure, temperature,
c                             xco2, and chemical
c                             potential increments use kelvins, bars and
c                             joules as units (if no mobile components
c                             enter two zeroes for each read).
      read (n1,*,err=998) vmax
      read (n1,*,err=998) vmin
      read (n1,*,err=998) dv
c                             read the default indices of the
c                             dependent, independent, and secondary
c                             independent intensive variables, p = 1,
c                             t = 2, and xco2 = 3, respectively.
      read (n1,*,err=998) (iv(i), i = 1, 5)
c                             check to make sure input requests are
c                             consistent:
      if (icopt.ne.0.and.icopt.ne.4) then
c                             first check iv(1):
         if (iv(1).eq.3.and.ifyn.eq.1) call error (110,r,i,'I')
         if (iv(1).eq.3.and.ifct.eq.1) then 
            if (icopt.ne.7.and.iv(2).ne.3) call error (111,r,i,'I')
         end if 

         if (vmin(iv(1)).ge.vmax(iv(1)).and.icopt.lt.5) 
     *                                 call error (112,r,i,'I') 
         if (vname(iv(1)).eq.blank) call error (116,dip,i,'I')
      end if
c                             now check iv(2):
      if (icopt.eq.1) then
         if (iv(2).eq.3.and.ifyn.eq.1) call error (110,r,i,'INPUT1')
         if (iv(2).eq.3.and.ifct.eq.1) call error (111,r,i,'INPUT1')
         if (vmin(iv(2)).ge.vmax(iv(2))) call error (112,r,i,'I')
         if (vname(iv(2)).eq.blank) call error (116,r,i,'INPUT1')
      end if
c                             if a chemical potential is specified as an
c                             independent variable (iv(1-3)), check if
c                             the variable is defined:
      kct = 0
      do i = 1, 3
         if (iv(i).gt.3) kct = kct + 1
      end do 
c                             identify the variable used to determine
c                             which phases lie on the left hand side
c                             of a reaction equation.
      if (icopt.eq.3) then
         ivfl = iv(1)
      else if (iv(1).eq.2.or.iv(2).eq.2) then
c                             choose T
         ivfl = 2
      else if (iv(1).eq.1.or.iv(2).eq.1) then
c                             no T, so choose P
         ivfl = 1
      else
c                             no P or T, choose independent V
         if (iv(2).ne.3) then
            ivfl = iv(2)
         else
            ivfl = iv(1)
         end if
      end if
c                             ok, now find out which variables are
c                             dummies and story the indexes of the
c                             non-dummy variables in jv.
      ipot = 0

      do i = 1, 5
c                             variables v(1) (p) and v(2) (t) are
c                             only dummies if idep is set.
         if ((iv(i).ne.idep.or.icopt.eq.7.or.icopt.eq.9).and.
     *       (iv(i).eq.1.or.iv(i).eq.2)) then
            ipot = ipot+1
            jv(ipot) = iv(i)
c                             variable v(3) is a dummy if ifyn = 1:
         else if ((iv(i).eq.3).and.(ifyn.eq.0)) then
            ipot = ipot+1
            jv(ipot) = iv(i)
c                             variables v(4) and v(4) are dummies if
c                             imyn = 1:
         else if (jmct.ne.0) then
            if (iv(i).eq.4) then
               ipot = ipot+1
               jv(ipot) = iv(i)
            else if (iv(i).eq.5.and.jmct.eq.2) then
               ipot = ipot+1
               jv(ipot) = iv(i)
            end if
         end if
      end do 
c                                 if dependent variable add to jv list, could
c                                 increment ipot, but maybe it's better not to.
      if (idep.ne.0) jv(ipot+1) = idep
c                                 set convergence criteria for routine univeq
      if (icopt.le.3) call concrt

      if (icopt.ne.0) close (n1)
c                                 open files requested in input
      call ftopen (n2name,n3name,n4name,n9name,jbulk,icp,icopt,jtest)
c                                 read auxilliary input for 2d fractionation
      if (icopt.eq.9) call rdain

      goto 999
c                                 archaic error trap
998   call error (27,r,i,n2name)

1000  format (a,a)
1010  format (a2,a5,a1)
1040  format ('**error ver191** FOPEN cannot find file:',/,a,/,
     *        'Try again (y/n)?')
1050  format ('Enter the problem definition file name (i.e. the file',
     *        ' created',/,'with BUILD), left justified: ')
1120  format (a,1x,i2)

999   end

      subroutine input2 (first)
c----------------------------------------------------------------------
c input2 reads the thermodynamic data file for most perplex programs, 
c a (the?) notable exception being frendly that calls the parallel 
c routine jnput2.
c----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'
 
      character*5 cname, zname, name*8,mnames(k16*k17)*8,names*8

      double precision twt(k5),ctot,dblk,cx, cst
 
      integer iff,idss,ifug,ifyn,isyn,iwt,icomp,
     *        istct,iphct,icp,
     *        icont,i,j,ict,im,k,ifer,inames, jphct, imak(k16)
 
      logical eof, good, first

      common/ cst10 /iff(2),idss(h5),ifug,ifyn,isyn/ cst8  /names(k1)
     *      / cst3  /ctot(k1)/ csta4 /cname(k5)
     *      / cst209 /iwt/ cst209a /zname/ cst6 /icomp,istct,iphct,icp
     *      / csta6 /name
     *      / cst314 /dblk(3,k5),cx(2),icont

      integer ifct,idfl
      common/ cst208 /ifct,idfl

      integer ic
      common/ cst42 /ic(k5)

      integer idh2o,idco2,ikind,icmpn,icout
      double precision therm,comp,tot
      common/ cst43 /therm(k4),comp(k0),tot,icout(k0),idh2o,idco2,
     *               ikind,icmpn

      character cmpnt*5, dname*40
      common/ csta5 /dname,cmpnt(k0)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      double precision vmax,vmin,dv
      common/ cst9  /vmax(l2),vmin(l2),dv(l2)

      integer ilam,idiso,lamin,idsin
      double precision tm,td
      common/ cst202 /tm(m7,m6),td(m8),ilam,idiso,lamin,idsin

      integer ids,isct,icp1,isat,io2
      common/ cst40 /ids(h5,h6),isct(h5),icp1,isat,io2

      double precision atwt
      common/ cst45 /atwt(k0)

      integer iemod,kmod
      logical smod
      double precision emod
      common/ cst319 /emod(k15,k10),smod(h9),iemod(k10),kmod

      integer jbulk
      double precision cblk
      common/ cst300 /cblk(k5),jbulk

      integer ixct,iexyn,ifact
      common/ cst37 /ixct,iexyn,ifact 

      character*8 exname,afname
      common/ cst36 /exname(h8),afname(2)

      integer ipoint,imyn
      common/ cst60 /ipoint,imyn

      integer jfct,jmct,jprct
      common/ cst307 /jfct,jmct,jprct

      integer imaf,idaf
      common/ cst33 /imaf(i6),idaf(i6)

      double precision mcomp
      character mknam*8
      integer nmak
      logical mksat
      common / cst333 /mcomp(k16,k0),nmak,mksat(k16),mknam(k16,k17)

      double precision mkcoef, mdqf
      integer mknum, mkind
      common / cst334 /mkcoef(k16,k17),mdqf(k16,k17),mkind(k16,k17),
     *                 mknum(k16)

      integer make
      common / cst335 /make(k10)

      integer isfp
      common/ cst303 /isfp(k10)

      integer ikp
      common/ cst61 /ikp(k1)

      double precision vnumu
      common/ cst44 /vnumu(i6,k10)
c-----------------------------------------------------------------------
c                               initialization for each data set
c                               for k10 endmembers
      do i = 1, k10
         make(i) = 0 
      end do
c                               for k1 phases:
      do i = 1, k1
         ikp(i) = 0
      end do 
c                               other counters and flags:
      do i = 1, h5
         isct(i) = 0
      end do 
c                               counters for bounds
      iphct = 0
      lamin = 0 
      idsin = 0 
      idfl = 0
c                               read data base header, do component
c                               transformations, read make definitions.
      call topn2 (0)
c                               general input data for main program

c                               reorder thermodynamic components
c                               if the saturated phase components are 
c                               present
      ict = 0 
c                               first component 1st
      do i = 1, icp
         if (cname(i).eq.cmpnt(idh2o)) then 
            ict = 1
            if (i.eq.1) exit 
            cname(i) = cname(1)

            do j = 1, 3
               cst = dblk(j,i)
               dblk(j,i) = dblk(j,1) 
               dblk(j,1) = cst
            end do 

            cname(1) = cmpnt(idh2o)
            ict = 1
            exit            
         end if 
      end do 
c                               now check second component
      do i = 1, icp
         if (cname(i).eq.cmpnt(idco2)) then 
c                               if i = 2, already in second position exit
            if (i.eq.ict+1) exit
            cname(i) = cname(ict+1)

            do j = 1, 3
               cst = dblk(j,i)
               dblk(j,i) = dblk(j,ict+1) 
               dblk(j,ict+1) = cst
            end do 

            cname(ict+1) = cmpnt(idco2)
            ict = 1
c                              used to re-execute loop? changed to exit 10/06
            exit            
         end if 
      end do 
c                              load the old cbulk array
      if (ifyn.ne.1) iphct = 2
c                               identify nonzero components.
c                               initialize icout(i) = 0
      do i = 1, icmpn
         icout(i) = 0
      end do

      do i = 1, icomp
         im = 0
         do j = 1, icmpn
            if (cname(i).eq.cmpnt(j)) then 
               twt(i) = atwt(j)
               ic(i) = j
               icout(j) = 1
               if (j.eq.idh2o) then
                  iff(1) = i
                  idfl = idfl + 1
               else if (j.eq.idco2) then
                  iff(2) = i
                  idfl = idfl + 1
               end if 
               im = 1
            end if 
         end do 
c                               write error message if a component
c                               was not found:
         if (im.eq.0) then 
            write (*,1230) cname(i), (cmpnt(k), k = 1, icmpn)
            write (*,1240)
            stop
         end if 
 
      end do 
c                                 this segment is to check if
c                                 a possible saturated phase component
c                                 has been made a mobile component,
c                                 if there is also a saturated phase
c                                 component idfl is the identity of the
c                                 mobile component otherwise idfl = 0.
      if (ifct.eq.1.and.idfl.eq.2) then
         if (zname.eq.cmpnt(idh2o)) then
            idfl = 1
         else if (zname.eq.cmpnt(idco2)) then
            idfl = 2
         end if
      else if (idfl.eq.1.and.ifct.eq.0) then
c                                 an ad-hoc correction for stefano's bug
c                                 nov 20, 2008
         ifct = 1
         if (zname.eq.cmpnt(idh2o)) then
            idfl = 1
            xco2 = 0 
         else if (zname.eq.cmpnt(idco2)) then
            idfl = 2
            xco2 = 1d0 
         end if
         vmin(3) = xco2
         vmax(3) = xco2
      else 
         idfl = 0
      end if
c                                 load atwts in updated order
      do i = 1, icomp
         atwt(i) = twt(i)
      end do 
c                                 convert weight to molar amounts
      if (jbulk.ne.0) then 

         if (iwt.eq.1) then 
            do i = 1, jbulk
               do j = 1, 3
                  dblk(j,i) = dblk(j,i)/atwt(i)
               end do 
            end do 
         end if 

         do i = 1, jbulk
            cblk(i) = dblk(1,i)
         end do   

      end if 
c                                 get composition vectors for entities
c                                 defined by a make definition:
      call makecp (inames,mnames,first)
c                                 loop to read reference phase data for
c                                 activity/fugacity variables
      ict = 0 

      if (ifact.gt.0) then
c                                 rewind and read 'til end of header
         call eohead (n2)

         good = .false.

         do

            call getphi (name,eof)

            if (eof) then 

               write (*,1000) (afname(i),i=1,jmct)
               write (*,1010)
               stop

            end if 
c                                 now look for a match with the 
c                                 reference phase names
            do i = 1, jmct

               if (name.eq.afname(i)) then 
c                                 got a match, count
                  iphct = iphct + 1

                  ict = ict + 1

                  idaf(i) = iphct
c                                 store thermodynamic parameters:
                  call loadit (iphct)
c                                 zero the component
                  vnumu(i,iphct) = 0d0

                  if (imaf(i).eq.2) then 
c                                 if some cretin chooses fugacity, prevent
c                                 gphase from calling the EoS.   
                     isfp(iphct) = 0 

                  else 
c                                 check for special component names
c                                 this is necessary because loadit 
c                                 will not set isfp if ifyn = 0.
                     if (name.eq.cmpnt(idh2o)) then 
                        isfp(iphct) = 1
                     else if (name.eq.cmpnt(idco2)) then
                        isfp(iphct) = 2 
                     end if
 
                  end if 
c                                 blank the name, this has two purposes,
c                                 it prevents problems if an entry is 
c                                 replicated in the data file, and flags
c                                 tagged entries 
                  afname(i) = ' '

                  if (ict.eq.jmct) good = .true.

                  exit 

               end if 

            end do 

            if (good) exit 

         end do 

      end if 
c                                 begin first read loop for data on
c                                 saturated components.
      if (isyn.ne.0.and.ifyn.ne.0) goto 40
c                                 read 'til end of header
      call eohead (n2)
c                                 loop to read real saturated
c                                 entities:
      ifer = 0

      do 

         call getphi (name,eof)

         if (eof) exit
 
         call chkphi (0,name,good)

         if (good) call sattst (ifer)

      end do 
c                                 loop to load made saturated entities
      do i = 1, nmak

         if (.not.mksat(i)) cycle
c                                 load make data 
         do j = 1, icmpn
            comp(j) = mcomp(i,j)
         end do 

         name = mknam(i,mknum(i)+1)
c                                 redundant check:
         call chkphi (2,name,good)
c                               
         if (.not.good) call error (57,comp(1),iphct,name)

         call sattst (ifer)

         make(iphct) = i
c                                 pointer used for iemod.
         imak(i) = iphct

      end do 
c                                 check that there is data for
c                                 every fluid component.
      if (ifyn.eq.0.and.ifer.ne.ifct) call error (36,r,i,'INPUT2')
c                                 check that there is one phase
c                                 for each saturation constraint
40    do i = 1, isat
         if (isct(i).lt.1) call error (15,r,i,cname(icp+i))
      end do
c                                 read data for the remaining
c                                 phases of appropriate composition.
      istct = iphct + 1
c                                 read till end of header
      call eohead (n2)
c                                 begin second read loop:

      do  
    
         call getphi (name,eof)

         if (eof) exit 
c                               check if valid phase:
         call chkphi (1,name,good)

         if (good) then 
c                               acceptable data, count the phase:
            iphct = iphct + 1
c                               for normalized composition:
            ctot(iphct) = tot
c                               store thermodynamic parameters:
            call loadit (iphct)
         end if 
      end do 

c                                loop to load made entities
      do i = 1, nmak

         if (mksat(i)) cycle
c                                load make data 
         do j = 1, icmpn
            comp(j) = mcomp(i,j)
         end do 

         name = mknam(i,mknum(i)+1)
c                                redundant check, but makes ctot.
         call chkphi (3,name,good)
c                               
         if (.not.good) call error (57,comp(1),iphct,name)

         iphct = iphct + 1
         ctot(iphct) = tot
         call loadit (iphct)

         make(iphct) = i
c                                pointer used for iemod.
         imak(i) = iphct

      end do 
c                                get/save data for makes, this
c                                data is saved in the arrays thermo
c                                and cp by loadit, but are not counted,
c                                i.e., the counters ipoint and iphct
c                                are reset. soload will then load the
c                                cp array over the values loaded here,
c                                but thermo should not be affected. gmake
c                                then gets the data using the array 
c                                mkind. the names array will also be 
c                                overwritten.
      jphct = iphct
c                                read header
      call eohead (n2)

      do 

         call getphi (name,eof)

         if (eof) exit

         do i = 1, inames

            if (name.ne.mnames(i)) cycle
c                                matched a name
            iphct = iphct + 1
c                                store thermodynamic parameters:
            call loadit (iphct)

         end do 

      end do 

      do i = 1, nmak
c                                remake pointer array for makes 
         do j = 1, mknum(i)
            do k = jphct + 1, iphct
               if (names(k).ne.mknam(i,j)) cycle
               mkind(i,j) = k
            end do
         end do 
      end do  
c                                reset ipoint counter, but do not 
c                                reset iphct, because the compositions
c                                of the make phases are necessary for
c                                chemical potential variables. 
      iphct = jphct 
      ipoint = jphct

      do 10 i = 1, nmak
c                                make an iemod flag for made
c                                endmembers:   
         do j = 1, mknum(i)
           if (iemod(mkind(i,j)).eq.0) goto 10 
         end do 

         iemod(imak(i)) = iemod(mkind(i,1))

10    end do 

1000  format ('**error ver007** at least one of the reference ',
     *        'phases:',/,5(a,1x))
1010  format ('needed to define an independent fugacity/activity ',
     *    'variable is missing from the',/,'thermodynamic data file',/)
1230  format ('**error ver013** ',a,' is an incorrect component'
     *       ,' name, valid names are:',/,12(1x,a))
1240  format ('check for upper/lower case matches or extra blanks',/)

      close (n2)

      end

      block data
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5  /p,t,xco2,u1,u2,tr,pr,r,ps

      double precision pp,tt,yy,rr
      common/ cst85 /pp,tt,yy,rr

      integer iff,idss,ifug,ifyn,isyn
      common / cst10 /iff(2),idss(h5),ifug,ifyn,isyn
    
      double precision ptx
      integer ipt2
      common/ cst32 /ptx(l5),ipt2

      double precision thermo,uf,us
      common/ cst1  /thermo(k4,k10),uf(2),us(h5)

      logical gflu,aflu,fluid,shear,lflu,volume
      common/ cxt20 /gflu,aflu,fluid(k5),shear,lflu,volume
c                                 interval limits conformal transformation
      integer intv
      double precision yint, yfrc
      common/ cst47 /yint(5,ms1,mst,h9),yfrc(4,ms1,mst,h9),intv(4)
c-----------------------------------------------------------------------
      data iff/2*0/,ipt2/0/
c
      data us, uf/ h5*0d0, 2*0d0/

      data r,rr/8.3144126d0,83.144126d0/

      data gflu/.false./

      data intv/1,2,4,1/

      end

      subroutine setvr0 (i,j)
c--------------------------------------------------------------------
c setvr1 computes nodal variables for node ij, three cases:

c icont = 1 -> independent variables are the 1st and 2nd potentials
c icont = 2 -> 1st independent variable is a composition variable,  
c              2nd independent variable is the 1st potential (iv1)
c icont = 3 -> independent variables are compositional variables

c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,j

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      double precision vmax,vmin,dv
      common/ cst9  /vmax(l2),vmin(l2),dv(l2)

      integer ipot,jv,iv1,iv2,iv3,iv4,iv5
      common / cst24 /ipot,jv(l2),iv1,iv2,iv3,iv4,iv5

      integer jvar
      double precision var,dvr,vmn,vmx
      common/ cxt18 /var(l3),dvr(2),vmn(l3),vmx(l3),jvar

      integer icont
      double precision dblk,cx
      common/ cst314 /dblk(3,k5),cx(2),icont
c----------------------------------------------------------------------

      if (icont.eq.1) then 

         v(iv1) = vmin(iv1) + (i-1)*dv(iv1)
         v(iv2) = vmin(iv2) + (j-1)*dv(iv2)
         call incdp0

      else if (icont.eq.2) then 

         v(iv1) = vmin(iv1) + (j-1)*dv(iv1)
         call incdep (iv1)

         cx(1) =  (i-1)*dvr(1)
         call setblk 

      else 

         cx(1) = (i-1) * dvr(1)
         cx(2) = (j-1) * dvr(2)
         call setblk

      end if 

      end

      subroutine setblk
c-----------------------------------------------------------------------
c for gridded minimization setblk computes the bulk composition
c and initializes the arrays for lpopt.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision ctotal,x2,x1

      integer i,j,icomp,istct,iphct,icp

      common/ cst6  /icomp,istct,iphct,icp

      double precision a,b,c
      common/ cst313 /a(k5,k1),b(k5),c(k1)

      integer jbulk
      double precision cblk
      common/ cst300 /cblk(k5),jbulk

      integer icont
      double precision dblk,cx
      common/ cst314 /dblk(3,k5),cx(2),icont

      logical usv
      common/ cst54 /usv
c-----------------------------------------------------------------------

      if (icont.eq.2) then
 
         x1 = 1d0 - cx(1)
         x2 = cx(1)
c                                 bulk composition computed as
c                                 c = (1-cx(1))*dblk(1) + cx(1)*dblk(2)
         do j = 1, jbulk
            cblk(j) = x1*dblk(1,j) + x2*dblk(2,j)
         end do 

      else

         x1 = cx(1)
         x2 = cx(2)
c                                 bulk composition computed as
c                                 c = dblk(1) + cx(1)*dblk(2) + cx(2)*dblk(3)  
         do j = 1, jbulk
            cblk(j) = dblk(1,j) + x1*dblk(2,j) + x2*dblk(3,j)
         end do 

      end if 
c                                 modify cblk here to change the 
c                                 composition before minimization.
      ctotal = 0d0

      if (usv) then

         j = jbulk

      else 

         j = icp 

      end if  
c                                 get total moles to compute mole fractions             
      do i = 1, j
         ctotal = ctotal + cblk(i)
      end do

      do i = 1, j 
         b(i) = cblk(i)/ctotal
      end do

      end 


      subroutine setvar 
c--------------------------------------------------------------------
c setvar initializes the variables for gridded minimization, three
c cases:

c icont = 1 -> independent variables are the 1st and 2nd potentials
c icont = 2 -> 1st independent variable is a composition variable,  
c              2nd independent variable is the 1st potential (iv1)
c icont = 3 -> independent variables are compositional variables
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,j

      double precision v,tr,pr,r,ps,vmax,vmin,dv,
     *                 rloopy,rloopx

      common/ cst5  /v(l2),tr,pr,r,ps
     *      / cst9  /vmax(l2),vmin(l2),dv(l2)

      integer ipot,jv,iv1,iv2,iv3,iv4,iv5
      common / cst24 /ipot,jv(l2),iv1,iv2,iv3,iv4,iv5

      integer jbulk
      double precision cblk
      common/ cst300 /cblk(k5),jbulk

      integer jlow,jlev,loopx,loopy,jinc
      common/ cst312 /jlow,jlev,loopx,loopy,jinc

      integer jvar
      double precision var,dvr,vmn,vmx
      common/ cxt18 /var(l3),dvr(2),vmn(l3),vmx(l3),jvar

      integer icont
      double precision dblk,cx
      common/ cst314 /dblk(3,k5),cx(2),icont

      integer isec,icopt,ifull,imsg,io3p
      common/ cst103 /isec,icopt,ifull,imsg,io3p
c----------------------------------------------------------------------

      rloopy = dfloat(loopy-1)
      rloopx = dfloat(loopx-1)
c                                 for 1d calculations
      if (rloopx.eq.0) rloopx = rloopy

      do i = 1, ipot
         v(jv(i)) = vmin(jv(i))
      end do

      call incdp0

      if (icopt.eq.10) then 
c                                using nodal coordinate system
         dvr(1) = 1

      else if (icont.eq.1) then 
c                                v(iv1) on x, v(iv2) on y
         dv(iv1) = (vmax(iv1) - vmin(iv1))/rloopx
         dvr(1) = dv(iv1)

         dv(iv2) = (vmax(iv2) - vmin(iv2))/rloopy
         dvr(2) = dv(iv2)

      else if (icont.eq.2) then 
c                               composition is on x, v(iv1) on y
         dvr(1) = 1d0/rloopx

         dv(iv1) = (vmax(iv1) - vmin(iv1))/rloopy
         dvr(2) = dv(iv1)

      else 
c                                compositions on both axes
         dvr(1) = 1d0/rloopx
         dvr(2) = 1d0/rloopy 
         cx(1) = 0d0
         cx(2) = 0d0

      end if 
c                                set the bulk composition:
      do j = 1, jbulk
         cblk(j) = dblk(1,j)
      end do 

      end 

      subroutine inipot 
c--------------------------------------------------------------------
c setvar initializes the independent potential variables to their 
c minimum values
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      double precision vmax,vmin,dv
      common/ cst9  /vmax(l2),vmin(l2),dv(l2)

      integer ipot,jv,iv
      common / cst24 /ipot,jv(l2),iv(l2)
c----------------------------------------------------------------------
c                                 initialize potentials
      do i = 1, ipot
         v(jv(i)) = vmin(jv(i))
      end do 
c                                 set dependent potential, if it exists
      call incdp0

      end 

      subroutine getcmp (jd,id,ids)
c-----------------------------------------------------------------------
c getcmp gets the composition of pseudocompund id, where:
c  if ids < 0, -ids points to the composition of a true compound in array cp
c  if ids > 0, id points to the composition of a solution defined in terms
c              on endmember fractions defined and saved by routine resub
c              in array zcoor.
c the composition is saved in arrays cp3 and x3, entry jd
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i,j,id,jd,ids
c                                 -------------------------------------
c                                 global variables:
      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      double precision cp
      common/ cst12 /cp(k5,k1)
c                                 bookkeeping variables
      integer ksmod, ksite, kmsol, knsp
      common/ cxt0  /ksmod(h9),ksite(h9),kmsol(h9,m4,mst),knsp(m4,h9)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)
c                                 working arrays
      double precision z, pa, p0a, x, w, y
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1)
c                                 single site solution coordinates:
      integer jend
      common/ cxt23 /jend(h9,k12)
c                                 refined compositions and solution 
c                                 pointer
      integer kkp,np,ncpd,ntot
      double precision cp3,ctot3
      common/ cxt15 /cp3(k5,k5),ctot3(k5),kkp(k5),np,ncpd,ntot

      integer ikp
      common/ cst61 /ikp(k1)
c----------------------------------------------------------------------

      kkp(jd) = ids
      ctot3(jd) = 0d0

      if (ids.lt.0) then 
c                                 simple compounds
         do i = 1, icomp
            cp3(i,jd) = cp(i,-ids)
         end do 
c                                 check if it's a solution endmember
         if (ikp(-ids).ne.0) call endcp (jd,-ids,ikp(-ids))

      else
c                                 solutions 
c                                 initialize
         do i = 1, icomp
            cp3(i,jd) = 0d0
         end do
c                                 get the x(i,j) coordinates for the
c                                 composition from the zcoor array,
c                                 this routine also saves a copy of the
c                                 xcoordinates in x3(jd,i,j)
         call getxz (jd,id,ids)
c                                 convert the x(i,j) coordinates to the
c                                 geometric y coordinates
         call xtoy (ids)

         if (ksmod(ids).le.6.or.ksmod(ids).ge.10) then
c                                 for solutions with no dependent endmembers
c                                 the y coordinates can be used 
c                                 to compute the composition
            do i = 1, mstot(ids)
               do j = 1, icomp
                  cp3(j,jd) = cp3(j,jd) + y(i) * cp(j,jend(ids,2+i))
               end do
            end do

         else 
c                                 get the p' coordinates (amounts of 
c                                 the independent disordered endmembers)     
            call getpp (ids) 

            do i = 1, lstot(ids)
               do j = 1, icomp 
                  cp3(j,jd) = cp3(j,jd) + p0a(i) * cp(j,jend(ids,2+i))
               end do 
            end do          

         end if 

      end if 

      do i = 1, icomp
         ctot3(jd) = ctot3(jd) + cp3(i,jd)
      end do 

      end 

      subroutine getpp (id)
c-----------------------------------------------------------------------
c getpp computes the amounts of the indepdendent edmembers of a reciprocal
c solution in terms of the disordered endmembers (i.e., the p coordinates
c corrected for the amounts of the ordered species if present [ksmod=8]).
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer id

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)
c----------------------------------------------------------------------
c                                  first convert the istot disordered
c                                  endmember coordinates to the 
c                                  kstot + nord p0 coordinates
      call y2p0 (id) 
c                                  decompose ordered species
      if (nord(id).gt.0) call p0dord (id)

      end

      subroutine inblnk (text,char)
c----------------------------------------------------------------------
c inblnk - scan text to last '/' or '\' and insert char after.
 
c     text - character string 
c----------------------------------------------------------------------
      implicit none

      integer i, nchar
 
      character text*(*), bitsy(400)*1, char*1 
c----------------------------------------------------------------------
      nchar = len(text) 
      read (text,1000) (bitsy(i), i = 1, nchar)
c                                 scan for blanks:

      do i = nchar,1,-1
c                                 this line may cause problems
c                                 on some operating systems that 
c                                 recognize the backslash as an escape
c                                 character.
         if (bitsy(i).eq.'/') goto 10
         bitsy(i+1) = bitsy(i)
      end do 

      i = 0

10    bitsy(i+1) = char

      write (text,1000) (bitsy(i), i = 1, nchar)
 
1000  format (400a1)
      end

      subroutine matchj (unnown,itis)
c----------------------------------------------------------------------
 
c matchj - subroutine to determine if the string unnown is a valid
c          solution or compound name.
 
c   itis = -id if compound
c   itis = ikp if solution 
c   itis = 0 if invalid
c----------------------------------------------------------------------
      implicit none

      integer i, itis
 
      character*10 unnown
 
      include 'perplex_parameters.h'
 
      integer isoct
      common/ cst79 /isoct

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      character*10 sname, names*8
      common/ cst8 /names(k1)/ csta7 /sname(h9)
c---------------------------------------------------------------------- 
 
      itis = 0

      do i = 1, isoct
         if (unnown.eq.sname(i)) then
             itis = i
             goto 99
         end if
      end do

      do i = 1, iphct
         if (unnown.eq.names(i)) then
            itis = -i
            goto 99
         end if
      end do 

99    end

      subroutine maktit 
c-----------------------------------------------------------------------
c create a title for graphics output, the title consists of the 
c calculation title + saturation hierarchy (provided one is 
c specified) and is the first two elements of title (csta8).
c if icopt = 1 or 3, also adds a blurb about reaction convention.

c title is max 3 lines, but four lines are written to be consistent
c with old plot file formats written by frendly, pt2curv etc.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      character*8 vname,xname,cname*5

      integer ivfl,i

      character*162 title
      common/ csta8 /title(4)

      common/ csta2  /xname(k5),vname(l2)
     *      / cst102 /ivfl
     *      / csta4 /cname(k5)

      integer isec,icopt,ifull,imsg,io3p
      common/ cst103 /isec,icopt,ifull,imsg,io3p

      integer iff,idss,ifug,ifyn,isyn
      common/ cst10  /iff(2),idss(h5),ifug,ifyn,isyn

      integer ids,isct,icp1,isat,io2
      common/ cst40 /ids(h5,h6),isct(h5),icp1,isat,io2

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  
c-----------------------------------------------------------------------
      do i = 2, 4
         title(i) = ' '
      end do                              
c                               saturated and buffered component names:
      if (isyn.eq.0) then 
         write (title(2),1070) (cname(i+icp), i= 1, isat)
      else 
         write (title(2),1000) ' '
      end if 
c                                 reaction convention
      if (icopt.eq.1.or.icopt.eq.3) write (title(3),1080) vname(ivfl)

      do i = 1, 3
         call deblnk (title(i))
      end do 

1000  format (a)
1070  format ('Component saturation hierarchy: ',7(a,1x))
1080  format ('Reaction equations are written with the high ',
     *         a,'assemblage to the right of the = sign')

      end

      subroutine rdain
c-----------------------------------------------------------------------
c a subprogram to read auxilliary input file for 2d fractionation 
c calculations, called by VERTEX and WERAMI
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer maxlay

      parameter (maxlay=6) 

      integer i,ierr

      double precision zlayer

      character*100 aname

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      character*100 n1name,cfname
      common/ cst228 /n1name,cfname

      integer ipot,jv,iv
      common/ cst24 /ipot,jv(l2),iv(l2)

      double precision vmax,vmin,dv
      common/ cst9  /vmax(l2),vmin(l2),dv(l2)

      integer gloopy,ilay,irep
      double precision a0,a1,a2,a3,b0,b1,b2,b3,c0,c1,c2,c3,dv1dz,
     *               zbox,iblk
      common/ cst66 /a0,a1,a2,a3,b0,b1,b2,b3,c0,c1,c2,c3,dv1dz,
     *               zbox,iblk(maxlay,k5),gloopy,ilay,irep(maxlay)
c-----------------------------------------------------------------------
c                                 create auxilliary input file (a+input file name)
      aname = n1name
c                                 look for input data from a file 
c                                 named by prefixing the problem definition 
c                                 file with 'a'
      call inblnk (aname,'a')
c
      open (n8,file=aname,status='old',iostat=ierr)

      if (ierr.ne.0) call error (51,zbox,gloopy,aname) 
c                                 set the number of independent variables
c                                 to 1 (the independent path variable must
c                                 be variable jv(1), and the dependent path
c                                 variable must be jv(2), the path variables
c                                 can only be pressure and temperature
      ipot = 1
c                                 in old versions this was the number of steps 
c                                 along the path, in 07 loopy is now set via
c                                 jlow/1dpath (above). gloopy is not used as a 
c                                 flag if top layer composition is to be refreshed
c                                 after each step (gloopy=999).
      read (n8,*) gloopy
c                                 thickness of a box in column
      read (n8,*) zbox 
c                                 gradient in variable jv(1) with z, jv(1)
c                                 is the independent variable, for subduction
c                                 this is logically pressure, i.e.,
c                                 dp(bar)/dz(m)
      read (n8,*) dv1dz 
c                                 now we need a path function for the dependent
c                                 variable, here we take a function defined in
c                                 terms of the absolute depth of the top of the
c                                 column (z0) and the relative depth (dz) within
c                                 the column
c                                 
c                                 v2 = a(z0)*dz^2 + b(z0)*dz + c(z0)

c                                 e.g., T(K) =  a(z0)*dz^2 + b(z0)*dz + c(z0)

c                                 where a(z0) = a0 + a1*z0 + a2*z0^2 + a3*z0^3 + ...
c                                 b(z0) = b0 + b1*z0 + b2*z0^2 + b3*z0^3 + ...
c                                 c(z0) = c0 + c1*z0 + c2*z0^2 + c3*z0^3 + ...
      read (n8,*) a0, a1, a2, a3
      read (n8,*) b0, b1, b2, b3
      read (n8,*) c0, c1, c2, c3
c                                 get the initial global composition array
c                                 consisting of ibox compositions defined 
c                                 in terms of icp components. this read
c                                 statement assumes that H2O an CO2 (if 
c                                 thermodynamic components) are the 1st and
c                                 2nd components (if present). 
      ilay = 0
      vmax(2) = 0d0
      vmin(2) = 0d0
c                                 number of nodes with appended composition
c                                 end of data indicated by zero 
      do 

         read (n8,*) zlayer

      if (zlayer.eq.0) exit 

         ilay = ilay + 1

         if (ilay.eq.maxlay) then 
            write (*,*) 'increase maxlay in routine FRAC2D'
            stop
         end if 

         read (n8,*) (iblk(ilay,i),i=1,icp)

         irep(ilay) = idint(zlayer/zbox)
c                                 set the y coodinate to depth below top
         vmin(2) = vmin(2) - irep(ilay)*zbox

      end do 

      close (n8)

      end 

      subroutine fr2dpt (p0,dz)
c----------------------------------------------------------------------
c subroutine to set p-t variables from i-j coordinates in 2d-fractionation
c calculations, called by VERTEX and WERAMI
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer maxlay

      parameter (maxlay=6) 

      double precision p0, z0, dz

      integer jlow,jlev,loopx,loopy,jinc
      common/ cst312 /jlow,jlev,loopx,loopy,jinc

      integer gloopy,ilay,irep
      double precision a0,a1,a2,a3,b0,b1,b2,b3,c0,c1,c2,c3,dv1dz,
     *               zbox,iblk
      common/ cst66 /a0,a1,a2,a3,b0,b1,b2,b3,c0,c1,c2,c3,dv1dz,
     *               zbox,iblk(maxlay,k5),gloopy,ilay,irep(maxlay)

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps
c----------------------------------------------------------------------
c                                 convert to depth at top of column
      z0 = p0/dv1dz
c                                 set the independent variable
      v(1) = p0 + dz * dv1dz   
c                                 set the dependent variable
      v(2) = (a0 + a1*z0 + a2*z0**2 + a3*z0**3)*dz**2 
     *     + (b0 + b1*z0 + b2*z0**2 + b3*z0**3)*dz 
     *     +  c0 + c1*z0 + c2*z0**2 + c3*z0**3                       
      end 
   