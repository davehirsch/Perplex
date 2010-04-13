      program build 
c----------------------------------------------------------------------
c                       ************************
c                       *                      *
c                       *  build.mar.10.1991   *
c                       *                      *
c                       ************************
c----------------------------------------------------------------------
c build is a fortran program which interactively creates the input
c file read from unit n1 by the vertex program.  build reads data
c from three sources n2, *, and n9.  the output file is written to
c unit n1
c-----------------------------------------------------------------------
c files (see vertex program documentation for additional information):
c
c  l.u.n  i/o       purpose
c  -----  ---   --------------------------------------------------------
c   n2     i    data file containing the names, compositional vectors,
c               and standard state thermodynamic data for potentially
c               stable condensed phases and fluid species.
c   n9     i    optional data file which contains data on the
c               compositional equations of state and subdivision schemes
c               to be used for solution phases.
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      character*8 fugact(3), phase, b1, b2, mnames(k16*k17)*8

      character*100 n3name, n4name, n9name, cfname, dpl, dpr, 
     *              dsol,  char6*6, opname

      character*5 mname(k5), nname(k5), oname(k5), qname(k0), pname(k5), 
     *            char5, uname(k0)

      character amount*6,nc(3)*2,text*256, dtext*200, new*3, blank*1, 
     *          stext*11, title*162, y*1

      character*10 blah, fname(h9), tname(i9)

      integer i, iind, ifugy, im, idum, ima, jsat(h5),
     *        ifct, ivct, iwt, jcth, icth, igood,
     *        jcmpn, j, ier,icopt,idep,
     *        ict, jvct, jc, ix, jst,ind, loopy, 
     *        loopx, ierr, idsol, isoct, kvct, inames

      logical eof, good, vertex, satflu, mobflu, oned, findph, bad

      double precision c(0:4)

      character*8 name
      common/ csta6 /name

      integer ids,isct,icp1,isat,io2
      common/ cst40 /ids(h5,h6),isct(h5),icp1,isat,io2

      integer jfct,jmct,jprct
      common/ cst307 /jfct,jmct,jprct

      integer imaf,idaf
      common/ cst33 /imaf(i6),idaf(i6)

      integer grid
      double precision rid 
      common/ cst327 /grid(5,2),rid(2)
 
      integer ictr, itrans
      double precision ctrans
      common/ cst207 /ctrans(k0,k5),ictr(k5),itrans

      integer ipoint,imyn
      common/ cst60 /ipoint,imyn

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      integer inew
      common/ cst140 /inew

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      integer io3,io4,io9
      common / cst41 /io3,io4,io9

      double precision vmax,vmin,dv
      common/ cst9  /vmax(l2),vmin(l2),dv(l2)

      character*8 names
      common/ cst8 /names(k1)

      character*8 vname, xname
      common/ csta2 /xname(k5),vname(l2)

      integer iff,idss,ifug,ifyn,isyn
      common/ cst10  /iff(2),idss(h5),ifug,ifyn,isyn

      double precision buf
      common/ cst112 /buf(5)

      integer ibuf,hu,hv,hw,hx   
      double precision dlnfo2,elag,gz,gy,gx
      common/ cst100 /dlnfo2,elag,gz,gy,gx,ibuf,hu,hv,hw,hx
  
      integer ixct,iexyn,ifact
      common/ cst37 /ixct,iexyn,ifact 

      character*8 exname,afname
      common/ cst36 /exname(h8),afname(2)

      integer idh2o,idco2,ikind,icmpn,icout
      double precision therm,comp,tot
      common/ cst43 /therm(k4),comp(k0),tot,icout(k0),idh2o,idco2,
     *               ikind,icmpn

      integer ic
      common/ cst42 /ic(k5)

      integer icont
      double precision dblk,cx
      common/ cst314 /dblk(3,k5),cx(2),icont

      character cmpnt*5, dname*40
      common/ csta5 /dname,cmpnt(k0)

      double precision wg,xmn,xmx,xnc
      integer iend,isub,imd,insp,ist,isp,isite,iterm,iord,istot,
     *        jstot,kstot
      common/ cst108 /wg(m1,m3),xmn(mst,msp),xmx(mst,msp),xnc(mst,msp),
     *      iend(m4),isub(m1,m2,2),imd(ms1,mst),insp(m4),ist(mst),
     *      isp(mst),isite,iterm,iord,istot,jstot,kstot

      integer jsmod
      double precision vlaar
      common/ cst221 /vlaar(m3,m4),jsmod

      integer ipot,jv,iv
      common/ cst24 /ipot,jv(l2),iv(l2)

      double precision mcomp
      character mknam*8
      integer nmak
      logical mksat
      common / cst333 /mcomp(k16,k0),nmak,mksat(k16),mknam(k16,k17)

      double precision mkcoef, mdqf
      integer mknum, mkind
      common / cst334 /mkcoef(k16,k17),mdqf(k16,k17),mkind(k16,k17),
     *                 mknum(k16)

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      data blank,fugact/' ','chem_pot','fugacity','activity'/

      data dpr,dpl,dsol/'pr            ','pl           ',
     *                  'solution_models.dat     '/
c-----------------------------------------------------------------------
      write (*,7020)
c                                     assign data files
c                                     get thermo data file choice
      call ftopen
c                                     get computational option file name
      write (*,1170) 
      read (*,'(a)') opname
      
      if (opname.eq.' ') opname = 'perplex_option.dat'

c                                     initialization:
      n3name = blank
      n4name = blank
      n9name = blank

      do i = 1, l2
         iv(i) = i
         vmin(i) = 0d0
         vmax(i) = 0d0
         dv(i) = 0d0 
      end do 
 
      phase ='  FLUID ' 
      amount = 'molar '
      dtext = ' '

      nc(1) = 'C0'
      nc(2) = 'C1'
      nc(3) = 'C2'
 
      iv(5) = 3
      idum = 0 
      iind = 0 
      ifyn = 0
      iexyn = 0
      ixct = 0
      ifugy = 0
      im = 0
      isat=0
      jmct = 0
      ifct = 0
      io3 = 1
      io4 = 1
      ivct = 2
      icont = 0
      char5 = 'b' 
      iphct = 0
      isoct = 0 
      iwt = 0
      jcth = 0 
      icth = 0 
      satflu = .false.
      mobflu = .false.
      vertex = .true.

      call redop1 (satflu,opname)

      do i = 0, 4
         c(i) = 0d0
      end do 
c                                 Read THERMODYNAMIC DATA file header
      call topn2 (3)

      jcmpn = icmpn
c                                 Component stuff first:
      do i = 1, icmpn
         uname(i) = cmpnt(i)
         qname(i) = cmpnt(i)
      end do
c                                 components of saturated phase:
      write (*,2030) phase
      read (*,2100) y

      if (y.eq.'y'.or.y.eq.'Y') then
c                                   write prompt
         write (*,2031) 
         write (*,1040) uname(idh2o),uname(idco2)
         write (*,2021)
c                                   write blurb
         write (*,1030)

         do 

            read (*,3000) char5 

            if (char5.ne.blank) then
c                                 check choice 
               call chknam (igood,jcmpn,1,good,char5,qname,uname)

               if (.not.good) cycle
c                                 count component
               ifct = ifct + 1  
               mname(ifct) = char5

               if (ifct.eq.2) then 
                  write (*,*)
                  exit 
               end if   
            else
c                                 blank input 
                exit 

            end if 
         end do 

         if (ifct.ne.0) then 
            ifyn = 1
            iv(3) = 3
            ivct = 3           
         end if 

      end if
c                                 saturated components
      write (*,2110)
      read (*,2100) y
 
      if (y.eq.'Y'.or.y.eq.'y') then

         call warn (15,r,i,'BUILD')
 
         write (*,2032) h5+1
         write (*,1040) (qname(I),I=1,JCMPN)
         write (*,2021)

         do

            read (*,3000) char5 

            if (char5.ne.blank) then
c                                 check if it's a saturated phase component:
               if (char5.eq.uname(idh2o).or.char5.eq.uname(idco2)) then 
                  if (ifct.gt.0) then 
                     call warn (7,r,i,'BUILD')
                     write (*,1000) 
                     read (*,2100) y
                     if (y.eq.'y'.or.y.eq.'Y') then
                        satflu = .true.
                        write (*,1070)
                     else 
                        write (*,1080) 
                        cycle
                     end if
                  else 
                     satflu = .true.
                  end if  
               end if 

               call chknam (igood,jcmpn,1,good,char5,qname,uname)

               if (.not.good) cycle
c                                  count component
               isat = isat + 1    
               jsat(isat) = igood
               nname(isat) = char5

            else 
c                                  blank input
               exit   
             
            end if 

         end do 
 
      end if
c                                 mobile components:
      write (*,2040)
      read (*,2100) y
 
      if (y.eq.'y'.or.y.eq.'Y') then
         
         do 
c                                 write prompt
            write (*,2050) 
            write (*,1040) (qname(i),i=1,jcmpn)
            read (*,3000) char5 
            if (char5.eq.blank) exit 
c                                 check if it's a saturated phase component,
c                                 or a saturated component:
            if (char5.eq.uname(idh2o).or.char5.eq.uname(idco2)) then
               if (ifct.gt.0.or.satflu) then  
                  call warn (5,r,i,'BUILD')
                  write (*,1000) 
                  read (*,2100) y
                  if (y.eq.'y'.or.y.eq.'Y') then
                     mobflu = .true.
                     write (*,1070)
                  else 
                     write (*,1080) 
                     cycle
                  end if
               else
                  mobflu = .true.
               end if 
            end if 

            call chknam (igood,jcmpn,1,good,char5,qname,uname)

            if (.not.good) cycle
c                                 count component
            jmct = jmct + 1    
            ivct = ivct + 1
            oname(jmct) = char5
c                                 ask if chemical potential, activity
c                                 or fugacity
            write (*,2060) char5
            do 
               read (*,*,iostat=ier) ima
               if (ier.eq.0) exit
            end do 
c                                 if fugacity/activity get the reference
c                                 phase
            if (ima.eq.2.or.ima.eq.3) then
                     
               call topn2 (2)

               do 
 
                  call getphi (name,eof)

                  if (eof) exit

                  if (findph(igood)) then
                     iphct = iphct + 1
                     names(iphct) = name
                  end if 

               end do 
c                                 names contains the list of candidates
               if (iphct.eq.0) then
                  write (*,2061) fugact(ima),char5
                  ima = 1
                  imaf(jmct) = ima
               else if (iphct.eq.1) then 
                  write (*,2062) names(1),char5,fugact(ima)
                  afname(jmct) = names(1)
                  icp = 1
               else 
                  write (*,2063) char5,fugact(ima)
                  write (*,2064) (names(i),i=1,iphct)
                  good = .false.
                  do 
                     read (*,2100) name
                     do i = 1, iphct
                        if (name.eq.names(i)) then
                           good = .true.
                           icp = i  
                           afname(jmct) = name
                           exit 
                        end if
                     end do 
                     if (good) exit
                     write (*,2310) name
                  end do 
               end if 
            else 
               ima = 1
               afname(jmct) = blank
            end if
c                                 now make the variable name
            if (ima.eq.2.or.ima.eq.3) then

               read (names(icp),2100) char6
               if (ima.eq.2) then 
                  write (vname(3+jmct),5000) 'f',char6
               else 
                  write (vname(3+jmct),5000) 'a',char6
               end if

               write (*,2065) char5,fugact(ima),vname(3+jmct)
c                                 if name is a special phase component 
c                                 set flag for eos request
               if (names(icp).eq.uname(idh2o).or.
     *             names(icp).eq.uname(idco2)) ifugy = ifugy + 1

            else 

               write (vname(3+jmct),5000) 'mu',char5                 
               write (*,2066) char5,vname(3+jmct)
c                                 if component is a special phase component 
c                                 set flag for eos request
               if (char5.eq.uname(idh2o).or.
     *             char5.eq.uname(idco2)) ifugy = ifugy + 1

            end if 
                
            imaf(jmct) = ima
            iphct = 0 
            if (jmct.eq.i6) exit  

         end do 
c                                 reset n2 
         call topn2 (2)
c                                 end of mobile component loop
      end if
 
      if (ifyn.eq.0) then 
c                                 no saturated fluid phase 
         iv(3)=4
         iv(4)=5

      else 
c                                 a saturated phase is present
c                                 but may have only one component
         iv(4)=4
         iv(5)=5

      end if
c                                 Thermodynamic components:
      write (*,2070) 
      write (*,1040) (qname(i),i=1,jcmpn)
      write (*,2021)

      icp = 0 
      char5 = 'b'

      do 

         read (*,3000) char5

         if (char5.ne.blank) then 
c                                 check if it's a saturated phase component,
c                                 mobile or a saturated component:
            if ((ifct.gt.0.or.satflu.or.mobflu).and.
     *          (char5.eq.uname(idh2o).or.char5.eq.uname(idco2))) then 
               call warn (5,r,i,'BUILD')
               write (*,1000) 
               read (*,2100) y
               if (y.eq.'y'.or.y.eq.'Y') then
                  write (*,1070)
               else 
                  write (*,1080) 
                  cycle
               end if 
            end if 

            call chknam (igood,jcmpn,1,good,char5,qname,uname)

            if (.not.good) cycle
c                                  good name, count and save index
            icp = icp + 1
            ic(icp) = igood
            cmpnt(icp) = char5    

         else 
c                                 blank input, check counter
            if (icp.eq.0) then 
               call warn (19,r,k5,'build')
               cycle
            end if

            exit 

         end if 

      end do 
c                                 check for fluid components 
c                                 in thermodynamic composition space.    
      do i = 1, icp
         if (cmpnt(i).eq.uname(idh2o).or.
     *       cmpnt(i).eq.uname(idco2)) then
             call warn (16,r,i,cmpnt(i))
             ifugy = ifugy + 1
         end if
      end do
c                                 
c                                 get fluid equation of state
      if (ifyn.ne.0.or.ifugy.ne.0.or.satflu) call rfluid (1,ifug)
c                                 eliminate composition variable
c                                 for saturated fluid if constrained
c                                 fugacity EoS is used:
c                                 probably should change iv(3)?
      if ((ifug.ge.7.and.ifug.le.9.or.ifug.eq.24).and.ifct.gt.1) 
     *   ivct = ivct - 1
c                                 jcmpn unused components are in qname
c                                 get the indices for chkphi
      do 33 i = 1, icmpn
         icout(i) = 0
         do j = 1, jcmpn
            if (qname(j).eq.uname(i)) goto 33
         end do 
         icout(i) = 1
33    continue 
c                                 also pad out ic array with saturated
c                                 component pointers for chkphi
      do i = 1, isat 
         ic(icp+i) = jsat(i)
      end do 
c                                 ====================================
c                                 next problem class and variable choice
c                                 and ranges
5003  write (*,1490)
      read (*,7010,iostat=ier) icopt
      call rerror (ier,*5003)

      if (icopt.lt.1.or.icopt.gt.5) icopt = 2
c                                 reorder for oned flag
      if (icopt.eq.3.or.icopt.eq.5) then 
         oned = .true.
         if (icopt.eq.5.and.icp.eq.1) call error (53,r,i,'BUILD')
c                                 fractionation from a file
         write (*,1090) 
         read (*,2100) y

         if (y.eq.'y'.or.y.eq.'Y') then 
            write (*,1100) (vname(iv(i)),i=1,ivct)
            write (*,1110)
            write (*,2010) 'coordinate','coor.dat'
            read (*,3040) cfname
            if (cfname.eq.blank) cfname = 'coor.dat'
            open (n8,file=cfname,iostat=ierr,status='old')
            if (ierr.ne.0) then 
               write (*,1140) cfname
            end if 
            close (n8)
            icopt = 10
            oned = .false.
            icont = 1
         end if 
      else 
         oned = .false.
      end if 

      if (icopt.gt.2.and.icopt.lt.10) icopt = icopt - 1
c                                 ====================================
c                                 ask if p = f(T) or vice versa for all
c                                 phase diagram calculations:
c                                
5102  idep = 0
      kvct = ivct 

      if (icopt.ne.3.and.icopt.ne.10) then 

         write (*,1050) vname(1),vname(2)
         if (icopt.eq.4) write (*,1160) 

         read (*,2100) y

         if (y.eq.'Y'.or.y.eq.'y') then
5101        write (*,1060) vname(1),vname(2),vname(2),vname(1)
            read (*,*,iostat=ier) idep
            call rerror (ier,*5101)
            if (idep.lt.1.or.idep.gt.2) goto 5102
c                                 reset the variable counters and flags:
c                                 depend changes the dependent variable
c                                 pointer to iv(ivct), iind points to the
c                                 position of the independent variable in 
c                                 the array v, idep to the dependent v.  
            call depend (ivct,idep,iind,iord,c,dtext)
c                                 kvct is the number of independent potentials
            kvct = ivct - 1

            if (icp.eq.1.and.kvct.eq.1) oned = .true.

         end if

      end if 


      if (oned) then 

         jvct = 1
         icont = 1
c                                 ======================================
c                                 for 1d calculations get the independent
c                                 variable, currently composition is not
c                                 allowed. 
         if (ivct.gt.1) then 
6026        write (*,1210)
            write (*,2140) (j,vname(iv(j)), j = 1, ivct) 
            if (ifct.eq.1.and.ifyn.eq.1) write (*,7150) vname(3) 
            read (*,*,iostat=ier) jc
            call rerror (ier,*6026)
         else 
            jc = 1
         end if 
c                                 get the minimum and maximum values for 
c                                 the path variable.                 
         call redvar (jc,1)
c                                 get sectioning variables values:
         do j = 1, ivct
            if (j.eq.jc) cycle 
            call redvar (j,2)
            vmax(iv(j)) = vmin(iv(j))
         end do     
c                                  put the independent variable in the 
c                                  "x" position (pointer iv(1)).
         ix = iv(1)
         iv(1) = iv(jc)
         iv(jc) = ix
c                                  set the icopt flag to its final value
         if (icopt.eq.2) then
c                                  1-d minimization
            icopt = 5

         else 
c                                  fractionation, also write the grid blurb
            icopt = 7
            write (*,3100) grid(4,1),grid(4,2)

         end if 


      else if (icopt.eq.3) then 
c                                  ==================================-
c                                  swash, not a lot to do
         icopt = 4   


      else if (icopt.eq.2) then
c                                  =====================
c                                  gridded minimization:
         icopt = 5
         jvct = 0
         icont = 1

         if (ivct.eq.1) then
c                                  there is only one potential variable, the
c                                  x variable must be composition, we assume
c                                  the user is not so stupid as to assign a 
c                                  potential dependency when he wants a composition 
c                                  diagram.  
            icont = 2

         else 
c                                  Select the x variable
6017        write (*,2111)
6013        write (*,2140) (j,vname(iv(j)), j = 1, ivct)
            if (icp.gt.1) write (*,1470) j
            if (ifct.eq.1.and.ifyn.eq.1) write (*,7150) vname(3) 
            if (icp.gt.1) write (*,1570)
            read (*,*,iostat=ier) jc
            call rerror (ier,*6013)
 
            if (jc.gt.ivct+1.or.jc.lt.1) then
               write (*,1150)
               goto 6017
            else if (jc.eq.ivct+1) then 
               icont = 2
            else 
               jvct = 1
               ix = iv(1)
               iv(1) = iv(jc)
               iv(jc) = ix
               call redvar (1,1)
            end if
         end if 

         if (ivct.eq.2.and.icont.eq.1) then 
c                                 there is no C variable and there 
c                                 are only 2 potentials, 
c                                 the y variable must be iv(2)     
            call redvar (2,1)
            jvct = ivct   

         else 
c                                 select the y variable 
            if (ivct.gt.1.or.icont.eq.2.and.icp.gt.2) then
               jst = 2
               if (icont.eq.2) jst = 1
6018           write (*,2130)
6015           write (*,2140) (j,vname(iv(j)), j = jst, ivct)
               if (icp.gt.2.and.icont.eq.2) write (*,1480) j
               write (*,*) ' '
               read (*,*,iostat=ier) jc
               call rerror (ier,*6015)
 
               if (jc.gt.ivct+1.or.jc.lt.jst) then
                  write (*,1150)
                  goto 6018
               else if (jc.eq.ivct+1) then
                  icont = 3
               end if

            else if (icont.eq.2) then
 
               jc = 1 

            else 

               jc = 2
 
            end if

            if (icont.lt.3) then 
               ind = 2
               if (icont.eq.2) ind = 1
               ix = iv(ind)
               iv(ind) = iv(jc)
               iv(jc) = ix
               jvct = jvct + 1
               call redvar (ind,1)
            end if  
         end if
c                                 get sectioning variables values:
         do j = jvct+1, ivct
            call redvar (j,2) 
            vmax(iv(j)) = vmin(iv(j))
         end do 

         write (*,3070)
c                                  inform the user of the grid settings:
         do j = 1, 2
            if (j.eq.1) then 
               stext = 'exploratory'
            else
               stext = 'auto-refine'
            end if 

            loopy = (grid(2,j)-1) * 2**(grid(3,j)-1) + 1
            loopx = (grid(1,j)-1) * 2**(grid(3,j)-1) + 1

            write (*,3080) stext,grid(3,j),grid(1,j),grid(2,j),loopx,
     *                     loopy

         end do 

         write (*,3090) opname

      else if (icopt.eq.1) then    
c                                  =========================
c                                  Normal computational mode
6003     write (*,1500)
         if (ivct.gt.1) write (*,1590)
         read (*,7010,iostat=ier) icopt

         call rerror (ier,*6003)

         if (icopt.gt.ivct) goto 6003

         if (icopt.lt.1.or.icopt.gt.2) icopt = 0

         if (icopt.eq.2) then 
            icopt = 1
         else if (icopt.eq.1) then 
            icopt = 3
         end if 

         if (icopt.eq.1) then
c                                  Select the x variable (IV(1)):
5017        write (*,2111)
 
5013        write (*,2140) (j,vname(iv(j)), j = 1, ivct)
            if (ifct.eq.1.and.ifyn.eq.1) write (*,7150) vname(3) 
            read (*,*,iostat=ier) jc
            call rerror (ier,*5013)
 
            if (jc.gt.ivct.or.jc.lt.1) then
               write (*,1150)
               goto 5017
            end if
 
            ix = iv(1)
            iv(1) = iv(jc)
            iv(jc) = ix
 
            call redvar (1,1)
c                                 select the y variable (iv(2)):
            if (ivct.gt.2) then
 
5018           write (*,2130)
5015           write (*,2140) (j,vname(iv(j)), j = 2, ivct)
               read (*,*,iostat=ier) jc
               call rerror (ier,*5015)
 
               if (jc.gt.ivct.or.jc.lt.2) then
                  write (*,1150)
                  goto 5018
               end if
 
            else
 
               jc=2
 
            end if
 
            ix = iv(2)
            iv(2) = iv(jc)
            iv(jc) = ix
 
            call redvar (2,1)
c                                 specify sectioning variables (iv(3)):
            do j = 3, ivct
               call redvar (j,2) 
               vmax(iv(j)) = vmin(iv(j))
            end do 

         else if (icopt.eq.3) then
c                                 select the y variable (iv(1)):
5027        write (*,2210)
 
5023        write (*,2140) (j,vname(iv(j)), j = 1, ivct)
            if (ifct.eq.1.and.ifyn.eq.1) write (*,7150) vname(3) 
            read (*,*,iostat=ier) jc
            call rerror (ier,*5023)
 
            if (jc.gt.ivct.or.jc.lt.1) then
               write (*,1150)
               goto 5027
            end if
 
            ix = iv(1)
            iv(1) = iv(jc)
            iv(jc) = ix
 
            call redvar (1,1)
c                                 specify sectioning variable (iv(2)):
            do j = 2, ivct
               call redvar (j,2) 
               vmax(iv(j)) = vmin(iv(j))
            end do 
         end if
      end if 
c                                 =====================================
c                                 check that X(O) is not 0 or 1
c                                 for fluid speciation routines:
      if (ifyn.eq.1.and.(ifug.ge.10.and.ifug.ne.13.and.ifug.ne.14
     *    .and.ifug.ne.15.and.ifug.ne.18)) then
         if (vmin(3).eq.vmax(3)) then
            if (vmin(3).lt.1d-6) then
               vmin(3) = 1d-6
               vmax(3) = vmin(3)
             else if (vmin(3).gt.0.999999d0) then
               vmin(3) = 0.999999d0
               vmax(3) = vmin(3)
             end if
          else 
             if (vmin(3).lt.1d-6) vmin(3) = 1d-6
             if (vmax(3).gt.0.999999d0) vmax(3) = 0.999999d0
          end if
      end if

      icth = icp + isat
      
      if (icopt.gt.4) then 
c                                 =========================
c                                 Compositional constraints, only
c                                 for constrained minimization
         if (isat.ne.0) then 
            write (*,1450)
            read (*,2100) y 
            if (y.ne.'y'.and.y.ne.'Y') then 
               jcth = icp
            else
               write (*,1460)
            end if 
         else 
            jcth = icp 
         end if 

         do i = 1, icp
            pname(i) = cmpnt(i)
         end do 

         do i = icp + 1, icp + isat
            pname(i) = nname(i-icp)
         end do

         if (jcth.ne.icp)  then 

            jcth = icp

            do i = icp + 1, icth
               write (*,1430) pname(i)
               read (*,2100) y
               if (y.eq.'Y'.or.y.eq.'y') then
                  jcth = jcth + 1
               else 
                  exit
               end if 
            end do 

         end if 
c                                ask if weight amounts, otherwise
c                                use molar amounts.
         write (*,1420) 
         read (*,2100) y
         if (y.eq.'Y'.or.y.eq.'y') then 
            iwt = 1
            amount = 'weight'
         end if 
c                                get the bulk composition:
         if (icont.eq.1) then 
c                                fixed bulk composition
5041        write (*,1390) amount
            write (*,1040) (pname(j), j = 1, jcth)
            write (*,1410) 
            read (*,*,iostat=ier) (dblk(1,j), j = 1, jcth)
            call rerror (ier,*5041)
         else 
c                                user must define compositional variables
            if (lopt(1)) then
c                                closed c-space
               if (icont.eq.2) then
                  write (*,1520) '   C = C0*(1-X(C1)) + C1*X(C1)'
                  write (*,1510) '   C = C0 + C1*X(C1)'
               else 
                  write (*,1540) 
     *           '   C = C0*(1-X(C1)-X(C2)) + C1*X(C1) + C2*X(C2)'
                  write (*,1510) '   C = C0 + C1*X(C1) + C2*X(C2)'
               end if
 
            else 
c                                 open c-space
               if (icont.eq.2) then
                  write (*,1520) '   C = C0 + C1*X(C1)'
                  write (*,1510) '   C = C0*(1-X(C1)) + C1*X(C1)'
               else 
                  write (*,1540) '   C = C0 + C1*X(C1) + C2*X(C2)'
                  write (*,1510) 
     *           '   C = C0*(1-X(C1)-X(C2)) + C1*X(C1) + C2*X(C2)'
               end if

            end if 
 
            do i = 1, icont
6041           write (*,1390) amount
               write (*,1040) (pname(j), j = 1, jcth)
               write (*,1530) nc(i)
               read (*,*,iostat=ier) (dblk(i,j), j = 1, jcth)
               call rerror (ier,*6041)
            end do 
         end if 
      end if 
c                                  ================================
c                                  print and plot file assignment 
c                                  and options:
      write (*,2010) 'print'
      read (*,3040) n3name

      if (n3name.ne.blank) then
         io3 = 0
         write (*,2010) 'plot'
         read (*,3040) n4name
         if (n4name.ne.blank) io4 = 0
      else
c                                if the user doesn't want print
c                                output, ask for plot file name
c                                by default.
         write (*,2011) 'plot','pl'
         read (*,3040) n4name
         if (n4name.eq.blank) n4name = dpl
         io4 = 0
      end if 
 
      if (io4.eq.0) then
c                                  warn about non-plottable results:
         if (icopt.eq.0.and.icp.ne.3) then
            call warn (37,r,i,'BUILD')
         else if (icopt.eq.3.and.icp.ne.2) then
            call warn (38,r,i,'BUILD')
         else if (jcth.gt.0.and.jcth.lt.icp) then
            call warn (39,r,i,'BUILD')
         end if 
      end if
c                                 warn the user he has not requested 
c                                 a plot file.
      if (jcth.ge.icp.and.io4.eq.1) call warn (89,r,k5,'build')
c                                 ======================================
c                                 now start phase data:
c                                 get composition vectors for entities
c                                 defined by a make definition:
      call makecp (inames,mnames,vertex)

      call eohead (n2)

      write (*,*) ' '
c                                 next get consistent real phases
      do 
 
         call getphi (name,eof)

         if (eof) exit
 
         call chkphi (0,name,good)

         if (.not.good) cycle
 
         iphct = iphct + 1

         call loadit (iphct)
 
      end do 
c                                 get all made phases
      do i = 1, nmak
         iphct = iphct + 1
         names(iphct) = mknam(i,mknum(i)+1)
      end do 
c                                 Excluded phases:
      write (*,2080)
      read (*,2100) y
 
      if (y.eq.'Y'.or.y.eq.'y') then
 
         write (*,2034)
         read (*,2100) y
 
         if (y.ne.'y'.and.y.ne.'Y') then
            name = 'b'
            write (*,2021)
 
            do while (ixct.lt.h8.and.name.ne.blank) 

               read (*,2100) name

               if (name.eq.blank) exit 

               good = .false.
 
               do i = 1, iphct
                  if (names(i).eq.name) then 
                     ixct = ixct + 1
                     exname(ixct) = name
                     good = .true.
                     exit 
                  end if 
               end do 

               if (ixct.eq.h8) then 
                  call warn (8,r,h8,'BUILD')
                  exit 
               end if 

               if (good) cycle

               write (*,1020) name
 
            end do
            
         else

            do i = 1, iphct
 
               write (*,1130) names(i)
               read (*,2100) y
 
               if (y.eq.'Y'.or.y.eq.'y') then
                  ixct = ixct + 1
                  if (ixct.gt.h8) then 
                     call warn (8,r,h8,'BUILD')
                     ixct = ixct - 1
                     exit
                  end if 
                  exname(ixct) = names(i)
              end if
 
           end do 
 
        end if
 
      end if
c                                read solution phases.
      write (*,2500)
      read (*,2100) y
 
      if (y.eq.'y'.or.y.eq.'Y') then
c                                 get the file containing the solution models
18       write (*,3010)
         read (*,3040) n9name
         if (n9name.eq.blank) n9name = dsol
         write (*,*) 
         open (n9,file=n9name,iostat=ierr,status='old')
 
         if (ierr.ne.0) then
c                                 system could not find the file
            write (*,3020) n9name
 
5033        write (*,7050)
            read (*,2100,iostat=ier) y
            call rerror (ier,*5033)
 
            if (y.ne.'Y'.and.y.ne.'y') goto 999
            goto 18
c                                 try again
         end if
 
         ict = 0
         ipoint = iphct

c         if (ifugy.eq.0) then
c                                 this was here to prevent users
c                                 from making h2o and/or co2 mobile
c                                 but using the components in a fluid.
c                                 changed as of 5/31/04, JADC.
c            b1 = uname(idh2o)
c            b2 = uname(idco2)
c         else
            b1 = ' '
            b2 = ' '
c         end if 
c                                 test file format
         read (n9,2100) new

         if (new.eq.'new'.or.new.eq.'NEW') then
            inew = 0 
         else if (new.eq.'007'.or.new.eq.'008') then
            inew = 1
         else 
            call warn (4,r,inew,new)
            goto 18
         end if 
c                                 read candidates:
110      call rmodel (blah,bad,vertex)
c                                 istot = 0 = eof
         if (.not.bad.and.istot.ne.0) then 
c                                 don't allow fluid models if 
c                                 the system is fluid saturated:
            if (jsmod.eq.0.and.ifyn.eq.1) goto 110
c                                 check for endmembers:
            call cmodel (im,idsol,blah,1,b1,b2,vertex)
            if (jstot.eq.0) goto 110
      
            ict = ict + 1
            if (ict.gt.i9) call error (24,r,i9,'build')

            tname(ict) = blah 

            goto 110
         end if 
c                                   we have the list, ask user for choices
         if (ict.eq.0) then
 
            write (*,7040)
 
         else
 
            write (*,2510)
            write (*,2520) (tname(i), i = 1, ict)
            write (*,*)

            blah = 'b'

            do 

               read (*,2100) blah
               if (blah.eq.blank) exit
c                                 check if same name entered twice
               do i = 1, isoct
                  if (blah.eq.fname(i)) cycle
               end do 
c                                 check if name in list
               good = .false.

               do i = 1, ict
                  if (blah.eq.tname(i)) then
                     isoct = isoct + 1
                     if (isoct.gt.h9) call error (25,r,h9,'BUILD') 
                     fname(isoct) = blah
                     good = .true.
                     exit 
                  end if 
               end do 

               if (good) cycle

               write (*,2310) blah

            end do  
         end if
      end if 
c                                  get title
      write (*,7070) 
      read (*,2100) title
c                               output options etc:
c                               thermodynamic data file is output by ftopen
c                               first.
      call mertxt (text,n3name,'print file, blank = no output')
      write (n1,3040) text
      call mertxt (text,n4name,'plot file, blank = no output')
      write (n1,3040) text
      call mertxt (text,n9name,
     *            'solution model file, blank = no solutions')
      write (n1,3040) text
c                               if special dependencies, put them in the
c                               title
      if (title.eq.' '.and.idep.ne.0) then 
         text = dtext
      else if (title.eq.' ') then
         text = 'Title Text'
      else if (idep.eq.0) then 
         text = title
      else
         call mertxt (text,title,dtext)
      end if 
c                                 title
      write (n1,3040) text
c                                 option file name
      call mertxt (text,opname,'computational option file')
      write (n1,3040) text 
c                                 computational mode:
      write (n1,1010) icopt,'calculation type: 0 - composition,',
     *                ' 1 - Schreinemakers,', 
     *                ' 3 - Mixed, 4 - gwash,',
     *                ' 5 - gridded min, 7 - 1d fract, 8 - gwash',
     *                ' 9 - 2d fract, 10 - 7 w/file input' 
c                                 coordinate file name if necessary
      if (icopt.eq.10) then 
         call mertxt (text,cfname,'coordinate file')
         write (n1,3040) text
      end if 

      write (n1,1010) idum,'unused place holder, post 06'
      write (n1,1010) idum,'unused place holder, post 06'
      write (n1,1010) idum,'unused place holder, post 06'
      write (n1,1010) idum,'unused place holder, post 06' 
      write (n1,1010) idum,'unused place holder, post 06'
      write (n1,1010) idum,'unused place holder, post 06'
      write (n1,1010) idum,'unused place holder, post 06'
      write (n1,1010) idum,'unused place holder, post 06'
      write (n1,1010) idum,'unused place holder, post 06'
      write (n1,1010) itrans,'number component transformations'
      write (n1,1010) icmpn,'number of components in the data base'
      if (itrans.gt.0) then
         do i = 1, itrans
            write (n1,1120) uname(ictr(i)), ictr(i)
            write (n1,1125) (ctrans(j,i), j = 1, icmpn)
         end do 
      end if 
      write (n1,1010) iwt,'component amounts, 0 - molar, 1 weight'
      write (n1,1010) idum,'unused place holder, post 06'
      write (n1,1010) idum,'unused place holder, post 06'
      write (n1,1010) idum,'unused place holder, post 05'
c                                 output saturated phase eos choice:
      if (ifyn.eq.0.and.ifugy.eq.0) then
         write (n1,1010) 0,'ifug EoS for saturated phase'
      else
         write (n1,1010) ifug,'ifug EoS for saturated phase'
            if (ifug.ge.7.and.ifug.le.12.and.ifug.ne.9.or.
     *          ifug.eq.14.or.ifug.eq.16.or.ifug.eq.17.or.
     *          ifug.eq.24.or.ifug.eq.24.or.ifug.eq.25) 
     *         write (n1,1330) ibuf, hu, dlnfo2, elag,
     *              'ibuf, ipro, choice dependent parameter, ln(ag)'
               if (ibuf.eq.5) write (n1,4020) buf,'a-e'
      end if

      if (oned) then 
         loopx = 1
      else 
         loopx = 2
      end if 

      write (n1,1010) loopx,'gridded minimization dimension (1 or 2)'
c                                 potential variable dependencies
      write (n1,1010) idep,'special dependencies: ',
     *                '0 - P and T independent, 1 - P(T), 2 - T(P)'
      write (n1,1560) c
c                                 output component data:
      write (n1,3060) 'thermodynamic component list'

      do i = 1, icp 
         if (i.gt.jcth) then 
            write (n1,3000) cmpnt(i),0,0.,0.,0.,'unconstrained'
         else 
            write (n1,3000) cmpnt(i),icont,(dblk(j,i),j=1,3),amount
         end if 
      end do 
      write (n1,3050) 'thermodynamic component list'

      write (n1,3060) 'saturated component list'

      do i = 1, isat

         if (i+icp.gt.jcth) then 
            write (n1,3000) nname(i),0,0.,0.,0.,'unconstrained'
         else 
            write (n1,3000) nname(i),icont,(dblk(j,i+icp),j=1,3),amount
         end if

      end do 

      write (n1,3050) 'saturated component list'

      write (n1,3060) 'saturated phase component list'
      do i = 1, ifct
         write (n1,2100) mname(i)
      end do 
      write (n1,3050) 'saturated phase component list'

      write (n1,3060) 
     *'independent potential/fugacity/activity list'
      do i = 1, jmct
         write (n1,2100) oname(i),vname(3+i),afname(i)
      end do 
      write (n1,3050) 'independent potential list'

      write (n1,3060) 'excluded phase list'
      do i = 1, ixct
         write (n1,2100) exname(i)
      end do 
      write (n1,3050) 'excluded phase list'

      write (n1,3060) 'solution phase list'
      do i = 1, isoct
         write (n1,2100) fname(i)
      end do 
      write (n1,3050) 'solution phase list'
c                                 output variable choices and values:
      write (n1,4020) (vmax(i), i= 1, l2),'max p, t, xco2, u1, u2'
      write (n1,4020) (vmin(i), i= 1, l2),'min p, t, xco2, u1, u2'
      write (n1,4020) (dv(i), i= 1, l2),'unused place holder post 06'

      if (oned) then 
         write (n1,1310) (iv(i), i = 1, l2),
     *  'index of the independent & sectioning variables'
      else if (icopt.ne.0) then 
         write (n1,1310) (iv(i), i = 1, l2),
     *  'indices of 1st & 2nd independent & sectioning variables'
      else
         write (n1,1310) (iv(i), i = 1, l2),
     *  'independent variables indices'
      end if 
c                                 get conditions for composition
c                                 diagrams:
      if (icopt.eq.0) then

         if (ifct.eq.1.and.ifyn.eq.1) write (*,7150) vname(3) 

         i = 0

         do 

            i = i + 1
 
            write (*,6020) (vname(iv(j)),j=1,ivct)
5036        write (*,6010) i
            read (*,*,iostat=ier) (vmin(iv(j)),j=1,ivct)
            call rerror (ier,*5036)
            if (vmin(iv(1))+vmin(iv(2)).eq.0d0) goto 99 
            write (n1,1340) vmin

         end do 
      end if
 
99    endfile (n1)
 
      close (n1)
 
999   stop

1000  format ('Are you sure you want to do this (y/n)?')
1010  format (i5,1x,a,a,a,a,a)
1020  format (a,' does not exist in the selected data base')
1030  format (/,'For C-O-H fluids it is only necessary to select ',
     *       'volatile species present in',/,'the solids of interest. ',
     *       'If the species listed here are H2O and CO2, then to',/,
     *       'constrain O2 chemical potential to be consistent with ',
     *       'C-O-H fluid speciation',/,'treat O2 as a saturated ',
     *       'component. Refer to the Perple_X Tutorial for details.',/)
1040  format (12(1x,a5))
1050  format (/,'The data base has ',a,' and ',a,' as default ',
     *       'independent potentials.',/,'Make one dependent on the ',
     *       'other, e.g., as along a geothermal gradient (y/n)? ')
1060  format (/,'Select dependent variable:',//,'  1 - ',a,' = f(',a,')'
     *       ,/,'  2 - ',a,' = f(',a,')',/)
1070  format ('ok, but dont say i didnt warn you.')
1080  format ('wise move, choose another component.')
1090  format (/,'Enter path coordinates from a file (Y/N)?')
1100  format (/,'In this mode VERTEX/WERAMI read path coordinates',
     *        'from a file',/,'the file must be formatted so that',
     *        ' the coordinates of each point',/,'are on a separate',
     *        ' line, the coordinates are, in order:',
     *        4x,5(a,2x))
1110  format (/)
1120  format (a,1x,i2,' component transformation')
1125  format (13(f6.2,1x))
1130  format ('Exclude ',a,' (Y/N)? ')
1140  format (/,'File: ',a,/,'Does not exist, you must create it',
     *        ' before running VERTEX.',/)
1150  format (/,'huh?',/)
1160  format (/,'Answer yes to specify a P-T path for phase ',
     *          'fractionation calculations.',/)
1170  format (/,'Enter the computational option file name',
     *       ', left justified, [default = perplex_option.dat]:',/,
     *       'See: www.perplex.ethz.ch/perplex_options.html')
1210  format ('Select the path variable for the calculation:',/)
1310  format (/,5(i2,1x),2x,a,/)
1330  format (i2,1x,i2,1x,g13.6,1x,g13.6,1x,a)
1340  format (5(g13.6,1x))
1390  format (/,'Enter ',a,' amounts of the components:')
1410  format ('for the bulk composition of interest:')
1420  format (/,'Specify component amounts by weight (Y/N)?')
1430  format ('Constrain component ',a,' (Y/N)?')
1450  format (/,'In this mode all thermodynamic components must be ',
     *         'constrained,',/'constrain ',
     *         'saturated components also (Y/N)?')
1460  format (/,'The next prompts are for the saturated',
     *        ' component constraints.',/,
     *        'Answering no at any point',
     *        ' completes the set of constraints.',/)
1470  format (5x,i1,' - Composition X(C1)* (user defined)')
1480  format (5x,i1,' - Composition X(C2) (user defined)')
1490  format (//,'Specify computational mode:',//,
     *    5x,'1 - Unconstrained minimization',/, 
     *    5x,'2 - Constrained minimization on a 2d grid [default]',/,
     *    5x,'3 - Constrained minimization on a 1d grid',/,
     *    5x,'4 - Output pseudocompound data',/,
     *    5x,'5 - Phase fractionation calculations',//,
     *        'Unconstrained optimization should be used for the ',
     *        'calculation of composition,',/,'mixed variable, and ',
     *        'Schreinemakers diagrams. Gridded minimization can be ',
     *        'used to',/,'construct phase diagram sections.',/)
1500  format (//,'Specify number of independent potential variables:',
     *         /,5x,'0 - Composition diagram [default]',/,
     *           5x,'1 - Mixed-variable diagram')
1510  format (/,'To compute bulk compositions as:',a,/,'change the ',
     *       'compuational option keyword closed_c_space.')
1520  format (/,'The bulk composition of the system will be computed ',
     *       'as:',/,a,/,
     *       'where X(C1) varies between 0 and 1, and C0 and C1 are ',
     *       'the compositions',/,'specified next.')
1530  format ('to define the composition ',a)
1540  format (/,'The bulk composition of the system will be computed ',
     *       'as:',/,a,/,
     *       'where X(C1) and X(C2) vary between 0 and 1, and C0, C1 ',
     *       'and C2 are the',/,'compositions specified next.')
1560  format (5(g12.6,1x),'Geothermal gradient polynomial coeffs.')
1570  format (/,'*X(C1) can not be selected as the y-axis variable',/)
1590  format (5x,'2 - Sections and Schreinemakers-type diagrams')
2000  format (/,'Do you want a ',a,' file (Y/N)?')
2010  format ('Enter the ',a,' file name, < 100 characters, left ',
     *        'justified (blank for none):')
2011  format ('Enter the ',a,' file name, < 100 characters, left ',
     *        'justified [default = ',a,'] (blank for default):')
2021  format ('Enter names, left justified, 1 per line,',
     *        ' and a blank to finish:')
2030  format (/,'Calculations with a saturated phase',
     *        ' (Y/N)?',/,'The phase is: ',a,/,
     *        'Its compositional variable is: Y(CO2), X(O), etc.')
2031  format (/,'Select the independent saturated phase components:')
2032  format (/,'Select < ',i1,' saturated components from the set:')
2034  format ('Do you want to be prompted for phases (Y/N)? ')
2040  format (/,'Use chemical potentials, activities or fugacities as',
     *        ' independent',/,'variables (Y/N)?')
2050  format (/,'Specify a component whose chemical potential, activi',
     *       'ty or fugacity is',/,'to be independent, ',
     *       'enter a blank to finish:')
2060  format (/,'Component ',a,' is to be characterized by:',//,
     *       5x,'1 - chemical potential [default]',/,
     *       5x,'2 - log10(fugacity)',/,
     *       5x,'3 - log10(activity)')
2061  format (/,'The data file contains no phase suitable to define ',
     *        'the',a,' of ',a,/,'This component will be ',
     *        'characterized by its chemical potential',/)
2062  format (/,'Phase ',a,' will be used to define ',a,' ',a)
2063  format (/,'Select the phase to be used to define ',a,' ',
     *        a,' from the following:',/)
2064  format (4(5x,9(a,1x),/))
2065  format (/,'The log10(',a,' ',a,') variable is named: ',a,/)
2066  format (/,'The chemical potential of ',a,' is named: ',a,/)
2070  format (/,'Select thermodynamic components from the set:')
2080  format (/,'Exclude phases (Y/N)?')
2100  format (3(a,1x))
2110  format (/,'Calculations with saturated components (Y/N)?')
2111  format (/,'Select x-axis variable:')
2130  format (/,'Select y-axis variable:')
2140  format (5x,I1,' - ',a)
2210  format (/,'Select vertical axis variable:')
2310  format (/,a,' is invalid. Check spelling, upper/lower',
     *        ' case matches. Try again:',/)
2500  format (/,'Do you want to treat solution phases (Y/N)? ')
2510  format (/,'Select phases from the following list, enter',
     *        ' 1 per line,',/,'left justified, enter a blank',
     *        ' to finish',/)
2520  format (6(2x,a))
3000  format (a,1x,i1,1x,3(g12.6,1x),a,' amount')
3010  format ('Enter solution model file name [default = ',
     *        'solution_models.dat]',/,'left justified, ',
     *        '< 100 characters: ')
3020  format (/,'**error ver191** FOPEN cannot find file:',/,a,/)
3040  format (a,10x,a)
3050  format ('end ',a,/)
3060  format (/,'begin ',a)
3070  format (/,'For gridded minimization, grid resolution is ',
     *        'determined by the number of levels',/,
     *        '(grid_levels) and the resolution at the lowest level ',
     *        'in the X- and Y-directions (x_nodes',/,
     *        'and y_nodes) these parameters are currently set for ',
     *        'the exploratory and autorefine cycles',/,
     *        'as follows:',//,'stage        grid_levels  xnodes  ',
     *        'ynodes    effective resolution')
3080  format (a,7x,i1,8x,i4,4x,i4,6x,i4,' x',i4,' nodes')
3090  format (/,'To change these options edit or create',
     *        ' the file ',a,/,'See: ',
     *        'www.perplex.ethz.ch/perplex_options.html#grid_parameters'
     *        ,/)
3100  format (/,'For phase fractionation calculations the number of ',
     *        'points computed along the path',/,
     *        'is determined by the 1d_path parameter. The values ',
     *        'for this parameter are currently',/,
     *        'set to ',i3,' and ',i3,' points for the exploratory ',
     *        'and autorefine cycles.')
4020  format (2(g11.5,1x),f10.8,1x,2(g11.5,1x),a)
5000  format (a,'_',a)
6020  format (/,'Specify values for:',/,(10x,5(a,2x)))
6010  format ('For calculation ',i2,', enter zeros to finish.')
7010  format (bn,i80)
7020  format (//,'NO is the default (blank) answer to all Y/N prompts',
     *        /)
7040  format (/,'The solution model file contains no',
     *               ' models valid for your problem.',/)
7050  format ('Try again (Y/N)? ')
7070  format (/,'Enter calculation title: ')
7150  format (/,'*Although only 1 component is specified for the ',
     *       'saturated fluid phase, the',/,' equation of state ',
     *       'permits fluid composition to vary through the variable ',
     *       a,'.',/,' For pure fluids, specify the appropriate ',
     *       'composition as a sectioning constraint.',/) 
      end
 
      subroutine ftopen
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'
 
      character*100 n2name, y*1, ddata*14, din*14, blank*1, text*140

      integer ierr

      data ddata, din, blank/'hp02ver.dat   ','in            ',' '/
c-----------------------------------------------------------------------
c                                 now the input data file to
c                                 be created for vertex.
10    write (*,1040)
      read (*,1020) n2name
      if (n2name.eq.blank) n2name = din
      open (n1,file=n2name,iostat=ierr,status='new')
      if (ierr.ne.0) then
         write (*,1050) n2name
         read (*,1020) y
         if (y.eq.'Y'.or.y.eq.'y') goto 10
         open (n1,file=n2name)
      end if
c                                 next the thermo data file
1     write (*,1000)
      read (*,1020) n2name
      if (n2name.eq.blank) n2name = ddata
      open (n2,file=n2name,iostat=ierr,status='old')
      if (ierr.ne.0) then
c                                 system could not find the file
         write (*,1010) n2name
         read (*,1020) y
         if (y.ne.'Y'.and.y.ne.'y') goto 999
         goto 1
c                                 try again
      end if
 
      call mertxt (text,n2name,'thermodynamic data file')
      write (n1,1020) text
 
      return
 
999   write (*,1060)
 
      stop
 
1000  format (/,'Enter thermodynamic data file name',
     *          ', left justified, [default = hp02ver.dat]:')
1010  format (/,'**warning ver191** FOPEN cannot find file ',/,a
     *         ,//,'try again (y/n)?')
1020  format (a,10x,a)
1040  format (/,'Enter the name of the problem definition file to be ',
     *       'creat',
     *       'ed,'/,'< 100 characters, left justified [default = in]:')
1050  format (/,'File:',/,a,/,'exists, do not overwrite it (y/n)? ')
1060  format (/,'O.K., I quit too.')
 
      end
 
      subroutine grxn (g)
c--------------------------------------------------------------------
c a dummy routine to allow rk to be linked with rlib.f
c--------------------------------------------------------------------
      double precision g
      g = g
      end

      subroutine depend (ivct,idep,iind,iord,c,dtext)
c---------------------------------------------------------------------------
c subroutine to reset variable flags and counters if one variable is
c made dependent on another; the routine also prompts and reads the 
c functional dependence
c---------------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      character xname*8, dtext*200, vname*8

      integer i, j1, j2, ivct, iind, idep, iord, ier

      double precision c(0:4)
 
      common/ csta2 /xname(k5),vname(l2)

      integer ipot,jv,iv
      common/ cst24 /ipot,jv(l2),iv(l2)
c                                 reset the variable counters and flags:
      ivct = ivct - 1

      if (idep.eq.1) then 
         iind = 2
         j1 = 1
         j2 = 2
         do i = 1, ivct 
            iv(i) = iv(i+1)
         end do 
      else 
         iind = 1
         j1 = 2
         j2 = 1
         do i = 2, ivct 
            iv(i) = iv(i+1)
         end do 
      end if 

      iv(ivct+1) = idep
c                                  now get the functional dependence
5103  write (*,1070) vname(j1),vname(j2)
      read (*,*,iostat=ier) iord
      call rerror (ier,*5103)
      if (iord.lt.0.or.iord.gt.5) goto 5103
      do i = 0, iord
5104     write (*,1080) i
         read (*,*,iostat=ier) c(i)
         call rerror (ier,*5104)
      end do 

      write (*,'(/)')
c                                  write a text version for the title
      write (dtext,1580) vname(idep),c(0),
     *                      (c(i),vname(iind),i,i=1,iord)
      call deblnk (dtext)

1080  format ('Enter c(',i2,')')
1070  format (/,'The dependence must be described by the polynomial',//,
     *        a,' = Sum ( c(i) * [',a,']^i, i = 0..n)',//,
     *       'Enter n (<5)')
1580  format (a,' = ',g12.6,4(' + ',g12.6,' * ',a,'^',i1))

      end 

      subroutine redvar (ind,iprompt)
c----------------------------------------------------------------------
c redvar interactively reads and checks values for the primary variable 
c indexed by ind. if num = 1 it reads only the vmin(ind) value, else both
c vmin(ind) and vmax(ind) are read.
c----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer ind, iprompt, icount, ier

      logical numbad

      double precision vmax,vmin,dv
      common/ cst9  /vmax(l2),vmin(l2),dv(l2)

      character*8 vname, xname
      common/ csta2 /xname(k5),vname(l2)

      integer ipot,jv,iv
      common/ cst24 /ipot,jv(l2),iv(l2)
c----------------------------------------------------------------------

      do 

         if (iprompt.eq.1) then 
            icount = 2
            write (*,1010) vname(iv(ind))
         else if (iprompt.eq.2) then
            icount = 1
            write (*,1020) vname(iv(ind))
         end if 

         if (icount.eq.1) then 
            read (*,*,iostat=ier) vmin(iv(ind))
         else 
            read (*,*,iostat=ier) vmin(iv(ind)),vmax(iv(ind))
         end if 

         if (ier.ne.0) then 
            write (*,1000)
            cycle 
         end if 

         if (numbad(1,ind)) cycle 

         if (icount.eq.2) then 
            if (numbad(2,ind)) cycle  
         end if 

         exit 
    
      end do 

1000  format (/,' Your input is incorrect, probably you are using ',
     *        'a character where',/,' you should be using a number ',
     *        'or vice versa, try again...',/)
1010  format (/,'Enter minimum and maximum values, respectively,',
     *        ' for: ',a)
1020  format (/,'Specify sectioning value for: ',a)

      end 

      logical function numbad (num,ind)
c----------------------------------------------------------------------
c numbad checks if a primary variable limit is reasonable, the variable
c is indexed by ind and is the lower (vmin) or upper (vmax) bound 
c depending on the value of num (1 or 2, respectively). 
c----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      character*1 y

      integer num, ind, jnd

      double precision value 

      integer imaf,idaf
      common/ cst33 /imaf(i6),idaf(i6)

      double precision vmax,vmin,dv
      common/ cst9  /vmax(l2),vmin(l2),dv(l2)

      character*8 vname, xname
      common/ csta2 /xname(k5),vname(l2)

      integer ipot,jv,iv
      common/ cst24 /ipot,jv(l2),iv(l2)
c----------------------------------------------------------------------

      numbad = .false.

      jnd = iv(ind) 

      if (num.eq.1) then 
         value = vmin(jnd)
      else
         value = vmax(jnd)
      end if 

      if (jnd.eq.1.or.jnd.eq.2) then 
c                                 pressure-temperauture
         if (value.le.0d0) then 
            write (*,1020) vname(jnd)
            numbad = .true.
         end if 

      else if (jnd.eq.3) then 
c                                 phase composition
         if (value.lt.0d0.or.value.gt.1d0) then
            write (*,1030) vname(jnd)
            numbad = .true.
         end if 
 
      else 
c                                 potential
         if (imaf(jnd-3).gt.1.and.value.gt.0d0) then
            write (*,1000) vname(jnd),vname(jnd)
            numbad = .true.
         end if 

      end if 

      if (numbad) then 
         read (*,1010) y
         if (y.ne.'y'.and.y.ne.'Y') numbad = .false.
      end if  

1000  format (/,'**WARNING** ',a,' is the base 10 log of activity/',
     *       'fugacity; therefore ',a,'> 0',/,'implies ',
     *       'supersaturation with respect to the reference phase.',/,
     *       'Specify new values (Y/N)?',/)
1010  format (a)
1020  format (/,'**WARNING** negative or zero values for: ',a,
     *       ' may destabilize',/,
     *       'calculations with some equations of state.',/,
     *       'Specify new values (Y/N)?',/)
1030  format (/,'**WARNING** ',a,'should be in the range 0-1.',/
     *       'Specify new values (Y/N)?',/)
      end 

      subroutine chknam (igood,jcmpn,iflu,good,char5,qname,uname)
c----------------------------------------------------------------------
c chknam - looks for a match between string char5 and the jcmpn strings in 
c array qname. if the match is found goof = .true, index is the index of
c the string in uname, and the string is eliminated from uname and jcmpn
c decremented. igood is the index of the equivalent string in the uname
c array.

c if iflu = 0, then chknam first matches fluid components in uname
c before eliminating the component from qname.
c----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      character*5 qname(k0), uname(k0), char5

      logical good

      integer igood, jcmpn, iflu, i, j

      integer idh2o,idco2,ikind,icmpn,icout
      double precision therm,comp,tot
      common/ cst43 /therm(k4),comp(k0),tot,icout(k0),idh2o,idco2,
     *               ikind,icmpn
c----------------------------------------------------------------------

      good = .false.

      if (iflu.eq.0.and.
     *    char5.ne.uname(idh2o).and.
     *    char5.ne.uname(idco2)) then
c                                 special check for saturated phase
c                                 components. 
         write (*,1000) char5
         return
      end if 

      do i = 1, jcmpn
         if (qname(i).eq.char5) then 
c                                  eliminate used components
            do j = i+1, jcmpn
                qname(j-1) = qname(j)
            end do 
            good = .true.
            exit 
         end if 
      end do

      if (good) then 
c                                 decrement qname counter
         jcmpn = jcmpn - 1
c                                 find index in uname array
         do i = 1, icmpn
            if (char5.ne.uname(i)) cycle 
            igood = i
            exit
         end do 

      else 

         write (*,1000) char5

      end if 

1000  format (/,a,' is invalid. Check spelling, upper/lower',
     *        ' case matches. Try again:',/)
      end 