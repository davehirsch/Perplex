
c   ctransf is a program to read a vertex thermo-data file and
c   rewrite the data in a new file with transformed components.  
c-----------------------------------------------------------------------
      write (*,1000)
c                               assign data files
      call ftopen
c                               Read THERMODYNAMIC DATA file (N2):
c                               read the data base header
      call topn2 (5)
c                               read and echo data cards with
c                               component conversion
35    call getph (*99)
      goto 35
      
1000  format (//,'NO is the default answer to all Y/N prompts',/)

99    end
 
      subroutine ftopen
c-----------------------------------------------------------------------
      implicit none



      include 'perplex_parameters.h'



      integer ierr
 
      character*100 n2name,yes*1
c-----------------------------------------------------------------------
c                                 first the thermo data file
1     write (*,1000)
      read (*,1020) n2name
      open (n2,file=n2name,iostat=ierr,status='old')
      if (ierr.ne.0) then
c                                 system could not find the file
         write (*,1010) n2name
         read (*,1020) yes
         if (yes.ne.'Y'.and.yes.ne.'y') goto 999
         goto 1
c                                 try again
      end if
 
      write (*,1070)

      open (n8,file='ctransf.dat')
 
      return
 
999   write (*,1060)
 
      stop
 
1000  format (/,'Enter thermodynamic data file name (e.g.',
     *         ' hp02ver.dat), left justified:')
1010  format (/,'**warning ver191** FOPEN cannot find file:',/,a
     *         ,//,'Try again (y/n)?')
1020  format (a)
1060  format (/,'O.K., then i quit too.',/)
1070  format (/,'Output will be written to file: ctransf.dat',/)
 
      end
  
      subroutine getph (*)
c------------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'
 
      character*8 name, oname, note*40, record*1



      integer ibase, i, it, j, jlam



      double precision ct, emod(k15)
 

      integer ilam,idiso,lamin,idsin

      double precision tm,td

      common/ cst202 /tm(m7,m6),td(m8),ilam,idiso,lamin,idsin



      integer ictr, itrans

      double precision ctrans

      common/ cst207 /ctrans(k0,k5),ictr(k5),itrans



      integer idh2o,idco2,ikind,icmpn,icout

      double precision therm,comp,tot

      common/ cst43 /therm(k4),comp(k0),tot,icout(k0),idh2o,idco2,

     *               ikind,icmpn


      save oname
 
      data oname/' '/
 
30    read (n2,1010,end=90,err=98) record

      if (record.eq.' ') then

c                          check for comments, i.e., data 

c                          with a blank 1st character 

         goto 30

      else  

         backspace (n2)

         read (n2,1020,end=90,err=98)

     *                     name, ibase, ikind, ilam, idiso, note

      end if


      write (n8,1020) name, ibase, ikind, ilam, idiso, note

      read (n2,*,err=98) (comp(i), i=1, icmpn), 

     *                   (therm(i), i=1, k14)
c                               do component transformation if
c                               itrans is not zero
      do i = 1, itrans
         it = ictr(i)
         if (comp(it).eq.0.d0) cycle
c                                ct is how much of the new
c                                component is in the phase.
         ct =  comp(it) / ctrans(it,i)
 
         do j = 1, icmpn
             comp(j) = comp(j) - ct * ctrans(j,i)
         end do 

 
         comp(it) = ct


      end do 

      write (n8,1040) (comp(i), i= 1, icmpn)
      write (n8,1050) (therm(i), i= 1, k14)

      if (ilam.ne.0) then
c                               determine number of transitions from
c                               flag ilam:
         jlam=ilam
         if (ilam.gt.3) jlam = ilam-3
         if (ilam.gt.6) jlam = ilam-6
         if (ilam.gt.9) jlam = ilam-9
 
         do i= 1, jlam
            read (n2,*,err=98) (tm(j,i), j = 1, m7-2)
            write (n8,1050) (tm(j,i), j = 1, m7-2)
         end do 

      end if
 
      if (idiso.ne.0) then 
         read (n2,*,err=98) td
         write (n8,1050) td
      end if 



      if (ikind.ne.0) then 

         read (n2,*,err=98) emod

         write (n8,1050) emod

      end if 
  
      oname = name
 
      return
 
90    return 1

98    if (oname.ne.' ') then
         call error (23,ct,i,oname)
      else
         call error (23,ct,i,'  none  ')
      end if


1010  format (a1)
1020  format (a,i2,i2,i2,i2,1x,a)
1040  format (12(f5.2,1X))
1050  format (5(g13.7,1X))
 
      end



      subroutine grxn (g)

c--------------------------------------------------------------------

c a dummy routine to allow rk to be linked with rlib.f

c--------------------------------------------------------------------

      implicit none

      double precision g

      g = g

      end