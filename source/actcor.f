      program actcor 
c----------------------------------------------------------------------
c                       ************************
c                       *                      *
c                       *  actcor.may.1989     *
c                       *                      *
c                       ************************
c----------------------------------------------------------------------
c a fortran program for making fixed activity corrections to the
c thermodynamic data file for vertex.  actcor creates a new data file 
c with the corrected data on unit n8
c-----------------------------------------------------------------------
c files (see vertex program documentation for additional information): 
c-----------------------------------------------------------------------
      implicit none

 

      include 'perplex_parameters.h'
                                                         
      character*8 blank8, name, y*1, test

      logical eof
           
      data blank8/'        '/
c-----------------------------------------------------------------------   
      write (*,1300) 
c                             open files
      call ftopen
c                             read and echo file header
      call topn2 (4)
      write (*,1010) 
c                             allow user to enter names:
      write (*,1030) 
      read (*,1000) y

      if (y.ne.'y'.and.y.ne.'Y') then 
c                             get the name:
100      write (*,1020) 
         read (*,1000) test

         if (test.eq.blank8) goto 999

         rewind n2
         call eohead (n2)

         do 
 
            call getphi (name,eof)
 
            if (eof) then 
               write (*,1050) test
               goto 100
            end if 

            if (name.eq.test) then
               call gotcha (name)
               goto 100
            end if 

         end do 

      else 
c                             read and modify individual entries  
         do 
 
            call getphi (name,eof)
 
            if (eof) exit

            write (*,1040) name
            read (*,1000) y
            if (y.eq.'y'.or.y.eq.'Y') call gotcha (name)   

         end do 
                     
      end if

1000  format (a) 
1010  format ('This program will create a new thermodynamic data',/,
     *        'file with (optionally) activity corrected entries.',/,
     *        'You must specify all phases that are to be included',/,
     *        'in the new data file (actcor.dat).',//)
1020  format ('Enter a phase to be included [<9 characters, blank to ',
     *        'finish]:')
1030  format ('Prompt for phases (y/n)?')
1040  format ('Include (y/n): ',a)
1050  format ('No such phase as: ',a)
1300  format (/,'NO is the default answer to all prompts',/)           

999   end 

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

      open (n8,file='actcor.dat')
 
      return
 
999   write (*,1060)
 
      stop
 
1000  format (/,'Enter the thermo data file name (e.g.',
     *         ' hp90ver.dat) < 140 characters, left justified: ')
1010  format (/,'**warning ver191** FOPEN cannot find file:',/,a
     *         ,//,'try again (y/n)? ')
1020  format (a)
1060  format (/,'O.K., then i am quitting also.',/)
1070  format (/,'Output will be written to file: actcor.dat',/)

      end

      subroutine grxn (g)

      implicit none
      double precision g
      end

      subroutine gotcha (name)

      implicit none


      include 'perplex_parameters.h'
                                                         
      character*8 blank8,name,y*1



      integer i,j,jlam



      double precision xmole,xmix,act



      character*5 cmpnt, dname*40

      common/ csta5 /dname,cmpnt(k0)

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



      double precision p,t,xco2,u1,u2,tr,pr,r,ps

      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps
                                                            
      data blank8/'        '/

      write (*,1060) name
      read (*,1000) y  
                       
      if (y.eq.'y'.or.y.eq.'Y') then 
c 
         write (*,1070) name
         read (*,1020) blank8
         write (*,1080) name
         write (*,2000) (cmpnt(i),i=1,icmpn)
         write (*,1040) (comp(i),i=1,icmpn) 
         write (*,1090)

         read (*,1000) y    

         if (y.eq.'y'.or.y.eq.'Y') then 
            write (*,1100) name,blank8

            read (*,*) xmole
               write (*,1110) name
               read (*,*) xmix

               act = xmole**xmix   

            else   
               write (*,1120) name
               read (*,*) act
            end if 

            write (*,1130) name,blank8,act 
            therm(1) = therm(1) + t * 8.314413 * dlog(act)
            therm(2) = therm(2) - 8.314413 * dlog(act)
            name = blank8

         end if 

         jlam = ilam

         if (ilam.gt.3) jlam = ilam - 3
         if (ilam.gt.6) jlam = ilam - 6
         if (ilam.gt.9) jlam = ilam - 9

         write (n8,1030) name,1,ikind,ilam,idiso
         write (n8,1040) (comp(i),i=1,icmpn) 
         write (n8,1050) (therm(i),i=1,18) 
         if (ilam.ne.0) then
            do 52 i = 1, jlam
52             write (n8,1050) (tm(j,i), j = 1, m7-2) 
         end if 
         if (idiso.ne.0) write (n8,1050) td

1000  format (a)
1020  format (a8,i2,i2,i2,i2)
1030  format (a8,i2,i2,i2,i2,a60)
1040  format (13(f5.2,1x))
1050  format (5(g13.7,1x))          
1060  format ('make an activity correction for ',a,' (y/n)?')
1070  format ('enter phase that ',a,' is a solution in:')
1080  format ('the stoichiometry of ',a,' is:')
1090  format (/,'ideal activity model (y/n)?')
1100  format ('enter mole fraction (x) of ',a,' in ',a,':')
1110  format ('activity of ',a,' will be computed as x**n',/,
     *        'enter number of mixing sites (n):')
1120  format ('enter activity of ',a,':')
1130  format (/,'activity of ',a,' in ',a,' is: ',g12.6)
2000  format (/,1x,13(a,1x),/,1x,13(a,1x))

      end

