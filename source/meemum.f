      program meemm       
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,ier,idead

      logical meemum, nodata, bulk

      character amount*6, yes*1

      integer itri(4),jtri(4),ijpt

      double precision ctot, wt(3) 

      integer iwt
      common/ cst209 /iwt

      double precision atwt
      common/ cst45 /atwt(k0)

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      double precision vmax,vmin,dv
      common/ cst9  /vmax(l2),vmin(l2),dv(l2)

      integer ipot,jv,iv
      common / cst24 /ipot,jv(l2),iv(l2)

      character*8 vname,xname
      common/ csta2  /xname(k5),vname(l2)

      character*5 cname
      common/ csta4 /cname(k5)

      integer jbulk
      double precision cblk
      common/ cst300 /cblk(k5),jbulk

      double precision a,b
      common/ cst313 /a(k5,k1),b(k5)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp
c                                 meemum flag used for output routines
c                                 to distunguish meemum from werami
      save meemum
      data meemum/.true./
c----------------------------------------------------------------------
c                                 initialization, read files etc. 
      call iniprp

      write (*,1000) 
      read (*,1050) yes
      if (yes.eq.'y'.or.yes.eq.'Y') then 
c                                 bulk is true, user enters composition and p-t conditions
         bulk = .true.
      else 
c                                 else user enters only p-t and composition read from input file.
         bulk = .false.
      end if 

c                                 iwt is set by input, it is only used below to determine
c                                 whether to convert weight compositions to molar. the 
c                                 computations are done solely in molar units. 
      amount = 'molar '
      if (iwt.eq.1) amount = 'weight'
c                                 computational loop
      do 
c                                 read potential variable values    
c                                 v(1) is P(bar), v(2) is T(K) the pointer jv used 
c                                 for general problems but can be eliminated for calculations 
c                                 simply as a f(P,T)       
         write (*,1070) (vname(jv(i)), i = 1, ipot)
         read (*,*,iostat=ier) (v(jv(i)), i = 1, ipot)
         if (ier.ne.0) cycle
         if (v(jv(1)).eq.0d0) exit 
          
         if (bulk) then 
c                                 load the composition into b, the component names are  
c                                 in cname, if iwt = 1 the composition is in mass fractions
c                                 otherwise in molar units. 
10          write (*,1060) amount
            write (*,1050) (cname(i),i=1,jbulk)
            read (*,*,iostat=ier) (cblk(i),i=1,jbulk)
            if (ier.ne.0) goto 10 
         
            if (iwt.eq.1) then 
c                                 convert mass to molar 
               do i = 1, jbulk
                  cblk(i) = cblk(i)/atwt(i)
               end do 
            end if
c                                 normalize the composition vector, this 
c                                 is necessary for reasons of stupidity (lpopt0). 
            ctot = 0d0

            do i = 1, icp
               ctot = ctot + cblk(i)
            end do 

            do i = 1, icp
               b(i) = cblk(i)/ctot
            end do

         end if 
c                                 lpopt does the minimization and outputs
c                                 the results to the print file.
         call lpopt0 (idead)

         if (idead.gt.0) then
 
            write (*,*) 'minimization failed'

         else 
c                                 compute derivative properties
            call getloc (itri,jtri,ijpt,wt,nodata,meemum)
c                                 print summary to LUN 6
            call calpr0 (6,meemum)

         end if 

      end do

1000  format (/,'Interactively enter bulk compositions (y/n)?',/,
     *          'If you answer no, MEEMUM uses the bulk composition',
     *         ' specified in the input file.',/)
1060  format (/,'Enter ',a,' amounts of the components:')
1050  format (12(a,1x))
1070  format (/,'Enter (zeroes to quit) ',7(a,1x))

      end 
    
      subroutine iniprp
c----------------------------------------------------------------------
c iniprp - read data files and initialization for mingee
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical first, output, vertex 

      character n4name*100

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      integer iemod,kmod
      logical smod
      double precision emod
      common/ cst319 /emod(k15,k10),smod(h9),iemod(k10),kmod
c----------------------------------------------------------------------- 
      first = .true.
      output = .false.
      vertex = .true.
c                                 elastic modulii flag
      kmod = 0 
c                                 -------------------------------------------
c                                 open statements for units n1-n5 and n9
c                                 are in subroutine input1
      call input1 (first,output,n4name)
c                                 for meemum turn auto_refine OFF
      iopt(6) = 0 
c                                 read thermodynamic data on unit n2:
      call input2 (first)
c                                 read data for solution phases on n9:
      call input9 (vertex,first,output)
c                                 call initlp to initialize arrays 
c                                 for optimization.
      call initlp     

      end
