      PROGRAM REWRIT 
C----------------------------------------------------------------------
C                       ************************
C                       *                      *
C                       *  REWRI.MAY.26.1997   *
C                       *                      *
C                       ************************
C----------------------------------------------------------------------
C REWRIT a FORTRAN program for rewriting or reformatting the N2 data
C file for the VERTEX program.  BUILD reads data from N2.  
C The output file is written to unit N1
C-----------------------------------------------------------------------
C FILES (see VERTEX program documentation for additional information): 
C-----------------------------------------------------------------------

 
      IMPLICIT none

      include 'perplex_parameters.h'
C                                                         
      CHARACTER*8 FNAME*100,CNAME*5,RECORD*132,
     *            DNAME*40,CMPNT*5,VNAME*8,XNAME*18,
     *            SIXTY*60,NAME 

      integer ilam,idiso,zi,zj,jlam,icomp,istct,iphct,icp,isoct,lamin,
     *        idsin,io3,io4,io9,ibase
      
      double precision td,v,pr,tr,r,ps,vmax,vmin,dv,tm
      
      DIMENSION FNAME(2)

      COMMON/ CST6 /ICOMP,ISTCT,IPHCT,ICP/ CST79 /ISOCT 
     *      / csta2 /xname(k5),vname(l2)
     *      / cst202 /tm(m7,m6),td(m8),ilam,idiso,lamin,idsin
     *      / cst5  /v(l2),tr,pr,r,ps/ cst9  /vmax(l2),vmin(l2),dv(l2)
     *      / CST41 /io3,io4,io9
     *      / CSTA4 /CNAME(k5)
     *      / CSTA5 /DNAME,CMPNT(k0) 

      integer idh2o,idco2,ikind,icmpn,icout
      double precision therm,comp,tot
      common/ cst43 /therm(k4),comp(k0),tot,icout(k0),idh2o,idco2,
     *               ikind,icmpn

      integer ixct,iexyn,ifact
      common/ cst37 /ixct,iexyn,ifact                       
C-----------------------------------------------------------------------   
C                               READ FILE HEADER
C-----------------------------------------------------------------------  
5000  format (a)
      write (6,*) 'input from file:'
      read (5,5000) fname(1)
      open (n2,file=fname(1),status='old')

      write (6,*) 'OUTPUT to file:'
      read (5,5000) fname(2)
      open  (n1,file=fname(2),status='unknown')
                                               

      icmpn = 12
4     read (n2,'(a132)') record
      write (n1,'(a132)') record
      if (record.eq.'end '.or.record.eq.'END ') goto 220
      goto 4

1000  FORMAT (3(A8,18X))
1010  FORMAT (I2,A78)
1020  FORMAT (A8,I2,I2,I2,I2)
1170  FORMAT (A)
1180  FORMAT (6(A5,1X)/6(A5,1X))                                        
C-----------------------------------------------------------------------   
C                             READ AND MODIFY INDIVIDUAL ENTRIES  
C----------------------------------------------------------------------- 
220   READ (N2,1030,END=999) NAME,IBASE,IKIND,ILAM,IDISO,SIXTY 
      WRITE (N1,1030) NAME,IBASE,IKIND,ILAM,IDISO,SIXTY 

      READ (N2,*) (COMP(ZI),ZI=1,ICMPN),THERM    


      WRITE (N1,1040) (COMP(ZI),ZI=1,ICMPN) 
      WRITE (N1,1050) THERM,0d0,0d0 
1040  FORMAT (12(F5.2,1X))
1050  FORMAT (5(D13.7,1X)) 
1060  FORMAT (A80)
      IF (ILAM.NE.0) THEN
C                               determine number of transitions from 
C                               flag ILAM:               
         JLAM=ILAM
         IF (ILAM.GT.3) JLAM=ILAM-3
         IF (ILAM.GT.6.and.ilam.lt.9) JLAM=ILAM-6  
         if (ilam.gt.9) jlam = ilam-9                                

         DO 51 ZI=1,JLAM
         READ (N2,*) (TM(ZJ,ZI),ZJ=1, m7 - 3) 
51       WRITE (N1,1050) (TM(ZJ,ZI),ZJ=1, m7 - 3),0. 
      END IF              
      IF (IDISO.NE.0) READ (N2,*) TD  
      IF (IDISO.NE.0) WRITE (N1,1050) TD  
      GOTO 220
C
999   STOP 
1030  FORMAT (A8,I2,I2,I2,I2,A60)
      END 
