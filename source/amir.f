
      subroutine DUMMY (ibef,p,t,rho,ks,gs,comps,ier,props,vp,vs)
c----------------------------------------------------------------------- 


      implicit none 
      include 'perplex_parameters.h'
c=========================================================================================
c        the second index identifies the endmember, the 1st index identifies the
c        property to be perturbed (1->energy, 2->volume, 3->Debye T or entropy.
c        on return ran could contain the absolute corrections (=rand*del).
c        k10 is 200, the actual number of endmembers is iphct
c        jphct is constant and greater than maxph.
c========================================================================================

      integer ier, i, j, maxph, ind,ibef

      real p,t

      parameter (maxph=41)
      double precision comps(5), rho, ks, gs,
     *                 vp, vs, props(maxph,8), tstop, tot

      double precision sprop,psys,psys1
      common/ cxt22 /sprop(i8,k5),psys(i8),psys1(i8)

      double precision pcomp
      common/ cst324 /pcomp(k5,k5)

      character pname*14
      common/ cxt21a /pname(k5)

      double precision a,b
      common/ cst313 /a(k5,k1),b(k5)

      integer jbulk
      double precision cblk
      common/ cst300 /cblk(k5),jbulk

      double precision atwt
      common/ cst45 /atwt(k0)

      integer kkp, np, ncpd, ntot
      double precision cp3, ctot3
      common/ cxt15 /cp3(k5,k5),ctot3(k5),kkp(k5),np,ncpd,ntot


      double precision pbar,tk,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /pbar,tk,xco2,u1,u2,tr,pr,r,ps

      integer junk
      double precision del, rand
      common/ cst321 /del(11,200),rand(12,200),junk

      integer ipoint, imyn
      common/ cst60  /ipoint,imyn

      integer iamir
      common/ cst53 /iamir
c----------------------------------------------------------------------

      if (ibef.eq.0) then
         iamir = 999
         call iniprp
         ibef = 1 
         junk = ipoint
      end if 


c1     write (*,*) 'enter t p(pa)'
c      read (*,*) t,p

      pbar = dble(p/1.e5)
      tk = dble(t+273.)
      tstop = 1073d0

      do i = 1,maxph
        do j = 1,8
          props(i,j) = 0.0
        enddo
      enddo   



c      write (*,*) 'enter comps'
c     read (*,*) (comps(i),i=1,5)


c     comps(1)=  2.90000000000000
c       comps(2)=   8.00000000000000
c       comps(3)=   35.1500000000000
c       comps(4)=   3.65000000000000
c       comps(5)=   49.9000000000000
c                                 assuming the composition is in weight
c                                 proportions create a molar composition
      tot = 0d0 

      do i = 1, 5
            cblk(i) = comps(i)/atwt(i)
            tot = tot + cblk(i)
      end do 

      do i = 1, 5
         b(i) = cblk(i)/tot
      end do 
c                                 do the minimization

         if (tk.lt.tstop) tk = tstop

         call lpopt0 (ier)

         tk = dble(t+273.15)

         if (ier.ne.0) then 
c            write (*,*) ' optimization failed'
             props(1,1) = -1.
         else 

            call getpar

c           call calpr0 (6)

            do i = 1, ntot

               if (pname(i).eq.'O(HP)'.or.pname(i).eq.'O(stx)') then
                  ind = 1
               else if (pname(i).eq.'Opx(HP)'.or.
     *                  pname(i).eq.'Opx(stx)') then
                  ind = 2
               else if (pname(i).eq.'Cpx(HP)'.or.
     *                  pname(i).eq.'Cpx(stx)') then
                  ind = 3
               else if (pname(i).eq.'Gt(HP)'.or.
     *                  pname(i).eq.'Gt(stx)') then
                  ind = 4
               else if (pname(i).eq.'an') then
                  ind = 5
               else if (pname(i).eq.'Sp(HP)'.or.
     *                  pname(i).eq.'Sp(stx)') then
                  ind = 6
               else if (pname(i).eq.'C2/c(stx)') then
                  ind = 7
               else if (pname(i).eq.'Wad(stx)') then
                  ind = 8
               else if (pname(i).eq.'Ring(stx)') then
                  ind = 9
               else if (pname(i).eq.'hCrd') then
                  ind = 10
               else if (pname(i).eq.'q') then
                  ind = 11
               else if (pname(i).eq.'coe') then
                  ind = 12
               else if (pname(i).eq.'stv') then
                  ind = 13
               else if (pname(i).eq.'trd') then
                  ind = 14
               else if (pname(i).eq.'crst') then
                  ind = 15
               else if (pname(i).eq.'and') then
                  ind = 16
               else if (pname(i).eq.'sil') then
                  ind = 17
               else if (pname(i).eq.'ky') then
                  ind = 18
               else if (pname(i).eq.'cor') then
                  ind = 19
               else if (pname(i).eq.'lime') then
                  ind = 20
               else if (pname(i).eq.'per') then
                  ind = 21
               else if (pname(i).eq.'wo') then
                  ind = 22
               else if (pname(i).eq.'aki') then
                  ind = 23
               else if (pname(i).eq.'perov') then
                  ind = 24
c                                                     real oddballs
               else if (pname(i).eq.'lrn') then 
                  ind = 25
               else if (pname(i).eq.'mont') then
                  ind = 26
               else if (pname(i).eq.'merw') then
                  ind = 27
               else if (pname(i).eq.'geh') then
                  ind = 28
               else if (pname(i).eq.'ak') then
                  ind = 29
               else if (pname(i).eq.'rnk') then
                  ind = 30
               else if (pname(i).eq.'cats') then
                  ind = 31
               else if (pname(i).eq.'pswo') then
                  ind = 32
               else if (pname(i).eq.'spr4') then
                  ind = 33
               else if (pname(i).eq.'spr7') then
                  ind = 34
               else if (pname(i).eq.'fspr') then
                  ind = 35
               else if (pname(i).eq.'Wus(fab)') then 
                  ind = 36
               else if (pname(i).eq.'Pv(fab)') then 
                  ind = 37
               else if (pname(i).eq.'Ppv(og)') then
                  ind = 38 
               else if (pname(i).eq.'Aki(fab)') then
                  ind = 39
               else if (pname(i).eq.'ca-pv') then
                  ind = 40
               else
                  ind = maxph 
                  write (189,*) ' i dunno this phase:',pname(i)
               end if 
c                                 shear modulus 1/bar
               props(ind,1) = sprop(5,i) 
c                                 bulk modulus 1/bar
               props(ind,2) = sprop(4,i)
c                                 volumetric mode, absolute fraction
               props(ind,3) = sprop(1,i)*sprop(16,i)/psys(1)
c                                 vp km/s 
               props(ind,4) = sprop(7,i)
c                                 vs km/s 
               props(ind,5) = sprop(8,i)
c                                 rho kg/m3 
               props(ind,6) =  sprop(17,i)/sprop(1,i)*1d2
c                                 molar mg/(mg+fe)
               tot = cp3(2,i) + cp3(3,i)
               if (tot.ne.0d0) then  
                  props(ind,7) = cp3(3,i)/tot
               else 
                  props(ind,7) = 0d0
               end if   
c                                 molar al2/(al2+si) 
               tot = cp3(4,i) + cp3(5,i)
               if (tot.ne.0d0) then  
                  props(ind,8) = cp3(4,i)/tot
               else 
                  props(ind,8) = 0d0
               end if

               do j = 1, 8
                  if (props(ind,j).lt.0d0) ier = 1
               end do 

               do j = 1, 6
                  if (props(ind,j).eq.0d0) ier = 1
               end do 

c              write (*,1010) pname(i),(props(ind,j),j=1,8)

            end do 

            rho = psys(10)
            ks = psys(4)
            gs = psys(5)
            vp = psys(7)
            vs = psys(8)

c           write (*,1000) rho, ks, gs, vp, vs

         end if 
 
1000  format (/,' rho (kg/m3) = ',g14.6,/,
     *          ' adiabatic bulk modulus (bar) = ',g14.6,/,
     *          ' adiabatic shear modulus (bar) = ',g14.6,/
     *          ' vp (km/s) = ',g14.6,/,' vs (km/s)  = ',g14.6,/)
1010  format (/,' stable phase  shear modulus  bulk modulus  ',
     *          ' volume mode        vp            vs           rho',/,
     *          12(1x,a,2x,8(g14.6,1x),/))

c     goto 1

      end 
    
      subroutine iniprp
c----------------------------------------------------------------------
c iniprp - read data files and initialization for mingee
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical vertex, output 

      character n4name*100

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      integer iemod,kmod
      logical smod
      double precision emod
      common/ cst319 /emod(k15,k10),smod(h9),iemod(k10),kmod
c----------------------------------------------------------------------- 
      vertex = .true.
c                              elastic modulii flag
      kmod = 0 
c                              get runtime parameters
      call redop1
c                            -------------------------------------------
c                              open statements for units n1-n5 and n9
c                              are in subroutine input1
      call input1 (vertex,n4name)
c                              read thermodynamic data on unit n2:
      call input2 (vertex)
c                              read data for solution phases on n9:
      call input9 (vertex,output)
c                              call initlp to initialize arrays 
c                              for optimization.
      call initlp     

      end

      subroutine getpar 
c-----------------------------------------------------------------------
c getpar loads the local version of the parameters and compositional 
c coordinates of an assemblage. 

c    ncpdg  -> ncpd
c    npg    -> np
c    idbulk -> kkp
c    xcoor  -> x3

c also sets flags (could set a solvus flag):

c    aflu  -> fluid present
c    fluid -> the phase is a fluid (indexed)

c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,j,ids,iwarn

      double precision mols, root, vol, chi, pgeo(i8),
     *                 pgeo1(i8), chi1

      logical ok, sick(i8), nodata, ssick
c                                 -------------------------------------
c                                 global variables

c                                 x-coordinates for the assemblage solutions
      double precision x3
      common/ cxt16 /x3(k21,mst,msp)
c                                 
      integer ixp,ifp
      double precision sxs,exces
      common/ cst304 /sxs(k13),exces(m3,k1),ixp(k1),ifp(k1)
c                                 composition and model flags
c                                 for final adaptive solution
      integer kkp, np, ncpd, ntot
      double precision cp3, ctot3
      common/ cxt15 /cp3(k5,k5),ctot3(k5),kkp(k5),np,ncpd,ntot
c                                 x coordinate description
      integer istg, ispg, imlt, imdg
      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)
c                                 bookkeeping variables
      integer ksmod, ksite, kmsol, knsp
      common/ cxt0  /ksmod(h9),ksite(h9),kmsol(h9,m4,mst),knsp(h9,m4)

      logical gflu,aflu,fluid,shear,lflu,volume
      common/ cxt20 /gflu,aflu,fluid(k5),shear,lflu,volume

      double precision gtot,fbulk,gtot1,fbulk1
      common/ cxt81 /gtot,fbulk(k5),gtot1,fbulk1(k5)

      double precision cp
      common/ cst12 /cp(k5,k1)

      character pname*14
      common/ cxt21a /pname(k5)

      character names*8, fname*10
      common/ cst8  /names(k1)/ csta7 /fname(h9)

      double precision atwt
      common/ cst45 /atwt(k0)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      double precision pcomp
      common/ cst324 /pcomp(k5,k5)

      integer iopt
      double precision nopt
      common / opts /nopt(i10),iopt(i10)

      double precision props,psys,psys1
      common/ cxt22 /props(i8,k5),psys(i8),psys1(i8)

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps
c                                 molar amounts (b)
      integer ipvt
      double precision a,b
      common/ cst301 /a(k5,k5),b(k5),ipvt(k5)

      save iwarn
      data iwarn/0/
c----------------------------------------------------------------------
      aflu = .false.
      shear = .true.
      volume = .true.
      nodata = .false.
      ssick = .false.
c                                 flag for bulk bad bulk properties
      do i = 1, i8
         sick(i) = .false.
      end do 
c                                 weighting scheme for seismic velocity
c                                 chi = 1 -> voigt
c                                 0 -> reuss
      chi = nopt(6)
      chi1 = 1d0 - chi

      do i = 1, ntot

         ids = kkp(i)

         if (i.le.np) then 

            if (ksmod(ids).eq.0.or.ksmod(ids).eq.26) then 
               aflu = .true.
               fluid(i) = .true.
            else 
               fluid(i) = .false.
            end if 

         else 

            if (ifp(-ids).eq.1) then 
               aflu = .true.
               fluid(i) = .true.
            else 
               fluid(i) = .false.
            end if

         end if
c                                 molar amounts
         props(16,i) = b(i)
c                                 convert x3 to y for calls to gsol            
         if (ids.gt.0) call x3toy (i,ids)
c                                 get shear moduli
         if (.not.fluid(i)) then 
            call moduli (ids,props(5,i),ok)
            if (.not.ok) shear = .false.  
         else
            props(5,i) = 0d0    
         end if      
c                                 get thermodynamic props  
         call dphase (ids,i)         
c                                 gruneisen parameter
         props(3,i) = props(1,i)/
     *                      (props(12,i)*props(14,i)/props(13,i) - 
     *                      v(2)*props(13,i)*props(1,i))
c                                 adiabatic bulk modulus
         props(4,i) = (1d0 + v(2)*props(13,i)
     *                           *props(3,i))/props(14,i)
c                                 like this so sick could be used to
c                                 say which property is bad
         if (props(1,i).le.0d0) sick(i) = .true.
         if (props(3,i).le.0d0) sick(i) = .true.
         if (props(4,i).le.0d0) sick(i) = .true.
         if (props(12,i).le.0d0) sick(i) = .true.
         if (props(13,i).le.0d0) sick(i) = .true.
         if (props(14,i).le.0d0) sick(i) = .true.
         if (sick(i)) volume = .false.
         if (sick(i).and.(.not.fluid(i))) ssick = .true.

      end do 
c                                 total weight of assemblage
      gtot = 0d0
      gtot1 = 0d0

      do i = 1, i8
         psys(i) = 0d0
         psys1(i) = 0d0
         pgeo(i) = 0d0 
         pgeo1(i) = 0d0
      end do 

      do i = 1, icomp
c                                 molar ammounts
         fbulk(i) = 0d0
         fbulk1(i) = 0d0

      end do

      do i = 1, ntot

         ids = kkp(i)
c                                 get atomic weight
         props(17,i) = 0d0

         do j = 1, icomp
c                                 formula weight
            props(17,i) = props(17,i) + cp3(j,i) * atwt(j) 
c                                 molar amounts of the components
            mols = props(16,i)*cp3(j,i)
c                                 mass of the components
            fbulk(j) = fbulk(j) + mols
            gtot = gtot + mols

            if (.not.fluid(i)) then 
               fbulk1(j) = fbulk1(j) + mols
               gtot1 = gtot1 + mols
            end if 
c                                 molar or weight phase compositions
            pcomp(j,i) = cp3(j,i)

         end do  
c                                 density 
         props(10,i) = props(17,i)/props(1,i)*1d2  

         if (.not.sick(i)) then 
c                                 sound velocity (km/s), for phases 
c                                 with h&p landau transitions (e.g. qtz)
c                                 the bulk modulus may be negative.
            root = props(4,i)*1d5/props(10,i)
            props(6,i) = dsqrt(root)/1d3
c                                 s-wave velocity
            root = (props(4,i)+4d0*props(5,i)/3d0)*1d5/props(10,i)
            if (root.ge.0d0) then 
               props(7,i) = dsqrt(root)/1d3
            else
               props(7,i) = nopt(7)
               shear = .false.
            end if 
c                                 p-wave velocity
            root = props(5,i)*1d5/props(10,i)
            if (root.ge.0d0) then 
               props(8,i) = dsqrt(root)/1d3
            else 
               props(8,i) = nopt(7)
               shear = .false.
            end if 
c                                 vp/vs
            if (props(7,i).ne.0d0) then 
               props(9,i) = props(8,i)/props(7,i)
            else
               props(9,i) = nopt(7)
            end if 

         else 
            do j = 3, 9
               props(j,i) = nopt(7)
            end do 
         end if 

         if (iopt(2).eq.1) then 
c                                 convert molar phase composition to 
c                                 mass % composition:
            do j = 1, icomp
               pcomp(j,i) = pcomp(j,i)/props(17,i)*1d2
            end do 
         end if      
c                                 make a name:
         if (ids.lt.0) then
c                                 simple compound:
            pname(i) = names(-ids)

         else  
c                                 solution phases:
            pname(i) = fname(ids)
         end if 

c                                 -------------------------------------
c                                 system properties:
c                                 vol of phase per mole of system
         vol = props(16,i)*props(1,i)
c                                 system molar volume
         psys(1)  = psys(1)  + vol
c                                 molar enthalpy
         psys(2)  = psys(2)  + props(2,i)*props(16,i) 
c                                 molar heat capacity
         psys(12) = psys(12) + props(12,i)*props(16,i) 
c                                 expansivity
         psys(13) = psys(13) + props(13,i)*vol 
c                                 compressibility
         psys(14) = psys(14) + props(14,i)*vol
c                                 molar entropy
         psys(15) = psys(15) + props(15,i)*props(16,i) 
c                                 moles of assemblage
         psys(16) = psys(16) + props(16,i)
c                                 mass of assemblage 
         psys(17) = psys(17) + props(17,i)*props(16,i)
       
         if (volume) then 
c                                 gruneisen
            psys(3) = psys(3) + vol*props(3,i)*chi
            pgeo(3) = pgeo(3) + vol/props(3,i)
c                                 Aggregate Bulk Modulus                                
            psys(4) = psys(4) + vol*props(4,i)*chi
            pgeo(4) = pgeo(4) + vol/props(4,i)
c                                 Aggregate Sound
            psys(6) = psys(6) + vol*props(6,i)*chi
            pgeo(6) = pgeo(6) + vol/props(6,i)

            if (shear) then 
c                                 Arithmetic means for mu and Vs
               psys(5) = psys(5) + vol*props(5,i)*chi
c                                 Aggregate Vp                                
               psys(7) = psys(7) + vol*props(7,i)*chi
               pgeo(7) = pgeo(7) + vol/props(7,i)

               psys(8) = psys(8) + vol*props(8,i)*chi
c                                 Aggregate Shear Modulus, only if
c                                 shear mod is available for all cpds.
               pgeo(5) = pgeo(5) + vol/props(5,i)
c                                 Aggregate Vs                               
               pgeo(8) = pgeo(8) + vol/props(8,i)
            end if 
         end if


         if (aflu) then 
c                                 assemblage includes fluid
            if (.not.fluid(i)) then
c                                 get total without fluid
               psys1(1)  = psys1(1) + vol 
c                                 molar enthalpy
               psys1(2)  = psys1(2) + props(2,i)*props(16,i) 
c                                 molar heat capacity
               psys1(12) = psys1(12) + props(12,i)*props(16,i) 
c                                 molar expansivity
               psys1(13) = psys1(13) + props(13,i)*vol 
c                                 molar compressibility
               psys1(14) = psys1(14) + props(14,i)*vol 
c                                 molar entropy
               psys1(15) = psys1(15) + props(15,i)*props(16,i) 
c                                 total number of moles of phases
               psys1(16) = psys1(16) + props(16,i)

               psys1(17) = psys1(17) + props(16,i)*props(17,i)

               if (.not.ssick) then 
c                                 gruneisen
                  psys1(3) = psys1(3) + vol*props(3,i)*chi
                  pgeo1(3) = pgeo1(3) + vol/props(3,i)
c                                 Aggregate Bulk Modulus                                
                  psys1(4) = psys1(4) + vol*props(4,i)*chi
                  pgeo1(4) = pgeo1(4) + vol/props(4,i)
c                                 Aggregate Sound
                  psys1(6) = psys1(6) + vol*props(6,i)*chi
                  pgeo1(6) = pgeo1(6) + vol/props(6,i)


                  if (shear) then 
c                                 Aggregate Vp                                
                     psys1(7) = psys1(7) + vol*props(7,i)*chi
                     pgeo1(7) = pgeo1(7) + vol/props(7,i)
c                                 Arithmetic means for mu and Vs
                     psys1(5) = psys1(5) + vol*props(5,i)*chi
                     psys1(8) = psys1(8) + vol*props(8,i)*chi
c                                 Aggregate Shear Modulus, only
c                                 if shear mod is available for all cpds.
                     pgeo1(5) = pgeo1(5) + vol/props(5,i)
c                                 Aggregate Vs                               
                     pgeo1(8) = pgeo1(8) + vol/props(8,i)
                  end if 
               end if

            end if 
         end if 
      end do 
c                                 correct for proportional wt 
c                                 on intensive properties (alpha, beta). 
      do i = 3, 8
         if ((i.eq.3.or.i.eq.4.or.i.eq.6.and.volume).or.
     *                             shear.and.volume) then
            psys(i) = psys(i)/psys(1)
            pgeo(i) = pgeo(i)/psys(1)
            if (psys1(1).ne.0d0) then  
               psys1(i) = psys1(i)/psys1(1)
               pgeo1(i) = pgeo1(i)/psys1(1)
            end if 

         end if 
      end do 

      psys(13) = psys(13)/psys(1)
      psys(14) = psys(14)/psys(1)

      if (psys1(1).ne.0d0) then 
         psys1(13) = psys1(13)/psys1(1)
         psys1(14) = psys1(14)/psys1(1)   
      end if 
c                                 aggregate properties
      if (psys1(1).ne.0d0) psys1(10) = psys1(17)/psys1(1)

      psys(3) = psys(3) + chi1/pgeo(3)
c                                 bulk modulus
      psys(4) = psys(4) + chi1/pgeo(4)
c                                 density, kg/m3
      psys(10) = psys(17)/psys(1)*1d2
c                                 sound velocity
      if (psys(4).gt.0d0) psys(6) = dsqrt(psys(4)*1d5/psys(10))/1d3

      if (volume.and.shear.and.pgeo(5).gt.0d0) then 
c                                 there are shear moduli for every phase:
c                                 shear modulus
         psys(5) = psys(5) + chi1/pgeo(5)
c                                 p-wave velocity
         psys(7) = dsqrt((psys(4)+4d0*psys(5)/3d0)*1d5/psys(10))/1d3
c                                 s-wave velocity
         psys(8) = dsqrt(psys(5)*1d5/psys(10))/1d3
         psys(9) = psys(8)/psys(7)

      else if (volume.and.shear.and.chi.gt.0) then 
c                                 shear moduli for everything but fluid
c                                 then don't add harmonic mean:
c                                 shear modulus
         psys(5) = psys(5)/chi
c                                 p-wave velocity
         psys(7) = dsqrt((psys(4)+4d0*psys(5)/3d0)*1d5/psys(10))/1d3
c                                 s-wave velocity
         psys(8) = dsqrt(psys(5)*1d5/psys(10))/1d3
         psys(9) = psys(8)/psys(7)
      else 
         do i = 3, 9
            psys(i) = nopt(7)
         end do 
      end if 

      if (aflu) then 
c                                 assemblage contains fluid, compute
c                                 fluid absent properties:
c                                 molar volume
         psys1(3) = psys1(3) + chi1/pgeo1(3)
         psys1(4) = psys1(4) + chi1/pgeo1(4) 

         psys1(10) = psys1(17)/psys1(1)*1d2
c                                 sound velocity
         psys1(6) = dsqrt(psys1(4)*1d5/psys1(10))/1d3

         if (shear.and.pgeo1(5).gt.0d0) then 
c                                 shear modulus
            psys1(5) = psys1(5) + chi1/pgeo1(5)
c                                 p-wave velocity
            psys1(7) = dsqrt((psys1(4)+4d0*psys1(5)/3d0)*1d5/psys1(10))
     *                 /1d3
c                                 s-wave velocity
            psys1(8) = dsqrt(psys1(5)*1d5/psys1(10))/1d3

            psys1(9) = psys1(7)/psys1(8)

         else if (shear.and.chi.gt.0) then 
c                                 shear moduli for everything but fluid
c                                 then don't add harmonic mean:
c                                 shear modulus
            psys1(5) = psys1(5)/chi
c                                 p-wave velocity
            psys1(7) = dsqrt((psys1(4)+4d0*psys1(5)/3d0)*1d5
     *                                       /psys1(10))/1d3
c                                 s-wave velocity
            psys1(8) = dsqrt(psys1(5)*1d5/psys1(10))/1d3
            psys1(9) = psys1(8)/psys1(7)

         else 

            psys1(3) = nopt(7)
            psys1(4) = nopt(7)
            psys1(5) = nopt(7)
            psys1(6) = nopt(7)
            psys1(7) = nopt(7)
            psys1(8) = nopt(7)
            psys1(9) = nopt(7)

         end if 
         
      end if 

      if (iwarn.lt.100) then 
         if (.not.volume.or..not.shear.and.(iwarn.lt.101)) then

            iwarn = iwarn + 1

            if (.not.shear.and.volume) then
                write (*,1000) v(1),v(2)
             else if (.not.volume) then 
               do i = 1, ntot
                  if (sick(i)) write (*,1010) v(1),v(2),pname(i)
               end do 
            end if 

            if (iwarn.eq.100) write (*,1020) 
                 
         end if 
      end if 

1000  format (/,'**WARNING** at T(K)=',g12.4,' P(bar)=',g12.4,1x,
     *        'aggregate seismic properties',/,'cannot be ',
     *        'computed because',
     *        ' of a missing/invalid shear modulus.',/)
1010  format (/,'**WARNING** at T(K)=',g12.4,' P(bar)=',g12.4,1x,
     *        'aggregate seismic properties ',/,'cannot be ',
     *        'computed because of missing/invalid properties ',/,
     *        'for phase:',a,/)
1020  format ('This warning will not be repeated for future instances',
     *        ' of the problem.',/)
99    end 

      subroutine dphase (id,jd)
c-----------------------------------------------------------------------
c dphase sets the props of the phase identified by id as 
c computed by centered finite differences from the Gibbs energy
c as stored in props(i8,jd)

c the difference increments are

c dt0, dp0 for 1st order derivatives (entropy,volume and enthalpy)
c dt1, dp1 for 2nd order derivatives (heat capacity, expansivity*, 
c          compressibility*)
c dt2, dp2 for 3rd order derivatives (alphat, alphap = betat, betap
c          alphtt), and the fourth order derivative alphtt, this 
c          derivative was used by werami prior to 5/31/04. 

c *expansivity (alpha) as returned here is 1/v*dv/dt
c *compressibility (beta) as returned here is -1/v*dv/dp

c corrected to check for negative pressure, in which case 
c forward differences are used. june 22, 2004, JADC.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer id,jd

      double precision dt0,dt1,dt2,g0,g1a, g2a, dg, ss,
     *                 dp0,dp1,dp2,e,alpha,v,ginc,beta,cp,s

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      double precision props,psys,psys1
      common/ cxt22 /props(i8,k5),psys(i8),psys1(i8)

      integer iopt
      double precision nopt
      common / opts /nopt(i10),iopt(i10)

      save dt0,dt1,dt2,dp0,dp1,dp2
      data dt0,dt1,dt2,dp0,dp1,dp2/0.5d0,5.,50.,5.,50.,500./
c----------------------------------------------------------------------

      dp0 = 0.5d-3 * p 
      dp1 = 0.5d-2 * p
      dp2 = 0.5d-1 * p
            
      g0 = ginc(0d0,0d0,id)
c                                 straight derivatives:
c                                 first order
      if (p-dp0.le.0d0) then 

         v = (ginc(0d0,dp0,id) - ginc(0d0,0d0,id))/dp0
         if (v.lt.0d0.or.dabs(v).gt.1d9)  
c                                 expand increment if invalid v
     *   v = (ginc(0d0,dp1,id) - ginc(0d0,0d0,id))/dp1
         if (v.lt.0d0.or.dabs(v).gt.1d9)  
c                                 expand increment more if invalid v
     *   v = (ginc(0d0,dp2,id) - ginc(0d0,0d0,id))/dp2

      else 

         v = (ginc(0d0,dp0,id) - ginc(0d0,-dp0,id))/dp0/2d0
         if ((v.lt.0d0.or.dabs(v).gt.1d9).and.p-dp1.gt.0d0)  
c                                 expand increment if invalid v
     *   v = (ginc(0d0,dp1,id) - ginc(0d0,-dp1,id))/dp1/2d0
         if ((v.lt.0d0.or.dabs(v).gt.1d9).and.p-dp2.gt.0d0)  
c                                 expand increment more if invalid v
     *   v = (ginc(0d0,dp2,id) - ginc(0d0,-dp2,id))/dp2/2d0

      end if 

      s = (ginc(-dt0,0d0,id) - ginc(dt0,0d0,id))/dt0/2d0
c                                 this crap is necessary because 
c                                 optimization or my bad programming
c                                 corrupts ginc with compaq visual fortran.
      g1a = ginc(-dt0,0d0,id)
      g2a =  ginc(dt0,0d0,id)
      dg = g1a-g2a
      ss = dg/dt0/2d0
      s = ss
c
c 
c     write (*,*) s, ss, dg, ginc(-dt0,0d0,id) - ginc(dt0,0d0,id), dt0

      e = g0 + t * s
c                                 second order
      cp = -t*(ginc(dt1,0d0,id) + ginc(-dt1,0d0,id) - 2d0*g0)/dt1/dt1
      if (cp.lt.0d0.or.dabs(cp).gt.1d9)  
c                                 expand increment if invalid cp
     *   cp = -t*(ginc(dt2,0d0,id) + ginc(-dt2,0d0,id) - 2d0*g0)/dt2/dt2
      if (cp.lt.0d0.or.dabs(cp).gt.1d9)  
c                                 shrink increment if invalid cp
     *   cp = -t*(ginc(dt0,0d0,id) + ginc(-dt0,0d0,id) - 2d0*g0)/dt0/dt0
   
      if (p-dp1.le.0d0) then 
c                                 use forward difference at small p's
         beta = (ginc(0d0,2d0*dp1,id) + g0 - 2d0*ginc(0d0,dp1,id))
     *          /dp1/dp1
         if (dabs(beta).gt.v)
c                                 expand increment if invalid beta
     *   beta = (ginc(0d0,2d0*dp2,id) + g0 - 2d0*ginc(0d0,dp2,id))
     *          /dp2/dp2                                 
         if (dabs(beta).gt.v)
c                                 shrink increment if invalid beta
     *   beta = (ginc(0d0,2d0*dp0,id) + g0 - 2d0*ginc(0d0,dp0,id))
     *          /dp0/dp0   

         alpha = ( ginc( dt1,dp1,id) - ginc( dt1,0d0,id)
     *            -ginc(-dt1,dp1,id) + ginc(-dt1,0d0,id))/dp1/dt1/2d0
         if (dabs(alpha).gt.v.or.alpha.lt.0d0)
c                                 expand increment if invalid alpha
     *   alpha = ( ginc( dt2,dp2,id) - ginc( dt2,0d0,id)
     *            -ginc(-dt2,dp2,id) + ginc(-dt2,0d0,id))/dp2/dt2/2d0
         if (dabs(alpha).gt.v.or.alpha.lt.0d0)
c                                 shrink increment if invalid alpha
     *   alpha = ( ginc( dt0,dp0,id) - ginc( dt0,0d0,id)
     *            -ginc(-dt0,dp0,id) + ginc(-dt0,0d0,id))/dp0/dt0/2d0

      else 
         beta = (ginc(0d0,dp1,id) + ginc(0d0,-dp1,id) - 2d0*g0)/dp1/dp1
         if (dabs(beta).gt.v.and.p-dp2.gt.0d0)
c                                 expand increment if invalid beta
     *   beta = (ginc(0d0,dp2,id) + ginc(0d0,-dp2,id) - 2d0*g0)/dp2/dp2                                
         if (dabs(beta).gt.v)
c                                 shrink increment if invalid beta
     *   beta = (ginc(0d0,dp0,id) + ginc(0d0,-dp0,id) - 2d0*g0)/dp0/dp0

         alpha = ( ginc( dt1,dp1,id) - ginc( dt1,-dp1,id)
     *            -ginc(-dt1,dp1,id) + ginc(-dt1,-dp1,id))/dp1/dt1/4d0
         if (dabs(alpha).gt.v.or.alpha.lt.0d0)
c                                 expand increment if invalid alpha
     *   alpha = ( ginc( dt2,dp2,id) - ginc( dt2,-dp2,id)
     *            -ginc(-dt2,dp2,id) + ginc(-dt2,-dp2,id))/dp2/dt2/4d0
         if (dabs(alpha).gt.v.or.alpha.lt.0d0)
c                                 shrink increment if invalid alpha
     *   alpha = ( ginc( dt0,dp0,id) - ginc( dt0,-dp0,id)
     *            -ginc(-dt0,dp0,id) + ginc(-dt0,-dp0,id))/dp0/dt2/4d0      
     
      end if 
c                                 convert to their normal forms:
      beta = -beta/v
      alpha = alpha/v
c                                 
      if (beta.gt.v.or.beta.lt.0d0) beta = nopt(7)
      if (alpha.gt.v.or.alpha.lt.0d0) alpha = nopt(7)
      if (cp.gt.1d9.or.cp.lt.0d0) cp = nopt(7)

      props(1,jd) = v 
      props(2,jd) = e
      props(12,jd) = cp
      props(13,jd) = alpha
      props(14,jd) = beta
      props(15,jd) = s

      end 

      double precision function ginc (dt,dp,id)
c-----------------------------------------------------------------------
      implicit none

      double precision dt,dp,p,t,xco2,u1,u2,tr,pr,r,ps,gee

      double precision gsol

      integer id

      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      p = p + dp 
      t = t + dt 

      call tset
      gee = gsol(id)

      p = p - dp 
      t = t - dt

      ginc = gee 

      end 

      subroutine moduli (ids,mu,ok) 
c-----------------------------------------------------------------------
c subroutine moduli determines shear moduli (mods) for entity ids, returns
c ok = false if moduli are unavailable.

c jmod is the number of cpds for which it is possible to calculate coeffs.
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      double precision mu, pmu

      integer i, ids, id

      logical ok
c                                 bookkeeping variables
      integer ksmod, ksite, kmsol, knsp
      common/ cxt0  /ksmod(h9),ksite(h9),kmsol(h9,m4,mst),knsp(h9,m4)

      integer ndep
      double precision y2pg
      common/ cxt4  /y2pg(m15,m4,h9)

      integer jend
      common/ cxt23 /jend(h9,k12)

      integer iemod,kmod
      logical smod
      double precision emod
      common/ cst319 /emod(k15,k10),smod(h9),iemod(k10),kmod

      double precision z, pa, p0a, x, w, y
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1)

      integer iopt
      double precision nopt
      common / opts /nopt(i10),iopt(i10)
c-----------------------------------------------------------------------
  
      ok = .true.

      mu = nopt(7)

      if (ids.le.0) then 

         if (iemod(-ids).ne.0) then

            call shearm (mu,-ids)

         else

            ok = .false.

         end if 

      else 

         if (smod(ids)) then 

            if (ksmod(ids).ne.7.and.ksmod(ids).ne.8) then
c                                 for solutions with no dependent endmembers
c                                 the y coordinates can be used to compute 
c                                 the composition
               do i = 1, mstot(ids)

                  id = jend(ids,2+i)

                  call shearm (pmu,id)

                  mu = mu + y(i) * pmu

               end do

            else 
c                                 get the p' coordinates (amounts of 
c                                 the independent disordered endmembers)     
               call getpp (ids) 

               do i = 1, lstot(ids)

                  id = jend(ids,2+i)

                  call shearm (pmu,id)

                  mu = mu + z(i) * pmu

               end do
 
            end if 

            if (mu.lt.0d0) then 
               mu = nopt(7)
               ok = .false.
            end if 
       
         else

            ok = .false.

         end if 

      end if  

      end

      subroutine shearm (mu,id)
c-----------------------------------------------------------------------
c shearm returns a linear model for the adiabatic shear modulus
c relative to the current pressure and temperature.

c three cases:

c make(id) = non-zero, use make definition to compute shear modulus.

c iemod = 1, linear model is input

c iemod = 2, shear modulus is known as a function of (V,T), then
c computed by centered finite differences.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer id

      double precision smu,mu,g,ginc

      common/ cst323 /smu

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      integer make
      common / cst335 /make(k10)

      integer iemod,kmod
      logical smod
      double precision emod
      common/ cst319 /emod(k15,k10),smod(h9),iemod(k10),kmod
c-----------------------------------------------------------------------

      if (make(id).ne.0) then 

         call makmod (id,mu)

      else if (iemod(id).eq.1) then 

         mu  = emod(1,id) + (p-pr)*emod(2,id) + (t-tr)*emod(3,id)

      else if (iemod(id).eq.2) then 
c                                 by calling ginc a call to
c                                 stixrudes EoS for the adiabatic
c                                 shear modulus is implicit (cst323)
         g = ginc(0d0,0d0,-id)
         mu = smu

      end if          

      end 

      subroutine makmod (id,mu)
c-----------------------------------------------------------------------
c gmake computes and sums the component g's for a make definition.
c the component g's may be calculated redundantly because gmake is
c called by gcpd, which in turn may be called by routines that call
c for a single g (e.g., gphase). 
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i, j, id, jd

      double precision mu, pmu

      double precision mkcoef, mdqf

      integer mknum, mkind
      common / cst334 /mkcoef(k16,k17),mdqf(k16,k17),mkind(k16,k17),
     *                 mknum(k16)

      integer make
      common / cst335 /make(k10)
c-----------------------------------------------------------------------

      jd = make(id)

      mu = 0d0
 
c                                compute the sum of the component g's
      do i = 1, mknum(jd)

         call shearm (pmu,mkind(jd,i))

         do j = 1, 3
            mu = mu + mkcoef(jd,i) * pmu
         end do 

      end do 

      end 

      subroutine calpr0 (n3)
c----------------------------------------------------------------------
c calpr0 - output properties of an assemblage
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'
 
      character cprop*18

      integer i,j,l,n3
c                                 -------------------------------------
c                                 global variables
      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      character cname*5
      common/ csta4  /cname(k5) 

      double precision gtot,fbulk,gtot1,fbulk1
      common/ cxt81 /gtot,fbulk(k5),gtot1,fbulk1(k5)

      double precision props,psys,psys1
      common/ cxt22 /props(i8,k5),psys(i8),psys1(i8)

      double precision pcomp
      common/ cst324 /pcomp(k5,k5)

      character pname*14
      common/ cxt21a /pname(k5)

      integer iopt
      double precision nopt
      common / opts /nopt(i10),iopt(i10)

      integer kkp, np, ncpd, ntot
      double precision cp3, ctot3
      common/ cxt15 /cp3(k5,k5),ctot3(k5),kkp(k5),np,ncpd,ntot

      double precision atwt
      common/ cst45 /atwt(k0)

      logical gflu,aflu,fluid,shear,lflu,volume
      common/ cxt20 /gflu,aflu,fluid(k5),shear,lflu,volume

      integer jbulk
      double precision cblk
      common/ cst300 /cblk(k5),jbulk

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      integer ipot,jv,iv
      common / cst24 /ipot,jv(l2),iv(l2)

      character*8 vname,xname
      common/ csta2  /xname(k5),vname(l2)
c---------------------------------------------------------------------- 
                                    
      write (n3,1000)  
      write (n3,1120) (vname(jv(i)),v(jv(i)), i = 1, ipot)
      write (n3,1120) (vname(jv(i)),v(jv(i)), i = 3, ipot)

      if (iopt(2).eq.0) then 
         cprop = 'molar  proportions'
      else
         cprop = 'weight percentages'
      end if

      write (n3,1020) cprop, (cname(i), i = 1, icomp)

      do i = 1, ntot

         write (n3,1030) pname(i), 
c                                 weight %
     *                   props(17,i)*props(16,i)/psys(17)*1d2,
c                                 vol %
     *                   props(1,i)*props(16,i)/psys(1)*1d2,
c                                 mol %
     *                   props(16,i)/psys(16)*1d2,
c                                 mol
     *                   props(16,i),
c                                 molar or weight composition
     *                   (pcomp(l,i), l = 1, icomp)
      end do 

      write (n3,1160)
c                                 phase/system summary, normal thermo:
         do i = 1, ntot
c                                 N, H, V, Cp, alpha, beta
            write (n3,1170) pname(i),props(17,i),props(2,i),props(15,i),
     *                      props(1,i),(props(j,i),j=12,14),
     *                      props(10,i)
         end do

         write (n3,1170) 'System        ',psys(17),psys(2),
     *                    psys(15),psys(1),psys(12),psys(13),psys(14)

         if (aflu) 
     *   write (n3,1170) 'System - Fluid',psys1(17),psys1(2),psys1(15),
     *                    psys1(1),psys1(12),psys1(13),psys1(14)

         write (n3,1190)
c                                 phase/system summary, seismic:
         do i = 1, ntot
            write (n3,1200) pname(i), (props(j,i), j = 3, 8) 
         end do

      if (.not.aflu) then 
c                                 no fluid is present:
         write (n3,1210)

         write (n3,1040)

         do i = 1, icomp

            write (n3,1110) cname(i),fbulk(i), fbulk(i)/gtot*1d2,
     *                      fbulk(i)*atwt(i)/psys(17)*1d2                            
         end do

         write (n3,1220)

         write (n3,1060) psys(10), 
c                                 enthalpy, specific enthalpy
     *                   psys(2)/psys(1)*1d5/psys(10),
     *                   psys(2)/psys(1)*1d5, 
c                                 entropy, specific entropy 
     *                   psys(15)/psys(1)*1d5/psys(10),
     *                   psys(15)/psys(1)*1d5,
c                                 cp, specific cp 
     *                   psys(12)/psys(1)*1d5/psys(10),
     *                   psys(12)/psys(1)*1d5,(psys(j),j= 3, 8)

      else 
c                                 fluid is present
         write (n3,1210)

         write (n3,1080)

         do i = 1, icomp

            write (n3,1110) cname(i),fbulk(i),fbulk(i)/gtot*1d2,
     *               fbulk(i)*atwt(i)/psys(17)*1d2,fbulk1(i)/gtot1*1d2,
     *               fbulk1(i)*atwt(i)/psys1(17)*1d2          
         end do

         write (n3,1220)
c                                 true bulk properties:
         write (n3,1060) psys(10), 
c                                 enthalpy, specific enthalpy
     *                   psys(2)/psys(1)*1d5/psys(10),
     *                   psys(2)/psys(1)*1d5, 
c                                 entropy, specific entropy 
     *                   psys(15)/psys(1)*1d5/psys(10),
     *                   psys(15)/psys(1)*1d5,
c                                 cp, specific cp 
     *                   psys(12)/psys(1)*1d5/psys(10),
     *                   psys(12)/psys(1)*1d5,(psys(j),j= 3, 8)
c                                 solid only bulk properties:
         write (n3,1100) psys1(10), 
c                                 enthalpy, specific enthalpy
     *                   psys1(2)/psys1(1)*1d5/psys1(10),
     *                   psys1(2)/psys1(1)*1d5, 
c                                 entropy, specific entropy 
     *                   psys1(15)/psys1(1)*1d5/psys1(10),
     *                   psys1(15)/psys1(1)*1d5,
c                                 cp, specific cp 
     *                   psys1(12)/psys1(1)*1d5/psys(10),
     *                   psys1(12)/psys1(1)*1d5,(psys1(j),j= 3, 8)

         write (n3,1230) 

      end if 

      write (n3,1070) 2, jbulk - ntot + 2 

1000  format (/,40('-'),//,'Stable phases at:')
1020  format (/,'Phase Compositions (',a,'):',
     *        /,19x,'wt %',6x,'vol %',5x,'mol %',5x,'mol  ',
     *          5x,20(1x,a,3x))
1030  format (1x,a,3x,3(f6.2,4x),g9.3,1x,20(f7.3,2x))
1040  format (/,14x,'mol',7x,'mol %',6x,'wt %')
1060  format (/,'Density (kg/m3) = ',f7.1,/, 
     *          'Enthalpy (J/kg) = ',g12.6,/,
     *          'Specific Enthalpy (J/m3) = ',g12.6,/,
     *          'Entropy (J/K/kg) = ',g12.6,/,
     *          'Specific Entropy (J/K/m3) = ',g12.6,/,
     *          'Heat Capacity (J/K/kg) = ',g12.6,/,
     *          'Specific Heat Capacity (J/K/m3) = ',g12.6,/,
     *          'Aggregate Grueneisen Ratio = ',g12.6,/,
     *          'Aggregate Adiabatic Bulk Modulus (bar) = ',g12.6,
     *        /,'Aggregate Shear Modulus (bar) = ',g12.6,/,
     *          'Aggregate Sound Velocity (km/s) = ',g12.6,/,
     *          'Aggregate P-Wave Velocity (km/s) = ',g12.6,/,
     *          'Aggregate S-Wave Velocity (km/s) = ',g12.6,/)
1070  format ('Variance (c-p+',i1,') = ',i2,//,40('-'),/)
1080  format (/,16x,'Complete Assemblage',15x,'Solid+Melt Only',
     *        /,14x,'mol',7x,' mol %',6x,'wt %',9x,' mol %',6x,'wt %')
1100  format (/,'Solid Density (kg/m3) = ',f7.1,/, 
     *          'Solid Enthalpy (J/kg) = ',g12.6,/,
     *          'Solid Secific Enthalpy (J/m3) (2) = ',g12.6,/,
     *          'Solid Entropy (J/K/kg) = ',g12.6,/,
     *          'Solid Specific Entropy (J/K/m3) = ',g12.6,/,
     *          'Solid Heat Capacity (J/K/kg) (1) = ',g12.6,/,
     *          'Solid Specific Heat Capacity (J/K/m3) (1) = ',g12.6,/,
     *          'Aggregate Solid Melt Grueneisen Ratio = ',g12.6,/
     *         ,'Aggregate Solid Adiabatic Bulk Modulus (bar) (2,3) = ',
     *  g12.6,/,'Aggregate Solid Shear Modulus (bar) = ',g12.6,/
     *         ,'Aggregate Solid Sound Velocity (km/s) = ',g12.6,/
     *         ,'Aggregate Solid P-Wave Velocity (km/s) = ',
     *  g12.6,/,'Aggregate Solid S-Wave Velocity (km/s) = ',g12.6,/)
1110  format (1x,a8,2x,f8.3,5x,2(f6.2,4x),5x,2(f6.2,4x))
1120  format (29x,a8,' = ',g12.6)
1160  format (/,'Molar Properties:'
     *        /,20x,'N(g)',8x,'H(J)',6x,'S(J/K)',6x,'V(J/bar)',6x,
     *         'Cp(J/K)'
     *         ,6x,'Alpha(1/K)',2x,'Beta(1/bar)',2x,'Density(kg/m3)')
1170  format (1x,a,1x,f9.2,3x,13(g12.5,1x))
1190  format (/,'Seismic (Molar) Properties:'
     *        /,17x,'Gruneisen',7x,'Ks(bar)',7x,'Mu(bar)',
     *        4x,'V0(km/s)',5x,'Vp(km/s)',5x,'Vs(km/s)')
1200  format (1x,a,3x,12(g12.5,1x))
1210  format (/,'Bulk Composition:')
1220  format (/,'Bulk Properties:')
1230  format (/,'N.B.: Aggregate properties represent the entire stable'
     *         ,' assemblage.',/,'Solid aggregate properties represent '
     *         ,'solid and melt properties,',/,'but do not include '
     *         ,'molecular fluid properties.',/)

      end 

      subroutine x3toy (id,ids)
c----------------------------------------------------------------------
c subroutine to convert geometric reciprocal solution compositions (x3(id,i,j))
c to geometric endmember fractions (y) for solution model ids. replicate 
c of subroutine xtoy, but for the x3 array (used only by getpar).
c----------------------------------------------------------------------
      implicit none 

      include 'perplex_parameters.h'
c                                 -------------------------------------
c                                 local variables:
      integer ids, l, m, ld, id
c                                 -------------------------------------
c                                 global variables:
c                                 bookkeeping variables
      integer ksmod, ksite, kmsol, knsp
      common/ cxt0  /ksmod(h9),ksite(h9),kmsol(h9,m4,mst),knsp(h9,m4)
c                                 working arrays
      double precision z, pa, p0a, x, w, y
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1)
c                                 x coordinate description
      integer istg, ispg, imlt, imdg
      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)

      double precision x3
      common/ cxt16 /x3(k21,mst,msp)
c----------------------------------------------------------------------

      do l = 1, mstot(ids)
c                                 the endmembers may have been
c                                 rearranged from the original order,
c                                 use knsp(ids,l) to assure correct
c                                 indexing
         ld = knsp(ids,l) 

         y(ld) = 1d0

         do m = 1, istg(ids)
            y(ld) = y(ld)*x3(id,m,kmsol(ids,ld,m))
         end do

      end do   

      end 
