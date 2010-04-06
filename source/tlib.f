c----------------------------------------------------------------------
 
c TLIB - a library of subprograms called by the PERPLEX programs.

c Copyright (c) 1998 by James A. D. Connolly, Institute for Mineralogy
c & Petrography, Swiss Federal Insitute of Technology, CH-8092 Zurich,
c SWITZERLAND. All rights reserved.
 
c Please do not distribute this source.

c----------------------------------------------------------------------

      subroutine redop1 (output,opname)
c----------------------------------------------------------------------
c redop1 - redop1 looks for the perplex_option.dat file, if it finds
c the option file it scans the file for keywords and sets options
c accordingly.

c option variables - keyword associations

c lopt(1)  - closed_c_space -> T = closed compositional variables
c lopt(2)  - set in getprp -> T = cumulative modes 
c lopt(3)  - hard_limits -> T = on
c lopt(4)  - helffrich murnaghan correction
c lopt(5)  - site_check -> T = reject invalid site fractions
c lopt(6)  - melt_is_fluid -> T = classify melts as fluids in output

c nopt(20) - T_melt - kill melt endmembers at T < nopt(20)
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ier, i, loopx, loopy

      logical output

      character*3 key*22, val, valu(i10), nval1*12, nval2*12,
     *            nval3*12,opname*100,strg*40,strg1*40

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      double precision bp1,bm1,bpm,lbpm
      common/ cst46 /bp1,bm1,bpm,lbpm

      integer grid
      double precision rid 
      common/ cst327 /grid(5,2),rid(2)

      integer jtest,jpot
      common/ debug /jtest,jpot

      integer isec,icopt,ifull,imsg,io3p
      common/ cst103 /isec,icopt,ifull,imsg,io3p
c----------------------------------------------------------------------
c                                 default option values:
c                                 nopt(1) and nopt(2) are flags used by reader
      nopt(1) = 0d0
c                                 closed or open compositional space
      lopt(1) = .true.
c                                 Anderson-Gruneisen correction
      lopt(4) = .true.
c                                 reject invalid site fractions
      lopt(5) = .true.
c                                 melt_is_fluid
      lopt(6) = .false.
c                                 minimum replicate label distance
      nopt(4) = 0.025
c
      iopt(2) = 0 
      valu(2) = 'mol'
      iopt(3) = 0 
      valu(3) = 'vol'
c                                 interpolation keyword
      iopt(4) = 0
      valu(4) = 'on '
c                                 extrapolation keyword
      iopt(5) = 0
      valu(5) = 'off'
c                                 vrh_weighting keyword
      nopt(6) = 0.5d0
c                                 bad_number keyword
      nopt(7) = 0d0
c                                 zero_mode (<0 off)
      nopt(9) = 1d-6
c                                 tolerance below which a component is considered to 
c                                 be zero during fractionation
      nopt(11) = 1d-6
c                                 number of iterations for
c                                 pseudocompound refinement
      iopt(10) = 4
c                                 number of increments used
c                                 in refinement
      iopt(11) = 3
c                                 max number of points to 
c                                 be refined in addition to 
c                                 active points
      iopt(12) = 10
c                                 quench temperature (K)
      nopt(12) = 0d0
c                                 initial resolution for adaptive 
c                                 refinement
      nopt(13) = 0.1d0
c                                 solvus_tolerance keyword
      nopt(8) = nopt(13)
c                                 reach for adaptive refinement
      nopt(14) = 0.9d0
c                                 perturbation to eliminate pseudocompound
c                                 degeneracies
      nopt(15) = 5d-3
c                                 poisson ratio to be use for calculating
c                                 missing shear moduli
      valu(15) = 'on '
      nopt(16) = 0.35d0
      iopt(16) = 1
c                                 stretch factor (b-1) for conformal
c                                 subdivision
      bm1 = 0.0164d0
c                                 subdivision model, 0 - solution model
c                                 1 - cartesian, 2 - stretch
      iopt(13) = 0 
      valu(13)  = 'off'
c                                 autorefine, 2 - automatic, 1 - manual, 0 - no
      iopt(6) = 1
      valu(6) = 'man'
c                                 increase in resolution for adaptive minimization 
      nopt(17) = 3d0 
c                                 increase in resolution for composition and mixed variable diagram calculations
      nopt(18) = 1d1
c                                 increase in resolution for Schreinemakers diagram calculations   
      nopt(19) = 3d0 
c                                 T_melt cutoff 
      nopt(20) = 873d0
c                                 fractional slop on autorefine limit
      nopt(3) = 1d-2
c                                 hard_limits for solution model refinement
      valu(16) = 'off'
      lopt(3) = .false.
c                                 compare local and max disorder state for o-d models
      valu(17) = 'on '
      iopt(17) = 1
c                                 assume linear boundaries within a cell during gridded minimization
      valu(18) = 'on '
      iopt(18) = 1
c                                 -------------------------------------
c                                 werami output options:

c                                 cumulative_modes
c     valu(15) = 'off'
c     iopt(15) = 0
c                                 -------------------------------------
c                                 for gridded minimization:
c                                 # nodes in i direction
      grid(1,1) = 20 
      grid(1,2) = 40
c                                 # nodes in j direction
      grid(2,1) = 20 
      grid(2,2) = 40
c                                 # of levels
      grid(3,1) = 1
      grid(3,2) = 4
c                                 1d fractionation path
      grid(4,1) = 20 
      grid(4,2) = 150
c                                 -------------------------------------
c                                 for schreinemakers etc:
c                                 max variance 
      grid(5,1) = 1
      grid(5,2) = 99
c                                 default increment (relative)
      rid(1) = 0.1d0
      rid(2) = 0.025d0
c                                 reaction format
      ifull = 0 
      valu(7) = 'min'
c                                 reaction list 
      jtest = 0 
      valu(9) = 'off'
c                                 console msgs
      imsg = 0
      valu(8) = 'on '
c                                 efficiency level
      isec = 3
c                                 short print
      io3p = 1
      valu(10) = 'on '
c                                 print dependent potentials
      jpot = 1
      valu(11) = 'off'
c                                 -------------------------------------
c                                 look for file
      open (n8, file = opname, iostat = ier, status = 'old')
c                                 no option file, set defaults
      if (ier.ne.0) then 
         write (*,1120) opname
      else 
         write (*,1130) opname
      end if
c                                 read cards to end of 
c                                 option file
      do while (ier.eq.0) 

         call redcd1 (ier,key,val,nval1,nval2,nval3,strg,strg1)
c                                 ier ne 0 bad record or eof
         if (ier.ne.0) exit 
c                                 if here we should have a keyword and
c                                 value
         if (key.eq.'composition') then 
c                                 phase composition key
            if (val.eq.'wt') then
               iopt(2) = 1
            else 
               iopt(2) = 0
            end if 
            valu(2) = val
         else if (key.eq.'proportions') then 
c                                 phase proportion key
            if (val.eq.'wt') then
               iopt(3) = 1
            else if (val.eq.'mol') then 
               iopt(3) = 2 
            else 
c                                 volume is default
               iopt(3) = 0
            end if 
            valu(3) = val
         else if (key.eq.'interpolation') then 
c                                 interpolation key
            if (val.eq.'off') iopt(4) = 0
            valu(4) = val
            if (val.eq.'on ') read (nval1,*) iopt(4)
         else if (key.eq.'extrapolation') then 
c                                 extrapolation key
            if (val.eq.'on ') then
               iopt(5) = 0
            else 
               iopt(5) = 1
            end if 
            valu(5) = val 
         else if (key.eq.'vrh_weighting') then 
c                                 vrh weighting key
            read (strg,*) nopt(6)
         else if (key.eq.'bad_number') then 
c                                 bad number key
            read (strg,*) nopt(7)
         else if (key.eq.'solvus_tolerance') then 

            read (strg,*) nopt(8)
         else if (key.eq.'zero_bulk') then
c                                 zero_bulk key
            read (strg,*) nopt(11)

         else if (key.eq.'zero_mode') then
c                                 zero_mode key
            read (strg,*) nopt(9)
         else if (key.eq.'iteration') then
c                                 max iteration key
            read (strg,*)  iopt(10)
            read (nval1,*) iopt(11)
            read (nval2,*) iopt(12)

         else if (key.eq.'initial_resolution') then
c                                 initial_resolution key 
            read (strg,*) nopt(13)

         else if (key.eq.'reach_factor') then 
c                                 reach_factor key
            read (strg,*) nopt(14)

         else if (key.eq.'stretch_factor') then
c                                 stretch_factor key = b - 1       
            read (strg,*) bm1

         else if (key.eq.'subdivision_override') then 
c                                 subdivision overide key
            valu(13) = val

            if (val.eq.'lin') then
               iopt(13) = 1
            else if (val.eq.'str') then
               iopt(13) = 2
            else 
               iopt(13) = 0 
            end if 
         else if (key.eq.'auto_refine') then
c                                 autorefine
            valu(6) = val

            if (val.eq.'off') then
               iopt(6) = 0
            else if (val.eq.'aut') then
               iopt(6) = 2
            end if 

            if (nval1.gt.'0') then 
               write (*,1100) nval1
               stop
            end if 

         else if (key.eq.'auto_refine_factor_I') then
c   
            read (strg,*) nopt(17)

         else if (key.eq.'auto_refine_factor_II') then
c   
            read (strg,*) nopt(18)

         else if (key.eq.'auto_refine_factor_III') then
c   
            read (strg,*) nopt(19)

         else if (key.eq.'auto_refine_slop') then
c   
            read (strg,*)  nopt(3)

         else if (key.eq.'x_nodes') then
c                                 number of x nodes at level 1 before autorefine
            read (strg,*) grid(1,1)
c                                 number of x nodes for autorefine
            read (nval1,*) grid(1,2)

         else if (key.eq.'y_nodes') then
c                                 number of y nodes at level 1
            read (strg,*) grid(2,1)
c                                 number of y nodes for autorefine
            read (nval1,*) grid(2,2)

         else if (key.eq.'grid_levels') then 
c                                 number of grid levels before autorefine
            read (strg,*) grid(3,1)
c                                 number of grid levels for autorefine
            read (nval1,*) grid(3,2)

         else if (key.eq.'1d_path') then 
c                                 number of grid points for 1d path before autorefine
            read (strg,*) grid(4,1)
c                                 number of grid points for 1d path for autorefine
            read (nval1,*) grid(4,2)

         else if (key.eq.'variance') then 
c                                 max variance of traced equilibria before autorefine
            read (strg,*) grid(5,1)
c                                 max variance of traced equilibria for autorefine
            read (nval1,*) grid(5,2)      

         else if (key.eq.'increment') then 
c                                 default exploratory relative increment    
            read (strg,*) rid(1)  
c                                 default autorefine relative increment
            read (nval1,*) rid(2)

         else if (key.eq.'reaction_format') then 

            valu(7) = val

            if (val.eq.'ful') then 
               ifull = 1
            else if (val.eq.'sto') then 
               ifull = 2
            else if (val.eq.'S+V') then 
               ifull = 3
            else if (val.eq.'eve') then
               ifull = 4
            else 
               valu(7) = 'min'
            end if 
  
         else if (key.eq.'console_messages') then 
            
            if (val.eq.'off') then 
               valu(8) = val
               imsg = 1
            else 
               valu(8) = 'on '
            end if 

         else if (key.eq.'reaction_list') then

            if (val.eq.'on ') then 
               valu(9) = val
               jtest = 3 
            else 
               valu(9) = 'off'
            end if

         else if (key.eq.'efficiency') then 

            read (strg,*) isec

            if (isec.lt.1.or.isec.gt.5) isec = 3 

         else if (key.eq.'short_print') then 

            if (val.eq.'off') then 
               io3p = 0
               valu(10) = 'off'
            end if 

         else if (key.eq.'dependent_potentials') then 

            if (val.eq.'on ') then 
               jpot = 0
               valu(11) = 'on '
            end if

         else if (key.eq.'replicate_label') then 
c                                 replicate lable tolerance for PSSECT
c                                 replicate fields are labled only if they
c                                 are further apart than the normalized distance 
c                                 specified by this keyword.   
            read (strg,*) nopt(4)

         else if (key.eq.'hard_limits') then 

            if (val.eq.'on ') then 
               lopt(3) = .true.
            else 
               valu(16) = 'off'
            end if

         else if (key.eq.'T_stop') then 
c                                 equilibrium cutoff T (K)    
            read (strg,*) nopt(12)

         else if (key.eq.'T_melt') then 
c                                 cutoff T (K) for melt endmember stability    
            read (strg,*) nopt(20)

         else if (key.eq.'order_check') then 
c                                 compare local and max disorder state for o-d models
            if (val.eq.'off') then 
               iopt(17) = 0
               valu(17) = 'off'
            end if 

         else if (key.eq.'linear_model') then   
c                                 assume linear boundaries within a cell during gridded minimization
            if (val.eq.'off') then 
               iopt(18) = 0
               valu(18) = 'off'
            end if 

         else if (key.eq.'closed_c_space') then
 
            if (val.ne.'T') lopt(1) = .false. 

         else if (key.eq.'Anderson-Gruneisen') then

            if (val.eq.'F') lopt(4) = .false.

         else if (key.eq.'site_check') then 

            if (val.eq.'F') lopt(5) = .false.

         else if (key.eq.'melt_is_fluid') then 

            if (val.eq.'T') lopt(6) = .true.

         else if (key.eq.'pc_perturbation') then
c                                 perturbation to eliminate pseudocompound degeneracies  
            read (strg,*) nopt(15)

         else if (key.eq.'poisson_ratio') then 
c                                 handle missing shear moduli
            if (val.eq.'on ') then
               read (nval1,*) nopt(16)
            else if (val.eq.'off') then 
               valu(15) = val
               iopt(16) = 0
            else if (val.eq.'all') then 
               read (nval1,*) nopt(16)
               valu(15) = val
               iopt(16) = 2
            end if             

         else if (key.ne.'|') then 

            write (*,1110) key

         end if 

      end do 
                
      close (n8)
c                                 -------------------------------------
c                                 error traps:

      if (iopt(10).gt.0) then 
         iopt(5) = 0
         valu(5) = 'off'
      end if 

      if (iopt(11).lt.3) then 
         write (*,1040)
         iopt(11) = 3
      end if 
c                                 reach factor
      if (2d0*nopt(14).gt.dfloat(iopt(11)).or.nopt(14).le.0d0) then
         write (*,1030)
         nopt(14) = 0.9d0 
      end if 
c                                 initial resolution
      if (nopt(13).ge.1d0.or.nopt(13).lt.0d0) then 
         write (*,1050)
         nopt(13) = 0.1d0
      end if 
c                                 stretching parameters
      if (bm1.lt.0d0) then 
         write (*,1060)
         bm1 = 0.0164
      end if 
c                                 auto-refine factor I
      if (nopt(17).lt.1d0) then 
         nopt(17) = 3d0
         write (*,1070) nopt(17)
      end if 
c                                 auto-refine factor II
      if (nopt(18).lt.1d0) then 
         nopt(18) = 1d1
         write (*,1070) nopt(18)
      end if 
c                                 auto-refine factor III
      if (nopt(19).lt.1d0) then 
         nopt(19) = 3d0
         write (*,1070) nopt(19)
      end if 
c                                 auto-refine slop
      if (nopt(3).lt.0d0.or.nopt(3).gt.1d0) then 
         nopt(3) = 0.01d0
         write (*,1080)
      end if 
c                                 grid parameters
      do i = 1, 2

         if (grid(3,i).le.0.or.grid(3,i).gt.l8) grid(3,i) = 4

         loopy = (grid(2,i)-1) * 2**(grid(3,i)-1) + 1
         loopx = (grid(1,i)-1) * 2**(grid(3,i)-1) + 1

         if (loopy.gt.l7) then 
            call warn (92,nopt(1),loopy,'y_node')
            grid(2,i) = (l7 - 1)/2**(grid(3,i)-1) + 1
         end if 

         if (loopx.gt.l7) then 
            call warn (92,nopt(1),loopx,'x_node')
            grid(1,i) = (l7 - 1)/2**(grid(3,i)-1) + 1
         end if 

         if (grid(4,i).gt.l7) then 
            call warn (92,nopt(1),loopx,'1dpath')
            grid(4,i) = l7 - 1
         end if  

         if (grid(5,i).lt.1) then 
            call warn (113,rid(i),grid(5,i),'VARIAN')
            grid(5,i) = 1
         end if 

         if (rid(i).lt.1d-2) then 
            call warn (114,rid(i),i,'INPUT1')
         end if 

      end do
c                                 stretching
      bp1 = 2d0 + bm1
      bpm = bp1/bm1
      lbpm = dlog(bpm)
c                                 compositional resolution = xinc*nopt(10)
      nopt(10) = 2d0*nopt(13)*nopt(14)*
     *           (nopt(14)/dfloat(iopt(11)))**iopt(10)
c                                 --------------------------------------
c                                 output
      if (output) then 
      
         write (*,1000) 

         write (*,1010) valu(11),nopt(7),
     *                  nopt(8),nopt(9),nopt(11),
     *                  iopt(10),iopt(11),iopt(12),k19,nopt(13)

         write (*,1011) bm1,
     *                  nopt(14),valu(13),valu(6),nopt(17),nopt(18),
     *                  nopt(19),nopt(3),valu(16),nopt(12),nopt(20),
     *                  valu(17),
     *                  lopt(1),nopt(15),lopt(4),lopt(5),lopt(6)
         
         write (*,1015) grid(1,1),grid(1,2),l7,
     *                  (grid(1,1)-1) * 2**(grid(3,1)-1) + 1,
     *                  (grid(1,2)-1) * 2**(grid(3,2)-1) + 1,
     *                  grid(2,1),grid(2,2),l7,
     *                  (grid(2,1)-1) * 2**(grid(3,1)-1) + 1,
     *                  (grid(2,2)-1) * 2**(grid(3,2)-1) + 1,
     *                  grid(3,1),grid(3,2),l8,
     *                  grid(4,1),grid(4,2),l7,
     *                  valu(18),
     *                  grid(5,1),grid(5,2),rid(1),rid(2), 
     *                  isec,valu(7),valu(9),valu(8),valu(10)
c                                 werami/meemum output options
         write (*,1016) (valu(i),i=2,5),nopt(6),valu(15),nopt(16)
c                                 resolution blurb
         if (nopt(13).gt.0d0) write (*,1090) 
     *                  'Exploratory',nopt(10)*1d2,
     *                  nopt(13)*1d2,nopt(13)*1d2,
     *                  'Auto-refine',nopt(10)*1d2/nopt(17),
     *                  nopt(13)/nopt(18)*1d2,
     *                  nopt(13)/nopt(19)*1d2
c                                 pssect/psvdraw options
         write (*,1017) nopt(4)

         write (*,1020) 

      end if
c                                 -------------------------------------
c                                 recalculate parameters:
c                                 proportionality constant for shear modulus
      nopt(16) = 1.5d0*(1d0-2d0*nopt(16))/(1d0+nopt(16))
c                                 the value of nopt(10) is multiplied to 
c                                 give the relaxation bounds.
      nopt(10) = 1d1*nopt(10)
      if (nopt(10).gt.nopt(13)) nopt(10) = nopt(13)

1000  format (/,'Perple_X options are currently set as:',//,
     *      '    Keyword:               Value:     Permitted values ',
     *          '[default]:')
1010  format (4x,'dependent_potentials   ',a3,8x,'off [on]',/,
     *        4x,'bad_number             ',f7.0,4x,'[0.0]',/,
     *        4x,'solvus_tolerance       ',f5.3,6x,
     *       '0->1 [0.05], 0 => p=c pseudocompounds, 1 => homogenize'
     *     ,/,4x,'zero_mode              ',e7.1,4x,
     *           '0->1 [1e-6], < 0 => off'
     *     ,/,4x,'zero_bulk              ',e7.1,4x,
     *           '0->1 [1e-6], < 0 => off'/,
     *        4x,'adaptive optimization ',/,
     *        4x,'   iteration limit     ',i2,9x,'>0 [3], value 1 of ',
     *                                          'iteration keyword',/,
     *        4x,'   refinement factor   ',i2,9x,
     *           '>2*r [5], value 2 of ',
     *                                          'iteration keyword',/,
     *        4x,'   refinement points   ',i2,9x,'0->',i3,
     *           ' [10], value 3'
     *                                       ,' of iteration keyword',/,
     *        4x,'initial_resolution     ',f4.2,7x,
     *           '0->1 [0.1], 0 => off')
1011  format (4x,'stretch_factor         ',f5.3,6x,'>0 [0.0164]',/,
     *        4x,'reach_factor           ',f3.1,8x,'>0.5 [0.9]',/,
     *        4x,'subdivision_override   ',a3,8x,'[off] lin str',/,
     *        4x,'auto_refine            ',a3,8x,'off [manual] auto',/,
     *        4x,'auto_refine_factor_I   ',f4.2,7x,'>=1 [3]',/
     *        4x,'auto_refine_factor_II  ',f4.2,7x,'>=1 [10]',/,
     *        4x,'auto_refine_factor_III ',f4.2,7x,'>=1 [3]',/,
     *        4x,'auto_refine_slop       ',f5.3,6x,'[0.01]'/,
     *        4x,'hard_limits            ',a3,8x,'[off] on',/,
     *        4x,'T_stop (K)             ',f6.1,5x,'[0]',/,
     *        4x,'T_melt (K)             ',f6.1,5x,'[873]',/,
     *        4x,'order_check            ',a3,8x,'off [on]',/,
     *        4x,'closed_c_space         ',l1,10x,'F [T]',/,
     *        4x,'pc_perturbation        ',f6.4,5x,'[5d-3]',/,
     *        4x,'Anderson-Gruneisen     ',l1,10x,'[T] F',/,
     *        4x,'site_check             ',l1,10x,'[T] F',/,
     *        4x,'melt_is_fluid          ',l1,10x,'[F] T',/)
1015  format (4x,'Gridded minimization parameters:',/,
     *        4x,'  x_nodes             ',i3,' /',i3,4x,'[20/40], >0, '
     *          ,'<',i4,'; effective x-resolution ',i4,' /',i4
     *          ,' nodes',/
     *        4x,'  y_nodes             ',i3,' /',i3,4x,'[20/40], >0, '
     *          ,'<',i4,'; effective y-resolution ',i4,' /',i4,
     *           ' nodes',/
     *        4x,'  grid_levels           ',i1,' /',i2,5x,'[1/4], >0, '
     *          ,'<',i2,/,
     *        4x,'  1d_path             ',i3,' /',i3,4x,
     *           '[20/150], >0, <',i4,/,
     *        4x,'  linear_model           ',a3,6x,'off [on]',//,
     *        4x,'Parameters for Schreinemakers and Mixed-variable ',
     *           'diagram calculations:',/,
     *        4x,'  variance             ',i2,' /',i2,5x,
     *           '[1/99], >0, maximum true variance',/,
     *        4x,'  increment         ',f5.3,' /',f5.3,2x,
     *           '[0.1/0.025], ',
     *           'default search/trace variable increment',/,
     *        4x,'  efficiency             ',i1,8x,'[3] >0 < 6',/,       
     *        4x,'  reaction_format      ',a3,8x,'[min] ',
     *           'full stoichiometry S+V everything',/,
     *        4x,'  reaction_list        ',a3,8x,'[off] on',/,
     *        4x,'  console_messages     ',a3,8x,'[on] off',/,
     *        4x,'  short_print_file     ',a3,8x,'[on] off',/) 
1016  format (4x,'WERAMI/MEEMUM output options:',/,
     *        4x,'  compositions           ',a3,8x,'wt  [mol]',/,
     *        4x,'  proportions            ',a3,8x,'wt  [vol] mol',/,
     *        4x,'  interpolation          ',a3,8x,'off [on ]',/,
     *        4x,'  extrapolation          ',a3,8x,'on  [off]',/,
     *        4x,'  vrh_weighting          ',f3.1,8x,'0->1 [0.5]',/,
     *        4x,'  poisson_ratio          ',a3,8x,'off [on ] all, ',
     *        'Poisson ratio = ',f4.2,/)
1017  format (4x,'PSSECT/PSVDRAW plot options:',/,
     *        4x,'  replicate_labels       ',f5.3,6x,'0->1 [0.025]',/)
1020  format ('To change these options edit or create ',
     *        'a computational option file',/,'See: ',
     *        'www.perplex.ethz.ch/perplex_options.html',/)
1030  format (/,'Warning: the reach_factor cannot exceed value2 of the',
     *         ' iteration keyword.',/,'the keyword will be',
     *         ' assigned its default value.',/)
1040  format (/,'Warning: value2 of the iteration keyword must be ',
     *         ' >= 3',/,'value2 will be',
     *         ' assigned its default value.',/)
1050  format (/,'Warning: initial_resolution keyword must be ',
     *         '< 1',/,'the keyword will be',
     *         ' assigned its default value.',/)
1060  format (/,'Warning: the stretch_factor cannot be less than zero',
     *        /,' the keyword will be  assigned its default value',/)
1070  format (/,'Warning: auto_refine_factors must be ',
     *         '> 1',/,'the keyword will be',
     *         ' assigned its default value (',i2,').',/)
1080  format (/,'Warning: auto_refine_slop must be ',
     *         ' >=0 and <1',/,'the keyword will be',
     *         ' assigned its default value.',/)
1090  format (
     *       'Worst case (Cartesian) compositional resolution (mol %)'
     *       ,/,3x,
     *       'Stage         Adaptive   Non-Adaptive*',
     *       '   Non-Adaptive**',/,
     *       (3x,a11,3x,f7.3,6x,f7.3,9x,f7.3),/,
     *       (3x,a11,3x,f7.3,6x,f7.3,9x,f7.3),/,
     *       '* - composition/mixed variable diagrams, ** - ',
     *        'Schreinemakers diagrams',/)
1100  format (/,'Error: data (',a
     *       ,') follows the auto_refine keyword value '
     *       ,'most probably',/,'you are using an obsolete copy of ',
     *        'perplex_option.dat recopy or edit the file.',/)
1110  format (/,'Warning: unrecognized option text: ',a,/,
     *       'If the text is intentional, check spelling and case.',/) 
1120  format (/,'Warning: the Perple_X option file: ',a,/,
     *       'was not found, default option values will be used.',/) 
1130  format (/,'Reading computational options from: ',a)

      end 

      subroutine redcd1 (ier,key,val,nval1,nval2,nval3,strg,strg1)
c----------------------------------------------------------------------
c this routine seeks a card containing a keyword and as many as 
c six values (char variables nval, nval1...), the first 3 letters of nval
c is also returned in val, if the second value is longer than 12 characters
c it is also saved in the character variable strg*80
c----------------------------------------------------------------------    
      implicit none

      include 'perplex_parameters.h'

      integer len, ier, iscan, i, iscnlt, ibeg, iend, ist

      character chars(240)*1, card*240, key*22, val*3,
     *          nval1*12, nval2*12, nval3*12, strg*40, strg1*40

      integer  iop0 
      common / basic /iop0

      logical  debug
      common / debugblk / debug

      ier = 0 
      key = ' '

      if (debug) PRINT *,'IN REDCD1'
      do 

         read (n8,'(a)',end=90) card

         if (debug) PRINT *, card

         if (card.ne.' ') then 
            read (card,'(240a)') chars
c                                 find end of data marker '|'
            len = iscan (1,240,chars,'|') - 1

            if (len.eq.0) cycle 
c                                 find a non blank character
            ibeg = iscnlt (1,len,chars,' ')

            exit 

         end if 

      end do 
c                                 find end of keyword 
      iend = ibeg + 1
      iend = iscan (iend,240,chars,' ') - 1
c                                 load chars into key
      write (key,'(22a1)') (chars(i), i = ibeg, iend)

      iend = iend + 1
c                                 now locate the value:
      ibeg = iscnlt (iend,len,chars,' ')
c                                 now find trailing blank
      iend = iscan (ibeg,240,chars,' ') 
c                                 save longer versions (only on first value)
c                                 this is done in case it's long text or 
c                                 several numbers on certain options. 
      strg = ' '
      strg1 = ' '
      if (iend-ibeg.gt.39) iend = ibeg+39
      write (strg,'(40a1)') (chars(i), i = ibeg, iend)
      write (strg1,'(40a1)') (chars(i), i = ibeg, ibeg+39)

      write (val,'(3a1)') (chars(i), i = ibeg, ibeg + 2)
c                                 look for a second value
      ist = iscan (ibeg,240,chars,' ')
c                                 find a non blank character
      ibeg = iscnlt (ist,len,chars,' ')
c                                 if no blank exit
      if (ibeg.gt.len) then 
         nval1 = '0'
      else 
         iend = iscan (ibeg,len,chars,' ')
         if (iend-ibeg.gt.11) iend = ibeg + 11 
         nval1 = ' '
         write (nval1,'(12a1)') (chars(i), i = ibeg, iend)
      end if 
c                                 look for a third value
      ist = iscan (ibeg,240,chars,' ')
c                                 find a non blank character
      ibeg = iscnlt (ist,len,chars,' ')
c                                 if no blank exit
      if (ibeg.gt.len) then 
         nval2 = '0'
      else 
         iend = iscan (ibeg,len,chars,' ')
         if (iend-ibeg.gt.11) iend = ibeg + 11 
         nval2 = ' '
         write (nval2,'(12a1)') (chars(i), i = ibeg, iend)
      end if 
c                                 look for a fourth value
      ist = iscan (ibeg,240,chars,' ')
c                                 find a non blank character
      ibeg = iscnlt (ist,len,chars,' ')
c                                 if no blank exit
      if (ibeg.gt.len) then 
         nval3 = '0'
      else 
         iend = iscan (ibeg,len,chars,' ')
         if (iend-ibeg.gt.11) iend = ibeg + 11 
         nval3 = ' '
         write (nval3,'(12a1)') (chars(i), i = ibeg, iend)
      end if 

      goto 99

90    ier = 1

99    end


      subroutine rdstrg (lun,nstrg,string,eof)
c----------------------------------------------------------------------
c rdstrg - read 240 column card images from unit lun until a non-blank
c (i.e., with data other than comments) record, then read up to three
c strings from the record. on output nstrg is the number of strings read
c from the record. 
c----------------------------------------------------------------------    
      implicit none

      include 'perplex_parameters.h'

      integer len, lun, iscan, i, iscnlt, ibeg, iend, ier, nstrg, imax

      logical eof

      character chars(240)*1, card*240, string(3)*8

      eof = .false.

      do 
c                                 read cards till a non-blank
c                                 card.
         read (lun,'(a)',iostat=ier) card

         if (ier.ne.0) then
c                                 error on read = eof
            eof = .true.

            return

         else if (card.ne.' ') then 

            read (card,'(240a)') chars
c                                 find end of data marker '|'
            len = iscan (1,240,chars,'|') - 1

            if (len.eq.0) cycle 
c                                 find a non blank character
            ibeg = iscnlt (1,len,chars,' ')

            exit 

         end if 

      end do 
c 
c                                 we have a non-blank card
      nstrg = 1

      do 
c                                 find the end of the string
         iend = iscan (ibeg,240,chars,' ') - 1

         if (iend-ibeg.gt.7) then

            imax = ibeg + 7 

         else

            imax = iend
 
         end if 
c                                 load chars into string
         write (string(nstrg),'(8a1)') (chars(i), i = ibeg, imax)
c                                 find the next string
         ibeg = iscnlt (iend+1,len,chars,' ')

         if (ibeg.gt.len.or.nstrg.eq.3) return
 
         nstrg = nstrg + 1

      end do 

      end

      subroutine eohead (n)
c----------------------------------------------------------------------
c eohead reads cards from n until an 'END ' or 'end ' is found in
c the first 3 columns
c----------------------------------------------------------------------
      implicit none

      integer n

      character tag*4

      rewind n

10    read (n,'(a)',end=20) tag
      if (tag.eq.'end'.or.tag.eq.'END') goto 99
      goto 10
20    call error (37,1d0,n,'EOHEAD')

99    end
 
      subroutine topn2 (iopt)
c----------------------------------------------------------------------
c topn2 reads the header of the thermodynamic data file, if iopt = 1
c then data base choice is known (idbase), else if 1 asks console for a 
c choice else if > 3 echos data except components and atwts to n8
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'
 
      character*5 tcname, tag*11, n2name*100, vname*8, xname*8, 
     *            string*140

      integer icod, ig(3), iopt, i, j

      double precision sum

      double precision delt,dtol,utol,ptol
      common/ cst87 /delt(l2),dtol,utol,ptol

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5  /p,t,xco2,u1,u2,tr,pr,r,ps

      integer io3,io4,io9
      common/ cst41 /io3,io4,io9

      double precision ctrans
      integer ictr,itrans
      common/ cst207 /ctrans(k0,k5),ictr(k5),itrans

      integer iff,idss,ifug,ifyn,isyn
      common/ cst10 /iff(2),idss(h5),ifug,ifyn,isyn

      integer idh2o,idco2,ikind,icmpn,icout
      double precision therm,comp,tot
      common/ cst43 /therm(k4),comp(k0),tot,icout(k0),idh2o,idco2,
     *               ikind,icmpn
  
      common/ csta2 /xname(k5),vname(l2)/ cst41a /n2name

      character cmpnt*5, dname*40
      common/ csta5 /dname,cmpnt(k0)

      common/ csta9 /tcname(k5)

      double precision atwt
      common/ cst45 /atwt(k0)
c-----------------------------------------------------------------------
      rewind n2
c                               a program other than vertex is reading
c                               the header if iopt = 1
      if (iopt.eq.1) itrans = 0
c                               read the number of data bases
c                               obsolete since always 1:
      read (n2,*,err=90) i
c                               read the extrinsic variable names:
      read (n2,1000,err=90) (vname(i), i = 1, 3)
      read (n2,1010,err=90) (ig(i), i = 1, 3)

      if (iopt.lt.4) then 
         if ((ifug.ge.10.and.ifug.le.12).or.
     *       ifug.eq.15.or.ifug.eq.17.or.ifug.eq.18) then 
            vname(3) = ' X(O) ' 
         else if (ifug.eq.25) then 
            vname(3) = 'Y(CO2)*'
         else if (ifug.eq.14.or.ifug.eq.13) then 
            vname(3) = 'X(H2)'
         end if 
      end if  
c                               read delt the finite difference
c                               increments for v1, v2, and v3. and dtol
c                               utol and ptol, the critical tolerances
c                               (in energy units) for determination of
c                               the stability of divariant, univariant
c                               and invariant equilibria or reactions.
      read (n2,*,err=90) delt, dtol
c                               dtol must be less than zero
c                               utol must be smaller than -utol
c                               ptol must be > 2*-dtol
      utol = -dtol/1d1
      ptol = -dtol*3d0
c                               read obsolete data base code and title
c                               and the data base reference state
c                               conditions consistent with v1 and v2.
      read (n2,*,err=90) icod,pr,tr
      read (n2,1170,err=90) dname   
c                               read number of data base comps 
      read (n2,*,err=90) icmpn
c                               read component names.
      read (n2,1180,err=90) (cmpnt(i), i = 1, icmpn)
c                               component atomic wts.
      read (n2,*,err=90) (atwt(i), i = 1, icmpn)
c                               read pointers for water and co2 cmpnt.  
c                               if data file doesn't contain water
c                               or co2 dummy values
c                               must be included in the file (ne. 0).
      read (n2,*,end=90) idh2o, idco2

      if (iopt.eq.3.or.iopt.eq.5) then 
c                                 get transformations if build or ctransf
         call gettrn (iopt)

      else 
c                                 substitute transformed component names
c                                 and compute the new formula wieghts
         do i = 1, itrans

            cmpnt(ictr(i)) = tcname(i)

            sum = 0d0

            do j = 1, icmpn
               sum = sum + ctrans(j,i) * atwt(j)
            end do
 
            atwt(ictr(i)) = sum

         end do 

      end if 

      if (iopt.gt.3) then 
c                                 echo formatted header data for ctransf/actcor:
         write (n8,9010) i
         write (n8,1000) (vname(i), i=1, 3)
         write (n8,1010) (ig(i), i=1, 3)
         write (n8,9020) delt, dtol, utol, ptol
         write (n8,9030) icod, pr, tr
         write (n8,1170) dname
         write (n8,9010) icmpn
         write (n8,1180) (cmpnt(i), i = 1, icmpn)
         write (n8,1030) (atwt(i), i = 1, icmpn)
         write (n8,9010) idh2o, idco2

      end if 
c                                 read and echo unformatted comments and make data                            
      do 

         read (n2,'(a)',end=90) string
         read (string,'(a)') tag
         if (iopt.gt.3) write (n8,'(a)') string

         if (tag.eq.'begin_makes') then
 
            call rmakes (iopt)

            cycle 

         else if (tag.ne.'end'.and.tag.ne.'END') then     

            cycle

         else 

            exit       

         end if

      end do  

      goto 99

90    call error (21,r,i,n2name)

1000  format (3(a8,18x))
1010  format (i2)
1030  format (6(g12.6,1x))
1170  format (a40)
1180  format (6(a5,1X)/6(a5,1X))
9010  format (10(i2,1x))
9020  format (8(g9.2,1x))
9030  format (i2,1x,8(g12.6,1x))

99    t = tr
      p = pr

      end 

      subroutine gettrn (iopt)
c----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i,iopt,ict,ier
 
      character*5 pname, rname, y*1

      double precision sum

      integer ictr, itrans
      double precision ctrans
      common/ cst207 /ctrans(k0,k5),ictr(k5),itrans

      character*5 cmpnt,dname*40
      common/ csta5 /dname,cmpnt(k0)

      double precision atwt
      common/ cst45 /atwt(k0)

      integer idh2o,idco2,ikind,icmpn,icout
      double precision therm,comp,tot
      common/ cst43 /therm(k4),comp(k0),tot,icout(k0),idh2o,idco2,
     *               ikind,icmpn
c-----------------------------------------------------------------------
c                                 recombine components:
10    write (*,1030)
      write (*,1040) (cmpnt(i), i = 1, icmpn)
      write (*,1050)
      read (*,3000) y
      if (y.ne.'Y'.and.y.ne.'y') goto 99

      write (*,1060)
      read (*,3000) pname
      if (pname.eq.'     ') goto 99
c                                 get the identity of the real comp
c                                 to be replaced.
50    write (*,1070) pname
      read (*,3000) rname

      do i = 1, icmpn
         if (rname.eq.cmpnt(i)) then
            if (iopt.eq.3) then
               if (i.eq.idh2o.or.i.eq.idco2) then
c                                 don't allow build users 
c                                 to transform saturated
c                                 phase components
                   call warn (14,atwt(1),i,cmpnt(i))
                   goto 60
                end if 
             end if

             icout(1) = i
             goto 70
         end if
      end do 
 
60    write (*,1080) 
      write (*,1040) (cmpnt(i), i = 1, icmpn)
      goto 50
c                                 get the identities of the other 
c                                 components in the new component:
70    itrans = itrans + 1
      ict = 1
      if (itrans.gt.k5) call error (999,atwt(1),ict,'GETTRN')
       
      write (*,4050) k5-1,pname
30    read (*,3000) rname
      if (rname.eq.'     ') goto 80

      do i = 1, icmpn
         if (rname.eq.cmpnt(i)) then 
            ict = ict + 1
            icout(ict) = i
            goto 30
         end if
      end do 
c                                 no match, try again message
      write (*,2300)
      goto 30
c                                 get the component stoichiometries:
80    write (*,4030) (cmpnt(icout(i)),i=1,ict)
      write (*,4040) pname
90    read (*,*,iostat=ier) (ctrans(icout(i),itrans), i= 1, ict)
      call rerror (ier,*90)
 
      write (*,1100) pname,(ctrans(icout(i),itrans),
     *                      cmpnt(icout(i)), i = 1, ict)
      write (*,1110)
      read (*,3000) y
 
      if (y.eq.'y'.or.y.eq.'Y') then
         sum = 0d0
         do i = 1, ict
            sum = sum + ctrans(icout(i),itrans) * atwt(icout(i))
         end do 
         atwt(icout(1)) = sum
         cmpnt(icout(1)) = pname
         ictr(itrans) = icout(1)
      else
         itrans = itrans - 1
         write (*,*) 'Try again.'
      end if

      goto 10

1030  format (/,'The current data base components are:')
1040  format (12(1x,a5))
1050  format ('Transform them (Y/N)? ')
1060  format ('Enter new component name, < 6 characters,',
     *          ' left justified: ')
1070  format ('Enter old component to be replaced',
     *          ' with ',a5,': ')
1080  format ('Select the component from the set: ')
1100  format (1x,a5,' = ',6(f6.2,1x,a5),/,9x,6(f6.2,1x,a5))
1110  format ('Is this correct (Y/N)? ')
2300  format (/,'You made a mistake, try again.',/
     *          'Check spelling and upper/lower case matches.',/)
3000  format (a)
4030  format ('Enter stoichiometric coefficients of:',/,
     *        2x,12(a5,1x))
4040  format ('in ',a5,' (in above order): ')
4050  format ('Enter other components (< ',i2,') in ',a,' 1 per',
     *        ' line, <cr> to finish:')

99    end

      subroutine rerror (ier,*)
c---------------------------------------------------------------------
c rerror - routine to check for errors during list directed i/o
 
      implicit none

      integer ier
c---------------------------------------------------------------------
 
      if (ier.eq.0) then
         return
      else
         write (*,1000)
         ier = 0
         return 1
      end if
 
1000  format (/,' Your input is incorrect, probably you are using ',
     *        'a character where',/,' you should be using a number ',
     *        'or vice versa, try again...',/)
 
      end

      subroutine error (ier,realv,int,char)
c---------------------------------------------------------------------
c write error messages and terminate execution
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ier, int
 
      character char*(*)

      double precision realv

      PRINT *,"In error.",ier, realv, int, char
      if (ier.eq.1.or.ier.eq.2) then 
         write (*,1) char,int
      else if (ier.eq.3) then 
         write (*,3)
      else if (ier.eq.4) then
         write (*,4) char 
      else if (ier.eq.5) then 
         write (*,5) int,char,j3
      else if (ier.eq.6) then 
         write (*,6) char
      else if (ier.eq.7) then 
         write (*,7) char
      else if (ier.eq.13) then
         write (*,13) h8
      else if (ier.eq.14) then
         write (*,14) char
      else if (ier.eq.15) then
         write (*,15) char
      else if (ier.eq.16) then
         write (*,16) h5
      else if (ier.eq.17) then
         write (*,17) int
      else if (ier.eq.18) then
         write (*,18) char
      else if (ier.eq.19) then
         write (*,19) char
      else if (ier.eq.20) then
         write (*,20) int, char
      else if (ier.eq.21) then
         write (*,21) char 
      else if (ier.eq.22) then
         write (*,22) int, char
      else if (ier.eq.23) then
         write (*,23) char
      else if (ier.eq.24) then
         write (*,24) int
      else if (ier.eq.25) then
         write (*,25) h9
      else if (ier.eq.26) then
         write (*,26) int, char
      else if (ier.eq.27) then
         write (*,27) 
      else if (ier.eq.28) then 
         write (*,28) int, char
      else if (ier.eq.29) then 
         write (*,29) int, char
      else if (ier.eq.30) then
         write (*,30) int,char
      else if (ier.eq.31) then
         write (*,31) char
      else if (ier.eq.32) then 
         write (*,32)
      else if (ier.eq.33) then 
         write (*,33) char, int
      else if (ier.eq.34) then
         write (*,34)
      else if (ier.eq.35) then
         write (*,35)
      else if (ier.eq.36) then
         write (*,36)
      else if (ier.eq.37) then
         write (*,37) int
      else if (ier.eq.38) then 
         write (*,38) 
      else if (ier.eq.39) then
         write (*,39) int
      else if (ier.eq.40) then
         write (*,40) int, char
      else if (ier.eq.41) then
         write (*,41) int, char
      else if (ier.eq.42) then
         write (*,42)
      else if (ier.eq.43) then
         write (*,43) char
      else if (ier.eq.44) then 
         write (*,44) 
      else if (ier.eq.45) then 
         write (*,45) 
      else if (ier.eq.46) then 
         write (*,46) realv, char
      else if (ier.eq.47) then 
         write (*,47) char
      else if (ier.eq.48) then 
         write (*,48) char,int
      else if (ier.eq.49) then 
         write (*,49) char,int
      else if (ier.eq.50) then 
         write (*,50) realv, char, int
      else if (ier.eq.51) then
         write (*,51) char
      else if (ier.eq.52) then 
         write (*,52) h9
      else if (ier.eq.53) then 
         write (*,53) 
      else if (ier.eq.54) then 
         write (*,54)
      else if (ier.eq.55) then 
         write (*,55) k16
      else if (ier.eq.56) then 
         write (*,56) k17
      else if (ier.eq.57) then 
         write (*,57) char
      else if (ier.eq.58) then 
         write (*,58) k21, char
      else if (ier.eq.59) then 
         write (*,59) k20, char
      else if (ier.eq.60) then 
         write (*,60) k22, char
      else if (ier.eq.61) then 
         write (*,61) k18, char
      else if (ier.eq.62) then
         write (*,62) char, int, realv
      else if (ier.eq.63) then 
         write (*,63) 
      else if (ier.eq.64) then
         write (*,64) char
      else if (ier.eq.65) then
         write (*,65) char
      else if (ier.eq.66) then
         write (*,66) char
      else if (ier.eq.67) then
         write (*,67) char
      else if (ier.eq.89) then
         write (*,89) 
      else if (ier.eq.90) then
         write (*,90) l6
      else if (ier.eq.106) then
         write (*,106) char
      else if (ier.eq.107) then
         write (*,107) int
      else if (ier.eq.108) then
         write (*,108) int
      else if (ier.eq.109) then
         write (*,109) int
      else if (ier.eq.110) then
         write (*,110)
      else if (ier.eq.111) then
         write (*,111)
      else if (ier.eq.112) then
         write (*,112)
      else if (ier.eq.116) then
         write (*,116)
      else if (ier.eq.117) then
         write (*,117)
      else if (ier.eq.118) then
         write (*,118)
      else if (ier.eq.120) then
         write (*,120) char
      else if (ier.eq.122) then
         write (*,122) char
      else if (ier.eq.125) then 
         write (*,125) realv, char
      else if (ier.eq.169) then
         write (*,169) int
      else if (ier.eq.180) then
         write (*,180) int,char
      else if (ier.eq.181) then
         write (*,181) int
      else if (ier.eq.182) then
         write (*,182) k2
      else if (ier.eq.183) then
         write (*,183) k2,char
      else if (ier.eq.197) then
         write (*,197)
      else if (ier.eq.200) then
         write (*,200)
      else if (ier.eq.204) then
         write (*,204) int
      else if (ier.eq.206) then
         write (*,206) int
      else if (ier.eq.207) then
         write (*,207) realv,char
      else if (ier.eq.208) then
         write (*,208) char
      else if (ier.eq.227) then
         write (*,227) char, int
      else if (ier.eq.228) then 
         write (*,228) char, realv, int
      else if (ier.eq.279) then
         write (*,279) int
      else if (ier.eq.323) then
         write (*,323)
      else
         write (*,999) ier,realv,int,char
      end if
 
      stop

1     format (/'**error ver001** increase parameter ',a,' to ',i7,' in',
     *       ' perplex_parameters.h and recompile Perple_X',/)
3     format (/,'**error ver003** the solution model file appears to ',
     *        'be in a format that is no longer supported.',/,
     *        'copy the current version from: ',
     *        'www.perplex.ethz.ch/datafiles/solut_08.dat',/)
4     format (/,'**error ver004** you must use ',a,' to analyze this ',
     *        'type of calculation.',/)
5     format (/,'**error ver005** too many ordered species (',i2,') in',
     *        ' solution model ',a,/,'increase dimension j3 (',i2,')',/)
6     format (/,'**error ver006** fractionation path coordinate file: '
     *          ,a,/,'does not exist.',/)
7     format (/,'**error ver007** reference phase ',a,' is not in the ',
     *          'thermodynamic data file.',/)
13    format ('**error ver013** too many excluded phases, ',
     *        'increase dimension h8 (',i3,')')
14    format ('**error ver014** programming error, routine ',a)
15    format (/,'**error ver015** missing composant for: ',a,/)
16    format (/,'**error ver016** too many saturated components, ',
     *        'increase dimension h5 (',i2,')')
17    format (/,'**error ver017** too many composants for a saturation',
     *        ' constraint increase dimension h6 (',i3,')')
18    format (/,'**error ver018** ',a,' is defined as a saturated ',
     *        'phase component in the thermodynamic data file.')
19    format (/,'**error ver019** probable cause missing composant,',
     *        ' executing routine ',a)
20    format (/,'**error ver020** error ',i2,' reading solution model',
     *        ' file.',/,'   Reading model: ',a,' Check format.',/)
21    format (/,'**error ver021**error reading ',
     *        'header section of',/,'thermodynamic data ',
     *        'file:',/,a,/,'Check formatting',/)
22    format (/,'**error ver022** too many divariant assemblages, ',
     *        'increase dimension j9 (',i8,') routine: ',a)
23    format (/,'**error ver023**error reading',
     *        ' thermodynamic data file.',/,' Last data read',
     *        ' without error was for: ',a,' Check formatting.',/)
24    format (/,'**error ver024** too many solution models in',
     *        ' solution model file',/,' increase parameter i9 (',
     *        i3,')')
25    format (/,'**error ver025** too many solution models ',
     *        ' increase parameter h9 (',i3,')')
26    format (/,'**error ver026** the number of fixed components (',
     *        i2,') in ',a,/,' is >= the number of components ',/)
27    format (/,'**error ver027** Error reading the computational',
     *        ' option file.',//,'Probable cause: You are using an ',
     *        'input file created by an out-of-date',/,
     *        '                version of BUILD, or you have',
     *        ' incorrectly edited the',/'                input file',
     *        ' created by BUILD',/)
28    format (/,'**error ver028** invalid buffer choice (',i3,') in',
     *          ' routine: ',a,/)
29    format (/,'**error ver029** unknown term type ',i6,' for',
     *          ' solution model: ',a,/)
30    format (/,'**error ver030** the number of mixing sites ',i2,
     *          ' is < the number of independent sites',/,' for',
     *          ' solution model: ',a,/)
31    format (/,'**error ver031** erroneous solution model',
     *          ' parameter for: ',a,/)
32    format (/,'**error ver032** stability field calculations (',
     *          'option 2) are disabled in this version of PERPLEX',/)
33    format (/,'**error ver033** expression with too many terms in ',a
     *       ,/,'increase m0 or j6 to',i2,'.',/)
34    format (/,'**error ver034** vmax is lt vmin, check input.')
35    format (/,'**error ver035** dv is lt 0, check input.')
36    format (/,'**error ver036** missing composant for the saturated',
     *        ' phase,',/,'you have probably excluded either H2O or',
     *        ' CO2,',/,'or a composant is duplicated in the',
     *        ' thermodynamic data file',/)
37    format (/,'**error ver037** no end marker in header',/,
     *        'section of thermodynamic data file unit ',i2,/)
38    format (/,'**error ver038** you have configured a ',
     *       'problem with only one independent variable.',/,
     *       'This case cannot be handled by constrained minimization',
     *       ' use the unconstrained computational mode.'/)
39    format (/,'**error ver039** too many end-members, ',
     *        'increase dimension k12 (',i2,') Routine: ',a)
40    format (/,'**error ver040** too many compositional coordinates, ',
     *        'increase dimension k13 (',i7,')  Routine: ',a)
41    format (/,'**error ver041** too many pseudocompounds, ',
     *        'increase dimension k1 or k21 (',i7,') Routine: ',a)
42    format (/,'**error ver042** the possible phases of the system do'
     *         ,' not span the specified bulk composition.',/,
     *        'To avoid this problem add phases or modify the bulk ',
     *        'composition.',/)
43    format (/,'**error ver043** you cannot simultaneously treat: ',
     *          a,/,'as a thermodynamic solution and as a saturated',
     *          ' phase.',/)
44    format (/,'**error ver044** too many saturated phase components.'
     *        /)
45    format (/,'**error ver045** too many mobile components.'/)
46    format (/,'**error ver046** temperature (',g12.6,' K) is out of ',
     *        'range for the equation of state (',a,'),',/,
     *        'most likely this problem can be corrected by setting ',
     *        'the Anderson_Gruneisen keyword to TRUE',/,
     *        'in the Perple_X option file.')
47    format (/,'**error ver047** solution model ',a,' is incorrectly ',
     *        'formatted (van Laar).',/)
48    format (/,'**error ver048** too many terms in solution model ',a,
     *        ' increase parameter m1 (',i2,').',/)
49    format (/,'**error ver049** the order of solution model ',a,
     *        ' is too high, increase parameter m2 (',i2,').',/)
50    format (/,'**error ver050** requested compositional resolution ',
     *          '(',f6.0,') for a component in',
     *          'solution model: ',a,/,'exceeds 1/MRES (MRES=',i4,') ',
     *          'reduce requested resolution or inrease parameter',/,
     *          'MRES in routine CARTES',/)
51    format (/,'**error ver051** DUMMY1 could not find the auxilliary'
     *         ,' input file:',/,a,/,'required for open system model ',
     *          'computations (ICOPT=9).',/)
52    format (/,'**error ver052** too many solution models in your'
     *         ,' calculation',/,'reduce the number of models or ',
     *          'increase parameter H9 (',i2,').',/)
53    format (/,'**error ver053** phase fractionation calculations '
     *         ,'require >1 thermodynamic component.',/)
54    format (/,'**error ver054** unanticipated condition, probable ',
     *          'cause is incorrect ordering of',/,'endmembers in the',
     *          ' solution model, which leads to inconsistent site ',
     *          'occupancies',/)
55    format (/,'**error ver055** too many make definitions,'
     *         ,'increase parameter K16 (',i2,') and recompile.',/)
56    format (/,'**error ver056** too many phases in a make definition'
     *         ,', increase parameter K17 (',i2,') and recompile.',/)
57    format (/,'**error ver057** failed on an accepted make definition'
     *         ,' for ',a,/,'routine INPUT2'/)
58    format (/,'**error ver058** too many pseudocompounds generated ',
     *     'by refinement, increase dimension k21 (',i8,') routine: ',a)
59    format (/,'**error ver059** too many coordinates generated by ',
     *        'refinement, increase dimension k20 (',i8,') routine: ',a)
60    format (/,'**error ver060** too many coordinates generated by ',
     *        'refinement, increase dimension k22 (',i8,') routine: ',a)
61    format (/,'**error ver061** too many solution coordinates, ',
     *        'increase dimension k18 (',i8,') routine: ',a)
62    format (/,'**error ver062** solution model ',a,' specifies non-',
     *          'Cartesian subdivision (',i1,')',/,' and must be refor',
     *          'mulated for adapative minimization, but VERTEX cannot',
     *        /,'do the reformulation because the initial_reolution ',
     *          'keyword specified in',/,'perplex_option.dat (',f5.2,
     *          ') is invalid',/)
63    format (/,'**error ver063** inconsistent auto-refine data.',
     *        ' Suppress or reinitialize auto-refinement.',/) 
64    format (/,'**error ver064** PSVDRAW plots only ',
     *          'binary mixed-variable and ',/,
     *          'ternary composition diagrams (',a,').',/)
65    format (/,'**error ver065** dimensioning error (',a,').',/)
66    format (/,'**error ver066** invalid format, most probably this',
     *          'result should be plotted with PSSECT (',a,').',/)
67    format (/,'**error ver067** file ',a,/,
     *        'is not formatted correctly for PSVDRAW.',/)
89    format (/,'**error ver089** SMPLX programming error. Change ',
     *        'minimnization method.',/)
90    format (/,'**error ver090** SMPLX failed to converge within ', 
     *        i6,' iterations.',/,'Probable cause: the possible ',
     *        'phases do not span the systems composition',/,
     *        'To avoid this problem add phases or modify the bulk ',
     *        'composition.',/,'Alternatively, although less ',
     *        'probably, increasing parameter L6 in perplex_',
     *        'parameters.h',/,
     *        'and recompiling VERTEX permit SMPLEX to converge.',/)
106   format (/,'**error ver106** programming error in ',a)
107   format (/,'**error ver107** the assemblage you input is ',
     *        'metastable (ier=',i3,').')
108   format (/,'**error ver108** the assemblage you input is ',
     *        'degenerate (ier=',i3,').')
109   format (/,'**error ver109** the assemblage you input does not '
     *       ,'span the specified bulk composition (ier=',i3,').')
110   format (/,'**error ver110** you have requested a calculation ',
     *        'with the composition ',/,'of a saturated phase as a ',
     *        'variable, but you have not defined its composition')
111   format (/,'**error ver111** you have requested a calculation ',
     *        'with the composition ',/,'of a saturated phase as a ',
     *        'variable, but the phase has only one component')
112   format (/,'**error ver112** the maximum value of an independent '
     *       ,'variable',/,'is less than or equal to the minimum value')
116   format (/,'**error ver116** an independent variable, or at least'
     *       ,' its name, is undefined')
117   format (/,'**error ver117** vmax(iv(3)) ne vmin(iv(3) but no ',
     *        'sectioning variable v(iv(3)) is defined')
118   format (/,'**error ver118** the default increment of the ',
     *        'sectioning variable will result ',/,
     *        'in the generation of more ',
     *        'than 10 sections, to avoid this',/,' error increase ',
     *        'the increment or modify this test')
120   format (/,'**error ver120** file:',/,a,/,
     *        'could not be opened, check that it exists or that it is',
     *        ' not in use by another program.',/) 
122   format ('**error ver120** plot file: ',a,/,'was not found, ',
     *        'you must generate it with VERTEX.')
125   format (/,'**error ver125** a site fraction (',g8.2,') is out',
     *          ' of range for : ',a,/,'   The configurational',
     *          ' entropy model is probably incorrect.',/)
169   format (/,'**error ver169** cart, imod=',i2,' is an invalid ',
     *          'request')
180   format (/,'**error ver180** too many pseudocompounds,increase ',
     *           'parameter k1 or k21 or k10 (',i7,') Routine: ',a)
181   format (/,'**error ver181** too many reactions,',
     *           ' increase dimension k2 (',i6,')')
182   format (/,'**error ver182** too many invariant points,',
     *           ' increase parameter k2 (',i6,')')
183   format (/,'**error ver183** too many assemblages; increase ',
     *        ' parameter k2 (',i6,'), routine ',a)
197   format (/,'**error ver197** to many components, increase',
     *        ' parameter k5',/)
200   format (/,'**error ver200** you are trying to use a fluid ',
     *        'equation of state for an invalid component',/)
204   format ('**error ver204** too many stable assemblages, i',
     *        'ncrease dimension j9 (',i8,')',/)
206   format ('**error ver206** too many univariant assemblages ',
     *        'intersect the edges of the diagram, i',
     *        'ncrease dimension k2 (',i6,')',/)
207   format (/,'**error ver207** the value of the stretching ',
     *          ' parameter (',g13.6,')',/,'for solution ',a,
     *          ' is invalid (<1) for transform subdivision,',/,
     *          'check section 4 of PERPLEX documentation.',/)
208   format (/,'**error ver208** too many phases on one side of a',/
     *        ' reaction.',/,'Do not use the full reaction',
     *        ' equation option (',a,').')
227   format (/,'**error ver227** in solution model ',a,' a DQF ',
     *          'correction is specified for endmember: ',i2,/,
     *          'DQF corrections can only be made on the idependent ',
     *          'endmembers of a solution model',/)
228   format (/,'**error ver228** in solution model ',a,' negative ',
     *          'composition (',g12.6,') for component ',i2,/,
     *          'probable cause is an incorrect stoichiometric ',
     *          'definition for a dependent endmember.',/)
279   format (/,'**error ver279** you have got a problem ',i3)
323   format (/,'**error ver323** prime9, imd(i)=0 is the only',/,
     *        'subdivision scheme permitted for this version' )
999   format (/,'**error vertex** unspecified error ier=',i3,
     *        ' real=',g13.6,' i=',i7,' char=',a6)
      end

      subroutine warn (ier,realv,int,char)
c---------------------------------------------------------------------
c write warning messages and terminate execution
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'
 
      integer ier,int

      double precision realv

      integer io3,io4,io9
      common/ cst41 /io3,io4,io9
 
      character char*(*)

      if (ier.eq.1) then 
         write (*,1) 
      else if (ier.eq.5) then
         write (*,5) 
      else if (ier.eq.6) then
         write (*,6) 
      else if (ier.eq.7) then
         write (*,7) 
      else if (ier.eq.8) then
         write (*,8) h8
      else if (ier.eq.9) then
         write (*,9) char
      else if (ier.eq.10) then
         write (*,10) int, realv, char
      else if (ier.eq.11) then
         write (*,11) char
      else if (ier.eq.12) then
         write (*,12) char
      else if (ier.eq.13) then
         write (*,13) char
      else if (ier.eq.14) then
         write (*,14) char
      else if (ier.eq.15) then
         write (*,15)
      else if (ier.eq.16) then
         write (*,16)
      else if (ier.eq.17) then
         write (*,17)
      else if (ier.eq.18) then
         write (*,18) realv
      else if (ier.eq.19) then
         write (*,19) 
      else if (ier.eq.20) then
         write (*,20)
      else if (ier.eq.21) then
         write (*,21) realv, char
      else if (ier.eq.22) then
         write (*,22) realv, char
      else if (ier.eq.23) then
         write (*,23) char     
      else if (ier.eq.24) then
         write (*,24) realv
      else if (ier.eq.25) then 
         write (*,25) int, char
      else if (ier.eq.26) then 
         write (*,26) char
      else if (ier.eq.28) then
         write (*,28)
      else if (ier.eq.29) then
         write (*,29) char
      else if (ier.eq.30) then
         write (*,30) char
      else if (ier.eq.31) then 
         write (*,31)
      else if (ier.eq.32) then
         write (*,32) char
      else if (ier.eq.33) then
         write (*,33) char
      else if (ier.eq.34) then
         write (*,34) char
      else if (ier.eq.35) then
         write (*,35) char,realv
      else if (ier.eq.36) then 
         write (*,36) realv, char 
      else if (ier.eq.37) then
         write (*,37)  
      else if (ier.eq.38) then
         write (*,38) 
      else if (ier.eq.39) then
         write (*,39) 
      else if (ier.eq.40) then
         write (*,40) 
      else if (ier.eq.41) then
         write (*,41) char
      else if (ier.eq.42) then
         write (*,42)     
      else if (ier.eq.43) then
         write (*,43) int, char
      else if (ier.eq.44) then
         write (*,44) char
      else if (ier.eq.45) then
         write (*,45) char
      else if (ier.eq.46) then 
         write (*,46) realv, char
      else if (ier.eq.47) then
         write (*,47) int, realv
      else if (ier.eq.48) then 
         write (*,48) 
      else if (ier.eq.50) then
         write (*,50) char
      else if (ier.eq.51) then 
         write (*,51) char
      else if (ier.eq.52) then 
         write (*,52) char
      else if (ier.eq.58) then
         write (*,58)
      else if (ier.eq.60) then
         write (*,60) char
      else if (ier.eq.63) then
         write (*,63)
      else if (ier.eq.73) then
         write (*,73) char, realv, int
      else if (ier.eq.74) then
         write (*,74)
      else if (ier.eq.79) then
         write (*,79) char
      else if (ier.eq.87) then
         write (*,87)
      else if (ier.eq.88) then
         write (*,88)
      else if (ier.eq.89) then
         write (*,89)
      else if (ier.eq.90) then
         write (*,90) 
      else if (ier.eq.91) then
         write (*,91)
      else if (ier.eq.92) then 
         write (*,92) int, l7, char
      else if (ier.eq.106) then
         write (*,106) char
      else if (ier.eq.108) then
         write (*,108)
      else if (ier.eq.109) then
         write (*,109)
      else if (ier.eq.113) then
         write (*,113) int
      else if (ier.eq.114) then
         write (*,114)
      else if (ier.eq.172) then
         write (*,172) 
      else if (ier.eq.173) then
         write (*,173) 
      else if (ier.eq.175) then
         write (*,175) char,ier,realv
      else if (ier.eq.176) then
         write (*,176) char
      else if (ier.eq.190) then
         write (*,190) int
      else if (ier.eq.205) then
         write (*,205) int
         write (*,900)
      else if (ier.eq.207) then
         write (*,207) int
      else if (ier.eq.208) then
         write (*,208) int
      else if (ier.eq.209) then
         write (*,209) int
      else
         write (*,999) ier, char, realv, int
      end if

1     format (/,'**warning ver001** the amount of a phase is < 0, this',
     *        ' usually indicates that',/,'the specified amount of a ',
     *        'saturated component is inadequate to saturate the system'
     *        ,/)
5     format (/,'**warning ver005** fluid components are specified',
     *        ' as thermodynamic AND as either',/,'saturated phase',   
     *      ' or saturated components; almost certainly a BAD idea.',/)
6     format (/,'**warning ver006** fluid components are specified',
     *        ' as both thermodynamic AND',/,'saturated ',    
     *        'components; almost certainly a BAD idea.',/)
7     format (/,'**warning ver007** fluid components are specified as'
     *       ,' a saturated phase component',/,'AND as a thermodynamic', 
     *        'or saturated component; almost certainly a BAD idea.',/)
8     format ('**warning ver08** too exclude more phases ',
     *        'increase paramter h8 (',i3,')')
9     format ('**warning ver009** unable to deconstruct transition,'
     *       ,/,'data for ',a,' will not be output.')
10    format (/,'**warning ver010** not able to traverse ',
     *          'the entire',/,'extent of equilibrium (',i6,')',
     *          ' at v(3)=',g12.6,/,
     *          'this error can usually be avoided by increasing the ',
     *          'finite',/,'difference increment delt(iv(1)) or delv',
     *          '(iv(2)), as defined on card 6 of',/,'the file on n2.',
     *          ' In routine:',a,/)
11    format (/,'**warning ver011** ',a,' has > 1',
     *          ' transition with dp/dT ne 0 and may not be treated ',/,
     *          ' correctly')
12    format (/,'**warning ver012** ',a,' has a transition ',
     *          ' with dp/dT < 0 and may not be treated ',/,
     *          ' correctly')
13    format (/,'**warning ver013** ',a,' has null or negative ',
     *          'composition and will be ',/,'rejected from the ',
     *          'composition space.')
14    format (/,'**warning ver014** You can not redefine the ',
     *          'saturated phase component:',a,/,'To circumvent this ',
     *          'restriction use CTRANSF to make a data base with the',/
     *         ,'the desired component transformations',/)
15    format (/,'**warning ver015** if you select > 1 saturated ',
     *          'component, then the order you',/,'enter the ',
     *          'components determines the saturation heirarchy and may'
     *          ,' effect your',/,'results (see Connolly 1990).',/)
16    format (/,'**warning ver016** you are going to treat a saturate',
     *        'd (fluid) phase component',/,'as a thermodynamic ',
     *        'component, this may not be what you want to do.',/)
17    format (/,'**warning ver017** you gotta be kidding, either ',
     *          ' 1 or 2 components, try again:',/)
18    format (/,'**warning ver018** the value of the default dependen',
     *         't variable (',g14.6,') for the following',/,
     *         'equilibrium was inconsistent with the an earlier ',
     *         'determination of the invariant condition',/,
     *         'and will be reset. This may cause the curve to ',
     *         'look kinked near the invariant point',/)
19    format (/,'**warning ver019** you must specify at least ',
     *        'one thermodynamic component, try again',/)
20    format ('**warning ver020** sfol2')
21    format ('**warning ver021** xmax (',g12.6,') > 1 for '
     *         ,' solution model ',a,/,' xmax will be reset to 1',
     *        /,' see documentation, section 4.0.')
22    format ('**warning ver022** xmin (',g12.6,') < 0 for '
     *         ,' solution model ',a,/,' xmin will be reset to 1',
     *        /,' see documentation, section 4.0.')
23    format ('**warning ver023** xmin > xmax for solution ',a,/,
     *        'xmin will be set to xmax NO PSEUDOCOMPOUNDS WILL BE',
     *        ' GENERATED.',/,'see documentation, section 4.0',/)
24    format (/,'**warning ver024** wway, increment refined out of',
     *          ' range (',g8.1,')',/,' before the stable',
     *          ' extension of the equilibria was located')
25    format ('**warning ver025** ',i1,' endmembers for ',a,
     *          ' The solution will not be considered.')
26    format ('**warning ver026** only one endmember for ',a,
     *          ' The solution will not be considered.')
28    format (/,'**warning ver028** minimization failed, ill-',
     *        'conditioned?',/)
29    format ('**warning ver029** programming error, routine ',a,/)
30    format (/,'**warning ver030** Because of missing endmembers, ',
     *        'or that the',/,'subdivision',
     *        ' scheme specified for solution model ',a,/,'is too',
     *        ' restrictive, there are no valid compositions for', 
     *        ' this model.',/)
31    format (/,'**warning ver031** this choice is disabled because ',
     *        'the dependent_potentials',/,'keyword is missing or off',
     *        ' in perplex_option.dat, to use this choice set the',/,
     *        'keyword to on and re-run VERTEX.',/)
32    format ('**warning ver032** fixed activity option requested',
     *          ' for ',a,/,'This option is disabled, the',
     *          ' solution will not be considered.')
33    format ('**warning ver033** missing endmembers for ',a,/,
     *        'The model may be recast in > 1 way for',
     *        ' the endmember subset.',/,'To control this choice',
     *        ' eliminate undesired endmembers.')
34    format ('**warning ver034** ',a,' could not be recast as',
     *          ' a simpler model.',/,'The solution will not be',
     *          ' considered. Add the missing endmembers or eliminate'
     *          ,/,'additional endmembers to allow this model.',/)
35    format (/,'**warning ver035** ',a,' is only for pure fluids',
     *        /,' XCO2 will be reset to: ',f4.2,/)
36    format ('**warning ver021** xinc (',g12.6,') < 0 for'
     *         ,' solution model ',a,/,'xinc will be reset to 1.'
     *         ,' see documentation, section 4.0. ',/)
37    format (/,'**warning ver37** you will not be able to plot the ',
     *       'results of this',/,'calculation with PSVDRAW. PSVDRAW ',
     *       'only plots ternary composition diagrams.',/)
38    format (/,'**warning ver38** you will not be able to plot the ',
     *       'results of this',/,'calculation with PSVDRAW. PSVDRAW ',
     *       'only plots mixed-variable diagrams for',/,'a binary ',
     *       'system with one independent potential variable.',/)
39    format (/,'**warning ver39** PSVDRAW will plot the results of ',
     *       'this calculation as a',/,'projected section, such plots ',
     *       'may be difficult to interpret. To plot',/,
     *       'pseudosections as a an explicit function of a systems ',
     *       'composition use',/, 'gridded minimization.',/)
40    format (/,'**warning ver040** you have configured a ',
     *       'problem with only one independent variable.',/)
41    format (/,'**warning ver041** icky pseudocompound names'
     *       ,' for solution model: ',a,/,'refer to pseudocompound_'
     *       ,'glossary.dat file for pseudocompound definitions.',/)
42    format (/,'**warning ver042** the possible phases of the system',
     *        ' do not span the specified bulk composition.',/,3x,
     *        'To avoid this problem add phases or modify the bulk ',
     *        'composition.',/)
43    format (/,'**warning ver043** ',i2,' solutions referenced ',
     *          'in your input',/,'were not found in the solution ',
     *          'model file, routine:',a,/)
44    format ('**warning ver044** a solution model has destabilized',
     *        ' the endmember: ',a,' (iend=2).')
45    format (/,'**warning ver045** the entity involves ',
     *        ' phases (',a,' )',/,'described by a nonlinear EoS (see',
     *        ' program documentation Eq 2.2)',/,
     *        ' NO OUTPUT WILL BE GENERATED FOR THIS ENTITY.',/)
46    format (/,'**warning ver046** temperature (',g12.6,' K) is out',
     *        ' of range for the equation of state (',a,').',/,
     *        'most likely this problem can be corrected by setting ',
     *        'the Anderson_Gruneisen keyword to TRUE',/,
     *        'in the Perple_X option file.')
47    format (/,'**warning ver047** invariant point ',i6,' could ',
     *        ' not be located within',/,' the specified tolerance ',
     *        ' (PTOL= ',g12.6,' ) reset PTOL to avoid this problem.',/)
48    format (/,'**warning ver048** fluid phase pseudocompound data ',
     *         'does not include',/,' volumetric properties (SWASH).',/)
c49    format (/,'**warning ver049** some pseudocompound data has not',
c     *          ' been output because',/' the bulk modulus pressure ',
c     *          'derivative is not constant for all endmembers ',/,
c     *          ' (SWASH, see program documentation Eq 2.3)',/)
50    format (/,'**warning ver050** reformulating reciprocal ',
     *          'solution: ',a,' because of missing endmembers. ',
     *        /,'(reformulation can be controlled explicitly ',
     *          'by excluding additional endmembers).',/)
51    format (/,'**warning ver051** cannot make ',a,' because of ',
     *          'missing data or'
     *       ,/,'an invalid definition in the thermodynamic data file.'
     *       ,/)
52    format (/,'**warning ver052** rejecting ',a,'; excluded or '
     *       ,'invalid composition.',/)
58    format (/,'**warning ver058** wway, the equilibrium of the '
     *         ,'following reaction',/,' is inconsistent with the ',
     *          'invariant equilibrium.',/)
60    format (/,'**warning ver060** non-fatal programming error ',
     *          'routine:',a,/)
63    format (/,'**warning ver063** wway, invariant point on an edge?',
     *        /)
73    format (/,'**warning ver073** an invariant point has been',
     *          ' skipped, routine:',a,/,' This problem typically',
     *          ' occurs because two phases in the thermodynamic data',
     *          ' file have identical properties.',/,' Otherwise ',
     *          ' decreasing DTOL (',g9.3,') in the thermodynamic', 
     *          ' data file (Doc Sect 3, Card 6) for  variable ',i1,/,
     *          ' may eliminate this problem',/)
74    format (/,'**warning ver074** no new equilibria identified,',
     *          ' if degenerate segments have',/,' been skipped',
     *          ' increase the computational reliability level.',/)
79    format (/,'**warning ver079** univeq failed on an edge for ',
     *          'the following equilibrium.',/,' Probable cause is ',
     *          'extreme independent variable limits (e.g., xco2=0)',/
     *          ' or poor convergence criteria ',
     *          'in the thermodynamic data file. In routine:',a,/)
87    format (/,'**warning ver087** wway-univeq did not converge ',
     *          'when div was refined',/)
88    format (/,'**warning ver088** SMPLX converged to a non-unique ',
     *        'solution.',/,3x,'Probable cause: system composition ',
     *        'coincides with that of ',//,3x,'a compound or a ',
     *        'tieline between compounds.',//,3x,'This may lead to ',
     *        'erratic results.',//,3x,'To avoid this problem ',
     *        'perturb the systems composition.',/)
89    format (//,'**warning ver089** BUILD you did not request',
     *        'plot file output.',/,' You will not be able to process',
     *        ' the results of the requested calculation.',//)
90    format (/,'**warning ver090** SMPLX failed to converge. ',/,
     *        3x,'Probable cause: the possible ',
     *        'phases do not span the systems composition',/,3x,
     *        'To avoid this problem add phases or modify the bulk ',
     *        'composition.',/,3x,'Alternatively, although less ',
     *        'probably, increasing parameter L6 in perplex_',
     *        'parameters.h',/,3x,
     *        'and recompiling VERTEX permit SMPLEX to converge.',/)
91    format (/,'**warning ver091** SMPLX programming error. Change ',
     *        'minimnization method',/)
92    format (/,'**warning ver092** you have requested ',i4,
     *        ' grid points. Current',/,'dimensioning is for ',
     *        i4,' points. To obtain the requested resolution',/,
     *        'increase parameter L7 and recompile; or reduce the ',
     *        'required resolution via',/,'the ',a,' keyword in ',
     *        'perplex_option.dat',/)
106   format ('**warning ver106** programming error in ',a)
108   format (/,'**warning ver108** wway, a phase field with the '
     *         ,'following',/,' reaction is stable on both ',
     *          'sides of an invariant point',/,' this error can ',
     *          'usually be avoided by increasing the finite ',/,
     *          ' difference increment delt(iv(1)) or delv',
     *          '(iv(2)), defined',
     *           /,' on card 6 of the thermodynamic data file',/)
109   format (/,'**warning ver109** you may ',
     *        'have assigned a mobile component as an independent ',/,
     *        ' variable without defining the component',/)
113   format (/,'**warning ver113** maximum variance for equilibrium',
     *        ' tracing [the variance keyword in ',/,
     *        'perplex_option.dat] must be > 0, but is ',i2,
     *        '. Set to 1 for the current calculation',/)
114   format (/,'**warning ver114** the default increment of an ',
     *        'independent variable is <',/,' 1 percent of ',
     *        'its range, this is usually inefficient',/)
172   format (/,'**warning ver172** you cannot use this equation of',
     *          ' state with Y(CO2)',/,' as an indepedent variable, ',
     *          ' pick another one:',/)
173   format (/,'**warning ver173** invalid buffer choice ',/)
175   format (/,'**warning ver175** speciation routine ',a,' did',
     *          ' not converge ',/,' possibly due to graphite super-',
     *          'saturation. ier = ',i1,' real = ',g16.5,/)
176   format ('**warning ver176** fluid equation of state routine ',
     *        a,' did not converge.',/,' If execution continues this ',
     *        'may lead to incorrect results.',/,' To avoid this ',
     *        'problem choose a different equation of state.')
190   format (/,'**warning ver190** SMPLX failed to converge within ', 
     *        i6,' iterations.',/,3x,'Probable cause: the possible ',
     *        'phases do not span the systems composition',/,3x,
     *        'To avoid this problem add phases or modify the bulk ',
     *        'composition.',/,3x,'Alternatively, although less ',
     *        'probably, increasing parameter L6 in perplex_',
     *        'parameters.h',/,3x,
     *        'and recompiling VERTEX permit SMPLEX to converge.',/)
205   format (/,'**error ver205** too many new phase assemblages, ',
     *        'found by routine newhld',/,' increase dimension j9 (',
     *        i8,')',/)
207   format ('**warning ver207** the assemblage you input is ',
     *        'metastable (ier=',i3,').')
208   format ('**warning ver208** the assemblage you input is ',
     *        'degenerate (ier=',i3,').')
209   format ('**warning ver209** the assemblage you input does not'
     *       ,' span the specified bulk composition (ier=',i3,').')
900   format (' the calculation may be incomplete !!!!',/)
999   format (/,'**warning unspecified** ier =',i3,' routine ',a6
     *         ,' r = ',g12.6,' int = ',i9,/)
      end

      subroutine rmakes (iopt)
c----------------------------------------------------------------------
c rmakes is called by topn2 to read make definitions of thermodynamic
c entities, these entities are defined as a linear combination of 
c exisiting entities as defined in the thermodynamic file, with an 
c optional pressure/temperature dependent DQF correction. The format
c assumes data on one line of less than 240 characters, the expected format
c is

c name = num_1 * name_1 + num_2 * name_2 ....+ num_int * name_int
c dnum_1 dnum_2 dnum_3

c where i num_j is a number or fraction (i.e., two numbers separated by a 
c '/') and name_j is the name of the int existing entities. 
c and the dqf correction to the entity 'name' is
c Gdqf(J/mol) = dnum_1 + T[K]*dnum_2 + P[bar]*dnum_3

c end_of_data is either a "|" or the end of the record.

c make definitions are preceeded by the keyword:

c begin_makes 

c and truncated by the keyword:

c end_makes

c if iopt > 3, data is echoed to LUN n8 (for ctransf/actcor).
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ibeg, iend, len, ier, iscan, i, nreact, iopt

      double precision rnum

      integer io3,io4,io9
      common / cst41 /io3,io4,io9

      character chars(240)*1, tname*8, name*8, rec*240, tag*3

      double precision mcomp
      character mknam*8
      integer nmak
      logical mksat
      common / cst333 /mcomp(k16,k0),nmak,mksat(k16),mknam(k16,k17)

      double precision mkcoef, mdqf
      integer mknum, mkind
      common / cst334 /mkcoef(k16,k17),mdqf(k16,k17),mkind(k16,k17),
     *                 mknum(k16)

      integer ixct,iexyn,ifact
      common/ cst37 /ixct,iexyn,ifact 

      character*8 exname,afname
      common/ cst36 /exname(h8),afname(2)
c----------------------------------------------------------------------

      call readcd (n2,len,ier,chars)
      if (ier.ne.0) goto 90 
c                                 echo data for ctransf/actcor
      if (iopt.gt.3) write (n8,'(240a1)') (chars(i),i=1,len)

      iend = len
      nmak = 0 

      write (rec,'(240a1)') chars
      read (rec,'(a3)') tag

      do while (tag.ne.'end')   

         nmak = nmak + 1
         if (nmak.gt.k16) call error (55,mkcoef(1,1),nmak,'RMAKES')
c                                 get first name
         ibeg = 1
         call readnm (ibeg,iend,len,ier,tname,chars)
         if (ier.ne.0) goto 90
c                                 find start of data marker '='
         ibeg = iscan (1,len,chars,'=') + 1
c                                 the rest of the data should
c                                 consist of coefficients followed
c                                 by names
         nreact = 0 

         do while (ibeg.lt.len) 
c                                 find the number
            call readfr (rnum,ibeg,iend,len,ier,chars)
            if (ier.eq.2) then 
c                                 ier = 2 = a read error
               goto 90
            else if (ier.eq.1) then 
c                                 ier = 1, end-of-definition
               exit 
            end if 
c                                 find the name
            call readnm (ibeg,iend,len,ier,name,chars)
            if (ier.ne.0) goto 90

            nreact = nreact + 1
            if (nreact.gt.k17) call error (56,mkcoef(1,1),nmak,'RMAKES')

            mkcoef(nmak,nreact) = rnum 
            mknam(nmak,nreact) = name
           
         end do

         if (nreact+1.gt.k17) call error (56,mkcoef(1,1),nmak,'RMAKES')
         mknam(nmak,nreact+1) = tname
         mknum(nmak) = nreact
c                                 now the dqf
         call readcd (n2,len,ier,chars)
         if (ier.ne.0) goto 90
c                                 echo data for ctransf/actcor 
         if (iopt.gt.3) write (n8,'(240a1)') (chars(i),i=1,len)

         write (rec,'(240a1)') chars
         read (rec,*) mdqf(nmak,1),mdqf(nmak,2),mdqf(nmak,3)
c                                 start next make definition
         call readcd (n2,len,ier,chars)
         write (rec,'(240a1)') chars
         read (rec,'(a3)') tag
c                                 echo data for ctransf/actcor
         if (iopt.gt.3) write (n8,'(240a1)') (chars(i),i=1,len)

c                                 reject excluded makes
         do i = 1, ixct
            if (tname.eq.exname(i)) then 
               nreact = nreact - 1
               exit 
            end if 
         end do 

      end do 

      goto 99

90    write (*,1000) (chars(i),i=1,len)
      stop
      
1000  format (/,'**error ver200** READMK bad make definition in the',
     *        ' thermodynamic data file',/,'currently reading: ',/
     *        ,240a)

99    end 

      subroutine readnm (ibeg,iend,len,ier,name,chars)
c----------------------------------------------------------------------
c readnm looks for the first word in a record chars, ibeg is the index
c of the 1st letter, iend is the index of the last letter.
c-----------------------------------------------------------------------
      implicit none

      integer ibeg, iend, len, iscan, iscnlt, ier, i, imax

      character chars(240)*1, name*8

      ier = 0 
c                                 find start of name
      ibeg = iscnlt (ibeg,len,chars,' ') 
c                                 find next blank
      iend = iscan (ibeg,len,chars,' ') - 1

      imax = iend - ibeg
c                                 initialize to be safe:
      name = '        '

      if (imax.le.7) then

         write (name,'(8a1)') (chars(i),i=ibeg,iend)

      else 
c                                 can't be a valid name, save it
c                                 anyway in case it's a tag
         write (name,'(8a1)') (chars(i),i=ibeg,ibeg+7)
         ier = 4

      end if 

      ibeg = iend + 1

      end 

      subroutine readcd (nloc,len,ier,chars)
c----------------------------------------------------------------------
c readcd - read 240 column card image from unit 9, strip out unwanted
c characters. ier = 1 no card found.
c----------------------------------------------------------------------    
      implicit none

      integer len, ier, iscan, ict, i, iscnlt, ibeg, nloc

      character chars(240)*1, card*240

      ier = 0 

      ibeg = 0
  
      len = 0 

      card = ' '

      do while (ibeg.ge.len) 

         read (nloc,'(a)',end=90) card

         if (card.ne.' ') then 
            read (card,'(240a)') chars
c                                 find end of data marker '|'
            len = iscan (1,240,chars,'|') - 1
c                                 find a non blank character
            ibeg = iscnlt (1,len,chars,' ')
         end if 

      end do 

      ict = 1

      do i = 2, len 
c                                 strip out '+' and '*' chars
         if (chars(i).eq.'+'.or.chars(i).eq.'*') chars(i) = ' '
c                                 eliminate blanks after '/' and '-'
c                                 and double blanks
         if ((chars(ict).eq.'/'.and.chars(i  ).ne.' ') .or. 
     *       (chars(ict).eq.'-'.and.chars(i  ).ne.' ') .or.
     *       (chars(ict).eq.' '.and.chars(i  ).ne.' ') .or.
     *       (chars(ict).ne.'-'.and.chars(ict).ne.'/'.and.
     *        chars(ict).ne.' ') ) then
             ict = ict + 1
             chars(ict) = chars(i)
         end if

      end do 

      len = ict

      goto 99

90    ier = 3

99    end

      subroutine readfr (rnum,ibeg,iend,len,ier,chars)
c----------------------------------------------------------------------
c readfr looks for a number or two numbers separated by a backslash / in
c that array chars(iend:ibeg), the latter case is interpreted as a ratio. 
c the result is returned as num
c-----------------------------------------------------------------------
      implicit none

      double precision rnum, rnum1 

      integer ibeg, iend, len, iback, ier, iscan, iscnlt, i

      character chars(240)*1, num*30

      ier = 0 
c                                 now find start of a number
      ibeg = iscnlt (ibeg,len,chars,' ')  
c                                 find backslash
      iback = iscan (ibeg,len,chars,'/') - 1
c                                 find next blank
      iend = iscan (ibeg,len,chars,' ') - 1
c                                 three cases:
      if (iend.ge.len) then

         ier = 1
         goto 99 

      else if (iback.gt.iend) then
c                                 no fraction
         if (iend-ibeg+1.gt.30) goto 90
c                                 first constant
         write (num,'(30a)') (chars(i),i=ibeg,iend)
         read (num,*,err=90) rnum

      else 
c                                 fraction write numerator
         if (iback+1-ibeg.gt.30) goto 90
c                                 first number
         write (num,'(30a)') (chars(i),i=ibeg,iback)       
         read (num,*,err=90) rnum
c                                 second number 

         if (iend-iback-1.gt.30) goto 90
         write (num,'(30a)') (chars(i),i=iback+2,iend)      
         read (num,*,err=90) rnum1

         rnum = rnum/rnum1

      end if 

      ibeg = iend + 1

      goto 99

90    ier = 2

99    end 

      integer function iscan (ibeg,iend,chars,char)
c----------------------------------------------------------------------
c iscan finds the first occurence of char in chars(ibeg..iend)
c----------------------------------------------------------------------
      implicit none

      character chars(240)*1, char*1

      integer ibeg, iend
c----------------------------------------------------------------------
      do iscan = ibeg, iend

         if (chars(iscan).eq.char) exit

      end do 

      end 

      integer function iscnlt (ibeg,iend,chars,char)
c----------------------------------------------------------------------
c iscan finds the first occurence of a character in chars(ibeg..iend) that
c is greater than char. assuming ascii collating sequence +/- < 0 < a
c----------------------------------------------------------------------
      implicit none

      character chars(240)*1, char*1

      integer ibeg, iend


      do iscnlt = ibeg, iend

         if (chars(iscnlt).gt.char) exit

      end do 

      end 

      subroutine getphi (name,eof)
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'
 
      integer ibase, i, it, j, jlam

      double precision ct

      logical eof

      character*8 name, oname, record*1
 
      double precision emodu
      common/ cst318 /emodu(k15)

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

      integer iamir
      common/ cst53 /iamir
 
      double precision delta
      common/ cst325 /delta(11)

      save oname
 
      data oname/' '/
c----------------------------------------------------------------------
      eof = .false.
 
30    read (n2,1000,end=90,err=98) record

      if (record.eq.' ') then
c                          check for comments, i.e., data 
c                          with a blank 1st character 
         goto 30
      else  
         backspace (n2)
         read (n2,1000,end=90,err=98)
     *                     name, ibase, ikind, ilam, idiso
      end if 

      read (n2,*,err=98) (comp(i), i= 1, icmpn), 
     *                   (therm(i), i= 1, k14)
c                               do component transformation if
c                               itrans is not zero
      if (itrans.gt.0) then
 
         do i = 1, itrans
            it = ictr(i)
            if (comp(it).ne.0d0.and.ctrans(it,i).ne.0d0) then
c                                ct is how much of the new
c                                component is in the phase.
               ct =  comp(it) / ctrans(it,i)
 
               do j = 1, icmpn
                  comp(j) = comp(j) - ct * ctrans(j,i)
               end do 
 
               comp(it) = ct
            end if 
         end do
      end if
 
      if (ilam.ne.0) then
c                                 determine number of transitions from
c                                 flag ilam:
         jlam = ilam

         if (ilam.gt.9) then 
            jlam = ilam - 9
         else if (ilam.gt.6) then
            jlam = ilam - 6
         else if (ilam.gt.3) then
            jlam = ilam - 3
         end if 
 
         do i= 1, jlam
            read (n2,*,err=98) (tm(j,i), j = 1, m7 - 2)
         end do 

      end if
 
      if (idiso.ne.0) read (n2,*,err=98) td

      if (ikind.ne.0) read (n2,*,err=98) emodu
c                                 read uncertainties for MC calculations
      if (iamir.eq.999) read (n2,*,err=98) delta
  
      oname = name
 
      return
 
90    eof = .true.

      goto 99

98    if (oname.ne.' ') then
         call error (23,ct,i,oname)
       else
         call error (23,ct,i,'  none  ')
       end if

1000  format (a,i2,i2,i2,i2)
 
99    end

      subroutine deblnk (text)
c----------------------------------------------------------------------
c deblnk - scan text and delete multiple blank characters, strip
c out sequential + - or - + operators, trailing -/+ operators, 
c leading blanks, and leading + operators.
 
c     text - character string 

c no test is made for a blank string or a string of "+" signs.
c----------------------------------------------------------------------
      implicit none 

      integer i, ict, nchar

      logical strip
 
      character text*(*), bitsy(400)*1 
 
      nchar = len(text) 

      read (text,1000) (bitsy(i), i = 1, nchar)
c                                find last non-blank
      ict = 1 
      
      do i = 1, nchar
         if (bitsy(i).gt.' ') ict = i
      end do  

      nchar = ict          
c                                 kill any trailing +/-
      if (bitsy(nchar).eq.'+'.or.bitsy(nchar).eq.'-') nchar = nchar - 1
         
c                                 scan for first non blank/+ character:
      ict = 0 
      
      do i = 1, nchar
         if (bitsy(i).eq.' '.or.bitsy(i).eq.'+') cycle
         ict = i
         exit 
      end do 
c                                 shift everything right
      if (ict.gt.1) then 

         ict = ict - 1
         
         do i = ict+1, nchar
            bitsy(i-ict) = bitsy(i)
         end do 

         nchar = nchar - ict

      end if 

      ict = 1
      
      do i = 2, nchar
c                                 strip out double blanks
         if (bitsy(i).eq.' '.and.bitsy(i-1).eq.' ') cycle 
         ict = ict + 1
         bitsy(ict) = bitsy(i)

      end do

      nchar = ict      
c                                 strip put + - and - + strings
      strip = .false.

      do i = 1, nchar - 2

         if (bitsy(i).eq.'+'.and.bitsy(i+1).eq.'-'.or.
     *       bitsy(i).eq.'-'.and.bitsy(i+1).eq.'+') then

             bitsy(i) = '-'
             bitsy(i+1) = ' '
             strip = .true.

         else if (bitsy(i).eq.'+'.and.bitsy(i+2).eq.'-'.or.
     *            bitsy(i).eq.'-'.and.bitsy(i+2).eq.'+') then

             bitsy(i) = '-'
             bitsy(i+2) = ' '
             strip = .true.

         end if 

      end do 

      if (strip) then 
c                                 strip out new double blanks
         ict = 1

         do i = 2, nchar

            if (bitsy(i).eq.' '.and.bitsy(i-1).eq.' ') cycle 
            ict = ict + 1
            bitsy(ict) = bitsy(i)

         end do 

      end if 

      write (text,1000) (bitsy(i), i = 1, ict),
     *                  (' ',i = ict+1, len(text))
 
1000  format (400a1)

      end

      subroutine getnam (name,ids)
c----------------------------------------------------------------------
c subroutine to retrieve phase name corresponding to index ids
c----------------------------------------------------------------------
      implicit none 

      include 'perplex_parameters.h'

      integer ids

      character names*8, fname*10, name*14
      common/ cst8  /names(k1)/ csta7 /fname(h9)

      if (ids.lt.0) then
c                                 simple compound:
         name = names(-ids)

      else  
c                                 solution phases:
         name = fname(ids)

      end if 
      
      end 

