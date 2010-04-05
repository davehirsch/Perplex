      implicit none
        
      character*100 fname,tname

      include 'perplex_parameters.h'

      double precision a(5*l5,40), b(5*l5,40), x, y, ymax, ymin,
     *                 xmax, xmin, z
       
        integer itic(40),i,icod,j

        write (*,*) 'in file?'
        read (*,1000) fname
        open (10,file=fname)
        write (*,*) 'out file?'
        read (*,1000) tname
        open (11,file=tname)

        do i = 1, 40
           itic(i) = 0
        end do
        xmax = -99d33
        xmin = 99d33
        ymin = 99d33
        ymax = -99d33
       
10      read (10,*,end=20) icod, x, y, z
        icod = icod + 1

        if (x.lt.xmin) xmin = x 
        if (x.gt.xmax) xmax = x
        if (y.lt.ymin) ymin = y
        if (y.gt.ymax) ymax = y
 
        itic(icod) = itic(icod) + 1

        if (itic(icod).gt.5*l5) then 
           write (*,*) '** error ** too many points,',
     *                 ' increase parameter l5.'
           goto 99
        end if 

        a(itic(icod),icod) = x
        b(itic(icod),icod) = y
        goto 10 

20    write (11,2020) 0.,' ',' ',' ',' ',xmax,xmin,ymax,ymin,
     *                '   X   ','   Y   '
      do i = 1, 30
         if (itic(i).ne.0) then

            write (11,2010) itic(i)*2,1,i,1,1,1,1,1,1,0e0
            write (11,*) (a(j,i),b(j,i),j=1,itic(i))

         end if
      end do 

2020  format ('1',/,'0 0 0',/,'0 0 0 0 0 0 ',/,g9.1,1x,a162,3(/,a162),/,
     *        '2 1 2 0 0',/,'0 0 0 0. 0. 0. 0. 0.',/,
     *        4(g12.6,1x),/,a,/,a)
2010  format (i5,1x,8(i3,1x),/,g12.6)
1000  format (a72)
99    end 
