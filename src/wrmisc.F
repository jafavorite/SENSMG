      subroutine wrmisc
c write misc input files
      use COMS
      implicit none
      integer i,j,k2,ng
      character mfile*20,ofile*20,sfile*20
c     character ce*100
c
      call set_ebins(neg,iug2,iuo,sens_data,libname,ebins)
      k2=0
      do i=1,nm
c write misc input file.
        write(mfile,'("misc/m",i3.3,"_misc.inp")')i
        write(ofile,'("m",i3.3,"_misc.out")')i
        write(sfile,'("m",i3.3,"_misc.src")')i
        open(2,file=mfile,status='unknown')
        write(2,'("output = ",a)')trim(ofile)
        write(2,'("srcout = ",a)')trim(sfile)
        write(2,'("decayfile = endf7.dk.xml")')
        write(2,'("abundfile = nist.na.xml")')
        write(2,'("pelib = endf7")')
        write(2,'("format = v")')
        write(2,'("particle = 1")')
        write(2,'("matspec = ")')
        do j=1,ncb(i)
          k2=k2+1
          write(2,'(12x,i7,1pe15.7)')nint(blk(1,k2)),-blk(2,k2)
        end do ! j
        write(2,'("density = ",0pf15.8)')-rho(i)
        write(2,'("ebins = ")')
        write(2,'(7x,1p5e12.5)')ebins(1,neg),(ebins(2,ng),ng=neg,1,-1)
        close(2)
      end do ! i
c
      return
      end subroutine
