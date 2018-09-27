      subroutine wrmisc
c write misc input files
      use COMS
      implicit none
      integer i,j,k2,ng
      character mfile*20,ofile*20,sfile*20
c     character ce*100
c
      k2=0
      do i=1,nm
c write misc input file.
        write(mfile,'("misc/m",i3.3,"_misc.inp")')i
        write(ofile,'("m",i3.3,"_misc.out")')i
        write(sfile,'("m",i3.3,"_misc.src")')i
        open(iun,file=mfile,status='unknown')
        write(iun,'("output = ",a)')trim(ofile)
        write(iun,'("srcout = ",a)')trim(sfile)
        write(iun,'("decayfile = endf7.dk.xml")')
        write(iun,'("abundfile = nist.na.xml")')
        write(iun,'("pelib = endf7")')
        write(iun,'("format = v")')
        write(iun,'("particle = 1")')
        write(iun,'("matspec = ")')
        do j=1,ncb(i)
          k2=k2+1
          write(iun,'(12x,i7,1pe15.7)')nint(blk(1,k2)),-blk(2,k2)
        end do ! j
        write(iun,'("density = ",0pf15.8)')-rho(i)
        write(iun,'("ebins = ")')
        write(iun,'(7x,1p5e12.5)')ebins(1,neg),(ebins(2,ng),ng=neg,1,-1)
        close(iun)
      end do ! i
c
      return
      end subroutine
