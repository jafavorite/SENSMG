      subroutine output3(sintr)
c output interface location derivatives for spheres and slabs.
      use COMS
      implicit none
      real*8 sintr(neg,0:nr,nz)
      integer i0,i1,i2,i,j,ng,nblock
c there are 30 radii in each block; can be changed here. 
c nblock+3 <= 999.
      parameter (nblock=30)
c
      j=1 ! this routine for spheres and slabs
      if(calc2d.eq.0)then ! will not be 1 in this routine
        i0=1
      else if(calc2d.eq.-1)then
        i0=0
      end if
c summary table first
      write(iur,'(" summary (sums over groups)")')
      do i=i0,nr
        write(iur,'(" radius",i4.4,1pe17.6)')i,sum(sintr(1:neg,i,j))
      end do ! i
c now the details
      i1=i0
      do while(i1.le.nr)
        i2=min(i1+nblock-1,nr)
        write(iur,'(1x,"group",2x,"e_lower",7x,"e_upper",7x,
     1   "e_average",5x)',advance='no')
        do i=i1,i2
          if(i.eq.i2)then
            write(iur,'("radius",i4.4,4x)',advance='yes')i
          else
            write(iur,'("radius",i4.4,4x)',advance='no')i
          end if
        end do ! i
        do ng=1,neg
          write(iur,'(i6,1p999e14.6)')ng,ebins(1:3,ng),
     1     (sintr(ng,i,j),i=i1,i2)
        end do ! ng
c 42x is 14*3 for the energy groups.
        write(iur,'(" total",42x,1p999e14.6)')
     1   (sum(sintr(1:neg,i,j)),i=i1,i2)
        i1=i2+1
      end do
c
      return
      end subroutine
