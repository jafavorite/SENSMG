      subroutine output3(sintr)
c output interface location derivatives for spheres and slabs.
      use COMS
      implicit none
      real*8 sintr(neg,0:nr,nz)
      integer i0,i,j,ng,n1,n2
c head is the heading " group  radius...." length should be
c 6+14*9999 = 139992. 566 is ok for nr.le.40.
      character head*510
c there are several instances of 9999e14.6 and 9999i14. the number
c of radii nr must be .le. 9999.
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
      head=" "
      write(head,'(1x,"group",2x,"e_lower",7x,"e_upper",7x,
     1 "e_average")')
      n1=6+3*14+3 ! i6,3e14.6, 3 spaces
      do i=i0,nr
        n2=n1+14
        write(head(n1:n2),'("radius",i4.4)')i
        n1=n2
      end do ! i
      write(iur,'(a)')trim(head)
      do ng=1,neg
        write(iur,'(i6,1p9999e14.6)')ng,ebins(1:3,ng),
     1   (sintr(ng,i,j),i=i0,nr)
      end do ! ng
c 42x is 14*3 for the energy groups.
      write(iur,'(" total",42x,1p9999e14.6)')
     1 (sum(sintr(1:neg,i,j)),i=i0,nr)
c
      return
      end subroutine
