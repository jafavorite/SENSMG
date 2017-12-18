      subroutine output4(sintr,sintz)
c output interface location derivatives for cylinders.
      use COMS
      implicit none
      real*8 sintr(neg,0:nr,nz),sintz(neg,nr,0:nz)
      integer i,j,ng,n1,n2
c head is the heading " group  radius...." length should be
c 6+14*9999 = 139992. 566 is ok for nr.le.40.
      character head*510
c there are several instances of 9999e14.6 and 9999i14. the number
c of radii nr must be .le. 9999.
c
      write(iur,'(" radial surfaces")')
c summary table first
      write(iur,'(" summary (sums over groups)")')
      head=" "
      n1=4+4+8+1 ! a4,i4.4, 3 spaces before, 5 after
      do i=1,nr
        n2=n1+14
        write(head(n1:n2),'("radius",i4.4)')i
        n1=n2
      end do ! i
      write(iur,'(a)')trim(head)
      do j=nz,1,-1
        write(iur,'("   axis",i4.4,3x,1p9999e14.6)')j,
     1   (sum(sintr(1:neg,i,j)),i=1,nr)
      end do ! j
c now the details
      do j=nz,1,-1
        write(iur,'(" axis",i4.4)')j
        head=" "
        write(head,'(1x,"group",2x,"e_lower",7x,"e_upper",7x,
     1   "e_average")')
        n1=6+3*14+3 ! i6,3e14.6, 3 spaces
        do i=1,nr
          n2=n1+14
          write(head(n1:n2),'("radius",i4.4)')i
          n1=n2
        end do ! i
        write(iur,'(a)')trim(head)
        do ng=1,neg
          write(iur,'(i6,1p9999e14.6)')ng,ebins(1:3,ng),
     1     (sintr(ng,i,j),i=1,nr)
        end do ! ng
c 42x is 14*3 for the energy groups.
        write(iur,'(" total",42x,1p9999e14.6)')
     1   (sum(sintr(1:neg,i,j)),i=1,nr)
      end do ! j
c
      write(iur,'(/," axial surfaces")')
c summary table first
      write(iur,'(" summary (sums over groups)")')
      head=" "
      n1=4+4+8+1 ! a4,i4.4, 3 spaces before, 5 after
      do i=1,nr
        n2=n1+14
        write(head(n1:n2),'("radius",i4.4)')i
        n1=n2
      end do ! i
      write(iur,'(a)')trim(head)
      do j=nz,0,-1
        write(iur,'("   axis",i4.4,3x,1p9999e14.6)')j,
     1   (sum(sintz(1:neg,i,j)),i=1,nr)
      end do ! j
c now the details
      do i=1,nr
        write(iur,'(" radius",i4.4)')i
        head=" "
        write(head,'(1x,"group",2x,"e_lower",7x,"e_upper",7x,
     1   "e_average")')
        n1=6+3*14+3 ! i6,3e14.6, 3 spaces
        do j=nz,0,-1
          n2=n1+14
          write(head(n1:n2),'("axis",i4.4)')j
          n1=n2
        end do ! j
        write(iur,'(a)')trim(head)
        do ng=1,neg
          write(iur,'(i6,1p9999e14.6)')ng,ebins(1:3,ng),
     1     (sintz(ng,i,j),j=nz,0,-1)
        end do ! ng
c 42x is 14*3 for the energy groups.
        write(iur,'(" total",42x,1p9999e14.6)')
     1   (sum(sintz(1:neg,i,j)),j=nz,0,-1)
      end do ! i
c
      return
      end subroutine
