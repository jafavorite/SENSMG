      subroutine output2(sden)
c output mass density sensitivities.
      use COMS
      implicit none
      real*8 sden(neg,nm)
      integer i,ng,n1,n2
c head is the heading " group  material...." length should be
c 6+14*9999 = 139992. 566 is ok for nm.le.40.
      character head*566
c there are several instances of 9999e14.6 and 9999i14. the number
c of materials nm must be .le. 9999.
c
c summary table first
      write(iur,'(" summary (sums over groups)")')
      do i=1,nm
        write(iur,'(" material",i4.4,1pe15.6)')i,sum(sden(1:neg,i))
      end do ! i
c now the details
      head=" "
      write(head,'(1x,"group",2x,"e_lower",7x,"e_upper",7x,
     1 "e_average")')
      n1=6+3*14+3 ! i6,3e14.6, 3 spaces
      do i=1,nm
        n2=n1+14
        write(head(n1:n2),'("material",i4.4)')i
        n1=n2
      end do ! i
      write(iur,'(a)')trim(head)
      do ng=1,neg
        write(iur,'(i6,1p9999e14.6)')ng,ebins(1:3,ng),
     1   (sden(ng,i),i=1,nm)
      end do ! ng
c 42x is 14*3 for the energy groups.
      write(iur,'(" total",42x,1p9999e14.6)')
     1 (sum(sden(1:neg,i)),i=1,nm)
c
      return
      end subroutine
