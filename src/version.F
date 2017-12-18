      subroutine version(iuo)
c write the version number
      implicit none
      integer iuo
      character(8), parameter :: versn = "4.2.2   "
      character(8), parameter :: datev = "20171214"
c
      write(*,'("code version ",a8," dated ",a8)')versn,datev
      write(iuo,'("code version ",a8," dated ",a8)')versn,datev
      return
      end subroutine
