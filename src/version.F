      subroutine version(iuo)
c write the version number
      implicit none
      integer iuo
      character(8), parameter :: versn = "5.0.7   "
      character(8), parameter :: datev = "20190716"
c
      write(*,'("code version ",a8," dated ",a8)')versn,datev
      write(iuo,'("code version ",a8," dated ",a8)')versn,datev
      return
      end subroutine
