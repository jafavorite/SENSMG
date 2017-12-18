      subroutine stoponerror
c stop program due to an input or runtime error.
      use COMS
      implicit none
      open(iut,file='stoponerror',status='unknown')
      close (iut)
      close (iuo) ! sensmg.log
      stop
      end subroutine
