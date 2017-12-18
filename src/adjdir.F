      subroutine adjdirmp(ndir,iadir)
c determine adjoint direction corresponding to forward direction
c for a 2-d cylinder when trnsol=1 (multiprocessor run).
      implicit none
      integer ndir,iadir(ndir)
      integer m,n4
c
      n4=ndir/4
c 1st quadrant to 4th
      do m=1,n4
        iadir(m)=(3-0)*n4+m
      end do ! m
c 2nd quadrant to 3rd
      do m=n4+1,2*n4
        iadir(m)=(2-1)*n4+m
      end do ! m
c 3rd quadrant to 2nd
      do m=2*n4+1,3*n4
        iadir(m)=(1-2)*n4+m
      end do ! m
c 4th quadrant to 1st
      do m=3*n4+1,4*n4
        iadir(m)=(0-3)*n4+m
      end do ! m
c
c     do m=1,ndir
c       write(*,'(2i6)')m,iadir(m)
c     end do ! m
c
      return
      end

