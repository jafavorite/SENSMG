      subroutine massin
c compute isotopic mass inventory for redoin/lnk3dnt to compare
c with partisn.
      use COMS
      implicit none
      real*8 massiso(niso)
      integer ir,jr,m,mt
c
      massiso(1:niso)=0.d0
      do ir=1,it
        do jr=1,jt
          if(idclnk3(ir,jr,1).eq.0)cycle
          do m=1,nmxsp
            mt=idclnk3(ir,jr,m)
            if(mt.eq.0)exit
            massiso(mt)=massiso(mt)+denlnk3(ir,jr,m)*vol(ir,jr)
c partisn unit is at/b-cm*cm^3=10^24 atoms
c    1       *atwt(mt)/avg
          end do ! m
        end do ! jr
      end do ! ir
      write(iuo,'("mass inventory",/,5x,"iso     mass (10^24 at)")')
      do m=1,niso
        write(iuo,'(i7,1pe17.9)')m,massiso(m)
      end do ! nm
c
      return
      end subroutine
