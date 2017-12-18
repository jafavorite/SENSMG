      subroutine calcmom(angf,momf)
c compute flux moments momf (in mesh centers) given angular fluxes angf
c (on mesh edges).
c for now this should only be called for spheres.
      use COMS
      implicit none
      real*8 plgndr ! function to return legendre polynomials
      real*8 angf(neg,ndir,nitm,njtm),momf(neg,0:nmom-1,nitm,njtm)
      real*8 wplg(0:isct,isn),angc
      integer i,j,l,m,ng
c
      if(calc2d.eq.1)then
        write(*,'("error. in calcmom to calculate moments from ",
     1   "angular fluxes for a cylinder. not set up for that.")')
        write(iuo,'("error. in calcmom to calculate moments from ",
     1   "angular fluxes for a cylinder. not set up for that.")')
        call stoponerror
      end if
c construct table of wgt*plegendre.
      do l=0,isct
        do m=1,isn
c         plg(l,m)=plgndr(l,0,dir(m))
          wplg(l,m)=wgt(m)*plgndr(l,0,dir(m))
c         write(*,'(2i3,1p20e14.6)')l,m,wgt(m),dir(m),wplg(l,m)
        end do ! m
      end do ! l
c
      momf(1:neg,0:nmom-1,1:nitm,1:njtm)=0.d0
      do i=1,it
        do j=1,1
          do ng=1,neg
            do l=0,isct
              do m=1,isn
                angc=(angf(ng,m,i,j)+angf(ng,m,i+1,j))*0.5d0
                momf(ng,l,i,j)=momf(ng,l,i,j)+wplg(l,m)*angc
              end do ! m
            end do ! l
          end do ! ng
        end do ! j
      end do ! i
c
      return
      end subroutine
