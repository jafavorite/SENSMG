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
c DEBUG MOMENTS
c     real*8 diff
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
c DEBUG MOMENTS
c         write(*,'(2i3,1p20e14.6)')l,m,wgt(m),dir(m),wplg(l,m),
c    1     plgndr(l,0,dir(m))
        end do ! m
      end do ! l
c
      momf(1:neg,0:nmom-1,1:nitm,1:njtm)=0.d0
      do i=1,it
        do j=1,1
          do ng=1,neg
            do l=0,isct
              do m=1,isn
                if(i.eq.1)then
                  angc=angf(ng,m,i,j)*0.5d0+angf(ng,m,i+1,j)*0.5d0
                else
                  angc=(angf(ng,m,i,j)+angf(ng,m,i+1,j))*0.5d0
                end if
                momf(ng,l,i,j)=momf(ng,l,i,j)+wplg(l,m)*angc
              end do ! m
            end do ! l
          end do ! ng
        end do ! j
      end do ! i
c DEBUG MOMENTS
c     do i=1,it
c       do j=1,1
c         do ng=1,neg
c           do l=0,isct
c             diff=(momf(ng,l,i,j)-fmom(ng,l,i,j))/fmom(ng,l,i,j)
c             if(abs(diff).le.1.d-6)then
c               write(*,'(3i6,1p3e16.8)')i,ng,l,momf(ng,l,i,j),
c    1           fmom(ng,l,i,j),diff
c             else
c               write(*,'(3i6,1p3e16.8," check")')i,ng,l,momf(ng,l,i,j),
c    1           fmom(ng,l,i,j),diff
c             end if
c           end do ! l
c         end do ! ng
c       end do ! j
c     end do ! i
c     i=1
c     j=1
c     do ng=1,neg
c       do m=1,isn
c         write(*,'(2i6,1p2e16.8)')ng,m,angf(ng,m,i,j),angf(ng,m,i+1,j)
c       end do ! m
c     end do ! ng
c     call stoponerror
c
      return
      end subroutine
