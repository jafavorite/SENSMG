      subroutine fissip(fiss,adjmom)
c construct <F0* adjmom, psi0> for input adjmom.
      use COMS
      implicit none
      real*8 fiss,adjmom(neg,0:nmom-1,nitm,njtm)
      real*8 st,chi1(1:neg),adj1(1:neg) !,kflux
      integer ir,jr,i,j,n,ng,mt !,lg
      real*8 fiss_omp(num_threads) !,kflux_omp(num_threads)
      integer mynum,omp_get_thread_num
c
      fiss=0.d0
      fiss_omp(1:num_threads)=0.d0
c     kflux=0.d0
c     kflux_omp(1:num_threads)=0.d0
      do ir=1,nr
        do jr=1,nz
          if(ifcel(ir,jr).eq.0)cycle
          mt=mat(ir,jr)
!$OMP PARALLEL
!$OMP& PRIVATE(i,j,ng,st,chi1,adj1,mynum)
!$OMP DO SCHEDULE(DYNAMIC)
          do i=iindex(ir-1)+1,iindex(ir)
            mynum=omp_get_thread_num()+1
            do j=jindex(jr-1)+1,jindex(jr)
              if(ilnk3dnt.gt.0)then
                call maclnk3(mt,i,j)
              end if
              do ng=1,lng
                st=nusigf(ng,mt)*fmom(ng,0,i,j)*dv(i,j)
c replace lg loop with dot_product in hopes of speed-up
                chi1(1:neg)=chi(1:neg,ng,mt)
                adj1(1:neg)=adjmom(1:neg,0,i,j)
                fiss_omp(mynum)=fiss_omp(mynum)+dot_product(chi1,adj1)
     1           *st
c               do lg=1,lng
c                 fiss_omp(mynum)=fiss_omp(mynum)+chi(lg,ng,mt)
c    1             *adjmom(lg,0,i,j)*st
c                 kflux_omp(mynum)=kflux_omp(mynum)+chi(lg,ng,mt)*st
c               end do ! lg
              end do ! ng
            end do ! j
c           if(calc2d.eq.1.and.neg.ge.79)write(*,'("fissip, ",
c    1       "finished i=",i5," of",i5)')i,it
          end do ! i
!$OMP END DO
!$OMP END PARALLEL
        end do ! jr
      end do ! ir
c
      do n=1,num_threads
        fiss=fiss+fiss_omp(n)
c       kflux=kflux+kflux_omp(n)
      end do ! n
c
      if(icalc.eq.1)then
        fiss=fiss/keff
      end if
c     write(*,'("<F psi>=",1pe15.7)')kflux
c     write(iuo,'("<F psi>=",1pe15.7)')kflux
c
      return
      end subroutine
