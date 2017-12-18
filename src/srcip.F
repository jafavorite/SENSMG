      subroutine srcip(fiss,adjmom)
c construct <adj, Q> for input adj (moments or angular fluxes).
      use COMS
      implicit none
      real*8 fiss,adjmom(neg,0:nmom-1,nitm,njtm)
      integer ir,jr,i,j,mt,n,ng
      real*8 fiss_omp(num_threads)
      integer mynum,omp_get_thread_num
c
      fiss=0.d0
      fiss_omp(1:num_threads)=0.d0
      do ir=1,nr
        do jr=1,nz
          mt=mat(ir,jr)
          if(ismat(mt).eq.0)cycle
!$OMP PARALLEL
!$OMP& PRIVATE(i,j,ng,mynum)
!$OMP DO SCHEDULE(DYNAMIC)
          do i=iindex(ir-1)+1,iindex(ir)
            mynum=omp_get_thread_num()+1
            do j=jindex(jr-1)+1,jindex(jr)
              do ng=1,lng
                fiss_omp(mynum)=fiss_omp(mynum)+adjmom(ng,0,i,j)
     1           *nsrcf(ng,mt)*dv(i,j)
              end do ! ng
            end do ! j
          end do ! i
!$OMP END DO
!$OMP END PARALLEL
        end do ! jr
      end do ! ir
c
      do n=1,num_threads
        fiss=fiss+fiss_omp(n)
      end do ! n
c
      return
      end subroutine
