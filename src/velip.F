      subroutine velip(fiss,rlip,adjmom,psiadj,neg0,ndir0,nitm0,njtm0)
c construct <1/vel adj, psi0> for input adj (moments or angular fluxes).
      use COMS
      implicit none
      integer neg0,ndir0,nitm0,njtm0
      real*8 fiss,rlip(0:isct),adjmom(neg,0:nmom-1,nitm,njtm),
     1 psiadj(neg0,ndir0,nitm0,njtm0)
      real*8 rldv,afregc,afadjc
      integer ir,jr,i,j,ng,k,ll,l,m,n
      real*8 fiss_omp(num_threads)
      integer mynum,omp_get_thread_num
c
      fiss=0.d0
      fiss_omp(1:num_threads)=0.d0
      do ir=1,nr
        do jr=1,nz
c         if(ifcel(ir,jr).eq.0)then
c           cycle
c         end if
!$OMP PARALLEL
!$OMP& PRIVATE(i,j,ng,ll,l,rldv,k,m,n,afregc,afadjc,mynum)
!$OMP DO SCHEDULE(DYNAMIC)
          do i=iindex(ir-1)+1,iindex(ir)
            mynum=omp_get_thread_num()+1
            do j=jindex(jr-1)+1,jindex(jr)
              do ng=1,lng
c estimate with moments.
                if(iangflux.eq.0)then
                  ll=-1
                  do l=0,isct
                    rldv=rlip(l)*dv(i,j)/vel(ng)
                    do k=0,l
                      ll=ll+1
                      fiss_omp(mynum)=fiss_omp(mynum)+rldv
     1                 *fmom(ng,ll,i,j)*adjmom(ng,ll,i,j)
                      if(calc2d.eq.0.or.calc2d.eq.-1)exit
                    end do ! k
                  end do ! l
c use angular fluxes.
                else if(iangflux.eq.1)then
                  rldv=dv(i,j)/vel(ng)
                  do m=1,isn
                    afregc=(afreg(ng,m,i,j)+afreg(ng,m,i+1,j))*0.5d0
                    n=isn+1-m ! adjoint direction
                    afadjc=(psiadj(ng,n,i,j)+psiadj(ng,n,i+1,j))*0.5d0
                    fiss=fiss+wgt(m)*afregc*afadjc*rldv
                  end do ! m
                end if
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
