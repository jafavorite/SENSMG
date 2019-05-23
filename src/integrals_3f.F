      subroutine integrals_3f(lam,phifor,phiadj,ipfis,dchig)
c computes the inner product of the scalar forward and adjoint fluxes
c with the induced fission importance, which is required for the 
c sensitivity of R2 to chi. 
c phifor and phiadj are input forward and adjoint moments.
c
      use COMS
      implicit none
      real*8 lam,phifor(neg,0:nmom-1,nitm,njtm),
     1 phiadj(neg,0:nmom-1,nitm,njtm),ipfis(it),
     2 dchig(neg,neg,nr,nz)
      real*8 st
      integer ir,jr,ng,lg,i,j
c
      dchig(1:neg,1:neg,1:nr,1:nz)=0.d0
c
      do ir=1,nr
        do jr=1,nz
          if(ifcel(ir,jr).eq.0)cycle
          do i=iindex(ir-1)+1,iindex(ir)
            do j=jindex(jr-1)+1,jindex(jr)
              do ng=1,neg
                st=lam*phiadj(ng,0,i,j)*ipfis(i)*dv(i,j)
                do lg=1,neg
                  dchig(ng,lg,ir,jr)=dchig(ng,lg,ir,jr)
     1             +st*phifor(lg,0,i,j)
                end do ! lg
              end do ! ng
            end do ! j
          end do ! i
        end do ! jr
      end do ! ir
c
      return
      end subroutine
