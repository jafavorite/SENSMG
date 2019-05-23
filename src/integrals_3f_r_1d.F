      subroutine integrals_3f_r_1d(lam,psifor,psiadj,ipfis_s,dchigr)
c computes the inner product of the scalar forward and adjoint fluxes
c with the induced fission importance, which is required for the 
c sensitivity of R2 to chi. 
c psifor and psiadj are input forward and adjoint moments.
c
      use COMS
      implicit none
      real*8 lam,psifor(neg,ndir,nitm,njtm),
     1 psiadj(neg,ndir,nitm,njtm),ipfis_s(0:nr),
     2 dchigr(neg,neg,0:nr,nz)
      real*8 fors(neg,0:nr,nz),adjs(neg,0:nr,nz),st,sa
      integer ir,jr,ng,lg,i,j,i0,m,n
c
      fors(1:neg,0:nr,1:nz)=0.d0
      adjs(1:neg,0:nr,1:nz)=0.d0
      dchigr(1:neg,1:neg,0:nr,1:nz)=0.d0
c
c get moments on coarse mesh surfaces.
c these are calculational moments (except omi).
      if(calc2d.eq.0)then ! will not be 1 in this routine
        i0=1
      else if(calc2d.eq.-1)then
        i0=0
      end if
      j=1 ! this is a 1d routine
      do ir=i0,nr
        i=iindex(ir)+1
        do jr=1,nz ! this will be only 1
          do ng=1,neg
            if(iangflux.eq.1)then
              do m=1,isn
                n=isn+1-m
                fors(ng,ir,jr)=fors(ng,ir,jr)+wgt(m)*psifor(ng,m,i,j)
                adjs(ng,ir,jr)=adjs(ng,ir,jr)+wgt(n)*psiadj(ng,m,i,j)
              end do ! m
            else if(iangflux.eq.0)then
c output moments are mesh centered. TODO check this
c             if(ir.lt.nr)then
c               fors(ng,ir,jr)=(phifor(ng,0,i-1,j)+phifor(ng,0,i,j))*0.5d0
c               adjs(ng,ir,jr)=(phiadj(ng,0,i-1,j)+phiadj(ng,0,i,j))*0.5d0
c             else
c               fors(ng,ir,jr)=phifor(ng,0,i-1,j)
c               adjs(ng,ir,jr)=phiadj(ng,0,i-1,j)
c             end if
            end if
          end do ! ng
        end do ! jr
      end do ! ir

      do ir=i0,nr
        if(calc2d.eq.0)then
          sa=4.d0*pi*r(ir)**2
        else if(calc2d.eq.-1)then
          sa=1.d0
        end if
        do jr=1,nz
          do ng=1,neg
            st=lam*adjs(ng,ir,jr)*sa
            do lg=1,neg
              dchigr(ng,lg,ir,jr)=dchigr(ng,lg,ir,jr)
     1         +st*fors(lg,ir,jr)*ipfis_s(ir)
            end do ! lg
          end do ! ng
        end do ! jr
      end do ! ir
c
      return
      end subroutine
