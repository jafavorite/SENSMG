      subroutine integrals_r_1d(is,lam,rlip,phifor,phiadj,psifor,psiadj,
     1 dtgr,dsgr,dfgr,dqgr,ddgr)
c do surface integrals for spheres and slabs.
c is=0/1 for inner products/forward only.
c phifor and phiadj are input forward and adjoint moments.
c psifor and psiadj are input forward and adjoint angular fluxes.
c
      use COMS
      implicit none
      integer is
      real*8 lam,rlip(0:isct),phifor(neg,0:nmom-1,nitm,njtm),
     1 phiadj(neg,0:nmom-1,nitm,njtm),
     2 psifor(neg,ndir,nitm,njtm),psiadj(neg,ndir,nitm,njtm),
     3 dtgr(neg,0:nr,nz),dsgr(neg,neg,0:isct,0:nr,nz),
     4 dfgr(neg,neg,0:nr,nz),dqgr(neg,0:nr,nz),ddgr(neg,0:nr,nz)
      real*8 omi(neg,0:nr,nz),fors(neg,0:isct,0:nr,nz),
     1 adjs(neg,0:isct,0:nr,nz),st,fmregc,fmadjc,sa,rlsa
      real*8 plgndr ! function to return legendre polynomials
      real*8 wplg(0:isct,isn)
      integer i0,ir,jr,ng,lg,lg1,lg2,i,j,l,m,m1,m2,n
c
      omi(1:neg,0:nr,1:nz)=0.d0
      fors(1:neg,0:isct,0:nr,1:nz)=0.d0
      adjs(1:neg,0:isct,0:nr,1:nz)=0.d0
      dtgr(1:neg,0:nr,1:nz)=0.d0
      dsgr(1:neg,1:neg,0:isct,0:nr,1:nz)=0.d0
      dfgr(1:neg,1:neg,0:nr,1:nz)=0.d0
      dqgr(1:neg,0:nr,1:nz)=0.d0
      ddgr(1:neg,0:nr,1:nz)=0.d0
c
c construct table of wgt*plegendre.
      if(iangflux.eq.1)then
        do l=0,isct
          do m=1,isn
c           plg(l,m)=plgndr(l,0,dir(m))
            wplg(l,m)=wgt(m)*plgndr(l,0,dir(m))
c           write(*,'(2i3,1p20e14.6)')l,m,wgt(m),dir(m),wplg(l,m)
          end do ! m
        end do ! l
      end if
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
                omi(ng,ir,jr)=omi(ng,ir,jr)+wgt(m)
     1           *psifor(ng,m,i,j)*psiadj(ng,n,i,j)
                do l=0,isct
                  fors(ng,l,ir,jr)=fors(ng,l,ir,jr)+wplg(l,m)
     1             *psifor(ng,m,i,j)
                  adjs(ng,l,ir,jr)=adjs(ng,l,ir,jr)+wplg(l,m)
     1             *psiadj(ng,m,i,j)
                end do ! l
              end do ! m
            else if(iangflux.eq.0)then
c output moments are mesh centered. TODO check this
              do l=0,isct
                if(ir.eq.0.and.i.eq.1)then ! for slabs; same for sphere?
                  fmregc=phifor(ng,l,i,j)
                  fmadjc=phiadj(ng,l,i,j)
                else if(ir.lt.nr)then
                  fmregc=(phifor(ng,l,i-1,j)+phifor(ng,l,i,j))*0.5d0 
                  fmadjc=(phiadj(ng,l,i-1,j)+phiadj(ng,l,i,j))*0.5d0 
                else
                  fmregc=phifor(ng,l,i-1,j)
                  fmadjc=phiadj(ng,l,i-1,j)
                end if
                fors(ng,l,ir,jr)=fmregc
                adjs(ng,l,ir,jr)=fmadjc
                omi(ng,ir,jr)=omi(ng,ir,jr)+rlip(l)*fmregc*fmadjc
              end do ! l
            end if
          end do ! ng
        end do ! jr
      end do ! ir
c
c adjoint moments (regular and generalized) here are already inner
c product moments.
      do ir=i0,nr
        if(calc2d.eq.0)then
          sa=4.d0*pi*r(ir)**2
        else if(calc2d.eq.-1)then
          sa=1.d0
        end if
        do jr=1,nz
          m1=mat(ir,jr)
          m2=mat(ir+1,jr)
          do ng=1,neg
            if(is.eq.0)then
              dtgr(ng,ir,jr)=omi(ng,ir,jr)*sa
c do this in the scattering loop to avoid looping over groups?
              if(ifcel(ir,jr).eq.1.or.ifcel(ir+1,jr).eq.1)then
                st=lam*adjs(ng,0,ir,jr)*sa
                do lg=1,neg
                  dfgr(ng,lg,ir,jr)=dfgr(ng,lg,ir,jr)+st
     1             *fors(lg,0,ir,jr)
                end do ! lg
              end if
              do l=0,isct
                rlsa=rlip(l)*sa
c this logic is here to avoid looping over groups. is it worth it?
c               if(ifcel(ir,jr).eq.1)then
c                 lg1=1
c                 lg2=neg
c               else
                  if(ilnk3dnt.eq.0)then
                    lg1=min(scgr(ng,l,m1,1),scgr(ng,l,m2,1))
                    lg2=max(scgr(ng,l,m1,2),scgr(ng,l,m2,2))
                  else if(ilnk3dnt.gt.0)then
                    lg1=1
                    lg2=neg
                  end if
c               end if
                st=adjs(ng,l,ir,jr)*rlsa
                do lg=lg1,lg2
                  dsgr(ng,lg,l,ir,jr)=dsgr(ng,lg,l,ir,jr)+st
     1             *fors(lg,l,ir,jr)
c                 if(l.eq.0.and.ifcel(ir,jr).eq.1)then
c                   dfgr(ng,lg,ir,jr)=lam*dsgr(ng,lg,0,ir,jr)
c                 end if
                end do ! lg
              end do ! l
              ddgr(ng,ir,jr)=ddgr(ng,ir,jr)+fors(ng,0,ir,jr)*sa
              dqgr(ng,ir,jr)=dqgr(ng,ir,jr)+adjs(ng,0,ir,jr)*sa
            else if(is.eq.1)then
c TODO forward only (if needed)
            end if
c is.eq.0 is above, inside the lg loop.
            if(is.eq.1.and.(ifcel(ir,jr).eq.1.or.ifcel(ir+1,jr).eq.1))
     1       then
            end if
          end do ! ng
        end do ! jr
      end do ! ir
c
      end subroutine
