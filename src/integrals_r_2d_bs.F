      subroutine integrals_r_2d_bs(lam,rlip,
     1 dtgr,dsgr,dfgr,dqgr,dtgz,dsgz,dfgz,dqgz)
c do surface integrals for cylinders. 
c for old partisn versions using bsleft/bsrite and bsbot/bstop.
c moments and angular fluxes are mesh-centered.
c
      use COMS
      implicit none
      real*8 lam,rlip(0:isct),
     1 dtgr(neg,0:nr,nz),dsgr(neg,neg,0:isct,0:nr,nz),
     2 dfgr(neg,neg,0:nr,nz),dqgr(neg,0:nr,nz),
     3 dtgz(neg,nr,0:nz),dsgz(neg,neg,0:isct,nr,0:nz),
     4 dfgz(neg,neg,nr,0:nz),dqgz(neg,nr,0:nz)
      real*8 omi(neg,0:max(nr,nz),max(it,jt)),
     1 fors(neg,1:nmom,0:max(nr,nz),max(it,jt)),
     2 adjs(neg,1:nmom,0:max(nr,nz),max(it,jt)),
     3 st,sa,rlsa,lsrf(2),lszf(2),lsra(2),lsza(2),
     4 mcon(nmom,ndir),di(ndir)
      integer ir,jr,ng,lg,lg1,lg2,i,j,k,ll,l,m,m1,m2,n,iadir(ndir)
      real*8 afregbs(neg,ndir,0:max(nr,nz),max(it,jt)),
     1 afadjbs(neg,ndir,0:max(nr,nz),max(it,jt))
c
c omi, fors, and adjs are used for radial and axial surfaces.
      omi(1:neg,0:max(nr,nz),1:max(it,jt))=0.d0
      fors(1:neg,1:nmom,0:max(nr,nz),1:max(it,jt))=0.d0
      adjs(1:neg,1:nmom,0:max(nr,nz),1:max(it,jt))=0.d0
c don't zero all values because most are calculated in integrals_r_2d first.
      if(iplane.gt.-1)then ! 0 is caught in rdctrl
        dtgr(1:neg,iplane,1:nz)=0.d0
        dsgr(1:neg,1:neg,0:isct,iplane,1:nz)=0.d0
        dfgr(1:neg,1:neg,iplane,1:nz)=0.d0
        dqgr(1:neg,iplane,1:nz)=0.d0
      end if
      if(jplane.gt.-1)then
        dtgz(1:neg,1:nr,jplane)=0.d0
        dsgz(1:neg,1:neg,0:isct,1:nr,jplane)=0.d0
        dfgz(1:neg,1:neg,1:nr,jplane)=0.d0
        dqgz(1:neg,1:nr,jplane)=0.d0
      end if
c
c get quadrature
      if(iplane.gt.-1)then
        call rddantbs(4,0,1,afregbs)
      else if(jplane.gt.-1)then
        call rddantbs(4,0,2,afregbs)
      end if
c
c these are angular fluxes in this routine.
c get adjoint directions.
      call adjdirmp(ndir,iadir)
      call momcon(isct,nmom,mcon,di)
c debug moments
c     write(*,'("integrals_r_2d_bs.")')
c     do m=1,ndir
c       write(*,'(i5,1p3e14.6)')m,di(m),eta(m),wgt(m)
c       write(*,'(i5,1p20e14.6)')m,(mcon(l,m),l=1,nmom)
c     end do ! m
c
      if(iplane.eq.-1)go to 500
c first flag is forward/adjoint; third is radial planes.
      call rddantbs(1,0,1,afregbs)
      call rddantbs(2,0,1,afadjbs)
c
c radial surfaces.
c get moments on radial mesh surfaces.
      ir=iplane
cc !$OMP PARALLEL
cc !$OMP& PRIVATE(j,ng,m,n)
cc !$OMP DO SCHEDULE(DYNAMIC)
      do j=1,jt ! axial segments on radii
        do ng=1,neg
          do m=1,ndir
            n=iadir(m)
cc !$OMP CRITICAL(FORSI)
            fors(ng,1:nmom,ir,j)=fors(ng,1:nmom,ir,j)
     1       +afregbs(ng,m,ir,j)*mcon(1:nmom,m)
            adjs(ng,1:nmom,ir,j)=adjs(ng,1:nmom,ir,j)
     1       +afadjbs(ng,n,ir,j)*mcon(1:nmom,n)
            omi(ng,ir,j)=omi(ng,ir,j)+wgt(m)*afregbs(ng,m,ir,j)
     1       *afadjbs(ng,n,ir,j)
cc !$OMP END CRITICAL(FORSI)
          end do ! m
        end do ! ng
c debug moments
c       write(47,'(i6,1p20e14.6)')j,fors(1,1:nmom,ir,j)
c       write(48,'(i6,1p20e14.6)')j,adjs(1,1:nmom,ir,j)
      end do ! j
cc !$OMP END DO
cc !$OMP END PARALLEL
c
c inner products on radial surfaces.
      lsrf(1:2)=0.d0 ! forward current on radial surface
      lsra(1:2)=0.d0 ! adjoint current on radial surface
      do jr=1,nz
        m1=mat(ir,jr)
        m2=mat(ir+1,jr)
        sa=2.d0*pi*r(ir)*(z(jr)-z(jr-1))/real(jints(jr))
cc !$OMP PARALLEL
cc !$OMP& PRIVATE(j,ng,m,k,st,lg,ll,l,rlsa,lg1,lg2)
cc !$OMP DO SCHEDULE(DYNAMIC)
        do j=jindex(jr-1)+1,jindex(jr)
          do ng=1,neg
            do m=1,ndir
              k=1
              if(di(m).lt.0.d0)k=2
              lsrf(k)=lsrf(k)+wgt(m)*abs(di(m))*afregbs(ng,m,ir,j)*sa
              lsra(k)=lsra(k)+wgt(m)*abs(di(m))*afadjbs(ng,m,ir,j)*sa
            end do ! m
            dqgr(ng,ir,jr)=dqgr(ng,ir,jr)+adjs(ng,1,ir,j)*sa
            dtgr(ng,ir,jr)=dtgr(ng,ir,jr)+omi(ng,ir,j)*sa
c do this in the scattering loop to avoid looping over groups?
            if(ifcel(ir,jr).eq.1.or.ifcel(ir+1,jr).eq.1)then
c zero'th moment has index 1 here.
              st=lam*adjs(ng,1,ir,j)*sa
              do lg=1,neg
                dfgr(ng,lg,ir,jr)=dfgr(ng,lg,ir,jr)+st
     1           *fors(lg,1,ir,j)
              end do ! lg
            end if
            ll=0
            do l=0,isct
              rlsa=rlip(l)*sa
              do k=0,l
                ll=ll+1
c this logic is here to avoid looping over groups. is it worth it?
c               if(ifcel(ir,jr).eq.1)then
c                 lg1=1
c                 lg2=neg
c               else
                  lg1=min(scgr(ng,l,m1,1),scgr(ng,l,m2,1))
                  lg2=max(scgr(ng,l,m1,2),scgr(ng,l,m2,2))
c               end if
                st=adjs(ng,ll,ir,j)*rlsa
                do lg=lg1,lg2
                  dsgr(ng,lg,l,ir,jr)=dsgr(ng,lg,l,ir,jr)+st
     1             *fors(lg,ll,ir,j)
c                 if(ll.eq.1.and.ifcel(ir,jr).eq.1)then
c                   dfgr(ng,lg,ir,jr)=lam*dsgr(ng,lg,0,ir,jr)
c                 end if
                end do ! lg
              end do ! k
            end do ! l
          end do ! ng
        end do ! j
cc !$OMP END DO
cc !$OMP END PARALLEL
      end do ! jr
      write(*,'("forward current on radial surface",i4,": j+, j-",
     1 1p2e14.6)')ir,lsrf(1:2)
      write(iuo,'("forward current on radial surface",i4,": j+, j-",
     1 1p2e14.6)')ir,lsrf(1:2)
      if(icalc.eq.0.and.ir.eq.nr)then
        write(*,'("adjoint current on radial surface",i4,": j+, j-",
     1   1p2e14.6," (includes source)")')ir,lsra(1:2)
        write(iuo,'("adjoint current on radial surface",i4,": j+, j-",
     1   1p2e14.6," (includes source)")')ir,lsra(1:2)
      else
        write(*,'("adjoint current on radial surface",i4,": j+, j-",
     1   1p2e14.6)')ir,lsra(1:2)
        write(iuo,'("adjoint current on radial surface",i4,": j+, j-",
     1   1p2e14.6)')ir,lsra(1:2)
      end if
c
c omi, fors, and adjs are used for radial and axial surfaces.
c note that the order of the last two indices is not (r,z); it is
c (this direction,the other one).
  500 if(jplane.eq.-1)go to 900
c first flag is forward/adjoint; third is axial planes.
      call rddantbs(1,0,2,afregbs)
      call rddantbs(2,0,2,afadjbs)
      omi(1:neg,0:max(nr,nz),1:max(it,jt))=0.d0
      fors(1:neg,1:nmom,0:max(nr,nz),1:max(it,jt))=0.d0
      adjs(1:neg,1:nmom,0:max(nr,nz),1:max(it,jt))=0.d0
c
c axial surfaces.
c get moments on axial mesh surfaces.
      jr=jplane
cc !$OMP PARALLEL
cc !$OMP& PRIVATE(i,ng,m,n)
cc !$OMP DO SCHEDULE(DYNAMIC)
      do i=1,it ! radial segments on planes
        do ng=1,neg
          do m=1,ndir
            n=iadir(m)
cc !$OMP CRITICAL(FORSJ)
            fors(ng,1:nmom,jr,i)=fors(ng,1:nmom,jr,i)
     1       +afregbs(ng,m,jr,i)*mcon(1:nmom,m)
            adjs(ng,1:nmom,jr,i)=adjs(ng,1:nmom,jr,i)
     1       +afadjbs(ng,n,jr,i)*mcon(1:nmom,n)
            omi(ng,jr,i)=omi(ng,jr,i)+wgt(m)*afregbs(ng,m,jr,i)
     1       *afadjbs(ng,n,jr,i)
cc !$OMP END CRITICAL(FORSJ)
          end do ! m
        end do ! ng
c debug moments
c       write(47,'(i6,1p20e14.6)')i,fors(1,1:nmom,jr,i)
c       write(48,'(i6,1p20e14.6)')i,adjs(1,1:nmom,jr,i)
      end do ! i
cc !$OMP END DO
cc !$OMP END PARALLEL
c
c inner products on axial surfaces.
      lszf(1:2)=0.d0 ! forward current on axial surface
      lsza(1:2)=0.d0 ! adjoint current on axial surface
      do ir=1,nr
        m1=mat(ir,jr)
        m2=mat(ir,jr+1)
cc !$OMP PARALLEL
cc !$OMP& PRIVATE(i,sa,ng,m,k,st,lg,ll,l,rlsa,lg1,lg2)
cc !$OMP DO SCHEDULE(DYNAMIC)
        do i=iindex(ir-1)+1,iindex(ir)
          sa=pi*(rfm(i)**2-rfm(i-1)**2)
          do ng=1,neg
            do m=1,ndir
              k=1
              if(eta(m).lt.0.d0)k=2
              lszf(k)=lszf(k)+wgt(m)*abs(eta(m))*afregbs(ng,m,jr,i)*sa
              lsza(k)=lsza(k)+wgt(m)*abs(eta(m))*afadjbs(ng,m,jr,i)*sa
            end do ! m
            dqgz(ng,ir,jr)=dqgz(ng,ir,jr)+adjs(ng,1,jr,i)*sa
            dtgz(ng,ir,jr)=dtgz(ng,ir,jr)+omi(ng,jr,i)*sa
c do this in the scattering loop to avoid looping over groups?
            if(ifcel(ir,jr).eq.1.or.ifcel(ir,jr+1).eq.1)then
c zero'th moment has index 1 here.
              st=lam*adjs(ng,1,jr,i)*sa
              do lg=1,neg
                dfgz(ng,lg,ir,jr)=dfgz(ng,lg,ir,jr)+st
     1           *fors(lg,1,jr,i)
              end do ! lg
            end if
            ll=0
            do l=0,isct
              rlsa=rlip(l)*sa
              do k=0,l
                ll=ll+1
c this logic is here to avoid looping over groups. is it worth it?
c               if(ifcel(ir,jr).eq.1)then
c                 lg1=1
c                 lg2=neg
c               else
                  lg1=min(scgr(ng,l,m1,1),scgr(ng,l,m2,1))
                  lg2=max(scgr(ng,l,m1,2),scgr(ng,l,m2,2))
c               end if
                st=adjs(ng,ll,jr,i)*rlsa
                do lg=lg1,lg2
                  dsgz(ng,lg,l,ir,jr)=dsgz(ng,lg,l,ir,jr)+st
     1             *fors(lg,ll,jr,i)
c                 if(ll.eq.1.and.ifcel(ir,jr).eq.1)then
c                   dfgz(ng,lg,ir,jr)=lam*dsgz(ng,lg,0,ir,jr)
c                 end if
                end do ! lg
              end do ! k
            end do ! l
          end do ! ng
        end do ! i
cc !$OMP END DO
cc !$OMP END PARALLEL
      end do ! ir
      write(*,'("forward current on axial surface",i4,": j+, j-",
     1 1p2e14.6)')jr,lszf(1:2)
      write(iuo,'("forward current on axial surface",i4,": j+, j-",
     1 1p2e14.6)')jr,lszf(1:2)
      if(icalc.eq.0.and.(jr.eq.0.or.jr.eq.nz))then
        write(*,'("adjoint current on axial surface",i4,": j+, j-",
     1   1p2e14.6," (includes source)")')jr,lsza(1:2)
        write(iuo,'("adjoint current on axial surface",i4,": j+, j-",
     1   1p2e14.6," (includes source)")')jr,lsza(1:2)
       else
        write(*,'("adjoint current on axial surface",i4,": j+, j-",
     1   1p2e14.6)')jr,lsza(1:2)
        write(iuo,'("adjoint current on axial surface",i4,": j+, j-",
     1   1p2e14.6)')jr,lsza(1:2)
       end if
c
  900 return
      end subroutine
