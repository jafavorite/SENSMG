      subroutine integrals_r_2d_as(lam,rlip,mp1,mp2,
     1 dtgr,dsgr,dfgr,dqgr,dtgz,dsgz,dfgz,dqgz)
c do surface integrals for cylinders. 
c for partisn 8_27 and later using asfluxx and asfluxy.
c angular fluxes are mesh-edge.
c
      use COMS
      implicit none
      real*8 lam,rlip(0:isct),
     1 dtgr(neg,0:nr,nz),dsgr(neg,neg,0:isct,0:nr,nz),
     2 dfgr(neg,neg,0:nr,nz),dqgr(neg,0:nr,nz),
     3 dtgz(neg,nr,0:nz),dsgz(neg,neg,0:isct,nr,0:nz),
     4 dfgz(neg,neg,nr,0:nz),dqgz(neg,nr,0:nz)
      integer mp1,mp2
c     real*8 omia(neg,0:max(nr,nz),max(it,jt)),
c    1 forsa(neg,1:nmom,0:max(nr,nz),max(it,jt)),
c    2 adjsa(neg,1:nmom,0:max(nr,nz),max(it,jt))
c     real*8 omia(neg,0:mp1,mp2),
c    1 forsa(neg,1:nmom,0:mp1,mp2),
c    2 adjsa(neg,1:nmom,0:mp1,mp2)
       real*8
     3 st,sa,rlsa,lsrf(2),lszf(2),lsra(2),lsza(2),
     4 mcon(nmom,ndir),di(ndir)
      integer ir,jr,ng,lg,lg1,lg2,i,j,k,ll,l,m,m1,m2,n,iadir(ndir)
c     real*8 afregas(neg,ndir,0:max(nr,nz),max(it,jt)),
c    1 afadjas(neg,ndir,0:max(nr,nz),max(it,jt))
c     real*8 afregas(neg,ndir,0:mp1,mp2),
c    1 afadjas(neg,ndir,0:mp1,mp2)
c
c omia, forsa, and adjsa are used for radial and axial surfaces.
c DEBUG
c     write(*,'("integrals_r_2d_as",8i6)')nr,nz,mp1,it,jt,nitm,njtm,mp2
c TEST
      call allocate_arrays3
c TEST
      omia(1:neg,0:mp1,1:mp2)=0.d0
      forsa(1:neg,1:nmom,0:mp1,1:mp2)=0.d0
      adjsa(1:neg,1:nmom,0:mp1,1:mp2)=0.d0
      dtgr(1:neg,0:nr,1:nz)=0.d0
      dsgr(1:neg,1:neg,0:isct,0:nr,1:nz)=0.d0
      dfgr(1:neg,1:neg,0:nr,1:nz)=0.d0
      dqgr(1:neg,0:nr,1:nz)=0.d0
      dtgz(1:neg,1:nr,0:nz)=0.d0
      dsgz(1:neg,1:neg,0:isct,1:nr,0:nz)=0.d0
      dfgz(1:neg,1:neg,1:nr,0:nz)=0.d0
      dqgz(1:neg,1:nr,0:nz)=0.d0
c
c get quadrature
      call rddantas(4,0,1,afregas)
c
c these are angular fluxes in this routine.
c get adjoint directions.
      call adjdirmp(ndir,iadir)
      call momcon(isct,nmom,mcon,di)
c debug moments
c     write(*,'("integrals_r_2d_as.")')
c     do m=1,ndir
c       write(*,'(i5,1p3e14.6)')m,di(m),eta(m),wgt(m)
c       write(*,'(i5,1p20e14.6)')m,(mcon(l,m),l=1,nmom)
c     end do ! m
c
c first flag is forward/adjoint; third is radial planes.
      call rddantas(1,0,1,afregas)
      call rddantas(2,0,1,afadjas)
c
c radial surfaces.
c get moments on radial mesh surfaces.
      do ir=1,nr
cc !$OMP PARALLEL
cc !$OMP& PRIVATE(j,ng,m,n)
cc !$OMP DO SCHEDULE(DYNAMIC)
      do j=1,jt ! axial segments on radii
        do ng=1,neg
          do m=1,ndir
            n=iadir(m)
cc !$OMP CRITICAL(FORSI)
            forsa(ng,1:nmom,ir,j)=forsa(ng,1:nmom,ir,j)
     1       +afregas(ng,m,ir,j)*mcon(1:nmom,m)
            adjsa(ng,1:nmom,ir,j)=adjsa(ng,1:nmom,ir,j)
     1       +afadjas(ng,n,ir,j)*mcon(1:nmom,n)
            omia(ng,ir,j)=omia(ng,ir,j)+wgt(m)*afregas(ng,m,ir,j)
     1       *afadjas(ng,n,ir,j)
cc !$OMP END CRITICAL(FORSI)
          end do ! m
        end do ! ng
c debug moments
c       write(47,'(i6,1p20e14.6)')j,forsa(1,1:nmom,ir,j)
c       write(48,'(i6,1p20e14.6)')j,adjsa(1,1:nmom,ir,j)
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
              lsrf(k)=lsrf(k)+wgt(m)*abs(di(m))*afregas(ng,m,ir,j)*sa
              lsra(k)=lsra(k)+wgt(m)*abs(di(m))*afadjas(ng,m,ir,j)*sa
            end do ! m
            dqgr(ng,ir,jr)=dqgr(ng,ir,jr)+adjsa(ng,1,ir,j)*sa
            dtgr(ng,ir,jr)=dtgr(ng,ir,jr)+omia(ng,ir,j)*sa
c do this in the scattering loop to avoid looping over groups?
            if(ifcel(ir,jr).eq.1.or.ifcel(ir+1,jr).eq.1)then
c zero'th moment has index 1 here.
              st=lam*adjsa(ng,1,ir,j)*sa
              do lg=1,neg
                dfgr(ng,lg,ir,jr)=dfgr(ng,lg,ir,jr)+st
     1           *forsa(lg,1,ir,j)
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
                st=adjsa(ng,ll,ir,j)*rlsa
                do lg=lg1,lg2
                  dsgr(ng,lg,l,ir,jr)=dsgr(ng,lg,l,ir,jr)+st
     1             *forsa(lg,ll,ir,j)
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
      if((icalc.eq.0.or.icalc.eq.3.or.icalc.eq.4).and.ir.eq.nr)then
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
      end do ! ir
c
c omia, forsa, and adjsa are used for radial and axial surfaces.
c note that the order of the last two indices is not (r,z); it is
c (this direction,the other one).
c
c first flag is forward/adjoint; third is axial planes.
      call rddantas(1,0,2,afregas)
      call rddantas(2,0,2,afadjas)
      omia(1:neg,0:mp1,1:mp2)=0.d0
      forsa(1:neg,1:nmom,0:mp1,1:mp2)=0.d0
      adjsa(1:neg,1:nmom,0:mp1,1:mp2)=0.d0
c
c axial surfaces.
c get moments on axial mesh surfaces.
      do jr=0,nz
cc !$OMP PARALLEL
cc !$OMP& PRIVATE(i,ng,m,n)
cc !$OMP DO SCHEDULE(DYNAMIC)
      do i=1,it ! radial segments on planes
        do ng=1,neg
          do m=1,ndir
            n=iadir(m)
cc !$OMP CRITICAL(FORSJ)
            forsa(ng,1:nmom,jr,i)=forsa(ng,1:nmom,jr,i)
     1       +afregas(ng,m,jr,i)*mcon(1:nmom,m)
            adjsa(ng,1:nmom,jr,i)=adjsa(ng,1:nmom,jr,i)
     1       +afadjas(ng,n,jr,i)*mcon(1:nmom,n)
            omia(ng,jr,i)=omia(ng,jr,i)+wgt(m)*afregas(ng,m,jr,i)
     1       *afadjas(ng,n,jr,i)
cc !$OMP END CRITICAL(FORSJ)
          end do ! m
        end do ! ng
c debug moments
c       write(47,'(i6,1p20e14.6)')i,forsa(1,1:nmom,jr,i)
c       write(48,'(i6,1p20e14.6)')i,adjsa(1,1:nmom,jr,i)
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
              lszf(k)=lszf(k)+wgt(m)*abs(eta(m))*afregas(ng,m,jr,i)*sa
              lsza(k)=lsza(k)+wgt(m)*abs(eta(m))*afadjas(ng,m,jr,i)*sa
            end do ! m
            dqgz(ng,ir,jr)=dqgz(ng,ir,jr)+adjsa(ng,1,jr,i)*sa
            dtgz(ng,ir,jr)=dtgz(ng,ir,jr)+omia(ng,jr,i)*sa
c do this in the scattering loop to avoid looping over groups?
            if(ifcel(ir,jr).eq.1.or.ifcel(ir,jr+1).eq.1)then
c zero'th moment has index 1 here.
              st=lam*adjsa(ng,1,jr,i)*sa
              do lg=1,neg
                dfgz(ng,lg,ir,jr)=dfgz(ng,lg,ir,jr)+st
     1           *forsa(lg,1,jr,i)
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
                st=adjsa(ng,ll,jr,i)*rlsa
                do lg=lg1,lg2
                  dsgz(ng,lg,l,ir,jr)=dsgz(ng,lg,l,ir,jr)+st
     1             *forsa(lg,ll,jr,i)
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
      if((icalc.eq.0.or.icalc.eq.3.or.icalc.eq.4).and.
     1 (jr.eq.0.or.jr.eq.nz))then
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
      end do ! jr
c
      return
      end subroutine
