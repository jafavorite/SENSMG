      subroutine integrals_r_2d(is,lam,rlip,phifor,phiadj,psifor,psiadj,
     1 neg0,ndir0,nitm0,njtm0,dtgr,dsgr,dfgr,dqgr,ddgr,dtgz,dsgz,dfgz,
     2 dqgz,ddgz)
c do surface integrals for cylinders. 
c is=0/1 for inner products/forward only.
c phifor and phiadj are input forward and adjoint moments.
c psifor and psiadj are input forward and adjoint angular fluxes.
c moments and angular fluxes are mesh-centered.
c
      use COMS
      implicit none
      integer is,neg0,ndir0,nitm0,njtm0
      real*8 lam,rlip(0:isct),phifor(neg,0:nmom-1,nitm,njtm),
     1 phiadj(neg,0:nmom-1,nitm,njtm),
     2 psifor(neg0,ndir0,nitm0,njtm0),psiadj(neg0,ndir0,nitm0,njtm0),
     3 dtgr(neg,0:nr,nz),dsgr(neg,neg,0:isct,0:nr,nz),
     4 dfgr(neg,neg,0:nr,nz),dqgr(neg,0:nr,nz),ddgr(neg,0:nr,nz),
     5 dtgz(neg,nr,0:nz),dsgz(neg,neg,0:isct,nr,0:nz),
     6 dfgz(neg,neg,nr,0:nz),dqgz(neg,nr,0:nz),ddgz(neg,nr,0:nz)
      real*8 omi(neg,0:max(nr,nz),max(it,jt)),
     1 fors(neg,1:nmom,0:max(nr,nz),max(it,jt)),
     2 adjs(neg,1:nmom,0:max(nr,nz),max(it,jt)),
     3 st,fmregc,fmadjc,sa,rlsa,mcon(nmom,ndir),di(ndir)
      integer ir,jr,ng,lg,lg1,lg2,i,j,k,ll,l,m,m1,m2,n,iadir(ndir)
c
c omi, fors, and adjs are used for radial and axial surfaces.
      omi(1:neg,0:max(nr,nz),1:max(it,jt))=0.d0
      fors(1:neg,1:nmom,0:max(nr,nz),1:max(it,jt))=0.d0
      adjs(1:neg,1:nmom,0:max(nr,nz),1:max(it,jt))=0.d0
      dtgr(1:neg,0:nr,1:nz)=0.d0
      dsgr(1:neg,1:neg,0:isct,0:nr,1:nz)=0.d0
      dfgr(1:neg,1:neg,0:nr,1:nz)=0.d0
      dqgr(1:neg,0:nr,1:nz)=0.d0
      ddgr(1:neg,0:nr,1:nz)=0.d0
      dtgz(1:neg,1:nr,0:nz)=0.d0
      dsgz(1:neg,1:neg,0:isct,1:nr,0:nz)=0.d0
      dfgz(1:neg,1:neg,1:nr,0:nz)=0.d0
      dqgz(1:neg,1:nr,0:nz)=0.d0
      ddgz(1:neg,1:nr,0:nz)=0.d0
c
c this section is only needed if flux moments are not used and
c angular fluxes are used.
      if(iangflux.eq.1)then
        call adjdirmp(ndir,iadir)
        call momcon(isct,nmom,mcon,di)
      end if ! iangflux.eq.1
c
c radial surfaces.
c get moments on radial mesh surfaces.
      do ir=1,nr
        i=iindex(ir)+1
        if(iangflux.eq.1)then
          do j=1,jt ! axial segments on radii
            do ng=1,neg
              do m=1,ndir
                n=iadir(m)
                if(ir.lt.nr)then
                  fmregc=(psifor(ng,m,i-1,j)+psifor(ng,m,i,j))*0.5d0
                  fmadjc=(psiadj(ng,n,i-1,j)+psiadj(ng,n,i,j))*0.5d0
                else
                  fmregc=psifor(ng,m,i-1,j)
                  fmadjc=psiadj(ng,n,i-1,j)
                end if
                fors(ng,1:nmom,ir,j)=fors(ng,1:nmom,ir,j)
     1           +fmregc*mcon(1:nmom,m)
                adjs(ng,1:nmom,ir,j)=adjs(ng,1:nmom,ir,j)
     1           +fmadjc*mcon(1:nmom,n)
                if((icalc.eq.1.or.icalc.eq.2).and.
     1            ir.eq.nr)then
c <psi*, psi>=0 on boundary for ev
                  omi(ng,ir,j)=0.d0
                else
                  omi(ng,ir,j)=omi(ng,ir,j)+wgt(m)*fmregc*fmadjc
                end if
              end do ! m
            end do ! ng
c debug moments
c           write(37,'(i6,1p20e14.6)')j,fors(1,1:nmom,nr,j)
c           write(38,'(i6,1p20e14.6)')j,adjs(1,1:nmom,nr,j)
          end do ! j
        else if(iangflux.eq.0)then
          if(ir.lt.nr)then
            fors(1:neg,1:nmom,ir,1:jt)=
     1       (phifor(1:neg,0:nmom-1,i-1,1:jt)
     2       +phifor(1:neg,0:nmom-1,i,1:jt))*0.5d0
            adjs(1:neg,1:nmom,ir,1:jt)=
     1       (phiadj(1:neg,0:nmom-1,i-1,1:jt)
     2       +phiadj(1:neg,0:nmom-1,i,1:jt))*0.5d0
          else
            fors(1:neg,1:nmom,ir,1:jt)=
     1       phifor(1:neg,0:nmom-1,i-1,1:jt)
            adjs(1:neg,1:nmom,ir,1:jt)=
     1       phiadj(1:neg,0:nmom-1,i-1,1:jt)
          end if
        end if
      end do ! ir
c
c inner products on radial surfaces.
      do ir=1,nr
        do jr=1,nz
          m1=mat(ir,jr)
          m2=mat(ir+1,jr)
          sa=2.d0*pi*r(ir)*(z(jr)-z(jr-1))/real(jints(jr),8)
          do j=jindex(jr-1)+1,jindex(jr)
            do ng=1,neg
              if(is.eq.0)then
                ddgr(ng,ir,jr)=ddgr(ng,ir,jr)+fors(ng,1,ir,j)*sa
                dqgr(ng,ir,jr)=dqgr(ng,ir,jr)+adjs(ng,1,ir,j)*sa
c estimate with moments.
                if(iangflux.eq.0)then
                  ll=0
                  do l=0,isct
                    rlsa=rlip(l)*sa
                    do k=0,l
                      ll=ll+1
                      dtgr(ng,ir,jr)=dtgr(ng,ir,jr)+rlsa
     1                 *fors(ng,ll,ir,j)*adjs(ng,ll,ir,j)
                    end do ! k
                  end do ! l
c use angular fluxes.
                else if(iangflux.eq.1)then
                  dtgr(ng,ir,jr)=dtgr(ng,ir,jr)+omi(ng,ir,j)*sa
                end if
c do this in the scattering loop to avoid looping over groups?
                if(ifcel(ir,jr).eq.1.or.ifcel(ir+1,jr).eq.1)then
c zero'th moment has index 1 here.
                  st=lam*adjs(ng,1,ir,j)*sa
                  do lg=1,neg
                    dfgr(ng,lg,ir,jr)=dfgr(ng,lg,ir,jr)+st
     1               *fors(lg,1,ir,j)
                  end do ! lg
                end if
                ll=0
                do l=0,isct
                  rlsa=rlip(l)*sa
                  do k=0,l
                    ll=ll+1
c this logic is here to avoid looping over groups. is it worth it?
c                   if(ifcel(ir,jr).eq.1)then
c                     lg1=1
c                     lg2=neg
c                   else
                      if(ilnk3dnt.eq.0)then
                        lg1=min(scgr(ng,l,m1,1),scgr(ng,l,m2,1))
                        lg2=max(scgr(ng,l,m1,2),scgr(ng,l,m2,2))
                      else if(ilnk3dnt.gt.0)then
                        lg1=1
                        lg2=neg
                      end if
c                   end if
                    st=adjs(ng,ll,ir,j)*rlsa
                    do lg=lg1,lg2
                      dsgr(ng,lg,l,ir,jr)=dsgr(ng,lg,l,ir,jr)+st
     1                 *fors(lg,ll,ir,j)
c                     if(ll.eq.1.and.ifcel(ir,jr).eq.1)then
c                       dfgr(ng,lg,ir,jr)=lam*dsgr(ng,lg,0,ir,jr)
c                     end if
                    end do ! lg
                  end do ! k
                end do ! l
              else if(is.eq.1)then
c TODO forward only (if needed)
              end if
c is.eq.0 is above, inside the lg loop.
              if(is.eq.1.and.(ifcel(ir,jr).eq.1.or.ifcel(ir+1,jr).eq.1))
     1         then
              end if
            end do ! ng
          end do ! j
        end do ! jr
      end do ! ir
c
c omi, fors, and adjs are used for radial and axial surfaces.
c note that the order of the last two indices is not (r,z); it is
c (this direction,the other one).
      omi(1:neg,0:max(nr,nz),1:max(it,jt))=0.d0
      fors(1:neg,1:nmom,0:max(nr,nz),1:max(it,jt))=0.d0
      adjs(1:neg,1:nmom,0:max(nr,nz),1:max(it,jt))=0.d0
c
c axial surfaces.
c get moments on axial mesh surfaces.
      do jr=0,nz
        j=jindex(jr)+1
        if(iangflux.eq.1)then
          do i=1,it ! radial segments on planes
            do ng=1,neg
              do m=1,ndir
                n=iadir(m)
                if(jr.gt.0.and.jr.lt.nz)then
                  fmregc=(psifor(ng,m,i,j-1)+psifor(ng,m,i,j))*0.5d0
                  fmadjc=(psiadj(ng,n,i,j-1)+psiadj(ng,n,i,j))*0.5d0
                else if(jr.eq.0)then
                  fmregc=psifor(ng,m,i,1)
                  fmadjc=psiadj(ng,n,i,1)
                else if(jr.eq.nz)then
                  fmregc=psifor(ng,m,i,j-1)
                  fmadjc=psiadj(ng,n,i,j-1)
                end if
                fors(ng,1:nmom,jr,i)=fors(ng,1:nmom,jr,i)
     1           +fmregc*mcon(1:nmom,m)
                adjs(ng,1:nmom,jr,i)=adjs(ng,1:nmom,jr,i)
     1           +fmadjc*mcon(1:nmom,n)
                if((icalc.eq.1.or.icalc.eq.2).and.
     1            (jr.eq.0.or.jr.eq.nz))then
c <psi*, psi>=0 on boundary for ev
                  omi(ng,jr,i)=0.d0
                else
                  omi(ng,jr,i)=omi(ng,jr,i)+wgt(m)*fmregc*fmadjc
                end if
              end do ! m
            end do ! ng
          end do ! i
        else if(iangflux.eq.0)then
          if(jr.gt.0.and.jr.lt.nz)then
            fors(1:neg,1:nmom,jr,1:it)=
     1       (phifor(1:neg,0:nmom-1,1:it,j-1)
     2       +phifor(1:neg,0:nmom-1,1:it,j))*0.5d0
            adjs(1:neg,1:nmom,jr,1:it)=
     1       (phiadj(1:neg,0:nmom-1,1:it,j-1)
     2       +phiadj(1:neg,0:nmom-1,1:it,j))*0.5d0
          else if(jr.eq.0)then
            fors(1:neg,1:nmom,jr,1:it)=
     1       phifor(1:neg,0:nmom-1,1:it,1)
            adjs(1:neg,1:nmom,jr,1:it)=
     1       phiadj(1:neg,0:nmom-1,1:it,1)
          else if(jr.eq.nz)then
            fors(1:neg,1:nmom,jr,1:it)=
     1       phifor(1:neg,0:nmom-1,1:it,j-1)
            adjs(1:neg,1:nmom,jr,1:it)=
     1       phiadj(1:neg,0:nmom-1,1:it,j-1)
          end if
        end if
      end do ! jr
c
c inner products on axial surfaces.
      do jr=0,nz
        do ir=1,nr
          m1=mat(ir,jr)
          m2=mat(ir,jr+1)
          do i=iindex(ir-1)+1,iindex(ir)
            sa=pi*(rfm(i)**2-rfm(i-1)**2)
            do ng=1,neg
              ddgz(ng,ir,jr)=ddgz(ng,ir,jr)+fors(ng,1,jr,i)*sa
              dqgz(ng,ir,jr)=dqgz(ng,ir,jr)+adjs(ng,1,jr,i)*sa
              if(is.eq.0)then
c estimate with moments.
                if(iangflux.eq.0)then
                  ll=0
                  do l=0,isct
                    rlsa=rlip(l)*sa
                    do k=0,l
                      ll=ll+1
                      dtgz(ng,ir,jr)=dtgz(ng,ir,jr)+rlsa
     1                 *fors(ng,ll,jr,i)*adjs(ng,ll,jr,i)
                    end do ! k
                  end do ! l
c use angular fluxes.
                else if(iangflux.eq.1)then
                  dtgz(ng,ir,jr)=dtgz(ng,ir,jr)+omi(ng,jr,i)*sa
                end if
c do this in the scattering loop to avoid looping over groups?
                if(ifcel(ir,jr).eq.1.or.ifcel(ir,jr+1).eq.1)then
c zero'th moment has index 1 here.
                  st=lam*adjs(ng,1,jr,i)*sa
                  do lg=1,neg
                    dfgz(ng,lg,ir,jr)=dfgz(ng,lg,ir,jr)+st
     1               *fors(lg,1,jr,i)
                  end do ! lg
                end if
                ll=0
                do l=0,isct
                  rlsa=rlip(l)*sa
                  do k=0,l
                    ll=ll+1
c this logic is here to avoid looping over groups. is it worth it?
c                   if(ifcel(ir,jr).eq.1)then
c                     lg1=1
c                     lg2=neg
c                   else
                      if(ilnk3dnt.eq.0)then
                        lg1=min(scgr(ng,l,m1,1),scgr(ng,l,m2,1))
                        lg2=max(scgr(ng,l,m1,2),scgr(ng,l,m2,2))
                      else if(ilnk3dnt.eq.1)then
                        lg1=1
                        lg2=neg
                      end if
c                   end if
                    st=adjs(ng,ll,jr,i)*rlsa
                    do lg=lg1,lg2
                      dsgz(ng,lg,l,ir,jr)=dsgz(ng,lg,l,ir,jr)+st
     1                 *fors(lg,ll,jr,i)
c                     if(ll.eq.1.and.ifcel(ir,jr).eq.1)then
c                       dfgz(ng,lg,ir,jr)=lam*dsgz(ng,lg,0,ir,jr)
c                     end if
                    end do ! lg
                  end do ! k
                end do ! l
              else if(is.eq.1)then
c TODO forward only (if needed)
              end if
c is.eq.0 is above, inside the lg loop.
              if(is.eq.1.and.(ifcel(ir,jr).eq.1.or.ifcel(ir,jr+1).eq.1))
     1         then
              end if
            end do ! ng
          end do ! j
        end do ! ir
      end do ! jr
c
      return
      end subroutine
