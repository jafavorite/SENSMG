      subroutine temp_integrals_mom(is,isct1,nmom1,lam,psifor,psiadj,
     1 neg0,ndir0,nitm0,njtm0,dtg,dsg,dfg,dqg)
c do the integrals. 
c is=0/1 for inner products/forward only.
c phifor and phiadj are forward and adjoint moments to be computed.
c psifor and psiadj are input forward and adjoint angular fluxes.
c
      use COMS
      implicit none
      integer is,isct1,nmom1,neg0,ndir0,nitm0,njtm0
      real*8 lam,
     1 psifor(neg0,ndir0,nitm0,njtm0),psiadj(neg0,ndir0,nitm0,njtm0),
     2 dtg(neg,nr,nz),dsg(neg,neg,0:isct,nr,nz),
     3 dfg(neg,neg,nr,nz),dqg(neg,nr,nz)
      real*8 st,afregc,afadjc,rldv,
     1 das(neg,nr,nz),dfs(neg,nr,nz),dis(neg,nr,nz),
     2 dos(neg,nr,nz),dss(neg,nr,nz),dastot,dfstot,distot,
     3 dostot,dsstot,phifor(neg,0:nmom1-1,it,jt),
     4 phiadj(neg,0:nmom1-1,it,jt),
     5 mcon(nmom1,ndir),di(ndir),rlip(0:isct1)
      integer ir,jr,ng,lg,lg1,lg2,i,j,k,ll,l,m,mt,n,iadir(ndir)
      real*8 temps
      integer mg,n1,n2,k1,k2,i2lp1,ifissdata
      character line*130
c
      dtg(1:neg,1:nr,1:nz)=0.d0
      dsg(1:neg,1:neg,0:isct,1:nr,1:nz)=0.d0
      dfg(1:neg,1:neg,1:nr,1:nz)=0.d0
      dqg(1:neg,1:nr,1:nz)=0.d0
      das(1:neg,1:nr,1:nz)=0.d0
      dfs(1:neg,1:nr,1:nz)=0.d0
      dis(1:neg,1:nr,1:nz)=0.d0
      dos(1:neg,1:nr,1:nz)=0.d0
      dss(1:neg,1:nr,1:nz)=0.d0
c
c adjoint moments (regular and generalized) are always calculational
c moments. the sign of rlip is used to convert to inner product moments.
      do l=0,isct1
        rlip(l)=real(2*l+1)
        if(mod(l,2).ne.0)then
          rlip(l)=-rlip(l)
        end if
c       write(*,'(0p2f8.3)')rlip(l)
      end do ! l
c
      if(iangflux.eq.1.and.calc2d.eq.1)then
        write(*,'("computing moments in temp_integrals_mom for isct1="
     1   ,i3,"...")')isct1
        call adjdirmp(ndir,iadir)
        call momcon(isct1,nmom1,mcon,di)
c
c compute computational moments from angular fluxes.
        phifor(1:neg,0:nmom1-1,1:it,1:jt)=0.d0
        phiadj(1:neg,0:nmom1-1,1:it,1:jt)=0.d0
        do i=1,it
          do j=1,jt
            do ng=1,neg
              do m=1,ndir
                phifor(ng,0:nmom1-1,i,j)=phifor(ng,0:nmom1-1,i,j)
     1           +psifor(ng,m,i,j)*mcon(1:nmom1,m)
                phiadj(ng,0:nmom1-1,i,j)=phiadj(ng,0:nmom1-1,i,j)
     1           +psiadj(ng,m,i,j)*mcon(1:nmom1,m)
              end do ! m
            end do ! ng
          end do ! j
        end do ! i
        iangflux=0 ! now use the moments
        write(*,'("done.")')
      end if ! iangflux.eq.1
c
c adjoint moments (regular and generalized) are always calculational
c moments. the sign of rlip is used to convert to inner product moments.
      do ir=1,nr
        do jr=1,nz
          mt=mat(ir,jr)
          do i=iindex(ir-1)+1,iindex(ir)
            do j=jindex(jr-1)+1,jindex(jr)
              do ng=1,neg
                if(is.eq.0)then
c estimate with moments.
                  if(iangflux.eq.0)then
                    ll=-1
                    do l=0,isct1
                      rldv=rlip(l)*dv(i,j)
                      do k=0,l
                        ll=ll+1
                        dtg(ng,ir,jr)=dtg(ng,ir,jr)+rldv
     1                   *phifor(ng,ll,i,j)*phiadj(ng,ll,i,j)
                        if(calc2d.eq.0.or.calc2d.eq.-1)exit
                      end do ! k
                    end do ! l
c use angular fluxes.
                  else if(iangflux.eq.1)then
                    if(calc2d.eq.0.or.calc2d.eq.-1)then
                      do m=1,isn
                        n=isn+1-m ! adjoint direction
                        afregc=(psifor(ng,m,i,j)+psifor(ng,m,i+1,j))
     1                   *0.5d0
                        afadjc=(psiadj(ng,n,i,j)+psiadj(ng,n,i+1,j))
     1                   *0.5d0
                        dtg(ng,ir,jr)=dtg(ng,ir,jr)+wgt(m)*afregc*afadjc
     1                   *dv(i,j)
                      end do ! m
                    else if(calc2d.eq.1)then
                      do m=1,ndir
                        n=iadir(m)
                        dtg(ng,ir,jr)=dtg(ng,ir,jr)+wgt(m)
     1                   *psifor(ng,m,i,j)*psiadj(ng,n,i,j)*dv(i,j)
                      end do ! m
                    end if
                  end if
c do this in the scattering loop to avoid looping over groups?
                  if(ifcel(ir,jr).eq.1)then
                    st=lam*phiadj(ng,0,i,j)*dv(i,j)
                    do lg=1,neg
                      dfg(ng,lg,ir,jr)=dfg(ng,lg,ir,jr)+st
     1                 *phifor(lg,0,i,j)
                    end do ! lg
                  end if
                  ll=-1
                  do l=0,isct
                    rldv=rlip(l)*dv(i,j)
                    do k=0,l
                      ll=ll+1
c this logic is here to avoid looping over groups. is it worth it?
c                     if(ifcel(ir,jr).eq.1)then
c                       lg1=1
c                       lg2=neg
c                     else
                        lg1=scgr(ng,l,mt,1)
                        lg2=scgr(ng,l,mt,2)
c                     end if
                      st=phiadj(ng,ll,i,j)*rldv
                      do lg=lg1,lg2
                        dsg(ng,lg,l,ir,jr)=dsg(ng,lg,l,ir,jr)+st
     1                   *phifor(lg,ll,i,j)
c                       if(ll.eq.0.and.ifcel(ir,jr).eq.1)then
c                         dfg(ng,lg,ir,jr)=lam*dsg(ng,lg,0,ir,jr)
c                       end if
                      end do ! lg
                      if(calc2d.eq.0.or.calc2d.eq.-1)exit
                    end do ! k
                  end do ! l
                  dqg(ng,ir,jr)=dqg(ng,ir,jr)+phiadj(ng,0,i,j)*dv(i,j)
                else if(is.eq.1)then
                  st=phifor(ng,0,i,j)*dv(i,j)
                  dtg(ng,ir,jr)=dtg(ng,ir,jr)+st
                  das(ng,ir,jr)=das(ng,ir,jr)+st*siga(ng,mt)
c                 do lg=1,neg
                  do lg=scgr(ng,0,mt,1),scgr(ng,0,mt,2)
c for reaction rates, only l=0 constributes to scattering.
                    st=phifor(lg,0,i,j)*dv(i,j)
                    dsg(ng,lg,0,ir,jr)=dsg(ng,lg,0,ir,jr)+st
                    st=sigs(ng,lg,0,mt)*st
                    if(lg.ne.ng)then
                      dis(ng,ir,jr)=dis(ng,ir,jr)+st
                      dos(lg,ir,jr)=dos(lg,ir,jr)+st
                    else
                      dss(ng,ir,jr)=dss(ng,ir,jr)+st
                    end if
                  end do ! lg
                end if
c is.eq.0 is above, inside the lg loop.
                if(is.eq.1.and.ifcel(ir,jr).eq.1)then
                  do lg=1,lng
                    st=lam*phifor(lg,0,i,j)*dv(i,j)
                    dfg(ng,lg,ir,jr)=dfg(ng,lg,ir,jr)+st
                    dfs(ng,ir,jr)=dfs(ng,ir,jr)+st*chi(ng,lg,mt)
     1               *nusigf(lg,mt)
                  end do ! lg
                end if
              end do ! ng
            end do ! j
            if(calc2d.eq.1.and.neg.gt.79)write(*,'("integrals, ",
     1       "finished i=",i5," of",i5)')i,it
          end do ! i
        end do ! jr
      end do ! ir
c
c holdover from inverse; might be useful.
c absorption inner product = total - (outscattering+selfscattering)
c     dagc2(1:neg,1:nr,1)=dtgc(1:neg,1:nr,1)-dsogc(1:neg,1:nr,1)
c    1 -dssgc(1:neg,1:nr,1)
c
c reaction rates to compare with partisn balance tables.
      if(is.eq.1)then
        dastot=0.d0
        dfstot=0.d0
        distot=0.d0
        dostot=0.d0
        dsstot=0.d0
        do ir=1,nr
          do jr=1,nz
            dastot=dastot+sum(das(1:neg,ir,jr))
            dfstot=dfstot+sum(dfs(1:neg,ir,jr))
            distot=distot+sum(dis(1:neg,ir,jr))
            dostot=dostot+sum(dos(1:neg,ir,jr))
            dsstot=dsstot+sum(dss(1:neg,ir,jr))
          end do ! jr
        end do ! ir
        write(*,'(/,"partisn balance table (in part)")')
        write(*,'(8x,"fiss source",4x,"absorption",5x,"in scatter",
     1   5x,"self scatter",3x,"out scatter")')
        do ng=1,neg
          write(*,'(i6,1p10e15.7)')ng,sum(dfs(ng,1:nr,1:nz)),
     1     sum(das(ng,1:nr,1:nz)),sum(dis(ng,1:nr,1:nz)),
     2     sum(dss(ng,1:nr,1:nz)),sum(dos(ng,1:nr,1:nz))
        end do ! ng
        write(*,'("total ",1p10e15.7)')dfstot,dastot,distot,dsstot,
     1   dostot
        write(*,'("fission source")')
        do ng=1,neg
          do jr=1,nz
            write(*,'(i6,1p10e14.6)')ng,dfs(ng,1:nr,jr)
          end do ! jr
        end do ! ng
        do jr=1,nz
          write(*,'("total ",1p10e14.6)')(sum(dfs(1:neg,i,jr)),i=1,nr)
        end do ! jr
        write(*,'("inscattering")')
        do ng=1,neg
          do jr=1,nz
            write(*,'(i6,1p10e14.6)')ng,dis(ng,1:nr,jr)
          end do ! jr
        end do ! ng
        do jr=1,nz
          write(*,'("total ",1p10e14.6)')(sum(dis(1:neg,i,jr)),i=1,nr)
        end do ! jr
        write(*,'("outscattering")')
        do ng=1,neg
          do jr=1,nz
            write(*,'(i6,1p10e14.6)')ng,dos(ng,1:nr,jr)
          end do ! jr
        end do ! ng
        do jr=1,nz
          write(*,'("total ",1p10e14.6)')(sum(dos(1:neg,i,jr)),i=1,nr)
        end do ! jr
        write(*,'("absorption")')
        do ng=1,neg
          do jr=1,nz
            write(*,'(i6,1p10e14.6)')ng,das(ng,1:nr,jr)
          end do ! jr
        end do ! ng
        do jr=1,nz
          write(*,'("total ",1p10e14.6)')(sum(das(1:neg,i,jr)),i=1,nr)
        end do ! jr
      end if
c check if the cross sections are still correct.
      go to 900
      temps=2.5301d-08
      i2lp1=0
      ifissdata=1
c for materials
      n1=1
      n2=nm
c for isotopes
c     n1=nm+1
c     n2=nm+nel+nrrx
      do ng=1,neg
        write(10,'(" *",/,1x,56("*"),"  g r o u p",i4,1x,56("*"))')ng
        write(10,'(" *",/," *",44x,"...principal cross sections...")')
        write(10,'(" *",/," *",26x,"zone",8x,"temperature",5x,"chi",
     1   9x,"nu*fission",6x,"total",
     2   7x,"absorption",/," *",24x,"no.  name")')
        do n=n1,n2
          line=" *   "
c         write(line(27:36),'(i2,"  zone",i2.2)')n,n
          write(line(27:36),'(i2,"  zone",i2.2)')n,0
          write(line(37:106),'(1p5e14.4)')temps,
     1     chi(ng,ng,n),nusigf(ng,n),sigt(ng,n),siga(ng,n)
          write(10,'(a)')trim(line)
        end do ! n
        if(ifissdata.eq.0)then
          write(10,'(" *",/," *",51x,"...fission matrices...",
     1     /," *",/," * zone",1x,"temperature",
     2     1x,"fission matrix ( gprime -->",i4," )")')ng
        else if(ifissdata.eq.1)then
          write(10,'(" *",/," *",39x,"...energy dependent fission ",
     1     "chi...",/," *",/," *   zone",7x,"temperature",
     2     4x,"chi ( gprime -->",i4," )")')ng
        end if
        do n=n1,n2
          line=" *"
          write(line(3:8),'(i6)')n
          if(sum(chi(ng,1:neg,n)).eq.0.)then
            line(21:36)='all entries zero'
            write(10,'(a)')trim(line)
            cycle
          end if
          write(line(14:27),'(1pe14.4)')temps
          mg=1
  520     do k=28,98,14
c this works for the special 2-mat cross sections (SU).
            if(ifissdata.eq.0)then
              write(line(k:k+13),'(1pe14.4)')nusigf(mg,n)*chi(ng,mg,n)
            else if(ifissdata.eq.1)then
              write(line(k:k+13),'(1pe14.4)')chi(ng,mg,n)
            end if
            mg=mg+1
            if(mg.eq.neg+1)go to 530
          end do ! k
  530     write(10,'(a)')trim(line)
          line=" *   "
          if(mg.lt.neg)go to 520
        end do ! n
        if(i2lp1.eq.0)then
          write(10,'(" *",/," *",49x,"...scattering matrices...",/,
     1     " *",52x,"(2l+1 not included)",/," *")')
        else
          write(10,'(" *",/," *",49x,"...scattering matrices...",/,
     1     " *",52x,"(2l+1 included)",/," *")')
        end if
        write(10,'(" *   zone   order  first grp",
     1   "  temperature    cross sections")')
        do n=n1,n2
          do l=0,isct
            mg=scgr(ng,l,n,2)
            if(mg.eq.0)cycle
            line=" *   "
            write(line(6:41),'(i3,6x,i1,4x,i3,5x,1pe14.4)')n,l,
     1       scgr(ng,l,n,2),temps
  620       do k=45,115,14
              write(line(k:k+13),'(1pe14.4)')sigs(ng,mg,l,n)
              mg=mg-1
              if(mg.eq.0.or.mg+1.eq.scgr(ng,l,n,1))go to 630
            end do ! k
  630       write(10,'(a)')trim(line)
            line=" *   "
            if(mg.ne.0.and.mg+1.ne.scgr(ng,l,n,1))go to 620
          end do ! l
        end do ! n
      end do ! ng
      close(10)
  900 continue
c check to see if forward moments are correct; compare with fluxp=2
      go to 950
      do ng=1,neg
        write(10,'(1x,128("*"),/," *",/,
     1   " *",41x,"flux components for group,",
     1   i4," and plane   1",/," *",/,1x,128("*"))')ng
        do l=0,isct
          write(10,'(" *",/," * component number",i3,
     1     /," *")')l+1
          k2=0
          do i=1,it,8
            k1=k2+1
            k2=min(it,k1+7)
            write(10,'(" *",6x,8("   r mesh",i5))')
     1       (k,k=k1,k2)
            write(10,'(" *")')
            write(10,'(" *",i6,8(1pe14.6))')1,(phifor(ng,l,k,1),
     1       k=k1,k2)
            write(10,'(" *")')
          end do ! i
        end do ! l
      end do ! ng
  950 continue
c
      return
      end subroutine
