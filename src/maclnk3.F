      subroutine maclnk3(m0,ir,jr)
c return macroscopic cross sections for fine mesh ir,jr using idclnk3
c and denlnk3. m0 is 0 or 1 for negative or positive side of the
c interface, or it can be any other index. remember cross sections
c after index nm are for isotopes.
      use COMS
      implicit none
      integer m0,ir,jr
      integer lg,ng,l,n,iso
      real*8 chiden
c     real*8 chitot(2)
c
      sigt(1:neg,m0)=0.d0
      chi(1:neg,1:neg,m0)=0.d0
      nusigf(1:neg,m0)=0.d0
      sigs(1:neg,1:neg,0:isct,m0)=0.d0
c scgr(ng,l,m0,1/2) are the top and bottom groups scattering into ng.
      scgr(1:neg,0:isct,m0,1)=neg
      scgr(1:neg,0:isct,m0,2)=1
      if(ir.lt.1.or.ir.gt.it.or.jr.lt.1.or.jr.gt.jt)then
        go to 990
      end if
c     chitot(1:2)=0.d0
      chiden=0.d0
      do n=1,nmxsp
        if(idclnk3(ir,jr,n).eq.0)exit
        iso=nm+idclnk3(ir,jr,n)
        sigt(1:neg,m0)=sigt(1:neg,m0)+sigt(1:neg,iso)
     1   *denlnk3(ir,jr,n)
        nusigf(1:neg,m0)=nusigf(1:neg,m0)+nusigf(1:neg,iso)
     1   *denlnk3(ir,jr,n)
        sigs(1:neg,1:neg,0:isct,m0)=sigs(1:neg,1:neg,0:isct,m0)
     1   +sigs(1:neg,1:neg,0:isct,iso)*denlnk3(ir,jr,n)
        do ng=1,neg
          do l=0,isct
            scgr(ng,l,m0,1)=min(scgr(ng,l,m0,1),scgr(ng,l,iso,1))
            scgr(ng,l,m0,2)=max(scgr(ng,l,m0,2),scgr(ng,l,iso,2))
          end do ! l
          if(fissdata.eq.0.or.fissdata.eq.1)then
            chi(1:neg,ng,m0)=chi(1:neg,ng,m0)+chi(1:neg,ng,iso)
     1       *denlnk3(ir,jr,n)*nusigf(ng,iso)
          else if(fissdata.eq.2)then
            do lg=1,neg
              chi(1:neg,ng,m0)=chi(1:neg,ng,m0)+chi(1:neg,lg,iso)
     1         *denlnk3(ir,jr,n)*nusigf(lg,iso)*wgtfct(lg,iso)
            end do ! lg
            chiden=chiden+denlnk3(ir,jr,n)*nusigf(ng,iso)
     1       *wgtfct(ng,iso)
          end if
        end do ! ng
      end do ! n
      do ng=1,neg
        if(fissdata.eq.0.or.fissdata.eq.1)then
          if(nusigf(ng,m0).ne.0.d0)then
            chi(1:neg,ng,m0)=chi(1:neg,ng,m0)/nusigf(ng,m0)
          else
            chi(1:neg,ng,m0)=0.d0
          end if
        else if(fissdata.eq.2)then
          if(chiden.ne.0.d0)then
            chi(1:neg,ng,m0)=chi(1:neg,ng,m0)/chiden
          else
            chi(1:neg,ng,m0)=0.d0
          end if
        end if
      end do ! ng
c check chi matrix normalization
c     if(iter.eq.0.and.ictrl.eq.4.and.ir.eq.1.and.jr.eq.1)then
c       write(iuo,'("chi matrix normalization (fissdata=",i1,")",
c    1   /,"  mat/iso  ng sum-out-to-ng  sum-ng-to-out")')fissdata
c       do n=m0,m0
c         if(sum(chi(1:neg,1:neg,n)).eq.0.d0)cycle
c         do ng=1,neg
c           if(fissdata.eq.0.or.fissdata.eq.1)then
c             write(iuo,'(2i6,1p2e15.7)')n,ng,sum(chi(ng,1:neg,n)),
c    1         sum(chi(1:neg,ng,n))
c             chitot(1)=chitot(1)+sum(chi(ng,1:neg,n))
c             chitot(2)=chitot(2)+sum(chi(1:neg,ng,n))
c           else if(fissdata.eq.2)then
c for fissdata=2, the chi vector is stored on the diagonal of the 
c chi matrix; see rdmacrxs.
c             write(iuo,'(2i6,1p2e15.7)')n,ng,chi(ng,ng,n),
c    1         chi(ng,ng,n)
c             chitot(1)=chitot(1)+chi(ng,ng,n)
c             chitot(2)=chitot(2)+chi(ng,ng,n)
c           end if
c         end do ! ng
c         write(iuo,'(i6,"   tot",1p2e15.7)')n,chitot(1:2)
c       end do ! n
c     end if
c     if(ir.eq.it.or.ir.eq.758.or.ir.eq.759.or.ir.eq.760)then
c       write(*,'(2i4,1pe14.6,2i4,e14.6)')ir,m0,sigt(12,m0),
c    1   scgr(12,0,m0,1),scgr(12,0,m0,2),sigs(12,12,0,m0)
c     end if
c     write(99,'(3i5,1p100e14.6)')m0,ir,jr,sigt(1:neg,m0)
  990 return
      end subroutine
