      subroutine rxnrateedit(rxnrate,rxnratet)
      use COMS
c calculate reaction rates and ratios.
      implicit none
      real*8 rxnrate(neg,nrrx),rxnratet(nrrx+nflux)
      real*8 fluxt,st
      integer i,j,k,ir,jr,k2,n,ng
c
      fluxt=0.d0
      rxnrate(1:neg,1:nrrx)=0.d0
      rxnratet(1:nrrx)=0.d0
      if(ilnk3dnt.eq.0)then
        do k=1,nedpoints
          ir=edpoints(k)
          jr=edpoints(nedpoints+k)
          do i=iindex(ir-1)+1,iindex(ir)
            do j=jindex(jr-1)+1,jindex(jr)
              do ng=1,neg
                st=fmom(ng,0,i,j)*dv(i,j)
                fluxt=fluxt+st
                do n=1,nrrx
                  rxnrate(ng,n)=rxnrate(ng,n)+st*rrxs(ng,n)
                  rxnratet(n)=rxnratet(n)+st*rrxs(ng,n)
                end do ! n
              end do ! ng
            end do ! j
          end do ! i
        end do ! k
      else if(ilnk3dnt.gt.0)then
        do k=1,nedpoints
          i=edpoints(k)
          j=edpoints(nedpoints+k)
          do ng=1,neg
            st=fmom(ng,0,i,j)*dv(i,j)
            fluxt=fluxt+st
            do n=1,nrrx
              rxnrate(ng,n)=rxnrate(ng,n)+st*rrxs(ng,n)
              rxnratet(n)=rxnratet(n)+st*rrxs(ng,n)
            end do ! n
          end do ! ng
        end do ! k
      end if
c
      if(nflux.eq.1)then
        if(icalc.eq.0)then ! icalc=3 or 4 not allowed.
          rxnratet(nrrx+nflux)=1.d0
          rrxs(1:neg,nrrx+nflux)=0.d0
        else if(icalc.eq.1.or.icalc.eq.2)then
          rxnratet(nrrx+nflux)=fluxt
          rrxs(1:neg,nrrx+nflux)=1.d0
        end if
      end if
c
      write(*,'(/,"group-sum reaction rates")')
      write(iuo,'(/,"group-sum reaction rates")')
      do n=1,nrrx
        k2=nel+n
        write(*,'(i3,i7,2x,a,2x,i5,1pe14.6)')n,nm+n,
     1   trim(zaidfull(k2)),irrx(2,n),rxnratet(n)
        write(iuo,'(i3,i7,2x,a,2x,i5,1pe14.6)')n,nm+n,
     1   trim(zaidfull(k2)),irrx(2,n),rxnratet(n)
      end do ! n
      write(*,'(12x,"flux",14x,1pe14.6)')fluxt
      write(iuo,'(12x,"flux",14x,1pe14.6)')fluxt
      write(*,'("reaction-rate ratios")')
      write(iuo,'("reaction-rate ratios")')
      do n=1,nrrr
        rxnratio(n)=rxnratet(irri(1,n))/rxnratet(irri(2,n))
        write(*,'(i3,i3,"/",i3,1pe14.6)')n,irri(1,n),irri(2,n),
     1   rxnratio(n)
        write(iuo,'(i3,i3,"/",i3,1pe14.6)')n,irri(1,n),irri(2,n),
     1   rxnratio(n)
      end do ! n
c
      return
      end subroutine
