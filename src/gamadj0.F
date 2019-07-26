      subroutine gamadj0(niso0)
      use COMS
c set up generalized adjoint source in zero'th iteration.
      implicit none
      integer niso0 ! number of isotopes in snxedt
      real*8 rxnratet(nrrx+nflux),gams(neg,0:nmom-1,it,jt),xfiss
c not needed in this routine but needed to pass to rxnrateedit.
c not sure how to do this efficiently.
      real*8 rxnrate(neg,nrrx)
      integer i,j,k,iconv,ir,jr,n,ng,ncalcs,ig
c ig=1/2 solve gamma* for reactions/solve gamma* for reaction
c differences
      if(icalc.eq.0.or.icalc.eq.3.or.icalc.eq.4)then
        ig=2
      else if(icalc.eq.1.or.icalc.eq.2)then
        ig=2
      end if
c
c solve gamma* for reactions
      if(ig.eq.1)then
        ncalcs=nrrx
c solve gamma* for reaction differences
      else if(ig.eq.2)then
        ncalcs=nrrr
      end if
c
c edit reaction rates.
      if(nrrx.gt.0)then
        call rxnrateedit(rxnrate,rxnratet)
      end if
c
      iconv=0
      do n=1,ncalcs
        gams(1:neg,0:nmom-1,1:it,1:jt)=0.d0
        if(ilnk3dnt.eq.0)then
          do k=1,nedpoints
            jr=edpoints(k)/nr+1
            if(mod(edpoints(k),nr).eq.0)jr=jr-1
            ir=edpoints(k)-(jr-1)*nr
            do i=iindex(ir-1)+1,iindex(ir)
              do j=jindex(jr-1)+1,jindex(jr)
                do ng=1,neg
                  if(ig.eq.1)then
                    gams(ng,0,i,j)=rrxs(ng,n)/rxnratet(n)
                  else if(ig.eq.2)then
c for source problem, this is
c Stacey's Eq. (1.1.15) divided by G3 (G3 is the ratio)
                    gams(ng,0,i,j)=
     1               rrxs(ng,irri(1,n))/rxnratet(irri(1,n))
     2               -rrxs(ng,irri(2,n))/rxnratet(irri(2,n))
                  end if
                end do ! ng
              end do ! j
            end do ! i
          end do ! k
        else if(ilnk3dnt.gt.0)then
          do k=1,nedpoints
c n=(j-1)it + i
            j=edpoints(k)/it+1
            if(mod(edpoints(k),it).eq.0)j=j-1
            i=edpoints(k)-(j-1)*it
            do ng=1,neg
              if(ig.eq.1)then
                gams(ng,0,i,j)=rrxs(ng,n)/rxnratet(n)
              else if(ig.eq.2)then
c for source problem, this is
c Stacey's Eq. (1.1.15) divided by G3 (G3 is the ratio)
                gams(ng,0,i,j)=
     1           rrxs(ng,irri(1,n))/rxnratet(irri(1,n))
     2           -rrxs(ng,irri(2,n))/rxnratet(irri(2,n))
              end if
            end do ! ng
          end do ! k
        end if
c check <S0*, psi0>=0 (or 1 for absolute reaction rate
c in fixed-source)
        xfiss=0.d0
        do ir=1,nr
          do jr=1,nz
            do i=iindex(ir-1)+1,iindex(ir)
              do j=jindex(jr-1)+1,jindex(jr)
                do ng=1,neg
                  xfiss=xfiss+gams(ng,0,i,j)*fmom(ng,0,i,j)*dv(i,j)
                end do ! ng
              end do ! j
            end do ! i
          end do ! jr
        end do ! ir
        if(n.eq.1)then
          write(*,'("")')
          write(iuo,'("")')
        end if
        if((icalc.eq.0.or.icalc.eq.3.or.icalc.eq.4).and.
     1   irri(2,n).eq.nrrx+nflux)then
          write(*,'("check <S0* psi0>=1")')
          write(iuo,'("check <S0* psi0>=1")')
        else
          write(*,'("check <S0* psi0>=0")')
          write(iuo,'("check <S0* psi0>=0")')
        end if
        write(*,'(i3," <S0* psi0>=",1pe14.6)')n,xfiss
        write(iuo,'(i3," <S0* psi0>=",1pe14.6)')n,xfiss
        call wrdantga(n,iconv,ig,niso0,gams)
      end do ! n
c
      if(ncalcs.eq.0)then
        open(iut,file='stopconverged',status='unknown')
        close (iut)
      end if
      return
      end subroutine
