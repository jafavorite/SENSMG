      subroutine rdxsmg0(ix,niso0,nedt0)
c read partisn output file ???/snxedt for microscopic cross sections.
c this scheme gives barns.
      use COMS
      implicit none
c ix=1/2 xs1/snxedt / for/snxedt
      integer ix,niso0,nedt0
      real*8 sigt0(neg,nrrx),nusigf0(neg,nrrx),chi0(neg,neg,nrrx),
     1 siga0(neg,nrrx)
      integer i,n,i2lp1,ierr  ,ng
c
      ierr=0
      do i=1,nrrx
        if(irrx(2,i).gt.nedt0)then
          ierr=ierr+1
        end if
      end do ! i
      if(ierr.gt.0)then
        write(*,'("error.",i3," specified reaction numbers do not ",
     1   "exist in the library. library nedt0=",i4)')ierr,nedt0
        call stoponerror
      end if
c zero the cross sections.
      sigf(1:neg,0:nxs)=0.d0
      sigc(1:neg,0:nxs)=0.d0
      sige(1:neg,0:nxs)=0.d0
      sigi(1:neg,0:nxs)=0.d0
c this routine only needs to return edit cross sections and sigf,
c sigc, sige, and sigi.
c however it is interesting to compare all with xs1_macrxs.
      chi0(1:neg,1:neg,1:nrrx)=0.d0
      nusigf0(1:neg,1:nrrx)=0.d0
      sigt0(1:neg,1:nrrx)=0.d0
      siga0(1:neg,1:nrrx)=0.d0
      i2lp1=0
c
      call rdsnxedt(ix,niso0,nedt0,sigt0,nusigf0,chi0,siga0)
c convert partisn_7 edit names
      cedits(0)="flux"
      do n=1,nedt0
        if(cedits(n)(1:5).eq."mend1")then
          cedits(n)="(n,n)"
        else if(cedits(n)(1:5).eq."mend2")then
          cedits(n)="(n,n')"
        else if(cedits(n)(1:5).eq."mend3")then
          cedits(n)="(n,2n)"
        else if(cedits(n)(1:5).eq."mend4")then
          cedits(n)="(n,3n)"
        else if(cedits(n)(1:5).eq."mend5")then
          cedits(n)="(n,g)"
        else if(cedits(n)(1:5).eq."mend6")then
          cedits(n)="(n,a)"
        else if(cedits(n)(1:5).eq."mend7")then
          cedits(n)="(n,p)"
        else if(cedits(n)(1:5).eq."mend8")then
          cedits(n)="(n,d)"
        else if(cedits(n)(1:5).eq."mend9")then
          cedits(n)="(n,t)"
        else if(cedits(n)(1:6).eq."mend10")then
          cedits(n)="(n,3He)"
        else if(cedits(n)(1:6).eq."mend12")then
          cedits(n)="chi_pr"
        else if(cedits(n)(1:6).eq."mend13")then
          cedits(n)="chi_tot"
        end if
      end do ! n
c
c compare
c     do n=1,nrrx
c       do ng=1,neg
c         write(*,'(2i4,1p3e14.6)')n,ng,rrxs(ng,n),sigf(ng,2)
c         do mg=1,neg
c           write(*,'(3i4,1p3e14.6)')n,ng,mg,
c    1       chi0(mg,ng,n),chi(mg,ng,nel+n),
c    2       chi0(mg,ng,n)/chi(mg,ng,nel+n)
c         end do ! mg
c       end do ! ng
c     end do ! n
c check sigf
c     do n=1,nel
c       do ng=1,neg
c         write(*,'(2i4,1p30e30.22)')n,ng,nusigf(ng,nm+n),
c    1     sigf(ng,nm+n),nusigf(ng,nm+n)/sigf(ng,nm+n),sigc(ng,nm+n)
c       end do ! ng
c     end do ! n
c     call stoponerror
c check elas and inelas.
c sige(ng)+sigi(ng).le.sum(sigs(1:neg,ng,0,nm+n))
c     do n=1,nel
c       do ng=1,neg
c         write(*,'(2i4,1p30e17.9)')n,ng,sige(ng,nm+n),sigi(ng,nm+n),
c    1     sige(ng,nm+n)+sigi(ng,nm+n),
c    2     sum(sigs(ng,1:neg,0,nm+n)),sum(sigs(1:neg,ng,0,nm+n))
c       end do ! ng
c     end do ! n
c
      return
      end subroutine
