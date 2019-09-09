      subroutine atbcm(iw,na)
c compute atoms/b-cm
      use COMS
      implicit none
      integer iw,na
      integer ia,i,j,k2,k3
      character line*17
c
c sometimes just print
      if(iw.eq.1)go to 100
      k2=0
      do i=1,nm
        k3=k2
        if(blk(2,k2+1).gt.0.d0.and.rho(i).eq.0.d0)then
c ilnk3dnt should be 1 in this case, all atom densities 1.
          ia=0
          blk(3,k3+1:k3+ncb(i))=blk(2,k3+1:k3+ncb(i))
          blk(2,k3+1:k3+ncb(i))=0.d0
        else if(blk(2,k2+1).lt.0.d0.and.rho(i).lt.0.d0)then
c given wgt. fracs and mass density
          ia=1
          rho(i)=-rho(i)
          rhoa(i)=0.d0
          blk(2,k3+1:k3+ncb(i))=-blk(2,k3+1:k3+ncb(i))
        else if(blk(2,k2+1).gt.0.d0.and.rho(i).gt.0.d0)then
c given at. fracs and atom density
          ia=2
          rhoa(i)=rho(i)
          rho(i)=0.d0
        else if(blk(2,k2+1).lt.0.d0.and.rho(i).gt.0.d0)then
c given wgt. fracs and atom density
          ia=3
          rhoa(i)=rho(i)
          rho(i)=0.d0
          blk(2,k3+1:k3+ncb(i))=-blk(2,k3+1:k3+ncb(i))
        else if(blk(2,k2+1).gt.0.d0.and.rho(i).lt.0.d0)then
c given at. fracs and mass density
          ia=4
          rho(i)=-rho(i)
          rhoa(i)=0.d0
        end if
   10   do j=1,ncb(i)
          k2=k2+1
          if(abs(ia).eq.1)then
            blk(3,k2)=blk(2,k2)/atwt(k2)*avg*rho(i) 
            if(ia.eq.1)rhoa(i)=rhoa(i)+blk(3,k2)
          else if(abs(ia).eq.2)then
            blk(3,k2)=blk(2,k2)*rhoa(i)
            blk(2,k2)=atwt(k2)*blk(3,k2)/avg
            if(ia.eq.2)rho(i)=rho(i)+blk(2,k2)
          else if(ia.eq.3)then
            rho(i)=rho(i)+blk(2,k2)*avg/atwt(k2)
          else if(ia.eq.4)then
            rhoa(i)=rhoa(i)+blk(2,k2)*atwt(k2)
          end if
        end do ! j
        if(abs(ia).eq.2)then
          blk(2,k3+1:k3+ncb(i))=blk(2,k3+1:k3+ncb(i))/rho(i)
        else if(ia.eq.3)then
          rho(i)=rhoa(i)/rho(i)
          k2=k3
          ia=-1
          go to 10
        else if(ia.eq.4)then
          rhoa(i)=rho(i)*avg/rhoa(i)
          k2=k3
          ia=-2
          go to 10
        end if
      end do ! i
c
c come here to just print
  100 continue
      write(iuo,'(/,3x,"mat    iso   zaid",9x,"at.wt.",11x,"at.dens.",
     1 9x,"wt.frac.",9x,"at.dens./mass_dens.")')
      k2=0
      do i=1,nm
        do j=1,ncb(i)
          k2=k2+1
          if(ilnk3dnt.eq.0)then
            write(iuo,'(2i6,4x,a11,1p4e17.9)')i,j,
     1       trim(adjustl(zaidfull(k2))),atwt(k2),blk(3,k2),
     2       blk(2,k2),blk(3,k2)/rho(i)
          else if(ilnk3dnt.eq.1)then
            write(iuo,'(2i6,4x,a11,1p2e17.9)')i,j,
     1       trim(adjustl(zaidfull(k2))),atwt(k2),blk(3,k2)
          else if(ilnk3dnt.eq.2)then
            write(iuo,'(2i6,4x,a11,1pe17.9,34x,e17.9)')i,j,
     1       trim(adjustl(zaidfull(k2))),atwt(k2),blk(3,k2)
          end if
        end do ! j
      end do ! i
      line="Avogadro's number"
      write(iuo,'(/,2x,a,1pe22.15)')line,avg
      if(ilnk3dnt.eq.0)then
        write(iuo,'(3x,"mat  atom_density",5x,"mass_density")')
        do i=1,nm
          write(iuo,'(i6,1p2e17.9)')i,rhoa(i),rho(i)
        end do ! nm
      else if(ilnk3dnt.eq.1)then
        write(iuo,'(3x,"no material densities printed.")')
      else if(ilnk3dnt.eq.2)then
        write(iuo,'(3x,"mat  atom_density/mass_density")')
        do i=1,nm
          write(iuo,'(i6,1pe17.9)')i,rhoa(i)
        end do ! nm
      end if
c
      if(na.ne.niso)then
        write(*,'("error. atomic weights missing.")')
        write(iuo,'("error. atomic weights missing.")')
        call stoponerror
      end if
      return
      end subroutine
