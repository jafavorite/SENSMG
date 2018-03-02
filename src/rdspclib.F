      subroutine rdspclib
c read "special" library for atomic weights
c it reads a file containing "diagnostics from rdgendir"
c named spc_gendir
      use COMS
      implicit none
      integer i,j,k2,l0,ls,lz,l,izaid0(niso),ios,na
      character sens_data*120,line0*120,linez*120,
     1 zl(niso)*24
c
      if(ictrl.eq.2)write(iuo,'(/,"diagnostics from rdspclib")')
c get full zaid for transport and reaction-rate isotopes.
c start with the digits.
      do i=1,nel
        izaid0(i)=nint(blk(1,i))
        write(zaidfull(i),'(i5)')izaid0(i)
      end do ! i
      do i=1,nrrx
        izaid0(nel+i)=irrx(1,i)
        write(zaidfull(nel+i),'(i5)')izaid0(nel+i)
      end do ! i 
c
      atwt(1:niso)=0.d0
      zl(1:niso)=adjustl(zaidfull(1:niso))
      call getenv('SENS_DATA',sens_data)
      l=len_trim(sens_data)
      if(l.eq.0)then
        write(*,'("error. no SENS_DATA set.")')
        call stoponerror
      end if
      if(sens_data(l:l).ne."/")sens_data(l+1:l+1)="/"
      sens_data=trim(sens_data)//"spc_gendir"
      if(ictrl.eq.2)write(iuo,'(a)')trim(sens_data)
      open(iug,file=sens_data,status="old",iostat=ios)
      if(ios.ne.0)then
        write(*,'("error. trouble opening ",a)')trim(sens_data)
        call stoponerror
      end if
c this logic is the same as used in rdgendir.F
      na=0
   20 read(iug,'(a)',iostat=ios)line0
      if(na.eq.niso)go to 100
      if(ios.ne.0)go to 900
      lz=index(line0,"z=")
      if(lz.eq.0)go to 20
      lz=lz+2 ! start of zaid or space
      linez=adjustl(line0(lz:len_trim(line0)))
      ls=min(len_trim(linez),index(linez," ")) ! first space or end of line
      linez=linez(1:ls) ! this is the full zaid after "z="
c the logic here allows for repetition of isotopes in the list.
      do i=1,niso
        l0=len_trim(zl(i))
        if(index(linez,zl(i)(1:l0)).ne.0.and.atwt(i).eq.0.d0)then
          if(ictrl.eq.2)write(iuo,'(a)')trim(line0)
          if(izaid0(i).lt.10000)then
            zaidfull(i)=" "//trim(linez)
          else
            zaidfull(i)=trim(linez)
          end if
          atwt(i)=-1.d0 ! put in something for logic to work
          na=na+1
        end if
      end do ! i
      go to 20
c
c "zaids and atomic weights" line already read
  100 do i=1,niso
        read(iug,'(a)',iostat=ios)line0
        lz=len_trim(line0)
        do j=lz,1,-1
          if(line0(j:j).eq." ")then
            ls=j+1
            exit
          end if
        end do ! j
        read(line0(ls:lz),*)atwt(i)
      end do ! i
      close(iug)
c
      if(ictrl.eq.2)then
        write(iuo,'("zaids and atomic weights")')
        do i=1,niso
          write(iuo,'(a,2x,f21.17)')trim(zaidfull(i)),atwt(i)
        end do ! i
      end if
c
c compute atoms/b-cm
      if(ictrl.eq.2)write(iuo,'(/,"  mat   iso zaid",9x,
     1 "at.wt.",11x,"at.dens.",9x,"at.dens./mass dens.")')
      k2=0
      do i=1,nm
        do j=1,ncb(i)
          k2=k2+1
          blk(3,k2)=blk(2,k2)/atwt(k2)*avg*rho(i) 
          if(ictrl.eq.2)write(iuo,'(2i5,2x,a,1p3e17.9)')i,j,
     1     trim(zaidfull(k2)),atwt(k2),blk(3,k2),blk(3,k2)/rho(i)
        end do ! j
      end do ! i
      if(na.ne.niso)then
        write(*,'("error. atomic weights missing.")')
        write(iuo,'("error. atomic weights missing.")')
        call stoponerror
      end if
      if(ictrl.eq.2)write(iuo,'("end rdgendir",/)')
      return
c
c eof error on sens_data
  900 write(*,'("error. unexpected eof on ",a)')trim(sens_data)
      write(*,'("status of atomic weights:")')
      do i=1,niso
        write(*,'(a,2x,f21.17)')trim(zaidfull(i)),atwt(i)
      end do ! i
      call stoponerror
      end subroutine
