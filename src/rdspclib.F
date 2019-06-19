      subroutine rdspclib(iwrite)
c read "special" library for atomic weights
c it reads a file containing "diagnostics from rdgendir"
c named spc_gendir
      use COMS
      implicit none
      integer iwrite
      integer i,j,l0,ls,lz,itmp,izaid0(niso),ios,na,ng
      character sens_dataf*120,line0*120,linez*120,zl(niso)*24
c
      if(iwrite.eq.1)write(iuo,'(/,"diagnostics from rdspclib")')
c get full zaid for transport and reaction-rate isotopes.
c start with the digits.
      do i=1,nel
        izaid0(i)=nint(blk(1,i))
        write(zaidfull(i),'(i7)')izaid0(i)
      end do ! i
      do i=1,nrrx
        izaid0(nel+i)=irrx(1,i)
        write(zaidfull(nel+i),'(i7)')izaid0(nel+i)
      end do ! i 
c
      atwt(1:niso)=0.d0
c ndi weight function will be read.
      wgtfct(1:neg,1:niso)=0.d0
      zl(1:niso)=adjustl(zaidfull(1:niso))
      sens_dataf=trim(sens_data)//"spc_gendir"
      if(iwrite.eq.1)write(iuo,'(a)')trim(sens_dataf)
      open(iug,file=sens_dataf,status="old",iostat=ios)
      if(ios.ne.0)then
        write(*,'("error. trouble opening ",a)')trim(sens_dataf)
        call stoponerror
      end if
c this logic is the same as used in rdgendir.F
      na=0
   20 read(iug,'(a)',iostat=ios)line0
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
          if(iwrite.eq.1)write(iuo,'(a)')trim(line0)
          if(izaid0(i).lt.10000)then
            zaidfull(i)=" "//trim(linez)
          else
            zaidfull(i)=trim(linez)
          end if
          atwt(i)=-1.d0 ! put in something for logic to work
          na=na+1
        end if
      end do ! i
      if(na.lt.niso)go to 20
c
      read(iug,'(a)',iostat=ios)line0 ! "zaids and atomic weights"
      do i=1,niso
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
c
   30 read(iug,'(a)',iostat=ios)line0
      if(ios.ne.0)then
        write(*,'("warning. spectrum weight functions do not exist ",
     1   "in ",a,". using 1.")')trim(sens_dataf)
        write(iuo,'("warning. spectrum weight functions do not exist ",
     1   "in ",a,". using 1.")')trim(sens_dataf)
        wgtfct(1:neg,1:niso)=1.d0
        go to 40
      end if
      if(line0(1:24).ne."spectrum weight function")go to 30
      do i=1,niso
        read(iug,'(a)',iostat=ios)linez
        if(ios.ne.0)go to 40
        do j=1,niso
          l0=len_trim(zl(j))
          if(index(linez,zl(j)(1:l0)).ne.0)then
            do ng=1,neg
              read(iug,*)itmp,wgtfct(ng,j)
            end do ! ng
            exit
          end if
        end do ! j
      end do ! i
   40 close(iug)
c
      if(iwrite.eq.1)then
        write(iuo,'("zaids and atomic weights")')
        do i=1,niso
          write(iuo,'(a,2x,f21.17)')trim(zaidfull(i)),atwt(i)
        end do ! i
      end if
c
c compute atoms/b-cm
      call atbcm(iwrite,na)
      if(iwrite.eq.1)write(iuo,'("end rdspclib",/)')
      return
c
c eof error on sens_dataf
  900 write(*,'("error. unexpected eof on ",a)')trim(sens_dataf)
      write(*,'("status of atomic weights:")')
      do i=1,niso
        write(*,'(a,2x,f21.17)')trim(zaidfull(i)),atwt(i)
      end do ! i
      call stoponerror
      end subroutine
