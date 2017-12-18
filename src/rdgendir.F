      subroutine rdgendir
c read gendir for atomic weights
c this version sets the full zaid based on the first encounter in the
c library.
      use COMS
      implicit none
      integer i,j,k2,l0,ls,lz,l,izaid0(niso),ios,na
      character ndi_gendir_path*120,line0*120,linez*120,line*120,ct*120,
     1 zl(niso)*24,azaid*8,cl*10
c
      if(ictrl.eq.2)write(iuo,'(/,"diagnostics from rdgendir")')
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
      call getenv('NDI_GENDIR_PATH',ndi_gendir_path)
      if(ictrl.eq.2)write(iuo,'(a)')trim(ndi_gendir_path)
      if(len_trim(ndi_gendir_path).eq.0)then
        write(*,'("error. no NDI_GENDIR_PATH set.")')
        call stoponerror
      end if
      open(iug,file=ndi_gendir_path,status="old",iostat=ios)
      if(ios.ne.0)then
        write(*,'("error. trouble opening ",a)')trim(ndi_gendir_path)
        call stoponerror
      end if
c look for first instance of libname. use mendf71x for special cases.
c kynea3 is handled in rdbxslib; this is left as an example.
      if(libname(1:6).eq."kynea3")then
        cl="l=mendf71x"
      else
        cl="l="//adjustl(libname)
      end if
      l0=len_trim(cl)
   10 read(iug,'(a)',iostat=ios)line0
      if(ios.ne.0)go to 900
      if(index(line0,cl(1:l0)).eq.0)go to 10
      backspace(iug)
      na=0
c now we are at the first isotope with l=<libname>.
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
c TODO ensure that l=<libname> is in the block before "end"
      do i=1,niso
        l0=len_trim(zl(i))
        if(index(linez,zl(i)(1:l0)).ne.0.and.atwt(i).eq.0.d0)then
c this is a crappy way of bringing in special cross sections. redo.
          if(index(ndi_gendir_path,"/mt71x_h_sab/gendir").ne.0)then
            if(trim(zl(i)).eq."1001")then
              if(linez(6:10).ne."001nm")then
c             if(linez(6:10).ne."002nm")then
                cycle
              else
                write(*,'("warning. 1001 with poly is being used.")')
                write(iuo,'("warning. 1001 with poly is being used.")')
c               write(*,'("warning. 1001 with lwtr is being used.")')
c               write(iuo,'("warning. 1001 with lwtr is being used.")')
              end if
            end if
          end if
          if(ictrl.eq.2)write(iuo,'(a)')trim(line0)
          if(izaid0(i).lt.10000)then
            zaidfull(i)=" "//trim(linez)
          else
            zaidfull(i)=trim(linez)
          end if
          backspace(iug) ! reread "z=" line to check for "aw="
   30     read(iug,'(a)',iostat=ios)line
          if(ios.ne.0)go to 900
          if(index(line,"aw=").eq.0)go to 30
          l=index(line,"aw=")+3
          ct=" "
          ct(1:120)=adjustl(line(l:len_trim(line))) ! should remove blank spaces in front
          l=index(ct," ")-1 ! the first space
          read(ct(1:l),*)atwt(i)
          na=na+1
        end if
      end do ! i
      go to 20
c
  100 close(iug)
c
c convert kynea3 zaids.
c again, this is left as an example for some other special case.
      if(libname(1:6).eq."kynea3")then
        if(ictrl.eq.2)write(iuo,'("kynea3 zaids")')
        do i=1,niso
          call kynea3_iso(iuo,0,izaid0(i),azaid)
          zaidfull(i)=azaid
          if(ictrl.eq.2)write(iuo,'(i3,i8,2x,a)')i,izaid0(i),zaidfull(i)
        end do ! i
      end if
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
c eof error on ndi_gendir_path
  900 write(*,'("error. unexpected eof on ",a)')trim(ndi_gendir_path)
      write(*,'("status of atomic weights:")')
      do i=1,niso
        write(*,'(a,2x,f21.17)')trim(zaidfull(i)),atwt(i)
      end do ! i
      call stoponerror
      end subroutine
