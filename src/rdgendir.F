      subroutine rdgendir(iwrite)
c read gendir for atomic weights
c this version sets the full zaid based on the first encounter in the
c library.
      use COMS
      implicit none
      integer iwrite
      integer i,l0,ls,lz,l,izaid0(niso),ios,na,nl,ng_ndi
      character ndi_gendir_path*120,line0*120,linez*120,line*120,ct*120,
     1 zl(niso)*24,azaid*8,cl*10,f_ndi*120
c
      if(iwrite.eq.1)write(iuo,'(/,"diagnostics from rdgendir")')
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
      call getenv('NDI_GENDIR_PATH',ndi_gendir_path)
      if(len_trim(ndi_gendir_path).eq.0)then
        write(*,'("error. no NDI_GENDIR_PATH set.")')
        call stoponerror
      end if
      if(iwrite.eq.1)write(iuo,'(a)')trim(ndi_gendir_path)
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
        nl=0 ! count lines read
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
          if(iwrite.eq.1)write(iuo,'(a)')trim(line0)
          if(izaid0(i).lt.10000)then
            zaidfull(i)=" "//trim(linez)
          else
            zaidfull(i)=trim(linez)
          end if
          backspace(iug) ! reread "z=" line to check for "aw="
   30     read(iug,'(a)',iostat=ios)line
          if(ios.ne.0)go to 900
          nl=nl+1
c look for aw=, ng=, f=.
c no else here because these can be on the same line.
          if(index(line,"aw=").ne.0)then
            l=index(line,"aw=")+3
            ct=" "
            ct(1:120)=adjustl(line(l:len_trim(line))) ! should remove blank spaces in front
            l=index(ct," ")-1 ! the first space
            read(ct(1:l),*)atwt(i)
          end if
          if(index(line,"ng=").ne.0)then
            l=index(line,"ng=")+3
            ct=" "
            ct(1:120)=adjustl(line(l:len_trim(line))) ! should remove blank spaces in front
            l=index(ct," ")-1 ! the first space
            read(ct(1:l),*)ng_ndi
          end if
          if(index(line,"f=").ne.0)then
            l=index(line,"f=")+2
            ct=" "
            ct(1:120)=adjustl(line(l:len_trim(line))) ! should remove blank spaces in front
            l=index(ct," ")-1 ! the first space
            f_ndi=ct(1:l)
          end if
          if(index(line," end").ne.0.or.index(line,"end").eq.1)then
            call rdspecwgt(i,ng_ndi,f_ndi)
            do l=1,nl-1
              backspace(iug)
            end do ! l
            na=na+1
            cycle
          end if
          go to 30
        end if
      end do ! i
      if(na.lt.niso)go to 20
c
      close(iug)
c
c convert kynea3 zaids.
c again, this is left as an example for some other special case.
      if(libname(1:6).eq."kynea3")then
        if(iwrite.eq.1)write(iuo,'("kynea3 zaids")')
        do i=1,niso
          call kynea3_iso(iuo,0,izaid0(i),azaid)
          zaidfull(i)=azaid
          if(iwrite.eq.1)write(iuo,'(i3,i8,2x,a)')i,izaid0(i),
     1     zaidfull(i)
        end do ! i
      end if
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
      if(iwrite.eq.1)write(iuo,'("end rdgendir",/)')
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

      subroutine rdspecwgt(iso,ng_ndi,f_ndi)
c read ndi files for the spectrum weighting function. only needed
c for fissdata=2 (chi vector), but fissdata can be changed, so do
c this anyway.
      use COMS
      implicit none
      integer iso,ng_ndi
      character f_ndi*120
      real*8 e_bounds(ng_ndi+1),wgts(ng_ndi)
      character line*120
      integer ng,lg,ios
c
      open(iug1,file=f_ndi,status="old",iostat=ios)
      if(ios.ne.0)then
        write(*,'("error. trouble opening ",a)')trim(f_ndi)
        call stoponerror
      end if
   10 read(iug1,'(a)',iostat=ios)line
      if(index(line,"e_bounds").eq.0)go to 10
      if(index(line,"e_bounds_0").ne.0)go to 10
      read(iug1,*)e_bounds(1:ng_ndi+1)
   20 read(iug1,'(a)',iostat=ios)line
      if(index(line,"wgts").eq.0)go to 20
      read(iug1,*)wgts(1:ng_ndi)
      close(iug1)
c
c collapse by simple sum and linear interpolation.
c ebins(1/2/3,ng) is min/max/avg.
      if(ng_ndi.gt.neg)then
        lg=1
        do ng=1,ng_ndi
          if(e_bounds(ng).le.ebins(2,lg).and.
     1     e_bounds(ng+1).ge.ebins(1,lg))then
            wgtfct(lg,iso)=wgtfct(lg,iso)+wgts(ng)
c           write(*,'(2i4,1p1e24.16,2e14.6,a)')ng,lg,e_bounds(ng),
c    1       ebins(1:2,lg)," all"
          else if(e_bounds(ng).gt.ebins(2,lg).and.
     1     e_bounds(ng+1).lt.ebins(2,lg))then
            wgtfct(lg,iso)=wgtfct(lg,iso)+wgts(ng)
     1       *(ebins(2,lg)-e_bounds(ng+1))/(e_bounds(ng)-e_bounds(ng+1))
c           write(*,'(2i4,1p1e24.16,2e14.6,a)')ng,lg,e_bounds(ng),
c    1       ebins(1:2,lg)," top"
          else if(e_bounds(ng).gt.ebins(1,lg).and.
     1     e_bounds(ng+1).lt.ebins(1,lg))then
            wgtfct(lg,iso)=wgtfct(lg,iso)+wgts(ng)
     1       *(e_bounds(ng)-ebins(1,lg))/(e_bounds(ng)-e_bounds(ng+1))
c           write(*,'(2i4,1p1e24.16,2e14.6,a)')ng,lg,e_bounds(ng),
c    1       ebins(1:2,lg)," bot"
          else
c           write(*,'(2i4,1p1e24.16,2e14.6,a)')ng,lg,e_bounds(ng),
c    1       ebins(1:2,lg)," outside"
          end if
          if(e_bounds(ng+1).le.ebins(1,lg))lg=lg+1
          if(lg.gt.neg)exit
        end do ! ng
c TODO check this
      else if(ng_ndi.eq.neg)then
        wgtfct(1:neg,iso)=wgts(1:ng_ndi)
      end if
      return
      end subroutine
