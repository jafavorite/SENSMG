      subroutine rdgendir
c read gendir for atomic weights
c this version sets the full zaid based on the first encounter in the
c library.
      use COMS
      implicit none
      integer i,l0,ls,lz,l,izaid0(niso),ios,na,nl,neg_ndi
      character ndi_gendir_path*120,line0*120,linez*120,line*120,ct*120,
     1 zl(niso)*24,azaid*8,cl*10,f_ndi*120
c
      write(iuo,'(/,"diagnostics from rdgendir")')
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
c ndi weight function will be read from ndi tables.
      wgtfct(1:neg,1:niso)=0.d0
      zl(1:niso)=adjustl(zaidfull(1:niso))
      call getenv('NDI_GENDIR_PATH',ndi_gendir_path)
      if(len_trim(ndi_gendir_path).eq.0)then
        write(*,'("error. no NDI_GENDIR_PATH set.")')
        call stoponerror
      end if
      write(iuo,'(a)')trim(ndi_gendir_path)
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
          write(iuo,'(a)')trim(line0)
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
            read(ct(1:l),*)neg_ndi
          end if
          if(index(line,"f=").ne.0)then
            l=index(line,"f=")+2
            ct=" "
            ct(1:120)=adjustl(line(l:len_trim(line))) ! should remove blank spaces in front
            l=index(ct," ")-1 ! the first space
            f_ndi=ct(1:l)
          end if
          if(index(line," end").ne.0.or.index(line,"end").eq.1)then
            call rdnditable(i,neg_ndi,f_ndi)
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
        write(iuo,'("kynea3 zaids")')
        do i=1,niso
          call kynea3_iso(iuo,0,izaid0(i),azaid)
          zaidfull(i)=azaid
          write(iuo,'(i3,i8,2x,a)')i,izaid0(i),zaidfull(i)
        end do ! i
      end if
c
      write(iuo,'("zaids and atomic weights")')
      do i=1,niso
        write(iuo,'(a,2x,f21.17)')trim(zaidfull(i)),atwt(i)
      end do ! i
c
c compute atoms/b-cm
      call atbcm(0,na)
      write(iuo,'("end rdgendir",/)')
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

      subroutine rdnditable(iso,neg_ndi,f_ndi)
c read ndi file for the spectrum weighting function. only needed
c for fissdata=2 (chi vector), but fissdata can be changed, so do
c this anyway.
c read scattering of order ISCT+1 for diagonal and bhs transport
c correction.
      use COMS
      implicit none
      integer iso,neg_ndi
      character f_ndi*120
      real*8 e_bounds(neg_ndi+1),wgts(neg_ndi),scat(neg_ndi,neg_ndi),tol
      parameter (tol=1.d-5)
      character line*120
      integer itrcorl,iso1,l,lg,lg1,lg2,ng,ng1,pn,pn_order,ios,
     1 ibds(neg+1)
c
      open(iug1,file=f_ndi,status="old",iostat=ios)
      if(ios.ne.0)then
        write(*,'("error in rdnditable. trouble opening ",a)')
     1   trim(f_ndi)
        write(iuo,'("error in rdnditable. trouble opening ",a)')
     1   trim(f_ndi)
        call stoponerror
      end if
      itrcorl=0
   10 read(iug1,'(a)',iostat=ios)line
      if(ios.ne.0)then
        write(*,'("error in rdnditable. e_bounds not found in ",a)')
     1   trim(f_ndi)
        write(iuo,'("error in rdnditable. e_bounds not found in ",a)')
     1   trim(f_ndi)
        go to 900
      end if
      if(index(line,"e_bounds").eq.0)go to 10
      if(index(line,"e_bounds_0").ne.0)go to 10
      read(iug1,*)e_bounds(1:neg_ndi+1)
   20 read(iug1,'(a)',iostat=ios)line
      if(ios.ne.0)then
        write(*,'("error in rdnditable. wgts not found in ",a)')
     1   trim(f_ndi)
        write(iuo,'("error in rdnditable. wgts not found in ",a)')
     1   trim(f_ndi)
        go to 900
      end if
      if(index(line,"wgts").eq.0)go to 20
      read(iug1,*)wgts(1:neg_ndi)
      if(itrcor.eq.1.or.itrcor.eq.2)then
        itrcorl=1
   30   read(iug1,'(a)',iostat=ios)line
        if(ios.ne.0)then
          write(*,'("error in rdnditable. pn_order not found in ",a)')
     1   trim(f_ndi)
          write(iuo,'("error in rdnditable. pn_order not found in ",a)')
     1   trim(f_ndi)
          go to 900
        end if
        if(index(line,"pn_order").eq.0)go to 30
        read(iug1,*)pn_order
c pn_order is indexed starting with 1, not 0
        if(pn_order-1.le.isct)then
          write(*,'("warning. no transport correction for nuclide ",
     1     a,".")')trim(zaidfull(iso))
          write(iuo,'("warning. no transport correction for nuclide ",
     1     a,".")')trim(zaidfull(iso))
          itrcorl=0
        end if
      end if ! itrcor
      if(itrcorl.eq.1)then
   40   read(iug1,'(a)',iostat=ios)line
        if(ios.ne.0)then
          write(*,'("error in rdnditable. pn_full not found in ",a)')
     1   trim(f_ndi)
          write(iuo,'("error in rdnditable. pn_full not found in ",a)')
     1   trim(f_ndi)
          go to 900
        end if
        if(index(line,"pn_full").eq.0)go to 40
        do l=0,isct+1
          read(iug1,*)pn ! should equal l
          read(iug1,*)((scat(ng,lg),ng=1,neg_ndi),lg=1,neg_ndi)
        end do ! l
      end if ! itrcorl
c
      close(iug1)
c
c find coarse-group boundaries in fine-group structure.
c this assumes the fine- and coarse-group boundaries coincide to
c within parameter tol (relative).
c ebins(1/2/3,ng) is min/max/avg.
c lg is coarse group, ng is fine group.
c if there is no collapse, ibds is not used, but still check
c consistency.
      ibds(1:neg+1)=0
      lg=1
      do ng=1,neg_ndi+1
        if(lg.le.neg)then
c         write(*,'(i6,1pe19.11,i6,e19.11)')ng,e_bounds(ng),lg,
c    1     ebins(2,lg)
          if(abs(e_bounds(ng)/ebins(2,lg)-1.d0).le.tol)then
            ibds(lg)=ng
            lg=lg+1
          end if
        else
c         write(*,'(i6,1pe19.11,i6,e19.11)')ng,e_bounds(ng),lg,
c    1     ebins(1,lg-1)
          if(abs(e_bounds(ng)/ebins(1,lg-1)-1.d0).le.tol)then
            ibds(lg)=ng
            exit
          end if
        end if
      end do ! ng
c check boundaries.
      do lg=1,neg+1
        if(ibds(lg).eq.0)then
          write(*,'("error in group collapse. boundaries are not ",
     1     "equal to within ",1pe8.2,"%. coarse group",i6,".")')
     2     tol*100.d0,lg
          write(iuo,'("error in group collapse. boundaries are not ",
     1     "equal to within ",1pe8.2,"%. coarse group",i6,".")')
     2     tol*100.d0,lg
          call stoponerror
        end if
c debug print to diagnose weight function collapse.
c       if(lg.le.neg)then
c         write(*,'(i6,1pe19.11,i6,2e19.11)')lg,ebins(2,lg),ibds(lg),
c    1     e_bounds(ibds(lg)),ebins(2,lg)-e_bounds(ibds(lg))
c       else
c         write(*,'(i6,1pe19.11,i6,2e19.11)')lg,ebins(1,neg),ibds(lg),
c    1     e_bounds(ibds(lg)),ebins(1,neg)-e_bounds(ibds(lg))
c       end if
      end do ! lg
c
c collapse weight function by simple sum.
c cross sections except scattering (below) appear to be collapsed
c according to this scheme as well.
c lg is coarse group, ng is fine group.
      if(neg_ndi.gt.neg)then
        lg=0
        do ng=ibds(1),ibds(neg+1)-1
          if(ng.eq.ibds(lg+1))then
            lg=lg+1
          end if
          if(ng.ge.ibds(lg))then
            wgtfct(lg,iso)=wgtfct(lg,iso)+wgts(ng)
c debug print to diagnose weight function collapse.
c           write(*,'(3i6,1p5e19.11)')iso,ng,lg,e_bounds(ng),wgts(ng),
c    1       ebins(1,lg),ebins(2,lg),wgtfct(lg,iso)
          end if
        end do ! ng
c TODO check this
      else if(neg_ndi.eq.neg)then
        wgtfct(1:neg,iso)=wgts(1:neg_ndi)
      end if
c
c collapse scattering if there is a transport correction (diag or
c bhs). collapse by simple sum as above. however, the bottom
c group is treated differently.
c lg1 is coarse to group, ng1 is fine to group.
c lg is coarse from group, ng is fine from group.
      iso1=nm+iso
      if(itrcorl.eq.1)then
        if(neg_ndi.gt.neg)then
          lg2=0
c in bottom coarse group, add scattering into all bottom fine groups
c (thus neg_ndi here, cap of neg in lg1= below).
c in top coarse group, add scattering into all top fine groups
c (thus 1 here).
          do ng1=1,neg_ndi ! go through all fine groups
            lg1=max(lg2,1) ! if lg2 is 0, set lg1=1
            if(ng1.eq.ibds(lg2+1))then
              lg2=min(lg2+1,neg)
              lg1=lg2
            end if
            if(ng1.ge.ibds(lg1))then
c debug print to diagnose weight function collapse.
c             write(*,'(3i6,1p5e19.11)')iso,ng1,lg1,e_bounds(ng1),
c    1         ebins(1,lg1),ebins(2,lg1)
              lg=0
              do ng=ibds(1),ibds(neg+1)-1
                if(ng.eq.ibds(lg+1))then
                  lg=lg+1
                end if
                if(ng.ge.ibds(lg))then
                  sigsp1(lg1,lg,iso1)=sigsp1(lg1,lg,iso1)+wgts(ng)
     1             *scat(ng1,ng)
c debug print to diagnose weight function collapse.
c                 write(*,'(5i6,1p8e19.11)')iso,ng,lg,ng1,lg1,
c    1             e_bounds(ng),e_bounds(ng1),wgts(ng)*scat(ng1,ng),
c    2             ebins(1:2,lg),ebins(1:2,lg1),sigsp1(lg1,lg,iso1)
                end if
              end do ! ng
            end if
          end do ! ng1
          do lg=1,neg
            sigsp1(1:neg,lg,iso1)=sigsp1(1:neg,lg,iso1)/wgtfct(lg,iso)
          end do ! lg
c TODO check this
        else if(neg_ndi.eq.neg)then
          sigsp1(1:neg,1:neg,iso1)=scat(1:neg_ndi,1:neg_ndi)
        end if ! neg_ndi
      end if ! itrcorl
      return
c
c eof error on file f_ndi.
  900 write(*,'("possible unexpected order. generalize read in ",
     1 "subroutine rdnditable.")')
      write(iuo,'("possible unexpected order. generalize read in ",
     1 "subroutine rdnditable.")')
      call stoponerror
      end subroutine
