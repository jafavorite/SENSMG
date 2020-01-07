      subroutine rdbxslib(niso0,nisoaw,nedt0)
c read bxslib for atomic weights
      use COMS
      implicit none
      integer niso0,nisoaw,nedt0
      character hname*8,huse(2)*8
      integer nwd
      parameter (nwd=9)
      integer*4 ivers,ia(nwd)
      character hsetid(12)*8,hisonm(niso0)*8,edname(nedt0)*8,
     1 atisoc(nisoaw)*8
      real*8 chin(neg),emax(neg),emin,atisor(2*nisoaw),k3a_7(145:151)
      integer nspi(niso0)
      integer i,iso,iso_last,l0,j,n,izaid0(niso),na,ni,n8,idum
      character afile*12,zl(niso)*24,azaid*8,k3c_7(145:151)*8,ctmp*11
      logical lf
c the last 7 isotopes in kynea3 don't have atomic weights.
c these data cards are repeated in wrdantnm.F
      data k3c_7/"H3      ", "He      ", "He3     ", "He4     ", 
     1 "Li      ", "Li6     ", "Li7     "/
      data k3a_7/3.0155005623387798d0, 4.0026018835759300d0,
     1 3.0149316753329298d0, 4.0026032366859798d0,
     2 6.9409329654160200d0, 6.0150722885135997d0,
     3 7.0160027500297302d0/
c
      write(iuo,'(/,"diagnostics from rdbxslib")')
      write(iuo,'("libname=",a)')trim(libname)
c get full zaid for transport and reaction-rate isotopes.
      do i=1,nel
        izaid0(i)=nint(blk(1,i))
      end do ! i
      do i=1,nrrx
        izaid0(nel+i)=irrx(1,i)
      end do ! i 
c
      if(libname(1:6).eq."kynea3")then
c convert kynea3 zaids.
        write(iuo,'("kynea3 zaids")')
c get full zaid for transport and reaction-rate isotopes.
        do i=1,niso
          call kynea3_iso(iuo,0,izaid0(i),azaid)
          zaidfull(i)=azaid
          write(iuo,'(i3,i8,2x,a)')i,izaid0(i),trim(zaidfull(i))
        end do ! i
      else if(libname(1:5).eq."scale")then
c convert scale zaids.
        write(iuo,'("scale zaids")')
c get full zaid for transport and reaction-rate isotopes.
        do i=1,niso
          call scale_iso(iuo,0,izaid0(i),azaid)
          zaidfull(i)=azaid
          write(iuo,'(i3,i8,2x,a)')i,izaid0(i),trim(zaidfull(i))
        end do ! i
      end if
c
      atwt(1:niso)=0.d0
c there is no ndi weight function to read.
      wgtfct(1:neg,1:niso)=1.d0
      zl(1:niso)=adjustl(zaidfull(1:niso))
      afile="for/"//libname
c file is opened as unformatted because it has record marks.
      open(iug,file=afile,status='old',position='rewind',
     1 form='unformatted',access='sequential',
     2 err=900)
      read(iug,err=910)hname,(huse(i),i=1,2),ivers
      read(iug,err=910)(ia(i),i=1,nwd)
      if(idbgw.eq.1)then
        write(41,'(a)')trim(afile)
        write(41,'(10i8)')(ia(i),i=1,nwd)
      endif
      if(ia(1).ne.neg)go to 795
      if(ia(2).ne.niso0)go to 800
      if(ia(4)-3.ne.nedt0)go to 805
      if(ia(7).ne.nisoaw)go to 810
c
      read(iug,err=910)(hsetid(i),i=1,12),(hisonm(i),i=1,ia(2)),
     1 (edname(i),i=1,ia(4)-3),(chin(n),n=1,ia(1)),(vel(n),n=1,ia(1)),
     2 (emax(n),n=1,ia(1)),emin,(nspi(i),i=1,ia(2))
c ia(8) should equal sum(nspi(1:ia(2)))
c this is how to find out what edits and isotopes are in the file.
c     write(*,'("hsetid")')
c     do i=1,12
c       write(*,'(i3,2x,a8)')i,hsetid(i)
c     end do ! i
c     write(*,'("hisonm")')
c     do i=1,ia(2)
c       write(*,'(i3,2x,a)')i,trim(hisonm(i))
c     end do ! i
c if no diagonal or bhs transport correction, no need to read
c everything.
      if(itrcor.ne.1.and.itrcor.ne.2)then
        sigsp1(1:neg,1:neg,0:nxs)=0.d0
        do i=1,ia(8)
          read(iug,err=910)idum
        end do ! i
      else
c read isct+1 scattering for diagonal or bhs transport correction.
c the logic here allows for repetition of isotopes in the list.
        n8=0
        na=0
        do n=1,ia(2)
          lf=.false.
          if(n8.ne.sum(nspi(1:n-1)))then
            write(*,'("error in rdbxslib.",3i6,1x,a)')n,n8,
     1       sum(nspi(1:n-1)),trim(hisonm(n))
            write(iuo,'("error in rdbxslib.",3i6,1x,a)')n,n8,
     1       sum(nspi(1:n-1)),trim(hisonm(n))
          end if
          if(na.eq.niso)go to 50
          do i=1,niso
            iso=nm+i
            l0=len_trim(zl(i))
c this test may need to be modified for a bxslib that is not kynea3.
            if(trim(hisonm(n)).eq.zl(i)(1:l0).and.
     1       sigsp1(1,1,iso).eq.0.d0)then
              if(nspi(n)-1.le.isct)then
c write warning and advance through this nuclide if not done already.
                if(.not.lf)then
                  write(*,'("warning. no transport correction for ",
     1             "isotope ",a,".")')trim(zaidfull(i))
                  write(iuo,'("warning. no transport correction for ",
     1             "isotope ",a,".")')trim(zaidfull(i))
                  do j=1,nspi(n)
                    read(iug,err=910)idum
                    n8=n8+1
                  end do ! j
                  lf=.true.
                end if
              else
c read cross sections or assign if repeated nuclide
                if(.not.lf)then
c                 write(*,'("read cross sections ",a)')trim(hisonm(n))
                  do j=0,isct
                    read(iug,err=910)idum
                    n8=n8+1
c                   write(*,'(4i6)')n,i,j,n8
                  end do ! j
                  call rdkynea3(iug,isct+1,neg,ia(3),ia(4),ia(5),ia(9),
     1             sigsp1(1:neg,1:neg,iso))
                  n8=n8+1
c                 write(*,'(4i6," isct+1")')n,i,isct+1,n8
                  do j=isct+2,nspi(n)-1
                    read(iug,err=910)idum
                    n8=n8+1
c                   write(*,'(4i6)')n,i,j,n8
                  end do ! j
                  lf=.true.
                else
c                 write(*,'("assign cross sections ",a)')trim(hisonm(n))
                  sigsp1(1:neg,1:neg,iso)=sigsp1(1:neg,1:neg,iso_last)
                end if
              end if
              na=na+1
            else if(i.eq.niso.and..not.lf)then
c just advance through this nuclide.
              do j=1,nspi(n)
                read(iug,err=910)idum
                n8=n8+1
c               write(*,'(4i6," advancing")')n,i,j,n8
              end do ! j
              lf=.true.
            end if
            iso_last=iso
          end do ! i
        end do ! n
c advance to atomic weights
   50   do j=1,ia(8)-n8
          read(iug,err=910)idum
          n8=n8+1
        end do ! j
      end if ! itrcor
c     write(*,'("n8,ia(8)",2i6)')n8,ia(8)
c
c now we are at the isotope label/atomic weight pair set.
c the logic here allows for repetition of isotopes in the list.
      read(iug,err=910)(atisor(n),n=1,2*ia(7))
      na=0
      ni=0
      do n=1,2*ia(7),2
        ni=ni+1
        if(na.eq.niso)go to 100
        atisoc(ni)=transfer(atisor(n),hisonm(1))
c       write(*,'(a8,1pe15.6)')atisoc(ni),atisor(n+1)
        do i=1,niso
          l0=len_trim(zl(i))
c this test may need to be modified for a bxslib that is not kynea3.
          if(trim(atisoc(ni)).eq.zl(i)(1:l0).and.atwt(i).eq.0.d0)then
            atwt(i)=atisor(n+1)
            na=na+1
          end if
        end do ! i
      end do ! n
c the last 7 isotopes in kynea3 don't have atomic weights.
      if(libname(1:6).eq."kynea3")then
        do ni=145,151
          if(na.eq.niso)go to 100
          do i=1,niso
            l0=len_trim(zl(i))
            if(trim(k3c_7(ni)).eq.zl(i)(1:l0).and.atwt(i).eq.0.d0)then
              atwt(i)=k3a_7(ni)
              na=na+1
            end if
          end do ! i
        end do ! ni
      end if
c
  100 close(iug)
c
c don't convert kynea3 and scale energies because they are not passed
c out from here.
c convert kynea3 and scale velocities from cm/(0.01 s) to cm/s.
      if(libname(1:6).eq."kynea3".or.libname(1:5).eq."scale")then
        vel(1:neg)=vel(1:neg)*0.01d0
      end if
c
      write(iuo,'("zaids and atomic weights")')
      do i=1,niso
        ctmp=adjustl(zaidfull(i))
        write(iuo,'(a11,2x,f21.17)')ctmp,atwt(i)
      end do ! i
c
c compute atoms/b-cm
      call atbcm(0,na)
      write(iuo,'("end rdbxslib",/)')
      return
c
c eof error on bxslib
  795 write(*,'("error. ia(1).ne.neg in ",a,"; ia(1)=",i8)')
     1 trim(afile),ia(1)
      write(*,'("neg=",i8)')neg
      call stoponerror
  800 write(*,'("error. ia(2).ne.niso0 in ",a,"; ia(2)=",i8)')
     1 trim(afile),ia(2)
      write(*,'("niso0=",i8)')niso0
      call stoponerror
  805 write(*,'("error. ia(4)-3.ne.nedt0 in ",a,"; ia(4)=",i8)')
     1 trim(afile),ia(4)
      write(*,'("nedt0=",i8)')nedt0
      call stoponerror
  810 write(*,'("error. ia(7).ne.nisoaw in ",a,"; ia(7)=",i8)')
     1 trim(afile),ia(7)
      write(*,'("nisoaw=",i8)')nisoaw
      call stoponerror
  900 write(*,'("error opening ",a,".")')trim(afile)
      write(iuo,'("error opening ",a,".")')trim(afile)
      call stoponerror
  910 write(*,'("error reading ",a,".")')trim(afile)
      write(iuo,'("error reading ",a,".")')trim(afile)
      write(*,'("status of atomic weights:")')
      do i=1,niso
        write(*,'(a,2x,f21.17)')trim(zaidfull(i)),atwt(i)
      end do ! i
      call stoponerror
      end subroutine
c
      subroutine kynea3_iso(iuo,iw,izaid0,azaid)
c kynea3 uses alphanumeric text, not zaids.
      implicit none
c iw=0/1 return azaid/izaid0
      integer iuo,iw,izaid0
      character azaid*8
      integer tablength
      parameter (tablength=151)
      integer iname(tablength),i
      character cname(tablength)*8
c
c duplicates:
c H1(H2O) 101001  | H1          1001
c Be9(th) 104009  | Be9         4009
c C12       6012  | C12(gph)  106012
c Am242    95242  | Am242m    195242
c H1 is H1(CH2) from Vitamin-B6.
      data cname/
     1 "H       ", "H1(H2O) ", "H1      ", "H2      ", "Be      ",
     2 "Be9(th) ", "Be9     ", "B       ", "B10     ", "B11     ",
     3 "C       ", "C12     ", "C12(gph)", "N       ", "N14     ",
     4 "N15     ", "O       ", "O16     ", "O17     ", "F       ",
     5 "F19     ", "Na      ", "Na23    ", "Mg      ", "Al      ",
     6 "Al27    ", "Si      ", "P       ", "P31     ", "S       ",
     7 "S32     ", "Cl      ", "K       ", "Ca      ", "Ti      ",
     8 "V       ", "Cr      ", "Cr50    ", "Cr52    ", "Cr53    ",
     9 "Cr54    ", "Mn      ", "Mn55    ", "Fe      ", "Fe54    ",
     1 "Fe56    ", "Fe57    ", "Fe58    ", "Co      ", "Co59    ",
     2 "Ni      ", "Ni58    ", "Ni60    ", "Ni61    ", "Ni62    ",
     3 "Ni64    ", "Cu      ", "Cu63    ", "Cu65    ", "Ga      ",
     4 "Y       ", "Y89     ", "Zr      ", "Nb      ", "Nb93    ",
     5 "Mo      ", "Ag      ", "Ag107   ", "Ag109   ", "Cd      ",
     6 "In      ", "Sn      ", "Ba      ", "Ba138   ", "Eu      ",
     7 "Eu151   ", "Eu152   ", "Eu153   ", "Eu154   ", "Eu155   ",
     8 "Hf      ", "Hf174   ", "Hf176   ", "Hf177   ", "Hf178   ",
     9 "Hf179   ", "Hf180   ", "Ta      ", "Ta181   ", "Ta182   ",
     1 "W       ", "W182    ", "W183    ", "W184    ", "W186    ",
     2 "Re      ", "Re185   ", "Re187   ", "Au      ", "Au197   ",
     3 "Pb      ", "Pb206   ", "Pb207   ", "Pb208   ", "Bi      ",
     4 "Bi209   ", "Th      ", "Th230   ", "Th232   ", "Pa      ",
     5 "Pa231   ", "Pa233   ", "U       ", "U232    ", "U233    ",
     6 "U234    ", "U235    ", "U236    ", "U237    ", "U238    ",
     7 "Np237   ", "Np238   ", "Np239   ", "Pu236   ", "Pu237   ",
     8 "Pu238   ", "Pu239   ", "Pu240   ", "Pu241   ", "Pu242   ",
     9 "Pu243   ", "Pu244   ", "Am241   ", "Am242   ", "Am242m  ",
     1 "Am243   ", "Cm241   ", "Cm242   ", "Cm243   ", "Cm244   ",
     2 "Cm245   ", "Cm246   ", "Cm247   ", "Cm248   ", "H3      ",
     3 "He      ", "He3     ", "He4     ", "Li      ", "Li6     ",
     4 "Li7     "/
      data iname/
     1 001000, 101001, 001001, 001002, 004000,
     2 104009, 004009, 005000, 005010, 005011,
     3 006000, 006012, 106012, 007000, 007014,
     4 007015, 008000, 008016, 008017, 009000,
     5 009019, 011000, 011023, 012000, 013000,
     6 013027, 014000, 015000, 015031, 016000,
     7 016032, 017000, 019000, 020000, 022000,
     8 023000, 024000, 024050, 024052, 024053,
     9 024054, 025000, 025055, 026000, 026054,
     1 026056, 026057, 026058, 027000, 027059,
     2 028000, 028058, 028060, 028061, 028062,
     3 028064, 029000, 029063, 029065, 031000,
     4 039000, 039089, 040000, 041000, 041093,
     5 042000, 047000, 047107, 047109, 048000,
     6 049000, 050000, 056000, 056138, 063000,
     7 063151, 063152, 063153, 063154, 063155,
     8 072000, 072174, 072176, 072177, 072178,
     9 072179, 072180, 073000, 073181, 073182,
     1 074000, 074182, 074183, 074184, 074186,
     2 075000, 075185, 075187, 079000, 079197,
     3 082000, 082206, 082207, 082208, 083000,
     4 083209, 090000, 090230, 090232, 091000,
     5 091231, 091233, 092000, 092232, 092233,
     6 092234, 092235, 092236, 092237, 092238,
     7 093237, 093238, 093239, 094236, 094237,
     8 094238, 094239, 094240, 094241, 094242,
     9 094243, 094244, 095241, 095242, 195242,
     1 095243, 096241, 096242, 096243, 096244,
     2 096245, 096246, 096247, 096248, 001003,
     3 002000, 002003, 002004, 003000, 003006,
     4 003007/
c
      if(iw.eq.0)then
        azaid="        "
      else if(iw.eq.1)then
        izaid0=0
      end if
      do i=1,tablength
        if(iw.eq.0)then
          if(izaid0.eq.iname(i))then
            azaid=cname(i)
            if(iname(i).eq.101001)then
              write(*,'("warning. H1(H2O) is being used.")')
              write(iuo,'("warning. H1(H2O) is being used.")')
            else if(iname(i).eq.104009)then
              write(*,'("warning. Be9(th) is being used.")')
              write(iuo,'("warning. Be9(th) is being used.")')
            else if(iname(i).eq.106012)then
              write(*,'("warning. C(gph) is being used.")')
              write(iuo,'("warning. C(gph) is being used.")')
            else if(iname(i).eq.195242)then
              write(*,'("warning. Am242m is being used.")')
              write(iuo,'("warning. Am242m is being used.")')
            end if
            go to 999
          end if
        else if(iw.eq.1)then
          if(azaid.eq.cname(i))then
            izaid0=iname(i)
            if(iname(i).eq.101001)then
              write(*,'("warning. H1(H2O) is being used.")')
              write(iuo,'("warning. H1(H2O) is being used.")')
            else if(iname(i).eq.104009)then
              write(*,'("warning. Be9(th) is being used.")')
              write(iuo,'("warning. Be9(th) is being used.")')
            else if(iname(i).eq.106012)then
              write(*,'("warning. C(gph) is being used.")')
              write(iuo,'("warning. C(gph) is being used.")')
            else if(iname(i).eq.195242)then
              write(*,'("warning. Am242m is being used.")')
              write(iuo,'("warning. Am242m is being used.")')
            end if
            go to 999
          end if
        end if
      end do ! i
      if(iw.eq.0)then
        write(*,'("error. no name for zaid",i7,
     1   " in subroutine kynea3_iso.")')izaid0
      else if(iw.eq.1)then
        write(*,'("error. no zaid for name",a8,
     1   " in subroutine kynea3_iso.")')azaid
      end if
      call stoponerror
c
  999 return
      end subroutine
c
      subroutine scale_iso(iuo,iw,izaid0,azaid)
c for scale, we have to know the order of the isotopes in the table:
c iso1 C12
c iso2 O16
c iso3 H1
c iso4 Pu239
c iso5 Pu240
      implicit none
c iw is unused
      integer iuo,iw,izaid0
      character azaid*8
      integer tablength
      parameter (tablength=5)
      integer iname(tablength),i
      character cname(tablength)*8
c
c this logic has no way of dealing with duplicates:
c H1(H2O) and H1, Be9(th) and Be9, C12 and C12(gph), Am242 and Am242m
c defaults are H1, Be9, C12, and Am242
c H1 is supposed to be H1(CH2) from Vitamin-B6.
      data cname/
     1 "C12     ", "O16     ", "H1      ", "Pu239   ", "Pu240   "/
c    1 "iso1    ", "iso2    ", "iso3    ", "iso4    ", "iso5    "/
      data iname/
     1 006000, 008016, 001001, 094239, 094240/
c
      azaid="        "
      do i=1,tablength
        if(izaid0.eq.iname(i))then
          azaid=cname(i)
          return
        end if
      end do ! i
      write(*,'("error. no name for zaid",i7,
     1 " in subroutine scale_iso.")')izaid0
      call stoponerror
c
      return
      end subroutine

      subroutine rdkynea3(iug,pn,neg,ihm,iht,ihs,i2lp1,scat)
      implicit none
      integer iug,pn,neg,ihm,iht,ihs,i2lp1
      real*8 scat(neg,neg)
      integer i,ng,lg
      real*8 c(ihm,neg),rl
      read(iug)((c(i,ng),i=1,ihm),ng=1,neg)
      if(i2lp1.eq.1)then
        rl=1.d0/real(2*pn+1,8)
        c(1:ihm,1:neg)=c(1:ihm,1:neg)*rl
      end if
      do ng=1,neg
        do lg=ihs+(neg-1),ihs-(neg-1),-1
          if(ng-(lg-ihs).lt.1.or.ng-(lg-ihs).gt.neg.or.
     1     lg.le.iht.or.lg.gt.ihm)cycle
          scat(ng,ng-(lg-ihs))=c(lg,ng)
        end do !lg
      end do ! ng
      return
      end subroutine
