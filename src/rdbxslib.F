      subroutine rdbxslib(niso0,nisoaw,nedt0)
c read bxslib for atomic weights
      use COMS
      implicit none
      integer niso0,nisoaw,nedt0
      character hname*8,huse(2)*8
      integer*4 ivers,ia(9)
      character hsetid(12)*8,hisonm(niso0)*8,edname(nedt0)*8,
     1 atisoc(nisoaw)*8
      real*8 chin(neg),emax(neg),emin,atisor(2*nisoaw),k3a_7(145:151)
      integer nspi(niso0)
      integer i,j,k2,l0,n,izaid0(niso),na,ni,idum
      character afile*12,zl(niso)*24,azaid*8,k3c_7(145:151)*8,ctmp*11
c the last 7 isotopes in kynea3 don't have atomic weights.
c these data cards are repeated in wrdantnm.F
      data k3c_7/"H3      ", "He      ", "He3     ", "He4     ", 
     1 "Li      ", "Li6     ", "Li7     "/
      data k3a_7/3.0155005623387798d0, 4.0026018835759300d0,
     1 3.0149316753329298d0, 4.0026032366859798d0,
     2 6.9409329654160200d0, 6.0150722885135997d0,
     3 7.0160027500297302d0/
c
      if(ictrl.eq.2)then
        write(iuo,'(/,"diagnostics from rdbxslib")')
        write(iuo,'("libname=",a)')trim(libname)
      end if
c get full zaid for transport and reaction-rate isotopes.
      do i=1,nel
        izaid0(i)=nint(blk(1,i))
      end do ! i
      do i=1,nrrx
        izaid0(nel+i)=irrx(1,i)
      end do ! i 
c
c convert kynea3 zaids.
      if(libname(1:6).eq."kynea3")then
        if(ictrl.eq.2)write(iuo,'("kynea3 zaids")')
c get full zaid for transport and reaction-rate isotopes.
        do i=1,niso
          call kynea3_iso(iuo,0,izaid0(i),azaid)
          zaidfull(i)=azaid
          if(ictrl.eq.2)write(iuo,'(i3,i8,2x,a)')i,izaid0(i),
     1     trim(zaidfull(i))
        end do ! i
      end if
c
      atwt(1:niso)=0.d0
      zl(1:niso)=adjustl(zaidfull(1:niso))
      afile="for/"//libname
c file is opened as unformatted because it has record marks.
      open(iug,file=afile,status='old',position='rewind',
     1 form='unformatted',access='sequential',
     2 err=900)
      read(iug,err=910)hname,(huse(i),i=1,2),ivers
      read(iug,err=910)(ia(i),i=1,9)
      if(idbgw.eq.1)then
        write(41,'(a)')trim(afile)
        write(41,'(10i8)')(ia(i),i=1,9)
      endif
      if(ia(1).ne.neg)go to 795
      if(ia(2).ne.niso0)go to 800
      if(ia(4)-3.ne.nedt0)go to 805
      if(ia(7).ne.nisoaw)go to 810
c
      read(iug,err=910)(hsetid(i),i=1,12),(hisonm(i),i=1,ia(2)),
     1 (edname(i),i=1,ia(4)-3),(chin(n),n=1,ia(1)),(vel(n),n=1,ia(1)),
     2 (emax(n),n=1,ia(1)),emin,(nspi(i),i=1,ia(2))
c     write(*,'("hsetid")')
c     do i=1,12
c       write(*,'(i3,2x,a8)')i,hsetid(i)
c     end do ! i
c     write(*,'("hisonm")')
c     do i=1,ia(2)
c       write(*,'(i3,2x,a)')i,trim(hisonm(i))
c     end do ! i
      do i=1,ia(8)
        read(iug,err=910)idum
      end do ! i
c
c now we are at the isotope label/atomic weight pair set
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
c convert kynea3 velocities from cm/(0.01 s) to cm/s.
c     if(libname(1:6).eq."kynea3")then
c       vel(1:neg)=vel(1:neg)*0.01d0
c     end if
c
      if(ictrl.eq.2)then
        write(iuo,'("zaids and atomic weights")')
        do i=1,niso
          ctmp=adjustl(zaidfull(i))
          write(iuo,'(a11,2x,f21.17)')ctmp,atwt(i)
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
          if(ictrl.eq.2)write(iuo,'(2i5,2x,a11,1p3e17.9)')i,j,
     1     trim(zaidfull(k2)),atwt(k2),blk(3,k2),blk(3,k2)/rho(i)
        end do ! j
      end do ! i
      if(na.ne.niso)then
        write(*,'("error. atomic weights missing.")')
        write(iuo,'("error. atomic weights missing.")')
        call stoponerror
      end if
      if(ictrl.eq.2)write(iuo,'("end rdbxslib",/)')
      return
c
c eof error on bxslib
  795 write(*,'("error. ia(1).ne.neg in ",a,"; ia(1)=",i8)')
     1 trim(afile),ia(1)
      call stoponerror
  800 write(*,'("error. ia(2).ne.niso0 in ",a,"; ia(2)=",i8)')
     1 trim(afile),ia(2)
      call stoponerror
  805 write(*,'("error. ia(4)-3.ne.nedt0 in ",a,"; ia(4)=",i8)')
     1 trim(afile),ia(4)
      call stoponerror
  810 write(*,'("error. ia(7).ne.nisoaw in ",a,"; ia(7)=",i8)')
     1 trim(afile),ia(7)
      call stoponerror
  900 write(*,'("error opening ",a,".")')trim(afile)
      call stoponerror
  910 write(*,'("error reading ",a,".")')trim(afile)
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
            return
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
            return
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
      return
      end subroutine
