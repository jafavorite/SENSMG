      subroutine wrdanttd(niso0)
c write ftd_inp, forward time-dependent partisn input file
c for subcritical neutron problem.  this routine assumes nm < 1,000,000.
c this routine is for computing the rise to the asymptotic Feynman Y.
c it is not yet called by the rest of the code, but is here to
c capture the knowledge gained in a prior project.
      use COMS
      implicit none
      integer niso0 ! number of isotopes in snxedt
c mt is materials for partisn input; =nm+nrrx for forward
      integer i,j,k2,l,iitm,il,l1,l2,l3,nf,mt,i_srcacc,ios
      character pfile*11,title*7,ct*26,k3c_7(7)*8,line*80,ctmp1*20
      real*8 k3a_7(7)
      logical lk
c the last 7 isotopes in kynea3 don't have atomic weights.
c these data cards are repeated in rdbxslib.F
      data k3c_7/"H3      ", "He      ", "He3     ", "He4     ",
     1 "Li      ", "Li6     ", "Li7     "/
      data k3a_7/3.0155005623387798d0, 4.0026018835759300d0,
     1 3.0149316753329298d0, 4.0026032366859798d0,
     2 6.9409329654160200d0, 6.0150722885135997d0,
     3 7.0160027500297302d0/
c
      iitm=999   ! if iitm.gt.0 it will be used for forward and adjoint.
      i_srcacc=0 ! 0/1 for default/srcacc=no
c
      if(epsi.gt.1.001d-6)write(*,'("****",/,
     1 "**** warning. epsi=",1pe9.2,/"****")')epsi
      if(isrcacc_no.eq.1.or.isrcacc_no.eq.3)then
        i_srcacc=1
        write(*,'(/,"comment. no source acceleration for ",
     1   "time-dependent calculation.")')
        write(iuo,'(/,"comment. no source acceleration for ",
     1   "time-dependent calculation.")')
      end if
c
      nf=iup
      pfile="ftd/ftd_inp"
      title="forward"
      open(nf,file=pfile,status='unknown')
      write(nf,'(5i6)')3,0,0,0,0
      write(nf,'(a)')trim(id)
      write(nf,'("input: ",a)')trim(ifile)
      write(nf,'(1x,a7," input file, time-dependent")')title
c
      if(ilnk3dnt.eq.0)then 
        mt=nm+nrrx
      else if(ilnk3dnt.eq.1)then 
        mt=nel
      else if(ilnk3dnt.eq.2)then 
        mt=nm+nrrx
      end if
      call block1(nf,mt,niso0)
      if(ilnk3dnt.eq.0)then
        call block2(nf)
      end if
c not really adjoint but use macrxs
      call block3(nf,1,0)
c
      write(nf,'("/ * * * * block iv * * * *")')
      if(ilnk3dnt.eq.0)then
        write(nf,'(" matspec=wtfrac")')
        write(nf,'(" matls=")')
        k2=0
        do i=1,nm
          write(nf,'(8x,"m",i6.6)')i
          do j=1,ncb(i)
            k2=k2+1
            ct=" "
            ct(2:25)=zaidfull(k2)(1:24)
            l=scan(ct,abc123)-1
            write(ct(l:l),'(1h")')
            l=len_trim(ct)+1
            write(ct(l:l),'(1h")')
c this is clunky but partisn 8 can't read unattached semi-colons.
            if(j.lt.ncb(i))then
              write(nf,'(12x,a,1pe15.7)')trim(ct),blk(2,k2)
            else
              write(nf,'(12x,a,1pe15.7,";")')trim(ct),blk(2,k2)
            end if
          end do ! j
        end do ! i
c
c the last 7 isotopes in kynea3 don't have atomic weights.
c if any are present, write atwt for them all.
        if(libname(1:6).eq."kynea3")then
          lk=.false.
          do i=1,niso
            do j=1,7
              if(trim(zaidfull(i)).eq.trim(k3c_7(j)))then
                write(nf,'(" atwt=")')
                write(nf,'((5x,2(a8,0pf20.16,2x)))')(k3c_7(l),k3a_7(l),
     1           l=1,7)
                lk=.true.
                exit
              end if
            end do ! j
            if(lk)exit
          end do ! i
        else if(libname(1:3).eq.'b96')then
c this is an example. do for real if needed.
          write(nf,'(" atwt=153 1.007800  159  15.9905")')
          write(nf,'("      107 234.0405  161 235.0442")')
          write(nf,'("      109 236.0458  162 238.0510")')
          write(nf,'("       40  18.9982  131  26.9818")')
          write(nf,'("      215  54.9380   98  28.0853")')
          write(nf,'("       33  62.9296              ")')
        end if
        write(nf,'(" assign=")')
        do i=1,nm
c         write(nf,'(8x,"zone",i6.6,"   m",i6.6,f14.8,";")')i,i,rho(i)
          write(nf,'(8x,"zn",i6.6,"   m",i6.6,f14.8,";")')i,i,rho(i)
        end do ! i
      else if(ilnk3dnt.eq.1)then
        write(nf,'(3x,"matspec=",/,6x,"atdens")')
        write(nf,'(3x,"matls=")')
        do i=1,nel
          ctmp1='"'//trim(adjustl(zaidfull(i)))//'"'
          write(nf,'(6x,"m",i6.6,2x,a,"  1.;")')i,trim(ctmp1)
        end do ! i 
        write(nf,'(3x,"assign=",/,6x,"matls")')
      else if(ilnk3dnt.eq.2)then
c open redoin input file and advance to block 4.
        open(iui,file=ifile,status='old',iostat=ios)
   10   read(iui,'(a)',iostat=ios)line
        if(line(1:15).ne."&block_4_arrays")go to 10
        il=0
        ctmp1=" "
        do
          read(iui,'(a)',iostat=ios)line
          if(index(line,"/").ne.0)then
            if(il.eq.3)write(nf,'(6x,";")')
            exit
          else if(index(line,"matspec=").ne.0)then
            if(il.eq.3)write(nf,'(6x,";")')
            il=1
            write(nf,'(a)')trim(line) ! matspec
          else if(index(line,"assign=").ne.0)then
            if(il.eq.3)write(nf,'(6x,";")')
            il=2
            if(index(line,"matls").ne.0)then ! assign=matls
              l1=index(line,'"')
              l2=index(line(l1+1:l1+20),'"')+l1
              line(l1:l1)=" "
              line(l2:l2)=" "
            end if
            write(nf,'(a)')trim(line) ! assign
          else if(index(line,"matls=").ne.0)then
            il=3
            write(nf,'(a)')trim(line) ! matls
          else if(il.eq.1)then
            write(nf,'(a)')trim(line)
          else if(il.eq.2)then
            if(index(line,"matls").ne.0)then ! assign=matls
              l1=index(line,'"')
              l2=index(line(l1+1:l1+20),'"')+l1
              line(l1:l1)=" "
              line(l2:l2)=" "
              write(nf,'(a)')trim(line)
            else
              write(nf,'(a,a)')trim(line),";"
            end if
          else if(il.eq.3)then
            l1=index(line,'"')
            l2=index(line(l1+1:l1+20),'"')+l1
            l3=l2-l1+1
            if(ctmp1.eq." ")then ! 1st material
              ctmp1(1:l3)=line(l1:l2)
              write(nf,'(6x,a)')trim(ctmp1)
              write(nf,'(14x,a)')
     1         adjustl(trim(line(l2+1:len_trim(line))))
            else if(ctmp1(1:l3).eq.line(l1:l2))then ! same material
              write(nf,'(14x,a)')
     1         adjustl(trim(line(l2+1:len_trim(line))))
            else ! new material
              ctmp1(1:l3)=line(l1:l2)
              write(nf,'(6x,";")')
              if(index(line,"matspec=").ne.0.or.
     1         index(line,"assign=").ne.0)then
                backspace(iui)
                cycle
              end if
              write(nf,'(6x,a)')trim(ctmp1)
              write(nf,'(14x,a)')
     1         adjustl(trim(line(l2+1:len_trim(line))))
            end if
          end if
        end do
        close(iui)
      end if ! ilnk3dnt
      write(nf,'(" t")')
c
      write(nf,'("/ * * * * block v * * * *")')
      if(icalc.eq.0.or.icalc.eq.3.or.icalc.eq.4)then
        write(nf,'(" ievt=",i1," isct=",i3)')0,isct
      else if(icalc.eq.1.or.icalc.eq.2)then
        write(nf,'(" ievt=",i1," isct=",i3)')icalc,isct
      end if
c JAF -- testing
c     write(iuo,'("warning. reflecting boundaries.")')
c     write(nf,'(" ibl=1 ibr=1 / reflective")')
c JAF -- testing
      write(nf,'(" epsi=",1pe9.2)')epsi
      if(cellsol.eq.-1)then
        write(nf,'(" balp=",i1)')1
      else if(cellsol.eq.0)then
        write(nf,'(" balp=",i1)')1
        write(nf,'(" npeg=2 cellsol=0")')
      else if(cellsol.eq.1)then
        write(nf,'(" balp=",i1)')0
        write(nf,'(" cellsol=1")')
      end if
      if(icalc.eq.1.or.icalc.eq.2)then
        write(nf,'(" norm=1.0")')
      end if
      if(icalc.eq.2)then
        write(nf,'(" aflxfrm=",i1)')aflxfrm
      end if
      write(nf,'(" rmflux=1")')
      if(iangflux.eq.0.or.cellsol.eq.1)then
        write(nf,'(" raflux=0")')
      else if(iangflux.eq.1.and.(cellsol.eq.-1.or.cellsol.eq.0))then
        write(nf,'(" raflux=1")')
      end if
      if(idbgw.eq.1)write(nf,'(" fluxp=2")')
      if(itrcor.eq.1)then
        write(nf,'(" trcor=diag")')
      else if(itrcor.eq.2)then
        write(nf,'(" trcor=bhs")')
      else if(itrcor.eq.3)then
        write(nf,'(" trcor=cesaro")')
      end if
c print spherical harmonics (search for "Spherical harmonics",
c compare with last debug print in subroutine mcon.)
c     write(nf,'(" extras=17r0 1 f0")')
c debug moments
c     write(nf,'(" fluxp=2")')
c debug lnk3dnt
c     write(nf,'(" wrlnk3d=1")')
      if(iitm.gt.0)write(nf,'(" iitm=",i3)')iitm
      write(nf,'(" iitl=0 oitm=9999")')
      if(i_srcacc.eq.1)write(nf,'(" srcacc=no")')
      if(nofxup.eq.1)write(nf,'(" nofxup=1")')
c get surface angular fluxes.
c iplane, jplane are used for partisn versions before 8_27.
c DEBUG MOMENTS (2 lines)
c     write(nf,'(" angp=1 fluxp=2")')
c     write(nf,'(" initaf=2")')
      if(calc2d.eq.1)then
        if(iplane.gt.-1)then ! 0 is caught in rdctrl
          write(nf,'(" asleft=",i6," asrite=",i6)')iindex(iplane)+1,
     1     iindex(iplane)+1
        else
          if(iaflux.eq.1)then
            write(nf,'(" afluxx=")')
            write(nf,'(8x,8i6)')(iindex(i)+1,i=1,nr)
          else
            write(*,'("warning. no asleft/asrite or afluxx.")')
            write(iuo,'("warning. no asleft/asrite or afluxx.")')
          end if
        end if
        if(jplane.gt.-1)then
          write(nf,'(" asbott=",i6," astop= ",i6)')jindex(jplane)+1,
     1     jindex(jplane)+1
        else
          if(iaflux.eq.1)then
            write(nf,'(" afluxy=")')
            write(nf,'(8x,8i6)')(jindex(j)+1,j=0,nz)
          else
            write(*,'("warning. no asbott/astop or afluxy.")')
            write(iuo,'("warning. no asbott/astop or afluxyx.")')
          end if
        end if
      end if
c
      write(nf,'(" ith=0 xsectp=0")')
c forward time-dependent:
c kynea3 velocities are in cm/us.
c ndi velocities are in cm/sh.
c 1 sh=10^-8 s, 1 us=100 sh
c time-dependent angular fluxes. check how this works.
c     write(nf,'(" angflx=1")')
      write(nf,'(" timedep=1 timeq=1 itmdis=0")')
c maximum time=4096 us, initial step size 0.1 us
      if(libname(1:6).eq."kynea3")then
        write(nf,'(" t0=0. ts=4096.")')
c       write(nf,'(" delti=0.10 deltmax=100.")')
c       write(nf,'(" delti=0.10 deltmax=1.")')
        write(nf,'(" delti=0.001 deltmax=100.")')
c default behaviour but specify for clarity.
        write(nf,'(" stimes= 0. 4096.")')
        write(nf,'(" samp=   1. 1.")')
c John uses this
c       write(nf,'(" deltmin=1 efact=1")')
      else
        write(nf,'(" t0=0. ts=409600.")')
c       write(nf,'(" delti=10. deltmax=10000.")')
c       write(nf,'(" delti=10. deltmax=100.")')
        write(nf,'(" delti=0.1 deltmax=10000.")')
c default behaviour but specify for clarity.
        write(nf,'(" stimes= 0. 409600.")')
        write(nf,'(" samp=   1. 1.")')
      end if
      call fixedsrc(nf,0)
      write(nf,'(" t")')
c
      close(nf)
c
      return
      end subroutine
