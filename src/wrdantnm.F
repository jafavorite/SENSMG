      subroutine wrdantnm(ith,niso0)
c write for_inp and adj_inp, forward and adjoint partisn input files
c for neutron problem.  this routine assumes nm < 1,000,000.  
      use COMS
      implicit none
      integer ith ! =0/1 forward/adjoint
      integer niso0 ! number of isotopes in snxedt
c mt is materials for partisn input; =nm+nrrx for forward
c n0 is niso to use
      integer i,j,k2,l,iitm,il,l1,l2,l3,nf,mt,n0,i_srcacc,npn,ios
      character pfile*11,title*7,ct*26,k3c_7(7)*8,line*80,ctmp1*20
      real*8 k3a_7(7)
      logical lk
      character edit7(18)*8,edit8(20)*8,edit(20)*8
      data edit7/ "chi", "nusigf", "total", "abs", "chi", "mend1",
     1 "mend2", "mend3", "mend4", "mend5", "mend6", "mend7",
     2 "mend8", "mend9", "mend10", "n-fiss", "mend12", "mend13"/
      data edit8/ "chi", "nusigf", "total", "abs", "chi", "(n,n)",
     1 "(n,n')", "(n,2n)", "(n,3n)", "(n,g)", "(n,p)", "(n,a)",
     2 "(n,f)", "(n,n')f", "(n,2n)f", "(n,F)", "chi_pr", "chi_tot",
     3 "(n,d)", "(n,t)"/
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
      if(libname(1:6).eq."kynea3")then
        npn=4
        edit(1:npn)=edit7(1:npn)
      else if(libname(1:5).eq."scale")then
c check this
        npn=4
        edit(1:npn)=edit7(1:npn)
      else
c if the supported version isn't here, don't do edits.
        if(partisnver(1:1).eq."5".or.partisnver(1:1).eq."7")then
          npn=18 ! partisn 5_97 and 7_72
          edit(1:npn)=edit7(1:npn)
        else if(partisnver(1:4).eq."8_27".or.partisnver(1:4).eq."8_29"
     1   .or.partisnver(1:4).eq."8_xx")then
          npn=20 ! partisn 8_27 and later
          edit(1:npn)=edit8(1:npn)
c be ready for changes in partisn.
        else if(partisnver(1:1).eq."8")then
          npn=20 ! partisn 8_??
          edit(1:npn)=edit8(1:npn)
c unsupported partisn.
        else
          npn=0
        end if
      end if
c
      if(epsi.gt.1.001d-6)write(*,'("****",/,
     1 "**** warning. epsi=",1pe9.2,/"****")')epsi
      if(ith.eq.0.and.(isrcacc_no.eq.1.or.isrcacc_no.eq.3))then
        i_srcacc=1
        write(*,'(/,"comment. no source acceleration for forward ",
     1   "calculation.")')
        write(iuo,'(/,"comment. no source acceleration for forward ",
     1   "calculation.")')
      else if(ith.eq.1.and.(isrcacc_no.eq.2.or.isrcacc_no.eq.3))then
        i_srcacc=1
        write(*,'(/,"comment. no source acceleration for adjoint ",
     1   "calculation.")')
        write(iuo,'(/,"comment. no source acceleration for adjoint ",
     1   "calculation.")')
      end if
c
      nf=iup
      if(ith.eq.0)then
        pfile="for/for_inp"
        title="forward"
        if(libname(1:6).eq.'kynea3')then
          n0=151
        else if(libname(1:5).eq.'scale')then
c check this
          n0=5
        else if(libname(1:3).eq.'b96')then
          n0=225
        else if(libname(1:7).eq.'special')then
c assume mendf71x
          n0=17570
        else
          n0=0
        end if
      else if(ith.eq.1)then
        pfile="adj/adj_inp"
        title="adjoint"
        n0=niso0
      end if
      open(nf,file=pfile,status='unknown')
      write(nf,'(5i6)')3,0,0,0,0
      write(nf,'(a)')trim(id)
      write(nf,'("input: ",a)')trim(ifile)
      if(icalc.eq.0.or.icalc.eq.3.or.icalc.eq.4)then
        write(nf,'(1x,a7," input file, fixed-source")')title
      else if(icalc.eq.1)then
        write(nf,'(1x,a7," input file, keff")')title
      else if(icalc.eq.2)then
        write(nf,'(1x,a7," input file, alpha")')title
      end if
c
      if(ilnk3dnt.eq.0)then 
        mt=nm+nrrx
      else if(ilnk3dnt.eq.1)then 
        mt=nel
      else if(ilnk3dnt.eq.2)then 
        mt=nm+nrrx
      end if
      call block1(nf,mt,n0)
      if(ilnk3dnt.eq.0)then
        call block2(nf)
      end if
      call block3(nf,ith,0)
c
      write(nf,'("/ * * * * block iv * * * *")')
      if(ilnk3dnt.eq.0)then
        write(nf,'(" matspec=wtfrac")')
        if(ith.eq.0)then
c materials are given as weight fractions; edit isotopes
c are given as atdens=1.
          do i=2,nm
            write(nf,'("         wtfrac")')
          end do ! i
          do i=1,nrrx
            write(nf,'("         atdens")')
          end do ! i
        end if
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
c edit isotopes
        if(ith.eq.0)then
          if(nrrx.gt.0)write(nf,'("/ edit isotopes")')
          do i=1,nrrx
            ct=" "
            ct(2:25)=zaidfull(nel+i)(1:24)
            l=scan(ct,abc123)-1
            write(ct(l:l),'(1h")')
            l=len_trim(ct)+1
            write(ct(l:l),'(1h")')
            write(nf,'(8x,"m",i6.6,2x,a,"  1.;")')nm+i,trim(ct)
          end do ! i
        end if
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
        write(nf,'(" ievt=",i1," isct=",i1)')0,isct
      else if(icalc.eq.1.or.icalc.eq.2)then
        write(nf,'(" ievt=",i1," isct=",i1)')icalc,isct
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
c this is not needed
        if(partisnver(1:4).eq."8_xx".and.ith.eq.0)then
          write(nf,'(" rstream=1")')
        end if
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
c different solve parameters possible for forward and adjoint files 
c forward:
      if(ith.eq.0)then
        write(nf,'(" ith=0 xsectp=2")')
c     
c adjoint:
      else if(ith.eq.1)then
        write(nf,'(" ith=1 xsectp=0")')
      end if
      if(icalc.eq.0.or.icalc.eq.3.or.icalc.eq.4)then
        call fixedsrc(nf,ith)
      end if
      write(nf,'(" t")')
c
c no edits for unsupported partisn.
c no edits for adjoint for eigenvalues.
c no edits for redoin/lnk3dnt.
c no edits for too many materials.
      if(npn.eq.0)then
        write(*,'("comment. no edits for unsupported partisn.")')
        write(iuo,'("comment. no edits for unsupported partisn.")')
        go to 100
      end if
      if(ith.eq.1.and.(icalc.eq.1.or.icalc.eq.2))then
        write(*,'("comment. no edits for adjoint for eigenvalues.")')
        write(iuo,'("comment. no edits for adjoint for eigenvalues.")')
        go to 100
      end if
      if(ilnk3dnt.gt.0)then
        write(*,'("comment. no edits for redoin/lnk3dnt.")')
        write(iuo,'("comment. no edits for redoin/lnk3dnt.")')
        go to 100
      end if
c max determined with one problem using partisn 8.31. may be smaller.
      if(nrrx.eq.0.and.nm.gt.460)then
        write(*,'("comment. no edits--too many materials.")')
        write(iuo,'("comment. no edits--too many materials.")')
        go to 100
      end if
c this line for nse01 problem
c     write(*,'("warning. skipping edit in wrdantnm.")')
c     write(iuo,'("warning. skipping edit in wrdantnm.")')
c     go to 100
      write(nf,'("/ * * * * block vi * * * *")')
      write(nf,'(" ajed=",i1)')ith
c     write(nf,'(" zned=1 igrped=1")')
c     write(nf,'(" rsfe= ",i3,"r 1.")')neg
c get reaction rates instead (comment above two lines).
      if(ith.eq.0)then
        write(nf,'(" zned=1 igrped=0 / resdnt=1")')
        write(nf,'(" edxs=")')
        do i=1,npn
          ct=" "
          ct(2:9)=edit(i)(1:8)
          write(ct(1:1),'(1h")')
          l=len_trim(ct)+1
          write(ct(l:l),'(1h")')
          write(nf,'(7x,a)')trim(ct)
        end do ! i
        write(nf,'(" edmats=")')
        if(nrrx.gt.0)then
          do i=1,nrrx
            write(nf,'(8x,"m",i6.6)')nm+i
          end do ! i
        else
          do i=1,nm
            write(nf,'(8x,"m",i6.6)')i
          end do ! i
        end if
        write(nf,'(" rsfnam=",1h","flux",1h")')
        write(nf,'(" rsfe=",i4,"r 1.")')neg
      else if(ith.eq.1)then
        call adjedit(nf)
      end if
c
      write(nf,'(" t")')
c
  100 close(nf)
c
      return
      end subroutine
c
      subroutine adjedit(nf)
c get <psi* q>. one response fct for each coarse mesh with a neut source.
      use COMS
      implicit none
      integer nf
      integer i,j,k,n1,ng,isi(max(nr,nz))
c
      write(nf,'(" zned=1 igrped=0")')
c     write(nf,'(" rsfnam=",1h","<psi* q>",1h")')
      write(nf,'(" rsfe=")')
      do i=1,nr
        do j=1,nz
          if(ismat(mat(i,j)).eq.0)cycle
          do ng=1,neg,3
            n1=min(neg,ng+2)
            write(nf,'(1x,1p3e16.9)')(nsrcf(k,mat(i,j)),k=ng,n1)
          end do ! ng
          write(nf,'("  ;")')
        end do ! j
      end do ! i
      write(nf,'(" rsfx=")')
      do i=1,nr
        do j=1,nz
          isi(1:nr)=0
          if(ismat(mat(i,j)).eq.0)cycle
          isi(i)=1
          write(nf,'(8(i6,"r ",i1,"."))')(iints(k),isi(k),k=1,nr)
          write(nf,'(7x,";")')
        end do ! j
      end do ! i
      if(calc2d.eq.1)then
        write(nf,'(" rsfy=")')
        do i=1,nr
          do j=1,nz
            isi(1:nz)=0
            if(ismat(mat(i,j)).eq.0)cycle
            isi(j)=1
            write(nf,'(8(i6,"r ",i1,"."))')(jints(k),isi(k),k=1,nz)
            write(nf,'(7x,";")')
          end do ! j
        end do ! i
      end if
      return
      end subroutine
