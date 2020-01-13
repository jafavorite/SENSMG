      subroutine wrdant2o(ith,nnn,niso0,gams)
c write for_inp and adj_inp, forward and adjoint partisn input files
c for neutron problem.  this routine assumes nm < 1,000,000.  
      use COMS
      implicit none
      integer niso0 ! number of isotopes in snxedt
      integer ith ! =0/1 forward/adjoint
      integer nnn
      real*8 gams(neg,0:nmom-1,it,jt)
      integer nf,i,ic,j,k,l,n,ng,nmq,ifido,insors,ios,mt
c ensure pfile has same length in wrdant2o, wrfixsrc, and wrdantpgga.
      character pfile*14,title*7,ct*26,k3c_7(7)*8
      integer k2,iitm,i_srcacc
      logical lk
      real*8 k3a_7(7)
          
c insors=0/1 source in partisn input/source in fixsrc file
c   (partisn's insors)
      insors=1
c
      ifido=-1   ! -1/0 for 4e18/6e12
c
      iitm=999   ! if iitm.gt.0 it will be used for forward and adjoint.
      i_srcacc=0 ! 0/1 for default/srcacc=no
c
      if(epsi.gt.1.001d-6)write(*,'("****",/,
     1 "**** warning.  epsi=",1pe9.2,/"****")')epsi
      if(ith.eq.0.and.(isrcacc_no.eq.1.or.isrcacc_no.eq.3))then
        i_srcacc=1
c        write(*,'(/,"comment. no source acceleration for forward ",
c     1   "calculation.")')
c        write(iuo,'(/,"comment. no source acceleration for forward ",
c     1   "calculation.")')
      else if(ith.eq.1.and.(isrcacc_no.eq.2.or.isrcacc_no.eq.3))then
        i_srcacc=1
c        write(*,'(/,"comment. no source acceleration for adjoint ",
c     1   "calculation.")')
c        write(iuo,'(/,"comment. no source acceleration for adjoint ",
c     1   "calculation.")')
      end if     
c
      nf=iup
      if(ith.eq.0)then
        write(pfile,'("fo2/",i5.5,"_inp")')nnn
        title="forward"
      else if(ith.eq.1)then
        write(pfile,'("ad2/",i5.5,"_inp")')nnn
        title="adjoint"
      end if
      open(nf,file=pfile,status='new',iostat=ios)
      if(ios.ne.0)then
        return ! it exists if use_existing = yes
      end if
      write(nf,'(5i6)')3,0,0,0,0
      write(nf,'(a)')trim(id)
      write(nf,'("input: ",a)')trim(ifile)
c icalc=-1 or we wouldn't be here.
      write(nf,'(1x,a7," input file ",i5.5,", lkg, second order")')
     1 title,nnn
c
      mt=nm+nrrx
      call block1(nf,mt,niso0)
      call block2(nf)
      call block3(nf,ith,4)      
c
      write(nf,'("/ * * * * block iv * * * *")')
c     if(ilnk3dnt.eq.0)then
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
!     else if(ilnk3dnt.eq.1)then
!       write(nf,'(3x,"matspec=",/,6x,"atdens")')
!       write(nf,'(3x,"matls=")')
!       do i=1,nel
!         ctmp1='"'//trim(adjustl(zaidfull(i)))//'"'
!         write(nf,'(6x,"m",i6.6,2x,a,"  1.;")')i,trim(ctmp1)
!       end do ! i 
!       write(nf,'(3x,"assign=",/,6x,"matls")')
!     else if(ilnk3dnt.eq.2)then
!c open redoin input file and advance to block 4.
!       open(iui,file=ifile,status='old',iostat=ios)
!  10   read(iui,'(a)',iostat=ios)line
!       if(line(1:15).ne."&block_4_arrays")go to 10
!       il=0
!       ctmp1=" "
!       do
!         read(iui,'(a)',iostat=ios)line
!         if(index(line,"/").ne.0)then
!           if(il.eq.3)write(nf,'(6x,";")')
!           exit
!         else if(index(line,"matspec=").ne.0)then
!           if(il.eq.3)write(nf,'(6x,";")')
!           il=1
!           write(nf,'(a)')trim(line) ! matspec
!         else if(index(line,"assign=").ne.0)then
!           if(il.eq.3)write(nf,'(6x,";")')
!           il=2
!           if(index(line,"matls").ne.0)then ! assign=matls
!             l1=index(line,'"')
!             l2=index(line(l1+1:l1+20),'"')+l1
!             line(l1:l1)=" "
!             line(l2:l2)=" "
!           end if
!           write(nf,'(a)')trim(line) ! assign
!         else if(index(line,"matls=").ne.0)then
!           il=3
!           write(nf,'(a)')trim(line) ! matls
!         else if(il.eq.1)then
!           write(nf,'(a)')trim(line)
!         else if(il.eq.2)then
!           if(index(line,"matls").ne.0)then ! assign=matls
!             l1=index(line,'"')
!             l2=index(line(l1+1:l1+20),'"')+l1
!             line(l1:l1)=" "
!             line(l2:l2)=" "
!             write(nf,'(a)')trim(line)
!           else
!             write(nf,'(a,a)')trim(line),";"
!           end if
!         else if(il.eq.3)then
!           l1=index(line,'"')
!           l2=index(line(l1+1:l1+20),'"')+l1
!           l3=l2-l1+1
!           if(ctmp1.eq." ")then ! 1st material
!             ctmp1(1:l3)=line(l1:l2)
!             write(nf,'(6x,a)')trim(ctmp1)
!             write(nf,'(14x,a)')
!    1         adjustl(trim(line(l2+1:len_trim(line))))
!           else if(ctmp1(1:l3).eq.line(l1:l2))then ! same material
!             write(nf,'(14x,a)')
!    1         adjustl(trim(line(l2+1:len_trim(line))))
!           else ! new material
!             ctmp1(1:l3)=line(l1:l2)
!             write(nf,'(6x,";")')
!             if(index(line,"matspec=").ne.0.or.
!    1         index(line,"assign=").ne.0)then
!               backspace(iui)
!               cycle
!             end if
!             write(nf,'(6x,a)')trim(ctmp1)
!             write(nf,'(14x,a)')
!    1         adjustl(trim(line(l2+1:len_trim(line))))
!           end if
!         end if
!       end do
!       close(iui)
!     end if ! ilnk3dnt
      write(nf,'(" t")')
c
      write(nf,'("/ * * * * block v * * * *")')
      if(icalc.eq.0.or.icalc.eq.3.or.icalc.eq.4)then
        write(nf,'(" ievt=",i1," isct=",i3)')0,isct
      else if(icalc.eq.1.or.icalc.eq.2)then
        write(nf,'(" ievt=",i1," isct=",i3)')icalc,isct
      end if
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
c print spherical harmonics (search for "Spherical harmonics",
c compare with last debug print in subroutine mcon.)
c     write(nf,'(" extras=17r0 1 f0")')
c debug moments
c     write(nf,'(" fluxp=2")')
      if(iitm.gt.0)write(nf,'(" iitm=",i3)')iitm
      write(nf,'(" iitl=0 oitm=9999")')
      if(i_srcacc.eq.1)write(nf,'(" srcacc=no")')
      if(nofxup.eq.1)write(nf,'(" nofxup=1")')
c get surface angular fluxes.
c iplane, jplane are used for partisn versions before 8_27.
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
        write(nf,'(" ith=0")')
c adjoint:
      else if(ith.eq.1)then
        write(nf,'(" ith=1")')
      end if
c      write(nf,'("/ line energies:",8f7.4)')(egpg(ng),ng=1,neg)
c source
      nmq=nmom-1
      if(insors.eq.1)then
        write(nf,'(" insors=1")')
        call wrfixsrc2o(pfile,nmq,ic,gams)
      else
        write(nf,'(" sourcf=")')
        do l=0,nmq
          if(nmq.gt.0)write(nf,'("/ order",i4)')l
          do ng=1,neg
            write(nf,'("/ group",i4)')ng
            do j=1,jt
c             write(nf,'(1x,1p5e14.6)')(gams(ng,l,i,j),i=1,it)
c             write(nf,'(7x,";")')
c this is clunky but partisn 8 can't read unattached semi-colons.
              do i=1,it,5
                k=min(it,i+4)
                if(k.ne.it)then
                  write(nf,'(1x,1p5e14.6)')(gams(ng,l,n,j),n=i,k)
                else
                  if(k-i.eq.0)then
                    write(nf,'(1x,1p1e14.6,";")')(gams(ng,l,n,j),n=i,k)
                  else if(k-i.eq.1)then
                    write(nf,'(1x,1p2e14.6,";")')(gams(ng,l,n,j),n=i,k)
                  else if(k-i.eq.2)then
                    write(nf,'(1x,1p3e14.6,";")')(gams(ng,l,n,j),n=i,k)
                  else if(k-i.eq.3)then
                    write(nf,'(1x,1p4e14.6,";")')(gams(ng,l,n,j),n=i,k)
                  else if(k-i.eq.4)then
                    write(nf,'(1x,1p5e14.6,";")')(gams(ng,l,n,j),n=i,k)
                  end if
                end if
              end do ! i
            end do ! j
          end do ! ng
        end do ! l
      end if ! insors
c
      if(nofxup.eq.1)write(nf,'(" nofxup=1 srcacc=no")')
c get surface angular fluxes.
c iplane, jplane are used for partisn versions before 8_27.
      if(calc2d.eq.1)then
        if(iplane.gt.-1)then ! 0 is caught in rdctrl
cccc      write(nf,'(" asleft=",i6," asrite=",i6)')iindex(iplane)+1,
cccc 1     iindex(iplane)+1
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
cccc      write(nf,'(" asbott=",i6," astop= ",i6)')jindex(jplane)+1,
cccc 1     jindex(jplane)+1
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
      write(nf,'(" t")')
c
c don't take time for these edits.
c     write(nf,'("/ * * * * block vi * * * *")')
c     write(nf,'(" zned=1 igrped=1 byvolp=0")')
c     if(ith.eq.0)then
c       write(nf,'(" ajed=0")')
c     else if(ith.eq.1)then
c       write(nf,'(" ajed=1")')
c     end if
c     write(nf,'(" rsfe=",i3,"r 1. ")')neg
c     write(nf,'(" t")')
c
      close(nf)
c 
      return
      end 
