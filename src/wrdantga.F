      subroutine wrdantga(nrx,iconv,ig,niso0,gams)
c write a??_inp, generalized adjoint partisn input files
c for neutron problem.  this routine assumes nm < 10000.  
c this routine handles feyny or sm2 first-order sensitivity inputs.
      use COMS
      implicit none
      integer niso0 ! number of isotopes in snxedt
      integer iconv ! flag for convergence
      integer ig ! flag for what the adjoint source is
      real*8 gams(neg,0:nmom-1,it,jt)
c for ig=1/2, nrx is the isotope and reaction index
c for ig=3, nrx tells which calculation to do
      integer nrx
      integer i,j,k,k2,l,n,iitm,mt,nf,ng,i_srcacc,insors,nmq
      character pfile*11,ct*26
      iitm=999   ! if iitm.gt.0 it will be used for forward and adjoint.
      i_srcacc=1 ! 0/1 for default/srcacc=no
c insors=0/1 source in partisn input/source in fixsrc file
c   (partisn's insors)
      insors=1
c
c these warnings are only needed once.
      if(nrx.eq.1.and.iter.eq.0)then
        if(epsi.gt.1.001d-6)then
          write(*,'("****",/,"**** warning.  epsi=",1pe9.2,/"****")')
     1     epsi
          write(iuo,'("****",/,"**** warning.  epsi=",1pe9.2,/"****")')
     1     epsi
        end if
        write(*,'(/,"comment. no source acceleration for ",
     1   "generalized adjoint calculation.")')
        write(iuo,'(/,"comment. no source acceleration for ",
     1   "generalized adjoint calculation.")')
      end if
c
      nf=iup
      if(ig.eq.1.or.ig.eq.2)then
        write(pfile,'("a",i2.2,"/a",i2.2,"_inp")')nrx,nrx
      else if(ig.eq.3)then
        if(nrx.eq.1.or.nrx.eq.2)then
          write(pfile,'("smf/",i2.2,"_inp")')nrx
        else if(nrx.eq.3)then
          write(pfile,'("sma/",i2.2,"_inp")')nrx
        end if
      end if
      open(nf,file=pfile,status='unknown')
c
      if(iconv.eq.1)then
        write(nf,'("generalized adjoint converged")')
        go to 100
      end if
c
      if(calc2d.eq.0.or.calc2d.eq.-1)then
        if(ig.eq.1)then
          write(nf,'("     4     0     0")')
        else if(ig.eq.2)then
          write(nf,'("     5     0     0")')
        else if(ig.eq.3)then
          write(nf,'("     4     0     0")')
        end if
      else if(calc2d.eq.1)then
        if(ig.eq.1)then
          write(nf,'("     4     0     0")')
        else if(ig.eq.2)then
          write(nf,'("     5     0     0")')
        else if(ig.eq.3)then
          write(nf,'("     4     0     0")')
        end if
      end if
      write(nf,'(a)')trim(ifile)
      write(nf,'(a)')trim(id)
      if(ig.eq.1)then
        write(nf,'(" generalized input file, z=",i6," r=",i2)')
     1   irrx(1:2,nrx)
      else if(ig.eq.2)then
        write(nf,'(" generalized input file, ",/,
     1   " numer z=",i6," r=",i2,"; denom z=",i6," r=",i2)')
     2   irrx(1:2,irri(1,nrx)),irrx(1:2,irri(2,nrx))
      else if(ig.eq.3)then
        if(nrx.eq.1)then
          write(nf,'(" subcrit. mult. sens., 1st forward")')
        else if(nrx.eq.2)then
          write(nf,'(" subcrit. mult. sens., 2nd forward")')
        else if(nrx.eq.3)then
          write(nf,'(" subcrit. mult. sens., adjoint")')
        end if
      end if
      if(icalc.eq.0.or.icalc.eq.3.or.icalc.eq.4)then
        write(nf,'(" (A*-F*)G*=S*")')
      else if(icalc.eq.1.or.icalc.eq.2)then
        if(iter.eq.0)then
          write(nf,'(" it.",i3,", (A*)G0*=S*")')iter
        else
          write(nf,'(" it.",i3,", A*[G(n)*]=F*/k[G(n-1)]*")')iter
        end if
      end if
c
c need to use same nm as forward
      mt=nm+nrrx
      call block1(nf,mt,niso0)
      call block2(nf)
c
      write(nf,'("/ * * * * block iii * * * *")')
c use macrxs from forward run.
      write(nf,'(" lib=macrxs")')
      write(nf,'(" lng=",i3)')lng
      write(nf,'(" t")')
c
      write(nf,'("/ * * * * block iv * * * *")')
      write(nf,'(" matspec=wtfrac")')
      write(nf,'(" matls=")')
      k2=0
      do i=1,nm
        write(nf,'(8x,"m",i4.4)')i
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
      write(nf,'(" assign=")')
      do i=1,nm
        write(nf,'(8x,"zone",i4.4,"   m",i4.4,f14.8,";")')i,i,rho(i)
      end do ! i
      write(nf,'(" t")')
c
      write(nf,'("/ * * * * block v * * * *")')
      write(nf,'(" ievt=0 isct=",i1)')isct
      if(icalc.eq.1.or.icalc.eq.2)then
        write(nf,'(" nosigf=1 / inhibit fission multiplication")')
      end if
      if(icalc.eq.2)then
        write(nf,'(" anorm=",1pe15.7)')alpha
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
      write(nf,'(" rmflux=1")')
c originally this routine did not have iangflux for calc2d=1; why?
      if(iangflux.eq.0.or.cellsol.eq.1)then
        write(nf,'(" raflux=0")')
      else if(iangflux.eq.1.and.(cellsol.eq.-1.or.cellsol.eq.0))then
        write(nf,'(" raflux=1")')
      end if
      if(idbgw.eq.1)write(nf,'(" fluxp=1")')
c     write(nf,'(" fluxp=2")')
      if(iitm.gt.0)write(nf,'(" iitm=",i3)')iitm
      write(nf,'(" iitl=0 oitm=9999")')
      if(i_srcacc.eq.1)write(nf,'(" srcacc=no")')
c
      if((ig.eq.1.or.ig.eq.2).or.(ig.eq.3.and.nrx.eq.3))then
        write(nf,'(" ith=1 xsectp=0")')
      else if(ig.eq.3.and.(nrx.eq.1.or.nrx.eq.2))then
        write(nf,'(" ith=0 xsectp=0")')
      end if
c ig.eq.2 has negative sources, requires nofxup=1
      if(ig.eq.2.or.nofxup.eq.1)then
        write(nf,'(" nofxup=1")')
      end if
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
          write(nf,'(" asbott=",i6," astop=",i6)')jindex(jplane)+1,
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
c source
      nmq=0 ! this can be changed in the future if needed.
      if(insors.eq.1)then
        write(nf,'(" insors=1")')
        call wrfixsrc(ig,nrx,nmq,gams)
      else
        write(nf,'(" sourcf=")')
        do l=0,nmq
          if(nmq.gt.0)write(nf,'("/ order",i4)')l
          do ng=1,lng
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
      write(nf,'(" t")')
c
c get <psi* q> for fixed-source problems.
      if(icalc.eq.0)then
        write(nf,'("/ * * * * block vi * * * *")')
        write(nf,'(" ajed=1")')
c       write(nf,'(" zned=1 igrped=1")')
c       write(nf,'(" rsfe=",i3,"r 1.")')neg
        call adjedit(nf)
        write(nf,'(" t")')
      end if
c
  100 close(nf)
      return
      end subroutine
