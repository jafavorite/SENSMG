      subroutine wrdantxs
c write xs1_inp, partisn input file that will only give microscopic
c cross sections.  this scheme gives barns.
c this routine assumes niso < 10000.  
c one coarse mesh and zone for each isotope in the problem,
c including edit isotopes.
      use COMS
      implicit none
      integer i,l,mt,nf,n0
      character afile*11,ct*26
      nf=iup ! output file number
c
      if(libname(1:6).eq.'kynea3')then
        n0=151
      else if(libname(1:3).eq.'b96')then
        n0=225
      else if(libname(1:7).eq.'special')then
c assume mendf71x
        n0=17570
      else
        n0=0
      end if
c recall niso=nel+nrrx
      mt=niso
      write(afile,'("xs1/xs1_inp")')
      open(nf,file=afile,status='unknown')
      write(nf,'("     4     0     0")')
      write(nf,'(a)')trim(ifile)
      write(nf,'(a)')trim(id)
      write(nf,'(" cross section input file; used only for ",/,
     1 " obtaining cross sections from macrxs.")')
c
      write(nf,'("/ * * * * block i * * * *")')
      write(nf,'(" igeom=sphere isn=",i3," ngroup=",i4)')2,neg
c     write(nf,'("/maxscm=20000000 maxlcm=9000000")')
      write(nf,'(" niso=",i6," mt=",i4," nzone=",i4)')n0,mt,mt 
      write(nf,'(" im=",i4," it=",i5)')mt,4*mt
c for xsectp to take effect (to print the cross sections), comment this
      write(nf,'(" nosolv=1")')
      write(nf,'(" t")')
c
      write(nf,'("/ * * * * block ii * * * *")')
      write(nf,'(" xmesh=  0.")')
      do i=1,mt
        write(nf,'(6x,0pf6.1)')real(i)
      end do ! i
      write(nf,'(" xints= ",i4,"r 4")')mt
      write(nf,'(" zones=")')
      do i=1,mt
        write(nf,'(4x,i4)')i
      end do ! i
      write(nf,'(" t")')
c
      write(nf,'("/ * * * * block iii * * * *")')
      if(libname(1:6).ne.'kynea3'.and.libname(1:3).ne.'b96'.and.
     1 libname(1:7).ne.'special')then
        write(nf,'(" lib=ndilib fissneut=1 / total nu-bar")')
c JAF -- testing
c       write(*,'("WARNING! USING PROMPT NU!")')
c       write(nf,'(" lib=ndilib fissneut=0 / prompt nu-bar")')
c JAF -- testing
        if(fissdata.eq.1)then
          if(partisn(1:1).eq."5")then
            write(nf,'(" nochimat=0 / chi matrix")')
          else
            write(nf,'(" fissdata=1 / chi matrix")')
          end if
        else if(fissdata.eq.2)then
          if(partisn(1:1).eq."5")then
            write(nf,'(" nochimat=1 / chi vector")')
          else
            write(nf,'(" fissdata=2 / chi vector")')
          end if
        else if(fissdata.eq.0)then
          if(partisn(1:1).ne."5")then
            write(nf,'(" fissdata=0 / fission transfer matrix")')
          end if
        end if
c       if(libname(1:5).ne."mt71x")then
          write(nf,'(" libname=",a)')trim(libname)
c       end if
      else if(libname(1:6).eq."kynea3")then
        write(nf,'(" lib=bxslib fissneut=1 / total nu-bar")')
        write(nf,'(" libname=",a)')trim(libname)
        if(partisn(1:1).eq."5")then
          write(nf,'(" nochimat=1 / chi vector")')
        else
          write(nf,'(" fissdata=2 / chi vector")')
        end if
      else if(libname(1:7).eq."special")then
        write(nf,'("/ this file does not actually get run")')
        write(nf,'(" lib=macrxs")')
      end if
c     write(nf,'(" grpstr=default ggrpstr=default")')
      write(nf,'(" lng=",i3)')lng
      write(nf,'(" t")')
c
      write(nf,'("/ * * * * block iv * * * *")')
      write(nf,'(" matspec=atdens")')
      write(nf,'(" matls=")')
      do i=1,niso
        ct=" "
        ct(2:25)=zaidfull(i)(1:24)
        l=scan(ct,abc123)-1
        write(ct(l:l),'(1h")')
        l=len_trim(ct)+1
        write(ct(l:l),'(1h")')
        write(nf,'(8x,"m",i4.4,2x,a,"  1.;")')i,trim(ct)
      end do ! i
      write(nf,'(" assign=")')
      do i=1,niso
        write(nf,'(8x,"zone",i4.4,"   m",i4.4,"  1.;")')i,i
      end do ! i
      write(nf,'(" t")')
c
      write(nf,'("/ * * * * block v * * * *")')
      write(nf,'(" ievt=0 isct=",i1)')isct
      write(nf,'(" iitm=0 kitm=0 iitl=0 oitm=2")')
      write(nf,'(" npeg=2")')
      write(nf,'(" ith=0 xsectp=2 epsi=10.")')
      write(nf,'(" sirite= ",i4,"r 1.")')neg
      write(nf,'(" t")')
c
      close(nf)
c
      return
      end subroutine
