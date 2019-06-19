      subroutine wrdantxs
c write xs1_inp, partisn input file that will only give microscopic
c cross sections.  this scheme gives barns.
c this routine assumes niso < 1,000,000.
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
      else if(libname(1:5).eq.'scale')then
        n0=5
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
      write(nf,'(5i6)')4,0,0,0,0
      write(nf,'(a)')trim(id)
      write(nf,'("input: ",a)')trim(ifile)
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
        write(nf,'(6x,0pf6.1)')real(i,8)
      end do ! i
      write(nf,'(" xints= ",i4,"r 4")')mt
      write(nf,'(" zones=")')
      do i=1,mt
        write(nf,'(4x,i6)')i
      end do ! i
      write(nf,'(" t")')
c
      call block3(nf,0,1)
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
        write(nf,'(8x,"m",i6.6,2x,a,"  1.;")')i,trim(ct)
      end do ! i
      write(nf,'(" assign=")')
      do i=1,niso
        write(nf,'(8x,"zn",i6.6,"   m",i6.6,"  1.;")')i,i
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
