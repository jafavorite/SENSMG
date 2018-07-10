      subroutine rdctrl(ncbc)
c read file "control" for the control integer ictrl.
      use COMS
      implicit none
      character ncbc*120
      integer ierr,ios
c
c----------------------------------------------------
c initialize values
      imisc=1 ! 0/1 sources4c/misc
      ialphan=1 ! 0/1 ignore/use (alpha,n) sources
      ictrl=0
      iter=-1
      isn=0
      isct=0
      neg=30
      iver=0
      ichinorm=1
      isrcacc_no=0
      iangflux=0
      aflxfrm=0
      idbgw=0
      ilkg=1
      nofxup=1
      nag=100
      epsi=1.d-6
      epsig=1.d-5
      partisn="7"
c
c read control file
      open(iuc,file='control',status='old',iostat=ios)
c
      ierr=0
c
      read(iuc,*)ictrl,iter
c no use for -1 yet, but here it is.
      if(ictrl.eq.-1) then
         close(iuc)
         return
      else if(ictrl.eq.1.or.ictrl.eq.2.or.ictrl.eq.3.or.ictrl.eq.4.or.
     1 ictrl.eq.5)then
         go to 100
      else
         write(*,'("error. ictrl.ne.-1, 1, 2, 3, 4, or 5; it is",
     1    i5)')ictrl
         call stoponerror
      end if
c     
c read other things here
  100 continue
      read(iuc,*)isn,isct,neg,ichinorm,isrcacc_no,fissdata,iangflux,
     1 aflxfrm,iver,idbgw
      read(iuc,*)nm, nel, nr, nz, nrrr
      read(iuc,'(a120)')ncbc
      read(iuc,*)imisc,ialphan,nag
      read(iuc,*)nofxup
      read(iuc,*)epsi,epsig
      read(iuc,*)ifile
      read(iuc,*)partisn
      read(iuc,*)iplane,jplane
c use defaults for now
c     read(iuc,*)ilkg
      close(iuc)
c
      if(ilkg.lt.0.or.ilkg.gt.1)then
        write(*,'("error. ilkg must be 0 or 1; it is",i5)')ilkg
        ierr=1
      end if
      if(nofxup.lt.0.or.nofxup.gt.1)then
        write(*,'("error. nofxup must be 0 or 1; it is",i6)')nofxup
        ierr=1
      end if
      if(ichinorm.lt.0.or.ichinorm.gt.2)then
        write(*,'("error. ichinorm must be 0, 1, or 2; it is",i5)')
     1   ichinorm
        ierr=1
      end if
      if(iangflux.lt.0.or.iangflux.gt.1)then
        write(*,'("error. iangflux must be 0 or 1; it is",i5)')iangflux
        ierr=1
      end if
      if(aflxfrm.lt.0.or.aflxfrm.gt.1)then
        write(*,'("error. aflxfrm must be 0 or 1; it is",i5)')aflxfrm
        ierr=1
      end if
      if(aflxfrm.eq.1.and.partisn(1:1).eq."5")then
        write(*,'("warning. aflxfrm is not used in partisn ",
     1   "version 5; resetting to 0.")')aflxfrm
        write(iuo,'("warning. aflxfrm is not used in partisn ",
     1   "version 5; resetting to 0.")')aflxfrm
        aflxfrm=0
      end if
      if(iver.lt.0.or.iver.gt.2)then
        write(*,'("error. iver must be 0, 1, or 2; it is",i5)')iver
        ierr=1
      end if
      if(nag.lt.1.or.nag.gt.9999)then
        write(*,'("error. nag must be between 1 and 9999; it is",i5)')
     1   nag
        ierr=1
      end if
      if(iplane.ne.-1.and.iplane.lt.1.and.iplane.gt.nr)then
        write(*,'("error. iplane must be -1 or between 1 and nr:",i5)')
     1   iplane
        ierr=1
      end if
      if(jplane.lt.-1.and.jplane.gt.nz)then
        write(*,'("error. jplane must be -1 or between 0 and nz:",i5)')
     1   jplane
        ierr=1
      end if
      if(ierr.eq.1)call stoponerror
      return
      end subroutine
