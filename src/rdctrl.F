      subroutine rdctrl(ic,iwrite,iplotg,iwrxsecs,iwrsensmg)
c read file "control" for the control integer ictrl.
      use COMS
      implicit none
c ic=0/1 read all except ncb/read ncb
      integer ic,iwrite,iplotg,iwrxsecs,iwrsensmg
      integer i,ierr,ios
      character line*1
c
c----------------------------------------------------
c initialize values
      if(ic.eq.0)then
        ictrl=0
        iter=-1
        iver=0
        ichinorm=1
        ilnk3dnt=0
        i_2nd_order=1
        i_xsecs=1
        isrcacc_no=0
        iangflux=0
        aflxfrm=0
        idbgw=0
        ilkg=1
        nofxup=1
        epsi=1.d-6
        epsig=1.d-5
        iplotg=0
        iwrxsecs=0
        iwrsensmg=0
      end if
c
c read control file
      open(iuc,file='control',status='old',iostat=ios)
c
      ierr=0
c
      if(ic.eq.0)then
        read(iuc,*)ictrl,iwrite,iter
c no use for -1 yet, but here it is.
        if(ictrl.eq.-1) then
           close(iuc)
           return
        else if(ictrl.eq.1.or.ictrl.eq.2.or.ictrl.eq.3.or.ictrl.eq.4.or.
     1   ictrl.eq.5)then
           go to 100
        else
           write(*,'("error. ictrl.ne.-1, 1, 2, 3, 4, or 5; it is",
     1      i5)')ictrl
           call stoponerror
        end if
c     
c read other things here
  100   continue
        read(iuc,*)isn,isct,neg,ichinorm,isrcacc_no,fissdata,chieff,
     1   fissneut,iangflux,aflxfrm,iver,idbgw
        read(iuc,*)nm,nel,nr,nz,nrrr,nedpoints
        read(iuc,*)ilnk3dnt,iplotg,iwrxsecs,iwrsensmg
        read(iuc,'(a)')line ! ncb; will return to read this from rdmdl.
        read(iuc,*)imisc,ialphan,nag
        read(iuc,*)nofxup,itrcor,i_2nd_order,i_xsecs
        read(iuc,*)epsi,epsig
        read(iuc,*)ifile
        read(iuc,*)partisnver
        read(iuc,*)iplane,jplane
c use defaults for now
c       read(iuc,*)ilkg
        close(iuc)
      else if(ic.eq.1)then
        read(iuc,*)line
        read(iuc,*)line
        read(iuc,*)line
        read(iuc,*)line
        read(iuc,*)(ncb(i),i=1,nm)
        close(iuc)
        return
      end if
c
      if(ilkg.lt.0.or.ilkg.gt.1)then
        write(*,'("error. ilkg must be 0 or 1; it is",i5)')ilkg
        ierr=1
      end if
      if(nofxup.lt.0.or.nofxup.gt.1)then
        write(*,'("error. nofxup must be 0 or 1; it is",i6)')nofxup
        ierr=1
      end if
c this is checked in sensmg.py
c     if(itrcor.lt.0.or.itrcor.gt.3)then
c       write(*,'("error. itrcor must be 0, 1, 2, or 3; it is",i6)')
c    1   itrcor
c       ierr=1
c     end if
c this is checked in sensmg.py
c     if(i_2nd_order.lt.0.or.i_2nd_order.gt.1)then
c       write(*,'("error. i_2nd_order must be 0 or 1; it is",i6)')
c    1   i_2nd_order
c       ierr=1
c     end if
c this is checked in sensmg.py
c     if(i_xsecs.lt.0.or.i_xsecs.gt.1)then
c       write(*,'("error. i_xsecs must be 0 or 1; it is",i6)')i_xsecs
c       ierr=1
c     end if
c this is checked in sensmg.py
c     if(ichinorm.lt.0.or.ichinorm.gt.2)then
c       write(*,'("error. ichinorm must be 0, 1, or 2; it is",i5)')
c    1   ichinorm
c       ierr=1
c     end if
      if(iangflux.lt.0.or.iangflux.gt.1)then
        write(*,'("error. iangflux must be 0 or 1; it is",i5)')iangflux
        ierr=1
      end if
c this is checked in sensmg.py
c     if(aflxfrm.lt.0.or.aflxfrm.gt.1)then
c       write(*,'("error. aflxfrm must be 0 or 1; it is",i5)')aflxfrm
c       ierr=1
c     end if
      if(aflxfrm.eq.1.and.partisnver(1:1).eq."5")then
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
