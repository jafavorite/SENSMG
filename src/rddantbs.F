      subroutine rddantbs(iw,ib,ix,fluxbs)
c read binary bs files from partisn.
c iw=1/2/3/4 read forward/adjoint/generalized adjoint/
c            forward for quadrature only.
c ib=index of generalized adjoint
c ix=1/2 r plane/z plane
c
      use COMS
      implicit none
      integer iw,ib,ix
      real*8 fluxbs(neg,ndir,0:max(nr,nz),max(it,jt))
      real*8 d1,d2,d3,xi(ndir),ybnd(jt+1),zbnd(2)
      character hname*8,huse(2)*8
      integer nwd
      parameter (nwd=7)
      integer*4 ivers,ia(nwd)
      integer iob,ipl,nf,ng,m,n1,n2,n3,i,mm,m1,m2,j,j2,mord(ndir/2)
      character afile*17,pre*8
c        
c note that the raflxm file does not contain the number of
c directions; ndir is needed separately.
c
      if(iw.eq.1.or.iw.eq.4)then
c read forward.
        pre="for/"
      else if(iw.eq.2)then
c read adjoint.
        pre="adj/"
      else if(iw.eq.3)then
c read generalized adjoint from the last iteration. 
c TODO update
        write(pre,'("a",i2.2,"/a",i2.2,"_")')ib,ib
      end if
c
      fluxbs(1:neg,1:ndir,0:max(nr,nz),1:max(it,jt))=0.d0
      do nf=1,2 ! first do left/bot, then rite/top
        iob=0
        if(ix.eq.1)then ! radial surface
          ipl=iplane
          j2=jt
          if(nf.eq.1)then
            afile=trim(pre)//'bsleft-000000'
            m1=ndir/4+1
            m2=3*ndir/4+1
          else if(nf.eq.2)then
            m1=1
            m2=ndir/2+1
            if(ipl.eq.nr)iob=1
            afile=trim(pre)//'bsrite-000000'
          end if
        else if(ix.eq.2)then ! axial surface
          ipl=jplane
          j2=it
          if(nf.eq.1)then
            m1=ndir/2+1
            m2=3*ndir/4+1
            if(ipl.eq.0)iob=2
            afile=trim(pre)//'bsbot-000000'
          else if(nf.eq.2)then
            m1=1
            m2=ndir/4+1
            if(ipl.eq.nz)iob=3
            afile=trim(pre)//'bstop-000000'
          end if
        end if
        mord(1)=m1
        do m=2,ndir/2
          if(m.eq.ndir/4+1)then
            mord(m)=m2
          else
            mord(m)=mord(m-1)+1
          end if
        end do ! m
c       do m=1,ndir/2
c         write(*,'(2i5)')m,mord(m)
c       end do ! m
c
c put in adjoint source for leakage problem.
c this assumes leakage only; for rxn rate, modify.
c this assumes total leakage.
        if(iob.gt.0)then
          if(iw.eq.2.and.(icalc.eq.0.or.icalc.eq.3.or.icalc.eq.4))then
            fluxbs(1:neg,mord(1:ndir/2),ipl,1:j2)=1.d0
          end if
          cycle
        end if
c
        if(iw.eq.4)then
          write(*,'("reading ",a," for quadrature...")')trim(afile)
          write(iuo,'("reading ",a," for quadrature...")')trim(afile)
        else
          write(*,'("reading ",a,"...")')trim(afile)
          write(iuo,'("reading ",a,"...")')trim(afile)
        end if
c
        open(iua,file=afile,status='old',position='rewind',
     1   form='unformatted',access='sequential',
     2   err=900)
        read(iua,err=910)hname,(huse(i),i=1,2),ivers
        read(iua,err=910)(ia(i),i=1,nwd)
        if(idbgw.eq.1)then
          write(41,'(a)')trim(afile)
          write(41,'(10i8)')(ia(i),i=1,nwd)
        endif
        if((calc2d.eq.0.or.calc2d.eq.-1).and.ia(1).ne.1)go to 795
        if(calc2d.eq.1.and.ia(1).ne.2)go to 797
        if(ia(2).ne.neg)go to 800
        if(ia(3).ne.it)go to 805
        if(ia(4).ne.jt)go to 810
        if(ia(5).ne.1)go to 815
        read(iua,err=910)d1,d2,d3
        if(idbgw.eq.1)then
          write(41,'(1p10e14.6)')d1,d2,d3
        endif
c
        read(iua,err=910)(ybnd(j),j=1,j2+1)
        read(iua,err=910)(zbnd(j),j=1,2) ! 2 zeroes
c
        if(iw.eq.1)then
          n1=1
          n2=ia(2)
          n3=1
        else if(iw.eq.2.or.iw.eq.3)then
c note that the adjoint fluxes are printed in reverse order and read in
c reverse order, so they are in forward order for the rest of the code.
          n1=ia(2)
          n2=1
          n3=-1
        else if(iw.eq.4)then
c only do enough to get quadrature.
          n1=1
          n2=1
          n3=1
        end if
c
        if(ia(1).eq.1)then
c TODO check this for one-d. normally would not be needed.
          do ng=n1,n2,n3
            do m=1,ndir/2
              read(iua,err=910)fluxbs(ng,m,ipl,1)
            end do ! m
          end do ! ng
        elseif(ia(1).eq.2)then
          do ng=n1,n2,n3
            read(iua,err=910)mm
            if(iw.eq.4)then
              read(iua,err=910)(wgt(m),m=1,mm)
              read(iua,err=910)(dir(m),m=1,mm)
              read(iua,err=910)(eta(m),m=1,mm)
              xi(1:mm)=sqrt(1.d0-dir(1:mm)**2-eta(1:mm)**2)
              do m=mm+1,ndir
                wgt(m)=wgt(m-mm)
                dir(m)=dir(m-mm) ! not correct but fills in the array (actually maybe it is correct)
                eta(m)=eta(m-mm)
                xi(m)=xi(m-mm)
              end do ! m
            else
              read(iua,err=910)(d1,m=1,mm)
              read(iua,err=910)(d1,m=1,mm)
              read(iua,err=910)(d1,m=1,mm)
            end if
            do j=1,j2
              read(iua,err=910)(fluxbs(ng,mord(m),ipl,j),m=1,ndir/2)
            end do ! j
          end do ! ng
        end if
        close(iua)
        write(*,'("done.")')
c       write(*,'("rddantbs.")')
c       do m=1,mm
c         write(*,'(i5,1p4e14.6)')m,dir(m),eta(m),xi(m),wgt(m)
c       end do ! m
        if(iw.eq.4)exit
      end do ! nf
c     if(calc2d.eq.0)then
c       do m=1,ndir
c         if(calc2d.eq.0)then
c           write(41,'(i4,1p20e14.6)')m,wgt(m),dir(m)
c         elseif(calc2d.eq.1)then
c           write(41,'(i4,1p20e14.6)')m,wgt(m),dir(m),eta(m),xi(m)
c         end if
c       end do ! m
c     end if
c
c write fluxes to a file for debugging
c redo all this for this file
c     if(idbgw.eq.0)go to 999
      go to 999
      if(calc2d.eq.1)go to 999
      if(iw.eq.2)go to 999
      write(*,'("writing one-d for. ang. fluxes for debugging...")')
      n1=0
      n2=0
      do ng=1,ia(2)
        do m=1,ndir 
          write(41,'("angular flux for group=",i4," mu=",i4," mu,wgt=",
     1     1p3e14.6)')ng,m,dir(m),wgt(m),dir(m)*wgt(m)
c         do i=1,ia(3)+1
          do i=ia(3)+1,ia(3)+1
            write(41,'(i6,1p2e14.6)')i,afreg(ng,m,i,1),afadj(ng,m,i,1)
            if(afreg(ng,m,i,1).lt.0..and.n1.eq.0)then
              write(*,'("negative regular flux",3i6)')ng,m,i
              n1=1 
            end if
            if(afadj(ng,m,i,1).lt.0..and.n2.eq.0)then 
              write(*,'("negative adjoint flux",3i6)')ng,m,i
              n2=1
            end if 
          end do ! i
        end do ! m 
      end do ! ng 
      go to 999 
c
  795 write(*,'("error. ia(1).ne.1 in ",a,"; ia(1)=",i8)')
     1 trim(afile),ia(1)
      call stoponerror
  797 write(*,'("error. ia(1).ne.2 in ",a,"; ia(1)=",i8)')
     1 trim(afile),ia(1)
      call stoponerror
  800 write(*,'("error. ia(2).ne.neg in ",a,"; ia(2)=",i8)')
     1 trim(afile),ia(2)
      write(*,'("neg=",i8)')neg
      call stoponerror
  805 write(*,'("error. ia(3).ne.it in ",a,"; ia(3)=",i8)')
     1 trim(afile),ia(3)
      write(*,'("it=",i8)')it
      call stoponerror
  810 write(*,'("error. ia(4).ne.jt in ",a,"; ia(4)=",i8)')
     1 trim(afile),ia(4)
      write(*,'("jt=",i8)')jt
      call stoponerror
  815 write(*,'("error. ia(5).ne.1 in ",a,"; ia(5)=",i8)')
     1 trim(afile),ia(5)
      call stoponerror
  900 write(*,'("error opening ",a,".")')trim(afile)
      write(iuo,'("error opening ",a,".")')trim(afile)
      call stoponerror
  910 write(*,'("error reading ",a,".")')trim(afile)
      write(iuo,'("error reading ",a,".")')trim(afile)
      call stoponerror
c
  999 return
      end subroutine
