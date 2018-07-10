      subroutine rddantas(iw,ib,ix,fluxas)
c read binary as files from partisn.
c iw=1/2/3/4 read forward/adjoint/generalized adjoint/
c            forward for quadrature only.
c ib=index of generalized adjoint
c ix=1/2 r plane/z plane
c
      use COMS
      implicit none
      integer iw,ib,ix
      real*8 fluxas(neg,ndir,0:max(nr,nz),max(it,jt))
      real*8 d1,d2,xi(ndir),ybnd(jt+1),zbnd(2)
      character hname*8,huse(2)*8
      integer*4 ivers,ia(10)
      integer ipl,ng,m,n,n1,n2,n3,i,i1,mm,m1,j,j2,noct,mord(ndir/2)
      character afile*17,pre*4
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
c TODO update if this has to be used; hopefully will not.
        write(afile,'("a",i2.2,"_aaflxm")')ib
      end if
c
      fluxas(1:neg,1:ndir,0:max(nr,nz),1:max(it,jt))=0.d0
      if(ix.eq.1)then ! radial surface
        j2=jt
        afile=pre//'asfluxx'
        i1=0 ! for axial surfaces, index starts at 0
      else if(ix.eq.2)then ! axial surface
        j2=it
        afile=pre//'asfluxy'
        i1=1 ! for radial surfaces, index starts at 0
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
     1 form='unformatted',access='sequential',
     2 err=900)
      read(iua,err=910)hname,(huse(i),i=1,2),ivers
      read(iua,err=910)(ia(i),i=1,7)
      if(idbgw.eq.1)then
        write(41,'(a)')trim(afile)
        write(41,'(10i8)')(ia(i),i=1,7)
      endif
      if((calc2d.eq.0.or.calc2d.eq.-1).and.ia(1).ne.1)go to 795
      if(calc2d.eq.1.and.ia(1).ne.2)go to 797
      if(ia(2).ne.neg)go to 800
      if(ia(3).ne.it)go to 805
      if(ia(4).ne.jt)go to 810
      if(ia(5).ne.1)go to 815
      if(ix.eq.1.and.ia(7).ne.nr)go to 820
      if(ix.eq.2.and.ia(7).ne.nz+1)go to 825
c     read(iua,err=910)d1,d2,d3
c     if(idbgw.eq.1)then
c       write(41,'(1p10e14.6)')d1,d2,d3
c     endif
c     noct=2**ia(1)
      if(ia(1).eq.1)then
        noct=2
      else if(ia(1).eq.2)then
        noct=4
      end if
c     write(*,'("noct=",i6)')noct
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
            read(iua,err=910)fluxas(ng,m,ipl,1)
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
c         if(ng.eq.n1)write(*,'("mm,ndir,noct=",3i6)')mm,ndir,noct
          do i=1-i1,ia(7)-i1
            read(iua,err=910)ipl ! plane index
            read(iua,err=910)d2  ! plane location
c           write(*,'("ipl,mm,ndir,d2=",3i6,1pe14.6)')ipl,mm,ndir,d2
            do n=1,noct ! for cylinders, noct=4
              if(n.eq.1)then
                m1=1
              else if(n.eq.2)then
                m1=ndir/4+1
              else if(n.eq.3)then
                m1=ndir/2+1
              else if(n.eq.4)then
                m1=3*ndir/4+1
              end if
              mord(1)=m1
              do m=2,mm ! mm=ndir/4
                mord(m)=mord(m-1)+1
              end do ! m
c             write(*,'(3i6)')ng,i,n
              do j=1,j2
                read(iua,err=910)(fluxas(ng,mord(m),i,j),m=1,mm)
              end do ! j
            end do ! n
          end do ! i
        end do ! ng
c put in adjoint source for leakage problem.
c this assumes leakage only; for rxn rate, modify.
c this assumes total leakage.
        if(iw.eq.2.and.(icalc.eq.0.or.icalc.eq.3.or.icalc.eq.4))then
          if(ix.eq.1)then
            m=1
            m1=ndir/4
            fluxas(1:neg,m:m1,nr,1:j2)=1.d0
            m=ndir/2+1
            m1=3*ndir/4
            fluxas(1:neg,m:m1,nr,1:j2)=1.d0
          else if (ix.eq.2)then
            m=ndir/2+1
            m1=ndir
            fluxas(1:neg,m:m1,0,1:j2)=1.d0
            m=1
            m1=ndir/2
            fluxas(1:neg,m:m1,nz,1:j2)=1.d0
          end if
        end if
      end if
      close(iua)
      write(*,'("done.")')
c     write(*,'("rddantas.")')
c     do m=1,mm
c       write(*,'(i5,1p4e14.6)')m,dir(m),eta(m),xi(m),wgt(m)
c     end do ! m
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
      call stoponerror
  805 write(*,'("error. ia(3).ne.it in ",a,"; ia(3)=",i8)')
     1 trim(afile),ia(3)
      call stoponerror
  810 write(*,'("error. ia(4).ne.jt in ",a,"; ia(4)=",i8)')
     1 trim(afile),ia(4)
      call stoponerror
  815 write(*,'("error. ia(5).ne.1 in ",a,"; ia(5)=",i8)')
     1 trim(afile),ia(5)
      call stoponerror
  820 write(*,'("error. ia(7).ne.nr in ",a,"; ia(7)=",i8)')
     1 trim(afile),ia(7)
      call stoponerror
  825 write(*,'("error. ia(7).ne.nz+1 in ",a,"; ia(7)=",i8)')
     1 trim(afile),ia(7)
      call stoponerror
  900 write(*,'("error opening ",a,".")')trim(afile)
      call stoponerror
  910 write(*,'("error reading ",a,".")')trim(afile)
      call stoponerror
c
  999 return
      end subroutine
