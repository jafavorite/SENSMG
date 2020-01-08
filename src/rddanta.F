      subroutine rddanta(iw,ib)
c read binary raflxm and aaflxm files from partisn.
c in one-d, afreg=psi on edges, afadj=psi* on edges.
c in two-d, afreg=psi in centers, afadj=psi* in centers.
c the one-d format only works for partisn after may 2004.
c iw=1/2/3/4/5/6 read forward/adjoint/generalized adjoint/
c                forward for quadrature only/
c                feyny or sm2 forward/feyny or sm2 adjoint.
c ib=index of generalized adjoint
c
      use COMS
      implicit none
      integer iw,ib
      real*8 d1,d2,d3,d4,xi(ndir)
      character hname*8,huse(2)*8
      integer nwd
      parameter (nwd=8)
      integer*4 ivers,ia(nwd)
      integer ng,m,n1,n2,n3,i,mm,m2,n,m1,j
      character afile*20
c        
c note that the raflxm file does not contain the number of
c directions; ndir is needed separately.
c
      if(iw.eq.1.or.iw.eq.4)then
c read raflxm.
        afile='for/raflxm'
      else if(iw.eq.2)then
c read aaflxm.
        afile='adj/aaflxm'
      else if(iw.eq.3)then
c read generalized adjoint from the last iteration. 
c afgad is used here because it is declared in common;
c it will be afgadn.
        write(afile,'("a",i2.2,"/aaflxm")')ib
      else if(iw.eq.5)then
c read raflxm.
        write(afile,'("smf/",i2.2,"_raflxm")')ib
      else if(iw.eq.6)then
c read aaflxm.
        write(afile,'("sma/",i2.2,"_aaflxm")')ib
      end if
      if(iw.eq.4)then
        write(*,'("reading ",a," for quadrature...")')trim(afile)
        write(iuo,'("reading ",a," for quadrature...")')trim(afile)
      else
        write(*,'("reading ",a,"...")')trim(afile)
        write(iuo,'("reading ",a,"...")')trim(afile)
      end if
c
      open(iua,file=afile,status='old',position='rewind',
#if ( defined INTEL )
     1 form='binary',
#elif ( defined GNU )
     1 access='stream', ! untested capability
#endif
     2 err=900)
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
      if((calc2d.eq.0.or.calc2d.eq.-1).and.ia(6).ne.2)go to 820
      if(calc2d.eq.1.and.ia(6).ne.4)go to 822
c
      if(iw.eq.1.or.iw.eq.5)then
        n1=1
        n2=ia(2)
        n3=1
      else if(iw.eq.2.or.iw.eq.3.or.iw.eq.6)then
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
        do ng=n1,n2,n3
          do m=1,ndir
            read(iua,err=910)dir(m),d1,d2,wgt(m),d3,d4
            if(iw.eq.1.or.iw.eq.5)then
              read(iua,err=910)(afreg(ng,m,i,1),i=1,ia(3)+1)
            else if(iw.eq.2.or.iw.eq.6)then
              read(iua,err=910)(afadj(ng,m,i,1),i=1,ia(3)+1)
            else if(iw.eq.3)then
              read(iua,err=910)(afgad(ng,m,i,1),i=1,ia(3)+1)
            end if
          end do ! m
        end do ! ng
      elseif(ia(1).eq.2)then
c       open(42,file="scalar",status="unknown")
        do ng=n1,n2,n3
          read(iua,err=910)mm ! should be ndir/4
          read(iua,err=910)(wgt(m),m=1,mm)
          read(iua,err=910)(dir(m),m=1,mm)
          read(iua,err=910)(eta(m),m=1,mm)
          read(iua,err=910)(xi(m),m=1,mm)
          m2=0
          do n=1,ia(6)
            m1=m2+1
            m2=mm*n
            if(iw.eq.1.or.iw.eq.5)then
              read(iua,err=910)(((afreg(ng,m,i,j),i=1,ia(3)),j=1,ia(4)),
     1         m=m1,m2)
            else if(iw.eq.2.or.iw.eq.6)then
              read(iua,err=910)(((afadj(ng,m,i,j),i=1,ia(3)),j=1,ia(4)),
     1         m=m1,m2)
            else if(iw.eq.3)then
              read(iua,err=910)(((afgad(ng,m,i,j),i=1,ia(3)),j=1,ia(4)),
     1         m=m1,m2)
            end if
          end do ! n
        end do ! ng
        do m=mm+1,ndir
          wgt(m)=wgt(m-mm)
          dir(m)=dir(m-mm) ! not correct but fills in the array (actually maybe it is correct)
          eta(m)=eta(m-mm)
          xi(m)=xi(m-mm)
        end do ! m
c       do ng=1,npl
c         do m=1,ndir
c           do j=1,jt
cc            do i=1,it
c               fmom(ng,0,it,j)=fmom(ng,0,it,j)+wgt(m)*afreg(ng,m,it,j)
cc            end do ! i
c           end do ! j
c         end do ! m
c       end do ! ng
c       do j=1,jt
c         write(42,'(1pe14.6,3x,e14.6,3x,e14.6,3x,e14.6,3x,e14.6)')
c    1     zfm(j),(fmom(ng,0,it,j),ng=1,npl)
c       end do ! j
c       stop
      end if
      close(iua)
      write(*,'("done.")')
c     if(calc2d.eq.0.or.iw.eq.4)then
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
c     if(idbgw.eq.0)go to 999
      go to 999
      if(calc2d.eq.1)go to 999
c     if(iw.ne.2)go to 999
      write(*,'("writing one-d for. ang. fluxes for debugging...")')
      write(*,'("file ",a)')trim(afile)
      write(41,'("file ",a)')trim(afile)
      n1=0
      n2=0
      do ng=1,ia(2)
        do m=1,ndir 
          write(41,'("angular flux for group=",i4," mu=",i4," mu,wgt=",
     1     1p3e14.6)')ng,m,dir(m),wgt(m),dir(m)*wgt(m)
          do i=1,ia(3)+1
c         do i=ia(3)+1,ia(3)+1
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
  820 write(*,'("error. ia(6).ne.2 in ",a,"; ia(6)=",i8)')
     1 trim(afile),ia(6)
      call stoponerror
  822 write(*,'("error. ia(6).ne.4 in ",a,"; ia(6)=",i8)')
     1 trim(afile),ia(6)
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
c
      subroutine rddants
c read binary sncons.
c
      use COMS
      implicit none
      character hname*8,huse(2)*8
      integer i,m,nwd
      parameter (nwd=4)
      integer*4 ivers,ia(nwd)
      character afile*8
c        
      afile='sncons'
      write(*,'("reading ",a,"...")')trim(afile)
c
      open(iua,file=afile,status='old',position='rewind',
#if ( defined INTEL )
     1 form='binary',
#elif ( defined GNU )
     1 access='stream', ! untested capability
#endif
     2 err=900)
      read(iua,err=910)hname,(huse(i),i=1,2),ivers
      read(iua,err=910)(ia(i),i=1,nwd)
      write(*,'(11i8)')(ia(i),i=1,nwd)
      if(idbgw.eq.1)then
        write(41,'(a)')trim(afile)
        write(41,'(11i8)')(ia(i),i=1,nwd)
      endif
      if((calc2d.eq.0.or.calc2d.eq.-1).and.ia(1).ne.1)go to 795
      if(calc2d.eq.1.and.ia(1).ne.2)go to 797
      if(ia(2).ne.ndir)go to 800
c
      if(ia(1).eq.1)then
        read(iua,err=910)(wgt(m),m=1,ia(2)),(dir(m),m=1,ia(2))
      elseif(ia(1).eq.2)then
c TODO check this for cylinders
        read(iua,err=910)(wgt(m),m=1,ia(2)),(dir(m),m=1,ia(2)),
     1   (eta(m),m=1,ia(2))
      end if
      close(iua)
      write(*,'("done.")')
c
c     do m=1,ndir
c       if(calc2d.eq.0.or.calc2d.eq.-1)then
c         write(42,'(i4,1p20e14.6)')m,wgt(m),dir(m)
c       elseif(calc2d.eq.1)then
c         write(42,'(i4,1p20e14.6)')m,wgt(m),dir(m),eta(m)
c       end if
c     end do ! m
      go to 999
c
  795 write(*,'("error. ia(1).ne.1 in ",a,"; ia(1)=",i8)')
     1 trim(afile),ia(1)
      call stoponerror
  797 write(*,'("error. ia(1).ne.2 in ",a,"; ia(1)=",i8)')
     1 trim(afile),ia(1)
      call stoponerror
  800 write(*,'("error. ia(2).ne.ndir in ",a,"; ia(2)=",i8)')
     1 trim(afile),ia(2)
      write(*,'("ndir=",i8)')ndir
      call stoponerror
  900 write(*,'("error opening ",a,".")')trim(afile)
      call stoponerror
  910 write(*,'("error reading ",a,".")')trim(afile)
      call stoponerror
c
  999 return
      end subroutine
