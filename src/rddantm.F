      subroutine rddantm(iw,ib)
c read binary rmflux and amflux files from partisn.
c the moments are in mesh centers.
c iw=1/2/3/5/6/7 read forward/adjoint/generalized adjoint/
c                feyny or sm2 forward/feyny or sm2 adjoint/
c                forward stream
c ib=index of generalized adjoint
c
      use COMS
      implicit none
      integer iw,ib
      character hname*8,huse(2)*8
      integer ng,n1,n2,n3,i,j,k,l,nwd
      parameter (nwd=11)
      integer*4 ivers,ia(nwd)
      character afile*20
c        
      if(iw.eq.1)then
c read rmflux.
        afile='for/rmflux'
      else if(iw.eq.2)then
c read amflux.
        afile='adj/amflux'
      else if(iw.eq.3)then
c read generalized adjoint from the last iteration. 
c gmom is used here because it is declared in common;
c it will be gmomn.
        write(afile,'("a",i2.2,"/amflux")')ib
      else if(iw.eq.5)then
c read raflxm.
        write(afile,'("smf/",i2.2,"_rmflux")')ib
      else if(iw.eq.6)then
c read aaflxm.
        write(afile,'("sma/",i2.2,"_amflux")')ib
      else if(iw.eq.7)then
c read stream.
        afile='for/stream'
      end if
      write(*,'("reading ",a,"...")')trim(afile)
      write(iuo,'("reading ",a,"...")')trim(afile)
c
      open(ium,file=afile,status='old',position='rewind',
#if ( defined INTEL )
     1 form='binary',
#elif ( defined GNU )
     1 access='stream', ! untested capability
#endif
     2 err=900)
      read(ium,err=910)hname,(huse(i),i=1,2),ivers
      read(ium,err=910)(ia(i),i=1,nwd)
      if(idbgw.eq.1)then
        write(41,'(a)')trim(afile)
        write(41,'(11i8)')(ia(i),i=1,nwd)
      endif
      if((calc2d.eq.0.or.calc2d.eq.-1).and.ia(1).ne.1)go to 795
      if(calc2d.eq.1.and.ia(1).ne.2)go to 797
      if(ia(2).ne.neg)go to 800
      if(ia(3).ne.it)go to 805
      if(ia(4).ne.jt)go to 810
      if(ia(5).ne.1)go to 815
      if(ia(6).ne.nmom)go to 845
      if(ia(10).ne.0)go to 850
      if(ia(11).ne.0)go to 855
c
      if(iw.eq.1.or.iw.eq.5.or.iw.eq.7)then
        n1=1
        n2=ia(2)
        n3=1
      else if(iw.eq.2.or.iw.eq.3.or.iw.eq.6)then
c note that the adjoint fluxes are printed in reverse order and read in
c reverse order, so they are in forward order for the rest of the code.
        n1=ia(2)
        n2=1
        n3=-1
      end if
c
      do ng=n1,n2,n3
        if(iw.eq.1.or.iw.eq.5.or.iw.eq.7)then
          read(ium,err=910)(((fmom(ng,0,i,j),i=1,ia(3)),j=1,ia(4)),
     1     k=1,ia(5))
          if(ia(6).gt.1)then
            if(cellsol.ne.1)then
              read(ium,err=910)((((fmom(ng,l,i,j),i=1,ia(3)),j=1,ia(4)),
     1         k=1,ia(5)),l=1,ia(6)-1)
            else if(cellsol.eq.1)then
              read(ium,err=910)((((fmom(ng,l,i,j),l=1,ia(6)-1),
     1         i=1,ia(3)),j=1,ia(4)),k=1,ia(5))
            end if
          end if
        elseif(iw.eq.2.or.iw.eq.6)then
          read(ium,err=910)(((amom(ng,0,i,j),i=1,ia(3)),j=1,ia(4)),
     1     k=1,ia(5))
          if(ia(6).gt.1)then
            if(cellsol.ne.1)then
              read(ium,err=910)((((amom(ng,l,i,j),i=1,ia(3)),j=1,ia(4)),
     1         k=1,ia(5)),l=1,ia(6)-1)
            else if(cellsol.eq.1)then
              read(ium,err=910)((((amom(ng,l,i,j),l=1,ia(6)-1),
     1         i=1,ia(3)),j=1,ia(4)),k=1,ia(5))
            end if
          end if
        elseif(iw.eq.3)then
          read(ium,err=910)(((gmom(ng,0,i,j),i=1,ia(3)),j=1,ia(4)),
     1     k=1,ia(5))
          if(ia(6).gt.1)then
            if(cellsol.ne.1)then
              read(ium,err=910)((((gmom(ng,l,i,j),i=1,ia(3)),j=1,ia(4)),
     1         k=1,ia(5)),l=1,ia(6)-1)
            else if(cellsol.eq.1)then
              read(ium,err=910)((((gmom(ng,l,i,j),l=1,ia(6)-1),
     1         i=1,ia(3)),j=1,ia(4)),k=1,ia(5))
            end if
          end if
        end if
      end do ! ng
      close(ium)
      write(*,'("done.")')
c
c write moments to a file for debugging
c     if(idbgw.eq.0)go to 999
      go to 999
      write(*,'("writing moments for debugging...")')
      do ng=n1,n2,n3
        do l=0,0 !ia(6)-1
          do i=1,ia(3)
            do j=1,ia(4)
              if(iw.eq.7)then
                write(32,'("ng: ",i3," l: ",i3," i: ",i5," j: ",i5,
     1           " fmom: ",1pe13.6)')ng,l,i,j,fmom(ng,l,i,j)
              else if(iw.eq.99)then
                write(32,'("ng: ",i3," l: ",i3," i: ",i5," j: ",i5,
     1           " amom: ",1pe13.6)')ng,l,i,j,amom(ng,l,i,j)
              end if
            end do ! j
          end do ! i
        end do ! l
      end do ! ng
c     call stoponerror
c     if(iw.eq.1)then
c     do ng=n1,n2,n3
c       do l=0,ia(6)-1
c         do i=1,ia(3)
c           do j=1,ia(4)
c               if(abs(fmom(ng,l,i,j)/1.102048d-07-1.d0).lt.1.d-5)
c               if(abs(fmom(ng,l,i,j)/1.102273E-07-1.d0).lt.1.d-5)
c    1            write(32,'("ng: ",i3," l: ",i3," i: ",i5," j: ",i5,
c    2            " fmom: ",1pe13.6)')ng,l,i,j,fmom(ng,l,i,j)
c           end do ! j
c         end do ! i
c       end do ! l
c     end do ! ng
c     end if
c     call stoponerror
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
  845 write(*,'("error. ia(6).ne.nmom in ",a,"; ia(6)=",i8)')
     1 trim(afile),ia(6)
      write(*,'("nmom=",i8)')nmom
      call stoponerror
  850 write(*,'("error. ia(10).ne.0 in ",a,"; ia(10)=",i8)')
     1 trim(afile),ia(10)
      call stoponerror
  855 write(*,'("error. ia(11).ne.0 in ",a,"; ia(11)=",i8)')
     1 trim(afile),ia(11)
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
