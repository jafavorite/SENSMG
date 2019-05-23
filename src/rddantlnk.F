      subroutine rddantlnk_g
c read binary lnk3dnt file from partisn for geometry.
      use COMS
      implicit none
      character hname*8,huse(2)*8
      integer i,nwd
      parameter (nwd=27)
      integer*4 ivers,ia(nwd)
      character afile*20
      real*8 dum1
      integer idum1
      logical lex,lod
c        
c read from for directory.
c     afile='for/lnk3dnt'
c read from head directory.
      afile='lnk3dnt'
      write(*,'("reading ",a," for geometry...")')trim(afile)
      write(iuo,'("reading ",a," for geometry...")')trim(afile)
c
      inquire(file=afile,exist=lex,opened=lod)
      if(lod)then
        write(*,'("error. ",a," is already open.")')
        write(iuo,'("error. ",a," is already open.")')
        call stoponerror
      end if
      open(iul,file=afile,status='old',position='rewind',
#if ( defined INTEL )
     1 form='binary',
#elif ( defined GNU )
     1 access='stream', ! untested capability
#endif
     2 err=900)
      read(iul,err=910)hname,(huse(i),i=1,2),ivers
      read(iul,err=910)(ia(i),i=1,nwd)
      if(idbgw.eq.1)then
        write(41,'(a)')trim(afile)
        write(41,'(9i8)')(ia(i),i=1,nwd)
      endif
c     write(*,'("ivers=",i1)')ivers
c     write(*,'(9i8)')(ia(i),i=1,nwd)
      if(calc2d.eq.-1.and.ia(1).ne.1)go to 795
      if(calc2d.eq.0.and.ia(1).ne.3)go to 800
      if(calc2d.eq.1.and.ia(1).ne.7)go to 805
c at this point, nm is the number of materials.
c depending on the type of lnk3dnt setup, it may go to 1.
      if(ia(2).ne.nm)go to 810
      if(ia(5).ne.nr)go to 815
      if(ilkg.eq.1.and.ia(6).ne.nz)go to 820
      if(ilkg.eq.0.and.ia(6).ne.nz+2)go to 825
      if(ia(7).ne.1)go to 830
      if(ia(8).ne.it)go to 835
      if(ia(9).ne.jt)go to 840
      if(ia(10).ne.1)go to 845
c mixing instructions
      if(ia(24).gt.nel)go to 850
c regions per cell
      if(ia(25).ne.0)go to 855
c presence of block levels
      if(ia(27).ne.0)go to 860
c
c will be used to dimension idclnk3 and denlnk3 later.
      nmxsp=ia(24)
c coarse mesh boundaries, number of equally spaced fine meshes
      if(ia(1).eq.1.or.ia(1).eq.3)then
        read(iul,err=910)(r(i),i=0,ia(5)),(iints(i),i=1,ia(5))
        jints(1)=1
c       write(*,'(1p10e15.7)')(r(i),i=0,ia(5))
c       write(*,'(10i15)')(iints(i),i=1,ia(5))
      else if(ia(7).eq.1)then
        read(iul,err=910)(r(i),i=0,ia(5)),(z(i),i=0,ia(6)),
     1   (iints(i),i=1,ia(5)),(jints(i),i=1,ia(6))
c       write(*,'(1p10e15.7)')(r(i),i=0,ia(5))
c       write(*,'(10i15)')(iints(i),i=1,ia(5))
c       write(*,'(1p10e15.7)')(z(i),i=0,ia(6))
c       write(*,'(10i15)')(jints(i),i=1,ia(6))
      end if
c fine mesh interval widths if ivers.eq.5
      if(ivers.eq.5)then
        read(iul,err=910)idum1
        read(iul,err=910)idum1
        read(iul,err=910)idum1
        read(iul,err=910)(dum1,i=1,ia(8))
        read(iul,err=910)(dum1,i=1,ia(9))
        read(iul,err=910)(dum1,i=1,ia(10))
      end if
c     close(iul) ! close in rddantlnk_m
      write(*,'("done.")')
c
c write geometry to a file for debugging
c     if(idbgw.eq.0)go to 999
      go to 999
      write(*,'("writing geomtry for debugging...")')
c     call stoponerror
      go to 999 
c
  795 write(*,'("error. ia(1).ne.1 in ",a,"; ia(1)=",i8)')
     1 trim(afile),ia(1)
      go to 990
  800 write(*,'("error. ia(1).ne.3 in ",a,"; ia(1)=",i8)')
     1 trim(afile),ia(1)
      go to 990
  805 write(*,'("error. ia(1).ne.7 in ",a,"; ia(1)=",i8)')
     1 trim(afile),ia(1)
      go to 990
  810 write(*,'("error. ia(2).ne.nm in ",a,"; ia(2)=",i8)')
     1 trim(afile),ia(2)
      write(*,'("nm=",i8)')nm
      go to 990
  815 write(*,'("error. ia(5).ne.nr in ",a,"; ia(5)=",i8)')
     1 trim(afile),ia(5)
      write(*,'("nr=",i8)')nr
      go to 990
  820 write(*,'("error. ia(6).ne.nz in ",a,"; ia(6)=",i8)')
     1 trim(afile),ia(6)
      write(*,'("nz=",i8)')nz
      go to 990
  825 write(*,'("error. ia(6).ne.nz+2 in ",a,"; ia(6)=",i8)')
     1 trim(afile),ia(6)
      write(*,'("nz+2=",i8)')nz+2
      go to 990
  830 write(*,'("error. ia(7).ne.1 in ",a,"; ia(7)=",i8)')
     1 trim(afile),ia(7)
      go to 990
  835 write(*,'("error. ia(8).ne.it in ",a,"; ia(8)=",i8)')
     1 trim(afile),ia(8)
      write(*,'("it=",i8)')it
      go to 990
  840 write(*,'("error. ia(9).ne.jt in ",a,"; ia(9)=",i8)')
     1 trim(afile),ia(9)
      write(*,'("jt=",i8)')jt
      go to 990
  845 write(*,'("error. ia(10).ne.1 in ",a,"; ia(10)=",i8)')
     1 trim(afile),ia(10)
      go to 990
  850 write(*,'("error. ia(24).gt.nel in ",a,"; ia(24)=",i8)')
     1 trim(afile),ia(24)
      write(*,'("nel=",i8)')nel
      go to 990
  852 write(*,'("error. ia(24).ne.nmxsp in ",a,"; ia(24)=",i8)')
     1 trim(afile),ia(24)
      write(*,'("nmxsp=",i8)')nmxsp
      go to 990
  855 write(*,'("error. ia(25).ne.0 in ",a,"; ia(25)=",i8)')
     1 trim(afile),ia(25)
      go to 990
  860 write(*,'("error. ia(27).ne.0 in ",a,"; ia(27)=",i8)')
     1 trim(afile),ia(27)
      go to 990
  900 write(*,'("error opening ",a,".")')trim(afile)
      write(iuo,'("error opening ",a,".")')trim(afile)
      call stoponerror
  910 write(*,'("error reading ",a,".")')trim(afile)
      write(iuo,'("error reading ",a,".")')trim(afile)
      call stoponerror
c
  990 write(*,'("ensure lnk3dnt file corresponds to input ",a,
     1 ".")')trim(ifile)
      write(iuo,'("ensure lnk3dnt file corresponds to input ",a,
     1 ".")')trim(ifile)
      call stoponerror
  999 return
      end subroutine
c
      subroutine rddantlnk_m
c read binary lnk3dnt file from partisn for materials.
c it is already open and positioned.
      use COMS
      implicit none
      integer i,j,n
      character afile*20
      logical lex,lod
c        
c read from for directory.
c     afile='for/lnk3dnt'
c read from head directory.
      afile='lnk3dnt'
      write(*,'("reading ",a," for materials...")')trim(afile)
      write(iuo,'("reading ",a," for materials...")')trim(afile)
c
      inquire(file=afile,exist=lex,opened=lod)
      if(.not.lod)then
        write(*,'("error. ",a," is not already open.")')
        write(iuo,'("error. ",a," is not already open.")')
        call stoponerror
      end if
c
c mixing arrays
c macroscopic material number
      read(iul,err=910)(((idclnk3(i,j,n),i=1,it),j=1,jt),n=1,nmxsp)
c 1 means void so decrement by 1 to recover partisn input, where
c 0 means void.
      idclnk3(1:it,1:jt,1:nmxsp)=idclnk3(1:it,1:jt,1:nmxsp)-1
c material density
      read(iul,err=910)(((denlnk3(i,j,n),i=1,it),j=1,jt),n=1,nmxsp)
c     write(*,'("idclnk3")')
c     write(*,'(10i15)')(((idclnk3(i,j,n),i=1,it),j=1,jt),
c    1 n=1,nmxsp)
c     write(*,'("denlnk3")')
c     write(*,'(1p10e15.7)')(((denlnk3(i,j,n),i=1,it),j=1,jt),
c    1 n=1,nmxsp)
c     do i=1,it
c       do j=1,jt
c         do n=1,nmxsp
c           if(idclnk3(i,j,n).gt.0)then
c             write(*,'(3i6,i8,1pe14.6)')i,j,n,
c    1         nint(blk(1,idclnk3(i,j,n))),denlnk3(i,j,n)
c           end if
c         end do ! n
c       end do ! j
c     end do ! i
c cell temperatures
c if the file stays open, we will not know ia(23)
c     if(ivers.eq.4.and.ia(23).ge.1)then
c       read(iul,err=910)((dum1,i=1,ia(8)),j=1,ia(9))
c     end if
      close(iul)
      write(*,'("done.")')
c
c write geometry to a file for debugging
c     if(idbgw.eq.0)go to 999
      go to 999
      write(*,'("writing geomtry for debugging...")')
c     call stoponerror
      go to 999 
c
  900 write(*,'("error opening ",a,".")')trim(afile)
      write(iuo,'("error opening ",a,".")')trim(afile)
      call stoponerror
  910 write(*,'("error reading ",a,".")')trim(afile)
      write(iuo,'("error reading ",a,".")')trim(afile)
      call stoponerror
c
  999 return
      end subroutine
