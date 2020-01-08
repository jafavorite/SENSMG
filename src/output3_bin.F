      subroutine output3_bin(sintr,sintz)
c write binary sensrx file.
      use COMS
      implicit none
      real*8 sintr(1:neg,0:nr,1:nz),sintz(1:neg,1:nr,0:nz)
      integer i,j,ng
      character afile*20
c
      write(afile,'("sensrx")')
      write(*,'("writing ",a,"...")')trim(afile)
      write(iuo,'("writing ",a,"...")')trim(afile)
c
c the file is still open.
      write(iub)(((sintr(ng,i,j),ng=1,neg),i=0,nr),j=1,nz)
      if(calc2d.eq.1)then
        write(iub)(((sintz(ng,i,j),ng=1,neg),i=1,nr),j=0,nz)
      end if
      close(iub)
c
      return
      end subroutine
c
      subroutine rdoutput3_bin(sintr,sintz)
c read binary sensrx file.
      use COMS
      implicit none
      real*8 sintr(1:neg,0:nr,1:nz),sintz(1:neg,1:nr,0:nz)
      integer i,j,ng
      character afile*20
c
      write(afile,'("sensrx")')
      write(*,'("reading ",a,"...")')trim(afile)
      write(iuo,'("reading ",a,"...")')trim(afile)
c
c the file is still open.
      read(iub)(((sintr(ng,i,j),ng=1,neg),i=0,nr),j=1,nz)
      if(calc2d.eq.1)then
        read(iub)(((sintz(ng,i,j),ng=1,neg),i=1,nr),j=0,nz)
      else
        sintz(1:neg,1:nr,0:nz)=0.d0
      end if
      close(iub)
      write(*,'("done.")')
      go to 999
c
  910 write(*,'("error reading ",a,".")')trim(afile)
      write(iuo,'("error reading ",a,".")')trim(afile)
      call stoponerror
c
  999 return
      end subroutine
