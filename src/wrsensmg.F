      subroutine wrsensmg
c write a sensmg input file from lnk3dnt data.
c called from rddantlnk.
      use COMS
      implicit none
      integer i,j,map(it,jt),n,nm1,ios,i_wgt
      character bfile*14,line*12
      real*8 rh(nmxsp)
      real(R8KIND), allocatable, dimension(:) :: rho1
c
c i_wgt=0/1 use atom densities/use weight densities
      i_wgt=1
c this is done in sensmg.py
c     if(ilnk3dnt.eq.0)then
c       write(*,'("error. wrsensmg capability is only for ",
c    1   "redoin/lnk3dnt.")')
c       write(iuo,'("error. wrsensmg capability is only for ",
c    1   "redoin/lnk3dnt.")')
c       call stoponerror
c     end if
c first get number of zones/materials and material map.
      map(1:it,1:jt)=0
      nm1=0
      do i=1,it
        do j=1,jt
          if(idclnk3(i,j,1).gt.0)then
            nm1=nm1+1
            map(i,j)=nm1
          end if
        end do ! j
      end do ! i
      allocate(rho1(nm1),stat=ios)
      bfile="tmp_sensmg_inp"
      write(*,'("writing ",a,"...")')trim(bfile)
      write(iuo,'("writing ",a,"...")')trim(bfile)
      open(iul1,file=bfile,status='unknown',position='rewind')
      write(iul1,'(a)')trim(id)
      if(calc2d.eq.-1)then
        line="slab"
      else if(calc2d.eq.0)then
        line="sphere"
      else if(calc2d.eq.1)then
        line="cyl"
      end if
      if(icalc.eq.1)then
        line=trim(line)//" keff"
      else if(icalc.eq.2)then
        line=trim(line)//" alpha"
      end if
      write(iul1,'(a)')trim(line)
      write(iul1,'(a)')trim(libname)
      write(iul1,'(i6," / no. of materials")')nm1
c mats
      rho1(1:nm1)=0.d0
      nm1=0
      do i=1,it
        do j=1,jt
          if(idclnk3(i,j,1).eq.0)cycle
          nm1=nm1+1
          write(iul1,'(i6)',advance='no')nm1
          do n=1,nmxsp
            if(idclnk3(i,j,n).eq.0)exit
            if(i_wgt.eq.0)then
              rh(n)=denlnk3(i,j,n)
            else if(i_wgt.eq.1)then
c if requested, convert atom density to mass density
              rh(n)=denlnk3(i,j,n)*atwt(idclnk3(i,j,n))/avg
            end if
            rho1(nm1)=rho1(nm1)+rh(n)
          end do ! n
          if(i_wgt.eq.1)then
            rho1(nm1)=-rho1(nm1)
          end if
          do n=1,nmxsp
            if(idclnk3(i,j,n).eq.0)exit
            rh(n)=rh(n)/rho1(nm1)
            write(iul1,'(i8,1pe16.8)',advance='no')
     1       nint(blk(1,idclnk3(i,j,n))),rh(n)
          end do ! n
          write(iul1,'(" /")',advance='yes')
        end do ! j
      end do ! i
c densities
      write(iul1,'(1p9999e16.8)')rho1(1:nm1)
c radii and heights
      if(calc2d.eq.-1.or.calc2d.eq.0)then
        write(iul1,'(i6," / no. of regions")')it
      else if(calc2d.eq.1)then
        write(iul1,'(2i6," / no. of cylinders, no. of heights")')it,jt
      end if
      if(calc2d.eq.-1)then
        write(iul1,'(1p9999e14.6)')r(0:it)
      else if(calc2d.eq.0.or.calc2d.eq.1)then
        write(iul1,'(1p9999e14.6)')r(1:it)
      end if
      if(calc2d.eq.1)then
        write(iul1,'(1p9999e14.6)')z(0:jt)
      end if
c material map
      do j=jt,1,-1
        write(iul1,'(9999i6)')(map(i,j),i=1,it)
      end do ! j
c reaction rates
      write(iul1,'(i6," / number of edit points")')nedpoints
      if(nedpoints.gt.0)then
        write(iul1,'(999999i6)')(edpoints(i),i=1,nedpoints)
      end if
      write(iul1,'(i6," / number of reaction-rate ratios")')nrrr
      do i=1,nrrr
        if(irrr(3,i).gt.0)then
          write(iul1,'(i8,i4,i8,i4," /")')irrr(1:4,i)
        else
          write(iul1,'(i8,i4,12x," /")')irrr(1:2,i)
        end if
      end do ! i
c     write(iul1,'(/,"to use this input, set drdzp=100.d0 ",
c    1 "and mdrdz=1 in rdmdl.")')
      close(iul1)
      deallocate(rho1,stat=ios)
      go to 999
c
  900 write(*,'("error opening ",a,".")')trim(bfile)
      write(iuo,'("error opening ",a,".")')trim(bfile)
      call stoponerror
  910 write(*,'("error writing ",a,".")')trim(bfile)
      write(iuo,'("error writing ",a,".")')trim(bfile)
      call stoponerror
c
  999 return
      end subroutine
