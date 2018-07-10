      subroutine set_deteff(neg,iug2,iuo,sens_data,detname,deteff,pi,
     1 lwrite)
c read detector response function file.
c to add a new file, do three things:
c 1. increase nf.
c 2. add the file name in the data statement.
c 3. add logic to set the name correctly.
      implicit none
      integer neg,iug2,iuo
      character sens_data*220,detname*8
      real*8 deteff(neg),pi
      logical lwrite
      integer l,ng,n
      real*8 x1,x2,y1,y2,d,geomeff
      character gfile*120
      integer nf
c nf is the number of detector efficiency files.
c increase to add new files.
      parameter (nf=6)
      character det_eff_files(nf)*12
c add new group-structure files in this list.
      data det_eff_files/
     & "npod_030",
     & "npod_044",
     & "npod_070",
     & "npod_079",
     & "npod_130",
     & "npod_618" 
     & /
c
c if detector file is "none," set all efficiencies to 1 to use
c total leakge and return.
      if(detname(1:4).eq."none")then
        deteff(1:neg)=1.d0
        go to 999 ! return
c add logic here for new detector files.
      else if(detname(1:4).eq."npod")then
        gfile=trim(detname)//"_"
c location of object center wrt to npod dimensions
        x2=43.18d0/2.d0
        x1=-x2
        y2=(42.243375d0+0.079375d0)/2.d0
        y1=-y2
c distance to object center
        d=50.d0
      else
        write(*,'("error. unknown detector name in set_deteff: ",a)')
     1   trim(detname)
        write(iuo,'("error. unknown detector name in set_deteff: ",a)')
     1   trim(detname)
        call stoponerror
      end if
      l=len_trim(gfile)
      write(gfile(l+1:l+3),'(i3.3)')neg
c ==============================================
      gfile=adjustl(gfile)
c
c check against the list of available files.
      do n=1,nf
        if(trim(adjustl(det_eff_files(n))).eq.
     1   trim(gfile))go to 10
      end do ! n
      go to 900 ! error
c
   10 gfile=trim(sens_data)//gfile
      open(iug2,file=gfile,status="old",position="rewind",err=910)
      do ng=1,neg
        read(iug2,*)l,deteff(ng)
      end do ! ng
      close(iug2)
c
c geometry efficiency.
c formula for rectangular detector, Eq. (3) in
c H.Gotoha and H.Yagi, "Solid angle subtended by a rectangular slit,"
c Nuclear Instruments and Methods, vol. 96, No. 3, 485-486 (1971);
c https://doi.org/10.1016/0029-554X(71)90624-0
      geomeff=atan(x2*y2/sqrt(x2**2+y2**2+d**2)/d)
     1 -atan(x1*y2/sqrt(x1**2+y2**2+d**2)/d)
     2 -atan(x2*y1/sqrt(x2**2+y1**2+d**2)/d)
     3 +atan(x1*y1/sqrt(x1**2+y1**2+d**2)/d)
      geomeff=geomeff/(4.d0*pi)
      deteff(1:neg)=deteff(1:neg)*geomeff
      if(lwrite)then
        write(iuo,'(/,"geometry efficiency for ",a)')trim(detname)
        write(iuo,'(2x,"x1=",1pe16.8,/,2x,"x2=",1pe16.8)')x1,x2
        write(iuo,'(2x,"y1=",1pe16.8,/,2x,"y2=",1pe16.8)')y1,y2
        write(iuo,'(2x," d=",1pe16.8)')d
        write(iuo,'(2x," geomeff=",1pe16.8,/)')geomeff
      end if
      go to 999
c
  900 write(*,'("error. ",a," does not exist. ",
     1 "available efficiencies are:")')trim(gfile)
      write(iuo,'("error. ",a," does not exist. ",
     1 "available efficiencies are:")')trim(gfile)
      write(*,'(5x,a,5x,a)')
     1 (adjustr(det_eff_files(n)),n=1,nf)
      write(iuo,'(5x,a,5x,a)')
     1 (adjustr(det_eff_files(n)),n=1,nf)
      write(*,'("to add a detector efficiency file, see source ",
     1 "subroutine set_deteff.F.")')
      write(iuo,'("to add a detector efficiency file, see source ",
     1 "subroutine set_deteff.F.")')
      call stoponerror
  910 write(*,'("error opening ",a,".")')trim(gfile)
      write(iuo,'("error opening ",a,".")')trim(gfile)
      call stoponerror
c
  999 return
      end subroutine
