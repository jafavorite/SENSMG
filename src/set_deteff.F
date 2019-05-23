      subroutine set_deteff(iwrite,neg,iug2,iuo,sens_data,detname,
     1 deteff,pi)
c read detector response function file.
c to add a new file, do three things:
c 1. increase nf.
c 2. add the file name in the data statement.
c 3. add logic to set the name correctly.
      implicit none
      integer iwrite,neg,iug2,iuo
      character sens_data*220,detname*8
      real*8 deteff(neg),pi
      integer l,ng,n
      real*8 x1,x2,y1,y2,d,geomeff
      character gfile*120
      integer nf
c nf is the number of detector efficiency files.
c increase to add new files.
      parameter (nf=8)
      character det_eff_files(nf)*12
c DEBUG_ALEX use Alex's efficiency
c     real*8 eff_alex(30)
c DEBUG_ALEX
c add new group-structure files in this list.
      data det_eff_files/
     & "npod_001",
     & "npod_002",
     & "npod_030",
     & "npod_044",
     & "npod_070",
     & "npod_079",
     & "npod_130",
     & "npod_618" 
     & /
c DEBUG_ALEX use Alex's efficiency
c     data eff_alex/2.071658040154163941d-03, 2.224700284734928805d-03,
c    1 2.335476799370791740d-03, 2.572063658984820138d-03,
c    1 3.037012310960913661d-03, 3.866870809450288164d-03,
c    1 5.172366429192214177d-03, 6.575874221297642319d-03,
c    1 7.703164521729598951d-03, 8.896692533892800944d-03,
c    1 9.953557489324679763d-03, 1.056121829505005105d-02,
c    1 1.263086628472844149d-02, 1.376168402911592467d-02,
c    1 1.450645037975756280d-02, 1.496835321426106917d-02,
c    1 1.550853156274624098d-02, 1.602818233321747290d-02,
c    1 1.656034460103910505d-02, 1.716626602999230325d-02,
c    1 1.773836341732824839d-02, 1.826678654303908259d-02,
c    1 1.871372676899101667d-02, 1.992197742623653955d-02,
c    1 2.045785944386214333d-02, 2.081965371217949987d-02,
c    1 2.022575534802075609d-02, 1.605680317317980527d-02,
c    1 1.286481107734700627d-03, 3.082674141738126793d-07/
c DEBUG_ALEX
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
c H. Gotoha and H. Yagi, "Solid angle subtended by a rectangular slit,"
c Nuclear Instruments and Methods, vol. 96, No. 3, 485-486 (1971);
c https://doi.org/10.1016/0029-554X(71)90624-0
      geomeff=atan(x2*y2/sqrt(x2**2+y2**2+d**2)/d)
     1 -atan(x1*y2/sqrt(x1**2+y2**2+d**2)/d)
     2 -atan(x2*y1/sqrt(x2**2+y1**2+d**2)/d)
     3 +atan(x1*y1/sqrt(x1**2+y1**2+d**2)/d)
      geomeff=geomeff/(4.d0*pi)
      deteff(1:neg)=deteff(1:neg)*geomeff
      if(iwrite.eq.1)then
        write(iuo,'("geometry efficiency for ",a)')trim(detname)
        write(iuo,'(2x,"x1=",1pe16.8,/,2x,"x2=",1pe16.8)')x1,x2
        write(iuo,'(2x,"y1=",1pe16.8,/,2x,"y2=",1pe16.8)')y1,y2
        write(iuo,'(2x," d=",1pe16.8)')d
        write(iuo,'(2x,"geomeff=",1pe16.8,/)')geomeff
      end if
c DEBUG_ALEX use Alex's efficiency
c     write(*,'("warning. using alex`s efficiency.")')
c     write(iuo,'("warning. using alex`s efficiency.")')
c     deteff(1:30)=eff_alex(1:30)
c DEBUG_ALEX
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
  910 write(*,'("set_deteff. error opening ",a,".")')trim(gfile)
      write(iuo,'("set_deteff. error opening ",a,".")')trim(gfile)
      call stoponerror
c
  999 return
      end subroutine
