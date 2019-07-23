      subroutine set_ebins(neg,iug2,iuo,sens_data,libname,ebins)
c read group-structure files for fixed-source problems.
c to add a new file, do three things:
c 1. increase nf.
c 2. add the file name in the data statement.
c 3. add logic to set the name correctly.
      implicit none
      integer neg,iug2,iuo
      character sens_data*220,libname*8
      real*8 ebins(3,neg)
      integer l,ng,n
      character gfile*120
      integer nf
c nf is the number of group-structure files.
c increase to add new files.
      parameter (nf=29)
      character group_structure_files(nf)*12
c add new group-structure files in this list.
      data group_structure_files/
     & "618_lanl",
     & "301_lanl",
     & "250_lanl",
     & "141_lanl",
     & "133_lanl",
     & "091_lanl",
     & "075_lanl",
     & "070_lanl",
     & "049_lanl",
     & "047_lanl",
     & "030_lanl",
     & "021_lanl",
     & "016_lanl",
     & "014_lanl",
     & "013_lanl",
     & "012_lanl",
     & "010_lanl",
     & "009_lanl",
     & "008_lanl",
     & "007_lanl",
     & "006_lanl",
     & "005_lanl",
     & "004_lanl",
     & "003_lanl",
     & "002_lanl",
     & "001_lanl",
     & "130_lanl",
     & "079_kynea3",
     & "044_scale"
c    & "016_lanl2",
c    & "012_lanl2",
c    & "002_lanl2",
c    & "087_llnl",
c    & "017_llnl",
c    & "016_llnl",
c    & "105_awe",
c    & "016_hr",
     & /
c group structures 087_llnl, 017_llnl, 016_llnl, 105_awe, and 016_hr
c exist, but because there is no logic to use these in partisn, they
c are not allowed in sensmg. same with 016_lanl2, 012_lanl2, and
c 002_lanl2.
c
      write(gfile,'(i3.3)')neg
c add logic here for new group-structure files.
      if(libname(1:6).eq."kynea3")then
        gfile=trim(gfile)//"_kynea3"
      else if(libname(1:5).eq."scale")then
        gfile=trim(gfile)//"_scale"
      else
        gfile=trim(gfile)//"_lanl"
      end if
c ==============================================
      gfile=adjustl(gfile)
c
c check against the list of available files.
      do n=1,nf
        if(trim(adjustl(group_structure_files(n))).eq.
     1   trim(gfile))go to 10
      end do ! n
      go to 900
c
   10 gfile=trim(sens_data)//gfile
      open(iug2,file=gfile,status="old",position="rewind",err=910)
c ebins(1/2/3,ng) is min/max/avg. read bounds as max(1:neg) and min(neg)
      do ng=1,neg
        read(iug2,*)l,ebins(2,ng)
      end do ! ng
      read(iug2,*)l,ebins(1,neg)
c populate ebins(1)
      ebins(1,1:neg-1)=ebins(2,2:neg)
c populate ebins(3)
      ebins(3,1:neg)=0.5d0*(ebins(1,1:neg)+ebins(2,1:neg))
      close(iug2)
      go to 999
c
  900 write(*,'("error. ",a," does not exist. ",
     1 "available group structures are:")')trim(gfile)
      write(iuo,'("error. ",a," does not exist. ",
     1 "available group structures are:")')trim(gfile)
      write(*,'(5x,a,5x,a)')
     1 (adjustr(group_structure_files(n)),n=1,nf)
      write(iuo,'(5x,a,5x,a)')
     1 (adjustr(group_structure_files(n)),n=1,nf)
      write(*,'("to add a group structure file, see source ",
     1 "subroutine set_ebins.F.")')
      write(iuo,'("to add a group structure file, see source ",
     1 "subroutine set_ebins.F.")')
      call stoponerror
  910 write(*,'("set_ebins. error opening ",a,".")')trim(gfile)
      write(iuo,'("set_ebins. error opening ",a,".")')trim(gfile)
      call stoponerror
c
  999 return
      end subroutine
