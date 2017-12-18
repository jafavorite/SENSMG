      program make_template
c in ?.sensmg.log, replace executable paths with /usr/projects/...
c compile with
c   ifort -o make_template -r8 -i8 -check bounds -warn unused -traceback make_template.F
      implicit none
      integer nfiles
      parameter (nfiles=8)
c iflag=1/2/3 convert original files/convert test files/named file
      integer iflag
      integer n0,nf,l1,l2,ios
      character files(nfiles)*120,ifile*130,ofile*130,line*9999,
     1 ctmp*120,cflag*1
      data files/
     1 "A.sensmg.log",
     1 "B.sensmg.log",
     1 "C.sensmg.log",
     1 "D.sensmg.log",
     1 "E.sensmg.log",
     1 "F.sensmg.log",
     1 "G.sensmg.log",
     1 "H.sensmg.log"/
c
      call get_command_argument(1, cflag)
      read(cflag,'(i1)')iflag
      if(iflag.lt.1.or.iflag.gt.3)then
        write(*,'("in ?.sensmg.log, replace executable paths with ",
     1   "/usr/projects/...")')
        write(*,'("enter iflag=1/2/3 ",
     1   "convert original files/convert test files/convert one file ",
     2   "to be named")')
        stop
      end if
c
      if(iflag.eq.1.or.iflag.eq.2)then
        n0=nfiles
      else if(iflag.eq.3)then
        n0=1
      end if
      do nf=1,n0
        if(iflag.eq.1)then
          ifile=trim(files(nf))//".org"
        else if(iflag.eq.2)then
          ifile=trim(files(nf))
        else if(iflag.eq.3)then
          write(*,'("enter file name: ")')
          read(*,*)ifile
        end if
        open(unit=2,file=ifile,status="old",IOSTAT=ios) 
        if(ios.ne.0)then
          write(*,'("*** error. could not open ",a)')trim(ifile)
          cycle
        end if 
        if(iflag.eq.1.or.iflag.eq.2)then
          ofile=trim(files(nf))//".tmp"
        else if(iflag.eq.3)then
          ofile=trim(ifile)//".tmp"
        end if
        write(*,'(a)')trim(ifile)//" --> "//trim(ofile)
        open(unit=3,file=ofile,status="unknown") 
c
   10   read(2,'(a)',IOSTAT=ios)line
        if(ios.ne.0)go to 20
c they are all sensmg.log but never mind
        if(index(ifile,"sensmg.log").ne.0)then
          if(index(line,"COMMAND LINE").ne.0)then
            l1=index(line,"=")+1
            l2=index(line,".py ")+3
            ctmp=line(l2:len_trim(line))
            line=line(1:l1-1)//
     1"/usr/projects/data/nuclear/working/sensitivities/bin/sensmg.py"
     2//trim(ctmp)
          else if(index(line,"SENS_PART_CODE").ne.0)then
            l1=index(line,"=")+1
            line=line(1:l1-1)//
     1"/usr/projects/data/nuclear/working/sensitivities/bin/sensmg"
          else if(index(line,"CONTROL SCRIPT").ne.0)then
            l1=index(line,"=")+1
            line=line(1:l1-1)//
     1"/usr/projects/data/nuclear/working/sensitivities/bin/sensmg.py"
          else if(index(line,"SOURCES").ne.0.and.
     1     index(line,"SOURCES_DIR").eq.0)then
            l1=index(line,"=")+1
            line=line(1:l1-1)//
     1"/usr/projects/data/nuclear/working/sensitivities/sources4c/"//
     2"bin/sources4c.jaf"
          end if
        end if
        write(3,'(a)')trim(line)
        go to 10
c
   20 close(2)
      close(3)
      end do ! nf
c
      end program
