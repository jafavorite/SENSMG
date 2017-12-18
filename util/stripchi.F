      program stripchi
c strip chi out of outputs.
c compile with
c   ifort -o stripchi -r8 -i8 -check bounds -warn unused -traceback stripchi.F
      implicit none
      integer nfiles
      parameter (nfiles=13)
c iflag=1/2/3 convert original files/convert test files/named file
      integer iflag
      integer ngroups(nfiles)
      integer n0,nf,ng,ios
      character files(nfiles)*120,ifile*130,ofile*130,line*9999,cflag*1
      data files/
     1 "A.sens_k_x",
     1 "B.sens_k_x",
     1 "C.sens_k_x","C.sens_rr_x",
     1 "D.sens_k_x","D.sens_rr_x",
     1 "E.sens_k_x","E.sens_rr_x",
     1 "F.sens_a_x","F.sens_rr_x",
     1 "G.sens_l_x","G.sens_rr_x",
     1 "H.sens_l_x"/
      data ngroups/79,30,70,70,30,30,30,30,30,30,30,30,30/
     
c
      call get_command_argument(1, cflag)
      read(cflag,'(i1)')iflag
      if(iflag.lt.1.or.iflag.gt.3)then
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
        if(line(2:15).eq."group  e_lower")then
          write(3,'(a)')trim(line)
          do ng=1,ngroups(nf)+1
            read(2,'(a)',IOSTAT=ios)line
            line(106:118)="             "
            write(3,'(a)')trim(line)
          end do ! ng
        else if(line(1:23).eq."fission chi sensitivity")then
          write(3,'(a)')trim(line)
          read(2,'(a)',IOSTAT=ios)line
          write(3,'(a)')trim(line)
          do ng=1,ngroups(nf)
            read(2,'(a)',IOSTAT=ios)line
            write(3,'(a)')line(1:6)
          end do ! ng
        else
          write(3,'(a)')trim(line)
        end if
        go to 10
c
   20 close(2)
      close(3)
      end do ! nf
c
      end program
