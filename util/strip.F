      program strip
c in n.nnnnnnE-xx, replace last iq n's with blanks and replace
c executable paths with blanks
c compile with
c   ifort -o strip -r8 -i8 -check bounds -warn unused -traceback strip.F
      implicit none
      integer nfiles,iq
      parameter (nfiles=39,iq=3)
c iflag=1/2/3 convert original files/convert test files/named file
      integer iflag
      integer n0,nf,i,j,l,ios
      character files(nfiles)*120,ifile*130,ofile*130,line*9999,cflag*1
      data files/
     1 "A.sens_k_r","A.sens_k_x","A.sensmg.log",
     1 "B.sens_k_r","B.sens_k_x","B.sensmg.log",
     1 "C.sens_k_r","C.sens_k_x","C.sens_rr_r","C.sens_rr_x",
     1  "C.sensmg.log",
     1 "D.sens_k_r","D.sens_k_x","D.sens_rr_r","D.sens_rr_x",
     1  "D.sensmg.log",
     1 "E.sens_k_r","E.sens_k_x","E.sens_rr_r","E.sens_rr_x",
     1  "E.sensmg.log",
     1 "F.sens_a_r","F.sens_a_x","F.sens_rr_r","F.sens_rr_x",
     1  "F.sensmg.log",
     1 "G.sens_l_r","G.sens_l_x","G.sens_rr_r","G.sens_rr_x",
     1  "G.sensmg.log",
     1 "H.sens_l_r","H.sens_l_x","H.sensmg.log",
     1 "I.sens_l_r","I.sens_l_x","I.sens_s_r","I.sens_s_x",
     1  "I.sensmg.log"/
     
c
      call get_command_argument(1, cflag)
      read(cflag,'(i1)')iflag
      if(iflag.lt.1.or.iflag.gt.3)then
        write(*,'("in specified files, in `n.nnnnnnE-xx`, ",
     1   "replace last ",i1," n`s with blanks.")')iq
        write(*,'("the number can be changed (iq).")')
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
        if(index(ifile,"sensmg.log").ne.0)then
          if(index(line,"COMMAND LINE").ne.0)then
            l=index(line,"=")
            do i=l+1,index(line,".py ")+2
              line(i:i)=" "
            end do ! i
          else if(index(line,"SENSMG_EXE").ne.0.or.
     1     index(line,"CONTROL SCRIPT").ne.0)then
            l=index(line,"=")
            do i=l+1,len_trim(line)
              line(i:i)=" "
            end do ! i
          end if
        end if
        do i=1,len_trim(line)-3
          if(line(i:i).ne."E")cycle
          l=index("-+",line(i+1:i+1))
          if(l.eq.0)cycle
          l=1
          do j=1,iq
            l=l*index("0123456789",line(i-j:i-j))
          end do ! j
          if(l.eq.0)cycle
          do j=1,iq
            line(i-j:i-j)=" "
          end do ! j
        end do ! i
        write(3,'(a)')trim(line)
        go to 10
c
   20 close(2)
      close(3)
      end do ! nf
c
      end program
