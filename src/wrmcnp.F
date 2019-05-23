      subroutine wrmcnp
c write an mcnp input file inpi -- at this point, only for plotting.
c a plot command file com is also written.
c for redoin/lnk3dnt, the geometry will be correct, but the materials
c may not be.
      use COMS
      implicit none
      integer i,j,j1,k2
      character bfile*4,line*80
c
c this is done in sensmg.py
c     if(ilnk3dnt.eq.1)then
c       write(*,'("error. cannot make mcnp file from redoin/lnk3dnt",
c    1   " files. first use -wrsensmg yes.")')
c       write(iuo,'("error. cannot make mcnp file from redoin/lnk3dnt",
c    1   " files. first use -wrsensmg yes.")')
c       call stoponerror
c     end if
c
      bfile="inpi"
      open(iumc, file=bfile,status='unknown')
      write(iumc,'(a)')trim(id)
      write(iumc,'("c for plotting only.")')
      write(iumc,'("c ilkg=",i1)')ilkg
c
c cells.
      if(calc2d.eq.0)then
        do i=1,lastr
          write(line,'(i6,i7,f14.8,2i8,14x,"imp:n,p=1")')1000+i,
     1     mat(i,1),-rho(mat(i,1)),-(1000+i),(1000+i-1)
          if(i.eq.1)line(36:43)="        "
          if(mat(i,1).eq.0)line(14:27)="              "
          write(iumc,'(a)')trim(line)
        end do ! i
        write(iumc,'(i6,i7,i22,i8,14x,"imp:n,p=1")')999998,0,-999999,
     1   1000+lastr
        write(iumc,'(i6,i7,i22,22x,"imp:n,p=0")')999999,0,999999
      else if(calc2d.eq.1)then
c this logic is from wrdantnm.
        if(ilkg.eq.0)then
c         j1=0 ! check this
        elseif(ilkg.eq.1)then
          j1=1
        end if
        do j=j1,nz
          do i=1,lastr
            write(line,'(i6,i7,f14.8,4i8,2x,"imp:n,p=1")')
     1       j*1000+i,mat(i,j),-rho(mat(i,j)),-(1000+i),(1000+i-1),
     2       -(2000+j),(2000+j-1)
            if(i.eq.1)line(36:43)="        "
            if(mat(i,j).eq.0)line(14:27)="              "
            write(iumc,'(a)')trim(line)
          end do ! i
        end do ! j
        write(iumc,'(i6,i7,i14,2x,"(",i8,":",i8,":",i8,")",4x,
     1   "imp:n,p=1")')999998,0,-999999,1000+lastr,-(2000+j1-1),
     1   (2000+nz)
        write(iumc,'(i6,i7,i14,34x,"imp:n,p=0")')999999,0,999999
      end if
c
c surfaces.
      write(iumc,'(" ")')
      if(calc2d.eq.0)then
        do i=1,lastr
          write(iumc,'(i6,"  so  ",1pe15.8)')1000+i,r(i)
        end do ! i
      elseif(calc2d.eq.1)then
        do i=1,lastr
          write(iumc,'(i6,"  cz  ",1pe15.8)')1000+i,r(i)
        end do ! i
        if(ilkg.eq.0)then
c         j1=0 ! check this
        elseif(ilkg.eq.1)then
          j1=0
        end if
        do j=j1,nz
          write(iumc,'(i6,"  pz  ",1pe15.8)')2000+j,z(j)
        end do ! j
      end if
      write(iumc,'(i6,"  so  ",1pe7.0)')999999,1.d+6
c
c data.
      write(iumc,'(/,"mode p")')
      write(iumc,'("nps 1e9")')
      write(iumc,'("prdmp j 1e8")')
      write(iumc,'("rand gen=2 seed=1000001")')
      write(iumc,'("notrn $ to match ray-tracing")')
      write(iumc,'("phys:p j j 1 j j j j $ no coherent scattering")')
      write(iumc,'("sdef")')
      k2=0
      do i=1,nm
        write(iumc,'("m",i6.6)')i
        do j=1,ncb(i)
          k2=k2+1
          write(iumc,'(i14,1pe14.6)')nint(blk(1,k2)),-blk(2,k2)
        end do ! j
      end do ! i
      write(iumc,'("c")')
      write(iumc,'("print -30",/)')
      write(iumc,'("end of input")')
      close(iumc)
c
c com file for convenience. plot x-z.
      if(ilkg.eq.0)then
c       j1=0 ! check this
      elseif(ilkg.eq.1)then
        j1=0
      end if
      bfile="com"
      open(iumc, file=bfile,status='unknown')
      write(iumc,'("ba 1 0 0 0 0 1 scale 1 &")')
      write(iumc,'("la 0 1 cel &")')
      write(iumc,'("or ",1p3e14.6," &")')0.d0,0.d0,(z(j1)+z(nz))/2.d0
      write(iumc,'("ext ",1pe14.6," &")')1.1d0*max(r(lastr),
     1 (z(j1)+z(nz))/2.d0)
      write(iumc,'("interact")')
      close(iumc)
c
      write(*,'("mcnp plot files created. use",/,
     1 "  mcnp ip i=inpi com=com")')
      write(iuo,'("mcnp plot files created. use",/,
     1 "  mcnp ip i=inpi com=com")')
      return
      end
