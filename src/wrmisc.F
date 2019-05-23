      subroutine wrmisc
c write misc input files
      use COMS
      implicit none
      integer i,ir,imod,j,jr,k2,n,ng,nm1
      character mfile*22,ofile*22,sfile*22
      real*8 rho1,rh(nmxsp)
c     character ce*100
c
      if(ilnk3dnt.eq.0)then
        nm1=nm
      else if(ilnk3dnt.gt.0)then
        nm1=nr*nz
      end if
      k2=0
      jr=0
      do i=1,nm1
        if(ilnk3dnt.gt.0)then
          ir=int(real(i,8)/real(nz,8))+1
          jr=jr+1
          imod=mod(i,nz)
          if(imod.eq.0)then
            ir=ir-1
          else if(imod.eq.1)then
            jr=1
          end if
          if(idclnk3(ir,jr,1).eq.0)cycle
        end if
c write misc input file.
        write(mfile,'("misc/m",i6.6,"_misc.inp")')i
        write(ofile,'("m",i6.6,"_misc.out")')i
        write(sfile,'("m",i6.6,"_misc.src")')i
        open(iun,file=mfile,status='unknown')
        write(iun,'("output = ",a)')trim(ofile)
        write(iun,'("srcout = ",a)')trim(sfile)
        write(iun,'("decayfile = endf7.dk.xml")')
        write(iun,'("abundfile = nist.na.xml")')
        write(iun,'("pelib = endf7")')
        write(iun,'("format = v")')
        write(iun,'("particle = 1")')
        write(iun,'("matspec = ")')
        if(ilnk3dnt.eq.0)then
          do j=1,ncb(i)
            k2=k2+1
            write(iun,'(12x,i7,1pe15.7)')nint(blk(1,k2)),-blk(2,k2)
          end do ! j
          write(iun,'("density = ",0pf15.8)')-rho(i)
        else if(ilnk3dnt.gt.0)then
          rho1=0.d0
          do n=1,nmxsp
            if(idclnk3(ir,jr,n).eq.0)exit
            rh(n)=denlnk3(ir,jr,n)*atwt(idclnk3(ir,jr,n))/avg
            rho1=rho1+rh(n)
          end do ! n
          do n=1,nmxsp
            k2=idclnk3(ir,jr,n)
            if(k2.eq.0)exit
            rh(n)=rh(n)/rho1
            write(iun,'(12x,i7,1pe15.7)')nint(blk(1,k2)),-rh(n)
          end do ! n
          write(iun,'("density = ",0pf15.8)')-rho1
        end if
        write(iun,'("ebins = ")')
        write(iun,'(7x,1p5e12.5)')ebins(1,neg),(ebins(2,ng),ng=neg,1,-1)
        close(iun)
      end do ! i
c
      return
      end subroutine
