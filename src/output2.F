      subroutine output2(nm1,sden)
c output mass density sensitivities.
      use COMS
      implicit none
      integer nm1
      real*8 sden(neg,nm1)
      integer i1,i2,i,ng,nblock
c there are 30 materials in each block; can be changed here.
c nblock+3 <= 999.
      parameter (nblock=30)
c
c summary table first
      write(iur,'(" summary (sums over groups)")')
      do i=1,nm1
        write(iur,'(" mat",i6.6,1pe15.6)')i,sum(sden(1:neg,i))
      end do ! i
c now the details
      i1=1
      do while(i1.le.nm1)
        i2=min(i1+nblock-1,nm1)
        write(iur,'(1x,"group",2x,"e_lower",7x,"e_upper",7x,
     1   "e_average",5x)',advance='no')
        do i=i1,i2
          if(i.eq.i2)then
            write(iur,'("mat",i6.6,5x)',advance='yes')i
          else
            write(iur,'("mat",i6.6,5x)',advance='no')i
          end if
        end do ! i
        do ng=1,neg
          write(iur,'(i6,1p999e14.6)')ng,ebins(1:3,ng),
     1     (sden(ng,i),i=i1,i2)
        end do ! ng
c 42x is 14*3 for the energy groups.
        write(iur,'(" total",42x,1p999e14.6)')
     1   (sum(sden(1:neg,i)),i=i1,i2)
        i1=i2+1
      end do
c
      return
      end subroutine
