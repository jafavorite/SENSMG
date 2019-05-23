      subroutine output4(sintr,sintz)
c output interface location derivatives for cylinders.
      use COMS
      implicit none
      real*8 sintr(neg,0:nr,nz),sintz(neg,nr,0:nz)
      integer i1,i2,i,j1,j2,j,ng,nblock
c there are 30 radii and 30 heights in each block; can be changed here.
c nblock+3 <= 999.
      parameter (nblock=30)
c
      write(iur,'(" radial surfaces")')
c summary table first
      write(iur,'(" summary (sums over groups)")')
      i1=1
      do while(i1.le.nr)
        i2=min(i1+nblock-1,nr)
        do i=i1,i2
          if(i.eq.i1)then
            write(iur,'(16x)',advance='no')
          end if
          if(i.eq.i2)then
            write(iur,'("radius",i4.4,4x)',advance='yes')i
          else
            write(iur,'("radius",i4.4,4x)',advance='no')i
          end if
        end do ! i
        do j=nz,1,-1
          write(iur,'("   axis",i4.4,3x,1p999e14.6)')j,
     1     (sum(sintr(1:neg,i,j)),i=i1,i2)
        end do ! j
        i1=i2+1
      end do ! i1
c now the details
      do j=nz,1,-1
        write(iur,'(" axis",i4.4)')j
        i1=1
        do while(i1.le.nr)
          write(iur,'(1x,"group",2x,"e_lower",7x,"e_upper",7x,
     1     "e_average",5x)',advance='no')
          i2=min(i1+nblock-1,nr)
          do i=i1,i2
            if(i.eq.i2)then
              write(iur,'("radius",i4.4,4x)',advance='yes')i
            else
              write(iur,'("radius",i4.4,4x)',advance='no')i
            end if
          end do ! i
          do ng=1,neg
            write(iur,'(i6,1p999e14.6)')ng,ebins(1:3,ng),
     1       (sintr(ng,i,j),i=i1,i2)
          end do ! ng
c 42x is 14*3 for the energy groups.
          write(iur,'(" total",42x,1p999e14.6)')
     1     (sum(sintr(1:neg,i,j)),i=i1,i2)
          i1=i2+1
        end do ! i1
      end do ! j
c
      write(iur,'(/," axial surfaces")')
c summary table first
      write(iur,'(" summary (sums over groups)")')
      i1=1
      do while(i1.le.nr)
        i2=min(i1+nblock-1,nr)
        do i=i1,i2
          if(i.eq.i1)then
            write(iur,'(16x)',advance='no')
          end if
          if(i.eq.i2)then
            write(iur,'("radius",i4.4,4x)',advance='yes')i
          else
            write(iur,'("radius",i4.4,4x)',advance='no')i
          end if
        end do ! i
        do j=nz,0,-1
          write(iur,'("   axis",i4.4,3x,1p999e14.6)')j,
     1     (sum(sintz(1:neg,i,j)),i=i1,i2)
        end do ! j
        i1=i2+1
      end do ! i1
c now the details
      do i=1,nr
        write(iur,'(" radius",i4.4)')i
        j2=nz
        do while(j2.ge.0)
          write(iur,'(1x,"group",2x,"e_lower",7x,"e_upper",7x,
     1     "e_average",5x)',advance='no')
          j1=max(j2-nblock+1,0)
          do j=j2,j1,-1
            if(j.eq.j1)then
              write(iur,'("axis",i4.4,6x)',advance='yes')j
            else
              write(iur,'("axis",i4.4,6x)',advance='no')j
            end if
          end do ! j
          do ng=1,neg
            write(iur,'(i6,1p999e14.6)')ng,ebins(1:3,ng),
     1       (sintz(ng,i,j),j=j2,j1,-1)
          end do ! ng
c 42x is 14*3 for the energy groups.
          write(iur,'(" total",42x,1p999e14.6)')
     1     (sum(sintz(1:neg,i,j)),j=j2,j1,-1)
          j2=j1-1
        end do ! j1
      end do ! i
c
      return
      end subroutine
