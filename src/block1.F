      subroutine block1(nf,mt,n0)
c write block1 of a partisn input file.
      use COMS
      implicit none
c nf is file number to write to
c mt is materials for partisn input; =nm+nrrx for forward
c  and generalized adjoint
c n0 is niso to use
      integer nf,mt,n0
c
      write(nf,'("/ * * * * block i * * * *")')
      if(calc2d.eq.0)then
        write(nf,'(" igeom=sphere isn=",i4," ngroup=",i5)')isn,neg
c       write(nf,'("/maxscm=20000000 maxlcm=9000000")')
      else if(calc2d.eq.-1)then
        write(nf,'(" igeom=slab isn=",i4," ngroup=",i5)')isn,neg
      else if(calc2d.eq.1)then
        write(nf,'(" igeom=r-z isn=",i4," ngroup=",i5)')isn,neg
c       write(nf,'("/maxscm=60000000 maxlcm=6000000")')
        write(nf,'(" maxscm=600000000")')
      end if
      write(nf,'(" niso=",i6," mt=",i4," nzone=",i4)')n0,mt,nm
      write(nf,'(" im=",i3," it=",i6)')lastr,it
      if(calc2d.eq.1)then
        if(ilkg.eq.1)then
          write(nf,'(" jm=",i3," jt=",i6)')nz,jt
        elseif(ilkg.eq.0)then
          write(nf,'(" jm=",i3," jt=",i6)')nz+2,jt
        end if
        write(nf,'(" iquad=6")')
c       write(nf,'(" iquad=4")') ! Galerkin
      end if
      write(nf,'(" t")')
c
      return
      end subroutine
c
      subroutine block2(nf)
c write block2 of a partisn input file.
      use COMS
      implicit none
      integer nf ! file number to write to
      integer i,j
c
      write(nf,'("/ * * * * block ii * * * *")')
      write(nf,'(" xmesh=",1pe16.8)')r(0)
      do i=1,lastr
        write(nf,'(7x,1pe16.8)')r(i)
      end do ! i
      if(calc2d.eq.1)then
        if(ilkg.eq.1)then
          write(nf,'(" ymesh=",1pe16.8)')z(0)
          do j=1,nz
            write(nf,'(7x,1pe16.8)')z(j)
          end do ! j
        elseif(ilkg.eq.0)then
c         write(nf,'(" ymesh=",1pe15.8)')z(0)
c         do j=firstz,lastz
c           write(nf,'(6x,1pe16.8)')z(j)
c         end do ! j
        end if
      end if
      write(nf,'(" xints=")')
      write(nf,'(10i6)')(iints(i),i=1,lastr)
      if(calc2d.eq.1)then
        write(nf,'(" yints=")')
        if(ilkg.eq.1)then
          write(nf,'(10i6)')(jints(i),i=1,nz)
        elseif(ilkg.eq.0)then
c         write(nf,'(10i6)')(jints(i),i=1,nz)
        end if
      end if
      if(calc2d.eq.0.or.calc2d.eq.-1)then
        write(nf,'(" zones=")')
        do i=1,lastr
          write(nf,'(4x,i6)')mat(i,1)
        end do ! i
      else if(calc2d.eq.1)then
        if(ilkg.eq.1)then
          write(nf,'(" zones=")')
          do j=1,nz
            write(nf,'(8x,20i3)')(mat(i,j),i=1,lastr)
            write(nf,'(8x,";")')
          end do ! j
        elseif(ilkg.eq.0)then
c         write(nf,'(" zones=",60i3)')(mat(i,firstz+1),i=1,lastr)
c         write(nf,'(";")')
c         do j=firstz+2,lastz
c           write(nf,'(8x,60i3)')(mat(i,j),i=1,lastr)
c           write(nf,'(";")')
c         end do ! j
        end if
      end if
      write(nf,'(" t")')
c
      return
      end subroutine
