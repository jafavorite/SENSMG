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
      if(ilnk3dnt.eq.0)then
        write(nf,'(" niso=",i6," mt=",i6," nzone=",i6)')n0,mt,nm
      else if(ilnk3dnt.eq.1)then
        write(nf,'(" niso=",i6," mt=",i6," nzone=",i6)')n0,mt,mt
      else if(ilnk3dnt.eq.2)then
        write(nf,'(" niso=",i6," mt=",i6," nzone=",i6)')n0,mt,nm
      end if
      write(nf,'(" im=",i6," it=",i6)')lastr,it
      if(calc2d.eq.1)then
        if(ilkg.eq.1)then
          write(nf,'(" jm=",i6," jt=",i6)')nz,jt
        elseif(ilkg.eq.0)then
          write(nf,'(" jm=",i6," jt=",i6)')nz+2,jt
        end if
        write(nf,'(" iquad=6")')
c       write(nf,'(" iquad=4")') ! Galerkin
      end if
c use lnk3dnt for geometry and materials
      if(ilnk3dnt.gt.0)then
        write(nf,'(" fmmix=1")')
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
            write(nf,'(8x,10i6)')(mat(i,j),i=1,lastr)
            write(nf,'(8x,";")')
          end do ! j
        elseif(ilkg.eq.0)then
c         write(nf,'(" zones=",60i6)')(mat(i,firstz+1),i=1,lastr)
c         write(nf,'(";")')
c         do j=firstz+2,lastz
c           write(nf,'(8x,10i6)')(mat(i,j),i=1,lastr)
c           write(nf,'(";")')
c         end do ! j
        end if
      end if
      write(nf,'(" t")')
c
      return
      end subroutine
c
      subroutine block3(nf,ith,iw)
c write block3 of a partisn input file.
      use COMS
      implicit none
c nf is the file number to write to
c ith=0/1 forward/adjoint
c iw=0/1/2/3 called from wrdantnm/wrdantxs/wrdantga/outputxs
      integer nf,ith,iw
      integer ng
c
      if(iw.ne.3)then
        write(nf,'("/ * * * * block iii * * * *")')
      end if
      if(ith.eq.0)then
c ndilib
        if(libname(1:6).ne.'kynea3'.and.libname(1:5).ne.'scale'.and.
     1   libname(1:3).ne.'b96'.and.libname(1:7).ne.'special')then
          write(nf,'(" lib=ndilib")')
c         if(libname(1:5).ne."mt71x")then
            write(nf,'(" libname=",a)')trim(libname)
c         end if
          if(fissneut.eq.1)then
            write(nf,'(" fissneut=1 / total nu-bar")')
          else if(fissneut.eq.0)then
            write(nf,'(" fissneut=0 / prompt nu-bar")')
          end if
          if(fissdata.eq.1)then
            if(partisnver(1:1).eq."5")then
              write(nf,'(" nochimat=0 / chi matrix")')
            else
              write(nf,'(" fissdata=1 / chi matrix")')
            end if
          else if(fissdata.eq.2)then
            if(partisnver(1:1).eq."5")then
              write(nf,'(" nochimat=1 / chi vector")')
            else
              write(nf,'(" fissdata=2 / chi vector")')
            end if
          else if(fissdata.eq.0)then
            if(partisnver(1:1).ne."5")then
              write(nf,'(" fissdata=0 / fission transfer matrix")')
            end if
          end if
        else if(libname(1:6).eq."kynea3".or.libname(1:5).eq."scale")then
          write(nf,'(" lib=bxslib")')
          write(nf,'(" libname=",a)')trim(libname)
          if(partisnver(1:1).eq."5")then
            write(nf,'(" nochimat=1 / chi vector")')
          else
            write(nf,'(" fissdata=2 / chi vector")')
          end if
        else if(libname(1:3).eq."b96")then
          if(libname(5:7).eq.'bxs'.or.libname(4:6).eq.'bxs')then
            write(nf,'(" lib=bxslib ifido=1 ititl=1 maxord=5 i2lp1=1")')
          else
            write(nf,'(" lib=xslib ifido=1 ititl=1 maxord=5 i2lp1=1")')
          end if
          write(nf,'(" libname=",a)')trim(libname)
          if(iw.ne.3)then
            if(libname(1:4).eq.'b96t')then
              write(nf,'(" iht=3 ihs=8 ihm=74")')
            else
              write(nf,'(" iht=3 ihs=4 ihm=70")')
            end if
          end if
        else if(libname(1:7).eq."special")then
          if(iw.eq.1)then
            write(nf,'("/ this file does not actually get run")')
          else if(iw.eq.3)then
            write(nf,'(" special nuclear data library.")')
          end if
          write(nf,'(" lib=macrxs")')
        end if
        if(iw.ne.3)then
          if(lng.eq.30)then
c           write(nf,'(" grpstr=default ggrpstr=default")')
          elseif(lng.eq.2)then
c           write(nf,'(" chivec=1. 0.")')
          elseif(lng.eq.47.and.libname(1:3).eq."b96")then ! bugle
          write(nf,'(" chivec=")')
          write(nf,'("  1.7061E-06 1.3729E-05 7.2920E-05 2.9617E-04")')
          write(nf,'("  7.8994E-04 2.0118E-03 4.7785E-03 1.0661E-02")')
          write(nf,'("  1.9656E-02 2.5942E-02 3.0160E-02 3.3180E-02")')
          write(nf,'("  3.4228E-02 3.5417E-02 3.9323E-02 4.4854E-02")')
          write(nf,'("  5.0259E-02 5.5792E-02 5.8909E-02 5.9571E-02")')
          write(nf,'("  5.9478E-02 5.8403E-02 5.5965E-02 5.2404E-02")')
          write(nf,'("  4.7284E-02 3.9306E-02 3.1774E-02 2.5309E-02")')
          write(nf,'("  2.0978E-02 1.8808E-02 1.7564E-02 1.6839E-02")')
          write(nf,'("  1.5119E-02 1.1764E-02 8.1157E-03 5.5876E-03")')
          write(nf,'("  3.5931E-03 2.0584E-03 1.4148E-03 9.3733E-04")')
          write(nf,'("  5.5136E-04 3.1573E-04 2.0916E-04 1.3162E-04")')
          write(nf,'("  9.0458E-05 5.7090E-05 2.5183E-05")')
          write(nf,'("  20r 0.;")')
          end if
c ebins(1/2/3,ng) is min/max/avg.
          if(iw.ne.1)then
            write(nf,'(" ebound=")')
            write(nf,'(1p5e14.6)')(ebins(2,ng),ng=1,neg),ebins(1,neg)
          end if
        end if
      else if(ith.eq.1)then
c use macrxs from forward run.
        write(nf,'(" lib=macrxs")')
      end if
      if(iw.ne.3)then
        write(nf,'(" lng=",i3)')lng
        write(nf,'(" t")')
      end if
c
      return
      end subroutine
