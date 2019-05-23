      subroutine asymptote(niso0)
c read partisn outputs and compute the Feynman Y asymptote
c using flux moment files.
c write partisn sources for sensitivities.
c equation numbers in this routine refer to
c John Mattingly, "Computation of Neutron Multiplicity Statistics Using
c Deterministic Transport," IEEE Transactions on Nuclear Science,
c vol. 59, no. 2, 314-322 (2012)
c use chivec until we figure out how to do this with a matrix.
      use COMS
      implicit none
      integer niso0 ! number of isotopes in snxedt
      real*8 ipsrc(it),ipfis(it),ipsrc_s(0:nr),ipfis_s(0:nr),
     1 twossrc,twosfis,r2,csrc,cfis
      real*8 mom2src(neg,0:nm),mom2fis(neg,0:nm)
c sources for sensitivities
      real*8 fsrc1(0:neg,it),fsrc2(0:neg,it),asrc(neg,it),
     1 gams(neg,0:nmom-1,it,jt)
      integer neg0,ndir0,nitm0,njtm0
      integer i,ir,i0,j,jr,nf,ng,m,m1,m2
c DEBUG_ALEX use this for feynman y central difference
      integer ia(4)
      real*8 chi_alex(30),tmp(2)
      data chi_alex/
     1 4.524071d-05,9.441084d-05,2.867205d-04,1.423820d-03,7.233329d-03,
     1 2.139274d-02,1.200992d-01,1.014603d-01,1.178932d-01,1.206688d-01,
     1 1.111523d-01,1.718359d-01,1.076289d-01,5.970784d-02,3.028628d-02,
     1 2.219630d-02,5.100366d-03,1.160542d-03,2.578353d-04,5.836856d-05,
     1 1.345363d-05,3.133590d-06,7.489831d-07,1.851863d-07,4.795220d-08,
     1 1.156804d-08,3.022614d-09,9.545161d-10,3.267145d-10,1.840000d-10/
c DEBUG_ALEX
c
      ipsrc(1:it)=0.d0
      ipfis(1:it)=0.d0
      ipsrc_s(0:nr)=0.d0
      ipfis_s(0:nr)=0.d0
c twossrc is the source excess variance, Eq. (6).
c twosfis is the fission excess variance, Eq. (8).
      twossrc=0.d0
      twosfis=0.d0
c
c fsrc1 is the source for phi_s,1
c fsrc2 is the source for phi_1
c asrc is the source for psi*_2
      fsrc1(0:neg,1:it)=0.d0
      fsrc2(0:neg,1:it)=0.d0
      asrc(1:neg,1:it)=0.d0
      gams(1:neg,0:nmom-1,1:it,1:jt)=0.d0
c
c read partisn fluxes.
      call rddantm(1,0) ! forward ss
      call rddantm(2,0) ! adjoint ss
c read material cross sections.
      call rdxsmg
c
c this is to dimension things correctly in subroutine integrals.
c for 2d problems with iangflux=0, the angular fluxes are dimensioned
c with 1's in coms.
      if(iangflux.eq.1)then
        neg0=neg
        ndir0=ndir
        nitm0=nitm
        njtm0=njtm
      else if(iangflux.eq.0)then
        neg0=1
        ndir0=1
        nitm0=1
        njtm0=1
      end if
c
c compute multiplicity moment ratios. 
      call multmom(mom2src,mom2fis)
c
c compute I0(r), Eq. (7)
      j=1
      nf=0
      do ir=1,nr
        m1=mat(ir,j)
        if(ismat(m1).eq.0)then
          nf=nf+iints(ir)
          cycle
        end if
        do i=1,iints(ir)
          nf=nf+1
          do ng=1,neg
            ipsrc(nf)=ipsrc(nf)+chisrc(ng,m1)*amom(ng,0,nf,j)
          end do ! ng
        end do ! i
      end do ! ir
c     write(*,'("chisrc",1p2e17.9)')chisrc(1:2,m1)
c
c compute I(r), Eq. (9)
c DEBUG_ALEX use this for feynman y central difference
c divide by rho because the macrxs file already includes it,
c but it was used in rdxsmg. so remove the extra factor of rho.
c     write(*,'("warning. wrong chivec, modify nusigf.")')
c     write(iuo,'("warning. wrong chivec, modify nusigf.")')
c     do i=1,nm
c       if(rho(i).eq.0.d0)cycle
c       nusigf(1:neg,i)=nusigf(1:neg,i)/rho(i)
c     end do ! i
c     open(iup,file="chibin",status='old',position='rewind',
c    1 form='unformatted',access='stream')
c     read(iup)(ia(i),i=1,4)
c     do ng=1,neg
c       read(iup)(chivec(ng,i),i=1,nxs)
c     end do
c     close(iup)
c     do ng=1,neg
c       write(*,'("asymptote chivec,nusigf",i4,1p2e17.9)')ng,
c    1   chivec(ng,1),nusigf(ng,1)
c       write(iuo,'("asymptote chivec,nusigf",i4,1p2e17.9)')ng,
c    1   chivec(ng,1),nusigf(ng,1)
c     end do ! ng
c
c     if(neg.eq.30)then
c       chivec(1:neg,1)=chi_alex(1:neg)
c     else if(neg.eq.1)then
c       chivec(1:neg,1)=1.d0
c     end if
c DEBUG_ALEX
c
      j=1
      nf=0
      do ir=1,nr
        if(ifcel(ir,j).eq.0)then
          nf=nf+iints(ir)
          cycle
        end if
        m1=mat(ir,j)
        do i=1,iints(ir)
          nf=nf+1
          do ng=1,neg
            ipfis(nf)=ipfis(nf)+chivec(ng,m1)*amom(ng,0,nf,j)
          end do ! ng
        end do ! i
      end do ! ir
c
c I0 and I on surfaces
      if(calc2d.eq.0)then ! will not be 1 yet in this routine
        i0=1
      else if(calc2d.eq.-1)then
        i0=0
      end if
      j=1 ! this is a 1d routine
      do ir=i0,nr
        i=iindex(ir)+1
        do jr=1,nz ! this will be only 1 for now
          m1=mat(ir,jr)
          m2=mat(ir+1,jr)
          do ng=1,neg
            do m=1,isn
              ipsrc_s(ir)=ipsrc_s(ir)+wgt(m)*afadj(ng,m,i,j)
     1         *(chisrc(ng,m1)-chisrc(ng,m2))
              ipfis_s(ir)=ipfis_s(ir)+wgt(m)*afadj(ng,m,i,j)
     1         *(chivec(ng,m1)-chivec(ng,m2))
            end do ! m
          end do ! ng
        end do ! jr
      end do ! ir
c     write(iuo,'("asy src",1p2e17.9)')ipsrc_s(0)**2,ipsrc_s(nr)**2
c     write(iuo,'("asy fis",1p2e17.9)')ipfis_s(0)**2,ipfis_s(nr)**2
c
c compute 2S0 and 2S [Eqs. (6) and (8)]. note Eq. (8) has a typo,
c there should be a phi(r,E,mu) in the integrand.
      j=1
      nf=0
      tmp(1:2)=0.d0
      do ir=1,nr
        m1=mat(ir,j)
        if(sum(mom2src(1:neg,m1)).eq.0.d0.and.
     1   sum(mom2fis(1:neg,m1)).eq.0.d0)then
          nf=nf+iints(ir)
          cycle
        end if
        do i=1,iints(ir)
          nf=nf+1
          csrc=ipsrc(nf)*ipsrc(nf)*dv(nf,j)
          cfis=ipfis(nf)*ipfis(nf)*dv(nf,j)
          tmp(1)=tmp(1)+cfis
          do ng=1,neg
            tmp(2)=tmp(2)+cfis*fmom(ng,0,nf,j)
            twossrc=twossrc+csrc*mom2src(ng,m1)
            twosfis=twosfis+cfis*mom2fis(ng,m1)*fmom(ng,0,nf,j)
            fsrc1(0,nf)=fsrc1(0,nf)+mom2src(ng,m1)
            fsrc2(0,nf)=fsrc2(0,nf)+mom2fis(ng,m1)
     1       *fmom(ng,0,nf,j)
            asrc(ng,nf)=mom2fis(ng,m1)*ipfis(nf)*ipfis(nf)
          end do ! ng
          fsrc1(1:neg,nf)=fsrc1(0,nf)*ipsrc(nf)*chisrc(1:neg,m1)
          fsrc2(1:neg,nf)=fsrc2(0,nf)*ipfis(nf)*chivec(1:neg,m1)
        end do ! i
      end do ! ir
c
      r2=twossrc+twosfis
      if(lkg.gt.0.d0)then
        yasym=r2/lkg
        sm2=0.5d0*r2/lkg**2
      else
        yasym=0.d0
        sm2=0.d0
      end if
c
      write(*,'(/,"results for subcritical multiplication")')
      write(iuo,'(/,"results for subcritical multiplication")')
      write(*,'(2x,"twossrc, twosfis, r2, lkg",1p4e17.9)')twossrc,
     1 twosfis,r2,lkg
      write(iuo,'(2x,"twossrc, twosfis, r2, lkg",1p4e17.9)')twossrc,
     1 twosfis,r2,lkg
      write(*,'(2x,"feynman_y_asymptote, sm2",1p2e17.9)')yasym,sm2
      write(iuo,'(2x,"feynman_y_asymptote, sm2",1p2e17.9)')yasym,sm2
c     write(*,'(2x,"feynman_y_asymptote, sm2, lkg, r2",1p4e17.9)')
c    1 yasym,sm2,lkg,r2
c     write(iuo,'(2x,"feynman_y_asymptote, sm2, lkg, r2",1p4e17.9)')
c    1 yasym,sm2,lkg,r2
c     write(*,'(2x,"integrals, phi*^2 and phi*^2 phi",1p2e17.9)')
c    1 tmp(1:2)
c     write(iuo,'(2x,"integrals, phi*^2 and phi*^2 phi",1p2e17.9)')
c    1 tmp(1:2)
c to compare with John's NSE results
c     write(*,'("mean and variance in 4096 us")')
c     write(*,'(1p2e14.6)')lkg*4096.d-6,(yasym+1.d0)*lkg*4096.d-6
c
c write partisn inputs for sensitivities.
c combine the two forward sources. separate if induced fission and
c spontaneous fission sources are needed separately.
      j=1
      fsrc1(1:neg,1:it)=fsrc1(1:neg,1:it)+fsrc2(1:neg,1:it)
      gams(1:neg,0,1:it,j)=fsrc1(1:neg,1:it)
      call wrdantga(1,0,3,niso0,gams)
c     gams(1:neg,0,1:it,j)=fsrc2(1:neg,1:it)
c     call wrdantga(2,0,3,niso0,gams)
      gams(1:neg,0,1:it,j)=asrc(1:neg,1:it)
      call wrdantga(3,0,3,niso0,gams)
c
c write other data needed for sensitivities.
      call wrsmbin(r2,ipsrc,ipfis,ipsrc_s,ipfis_s)
c
      return
      end
