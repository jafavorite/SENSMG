      subroutine asymptote(niso0)
c read partisn outputs and compute the Feynman Y asymptote
c using flux moment files.
c write partisn sources for sensitivities.
c equation numbers in this routine refer to
c John Mattingly, "Computation of Neutron Multiplicity Statistics Using
c Deterministic Transport," IEEE Transactions on Nuclear Science,
c vol. 59, no. 2, 314-322 (2012)
      use COMS
      implicit none
      integer niso0 ! number of isotopes in snxedt
      real*8 chivec(neg,nm),ipsrc(it),ipfis(it),
     1 twossrc,twosfis,r2,csrc,cfis,sum1
      real*8 mrsrc(0:nm),mrfis(0:nm)
c sources for sensitivities
      real*8 fsrc1(neg,it),fsrc2(neg,it),asrc(neg,it),
     1 gams(neg,0:nmom-1,it,jt)
      integer neg0,ndir0,nitm0,njtm0
      integer ir,nf,lg,ng,i,m1
c
      ipsrc(1:it)=0.d0
      ipfis(1:it)=0.d0
c twossrc is the source excess variance, Eq. (6).
c twosfis is the fission excess variance, Eq. (8).
      twossrc=0.d0
      twosfis=0.d0
c
c fsrc1 is the source for phi_s,1
c fsrc2 is the source for phi_1
c asrc is the source for psi*_2
      fsrc1(1:neg,1:it)=0.d0
      fsrc2(1:neg,1:it)=0.d0
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
      call multmom(mrsrc,mrfis)
c
c compute I0(r), Eq. (7)
      nf=0
      do ir=1,nr
        m1=mat(ir,1)
        if(ismat(m1).eq.0)then
          nf=nf+iints(ir)
          cycle
        end if
        do i=1,iints(ir)
          nf=nf+1
          do ng=1,neg
            ipsrc(nf)=ipsrc(nf)+chisrc(ng,m1)*amom(ng,0,nf,1)
          end do ! ng
        end do ! i
      end do ! ir
c
c compute I(r), Eq. (9)
c first get chivec. if chi is already a vector, this logic keeps that.
      chivec(1:neg,1:nm)=0.d0
      do i=1,nm
        sum1=0.d0
        do ng=1,lng
          do lg=1,lng
            chivec(ng,i)=chivec(ng,i)+chi(ng,lg,i)
          end do ! lg
          sum1=sum1+chivec(ng,i)
        end do ! ng
        if(sum1.gt.0.d0)then
          chivec(1:neg,i)=chivec(1:neg,i)/sum1
        else
          chivec(1:neg,i)=0.d0
        end if
      end do ! i
      nf=0
      do ir=1,nr
        if(ifcel(ir,1).eq.0)then
          nf=nf+iints(ir)
          cycle
        end if
        m1=mat(ir,1)
        do i=1,iints(ir)
          nf=nf+1
          do ng=1,neg
            ipfis(nf)=ipfis(nf)+chivec(ng,m1)*amom(ng,0,nf,1)
          end do ! ng
        end do ! i
      end do ! ir
c
c compute 2S0 and 2S [Eqs. (6) and (8)]. note Eq. (8) has a typo,
c there should be a phi(r,E,mu) in the integrand.
      nf=0
      do ir=1,nr
        m1=mat(ir,1)
        if(mrsrc(m1).eq.0.d0.and.mrfis(m1).eq.0.d0)then
          nf=nf+iints(ir)
          cycle
        end if
        do i=1,iints(ir)
          nf=nf+1
          csrc=mrsrc(m1)*ipsrc(nf)*ipsrc(nf)*dv(nf,1)
          cfis=mrfis(m1)*ipfis(nf)*ipfis(nf)*dv(nf,1)
          do ng=1,neg
            twossrc=twossrc+csrc*nsrcf(ng,m1)
            twosfis=twosfis+cfis*nusigf(ng,m1)*fmom(ng,0,nf,1)
            fsrc1(ng,nf)=fsrc1(ng,nf)+mrsrc(m1)*nsrcf(ng,m1)
            fsrc2(ng,nf)=fsrc2(ng,nf)+mrfis(m1)*nusigf(ng,m1)
     1       *fmom(ng,0,nf,1)
            asrc(ng,nf)=mrfis(m1)*ipfis(nf)*ipfis(nf)*nusigf(ng,m1)
          end do ! ng
          fsrc1(1:neg,nf)=fsrc1(1:neg,nf)*ipsrc(nf)
          fsrc2(1:neg,nf)=fsrc2(1:neg,nf)*ipfis(nf)*chivec(1:neg,m1)
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
      write(*,'(/,"twossrc,twosfis,r2,lkg=",1p4e14.6)')twossrc,
     1 twosfis,r2,lkg
      write(iuo,'(/,"twossrc,twosfis,r2,lkg=",1p4e14.6)')twossrc,
     1 twosfis,r2,lkg
c to compare with John's NSE results
c     write(*,'("mean and variance in 4096 us")')
c     write(*,'(1p2e14.6)')lkg*4096.d-6,(yasym+1.d0)*lkg*4096.d-6
c write sm2
      write(*,'("feynman_y_asymptote, sm2",1p2e14.6)')yasym,sm2
      write(iuo,'("feynman_y_asymptote, sm2",1p2e14.6)')yasym,sm2
c
c write partisn inputs for sensitivities.
c combine the two forward sources. separate if induced fission and
c spontaneous fission sources are needed separately.
      fsrc1(1:neg,1:it)=fsrc1(1:neg,1:it)+fsrc2(1:neg,1:it)
      gams(1:neg,0,1:it,1)=fsrc1(1:neg,1:it)
      call wrdantga(1,0,3,niso0,gams)
c     gams(1:neg,0,1:it,1)=fsrc2(1:neg,1:it)
c     call wrdantga(2,0,3,niso0,gams)
      gams(1:neg,0,1:it,1)=asrc(1:neg,1:it)
      call wrdantga(3,0,3,niso0,gams)
c
c write other data needed for sensitivities.
      call wrsmbin(r2,mrsrc,mrfis,ipsrc,ipfis,chivec)
c
      return
      end
