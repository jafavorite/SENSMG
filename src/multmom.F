      subroutine multmom(mom2src,mom2fis)
c compute multiplicity moments and ratios for materials.
c includes subroutine iso_2nd_mom to return the spontaneous- and
c induced-fission second moment for an input isotope.
      use COMS
      implicit none
      real*8 mom2src(neg,0:nm),mom2fis(neg,0:nm)
      logical lwrite
      real*8 mmsrc2(neg),mmfis2(neg)
      integer i,iso,iz,k2,n,ng
c
c no_sigf is set in rdsnxedt, 0 if sigf is present in the edit cross
c sections, 1 if sigf is not present.
c reset here (uncomment) if desired to use [nu(nu-1)-bar/nu-bar]
c     no_sigf=1
      if(no_sigf.eq.1.and.ictrl.eq.4)then
        write(*,'("warning. using [nu(nu-1)-bar/nu-bar] for ",
     1   "induced fission.")')
        write(iuo,'("warning. using [nu(nu-1)-bar/nu-bar] for ",
     1   "induced fission.")')
      end if
c sum of {[nu(nu-1)]-bar/nu-bar}_i*q_i, spontaneous fission
      mom2src(1:neg,0:nm)=0.d0
c sum of [nu(nu-1)]-bar]_i*sigf_i, induced fission
      mom2fis(1:neg,0:nm)=0.d0
c
      if(ictrl.eq.4)then
        lwrite=.true.
        write(iuo,'(/,"neutron multiplicity data in iso_2nd_mom")')
      else
        lwrite=.false.
      end if
      k2=0
      do n=1,nm
        k2=k2+ncb(n-1)
        do i=1,ncb(n)
          iso=nm+k2+i ! nm is the offset to the isotope cross sections
          iz=nint(blk(1,k2+i))
          call iso_2nd_mom(lwrite,iuo,iz,no_sigf,neg,ebins,mmsrc2,
     1     mmfis2)
          mom2src(1:neg,n)=mom2src(1:neg,n)+sfiso(1:neg,k2+i)
     1     *mmsrc2(1:neg)
          if(no_sigf.eq.0)then
            mom2fis(1:neg,n)=mom2fis(1:neg,n)+blk(3,k2+i)
     1       *sigf(1:neg,iso)*mmfis2(1:neg)
          else if(no_sigf.eq.1)then
c this capability is mainly for kynea3 but there may be other libraries
c with no sigf in the edits.
            mom2fis(1:neg,n)=mom2fis(1:neg,n)+blk(3,k2+i)
     1       *nusigf(1:neg,iso)*mmfis2(1:neg)
          end if
        end do ! i
      end do ! n
      lwrite=.false.
c
      if(ictrl.eq.4)then
        write(iuo,'(/,"spontaneous-fission second moment ratio * q")')
        do n=1,nm
          write(iuo,'("  material",i7)')n
          if(sum(mom2src(1:neg,n)).gt.0.d0)then
            do ng=1,neg
              write(iuo,'(i6,1pe16.8)')ng,mom2src(ng,n)
            end do ! ng
          else
            write(iuo,'("  all entries are zero.")')
          end if
        end do ! n
        write(iuo,'("induced-fission second moment * sigf")')
        do n=1,nm
          write(iuo,'("  material",i7)')n
          if(sum(mom2fis(1:neg,n)).gt.0.d0)then
            do ng=1,neg
              write(iuo,'(i6,1pe16.8)')ng,mom2fis(ng,n)
            end do ! ng
          else
            write(iuo,'("  all entries are zero.")')
          end if
        end do ! n
      end if
      return
      end
c
      subroutine iso_2nd_mom(lwrite,iuo,iz,no_sigf,neg,ebins,mmsrc2,
     1 mmfis2)
c return spontaneous- and induced-fission multiplicity second moments
c for isotope iz.
      implicit none
      logical lwrite
      integer iuo,iz,no_sigf,neg
      real*8 ebins(3,neg),mmsrc2(neg),mmfis2(neg)
c nsrcp,nfis_bh,nfis_zh are the number of isotopes in the source table,
c B&H fission table, and Z&H fission table, respectively.
      integer nsrcp,nfis_bh,nfis_zh
      parameter (nsrcp=6,nfis_bh=4,nfis_zh=5)
      real*8 mmsrc_bh(nsrcp,2),mmsrc_va(nsrcp,2),mmfis_bh(nfis_bh,2),
     1 mmfis_zh(nfis_zh,6)
      integer j,n,ng,ipnusrc(nsrcp),ipnufis_bh(nfis_bh),
     1 ipnufis_zh(nfis_zh)
      integer iva,izh
c for spontaneous fission, iva=1/0 use Valentine/use Boldeman & Hines
      iva=1
c for induced fission, izh=1/0 use energy-dependent Zucker & Holden
c                              /use Boldeman & Hines
      izh=1
c mmsrc data is 1st moment, 2nd moment
c 1st moment is nu-bar=sum( P(nu) * nu )
c 2nd moment is [nu(nu-1)]-bar=sum( P(nu) * nu * (nu-1) )
c the moments were computed in an Excel spreadsheet from P(nu)
c from B&H or Val.
c B&H:
c J. W. Boldeman and M. G. Hines, "Prompt Neutron Emission Probabilities
c Following Spontaneous and Thermal Neutron Fission," Nuclear Science
c and Engineering, vol. 91, 114-116 (1985).
c Val:
c Timothy E. Valentine, "MCNP-DSP Users Manual," Oak Ridge National
c Laboratory report ORNL/TM-13334, R2, Appendix B (January 2001);
c http://dx.doi.org/10.2172/777654. 
c use Valentine for isotopes that B&H don't include.
      data ipnusrc/94240,94242,98252,92238,96242,96244/
      data (mmsrc_bh(1,n),n=1,2)/2.1563d0, 3.8242d0/ ! pu240, src, B&H
      data (mmsrc_bh(2,n),n=1,2)/2.1450d0, 3.7942d0/ ! pu242, src, B&H
      data (mmsrc_bh(3,n),n=1,2)/3.7570d0,11.9688d0/ ! cf252, src, B&H
      data (mmsrc_bh(4,n),n=1,2)/1.9900002d0, 2.8742568d0/ ! u238, src, Val
      data (mmsrc_bh(5,n),n=1,2)/2.5400000d0, 5.1321158d0/ ! cm242, src, Val
      data (mmsrc_bh(6,n),n=1,2)/2.7200005d0, 5.9387782d0/ ! cm244, src, Val
      data (mmsrc_va(1,n),n=1,2)/2.1540000d0, 3.7889394d0/ ! pu240, src, Val
      data (mmsrc_va(2,n),n=1,2)/2.1489998d0, 3.8087442d0/ ! pu242, src, Val
      data (mmsrc_va(3,n),n=1,2)/3.77269d0, 12.04224d0/    ! cf252, src, Val
      data (mmsrc_va(4,n),n=1,2)/1.9900002d0, 2.8742568d0/ ! u238, src, Val
      data (mmsrc_va(5,n),n=1,2)/2.5400000d0, 5.1321158d0/ ! cm242, src, Val
      data (mmsrc_va(6,n),n=1,2)/2.7200005d0, 5.9387782d0/ ! cm244, src, Val
c mmfis data from B&H is also 1st moment, 2nd moment; no energy dep.
      data ipnufis_bh/92233,92235,94239,94241/
      data (mmfis_bh(1,n),n=1,2)/2.4801d0, 4.8854d0/ ! u233, fis, B&H
      data (mmfis_bh(2,n),n=1,2)/2.4055d0, 4.6260d0/ ! u235, fis, B&H
      data (mmfis_bh(3,n),n=1,2)/2.8794d0, 6.7728d0/ ! pu239, fis, B&H
      data (mmfis_bh(4,n),n=1,2)/2.9237d0, 7.0150d0/ ! pu241, fis, B&H
c mmfis data from Z&H is a quadratic fit of nu-bar and nu(nu-1)-bar with
c respect to energy. 1st 3 are for nu-bar, 2nd three for nu(nu-1)-bar.
c the fit was done in an Excel spreadsheet. energy in MeV.
c Z&H:
c M. S. Zucker and N. E Holden, "Energy Dependence of the Neutron
c Multiplicity P_nu in Fast Neutron Induced Fission of 235,238U and
c 239Pu," Brookhaven National Laboratory report BNL-38491 (1986).
c use B&H for isotopes that Z&H don't include.
      data ipnufis_zh/92235,92238,94239,92233,94241/
      data (mmfis_zh(1,n),n=1,6)/3.729018d-03, 1.150583d-01,
     1 2.396649d+00, 3.492494d-02, 4.659257d-01, 4.563657d+00/ ! u235, fis, Z&H
      data (mmfis_zh(2,n),n=1,6)/-5.827506d-10, 1.551850d-01,
     1 2.275378d+00, 2.675155d-02, 5.397033d-01, 4.334544d+00/ ! u238, fis, Z&H
      data (mmfis_zh(3,n),n=1,6)/3.683439d-04, 1.493042d-01,
     1 2.867202d+00, 2.324732d-02, 7.475031d-01, 6.688345d+00/ ! pu239, fis, Z&H
      data (mmfis_zh(4,n),n=1,6)/0.d0, 0.d0, 2.4801d0,
     1 0.d0, 0.d0, 4.8854d0/                                   ! u233, fis, B&H
      data (mmfis_zh(5,n),n=1,6)/0.d0, 0.d0, 2.9237d0,
     1 0.d0, 0.d0, 7.0150d0/                                   ! pu241, fis, B&H
c
c print everything
c     write(iuo,'(/,"neutron multiplicity data in iso_2nd_mom")')
c     if(iva.eq.1.and.neg.gt.1)then
c       write(iuo,'(2x,"spontaneous fission, Valentine")')
c       do j=1,nsrcp
c         write(iuo,'(i9,0p2f11.7)')ipnusrc(j),mmsrc_va(j,1),
c    1     mmsrc_va(j,2)
c       end do ! i
c     else
c       write(iuo,'(2x,"spontaneous fission, Boldeman & Hines")')
c       do j=1,nsrcp
c         write(iuo,'(i9,0p2f11.7)')ipnusrc(j),mmsrc_bh(j,1),
c    1     mmsrc_bh(j,2)
c       end do ! i
c       write(iuo,'(4x,"(92238, 96242, and 96244 are from Valentine)")')
c     end if
c     if(izh.eq.1.and.neg.gt.1)then
c       write(iuo,'(2x,"induced fission, Zucker & Holden")')
c       do j=1,nfis_zh
c         write(iuo,'(i9,1p6e14.6)')ipnufis_zh(j),(mmfis_zh(j,n),n=1,6)
c       end do ! j
c       write(iuo,'(4x,"(92233 and 94241 are from Boldeman & Hines)")')
c     else
c       write(iuo,'(2x,"induced fission, Boldeman & Hines")')
c       do j=1,nfis_bh
c         write(iuo,'(i9,0p2f11.7)')ipnufis_bh(j),mmfis_bh(j,1),
c    1     mmfis_bh(j,2)
c       end do ! j
c     end if
c
c spontaneous-fission
      mmsrc2(1:neg)=0.d0
c Valentine
      if(iva.eq.1.and.neg.gt.1)then
        do j=1,nsrcp
          if(iz.eq.ipnusrc(j))then
            mmsrc2(1:neg)=mmsrc_va(j,2)/mmsrc_va(j,1)
            if(lwrite)then
              write(iuo,'(i9,0p2f11.7,a)')ipnusrc(j),mmsrc_va(j,1),
     1         mmsrc_va(j,2)," spontaneous fission, Valentine"
            end if
            exit
          end if
        end do ! j
      else
c Boldeman & Hines
        do j=1,nsrcp
          if(iz.eq.ipnusrc(j))then
            mmsrc2(1:neg)=mmsrc_bh(j,2)/mmsrc_bh(j,1)
            if(lwrite)then
              write(iuo,'(i9,0p2f11.7,a)',advance='no')ipnusrc(j),
     1         mmsrc_bh(j,1),mmsrc_bh(j,2)," spontaneous fission, "
              if(iz.eq.92238.or.iz.eq.96242.or.iz.eq.96244)then
                write(iuo,'(a)')"Valentine"
              else
                write(iuo,'(a)')"Boldeman & Hines"
              end if
            end if
            exit
          end if
        end do ! j
      end if
c
c induced-fission
      mmfis2(1:neg)=0.d0
c Zucker & Holden
      if(izh.eq.1.and.neg.gt.1)then
        do j=1,nfis_zh
          if(iz.eq.ipnufis_zh(j))then
            do ng=1,neg
              mmfis2(ng)=mmfis_zh(j,4)*ebins(3,ng)**2
     1         +mmfis_zh(j,5)*ebins(3,ng)+mmfis_zh(j,6)
              if(no_sigf.eq.1)then
                mmfis2(ng)=mmfis2(ng)/(mmfis_zh(j,1)*ebins(3,ng)**2
     1           +mmfis_zh(j,2)*ebins(3,ng)+mmfis_zh(j,3))
              end if
            end do ! ng
            if(lwrite)then
              write(iuo,'(i9,1p6e14.6,a)',advance='no')ipnufis_zh(j),
     1         (mmfis_zh(j,n),n=1,6)," induced fission, "
              if(iz.eq.92233.or.iz.eq.94241)then
                write(iuo,'(a)')"Boldeman & Hines"
              else
                write(iuo,'(a)')"Zucker & Holden"
              end if
            end if
            exit
          end if
        end do ! j
      else
c Boldeman & Hines
        do j=1,nfis_bh
          if(iz.eq.ipnufis_bh(j))then
            mmfis2(1:neg)=mmfis_bh(j,2)
            if(no_sigf.eq.1)then
              mmfis2(1:neg)=mmfis2(1:neg)/mmfis_bh(j,1)
            end if
            if(lwrite)then
              write(iuo,'(i9,0p2f11.7,a)')ipnufis_bh(j),mmfis_bh(j,1),
     1         mmfis_bh(j,2)," induced fission, Boldeman & Hines"
            end if
            exit
          end if
        end do ! j
      end if
c
      return
      end subroutine
