      subroutine multmom(mom2src,mom2fis)
c compute multiplicity moments and ratios for materials.
c includes subroutine iso_2nd_mom to return the spontaneous- and
c induced-fission second moment for an input isotope.
      use COMS
      implicit none
      real*8 mom2src(neg,0:nm),mom2fis(neg,0:nm)
      real*8 mmsrc2(neg),mmfis2(neg)
      integer i,iso,iz,k2,n,ng
c
c no_sigf is set in rdsnxedt, 0 if sigf is present in the edit cross
c sections, 1 if sigf is not present.
c reset here (uncomment) if desired to use [nu(nu-1)-bar/nu-bar]
c     no_sigf=1
      if(no_sigf.eq.1)then
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
      k2=0
      do n=1,nm
        k2=k2+ncb(n-1)
        do i=1,ncb(n)
          iso=nm+k2+i ! nm is the offset to the isotope cross sections
          iz=nint(blk(1,k2+i))
          call iso_2nd_mom(iz,no_sigf,neg,mmsrc2,mmfis2)
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
      subroutine iso_2nd_mom(iz,no_sigf,neg,mmsrc2,mmfis2)
c return spontaneous- and induced-fission multiplicity second moments
c for isotope iz.
c this can be generalized to be energy-dependent.
      implicit none
      integer iz,no_sigf,neg
      real*8 mmsrc2(neg),mmfis2(neg)
c nsrcp,nfis are the number of isotopes in the source table and
c the fission table, respectively.
      integer nsrcp,nfis
      parameter (nsrcp=3,nfis=4)
      real*8 mmsrc(nsrcp,2),mmfis(nfis,2)
      integer j,n,ipnusrc(nsrcp),ipnufis(nfis)
c mmsrc data is 1st moment, 2nd moment
c 1st moment is nu-bar=sum( P(nu) * nu )
c 2nd moment is [nu(nu-1)]-bar=sum( P(nu) * nu * (nu-1) )
c P(nu) are from Boldeman and Hines; the moments were computed in a 
c spreadsheet.
c J. W. Boldeman and M. G. Hines, "Prompt Neutron Emission Probabilities
c Following Spontaneous and Thermal Neutron Fission," Nuclear Science
c and Engineering, vol. 91, 114-116 (1985).
      data ipnusrc/94240,94242,98252/
      data (mmsrc(1,n),n=1,2)/2.1563d0, 3.8242d0/ ! pu240, src
      data (mmsrc(2,n),n=1,2)/2.1450d0, 3.7942d0/ ! pu242, src
      data (mmsrc(3,n),n=1,2)/3.7570d0,11.9688d0/ ! cf252, src
      data ipnufis/92233,92235,94239,94241/
      data (mmfis(1,n),n=1,2)/2.4801d0, 4.8854d0/ ! u233, fis
      data (mmfis(2,n),n=1,2)/2.4055d0, 4.6260d0/ ! u235, fis
      data (mmfis(3,n),n=1,2)/2.8794d0, 6.7728d0/ ! pu239, fis
      data (mmfis(4,n),n=1,2)/2.9237d0, 7.0150d0/ ! pu241, fis
c
c print table entries
c     do j=1,nsrcp
c       write(*,'(i9,0p2f8.4)')ipnusrc(j),mmsrc(j,1),mmsrc(j,2)
c     end do ! i
c     do j=1,nfis
c       write(*,'(i9,0p2f8.4)')ipnufis(j),mmfis(j,1),mmfis(j,2)
c     end do ! j
c
c spontaneous-fission
      mmsrc2(1:neg)=0.d0
      do j=1,nsrcp
        if(iz.eq.ipnusrc(j))then
          mmsrc2(1:neg)=mmsrc(j,2)/mmsrc(j,1)
          exit
        end if
      end do ! j
c
c induced-fission
      mmfis2(1:neg)=0.d0
      do j=1,nfis
        if(iz.eq.ipnufis(j))then
          if(no_sigf.eq.0)then
            mmfis2(1:neg)=mmfis(j,2)
          else if(no_sigf.eq.1)then
            mmfis2(1:neg)=mmfis(j,2)/mmfis(j,1)
          end if
          exit
        end if
      end do ! j
      return
      end subroutine
