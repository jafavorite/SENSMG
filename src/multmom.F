      subroutine multmom(mrsrc,mrfis)
c compute multiplicity moments and ratios.
      use COMS
      implicit none
      real*8 mrsrc(0:nm),mrfis(0:nm)
c nsrcp,nfis are the number of isotopes in the source table and
c the fission table, respectively.
      integer nsrcp,nfis
      parameter (nsrcp=3,nfis=4)
      real*8 mmsrc(nsrcp,2),mmfis(nfis,2),msrc(2)
      integer i,j,k2,n,ipnusrc(nsrcp),ipnufis(nfis)
c mmsrc data is 1st moment, 2nd moment
c 1st moment is nu-bar=sum( P(nu) * nu )
c 2nd moment is [nu(nu-1)]-bar=sum( P(nu) * nu * (nu-1) )
c P(nu) are from Boldeman and Hines; the moments were computed in a 
c spreadsheet.
c J. W. Boldeman and M. G. Hines, "Prompt Neutron Emission Probabilities
c Following Spontaneous and Thermal Neutron Fission," Nuclear Science and
c Engineering, vol. 91, 114-116 (1985).
      data ipnusrc/94240,94242,98252/
      data (mmsrc(1,n),n=1,2)/2.1563d0, 3.8242d0/ ! pu240, src
      data (mmsrc(2,n),n=1,2)/2.1450d0, 3.7942d0/ ! pu242, src
      data (mmsrc(3,n),n=1,2)/3.7570d0,11.9688d0/ ! cf252, src
      data ipnufis/92233,92235,94239,94241/
      data (mmfis(1,n),n=1,2)/2.4801d0,4.8854d0/ ! u233, fis
      data (mmfis(2,n),n=1,2)/2.4055d0,4.6260d0/ ! u235, fis
      data (mmfis(3,n),n=1,2)/2.8794d0,6.7728d0/ ! pu239, fis
      data (mmfis(4,n),n=1,2)/2.9237d0,7.0150d0/ ! pu241, fis
c
c print table entries
c     do i=1,nsrcp
c       write(*,'(i9,0p2f8.4)')ipnusrc(i),mmsrc(i,1),mmsrc(i,2)
c     end do ! i
c     do i=1,nfis
c       write(*,'(i9,0p2f8.4)')ipnufis(i),mmfis(i,1),mmfis(i,2)
c     end do ! i
c
c multiplicity ratios [nu(nu-1)]-bar / nu-bar
      mrsrc(0:nm)=0.d0
      mrfis(0:nm)=0.d0
c weight source isotopes by neutrons/s/cm**3. the fission
c sources would be normalized if the numerator and denominator were
c needed separately.
      k2=-ncb(1)
      do n=1,nm
        k2=k2+ncb(n)
        if(ismat(n).eq.0)cycle
        msrc(1:2)=0.d0
        do i=1,ncb(n)
          do j=1,nsrcp
            if(nint(blk(1,k2+i)).eq.ipnusrc(j))then
              msrc(1:2)=msrc(1:2)+mmsrc(j,1:2)*sum(sfiso(1:neg,k2+i))
              exit
            end if
          end do ! j
        end do ! i
        mrsrc(n)=msrc(2)/msrc(1)
      end do ! n
c TODO these values are for berp ball test problem.
c do for general mixtures.
      mrfis(1)=mmfis(3,2)/mmfis(3,1)
c
      write(iuo,'("source and fission moment ratios",0p2f8.5)')
     1 mrsrc(1),mrfis(1)
      return
      end
