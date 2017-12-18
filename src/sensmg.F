      program sensmg
c compute first-order sensitivities of a PARTISN calculation.
c
c this code is used with a python script, sensmg.py
c
      use COMS
      implicit none
      character ncbc*120
      integer niso0,nedt0
c
c read file "control" to see what to do.
      call rdctrl(ncbc)
      call rdmdl(ncbc)
      if(ictrl.eq.1)then
c write misc and/or sources4c input
        if(icalc.eq.0)then
c imisc=1, ialphan=1: misc and sources4c
c imisc=1, ialphan=0: misc
c imisc=0, ialphan=0 or 1: sources4c
          if(imisc.eq.1)then
            call wrmisc
          end if
          if(imisc.eq.0.or.ialphan.eq.1)then
            call wrsources
          end if
        end if
c
      elseif(ictrl.eq.2)then
c read misc and/or sources output, write partisn forward inputs
        if(icalc.eq.0)then
          if(imisc.eq.1)then
            call rdmisc
          end if
          if(imisc.eq.0.or.ialphan.eq.1)then
            call rdsources
          end if
        end if
        call wrdantnm(0,0)
        call wrdantxs
c
      elseif(ictrl.eq.3)then
c read cross sections, write partisn adjoint inputs
        if(icalc.eq.0)then
          if(imisc.eq.1)then
            call rdmisc
          end if
          if(imisc.eq.0.or.ialphan.eq.1)then
            call rdsources
          end if
        end if
        call seteff
c read forward flux moments and keff or alpha or leakage
        call rddantm(1,0)
        if(iangflux.eq.1)call rddanta(1,0)
        call rddantk
c read transport cross sections for all materials and isotopes
        call rdxsmg
c read reaction-rate edit cross sections
        call rdsnxedth(1,iue,niso0,nedt0)
        call rdxsmg0(1,niso0,nedt0)
        if(iter.eq.0)then
c write regular adjoint input and generalized adjoint input
          call wrdantnm(1,niso0)
          call gamadj0(niso0)
        else
c read regular adjoint flux moments
          call rddantm(2,0)
          if(iangflux.eq.1)call rddanta(2,0)
c check convergence of generalized adjoints
          call gamadj(niso0)
        end if
      end if ! ictrl
c
c ictrl is set to 4 in gamadj when adjoints are converged.
      if(ictrl.eq.4)then
        call calcsens
      end if
c
  800 continue
      close(iuo) ! sensmg.log
      end program
