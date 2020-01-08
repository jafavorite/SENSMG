      program sensmg
c compute first-order sensitivities of a PARTISN calculation.
c
c this code is used with a python script, sensmg.py
c
      use COMS
      implicit none
      integer niso0,nedt0,iwrite,iplotg,iwrxsecs,iwrsensmg
c
      niso0=0
      nedt0=0
c read file "control" to see what to do.
c (there is another call in rdmdl.)
      call rdctrl(0,iwrite,iplotg,iwrxsecs,iwrsensmg)
      call rdmdl(iwrite)
      if(iwrsensmg.eq.1)then
        call wrsensmg
        if(iplotg.eq.0.and.iwrxsecs.eq.0)then
          call stoponerror
        end if
      end if
      if(ictrl.eq.1)then
c write misc and/or sources4c input
        if(icalc.eq.0.or.icalc.eq.3.or.icalc.eq.4)then
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
        if(icalc.eq.0.or.icalc.eq.3.or.icalc.eq.4)then
          if(imisc.eq.1)then
            call rdmisc
          end if
          if(imisc.eq.0.or.ialphan.eq.1)then
            call rdsources
          end if
        end if
        call wrdantnm(0,0)
c test lnk3dnt (1 line)
c       call wrdantlnk ; call stoponerror
c write mcnp input for plotting geometry
        if(iplotg.eq.1)then
          call wrmcnp
          if(iwrxsecs.eq.0)then
            call stoponerror
          end if
        end if
        call wrdantxs
c
      elseif(ictrl.eq.3)then
c read cross sections, write partisn adjoint inputs
        if(icalc.eq.0.or.icalc.eq.3.or.icalc.eq.4)then
          if(imisc.eq.1)then
            call rdmisc
          end if
          if(imisc.eq.0.or.ialphan.eq.1)then
            call rdsources
          end if
        end if
c read forward flux moments and keff or alpha or leakage
        call rddantm(1,0)
        if(iangflux.eq.1)call rddanta(1,0)
        call rddantk
c read transport cross sections for all materials and isotopes
        call rdxsmg
c read reaction-rate edit cross sections
        call rdsnxedth(1,iuo,iue,niso0,nedt0)
        call rdxsmg0(1,niso0,nedt0)
        if(iwrxsecs.eq.1)then
          call write_xsecs
          call stoponerror
        end if
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
c
c do feyny or sm2, if called for.
c sources for sensitivity equations are computed in asymptote.
        if(icalc.eq.3.or.icalc.eq.4)then
          if(niso0.eq.0)then
            call rdsnxedth(1,iuo,iue,niso0,nedt0)
          end if
          call asymptote(niso0)
        end if
c
      else if(ictrl.eq.5)then
c for feyny or sm2. sensitivity of R1 is needed again. 
c it was written to binary.
        call calcsens_sm
      end if
c
c 2nd-order sources
      if(ictrl.eq.6)then
c read data for computing 2nd-order sources    
c read cross sections
        if(icalc.eq.0.or.icalc.eq.3.or.icalc.eq.4)then
          if(imisc.eq.1)then
            call rdmisc
          end if
          if(imisc.eq.0.or.ialphan.eq.1)then
            call rdsources
          end if
        end if
c read forward flux moments and keff or alpha or leakage
        call rddantm(1,0)
        if(iangflux.eq.1)call rddanta(1,0)
        call rddantk
c read transport cross sections for all materials and isotopes
        call rdxsmg
c read reaction-rate edit cross sections
        call rdsnxedth(1,iuo,iue,niso0,nedt0)
        call rdxsmg0(1,niso0,nedt0)
c read regular adjoint flux moments
        call rddantm(2,0)
        if(iangflux.eq.1)call rddanta(2,0)
c check convergence of generalized adjoints
        call gamadj(niso0)
c compute 2nd-order sources and write partisn forward/adjoint inputs
        call src2o(nrrr+1,niso0)
      end if
c      
c second-order sensitivities.
      if(ictrl.eq.7)then
        call calcsens2
      end if
c      
  800 continue
      close(iuo) ! sensmg.log
      end program
