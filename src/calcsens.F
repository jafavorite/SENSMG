      subroutine calcsens
c compute sensitivities.
      use COMS
      implicit none
      real*8 dtg(neg,nr,nz),dsg(neg,neg,0:isct,nr,nz),
     1 dfg(neg,neg,nr,nz),dqg(neg,nr,nz),
     2 sden(0:nrrr,neg,2*nel),stot(0:nrrr,neg,2*nel),
     3 sabs(0:nrrr,neg,2*nel),scap(0:nrrr,neg,2*nel),
     4 sfis(0:nrrr,neg,2*nel),sins(0:nrrr,neg,2*nel),
     5 snu(0:nrrr,neg,neg,2*nel),schi(0:nrrr,neg,neg,2*nel),
     6 ssct(0:nrrr,neg,neg,0:isct,2*nel),sscti(0:nrrr,neg,0:isct,2*nel),
     7 sscto(0:nrrr,neg,0:isct,2*nel),sscts(0:nrrr,neg,0:isct,2*nel),
     8 ssctt(0:nrrr,neg,2*nel),st,fissd,lam,rlip(0:isct),
     9 rxnrate(neg,nrrx),rxnratet(nrrx+nflux),
     1 sdir(nrrr,neg,2),ediv(neg),xfiss,
     2 sensverify(neg+1,neg+1,0:isct,nel),tmp(5)
      integer neg0,ndir0,nitm0,njtm0
      integer iconv,i,ir,iso,j,jr,k2,l,lg,m,n,ng,niso0,nedt0,ix,
     1 iz,izlist(nel),nelu,nsen,ienergyout,icalc_org
      character sfile*20,rfile*20
c temp. use more moments.
c     integer isct1,nmom1
c
c feyny and sm2, icalc=3 and icalc=4, will go through all of this. set
c icalc=0 but remember what it was.
      icalc_org=icalc
      if(icalc.eq.3.or.icalc.eq.4)then
        icalc=0
      end if
c
c ienergyout=0/1 no division/divide by ln(E_{g})-ln(E_{g-1})
      ienergyout=0
c for debugging.
      tmp(1:5)=0.d0
c
c iter can be 0 when using existing flux output files.
      if(iter.eq.0)then
c read sources
        if(icalc.eq.0)then
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
c read transport cross sections for all isotopes
        call rdxsmg
c read reaction-rate edit cross sections
        call rdsnxedth(1,iue,niso0,nedt0)
        call rdxsmg0(1,niso0,nedt0)
c read regular adjoint flux moments
        call rddantm(2,0)
        if(iangflux.eq.1)call rddanta(2,0)
c in case unconverged adjoints are used. treat them as converged.
        open(iut,file='stopconverged',status='unknown')
        close (iut)
      end if
c read pdata file written by Favorite's version of sources4c
c for this material.
      if(icalc.eq.0.and.(imisc.eq.0.or.ialphan.eq.1))then
        call rdsources_pdata
      end if
c
      if(icalc.eq.0.or.icalc.eq.2)then
        lam=1.d0
      else if(icalc.eq.1)then
        lam=1.d0/keff
      end if
c adjoint moments (regular and generalized) are always calculational
c moments. the sign of rlip is used to convert to inner product moments.
      do l=0,isct
        rlip(l)=real(2*l+1)
        if(mod(l,2).ne.0)then
          rlip(l)=-rlip(l)
        end if
c       write(*,'(0p2f8.3)')rlip(l)
      end do ! l
c
      sden(0:nrrr,1:neg,1:2*nel)=0.d0
      stot(0:nrrr,1:neg,1:2*nel)=0.d0
      sabs(0:nrrr,1:neg,1:2*nel)=0.d0
      scap(0:nrrr,1:neg,1:2*nel)=0.d0
      sfis(0:nrrr,1:neg,1:2*nel)=0.d0
      sins(0:nrrr,1:neg,1:2*nel)=0.d0
      snu(0:nrrr,1:neg,1:neg,1:2*nel)=0.d0
      schi(0:nrrr,1:neg,1:neg,1:2*nel)=0.d0
      ssct(0:nrrr,1:neg,1:neg,0:isct,1:2*nel)=0.d0
      ssctt(0:nrrr,1:neg,1:2*nel)=0.d0
      sscti(0:nrrr,1:neg,0:isct,1:2*nel)=0.d0
      sscto(0:nrrr,1:neg,0:isct,1:2*nel)=0.d0
      sscts(0:nrrr,1:neg,0:isct,1:2*nel)=0.d0
      sdir(1:nrrr,1:neg,1:2)=0.d0
c
c maybe divide by lethargy bins.
      if(ienergyout.eq.1)then
        do ng=1,neg
          ediv(ng)=1.d0/log(ebins(2,ng)/ebins(1,ng)) 
        end do ! ng
      end if
c
c edit reaction rates.
      if(nrrx.gt.0)then
        call rxnrateedit(rxnrate,rxnratet)
      end if
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
c get reaction rates to compare with partisn balance tables
c for testing.
c     call integrals(1,lam,rlip,fmom,amom,afreg,afadj,neg0,
c    1 ndir0,nitm0,njtm0,dtg,dsg,dfg,dqg)
c
c compute and write sensitivities.
c n=0 is for keff, alpha, or leakage sensitivities,
c n > 0 is for reaction-rate sensitivities.
      do n=0,nrrr
        if(n.eq.0)then
          if(icalc.eq.0)then
            sfile="sens_l_x"
          else if(icalc.eq.1)then
            sfile="sens_k_x"
          else if(icalc.eq.2)then
            sfile="sens_a_x"
          end if
        else if(n.gt.0)then
          sfile="sens_rr_x"
        end if
        if(n.le.1.and.ictrl.eq.4)then
          write(*,'(/,"writing sensitivities to file ",a,".")')
     1     trim(sfile)
          write(iuo,'(/,"writing sensitivities to file ",a,".")')
     1     trim(sfile)
          open(ius,file=sfile,status='unknown')
          write(ius,'(a)')trim(id)
          write(ius,'("relative sensitivities (%/%), first order")')
        end if
c
        if(n.eq.0)then
          if(icalc.eq.0)then
            call srcip(xfiss,amom)
            if(ictrl.eq.4)then
              write(iuo,'(" <psi*, Q> using amom:",1pe14.6)')xfiss
              write(iuo,'(" <psi*, Q>/lkg - 1:",1pe14.6)')xfiss/lkg-1.d0
            end if
            fissd=1.d0/lkg
          else if(icalc.eq.1)then
            call fissip(fissd,amom)
            if(ictrl.eq.4)then
              write(iuo,'(" <psi*, F psi> using fmom and amom:",
     1         1pe14.6)')fissd
            end if
            fissd=1.d0/fissd
          else if(icalc.eq.2)then
            call velip(fissd,rlip,amom,afadj,neg0,ndir0,nitm0,njtm0)
            if(ictrl.eq.4)then
              write(iuo,'(" <psi*, 1/vel psi> using fmom and amom:",
     1         1pe14.6)')fissd
            end if
            fissd=1.d0/(fissd*alpha)
          end if
        else
c set this to pass in to other routines; is not used for n.ne.0
          fissd=1.d0
          if(icalc.eq.0)then
            call rddantm(3,n)
            if(iangflux.eq.1)then
              call rddanta(3,n)
            end if
            call srcip(xfiss,gmom)
            if(ictrl.eq.4)then
              write(iuo,'(" <psi*, Q> using gmom:",1pe14.6)')xfiss
              if(irri(2,n).le.nrrx)then
                write(iuo,'(" compare with 0.")')
              else
                write(iuo,'(" <psi*, Q> - 1:",1pe14.6)')xfiss-1.d0
              end if
            end if
          else if(icalc.eq.1.or.icalc.eq.2)then
            if(iangflux.eq.0)then
              call rdgmom(n,iconv)
            else if(iangflux.eq.1)then
              call rdafgad(n,iconv)
c construct gmom
              call calcmom(afgad,gmom)
            end if
c debug print
c           if(n.eq.1.and.ictrl.eq.4)then
c             write(iuo,'("debug print, iteration",i3)')iter
c             l=0
c             do i=1,it
c               write(101,'(i6,1p100e14.6)')i,0.5d0*(rfm(i)+rfm(i-1)),
c    1           gmom(1:neg,l,i,1),sum(gmom(1:neg,l,i,1))
c             end do ! i
c           end if
            if(iconv.ne.1.and.ictrl.eq.4)then
              write(*,'("warning. adjoint",i3," is not converged. ",
     1         "continuing.")')n
              write(iuo,'(/,"warning. adjoint",i3," is not converged. ",
     1         "continuing.")')n
              write(ius,'(/,"warning. adjoint",i3," is not converged. ",
     1         "continuing.")')n
            end if
          end if
        end if
c
c call integrals with 0/1 for inner products/forward only.
        if(n.eq.0)then
          call integrals(0,lam,rlip,fmom,amom,afreg,afadj,
     1     neg0,ndir0,nitm0,njtm0,dtg,dsg,dfg,dqg)
c temp. use more moments. comment call integrals (preceding line).
c         isct1=17 !isct ! test that my moments are correct.
c         write(*,'("for testing. calling temp_integrals_mom with ",
c    1     "isct1=",i3)')isct1
c         write(iuo,'("for testing. calling temp_integrals_mom with ",
c    1     "isct1=",i3)')isct1
c         if(calc2d.eq.0.or.calc2d.eq.-1)then
c           nmom1=isct1+1
c         elseif(calc2d.eq.1)then
c           nmom1=(isct1+1)*(isct1+2)/2 ! for iquad=5 and 6
c         end if
c         call temp_integrals_mom(0,isct1,nmom1,lam,afreg,afadj,neg0,
c    1     ndir0,nitm0,njtm0,dtg,dsg,dfg,dqg)
c debug
c         do ng=1,neg
c           write(10,'("dtg",i6,1pe15.3)')ng,dtg(ng,1,1)
c         end do ! ng
c         do ng=1,neg
c           do lg=1,neg
c             write(10,'("dsg",2i6,1pe15.3)')ng,lg,dsg(ng,lg,0,1,1)
c           end do ! lg
c         end do ! ng
c         do ng=1,neg
c           do lg=1,neg
c             write(10,'("dfg",2i6,1pe15.3)')ng,lg,dfg(ng,lg,1,1)
c           end do ! lg
c         end do ! ng
c         close(10)
c         call stoponerror
c         do ng=1,neg
c           write(10,'("sigc",i6,1p200e15.6)')ng,(sigc(ng,i),i=1,nxs)
c         end do ! ng
        else
          call integrals(0,lam,rlip,fmom,gmom,afreg,afgad,
     1     neg0,ndir0,nitm0,njtm0,dtg,dsg,dfg,dqg)
        end if
c
c cycle through transport isotopes for the indirect effect.
c nm (number of materials) is the offset into the isotope part
c of the cross section arrays. blk(3,k2) is atoms/b-cm.
c this logic catches a material used in different regions.
c dtg = <G*, psi> by energy, coarse mesh
c dsg = <G*, (2l+1)psi> by energy out, energy in, moment, coarse mesh
c dfg = <G*, psi/keff> by energy out, energy in, coarse mesh
c dqg = <G*> by energy, coarse mesh
        k2=0
        do i=1,nm
          do j=1,ncb(i)
            k2=k2+1
            iso=nm+k2 ! nm is the offset to the isotope cross sections
            do ir=1,nr
              do jr=1,nz
                if(mat(ir,jr).ne.i)cycle
                do ng=1,neg
                  stot(n,ng,k2)=stot(n,ng,k2)-dtg(ng,ir,jr)*sigt(ng,iso)
     1             *blk(3,k2)
                  sabs(n,ng,k2)=sabs(n,ng,k2)-dtg(ng,ir,jr)*siga(ng,iso)
     1             *blk(3,k2)
                  scap(n,ng,k2)=scap(n,ng,k2)-dtg(ng,ir,jr)*sigc(ng,iso)
     1             *blk(3,k2)
                  sfis(n,ng,k2)=sfis(n,ng,k2)-dtg(ng,ir,jr)*sigf(ng,iso)
     1             *blk(3,k2)
c contribution of the inhomogeneous source. the spectrum is an
c approximation but the total is correct.
c sfiso is the absolute spont. fiss. source rate density.
c saniso is the (alpha,n) spectrum for k2 times dQ/dN.
                  if(icalc.eq.0.and.ismat(i).eq.1)then
                    sins(n,ng,k2)=sins(n,ng,k2)+dqg(ng,ir,jr)
     1               *(sfiso(ng,k2)+saniso(ng,k2)*blk(3,k2))
#if ( defined RPSD18 )
                    if(nint(blk(1,k2)).eq.94242)then
                      tmp(1)=tmp(1)+dqg(ng,ir,jr)*sfiso(ng,k2)
                      tmp(2)=tmp(2)+dqg(ng,ir,jr)*saniso(ng,k2)
     1                 *blk(3,k2)
                    end if
#endif
                  end if
                  do lg=1,neg
                    st=dfg(lg,ng,ir,jr)*chi(lg,ng,iso)*nusigf(ng,iso)
     1               *blk(3,k2)
                    snu(n,lg,ng,k2)=snu(n,lg,ng,k2)+st
                    sfis(n,ng,k2)=sfis(n,ng,k2)+st
                    do l=0,isct
                      st=dsg(ng,lg,l,ir,jr)*sigs(ng,lg,l,iso)*blk(3,k2)
                      ssctt(n,lg,k2)=ssctt(n,lg,k2)+st
c sigs contributes to sigt for l=0 and outgoing groups.
                      if(l.eq.0)then
                        st=st-dtg(lg,ir,jr)*sigs(ng,lg,l,iso)*blk(3,k2)
                      end if
                      ssct(n,ng,lg,l,k2)=ssct(n,ng,lg,l,k2)+st
                      if(lg.ne.ng)then
                        sscti(n,ng,l,k2)=sscti(n,ng,l,k2)+st
                        sscto(n,lg,l,k2)=sscto(n,lg,l,k2)+st
                      else
                        sscts(n,ng,l,k2)=sscts(n,ng,l,k2)+st
                      end if
                    end do ! l
                  end do ! lg
                end do ! ng
              end do ! jr
            end do ! ir
          end do ! j
        end do ! i
        do k2=1,nel
          do ng=1,neg
c the terms already have the correct sign.
c sins is 0 except for fixed-source problems.
            sden(n,ng,k2)=sins(n,ng,k2)+sum(snu(n,1:neg,ng,k2))
     1       +stot(n,ng,k2)+ssctt(n,ng,k2)
          end do ! ng
c debug print
c         write(*,'(2i4,1p5e14.6)')n,k2,sum(sden(n,1:neg,k2)),
c    1     sum(stot(n,1:neg,k2)),sum(ssctt(n,1:neg,k2)),
c    2     sum(snu(n,1:neg,1:neg,k2)),sum(sins(n,1:neg,k2))
        end do ! k2
c
c calculate the direct effect.
c this assumes that the reaction is not a transport cross section.
c therefore the totals (direct+indirect) may not be correct.
c check this.
        if(n.gt.0)then
          j=irri(1,n)
          sdir(n,1:neg,1)=rxnrate(1:neg,j)/rxnratet(j)
          j=irri(2,n)
          if(j.le.nrrx)then
            sdir(n,1:neg,2)=-rxnrate(1:neg,j)/rxnratet(j)
          end if
        end if
c
c compute totals by isotope
        nelu=0 ! number of unique isotopes, not counting edits only
        do i=1,nel
          iz=nint(blk(1,i))
          if(i.eq.1)then
            izlist(i)=1
            nelu=nelu+1
          else
            do j=1,nelu
              if(iz.eq.nint(blk(1,izlist(j))))go to 100
            end do ! j 
            nelu=nelu+1
            izlist(nelu)=i
          end if
          k2=nel+nelu ! isotope totals are offset by nel
c add direct effect first. density and total get it.
c then check reactions.
          if(n.gt.0)then
            do j=1,2
              if(iz.eq.irrx(1,irri(j,n)))then
c if chi is allowed as a reaction, test against that for sden here.
                sden(n,1:neg,k2)=sden(n,1:neg,k2)+sdir(n,1:neg,j)
c this test is because 2 is nusigf which may be nu.
c comment this test if nusigf is sigf.
c               if(irrx(2,irri(j,n)).ne.2)then
                  stot(n,1:neg,k2)=stot(n,1:neg,k2)+sdir(n,1:neg,j)
c               end if
                if(irrx(2,irri(j,n)).eq.2)then
c if nusigf is considered sigf, use this.
                  sfis(n,1:neg,k2)=sfis(n,1:neg,k2)+sdir(n,1:neg,j)
c if nusigf is considered nu, figure out what to do with this:
c                 snu(n,1:neg,1:neg,k2)=snu(n,1:neg,1:neg,k2)
c    1            +sdir(n,1:neg,j)
                else if(irrx(2,irri(j,n)).eq.4)then
                  sabs(n,1:neg,k2)=sabs(n,1:neg,k2)+sdir(n,1:neg,j)
                else if(irrx(2,irri(j,n)).eq.10)then
                  scap(n,1:neg,k2)=scap(n,1:neg,k2)+sdir(n,1:neg,j)
                else if(irrx(2,irri(j,n)).eq.16)then
                  sfis(n,1:neg,k2)=sfis(n,1:neg,k2)+sdir(n,1:neg,j)
                else if(irrx(2,irri(j,n)).ne.3)then
                  write(*,'("*** error. edit reaction is not in the ",
     1             "list for direct effect in calcsens. rxn=",i5)')
     2             irrx(2,irri(j,n))
                  write(iuo,'("*** error. edit reaction is not in the ",
     1             "list for direct effect in calcsens. rxn=",i5)')
     2             irrx(2,irri(j,n))
                end if
              end if
            end do ! j
          end if
          do j=i,nel
            if(nint(blk(1,j)).ne.iz)cycle
            sden(n,1:neg,k2)=sden(n,1:neg,k2)+sden(n,1:neg,j)
            stot(n,1:neg,k2)=stot(n,1:neg,k2)+stot(n,1:neg,j)
            sabs(n,1:neg,k2)=sabs(n,1:neg,k2)+sabs(n,1:neg,j)
            scap(n,1:neg,k2)=scap(n,1:neg,k2)+scap(n,1:neg,j)
            sfis(n,1:neg,k2)=sfis(n,1:neg,k2)+sfis(n,1:neg,j)
            sins(n,1:neg,k2)=sins(n,1:neg,k2)+sins(n,1:neg,j)
            snu(n,1:neg,1:neg,k2)=snu(n,1:neg,1:neg,k2)
     1       +snu(n,1:neg,1:neg,j)
            ssct(n,1:neg,1:neg,0:isct,k2)=ssct(n,1:neg,1:neg,0:isct,k2)
     1       +ssct(n,1:neg,1:neg,0:isct,j)
            ssctt(n,1:neg,k2)=ssctt(n,1:neg,k2)+ssctt(n,1:neg,j)
            sscti(n,1:neg,0:isct,k2)=sscti(n,1:neg,0:isct,k2)
     1       +sscti(n,1:neg,0:isct,j)
            sscto(n,1:neg,0:isct,k2)=sscto(n,1:neg,0:isct,k2)
     1       +sscto(n,1:neg,0:isct,j)
            sscts(n,1:neg,0:isct,k2)=sscts(n,1:neg,0:isct,k2)
     1       +sscts(n,1:neg,0:isct,j)
          end do ! j
  100     continue
        end do ! i
c there is a block of sensitivities for all the isotopes (1:nel)
c and another block for sums over unique isotopes (nel+1:nelu)
        nsen=nel+nelu
c
c for leakage sensitivities, divide by leakage
c for keff sensitivities, divide by <psi*, F psi>
c for alpha sensitivities, divide by <psi*, 1/vel psi>*alpha
        if(n.eq.0)then
          sden(n,1:neg,1:nsen)=sden(n,1:neg,1:nsen)*fissd
          stot(n,1:neg,1:nsen)=stot(n,1:neg,1:nsen)*fissd
          sabs(n,1:neg,1:nsen)=sabs(n,1:neg,1:nsen)*fissd
          scap(n,1:neg,1:nsen)=scap(n,1:neg,1:nsen)*fissd
          sfis(n,1:neg,1:nsen)=sfis(n,1:neg,1:nsen)*fissd
          sins(n,1:neg,1:nsen)=sins(n,1:neg,1:nsen)*fissd
          snu(n,1:neg,1:neg,1:nsen)=snu(n,1:neg,1:neg,1:nsen)*fissd
          ssct(n,1:neg,1:neg,0:isct,1:nsen)=
     1     ssct(n,1:neg,1:neg,0:isct,1:nsen)*fissd
          sscti(n,1:neg,0:isct,1:nsen)=
     1     sscti(n,1:neg,0:isct,1:nsen)*fissd
          sscto(n,1:neg,0:isct,1:nsen)=
     1     sscto(n,1:neg,0:isct,1:nsen)*fissd
          sscts(n,1:neg,0:isct,1:nsen)=
     1     sscts(n,1:neg,0:isct,1:nsen)*fissd
          ssctt(n,1:neg,1:nsen)=ssctt(n,1:neg,1:nsen)*fissd
#if ( defined RPSD18 )
          tmp(1:2)=tmp(1:2)*fissd
#endif
        end if
c
c normalize chi as requested.
c ichinorm = 0/1/2 unnormalized/full norm/partial norm
c the sum of the elements of the chi matrix (fissdata=0 or 1)
c is neg, so the normalization is 1/neg, not 1.
        if(ichinorm.eq.0)then
          schi(n,1:neg,1:neg,1:nsen)=snu(n,1:neg,1:neg,1:nsen)
        else if(ichinorm.eq.1.or.ichinorm.eq.2)then
          do i=1,nsen
            st=sum(snu(n,1:neg,1:neg,i))
            if(st.eq.0.d0)then
              cycle
            end if
            if(i.le.nel)then
              iso=nm+i
            else
              iso=nm+izlist(i-nel)
            end if
            do ng=1,neg
c full normalization.
              if(fissdata.eq.2)then ! chi vector
                schi(n,ng,ng,i)=sum(snu(n,ng,1:neg,i))
     1           -st*chi(ng,ng,iso)
              else
                st=sum(snu(n,1:neg,ng,i))
                schi(n,1:neg,ng,i)=snu(n,1:neg,ng,i)
     1           -st*chi(1:neg,ng,iso)
              end if
c one more division for partial normalization.
c needs verification.
              if(ichinorm.eq.2)then
                if(fissdata.eq.2)then ! chi vector
                  schi(n,ng,1:neg,i)=schi(n,ng,1:neg,i)
     1             /(1.d0-chi(ng,ng,iso))
                else
                  schi(n,1:neg,ng,i)=schi(n,1:neg,ng,i)
     1             /(1.d0-chi(1:neg,ng,iso))
                end if
              end if
            end do ! ng
          end do ! i
        end if
c
c output
c indirect effect, by isotope and material.
        if(n.eq.0)then
          if(icalc.eq.0)then
            write(ius,'(/," leakage sensitivities",10x,"leakage=",
     1       1pe15.7)')lkg
          else if(icalc.eq.1)then
            write(ius,'(/," keff sensitivities",10x,"keff=",0pf12.8)')
     1       keff
          else if(icalc.eq.2)then
            write(ius,'(/," alpha sensitivities",10x,"alpha=",
     1       1pe15.7)')alpha
          end if
        else
          call outputh(ius,n,irrx(1:2,irri(1,n)),irrx(1:2,irri(2,n)),
     1     icalc,rxnratio(n))
          write(ius,'(" indirect effect")')
        end if
c k2 and m are used to keep track of materials for writing the header.
        k2=0
        m=1
        do i=1,nel
          k2=k2+1
          write(ius,'(/,1x,a24," isotope",i4," in material",i4)')
     1     zaidfull(i),k2,m
#if ( defined RPSD18 )
c test
c         write(*,'(2i4,1pe14.6)')n,i,sum(sins(n,1:neg,i))
          if(nint(blk(1,k2)).eq.94242)then
            write(ius,'(" source: fiss,an=",1p2e14.6)')tmp(1),tmp(2)
          end if
#endif
          call output1(sden(n,1:neg,i),stot(n,1:neg,i),
     1     sabs(n,1:neg,i),scap(n,1:neg,i),
     2     schi(n,1:neg,1:neg,i),snu(n,1:neg,1:neg,i),
     3     sfis(n,1:neg,i),ssct(n,1:neg,1:neg,0:isct,i),
     4     sscti(n,1:neg,0:isct,i),sscts(n,1:neg,0:isct,i),
     5     sscto(n,1:neg,0:isct,i),sins(n,1:neg,i),
     6     ssctt(n,1:neg,i))
          if(k2.eq.ncb(m))then
            k2=0
            m=m+1
          end if
        end do ! i
c
c direct effect (for reaction-rate ratios only)
        if(n.gt.0)then
          call outputh(ius,n,irrx(1:2,irri(1,n)),irrx(1:2,irri(2,n)),
     1     icalc,rxnratio(n))
          write(ius,'(" direct effect")')
          write(ius,'(/,1x,a24," isotope in numerator")')
     1     zaidfull(nel+irri(1,n))
          if(irri(2,n).le.nrrx)then
            write(ius,'(1x,a24," isotope in denominator")')
     1       zaidfull(nel+irri(2,n))
          else
            write(ius,'(1x,"none",20x," isotope in denominator")')
          end if
c summary table first. these should be 1. and -1.
          write(ius,'(" summary (sums over groups)")')
          write(ius,'(a8,3x,1pe16.6)')cedits(irrx(2,irri(1,n))),
     1     sum(sdir(n,1:neg,1))
          write(ius,'(a8,3x,1pe16.6)')cedits(irrx(2,irri(2,n))),
     1     sum(sdir(n,1:neg,2))
c now the details
          write(ius,'(1x,"group",2x,"e_lower",7x,"e_upper",7x,
     1     "e_average",5x,a8,6x,a8)')cedits(irrx(2,irri(1,n))),
     2     cedits(irrx(2,irri(2,n)))
          do ng=1,neg
            write(ius,'(i6,1p5e14.6)')ng,ebins(1:3,ng),
     1       sdir(n,ng,1),sdir(n,ng,2)
          end do ! ng
c 42x is 14*3 for the energy groups.
          write(ius,'(" total",42x,1p2e14.6)')sum(sdir(n,1:neg,1)),
     1     sum(sdir(n,1:neg,2))
        end if
c
c isotope totals including indirect+direct
c for keff or alpha (n.eq.0), this is not needed if nelu.eq.nel
c (it is done anyway).
        if(n.eq.0)then
          write(ius,'(/" isotope totals")')
        else
          call outputh(ius,n,irrx(1:2,irri(1,n)),irrx(1:2,irri(2,n)),
     1     icalc,rxnratio(n))
          write(ius,'(" isotope totals, indirect + direct")')
        end if
c
        do i=1,nelu
c the zaid for this isotope is nint(blk(1,izlist(i)))
          ix=nel+i
          write(ius,'(/,1x,a24," isotope",i4)')zaidfull(izlist(i)),i
          call output1(sden(n,1:neg,ix),stot(n,1:neg,ix),
     1     sabs(n,1:neg,ix),scap(n,1:neg,ix),
     2     schi(n,1:neg,1:neg,ix),snu(n,1:neg,1:neg,ix),
     3     sfis(n,1:neg,ix),ssct(n,1:neg,1:neg,0:isct,ix),
     4     sscti(n,1:neg,0:isct,ix),sscts(n,1:neg,0:isct,ix),
     5     sscto(n,1:neg,0:isct,ix),sins(n,1:neg,ix),
     6     ssctt(n,1:neg,ix))
        end do ! i
c
        if(n.eq.0.or.n.eq.nrrr)then
          close(ius)
        end if
c
c for feyny or sm2, write sensitivities to binary file.
        if(icalc_org.eq.3.or.icalc_org.eq.4)then
          call output1_bin(sden(n,1:neg,1:nel),stot(n,1:neg,1:nel),
     1     sabs(n,1:neg,1:nel),scap(n,1:neg,1:nel),
     2     schi(n,1:neg,1:neg,1:nel),snu(n,1:neg,1:neg,1:nel),
     3     sfis(n,1:neg,1:nel),ssct(n,1:neg,1:neg,0:isct,1:nel),
     4     sscti(n,1:neg,0:isct,1:nel),sscts(n,1:neg,0:isct,1:nel),
     5     sscto(n,1:neg,0:isct,1:nel),sins(n,1:neg,1:nel),
     6     ssctt(n,1:neg,1:nel))
        end if
c
c density and geometry sensitivities.
        if(n.eq.0)then
          if(icalc.eq.0)then
            rfile="sens_l_r"
          else if(icalc.eq.1)then
            rfile="sens_k_r"
          else if(icalc.eq.2)then
            rfile="sens_a_r"
          end if
        else if(n.gt.0)then
          rfile="sens_rr_r"
        end if
        if(n.le.1)then
          open(iur,file=rfile,status='unknown')
          write(iur,'(a)')trim(id)
          write(iur,'("relative sensitivities (%/%) to densities, ",
     1     "first order (constant volume)")')
          write(iur,'("first derivatives with respect to ",
     1     "interface locations (constant density)")')
        end if
c TODO make sure this works for rxn rates; remove this if
        if(n.eq.0)then
          call calcsens_d(n,ediv,fissd,dtg,dsg,dfg,dqg,rfile,
     1     rxnrate,rxnratet)
          call calcsens_r(n,ediv,fissd,lam,rlip,rfile,
     1     rxnrate,rxnratet)
        else
          write(iur,'(/," reaction-rate ratio capability to be added")')
        end if
        if(n.eq.0.or.n.eq.nrrr)then
          close(iur)
        end if
c
      end do ! n
c
      if(nrrr.eq.0.and.ictrl.eq.4)then
        write(*,'("no sens_rr file because there are no reaction ",
     1   "rates specified.")')
        write(iuo,'("no sens_rr file because there are no reaction ",
     1   "rates specified.")')
      end if
c
c debug using direct pert
c     call wrmacrsx
      if(iver.eq.1.and.ictrl.eq.4)then
        sensverify(1:neg+1,1:neg+1,0:isct,1:nel)=0.d0
c set to the sens that we want to verify and pass to wrodninp. then the
c cross section perts can be set using (desired pert)/sensverify.
        n=0 ! 0/>0 keff or alpha/rxn-rate ratio
        lg=1
        l=0
        do i=1,nel
c for scat
c         do l=0,isct
            do ng=1,neg
c for tot
c             sensverify(ng,lg,l,i)=stot(n,ng,i)
c for dens
              sensverify(ng,lg,l,i)=sden(n,ng,i)
c for fiss
c             sensverify(ng,lg,l,i)=sfis(n,ng,i)
c for nu
c             sensverify(ng,lg,l,i)=sum(snu(n,1:neg,ng,i))
c for cap
c             sensverify(ng,lg,l,i)=scap(n,ng,i)
c for scat
c             do lg=1,neg
c               sensverify(ng,lg,l,i)=ssct(n,ng,lg,l,i)
c               sensverify(neg+1,lg,l,i)=sum(ssct(n,1:neg,lg,l,i))
c             end do ! lg
            end do ! ng
c for tot
c           sensverify(neg+1,lg,l,i)=sum(stot(n,1:neg,i))
c for dens
            sensverify(neg+1,lg,l,i)=sum(sden(n,1:neg,i))
c for fiss
c           sensverify(neg+1,lg,l,i)=sum(sfis(n,1:neg,i))
c for nu
c           sensverify(neg+1,lg,l,i)=sum(snu(n,1:neg,1:neg,i))
c for cap
c           sensverify(neg+1,lg,l,i)=sum(scap(n,1:neg,i))
c for scat
c         end do ! l
        end do ! i
c write odninp file for verification.
        call wrodninp(sensverify)
      else if(iver.eq.2)then
        call wrodninp2
      end if
c
c reset icalc.
      icalc=icalc_org
c
      return
      end subroutine
