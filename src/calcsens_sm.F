      subroutine calcsens_sm
c compute sensitivities for feyny or sm2. this routine computes
c relative sensitivities of the second moment, R2. relative
c sensitivities of the first moment R1 (the mean count rate)
c are stored in a binary file, senssm.
      use COMS
      implicit none
      integer ncon
      parameter (ncon=4)
      real*8 dtg(neg,nr,nz),dsg(neg,neg,0:isct,nr,nz),
     1 dfg(neg,neg,nr,nz),dqg(neg,nr,nz),
     2 sden(neg,2*nel),stot(neg,2*nel),
     3 sabs(neg,2*nel),scap(neg,2*nel),
     4 sfis(neg,2*nel),sins(neg,2*nel),
     5 selas(neg,2*nel),sinel(neg,2*nel),
     6 snu(neg,neg,2*nel),schi(neg,neg,2*nel),
     7 ssct(neg,neg,0:isct,2*nel),sscti(neg,0:isct,2*nel),
     8 sscto(neg,0:isct,2*nel),sscts(neg,0:isct,2*nel),
     9 ssctt(neg,2*nel),ssctp1(neg,neg,2*nel),
     1 st,st1,fissd,lam,rlip(0:isct),
     2 ediv(neg),mom2src(neg,0:nm),mom2fis(neg,0:nm),
     3 r2,mmsrc2(neg),mmfis2(neg),ipsrc(it),ipfis(it),
     4 ipsrc_s(0:nr),ipfis_s(0:nr),sden_m(neg,nm),
     5 nfratio,nfdenom,chifacx,chifacd,dchig(neg,neg,nr,nz),trcor
     6 ,tmp(0:5,2*nel)
      integer neg0,ndir0,nitm0,njtm0
      integer i,ir,iso,j,jr,k2,l,lg,m,n,ng,niso0,nedt0,ix,iz,
     1 izlist(nel),nelu,nel2,nsen,ienergyout,iw(2),ib(2),i_no_combine
      integer nm1 ! TODO lnk3dnt with feyny/sm2
      character sfile*20,rfile*20
c
c i_no_combine=0/1 combine as normal, to get feyny or sm2 sensitivities/
c                  do not combine, to get r2 sensitivities
      i_no_combine=0
c
c ienergyout=0/1 no division/divide by ln(E_{g})-ln(E_{g-1})
      ienergyout=0
c for debugging.
      tmp(0:5,1:2*nel)=0.d0
c
c iter can be 0 when using existing flux output files.
c     write(*,'("calcsens_sm iter,ictrl=",2i6)')iter,ictrl
      if(iter.eq.0)then
c read sources
        if(icalc.eq.0.or.icalc.eq.3.or.icalc.eq.4)then
          if(imisc.eq.1)then
            call rdmisc
          end if
          if(imisc.eq.0.or.ialphan.eq.1)then
            call rdsources
          end if
        end if
c read transport cross sections for all isotopes
        call rdxsmg
c read reaction-rate edit cross sections
        call rdsnxedth(1,iue,niso0,nedt0)
        call rdxsmg0(1,niso0,nedt0)
c read data from asymptote calculation
        call rdsmbin(r2,ipsrc,ipfis,ipsrc_s,ipfis_s)
      end if
c read pdata file written by Favorite's version of sources4c
c for this material.
      if((icalc.eq.0.or.icalc.eq.3.or.icalc.eq.4).and.
     1 (imisc.eq.0.or.ialphan.eq.1))then
        call rdsources_pdata
      end if
c
      if(icalc.eq.0.or.icalc.eq.3.or.icalc.eq.4)then
        lam=1.d0
      end if
c adjoint moments (regular and generalized) are always calculational
c moments. the sign of rlip is used to convert to inner product moments.
      do l=0,isct
        rlip(l)=real(2*l+1,8)
        if(mod(l,2).ne.0)then
          rlip(l)=-rlip(l)
        end if
c       write(*,'(0p2f8.3)')rlip(l)
      end do ! l
c
      nel2=2*nel
      sden(1:neg,1:nel2)=0.d0
      stot(1:neg,1:nel2)=0.d0
      sabs(1:neg,1:nel2)=0.d0
      scap(1:neg,1:nel2)=0.d0
      sfis(1:neg,1:nel2)=0.d0
      sins(1:neg,1:nel2)=0.d0
      selas(1:neg,1:nel2)=0.d0
      sinel(1:neg,1:nel2)=0.d0
      snu(1:neg,1:neg,1:nel2)=0.d0
      schi(1:neg,1:neg,1:nel2)=0.d0
      ssct(1:neg,1:neg,0:isct,1:nel2)=0.d0
      ssctp1(1:neg,1:neg,1:nel2)=0.d0
      ssctt(1:neg,1:nel2)=0.d0
      sscti(1:neg,0:isct,1:nel2)=0.d0
      sscto(1:neg,0:isct,1:nel2)=0.d0
      sscts(1:neg,0:isct,1:nel2)=0.d0
      sden_m(1:neg,1:nm)=0.d0
c
c compute multiplicity moment ratios for all materials.
      call multmom(mom2src,mom2fis)
c     write(*,'("mom2",1p4e17.9)')mom2src(1,0:1),mom2fis(1,0:1)
c
c maybe divide by lethargy bins.
      if(ienergyout.eq.1)then
        do ng=1,neg
          ediv(ng)=1.d0/log(ebins(2,ng)/ebins(1,ng)) 
        end do ! ng
      end if
c
      fissd=1.d0/r2
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
      write(*,'("calcsens_sm. feynman_y_asymptote, sm2, lkg, r2",
     1 1p4e17.9)')yasym,sm2,lkg,r2
      write(iuo,'("calcsens_sm. feynman_y_asymptote, sm2, lkg, r2",
     1 1p4e17.9)')yasym,sm2,lkg,r2
c if cross section sensitivities are not wanted, skip all this
      if(i_xsecs.eq.0)go to 400
c otherwise compute and write sensitivities.
c  0  psi2*,psi and Q2*,psi
c  1  psi1*,(phi1+phis,1)
cc 2  psi1*,phis,1 -- unused because the sum of these is done as n=1
c  3  Qs2*
c
      if(icalc.eq.0.or.icalc.eq.3.or.icalc.eq.4)then
        sfile="sens_s_x"
      end if
      write(*,'(/,"writing sensitivities to file ",a,".")')
     1 trim(sfile)
      write(iuo,'(/,"writing sensitivities to file ",a,".")')
     1 trim(sfile)
      open(ius,file=sfile,status='unknown')
      write(ius,'(a)')trim(id)
      write(ius,'("relative sensitivities (%/%), first order")')
      do n=0,ncon
c call integrals with 0/1 for inner products/forward only.
        if(n.eq.0)then
c 1st-moment forward
          iw(1)=1
          ib(1)=0
c 2nd-moment adjoint
          iw(2)=6
          ib(2)=3
        else if(n.eq.1)then
c <I_1**2,forward> inner product
c
c 1st-moment forward flux is already loaded and importance is
c in memory. need to dimension the importance squared such that 
c the integrals function can treat it like the angular and
c moment adjoint fluxes.
c
c importance has the same value for each group/angle and one 
c expansion moment for a given fine mesh cell, so to integrate 
c properly, only give amom a nonzero value for one group and 
c one expansion moment.
c
c don't need to redefine the angular adjoint flux because, for
c n=1, only the flux moments are used.
          amom(1:neg,1:isct,1:it,1:jt)=0.d0
          do i=1,it
            amom(1:neg,0,i,1)=ipfis(i)*ipfis(i)
          end do ! i
        else if(n.eq.2)then
c 2nd-moment forward 1+2
          iw(1)=5
          ib(1)=1
c 1st-moment adjoint
          iw(2)=2
          ib(2)=0
        else if(n.eq.3)then
c I_{1,sf}**2 volume integral
c
c similar to n.eq.1 but without a forward flux and for the 
c spontaneous fission importance. required for the sensitivity
c of the Feynman Y to density.
          fmom(1:neg,1:isct,1:it,1:jt)=0.d0
          amom(1:neg,1:isct,1:it,1:jt)=0.d0
          do i=1,it
            amom(1:neg,0,i,1)=ipsrc(i)*ipsrc(i)
            fmom(1:neg,0,i,1)=1.d0
          end do ! i
        else if(n.eq.4)then
c <I_1,m1-forward,m1-adjoint> inner product.
c required for sensitivity of R2/Y/Sm2 to chi.
c contribution of the m2-adjoint source term to this sensitivity.
c
c 1st-moment forward
          iw(1)=1
          ib(1)=0
c 1st-moment adjoint
          iw(2)=2
          ib(2)=0
        end if ! n
        if(n.eq.0.or.n.eq.2.or.n.eq.4)then
          call rddantm(iw(1),ib(1))
          if(iangflux.eq.1)call rddanta(iw(1),ib(1))
          call rddantm(iw(2),ib(2))
          if(iangflux.eq.1)call rddanta(iw(2),ib(2))
        end if
c is there an inner product that diagnoses the R2 fluxes?
c       call srcip(xfiss,amom)
c       write(iuo,'(" <psi*, Q> using amom:",1pe14.6)')xfiss
c       write(iuo,'(" <psi*, Q>/lkg - 1:",1pe14.6)')xfiss/lkg-1.d0
        if(n.eq.4)then
          call integrals_3f(lam,fmom,amom,ipfis,dchig)
        else
          call integrals(0,lam,rlip,fmom,amom,afreg,afadj,
     1     neg0,ndir0,nitm0,njtm0,dtg,dsg,dfg,dqg)
        end if ! n.eq.4
c
        if(n.eq.2)then
c no dqg because this is only for the forward flux.
          dtg(1:neg,1:nr,1:nz)=2.d0*dtg(1:neg,1:nr,1:nz)
          dsg(1:neg,1:neg,0:isct,1:nr,1:nz)=2.d0
     1     *dsg(1:neg,1:neg,0:isct,1:nr,1:nz)
          dfg(1:neg,1:neg,1:nr,1:nz)=2.d0*dfg(1:neg,1:neg,1:nr,1:nz)
        end if ! n.eq.2
c
c cesaro transport correction
        if(itrcor.eq.3)then
          do l=1,isct
            trcor=real((isct+2)*(isct+1),8)
            trcor=real((isct+2-l)*(isct+1-l),8)/trcor
            dsg(1:neg,1:neg,l,1:nr,1:nz)=trcor
     1       *dsg(1:neg,1:neg,l,1:nr,1:nz)
          end do ! l
        end if
c debug
c       do ng=1,neg
c         write(10,'("dtg",i6,1pe15.3)')ng,dtg(ng,1,1)
c       end do ! ng
c       do ng=1,neg
c         do lg=1,neg
c           write(10,'("dsg",2i6,1pe15.3)')ng,lg,dsg(ng,lg,0,1,1)
c         end do ! lg
c       end do ! ng
c       do ng=1,neg
c         do lg=1,neg
c           write(10,'("dfg",2i6,1pe15.3)')ng,lg,dfg(ng,lg,1,1)
c         end do ! lg
c       end do ! ng
c       close(10)
c       call stoponerror
c       do ng=1,neg
c         write(10,'("sigc",i6,1p200e15.6)')ng,(sigc(ng,i),i=1,nxs)
c       end do ! ng
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
          if(fissdata.eq.2.and.ilnk3dnt.eq.0)then
            call nusigf_denom(i,nfdenom)
          end if
          do j=1,ncb(i)
            k2=k2+1
            iso=nm+k2 ! nm is the offset to the isotope cross sections
            iz=nint(blk(1,k2))
            call iso_2nd_mom(.false.,iuo,iz,no_sigf,neg,ebins,mmsrc2,
     1       mmfis2)
c           write(*,'(i14,1p2e17.9)')iz,mmsrc2(1),mmfis2(1)
c           write(*,'("mom2",1p4e17.9)')mom2src(1,0:1),mom2fis(1,0:1)
            if(fissdata.eq.2)then
              call nusigf_ratio(iso,k2,nfdenom,nfratio)
            end if
            do ir=1,nr
              do jr=1,nz
                if(mat(ir,jr).ne.i)cycle
                do ng=1,neg
                  if(n.eq.0.or.n.eq.2)then
c contribution from total loss interaction
                    stot(ng,k2)=stot(ng,k2)-dtg(ng,ir,jr)*sigt(ng,iso)
     1               *blk(3,k2)
                    sabs(ng,k2)=sabs(ng,k2)-dtg(ng,ir,jr)*siga(ng,iso)
     1               *blk(3,k2)
                    scap(ng,k2)=scap(ng,k2)-dtg(ng,ir,jr)*sigc(ng,iso)
     1               *blk(3,k2)
                    sfis(ng,k2)=sfis(ng,k2)-dtg(ng,ir,jr)*sigf(ng,iso)
     1               *blk(3,k2)
                    selas(ng,k2)=selas(ng,k2)-dtg(ng,ir,jr)
     1               *sige(ng,iso)*blk(3,k2)
                    sinel(ng,k2)=sinel(ng,k2)-dtg(ng,ir,jr)
     1               *sigi(ng,iso)*blk(3,k2)
c                   tmp(n+1,k2)=tmp(n+1,k2)-dtg(ng,ir,jr)*sigf(ng,iso)
c    1               *blk(3,k2)
c contribution from fission source term.
                    do lg=1,neg
                      st=dfg(lg,ng,ir,jr)*chi(lg,ng,iso)*nusigf(ng,iso)
     1                 *blk(3,k2)
                      snu(lg,ng,k2)=snu(lg,ng,k2)+st
                      sfis(ng,k2)=sfis(ng,k2)+st
                      if(chieff.eq.1)then
c get relative sensitivity of chi(ng to lg) to nu(ng) and sigf(ng)
c (they are equal) and to isotope density
                        call chideriv(i,iso,k2,lg,ng,nfdenom,chifacx,
     1                   chifacd)
                        st1=st*chifacx
                        snu(lg,ng,k2)=snu(lg,ng,k2)+st1
                        sfis(ng,k2)=sfis(ng,k2)+st1
c subtract st1 because later sden will use the sum of snu and the
c derivative in st1 shouldn't be included.
                        sden(ng,k2)=sden(ng,k2)+st*chifacd-st1
                      end if
                      if(fissdata.eq.2)then
                        st=dfg(lg,ng,ir,jr)*chi(lg,ng,iso)*nusigf(ng,i)
     1                   *nfratio
                      end if
                      schi(lg,ng,k2)=schi(lg,ng,k2)+st
                      tmp(n,k2)=tmp(n,k2)+st
c contribution from scatter source term.
                      do l=0,isct
                        st=dsg(ng,lg,l,ir,jr)*sigs(ng,lg,l,iso)
     1                   *blk(3,k2)
                        ssctt(lg,k2)=ssctt(lg,k2)+st
c sigs contributes to sigt for l=0 and outgoing groups.
                        if(l.eq.0)then
                          st=st-dtg(lg,ir,jr)*sigs(ng,lg,l,iso)
     1                     *blk(3,k2)
                        end if
                        ssct(ng,lg,l,k2)=ssct(ng,lg,l,k2)+st
                        if(lg.ne.ng)then
                          sscti(ng,l,k2)=sscti(ng,l,k2)+st
                          sscto(lg,l,k2)=sscto(lg,l,k2)+st
                        else
                          sscts(ng,l,k2)=sscts(ng,l,k2)+st
                        end if
c diagonal and bhs transport corrections
                        if((itrcor.eq.1.and.lg.eq.ng).or.itrcor.eq.2)
     1                   then
                          ssctp1(ng,lg,k2)=ssctp1(ng,lg,k2)
     1                     -dsg(ng,lg,l,ir,jr)*sigsp1(ng,lg,iso)
     2                     *blk(3,k2)
                          if(l.eq.0)then
                            ssctp1(ng,lg,k2)=ssctp1(ng,lg,k2)
     1                       +dtg(lg,ir,jr)*sigsp1(ng,lg,iso)*blk(3,k2)
                          end if
                        end if
                      end do ! l
c sige(ng)+sigi(ng).le.sum(sigs(1:neg,ng,0,nm+n))
c selas+sinel should be .le. sum(sscto+sscts), but that can't happen
c using just the edits. not sure if these sensitivities are meaningful.
c this has the same form as out-scattering.
                      if(sigs(ng,lg,0,iso).ne.0.d0)then
                        selas(lg,k2)=selas(lg,k2)+dsg(ng,lg,0,ir,jr)
     1                   *sige(lg,iso)*blk(3,k2)
                        sinel(lg,k2)=sinel(lg,k2)+dsg(ng,lg,0,ir,jr)
     1                   *sigi(lg,iso)*blk(3,k2)
                      end if
                    end do ! lg
c contribution of the inhomogeneous source. the spectrum is an
c approximation but the total is correct.
c sfiso is the absolute spont. fiss. source rate density.
c saniso is the (alpha,n) spectrum for k2 times dQ/dN.
                    if(n.eq.0)then
                      if(ismat(i).eq.1)then
c ignore (alpha,n)?
                        sins(ng,k2)=sins(ng,k2)+dqg(ng,ir,jr)
     1                   *(sfiso(ng,k2)+saniso(ng,k2)*blk(3,k2))
                      end if
                    end if
c contribution from second-moment adjoint source.
c only incoming groups are used for dfg. this looks wrong but it is correct.
                  else if(n.eq.1)then
                    st=dfg(1,ng,ir,jr)*mmfis2(ng)*sigf(ng,iso)*blk(3,k2)
                    sfis(ng,k2)=sfis(ng,k2)+st
                    sden(ng,k2)=sden(ng,k2)+st
c                   tmp(5,k2)=tmp(5,k2)+dfg(ng,ng,ir,jr)*mmfis2(ng)
c    1               *sigf(ng,iso)*blk(3,k2)
                  else if(n.eq.3)then
                    sins(ng,k2)=sins(ng,k2)+dfg(1,ng,ir,jr)*mmsrc2(ng)
     1               *sfiso(ng,k2)
                  else if(n.eq.4)then
                    do lg=1,neg
c sum over forward energy groups.
                      if(fissdata.eq.0.or.fissdata.eq.1)then
c TODO check this
                        st=2.d0*dchig(lg,ng,ir,jr)*mom2fis(ng,i)
     1                   *chi(lg,ng,iso)
                      else if(fissdata.eq.2)then
                        st=2.d0*dchig(lg,ng,ir,jr)*mom2fis(ng,i)
     1                   *chi(lg,ng,iso)*nfratio
                      end if
c                     if(ng.eq.1.and.lg.eq.1)then
c                       write(*,'("i,k2",2i4)')i,k2
c                       write(*,'("st",1pe17.9)')st
c                       write(*,'("sum dfg",1pe17.9)')
c    1                   sum(dfg(1:neg,1:neg,1:nr,1:nz))
c                       write(*,'("mom2fis",1pe17.9)')
c    1                   sum(mom2fis(1:neg,i))
c                       write(*,'("chi",1pe17.9)')sum(chi(lg,1:neg,iso))
c                       write(*,'("dchig",1pe17.9)')dchig(lg,ng,1,1)
c                       write(*,'("nfratio",1pe17.9)')nfratio
c                     end if
                      schi(lg,ng,k2)=schi(lg,ng,k2)+st
                      tmp(n,k2)=tmp(n,k2)+st
                    end do ! lg
                  end if ! n
                end do ! ng
              end do ! jr
            end do ! ir
          end do ! j
        end do ! i
        do k2=1,nel
          do ng=1,neg
            if(n.eq.3)then
c the terms already have the correct sign.
              sden(ng,k2)=sden(ng,k2)+sins(ng,k2)
            else if(n.eq.2)then
              sden(ng,k2)=sden(ng,k2)+sum(snu(1:neg,ng,k2))
     1         +stot(ng,k2)+ssctt(ng,k2)
            end if
c debug print
c           write(*,'(2i4,1p5e14.6)')n,k2,sum(sden(1:neg,k2)),
c    1       sum(stot(1:neg,k2)),sum(ssctt(1:neg,k2)),
c    2       sum(snu(1:neg,1:neg,k2)),sum(sins(1:neg,k2))
          end do ! ng
        end do ! k2
c
      end do ! n
c     write(*,'("tmp",1p20e17.9)')tmp(0:4,1)
c     write(*,'("tmp",1p20e17.9)')tmp(0:4,2)
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
        do j=i,nel
          if(nint(blk(1,j)).ne.iz)cycle
          sden(1:neg,k2)=sden(1:neg,k2)+sden(1:neg,j)
          stot(1:neg,k2)=stot(1:neg,k2)+stot(1:neg,j)
          sabs(1:neg,k2)=sabs(1:neg,k2)+sabs(1:neg,j)
          scap(1:neg,k2)=scap(1:neg,k2)+scap(1:neg,j)
          sfis(1:neg,k2)=sfis(1:neg,k2)+sfis(1:neg,j)
          sins(1:neg,k2)=sins(1:neg,k2)+sins(1:neg,j)
          selas(1:neg,k2)=selas(1:neg,k2)+selas(1:neg,j)
          sinel(1:neg,k2)=sinel(1:neg,k2)+sinel(1:neg,j)
          snu(1:neg,1:neg,k2)=snu(1:neg,1:neg,k2)
     1     +snu(1:neg,1:neg,j)
          schi(1:neg,1:neg,k2)=schi(1:neg,1:neg,k2)
     1     +schi(1:neg,1:neg,j)
          ssct(1:neg,1:neg,0:isct,k2)=ssct(1:neg,1:neg,0:isct,k2)
     1     +ssct(1:neg,1:neg,0:isct,j)
          ssctp1(1:neg,1:neg,k2)=ssctp1(1:neg,1:neg,k2)
     1     +ssctp1(1:neg,1:neg,j)
          ssctt(1:neg,k2)=ssctt(1:neg,k2)+ssctt(1:neg,j)
          sscti(1:neg,0:isct,k2)=sscti(1:neg,0:isct,k2)
     1     +sscti(1:neg,0:isct,j)
          sscto(1:neg,0:isct,k2)=sscto(1:neg,0:isct,k2)
     1     +sscto(1:neg,0:isct,j)
          sscts(1:neg,0:isct,k2)=sscts(1:neg,0:isct,k2)
     1     +sscts(1:neg,0:isct,j)
        end do ! j
  100   continue
      end do ! i
c there is a block of sensitivities for all the isotopes (1:nel)
c and another block for sums over unique isotopes (nel+1:nelu)
      nsen=nel+nelu
c
c for R2 sensitivities, divide by R2
      sden(1:neg,1:nsen)=sden(1:neg,1:nsen)*fissd
      stot(1:neg,1:nsen)=stot(1:neg,1:nsen)*fissd
      sabs(1:neg,1:nsen)=sabs(1:neg,1:nsen)*fissd
      scap(1:neg,1:nsen)=scap(1:neg,1:nsen)*fissd
      sfis(1:neg,1:nsen)=sfis(1:neg,1:nsen)*fissd
      sins(1:neg,1:nsen)=sins(1:neg,1:nsen)*fissd
      selas(1:neg,1:nsen)=selas(1:neg,1:nsen)*fissd
      sinel(1:neg,1:nsen)=sinel(1:neg,1:nsen)*fissd
      snu(1:neg,1:neg,1:nsen)=snu(1:neg,1:neg,1:nsen)*fissd
      schi(1:neg,1:neg,1:nsen)=schi(1:neg,1:neg,1:nsen)*fissd
      ssct(1:neg,1:neg,0:isct,1:nsen)=
     1 ssct(1:neg,1:neg,0:isct,1:nsen)*fissd
      ssctp1(1:neg,1:neg,1:nsen)=
     1 ssctp1(1:neg,1:neg,1:nsen)*fissd
      sscti(1:neg,0:isct,1:nsen)=
     1 sscti(1:neg,0:isct,1:nsen)*fissd
      sscto(1:neg,0:isct,1:nsen)=
     1 sscto(1:neg,0:isct,1:nsen)*fissd
      sscts(1:neg,0:isct,1:nsen)=
     1 sscts(1:neg,0:isct,1:nsen)*fissd
      ssctt(1:neg,1:nsen)=ssctt(1:neg,1:nsen)*fissd
c
c normalize chi as requested.
c ichinorm = 0/1/2 unnormalized/full norm/partial norm
      if(ichinorm.eq.1.or.ichinorm.eq.2)then
        call normchi(nsen,izlist,schi(1:neg,1:neg,1:nsen))
      end if
c
c combine with sensitivities of R1 to get feyny or sm2.
      if(i_no_combine.eq.0)then
        call smcombine1(sden,stot,sabs,scap,schi,snu,sfis,
     1   ssct,sscti,sscts,sscto,sins,selas,sinel,ssctt,ssctp1)
      end if
c the sensitivities are now relative sensitivities for feyny or sm2.
c
c output
c indirect effect, by isotope and material.
      if(i_no_combine.eq.0)then
        if(icalc.eq.3)then
          write(ius,'(/," feyny sensitivities",10x,"feyny=",
     1     1pe15.7)')yasym
        else if(icalc.eq.4)then
          write(ius,'(/," sm2 sensitivities",10x,"sm2=",
     1     1pe15.7)')sm2
        end if
      else if(i_no_combine.eq.1)then
        write(ius,'(/," r2 sensitivities",10x,"r2=",
     1   1pe15.7)')r2
      end if
c k2 and m are used to keep track of materials for writing the header.
      k2=0
      m=1
      do i=1,nel
        k2=k2+1
        write(ius,'(/,1x,a24," isotope",i4," in material",i7)')
     1   zaidfull(i),k2,m
        call output1(sden(1:neg,i),stot(1:neg,i),
     1   sabs(1:neg,i),scap(1:neg,i),
     2   schi(1:neg,1:neg,i),snu(1:neg,1:neg,i),
     3   sfis(1:neg,i),ssct(1:neg,1:neg,0:isct,i),
     4   sscti(1:neg,0:isct,i),sscts(1:neg,0:isct,i),
     5   sscto(1:neg,0:isct,i),sins(1:neg,i),
     6   selas(1:neg,i),sinel(1:neg,i),ssctt(1:neg,i),
     7   ssctp1(1:neg,1:neg,i))
c material density
        sden_m(1:neg,m)=sden_m(1:neg,m)+sden(1:neg,i)
        if(k2.eq.ncb(m))then
          k2=0
          m=m+1
        end if
      end do ! i
c
c isotope totals including indirect+direct
c for keff or alpha (n.eq.0), this is not needed if nelu.eq.nel
c (it is done anyway).
      write(ius,'(/" isotope totals")')
c
      do i=1,nelu
c the zaid for this isotope is nint(blk(1,izlist(i)))
        ix=nel+i
        write(ius,'(/,1x,a24," isotope",i4)')zaidfull(izlist(i)),i
        call output1(sden(1:neg,ix),stot(1:neg,ix),
     1   sabs(1:neg,ix),scap(1:neg,ix),
     2   schi(1:neg,1:neg,ix),snu(1:neg,1:neg,ix),
     3   sfis(1:neg,ix),ssct(1:neg,1:neg,0:isct,ix),
     4   sscti(1:neg,0:isct,ix),sscts(1:neg,0:isct,ix),
     5   sscto(1:neg,0:isct,ix),sins(1:neg,ix),
     6   selas(1:neg,ix),sinel(1:neg,ix),ssctt(1:neg,ix),
     7   ssctp1(1:neg,1:neg,ix))
      end do ! i
c
      close(ius)
c
c skip to here if cross section sensitivities are not wanted.
  400 rfile="sens_s_r"
      write(*,'(/,"writing sensitivities to file ",a,".")')trim(rfile)
      write(iuo,'(/,"writing sensitivities to file ",a,".")')trim(rfile)
c temporary (4 lines)
      write(*,'("warning. results for interfaces in ",a,
     1 " are preliminary.")')trim(rfile)
      write(iuo,'("warning. results for interfaces in ",a,
     1 " are preliminary.")')trim(rfile)
      open(iur,file=rfile,status='unknown')
      write(iur,'(a)')trim(id)
c temporary (1 line)
      write(iur,'("warning. results for interfaces in this file ",
     1 "are preliminary.")')
      write(iur,'("relative sensitivities (%/%) to densities, ",
     1 "first order (constant volume)")')
      write(iur,'("first derivatives with respect to ",
     1 "interface locations (constant density)")')
c output material density sensitivities, or compute them again.
c there needs to be a switch for this.
      nm1=nm
c either this
c     if(icalc.eq.3)then
c       write(iur,'(/," feyny sensitivities",10x,"feyny=",
c    1   1pe15.7,/)')yasym
c     else if(icalc.eq.4)then
c       write(iur,'(/," sm2 sensitivities",10x,"sm2=",1pe15.7,/)')sm2
c     end if
c     call output2(nm1,sden_m(1:neg,1:nm1))
c material densities are already combined and output. this is only here
c to open sensrx.
c     if(i_no_combine.eq.0)then
c       call smcombine2(nm1,sden_m(1:neg,1:nm1))
c     end if
c or this
      call calcsens_sm_d(ediv,fissd,lam,r2,ipsrc,ipfis,mom2src,mom2fis,
     1 rlip,rfile,nm1,i_no_combine)
c end
      call calcsens_sm_r(ediv,fissd,lam,r2,ipsrc,ipfis,ipsrc_s,ipfis_s,
     1 mom2src,mom2fis,rlip,rfile,i_no_combine)
      close(iur)
c
      return
      end subroutine
c
      subroutine smcombine1(sden,stot,sabs,scap,schi,snu,sfis,
     1 ssct,sscti,sscts,sscto,sins,selas,sinel,ssctt,ssctp1)
c combine r1 and r2 sensitivities. input variables are r2
c sensitivities. read r1 from a binary file.
      use COMS
      implicit none
      real*8 sden(neg,2*nel),stot(neg,2*nel),sabs(neg,2*nel),
     1 scap(neg,2*nel),sfis(neg,2*nel),schi(neg,neg,2*nel),
     2 snu(neg,neg,2*nel),ssct(neg,neg,0:isct,2*nel),
     3 sscti(neg,0:isct,2*nel),sscto(neg,0:isct,2*nel),
     4 sscts(neg,0:isct,2*nel),sins(neg,2*nel),selas(neg,2*nel),
     5 sinel(neg,2*nel),ssctt(neg,2*nel),ssctp1(neg,neg,2*nel)
      real*8 tt
      real*8 sden_tmp(neg,2*nel),stot_tmp(neg,2*nel),
     1 sabs_tmp(neg,2*nel),
     2 scap_tmp(neg,2*nel),sfis_tmp(neg,2*nel),schi_tmp(neg,neg,2*nel),
     3 snu_tmp(neg,neg,2*nel),ssct_tmp(neg,neg,0:isct,2*nel),
     4 sscti_tmp(neg,0:isct,2*nel),sscto_tmp(neg,0:isct,2*nel),
     5 sscts_tmp(neg,0:isct,2*nel),sins_tmp(neg,2*nel),
     6 selas_tmp(neg,2*nel),sinel_tmp(neg,2*nel),ssctt_tmp(neg,2*nel),
     7 ssctp1_tmp(neg,neg,2*nel)
      integer nel2
c
      call rdoutput1_bin(sden_tmp,stot_tmp,sabs_tmp,scap_tmp,
     1 schi_tmp,snu_tmp,sfis_tmp,
     2 ssct_tmp,sscti_tmp,sscts_tmp,sscto_tmp,sins_tmp,selas_tmp,
     3 sinel_tmp,ssctt_tmp,ssctp1_tmp)
c
      if(icalc.eq.3)then
       tt=1.d0 ! feyny
      else if(icalc.eq.4)then
       tt=2.d0 ! sm2
      end if
      nel2=2*nel
      sden(1:neg,1:nel2)=sden(1:neg,1:nel2)-tt*sden_tmp(1:neg,1:nel2)
      stot(1:neg,1:nel2)=stot(1:neg,1:nel2)-tt*stot_tmp(1:neg,1:nel2)
      sabs(1:neg,1:nel2)=sabs(1:neg,1:nel2)-tt*sabs_tmp(1:neg,1:nel2)
      scap(1:neg,1:nel2)=scap(1:neg,1:nel2)-tt*scap_tmp(1:neg,1:nel2)
      sfis(1:neg,1:nel2)=sfis(1:neg,1:nel2)-tt*sfis_tmp(1:neg,1:nel2)
      sins(1:neg,1:nel2)=sins(1:neg,1:nel2)-tt*sins_tmp(1:neg,1:nel2)
      selas(1:neg,1:nel2)=selas(1:neg,1:nel2)-tt*selas_tmp(1:neg,1:nel2)
      sinel(1:neg,1:nel2)=sinel(1:neg,1:nel2)-tt*sinel_tmp(1:neg,1:nel2)
      schi(1:neg,1:neg,1:nel2)=schi(1:neg,1:neg,1:nel2)
     1 -tt*schi_tmp(1:neg,1:neg,1:nel2)
      snu(1:neg,1:neg,1:nel2)=snu(1:neg,1:neg,1:nel2)
     1 -tt*snu_tmp(1:neg,1:neg,1:nel2)
      ssct(1:neg,1:neg,0:isct,1:nel2)=ssct(1:neg,1:neg,0:isct,1:nel2)
     1 -tt*ssct_tmp(1:neg,1:neg,0:isct,1:nel2)
      ssctp1(1:neg,1:neg,1:nel2)=ssctp1(1:neg,1:neg,1:nel2)
     1 -tt*ssctp1_tmp(1:neg,1:neg,1:nel2)
      sscti(1:neg,0:isct,1:nel2)=sscti(1:neg,0:isct,1:nel2)
     1 -tt*sscti_tmp(1:neg,0:isct,1:nel2)
      sscto(1:neg,0:isct,1:nel2)=sscto(1:neg,0:isct,1:nel2)
     1 -tt*sscto_tmp(1:neg,0:isct,1:nel2)
      sscts(1:neg,0:isct,1:nel2)=sscts(1:neg,0:isct,1:nel2)
     1 -tt*sscts_tmp(1:neg,0:isct,1:nel2)
      ssctt(1:neg,1:nel2)=ssctt(1:neg,1:nel2)-tt*ssctt_tmp(1:neg,1:nel2)
      return
      end subroutine
