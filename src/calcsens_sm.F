      subroutine calcsens_sm
c compute sensitivities for feyny or sm2. this routine computes
c relative sensitivities of the second moment, R2. relative
c sensitivities of the first moment R1 (the mean count rate)
c are stored in a binary file, senssm.
      use COMS
      implicit none
      integer ncon
      parameter (ncon=3)
      real*8 dtg(neg,nr,nz),dsg(neg,neg,0:isct,nr,nz),
     1 dfg(neg,neg,nr,nz),dqg(neg,nr,nz),
     2 sden(neg,2*nel),stot(neg,2*nel),
     3 sabs(neg,2*nel),scap(neg,2*nel),
     4 sfis(neg,2*nel),sins(neg,2*nel),
     5 snu(neg,neg,2*nel),schi(neg,neg,2*nel),
     6 ssct(neg,neg,0:isct,2*nel),sscti(neg,0:isct,2*nel),
     7 sscto(neg,0:isct,2*nel),sscts(neg,0:isct,2*nel),
     8 ssctt(neg,2*nel),st,fissd,lam,rlip(0:isct),
     1 ediv(neg),
     2 sensverify(neg+1,neg+1,0:isct,nel),tmp(5),
     3 r2,mrsrc(0:nm),mrfis(0:nm),ipsrc(it),ipfis(it),
     4 chivec(neg,nm)
      integer neg0,ndir0,nitm0,njtm0
      integer i,ir,iso,j,jr,k2,l,lg,m,n,ng,niso0,nedt0,ix,iz,
     1 izlist(nel),nelu,ns,nsen,ienergyout,icalc_org,iw(2),ib(2)
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
      write(*,'("calcsens_sm",2i6)')iter,ictrl
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
c read transport cross sections for all isotopes
        call rdxsmg
c read reaction-rate edit cross sections
        call rdsnxedth(1,iue,niso0,nedt0)
        call rdxsmg0(1,niso0,nedt0)
c read data from asymptote calculation
        call rdsmbin(r2,mrsrc,mrfis,ipsrc,ipfis,chivec)
      end if
c read pdata file written by Favorite's version of sources4c
c for this material.
      if(icalc.eq.0.and.(imisc.eq.0.or.ialphan.eq.1))then
        call rdsources_pdata
      end if
c
      if(icalc.eq.0)then
        lam=1.d0
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
      ns=2*nel
      sden(1:neg,1:ns)=0.d0
      stot(1:neg,1:ns)=0.d0
      sabs(1:neg,1:ns)=0.d0
      scap(1:neg,1:ns)=0.d0
      sfis(1:neg,1:ns)=0.d0
      sins(1:neg,1:ns)=0.d0
      snu(1:neg,1:neg,1:ns)=0.d0
      schi(1:neg,1:neg,1:ns)=0.d0
      ssct(1:neg,1:neg,0:isct,1:ns)=0.d0
      ssctt(1:neg,1:ns)=0.d0
      sscti(1:neg,0:isct,1:ns)=0.d0
      sscto(1:neg,0:isct,1:ns)=0.d0
      sscts(1:neg,0:isct,1:ns)=0.d0
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
c compute and write sensitivities.
c TODO n indicates these combinations:
c  0  psi2*,psi and Q2*,psi
c  1  psi1*,(phi1+phis,1)
cc 2  psi1*,phis,1
c  3  Qs2*
c
      write(*,'("calcsens_sm. feynman_y_asymptote, sm2, lkg, r2",
     1 1p4e14.6)')yasym,sm2,lkg,r2
      write(iuo,'("calcsens_sm. feynman_y_asymptote, sm2, lkg, r2",
     1 1p4e14.6)')yasym,sm2,lkg,r2
      do n=0,ncon
        if(n.eq.0)then
          if(icalc.eq.0)then
            sfile="sens_s_x"
          end if
          write(*,'(/,"writing sensitivities to file ",a,".")')
     1     trim(sfile)
          write(iuo,'(/,"writing sensitivities to file ",a,".")')
     1     trim(sfile)
          open(ius,file=sfile,status='unknown')
          write(ius,'(a)')trim(id)
          write(ius,'("relative sensitivities (%/%), first order")')
        end if
c
c call integrals with 0/1 for inner products/forward only.
        if(n.eq.0)then
c 1st-moment forward
          iw(1)=1
          ib(1)=0
c 2nd-moment adjoint
          iw(2)=6
          ib(2)=3
        else if(n.eq.1)then
c 2nd-moment forward 1
          iw(1)=5
          ib(1)=1
c 1st-moment adjoint
          iw(2)=2
          ib(2)=0
        else if(n.eq.2)then
          cycle ! for now, use the sum of the two forwards.
c 2nd-moment forward 2
          iw(1)=1
          ib(1)=0
c 1st-moment adjoint
          iw(2)=6
          ib(2)=3
        else if(n.eq.3)then
c 1st-moment forward
          iw(1)=1
          ib(1)=0
c 1st-moment adjoint
          iw(2)=2
          ib(2)=0
        end if
        call rddantm(iw(1),ib(1))
        if(iangflux.eq.1)call rddanta(iw(1),ib(1))
        call rddantm(iw(2),ib(2))
        if(iangflux.eq.1)call rddanta(iw(2),ib(2))
c is there an inner product that diagnoses the R2 fluxes?
c         call srcip(xfiss,amom)
c         write(iuo,'(" <psi*, Q> using amom:",1pe14.6)')xfiss
c         write(iuo,'(" <psi*, Q>/lkg - 1:",1pe14.6)')xfiss/lkg-1.d0
        call integrals(0,lam,rlip,fmom,amom,afreg,afadj,
     1   neg0,ndir0,nitm0,njtm0,dtg,dsg,dfg,dqg)
        if(n.eq.1.or.n.eq.2)then
c no dqg because this is only for the forward flux.
          dtg(1:neg,1:nr,1:nz)=2.d0*dtg(1:neg,1:nr,1:nz)
          dsg(1:neg,1:neg,0:isct,1:nr,1:nz)=2.d0*
     1     dsg(1:neg,1:neg,0:isct,1:nr,1:nz)
          dfg(1:neg,1:neg,1:nr,1:nz)=2.d0*dfg(1:neg,1:neg,1:nr,1:nz)
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
          do j=1,ncb(i)
            k2=k2+1
            iso=nm+k2 ! nm is the offset to the isotope cross sections
            do ir=1,nr
              do jr=1,nz
                if(mat(ir,jr).ne.i)cycle
                do ng=1,neg
                  stot(ng,k2)=stot(ng,k2)-dtg(ng,ir,jr)*sigt(ng,iso)
     1             *blk(3,k2)
                  sabs(ng,k2)=sabs(ng,k2)-dtg(ng,ir,jr)*siga(ng,iso)
     1             *blk(3,k2)
                  scap(ng,k2)=scap(ng,k2)-dtg(ng,ir,jr)*sigc(ng,iso)
     1             *blk(3,k2)
                  sfis(ng,k2)=sfis(ng,k2)-dtg(ng,ir,jr)*sigf(ng,iso)
     1             *blk(3,k2)
c contribution of the inhomogeneous source. the spectrum is an
c approximation but the total is correct.
c sfiso is the absolute spont. fiss. source rate density.
c saniso is the (alpha,n) spectrum for k2 times dQ/dN.
                  if(icalc.eq.0.and.ismat(i).eq.1)then
                    sins(ng,k2)=sins(ng,k2)+dqg(ng,ir,jr)
     1               *(sfiso(ng,k2)+saniso(ng,k2)*blk(3,k2))
                  end if
                  do lg=1,neg
                    st=dfg(lg,ng,ir,jr)*chi(lg,ng,iso)*nusigf(ng,iso)
     1               *blk(3,k2)
                    snu(lg,ng,k2)=snu(lg,ng,k2)+st
                    sfis(ng,k2)=sfis(ng,k2)+st
                    do l=0,isct
                      st=dsg(ng,lg,l,ir,jr)*sigs(ng,lg,l,iso)*blk(3,k2)
                      ssctt(lg,k2)=ssctt(lg,k2)+st
c sigs contributes to sigt for l=0 and outgoing groups.
                      if(l.eq.0)then
                        st=st-dtg(lg,ir,jr)*sigs(ng,lg,l,iso)*blk(3,k2)
                      end if
                      ssct(ng,lg,l,k2)=ssct(ng,lg,l,k2)+st
                      if(lg.ne.ng)then
                        sscti(ng,l,k2)=sscti(ng,l,k2)+st
                        sscto(lg,l,k2)=sscto(lg,l,k2)+st
                      else
                        sscts(ng,l,k2)=sscts(ng,l,k2)+st
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
            sden(ng,k2)=sins(ng,k2)+sum(snu(1:neg,ng,k2))
     1       +stot(ng,k2)+ssctt(ng,k2)
          end do ! ng
c debug print
c         write(*,'(2i4,1p5e14.6)')n,k2,sum(sden(1:neg,k2)),
c    1     sum(stot(1:neg,k2)),sum(ssctt(1:neg,k2)),
c    2     sum(snu(1:neg,1:neg,k2)),sum(sins(1:neg,k2))
        end do ! k2
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
            snu(1:neg,1:neg,k2)=snu(1:neg,1:neg,k2)
     1       +snu(1:neg,1:neg,j)
            ssct(1:neg,1:neg,0:isct,k2)=ssct(1:neg,1:neg,0:isct,k2)
     1       +ssct(1:neg,1:neg,0:isct,j)
            ssctt(1:neg,k2)=ssctt(1:neg,k2)+ssctt(1:neg,j)
            sscti(1:neg,0:isct,k2)=sscti(1:neg,0:isct,k2)
     1       +sscti(1:neg,0:isct,j)
            sscto(1:neg,0:isct,k2)=sscto(1:neg,0:isct,k2)
     1       +sscto(1:neg,0:isct,j)
            sscts(1:neg,0:isct,k2)=sscts(1:neg,0:isct,k2)
     1       +sscts(1:neg,0:isct,j)
          end do ! j
  100     continue
        end do ! i
c there is a block of sensitivities for all the isotopes (1:nel)
c and another block for sums over unique isotopes (nel+1:nelu)
        nsen=nel+nelu
c
      end do ! n
c
c for R2 sensitivities, divide by R2
      sden(1:neg,1:nsen)=sden(1:neg,1:nsen)*fissd
      stot(1:neg,1:nsen)=stot(1:neg,1:nsen)*fissd
      sabs(1:neg,1:nsen)=sabs(1:neg,1:nsen)*fissd
      scap(1:neg,1:nsen)=scap(1:neg,1:nsen)*fissd
      sfis(1:neg,1:nsen)=sfis(1:neg,1:nsen)*fissd
      sins(1:neg,1:nsen)=sins(1:neg,1:nsen)*fissd
      snu(1:neg,1:neg,1:nsen)=snu(1:neg,1:neg,1:nsen)*fissd
      ssct(1:neg,1:neg,0:isct,1:nsen)=
     1 ssct(1:neg,1:neg,0:isct,1:nsen)*fissd
      sscti(1:neg,0:isct,1:nsen)=
     1 sscti(1:neg,0:isct,1:nsen)*fissd
      sscto(1:neg,0:isct,1:nsen)=
     1 sscto(1:neg,0:isct,1:nsen)*fissd
      sscts(1:neg,0:isct,1:nsen)=
     1 sscts(1:neg,0:isct,1:nsen)*fissd
      ssctt(1:neg,1:nsen)=ssctt(1:neg,1:nsen)*fissd
c
c combine with sensitivities of R1.
      call smcombine(sden,stot,sabs,scap,schi,snu,sfis,
     1 ssct,sscti,sscts,sscto,sins,ssctt)
c the sensitivities are now relative sensitivities for feyny or sm2.
c
c normalize chi as requested.
c ichinorm = 0/1/2 unnormalized/full norm/partial norm
c the sum of the elements of the chi matrix (fissdata=0 or 1)
c is neg, so the normalization is 1/neg, not 1.
      if(ichinorm.eq.0)then
        schi(1:neg,1:neg,1:nsen)=snu(1:neg,1:neg,1:nsen)
      else if(ichinorm.eq.1.or.ichinorm.eq.2)then
        do i=1,nsen
          st=sum(snu(1:neg,1:neg,i))
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
              schi(ng,ng,i)=sum(snu(ng,1:neg,i))
     1         -st*chi(ng,ng,iso)
            else
              st=sum(snu(1:neg,ng,i))
              schi(1:neg,ng,i)=snu(1:neg,ng,i)
     1         -st*chi(1:neg,ng,iso)
            end if
c one more division for partial normalization.
c needs verification.
            if(ichinorm.eq.2)then
              if(fissdata.eq.2)then ! chi vector
                schi(ng,1:neg,i)=schi(ng,1:neg,i)
     1           /(1.d0-chi(ng,ng,iso))
              else
                schi(1:neg,ng,i)=schi(1:neg,ng,i)
     1           /(1.d0-chi(1:neg,ng,iso))
              end if
            end if
          end do ! ng
        end do ! i
      end if
c
c output
c indirect effect, by isotope and material.
      if(icalc_org.eq.3)then
        write(ius,'(/," feyny sensitivities",10x,"feyny=",
     1   1pe15.7)')yasym
      else if(icalc_org.eq.4)then
        write(ius,'(/," sm2 sensitivities",10x,"sm2=",
     1   1pe15.7)')sm2
      end if
c k2 and m are used to keep track of materials for writing the header.
      k2=0
      m=1
      do i=1,nel
        k2=k2+1
        write(ius,'(/,1x,a24," isotope",i4," in material",i4)')
     1   zaidfull(i),k2,m
        call output1(sden(1:neg,i),stot(1:neg,i),
     1   sabs(1:neg,i),scap(1:neg,i),
     2   schi(1:neg,1:neg,i),snu(1:neg,1:neg,i),
     3   sfis(1:neg,i),ssct(1:neg,1:neg,0:isct,i),
     4   sscti(1:neg,0:isct,i),sscts(1:neg,0:isct,i),
     5   sscto(1:neg,0:isct,i),sins(1:neg,i),
     6   ssctt(1:neg,i))
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
     6   ssctt(1:neg,ix))
      end do ! i
c
      close(ius)
c
      if(icalc.eq.0)then
        rfile="sens_s_r"
      end if
      open(iur,file=rfile,status='unknown')
      write(iur,'(a)')trim(id)
      write(iur,'("relative sensitivities (%/%) to densities, ",
     1 "first order (constant volume)")')
      write(iur,'("first derivatives with respect to ",
     1 "interface locations (constant density)")')
c TODO
c     call calcsens_sm_d(ediv,fissd,dtg,dsg,dfg,dqg,rfile)
c     call calcsens_sm_r(ediv,fissd,lam,rlip,rfile)
      close(iur)
c
c debug using direct pert
c     call wrmacrsx
      if(iver.eq.1)then
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
c             sensverify(ng,lg,l,i)=stot(ng,i)
c for dens
              sensverify(ng,lg,l,i)=sden(ng,i)
c for fiss
c             sensverify(ng,lg,l,i)=sfis(ng,i)
c for nu
c             sensverify(ng,lg,l,i)=sum(snu(1:neg,ng,i))
c for cap
c             sensverify(ng,lg,l,i)=scap(ng,i)
c for scat
c             do lg=1,neg
c               sensverify(ng,lg,l,i)=ssct(ng,lg,l,i)
c               sensverify(neg+1,lg,l,i)=sum(ssct(1:neg,lg,l,i))
c             end do ! lg
            end do ! ng
c for tot
c           sensverify(neg+1,lg,l,i)=sum(stot(1:neg,i))
c for dens
            sensverify(neg+1,lg,l,i)=sum(sden(1:neg,i))
c for fiss
c           sensverify(neg+1,lg,l,i)=sum(sfis(1:neg,i))
c for nu
c           sensverify(neg+1,lg,l,i)=sum(snu(1:neg,1:neg,i))
c for cap
c           sensverify(neg+1,lg,l,i)=sum(scap(1:neg,i))
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
c
      subroutine smcombine(sden,stot,sabs,scap,schi,snu,sfis,
     1 ssct,sscti,sscts,sscto,sins,ssctt)
c combine r1 and r2 sensitivities. input variables are r2
c sensitivities. read r1 from a binary file.
      use COMS
      implicit none
      real*8 sden(neg,2*nel),stot(neg,2*nel),sabs(neg,2*nel),
     1 scap(neg,2*nel),sfis(neg,2*nel),snu(neg,neg,2*nel),
     2 schi(neg,neg,2*nel),ssct(neg,neg,0:isct,2*nel),
     3 sscti(neg,0:isct,2*nel),sscto(neg,0:isct,2*nel),
     4 sscts(neg,0:isct,2*nel),sins(neg,2*nel),ssctt(neg,2*nel)
      real*8 tt
      real*8 sden_tmp(neg,2*nel),stot_tmp(neg,2*nel),
     1 sabs_tmp(neg,2*nel),
     1 scap_tmp(neg,2*nel),sfis_tmp(neg,2*nel),snu_tmp(neg,neg,2*nel),
     2 schi_tmp(neg,neg,2*nel),ssct_tmp(neg,neg,0:isct,2*nel),
     3 sscti_tmp(neg,0:isct,2*nel),sscto_tmp(neg,0:isct,2*nel),
     4 sscts_tmp(neg,0:isct,2*nel),sins_tmp(neg,2*nel),
     1 ssctt_tmp(neg,2*nel)
      integer ns
c
      call rdoutput1_bin(sden_tmp,stot_tmp,sabs_tmp,scap_tmp,
     1 schi_tmp,snu_tmp,sfis_tmp,
     2 ssct_tmp,sscti_tmp,sscts_tmp,sscto_tmp,sins_tmp,ssctt_tmp)
c
      if(icalc.eq.4)then
       tt=1.d0 ! feyny
      else if(icalc.eq.5)then
       tt=2.d0 ! sm2
      end if
      ns=2*nel
      sden(1:neg,1:ns)=sden(1:neg,1:ns)-tt*sden_tmp(1:neg,1:ns)
      stot(1:neg,1:ns)=stot(1:neg,1:ns)-tt*stot_tmp(1:neg,1:ns)
      sabs(1:neg,1:ns)=sabs(1:neg,1:ns)-tt*sabs_tmp(1:neg,1:ns)
      scap(1:neg,1:ns)=scap(1:neg,1:ns)-tt*scap_tmp(1:neg,1:ns)
      sfis(1:neg,1:ns)=sfis(1:neg,1:ns)-tt*sfis_tmp(1:neg,1:ns)
      sins(1:neg,1:ns)=sins(1:neg,1:ns)-tt*sins_tmp(1:neg,1:ns)
      snu(1:neg,1:neg,1:ns)=snu(1:neg,1:neg,1:ns)
     1 -tt*snu_tmp(1:neg,1:neg,1:ns)
      ssct(1:neg,1:neg,0:isct,1:ns)=ssct(1:neg,1:neg,0:isct,1:ns)
     1 -tt*ssct_tmp(1:neg,1:neg,0:isct,1:ns)
      sscti(1:neg,0:isct,1:ns)=sscti(1:neg,0:isct,1:ns)
     1 -tt*sscti_tmp(1:neg,0:isct,1:ns)
      sscto(1:neg,0:isct,1:ns)=sscto(1:neg,0:isct,1:ns)
     1 -tt*sscto_tmp(1:neg,0:isct,1:ns)
      sscts(1:neg,0:isct,1:ns)=sscts(1:neg,0:isct,1:ns)
     1 -tt*sscts_tmp(1:neg,0:isct,1:ns)
      ssctt(1:neg,1:ns)=ssctt(1:neg,1:ns)-tt*ssctt_tmp(1:neg,1:ns)
      return
      end subroutine
