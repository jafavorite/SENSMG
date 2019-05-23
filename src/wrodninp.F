      subroutine wrodninp(sensverify)
c write odninp, cross sections within input file.
c this routine is for perturbed cross sections but can be used
c for unperturbed also.
      use COMS
      implicit none
      real*8 sensverify(neg+1,neg+1,0:isct,nel)
      real*8 sigt0(neg,nxs),nusigf0(neg,nxs),siga0(neg,nxs),
     1 chi0(neg,neg,nxs),sigs0(neg,neg,0:isct,nxs),
     2 sigf0(neg,nxs),sigc0(neg,nxs),sige0(neg,nxs),sigi0(neg,nxs)
      real*8 xsa(neg),xsn(neg),xst(neg),xsf(neg),xsc(neg),xse(neg),
     1 xsi(neg),xss(neg,2*neg-1),adens,
     2 pertf,pertn,perta,perts,pertt,pertd,pertx,pertc,perte,perti,
     3 ppf,ppn,ppa,pps,ppx,ppc,ppe,ppi,ppt,pp,
     4 pptadd,pptaddi,ppsadd,ppsaddi,pert,unpert,
     5 diff,maxdif_sigt(2),maxdif_nusigf(2),maxdif_sigs(2),
     6 chinum(neg),chiden,nusigf0j(nel)
      integer k1,k2,iso,i,ia(nel),iitm,i_srcacc,j,kg,kh,kl,km,ki,l,lg,
     1 mt,ng,n0,np,igrp,igrp_in,imt,iis,il,i_use_chimat,i_use_wgtfct
      character pfile*120,title*7,ct*3,typel*100,ctmp*9,c8*8,edname(4)*7
      integer iht,ihs,ihm,ith,ifido
      logical np0done
      data edname/"(n,g)","n-fiss","(n,n)","(n,n') "/
      iitm=999   ! if iitm.gt.0 it will be used for forward and adjoint.
      i_srcacc=0 ! 0/1 for default/srcacc=no
      ifido=-1   ! -1/0 for 4e18/6e12
      i_use_chimat=0 ! 0/1 recompute material chi/use given material chi
      i_use_wgtfct=0 ! 0/1 use 1 for wgtfct/use wgtfct from NDI
c
      write(*,'("writing verification inputs.")')
      write(iuo,'("writing verification inputs.")')
      if(fissdata.eq.2)then
        write(*,'("assumes fissdata=2. (it does.)")')
        write(iuo,'("assumes fissdata=2. (it does.)")')
      else
        write(*,'("assumes fissdata=2.")')
        write(*,'("** warning. fissdata=",i1)')fissdata
        write(iuo,'("assumes fissdata=2.")')
        write(iuo,'("** warning. fissdata=",i1)')fissdata
      end if
      np0done=.false.
      k1=0
      do km=1,1 !nm
      do ki=1,ncb(km)
      k1=k1+1
c     if(ki.ne.44.and.ki.ne.46.and.ki.ne.49)cycle
      write(*,'(" isotope",i4," in material",i7)')ki,km
      write(iuo,'(" isotope",i4," in material",i7)')ki,km
      do kg=1,neg+1
      do kh=1,1 ! 1,neg for scattering; kg,kg for self; otherwise 1,1
      do kl=0,0 ! 0,isct for scattering, otherwise 0,0
c
      do np=-1,1
      if(np.eq.0.and.np0done)cycle
c set the group, mat, isotope, and moment to perturb and the pert fraction.
c igrp=neg+1 means all groups.
c for scattering, igrp is the outgoing group. igrp_in is the ingoing
c group.
c pert fractions must be positive.
      igrp=kg           ! group to perturb (out for scat)
      imt=km            ! material
      iis=ki            ! isotope in material imt
      il=kl             ! moment for scat
      igrp_in=kh        ! in group for scat
      pertf=0.000d0     ! fission
      pertn=0.000d0     ! nu only (does not affect total)
      perta=0.000d0     ! absorption
      pertc=0.000d0     ! capture
      perts=0.000d0     ! scattering
      pertt=0.000d0     ! total (nothing else)
      pertd=0.000d0     ! density
      pertx=0.000d0     ! chi (fission spectrum)
      perte=0.000d0     ! elastic scatter
      perti=0.000d0     ! inelastic scatter
c for scattering
c     if(sensverify(kh,kg,kl,k1).eq.0.d0)then
c for scattering
c     if(sensverify(kg,kh,kl,k1).eq.0.d0)then
c       write(*,'("sensverify=0. kh,kg,kl,k1=",4i5)')kh,kg,kl,k1
c       cycle
c     end if
      if(np.eq.-1)then
c       pertt=abs(0.001d0/sensverify(kg,1,0,k1))
c       pertt=min(0.1d0,pertt)
c       pert=pertt
c       perts=abs(0.001d0/sensverify(kh,kg,kl,k1))
c       perts=min(0.75d0,perts)
c       pert=perts
c note that changing density this way does not change the source.
c       pertd=abs(0.001d0/sensverify(kg,1,0,k1))
c       pertd=min(0.999d0,pertd)
c       pert=pertd
c       pertf=abs(0.001d0/sensverify(kg,1,0,k1))
c       pertf=min(0.999d0,pertf)
c       pert=pertf
c       pertn=abs(0.001d0/sensverify(kg,1,0,k1))
c       pertn=min(0.9d0,pertn)
c       pert=pertn
c       pertc=abs(0.01d0/sensverify(kg,1,0,k1))
c       pertc=min(0.999d0,pertc)
c       pert=pertc
        pertx=0.1d0 ! abs(0.001d0/sensverify(kg,1,0,k1))
        pertx=min(0.999d0,pertx)
        pert=pertx
c partisn can't read "(n,n') "; don't use these for now
c       perte=abs(0.01d0/sensverify(kg,1,0,k1))
c       perte=min(0.999d0,perte)
c       pert=perte
c       perti=abs(0.01d0/sensverify(kg,1,0,k1))
c       perti=min(0.999d0,perti)
c       pert=perti
      else if(np.ge.0)then
c       pertt=pert
c       perts=pert
c       pertd=pert
c       pertf=pert
c       pertn=pert
c       pertc=pert
        pertx=pert
c       perte=pert
c       perti=pert
      end if
   10 ct="   "
      unpert=0.d0
      pp=0.d0
      if(pertf.gt.0.d0)then
        ct="fis"
c       if(np.eq.-1)write(*,'("A",2i3,1pe15.7)')k1,kg,pertf
        write(ctmp,'(1pe8.2)')pertf
        read(ctmp,'(e8.2)')pertf
c       if(np.eq.-1)write(*,'("B",2i3,1pe15.7)')k1,kg,pertf
        pp=pertf
      end if
      if(pertn.gt.0.d0)then
        if(ct.ne."   ")write(*,'("warning in wrodninp. ",
     1   "multiple perturbations.")')
        ct="0nu"
        write(ctmp,'(1pe8.2)')pertn
        read(ctmp,'(e8.2)')pertn
        pp=pertn
      end if
      if(perta.gt.0.d0)then
        if(ct.ne."   ")write(*,'("warning in wrodninp. ",
     1   "multiple perturbations.")')
        ct="abs"
        write(ctmp,'(1pe8.2)')perta
        read(ctmp,'(e8.2)')perta
        pp=perta
      end if
      if(pertc.gt.0.d0)then
        if(ct.ne."   ")write(*,'("warning in wrodninp. ",
     1   "multiple perturbations.")')
        ct="cap"
        write(ctmp,'(1pe8.2)')pertc
        read(ctmp,'(e8.2)')pertc
        pp=pertc
      end if
      if(perts.gt.0.d0)then
        if(ct.ne."   ")write(*,'("warning in wrodninp. ",
     1   "multiple perturbations.")')
        ct="sct"
        write(ctmp,'(1pe8.2)')perts
        read(ctmp,'(e8.2)')perts
        pp=perts
      end if
      if(pertt.gt.0.d0)then
        if(ct.ne."   ")write(*,'("warning in wrodninp. ",
     1   "multiple perturbations.")')
        ct="tot"
        write(ctmp,'(1pe8.2)')pertt
        read(ctmp,'(e8.2)')pertt
        pp=pertt
      end if
      if(pertd.gt.0.d0)then
        if(ct.ne."   ")write(*,'("warning in wrodninp. ",
     1   "multiple perturbations.")')
        ct="den"
        write(ctmp,'(1pe8.2)')pertd
        read(ctmp,'(e8.2)')pertd
        pp=pertd
      end if
      if(pertx.gt.0.d0)then
        if(ct.ne."   ")write(*,'("warning in wrodninp. ",
     1   "multiple perturbations.")')
        ct="chi"
        write(ctmp,'(1pe8.2)')pertx
        read(ctmp,'(e8.2)')pertx
        pp=pertx
      end if
      if(perte.gt.0.d0)then
        if(ct.ne."   ")write(*,'("warning in wrodninp. ",
     1   "multiple perturbations.")')
        ct="ela"
        write(ctmp,'(1pe8.2)')perte
        read(ctmp,'(e8.2)')perte
        pp=perte
      end if
      if(perti.gt.0.d0)then
        if(ct.ne."   ")write(*,'("warning in wrodninp. ",
     1   "multiple perturbations.")')
        ct="ine"
        write(ctmp,'(1pe8.2)')perti
        read(ctmp,'(e8.2)')perti
        pp=perti
      end if
      if(np.eq.0)pp=0.d0
c
c zero the cross sections.
      nusigf0(1:neg,1:nxs)=0.d0
      sigt0(1:neg,1:nxs)=0.d0
      siga0(1:neg,1:nxs)=0.d0
      sigs0(1:neg,1:neg,0:isct,1:nxs)=0.d0
      chi0(1:neg,1:neg,1:nxs)=0.d0
      nusigf0j(1:nel)=0.d0
c     nusigf0(1:neg,nm+1:nm+niso)=nusigf(1:neg,nm+1:nm+niso)
c     sigt0(1:neg,nm+1:nm+niso)=sigt(1:neg,nm+1:nm+niso)
c     siga0(1:neg,nm+1:nm+niso)=siga(1:neg,nm+1:nm+niso)
c     sigs0(1:neg,1:neg,0:isct,nm+1:nm+niso)=
c    1 sigs(1:neg,1:neg,0:isct,nm+1:nm+niso)
c     sigf0(1:neg,nm+1:nm+niso)=sigf(1:neg,nm+1:nm+niso)
c     sigc0(1:neg,nm+1:nm+niso)=sigc(1:neg,nm+1:nm+niso)
      sigf0(1:neg,1:nxs)=0.d0
      sigc0(1:neg,1:nxs)=0.d0
      sige0(1:neg,1:nxs)=0.d0
      sigi0(1:neg,1:nxs)=0.d0
c
c set macroscopic material cross sections.
      k2=0
      do i=1,nm
        do j=1,ncb(i)
          k2=k2+1
          iso=nm+k2
          do ng=1,neg
            adens=blk(3,k2)
            if(np.ne.0.and.(ng.eq.igrp.or.igrp.eq.neg+1).and.i.eq.imt
     1       .and.j.eq.iis.and.pertd.gt.0.d0)then
              unpert=adens
              adens=adens*(1.d0+sign(pertd,real(np,8)))
            end if
            pptadd=0.d0
            pptaddi=0.d0
            ppsadd=0.d0
            ppsaddi=0.d0
            ppn=1.d0
            if(np.ne.0.and.(ng.eq.igrp.or.igrp.eq.neg+1).and.i.eq.imt
     1       .and.j.eq.iis.and.pertn.gt.0.d0)then
              ppn=1.d0+sign(pertn,real(np,8))
              if(sigf(ng,iso).eq.0.d0)then
                write(*,'("error. sig=0. probably skipped reading c",
     1           " in rdsnxedt.")')
                call stoponerror
              end if
              unpert=nusigf(ng,iso)/sigf(ng,iso)
            end if
            ppf=1.d0
            if(np.ne.0.and.(ng.eq.igrp.or.igrp.eq.neg+1).and.i.eq.imt
     1       .and.j.eq.iis.and.pertf.gt.0.d0)then
              ppf=1.d0+sign(pertf,real(np,8))
              pptadd=pptadd+sigf(ng,iso)*adens*(ppf-1.d0)
              pptaddi=pptaddi+sigf(ng,iso)*(ppf-1.d0)
              unpert=sigf(ng,iso)
            end if
c only one of ppn or ppf should be > 1.
            nusigf0(ng,i)=nusigf0(ng,i)+nusigf(ng,iso)*adens
     1       *ppn*ppf
            nusigf0(ng,iso)=nusigf(ng,iso)*ppn*ppf
            sigf0(ng,iso)=sigf(ng,iso)*ppn*ppf
            if(i_use_wgtfct.eq.0)then
              nusigf0j(k2)=nusigf0j(k2)+nusigf(ng,iso)*adens
     1         *ppn*ppf
            else if(i_use_wgtfct.eq.1)then
              nusigf0j(k2)=nusigf0j(k2)+nusigf(ng,iso)*adens
     1         *wgtfct(ng,k2)*ppn*ppf
            end if
c           if(np.eq.-1)write(*,'("C",4i3,1p8e15.7)')k1,kg,k2,ng,
c    1       pertf,ppf,nusigf0(ng,i),nusigf(ng,iso),adens,ppn,
c    2       sigf(ng,i),sigf(ng,iso)
            ppe=1.d0
            if(np.ne.0.and.(ng.eq.igrp.or.igrp.eq.neg+1).and.i.eq.imt
     1       .and.j.eq.iis.and.perte.gt.0.d0)then
              ppe=1.d0+sign(perte,real(np,8))
              pptadd=pptadd+sige(ng,iso)*adens*(ppe-1.d0)
              pptaddi=pptaddi+sige(ng,iso)*(ppe-1.d0)
              ppsadd=ppsadd+sige(ng,iso)*adens*(ppe-1.d0)
              ppsaddi=ppsaddi+sige(ng,iso)*(ppe-1.d0)
              unpert=sige(ng,iso)
            end if
            sige0(ng,iso)=sige(ng,iso)*ppe
            ppi=1.d0
            if(np.ne.0.and.(ng.eq.igrp.or.igrp.eq.neg+1).and.i.eq.imt
     1       .and.j.eq.iis.and.perti.gt.0.d0)then
              ppi=1.d0+sign(perti,real(np,8))
              pptadd=pptadd+sigi(ng,iso)*adens*(ppi-1.d0)
              pptaddi=pptaddi+sigi(ng,iso)*(ppi-1.d0)
              ppsadd=ppsadd+sigi(ng,iso)*adens*(ppi-1.d0)
              ppsaddi=ppsaddi+sigi(ng,iso)*(ppi-1.d0)
              unpert=sigi(ng,iso)
            end if
            sigi0(ng,iso)=sigi(ng,iso)*ppi
            do lg=1,neg
              do l=0,isct
                pps=1.d0
                if(np.ne.0.and.(ng.eq.igrp.or.igrp.eq.neg+1)
     1           .and.(lg.eq.igrp_in.or.igrp_in.eq.neg+1)
     2           .and.i.eq.imt.and.j.eq.iis.and.l.eq.il
     3           .and.perts.gt.0.d0)then
                  pps=1.d0+sign(perts,real(np,8))
c only zero'th order affects total.
                  if(l.eq.0)then
                    pptadd=pptadd+sigs(lg,ng,l,iso)*adens*(pps-1.d0)
                    pptaddi=pptaddi+sigs(lg,ng,l,iso)*(pps-1.d0)
                  end if
                  unpert=sigs(lg,ng,l,iso)
                end if
c TODO include ppsadd and ppsaddi here
                sigs0(lg,ng,l,i)=sigs0(lg,ng,l,i)+sigs(lg,ng,l,iso)
     1           *adens*pps
                sigs0(lg,ng,l,iso)=sigs(lg,ng,l,iso)*pps
              end do ! l
            end do ! lg
            ppa=1.d0
            if(np.ne.0.and.(ng.eq.igrp.or.igrp.eq.neg+1).and.i.eq.imt
     1       .and.j.eq.iis.and.perta.gt.0.d0)then
              ppa=1.d0+sign(perta,real(np,8))
              pptadd=pptadd+siga(ng,iso)*adens*(ppa-1.d0)
              pptaddi=pptaddi+siga(ng,iso)*(ppa-1.d0)
              unpert=siga(ng,iso)
            end if
            siga0(ng,i)=siga0(ng,i)+siga(ng,iso)*adens*ppa
            siga0(ng,iso)=siga(ng,iso)*ppa
            ppc=1.d0
            if(np.ne.0.and.(ng.eq.igrp.or.igrp.eq.neg+1).and.i.eq.imt
     1       .and.j.eq.iis.and.pertc.gt.0.d0)then
              ppc=1.d0+sign(pertc,real(np,8))
              pptadd=pptadd+sigc(ng,iso)*adens*(ppc-1.d0)
              pptaddi=pptaddi+sigc(ng,iso)*(ppc-1.d0)
              unpert=sigc(ng,iso)
            end if
            sigc0(ng,iso)=sigc(ng,iso)*ppc
            ppt=1.d0
            if(np.ne.0.and.(ng.eq.igrp.or.igrp.eq.neg+1).and.i.eq.imt
     1       .and.j.eq.iis.and.pertt.gt.0.d0)then
              ppt=1.d0+sign(pertt,real(np,8))
              unpert=sigt(ng,iso)
            end if
            sigt0(ng,i)=sigt0(ng,i)+sigt(ng,iso)*adens*ppt+pptadd
            sigt0(ng,iso)=sigt(ng,iso)*ppt+pptaddi
c check P0 scattering and total for isotope.
            if(sigt(ng,iso)*adens*ppt+pptadd.lt.
     1       sigs(ng,ng,0,iso)*adens*pps)then
              write(*,'("p0 self-scatter greater than total; ",
     1         "iso,ng,sigt,sigs=",2i4,1p2e14.6)')k2,ng,
     2         sigt(ng,iso)*adens*ppt+pptadd,sigs(ng,ng,0,iso)*adens*pps
            end if
          end do ! ng
        end do ! j
      end do ! i
c perturb chi with full normalization
      k2=0
      do i=1,nm
        do j=1,ncb(i)
          k2=k2+1
          iso=nm+k2
          chiden=0.d0
          do ng=1,neg
            ppx=1.d0
            if(np.ne.0.and.ng.eq.igrp.and.i.eq.imt
     1       .and.j.eq.iis.and.pertx.gt.0.d0)then
              ppx=1.d0+sign(pertx,real(np,8))
              unpert=chi(ng,ng,iso)
            end if
            chi0(ng,ng,iso)=chi(ng,ng,iso)*ppx
            chiden=chiden+chi0(ng,ng,iso)
          end do ! ng
          if(ichinorm.eq.1.and.chiden.gt.0.d0)then
            do ng=1,neg
              chi0(ng,ng,iso)=chi0(ng,ng,iso)/chiden
            end do ! ng
          end if
        end do ! j
      end do ! i
c
      if(np.ne.0.and.unpert.eq.0.d0)then
        write(*,'("unpert.eq.0. iis,imt,igrp,igrp_in,np",5i5)')iis,imt,
     1   igrp,igrp_in,np
        cycle
      end if
c set material chi
      k2=0
      do i=1,nm
        chinum(1:neg)=0.d0
        chiden=0.d0
        do j=1,ncb(i)
          k2=k2+1
          iso=nm+k2
          do ng=1,neg
            chinum(ng)=chinum(ng)+chi0(ng,ng,iso)*nusigf0j(k2)
          end do ! ng
c chiden is equal to sum(nusigf0(1:neg,i))
          chiden=chiden+nusigf0j(k2)
        end do ! j
        if(chiden.eq.0.d0)cycle
        do ng=1,neg
          if(i_use_chimat.eq.0)then
            chi0(ng,ng,i)=chinum(ng)/chiden
          else if(i_use_chimat.eq.1)then
            chi0(ng,ng,i)=chi(ng,ng,i)
          end if
        end do ! ng
      end do ! i
c check P0 scattering and total for material.
      do i=1,nm
        do ng=1,neg
          if(sigt0(ng,i).lt.sigs0(ng,ng,0,i))then
            write(*,'("p0 self-scatter greater than total; ",
     1       "ng,sigt,sigs=",i4,1p2e14.6)')ng,sigt0(ng,i),
     2       sigs0(ng,ng,0,i)
            if(pertt.gt.0.d0)then
              pertt=0.9d0*pertt
            else if(perts.gt.0.d0)then
              perts=0.9d0*perts
            end if
            go to 10 
          end if
        end do ! ng
      end do ! i
c
c size of cross section tables. can make these a function of
c upscattering. iht=3 always.
      iht=3
      ihs=3+neg
      ihm=3+2*neg-1
c
      if(np.eq.0)then
        igrp=0
        imt=0
        iis=0
        il=0
        igrp_in=0
        np0done=.true.
      end if
c pert_xxx_00_00_00_m
      if(ct.eq."tot".or.ct.eq."den".or.ct.eq."fis".or.ct.eq."0nu"
     1 .or.ct.eq."cap".or.ct.eq."abs".or.ct.eq."chi"
     2 .or.ct.eq."ela".or.ct.eq."ine")then
        if(neg.lt.100)then
          write(pfile,'("pert_",a3,"_",i2.2,"_",i2.2,"_",i2.2)')
     1     ct,imt,iis,igrp
          write(typel,'(" type,mat,iso,group: ",a3,3i3,
     1     1pe12.4,e19.10)')ct,imt,iis,igrp,sign(pp,real(np,8)),
     2     unpert
        else if(neg.lt.1000)then
          write(pfile,'("pert_",a3,"_",i2.2,"_",i2.2,"_",i3.3)')
     1     ct,imt,iis,igrp
          write(typel,'(" type,mat,iso,group: ",a3,3i4,
     1     1pe12.4,e19.10)')ct,imt,iis,igrp,sign(pp,real(np,8)),
     2     unpert
        end if
      else if(ct.eq."sct")then
        if(neg.lt.100)then
          write(pfile,'("pert_",a3,"_",i2.2,"_",i2.2,"_",i2.2,
     1     "_",i2.2,"_",i2.2)')ct,imt,iis,igrp,igrp_in,il
          write(typel,'(" type,mat,iso,group,ingrp,l: ",a3,5i3,
     1     1pe12.4,e19.10)')ct,imt,iis,igrp,igrp_in,il,
     2     sign(pp,real(np,8)),unpert
        else if(neg.lt.1000)then
          write(pfile,'("pert_",a3,"_",i2.2,"_",i2.2,"_",i3.3,
     1     "_",i3.3,"_",i2.2)')ct,imt,iis,igrp,igrp_in,il
          write(typel,'(" type,mat,iso,group,ingrp,l: ",a3,5i4,
     1     1pe12.4,e19.10)')ct,imt,iis,igrp,igrp_in,il,
     2     sign(pp,real(np,8)),unpert
        end if
      end if
      if(np.eq.-1)then
        pfile=pfile(1:len_trim(pfile))//"_m"
      else if(np.eq.0)then
        write(pfile,'("pert_",a3,"_0")')ct
      else if(np.eq.1)then
        pfile=pfile(1:len_trim(pfile))//"_p"
      else
        write(*,'("unexpected pert index in wrodninp.")')
        write(pfile,'("pert_x")')
      end if
      ith=0 ! this is for forward
      if(ith.eq.0)then
        if(isrcacc_no.eq.1.or.isrcacc_no.eq.3)i_srcacc=1
        title="forward"
      else if(ith.eq.1)then
        if(isrcacc_no.eq.2.or.isrcacc_no.eq.3)i_srcacc=1
        title="adjoint"
      end if
      open(iup,file=pfile,status='unknown')
      write(iup,'(5i6)')3,0,0,0,0
      write(iup,'(a)')trim(id)
      if(icalc.eq.0.or.icalc.eq.3.or.icalc.eq.4)then
        write(iup,'(1x,a7," input file, fixed-source")')title
      else if(icalc.eq.1)then
        write(iup,'(1x,a7," input file, keff")')title
      else if(icalc.eq.2)then
        write(iup,'(1x,a7," input file, alpha")')title
      end if
      write(iup,'(a)')trim(typel)
c
      mt=nm
      n0=nm
      call block1(iup,mt,n0)
      call block2(iup)
      write(iup,'("/ * * * * block iii * * * *")')
      write(iup,'(" lib=odninp")')
      write(iup,'(" iht=",i4," ihs=",i4," ihm=",i4)')iht,ihs,ihm
      write(iup,'(" ifido=",i2," ititl=1")')ifido
      write(iup,'(" maxord=",i2)')isct
      write(iup,'(" names= i01")')
      do i=2,nm
        write(iup,'("        i",i2.2)')i
      end do ! i
      write(iup,'(" lng=",i3)')lng
      if(icalc.eq.2)then
        write(iup,'(" vel=")')
        write(iup,'(1p,4e18.11)')(vel(ng),ng=1,neg)
      end if
c this applies to the whole problem. use if only one fissionable
c region. otherwise use chi in block v.
      write(iup,'(" chivec=")')
      write(iup,'(2x,1p4e19.11)')(chi0(ng,ng,1),ng=1,neg)
      write(iup,'(2x,";")')
      write(iup,'(" t")')
c
c write cross sections.
      do i=1,nm
        xsa(1:neg)=siga0(1:neg,i)
        xsn(1:neg)=nusigf0(1:neg,i)
        xst(1:neg)=sigt0(1:neg,i)
        do l=0,isct
          if(ifido.eq.0)then
            write(iup,'(1x,"siga",8x,"nusigf",6x,"sigt",
     1       9x,"mat,ord=",2i3)')i,l
          elseif(ifido.eq.-1)then
            write(iup,'(1x,"siga",14x,"nusigf",12x,"sigt",
     1       9x,"mat,ord=",2i3)')i,l
          end if
          xss(1:neg,1:2*neg-1)=0.d0
          do ng=1,neg
            do lg=(neg-1),-(neg-1),-1
              if(ng+lg.lt.1.or.ng+lg.gt.neg)cycle
              xss(ng,neg-lg)=sigs0(ng,ng+lg,l,i)
            end do !lg
          end do ! ng
          if(ifido.eq.0)then
            write(iup,'(1p,6e12.5 )')(xsa(ng),xsn(ng),
     1       xst(ng),(xss(ng,lg),lg=1,2*neg-1),ng=1,neg)
          elseif(ifido.eq.-1)then
            write(iup,'(1p,4e18.11)')(xsa(ng),xsn(ng),
     1       xst(ng),(xss(ng,lg),lg=1,2*neg-1),ng=1,neg)
          end if
          if(l.gt.0)cycle
          xsa(1:neg)=0.d0
          xsn(1:neg)=0.d0
          xst(1:neg)=0.d0
        end do ! l
      end do ! i
c
      write(iup,'("/ * * * * block iv * * * *")')
      write(iup,'(" matspec=atdens")')
      write(iup,'(" matls= m01   i01  1.;")')
      do i=2,nm
        write(iup,'("        ","m",i2.2,"   i",i2.2,"  1.;")')i,i
      end do ! i
      write(iup,'(" assign= zone01   m01  1.;")')
      do i=2,nm
        write(iup,'("         zone",i2.2,"   m",i2.2,"  1.;")')i,i
      end do ! i
      write(iup,'(" zonetemp= 2.5301E-08 ;")')
      do i=2,nm
        write(iup,'("           2.5301E-08 ;")')
      end do ! i
      write(iup,'(" t")')
c
      write(iup,'("/ * * * * block v * * * *")')
      if(icalc.eq.0.or.icalc.eq.3.or.icalc.eq.4)then
        write(iup,'(" ievt=",i1," isct=",i1)')0,isct
      else
        write(iup,'(" ievt=",i1," isct=",i1)')icalc,isct
      end if
c     if(np.eq.0)write(*,'("warning. reflecting boundaries.")')
c     write(iup,'(" ibl=1 ibr=1 / reflective")')
c special, 8_31, don't normalize chi
      if(ichinorm.eq.0.and.partisnver(1:4).eq."8_31")then
        write(iup,'(" extras=52r 0 1 f 0 / dont normalize chi")')
      end if
      write(iup,'(" epsi=",1pe9.2," balp=1")')1.d-8
      if(icalc.eq.1.or.icalc.eq.2)then
        write(iup,'(" norm=1.0")')
      end if
      write(iup,'(" npeg=2")')
      if(idbgw.eq.1)write(iup,'(" fluxp=2")')
c     write(iup,'(" fluxp=2")')
      write(iup,'(" rmflux=1 raflux=1")')
      if(iitm.gt.0)write(iup,'(" iitm=",i3)')iitm
      write(iup,'(" iitl=0 oitm=9999")')
      if(i_srcacc.eq.1)write(iup,'(" srcacc=no")')
      if(nofxup.eq.1)write(iup,'(" nofxup=1")')
c
c different solve parameters possible for forward and adjoint files 
c forward:
      if(ith.eq.0)then
        write(iup,'(" ith=0 xsectp=2")')
c     
c adjoint:
      else if(ith.eq.1)then
        write(iup,'(" ith=1 xsectp=2")')
      end if
c for test problems with more than one fissionable region, this is
c needed. otherwise use chivec so adjoint will have perturbed chi.
c     write(iup,'(" chi=")')
c     do i=1,nm
c       write(iup,'(2x,1p4e19.11)')(chi0(ng,ng,i),ng=1,neg)
c       write(iup,'(2x,";")')
c     end do ! i
      if(icalc.eq.0.or.icalc.eq.3.or.icalc.eq.4)then
        call fixedsrc(iup,ith)
      end if
      write(iup,'(" t")')
c no edit for adjoint.
      if(ith.eq.1)go to 100
      write(iup,'("/ * * * * block vi * * * *")')
      write(iup,'(" ajed=",i1)')ith
      write(iup,'(" zned=1 igrped=0")')
      write(iup,'(" rsfnam=")')
      do i=1,nrrx
        write(iup,'(2x,1h","m",i6.6,1h"," ;")')nm+i
      end do ! i
      write(iup,'(2x,1h","flux",1h"," ;")')
      write(iup,'(" rsfe=")')
      do i=1,nrrx
        write(iup,'(2x,1p5e14.6)')(rrxs(ng,i),ng=1,neg)
        write(iup,'(2x,";")')
      end do ! i
      write(iup,'(2xi3,"r 1. ;")')neg
      write(iup,'(" t")')
c
  100 close(iup)
c
c cross-section input file in directory xs1
c check for repeats
      ia(1:nel)=1
      do i=2,nel
        do j=1,i-1
          if(zaidfull(i).eq.zaidfull(j))then
            ia(i)=0
            exit
          end if
        end do ! j
      end do ! i
c recall niso=nel+nrrx
      n0=sum(ia(1:nel))
      mt=niso
c size of cross section tables. add 4 to accommodate edits.
c     iht=iht+4
c     ihs=ihs+4
c     ihm=ihm+4
c partisn can't read "(n,n') " for now
      iht=iht+2
      ihs=ihs+2
      ihm=ihm+2
c
c pxs1_xxx_00_00_00_m
      pfile(1:4)="pxs1"
      open(iup,file=pfile,status='unknown')
      write(iup,'(5i6)')5,0,0,0,0
      write(iup,'(a)')trim(id)
      write(iup,'("input: ",a)')trim(ifile)
      write(iup,'(" cross section input file; used only for ",/,
     1 " obtaining cross sections from macrxs.")')
      write(iup,'(a)')trim(typel)
c
      write(iup,'("/ * * * * block i * * * *")')
      write(iup,'(" igeom=sphere isn=",i3," ngroup=",i4)')2,neg
c     write(iup,'("/maxscm=20000000 maxlcm=9000000")')
      write(iup,'(" niso=",i6," mt=",i4," nzone=",i4)')n0,mt,mt
      write(iup,'(" im=",i4," it=",i5)')mt,4*mt
c for xsectp to take effect (to print the cross sections), comment this
      write(iup,'(" nosolv=1")')
      write(iup,'(" t")')
c
      write(iup,'("/ * * * * block ii * * * *")')
      write(iup,'(" xmesh=  0.")')
      do i=1,mt
        write(iup,'(6x,0pf6.1)')real(i,8)
      end do ! i
      write(iup,'(" xints= ",i4,"r 4")')mt
      write(iup,'(" zones=")')
      do i=1,mt
        write(iup,'(4x,i6)')i
      end do ! i
      write(iup,'(" t")')
c
      write(iup,'("/ * * * * block iii * * * *")')
      write(iup,'(" lib=odninp")')
      write(iup,'(" iht=",i4," ihs=",i4," ihm=",i4)')iht,ihs,ihm
      write(iup,'(" ifido=",i2," ititl=1")')ifido
      write(iup,'(" maxord=",i2)')isct
      write(iup,'(" names=")')
c limited to 8 characters.
      do i=1,nel
        if(ia(i).eq.0)cycle
        c8=adjustl(zaidfull(i))
        write(iup,'(8x,1h",a,1h")')trim(c8)
      end do ! i
c partisn can't read "(n,n') " for now
c     write(iup,'(" edname=",1x,1h",a,1h",1x,1h",a,1h",1x,1h",a,1h",
c    1 1x,1h",a,1h")')(trim(edname(i)),i=1,4)
      write(iup,'(" edname=",1x,1h",a,1h",
     1 1x,1h",a,1h")')(trim(edname(i)),i=1,2)
      write(iup,'(" lng=",i3)')lng
      if(icalc.eq.2)then
        write(iup,'(" vel=")')
        write(iup,'(1p,4e18.11)')(vel(ng),ng=1,neg)
      end if
      write(iup,'(" t")')
c
c write cross sections.
      do i=1,niso
        if(ia(i).eq.0)cycle
        iso=nm+i
        xsa(1:neg)=siga0(1:neg,iso)
        xsn(1:neg)=nusigf0(1:neg,iso)
        xst(1:neg)=sigt0(1:neg,iso)
        xsf(1:neg)=sigf0(1:neg,iso)
        xsc(1:neg)=sigc0(1:neg,iso)
        xse(1:neg)=sige0(1:neg,iso)
        xsi(1:neg)=sigi0(1:neg,iso)
        do l=0,isct
          write(iup,'(1x,"sigc sigf sige sigi siga nusigf sigt",
     1     9x,"iso,ord=",2i3)')iso,l
          xss(1:neg,1:2*neg-1)=0.d0
          do ng=1,neg
            do lg=(neg-1),-(neg-1),-1
              if(ng+lg.lt.1.or.ng+lg.gt.neg)cycle
              xss(ng,neg-lg)=sigs0(ng,ng+lg,l,iso)
            end do !lg
          end do ! ng
          if(ifido.eq.0)then
c partisn can't read "(n,n') " for now
c           write(iup,'(1p,6e12.5 )')(xsc(ng),xsf(ng),xse(ng),xsi(ng),
c    1       xsa(ng),xsn(ng),xst(ng),(xss(ng,lg),lg=1,2*neg-1),ng=1,neg)
            write(iup,'(1p,6e12.5 )')(xsc(ng),xsf(ng),
     1       xsa(ng),xsn(ng),xst(ng),(xss(ng,lg),lg=1,2*neg-1),ng=1,neg)
          elseif(ifido.eq.-1)then
c partisn can't read "(n,n') " for now
c           write(iup,'(1p,4e18.11)')(xsc(ng),xsf(ng),xse(ng),xsi(ng),
c    1       xsa(ng),xsn(ng),xst(ng),(xss(ng,lg),lg=1,2*neg-1),ng=1,neg)
            write(iup,'(1p,4e18.11)')(xsc(ng),xsf(ng),
     1       xsa(ng),xsn(ng),xst(ng),(xss(ng,lg),lg=1,2*neg-1),ng=1,neg)
          end if
          if(l.gt.0)cycle
          xsa(1:neg)=0.d0
          xsn(1:neg)=0.d0
          xst(1:neg)=0.d0
          xsf(1:neg)=0.d0
          xsc(1:neg)=0.d0
          xse(1:neg)=0.d0
          xsi(1:neg)=0.d0
        end do ! l
      end do ! i
c
      write(iup,'("/ * * * * block iv * * * *")')
      write(iup,'(" matspec=atdens")')
      write(iup,'(" matls=")')
      do i=1,niso
        c8=adjustl(zaidfull(i))
        write(iup,'(8x,"m",i2.2,2x,1h",a,1h","  1.;")')i,trim(c8)
      end do ! i
      write(iup,'(" assign= zone01   m01  1.;")')
      do i=2,niso
        write(iup,'("         zone",i2.2,"   m",i2.2,"  1.;")')i,i
      end do ! i
      write(iup,'(" zonetemp= 2.5301E-08 ;")')
      do i=2,niso
        write(iup,'("           2.5301E-08 ;")')
      end do ! i
      write(iup,'(" t")')
c
      write(iup,'("/ * * * * block v * * * *")')
      write(iup,'(" ievt=0 isct=",i1)')isct
      write(iup,'(" iitm=0 kitm=0 iitl=0 oitm=2")')
      write(iup,'(" npeg=2")')
      write(iup,'(" ith=0 xsectp=2 epsi=10.")')
      write(iup,'(" sirite= ",i4,"r 1.")')neg
      write(iup,'(" t")')
c
      close(iup)
c
c pchi_xxx_00_00_00_m
      pfile(1:4)="pchi"
      open(iup,file=pfile,status='unknown',position='rewind',
#if ( defined INTEL )
     1 form='unformatted',access='stream'
#elif ( defined GNU )
     1 access='stream'  ! untested capability
#endif
     2 )
      write(iup)neg,nm,nel,nrrx
      do ng=1,neg
        write(iup)(chi0(ng,ng,i),i=1,nxs)
      end do ! ng
      write(iup)((sigf0(ng,i),i=1,niso),ng=1,neg)
      close(iup)
c     write(*,'(a,2i6)')trim(pfile),k1,iis
c     do ng=1,neg
c       write(*,'(i4,1pe19.11)')ng,chi0(ng,ng,nm+k1)
c     end do ! ng
c
c write cross sections as matrices.
c use avg to skip; 0 to run.
      if(np.eq.999)then
        i=1
        iso=i+2
        write(*,'("nusigf")')
        write(*,'(1p9999e18.11)')nusigf0(1:neg,iso)
        write(*,'("sigf")')
        write(*,'(1p9999e18.11)')sigf0(1:neg,iso)
        write(*,'("nu")')
        write(*,'(1p9999e18.11)')nusigf0(1:neg,iso)/sigf0(1:neg,iso)
        write(*,'("sigt")')
        write(*,'(1p9999e18.11)')sigt0(1:neg,iso)
        write(*,'("sigc")')
        write(*,'(1p9999e18.11)')sigc0(1:neg,iso)
        write(*,'("sigc_tot")')
        write(*,'(1p9999e18.11)')(sigt0(ng,iso)-sigf0(ng,iso)-
     1   sum(sigs0(1:neg,ng,0,iso)),ng=1,neg)
c       do ng=1,neg
c         if(ng.eq.1)then
c           write(*,'(1p9999e18.11)')sigt0(ng,iso),(0.d0,lg=ng+1,neg)
c         else if(ng.eq.neg)then
c           write(*,'(1p9999e18.11)')(0.d0,lg=1,ng-1),sigt0(ng,iso)
c         else
c           write(*,'(1p9999e18.11)')(0.d0,lg=1,ng-1),sigt0(ng,iso),
c    1       (0.d0,lg=ng+1,neg)
c         end if
c       end do ! ng
        write(*,'("sigs")')
        do ng=1,neg
c format for spreadsheet
c         write(*,'(1p9999e18.11)')(sigs0(ng,lg,0,iso),lg=1,neg)
c format for simple_ace
          write(*,'(1p9999e18.11)')(sigs0(lg,ng,0,iso),lg=1,neg)
        end do ! ng
        write(*,'("chi")')
        do ng=1,neg
          write(*,'(1p9999e18.11)')chi0(ng,ng,iso)
        end do ! ng
      end if
c compare cross sections.
c use avg to skip; 0 to run.
      if(np.eq.999)then
        maxdif_sigt(1)=1.d0
        maxdif_sigt(2)=-1.d0
        maxdif_nusigf(1)=1.d0
        maxdif_nusigf(2)=-1.d0
        maxdif_sigs(1)=1.d0
        maxdif_sigs(2)=-1.d0
        write(*,'(/,"material cross section comparison")')
        do i=1,nm
          do ng=1,neg
            diff=(sigt0(ng,i)-sigt(ng,i))/sigt(ng,i)
            if(diff.lt.0.d0)then
              maxdif_sigt(1)=min(maxdif_sigt(1),diff)
            else if(diff.gt.0.d0)then
              maxdif_sigt(2)=max(maxdif_sigt(2),diff)
            end if
            write(*,'("  sigt",2i5,1p3e14.6)')i,ng,sigt0(ng,i),
     1       sigt(ng,i),diff
            if(abs(diff).gt.1.d-6)then
              write(*,'(" **** sigt ****")')
            end if
            if(nusigf(ng,i).ne.0.d0)then
              diff=(nusigf0(ng,i)-nusigf(ng,i))/nusigf(ng,i)
              if(diff.lt.0.d0)then
                maxdif_nusigf(1)=min(maxdif_nusigf(1),diff)
              else if(diff.gt.0.d0)then
                maxdif_nusigf(2)=max(maxdif_nusigf(2),diff)
              end if
              write(*,'("nusigf",2i5,1p3e14.6)')i,ng,nusigf0(ng,i),
     1         nusigf(ng,i),diff
              if(abs(diff).gt.1.d-6)then
                write(*,'(" **** nusigf ****")')
              end if
            else
              if(nusigf0(ng,i).ne.0.d0)then
                write(*,'(" **** nusigf0.ne.0, nusigf.eq.0 ****")')
              end if
            end if
            do lg=1,neg
              do l=0,isct
                if(sigs(lg,ng,l,i).ne.0.d0)then
                 diff=(sigs0(lg,ng,l,i)-sigs(lg,ng,l,i))/sigs(lg,ng,l,i)
                 if(diff.lt.0.d0)then
                   maxdif_sigs(1)=min(maxdif_sigs(1),diff)
                 else if(diff.gt.0.d0)then
                   maxdif_sigs(2)=max(maxdif_sigs(2),diff)
                 end if
                 write(*,'(4i5,1p3e14.6)')i,ng,lg,l,sigs0(lg,ng,l,i),
     1            sigs(lg,ng,l,i),diff
                 if(abs(diff).gt.1.d-6)then
                   write(*,'(" **** sigs ****")')
                 end if
                else
                  if(sigs0(lg,ng,l,i).ne.0.d0)then
                    write(*,'(" **** sigs0.ne.0, sigs.eq.0 ****")')
                  end if
                end if
              end do ! l
            end do ! lg
          end do ! ng
        end do ! i
        write(*,'("max diffs, neg and pos")')
        write(*,'("sigt  ",1p2e14.6)')maxdif_sigt(1:2)
        write(*,'("nusigf",1p2e14.6)')maxdif_nusigf(1:2)
        write(*,'("sigs  ",1p2e14.6)')maxdif_sigs(1:2)
      end if
      end do ! np
      end do ! kl
      end do ! kh
      end do ! kg
      end do ! ki
      end do ! km
c
      return
      end subroutine
