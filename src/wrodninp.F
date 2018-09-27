      subroutine wrodninp(sensverify)
c write odninp, cross sections within input file.
c this routine is for perturbed cross sections but can be used
c for unperturbed also.
      use COMS
      implicit none
      real*8 sensverify(neg+1,neg+1,0:isct,nel)
      real*8 sigt0(neg,nm),nusigf0(neg,nm),siga0(neg,nm),
     1 chi0(neg,neg,nm),sigs0(neg,neg,0:isct,nm)
      real*8 xsa(neg),xsn(neg),xst(neg),xss(neg,2*neg-1),adens,
     1 pertf,pertn,perta,perts,pertt,pertd,ppf,ppn,ppa,pps,ppt,pptadd,
     2 diff,maxdif_sigt(2),maxdif_nusigf(2),maxdif_sigs(2),pp,
     3 pertc,ppc,pert,chinum(neg),chiden,nusigf0j(nel),unpert
      integer k1,k2,iso,i,iitm,i_srcacc,j,kg,kh,kl,km,ki,l,lg,mt,
     1 nf,ng,n0,np,igrp,igrp_in,imt,iis,il
      character pfile*120,title*7,ct*3,typel*100,ctmp*9
      integer iht,ihs,ihm,ith,ifido
      logical np0done
      iitm=999   ! if iitm.gt.0 it will be used for forward and adjoint.
      i_srcacc=0 ! 0/1 for default/srcacc=no
      ifido=-1   ! -1/0 for 4e18/6e12
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
      do ki=1,2 !ncb(km)
      k1=k1+1
      write(*,'(" isotope",i4," in material",i4)')ki,km
      write(iuo,'(" isotope",i4," in material",i4)')ki,km
      do kg=1,neg+1
      do kh=1,1 !neg  ! for scattering
      do kl=0,0 !isct ! for scattering
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
      if(sensverify(kg,kh,kl,k1).eq.0.d0)cycle
      if(np.eq.-1)then
c       pertt=abs(0.001d0/sensverify(kg,1,0,k1))
c       pertt=min(0.1d0,pertt)
c       pert=pertt
c       perts=abs(0.001d0/sensverify(kg,kh,kl,k1))
c       pert=perts
c       pertd=abs(0.001d0/sensverify(kg,1,0,k1))
c       pertd=min(0.999d0,pertd)
c       pert=pertd
c       pertf=abs(0.001d0/sensverify(kg,1,0,k1))
c       pertf=min(0.999d0,pertf)
c       pert=pertf
        pertn=abs(0.001d0/sensverify(kg,1,0,k1))
        pertn=min(0.9d0,pertn)
        pert=pertn
c       pertc=abs(0.01d0/sensverify(kg,1,0,k1))
c       pertc=min(0.999d0,pertc)
c       pert=pertc
      else if(np.ge.0)then
c       pertt=pert
c       perts=pert
c       pertd=pert
c       pertf=pert
        pertn=pert
c       pertc=pert
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
      if(np.eq.0)pp=0.d0
c
c zero the cross sections.
      nusigf0(1:neg,1:nm)=0.d0
      sigt0(1:neg,1:nm)=0.d0
      siga0(1:neg,1:nm)=0.d0
      sigs0(1:neg,1:neg,0:isct,1:nm)=0.d0
      chi0(1:neg,1:neg,1:nm)=0.d0
      nusigf0j(1:nel)=0.d0
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
            ppn=1.d0
            if(np.ne.0.and.(ng.eq.igrp.or.igrp.eq.neg+1).and.i.eq.imt
     1       .and.j.eq.iis.and.pertn.gt.0.d0)then
              ppn=1.d0+sign(pertn,real(np,8))
              if(sigf(ng,iso).eq.0.d0)then
                write(*,'("error. sig=0. probably skip reading c in ",
     1           "rdsnxedt.")')
                call stoponerror
              end if
              unpert=nusigf(ng,iso)/sigf(ng,iso)
            end if
            ppf=1.d0
            if(np.ne.0.and.(ng.eq.igrp.or.igrp.eq.neg+1).and.i.eq.imt
     1       .and.j.eq.iis.and.pertf.gt.0.d0)then
              ppf=1.d0+sign(pertf,real(np,8))
              pptadd=pptadd+sigf(ng,iso)*adens*(ppf-1.d0)
              unpert=sigf(ng,iso)
            end if
c only one of ppn or ppf should be > 1.
            nusigf0(ng,i)=nusigf0(ng,i)+nusigf(ng,iso)*adens
     1       *ppn*ppf
            nusigf0j(k2)=nusigf0j(k2)+nusigf(ng,iso)*adens
     1       *ppn*ppf
c           if(np.eq.-1)write(*,'("C",4i3,1p8e15.7)')k1,kg,k2,ng,
c    1       pertf,ppf,nusigf0(ng,i),nusigf(ng,iso),adens,ppn,
c    2       sigf(ng,i),sigf(ng,iso)
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
                    pptadd=pptadd+sigs(lg,ng,l,iso)*adens
     1               *(pps-1.d0)
                  end if
                  unpert=sigs(lg,ng,l,iso)
                end if
                sigs0(lg,ng,l,i)=sigs0(lg,ng,l,i)+sigs(lg,ng,l,iso)
     1           *adens*pps
              end do ! l
            end do ! lg
            ppa=1.d0
            if(np.ne.0.and.(ng.eq.igrp.or.igrp.eq.neg+1).and.i.eq.imt
     1       .and.j.eq.iis.and.perta.gt.0.d0)then
              ppa=1.d0+sign(perta,real(np,8))
              pptadd=pptadd+siga(ng,iso)*adens*(ppa-1.d0)
              unpert=siga(ng,iso)
            end if
            siga0(ng,i)=siga0(ng,i)+siga(ng,iso)*adens*ppa
            ppc=1.d0
            if(np.ne.0.and.(ng.eq.igrp.or.igrp.eq.neg+1).and.i.eq.imt
     1       .and.j.eq.iis.and.pertc.gt.0.d0)then
              ppc=1.d0+sign(pertc,real(np,8))
              pptadd=pptadd+sigc(ng,iso)*adens*(ppc-1.d0)
              unpert=sigc(ng,iso)
            end if
            ppt=1.d0
            if(np.ne.0.and.(ng.eq.igrp.or.igrp.eq.neg+1).and.i.eq.imt
     1       .and.j.eq.iis.and.pertt.gt.0.d0)then
              ppt=1.d0+sign(pertt,real(np,8))
              unpert=sigt(ng,iso)
            end if
            sigt0(ng,i)=sigt0(ng,i)+sigt(ng,iso)*adens*ppt+pptadd
c check P0 scattering and total for isotope.
c           if(sigt(ng,iso)*adens*ppt+pptadd.lt.
c    1       sigs(ng,ng,0,iso)*adens*pps)then
c             write(*,'("p0 self-scatter greater than total; ",
c    1         "iso,ng,sigt,sigs=",2i4,1p2e14.6)')iso,ng,
c    2         sigt(ng,iso)*adens*ppt+pptadd,sigs(ng,ng,0,iso)*adens*pps
c           end if
          end do ! ng
        end do ! j
      end do ! i
c set chi
      k2=0
      do i=1,nm
        chinum(1:neg)=0.d0
        chiden=0.d0
        do j=1,ncb(i)
          k2=k2+1
          iso=nm+k2
          do ng=1,neg
            chinum(ng)=chinum(ng)+chi(ng,ng,iso)*nusigf0j(k2)
          end do ! ng
c chiden is equal to sum(nusigf0(1:neg,i))
          chiden=chiden+nusigf0j(k2)
        end do ! j
        if(chiden.eq.0.d0)cycle
        do ng=1,neg
          chi0(ng,ng,i)=chinum(ng)/chiden
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
c upscattering.  iht=3 always.
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
     1 .or.ct.eq."cap")then
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
      nf=iup
      open(nf,file=pfile,status='unknown')
      write(nf,'("     3     0     0")')
      write(nf,'(a)')trim(id)
      if(icalc.eq.0.or.icalc.eq.3.or.icalc.eq.4)then
        write(nf,'(1x,a7," input file, fixed-source")')title
      else if(icalc.eq.1)then
        write(nf,'(1x,a7," input file, keff")')title
      else if(icalc.eq.2)then
        write(nf,'(1x,a7," input file, alpha")')title
      end if
      write(nf,'(a)')trim(typel)
c
      mt=nm
      n0=nm
      call block1(nf,mt,n0)
      call block2(nf)
      write(nf,'("/ * * * * block iii * * * *")')
      write(nf,'(" lib=odninp")')
      write(nf,'(" iht=3 ihs=",i4," ihm=",i4)')ihs,ihm
      write(nf,'(" ifido=",i2," ititl=1")')ifido
      write(nf,'(" maxord=",i2)')isct
      write(nf,'(" names= i01")')
      do i=2,nm
        write(nf,'("        i",i2.2)')i
      end do ! i
      write(nf,'(" lng=",i3)')lng
      if(icalc.eq.2)then
        write(nf,'(" vel=")')
        write(nf,'(1p,4e18.11)')(vel(ng),ng=1,neg)
      end if
      write(nf,'(" t")')
c
c write cross sections.
      do i=1,nm
        xsa(1:neg)=siga0(1:neg,i)
        xsn(1:neg)=nusigf0(1:neg,i)
        xst(1:neg)=sigt0(1:neg,i)
        do l=0,isct
          if(ifido.eq.0)then
            write(nf,'(1x,"siga",8x,"nusigf",6x,"sigt",
     1       9x,"mat,ord=",2i3)')i,l
          elseif(ifido.eq.-1)then
            write(nf,'(1x,"siga",14x,"nusigf",12x,"sigt",
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
            write(nf,'(1p,6e12.5 )')(xsa(ng),xsn(ng),
     1       xst(ng),(xss(ng,lg),lg=1,2*neg-1),ng=1,neg)
          elseif(ifido.eq.-1)then
            write(nf,'(1p,4e18.11)')(xsa(ng),xsn(ng),
     1       xst(ng),(xss(ng,lg),lg=1,2*neg-1),ng=1,neg)
          end if
          if(l.gt.0)cycle
          xsa(1:neg)=0.d0
          xsn(1:neg)=0.d0
          xst(1:neg)=0.d0
        end do ! l
      end do ! i
c
      write(nf,'("/ * * * * block iv * * * *")')
      write(nf,'(" matspec=atdens")')
      write(nf,'(" matls= m01   i01  1.;")')
      do i=2,nm
        write(nf,'("        ","m",i2.2,"   i",i2.2,"  1.;")')i,i
      end do ! i
      write(nf,'(" assign= zone01   m01  1.;")')
      do i=2,nm
        write(nf,'("         zone",i2.2,"   m",i2.2,"  1.;")')i,i
      end do ! i
      write(nf,'(" zonetemp= 2.5301E-08 ;")')
      do i=2,nm
      write(nf,'("           2.5301E-08 ;")')
      end do ! i
      write(nf,'(" t")')
c
      write(nf,'("/ * * * * block v * * * *")')
      if(icalc.eq.0.or.icalc.eq.3.or.icalc.eq.4)then
        write(nf,'(" ievt=",i1," isct=",i1)')0,isct
      else
        write(nf,'(" ievt=",i1," isct=",i1)')icalc,isct
      end if
      write(nf,'(" epsi=",1pe9.2," balp=1")')1.d-8
      if(icalc.eq.1.or.icalc.eq.2)then
        write(nf,'(" norm=1.0")')
      end if
      write(nf,'(" npeg=2")')
      if(idbgw.eq.1)write(nf,'(" fluxp=2")')
c     write(nf,'(" fluxp=2")')
      if(iitm.gt.0)write(nf,'(" iitm=",i3)')iitm
      write(nf,'(" iitl=0 oitm=9999")')
      if(i_srcacc.eq.1)write(nf,'(" srcacc=no")')
      if(nofxup.eq.1)write(nf,'(" nofxup=1")')
c
c different solve parameters possible for forward and adjoint files 
c forward:
      if(ith.eq.0)then
        write(nf,'(" ith=0 xsectp=2")')
c     
c adjoint:
      else if(ith.eq.1)then
        write(nf,'(" ith=1 xsectp=0")')
      end if
      write(nf,'(" chi=")')
      do i=1,nm
        write(nf,'(2x,1p5e14.6)')(chi0(ng,ng,i),ng=1,neg)
        write(nf,'(2x,";")')
      end do ! i
      if(icalc.eq.0.or.icalc.eq.3.or.icalc.eq.4)then
        call fixedsrc(nf,ith)
      end if
      write(nf,'(" t")')
c no edit for adjoint.
      if(ith.eq.1)go to 100
      write(nf,'("/ * * * * block vi * * * *")')
      write(nf,'(" ajed=",i1)')ith
      write(nf,'(" zned=1 igrped=0")')
      write(nf,'(" rsfnam=")')
      do i=1,nrrx
        write(nf,'(2x,1h","m",i4.4,1h"," ;")')nm+i
      end do ! i
      write(nf,'(2x,1h","flux",1h"," ;")')
      write(nf,'(" rsfe=")')
      do i=1,nrrx
        write(nf,'(2x,1p5e14.6)')(rrxs(ng,i),ng=1,neg)
        write(nf,'(2x,";")')
      end do ! i
      write(nf,'(2xi3,"r 1. ;")')neg
      write(nf,'(" t")')
c
  100 close(nf)
c
c compare cross sections.
      cycle ! skip all this
      if(np.eq.0)then
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
