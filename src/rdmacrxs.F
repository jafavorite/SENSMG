      subroutine rdmacrxs(iw,temps)
c read binary macrxs file from partisn for macroscopic cross sections.
c to recover cross sections from "key start mac" in for_out,
c multiply values by material density and divide scattering by (2l+1)
      use COMS
      implicit none
c iw=1/2 "for/macrxs" (nm)/"xs1/macrxs" (niso)
c maxord is the largest scattering expansion ever expected.
c maxxs is the largest number of principal+edit cross sections.
c maxmult is the largest number of fission multiplicities.
      integer iw,nof,maxord,maxxs,maxmult,nwd
      parameter (maxord=8,maxxs=10,maxmult=20,nwd=20)
      real*8 emax(neg),emin,c(nxs,maxxs),scat(neg*nxs*maxord),
     1 fiss_wgt(nxs),fiss_mult(maxmult,neg,nxs),sf_spectra(neg,nxs),
     2 upgam,temps,eavg(neg),einv(neg),wgts(neg),wgts_bnd(neg+1),
     3 fiss_mat(neg,neg,nxs),sourcm(neg,nxs)
      character hname*8,huse(2)*8,hmat(nxs)*8,hed(maxxs)*8
      integer*4 ivers,ia(nwd),nedt,ngpb(maxord,nxs),ifsg(maxord,nxs)
      integer ng,lg,i,if0,j,l,n,ntab
      character afile*10
      real*8 rl(0:isct)
c        
c nof is an offset in the cross-section arrays.
c they have materials first, then isotopes (including edit isotopes).
c this logic reads materials and edit isotopes in for/macrxs,
c then transport isotopes and edit isotopes in xs1/macrxs.
c in the second read, the transport isotopes overwrite the first
c read of the edit isotopes.
      if(iw.eq.1)then
        afile="for/macrxs"
        nof=0
      else if(iw.eq.2)then
        afile="xs1/macrxs"
        nof=nm
      end if
      write(*,'("reading multigroup cross sections from ",a,"...")')
     1 trim(afile)
      write(iuo,'("reading multigroup cross sections from ",a,"...")')
     1 trim(afile)
      open(iux,file=afile,status='old',position='rewind',
     1 form='unformatted',access='sequential',
     2 err=900)
      read(iux,err=910)hname,(huse(i),i=1,2),ivers
      ia(1:nwd)=0
      read(iux,err=910)(ia(i),i=1,nwd)
      if(idbgw.eq.1)then
        write(41,'(a)')trim(afile)
        write(41,'(10i8)')(ia(i),i=1,nwd)
      endif
      maxup=ia(7)
      nedt=ia(4)+ia(9)
      if(ia(1).ne.neg)go to 795
      if(ilnk3dnt.eq.0)then
        if(iw.eq.1)then
          if(ia(2).ne.nm+nrrx)go to 800
        else if(iw.eq.2)then
          if(ia(2).ne.niso)go to 802
        end if
      else if(ilnk3dnt.eq.1)then
        if(iw.eq.1)then
          if(ia(2).ne.nel)go to 803
        else if(iw.eq.2)then
          if(ia(2).ne.niso)go to 802
        end if
      else if(ilnk3dnt.eq.2)then
        if(iw.eq.1)then
          if(ia(2).ne.nm+nrrx)go to 800
        else if(iw.eq.2)then
          if(ia(2).ne.niso)go to 802
        end if
      end if
      if(ia(3).gt.maxord)go to 805
      if(ia(6).ne.lng)go to 810
      if(nedt.gt.maxxs)go to 815
      if(ia(12).gt.maxmult)go to 820
c ntemp must be 1. not set up to read or use more than 1.
      if(ia(14).ne.1)go to 825
      write(iuo,'("  fissdata=",i2)')ia(13)
      if0=fissdata
c how nochimat maps to fissdata
      if(partisnver(1:1).eq."5")if0=if0-1
      if(ia(13).ne.if0)then
        write(*,'("** warning. fissdata from ",a," not equal to ",
     1   "fissdata set in code = ",i1,".",/,
     2   "  resetting to match.")')trim(afile),fissdata
        write(iuo,'("** warning. fissdata from ",a," not equal to ",
     1   "fissdata set in code = ",i1,".",/,
     2   "  resetting to match.")')trim(afile),fissdata
        fissdata=ia(13)
        if(partisnver(1:1).eq."5")fissdata=fissdata+1
      end if
c
c 2l+1 is always included in macrxs. remove it.
      do l=0,isct
        rl(l)=real(2*l+1,8)
      end do
c
      if(partisnver(1:1).ne."5")then
        read(iux,err=910)(hmat(i),i=1,ia(2)),(hed(j),j=1,nedt),
     1   (vel(n),n=1,ia(1)),(emax(n),n=1,ia(1)),emin,temps,
     2   (eavg(n),n=1,ia(1)),(einv(n),n=1,ia(1)),(wgts(n),n=1,ia(1)),
     3   (wgts_bnd(n),n=1,ia(1)+1)
        ebins(1,1:ia(1)-1)=emax(2:ia(1))
        ebins(1,ia(1))=emin
        ebins(2,1:ia(1))=emax(1:ia(1))
        ebins(3,1:ia(1))=eavg(1:ia(1))
      else
        read(iux,err=910)(hmat(i),i=1,ia(2)),(hed(j),j=1,nedt),
     1   (vel(n),n=1,ia(1)),(emax(n),n=1,ia(1)),emin,temps
        ebins(1,1:ia(1)-1)=emax(2:ia(1))
        ebins(1,ia(1))=emin
        ebins(2,1:ia(1))=emax(1:ia(1))
        do ng=1,neg
          ebins(3,ng)=0.5d0*(ebins(1,ng)+ebins(2,ng))
        end do ! ng
      end if
c kynea3 energies are in keV; convert to MeV.
c convert velocities from cm/(0.01 sh) to cm/sh.
      if(libname(1:6).eq."kynea3")then
        ebins(1:3,1:neg)=ebins(1:3,1:neg)*1.d-3
        vel(1:neg)=vel(1:neg)*0.01d0
      end if
c kynea3 and scale have no eavg.
      if(libname(1:6).eq."kynea3".or.libname(1:5).eq."scale")then
        do ng=1,neg
          ebins(3,ng)=0.5d0*(ebins(1,ng)+ebins(2,ng))
        end do ! ng
      end if
      do ng=1,ia(1)
        read(iux,err=910)((c(i,j),i=1,ia(2)),j=1,nedt)
        do n=1,ia(2)
          chivec(ng,nof+n)=c(n,1)
          nusigf(ng,nof+n)=c(n,2)
          sigt(ng,nof+n)=c(n,3)
          siga(ng,nof+n)=c(n,4)
        end do ! n
        if(ia(3).gt.0)then
          read(iux,err=910)((ngpb(l,j),l=1,ia(3)),j=1,ia(2)),
     1     ((ifsg(l,j),l=1,ia(3)),j=1,ia(2))
          ntab=sum(ngpb(1:ia(3),1:ia(2)))
          if(ilnk3dnt.eq.0)then
            if(iw.eq.1)then
              if(ntab.gt.maxord*neg*(nm+nrrx))go to 930
            else if(iw.eq.1)then
              if(ntab.gt.maxord*neg*niso)go to 932
            end if
          else if(ilnk3dnt.gt.0)then
            if(iw.eq.1)then
              if(ntab.gt.maxord*neg*nel)go to 933
            else if(iw.eq.1)then
              if(ntab.gt.maxord*neg*niso)go to 932
            end if
          end if
          read(iux,err=920)(scat(i),i=1,ntab)
c separate the bands.
c sigs(ng,lg,l,ir) means scattering from lg into ng.
c scgr(ng,l,ir,1/2) are the top and bottom groups scattering into ng.
          i=0
          do j=1,ia(2)
            do l=1,ia(3)
              if(l-1.le.isct)then
                scgr(ng,l-1,nof+j,1)=ifsg(l,j)-ngpb(l,j)+1
                scgr(ng,l-1,nof+j,2)=ifsg(l,j)
              end if
              do lg=ifsg(l,j),ifsg(l,j)-ngpb(l,j)+1,-1
                i=i+1
                if(l-1.le.isct)sigs(ng,lg,l-1,nof+j)=scat(i)
     1           /rl(l-1)
              end do ! lg
            end do ! l
          end do ! j
        end if
      end do ! ng
c
      if(ia(11).ne.0)read(iux,err=940)(fiss_wgt(nof+i),i=1,ia(2))
      if(ia(12).ne.0)read(iux,err=940)(((fiss_mult(i,ng,nof+n),
     1 i=1,ia(12)),ng=1,ia(1)),n=1,ia(2))
      if(ia(15).gt.0)read(iux,err=940)((sf_spectra(n,nof+i),n=1,ia(1)),
     1 i=1,ia(2))
      if(ia(17).gt.0)read(iux,err=940)((sourcm(n,nof+i),n=1,ia(1)),
     1 i=1,ia(2))
      read(iux,err=940)upgam
      if0=0
      if(partisnver(1:1).eq."5")then
        if(ia(13).eq.0)if0=1
      else
        if(ia(13).le.1)if0=1
      end if
      if(if0.eq.1)then
c chi(ng,lg,nm) is from lg into ng (like sigs).
c DEBUG. why does reading this take so long?
        write(*,'("  rdmacrxs, reading chi...")')
        read(iux,err=940)(((chi(ng,lg,nof+n),n=1,ia(2)),lg=1,ia(1)),
     1   ng=1,ia(1))
      end if
c     do lg=1,lng
c       write(41,'("group",i4)')lg
c       write(41,'(1p6e14.4)')((chi(lg,ng,nof+n),n=1,ia(2)),ng=1,ia(1))
c       write(41,'(1p30e14.4)')((chi(lg,ng,nof+7)),ng=1,ia(1))
c     end do ! lg
c     close(41)
      if(partisnver(1:1).ne."5".and.ia(13).eq.0)then
c fiss_mat(ng,lg,nm) is from lg into ng (like sigs).
c DEBUG. why does reading this take so long?
        write(*,'("  rdmacrxs, reading fiss_mat...")')
        read(iux,err=940)(((fiss_mat(ng,lg,nof+n),n=1,ia(2)),
     1   lg=1,ia(1)),ng=1,ia(1))
        write(*,'("  done.")')
c once I thought fiss_mat should be multiplied by rho. I don't know
c why. it does not seem to be true.
c       if(iw.eq.1)then
c nof=0, ia(2)=nm+nrrx for this case
c         do n=1,nm
c           fiss_mat(1:ia(1),1:ia(1),n)=fiss_mat(1:ia(1),1:ia(1),n)
c    1       *rho(n)
c         end do ! n
c       end if
c       do lg=1,lng
c         write(41,'("group",i4)')lg
c         write(41,'(1p6e14.4)')((fiss_mat(lg,ng,nof+n),n=1,ia(2)),
c    1     ng=1,ia(1))
c       end do ! lg
        go to 50
        write(*,'("** constructing chi matrix from fission ",
     1   "transfer matrix.")')
        write(iuo,'("** constructing chi matrix from fission ",
     1   "transfer matrix.")')
        chi(1:neg,1:neg,nof+1:nof+ia(2))=0.d0
        do n=1,ia(2)
          do ng=1,neg
            if(nusigf(ng,nof+n).gt.0.)chi(1:neg,ng,nof+n)=
     1       fiss_mat(1:neg,ng,nof+n)/nusigf(ng,nof+n)
          end do ! ng
        end do ! n
      end if
   50 if0=0
      if(partisnver(1:1).eq."5")then
        if(ia(13).eq.1)if0=1
      else
        if(ia(13).eq.2)if0=1
      end if
      if(if0.eq.1)then
c populate matrix with vector chi.
c from reading of principal cross sections, chivec(ng) is into 
c group ng from all groups. chi(ng,lg) is into group ng from
c group lg. set to chivec(ng) for all groups.
        do ng=1,neg
          do lg=1,neg
            chi(ng,lg,nof+1:nof+ia(2))=chivec(ng,nof+1:nof+ia(2))
          end do ! ng
        end do ! ng
      end if
      close(iux)
      write(*,'("done.")')
      goto 999
c
  795 write(*,'("error. ia(1).ne.neg in ",a,"; ia(1)=",i8)')
     1 trim(afile),ia(1)
      write(*,'("neg=",i8)')neg
      call stoponerror
  800 write(*,'("error. ia(2).ne.nm+nrrx in ",a,"; ia(2)=",i8)')
     1 trim(afile),ia(2)
      write(*,'("nm+nrrx=",i8)')nm+nrrx
      call stoponerror
  802 write(*,'("error. ia(2).ne.niso in ",a,"; ia(2)=",i8)')
     1 trim(afile),ia(2)
      write(*,'("niso=",i8)')niso
      call stoponerror
  803 write(*,'("error. ia(2).ne.nel in ",a,"; ia(2)=",i8)')
     1 trim(afile),ia(2)
      write(*,'("nel=",i8)')nel
      call stoponerror
  804 write(*,'("error. ia(2).ne.nm in ",a,"; ia(2)=",i8)')
     1 trim(afile),ia(2)
      write(*,'("nm=",i8)')nm
      call stoponerror
  805 write(*,'("error. ia(3).gt.maxord in ",a,"; ia(3)=",i8)')
     1 trim(afile),ia(3)
      write(*,'("maxord=",i8)')maxord
      write(*,'("increase maxord in rdmacrxs.")')
      call stoponerror
  810 write(*,'("error. ia(6).ne.lng in ",a,"; ia(6)=",i8)')
     1 trim(afile),ia(6)
      write(*,'("lng=",i8)')lng
      call stoponerror
  815 write(*,'("error. nedt.gt.maxxs in ",a,"; nedt=",i8)')
     1 trim(afile),nedt
      write(*,'("maxxs=",i8)')maxxs
      write(*,'("increase maxxs in rdmacrxs.")')
      call stoponerror
  820 write(*,'("error. ia(12).gt.maxmult in ",a,"; ia(12)=",i8)')
     1 trim(afile),ia(12)
      write(*,'("maxmult=",i8)')maxmult
      write(*,'("increase maxmult in rdmacrxs.")')
      call stoponerror
  825 write(*,'("error. ntemp.ne.1 in ",a,"; ntemp=",i8)')
     1 trim(afile),ia(14)
      call stoponerror
  900 write(*,'("error opening ",a,".")')trim(afile)
      call stoponerror
  910 write(*,'("error reading ",a,".")')trim(afile)
      call stoponerror
  920 write(*,'("error reading ",a," for scattering.")')trim(afile)
      call stoponerror
  930 write(*,'("error reading ",a," for scattering: ntab too big.",
     1 " ntab, maxord*neg*(nm+nrrx)=",2i6)')trim(afile),
     2 ntab,maxord*neg*(nm+nrrx)
      write(*,'("increase maxord in rdmacrxs?")')
      call stoponerror
  932 write(*,'("error reading ",a," for scattering: ntab too big.",
     1 " ntab, maxord*neg*niso=",2i6)')trim(afile),
     2 ntab,maxord*neg*niso
      write(*,'("increase maxord in rdmacrxs?")')
      call stoponerror
  933 write(*,'("error reading ",a," for scattering: ntab too big.",
     1 " ntab, maxord*neg*nel=",2i6)')trim(afile),
     2 ntab,maxord*neg*nel
      write(*,'("increase maxord in rdmacrxs?")')
      call stoponerror
  940 write(*,'("error reading ",a," for chi matrix.")')trim(afile)
      call stoponerror
c
  999 return
      end subroutine
