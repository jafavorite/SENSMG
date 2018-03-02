      subroutine rdsources
c read sources4c output files for neutron source.
c the source is in n/s-cm^3.
      use COMS
      implicit none
      real*8 nsrcfan(neg,nm),nsrcfsf(neg,nm),sf(nel),
     1 sant(nel),sans(nel),totns(neg,nm),st,siniso(nel),
     2 sanisol(nel),sfisol(nel),an_norm,sftot
      integer i,i1,i2,iz,j,k,k2,ios,ng,nant,nans,nsf,isponfis
      character ofile*20,line*120,l1*1,csf(nel)*5,cant(nel)*5,
     1 cans(nel)*5,cc*2,ct*5,cs*5
c in feyny, these are used in multmom. not sure if they will be needed here.
c     real*8 srcnuis(nel,nm)
c in feyny, these are used in asymptote. not sure if they will be needed here.
c     real*8 nsrcsf(nm)
c isponfis=0/1 keep/discard spontaneous fission source for testing
c (alpha,n) stuff
      isponfis=0
c
      if(ialphan.eq.0)then
        write(*,'("warning. (alpha,n) sources are ignored.")')
        write(iuo,'("warning. (alpha,n) sources are ignored.")')
      end if
      if(isponfis.eq.1)then
        write(*,'("warning. spontaneous fission sources are ",
     1   "ignored.")')
        write(iuo,'("warning. spontaneous fission sources are ",
     1   "ignored.")')
      end if
c
c totns is read but not used for anything.
      totns(1:neg,1:nm)=0.d0
      nsrcfan(1:neg,1:nm)=0.d0
      nsrcfsf(1:neg,1:nm)=0.d0
      siniso(1:nel)=0.d0
      sanisol(1:nel)=0.d0 ! l is for local (this routine)
      sfisol(1:nel)=0.d0 ! l is for local (this routine)
      k2=-ncb(1)
      if(ictrl.eq.2)then
        if(imisc.eq.0.and.ialphan.eq.0)then
          write(iuo,'(/,"spontaneous fission source from sources4c, ",
     1     "(alpha,n) source ignored")')
        else if(imisc.eq.0.and.ialphan.eq.1)then
          write(iuo,'(/,"spontaneous fission source from sources4c, ",
     1     "(alpha,n) source from sources4c")')
        else if(imisc.eq.1.and.ialphan.eq.1)then
          write(iuo,'(/,"spontaneous fission source from misc, ",
     1     "(alpha,n) source from sources4c")')
        end if
        write(iuo,'("output from sources4c")')
      end if
      do i=1,nm
        k2=k2+ncb(i)
        sant(1:nel)=0.d0
        sans(1:nel)=0.d0
        sf(1:nel)=0.d0
c spectrum file
        write(ofile,'("sources/m",i3.3,"_outp2")')i
        open(iun,file=ofile,status='old',iostat=ios)
        if(ios.ne.0)then
          write(*,'("error. no sources4c output file ",a)')trim(ofile)
          write(iuo,'("error. no sources4c output file ",a)')trim(ofile)
          call stoponerror
        end if
        read(iun,'(a)',iostat=ios)line
c an empty file because no sources.
        if(ios.ne.0)then
          if(ictrl.eq.2)then
            write(*,'(/,2x,"no sources4c neutron source in material",
     1       i3)')i
            write(iuo,'(/,2x,"no sources4c neutron source in material",
     1       i3)')i
          end if
          close(iun)
          cycle
        end if
c or, with mods by JAF, there is NaN in the file when no sources.
        do j=1,13
          read(iun,'(a)')line
        end do ! j
        read(iun,'(a)')line
c more digits in sources4c.
c       if(line(19:21).eq."NaN")then
        if(line(25:27).eq."NaN")then
          if(ictrl.eq.2)then
            write(*,'(/,2x,"no sources4c neutron source in material",
     1       i3)')i
            write(iuo,'(/,2x,"no sources4c neutron source in material",
     1       i3)')i
          end if
          close(iun)
          cycle
        end if
c a source material.
        ismat(i)=1
c (alpha,n) sources are not normalized correctly. the total in
c outp2 is right, but the spectrum is wrong. renormalize by the
c "normalized total".
c the file is already open to the first energy group.
        if(ialphan.eq.1)then
          do ng=1,neg
            read(iun,'(a)')line
          end do ! ng 
          read(iun,'(13x,e14.6)')an_norm
c identify an (alpha,n) source material.
          if(an_norm.gt.0.d0)then
            isan(i)=1
          end if
        end if
   20   read(iun,'(a)')line
        if(line(13:27).ne."Absolute Totals")go to 20
        read(iun,'(a)')line ! ==
        read(iun,'(a)')line ! heading
        read(iun,'(a)')line ! --
        do ng=1,neg
c read (alpha,n) and spont. fiss. sources and total
c more digits in sources4c.
c         read(iun,'(10x,2e11.3,11x,e11.3)')nsrcfan(ng,i),nsrcfsf(ng,i),
c    1     totns(ng,i)
          read(iun,'(13x,2e14.6,14x,e14.6)')nsrcfan(ng,i),nsrcfsf(ng,i),
     1     totns(ng,i)
c renormalize due to sources4c weirdness.
          if(isan(i).eq.1)nsrcfan(ng,i)=nsrcfan(ng,i)/an_norm
          if(nsrcfsf(ng,i).lt.0.d0)then
            nsrcfsf(ng,i)=0.d0
            totns(ng,i)=nsrcfan(ng,i)
          end if
          if(totns(ng,i).lt.0.d0)then
            totns(ng,i)=0.d0
          end if
c turn off spont. fiss. for testing (alpha,n) stuff.
          if(isponfis.eq.1)then
            nsrcfsf(ng,i)=0.d0
          end if
c nsrcf is used in wrdantnm and calcsens_d; nsrc is the group total.
          if(imisc.eq.0)then
            nsrcf(ng,i)=nsrcfsf(ng,i) ! s.f.
            nsrc(i)=nsrc(i)+nsrcfsf(ng,i)
          end if
          if(ialphan.eq.1)then
            nsrcf(ng,i)=nsrcf(ng,i)+nsrcfan(ng,i)
            nsrc(i)=nsrc(i)+nsrcfan(ng,i)
c           write(*,'("src4.",i4,1p3e15.7)')ng,nsrc(i),nsrcfan(ng,i),
c    1       nsrcfsf(ng,i)
          end if
        end do ! ng
        close(iun)
        sftot=sum(nsrcfsf(1:neg,i)) ! total spont. fiss.
c debug prints
c       write(*,'("src4.",1p30e15.7)')nsrc(i),sum(nsrcfan(1:neg,i)),
c    1   sum(nsrcfsf(1:neg,i)),an_norm
c       write(iuo,'("debug in rdsources. ",
c    1   "(alpha,n) source from sources4c")')
c       do ng=1,neg
c         write(iuo,'(1pe14.6)')nsrcfan(ng,i)
c       end do ! ng
c       write(iuo,'(1pe14.6)')sum(nsrcfan(1:neg,i))
c
c magnitude file
        write(ofile,'("sources/m",i3.3,"_tape6")')i
        open(iun,file=ofile,status='old',iostat=ios)
        if(ios.ne.0)then
          write(*,'("error. no sources4c tape6 file ",a)')trim(ofile)
          write(iuo,'("error. no sources4c tape6 file ",a)')trim(ofile)
          call stoponerror
        end if
        nant=0
        nans=0
c skip (alpha,n) sources if none
        if(isan(i).eq.0)go to 50
   30   read(iun,'(a1)',iostat=ios)l1
        if(l1.ne."1")go to 30
c (alpha,n) sources
   32   read(iun,'(a1)',iostat=ios)l1
        if(l1.ne."+")go to 32
   40   read(iun,'(a)',iostat=ios)line
        if(line(49:59).eq."all targets")go to 50
        if(line(1:1).eq."2")go to 52 ! if there is only one target
        if(line(2:7).eq."failed")go to 40
        if(line(8:12).ne."     ")then
          read(line,'(7x,a5,15x,a5)')ct,cs
          i1=0
          do k=1,nant
            if(ct.eq.cant(k))then
              i1=k
              exit
            end if
          end do ! k
          if(i1.eq.0)then
            nant=nant+1
            i1=nant
            cant(i1)=ct
          end if
          i2=0
          do k=1,nans
            if(cs.eq.cans(k))then
              i2=k
              exit
            end if
          end do ! k
          if(i2.eq.0)then
            nans=nans+1
            i2=nans
            cans(i2)=cs
          end if
          backspace(iun)
   42     read(iun,'(a)',iostat=ios)line
          if(ios.ne.0)then
            write(*,'("error. end of file ",a)')trim(ofile)
            call stoponerror
          end if
          if(line(56:61).ne."Total:")go to 42
c more digits in sources4c.
c         read(line,'(65x,e12.4)')st
          read(line,'(65x,e14.6)')st
          sans(i2)=sans(i2)+st
          read(iun,'(a)',iostat=ios)line ! blank
        else if(line(49:59).eq."this target")then
c more digits in sources4c.
c         read(line,'(61x,e16.4)')sant(i1)
          read(line,'(61x,e18.6)')sant(i1)
        end if
        go to 40
c
   50   read(iun,'(a1)',iostat=ios)l1
        if(l1.ne."2")go to 50
c spontaneous fission sources
   52   read(iun,'(a1)',iostat=ios)l1
        if(l1.ne."+")go to 52
        nsf=0
   60   read(iun,'(a)',iostat=ios)line
        if(ios.ne.0)go to 70
        if(line(2:7).eq."failed")go to 60
        if(line(1:1).eq." ")then
          nsf=nsf+1
c more digits in sources4c.
c         read(line,'(9x,a5,46x,e11.3)')csf(nsf),sf(nsf)
          read(line,'(9x,a5,46x,e14.6)')csf(nsf),sf(nsf)
          go to 60
        end if
   70   close(iun)
c
        do k=1,nant
          cc=cant(k)(1:2)
          call convz(cc,iz)
          read(cant(k)(3:5),'(i3)')i1
          iz=iz+i1
          do j=1,ncb(i)
c 6013 is written in wrsources for 6000 and 6012.
            if((iz.eq.nint(blk(1,k2+j))).or.(iz.eq.6013.and.
     1       (nint(blk(1,k2+j)).eq.6000.or.nint(blk(1,k2+j)).eq.6012)))
     2       then
              if(ialphan.eq.1)then
                sanisol(k2+j)=sant(k)
                siniso(k2+j)=siniso(k2+j)+sant(k)
              end if
            end if
          end do ! j
        end do ! k
        do k=1,nans
          cc=cans(k)(1:2)
          call convz(cc,iz)
          read(cans(k)(3:5),'(i3)')i1
          iz=iz+i1
          do j=1,ncb(i)
            if(iz.eq.nint(blk(1,k2+j)))then
              if(ialphan.eq.1)then
                sanisol(k2+j)=sans(k)
                siniso(k2+j)=siniso(k2+j)+sans(k)
              end if
            end if
          end do ! j
        end do ! k
        do k=1,nsf
          cc=csf(k)(1:2)
          call convz(cc,iz)
          read(csf(k)(3:5),'(i3)')i1
          iz=iz+i1
          do j=1,ncb(i)
            if(iz.eq.nint(blk(1,k2+j)))then
              if(imisc.eq.0)then
c turn off spont. fiss. for testing (alpha,n) stuff.
                if(isponfis.eq.1)then
                  sf(k)=0.d0
                end if
                sfisol(k2+j)=sf(k)
                siniso(k2+j)=siniso(k2+j)+sf(k)
c apply the overall s.f. spectrum to all isotopes until sources4c writes
c that info.
                if(isponfis.eq.0)then
                  sfiso(1:neg,k2+j)=sf(k)*nsrcfsf(1:neg,i)/sftot
                end if
              else if(imisc.eq.1)then
                sf(k)=0.d0 ! so that the printed total is right.
              end if
            end if
          end do ! j
        end do ! k
c
        if(ictrl.eq.2)then
          write(iuo,'(/,"  material",i3,/,
     1     "  isotope   (alpha,n)    spont.fiss.  total")')i
          do j=1,ncb(i)
            write(iuo,'(i10,1p3e13.5)')nint(blk(1,k2+j)),
     1       sanisol(k2+j),sfisol(k2+j),siniso(k2+j)
          end do ! j
          write(iuo,'(2x,"totals [(alpha,n) targets and sources ",
     1     "should be equal]",/,2x,"(alpha,n) targets",1pe17.5,/,
     2     2x,"(alpha,n) sources",e17.5,/,
     3     2x,"spont.fiss. sources",e15.5,/,2x,"total",e29.5)')
     4     sum(sant(1:nant)),sum(sans(1:nans)),sum(sf(1:nsf)),
     5     sum(sant(1:nant))+sum(sf(1:nsf))
          if(imisc.eq.1)then
            write(iuo,'("  in this calculation, spont. fiss. ",
     1       "source is from misc, not sources4c.")')
          end if
          if(ialphan.eq.0)then
            write(iuo,'("  the (alpha,n) sources are ignored in this ",
     1       "calculation.")')
          end if
c
c redo this if a total neutron source rate for each isotope is desired.
c         write(iuo,'(/,"  material",i3,/,
c    1     "  isotope   total_neut_src_(n/cm^3/s)")')i
c a total over isotopes makes no sense because (alpha,n) targets
c and alpha sources are both included, so don't bother.
c         do j=1,ncb(i)
c           write(iuo,'(i10,1pe13.5)')nint(blk(1,k2+j)),siniso(k2+j)
c         end do ! j
        end if
      end do ! i
c
c for Brian's Trans. Am. Nucl. Soc. 109 (2013)
c     siniso(1:nel)=0.d0 
c     nsrc(1:nm)=0.d0
c     nsrcf(1:neg,1:nm)=0.d0
c     call set_ebins(neg,iug2,iuo,libname,ebins)
c     do ng=1,neg
c       if(ebins(1,ng).lt.1.d0.and.ebins(2,ng).ge.1.d0)then
c         nsrcf(ng,1)=1.d0/vol(1,1)
c         nsrc(1)=nsrc(1)+nsrcf(ng,1)
c         exit
c       end if
c     end do ! ng
c
      if(sum(ismat(1:nm)).eq.0)then
        write(*,'("error. no neutron sources in rdsources.")')
        write(iuo,'("error. no neutron sources in rdsources.")')
        call stoponerror
      end if
c
c write mcnp source for mat 1
c     call set_ebins(neg,iug2,iuo,libname,ebins)
c     i=1
c     write(*,'(3x,1p2e13.5)')ebins(1,neg),0.d0
c     do ng=neg,1,-1
c       write(*,'(3x,1p2e13.5)')ebins(2,ng),nsrcf(ng,i)
c     end do ! ng
c     call stoponerror
c
      return
      end subroutine
c
      subroutine convz(cc,iz)
      implicit none
      integer iz
      character cc*2
      integer ncz,i
      parameter (ncz=109)
      character cz(ncz)*2

c these are the spontaneous fissioners in sources4c.
      data cz/
     & " h", "he", "li", "be", " b", " c", " n", " o", " f", "ne", "na",
     & "mg", "al", "si", " p", " s", "cl", "ar", " k", "ca", "sc", "ti",
     & " v", "cr", "mn", "fe", "co", "ni", "cu", "zn", "ga", "ge", "as",
     & "se", "br", "kr", "rb", "sr", " y", "zr", "nb", "mo", "tc", "ru",
     & "rh", "pd", "ag", "cd", "in", "sn", "sb", "te", " i", "xe", "cs",
     & "ba", "la", "hf", "ta", " w", "re", "os", "ir", "pt", "au", "hg",
     & "tl", "pb", "bi", "po", "at", "rn", "fr", "ra", "ac", "ce", "pr",
     & "nd", "pm", "sm", "eu", "gd", "tb", "dy", "ho", "er", "tm", "yb",
     & "lu", "th", "pa", " u", "np", "pu", "am", "cm", "bk", "cf", "es",
     & "fm", "md", "no", "lr", "rf", "ha", "sg", "ns", "hs", "mt"/
      iz=0
      do i=1,ncz
        if(cc.eq.cz(i))then
          iz=i*1000
          exit
        end if
      end do ! i
      if(iz.eq.0)then
        write(*,'("error finding element in convz. cc=",a2)')cc
      end if
c
      return
      end subroutine
