      subroutine rdmisc
c read misc output files for neutron source.
c the source is in n/s-cm^3.
c if misc is modified to handle repeated isotopes (lnk3dnt),
c then this routine needs to be updated.
      use COMS
      implicit none
      real*8 nsrcfsf(neg,nm),siniso(nel),sf(nel),sftot,st
      integer i,ic,ios,ng,itmp,j,k,k2,nsf,isf(nel)
      character ofile*22,line*120
c
c (alpha,n) not yet included in misc.
c     nsrcfan(1:neg,1:nm)=0.d0
      nsrcfsf(1:neg,1:nm)=0.d0
      siniso(1:nel)=0.d0
      chisrc(1:neg,0:nm)=0.d0
      if(ictrl.eq.2)then
        write(iuo,'(/,"spontaneous fission source from misc, ",
     1   "(alpha,n) source ignored in misc")')
        write(iuo,'("output from misc")')
      end if
      k2=0
      do i=1,nm
        k2=k2+ncb(i-1)
c       sant(1:nel)=0.d0
c       sans(1:nel)=0.d0
        sf(1:nel)=0.d0
c spectrum file
        write(ofile,'("misc/m",i6.6,"_misc.src")')i
        open(iun,file=ofile,status='old',iostat=ios)
        if(ios.ne.0)then
          write(*,'("error. no misc .src file ",a)')trim(ofile)
          write(iuo,'("error. no misc .src file ",a)')trim(ofile)
          call stoponerror
        end if
        read(iun,'(a)',iostat=ios)line
c if no sources, "no sources" written by control script
        if(line(1:10).eq."no sources")then
          if(ictrl.eq.2)then
            write(*,'(/,2x,"no misc neutron source in material",i7)')i
            write(iuo,'(/,2x,"no misc neutron source in material",i7)')i
          end if
          close(iun)
          cycle
        end if
c a source material.
        ismat(i)=1
        read(line(19:30),'(e12.5)')nsrc(i) ! will be overwritten
        read(iun,'(a)')line ! si, sp
        read(iun,'(a)')line ! h, d
        read(iun,'(a)')line ! extra bin
        read(iun,'(a)')line ! extra bin
        sftot=0.d0
        do ng=neg,1,-1
          read(iun,'(17x,e12.5)')nsrcfsf(ng,i)
          if(nsrcfsf(ng,i).lt.0.d0)nsrcfsf(ng,i)=0.d0
          sftot=sftot+nsrcfsf(ng,i)
        end do ! ng
        close(iun)
c nsrcf is used in wrdantnm; nsrc is the group total (read above).
c here these are normalized spectra, not absolute.
c convert spectra to n/s/cm^3 later on
c renormalize for consistency, to overcome roundoff errors from reading
        nsrcfsf(1:neg,i)=nsrcfsf(1:neg,i)/sftot
        nsrcf(1:neg,i)=nsrcfsf(1:neg,i)
c source fission spectrum for feyny and sm2.
c this uses spont. fiss. only; how to account for (alpha,n)?
        chisrc(1:neg,i)=nsrcfsf(1:neg,i)
        if(ictrl.eq.2.and.abs(sftot-1.d0).gt.1.d-5)then
          write(*,'("misc s.f. source spectrum not normalized for",
     1     " material",i7,"; sum=",0pf12.9)')i,sftot
          write(iuo,'("misc s.f. source spectrum not normalized for",
     1     " material",i7,"; sum=",0pf12.9)')i,sftot
        end if
c
c magnitude file by isotope
        write(ofile,'("misc/m",i6.6,"_misc.out")')i
        open(iun,file=ofile,status='old',iostat=ios)
        if(ios.ne.0)then
          write(*,'("error. no misc .out file ",a)')trim(ofile)
          write(iuo,'("error. no misc .out file ",a)')trim(ofile)
          call stoponerror
        end if
   30   read(iun,'(a)',iostat=ios)line
        if(line(1:30).ne."ISOMER PARTICLE EMISSION RATES")go to 30
        nsf=0
        ios=0
        do while(ios.eq.0)
          read(iun,'(a)',iostat=ios)line ! blank
          if(ios.ne.0)exit
          read(iun,'(a)')line
          ic=index(line,",")
          read(line(1:ic-1),*)itmp
          read(iun,'(a)')line ! blank
          read(iun,'(a)')line ! Particle
   50     read(iun,'(a)')line
          if(line(11:11).eq."n")then
            nsf=nsf+1
            isf(nsf)=itmp
            read(line,'(12x,e12.5)')sf(nsf)
          end if
          if(line(7:11).ne."Total")go to 50
        end do
        do k=1,nsf
          do j=1,ncb(i)
            if(isf(k).eq.nint(blk(1,k2+j)))then
c apply the overall s.f. spectrum to all isotopes until misc writes
c that info.
              sfiso(1:neg,k2+j)=sf(k)*nsrcfsf(1:neg,i)
              siniso(k2+j)=siniso(k2+j)+sf(k)
            end if
          end do ! j
        end do ! k
c convert spectrum to n/s/cm^3 and overwrite the total that was read
c (use sum of isotopes instead for conistency).
        nsrc(i)=sum(sf(1:nsf))
        nsrcfsf(1:neg,i)=nsrcfsf(1:neg,i)*nsrc(i)
        nsrcf(1:neg,i)=nsrcf(1:neg,i)*nsrc(i)
c
c write s.f. sources.
        if(ictrl.eq.2)then
          write(iuo,'(/,"  material",i7,/,
     1     "  isotope   (alpha,n)    spont.fiss.  total")')i
          do j=1,ncb(i)
            write(iuo,'(i10,1p3e13.5)')nint(blk(1,k2+j)),0.d0,
     1       sum(sfiso(1:neg,k2+j)),siniso(k2+j)
          end do ! j
          write(iuo,'(2x,"totals [(alpha,n) targets and sources ",
     1     "should be equal]",/,2x,"(alpha,n) targets",1pe17.5,/,
     2     2x,"(alpha,n) sources",e17.5,/,
     3     2x,"spont.fiss. sources",e15.5,/,2x,"total",e29.5)')
     4     0.d0,0.d0,sum(sf(1:nsf)),0.d0+sum(sf(1:nsf))
          st=(sum(nsrcfsf(1:neg,i))-sum(sf(1:nsf)))/sum(sf(1:nsf))
          write(iuo,'(2x,"compare misc spont. fiss. overall ",
     1     "total and sum of isotopes.",/,2x,
     2     "overall total=",1pe16.8,/,2x,
     3     "sum of isos=  ",1pe16.8,/,2x,"rel. diff.=   ",1pe16.8)')
     4     sum(nsrcfsf(1:neg,i)),sum(sf(1:nsf)),st
          write(iuo,'("  the (alpha,n) sources are ignored in misc ",
     1     "calculation.")')
        end if
c
      end do ! i
c
      if(sum(ismat(1:nm)).eq.0)then
        write(*,'("error. no neutron sources in rdmisc.")')
        write(iuo,'("error. no neutron sources in rdmisc.")')
        call stoponerror
      end if
c
      return
      end
