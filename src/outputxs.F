      subroutine write_xsecs
c output cross sections.
      use COMS
      implicit none
      integer i,j,iso,ir,jr,k2,k3,m,ng,nx
      character sfile*20
      sfile="xsecs"
      write(*,'(/,"writing cross sections to file ",a,".")')trim(sfile)
      open(iux2,file=sfile,status='unknown')
      write(iux2,'(a)')trim(id)
      write(iux2,'("material macroscopic and isotope microscopic ",
     1 "cross sections")')
      call block3(iux2,0,3)
c materials
      if(ilnk3dnt.eq.0)then
        k2=0
        do i=1,nm
c compute edit cross sections for material. never needed.
c cross sections are zeroed in rdxsmg0.
          do j=1,ncb(i)
            k2=k2+1
            iso=nm+k2
            sigf(1:neg,i)=sigf(1:neg,i)+sigf(1:neg,iso)*blk(3,k2)
            sigc(1:neg,i)=sigc(1:neg,i)+sigc(1:neg,iso)*blk(3,k2)
            sige(1:neg,i)=sige(1:neg,i)+sige(1:neg,iso)*blk(3,k2)
            sigi(1:neg,i)=sigi(1:neg,i)+sigi(1:neg,iso)*blk(3,k2)
            if(itrcor.eq.1.or.itrcor.eq.2)then
              sigsp1(1:neg,1:neg,i)=sigsp1(1:neg,1:neg,i)
     1         +sigsp1(1:neg,1:neg,iso)*blk(3,k2)
            end if
          end do ! j
          write(iux2,'(/,"material",i7,14x,"atom density",1pe19.11)')i,
     1     rhoa(i)
          call outputxs(iux2,neg,isct,ebins,sigt(1:neg,i),
     1     siga(1:neg,i),sigc(1:neg,i),chi(1:neg,1:neg,i),
     2     nusigf(1:neg,i),sigf(1:neg,i),
     3     sigs(1:neg,1:neg,0:isct,i),sige(1:neg,i),sigi(1:neg,i))
c write isct+1 scattering if transport correction (diag or bhs)
          if(itrcor.eq.1.or.itrcor.eq.2)then
            call outputsigsp1(iux2,neg,isct,sigsp1(1:neg,1:neg,i))
          end if
        end do ! i
      else if(ilnk3dnt.gt.0)then
        i=0
        k2=0
        do ir=1,nr
          do jr=1,nz
            if(idclnk3(ir,jr,1).eq.0)cycle
            k2=k2+1
            call maclnk3(i,ir,jr)
c compute edit cross sections for material. never needed.
c cross sections are zeroed in rdxsmg0.
            do nx=1,nmxsp
              iso=idclnk3(ir,jr,nx)
              if(iso.eq.0)exit
              iso=nm+iso
              sigf(1:neg,i)=sigf(1:neg,i)+sigf(1:neg,iso)
     1         *denlnk3(ir,jr,nx)
              sigc(1:neg,i)=sigc(1:neg,i)+sigc(1:neg,iso)
     1         *denlnk3(ir,jr,nx)
              sige(1:neg,i)=sige(1:neg,i)+sige(1:neg,iso)
     1         *denlnk3(ir,jr,nx)
              sigi(1:neg,i)=sigi(1:neg,i)+sigi(1:neg,iso)
     1         *denlnk3(ir,jr,nx)
            end do ! nx
            write(iux2,'(/,"material",i7)')k2
            call outputxs(iux2,neg,isct,ebins,sigt(1:neg,i),
     1       siga(1:neg,i),sigc(1:neg,i),chi(1:neg,1:neg,i),
     2       nusigf(1:neg,i),sigf(1:neg,i),
     3       sigs(1:neg,1:neg,0:isct,i),sige(1:neg,i),sigi(1:neg,i))
c write isct+1 scattering if transport correction (diag or bhs)
            if(itrcor.eq.1.or.itrcor.eq.2)then
              call outputsigsp1(iux2,neg,isct,sigsp1(1:neg,1:neg,i))
            end if
          end do ! jr
        end do ! ir
      end if
c isotopes
      k2=0
      m=1
      do i=1,nel
        k2=k2+1
        k3=nm+i
        write(iux2,'(/,1x,a24," isotope",i4," in material",i7,
     1   14x,"atom density",1pe19.11)')zaidfull(i),k2,m,blk(3,i)
        call outputxs(iux2,neg,isct,ebins,sigt(1:neg,k3),
     1   siga(1:neg,k3),sigc(1:neg,k3),chi(1:neg,1:neg,k3),
     2   nusigf(1:neg,k3),sigf(1:neg,k3),
     3   sigs(1:neg,1:neg,0:isct,k3),sige(1:neg,k3),sigi(1:neg,k3))
c write isct+1 scattering if transport correction (diag or bhs)
        if(itrcor.eq.1.or.itrcor.eq.2)then
          call outputsigsp1(iux2,neg,isct,sigsp1(1:neg,1:neg,k3))
        end if
c       if(sum(chi(1:neg,1:neg,k3)).ne.0.d0)then
          write(iux2,'("spectrum weight function")')
          do ng=1,neg
            write(iux2,'(i6,1pe19.11)')ng,wgtfct(ng,i)
          end do ! ng
c       end if
        if(k2.eq.ncb(m))then
          k2=0
          m=m+1
        end if
      end do ! i
c velocities
      write(iux2,'(1x,"group",2x,"e_lower",12x,"e_upper",12x,
     1 "e_average",10x,"velocity")')
      do ng=1,neg
        write(iux2,'(i6,1p9999e19.11)')ng,ebins(1:3,ng),vel(ng)
      end do ! ng
      close(iux2)
      return
      end subroutine
c
      subroutine outputxs(iux2,neg,isct,ebins,sigt,siga,sigc,chi,
     1 nusigf,sigf,sigs,sige,sigi)
c output cross sections.
      implicit none
      integer iux2,neg,isct
      real*8 ebins(3,neg),sigt(neg),siga(neg),sigc(neg),
     1 nusigf(neg),sigf(neg),chi(neg,neg),sigs(neg,neg,0:isct),
     2 sige(neg),sigi(neg)
      real*8 scti(neg,0:isct),scto(neg,0:isct),sctt(neg)
      integer l,lg,ng
c there are several instances of 9999e19.11 and 9999i19. the number
c of energy groups neg must be .le. 9999 and the scattering order
c isct must satisfy 13+3*isct.le.9999 or isct.lt.3300.
c
c calculate total scattering, in-scattering, and out-scattering.
      scti(1:neg,0:isct)=0.d0
      scto(1:neg,0:isct)=0.d0
      sctt(1:neg)=0.d0
      do ng=1,neg
        do lg=1,neg
          do l=0,isct
            sctt(ng)=sctt(ng)+sigs(ng,lg,l)
            if(lg.ne.ng)then
              scti(ng,l)=scti(ng,l)+sigs(ng,lg,l)
              scto(lg,l)=scto(lg,l)+sigs(ng,lg,l)
            end if
          end do ! l
        end do ! lg
      end do ! ng
c no summary, just the details
      write(iux2,'(1x,"group",2x,"e_lower",12x,"e_upper",12x,
     1 "e_average",10x,"total",14x,"abs",16x,"(n,g)",14x,
     2 "chi",16x,"nu*fiss",12x,"fiss",15x,"elastic",12x,"inelastic",
     3 10x)',advance='no')
      do l=0,isct
        if(l.le.9)then
          write(iux2,'("in-scat-",i1,10x,"self-scat-",i1,8x,
     1     "out-scat-",i1,9x)',advance='no')l,l,l
        else if(l.le.99)then
          write(iux2,'("in-scat-",i2,9x,"self-scat-",i2,7x,
     1     "out-scat-",i2,8x)',advance='no')l,l,l
        else
          write(iux2,'("in-scat-",i3,8x,"self-scat-",i3,6x,
     1     "out-scat-",i3,7x)',advance='no')l,l,l
        end if
      end do ! l
      write(iux2,'("sctt")',advance='yes')
      do ng=1,neg
        write(iux2,'(i6,1p9999e19.11)')ng,ebins(1:3,ng),
     1   sigt(ng),siga(ng),sigc(ng),sum(chi(ng,1:neg))/real(neg,8),
     2   nusigf(ng),sigf(ng),sige(ng),sigi(ng),
     4   (scti(ng,l),sigs(ng,ng,l),scto(ng,l),
     5   l=0,isct),sctt(ng)
      end do ! ng
      if(sum(chi(1:neg,1:neg)).ne.0.d0)then
        write(iux2,'("fission chi matrix",
     1   " (columns are from, rows are to)")')
        write(iux2,'(9999i19)')(ng,ng=1,neg)
        do ng=1,neg
          write(iux2,'(i6,1p9999e19.11)')ng,(chi(ng,lg),lg=1,neg)
        end do ! ng
      end if
      do l=0,isct
        write(iux2,'("scattering matrix, order=",i2,
     1   " (columns are outscattering, rows are inscattering)")')l
        write(iux2,'(9999i19)')(ng,ng=1,neg)
        do ng=1,neg
          write(iux2,'(i6,1p9999e19.11)')ng,(sigs(ng,lg,l),lg=1,neg)
        end do ! ng
      end do ! l
c
      return
      end subroutine
c
      subroutine outputsigsp1(iux2,neg,isct,sigsp1)
c output sigsp1 (for transport correction).
      implicit none
      integer iux2,neg,isct
      real*8 sigsp1(neg,neg)
      integer l,lg,ng
c
      l=isct+1
      write(iux2,'("scattering matrix, order=",i2,
     1 " (columns are outscattering, rows are inscattering)")')l
      write(iux2,'(9999i19)')(ng,ng=1,neg)
      do ng=1,neg
        write(iux2,'(i6,1p9999e19.11)')ng,(sigsp1(ng,lg),lg=1,neg)
      end do ! ng
      return
      end subroutine
