      subroutine rdsources_pdata
c read pdata file written by Favorite's version of sources4c
c equation references are to Favorite and Weidenbenner,
c "Sensitivities of a Response to the Composition of an (alpha,n)
c Neutron Source," RPSD-2018.
      use COMS
      implicit none
      integer max_alphas
      parameter (max_alphas=30) ! from sources4c
      real*8 atdens,tmp,tt1,tt2,tt3,term1(maxval(ncb(1:nm))),
     1 term2(maxval(ncb(1:nm))),term3(maxval(ncb(1:nm))),dqdn(nel),
     2 azm,alam,aq,tot_src,tot_tar,an_norm,santmp(neg),
     3 term1m,term2m,term3m,dqdnm,mgspec(neg,2)
      integer i,ic,j,k,k2,k3,kp,nele_stp,niso_tar,niso_src,ng,iz,ios,
     1 jzm,idt,idq,ifs(maxval(ncb(1:nm))),istop(maxval(ncb(1:nm)))
      character ofile*20,line*120
      real*8 tmp1(4,neg),tmp2(4)
c
      do i=1,nm
      if(i.eq.1)then
        k3=0
      else
        k3=k3+ncb(i-1)
      end if
      if(isan(i).eq.0)then
        if(ictrl.eq.3.or.ictrl.eq.4)then
          write(iuo,'(/,"dQ/dN for (alpha,n) sources",/,"  material",i3,
     1     /,"  no (alpha,n) sources in this material.")')i
        end if
        cycle
      end if
c
      atdens=sum(blk(3,k3+1:k3+ncb(i)))
c first derivatives
      term1(1:ncb(i))=0.d0 ! target
      term2(1:ncb(i))=0.d0 ! source
      term3(1:ncb(i))=0.d0 ! stopping
      istop(1:ncb(i))=0
      tmp1(1:4,1:neg)=0.d0
      tmp2(1:4)=0.d0
c
      write(ofile,'("sources/m",i3.3,"_pdata")')i
      open(iun,file=ofile,status='old',iostat=ios)
      if(ios.ne.0)then
        write(*,'("error. no sources4c pdata file ",a)')trim(ofile)
        write(iuo,'("error. no sources4c pdata file ",a)')trim(ofile)
        call stoponerror
      end if
c
c skip ahead to actual number of sources and targets
   10 read(iun,'(a)',iostat=ios)line
      if(ios.ne.0)then
        write(*,'("error in sources4c pdata file ",a)')trim(ofile)
        write(iuo,'("error in sources4c pdata file ",a)')trim(ofile)
        call stoponerror
      end if
      if(line(1:26).ne."number of targets, sources")go to 10
      read(line(27:38),'(2i6)')niso_tar,niso_src
c skip ahead to normalization
      do k=1,niso_tar+niso_src+1
        do j=1,neg+4
          read(iun,'(a)',iostat=ios)line
        end do ! j
      end do ! k
      read(iun,'(a)',iostat=ios)line ! blank
      read(iun,'(28x,e15.7)',iostat=ios)an_norm
      rewind (iun)

c read number of stopping elements and alpha groups
   20 read(iun,'(a)',iostat=ios)line
      if(line(1:18).ne."number of stopping")go to 20
      read(line(26:31),'(i6)')nele_stp
      read(iun,'(a)',iostat=ios)line ! input sources
      read(iun,'(a)',iostat=ios)line ! input targets
      read(iun,'(a)',iostat=ios)line ! note about indices
c
c read combinations of sources and targets
      do ic=1,niso_tar*niso_src
        read(iun,'(a)',iostat=ios)line ! blank
        read(iun,'(20x,i11)',iostat=ios)idt
        read(iun,'(a)',iostat=ios)line ! target at. frac.
        read(iun,'(20x,i11)',iostat=ios)idq
        read(iun,'(12x,e15.7)',iostat=ios)aq ! source atom density
        aq=aq*1.d-24
        read(iun,'(12x,e15.7)',iostat=ios)alam ! source decay const.
        read(iun,'(a)',iostat=ios)line ! blank
c tt1 is sum_l fal*Pi
        read(iun,'(18x,e15.7)',iostat=ios)tt1
c get term1 for target. first term in Eq. (28).
        k2=k3
        do j=1,ncb(i)
          k2=k2+1
c 6013 is written in wrsources for 6000 and 6012.
c 12025 and 12026 are written in wrsources for 12000.
          if(nint(blk(1,k2)).eq.idt.or.(idt.eq.6013.and.
     1     (nint(blk(1,k2)).eq.6000.or.nint(blk(1,k2)).eq.6012)).or.
     2     ((idt.eq.12025.or.idt.eq.12026).and.
     3     nint(blk(1,k2)).eq.12000))go to 50
        end do ! j
        write(*,'("error in pdata file, targets. isotope not ",
     1   "found in material.",/,"material",i4," isotope",i11)')i,idt
        write(iuo,'("error in pdata file, targets. isotope not ",
     1   "found in material.",/,"material",i4," isotope",i11)')i,idt
        call stoponerror
   50   term1(j)=term1(j)+alam*aq/blk(3,k2)*tt1
c get term2 for source. second term in Eq. (28).
        k2=k3
        do j=1,ncb(i)
          k2=k2+1
          if(nint(blk(1,k2)).eq.idq)go to 60
        end do ! j
        write(*,'("error in pdata file, sources. isotope not ",
     1   "found in material.",/,"material",i4," isotope",i11)')i,idq
        write(iuo,'("error in pdata file, sources. isotope not ",
     1   "found in material.",/,"material",i4," isotope",i11)')i,idq
        call stoponerror
   60   term2(j)=term2(j)+alam*tt1
c get term3 for stopping elements (all isotopes). third term in
c Eq. (28).
        read(iun,'(a)',iostat=ios)line ! blank
        read(iun,'(a)',iostat=ios)line ! title
        do k=1,nele_stp
          read(iun,'(4x,i4,2e15.7)')jzm,azm,tt2
          k2=k3
          do j=1,ncb(i)
            k2=k2+1
            if(int(blk(1,k2))/1000.eq.jzm)then
c tt2 is the sum over l in term 3 of Eq. (28), including the factor
c Ni/N, where i is for the target. every isotope in a stopping
c element has the same value for this term.
              term3(j)=term3(j)+alam*aq/atdens*tt2
            end if
          end do ! j
        end do ! k
c terms for second derivative are here.
c TODO
        read(iun,'(a)',iostat=ios)line ! blank
        read(iun,'(a)',iostat=ios)line ! title
        do k=1,nele_stp
          do kp=1,nele_stp
            read(iun,'(4x,i4,e15.7,4x,i4,2e15.7)')jzm,azm,jzm,azm,tt2
          end do ! kp
        end do ! k
      end do ! ic
c
      term1(1:ncb(i))=term1(1:ncb(i))*1.d24
      term2(1:ncb(i))=term2(1:ncb(i))*1.d24
      term3(1:ncb(i))=-term3(1:ncb(i))*1.d24
      dqdn(k3+1:k3+ncb(i))=term1(1:ncb(i))+term2(1:ncb(i))
     1 +term3(1:ncb(i))
c these are only for printing.
      term1m=0.d0
      term2m=0.d0
      term3m=0.d0
      k2=k3
      do j=1,ncb(i)
        k2=k2+1
        term1m=term1m+term1(j)*blk(3,k2)
        term2m=term2m+term2(j)*blk(3,k2)
        term3m=term3m+term3(j)*blk(3,k2)
c a nuclide cannot be a target and an alpha emitter.
c find out if stopping term is largest.
        tmp=abs(term1(j))+abs(term2(j))
        if(tmp.gt.0.d0.and.tmp.lt.abs(term3(j)))then
          istop(j)=1
        end if
      end do ! j
      if(ictrl.eq.3.or.ictrl.eq.4)then
        term1m=term1m/atdens
        term2m=term2m/atdens
        term3m=term3m/atdens
        dqdnm=term1m+term2m+term3m
        write(iuo,'(/,"dQ/dN for (alpha,n) sources",/,"  material",i3,
     1   /,"  isotope",3x,"target",9x,"alpha_src",6x,"stop.elem.",
     2   5x,"total")')i
        k2=k3
        do j=1,ncb(i)
          k2=k2+1
          line=" "
          write(line,'(i10,1p4e15.7)')nint(blk(1,k2)),
     1     term1(j),term2(j),term3(j),dqdn(k2)
          if(istop(k2).eq.1)then
            line=trim(line)//" overall spectrum used"
          end if
          write(iuo,'(a)')trim(line)
        end do ! j
        write(iuo,'(5x,"total",1p4e15.7," sum(dQ/dNj Nj/N)")')term1m,
     1   term2m,term3m,dqdnm
        term1m=term1m*atdens
        term2m=term2m*atdens
        term3m=term3m*atdens
        dqdnm=term1m+term2m+term3m
        write(iuo,'(3x,"N*total",1p4e15.7," sum(dQ/dNj Nj)")')term1m,
     1   term2m,term3m,dqdnm
      end if
c
c read spectra by isotope
  100 read(iun,'(a)',iostat=ios)line
      if(line(1:26).ne."number of targets, sources")go to 100
      tot_src=0.d0
      tot_tar=0.d0
      ifs(1:ncb(i))=0
      do k=1,niso_tar+niso_src
        read(iun,'(a)',iostat=ios)line ! blank
        read(iun,'(20x,i11)',iostat=ios)iz ! not last digit, should be 0
c       write(*,'(i11)')iz
        k2=k3
        do j=1,ncb(i)
          k2=k2+1
c 6013 is written in wrsources for 6000 and 6012.
c 12025 and 12026 are written in wrsources for 12000.
          if(nint(blk(1,k2)).eq.iz.or.(iz.eq.6013.and.
     1     (nint(blk(1,k2)).eq.6000.or.nint(blk(1,k2)).eq.6012)).or.
     2     ((iz.eq.12025.or.iz.eq.12026).and.
     3     nint(blk(1,k2)).eq.12000))go to 110
        end do ! j
        write(*,'("error in pdata file, spectrum block. isotope not ",
     1   "found in material.",/,"material",i4," isotope",i11)')i,iz
        write(iuo,'("error in pdata file, spectrum block. isotope not ",
     1   "found in material.",/,"material",i4," isotope",i11)')i,iz
        call stoponerror
  110   read(iun,'(a)',iostat=ios)line ! heading
        do ng=1,neg
          if(nint(blk(1,k2)).eq.12000.and.iz.eq.12025)then
            read(iun,'(6x,e15.7)')mgspec(ng,1)
          else if(nint(blk(1,k2)).eq.12000.and.iz.eq.12026)then
            read(iun,'(6x,e15.7)')mgspec(ng,2)
          else
            read(iun,'(21x,e15.7)')saniso(ng,k2)
          end if
        end do ! ng
        read(iun,'(6x,e15.7)',iostat=ios)tt3
        if(k.le.niso_tar)then
          tot_tar=tot_tar+tt3
        else
          tot_src=tot_src+tt3
        end if
c use overall spectrum instead if stopping term is large.
        if(istop(k2).eq.1)cycle
c handle elemental magnesium.
        if(nint(blk(1,k2)).eq.12000)then
          if(iz.eq.12025)cycle
          if(iz.eq.12026)then
            tt3=0.d0
            do ng=1,neg
              saniso(ng,k2)=mgspec(ng,1)+mgspec(ng,2)
              tt3=tt3+saniso(ng,k2)
            end do ! ng
            saniso(1:neg,k2)=saniso(1:neg,k2)/tt3
          end if 
        end if
        saniso(1:neg,k2)=saniso(1:neg,k2)*dqdn(k2)
#if ( defined RPSD18 )
c for rpsd 2018 paper--Carbon, Pu-242, and Be-9
        if(nint(blk(1,k2))/1000.eq.6)then
          write(*,'("rdsources_pdata",i10)')nint(blk(1,k2))
          tmp1(2,1:neg)=tmp1(2,1:neg)+saniso(1:neg,k2)*blk(3,k2)
          tmp2(2)=tmp2(2)+blk(3,k2)
          write(*,'("tmp2",1pe16.8)')tmp2(2)
        else if(nint(blk(1,k2)).eq.94242)then
c for Pu add spont. fiss. this logic works because there is only one.
          write(*,'("rdsources_pdata",i10)')nint(blk(1,k2))
          tmp1(3,1:neg)=sfiso(1:neg,k2)+saniso(1:neg,k2)*blk(3,k2)
          tmp2(3)=blk(3,k2)
          write(*,'("tmp2",1pe16.8)')tmp2(3)
        else if(nint(blk(1,k2)).eq.4009)then
          write(*,'("rdsources_pdata",i10)')nint(blk(1,k2))
          tmp1(4,1:neg)=tmp1(4,1:neg)+saniso(1:neg,k2)*blk(3,k2)
          tmp2(4)=blk(3,k2)
          write(*,'("tmp2",1pe16.8)')tmp2(4)
        end if
#endif
        if(ifs(j).eq.1)then
          write(*,'("warning. rdsources_pdata. ifs already 1. ",
     1     "i,k2,iz",2i4,i11)')i,k2,iz
          write(iuo,'("warning. rdsources_pdata. ifs already 1. ",
     1     "i,k2,iz",2i4,i11)')i,k2,iz
        end if
        ifs(j)=1 ! flag to idicate this isotope is done
      end do ! k
c read overall spectrum
      read(iun,'(a)',iostat=ios)line ! blank
      read(iun,'(a)',iostat=ios)line ! title
      read(iun,'(a)',iostat=ios)line ! headings
      do ng=1,neg
        read(iun,'(21x,e15.7)')santmp(ng)
      end do ! ng
      read(iun,'(6x,e15.7)',iostat=ios)tt3
      k2=k3
      do j=1,ncb(i)
        k2=k2+1
        if(ifs(j).eq.0)then
          saniso(1:neg,k2)=santmp(1:neg)*dqdn(k2)
          ifs(j)=1
#if ( defined RPSD18 )
c for rpsd 2018 paper--Mo
          if(nint(blk(1,k2))/1000.eq.42)then
            write(*,'("rdsources_pdata",i10)')nint(blk(1,k2))
            tmp1(1,1:neg)=tmp1(1,1:neg)+saniso(1:neg,k2)*blk(3,k2)
            tmp2(1)=tmp2(1)+blk(3,k2)
            write(*,'("tmp2",1pe16.8)')tmp2(1)
c for rpsd 2018 paper--Carbon, Pu-242, and Be-9
          else if(nint(blk(1,k2))/1000.eq.6)then
            write(*,'("rdsources_pdata",i10)')nint(blk(1,k2))
            tmp1(2,1:neg)=tmp1(2,1:neg)+saniso(1:neg,k2)*blk(3,k2)
            tmp2(2)=tmp2(2)+blk(3,k2)
            write(*,'("tmp2",1pe16.8)')tmp2(2)
          else if(nint(blk(1,k2)).eq.94242)then
c for Pu add spont. fiss. this logic works because there is only one.
            write(*,'("rdsources_pdata",i10)')nint(blk(1,k2))
            tmp1(3,1:neg)=sfiso(1:neg,k2)+saniso(1:neg,k2)*blk(3,k2)
            tmp2(3)=blk(3,k2)
            write(*,'("tmp2",1pe16.8)')tmp2(3)
          else if(nint(blk(1,k2)).eq.4009)then
            write(*,'("rdsources_pdata",i10)')nint(blk(1,k2))
            tmp1(4,1:neg)=tmp1(4,1:neg)+saniso(1:neg,k2)*blk(3,k2)
            tmp2(4)=blk(3,k2)
            write(*,'("tmp2",1pe16.8)')tmp2(4)
          end if
#endif
        end if
      end do ! j
c
      if(ictrl.eq.3.or.ictrl.eq.4)then
        write(iuo,'(/,"(alpha,n) sources",/,"material",i4,/,
     1   "total for targets",1pe15.7,2x,"normed",e15.7,/,
     2   "total for sources",e15.7,2x,"normed",e15.7,/,
     3   "overall",10x,e15.7,2x,"normed",e15.7)')i,tot_tar,
     4   tot_tar/an_norm,tot_src,tot_src/an_norm,tt3,tt3/an_norm
#if ( defined RPSD18 )
c for rpsd 2018 paper
        do j=1,4
          if(tmp2(j).eq.0.d0)cycle
          tmp1(j,1:neg)=tmp1(j,1:neg)/tmp2(j)
        end do ! j
        write(iuo,'(4x,14x,"Mo",15x,"C",11x,"Pu242",14x,"Be")')
        do ng=1,neg
          write(iuo,'(i4,1p4e16.8)')ng,(tmp1(j,ng),j=1,4)
        end do ! ng
#endif
c
c debug print
c       write(iuo,'(/,"dQ_g/dN",/,5x,"g",13x,"E",9999i14)')
c    1   (j,j=k3+1,k3+ncb(i))
c       do ng=1,neg
c         write(iuo,'(i6,1p9999e14.6)')ng,ebins(3,ng),
c    1     saniso(ng,k3+1:k3+ncb(i))
c       end do ! ng
c       write(iuo,'(6x,14x,1p9999e14.6)')
c    1   (sum(saniso(1:neg,j)),j=k3+1,k3+ncb(i))
      end if ! ictrl
c
  200 close(iun)
c debug print
c     k2=k3
c     do j=1,ncb(i)
c       k2=k2+1
c       write(*,'(i11)')nint(blk(1,k2))
c       do ng=1,neg
c         write(*,'(i6,1pe15.7)')ng,saniso(ng,k2)
c       end do ! ng
c     end do ! j
      end do ! i
      return
      end subroutine
