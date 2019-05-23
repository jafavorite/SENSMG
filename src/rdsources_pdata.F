      subroutine rdsources_pdata
c read pdata file written by Favorite's version of sources4c
c equation references (1st) are to Favorite and Weidenbenner,
c "Sensitivities of a Response to the Composition of an (alpha,n)
c Neutron Source," RPSD-2018; (2nd) are to Favorite's memo,
c "Second Derivative of an (alpha,n) Neutron Source with Respect
c to Constituent Isotope Densities," LA-UR-18-29134, Sept. 19, 2018.
c this routine requires sources4c.jaf6.for.
      use COMS
      implicit none
      integer max_alphas
      parameter (max_alphas=30) ! from sources4c
      real*8 tmp,tt1,tt2,tt3,term1(maxval(ncb(1:nm))),
     1 term2(maxval(ncb(1:nm))),term3(maxval(ncb(1:nm))),
     2 at,alam,aq,tot_an,an_norm,cnst2,cnst3,
     3 term1m,term2m,term3m,dqdnm,
     4 term12(maxval(ncb(1:nm)),maxval(ncb(1:nm))),
     5 term22(maxval(ncb(1:nm)),maxval(ncb(1:nm))),
     6 term32(maxval(ncb(1:nm)),maxval(ncb(1:nm))),
     7 term42(maxval(ncb(1:nm)),maxval(ncb(1:nm))),
     8 dqdn2,dqdn2m,at_to_use,espec(neg)
      real*8 saniso2(neg,nel,nel) ! will go to coms
      integer i,i1,i1i,i2i,i2,j,j2,jp,k,k2,k2p,k3,kp,ng,ntmp,
     1 nele_stp,niso_tar,niso_src,iz,ios,jzm,jzm2,idt,idq,
     2 itarget(maxval(ncb(1:nm))),isource(maxval(ncb(1:nm))),
     3 itarindex(maxval(ncb(1:nm))),isrcindex(maxval(ncb(1:nm)))
      character ofile*22,line*120,tape8*22
      real*8 tmp1(4,neg),tmp2(4),tmp3,tmp4(0:4,neg)
      real*8 tmp5(100,neg)
c
      do i=1,nm
      if(i.eq.1)then
        k3=0
      else
        k3=k3+ncb(i-1)
      end if
      if(isan(i).eq.0)then
        if(ictrl.eq.3.or.ictrl.eq.4)then
          write(iuo,'(/,"dQ/dN for (alpha,n) sources",/,"  material",i7,
     1     /,"  no (alpha,n) sources in this material.")')i
        end if
        cycle
      end if
c
      tot_an=0.d0
c first derivatives
      term1(1:ncb(i))=0.d0 ! target
      term2(1:ncb(i))=0.d0 ! source
      term3(1:ncb(i))=0.d0 ! stopping
      tmp1(1:4,1:neg)=0.d0
      tmp2(1:4)=0.d0
      tmp3=0.d0
      tmp4(0:4,1:neg)=0.d0
      saniso(1:neg,1:nel)=0.d0
c second derivatives
      term12(1:ncb(i),1:ncb(i))=0.d0 !
      term22(1:ncb(i),1:ncb(i))=0.d0 !
      term32(1:ncb(i),1:ncb(i))=0.d0 !
      term42(1:ncb(i),1:ncb(i))=0.d0 !
      itarget(1:ncb(i))=0
      isource(1:ncb(i))=0
      itarindex(1:ncb(i))=0
      isrcindex(1:ncb(i))=0
      saniso2(1:neg,1:nel,1:nel)=0.d0
c
      write(ofile,'("sources/m",i6.6,"_pdata")')i
      open(iun,file=ofile,status='old',iostat=ios)
      if(ios.ne.0)then
        write(*,'("error. no sources4c pdata file ",a)')trim(ofile)
        write(iuo,'("error. no sources4c pdata file ",a)')trim(ofile)
        call stoponerror
      end if
      write(tape8,'("sources/m",i6.6,"_tape8")')i
      open(iun8,file=tape8,status='old',iostat=ios)
      if(ios.ne.0)then
        write(*,'("error. no sources4c tape8 file ",a)')trim(tape8)
        write(iuo,'("error. no sources4c tape8 file ",a)')trim(tape8)
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
c read targets and sources
      do k=1,niso_tar
        read(iun,'(20x,i11)',iostat=ios)itarget(k)
      end do ! k
      do k=1,niso_src
        read(iun,'(20x,i11)',iostat=ios)isource(k)
      end do ! k
      rewind (iun)
c
c associate targets and sources with isotopes in ncb
c 6013 is written in wrsources for 6000 and 6012.
c 12025 and 12026 are written in wrsources for 12000.
c TODO: handle all similar cases.
      do j=1,ncb(i)
        iz=nint(blk(1,j))
        do j2=1,niso_tar
          if(iz.eq.itarget(j2).or.(itarget(j2).eq.6013.and.
     1     (iz.eq.6000.or.iz.eq.6012)).or.
     2     ((itarget(j2).eq.12025.or.itarget(j2).eq.12026).and.
     3     iz.eq.12000))then
            itarindex(j2)=j
          end if
        end do ! j2
        do j2=1,niso_src
          if(iz.eq.isource(j2))then
            isrcindex(j2)=j
          end if
        end do ! j2
      end do ! j
c     write(iuo,'(/,"targets")')
      do j2=1,niso_tar
c       write(iuo,'(i8,i4)')itarget(j2),itarindex(j2)
        if(itarindex(j2).eq.0)then
          write(*,'("error in pdata file. target isotope not ",
     1     "found in material.",/,"material",i7," isotope",i11)')i,
     2     itarget(j2)
          write(iuo,'("error in pdata file. target isotope not ",
     1     "found in material.",/,"material",i7," isotope",i11)')i,
     2     itarget(j2)
          call stoponerror
        end if
      end do ! j2
c     write(iuo,'(/,"sources")')
      do j2=1,niso_src
c       write(iuo,'(i8,i4)')isource(j2),isrcindex(j2)
        if(isrcindex(j2).eq.0)then
          write(*,'("error in pdata file. source isotope not ",
     1     "found in material.",/,"material",i7," isotope",i11)')i,
     2     isource(j2)
          write(iuo,'("error in pdata file. source isotope not ",
     1     "found in material.",/,"material",i7," isotope",i11)')i,
     2     isource(j2)
          call stoponerror
        end if
      end do ! j2
c
c read number of stopping elements
   20 read(iun,'(a)',iostat=ios)line
      if(line(1:18).ne."number of stopping")go to 20
      read(line(26:31),'(i6)')nele_stp
      read(iun,'(a)',iostat=ios)line ! input sources
      read(iun,'(a)',iostat=ios)line ! input targets
      read(iun,'(a)',iostat=ios)line ! note about indices
c
c advance tape8
      do ng=1,neg/8+13
        read(iun8,'(a)',iostat=ios)line
      end do ! ng
c
c read combinations of sources and targets
      ntmp=0
      do i1=1,niso_tar
        i1i=itarindex(i1)
      do i2=1,niso_src
        i2i=isrcindex(i2)
        ntmp=ntmp+1
        read(iun,'(a)',iostat=ios)line ! blank
        read(iun,'(20x,i11)',iostat=ios)idt ! target zaid
        read(iun,'(12x,e15.7)',iostat=ios)at ! target at. frac.
        at=at*rhoa(i)  ! at is now target atom density
        read(iun,'(20x,i11)',iostat=ios)idq ! source zaid
        read(iun,'(12x,e15.7)',iostat=ios)aq ! source atom density
        aq=aq*1.d-24 ! convert to atom/b-cm
        read(iun,'(12x,e15.7)',iostat=ios)alam ! source decay const.
        read(iun,'(a)',iostat=ios)line ! blank
c tt1 is the (alpha,n) source rate density for this source and target.
        read(iun,'(59x,e15.7)',iostat=ios)tt1
        tot_an=tot_an+tt1
c       write(iuo,'("*** tar,src ***",4i8)')i1i,idt,i2i,idq
c
c read spectrum for this source and target
   30   read(iun8,'(a)',iostat=ios)line
        if(line(2:6).ne."(a,n)")go to 30
        read(iun8,*)espec(1:neg)
        read(iun8,'(55x,e15.7)')an_norm
        espec(1:neg)=espec(1:neg)/an_norm
        tmp5(ntmp,1:neg)=espec(1:neg)
c
c get term1 for target. first term in Eq. (29) (1st).
c use blk(3,k2) for 6013 because we want the derivative with respect to
c carbon atom density, not C-13 atom density.
c same argument for components of magnesium (12000).
        if(idt.eq.6013.or.idt/1000.eq.12)then
          at_to_use=blk(3,k3+i1i)
        else
          at_to_use=at
        end if
        cnst2=tt1/at_to_use
        term1(i1i)=term1(i1i)+cnst2
        saniso(1:neg,k3+i1i)=saniso(1:neg,k3+i1i)
     1   +espec(1:neg)*cnst2
c get term2 for source. second term in Eq. (29) (1st).
        cnst2=tt1/aq
        term2(i2i)=term2(i2i)+cnst2
        saniso(1:neg,k3+i2i)=saniso(1:neg,k3+i2i)
     1   +espec(1:neg)*cnst2
c get term3 for stopping elements (all isotopes). third term in
c Eq. (29) (1st).
        read(iun,'(a)',iostat=ios)line ! blank
        read(iun,'(a)',iostat=ios)line ! title
        do k=1,nele_stp
          read(iun,'(4x,i4,2e15.7)')jzm,tmp,tt2
          tt2=tt2*alam/rhoa(i)*1.d24
          k2=k3
          do j=1,ncb(i)
            k2=k2+1
            iz=nint(blk(1,k2))
            if(iz/1000.eq.jzm)then
c j and k2 identify an isotope corresponding to this stopping element.
c tt2 is the sum over l in term 3 of Eq. (29) (1st), including the
c factor Ni/N from Eq. (22) (1st), where i is for the target. every 
c isotope in a stopping element has the same value for this term.
              cnst2=-aq*tt2
              term3(j)=term3(j)+cnst2
              saniso(1:neg,k2)=saniso(1:neg,k2)
     1         +espec(1:neg)*cnst2
c second derivatives, second and third terms in Eq. (24) (2nd).
c i1i is the present target, i2i the present source
              if(i_2nd_order.eq.1)then
c               write(iuo,'("2nd and 3rd",5i8,1p2e14.6)')jzm,iz,j,
c    1           i1i,i2i,aq/at_to_use*tt2,tt2
                cnst2=-aq/at_to_use*tt2
                term22(i1i,j)=term22(i1i,j)+cnst2
                term22(j,i1i)=term22(j,i1i)+cnst2
                saniso2(1:neg,k3+i1i,k2)=saniso2(1:neg,k3+i1i,k2)
     1           +espec(1:neg)*cnst2
                saniso2(1:neg,k2,k3+i1i)=saniso2(1:neg,k2,k3+i1i)
     1           +espec(1:neg)*cnst2
                cnst2=-tt2
                term32(i2i,j)=term32(i2i,j)+cnst2
                term32(j,i2i)=term32(j,i2i)+cnst2
                saniso2(1:neg,k3+i2i,k2)=saniso2(1:neg,k3+i2i,k2)
     1           +espec(1:neg)*cnst2
                saniso2(1:neg,k2,k3+i2i)=saniso2(1:neg,k2,k3+i2i)
     1           +espec(1:neg)*cnst2
              end if
            end if
          end do ! j
        end do ! k
c read terms for second derivative. fourth term in Eq. (24) (2nd).
        read(iun,'(a)',iostat=ios)line ! blank
        read(iun,'(a)',iostat=ios)line ! title
        cnst2=2.d24*aq*alam/rhoa(i)**2
        do k=1,nele_stp
          do kp=1,nele_stp
            read(iun,'(4x,i4,e15.7,4x,i4,2e15.7)')jzm,tmp,jzm2,tmp,tt3
            if(i_2nd_order.eq.1)then
              k2=k3
              do j=1,ncb(i)
                k2=k2+1
                k2p=k3
                do jp=1,ncb(i)
                  k2p=k2p+1
                  if(nint(blk(1,k2))/1000.eq.jzm.and.
     1             nint(blk(1,k2p))/1000.eq.jzm2)then
                    cnst3=cnst2*tt3
                    term42(j,jp)=term42(j,jp)+cnst3
                    saniso2(1:neg,k2,k2p)=saniso2(1:neg,k2,k2p)
     1               +espec(1:neg)*cnst3
                  end if
                end do ! jp
              end do ! j
            end if
          end do ! kp
        end do ! k
c
c second derivatives, first term in Eq. (24) (2nd).
        if(i_2nd_order.eq.1)then
          cnst2=tt1/at_to_use/aq
          term12(i1i,i2i)=term12(i1i,i2i)+cnst2
          term12(i2i,i1i)=term12(i2i,i1i)+cnst2
          saniso2(1:neg,k3+i1i,k3+i2i)=saniso2(1:neg,k3+i1i,k3+i2i)
     1      +espec(1:neg)*cnst2
          saniso2(1:neg,k3+i2i,k3+i1i)=saniso2(1:neg,k3+i2i,k3+i1i)
     1      +espec(1:neg)*cnst2
        end if
      end do ! i2
      end do ! i1
c
c first derivatives
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
      end do ! j
c
c print
      if(ictrl.eq.3.or.ictrl.eq.4)then
c first derivatives
        term1m=term1m/rhoa(i)
        term2m=term2m/rhoa(i)
        term3m=term3m/rhoa(i)
        dqdnm=term1m+term2m+term3m
        write(iuo,'(/,"dQ/dN for (alpha,n) sources",/,"  material",i7,
     1   /,"  isotope",3x,"target",9x,"alpha_src",6x,"stop.elem.",
     2   5x,"total")')i
        k2=k3
        do j=1,ncb(i)
          k2=k2+1
          write(iuo,'(i10,1p4e15.7)')nint(blk(1,k2)),
     1     term1(j),term2(j),term3(j),term1(j)+term2(j)+term3(j)
        end do ! j
        write(iuo,'(5x,"total",1p4e15.7," sum(dQ/dNj Nj/N)")')term1m,
     1   term2m,term3m,dqdnm
        term1m=term1m*rhoa(i)
        term2m=term2m*rhoa(i)
        term3m=term3m*rhoa(i)
        dqdnm=term1m+term2m+term3m
        write(iuo,'(3x,"N*total",1p4e15.7," sum(dQ/dNj Nj)")')term1m,
     1   term2m,term3m,dqdnm
        tmp=abs(term1m/tot_an-1.d0)
        if(tmp.gt.1.d-6)then
          write(*,'("warning. dQ/dN for targets not equal to Q; ",
     1     "diff=",2pf8.4,"%")')tmp
          write(iuo,'("warning. dQ/dN for targets not equal to Q; ",
     1     "diff=",2pf8.4,"%")')tmp
        end if
        tmp=abs(term2m/tot_an-1.d0)
        if(tmp.gt.1.d-6)then
          write(*,'("warning. dQ/dN for sources not equal to Q; ",
     1     "diff=",2pf8.4,"%")')tmp
          write(iuo,'("warning. dQ/dN for sources not equal to Q; ",
     1     "diff=",2pf8.4,"%")')tmp
        end if
        tmp=abs(term3m/tot_an+1.d0)
        if(tmp.gt.1.d-6)then
          write(*,'("warning. dQ/dN for stopping not equal to -Q; ",
     1     "diff=",2pf8.4,"%")')tmp
          write(iuo,'("warning. dQ/dN for stopping not equal to -Q; ",
     1     "diff=",2pf8.4,"%")')tmp
        end if
        tmp=abs(dqdnm/tot_an-1.d0)
        if(tmp.gt.1.d-6)then
          write(*,'("warning. total N*dQ/dN not equal to Q; ",
     1     "diff=",2pf8.4,"%")')tmp
          write(iuo,'("warning. total N*dQ/dN not equal to Q; ",
     1     "diff=",2pf8.4,"%")')tmp
        end if
c second derivatives
        if(i_2nd_order.eq.1)then
          dqdnm=dqdnm/rhoa(i) ! total 1st deriv wrt material atom dens
          dqdn2m=0.d0
          tmp4(0,1:neg)=0.d0
          write(iuo,'(/,"d^2Q/dNjdNj` for (alpha,n) sources")')
          k2=k3
          do j=1,ncb(i)
            k2=k2+1
            write(iuo,'("  material",i7,"  isotope",i10)')i,
     1       nint(blk(1,k2))
            k2p=k3
            do jp=1,ncb(i)
              k2p=k2p+1
              dqdn2=term12(k2,k2p)+term22(k2,k2p)+term32(k2,k2p)
     1         +term42(k2,k2p)
              write(iuo,'(i10,1p5e15.7)')nint(blk(1,k2p)),
     1         term12(k2,k2p),term22(k2,k2p),term32(k2,k2p),
     2         term42(k2,k2p),dqdn2
              dqdn2m=dqdn2m+dqdn2*blk(3,k2)*blk(3,k2p)
              tmp4(0,1:neg)=tmp4(0,1:neg)+saniso2(1:neg,k2,k2p)
     1         *blk(3,k2)*blk(3,k2p)
            end do ! jp
          end do ! j
          dqdn2m=dqdn2m/rhoa(i)**2
          tmp4(0,1:neg)=tmp4(0,1:neg)/rhoa(i)**2
          write(iuo,'("  material",i7,2x,"dQ/dN",10x,"d^2Q/dN^2",/,
     1     10x,1p2e15.7)')i,dqdnm,dqdn2m
          if(abs(dqdn2m).gt.1.d-6)then
            write(*,'("warning. total d^2Q/dN^2 not equal to 0; ",
     1       1pe12.4)')dqdn2m
            write(iuo,'("warning. total d^2Q/dN^2 not equal to 0; ",
     1       1pe12.4)')dqdn2m
          end if
          tmp=abs(dqdn2m/dqdnm)
          if(tmp.gt.1.d-6)then
            write(*,'("warning. total d^2Q/dN^2 not small compared ",
     1       "to dQ/dN; ",1pe12.4)')tmp
            write(iuo,'("warning. total d^2Q/dN^2 not small compared ",
     1       "to dQ/dN; ",1pe12.4)')tmp
          end if
        end if ! i_2nd_order
      end if
c
      close(iun)
      close(iun8)
c
c     if(ictrl.eq.3.or.ictrl.eq.4)then
c debug print
c       write(iuo,'(/,"dQ_g/dN",/,5x,"g",13x,"E",9999i14)')
c    1   (j,j=k3+1,k3+ncb(i))
c       do ng=1,neg
c         write(iuo,'(i6,1p9999e14.6)')ng,ebins(3,ng),
c    1     saniso(ng,k3+1:k3+ncb(i))
c       end do ! ng
c       write(iuo,'(6x,14x,1p9999e14.6)')
c    1   (sum(saniso(1:neg,j)),j=k3+1,k3+ncb(i))
c     end if ! ictrl
#if ( defined RPSD18 )
c debug print
      k2=k3
      do j=1,ncb(i)
        k2=k2+1
c for rpsd 2018 paper--Carbon, Pu-242, Be-9, and Mo
        if(nint(blk(1,k2))/1000.eq.6)then
          write(*,'("rdsources_pdata",i10)')nint(blk(1,k2))
          tmp1(2,1:neg)=tmp1(2,1:neg)+saniso(1:neg,k2)*blk(3,k2)
          tmp2(2)=tmp2(2)+blk(3,k2)
        else if(nint(blk(1,k2)).eq.94242)then
c for Pu add spont. fiss. this logic works because there is only one.
          write(*,'("rdsources_pdata",i10)')nint(blk(1,k2))
          tmp1(3,1:neg)=sfiso(1:neg,k2)+saniso(1:neg,k2)*blk(3,k2)
          tmp2(3)=blk(3,k2)
        else if(nint(blk(1,k2)).eq.4009)then
          write(*,'("rdsources_pdata",i10)')nint(blk(1,k2))
          tmp1(4,1:neg)=tmp1(4,1:neg)+saniso(1:neg,k2)*blk(3,k2)
          tmp2(4)=blk(3,k2)
        else if(nint(blk(1,k2))/1000.eq.42)then
          write(*,'("rdsources_pdata",i10)')nint(blk(1,k2))
          tmp1(1,1:neg)=tmp1(1,1:neg)+saniso(1:neg,k2)*blk(3,k2)
          tmp2(1)=tmp2(1)+blk(3,k2)
          if(i_2nd_order.eq.1)then
            do k2p=k3+14,k3+20 ! only Mo
              dqdn2=sum(saniso2(1:neg,k2,k2p))
              tmp3=tmp3+dqdn2*blk(3,k2)*blk(3,k2p)
              tmp4(1,1:neg)=tmp4(1,1:neg)+saniso2(1:neg,k2,k2p)
     1         *blk(3,k2)*blk(3,k2p)
            end do ! k2p
            if(k2.eq.k3+20)then
              tmp3=tmp3/tmp2(1)**2
              tmp4(1,1:neg)=tmp4(1,1:neg)/tmp2(1)**2
            end if
          end if
        end if
      end do ! j
      do j=1,4
        if(tmp2(j).eq.0.d0)cycle
        tmp1(j,1:neg)=tmp1(j,1:neg)/tmp2(j)
      end do ! j
      write(iuo,'(4x,14x,"Mo",15x,"C",11x,"Pu242",14x,"Be")')
      do ng=1,neg
        write(iuo,'(i4,1p4e16.8)')ng,(tmp1(j,ng),j=1,4)
      end do ! ng
      if(i_2nd_order.eq.1.and.nel.eq.53)then
        write(iuo,'("d^2Q/dN^2 for Mo",1pe16.8)')tmp3
        write(iuo,'(4x,14x,"Mo",15x,"C",11x,"Pu242",14x,"Be")')
        do ng=1,neg
          write(iuo,'(i4,1p4e16.8)')ng,tmp4(1,ng),
     1     saniso2(ng,3,3),saniso2(ng,51,51),saniso2(ng,36,36)
        end do ! ng
      else if(i_2nd_order.eq.1.and.nel.eq.35)then
        write(iuo,'("d^2Q/dN^2 for Mo",1pe16.8)')tmp3
        write(iuo,'(4x,14x,"Mo",15x,"C",11x,"Pu242",14x,"Be")')
        do ng=1,neg
          write(iuo,'(i4,1p4e16.8)')ng,saniso2(ng,8,8),
     1     saniso2(ng,2,2),saniso2(ng,33,33),saniso2(ng,20,20)
        end do ! ng
      else if(i_2nd_order.eq.1.and.nel.eq.6)then
        write(iuo,'(4x,14x,"01",14x,"04",14x,"05",12x,"0405",
     1   12x,"0105",13x,"mat")')
        do ng=1,neg
          write(iuo,'(i4,1p40e16.8)')ng,saniso2(ng,1,1),saniso2(ng,4,4),
     1     saniso2(ng,5,5),saniso2(ng,4,5),saniso2(ng,1,5),tmp4(0,ng)
        end do ! ng
      else if(i_2nd_order.eq.1.and.nel.eq.4)then
        write(iuo,'(4x,14x,"01",14x,"02",12x,"0102",13x,"mat")')
        do ng=1,neg
          write(iuo,'(i4,1p40e16.8)')ng,saniso2(ng,1,1),saniso2(ng,2,2),
     1     saniso2(ng,1,2),tmp4(0,ng)
        end do ! ng
      end if
#endif
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
