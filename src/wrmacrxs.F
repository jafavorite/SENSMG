      subroutine wrmacrxs
c write a macrxs file.
c this routine is for perturbed cross sections but can be used
c for unperturbed also.
      use COMS
      implicit none
      real*8 sigt0(neg,nm+nrrx),nusigf0(neg,nm+nrrx),
     1 chi0(neg,neg,nm+nrrx),siga0(neg,nm+nrrx),
     2 sigs0(neg,neg,0:isct,nm+nrrx),rl(0:isct)
      integer k2,i,j,jr,l,lg,ng,n,ntab
      character afile*10,bfile*10
c maxord is the largest scattering expansion ever expected.
c maxxs is the largest number of principal+edit cross sections.
c maxmult is the largest number of fission multiplicities.
      integer maxord,maxxs,maxmult,nwd
      parameter (maxord=7,maxxs=10,maxmult=20,nwd=20)
      real*8 emax(neg),emin,c(nxs,maxxs),scat(neg*nxs*maxord),
     1 fiss_wgt(nxs),fiss_mult(maxmult,neg,nxs),sf_spectra(neg,nxs),
     2 upgam,temps,eavg(neg),einv(neg),wgts(neg),wgts_bnd(neg+1),
     3 fiss_mat(neg,neg,nxs),sourcm(neg,nxs),tmp1
      character hname*8,huse(2)*8,hmat(nxs)*8,hed(maxxs)*8
      integer*4 ivers,ia(nwd),nedt,ngpb(maxord,nxs),ifsg(maxord,nxs)
c
c include 2l+1 in macrxs.
      do l=0,isct
        rl(l)=real(2*l+1,8)
      end do
c
c zero the cross sections.
      nusigf0(1:neg,1:nm+nrrx)=0.d0
      sigt0(1:neg,1:nm+nrrx)=0.d0
      siga0(1:neg,1:nm+nrrx)=0.d0
      sigs0(1:neg,1:neg,0:isct,1:nm+nrrx)=0.d0
c
c set chi and material cross sections.
c cross sections for nm+1:nm+nrrx will be copied from existing macrxs.
      chi0(1:neg,1:neg,1:nm+nrrx)=chi(1:neg,1:neg,1:nm+nrrx)
c     write(10,'(9999i14)')(ng,ng=1,neg)
c     do ng=1,neg
c       write(10,'(i6,1p9999e14.6)')ng,(chi0(ng,lg,1),lg=1,neg)
c     end do ! ng
c perturb one group and renormalize (full normalization)
      if(fissdata.eq.2)then
c       chi0(152,152,1:nm+nrrx)=chi0(152,152,1:nm+nrrx)*(1.d0-0.2d0)
c       chi0(152,152,1:nm+nrrx)=chi0(152,152,1:nm+nrrx)*(1.d0+0.2d0)
c       chi0(7,7,1:nm+nrrx)=chi0(7,7,1:nm+nrrx)*(1.d0-0.2d0)
        chi0(7,7,1:nm+nrrx)=chi0(7,7,1:nm+nrrx)*(1.d0+0.2d0)
        do i=1,nm+nrrx
          tmp1=sum(chi0(1:neg,1:neg,i))
          chi0(1:neg,1:neg,i)=chi0(1:neg,1:neg,i)/tmp1
        end do ! i
c       write(*,'("chi0(152,152)",1p100e14.6)')chi0(152,152,1:nm+nrrx)
c       write(*,'("chi0(7,7)",1p100e14.6)')chi0(7,7,1:nm+nrrx)
      else
c       chi0(152,1:neg,1:nm+nrrx)=chi0(152,1:neg,1:nm+nrrx)*(1.d0-0.2d0)
c       chi0(152,1:neg,1:nm+nrrx)=chi0(152,1:neg,1:nm+nrrx)*(1.d0+0.2d0)
c       chi0(152,152,1:nm+nrrx)=chi0(152,152,1:nm+nrrx)*(1.d0-0.9d0)
c       chi0(152,152,1:nm+nrrx)=chi0(152,152,1:nm+nrrx)*(1.d0+0.9d0)
c       chi0(7,1:neg,1:nm+nrrx)=chi0(7,1:neg,1:nm+nrrx)*(1.d0-0.2d0)
        chi0(7,1:neg,1:nm+nrrx)=chi0(7,1:neg,1:nm+nrrx)*(1.d0+0.2d0)
c       chi0(7,7,1:nm+nrrx)=chi0(7,7,1:nm+nrrx)*(1.d0-0.5d0)
c       chi0(7,7,1:nm+nrrx)=chi0(7,7,1:nm+nrrx)*(1.d0+0.5d0)
        do i=1,nm+nrrx
          do ng=1,neg
            tmp1=sum(chi0(1:neg,ng,i))
            chi0(1:neg,ng,i)=chi0(1:neg,ng,i)/tmp1
          end do ! ng
c         write(*,'("sum",1pe15.6)')sum(chi0(1:neg,1:neg,i))
        end do ! i
c       write(*,'("chi0(152,152)",1p100e14.6)')chi0(152,152,1:nm+nrrx)
c       write(*,'("chi0(7,7)",1p100e14.6)')chi0(7,7,1:nm+nrrx)
      end if
c
      k2=0
      jr=1
      do i=1,nm
        do j=1,ncb(i)
          k2=k2+1
          do ng=1,neg
            sigt0(ng,i)=sigt0(ng,i)+sigt(ng,nm+k2)*blk(3,k2)
            siga0(ng,i)=siga0(ng,i)+siga(ng,nm+k2)*blk(3,k2)
            nusigf0(ng,i)=nusigf0(ng,i)+nusigf(ng,nm+k2)*blk(3,k2)
            do lg=1,neg
              do l=0,isct
                sigs0(ng,lg,l,i)=sigs0(ng,lg,l,i)+sigs(ng,lg,l,nm+k2)
     1           *blk(3,k2)*rl(l)
              end do ! l
            end do ! lg
          end do ! ng
        end do ! j
      end do ! i
c
      afile="for/macrxs"
      bfile="tmp_macrxs"
      write(*,'("writing ",a,"...")')trim(bfile)
      write(iuo,'("writing ",a,"...")')trim(bfile)
      open(iux,file=afile,status='old',position='rewind',
     1 form='unformatted',access='sequential',
     2 err=900)
      open(iux2,file=bfile,status='unknown',position='rewind',
     1 form='unformatted',access='sequential')
      read(iux)hname,(huse(i),i=1,2),ivers
      write(iux2)hname,(huse(i),i=1,2),ivers
      ia(1:nwd)=0
      read(iux)(ia(i),i=1,nwd)
      write(iux2)(ia(i),i=1,nwd)
      nedt=ia(4)+ia(9)
      read(iux)(hmat(i),i=1,ia(2)),(hed(j),j=1,nedt),
     1 (vel(n),n=1,ia(1)),(emax(n),n=1,ia(1)),emin,temps,
     2 (eavg(n),n=1,ia(1)),(einv(n),n=1,ia(1)),(wgts(n),n=1,ia(1)),
     3 (wgts_bnd(n),n=1,ia(1)+1)
      write(iux2)(hmat(i),i=1,ia(2)),(hed(j),j=1,nedt),
     1 (vel(n),n=1,ia(1)),(emax(n),n=1,ia(1)),emin,temps,
     2 (eavg(n),n=1,ia(1)),(einv(n),n=1,ia(1)),(wgts(n),n=1,ia(1)),
     3 (wgts_bnd(n),n=1,ia(1)+1)
      do ng=1,ia(1)
        read(iux)((c(i,j),i=1,ia(2)),j=1,nedt)
        do n=1,nm
          c(n,1)=chi0(ng,ng,n)
          c(n,2)=nusigf0(ng,n)/rho(n)
          c(n,3)=sigt0(ng,n)/rho(n)
          c(n,4)=siga0(ng,n)/rho(n)
        end do ! n 
        write(iux2)((c(i,j),i=1,ia(2)),j=1,nedt)
        if(ia(3).gt.0)then
          read(iux)((ngpb(l,j),l=1,ia(3)),j=1,ia(2)),
     1     ((ifsg(l,j),l=1,ia(3)),j=1,ia(2))
          write(iux2)((ngpb(l,j),l=1,ia(3)),j=1,ia(2)),
     1     ((ifsg(l,j),l=1,ia(3)),j=1,ia(2))
          ntab=sum(ngpb(1:ia(3),1:ia(2)))
          read(iux)(scat(i),i=1,ntab)
c separate the bands.
c sigs(ng,lg,l,ir) means scattering from lg into ng.
c scgr(ng,l,ir,1/2) are the top and bottom groups scattering into ng.
          i=0
          do j=1,nm
            do l=1,ia(3)
              do lg=ifsg(l,j),ifsg(l,j)-ngpb(l,j)+1,-1
                i=i+1
                if(l-1.le.isct)scat(i)=sigs0(ng,lg,l-1,j)/rho(j)
              end do ! lg
            end do ! l
          end do ! j
          write(iux2)(scat(i),i=1,ntab)
        end if
      end do ! ng
c
      if(ia(11).ne.0)read(iux)(fiss_wgt(i),i=1,ia(2))
      if(ia(11).ne.0)write(iux2)(fiss_wgt(i),i=1,ia(2))
      if(ia(12).ne.0)read(iux)(((fiss_mult(i,ng,n),
     1 i=1,ia(12)),ng=1,ia(1)),n=1,ia(2))
      if(ia(12).ne.0)write(iux2)(((fiss_mult(i,ng,n),
     1 i=1,ia(12)),ng=1,ia(1)),n=1,ia(2))
      if(ia(15).gt.0)read(iux)((sf_spectra(n,i),n=1,ia(1)),
     1 i=1,ia(2))
      if(ia(15).gt.0)write(iux2)((sf_spectra(n,i),n=1,ia(1)),
     1 i=1,ia(2))
      if(ia(17).gt.0)read(iux)((sourcm(n,i),n=1,ia(1)),
     1 i=1,ia(2))
      if(ia(17).gt.0)write(iux2)((sourcm(n,i),n=1,ia(1)),
     1 i=1,ia(2))
      read(iux)upgam
      write(iux2)upgam
      if(ia(13).le.1)then
c chi(ng,lg,nm) is from lg into ng (like sigs).
        read(iux)(((chi(ng,lg,n),n=1,ia(2)),lg=1,ia(1)),
     1   ng=1,ia(1))
        write(iux2)(((chi0(ng,lg,n),n=1,ia(2)),lg=1,ia(1)),
     1   ng=1,ia(1))
      end if
      if(ia(13).eq.0)then
c fiss_mat(ng,lg,nm) is from lg into ng (like sigs).
        read(iux)(((fiss_mat(ng,lg,n),n=1,ia(2)),lg=1,ia(1)),
     1   ng=1,ia(1))
        do ng=1,ia(1)
          do lg=1,ia(1)
            do n=1,ia(2)
              fiss_mat(ng,lg,n)=chi0(ng,lg,n)*nusigf0(lg,n)
            end do ! n
          end do ! lg
        end do ! ng
        write(iux2)(((fiss_mat(ng,lg,n),n=1,ia(2)),lg=1,ia(1)),
     1   ng=1,ia(1))
      end if
c
      close(iux)
      close(iux2)
      go to 999
c
  900 write(*,'("error opening ",a,".")')trim(afile)
      write(iuo,'("error opening ",a,".")')trim(afile)
      call stoponerror
  999 return
      end subroutine
