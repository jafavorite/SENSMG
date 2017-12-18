      subroutine rdsnxedth(ix,iue,niso0,nedt0)
c read header of snxedt to set array dimensions.
      implicit none
      integer ix,iue,niso0,nedt0
      character hname*8,huse(2)*8
      integer*4 ivers,ia(20)
      integer i
      character afile*10
c
      if(ix.eq.1)then
        afile="xs1/snxedt"
      else if(ix.eq.2)then
        afile="for/snxedt"
      end if
      open(iue,file=afile,status='old',position='rewind',
     1 form='unformatted',access='sequential',
     2 err=900)
      read(iue,err=910)hname,(huse(i),i=1,2),ivers
      read(iue,err=910)(ia(i),i=1,10)
      niso0=ia(2)
      nedt0=ia(4)+ia(9)
      close(iue)
      return
c
  900 write(*,'("error opening ",a,".")')trim(afile)
      call stoponerror
  910 write(*,'("error reading ",a,".")')trim(afile)
      call stoponerror
      end subroutine
c
      subroutine rdsnxedt(ix,nx,ne,sigt0,nusigf0,chi0,siga0)
c read binary snxedt file from partisn for microscopic cross sections.
c only read principal cross sections.
c return sigf for all isotopes, edit reactions for reaction-rate
c isotopes
c
      use COMS
      implicit none
c nx is the number of isotopes on snxedt
c ne is the number of edit cross sections
      integer nx,ne,ix
      real*8 velm(neg),emax(neg),emin,c(nx,ne)
      character hname*8,huse(2)*8,hiso(nx)*24,hed(ne)*8
      integer*4 ivers,ia(10),nedt
      integer ng,i,ipf,ipc,j,n,ieliso(nel),irriso(nrrx),ierr
      character afile*10
      real*8 sigt0(neg,nrrx),nusigf0(neg,nrrx),chi0(neg,neg,nrrx),
     1 siga0(neg,nrrx)
c        
      if(ix.eq.1)then
        afile="xs1/snxedt"
      else if(ix.eq.2)then
        afile="for/snxedt"
      end if
      write(*,'("reading multigroup cross sections from ",a,"...")')
     1 trim(afile)
      write(iuo,'("reading multigroup cross sections from ",a,"...")')
     1 trim(afile)
      open(iue,file=afile,status='old',position='rewind',
     1 form='unformatted',access='sequential',
     2 err=900)
      read(iue,err=910)hname,(huse(i),i=1,2),ivers
      read(iue,err=910)(ia(i),i=1,10)
      if(idbgw.eq.1)then
        write(41,'(a)')trim(afile)
        write(41,'(10i8)')(ia(i),i=1,10)
      endif
      nedt=ia(4)+ia(9)
      if(ia(1).ne.neg)go to 795
      if(ia(2).ne.nx)go to 800
      if(ia(6).ne.lng)go to 810
      if(nedt.ne.ne)go to 815
      if(nedt.gt.maxedits)go to 820
c
      ieliso(1:nel)=0
      irriso(1:nrrx)=0
      ipf=0
      ipc=0
      ierr=0
      read(iue,err=910)(hiso(i),i=1,ia(2)),(hed(j),j=1,nedt),
     1 (velm(n),n=1,ia(1)),(emax(n),n=1,ia(1)),emin
c     write(*,'(a10,2x,"hiso")')afile
      do i=1,ia(2)
c       write(*,'(i6,2x,a)')i,hiso(i)
c this equality works because zaidfull and hiso are both length 24.
c ieliso is used to find sigf and sigc for transport isotopes.
        do n=1,nel
          if(adjustl(zaidfull(n)).eq.hiso(i))ieliso(n)=i
        end do ! n
c irriso is used to find cross sections for reaction rates.
        do n=1,nrrx
          if(adjustl(zaidfull(nel+n)).eq.hiso(i))irriso(n)=i
        end do ! n
      end do ! i
c     do n=1,nrrx
c       write(*,'(i6,2x,a)')irriso(n),zaidfull(nel+n)
c     end do ! n
      if(iter.eq.0.and.ictrl.eq.3)then
        write(*,'(/,2x,"available reactions (hed)")')
        write(iuo,'(/,2x,"available reactions (hed)")')
      end if
      do i=1,nedt
        cedits(i)=hed(i)
        if(hed(i)(1:6).eq."n-fiss".or.hed(i)(1:6).eq."(n,F) ")then
          ipf=i
        else if(hed(i)(1:5).eq."mend5".or.hed(i)(1:5).eq."(n,g)")then
          ipc=i
        end if
        if(iter.eq.0.and.ictrl.eq.3)then
          write(*,'(i5,2x,a)')i,trim(hed(i))
          write(iuo,'(i5,2x,a)')i,trim(hed(i))
        end if
      end do ! i
      do n=1,nel
        if(ieliso(n).gt.0)cycle
        write(*,'("error. no zaid found in snxedt: ",a)')
     1   trim(zaidfull(n))
        write(iuo,'("error. no zaid found in snxedt: ",a)')
     1   trim(zaidfull(n))
        ierr=1
      end do ! i
      do n=1,nrrx
        if(irriso(n).gt.0)cycle
        write(*,'("error. no zaid found in snxedt: ",a)')
     1   trim(zaidfull(nel+n))
        write(iuo,'("error. no zaid found in snxedt: ",a)')
     1   trim(zaidfull(nel+n))
        ierr=1
      end do ! i
      if(ipf.eq.0)then
        write(*,'("warning. no n-fiss or (n,F) found in snxedt.")')
        write(iuo,'("warning. no n-fiss or (n,F) found in snxedt.")')
      end if
      if(ipc.eq.0)then
        write(*,'("warning. no mend5 or (n,g) found in snxedt.")')
        write(iuo,'("warning. no mend5 or (n,g) found in snxedt.")')
      end if
      if(ierr.eq.1)call stoponerror
c
      do ng=1,ia(1)
        read(iue,err=910)((c(i,j),i=1,ia(2)),j=1,nedt)
c test what is in the file.
c       do i=1,ia(2)
c         do j=1,nedt
c           if(c(i,j).gt.0.d0)write(19+ix,'(3i6,1pe14.6)')ng,i,j,c(i,j)
c         end do ! j
c       end do ! i
c this is getting sigf and sigc for isotopes; nm is the offset into the
c isotope cross sections.
        do n=1,nel
          if(ipf.gt.0)then
            sigf(ng,nm+n)=c(ieliso(n),ipf)
          end if
          if(ipc.gt.0)then
            sigc(ng,nm+n)=c(ieliso(n),ipc)
          end if
        end do ! n
c this is getting reaction-rate cross sections. the sig0 cross
c sections aren't needed; they are gathered here for comparison
c later if desired.
        do n=1,nrrx
          chi0(ng,ng,n)=c(irriso(n),1)
          nusigf0(ng,n)=c(irriso(n),2)
          sigt0(ng,n)=c(irriso(n),3)
          siga0(ng,n)=c(irriso(n),4)
          rrxs(ng,n)=c(irriso(n),irrx(2,n))
        end do ! n
      end do ! ng
c
      close(iue)
      goto 999
c
  795 write(*,'("error. ia(1).ne.neg in ",a,"; ia(1)=",i8)')
     1 trim(afile),ia(1)
      call stoponerror
  800 write(*,'("error. ia(2).ne.nx in ",a,"; ia(2)=",i8)')
     1 trim(afile),ia(2)
      call stoponerror
  810 write(*,'("error. ia(6).ne.lng in ",a,"; ia(6)=",i8)')
     1 trim(afile),ia(6)
      call stoponerror
  815 write(*,'("error. nedt.ne.ne in ",a,"; nedt=",i8)')
     1 trim(afile),nedt
      call stoponerror
  820 write(*,'("error. nedt.gt.maxedits in ",a,"; nedt=",i8)')
     1 trim(afile),nedt
      write(*,'("increase maxedits in coms.F90.")')
      call stoponerror
  900 write(*,'("error opening ",a,".")')trim(afile)
      call stoponerror
  910 write(*,'("error reading ",a,".")')trim(afile)
      call stoponerror
  920 write(*,'("error reading ",a," for scattering.")')trim(afile)
      call stoponerror
  940 write(*,'("error reading ",a," for chi matrix.")')trim(afile)
      call stoponerror
c
  999 return
      end subroutine
