      subroutine output1_bin(sden,stot,sabs,scap,schi,snu,sfis,
     1 ssct,sscti,sscts,sscto,sins,selas,sinel,ssctt,ssctp1)
c write binary senslx file.
      use COMS
      implicit none
      real*8 sden(neg,2*nel),stot(neg,2*nel),sabs(neg,2*nel),
     1 scap(neg,2*nel),sfis(neg,2*nel),schi(neg,neg,2*nel),
     2 snu(neg,neg,2*nel),ssct(neg,neg,0:isct,2*nel),
     3 sscti(neg,0:isct,2*nel),sscto(neg,0:isct,2*nel),
     4 sscts(neg,0:isct,2*nel),sins(neg,2*nel),selas(neg,2*nel),
     5 sinel(neg,2*nel),ssctt(neg,2*nel),ssctp1(neg,neg,2*nel)
      character hname*8,huse(2)*8
      integer i,l,lg,ng,nel2,nwd
      parameter (nwd=4)
      integer*4 ivers,ia(nwd)
      character afile*20
c
      write(afile,'("senslx")')
      write(hname,'("senslx")')
      write(*,'("writing ",a,"...")')trim(afile)
      write(iuo,'("writing ",a,"...")')trim(afile)
c
c 5.97 opens this file with form='binary',access='sequential' but this
c still works.
      open(iub,file=afile,status='unknown',position='rewind',
#if ( defined INTEL )
     1 form='unformatted',access='stream'
#elif ( defined GNU )
     1 access='stream'  ! untested capability
#endif
     2 )
      if(calc2d.eq.0.or.calc2d.eq.-1)then
        ia(1)=1
      else if(calc2d.eq.1)then
        ia(1)=2
      end if
      ia(2)=neg
      nel2=2*nel
      ia(3)=nel2
      ia(4)=isct+1 ! indexed from 1
      ivers=0
      huse(1:2)=" "
      write(iub)hname,(huse(i),i=1,2),ivers
      write(iub)(ia(i),i=1,nwd)
      do i=1,nel2
        do ng=1,neg
          write(iub)sden(ng,i),stot(ng,i),sabs(ng,i),scap(ng,i),
     1     sfis(ng,i),selas(ng,i),sinel(ng,i),sins(ng,i),ssctt(ng,i)
          write(iub)(schi(ng,lg,i),lg=1,ia(2))
          write(iub)(snu(ng,lg,i),lg=1,ia(2))
          write(iub)(sscti(ng,l,i),l=0,ia(4)-1),
     1     (sscto(ng,l,i),l=0,ia(4)-1),(sscts(ng,l,i),l=0,ia(4)-1)
          write(iub)((ssct(ng,lg,l,i),lg=1,ia(2)),l=0,ia(4)-1)
          if(itrcor.eq.1.or.itrcor.eq.2)then
            write(iub)(ssctp1(ng,lg,i),lg=1,ia(2))
          end if
        end do ! ng
      end do ! i
      close(iub)
c
      return
      end subroutine
c
      subroutine rdoutput1_bin(sden,stot,sabs,scap,schi,snu,sfis,
     1 ssct,sscti,sscts,sscto,sins,selas,sinel,ssctt,ssctp1)
c read binary senslx file.
      use COMS
      implicit none
      real*8 sden(neg,2*nel),stot(neg,2*nel),sabs(neg,2*nel),
     1 scap(neg,2*nel),sfis(neg,2*nel),schi(neg,neg,2*nel),
     2 snu(neg,neg,2*nel),ssct(neg,neg,0:isct,2*nel),
     3 sscti(neg,0:isct,2*nel),sscto(neg,0:isct,2*nel),
     4 sscts(neg,0:isct,2*nel),sins(neg,2*nel),selas(neg,2*nel),
     5 sinel(neg,2*nel),ssctt(neg,2*nel),ssctp1(neg,neg,2*nel)
      character hname*8,huse(2)*8
      integer i,l,lg,ng,nel2,nwd
      parameter (nwd=4)
      integer*4 ivers,ia(nwd)
      character afile*20
c
      write(afile,'("senslx")')
      write(hname,'("senslx")')
      write(*,'("reading ",a,"...")')trim(afile)
      write(iuo,'("reading ",a,"...")')trim(afile)
c
c 5.97 opens this file with form='binary',access='sequential' but this
c still works.
      open(iub,file=afile,status='old',position='rewind',
#if ( defined INTEL )
     1 form='unformatted',access='stream',
#elif ( defined GNU )
     1 access='stream',  ! untested capability
#endif
     2 err=900)
      read(iub,err=910)hname,(huse(i),i=1,2),ivers
      read(iub,err=910)(ia(i),i=1,nwd)
      if(idbgw.eq.1)then
        write(41,'(a)')trim(afile)
        write(41,'(10i8)')(ia(i),i=1,nwd)
      endif
      if((calc2d.eq.0.or.calc2d.eq.-1).and.ia(1).ne.1)go to 795
      if(calc2d.eq.1.and.ia(1).ne.2)go to 797
      if(ia(2).ne.neg)go to 800
      nel2=2*nel
      if(ia(3).ne.nel2)go to 805
      if(ia(4).ne.isct+1)go to 810
c
      do i=1,nel2
        do ng=1,neg
          read(iub)sden(ng,i),stot(ng,i),sabs(ng,i),scap(ng,i),
     1     sfis(ng,i),selas(ng,i),sinel(ng,i),sins(ng,i),ssctt(ng,i)
          read(iub)(schi(ng,lg,i),lg=1,ia(2))
          read(iub)(snu(ng,lg,i),lg=1,ia(2))
          read(iub)(sscti(ng,l,i),l=0,ia(4)-1),
     1     (sscto(ng,l,i),l=0,ia(4)-1),(sscts(ng,l,i),l=0,ia(4)-1)
          read(iub)((ssct(ng,lg,l,i),lg=1,ia(2)),l=0,ia(4)-1)
          if(itrcor.eq.1.or.itrcor.eq.2)then
            read(iub)(ssctp1(ng,lg,i),lg=1,ia(2))
          end if
        end do ! ng
      end do ! i
      close(iub)
      write(*,'("done.")')
      go to 999
c
  795 write(*,'("error. ia(1).ne.1 in ",a,"; ia(1)=",i8)')
     1 trim(afile),ia(1)
      call stoponerror
  797 write(*,'("error. ia(1).ne.2 in ",a,"; ia(1)=",i8)')
     1 trim(afile),ia(1)
      call stoponerror
  800 write(*,'("error. ia(2).ne.neg in ",a,"; ia(2)=",i8)')
     1 trim(afile),ia(2)
      write(*,'("neg=",i8)')neg
      call stoponerror
  805 write(*,'("error. ia(3).ne.2*nel in ",a,"; ia(3)=",i8)')
     1 trim(afile),ia(3)
      write(*,'("2*nel=",i8)')2*nel
      call stoponerror
  810 write(*,'("error. ia(4).ne.isct+1 in ",a,"; ia(4)=",i8)')
     1 trim(afile),ia(4)
      write(*,'("isct+1=",i8)')isct+1
      call stoponerror
  900 write(*,'("error opening ",a,".")')trim(afile)
      call stoponerror
  910 write(*,'("error reading ",a,".")')trim(afile)
      call stoponerror
c
  999 return
      end subroutine
