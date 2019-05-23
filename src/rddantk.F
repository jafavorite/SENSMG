      subroutine rddantk
c read ascii output for keff, alpha, or leakage.
      use COMS
      implicit none
      real*8 keffadj,alphaadj,lkgadj,dum
      integer i,j,l,n,ios,nrspf
      character afile*11,line*130,resp(0:4)*5
      data resp/"lkg  ","keff ","alpha","feyny","sm2  "/
c
      afile="for/for_out"
      write(*,'("reading ",a,"...")')trim(afile)
      write(iuo,'("reading ",a,"...")')trim(afile)
      open(iuk,file=afile,status='old',position='rewind',err=900)
   10 read(iuk,'(a)',iostat=ios)line
      if(ios.ne.0)go to 910
      if(icalc.eq.0.or.icalc.eq.3.or.icalc.eq.4)then
        if(line(3:25).ne."key start balance table")go to 10
c right leakage for slab and sphere, total for cylinder.
        if(calc2d.eq.-1.or.calc2d.eq.0)then
          do i=1,neg+12
            read(iuk,'(a)',iostat=ios)line
          end do ! i
          do i=1,neg
            read(iuk,'(6x,e16.7)',iostat=ios)dum
            lkg=lkg+dum*deteff(i)
          end do ! i
        else if(calc2d.eq.1)then
          do i=1,5
            read(iuk,'(a)',iostat=ios)line
          end do ! i
          do i=1,neg
            read(iuk,'(112x,e18.7)',iostat=ios)dum
            lkg=lkg+dum*deteff(i)
          end do ! i
        end if
        if(icalc.eq.0)then
          write(*,'("leakage=",1pe14.7)')lkg
          write(iuo,'("leakage=",1pe14.7)')lkg
        else if(icalc.eq.3.or.icalc.eq.4)then
          write(*,'("count rate=",1pe14.7)')lkg
          write(iuo,'("count rate=",1pe14.7)')lkg
        end if
      else if(icalc.eq.1)then
        if(line(2:23).ne."summary integral-k-eff")go to 10
        read(line(24:41),*)keff
        if(keff.ge.0.01d0)then
          write(*,'("keff=",0pf12.8)')keff
          write(iuo,'("keff=",0pf12.8)')keff
        else
          write(*,'("keff=",1pe15.7)')keff
          write(iuo,'("keff=",1pe15.7)')keff
        end if
      else if(icalc.eq.2)then
        if(line(2:23).ne."summary integral-alpha")go to 10
        read(line(24:41),*)alpha
        read(iuk,'(a)',iostat=ios)line ! integral-neutrons
        write(*,'("alpha=",1pe15.7)')alpha
        write(iuo,'("alpha=",1pe15.7)')alpha
        if(alpha.gt.0.d0.and.aflxfrm.eq.1)then
          write(*,'("warning. aflxfrm=1 is recommended only for ",
     1     "subcritical problems.")')
          write(iuo,'("warning. aflxfrm=1 is recommended only for ",
     1     "subcritical problems.")')
        end if
      end if
      write(*,'("iteration=",i6)')iter
      close(iuk)
c
c read adjoint keff, alpha, or leakage
      if(iter.eq.1.or.(iter.eq.0.and.ictrl.eq.4))then
        afile="adj/adj_out"
        write(*,'("reading ",a,"...")')trim(afile)
        write(iuo,'("reading ",a,"...")')trim(afile)
        open(iuk,file=afile,status='old',position='rewind',err=900)
   20   read(iuk,'(a)',iostat=ios)line
        if(ios.ne.0)go to 910
        if(icalc.eq.0.or.icalc.eq.3.or.icalc.eq.4)then
          if(line(50:92).ne.
     1     "zone edit for the sum of the neutron groups")go to 20
          nrspf=0
          do i=1,nr
            do j=1,nz
              nrspf=nrspf+ismat(mat(i,j))
            end do ! j
          end do ! i
          lkgadj=0.d0
   30     read(iuk,'(a)',iostat=ios)line
          l=index(line,"rspf") ! the first instance
          if(l.eq.0)go to 30
          l=l-5
          n=min(7,nrspf)
          do i=1,nr*nz+3
            read(iuk,'(a)',iostat=ios)line
          end do ! i
          do i=1,n
            read(line(l:l+13),*)dum
            lkgadj=lkgadj+dum
            l=l+14
          end do ! i
          nrspf=nrspf-n
          if(nrspf.gt.0)go to 30
c
          if(icalc.eq.0)then
            write(*,'("forward_leakage=",1pe15.7,"  adjoint_leakage=",
     1       e15.7,"  adjoint/forward=",0pf12.8)')lkg,lkgadj,lkgadj/lkg
            write(iuo,'("forward_leakage=",1pe15.7,"  adjoint_leakage=",
     1       e15.7,"  adjoint/forward=",0pf12.8)')lkg,lkgadj,lkgadj/lkg
          else if(icalc.eq.3.or.icalc.eq.4)then
            write(*,'("forward_count_rate=",1pe15.7,
     1       " adjoint_count_rate=",
     2       e15.7,"  adjoint/forward=",0pf12.8)')lkg,lkgadj,lkgadj/lkg
            write(iuo,'("forward_count_rate=",1pe15.7,
     1       " adjoint_count_rate=",
     2       e15.7,"  adjoint/forward=",0pf12.8)')lkg,lkgadj,lkgadj/lkg
          end if
        else if(icalc.eq.1)then
          if(line(2:23).ne."summary integral-k-eff")go to 20
          read(line(24:41),*)keffadj
          write(*,'("forward_keff=",0pf12.8,"  adjoint_keff=",f12.8,
     1     "  adjoint/forward=",f12.8)')keff,keffadj,keffadj/keff
          write(iuo,'("forward_keff=",0pf12.8,"  adjoint_keff=",f12.8,
     1     "  adjoint/forward=",f12.8)')keff,keffadj,keffadj/keff
        else if(icalc.eq.2)then
          if(line(2:23).ne."summary integral-alpha")go to 20
          read(line(24:41),*)alphaadj
          write(*,'("forward_alpha=",1pe15.7,"  adjoint_alpha=",e15.7,
     1     "  adjoint/forward=",0pf12.8)')alpha,alphaadj,alphaadj/alpha
          write(iuo,'("forward_alpha=",1pe15.7,"  adjoint_alpha=",e15.7,
     1     "  adjoint/forward=",0pf12.8)')alpha,alphaadj,alphaadj/alpha
        end if
        close(iuk)
      end if
      go to 999
c
  900 write(*,'("error opening ",a,".")')trim(afile)
      write(iuo,'("error opening ",a,".")')trim(afile)
      call stoponerror
  910 write(*,'("error finding ",a," in ",a,".")')trim(resp(icalc)),
     1 trim(afile)
      write(iuo,'("error finding ",a," in ",a,".")')trim(resp(icalc)),
     1 trim(afile)
c no need to stop if it's the adjoint.
      if(trim(afile).eq."for/for_out")then
        call stoponerror
      end if
c
c reset icalc.
  999 return
      end subroutine
