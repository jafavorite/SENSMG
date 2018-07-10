      subroutine gamadj(niso0)
c set up generalized adjoint source in iterations after the zero'th.
      use COMS
      implicit none
      integer niso0 ! number of isotopes in snxedt
      real*8 gams(neg,0:nmom-1,it,jt),gmomn(neg,0:nmom-1,nitm,njtm),
     1 afgadn(neg,ndir,nitm,njtm),
     2 lam,fissd,xfiss,rat,dd,difmax,gamax,rlip(0:isct)
      integer neg0,ndir0,nitm0,njtm0
      integer i,j,iconv,nconv,ir,jr,l,m,mt,n,lg,ng,nu,ncalcs,ig,nif,njf,
     1 idifmax(4)
      integer mynum,omp_get_thread_num
      real*8 gams_omp(num_threads,neg,it,jt),difmax_omp(num_threads),
     1 gamax_omp(num_threads)
      integer idifmax_omp(num_threads,4)
      character cnote(0:1)*16,c1*1,c2*6
      data cnote/" no convergence.","    convergence."/
c ig=1/2 solve gamma* for reactions/solve gamma* for reaction
c differences
      if(icalc.eq.0.or.icalc.eq.3.or.icalc.eq.4)then
        ig=2
      else if(icalc.eq.1.or.icalc.eq.2)then
        ig=2
      end if
c
      if(epsig.gt.1.001d-5)write(*,'("****",/,
     1 "**** warning.  epsig=",1pe9.2,/"****")')epsig
c
c solve gamma* for reactions
      if(ig.eq.1)then
        ncalcs=nrrx
c solve gamma* for reaction differences
      else if(ig.eq.2)then
        ncalcs=nrrr
      end if
c
c for fixed-source, the adjoint is converged already.
      if(icalc.eq.0.or.icalc.eq.3.or.icalc.eq.4)then
        nconv=ncalcs
        go to 900
      end if
c this is to dimension things correctly in subroutine velip.
c for 2d problems with iangflux=0, the angular fluxes are dimensioned
c with 1's in coms.
      if(iangflux.eq.1)then
        neg0=neg
        ndir0=ndir
        nitm0=nitm
        njtm0=njtm
      else if(iangflux.eq.0)then
        neg0=1
        ndir0=1
        nitm0=1
        njtm0=1
      end if
c
c adjoint moments (regular and generalized) are always calculational
c moments. the sign of rlip is used to convert to inner product moments.
      if(icalc.eq.2)then
        do l=0,isct
          rlip(l)=real(2*l+1)
          if(mod(l,2).ne.0)then
            rlip(l)=-rlip(l)
          end if
c         write(*,'(0p2f8.3)')rlip(l)
        end do ! l
      end if
c
c get <F0* psi*0, psi0>=fissd
      if(icalc.eq.1)then
        call fissip(fissd,amom)
      else if(icalc.eq.2)then
        call velip(fissd,rlip,amom,afadj,neg0,ndir0,nitm0,njtm0)
      end if
c
      if(iangflux.eq.0)then
        c1="l"
        c2="ga_mom"
        nif=it
        njf=jt
      else if(iangflux.eq.1)then
        c1="m"
        c2="ga_af "
        nif=nitm
        njf=njtm
      end if

      if(icalc.eq.0.or.icalc.eq.2.or.icalc.eq.3.or.icalc.eq.4)then
        lam=1.d0
      else if(icalc.eq.1)then
        lam=1.d0/keff
      end if
      nconv=0
      do n=1,ncalcs
        iconv=0 ! convergence flag
c read latest computed generalized adjoint, gmomn or afgadn.
c read into gmom or afgad temporarily, then store in gmomn or afgadn.
        if(iangflux.eq.0)then
          call rddantm(3,n)
          gmomn(1:neg,0:nmom-1,1:it,1:jt)=gmom(1:neg,0:nmom-1,1:it,1:jt)
        else if(iangflux.eq.1)then
          call rddanta(3,n)
          afgadn(1:neg,1:ndir,1:nitm,1:njtm)=
     1     afgad(1:neg,1:ndir,1:nitm,1:njtm)
c construct gmomn
          call calcmom(afgadn,gmomn)
        end if
c
c remove fundamental mode contamination.
        if(icalc.eq.1)then
          call fissip(xfiss,gmomn)
        else if(icalc.eq.2)then
          call velip(xfiss,rlip,gmomn,afgadn,neg0,ndir0,nitm0,njtm0)
        end if
        rat=xfiss/fissd
        if(iangflux.eq.0)then
          gmomn(1:neg,0:nmom-1,1:it,1:jt)=
     1     gmomn(1:neg,0:nmom-1,1:it,1:jt)
     2     -rat*amom(1:neg,0:nmom-1,1:it,1:jt)
        else if(iangflux.eq.1)then
          afgadn(1:neg,1:ndir,1:nitm,1:njtm)=
     1     afgadn(1:neg,1:ndir,1:nitm,1:njtm)
     2     -rat*afadj(1:neg,1:ndir,1:nitm,1:njtm)
c construct gmomn
          call calcmom(afgadn,gmomn)
        end if
        write(iuo,'("iteration",i3," adjoint",i3,
     1   " xfiss,fissd,xfiss/fissd=",1p3e14.5)')iter,n,xfiss,fissd,rat
        if(iter.eq.1)then
          gmom(1:neg,0:nmom-1,1:it,1:jt)=
     1     gmomn(1:neg,0:nmom-1,1:it,1:jt)
          if(iangflux.eq.1)then
            afgad(1:neg,1:ndir,1:nitm,1:njtm)=
     1       afgadn(1:neg,1:ndir,1:nitm,1:njtm)
          end if
c read current sum of generalized adjoints, gmom.
        else if(iter.gt.1)then
          if(iangflux.eq.0)then
            call rdgmom(n,iconv)
          else if(iangflux.eq.1)then
            call rdafgad(n,iconv)
          end if
          if(iconv.eq.1)then
            nconv=nconv+1
            cycle ! do no more with this one
          end if
c check convergence.
          idifmax(1:4)=0
          difmax=0.d0
          idifmax_omp(1:num_threads,1:4)=0
          difmax_omp(1:num_threads)=0.d0
!$OMP PARALLEL
!$OMP& PRIVATE(i,j,ng,l,dd,m,mynum)
!$OMP DO SCHEDULE(DYNAMIC)
          do i=1,nif
            mynum=omp_get_thread_num()+1
            do j=1,njf
              do ng=1,neg
                if(iangflux.eq.0)then
                  do l=0,nmom-1
                    if(gmom(ng,l,i,j).ne.0.d0)then
                      dd=gmomn(ng,l,i,j)/gmom(ng,l,i,j)
                      if(abs(dd).gt.abs(difmax_omp(mynum)))then
                        difmax_omp(mynum)=dd
                        gamax_omp(mynum)=gmom(ng,l,i,j)
                        idifmax_omp(mynum,1)=ng
                        idifmax_omp(mynum,2)=l
                        idifmax_omp(mynum,3)=i
                        idifmax_omp(mynum,4)=j
                      end if
                    end if
                  end do ! l
                else if(iangflux.eq.1)then
                  do m=1,ndir
                    if(afgad(ng,m,i,j).ne.0.d0)then
                      dd=afgadn(ng,m,i,j)/afgad(ng,m,i,j)
                      if(abs(dd).gt.abs(difmax_omp(mynum)))then
                        difmax_omp(mynum)=dd
                        gamax_omp(mynum)=afgad(ng,m,i,j)
                        idifmax_omp(mynum,1)=ng
                        idifmax_omp(mynum,2)=m
                        idifmax_omp(mynum,3)=i
                        idifmax_omp(mynum,4)=j
                      end if
                    end if
                  end do ! l
                end if
              end do ! ng
            end do ! j
          end do ! i
!$OMP END DO
!$OMP END PARALLEL
          do nu=1,num_threads
            if(abs(difmax_omp(nu)).gt.abs(difmax))then
              difmax=difmax_omp(nu)
              gamax=gamax_omp(nu)
              idifmax(1:4)=idifmax_omp(nu,1:4)
            end if
          end do ! nu
c
          if(abs(difmax).le.epsig)then
            iconv=1
            nconv=nconv+1
          end if
          write(iuo,'("iteration",i3," adjoint",i3,
     1     a16," ng,",a1,",i,rel.dif.,",a6,"=",2i4,i6,1pe12.3,e14.5)')
     2     iter,n,cnote(iconv),c1,c2,idifmax(1:3),difmax,gamax
        end if
c
        if(iter.gt.1)then
          if(iangflux.eq.0)then
            gmom(1:neg,0:nmom-1,1:it,1:jt)=
     1       gmom(1:neg,0:nmom-1,1:it,1:jt)
     2       +gmomn(1:neg,0:nmom-1,1:it,1:jt)
          else if(iangflux.eq.1)then
            afgad(1:neg,1:ndir,1:nitm,1:njtm)=
     1       afgad(1:neg,1:ndir,1:nitm,1:njtm)
     2       +afgadn(1:neg,1:ndir,1:nitm,1:njtm)
c construct gmom
            call calcmom(afgad,gmom)
          end if
        end if
c check normalization
        if(icalc.eq.1)then
          call fissip(xfiss,gmom)
        else if(icalc.eq.2)then
          call velip(xfiss,rlip,gmom,afgad,neg0,ndir0,nitm0,njtm0)
        end if
        write(iuo,'("iteration",i3," adjoint",i3,
     1   " norm=",1pe14.5)')iter,n,xfiss
        if(iangflux.eq.0)then
          call wrgmom(n,iconv)
        else if(iangflux.eq.1)then
          call wrafgad(n,iconv)
        end if
c
        gams(1:neg,0:nmom-1,1:it,1:jt)=0.d0
        gams_omp(1:num_threads,1:neg,1:it,1:jt)=0.d0
        if(iconv.eq.0)then
          do ir=1,nr
            do jr=1,nz
              mt=mat(ir,jr)
!$OMP PARALLEL
!$OMP& PRIVATE(i,j,ng,lg,mynum)
!$OMP DO SCHEDULE(DYNAMIC)
              do i=iindex(ir-1)+1,iindex(ir)
                mynum=omp_get_thread_num()+1
                do j=jindex(jr-1)+1,jindex(jr)
                  do ng=1,neg
                    do lg=1,neg
                      gams_omp(mynum,ng,i,j)=gams_omp(mynum,ng,i,j)
     1                 +lam*gmomn(lg,0,i,j)*chi(lg,ng,mt)*nusigf(ng,mt)
                    end do ! lg
                  end do ! ng
                end do ! j
              end do ! i
!$OMP END DO
!$OMP END PARALLEL
            end do ! jr
          end do ! ir
        end if
c
        do nu=1,num_threads
          gams(1:neg,0,1:it,1:jt)=gams(1:neg,0,1:it,1:jt)
     1     +gams_omp(nu,1:neg,1:it,1:jt)
        end do ! nu
c
        call wrdantga(n,iconv,ig,niso0,gams)
c debug print
c     if(n.eq.1)then
c       write(iuo,'("debug print, iteration",i3)')iter
c       l=0
c       do i=1,it
c         write(100+iter,'(i6,1p100e14.6)')i,0.5d0*(rfm(i)+rfm(i-1)),
c    1     gmomn(1:neg,l,i,1),sum(gmomn(1:neg,l,i,1))
c         write(200+iter,'(i6,1p100e14.6)')i,0.5d0*(rfm(i)+rfm(i-1)),
c    1     gams(1:neg,l,i,1),sum(gams(1:neg,l,i,1))
c       end do ! i
c       call stoponerror
c     end if
      end do ! n
c
      write(*,'("iteration",i3,i3," of",i3," generalized ",
     1 "adjoints converged.")')iter,nconv,ncalcs
      write(iuo,'("iteration",i3,i3," of",i3," generalized ",
     1 "adjoints converged.")')iter,nconv,ncalcs
c
c if all gamma* files are converged
  900 if(nconv.eq.ncalcs)then
        ictrl=4
        open(iut,file='stopconverged',status='unknown')
        close (iut)
      end if
c
      return
      end subroutine
