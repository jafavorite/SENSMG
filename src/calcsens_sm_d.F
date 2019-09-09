      subroutine calcsens_sm_d(ediv,fissd,lam,r2,ipsrc,ipfis,mom2src,
     1 mom2fis,rlip,rfile,nm1,i_no_combine)
c compute sensitivities of Feynman Y or Sm2 to material densities.
      use COMS
      implicit none
      integer ncon
      parameter (ncon=3)
      integer nm1,i_no_combine
      real*8 ediv(neg),fissd,lam,r2,ipsrc(it),ipfis(it),
     1 mom2src(neg,0:nm),mom2fis(neg,0:nm),rlip(0:isct)
      character rfile*20
      real*8 dtg(neg,nr,nz),dsg(neg,neg,0:isct,nr,nz),
     1 dfg(neg,neg,nr,nz),dqg(neg,nr,nz),
     2 sden(neg,nm1),stot(neg,nm1),sfis(neg,nm1),sins(neg,nm1),
     3 ssctt(neg,nm1),trcor
      integer i,ir,iso,j,jr,l,lg,k2,ng,m1,m2,nm2,n,iw(2),ib(2)
      integer neg0,ndir0,nitm0,njtm0
c
c if ilnk3dnt=0, nm1=nm; if ilnk3dnt>0, nm1=nr*nz (all meshes)
      sden(1:neg,1:nm1)=0.d0
      stot(1:neg,1:nm1)=0.d0
      sfis(1:neg,1:nm1)=0.d0
      sins(1:neg,1:nm1)=0.d0
      ssctt(1:neg,1:nm1)=0.d0
c
c this is to dimension things correctly in subroutine integrals
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
c compute and write sensitivities.
      write(*,'(/,"writing sensitivities to file ",a,".")')
     1 trim(rfile)
      write(iuo,'(/,"writing sensitivities to file ",a,".")')
     1 trim(rfile)
c
c need material isct+1 scattering for transport corrections.
      if(itrcor.eq.1.or.itrcor.eq.2)then
        k2=0
        do i=1,nm
          do j=1,ncb(i)
            k2=k2+1
            iso=nm+k2
            sigsp1(1:neg,1:neg,i)=sigsp1(1:neg,1:neg,i)
     1       +sigsp1(1:neg,1:neg,iso)*blk(3,k2)
          end do ! j
        end do ! i
      end if
c
      do n=0,ncon
c call integrals with 0/1 for inner products/forward only.
        if(n.eq.0)then
c 1st-moment forward
          iw(1)=1
          ib(1)=0
c 2nd-moment adjoint
          iw(2)=6
          ib(2)=3
        else if(n.eq.1)then
c <I_1**2,forward> inner product
          amom(1:neg,0:isct,1:it,1:jt)=0.d0
          do i=1,it
            amom(1:neg,0,i,1)=ipfis(i)*ipfis(i)
          end do ! i
        else if(n.eq.2)then
c 2nd-moment forward 1+2
          iw(1)=5
          ib(1)=1
c 1st-moment adjoint
          iw(2)=2
          ib(2)=0
        else if(n.eq.3)then
c I_{1,sf}**2 volume integral
          fmom(1:neg,0:isct,1:it,1:jt)=0.d0
          amom(1:neg,0:isct,1:it,1:jt)=0.d0
          do i=1,it
            amom(1:neg,0,i,1)=ipsrc(i)*ipsrc(i)
            fmom(1:neg,0,i,1)=1.d0
          end do ! i
        end if ! n
        if(n.eq.0.or.n.eq.2) then
          call rddantm(iw(1),ib(1))
          if(iangflux.eq.1)call rddanta(iw(1),ib(1))
          call rddantm(iw(2),ib(2))
          if(iangflux.eq.1)call rddanta(iw(2),ib(2))
        end if
        call integrals(0,lam,rlip,fmom,amom,afreg,afadj,
     1   neg0,ndir0,nitm0,njtm0,dtg,dsg,dfg,dqg)
c
        if(n.eq.2)then
c no dqg because this is only for the forward flux.
          dtg(1:neg,1:nr,1:nz)=2.d0*dtg(1:neg,1:nr,1:nz)
          dsg(1:neg,1:neg,0:isct,1:nr,1:nz)=2.d0
     1     *dsg(1:neg,1:neg,0:isct,1:nr,1:nz)
          dfg(1:neg,1:neg,1:nr,1:nz)=2.d0*dfg(1:neg,1:neg,1:nr,1:nz)
        end if
c
c cesaro transport correction
        if(itrcor.eq.3)then
          do l=1,isct
            trcor=real((isct+2)*(isct+1),8)
            trcor=real((isct+2-l)*(isct+1-l),8)/trcor
            dsg(1:neg,1:neg,l,1:nr,1:nz)=trcor
     1       *dsg(1:neg,1:neg,l,1:nr,1:nz)
          end do ! l
        end if
c
c this logic catches a material used in different regions.
c if ilnk3dnt>0, ir and jr cycle through the materials.
c dtg = <G*, psi> by energy, coarse mesh
c dsg = <G*, (2l+1)psi> by energy out, energy in, moment, coarse mesh
c dfg = <G*, psi/keff> by energy out, energy in, coarse mesh
c dqg = <G*> by energy, coarse mesh
        if(ilnk3dnt.eq.0)then
          nm2=nm
        else if(ilnk3dnt.gt.0)then
          nm2=1
          m1=0
        end if
        m2=0
        do i=1,nm2
          do ir=1,nr
            do jr=1,nz
              if(ilnk3dnt.eq.0)then
                if(mat(ir,jr).ne.i)cycle
                m1=i
                m2=i
              else if(ilnk3dnt.gt.0)then
                m2=m2+1
                if(idclnk3(ir,jr,1).eq.0)cycle
                call maclnk3(m1,ir,jr)
              end if
              do ng=1,neg
                if(n.eq.0.or.n.eq.2)then
                  stot(ng,m2)=stot(ng,m2)-dtg(ng,ir,jr)*sigt(ng,m1)
c diagonal and bhs transport corrections
                  if(itrcor.eq.1)then
                    stot(ng,m2)=stot(ng,m2)+dtg(ng,ir,jr)
     1               *sigsp1(ng,ng,m1)
                  else if(itrcor.eq.2)then
                    stot(ng,m2)=stot(ng,m2)+dtg(ng,ir,jr)
     1               *sum(sigsp1(1:neg,ng,m1))
                  end if
                  do lg=1,neg
                    sfis(ng,m2)=sfis(ng,m2)+dfg(lg,ng,ir,jr)
     1               *chi(lg,ng,m1)*nusigf(ng,m1)
                    do l=0,isct
                      ssctt(lg,m2)=ssctt(lg,m2)+dsg(ng,lg,l,ir,jr)
     1                 *sigs(ng,lg,l,m1)
c diagonal and bhs transport corrections
                      if((itrcor.eq.1.and.lg.eq.ng).or.itrcor.eq.2)then
                        ssctt(lg,m2)=ssctt(lg,m2)-dsg(ng,lg,l,ir,jr)
     1                   *sigsp1(ng,lg,m1)
                      end if
                    end do ! l
                  end do ! lg
                  if(n.eq.0)then
c ignore (alpha,n)?
                    sins(ng,m2)=sins(ng,m2)+dqg(ng,ir,jr)*nsrcf(ng,m1)
                  end if
                else if(n.eq.1)then
                  sden(ng,m2)=sden(ng,m2)+dfg(1,ng,ir,jr)*mom2fis(ng,m1)
                else if(n.eq.3)then
                  sden(ng,m2)=sden(ng,m2)+dfg(1,ng,ir,jr)*mom2src(ng,m1)
                end if
              end do ! ng
            end do ! jr
          end do ! ir
        end do ! i
        do i=1,nm1
          do ng=1,neg
c the terms already have the correct sign.
c sins is 0 except for fixed-source problems.
            if(n.eq.0)then
              sden(ng,i)=sden(ng,i)+sins(ng,i)
            else if(n.eq.2)then
              sden(ng,i)=sden(ng,i)+sfis(ng,i)+stot(ng,i)+ssctt(ng,i)
            end if
          end do ! ng
        end do ! i
c
      end do ! n
c
c for R2 sensitivities, divide by R2
      sden(1:neg,1:nm1)=sden(1:neg,1:nm1)*fissd
c
c combine with sensitivities of R1 to get feyny or sm2.
c set i_no_combine at the top of calcsens_sm.
      if(i_no_combine.eq.0)then
        call smcombine2(nm1,sden)
      end if
c the sensitivities are now relative sensitivities for feyny or sm2.
c
c output
      if(i_no_combine.eq.0)then
        if(icalc.eq.3)then
          write(iur,'(/," feyny sensitivities",10x,"feyny=",
     1     1pe15.7,/)')yasym
        else if(icalc.eq.4)then
          write(iur,'(/," sm2 sensitivities",10x,"sm2=",1pe15.7,/)')sm2
        end if
      else if(i_no_combine.eq.1)then
        write(iur,'(/," r2 sensitivities",10x,"r2=",
     1   1pe15.7)')r2
      end if
      call output2(nm1,sden(1:neg,1:nm1))
c
      return
      end subroutine
c
      subroutine smcombine2(nm1,sden)
c combine r1 and r2 sensitivities. input variables are r2
c sensitivities. read r1 from a binary file.
      use COMS
      implicit none
      integer nm1
      real*8 sden(neg,nm1)
      real*8 tt
      real*8 sden_tmp(neg,nm1)
c
      call rdoutput2_bin(nm1,sden_tmp)
c
      if(icalc.eq.3)then
       tt=1.d0 ! feyny
      else if(icalc.eq.4)then
       tt=2.d0 ! sm2
      end if
      sden(1:neg,1:nm1)=sden(1:neg,1:nm1)-tt*sden_tmp(1:neg,1:nm1)
      return
      end subroutine
