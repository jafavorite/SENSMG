      subroutine calcsens_sm_r(ediv,fissd,lam,r2,ipsrc,ipfis,
     1 ipsrc_s,ipfis_s,mom2src,mom2fis,rlip,rfile,i_no_combine)
c compute sensitivities of Feynman Y or Sm2 to surface locations.
      use COMS
      implicit none
      real*8 ediv(neg),fissd,lam,r2,ipsrc(it),ipfis(it),ipsrc_s(0:nr),
     1 ipfis_s(0:nr),mom2src(neg,0:nm),mom2fis(neg,0:nm),rlip(0:isct)
      character rfile*20
      integer i_no_combine
      integer ncon
      parameter (ncon=4)
      real*8 dtgr(neg,0:nr,nz),dsgr(neg,neg,0:isct,0:nr,nz),
     1 dfgr(neg,neg,0:nr,nz),dqgr(neg,0:nr,nz),ddgr(neg,0:nr,nz),
     2 dtgz(neg,nr,0:nz),dsgz(neg,neg,0:isct,nr,0:nz),
     2 dfgz(neg,neg,nr,0:nz),dqgz(neg,nr,0:nz),ddgz(neg,nr,0:nz),
     3 sintr(neg,0:nr,nz),sintz(neg,nr,0:nz),
     4 stotr(neg,0:nr,nz),sfisr(neg,0:nr,nz),
     5 sinsr(neg,0:nr,nz),sscttr(neg,0:nr,nz),
     6 stotz(neg,nr,0:nz),sfisz(neg,nr,0:nz),
     7 sinsz(neg,nr,0:nz),sscttz(neg,nr,0:nz),
     8 dchigr(neg,neg,0:nr,nz),trcor
      integer i,j,i0,ir,jr,l,lg,m1,m2,ng,n,iw(2),ib(2)
      integer neg0,ndir0,nitm0,njtm0
c
      sintr(1:neg,0:nr,1:nz)=0.d0
      sintz(1:neg,1:nr,0:nz)=0.d0
      stotr(1:neg,0:nr,1:nz)=0.d0
      sfisr(1:neg,0:nr,1:nz)=0.d0
      sinsr(1:neg,0:nr,1:nz)=0.d0
      sscttr(1:neg,0:nr,1:nz)=0.d0
      stotz(1:neg,1:nr,0:nz)=0.d0
      sfisz(1:neg,1:nr,0:nz)=0.d0
      sinsz(1:neg,1:nr,0:nz)=0.d0
      sscttz(1:neg,1:nr,0:nz)=0.d0
      dchigr(1:neg,1:neg,0:nr,1:nz)=0.d0
c
c this is to dimension things correctly in subroutine integrals_r_2d
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
      write(*,'(/,"writing derivatives to file ",a,".")')
     1 trim(rfile)
      write(iuo,'(/,"writing derivatives to file ",a,".")')
     1 trim(rfile)
c
c compute and write sensitivities.
c  0  psi2*,psi and Q2*,psi
c  1  psi1*,(phi1+phis,1)
cc 2  psi1*,phis,1 -- unused because the sum of these is done as n=1
c  3  Qs2*
c
c call integrals with 0/1 for inner products/forward only.
      do n=0,ncon
        if(n.eq.0)then
c 1st-moment forward
          iw(1)=1
          ib(1)=0
c 2nd-moment adjoint
          iw(2)=6
          ib(2)=3
        else if(n.eq.1)then
c <I_1**2,forward> inner product
          amom(1:neg,1:isct,1:it,1:jt)=0.d0
          j=1 ! this is a 1d routine
          do i=1,it
            amom(1:neg,0,i,j)=ipfis(i)*ipfis(i)
          end do ! i
          if(calc2d.eq.0)then ! will not be 1 in this routine
            i0=1
          else if(calc2d.eq.-1)then
            i0=0
          end if
          do ir=i0,nr
            i=iindex(ir)+1
            afadj(1:neg,1:ndir,i,j)=ipfis_s(ir)*ipfis_s(ir)
          end do ! ir
c         write(iuo,'("fis",1p2e17.9)')afadj(1,1,1,j),afadj(1,1,it+1,j)
        else if(n.eq.2)then
c 2nd-moment forward 1+2
          iw(1)=5
          ib(1)=1
c 1st-moment adjoint
          iw(2)=2
          ib(2)=0
        else if(n.eq.3)then
c I_{1,sf}**2 surface integral
          fmom(1:neg,1:isct,1:it,1:jt)=0.d0
          amom(1:neg,1:isct,1:it,1:jt)=0.d0
          j=1 ! this is a 1d routine
          do i=1,it
            fmom(1:neg,0,i,j)=1.d0
            amom(1:neg,0,i,j)=ipsrc(i)*ipsrc(i)
          end do ! i
          afreg(1:neg,1:ndir,1:nitm,1:njtm)=1.d0
          if(calc2d.eq.0)then ! will not be 1 in this routine
            i0=1
          else if(calc2d.eq.-1)then
            i0=0
          end if
          do ir=i0,nr
            i=iindex(ir)+1
            afadj(1:neg,1:ndir,i,j)=ipsrc_s(ir)*ipsrc_s(ir)
          end do ! ir
c         write(iuo,'("src",1p2e17.9)')afadj(1,1,1,j),afadj(1,1,it+1,j)
        else if(n.eq.4)then
c this is experimental.
c 1st-moment forward
          iw(1)=1
          ib(1)=0
c 1st-moment adjoint
          iw(2)=2
          ib(2)=0
        end if ! n
        if(n.eq.0.or.n.eq.2.or.n.eq.4)then
          call rddantm(iw(1),ib(1))
          if(iangflux.eq.1)call rddanta(iw(1),ib(1))
          call rddantm(iw(2),ib(2))
          if(iangflux.eq.1)call rddanta(iw(2),ib(2))
        end if
        if(n.eq.4)then
c this is experimental.
          call integrals_3f_r_1d(lam,afreg,afadj,ipfis_s,dchigr)
        else
          if(calc2d.eq.0.or.calc2d.eq.-1)then
            call integrals_r_1d(0,lam,rlip,fmom,amom,afreg,afadj,
     1       dtgr,dsgr,dfgr,dqgr,ddgr)
          else if(calc2d.eq.1)then
c redundant -- if rplane is used, then zplane is used?
c for old partisn versions using bsleft/bsrite and bsbot/bstop.
            if(iplane.gt.-1.and.jplane.gt.-1)then
              call integrals_r_2d(0,lam,rlip,fmom,amom,afreg,afadj,neg0,
     1         ndir0,nitm0,njtm0,dtgr,dsgr,dfgr,dqgr,ddgr,dtgz,dsgz,
     2         dfgz,dqgz,ddgz)
              call integrals_r_2d_bs(0,lam,rlip,dtgr,dsgr,dfgr,dqgr,
     1         ddgr,dtgz,dsgz,dfgz,dqgz,ddgz)
            else
c for partisn 8_27 and later using asfluxx and asfluxy.
              if(iaflux.eq.1)then
                call integrals_r_2d_as(0,lam,rlip,max(nr,nz),max(it,jt),
     1           dtgr,dsgr,dfgr,dqgr,ddgr,dtgz,dsgz,dfgz,dqgz,ddgz)
              else
                write(*,'("warning. no asleft/asrite, asbott/astop ",
     1           "or afluxx, afluxy.")')
                write(iuo,'("warning. no asleft/asrite, asbott/astop ",
     1           "or afluxx, afluxy.")')
                call integrals_r_2d(0,lam,rlip,fmom,amom,afreg,afadj,
     1           neg0,ndir0,nitm0,njtm0,dtgr,dsgr,dfgr,dqgr,ddgr,dtgz,
     2           dsgz,dfgz,dqgz,ddgz)
              end if
            end if ! iplane,jplane
          end if ! calc2d
        end if ! n
c
        if(n.eq.2)then
c no dqg because this is only for the forward flux.
          dtgr(1:neg,0:nr,1:nz)=2.d0*dtgr(1:neg,0:nr,1:nz)
          dsgr(1:neg,1:neg,0:isct,0:nr,1:nz)=2.d0
     1     *dsgr(1:neg,1:neg,0:isct,0:nr,1:nz)
          dfgr(1:neg,1:neg,0:nr,1:nz)=2.d0*dfgr(1:neg,1:neg,0:nr,1:nz)
          dtgz(1:neg,1:nr,0:nz)=2.d0*dtgz(1:neg,1:nr,0:nz)
          dsgz(1:neg,1:neg,0:isct,1:nr,0:nz)=2.d0
     1     *dsgz(1:neg,1:neg,0:isct,1:nr,0:nz)
          dfgz(1:neg,1:neg,1:nr,0:nz)=2.d0*dfgz(1:neg,1:neg,1:nr,0:nz)
        end if
c
c experimental. cesaro correction
        if(itrcor.eq.3)then
          do l=1,isct
            trcor=real((isct+2)*(isct+1),8)
            trcor=real((isct+2-l)*(isct+1-l),8)/trcor
            dsgr(1:neg,1:neg,l,0:nr,1:nz)=trcor
     1       *dsgr(1:neg,1:neg,l,0:nr,1:nz)
            dsgz(1:neg,1:neg,l,1:nr,0:nz)=trcor
     1       *dsgz(1:neg,1:neg,1,1:nr,0:nz)
          end do ! l
        end if
c
c dtgr = <G*, psi> by energy, radial surfaces
c dsgr = <G*, (2l+1)psi> by energy out, energy in, moment, radial surfaces
c dfgr = <G*, psi/keff> by energy out, energy in, radial surfaces
c dqgr = <G*> by energy, radial surfaces
c dtgz = <G*, psi> by energy, axial surfaces
c dsgz = <G*, (2l+1)psi> by energy out, energy in, moment, axial surfaces
c dfgz = <G*, psi/keff> by energy out, energy in, axial surfaces
c dqgz = <G*> by energy, axial surfaces
c
      if(ilnk3dnt.gt.0)then
        m1=0
        m2=1
      end if
c radial surfaces.
      if(calc2d.eq.0.or.calc2d.eq.1)then ! sphere or cylinder
        i0=1
      else if(calc2d.eq.-1)then ! slab
        i0=0
      end if
      do ir=i0,nr
        do jr=1,nz
          if(ilnk3dnt.eq.0)then
            m1=mat(ir,jr)
            m2=mat(ir+1,jr)
          else if(ilnk3dnt.gt.0)then
            if(calc2d.eq.-1.or.calc2d.eq.0)then
              if(ir.eq.i0)then
                call maclnk3(m1,ir,jr)
              else
                sigt(1:neg,m1)=sigt(1:neg,m2)
                chi(1:neg,1:neg,m1)=chi(1:neg,1:neg,m2)
                nusigf(1:neg,m1)=nusigf(1:neg,m2)
                sigs(1:neg,1:neg,0:isct,m1)=sigs(1:neg,1:neg,0:isct,m2)
                scgr(1:neg,0:isct,m1,1:2)=scgr(1:neg,0:isct,m2,1:2)
              end if
            else if(calc2d.eq.1)then
              call maclnk3(m1,ir,jr)
            end if
            call maclnk3(m2,ir+1,jr)
          end if
          do ng=1,neg
            if(n.eq.0.or.n.eq.2)then
              stotr(ng,ir,jr)=stotr(ng,ir,jr)-dtgr(ng,ir,jr)
     1         *(sigt(ng,m1)-sigt(ng,m2))
              do lg=1,neg
                sfisr(ng,ir,jr)=sfisr(ng,ir,jr)+dfgr(lg,ng,ir,jr)
     1           *(chi(lg,ng,m1)*nusigf(ng,m1)
     2           -chi(lg,ng,m2)*nusigf(ng,m2))
                do l=0,isct
                  sscttr(lg,ir,jr)=sscttr(lg,ir,jr)
     1             +dsgr(ng,lg,l,ir,jr)
     2             *(sigs(ng,lg,l,m1)-sigs(ng,lg,l,m2))
                end do ! l
              end do ! lg
              if(n.eq.0)then
c ignore (alpha,n)?
                sinsr(ng,ir,jr)=sinsr(ng,ir,jr)+dqgr(ng,ir,jr)
     1           *(nsrcf(ng,m1)-nsrcf(ng,m2))
              end if
            else if(n.eq.1)then
              sfisr(ng,ir,jr)=sfisr(ng,ir,jr)+dfgr(1,ng,ir,jr)
     1         *(mom2fis(ng,m1)-mom2fis(ng,m2))
            else if(n.eq.3)then
              sintr(ng,ir,jr)=sintr(ng,ir,jr)+dfgr(1,ng,ir,jr)
     1         *(mom2src(ng,m1)-mom2src(ng,m2))
            else if(n.eq.4)then
c this is experimental. mom2fis is evaluated on the interface
c but for a material; what to do about that?
c             do lg=1,neg
c               sintr(ng,ir,jr)=sintr(ng,ir,jr)+2.d0
c    1           *dchigr(lg,ng,ir,jr)*mom2fis(ng,m1)
c    2           *(chi(lg,ng,m1)-chi(lg,ng,m2))
c             end do ! lg
            end if ! n
          end do ! ng
        end do ! jr
      end do ! ir
      do ir=i0,nr
        do jr=1,nz
          do ng=1,neg
c the terms already have the correct sign.
            if(n.eq.0)then
              sintr(ng,ir,jr)=sintr(ng,ir,jr)+sinsr(ng,ir,jr)
            else if(n.eq.2)then
              sintr(ng,ir,jr)=sintr(ng,ir,jr)+sfisr(ng,ir,jr)
     1         +stotr(ng,ir,jr)+sscttr(ng,ir,jr)
            end if
          end do ! ng
c debug print
c         write(*,'(3i4,1p5e14.6)')n,ir,jr,sum(sintr(1:neg,ir,jr)),
c    1     sum(sinsr(1:neg,ir,jr)),sum(sfisr(1:neg,ir,jr)),
c    2     sum(stotr(1:neg,ir,jr)),sum(sscttr(1:neg,ir,jr))
        end do ! jr
      end do ! ir
c
c axial surfaces.
      if(calc2d.eq.1)then
        do jr=0,nz
          do ir=1,nr
            if(ilnk3dnt.eq.0)then
              m1=mat(ir,jr)
              m2=mat(ir,jr+1)
            else if(ilnk3dnt.gt.0)then
              call maclnk3(m1,ir,jr)
              call maclnk3(m2,ir,jr+1)
            end if
            do ng=1,neg
              stotz(ng,ir,jr)=stotz(ng,ir,jr)-dtgz(ng,ir,jr)
     1         *(sigt(ng,m1)-sigt(ng,m2))
              if(icalc.eq.0)then
c ignore (alpha,n)?
                sinsz(ng,ir,jr)=sinsz(ng,ir,jr)+dqgz(ng,ir,jr)
     1           *(nsrcf(ng,m1)-nsrcf(ng,m2))
              end if
              do lg=1,neg
                sfisz(ng,ir,jr)=sfisz(ng,ir,jr)+dfgz(lg,ng,ir,jr)
     1           *(chi(lg,ng,m1)*nusigf(ng,m1)
     2           -chi(lg,ng,m2)*nusigf(ng,m2))
                do l=0,isct
                  sscttz(lg,ir,jr)=sscttz(lg,ir,jr)
     1             +dsgz(ng,lg,l,ir,jr)
     2             *(sigs(ng,lg,l,m1)-sigs(ng,lg,l,m2))
                end do ! l
              end do ! lg
            end do ! ng
          end do ! ir
        end do ! jr
        do jr=0,nz
          do ir=1,nr
            do ng=1,neg
c the terms already have the correct sign.
c sinsr is 0 except for fixed-source problems.
              sintz(ng,ir,jr)=sinsz(ng,ir,jr)+sfisz(ng,ir,jr)
     1         +stotz(ng,ir,jr)+sscttz(ng,ir,jr)
            end do ! ng
          end do ! ir
        end do ! jr
      end if
c
      end do ! n
c
c for R2 derivatives, divide by 1 (i.e. do nothing).
c     sintr(1:neg,0:nr,1:nz)=sintr(1:neg,0:nr,1:nz)*fissd
c     sintz(1:neg,1:nr,0:nz)=sintz(1:neg,1:nr,0:nz)*fissd
c
c combine with derivatives of R1 to get feyny or sm2.
c set i_no_combine at the top of calcsens_sm.
      if(i_no_combine.eq.0)then
        call smcombine3(r2,sintr,sintz)
      end if
c the derivatives are now derivatives for feyny or sm2.
c
c output
      if(i_no_combine.eq.0)then
        if(icalc.eq.3)then
          write(iur,'(/," feyny derivatives",10x,"feyny=",
     1     1pe15.7)')yasym
        else if(icalc.eq.4)then
          write(iur,'(/," sm2 derivatives",10x,"sm2=",1pe15.7)')sm2
        end if
      else if(i_no_combine.eq.1)then
        write(iur,'(/," r2 derivatives",10x,"r2=",
     1   1pe15.7)')r2
      end if
      if(calc2d.eq.0.or.calc2d.eq.-1)then
        call output3(sintr(1:neg,0:nr,1:nz))
      else if(calc2d.eq.1)then
        call output4(sintr(1:neg,0:nr,1:nz),sintz(1:neg,1:nr,0:nz))
      end if
c
      return
      end subroutine
c
      subroutine smcombine3(r2,sintr,sintz)
c combine r1 and r2 derivatives. input variables are r2
c derivatives. read r1 from a binary file.
      use COMS
      implicit none
      real*8 r2,sintr(1:neg,0:nr,1:nz),sintz(1:neg,1:nr,0:nz)
      real*8 tt,tu
      real*8 sintr_tmp(1:neg,0:nr,1:nz),sintz_tmp(1:neg,1:nr,0:nz)
c
      call rdoutput3_bin(sintr_tmp,sintz_tmp)
c
      if(icalc.eq.3)then
       tt=1.d0/lkg  ! feyny
      else if(icalc.eq.4)then
       tt=0.5d0/lkg ! sm2
      end if
      tu=r2/lkg**2
      sintr(1:neg,0:nr,1:nz)=tt*sintr(1:neg,0:nr,1:nz)
     1 -tu*sintr_tmp(1:neg,0:nr,1:nz)
      if(calc2d.eq.1)then
        sintz(1:neg,1:nr,0:nz)=tt*sintz(1:neg,1:nr,0:nz)
     1   -tu*sintz_tmp(1:neg,1:nr,0:nz)
      end if
      if(icalc.eq.4)then
        sintr(1:neg,0:nr,1:nz)=sintr(1:neg,0:nr,1:nz)/lkg
        if(calc2d.eq.1)then
          sintz(1:neg,1:nr,0:nz)=sintz(1:neg,1:nr,0:nz)/lkg
        end if
      end if
      return
      end subroutine

