      subroutine calcsens_r(nn,ediv,fissd,lam,rlip,rfile,
     1 rxnrate,rxnratet)
c compute sensitivities to surface locations.
      use COMS
      implicit none
      integer nn
      real*8 ediv(neg),fissd,lam,rlip(0:isct),
     1 rxnrate(neg,nrrx),rxnratet(nrrx+nflux)
      character rfile*20
      real*8 dtgr(neg,0:nr,nz),dsgr(neg,neg,0:isct,0:nr,nz),
     1 dfgr(neg,neg,0:nr,nz),dqgr(neg,0:nr,nz),ddgr(neg,0:nr,nz),
     2 dtgz(neg,nr,0:nz),dsgz(neg,neg,0:isct,nr,0:nz),
     2 dfgz(neg,neg,nr,0:nz),dqgz(neg,nr,0:nz),ddgz(neg,nr,0:nz),
     3 sintr(0:nrrr,neg,0:nr,nz),sintz(0:nrrr,neg,nr,0:nz),
     4 stotr(0:nrrr,neg,0:nr,nz),sfisr(0:nrrr,neg,0:nr,nz),
     5 sinsr(0:nrrr,neg,0:nr,nz),sscttr(0:nrrr,neg,0:nr,nz),
     6 stotz(0:nrrr,neg,nr,0:nz),sfisz(0:nrrr,neg,nr,0:nz),
     7 sinsz(0:nrrr,neg,nr,0:nz),sscttz(0:nrrr,neg,nr,0:nz),
     8 st(neg),trcor,trcord
      integer i0,ir,jr,i1,i2,j1,j2,k,l,lg,m1,m2,ng
      integer neg0,ndir0,nitm0,njtm0
c
      sintr(0:nrrr,1:neg,0:nr,1:nz)=0.d0
      sintz(0:nrrr,1:neg,1:nr,0:nz)=0.d0
      stotr(0:nrrr,1:neg,0:nr,1:nz)=0.d0
      sfisr(0:nrrr,1:neg,0:nr,1:nz)=0.d0
      sinsr(0:nrrr,1:neg,0:nr,1:nz)=0.d0
      sscttr(0:nrrr,1:neg,0:nr,1:nz)=0.d0
      stotz(0:nrrr,1:neg,1:nr,0:nz)=0.d0
      sfisz(0:nrrr,1:neg,1:nr,0:nz)=0.d0
      sinsz(0:nrrr,1:neg,1:nr,0:nz)=0.d0
      sscttz(0:nrrr,1:neg,1:nr,0:nz)=0.d0
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
c nn=0 is for keff, alpha, or leakage sensitivities,
c nn > 0 is for reaction-rate sensitivities.
      if(nn.le.1)then
        write(*,'(/,"writing derivatives to file ",a,".")')
     1   trim(rfile)
        write(iuo,'(/,"writing derivatives to file ",a,".")')
     1   trim(rfile)
      end if
c
      if(nn.eq.0)then
        if(icalc.eq.0.or.icalc.eq.3.or.icalc.eq.4)then
          fissd=1.d0
        else if(icalc.eq.1)then
          fissd=keff*fissd
        else if(icalc.eq.2)then
          fissd=alpha*fissd
        end if
      end if
c
c call integrals with 0/1 for inner products/forward only.
      if(nn.eq.0)then
        if(calc2d.eq.0.or.calc2d.eq.-1)then
          call integrals_r_1d(0,lam,rlip,fmom,amom,afreg,afadj,
     1     dtgr,dsgr,dfgr,dqgr,ddgr)
        else if(calc2d.eq.1)then
c redundant -- if rplane is used, then zplane is used?
c for old partisn versions using bsleft/bsrite and bsbot/bstop.
          if(iplane.gt.-1.and.jplane.gt.-1)then
            call integrals_r_2d(0,lam,rlip,fmom,amom,afreg,afadj,neg0,
     1       ndir0,nitm0,njtm0,dtgr,dsgr,dfgr,dqgr,ddgr,dtgz,dsgz,dfgz,
     2       dqgz,ddgz)
            call integrals_r_2d_bs(nn,lam,rlip,dtgr,dsgr,dfgr,dqgr,ddgr,
     1       dtgz,dsgz,dfgz,dqgz,ddgz)
          else
c for partisn 8_27 and later using asfluxx and asfluxy.
            if(iaflux.eq.1)then
              call integrals_r_2d_as(nn,lam,rlip,max(nr,nz),max(it,jt),
     1         dtgr,dsgr,dfgr,dqgr,ddgr,dtgz,dsgz,dfgz,dqgz,ddgz)
            else
              write(*,'("warning. no asleft/asrite, asbott/astop ",
     1         "or afluxx, afluxy.")')
              write(iuo,'("warning. no asleft/asrite, asbott/astop ",
     1         "or afluxx, afluxy.")')
              call integrals_r_2d(0,lam,rlip,fmom,amom,afreg,afadj,neg0,
     1         ndir0,nitm0,njtm0,dtgr,dsgr,dfgr,dqgr,ddgr,dtgz,dsgz,
     2         dfgz,dqgz,ddgz)
            end if
          end if
        end if
      else ! nn
        if(calc2d.eq.0.or.calc2d.eq.-1)then
          call integrals_r_1d(0,lam,rlip,fmom,gmom,afreg,afgad,
     1     dtgr,dsgr,dfgr,dqgr,ddgr)
        else if(calc2d.eq.1)then
c use this until asfluxx/asflux are written for generalized adjoints.
          call integrals_r_2d(0,lam,rlip,fmom,gmom,afreg,afgad,neg0,
     1     ndir0,nitm0,njtm0,dtgr,dsgr,dfgr,dqgr,ddgr,dtgz,dsgz,dfgz,
     2     dqgz,ddgz)
c redundant -- if rplane is used, then zplane is used?
c for old partisn versions using bsleft/bsrite and bsbot/bstop.
c         if(iplane.gt.-1.and.jplane.gt.-1)then
c           call integrals_r_2d(0,lam,rlip,fmom,gmom,afreg,afgad,neg0,
c    1       ndir0,nitm0,njtm0,dtgr,dsgr,dfgr,dqgr,ddgr,dtgz,dsgz,dfgz,
c    2       dqgz,ddgz)
c           call integrals_r_2d_bs(nn,lam,rlip,dtgr,dsgr,dfgr,dqgr,ddgr,
c    1       dtgz,dsgz,dfgz,dqgz,ddgz)
c         else
c for partisn 8_27 and later using asfluxx and asfluxy.
c           if(iaflux.eq.1)then
c             call integrals_r_2d_as(nn,lam,rlip,max(nr,nz),max(it,jt),
c    1         dtgr,dsgr,dfgr,dqgr,ddgr,dtgz,dsgz,dfgz,dqgz,ddgz)
c           else
c             write(*,'("warning. no asleft/asrite, asbott/astop ",
c    1         "or afluxx, afluxy.")')
c             write(iuo,'("warning. no asleft/asrite, asbott/astop ",
c    1         "or afluxx, afluxy.")')
c             call integrals_r_2d(0,lam,rlip,fmom,gmom,afreg,afgad,neg0,
c    1         ndir0,nitm0,njtm0,dtgr,dsgr,dfgr,dqgr,ddgr,dtgz,dsgz,
c    2         dfgz,dqgz,ddgz)
c           end if
c         end if ! iplane,jplane
        end if ! calc2d
      end if ! nn
c
c cesaro transport correction
      if(itrcor.eq.3)then
        trcord=real((isct+2)*(isct+1),8)
        do l=1,isct
          trcor=real((isct+2-l)*(isct+1-l),8)/trcord
          dsgr(1:neg,1:neg,l,0:nr,1:nz)=trcor
     1     *dsgr(1:neg,1:neg,l,0:nr,1:nz)
          dsgz(1:neg,1:neg,l,1:nr,0:nz)=trcor
     1     *dsgz(1:neg,1:neg,l,1:nr,0:nz)
        end do ! l
      end if
c
c dtgr = <G*, psi> by energy, radial surfaces
c dsgr = <G*, (2l+1)psi> by energy out, energy in, moment, radial surfaces
c dfgr = <G*, psi/keff> by energy out, energy in, radial surfaces
c dqgr = <G*> by energy, radial surfaces
c ddgr = <G> by energy, radial surfaces
c dtgz = <G*, psi> by energy, axial surfaces
c dsgz = <G*, (2l+1)psi> by energy out, energy in, moment, axial surfaces
c dfgz = <G*, psi/keff> by energy out, energy in, axial surfaces
c dqgz = <G*> by energy, axial surfaces
c ddgz = <G> by energy, axial surfaces
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
            stotr(nn,ng,ir,jr)=stotr(nn,ng,ir,jr)-dtgr(ng,ir,jr)
     1       *(sigt(ng,m1)-sigt(ng,m2))
c diagonal and bhs transport corrections
            if(itrcor.eq.1)then
              stotr(nn,ng,ir,jr)=stotr(nn,ng,ir,jr)+dtgr(ng,ir,jr)
     1         *(sigsp1(ng,ng,m1)-sigsp1(ng,ng,m2))
            else if(itrcor.eq.2)then
              stotr(nn,ng,ir,jr)=stotr(nn,ng,ir,jr)+dtgr(ng,ir,jr)
     1         *(sum(sigsp1(1:neg,ng,m1))-sum(sigsp1(1:neg,ng,m2)))
            end if
            if(icalc.eq.0.or.icalc.eq.3.or.icalc.eq.4)then
              sinsr(nn,ng,ir,jr)=sinsr(nn,ng,ir,jr)+dqgr(ng,ir,jr)
     1         *(nsrcf(ng,m1)-nsrcf(ng,m2))
            end if
            do lg=1,neg
              sfisr(nn,ng,ir,jr)=sfisr(nn,ng,ir,jr)+dfgr(lg,ng,ir,jr)
     1         *(chi(lg,ng,m1)*nusigf(ng,m1)
     2         -chi(lg,ng,m2)*nusigf(ng,m2))
              do l=0,isct
                sscttr(nn,lg,ir,jr)=sscttr(nn,lg,ir,jr)
     1           +dsgr(ng,lg,l,ir,jr)
     2           *(sigs(ng,lg,l,m1)-sigs(ng,lg,l,m2))
c diagonal and bhs transport corrections
                if((itrcor.eq.1.and.lg.eq.ng).or.itrcor.eq.2)then
                  sscttr(nn,lg,ir,jr)=sscttr(nn,lg,ir,jr)
     1             -dsgr(ng,lg,l,ir,jr)
     2             *(sigsp1(ng,lg,m1)-sigsp1(ng,lg,m2))
                end if
              end do ! l
            end do ! lg
          end do ! ng
        end do ! jr
      end do ! ir
      do ir=i0,nr
        do jr=1,nz
          do ng=1,neg
c the terms already have the correct sign.
c sinsr is 0 except for fixed-source problems.
            sintr(nn,ng,ir,jr)=sinsr(nn,ng,ir,jr)+sfisr(nn,ng,ir,jr)
     1       +stotr(nn,ng,ir,jr)+sscttr(nn,ng,ir,jr)
          end do ! ng
c debug print
c         write(*,'(2i4,1p5e14.6)')ir,jr,sum(sintr(nn,1:neg,ir,jr)),
c    1     sum(sinsr(nn,1:neg,ir,jr)),sum(sfisr(nn,1:neg,ir,jr)),
c    2     sum(stotr(nn,1:neg,ir,jr)),sum(sscttr(nn,1:neg,ir,jr))
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
              stotz(nn,ng,ir,jr)=stotz(nn,ng,ir,jr)-dtgz(ng,ir,jr)
     1         *(sigt(ng,m1)-sigt(ng,m2))
c diagonal and bhs transport corrections
              if(itrcor.eq.1)then
                stotz(nn,ng,ir,jr)=stotz(nn,ng,ir,jr)+dtgz(ng,ir,jr)
     1           *(sigsp1(ng,ng,m1)-sigsp1(ng,ng,m2))
              else if(itrcor.eq.2)then
                stotz(nn,ng,ir,jr)=stotz(nn,ng,ir,jr)+dtgz(ng,ir,jr)
     1           *(sum(sigsp1(1:neg,ng,m1))-sum(sigsp1(1:neg,ng,m2)))
              end if
              if(icalc.eq.0.or.icalc.eq.3.or.icalc.eq.4)then
                sinsz(nn,ng,ir,jr)=sinsz(nn,ng,ir,jr)+dqgz(ng,ir,jr)
     1           *(nsrcf(ng,m1)-nsrcf(ng,m2))
              end if
              do lg=1,neg
                sfisz(nn,ng,ir,jr)=sfisz(nn,ng,ir,jr)+dfgz(lg,ng,ir,jr)
     1           *(chi(lg,ng,m1)*nusigf(ng,m1)
     2           -chi(lg,ng,m2)*nusigf(ng,m2))
                do l=0,isct
                  sscttz(nn,lg,ir,jr)=sscttz(nn,lg,ir,jr)
     1             +dsgz(ng,lg,l,ir,jr)
     2             *(sigs(ng,lg,l,m1)-sigs(ng,lg,l,m2))
c diagonal and bhs transport corrections
                  if((itrcor.eq.1.and.lg.eq.ng).or.itrcor.eq.2)then
                    sscttz(nn,lg,ir,jr)=sscttz(nn,lg,ir,jr)
     1               -dsgz(ng,lg,l,ir,jr)
     2               *(sigsp1(ng,lg,m1)-sigsp1(ng,lg,m2))
                  end if
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
              sintz(nn,ng,ir,jr)=sinsz(nn,ng,ir,jr)+sfisz(nn,ng,ir,jr)
     1         +stotz(nn,ng,ir,jr)+sscttz(nn,ng,ir,jr)
            end do ! ng
          end do ! ir
        end do ! jr
      end if
c
c calculate the direct effect for rxn-rate ratios.
      if(nn.gt.0)then
c       if(ilnk3dnt.eq.0)then
          do k=1,nedpoints
c volume mesh indices
            ir=edpoints(k)
            jr=edpoints(nedpoints+k)
            st(1:neg)=rrxs(1:neg,irri(1,nn))/rxnratet(irri(1,nn))
     1       -rrxs(1:neg,irri(2,nn))/rxnratet(irri(2,nn))
c radial surfaces
            i1=ir-1
            i2=ir
            do ng=1,neg
              sintr(nn,ng,i1,jr)=sintr(nn,ng,i1,jr)-st(ng)
     1         *ddgr(ng,i1,jr)
              sintr(nn,ng,i2,jr)=sintr(nn,ng,i2,jr)+st(ng)
     1         *ddgr(ng,i2,jr)
            end do ! ng
c axial surfaces.
            if(calc2d.eq.1)then
              j1=jr-1
              j2=jr
              do ng=1,neg
                sintz(nn,ng,ir,j1)=sintz(nn,ng,ir,j1)-st(ng)
     1           *ddgz(ng,ir,j1)
                sintz(nn,ng,ir,j2)=sintz(nn,ng,ir,j2)+st(ng)
     1           *ddgz(ng,ir,j2)
              end do ! ng
            end if
          end do ! k
c       end if ! ilnk3dnt
      end if ! nn
c
c for leakage derivatives, divide by 1.
c for keff derivatives, divide by <psi*, F psi>/keff
c for alpha derivatives, divide by <psi*, 1/vel psi>
c for rxn rate derivatives, multiply by rxn-rate ratio
      if(nn.eq.0)then
        sintr(0,1:neg,0:nr,1:nz)=sintr(0,1:neg,0:nr,1:nz)*fissd
        sintz(0,1:neg,1:nr,0:nz)=sintz(0,1:neg,1:nr,0:nz)*fissd
      else
        do ng=1,neg
          sintr(nn,ng,0:nr,1:nz)=sintr(nn,ng,0:nr,1:nz)*
     1     rxnratio(nn)
          sintz(nn,ng,1:nr,0:nz)=sintz(nn,ng,1:nr,0:nz)*
     1     rxnratio(nn)
        end do ! ng
      end if
c
c output
      if(nn.eq.0)then
        if(icalc.eq.0)then
          write(iur,'(/," leakage derivatives",10x,"leakage=",
     1     1pe15.7,/)')lkg
        else if(icalc.eq.3.or.icalc.eq.4)then
          write(iur,'(/," count rate sensitivities",10x,
     1     "count rate=",1pe15.7)')lkg
        else if(icalc.eq.1)then
          write(iur,'(/," keff derivatives",10x,"keff=",0pf12.8,/)')
     1     keff
        else if(icalc.eq.2)then
          write(iur,'(/," alpha derivatives",10x,"alpha=",1pe15.7,/)')
     1     alpha
        end if
      else
        call outputh(iur,nn,irrx(1:2,irri(1,nn)),irrx(1:2,irri(2,nn)),
     1   icalc,rxnratio(nn))
      end if
      if(calc2d.eq.0.or.calc2d.eq.-1)then
        call output3(sintr(nn,1:neg,0:nr,1:nz))
      else if(calc2d.eq.1)then
        call output4(sintr(nn,1:neg,0:nr,1:nz),
     1   sintz(nn,1:neg,1:nr,0:nz))
      end if
c
c for feyny or sm2, write sensitivities to binary file.
      if(icalc.eq.3.or.icalc.eq.4)then
        call output3_bin(sintr(nn,1:neg,0:nr,1:nz),
     1   sintz(nn,1:neg,1:nr,0:nz))
      end if
c
      return
      end subroutine
