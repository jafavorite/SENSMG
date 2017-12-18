      subroutine calcsens_d(n,ediv,fissd,dtg,dsg,dfg,dqg,rfile,
     1 rxnrate,rxnratet)
c compute sensitivities to densities.
      use COMS
      implicit none
      integer n
      real*8 ediv(neg),fissd,dtg(neg,nr,nz),dsg(neg,neg,0:isct,nr,nz),
     1 dfg(neg,neg,nr,nz),dqg(neg,nr,nz),
     2 rxnrate(neg,nrrx),rxnratet(nrrx+nflux)
      character rfile*20
      real*8 sden(0:nrrr,neg,nm),stot(0:nrrr,neg,nm),
     1 sfis(0:nrrr,neg,nm),sins(0:nrrr,neg,nm),ssctt(0:nrrr,neg,nm)
      integer i,ir,jr,l,lg,ng
c
      sden(0:nrrr,1:neg,1:nm)=0.d0
      stot(0:nrrr,1:neg,1:nm)=0.d0
      sfis(0:nrrr,1:neg,1:nm)=0.d0
      sins(0:nrrr,1:neg,1:nm)=0.d0
      ssctt(0:nrrr,1:neg,1:nm)=0.d0
c
c compute and write sensitivities.
c n=0 is for keff, alpha, or leakage sensitivities,
c n > 0 is for reaction-rate sensitivities.
      if(n.le.1)then
        write(*,'(/,"writing sensitivities to file ",a,".")')
     1   trim(rfile)
        write(iuo,'(/,"writing sensitivities to file ",a,".")')
     1   trim(rfile)
      end if
c
c this logic catches a material used in different regions.
c dtg = <G*, psi> by energy, coarse mesh
c dsg = <G*, (2l+1)psi> by energy out, energy in, moment, coarse mesh
c dfg = <G*, psi/keff> by energy out, energy in, coarse mesh
c dqg = <G*> by energy, coarse mesh
      do i=1,nm
        do ir=1,nr
          do jr=1,nz
            if(mat(ir,jr).ne.i)cycle
            do ng=1,neg
              stot(n,ng,i)=stot(n,ng,i)-dtg(ng,ir,jr)*sigt(ng,i)
              if(icalc.eq.0)then
                sins(n,ng,i)=sins(n,ng,i)+dqg(ng,ir,jr)*nsrcf(ng,i)
              end if
              do lg=1,neg
                sfis(n,ng,i)=sfis(n,ng,i)+dfg(lg,ng,ir,jr)
     1           *chi(lg,ng,i)*nusigf(ng,i)
                do l=0,isct
                  ssctt(n,lg,i)=ssctt(n,lg,i)+dsg(ng,lg,l,ir,jr)
     1             *sigs(ng,lg,l,i)
                end do ! l
              end do ! lg
            end do ! ng
          end do ! jr
        end do ! ir
      end do ! i
      do i=1,nm
        do ng=1,neg
c the terms already have the correct sign.
c sins is 0 except for fixed-source problems.
          sden(n,ng,i)=sins(n,ng,i)+sfis(n,ng,i)+stot(n,ng,i)
     1     +ssctt(n,ng,i)
        end do ! ng
      end do ! i
c
c for leakage sensitivities, divide by leakage
c for keff sensitivities, divide by <psi*, F psi>
c for alpha sensitivities, divide by <psi*, 1/vel psi>*alpha
      if(n.eq.0)then
        sden(0,1:neg,1:nm)=sden(0,1:neg,1:nm)*fissd
      end if
c
c output
      if(n.eq.0)then
        if(icalc.eq.0)then
          write(iur,'(/," leakage sensitivities",10x,"leakage=",
     1     1pe15.7,/)')lkg
        else if(icalc.eq.1)then
          write(iur,'(/," keff sensitivities",10x,"keff=",0pf12.8,/)')
     1     keff
        else if(icalc.eq.2)then
          write(iur,'(/," alpha sensitivities",10x,"alpha=",1pe15.7,/)')
     1     alpha
        end if
      else
        call outputh(iur,n,irrx(1:2,irri(1,n)),irrx(1:2,irri(2,n)),
     1   icalc,rxnratio(n))
      end if
      call output2(sden(n,1:neg,1:nm))
c
      return
      end subroutine
