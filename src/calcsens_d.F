      subroutine calcsens_d(nn,ediv,fissd,dtg,dsg,dfg,dqg,rfile,
     1 rxnrate,rxnratet,nm1)
c compute sensitivities to material densities.
      use COMS
      implicit none
      integer nn,nm1
      real*8 ediv(neg),fissd,dtg(neg,nr,nz),dsg(neg,neg,0:isct,nr,nz),
     1 dfg(neg,neg,nr,nz),dqg(neg,nr,nz),
     2 rxnrate(neg,nrrx),rxnratet(nrrx+nflux)
      character rfile*20
      real*8 sden(neg,nm1),stot(neg,nm1),sfis(neg,nm1),sins(neg,nm1),
     1 ssctt(neg,nm1),tau_lkg
      integer i,ir,iso,j,jr,k2,l,lg,ng,m1,m2,nm2
c
c if ilnk3dnt=0, nm1=nm; if ilnk3dnt>0, nm1=nr*nz (all meshes)
      sden(1:neg,1:nm1)=0.d0
      stot(1:neg,1:nm1)=0.d0
      sfis(1:neg,1:nm1)=0.d0
      sins(1:neg,1:nm1)=0.d0
      ssctt(1:neg,1:nm1)=0.d0
c
c compute and write sensitivities.
c nn=0 is for keff, alpha, or leakage sensitivities,
c nn > 0 is for reaction-rate sensitivities.
      if(nn.le.1)then
        write(*,'(/,"writing sensitivities to file ",a,".")')
     1   trim(rfile)
        write(iuo,'(/,"writing sensitivities to file ",a,".")')
     1   trim(rfile)
      end if
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
              stot(ng,m2)=stot(ng,m2)-dtg(ng,ir,jr)*sigt(ng,m1)
c diagonal and bhs transport corrections
              if(itrcor.eq.1)then
                stot(ng,m2)=stot(ng,m2)+dtg(ng,ir,jr)*sigsp1(ng,ng,m1)
              else if(itrcor.eq.2)then
                stot(ng,m2)=stot(ng,m2)+dtg(ng,ir,jr)
     1           *sum(sigsp1(1:neg,ng,m1))
              end if
              if(icalc.eq.0.or.icalc.eq.3.or.icalc.eq.4)then
                sins(ng,m2)=sins(ng,m2)+dqg(ng,ir,jr)*nsrcf(ng,m1)
              end if
              do lg=1,neg
                sfis(ng,m2)=sfis(ng,m2)+dfg(lg,ng,ir,jr)
     1           *chi(lg,ng,m1)*nusigf(ng,m1)
                do l=0,isct
                  ssctt(lg,m2)=ssctt(lg,m2)+dsg(ng,lg,l,ir,jr)
     1             *sigs(ng,lg,l,m1)
c diagonal and bhs transport corrections
                  if((itrcor.eq.1.and.lg.eq.ng).or.itrcor.eq.2)then
                    ssctt(lg,m2)=ssctt(lg,m2)-dsg(ng,lg,l,ir,jr)
     1               *sigsp1(ng,lg,m1)
                  end if
                end do ! l
              end do ! lg
            end do ! ng
          end do ! jr
        end do ! ir
      end do ! i
      do i=1,nm1
        do ng=1,neg
c the terms already have the correct sign.
c sins is 0 except for fixed-source problems.
          sden(ng,i)=sins(ng,i)+sfis(ng,i)+stot(ng,i)+ssctt(ng,i)
        end do ! ng
      end do ! i
c
c for leakage sensitivities, divide by leakage
c for keff sensitivities, divide by <psi*, F psi>
c for alpha sensitivities, divide by <psi*, 1/vel psi>*alpha
      if(nn.eq.0)then
        sden(1:neg,1:nm1)=sden(1:neg,1:nm1)*fissd
      end if
c
c for alpha eigenvalue, compute 1/tau for leakage
      if(icalc.eq.2.and.nn.eq.0)then
        tau_lkg=alpha*(-1.d0+sum(sden(1:neg,1:nm1)))
        write(*,'(/,"1/tau,leakage=",1pe17.9)')tau_lkg
        write(iuo,'(/,"1/tau,leakage=",1pe17.9)')tau_lkg
      end if
c
c output
      if(nn.eq.0)then
        if(icalc.eq.0)then
          write(iur,'(/," leakage sensitivities",10x,"leakage=",
     1     1pe15.7,/)')lkg
        else if(icalc.eq.3.or.icalc.eq.4)then
          write(iur,'(/," count rate sensitivities",10x,
     1     "count rate=",1pe15.7)')lkg
        else if(icalc.eq.1)then
          write(iur,'(/," keff sensitivities",10x,"keff=",0pf12.8,/)')
     1     keff
        else if(icalc.eq.2)then
          write(iur,'(/," alpha sensitivities",10x,"alpha=",1pe15.7,/)')
     1     alpha
        end if
      else
        call outputh(iur,nn,irrx(1:2,irri(1,nn)),irrx(1:2,irri(2,nn)),
     1   icalc,rxnratio(nn))
      end if
      call output2(nm1,sden(1:neg,1:nm1))
c for feyny or sm2, write sensitivities to binary file.
      if(icalc.eq.3.or.icalc.eq.4)then
        call output2_bin(nm1,sden(1:neg,1:nm1))
      end if
c
      return
      end subroutine
