      subroutine integrals_2nd(is,rlip,phifor,phiadj,
     1 psifor,psiadj,phifor2,phiadj2,psifor2,psiadj2,neg0,ndir0,
     2 nitm0,njtm0,dtg,dsg,dfg)
c
c do the integrals for dtg, dfg, and dsg for 2nd-order sensitivities 
c  
c is=0 for dtg+dfg+dsg
c is=1 for dtg only
c is=2 for dfg and dtg only
c is=3 for dfg only
c is=4 for dsg and dtg only
c is=5 for dsg only
c
c phifor and phiadj are input forward and adjoint moments
c psifor and psiadj are input forward and adjoint angular fluxes
c phifor2 and phiadj2 are input 2nd-order forward and adjoint moments
c psifor2 and psiadj2 are input 2nd-order forward and adjoint angular
c         fluxes
c
      use COMS
      implicit none
      integer is,neg0,ndir0,nitm0,njtm0
      real*8 rlip(0:isct),
     1 phifor(neg,0:nmom-1,nitm,njtm),phiadj(neg,0:nmom-1,nitm,njtm),
     2 psifor(neg0,ndir0,nitm0,njtm0),psiadj(neg0,ndir0,nitm0,njtm0),     
     3 phifor2(neg,0:nmom-1,nitm,njtm),phiadj2(neg,0:nmom-1,nitm,njtm),
     4 psifor2(neg0,ndir0,nitm0,njtm0),psiadj2(neg0,ndir0,nitm0,njtm0), 
     5 dtg(neg,3,nr,nz),dsg(3,neg,neg,0:isct,nr,nz),dfg(3,neg,neg,nr,nz)
      real*8 st,afregc,afadjc,rldv
      integer ir,jr,ng,lg,i,j,k,ll,l,m,mt,n,lg1,lg2

      dtg(1:neg,1:3,1:nr,1:nz)=0.d0
      dsg(1:3,1:neg,1:neg,0:isct,1:nr,1:nz)=0.d0
      dfg(1:3,1:neg,1:neg,1:nr,1:nz)=0.d0
c
c adjoint moments (regular and generalized) are always calculational
c moments. the sign of rlip is used to convert to inner product moments.
c
      do ir=1,nr
        do jr=1,nz
          mt=mat(ir,jr)
          do i=iindex(ir-1)+1,iindex(ir)
            do j=jindex(jr-1)+1,jindex(jr)
              do ng=1,neg
                if(is.eq.0.or.is.eq.1.or.is.eq.2.or.is.eq.4)then
c estimate with moments.
                  if(iangflux.eq.0)then
                    ll=-1
                    do l=0,isct
                      rldv=rlip(l)*dv(i,j)
                      do k=0,l
                        ll=ll+1
                        dtg(ng,1,ir,jr)=dtg(ng,1,ir,jr)
     1                   +rldv*phifor(ng,ll,i,j)*phiadj(ng,ll,i,j)
                        dtg(ng,2,ir,jr)=dtg(ng,2,ir,jr)
     1                   +rldv*phifor2(ng,ll,i,j)*phiadj(ng,ll,i,j)
                        dtg(ng,3,ir,jr)=dtg(ng,3,ir,jr)
     1                   +rldv*phiadj2(ng,ll,i,j)*phifor(ng,ll,i,j)     
                        if(calc2d.eq.0.or.calc2d.eq.-1)exit
                      end do ! k
                    end do ! l
c use angular fluxes.
                  else if(iangflux.eq.1)then
                    if(calc2d.eq.0.or.calc2d.eq.-1)then
                      do m=1,isn
                        n=isn+1-m ! adjoint direction
                        afregc=(psifor(ng,m,i,j)+psifor(ng,m,i+1,j))
     1                   *0.5d0
                        afadjc=(psiadj(ng,n,i,j)+psiadj(ng,n,i+1,j))
     1                   *0.5d0
                        dtg(ng,1,ir,jr)=dtg(ng,1,ir,jr)
     1                   +wgt(m)*afregc*afadjc*dv(i,j)
                        afregc=(psifor2(ng,m,i,j)+psifor2(ng,m,i+1,j))
     1                   *0.5d0
                        afadjc=(psiadj(ng,n,i,j)+psiadj(ng,n,i+1,j))
     1                   *0.5d0
                        dtg(ng,2,ir,jr)=dtg(ng,2,ir,jr)
     1                   +wgt(m)*afregc*afadjc*dv(i,j)
                        afregc=(psifor(ng,m,i,j)+psifor(ng,m,i+1,j))
     1                   *0.5d0
                        afadjc=(psiadj2(ng,n,i,j)+psiadj2(ng,n,i+1,j))
     1                   *0.5d0
                        dtg(ng,3,ir,jr)=dtg(ng,3,ir,jr)
     1                   +wgt(m)*afregc*afadjc*dv(i,j)
                    end do ! m
                  end if
                end if
              end if
            end do ! ng
c
c dfg
c
              if(is.eq.0.or.is.eq.2.or.is.eq.3)then
                do ng=1,neg
                  do lg=1,neg
                    st=phiadj(ng,0,i,j)*dv(i,j)
                    dfg(1,ng,lg,ir,jr)=dfg(1,ng,lg,ir,jr)
     1               +st*phifor(lg,0,i,j)
c 
                    st=phiadj(ng,0,i,j)*dv(i,j)
                    dfg(2,ng,lg,ir,jr)=dfg(2,ng,lg,ir,jr)
     1               +st*phifor2(lg,0,i,j)
c 
                    st=phiadj2(ng,0,i,j)*dv(i,j)
                    dfg(3,ng,lg,ir,jr)=dfg(3,ng,lg,ir,jr)
     1               +st*phifor(lg,0,i,j)
                  end do !lg
                end do !ng
              end if
c
c dsg 
c
              if(is.eq.0.or.is.eq.4.or.is.eq.5)then
                do ng=1,neg
                  ll=-1
                  do l=0,isct
                    rldv=rlip(l)*dv(i,j)
                    do k=0,l
                      ll=ll+1
                      lg1=scgr(ng,l,mt,1)
                      lg2=scgr(ng,l,mt,2)
c 
                      st=phiadj(ng,ll,i,j)*rldv
                      do lg=lg1,lg2
                        dsg(1,ng,lg,l,ir,jr)=dsg(1,ng,lg,l,ir,jr)
     1                   +st*phifor(lg,ll,i,j)
                      end do ! lg
c
                      st=phiadj(ng,ll,i,j)*rldv
                      do lg=lg1,lg2
                        dsg(2,ng,lg,l,ir,jr)=dsg(2,ng,lg,l,ir,jr)
     1                   +st*phifor2(lg,ll,i,j)
                      end do ! lg
c
                      st=phiadj2(ng,ll,i,j)*rldv
                      do lg=lg1,lg2
                        dsg(3,ng,lg,l,ir,jr)=dsg(3,ng,lg,l,ir,jr)
     1                   +st*phifor(lg,ll,i,j)
                      end do ! lg
                      if(calc2d.eq.0.or.calc2d.eq.-1)exit
                    end do ! k
                  end do ! l
                end do ! ng
              end if
c
c dqg
c
            end do ! j
          end do ! i
        end do ! jr
      end do ! ir

      return
      end subroutine
