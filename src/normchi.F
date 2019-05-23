      subroutine normchi(nsen,izlist,schi)
c normalize chi as requested.
      use COMS
      implicit none
      integer nsen,izlist(nel)
      real*8 schi(neg,neg,nsen)
      integer i,iso,ng
      real*8 st
c ichinorm = 0/1/2 unnormalized/full norm/partial norm
c the sum of the elements of the chi matrix (fissdata=0 or 1)
c is neg, so the normalization is 1/neg, not 1.
      do i=1,nsen
        st=sum(schi(1:neg,1:neg,i))
        if(st.eq.0.d0)cycle
        if(i.le.nel)then
          iso=nm+i
        else
          iso=nm+izlist(i-nel)
        end if
        do ng=1,neg
c full normalization.
          if(fissdata.eq.2)then ! chi vector
            schi(ng,ng,i)=sum(schi(ng,1:neg,i))
     1       -st*chi(ng,ng,iso)
            if(neg.gt.1)then
              if(ng.eq.1)then
                schi(ng,ng+1:neg,i)=0.d0
              else if(ng.eq.neg)then
                schi(ng,1:ng-1,i)=0.d0
              else
                schi(ng,1:ng-1,i)=0.d0
                schi(ng,ng+1:neg,i)=0.d0
              end if
            end if
          else
            st=sum(schi(1:neg,ng,i))
            schi(1:neg,ng,i)=schi(1:neg,ng,i)
     1       -st*chi(1:neg,ng,iso)
          end if
c one more division for partial normalization.
c needs verification.
          if(ichinorm.eq.2)then
            if(fissdata.eq.2)then ! chi vector
              schi(ng,1:neg,i)=schi(ng,1:neg,i)
     1         /(1.d0-chi(ng,ng,iso))
            else
              schi(1:neg,ng,i)=schi(1:neg,ng,i)
     1         /(1.d0-chi(1:neg,ng,iso))
            end if
          end if
        end do ! ng
      end do ! i
      return
      end subroutine
