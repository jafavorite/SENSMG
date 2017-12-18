      subroutine momcon(isct1,nmom1,mcon,di)
c compute constants for angular moments in cylinders.
      use COMS
      implicit none
      integer isct1,nmom1
      real*8 mcon(nmom1,ndir),di(ndir)
      real*8 cosine(0:isct1,1:ndir),fact(0:isct1*2)
      integer i,j,k,ll,l,m,n4
      real*8 plgndr ! function to return legendre polynomials
      real*8 wplg2d(nmom1,ndir)
c
c calculate cosines for each ordinate.
      n4=ndir/4
      do m=1,n4
        di(m)=-dir(m)
        eta(m)=-eta(m)
      end do ! m
      do m=n4+1,2*n4
        di(m)=dir(m)
        eta(m)=-eta(m)
      end do ! m
      do m=2*n4+1,3*n4
        di(m)=-dir(m)
      end do ! m
      do m=3*n4+1,4*n4
        di(m)=dir(m)
      end do ! m
      cosine(0,1:ndir)=1.d0
      if(isct1.gt.0)then
        cosine(1,1:ndir)=eta(1:ndir)/sqrt(1.d0-di(1:ndir)**2)
      end if
c cos(k phi) recursion is Eq. (24) from
c <http://mathworld.wolfram.com/TrigonometricAdditionFormulas.html>.
      do k=2,isct1
        cosine(k,1:ndir)=2.d0*cosine(k-1,1:ndir)*cosine(1,1:ndir)
     1   -cosine(k-2,1:ndir)
      end do ! l
c
c make a table of factorials.
      fact(0)=1.d0
      do i=1,isct1*2
        fact(i)=i
        do j=i-1,1,-1
          fact(i)=fact(i)*j
        end do ! j
      end do ! i
c
c make table of weight * legendre polynomial and set up constants.
c mcon is the constant that goes in the equation for the moments,
c Eq. (8-17b) of partisn manual.  in the l,k sums, when k=0, mcon
c and cosine are not used [partisn manual, Eq. (8-16)].
      ll=0
      do l=0,isct1
        do k=0,l
          ll=ll+1
          do m=1,ndir
            wplg2d(ll,m)=wgt(m)*plgndr(l,k,di(m))
          end do ! m
          if(k.eq.0)then
            mcon(ll,1:ndir)=1.d0
          else
            mcon(ll,1:ndir)=sqrt(2.d0*fact(l-k)/fact(l+k))
     1       *cosine(k,1:ndir)
          end if
        end do ! k
      end do ! l
      mcon(1:nmom1,1:ndir)=mcon(1:nmom1,1:ndir)*wplg2d(1:nmom1,1:ndir)
      return
      end subroutine
