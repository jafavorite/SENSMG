      function plgndr(l,m,x)
c  adapted almost line for line from Press et al., "Numerical Recipes in
c  FORTRAN: The Art of Scientific Computing," 2nd. Ed., Chapter 6, p. 247.
      implicit none
      integer l,m
      real*8 plgndr,x
c computes the associated legendre polynomial p_l^m(x). m and l are
c integers satisfying 0<=m<=l, x lies in the range -1<=x<=1.
      integer i,ll
      real*8 fact,pll,pmm,pmmp1,somx2
      if(m.lt.0.or.m.gt.l.or.abs(x).gt.1.)print *, ! was pause
     1 'bad arguments in plgndr'
      pmm=1.d0
      if(m.gt.0)then
        somx2=sqrt((1.d0-x)*(1.d0+x))
        fact=1.d0
        do i=1,m
          pmm=-pmm*fact*somx2
          fact=fact+2.d0
        end do ! i
      end if
      if(l.eq.m)then
        plgndr=pmm
      else
        pmmp1=x*(2*m+1)*pmm
        if(l.eq.m+1)then
          plgndr=pmmp1
        else
          do ll=m+2,l
            pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
            pmm=pmmp1
            pmmp1=pll
          end do ! ll
          plgndr=pll
        end if
      end if
      return
      end function
