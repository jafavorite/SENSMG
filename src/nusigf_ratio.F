      subroutine nusigf_ratio(iso,k2,nfdenom,nfratio)
c returns the ratio of the group-sum of the isotopic 
c macroscopic nusigf to the group-sum of the material
c macroscopic nusigf.
      use COMS
      implicit none
      integer iso,k2
      real*8 nfdenom,nfratio
      integer ng
c
      nfratio=0.d0
      if(nfdenom.ne.0.d0)then
        do ng=1,neg
          nfratio=nfratio+nusigf(ng,iso)*wgtfct(ng,k2)
        end do ! ng
        nfratio=blk(3,k2)*nfratio/nfdenom
      end if
      return
      end subroutine

      subroutine nusigf_denom(i0,nfdenom)
c returns the group-sum of the material macroscopic nusigf multiplied
c by the spectrum weight function.
      use COMS
      implicit none
      integer i0
      real*8 nfdenom
      integer iso,j,k2,ng
c
      nfdenom=0.d0
      k2=sum(ncb(1:i0-1))
      do j=1,ncb(i0)
        k2=k2+1
        iso=nm+k2
        do ng=1,neg
          nfdenom=nfdenom+blk(3,k2)*nusigf(ng,iso)*wgtfct(ng,k2)
        end do ! ng
      end do ! j
      return
      end subroutine

      subroutine nusigf_denom_lnk3(ir,jr,nfdenom)
c returns the group-sum of the material macroscopic nusigf multiplied
c by the spectrum weight function.
      use COMS
      implicit none
      integer ir,jr
      real*8 nfdenom
      integer iso,k2,n,ng
c
      nfdenom=0.d0
      do n=1,nmxsp
        k2=idclnk3(ir,jr,n)
        if(k2.eq.0)exit
        iso=nm+k2
        do ng=1,neg
          nfdenom=nfdenom+denlnk3(ir,jr,n)*nusigf(ng,iso)*wgtfct(ng,k2)
        end do ! ng
      end do ! n
      return
      end subroutine
