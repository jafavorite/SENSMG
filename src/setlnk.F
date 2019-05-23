      subroutine setlnk
c set up lnk3dnt variables for lnk3dnt files from partisn and mcnp.
c partisn writes 1. for densities (when a sensmg for_inp file is used),
c so use mass densities on assign keyword.
c mcnp writes mass densities, so use 1. for densities on assign
c keyword.
      use COMS
      implicit none
      integer ir,jr,k2,m,mt,n,nmxsp_new,nmxsp0,nx
      real(R8KIND), allocatable, dimension(:,:,:) :: denlnk3_tmp
      integer, allocatable, dimension(:,:,:)      :: idclnk3_tmp
c
      allocate(idclnk3_tmp(it,jt,nmxsp))
      allocate(denlnk3_tmp(it,jt,nmxsp))
      idclnk3_tmp(1:it,1:jt,1:nmxsp)=idclnk3(1:it,1:jt,1:nmxsp)
      denlnk3_tmp(1:it,1:jt,1:nmxsp)=denlnk3(1:it,1:jt,1:nmxsp)
c
c find new nmxsp (max isotopes in a mesh)
c     write(*,'("setlnk. nmxsp=",i6)')nmxsp
      nmxsp_new=0
      do ir=1,it
        do jr=1,jt
          if(idclnk3(ir,jr,1).eq.0)cycle
          nmxsp0=0
          do m=1,nmxsp
            mt=idclnk3(ir,jr,m)
            if(mt.eq.0)exit
            k2=sum(ncb(1:mt-1))
            nmxsp0=nmxsp0+ncb(mt)
c debug print
c           do n=k2+1,k2+ncb(mt)
c             write(*,'(5i4,i10,1p2e34.26)')ir,jr,m,mt,n,nint(blk(1,n)),
c    1         denlnk3(ir,jr,m),blk(2,n)
c           end do ! n
          end do ! m
          nmxsp_new=max(nmxsp_new,nmxsp0)
        end do ! jr
      end do ! ir
c     write(*,'("setlnk. nmxsp_new=",i6)')nmxsp_new
      deallocate(idclnk3)
      deallocate(denlnk3)
      allocate(idclnk3(it,jt,nmxsp_new))
      allocate(denlnk3(it,jt,nmxsp_new))
c
      idclnk3(1:it,1:jt,1:nmxsp_new)=0
      denlnk3(1:it,1:jt,1:nmxsp_new)=0.d0
      do ir=1,it
        do jr=1,jt
          if(idclnk3_tmp(ir,jr,1).eq.0)cycle
          nx=0
          do m=1,nmxsp
            mt=idclnk3_tmp(ir,jr,m)
            if(mt.eq.0)exit
            k2=sum(ncb(1:mt-1))
            do n=k2+1,k2+ncb(mt)
              nx=nx+1
              idclnk3(ir,jr,nx)=n
c convert mass density to atom density
              denlnk3(ir,jr,nx)=denlnk3_tmp(ir,jr,m)*blk(2,n)
     1         *rho(mt)*avg/atwt(n)
            end do ! n
          end do ! m
        end do ! jr
      end do ! ir
c
      deallocate(idclnk3_tmp)
      deallocate(denlnk3_tmp)
      nmxsp=nmxsp_new
c
c debug print
c     do ir=1,it
c       do jr=1,jt
c         if(idclnk3(ir,jr,1).eq.0)cycle
c         do m=1,nmxsp
c           mt=idclnk3(ir,jr,m)
c           if(mt.eq.0)exit
c           write(*,'(3i4,i10,1p3e34.26)')ir,jr,m,nint(blk(1,mt)),
c    1       denlnk3(ir,jr,m),blk(2,mt),atwt(mt)
c         end do ! m
c       end do ! jr
c     end do ! ir
c     call stoponerror
c
      return
      end subroutine
