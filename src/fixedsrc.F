      subroutine fixedsrc(nf,ith)
c write fixed source to partisn input file.
      use COMS
      implicit none
      integer nf ! file number
      integer ith ! =0/1 forward/adjoint
      integer i,j,ng,n1
c
      if(ith.eq.0)then
c forward source.
        write(nf,'("/ total source rate density in neuts/s/cm^3 ",
     1   "by region:")')
        do j=1,nz
          write(nf,'("/",1p5e13.5)')(nsrc(mat(i,j)),i=1,nr)
        end do ! j
        write(nf,'("/ total source rate in neuts/s by region:")')
        do j=1,nz
          write(nf,'("/",1p5e13.5)')(nsrc(mat(i,j))*vol(i,j),i=1,nr)
        end do ! j
        write(nf,'(" sourcf=")')
        if(calc2d.eq.0.or.calc2d.eq.-1)then
          do ng=1,lng
            do i=1,nr,3
              n1=min(nr,i+2)
              write(nf,'(1x,3(i6,"r",1pe16.9))')(iints(j),
     1         nsrcf(ng,mat(j,1)),j=i,n1)
            end do ! i
            write(nf,'(7x,";")')
          end do ! ng
c write mcnp style sources
c         do i=1,nr
c           write(*,'("region",i6," mat",i6)')i,mat(i,1) 
c           write(*,'("total source in neuts/s:",1pe13.5)')
c    1       nsrc(mat(i,1))*vol(i,1)
c           write(*,'(7x,1p5e12.5)')0.d0,0.d0,(nsrcf(ng,mat(i,1)),
c    1       ng=neg,1,-1)
c         end do ! i
c         call stoponerror
#if ( defined RPSD18 )
c write source for rpsd2018
          write(iuo,'("neutron source for debugging.")')
          do ng=1,lng
            write(iuo,'(i4,1pe16.9)')ng,nsrcf(ng,1)
          end do ! ng
c         call stoponerror
#endif
        elseif(calc2d.eq.1)then
          do ng=1,lng
            write(nf,'("/ group: ",i3)')ng
            if(ilkg.eq.1)then
              do j=1,nz
                do i=1,nr
                  write(nf,'(("   ",i6,"r",1pe16.9))')iints(i),
     1             nsrcf(ng,mat(i,j))
                end do ! i
                write(nf,'("   ;",i6,"y",1x,"1;")')jints(j)-1
              end do ! j
            elseif(ilkg.eq.0)then
c do something here if this is used.
            end if
          end do ! ng
        end if
      else if(ith.eq.1)then
c adjoint source for leakage.
c right leakage for slab and sphere, total for cylinder.
        if(calc2d.eq.0.or.calc2d.eq.-1)then
          write(nf,'(" sirite=")')
          do ng=1,neg,3
            n1=min(neg,ng+2)
            write(nf,'(1x,1p3e16.9)')(deteff(j),j=ng,n1)
          end do ! ng
        else if(calc2d.eq.1)then
          write(nf,'(" sirite=",i4,"r 1.;",i6,"y 1;")')neg,jt-1
          write(nf,'(" sitop= ",i4,"r 1.;",i6,"y 1;")')neg,it-1
          write(nf,'(" sibott=",i4,"r 1.;",i6,"y 1;")')neg,it-1
        end if
      end if
      return
      end subroutine
