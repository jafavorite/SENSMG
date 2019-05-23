      subroutine wrdantlnk
c write a lnk3dnt file, converting a regular coarse-mesh, fine-mesh
c input to an all-coarse-mesh input.
c called from sensmg.
c run the problem. copy the redoin file to the working directory.
c in the redoin file, remove block 2 after "block_2" and block 2
c arrays after "block_2_arrays". in block 1, set mt and nzone to nel.
c in block 1, set im=it, jm=jt. in block 4, change matspec to "atdens".
c list each isotope separately with density 1. change assign to
c "matls".
      use COMS
      implicit none
      character hname*8,huse(2)*8
      integer i,j,ir,jr,k2,i1,j1,mt,n,nwd
      parameter (nwd=27)
      integer*4 ivers,ia(nwd)
      character bfile*11
c
      if(ilnk3dnt.gt.0)then
        write(*,'("error. wrdantlnk can only be called for ",
     1   "sensmg input, not redoin/lnk3dnt.")')
        call stoponerror
      end if
      bfile="tmp_lnk3dnt"
      write(*,'("writing ",a,"...")')trim(bfile)
      write(iuo,'("writing ",a,"...")')trim(bfile)
      open(iul,file=bfile,status='unknown',position='rewind',
     1 form='binary')
c
      hname="lnk3dnt"
      huse(1)="090538"
      huse(2)="011620"
      ivers=4
      write(iul,err=910)hname,(huse(i),i=1,2),ivers
      if(calc2d.eq.-1)then
        ia(1)=1
      else if(calc2d.eq.0)then
        ia(1)=3
      else if(calc2d.eq.1)then
        ia(1)=7
      end if
      ia(2)=niso ! nzone?
      ia(3:4)=0 ! unused
      ia(5)=it
      ia(6)=jt
      ia(7)=1
      ia(8)=it
      ia(9)=jt
      ia(10)=1
      ia(11:22)=0 ! unused
      ia(23)=0 ! itemp
      nmxsp=maxval(ncb(1:nm))
      ia(24)=nmxsp
      ia(25:26)=0 ! unused
      ia(27)=0 ! ilevel
      write(iul)(ia(i),i=1,nwd)
coarse mesh boundaries, number of equally spaced fine meshes
      if(ia(1).eq.1.or.ia(1).eq.3)then
        write(iul,err=910)(rfm(i),i=0,it),(1,i=1,it)
      else if(ia(7).eq.1)then
        write(iul,err=910)(rfm(i),i=0,it),(zfm(i),i=0,jt),(1,i=1,it),
     1   (1,i=1,jt)
      end if
      write(*,'("ivers=",i1)')ivers
      write(*,'(9i8)')(ia(i),i=1,nwd)
c fine mesh interval widths if ivers.eq.5
      if(ivers.eq.5)then
c TODO
      end if
c mixing arrays
c macroscopic material number
      deallocate(idclnk3)
      deallocate(denlnk3)
      allocate(idclnk3(ia(8),ia(9),ia(24)))
      allocate(denlnk3(ia(8),ia(9),ia(24)))
      idclnk3(1:ia(8),1:ia(9),1:ia(24))=0
      denlnk3(1:ia(8),1:ia(9),1:ia(24))=0.d0
      do ir=1,nr
        do jr=1,nz
          mt=mat(ir,jr)
          if(mt.eq.0)then
            cycle
          else
            k2=sum(ncb(1:mt-1))
          end if
          do n=1,ncb(mt)
            k2=k2+1
            do i1=iindex(ir-1)+1,iindex(ir)
              do j1=jindex(jr-1)+1,jindex(jr)
                idclnk3(i1,j1,n)=k2
                denlnk3(i1,j1,n)=blk(3,k2)
              end do ! j1
            end do ! i1
          end do ! k
        end do ! jr
      end do ! jz
c     do i=1,it
c       do j=1,jt
c         do n=1,nmxsp
c           if(idclnk3(i,j,n).gt.0)then
c             write(*,'(3i6,i8,1pe14.6)')i,j,n,
c    1         nint(blk(1,idclnk3(i,j,n))),denlnk3(i,j,n)
c           end if
c         end do ! n
c       end do ! j
c     end do ! i
c 1 means void so increment by 1 to recover partisn lnk3dnt.
      idclnk3(1:ia(8),1:ia(9),1:ia(24))=idclnk3(1:ia(8),1:ia(9),
     1 1:ia(24))+1
      write(iul,err=910)(((idclnk3(i,j,n),i=1,ia(8)),j=1,ia(9)),
     1 n=1,ia(24))
c material density
      write(iul,err=910)(((denlnk3(i,j,n),i=1,ia(8)),j=1,ia(9)),
     1 n=1,ia(24))
      write(*,'("done.")')
      close(iul)
      go to 999
c
  900 write(*,'("error opening ",a,".")')trim(bfile)
      write(iuo,'("error opening ",a,".")')trim(bfile)
      call stoponerror
  910 write(*,'("error writing ",a,".")')trim(bfile)
      write(iuo,'("error writing ",a,".")')trim(bfile)
      call stoponerror
c
  999 return
      end subroutine
