      subroutine wrsources
c write sources4c input files
      use COMS
      implicit none
      integer maxtarg
      parameter (maxtarg=23)
      real*8 suma(nel)
      integer i,j,k,l,m,ir,i1,ic,ne,ng,ielem(nel),nas,iant(maxtarg),
     1 jant(maxtarg),rra
      character sfile*20,aline(0:maxtarg)*80,rrc*80
c (alpha,n) target isotopes. 060000 and 060120 added to
c catch c-13 in elemental carbon or c-12.
c see sources4c data file tape3.
      data iant/
     1 030070, 040090, 050100, 050110, 050000, 060130,
     2 070140, 080170, 080180, 090190, 100210, 100220,
     3 110230, 120000, 120250, 120260, 130270, 140290,
     4 140300, 150310, 170370, 
     5 060000, 060120/
c we do not use the sources4c (alpha,n) sources or spont. fis. sources.
c we just list all the isotopes as sources in the input. we commented
c the error in sources4c when a source isotope is not actually a source
c (replaced stop with cycle).
c
c blk(3) is atoms/b-cm
      ic=1
      do i=1,nm
        if(i.eq.1)then
          l=1
        else
          l=l+ncb(i-1)
        end if
        ic=ic+ncb(i)
        ne=0 ! number of elements in i
        suma(1:nel)=0.d0
        do j=1,ncb(i)
c check to see if this element has appeared.
          do k=1,ne
            if(int(blk(1,l+j-1))/1000.eq.ielem(k))go to 20
          end do ! k
c a new element, search for others like it in i.
          ne=ne+1
          ielem(ne)=int(blk(1,l+j-1))/1000
          suma(ne)=blk(3,l+j-1)
          do k=j+1,ncb(i)
            if(int(blk(1,l+k-1))/1000.eq.ielem(ne))then
              suma(ne)=suma(ne)+blk(3,l+k-1)
            end if
          end do ! k
   20     continue
        end do ! j
c       do j=1,ne
c         write(*,'(2i4,1pe15.7)')i,ielem(j),suma(j)/rhoa(i)
c       end do ! j
c
c write sources4c input file.
        write(sfile,'("sources/m",i3.3,"_tape1")')i
        open(iun,file=sfile,status='unknown')
        write(iun,'(a," mat ",i2)')trim(id),i
        write(iun,'("1 2 -1")')
        write(iun,'(i3," 0")')ne
        do j=1,ne
          write(iun,'(5x,i3.3,1pe15.7)')ielem(j),suma(j)/rhoa(i)
        end do ! j
        write(iun,'(i4,1p2e12.5)')-neg,ebins(2,1),ebins(1,neg)
        write(iun,'(4x,1p5e12.5)')(ebins(2,ng),ng=1,neg)
        write(iun,'(i3)')ncb(i)
        do j=1,ncb(i)
          write(iun,'(4x,i7.7,1pe15.7)')int(blk(1,l+j-1))*10,
     1     blk(3,l+j-1)*1.d24
        end do ! j
c find (alpha,n) source isotopes
        nas=0
        do j=1,ncb(i)
c chlorine is the largest-z target
          if(int(blk(1,l+j-1))/1000.gt.17)cycle
          do k=1,maxtarg
            if(int(blk(1,l+j-1))*10.eq.iant(k))go to 100
          end do ! k
          cycle
  100     nas=nas+1
          jant(nas)=iant(k)
          if(jant(nas).eq.120000)then
c use 12025 and 12026 instead of 12000 because 12000 isn't in tape4.
            write(iuo,'("comment. using Mg-25 and Mg-26 instead of Mg",
     1       " in sources4c; material",i4,".")')i
            jant(nas)=120250
            write(aline(nas),'(4x,i7.7,1pe15.7)')jant(nas),
     1       blk(3,l+j-1)/rhoa(i)*0.1000d0
            nas=nas+1
            jant(nas)=120260
            write(aline(nas),'(4x,i7.7,1pe15.7)')jant(nas),
     1       blk(3,l+j-1)/rhoa(i)*0.1101d0
          else if(k.le.maxtarg-2)then
            write(aline(nas),'(4x,i7.7,1pe15.7)')int(blk(1,l+j-1))*10,
     1       blk(3,l+j-1)/rhoa(i)
          else
c c-13 in c or c-12. mult by atom fraction.
            write(iuo,'("comment. using C-13 instead of C or C-12",
     1       " in sources4c; material",i4,".")')i
            write(aline(nas),'(4x,i7.7,1pe15.7)')060130,
     1       blk(3,l+j-1)/rhoa(i)*0.0110d0
          end if
        end do ! j
        if(nas.eq.0)then ! need a dummy?
          write(aline(nas),'(4x,i7.7,1pe15.7)')0040090,1.d-24
        else
c sort the lines. seems to be a sources4c requirement for targets.
c heapsort, numerical recipes, p. 329.
          if(nas.eq.1)go to 130
          m=nas/2+1
          ir=nas
  110     continue
          if(m.gt.1)then
            m=m-1
            rra=jant(m)
            rrc=aline(m)
          else
            rra=jant(ir)
            rrc=aline(ir)
            jant(ir)=jant(1)
            aline(ir)=aline(1)
            ir=ir-1
            if(ir.eq.1)then
              jant(1)=rra
              aline(1)=rrc
              go to 130
            endif
          endif
          i1=m
          j=m+m
  120     if(j.le.ir)then
            if(j.lt.ir)then
              if(jant(j).lt.jant(j+1))j=j+1
            endif
            if(rra.lt.jant(j))then
              jant(i1)=jant(j)
              aline(i1)=aline(j)
              i1=j
              j=j+j
            else
              j=ir+1
            endif
            go to 120
          endif
          jant(i1)=rra
          aline(i1)=rrc
          go to 110
  130     continue
        end if
c nag = number of alpha-particle groups
        write(iun,'(i3,1x,i4)')nas,nag
        do j=1,nas
          write(iun,'(a)')trim(aline(j))
        end do ! j
        close(iun)
      end do ! i
c write mcnp style e bins
c     write(*,'(7x,1p5e12.5)')0.d0,ebins(1,neg),(ebins(2,ng),
c    1 ng=neg,1,-1)
c
c     do i=1,nm
c       l=lcm(i)
c       do j=1,ncb(i)
c         write(*,'(3i8,1pe15.7)')i,j,int(blk(1,l+j-1)),
c    1     blk(3,l+j-1)/rhoa(i)
c         write(*,'(3i8,1pe15.7)')i,j,int(blk(1,l+j-1)),
c    1     blk(3,l+j-1)*1.d24
c       end do ! j
c       write(*,'(i8,1pe15.7)')i,rhoa(i)
c     end do ! i
c
      return
      end subroutine
