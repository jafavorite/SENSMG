      subroutine rdmdl(ncbc)
c read the model.
      use COMS
      implicit none
      character ncbc*120
      real*8 drdzp,d0,sumw
      integer i,i0,j,j0,ia,ir,jr,k,n,n1,ios,mdrdz,mtr(nm),mt,
     1 iwe(10),nw,iwb,idum,lcm(nm),lcb(nm)
      logical lwrite
      character cprob(10)*80,auc*27,keyword*81
      parameter (auc='ABCDEFGHIJKLMNOPQRSTUVWXYZ ')
      integer omp_get_num_threads
c
c open log file.
      open(iuo,file='sensmg.log',status='unknown',position='append')
c open input file.
      open(iui,file=ifile,status='old',iostat=ios)
      if(ios.ne.0)then
        write(*,'("error. no input file sensmg_inp for sensmg.")')
        call stoponerror
      end if
c
c allocate arrays using nr, nz, nm, nel, nrrr
      call allocate_arrays_1
c
      read(iui,'(a)')id
c
      calc2d=0 ! default is spherical
c icalc = 0/1/2 fixed-source/keff/alpha
      icalc=-99
      keyword=" "
      read(iui,'(a81)')keyword
      keyword(81:81)=" "
      cprob(1:10)=" "
      nw=0
      iwb=0
      do i=1,81
        if(iwb.eq.0.and.keyword(i:i).eq.' ')cycle
        if(iwb.ne.0.and.keyword(i:i).ne.' '.and.keyword(i:i).ne.'/')
     1   cycle
        if(iwb.eq.0.and.keyword(i:i).ne.' '.and.keyword(i:i).ne.'/')
     1   iwb=i
        if(iwb.ne.0.and.(keyword(i:i).eq.' '.or.keyword(i:i).eq.'/'))
     1   then
          nw=nw+1
          iwe(nw)=i-iwb
          cprob(nw)(1:iwe(nw))=keyword(iwb:i-1)
c convert to lower case
          do j=1,iwe(nw)
            if(index(auc,cprob(nw)(j:j)).ne.0)cprob(nw)(j:j)=
     1       abc(index(auc,cprob(nw)(j:j)):index(auc,cprob(nw)(j:j)))
          end do ! j
          iwb=0
        end if
        if(keyword(i:i).eq.'/')go to 20 ! ignore comments
      end do ! i
   20 do j=1,nw
        if(cprob(j)(1:iwe(j)).eq.'sph'.or.
     1   cprob(j)(1:iwe(j)).eq.'sphere')then
          calc2d=0
        else if(cprob(j)(1:iwe(j)).eq.'cyl'.or.
     1   cprob(j)(1:iwe(j)).eq.'cylinder')then
          calc2d=1
        else if(cprob(j)(1:iwe(j)).eq.'slab')then
          calc2d=-1
        else if(cprob(j)(1:iwe(j)).eq.'keff')then
          if(icalc.eq.-99)then
            icalc=1
            ilkg=1
          else
            call version(iuo)
            write(*,'("error. use keff, alpha, or lkg keyword.")')
            write(iuo,'("error. use keff, alpha, or lkg keyword.")')
            call stoponerror
          end if
        else if(cprob(j)(1:iwe(j)).eq.'lkg')then
          if(icalc.eq.-99)then
            icalc=0
          else
            call version(iuo)
            write(*,'("error. use keff, alpha, or lkg keyword.")')
            write(iuo,'("error. use keff, alpha, or lkg keyword.")')
            call stoponerror
          end if
        else if(cprob(j)(1:iwe(j)).eq.'alpha')then
          if(icalc.eq.-99)then
            icalc=2
          else
            call version(iuo)
            write(*,'("error. use keff, alpha, or lkg keyword.")')
            write(iuo,'("error. use keff, alpha, or lkg keyword.")')
            call stoponerror
          end if
        else
          call version(iuo)
          write(*,'("error. bad keyword on second line: ",/,a80)')
     1     keyword
          write(iuo,'("error. bad keyword on second line: ",/,a80)')
     1     keyword
          call stoponerror
        end if
      end do ! j
c
c decide whether to write problem specs to output file.
      lwrite=.false.
      if(iter.lt.0.and.((icalc.eq.0.and.ictrl.eq.1).or.
     1 (icalc.gt.0.and.ictrl.eq.2)))then
        lwrite=.true.
      end if
      if(lwrite)then
        call version(iuo)
        write(iuo,'(a)')trim(id)
        do j=1,nw
          write(iuo,'(2x,a)')trim(cprob(j))
        end do
      end if
c
      if(icalc.eq.-99)then
        write(*,'("error. use keff, alpha, or lkg keyword on ",
     1   "second line.")')
        write(iuo,'("error. use keff, alpha, or lkg keyword on ",
     1   "second line.")')
        call stoponerror
      end if
c
      if(iangflux.eq.0)then
        if(calc2d.eq.0)then
          write(*,'("***",/,
     1     "*** warning. using flux moments for a sphere.",/,"***")')
          write(iuo,'("***",/,
     1     "*** warning. using flux moments for a sphere.",/,"***")')
        else if(calc2d.eq.-1)then
          write(*,'("***",/,
     1     "*** warning. using flux moments for a slab.",/,"***")')
          write(iuo,'("***",/,
     1     "*** warning. using flux moments for a slab.",/,"***")')
        end if
      end if
c
      read(iui,*,iostat=ios)libname
c
      read(iui,*,iostat=ios)idum ! nm read in control file
c on error, do something, but what? then return.
      if(ios.ne.0)then
        return
      end if
      if(lwrite)write(iuo,'("number of materials =",i4)')nm
c
c ncbc is a string read in rdctrcl.
c ncb(i)=number of isotopes in material i
      read(ncbc,*,iostat=ios)(ncb(i),i=1,nm)
      if(ios.ne.0)then
        write(*,'("error reading ncbc. increase its length in ",
     1   "rdctrl.F, rdmdl.F, and sensmg.F.")')
        call stoponerror
      end if
      blk=0.d0
      lcm(1)=1
      do i=2,nm
        lcm(i)=lcm(i-1)+ncb(i-1)
      end do ! i
      do i=1,nm
        read(iui,*)mtr(i),(blk(1,j),blk(2,j),j=lcm(i),lcm(i)+ncb(i)-1) 
c ensure all are positive or negative. preparation for allowing
c atom fractions; not fully implemented.
        if(blk(2,lcm(i)).lt.0.d0)then
          ia=-1
        else if(blk(2,lcm(i)).gt.0.d0)then
          ia=1
        else
          ia=0
        end if
        do j=lcm(i)+1,lcm(i)+ncb(i)-1
          if((ia.eq.-1.and.blk(2,j).gt.0.d0).or.
     1     (ia.eq.1.and.blk(2,j).lt.0.d0))then
            write(*,'("error. material",i5," has positives and ",
     1       "negatives.")')i
            write(iuo,'("error. material",i5," has positives and ",
     1       "negatives.")')i
            call stoponerror
          else if(ia.eq.0.or.blk(2,j).eq.0.d0)then
            write(*,'("error. material",i5," has a nuclide with ",
     1       "zero mass.")')i
            write(iuo,'("error. material",i5," has a nuclide with ",
     1       "zero mass.")')i
            call stoponerror
          end if
        end do ! j
c normalize; not sure why gray calc doesn't do this.
        sumw=0.d0
        do j=lcm(i),lcm(i)+ncb(i)-1
          sumw=sumw+abs(blk(2,j))
        end do ! j
        do j=lcm(i),lcm(i)+ncb(i)-1
          blk(2,j)=blk(2,j)/sumw
        end do ! j
      end do ! i
c
      read(iui,*)(rho(i),i=1,nm)
      if(calc2d.eq.0)then
        read(iui,*)idum ! nr read in control file
        nz=1
        r(0)=0.d0
        read(iui,*)(r(i),i=1,nr)
        read(iui,*)(mat(i,1),i=1,nr)
        mat(nr+1,1)=0
        read(iui,*)irrcmr
        irrcmz=1
      elseif(calc2d.eq.-1)then
        read(iui,*)idum ! nr read in control file
        nz=1
        read(iui,*)(r(i),i=0,nr)
        read(iui,*)(mat(i,1),i=1,nr)
        mat(0,1)=0
        mat(nr+1,1)=0
        read(iui,*)irrcmr
        irrcmz=1
      elseif(calc2d.eq.1)then
        read(iui,*)idum,idum ! nr,nz read in control file
        r(0)=0.d0
        read(iui,*)(r(i),i=1,nr)
        read(iui,*)(z(j),j=0,nz)
        do j=nz,1,-1
          read(iui,*,iostat=ios)(mat(i,j),i=1,nr)
          mat(nr+1,j)=0
        end do
        mat(1:nr+1,0)=0
        mat(1:nr+1,nz+1)=0
        read(iui,*)irrcmr,irrcmz
      end if
c
c read reaction rates to compute.
c 2*nrrr is the max possible number of unique reaction rates.
      nrrx=0
      nflux=0
      if(nrrr.gt.0)irrx(1:2,2*nrrr)=-1
      read(iui,*)idum ! nrrr ! number of reaction-rate ratios read in control file
      if(lwrite)then
        write(iuo,'("number of reaction-rate ratios =",i4)')nrrr
        if(nrrr.gt.0)write(iuo,'("ratios",/,4x,"index",
     1   6x,"numerator",9x,"denominator")')
      end if
      do i=1,nrrr
        irrr(1:4,i)=0
        read(iui,*,iostat=ios)(irrr(j,i),j=1,4)
        if(ios.ne.0)then
          write(*,'("error reading reaction-rate ratios. (end each ",
     1     "material composition line and each ratio line with /.)")')
          call stoponerror
        end if
        if(lwrite)write(iuo,'(5i9)')i,irrr(1:4,i)
        if(irrr(1,i).eq.irrr(3,i).and.irrr(2,i).eq.irrr(4,i))then
          write(*,'("error. equivalent numerator and denominator not ",
     1     "allowed:",/,5i9)')irrr(1:4,i)
          call stoponerror
        end if
        if(i.eq.1)then
          irrx(1:2,1)=irrr(1:2,i)
          irrx(1:2,2)=irrr(3:4,i)
          nrrx=2
        else
          n1=nrrx
          do j=1,n1
            if(irrr(1,i).eq.irrx(1,j).and.irrr(2,i).eq.irrx(2,j))then
              go to 80
            end if
          end do ! j
          nrrx=nrrx+1
          irrx(1:2,nrrx)=irrr(1:2,i)
   80     do j=1,n1
            if(irrr(3,i).eq.irrx(1,j).and.irrr(4,i).eq.irrx(2,j))then
              go to 81
            end if
          end do ! j
          nrrx=nrrx+1
          irrx(1:2,nrrx)=irrr(3:4,i)
   81     continue
        end if
c if there is no denominator, use the flux. here nflux takes
c the index that the flux has; later, the flux will be moved
c to the end of the list.
c (for fixed-source, will not divide by flux; see rxnrateedit.)
        if(irrx(1,nrrx).eq.0.and.irrx(2,nrrx).eq.0)then
          nflux=nrrx
        end if
      end do ! i
      close(iui)
c move flux "reaction rate" to the end of the list.
c reset nflux to 1 to indicate that it exists.
      if(nflux.ne.0)then
        nrrx=nrrx-1
        do i=nflux,nrrx
          irrx(1:2,i)=irrx(1:2,i+1)
        end do ! i
        nflux=1
        irrx(1:2,nrrx+nflux)=0
      end if
c match reaction-rate index with numerator and denominator
      do i=1,nrrr
        do j=1,nrrx+nflux
          if(irrr(1,i).eq.irrx(1,j).and.irrr(2,i).eq.irrx(2,j))then
            irri(1,i)=j ! numerator
          end if
          if(irrr(3,i).eq.irrx(1,j).and.irrr(4,i).eq.irrx(2,j))then
            irri(2,i)=j ! denominator
          end if
        end do ! j
      end do ! i
      if(lwrite)then
        if(nrrx.gt.0)write(iuo,'("unique reaction rates")')
        do i=1,nrrx+nflux
          write(iuo,'(3i9)')i,irrx(1:2,i)
        end do ! i
        if(nrrr.gt.0)write(iuo,'("indices of ratios")')
        do i=1,nrrr
          write(iuo,'(i9,i3,"/",i3)')i,irri(1,i),irri(2,i)
        end do ! i
      end if
      nxs=nm+nel+nrrx
      niso=nel+nrrx
      if(lwrite)then
        write(iuo,'(" nm,nel,nrrr,nrrx,nxs,niso=",6i4)')nm,nel,
     1   nrrr,nrrx,nxs,niso
      end if
c
c lcm(k)=beginning isotope in material k
c lcb(mt)=beginning isotope in material mt
      do i=1,nr
        do j=1,nz
          mt=mat(i,j)
          if(mt.eq.0)cycle
          do k=1,nm
            if(mt.eq.mtr(k))go to 210
          end do ! k
          go to 880 ! error trap
  210     lcb(mt)=lcm(k)
        end do ! j
      end do ! i
c
c calculate volume and mass.
c do this lastr bit for cylinders if needed.
      if(ilkg.eq.1)then
        lastr=nr
      elseif(ilkg.eq.0)then
c       lastr=nr+2
c       d0=0.1d0
c       r(nr+1)=detdisr(1)-d0
c       r(nr+2)=detdisr(1)+d0
      end if
      if(calc2d.eq.0)then
        do i=1,lastr
          vol(i,1)=pi43*(r(i)**3-r(i-1)**3)
          mass(i,1)=rho(mat(i,1))*vol(i,1)
        end do ! i
      else if(calc2d.eq.-1)then
        do i=1,lastr
          vol(i,1)=r(i)-r(i-1)
          mass(i,1)=rho(mat(i,1))*vol(i,1)
        end do ! i
      else if(calc2d.eq.1)then
        do i=1,lastr
          do j=1,nz
            vol(i,j)=pi*(r(i)**2-r(i-1)**2)*(z(j)-z(j-1))
            mass(i,j)=rho(mat(i,j))*vol(i,j)
          end do ! j
        end do ! i
      end if
c
c find cells with fission.
      ifcel(0:nr+1,0:nz+1)=0
      n=0
      do i=1,nr
        do j=1,nz
          mt=mat(i,j)
          if(mt.eq.0)cycle
          do k=lcb(mt),lcb(mt)+ncb(mt)-1
c           write(*,'(2i4,i8)')i,j,nint(blk(1,k))
            if(nint(blk(1,k)).eq.92233.or.nint(blk(1,k)).eq.92235
     1       .or.nint(blk(1,k)).eq.93237.or.nint(blk(1,k)).eq.94239
     2       .or.nint(blk(1,k)).eq.92238)ifcel(i,j)=1
          end do ! k
          n=n+ifcel(i,j)
        end do ! j
      end do ! i
      if(n.eq.0)then
        write(*,'("error. no fissionable material.",/,
     1   "  need 92233, 92235, 93237, 94239, or 92238.")')
        call stoponerror
      end if
c
      if(libname(1:6).eq."vitb6 ")then
        lng=199
      else
        lng=neg ! no gamma rays
      end if
c
c get the number of threads
!$OMP PARALLEL
      num_threads=omp_get_num_threads()
!$OMP END PARALLEL
      if(lwrite)then
        write(*,'("number of threads for inner products=",i4)')
     1   num_threads
        write(iuo,'("number of threads for inner products=",i4)')
     1   num_threads
      end if
      if(num_threads.ne.1)then
        write(*,'("  error. num_threads must be 1.",/,
     1   "  setenv OMP_NUM_THREADS 1")')
        write(iuo,'("  error. num_threads must be 1.",/,
     1   "  setenv OMP_NUM_THREADS 1")')
        call stoponerror
      end if
c
c set cellsol and iaflux
      cellsol=0
      iaflux=0
      if(partisn(1:1).eq."5")then
        cellsol=-1 ! cellsol is not an input for ver. 5
c new logic will be required when more partisn versions
c have this feature.
      else if(partisn(1:4).eq."8_27")then
        if(calc2d.eq.1)then
          if(iplane.eq.-1.and.jplane.eq.-1)then
            cellsol=1
            iaflux=1
          end if
        end if
      end if
c
c ictrl=1 is to write sources4c or misc input.
c only for icalc=0.
      if(ictrl.eq.1)then
        nitm=1
        njtm=1
        nmom=1
        ndir=1
c allocate arrays.
        call allocate_arrays
c
c zaids and atomic weights
        if(libname.eq."kynea3")then
c call with known niso,nisoaw,iht-3
          call rdbxslib(151,144,0)
        else
          call rdgendir
        end if
c
      else if(ictrl.ne.1)then
c decide on fine mesh intervals for partisn input.
        if(calc2d.eq.0.or.calc2d.eq.-1)then
          drdzp=0.005d0
          mdrdz=40
        else if(calc2d.eq.1)then
          drdzp=0.05d0
          mdrdz=20
c these are needed for bck01, it's large
          if(ifile(1:5).eq."bck01".and.id(1:34).eq.
     1     "brian's problem, inspired by TREAT")then
            write(*,'("warning. setting ints for bck01.")')
            write(iuo,'("warning. setting ints for bck01.")')
            drdzp=0.5d0
            mdrdz=10
          end if
        end if
c this is used for testing and studies.
c       write(*,'("drdzp=drdzp*0.1d0 in rdmdl.")')
c       write(iuo,'("drdzp=drdzp*0.1d0 in rdmdl.")')
c       drdzp=drdzp*0.1d0
c for the nse01 and berp problems:
c       if(icalc.eq.1)then
c         write(*,'("drdzp=drdzp*0.2d0 in rdmdl.")')
c         write(iuo,'("drdzp=drdzp*0.2d0 in rdmdl.")')
c         drdzp=drdzp*0.2d0
c       else if(icalc.eq.0)then
c         write(*,'("drdzp=drdzp*0.5d0 in rdmdl.")')
c         write(iuo,'("drdzp=drdzp*0.5d0 in rdmdl.")')
c         drdzp=drdzp*0.5d0
c       end if
c       write(*,'("drdzp=drdzp*2.d0 in rdmdl.")')
c       write(iuo,'("drdzp=drdzp*2.d0 in rdmdl.")')
c       drdzp=drdzp*2.d0
c       mdrdz=0
c for bckfs01
c       write(*,'("drdzp=drdzp*10.d0 in rdmdl.")')
c       write(iuo,'("drdzp=drdzp*10.d0 in rdmdl.")')
c       drdzp=drdzp*10.d0
c
        do i=1,nr
          iints(i)=max(mdrdz,nint((r(i)-r(i-1))/drdzp))
        end do ! i 
        if(ilkg.eq.0)then
          do i=nr+1,nr+2
            drdzp=0.05d0
            mdrdz=1
            iints(i)=max(mdrdz,nint((r(i)-r(i-1))/drdzp))
          end do ! i 
        endif
        if(calc2d.eq.0.or.calc2d.eq.-1)then
          jints(1)=1
        else if(calc2d.eq.1)then
          do i=1,nz
            jints(i)=max(mdrdz,nint((z(i)-z(i-1))/drdzp))
          end do ! i 
        end if
c
        it=0
        do i=1,lastr
          dr(i)=(r(i)-r(i-1))/real(iints(i))
          if(dr(i).le.0.)then
            write(*,'("error. radii are bad. r(",i2,")=",f8.4,
     1       "; r(",i2,")=",f8.4".")')i-1,r(i-1),i,r(i)
            call stoponerror
          end if
          it=it+iints(i)
        end do ! i
c define total number of fine meshes +1 for dynamic memory allocations.
        nitm=it+1
        if(calc2d.eq.0.or.calc2d.eq.-1)then
          jt=1
          njtm=1
        else if(calc2d.eq.1)then
          jt=0
          do i=1,nz
            dz(i)=(z(i)-z(i-1))/real(jints(i))
            if(dz(i).le.0.)then
              write(*,'("error. heights are bad. z(",i2,")=",f8.4,
     1         "; z(",i2,")=",f8.4".")')i-1,z(i-1),i,z(i)
              call stoponerror
            end if
            jt=jt+jints(i)
          end do ! i
          njtm=jt+1
        end if
c
c number of angles and moments
        if(calc2d.eq.0.or.calc2d.eq.-1)then
          ndir=isn
          nmom=isct+1
        elseif(calc2d.eq.1)then
c         ndir=isn*(isn+2)/2 ! for iquad=4 and 5
          ndir=isn**2        ! for iquad=6
c number of moments for two-d geometries (partisn manual, table 8.2).
          nmom=(isct+1)*(isct+2)/2 ! for iquad=5 and 6
c         nmom=isn*(isn+2)/2       ! for iquad=4
        end if
c
c use angular fluxes for cylinders (if commented)
        if(calc2d.eq.1)iangflux=0
c allocate arrays.
        call allocate_arrays
c
c zaids and atomic weights
        if(libname.eq."kynea3")then
c call with known niso,nisoaw,iht-3
          call rdbxslib(151,144,0)
        else
          call rdgendir
        end if
        if(ictrl.eq.2)then
          write(iuo,'(2x,"material atom densities")')
          j=1
          do i=1,nm
            write(iuo,'(i5,1pe17.9)')i,sum(blk(3,j:j-1+ncb(i)))
            j=j+ncb(i)
          end do ! nm
        end if
c
        iindex(0)=0
        jindex(0)=0
        do i=1,nr
          iindex(i)=iindex(i-1)+iints(i)
        end do ! i
        do j=1,nz
          jindex(j)=jindex(j-1)+jints(j)
        end do ! j
c
c calculate fine mesh surfaces and volumes
        rfm(0)=r(0)
        do ir=1,lastr
          do i=iindex(ir-1)+1,iindex(ir)
            i0=i-iindex(ir-1)
            rfm(i)=r(ir-1)+i0*dr(ir)
            if(calc2d.eq.0)then
              dv(i,1)=pi43*(rfm(i)**3-rfm(i-1)**3)
            elseif(calc2d.eq.-1)then
              dv(i,1)=rfm(i)-rfm(i-1)
            end if
            if(isn.eq.1)dv(i,1)=rfm(i)-rfm(i-1)
c           write(*,'(i6,1p2e14.6)')i,rfm(i),dv(i,1)
          end do ! i
        end do ! ir
c
        if(calc2d.eq.1)then
          zfm(0)=z(0)
          do jr=1,nz
            do j=jindex(jr-1)+1,jindex(jr)
              j0=j-jindex(jr-1)
              zfm(j)=z(jr-1)+j0*dz(jr)
            end do ! j
          end do ! jr
          do i=1,it
            do j=1,jt
              dv(i,j)=pi*(rfm(i)**2-rfm(i-1)**2)*(zfm(j)-zfm(j-1))
            end do ! j
          end do ! i
        end if
      end if ! ictrl
c
c initialize spontaneous fission source rate density
c ismat(0) is for mat=0 (void)
      ismat(0:nm)=0
      sfiso(1:neg,1:nel)=0.d0
      saniso(1:neg,1:nel)=0.d0
      nsrcf(1:neg,0:nm)=0.d0
      nsrc(0:nm)=0.d0
c isan identifies (alpha,n) source materials
      isan(0:nm)=0
c
      return
c
  880 write(*,'("error. material missing.")')
      call stoponerror
      end subroutine
