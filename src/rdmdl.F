      subroutine rdmdl(iwrite)
c read the model.
      use COMS
      implicit none
      integer iwrite
      real*8 drdzp,d0,sumw
      integer i,i0,j,j0,ia,ir,itmp,jr,k,k2,k2p,l,l1,l2,le,n,ios,mdrdz,
     1 mtr(nm),mt,m,nrad,iwe(10),nel1,nw,iwb,idum,lcm(nm),imatls_err,
     2 ierr,ipart
      character cprob(10)*80,auc*27,keyword*81,ctmp*20,matspec*20
      parameter (auc='ABCDEFGHIJKLMNOPQRSTUVWXYZ ')
      logical lex
      integer omp_get_num_threads
c
c open log file.
      open(iuo,file='sensmg.log',status='unknown',position='append')
      if(iwrite.eq.1)then
        call version(iuo)
      end if
c open input file.
      open(iui,file=ifile,status='old',iostat=ios)
      if(ios.ne.0)then
        write(*,'("error. no input file sensmg_inp for sensmg.")')
        call stoponerror
      end if
c
c allocate arrays using nr, nz, nm, nel, nrrr, nedpoints. sensmg.py
c reads them from the input file and writes them to the control file.
c they are read in rdctrl.
      call allocate_arrays_1
c
      calc2d=0 ! default is spherical
      nmxsp=1  ! read in rddantlnk_g if used
c initialize reaction-rate edits
      if(nrrr.gt.0)then
        irrr(1:4,1:nrrr)=0
        irrx(1:2,2*nrrr)=-1
      end if
c read regular input file.
      if(ilnk3dnt.eq.0)then
        read(iui,'(a)')id
c
        ierr=0
c icalc = 0/1/2/3/4 fixed-source/keff/alpha/feyny/sm2
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
     1     cycle
          if(iwb.eq.0.and.keyword(i:i).ne.' '.and.keyword(i:i).ne.'/')
     1     iwb=i
          if(iwb.ne.0.and.(keyword(i:i).eq.' '.or.keyword(i:i).eq.'/'))
     1     then
            nw=nw+1
            iwe(nw)=i-iwb
            cprob(nw)(1:iwe(nw))=keyword(iwb:i-1)
c convert to lower case
            do j=1,iwe(nw)
              if(index(auc,cprob(nw)(j:j)).ne.0)cprob(nw)(j:j)=
     1         abc(index(auc,cprob(nw)(j:j)):index(auc,cprob(nw)(j:j)))
            end do ! j
            iwb=0
          end if
          if(keyword(i:i).eq.'/')go to 20 ! ignore comments
        end do ! i
   20   do j=1,nw
          if(cprob(j)(1:iwe(j)).eq.'sph'.or.
     1     cprob(j)(1:iwe(j)).eq.'sphere')then
            calc2d=0
          else if(cprob(j)(1:iwe(j)).eq.'cyl'.or.
     1     cprob(j)(1:iwe(j)).eq.'cylinder')then
            calc2d=1
          else if(cprob(j)(1:iwe(j)).eq.'slab')then
            calc2d=-1
          else if(cprob(j)(1:iwe(j)).eq.'keff')then
            if(icalc.eq.-99)then
              icalc=1
              ilkg=1
            else
              ierr=1
            end if
          else if(cprob(j)(1:iwe(j)).eq.'lkg')then
            if(icalc.eq.-99)then
              icalc=0
            else
              ierr=1
            end if
          else if(cprob(j)(1:iwe(j)).eq.'alpha')then
            if(icalc.eq.-99)then
              icalc=2
              ilkg=1
            else
              ierr=1
            end if
          else if(cprob(j)(1:iwe(j)).eq.'feyny')then
            if(icalc.eq.-99)then
              icalc=3
              ilkg=1
            else
              ierr=1
            end if
          else if(cprob(j)(1:iwe(j)).eq.'sm2')then
            if(icalc.eq.-99)then
              icalc=4
              ilkg=1
            else
              ierr=1
            end if
          else
            ierr=2
          end if
        end do ! j
        if(ierr.eq.1)then
          write(*,'("error. use keff, alpha, lkg, feyny, or sm2 ",
     1     "keyword.")')
          write(iuo,'("error. use keff, alpha, lkg, feyny, or sm2 ",
     1     "keyword.")')
          call stoponerror
        else if(ierr.eq.2)then
          write(*,'("error. bad keyword on second line: ",/,a80)')
     1     keyword
          write(iuo,'("error. bad keyword on second line: ",/,a80)')
     1     keyword
          call stoponerror
        end if
c
        if(iwrite.eq.1)then
          write(iuo,'(a)')trim(id)
          do j=1,nw
            write(iuo,'(2x,a)')trim(cprob(j))
          end do
        end if
c
        if(icalc.eq.-99)then
          write(*,'("error. use keff, alpha, or lkg keyword on ",
     1     "second line.")')
          write(iuo,'("error. use keff, alpha, or lkg keyword on ",
     1     "second line.")')
          call stoponerror
        end if
c
        if(iangflux.eq.0)then
          if(calc2d.eq.0)then
            write(*,'("***",/,
     1       "*** warning. using flux moments for a sphere.",/,"***")')
            write(iuo,'("***",/,
     1       "*** warning. using flux moments for a sphere.",/,"***")')
          else if(calc2d.eq.-1)then
            write(*,'("***",/,
     1       "*** warning. using flux moments for a slab.",/,"***")')
            write(iuo,'("***",/,
     1       "*** warning. using flux moments for a slab.",/,"***")')
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
        if(iwrite.eq.1)write(iuo,'("number of materials =",i4)')nm
c
c control file is read again for ncb.
c ncb(i)=number of isotopes in material i
        itmp=0
        call rdctrl(1,itmp,itmp,itmp,itmp)
        ncb(0)=0
        blk=0.d0
        lcm(1)=1
        do i=2,nm
          lcm(i)=lcm(i-1)+ncb(i-1)
        end do ! i
        do i=1,nm
          read(iui,*)mtr(i),(blk(1,j),blk(2,j),j=lcm(i),lcm(i)+ncb(i)-1)
c ensure all are positive or negative.
          if(blk(2,lcm(i)).lt.0.d0)then
            ia=-1
          else if(blk(2,lcm(i)).gt.0.d0)then
            ia=1
          else
            ia=0
          end if
          do j=lcm(i)+1,lcm(i)+ncb(i)-1
            if((ia.eq.-1.and.blk(2,j).gt.0.d0).or.
     1       (ia.eq.1.and.blk(2,j).lt.0.d0))then
              write(*,'("error. material",i7," has positives and ",
     1         "negatives.")')i
              write(iuo,'("error. material",i7," has positives and ",
     1         "negatives.")')i
              call stoponerror
            else if(ia.eq.0.or.blk(2,j).eq.0.d0)then
              write(*,'("error. material",i7," has a nuclide with ",
     1         "zero mass.")')i
              write(iuo,'("error. material",i7," has a nuclide with ",
     1         "zero mass.")')i
              call stoponerror
            end if
          end do ! j
c normalize; not sure why gray calc doesn't do this.
          sumw=0.d0
          do j=lcm(i),lcm(i)+ncb(i)-1
            sumw=sumw+blk(2,j)
          end do ! j
          do j=lcm(i),lcm(i)+ncb(i)-1
            blk(2,j)=blk(2,j)/abs(sumw)
          end do ! j
        end do ! i
c
        rho(0)=0.d0
        read(iui,*)(rho(i),i=1,nm)
        do i=1,nm
          if(rho(i).eq.0.d0)then
            write(*,'("warning. density is 0 for material",i7,".")')i
            write(iuo,'("warning. density is 0 for material",i7,".")')i
          end if
        end do ! i
        if(calc2d.eq.0)then
          read(iui,*)idum ! nr read in control file
          nz=1
          r(0)=0.d0
          read(iui,*)(r(i),i=1,nr)
          read(iui,*)(mat(i,1),i=1,nr)
          mat(nr+1,1)=0
        elseif(calc2d.eq.-1)then
          read(iui,*)idum ! nr read in control file
          nz=1
          read(iui,*)(r(i),i=0,nr)
          read(iui,*)(mat(i,1),i=1,nr)
          mat(0,1)=0
          mat(nr+1,1)=0
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
        end if
c
c read reaction rates to compute.
c 2*nrrr is the max possible number of unique reaction rates.
        if((icalc.eq.3.or.icalc.eq.4).and.nrrr.gt.0)then
          write(*,'("error. reaction rates and feyny or sm2 ",
     1     "not allowed.")')
          write(iuo,'("error. reaction rates and feyny or sm2 ",
     1     "not allowed.")')
          call stoponerror
        end if
c number of edit points and number of reaction-rate ratios
c are read in control file. consistency of nedpoints and nrrr
c is checked in sensmg.py.
        read(iui,*)idum ! nedpoints
        if(iwrite.eq.1)then
          write(iuo,'("number of reaction-rate edit points =",i4)')
     1     nedpoints
        end if
c first half of edpoints array is temporarily the point indices
        if(nedpoints.gt.0)then
          read(iui,*)(edpoints(i),i=1,nedpoints)
        end if
        read(iui,*)idum ! nrrr
        if(iwrite.eq.1)then
          write(iuo,'("number of reaction-rate ratios =",i4)')nrrr
        end if
        do i=1,nrrr
          read(iui,*,iostat=ios)(irrr(j,i),j=1,4)
          if(ios.ne.0)then
            write(*,'("error reading reaction-rate ratios. (end each ",
     1       "material composition line and each ratio line with /.)")')
            call stoponerror
          end if
        end do ! i
c
c read redoin file.
      else if(ilnk3dnt.gt.0)then
        matspec=" "
        read(iui,'(a)')keyword
        read(iui,'(a)')id
        ios=0
        jt=1 ! default
        do while(ios.eq.0)
          read(iui,'(a)',iostat=ios)keyword
          le=len_trim(keyword)
          l=index(keyword,"=")
          if(index(keyword,"igeom=").ne.0)then
            if(index(keyword,"sph").ne.0)then
              calc2d=0
            else if(index(keyword,"slab").ne.0)then
              calc2d=-1
            else if(index(keyword,"r-z").ne.0)then
              calc2d=1
            else if(index(keyword,'"3"').ne.0)then
              calc2d=0
            else if(index(keyword,'"1"').ne.0)then
              calc2d=-1
            else if(index(keyword,'"7"').ne.0)then
              calc2d=1
            else
              write(*,'("error. only sph, slab, or r-z accepted.")')
              write(iuo,'("error. only sph, slab, or r-z accepted.")')
              call stoponerror
            end if
          else if(index(keyword,"it=").ne.0)then
            read(keyword(l+1:le),*)it
          else if(index(keyword,"jt=").ne.0)then
            read(keyword(l+1:le),*)jt
          else if(index(keyword,"glibname=").ne.0)then
            cycle
          else if(index(keyword,"libname=").ne.0)then
            libname=keyword(l+2:le-1)
          else if(index(keyword,"matls=").ne.0)then
            imatls_err=0
            ctmp=" "
            m=0
            nel1=0 ! nel may be reduced here for reaction rates
            do i=1,nel
              read(iui,'(a)',iostat=ios)keyword
              le=len_trim(keyword)
c first pair of quotation marks has material
              l1=index(keyword,'"')+1
              l2=index(keyword(l1:le),'"')+l1-2
c mNNNNNN is expected, maybe fewer Ns but not more.
              if(iwrite.eq.1.and.l2-l1+1.gt.7)then
                write(*,'("error. too many characters in matls ",
     1           "entry. mNNNNNN is the max.",/,a)')trim(keyword)
                write(iuo,'("error. too many characters in matls ",
     1           "entry. mNNNNNN is the max.",/,a)')trim(keyword)
                imatls_err=imatls_err+1
              end if
c if materials don't repeat, it is a new material.
              if(ctmp.ne.keyword(l1:l2))then
                ctmp=keyword(l1:l2)
                m=m+1
                read(keyword(l1+1:l2),'(i6)',iostat=ios)i0
                if(ios.ne.0)then
                  write(*,'("error reading matls. line is",/,a)')
     1             trim(keyword)
                  write(iuo,'("error reading matls. line is",/,a)')
     1             trim(keyword)
                  imatls_err=imatls_err+1
                end if
                if(iwrite.eq.1.and.i0.ne.m)then
                  write(*,'("error in ",a,". materials out of order",
     1             " or bad matls entry.")')trim(ifile)
                  write(iuo,'("error in ",a,". materials out of order",
     1             " or bad matls entry.")')trim(ifile)
                  imatls_err=imatls_err+1
                end if
                if(m.le.nm)then
                  lcm(m)=i ! beginning isotope of material m
                else
                  nel1=nel1+1 ! count isotopes only for reaction rates
                  cycle
                end if
              end if
c second pair of quotation marks has isotope
              l1=index(keyword(l2+2:le),'"')+l2+2
              l2=index(keyword(l1:le),'"')+l1-2
c maybe there is a zaid extension
              if(index(keyword(l1:l2),'.').ne.0)then
                l2=index(keyword(l1:l2),'.')+l1-2
              end if
              read(keyword(l1:l2),*)blk(1,i)
c third word is density
              l2=index(keyword(l1:le),'"')+l1-1
              l1=l2+1
              read(keyword(l1:le),*)blk(2,i)
            end do ! i
            nel=nel-nel1 ! remove materials only for reaction rates
            if(imatls_err.gt.0)call stoponerror
          else if(index(keyword,"matspec=").ne.0)then
c check if the value is on the same line
            if(index(keyword,"wtfrac").ne.0)then
              matspec="wtfrac"
            else if(index(keyword,"atfrac").ne.0)then
              matspec="atfrac"
            else if(index(keyword,"atdens").ne.0)then
              matspec="atdens"
            else
              read(iui,'(a)',iostat=ios)keyword
              le=len_trim(keyword)
c one pair of quotation marks
              l1=index(keyword,'"')+1
              l2=index(keyword(l1:le),'"')+l1-2
              matspec=keyword(l1:l2)
            end if
          else if(index(keyword,"assign=").ne.0)then
c check if the value is on the same line
            if(index(keyword,"matls").ne.0)then
              rho(1:nm)=0.d0
            else
              do i=1,nm
                read(iui,'(a)',iostat=ios)keyword
                le=len_trim(keyword)
c first pair of quotation marks has zone or "matls"
                l1=index(keyword,'"')+1
                l2=index(keyword(l1:le),'"')+l1-2
                if(keyword(l1:l2).eq."matls")then
                  rho(1:nm)=0.d0
                  exit
                end if
                ilnk3dnt=2 ! identifies mcnp-style lnk3dnt
c second pair of quotation marks has material
                l1=index(keyword(l2+2:le),'"')+l2+2
                l2=index(keyword(l1:le),'"')+l1-2
c third word is density
                l2=index(keyword(l1:le),'"')+l1-1
                l1=l2+1
                read(keyword(l1:le),*)rho(i)
              end do ! i
c diagnose when assign with materials should be assign="matls"
c for ilnk3dnt=1, nm will be set to 1 later. rho is 0.
              if(iwrite.eq.1)then
                if(abs(sum(rho(1:nm))-real(nm,8)).le.1.d-3
     1           .and.nm.eq.nel)then
                  write(*,'("****",/,"**** warning. possibly assign=",
     1             1h","matls",1h"," was intended.",/,"****")')
                  write(iuo,'("****",/,"**** warning. possibly assign=",
     1             1h","matls",1h"," was intended.",/,"****")')
                end if
              end if
            end if
          else if(index(keyword,"ievt=").ne.0)then
            read(keyword(l+1:l+3),'(i3)')icalc
            if(icalc.eq.0)then
              write(*,'("warning. ievt=0 (fixed source) in input. ",
     1         "changing to ievt=2 (alpha).")')
              write(iuo,'("warning. ievt=0 (fixed source) in input. ",
     1         "changing to ievt=2 (alpha).")')
              icalc=2
            end if
            ilkg=1
          else if(index(keyword,"edisos=").ne.0.and.nrrr.gt.0)then
            do i=1,nrrr
              irrr(1,i)=0
              irrr(3,i)=0
              read(iui,'(a)',iostat=ios)keyword
c remove quotation marks
              do j=1,len_trim(keyword)
                if(keyword(j:j).eq.'"')keyword(j:j)=" "
              end do ! j
              read(keyword,*,iostat=ios)(irrr(j,i),j=1,3,2)
            end do ! i
          else if(index(keyword,"edxs=").ne.0.and.nrrr.gt.0)then
            do i=1,nrrr
              irrr(2,i)=0
              irrr(4,i)=0
              read(iui,'(a)',iostat=ios)keyword
c remove quotation marks
              do j=1,len_trim(keyword)
                if(keyword(j:j).eq.'"')keyword(j:j)=" "
              end do ! j
c TODO allow edit names
              read(keyword,*,iostat=ios)(irrr(j,i),j=2,4,2)
            end do ! i
          else if(index(keyword,"points=").ne.0.and.nrrr.gt.0)then
c first half of edpoints array is temporarily the point indices
            read(iui,*)(edpoints(i),i=1,nedpoints)
          end if
        end do
c assume all materials have the same matspec. don't include isotopes
c in reaction rates. there should be a check for this, or do it right.
        if(ilnk3dnt.eq.2)then
          if(matspec.eq."wtfrac")then
            rho(1:nm)=-rho(1:nm)
            blk(2,1:nel)=-blk(2,1:nel)
          else
c mcnp writes mass densities in the lnk3dnt file
            rho(1:nm)=-rho(1:nm)
          end if
        end if
        if(iwrite.eq.1)then
          write(iuo,'(a)')trim(id)
          write(iuo,'("number of materials =",i4)')nm
          write(iuo,'("number of reaction-rate edit points =",i4)')
     1     nedpoints
          write(iuo,'("number of reaction-rate ratios =",i4)')nrrr
        end if
        do m=1,nm-1
          ncb(m)=lcm(m+1)-lcm(m)
        end do ! m
        ncb(nm)=nel+1-lcm(nm)
c normalize.
        if(ilnk3dnt.eq.2)then
          do m=1,nm
            sumw=0.d0
            do j=lcm(m),lcm(m)+ncb(m)-1
              sumw=sumw+blk(2,j)
            end do ! j
            do j=lcm(m),lcm(m)+ncb(m)-1
              blk(2,j)=blk(2,j)/abs(sumw)
            end do ! j
          end do ! nm
        end if
c call for geometry and nmxsp.
        call rddantlnk_g
c for "regular" redoin/lnk3dnt, ilnk3dnt=2 now. otherwise,
c there is one material. ncb=nel and all cells have that material.
        if(ilnk3dnt.eq.1)then
          nm=1
          ncb(1)=nel
          mat(1:nr,1:nz)=1
          mtr(1)=1
        end if
      end if ! lnk3dnt
      close(iui)
c
c set up reaction-rate edits
      nrrx=0
      nflux=0
      if(nrrr.gt.0)then
        call edits(iwrite,iuo,nrrr,irrr,irrx,irri,nrrx,nflux)
      end if
c
      nxs=nm+nel+nrrx
      niso=nel+nrrx
      if(iwrite.eq.1)then
        write(iuo,'(" nm,nel,nedpoints,nrrr,nrrx,nxs,niso=",7i7)')
     1   nm,nel,nedpoints,nrrr,nrrx,nxs,niso
      end if
c
c check if all materials are used. not sure how important this is.
      do i=1,nr
        do j=1,nz
          mt=mat(i,j)
          if(mt.eq.0)cycle
          do k=1,nm
            if(mt.eq.mtr(k))go to 210
          end do ! k
          write(*,'("error. material",i7," is unused.")')mtr(k)
          write(iuo,'("error. material",i7," is unused.")')mtr(k)
          call stoponerror
  210     continue
        end do ! j
      end do ! i
c
c do this lastr bit for cylinders if needed.
      if(ilkg.eq.1)then
        lastr=nr
      elseif(ilkg.eq.0)then
c       lastr=nr+2
c       d0=0.1d0
c       r(nr+1)=detdisr(1)-d0
c       r(nr+2)=detdisr(1)+d0
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
      if(iwrite.eq.1)then
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
      if(partisnver(1:1).eq."5")then
        cellsol=-1 ! cellsol is not an input for ver. 5
c new logic will be required when more partisn versions
c have this feature. this assumes minimum version 8 is 8_27_15.
      else if(partisnver(1:2).eq."8_")then
        if(calc2d.eq.1)then
          if(iplane.eq.-1.and.jplane.eq.-1)then
            cellsol=1
            iaflux=1
          end if
        end if
      end if
c
c set detector if needed.
c TODO enter this in input file.
      detname="none"
      if(icalc.eq.3.or.icalc.eq.4)then
        detname="npod"
      end if
      if(iwrite.eq.1)then
        write(iuo,'("detname is ",a,".")')trim(detname)
      end if
c
c get SENSDATA.
      call getenv('SENS_DATA',sens_data)
      l=len_trim(sens_data)
      if(l.eq.0)then
        write(*,'("error. no SENS_DATA set.")')
        write(iuo,'("error. no SENS_DATA set.")')
        call stoponerror
      end if
      if(sens_data(l:l).ne."/")sens_data(l+1:l+1)="/"
c
c ictrl=1 is to write sources4c or misc input.
c only for icalc=0, 3, or 4.
      if(ictrl.eq.1)then
c check for repeated isotope in a material and fixed-source.
c after repeated isotopes are allowed, this will not need to be
c done, and ilnk3dnt>0 will work with icalc=0.
        if(icalc.eq.0)then
          if(ilnk3dnt.eq.0)then
            k2=0
            do i=1,nm
              do j=1,ncb(i)
                k2=k2+1
                if(ncb(i).eq.1)cycle
                k2p=k2
                do k=j+1,ncb(i)
                  k2p=k2p+1
                  if(nint(blk(1,k2)).eq.nint(blk(1,k2p)))then
                    write(*,'("error. repeated isotope in a material",
     1               " does not work with fixed-source calculation.",
     2               /,"  misc and/or sources4c and sensmg have to be ",
     3               "updated.",/,"  material",i7,",",i8," repeats.")')
     4               i,nint(blk(1,k2))
                    write(iuo,'("error. repeated isotope in a material",
     1               " does not work with fixed-source calculation.",
     2               /,"  misc and/or sources4c and sensmg have to be ",
     3               "updated.",/,"  material",i7,",",i8," repeats.")')
     4               i,nint(blk(1,k2))
                    call stoponerror
                  end if
                end do ! k
              end do ! j
            end do ! i
          end if ! ilnk3dnt
        end if
        nitm=1
        njtm=1
        nmom=1
        ndir=1
c allocate arrays.
        call allocate_arrays
c
      else if(ictrl.ne.1)then
c for regular input, decide on fine mesh intervals for partisn input.
        if(ilnk3dnt.eq.0)then
          if(calc2d.eq.0.or.calc2d.eq.-1)then
c DEBUG MOMENTS (2nd 2 lines)
            drdzp=0.005d0
            mdrdz=40
c           drdzp=100.d0
c           mdrdz=4
c DEBUG feynman y (slab)
c           drdzp=0.0005d0
          else if(calc2d.eq.1)then
            drdzp=0.05d0
            mdrdz=20
c these are needed for bck01, it's large
            if(ifile(1:5).eq."bck01".and.id(1:34).eq.
     1       "brian's problem, inspired by TREAT")then
              write(*,'("warning. setting ints for bck01.")')
              write(iuo,'("warning. setting ints for bck01.")')
              drdzp=0.5d0
              mdrdz=10
            end if
          end if
c to use with tmp_sensmg_inp (2 lines)
c         drdzp=100.d0
c         mdrdz=1
c this is used for testing and studies.
c         write(*,'("drdzp=drdzp*0.1d0 in rdmdl.")')
c         write(iuo,'("drdzp=drdzp*0.1d0 in rdmdl.")')
c         drdzp=drdzp*0.1d0
c for the nse01 and berp problems:
c         if(icalc.eq.1)then
c           write(*,'("drdzp=drdzp*0.2d0 in rdmdl.")')
c           write(iuo,'("drdzp=drdzp*0.2d0 in rdmdl.")')
c           drdzp=drdzp*0.2d0
c         else if(icalc.eq.0)then
c           write(*,'("drdzp=drdzp*0.5d0 in rdmdl.")')
c           write(iuo,'("drdzp=drdzp*0.5d0 in rdmdl.")')
c           drdzp=drdzp*0.5d0
c         end if
c         write(*,'("drdzp=drdzp*2.d0 in rdmdl.")')
c         write(iuo,'("drdzp=drdzp*2.d0 in rdmdl.")')
c         drdzp=drdzp*2.d0
c         mdrdz=0
c for bckfs01
c         write(*,'("drdzp=drdzp*10.d0 in rdmdl.")')
c         write(iuo,'("drdzp=drdzp*10.d0 in rdmdl.")')
c         drdzp=drdzp*10.d0
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
            dr(i)=(r(i)-r(i-1))/real(iints(i),8)
            if(dr(i).le.0.)then
              write(*,'("error. radii are bad. r(",i2,")=",f8.4,
     1         "; r(",i2,")=",f8.4".")')i-1,r(i-1),i,r(i)
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
              dz(i)=(z(i)-z(i-1))/real(jints(i),8)
              if(dz(i).le.0.)then
                write(*,'("error. heights are bad. z(",i2,")=",f8.4,
     1           "; z(",i2,")=",f8.4".")')i-1,z(i-1),i,z(i)
                call stoponerror
              end if
              jt=jt+jints(i)
            end do ! i
            njtm=jt+1
          end if
c for lnk3dnt, it and jt have been read from redoin;
c iints and jints have been read from lnk3dnt.
        else if(ilnk3dnt.gt.0)then
          do i=1,lastr
            dr(i)=(r(i)-r(i-1))/real(iints(i),8)
          end do ! i
c define total number of fine meshes +1 for dynamic memory allocations.
          nitm=it+1
          if(calc2d.eq.0.or.calc2d.eq.-1)then
            njtm=1
          else if(calc2d.eq.1)then
            do i=1,nz
              dz(i)=(z(i)-z(i-1))/real(jints(i),8)
            end do ! i
            njtm=jt+1
          end if
        end if ! ilnk3dnt
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
c if it's a cylinder and ianglux=1, set cellsol to 0.
c for testing now, but we do want to use angular fluxes for cylinders.
        if(calc2d.eq.1.and.iangflux.eq.1)then
          cellsol=0
          iaflux=0
        end if
c       write(*,'("cellsol,iangflux,iaflux",3i4)')cellsol,iangflux,
c    1   iaflux
c       write(iuo,'("cellsol,iangflux,iaflux",3i4)')cellsol,iangflux,
c    1   iaflux
c allocate arrays.
        call allocate_arrays
      end if ! ictrl
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
      if(ictrl.ne.1)then
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
c if lnk3dnt, read idclnk3 and denlnk3. nmxsp is already read.
      if(ilnk3dnt.gt.0)then
        call rddantlnk_m
      end if
c assign idclnk3 to mat if needed. use the first fine mesh in each
c coarse mesh. afterward, everything should be transferred to regular
c variables now.
      if(ilnk3dnt.eq.2)then
        do i=1,nr
          do j=1,nz
            mat(i,j)=idclnk3(iindex(i-1)+1,jindex(j-1)+1,1)
          end do ! j
        end do ! i
c sensmg is not set up to do lnk3dnt files in which there is more
c than one fine mesh per coarse mesh. partisn writes those files,
c however. in that case, there is only one material in each fine
c mesh. look for this situation and reset ilnk3dnt.
        itmp=ilnk3dnt
        if(nmxsp.eq.1)then
          if(any(iints(1:nr).gt.1))then
            if(calc2d.eq.0.or.calc2d.eq.-1)then
              itmp=0
            else if(calc2d.eq.1)then
              if(any(jints(1:nz).gt.1))then
                itmp=0
              end if
            end if
          end if
        end if
        if(itmp.ne.ilnk3dnt)then
c check to ensure that all fine meshes in a coarse mesh have the same
c material.
          do ir=1,nr
            do jr=1,nz
              mt=mat(ir,jr)
              do i=iindex(ir-1)+1,iindex(ir)
                do j=jindex(jr-1)+1,jindex(jr)
                  if(idclnk3(i,j,1).ne.mt)then
                    write(*,'("error. fine mesh i,j=",2i7," is in ",
     1               "coarse mesh ir,jr=",2i7,".",/,"fine mesh ",
     2               "material is",i7,"; coarse mesh material is",
     3               i7,".")')i,j,ir,jr,mt,idclnk3(i,j,1)
                    write(iuo,'("error. fine mesh i,j=",2i7," is in ",
     1               "coarse mesh ir,jr=",2i7,".",/,"fine mesh ",
     2               "material is",i7,"; coarse mesh material is",
     3               i7,".")')i,j,ir,jr,mt,idclnk3(i,j,1)
                    call stoponerror
                  end if
                end do ! j
              end do ! i
            end do ! jr
          end do ! ir
          if(iwrite.eq.1)then
            write(*,'("warning. redoin/lnk3dnt are input, "
     1       "but regular partisn input will now be used.")')
            write(iuo,'("warning. redoin/lnk3dnt are input, "
     1       "but regular partisn input will now be used.")')
          end if
          ilnk3dnt=itmp
        end if
      end if
c
c reaction-rate edit points
      if(nedpoints.gt.0)then
        if(iwrite.eq.1)then
          write(iuo,'(/,"radial, axial coordinates of edit points")',
     1     advance='no')
        end if
        if(ilnk3dnt.eq.0)then
          nrad=nr
        else if(ilnk3dnt.gt.0)then
          nrad=it
        end if
        if(iwrite.eq.1)then
          if(ilnk3dnt.eq.0)then
            write(iuo,'(" (coarse meshes)")')
          else if(ilnk3dnt.gt.0)then
            write(iuo,'(" (fine meshes)")')
          end if
          write(iuo,'(4x,"index",5x,"n",5x,"r",5x,"z")')
        end if
        do k=1,nedpoints
c calculate and store radial, axial indices.
c n=(jr-1)nrad + ir
c first half of edpoints array is radial, second is axial.
          jr=edpoints(k)/nrad+1
          if(mod(edpoints(k),nrad).eq.0)jr=jr-1
          ir=edpoints(k)-(jr-1)*nrad
          if(iwrite.eq.1)then
            write(iuo,'(i9,3i6)')k,edpoints(k),ir,jr
          end if
          edpoints(k)=ir
          edpoints(nedpoints+k)=jr
        end do ! k
      end if
c
c set energy bins. needed in rdgendir if fissdata=2 or trcor=1 or 2.
      call set_ebins(neg,iug2,iuo,sens_data,libname,ebins)
c
c scattering isct+1 for diagonal or bhs transport correction
      sigsp1(1:neg,1:neg,0:nxs)=0.d0
c zaids, atomic weights, densities, spectral weight functions, and
c sigsp1 are computed once and stored in a binary file for
c quick access in iterations. since these calls will only happen
c once, write diagnostics without iwrite.
      inquire(file='sensaw',exist=lex)
      if(.not.lex)then
        if(trim(libname).eq."kynea3")then
c call with known niso,nisoaw,iht-3
          call rdbxslib(151,144,0)
        else if(trim(libname).eq."scale")then
c call with known niso,nisoaw,iht-3
c this library is only for testing. it has 5 isotopes and 6 edits.
          call rdbxslib(5,5,6)
        else if(libname(1:7).eq."special")then
          call rdspclib
        else
          call rdgendir
        end if
c debug for trcor (3 lines)
c       write(*,'("warning. setting sigsp1 to 0 in rdmdl.")')
c       write(iuo,'("warning. setting sigsp1 to 0 in rdmdl.")')
c       sigsp1(1:neg,1:neg,0:nxs)=0.d0
        call wrawbin
      else
        call rdawbin
      end if
c the redoin file has materials to be mixed, not just isotopes.
      if(ilnk3dnt.eq.2)then
        call setlnk
      end if
c
c calculate volume and mass.
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
c partisn's mass inventory
      if(iwrite.eq.1.and.ilnk3dnt.eq.1)then
        call massin
      end if
c
c find cells with fission.
      ifcel(0:nr+1,0:nz+1)=0
      n=0
      if(ilnk3dnt.eq.0)then
        do i=1,nr
          do j=1,nz
            mt=mat(i,j)
            if(mt.eq.0)cycle
            k2=sum(ncb(1:mt-1))
            do k=k2+1,k2+ncb(mt)
c             write(*,'(2i4,i8)')i,j,nint(blk(1,k))
              if(nint(blk(1,k)).eq.92233.or.nint(blk(1,k)).eq.92235
     1         .or.nint(blk(1,k)).eq.93237.or.nint(blk(1,k)).eq.94239
     2         .or.nint(blk(1,k)).eq.92238)ifcel(i,j)=1
            end do ! k
            n=n+ifcel(i,j)
          end do ! j
        end do ! i
      else if(ilnk3dnt.gt.0)then
        do ir=1,nr
          do jr=1,nz
            do i=iindex(ir-1)+1,iindex(ir)
              do j=jindex(jr-1)+1,jindex(jr)
                do k2=1,nmxsp
                  if(idclnk3(i,j,k2).eq.0)exit
                  k=nint(blk(1,idclnk3(i,j,k2)))
c                 write(*,'(2i4,i8)')i,j,k
                  if(k.eq.92233.or.k.eq.92235
     1             .or.k.eq.93237.or.k.eq.94239
     2             .or.k.eq.92238)ifcel(ir,jr)=1
                end do ! k2
                n=n+ifcel(ir,jr)
              end do ! j
            end do ! i
          end do ! jr
        end do ! ir
      end if
      if(n.eq.0)then
        write(*,'("warning. no fissionable material.",/,
     1   "  need 92233, 92235, 93237, 94239, or 92238.")')
        write(iuo,'("warning. no fissionable material.",/,
     1   "  need 92233, 92235, 93237, 94239, or 92238.")')
      end if
c
c put here to check existence of efficiency file. might as well
c load deteff as well. (sens_data is needed and deteff has to be
c allocated first.)
      call set_deteff(iwrite,neg,iug2,iuo,sens_data,detname,deteff,pi)
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
c use stream file?
c only for partisn version 8_31 or later
      lstream=.false.
      if(icalc.eq.2)then
        ipart=0
        l1=index(partisnver,"_")-1
        if(l1.ge.1)then
          read(partisnver(1:l1),*)ipart
        end if
        if(ipart.gt.8)then
          lstream=.true.
        else
          l1=l1+2
          l2=index(partisnver(l1:len_trim(partisnver)),"_")+l1-2
          read(partisnver(l1:l2),*)ipart
          if(ipart.ge.31)then
            lstream=.true.
          end if
        end if
      end if ! icalc
c it isn't needed
      lstream=.false.
c
      return
      end subroutine
c
      subroutine edits(iwrite,iuo,nrrr,irrr,irrx,irri,nrrx,nflux)
      implicit none
      integer iwrite,iuo,nrrr,irrr(4,nrrr),irrx(2,2*nrrr),irri(2,nrrr),
     1 nrrx,nflux
      integer i,j,n1
      if(iwrite.eq.1)then
        write(iuo,'("reaction-rate ratios",/,4x,"index",6x,"numerator",
     1   9x,"denominator")')
      end if
      do i=1,nrrr
        if(iwrite.eq.1)write(iuo,'(5i9)')i,irrr(1:4,i)
        if(irrr(1,i).eq.irrr(3,i).and.irrr(2,i).eq.irrr(4,i))then
          write(*,'("error. equivalent numerator and denominator ",
     1     "not allowed:",/,5i9)')irrr(1:4,i)
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
              go to 10
            end if
          end do ! j
          nrrx=nrrx+1
          irrx(1:2,nrrx)=irrr(1:2,i)
   10     do j=1,n1
            if(irrr(3,i).eq.irrx(1,j).and.irrr(4,i).eq.irrx(2,j))then
              go to 20
            end if
          end do ! j
          nrrx=nrrx+1
          irrx(1:2,nrrx)=irrr(3:4,i)
   20     continue
        end if
c if there is no denominator, use the flux. here nflux takes
c the index that the flux has; later, the flux will be moved
c to the end of the list.
c (for fixed-source, will not divide by flux; see rxnrateedit.)
        if(irrx(1,nrrx).eq.0.and.irrx(2,nrrx).eq.0)then
          nflux=nrrx
        end if
      end do ! i
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
      if(iwrite.eq.1)then
        if(nrrx.gt.0)write(iuo,'("unique reaction rates")')
        do i=1,nrrx+nflux
          write(iuo,'(3i9)')i,irrx(1:2,i)
        end do ! i
        if(nrrr.gt.0)write(iuo,'("indices of ratios")')
        do i=1,nrrr
          write(iuo,'(i9,i3,"/",i3)')i,irri(1,i),irri(2,i)
        end do ! i
      end if
      return
      end subroutine
