      subroutine rdxsmg
c read partisn output files for multigroup cross sections
      use COMS
      implicit none
      character line*130
      integer ng,mg,n,k,l,nu,i2lp1,ifissdata,n1,n2
      real*8 chitot(2)
c temperature only used in cross section write
      real*8 temps
#if defined ( SU )
      character bfile*7
      integer lg,ln
#endif
c
      maxup=-neg
c recall nxs=nm+nel+nrrx (this is the order of the cross sections).
c nm is the number of materials.
c nel is the number if isotopes in the materials.
c nrrx is the number of unique isotope/reaction combinations.
c zero the cross sections up to nm. load defaults into scgr.
      chi(1:neg,1:neg,0:nm)=0.d0
      nusigf(1:neg,0:nm)=0.d0
      sigt(1:neg,0:nm)=0.d0
      siga(1:neg,0:nm)=0.d0
      scgr(1:neg,0:isct,1:nm,1)=1
      scgr(1:neg,0:isct,1:nm,2)=neg
      scgr(1:neg,0:isct,0,   1)=neg
      scgr(1:neg,0:isct,0,   2)=1
      sigs(1:neg,1:neg,0:isct,0:nm)=0.d0
c read for_macrxs for material cross-sections
      call rdmacrxs(1,temps)
      do n=1,nm
        nusigf(1:neg,n)=nusigf(1:neg,n)*rho(n)
        sigt(1:neg,n)=sigt(1:neg,n)*rho(n)
        siga(1:neg,n)=siga(1:neg,n)*rho(n)
        sigs(1:neg,1:neg,0:isct,n)=sigs(1:neg,1:neg,0:isct,n)*rho(n)
      end do ! n
c write some cross sections for debugging.
c     do ng=1,neg
c       write(*,'(i6,1p3e18.10)')ng,sigt(ng,1),nusigf(ng,1),siga(ng,1)
c     end do ! ng
c     call stoponerror
c
c zero the cross sections from nm+1 to the end. load defaults into scgr.
      chi(1:neg,1:neg,nm+1:nxs)=0.d0
      nusigf(1:neg,nm+1:nxs)=0.d0
      sigt(1:neg,nm+1:nxs)=0.d0
      siga(1:neg,nm+1:nxs)=0.d0
      scgr(1:neg,0:isct,nm+1:nxs,1)=1
      scgr(1:neg,0:isct,nm+1:nxs,2)=neg
      sigs(1:neg,1:neg,0:isct,nm+1:nxs)=0.d0
c read xs1_macrxs for element and isotope cross-sections
c sigf and sigc are read in rdxsmg0/snxedt0
      call rdmacrxs(2,temps)
c write some cross sections for debugging.
c     nu=0
c     do n=nm+1,nxs
c       nu=nu+1
c       write(*,'(i6,1p3e18.10)')nu,sigt(467,n)*blk(3,nu),
c    1   nusigf(467,n)*blk(3,nu),siga(467,n)*blk(3,nu)
c     end do ! ng
c     call stoponerror
c
      nu=0 ! count groups with upscattering.
      do ng=1,neg
        do n=1,nxs
          if(scgr(ng,0,n,2).gt.ng)then
            nu=nu+1
            exit
          end if
        end do ! n
      end do ! ng
      if(maxup.eq.-neg)maxup=nu ! only if cross sections not read from macrxs
      if(nu.gt.0)then
        write(*,'(5x,"cross sections have p0 upscattering in",i4,
     1   " groups; maxup=",i3,".")')nu,maxup
        write(iuo,'(5x,"cross sections have p0 upscattering in",i4,
     1   " groups; maxup=",i3,".")')nu,maxup
      else
        write(*,'(5x,"no p0 upscattering in cross sections; maxup=",
     1   i3,".")')maxup
        write(iuo,'(5x,"no p0 upscattering in cross sections; maxup=",
     1   i3,".")')maxup
      end if
c
c check for zero velocities.
      do ng=1,neg
        if(vel(ng).eq.0.d0)then
          vel(ng)=1.d0
        endif
      end do ! ng
c
c check chi matrix normalization
      if(iter.eq.0.and.ictrl.eq.3)then
        write(iuo,'("chi matrix normalization (fissdata=",i1,")",
     1   /," mat/iso ng sum-out-to-ng  sum-ng-to-out")')fissdata
        do n=1,nm+nel
          if(sum(chi(1:neg,1:neg,n)).eq.0.d0)cycle
          chitot(1:2)=0.d0
          do ng=1,neg
            if(fissdata.ne.2)then
              write(iuo,'(2i5,1p2e15.7)')n,ng,sum(chi(ng,1:neg,n)),
     1         sum(chi(1:neg,ng,n))
              chitot(1)=chitot(1)+sum(chi(ng,1:neg,n))
              chitot(2)=chitot(2)+sum(chi(1:neg,ng,n))
            else
c for fissdata=2, the chi vector is stored on the diagonal of the 
c chi matrix; see rdmacrxs.
              write(iuo,'(2i5,1p2e15.7)')n,ng,chi(ng,ng,n),chi(ng,ng,n)
              chitot(1)=chitot(1)+chi(ng,ng,n)
              chitot(2)=chitot(2)+chi(ng,ng,n)
            end if
          end do ! ng
          write(iuo,'(i5,"  tot",1p2e15.7)')n,chitot(1:2)
        end do ! n
      end if
c
c write cross sections for debugging
c     do ng=1,neg
c carbon p0 self-scatter
c       write(100,'(i6,1pe15.7)')ng,sigs(ng,ng,0,nm+3)
c carbon total
c       write(100,'(i6,1pe15.7)')ng,sigt(ng,nm+3)
c material total
c       write(100,'(i6,1pe15.7)')ng,sigt(ng,1)
c     end do ! ng
c     close(100)
c     write(*,'("done writing debug xsecs.")')
c     call stoponerror
c
c check if the cross sections were read and stored properly.
      go to 990
      ifissdata=1
      i2lp1=0
c for materials
      n1=1
      n2=nm
c for isotopes
c     n1=nm+1
c     n2=nm+nel+nrrx
      do ng=1,neg
        write(10,'(" *",/,1x,56("*"),"  g r o u p",i4,1x,56("*"))')ng
        write(10,'(" *",/," *",44x,"...principal cross sections...")')
        write(10,'(" *",/," *",26x,"zone",8x,"temperature",5x,"chi",
     1   9x,"nu*fission",6x,"total",
     2   7x,"absorption",/," *",24x,"no.  name")')
        do n=n1,n2
          line=" *   "
c         write(line(27:36),'(i2,"  zone",i2.2)')n,n
          write(line(27:36),'(i2,"  zone",i2.2)')n,0
          write(line(37:106),'(1p5e14.4)')temps,
     1     chi(ng,ng,n),nusigf(ng,n),sigt(ng,n),siga(ng,n)
          write(10,'(a)')trim(line)
        end do ! n
        if(ifissdata.eq.0)then
          write(10,'(" *",/," *",51x,"...fission matrices...",
     1     /," *",/," * zone",1x,"temperature",
     2     1x,"fission matrix ( gprime -->",i4," )")')ng
        else if(ifissdata.eq.1)then
          write(10,'(" *",/," *",39x,"...energy dependent fission ",
     1     "chi...",/," *",/," *   zone",7x,"temperature",
     2     4x,"chi ( gprime -->",i4," )")')ng
        end if
        do n=n1,n2
          line=" *"
          write(line(3:8),'(i6)')n
          if(sum(chi(ng,1:neg,n)).eq.0.)then
            line(21:36)='all entries zero'
            write(10,'(a)')trim(line)
            cycle
          end if 
          write(line(14:27),'(1pe14.4)')temps
          mg=1
  520     do k=28,98,14
c this works for the special 2-mat cross sections (SU).
            if(ifissdata.eq.0)then
              write(line(k:k+13),'(1pe14.4)')nusigf(mg,n)*chi(ng,mg,n)
            else if(ifissdata.eq.1)then
              write(line(k:k+13),'(1pe14.4)')chi(ng,mg,n)
            end if
            mg=mg+1
            if(mg.eq.neg+1)go to 530
          end do ! k
  530     write(10,'(a)')trim(line)
          line=" *   "
          if(mg.lt.neg)go to 520
        end do ! n
        if(i2lp1.eq.0)then
          write(10,'(" *",/," *",49x,"...scattering matrices...",/,
     1     " *",52x,"(2l+1 not included)",/," *")')
        else
          write(10,'(" *",/," *",49x,"...scattering matrices...",/,
     1     " *",52x,"(2l+1 included)",/," *")')
        end if
        write(10,'(" *   zone   order  first grp",
     1   "  temperature    cross sections")')
        do n=n1,n2
          do l=0,isct
            mg=scgr(ng,l,n,2)
            if(mg.eq.0)cycle
            line=" *   "
            write(line(6:41),'(i3,6x,i1,4x,i3,5x,1pe14.4)')n,l,
     1       scgr(ng,l,n,2),temps
  620       do k=45,115,14
              write(line(k:k+13),'(1pe14.4)')sigs(ng,mg,l,n)
              mg=mg-1
              if(mg.eq.0.or.mg+1.eq.scgr(ng,l,n,1))go to 630
            end do ! k
  630       write(10,'(a)')trim(line)
            line=" *   "
            if(mg.ne.0.and.mg+1.ne.scgr(ng,l,n,1))go to 620
          end do ! l
        end do ! n
      end do ! ng
c
c check scgr.
      write(*,'(/,"scgr(ng,l,ir,1/2) are the top and bottom groups",
     1 " scattering into ng.")')
      if(isct.eq.0)then
        write(*,'(7x,1(5x,"l=",i1))')(l,l=0,isct)
      else if(isct.eq.1)then
        write(*,'(7x,2(5x,"l=",i1))')(l,l=0,isct)
      else if(isct.eq.2)then
        write(*,'(7x,3(5x,"l=",i1))')(l,l=0,isct)
      else if(isct.eq.3)then
        write(*,'(7x,4(5x,"l=",i1))')(l,l=0,isct)
      else if(isct.eq.4)then
        write(*,'(7x,5(5x,"l=",i1))')(l,l=0,isct)
      else if(isct.eq.5)then
        write(*,'(7x,6(5x,"l=",i1))')(l,l=0,isct)
      end if
      write(*,'(3x,"n  ng",5a8)')("   1   2",l=0,isct)
      do n=0,nm+1
        do ng=1,neg
          write(*,'(12i4)')n,ng,(scgr(ng,l,n,1),scgr(ng,l,n,2),l=0,isct)
        end do ! ng
      end do ! n
c
      close(10)
      call stoponerror ! debug stop
c
  990 write(*,'("done.")')
c
#if defined ( SU )
c write macroscopic cross sections for input to makemg.
      go to 992 ! skip this whole thing.
      if(neg.eq.30.and.lng.eq.30)then
        bfile=" "
        do ng=0,neg+1
          jr=1
          if(ng.gt.0)jr=8
          do l=n1,n2
            ln=1 ! small perturbations
            if((l.eq.1.and.(ng.le.1.or.ng.gt.30)).or.
     1         (l.eq.2.and.(ng.le.1.or.ng.gt.30)))ln=2 ! big perturbations
            if(ng.eq.31)ln=1
            do k=1,jr
              write(bfile,'("i",3i2.2)')ng,l,k
c             write(*,'(a10,i3)')bfile,ln
              open(iut0,file=bfile,status='unknown')
              if(ng.eq.0)then
                write(iut0,'(2i1)')l,0
                write(iut0,'("xs",2i1)')l,0
                write(iut0,'("unpert. mat",i2)')l
              else
                write(iut0,'(2i1)')l,k
                write(iut0,'("xs",2i1)')l,k
                write(iut0,'("mat",i2," group",i3," pert",i2,":",0pf8.1,
     1           "%")')l,ng,k,(prt(k,ln)-1.d0)*100.d0
              end if
              write(iut0,'(i2.2)')neg
c energies
              write(iut0,'("17 15 13.5 12 10. 7.79 6.07 3.68 2.865 ",
     1         "2.232",
     1         " 1.738 1.353 .823 .5 .303 .184 .0676 .0248 9.12E-03 ",
     2         "3.35E-03 1.235E-03 4.54E-04 1.67E-04 6.14E-05 2.26E-05",
     3         " 8.32E-06 3.06E-06 1.13E-06 4.14E-07 1.52E-07 1.39E-10"
     4         )')
c look at abs
              xst(1:neg)=sigt(1:neg,l)
              xsa(1:neg)=siga(1:neg,l)-nusigf(1:neg,l)/nubar(1:neg)
              if(ng.eq.0)then
                do lg=1,neg
                  if(xsa(lg).ge.0.)cycle
                  if(ng.eq.0)
     1             write(*,'("negative siga reset. lg,l,siga,nusigf,",
     2             "nu,sigf",2i3,1p20e14.6)')lg,l,siga(lg,l),
     3             nusigf(lg,l),nubar(lg),nusigf(lg,l)/nubar(lg),
     4             xsa(lg)
                   san=1.d-4 ! macroscopic u235
                  xst(lg)=xst(lg)-xsa(lg)+san
                  xsa(lg)=san
                end do ! lg
                sigt(1:neg,l)=xst(1:neg)
                siga(1:neg,l)=xsa(1:neg)
              end if
c
c total
c             write(iut0,'("sigt")')
              if(ng.gt.0.and.ng.le.neg)then
                do lg=ng,neg
                  xst(ng)=xst(ng)-sigs(lg,ng,0,l)+sigs(lg,ng,0,l)
     1             *prt(k,ln)
                end do ! lg
              else if(ng.eq.neg+1)then
                do mg=1,neg
                  do lg=mg,neg
                    xst(mg)=xst(mg)-sigs(lg,mg,0,l)+sigs(lg,mg,0,l)
     1               *prt(k,ln)
                  end do ! lg
                end do ! mg
              end if
              do lg=1,neg
                write(iut0,'(1pe20.12)')xst(lg)
              end do ! lg
c fission
c             write(iut0,'("sigf")')
              do lg=1,neg
                write(iut0,'(1pe20.12)')nusigf(lg,l)/nubar(lg)
              end do ! lg
c nubar
              if(l.eq.1)then
c               write(iut0,'("nubar")')
                do lg=1,neg
                  write(iut0,'(1pe20.12)')nubar(lg)
                end do ! lg
c chi
c               write(iut0,'("chi")')
                do lg=1,neg
                  write(iut0,'(1pe20.12)')chi(lg,lg,l)
                end do ! lg
              end if
c absorption
c             write(iut0,'("siga")')
              do lg=1,neg
                write(iut0,'(1pe20.12)')xsa(lg)
              end do ! lg
              write(iut0,'("n")') ! upscatter
              write(iut0,'("n")') ! anisotropic
              write(iut0,'("0")') ! max order
c scat, group 1 to 1
c scat, group 1 to 2
c etc.
              do lg=1,neg
c             write(iut0,'("sigs from group",i3)')lg
                do mg=lg,neg
                  if(lg.eq.ng.or.ng.eq.neg+1)then
                    write(iut0,'(1pe20.12)')sigs(mg,lg,0,l)*prt(k,ln)
                  else
                    write(iut0,'(1pe20.12)')sigs(mg,lg,0,l)
                  end if
                end do ! mg
              end do ! lg
              close(iut0)
c
c write partisn input files
              if(ng.eq.0.and.l.eq.1)cycle
              xstt(1:neg,1:nm)=sigt(1:neg,1:nm)
              xsat(1:neg,1:nm)=siga(1:neg,1:nm)
              sigt(1:neg,l)=xst(1:neg)
              siga(1:neg,l)=xsa(1:neg)
              call wrdantep(1,ng,l,k,ln)
              sigt(1:neg,1:nm)=xstt(1:neg,1:nm)
              siga(1:neg,1:nm)=xsat(1:neg,1:nm)
c              
            end do ! k
          end do ! l
        end do ! ng
        stop
      end if
  992 continue
#endif
c
  900 format(a130)
      return
      end subroutine
