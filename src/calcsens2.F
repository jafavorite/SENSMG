      subroutine calcsens2
c compute second-order sensitivities.
c
c equation for computing the 2nd-order sensitivities is
c Ref. 1: D. G. Cacuci, "Application of the Second-order Comprehensive
c Adjoint Sensitivity Analysis Methodology to Compute 1st- and 2nd-Order
c Sensitivities of Flux Functionals in a Multiplying System With
c Source," Nucl. Sci. Eng., 193, 555-600 (2019); 
chttps://doi.org/10.1080/00295639.2018.1553910.
c
c expressions derived from the equations in Ref. 1 can also be found 
c from the following references:
c Ref. 2: D. G. Cacuci, R. Fang, and J. A. Favorite,"Comprehensive
c Second-Order Adjoint Sensitivity Analysis Methodology (2nd-ASAM)
C Applied to a Subcritical Experimental Reactor Physics Benchmark:
c I. Effects of Imprecisely Known Microscopic Total and Capture 
c Cross Sections," Energies, 12, 21 (2019);
c https://doi.org/10.3390/en12214219.
c
c Ref. 3: R. Fang and D. G. Cacuci,"Comprehensive Second-Order Adjoint
c Sensitivity Analysis Methodology (2nd-ASAM) Applied to a Subcritical
c Experimental Reactor Physics Benchmark: II. Effects of Imprecisely
c Known Microscopic Scattering Cross Sections," Energies, 12, 21 (2019);
c https://doi.org/10.3390/en12214114.
c
c Ref. 4: D. G. Cacuci, R. Fang, J. A. Favorite, M. C. Badea, and
c F. Di Rocco, "Comprehensive Second-Order Adjoint Sensitivity Analysis
c Methodology (2nd-ASAM) Applied to a Subcritical Experimental Reactor
c Physics Benchmark: III. Effects of Imprecisely Known Microscopic 
c Fission Cross Sections and Average Number of Nuetrons per Fission,"
c Energies, 12, 21 (2019); https://doi.org/10.3390/en12214100.
c 
c note: the order of parameters are important. In this subroutine,
c "nf" is used to index to the correct files to read from for 2nd-order
c adjoint functions. Each file name in fo2/ad2 contains a number which      
c corresponds to a specific parameter. Refer to "src2o.F" for the 
c description of the order of the parameters.
c
      use COMS
      implicit none
      real*8 dtg(neg,3,nr,nz),dfg(3,neg,neg,nr,nz),
     1 dsg(3,neg,neg,0:isct,nr,nz),rlip(0:isct),
     2 s2tt(0:nrrr,neg*nel,neg*nel),
c s2tt means the matrix of 2nd-order sensitivities w.r.to sigt and sigt
c s2tf means the matrix of 2nd-order sensitivities w.r.to sigt and sigf
     3 s2tf(0:nrrr,neg*nel,neg*nel),
     4 s2ts(0:nrrr,neg*nel,neg*neg*nmom*nel),
     5 s2nut(0:nrrr,neg*nel,neg*nel), 
     6 s2tnu(0:nrrr,neg*nel,neg*nel),
     7 s2nunu(0:nrrr,neg*nel,neg*nel),
     8 s2nuf(0:nrrr,neg*nel,neg*nel),
     9 s2nus(0:nrrr,neg*nel,neg*neg*nmom*nel),               
     1 s2ft(0:nrrr,neg*nel,neg*nel),
     2 s2fs(0:nrrr,neg*nel,neg*neg*nmom*nel),
     3 s2fnu(0:nrrr,neg*nel,neg*nel),
     4 s2st(0:nrrr,neg*neg*nmom*nel,neg*nel),
     5 s2sf(0:nrrr,neg*neg*nmom*nel,neg*nel),
     6 s2snu(0:nrrr,neg*neg*nmom*nel,neg*nel),
     7 s2ff(0:nrrr,neg*nel,neg*nel), 
c ideally, s2ss should be defined as s2ss(0:nrrr,1:neg*neg*nmom*nel,
c 1:neg*neg*nmom*nel), but that will require a huge memory block. 
c instead it is defined as follows 
     8 s2ss(0:nrrr,neg*neg,neg*neg*nmom*nel,0:isct),
     9 s2ssdiag(0:nrrr,neg,neg,0:isct,nel)     
      real*8 st,fissd
      integer neg0,ndir0,nitm0,njtm0
      integer i,ir,iso,jr,k2,l,m,n,nf,ng,niso0,nedt0,lg,i1,j1,noff2,
     1 noff3,noff4,noff5,nel2,tj,m2,ng1,lg1,i2,j2,iso2,k3,j,l1,noff6,
     2 noff7,noff8,k4,m4
      logical okf,oka  
      character afile*70,sfilett*70,sfiletf*70,sfiletnu*70,sfilets*70,
     1 sfilenut*70,sfilenunu*70,sfilenuf*70,sfilenus*70,sfileft*70,
     2 sfilefs*70,sfilefnu*70,sfileff*70,sfilest*70,sfilesf*70,
     3 sfilesnu*70,sfiless*70,head*33
c
      head="relative sensitivities, 2nd order"
c read sources
      if(icalc.eq.0.or.icalc.eq.3.or.icalc.eq.4)then
        if(imisc.eq.1)then
          call rdmisc
        end if
        if(imisc.eq.0.or.ialphan.eq.1)then
          call rdsources
        end if
      end if
c read forward flux moments and leakage
      call rddantm(1,0)
      if(iangflux.eq.1)call rddanta(1,0)
      call rddantk
c read transport cross sections for all isotopes
      call rdxsmg
c read reaction-rate edit cross sections
      call rdsnxedth(1,iuo,iue,niso0,nedt0)
      call rdxsmg0(1,niso0,nedt0)      
c read regular adjoint flux moments
      call rddantm(2,0)
      if(iangflux.eq.1)call rddanta(2,0)
      open(iut,file='stopconverged',status='unknown')
      close (iut)
c
c adjoint moments (regular and generalized) are always calculational
c moments. the sign of rlip is used to convert to inner product moments.
      do l=0,isct
        rlip(l)=real(2*l+1)
        if(mod(l,2).ne.0)then
          rlip(l)=-rlip(l)
        end if
      end do ! l
c
c this is to dimension things correctly in subroutine integral_2nd.
      if(iangflux.eq.1)then
        neg0=neg
        ndir0=ndir
        nitm0=nitm
        njtm0=njtm
      else if(iangflux.eq.0)then
        neg0=1
        ndir0=1
        nitm0=1
        njtm0=1
      end if
c
c offsets
c the offsets are the same as defined in "src2o.F". the following offsets
c are used to help to find the correct file indexed by "nf".
c
c noff1 - not used in thie routine. 
c
c noff2 = neg*nel
c        used to offset between sigt, nu, sigf, and sigs:
c        "1*noff2" offsets to the first parameter of nu;
c        "2*noff2" offsets to the first parameter of sigf;
c        "3*noff2" offsets to the first parameter of sigs.     
c
c noff3 = neg*(k2-1)
c        used to offset between isotopes for the parameters of sigt, nu,
c        and sigf, respectively.
c
c noff4 = neg*neg*nmom*(k2-1)
c        used to offset between isotopes for the parameters of sigs
c
c noff5 = neg*neg*l
c        used to offset between the l-th (l=0,...,nmom-1) order of sigs     
c
c briefly, the parameter table consists of the following blocks (in order):
c       neg*nel - for sigt
c       neg*nel - for nu
c       neg*nel - for sigf
c       neg*neg*nmom*nel - for sigs  
c
      noff2=neg*nel 
c
c compute and write sensitivities.
c n=0 is for leakage sensitivities,
      do n=0,nrrr
c
        write(*,'(/,"initializing arrays ...")')
        nel2=nel
        s2tt(0:nrrr,1:neg*nel,1:neg*nel)=0.d0
        s2tf(0:nrrr,1:neg*nel,1:neg*nel)=0.d0
        s2tnu(0:nrrr,1:neg*nel,1:neg*nel)=0.d0
        s2ts(0:nrrr,1:neg*nel,1:neg*neg*nmom*nel)=0.d0
        s2nut(0:nrrr,1:neg*nel,1:neg*nel)=0.d0
        s2nunu(0:nrrr,1:neg*nel,1:neg*nel)=0.d0
        s2nuf(0:nrrr,1:neg*nel,1:neg*nel)=0.d0                
        s2nus(0:nrrr,1:neg*nel,1:neg*neg*nmom*nel)=0.d0
        s2ft(0:nrrr,1:neg*nel,1:neg*nel)=0.d0
        s2fs(0:nrrr,1:neg*nel,1:neg*neg*nmom*nel)=0.d0
        s2fnu(0:nrrr,1:neg*nel,1:neg*nel)=0.d0
        s2ff(0:nrrr,1:neg*nel,1:neg*nel)=0.d0
        s2st(0:nrrr,1:neg*neg*nmom*nel,1:neg*nel)=0.d0
        s2sf(0:nrrr,1:neg*neg*nmom*nel,1:neg*nel)=0.d0
        s2snu(0:nrrr,1:neg*neg*nmom*nel,1:neg*nel)=0.d0
        s2ss(0:nrrr,1:neg*neg,1:neg*neg*nmom*nel,0:isct)=0.d0
        s2ssdiag(0:nrrr,1:neg,1:neg,0:isct,1:nel2)=0.d0
c                                                   
        if(n.eq.0)then
          if(icalc.eq.0)then
c write 2nd-order sensitivities to separate files
            sfilett="sens2_l_sigt_sigt"
            sfiletf="sens2_l_sigt_sigf"
            sfiletnu="sens2_l_sigt_nu"
            sfilets="sens2_l_sigt_sigs"
            sfilenut="sens2_l_nu_sigt"
            sfilenunu="sens2_l_nu_nu"
            sfilenuf="sens2_l_nu_sigf"
            sfilenus="sens2_l_nu_sigs"
            sfileft="sens2_l_sigf_sigt"
            sfilefs="sens2_l_sigf_sigs"
            sfilefnu="sens2_l_sigf_nu"
            sfileff="sens2_l_sigf_sigf"
            sfilest="sens2_l_sigs_sigt"
            sfilesf="sens2_l_sigs_sigf"
            sfilesnu="sens2_l_sigs_nu"
            sfiless="sens2_l_unmixed_sigs_sigs"
          end if
        end if
        if(n.le.1)then
          write(*,'(/,"writing second-order sensitivities.")')
          write(iuo,'(/,"writing second-order sensitivities.")')
          open(ius2tt,file=sfilett,status='unknown')
          open(ius2tf,file=sfiletf,status='unknown')
          open(ius2tnu,file=sfiletnu,status='unknown')
          open(ius2ts,file=sfilets,status='unknown')
          open(ius2nut,file=sfilenut,status='unknown')
          open(ius2nunu,file=sfilenunu,status='unknown')
          open(ius2nuf,file=sfilenuf,status='unknown')
          open(ius2nus,file=sfilenus,status='unknown')
          open(ius2ft,file=sfileft,status='unknown')
          open(ius2fs,file=sfilefs,status='unknown')
          open(ius2fnu,file=sfilefnu,status='unknown')
          open(ius2ff,file=sfileff,status='unknown')
          open(ius2st,file=sfilest,status='unknown')
          open(ius2sf,file=sfilesf,status='unknown')
          open(ius2snu,file=sfilesnu,status='unknown')
          open(ius2ss,file=sfiless,status='unknown')
          write(ius2tt,'(a)')trim(id)
          write(ius2tf,'(a)')trim(id)
          write(ius2tnu,'(a)')trim(id)
          write(ius2ts,'(a)')trim(id)
          write(ius2nut,'(a)')trim(id)
          write(ius2nunu,'(a)')trim(id)
          write(ius2nuf,'(a)')trim(id)
          write(ius2nus,'(a)')trim(id)
          write(ius2ft,'(a)')trim(id)
          write(ius2fs,'(a)')trim(id)
          write(ius2fnu,'(a)')trim(id)
          write(ius2ff,'(a)')trim(id)
          write(ius2st,'(a)')trim(id)
          write(ius2sf,'(a)')trim(id)
          write(ius2snu,'(a)')trim(id)
          write(ius2ss,'(a)')trim(id)
          write(ius2tt,'(/,a,/,"sigt-sigt")')trim(head)
          write(ius2tf,'(/,a,/,"sigt-sigf")')trim(head)
          write(ius2tnu,'(/,a,/,"sigt-nu")')trim(head)
          write(ius2ts,'(/,a,/,"sigt-sigs")')trim(head)
          write(ius2nunu,'(/,a,/,"nu-nu")')trim(head)
          write(ius2nut,'(/,a,/,"nu-sigt")')trim(head)
          write(ius2nuf,'(/,a,/,"nu-sigf")')trim(head)
          write(ius2nus,'(/,a,/,"nu-sigs")')trim(head)
          write(ius2ft,'(/,a,/,"sigf-sigt")')trim(head)
          write(ius2fs,'(/,a,/,"sigf-sigs")')trim(head)
          write(ius2fnu,'(/,a,/,"sigf-nu")')trim(head)
          write(ius2ff,'(/,a,/,"sigf-sigf")')trim(head)
          write(ius2st,'(/,a,/,"sigs-sigt")')trim(head)
          write(ius2sf,'(/,a,/,"sigs-sigf")')trim(head)
          write(ius2snu,'(/,a,/,"sigs-nu")')trim(head)
          write(ius2ss,'(/,a,/,"sigs-sigs",
     1     " for unmixed 2nd-order sensitivities to the same scattering"
     2     " cross sections only,"/,"which are the diagonal components"
     3     " of the full matrix s2ss.")')trim(head)
        end if
c      
        k2=0
        do i1=1,nm
          do j1=1,ncb(i1)
c for testing
!        do i1=1,1
!          do j1=1,1
            k2=k2+1
            iso=nm+k2 ! nm is the offset to the isotope cross sections
c offset between isotopes
            noff3=neg*(k2-1)
c
c sigt
c
c loop over all the parameters in vector sigt=[sigt(1,1),sigt(2,1),...]
c indexed by tj,where tj=1,...,Jsigt.
            do ng=1,neg
              tj=neg*(k2-1)+ng 
c
c reading 2nd-order adjoint functions psi(2,g;1,j) and psi(2,g;2,j)
c computed from Eqs. (164) and (165), respectively. Use "nf" to find 
c the correct files to read from.  
c inquire to see if the files indexed by nf exist or not. nf indicates
c the "at most" forward and adjoint files. for sigt, all the files
c indexed by nf should exist. In some cases,such as for sigf or nu,  
c forward and adjoint files only exist for fissionable isotopes;  
c similarly, for sigs, files only exist for nonzero sigs. Therefore,it 
c is necessary to do the inquiry before computing the integrals.  
c
              nf=noff3+ng
              write(afile,'("fo2/",i5.5,"_rmflux")')nf
              inquire(file=afile, exist=okf)
              if(okf)then 
                call rddantm2o(1,nf,n)
                if(iangflux.eq.1)call rddanta2o(1,nf,n)
              endif
              write(afile,'("ad2/",i5.5,"_amflux")')nf
              inquire(file=afile, exist=oka)
              if(oka)then
                call rddantm2o(2,nf,n)
                if(iangflux.eq.1)call rddanta2o(2,nf,n)
              endif
c the subroutine "integrals_2nd" is used to compute the integrals of 
c dtg, dfg, and dsg, for each parameter tj
              if(okf.or.oka)then
                call integrals_2nd(0,rlip,fmom,amom,
     1           afreg,afadj,fmom2,amom2,afreg2,afadj2,neg0,ndir0,
     2           nitm0,njtm0,dtg,dsg,dfg)
c to compute the mixed second derivatives of leakage response to                    
c parameters indexed by tj and m2, keeping tj unchanged, loop over  
c the parameters indexed by m2.
                k3=0
                do i2=1,nm
                  do j2=1,ncb(i2)
                    k3=k3+1
c isotope corresponds to parameters indexed by m2
                    iso2=nm+k3 
                    noff5=neg*neg*nmom*(k3-1)
                    do ir=1,nr
                      do jr=1,nz
                        if(mat(ir,jr).ne.i2)cycle
                        do ng1=1,neg
                          m2=neg*(k3-1)+ng1                    
c Eq. (158) to sigt-sigt
c where the first sigt indexed by tj, and the second sigt indexed by m2,
c with m2=1,...,Jsigt;
c for s2tt (i.e., 2nd-order sensitivities with respect to sigt and sigt),
c Eq. (158) is the only contribution, therefore, this is the final result.
                          s2tt(n,tj,m2)=s2tt(n,tj,m2)
     1                     -sum(dtg(ng1,2:3,ir,jr))*blk(3,k3)
     2                     *sigt(ng,iso)*sigt(ng1,iso2)
c Eq. (158) to sigt-sigf
c where m2=1,...,Jsigf
c s2tf comprises 2 contributions: Eq158tf+Eq160tf
c here is the 1st contribution resulted from Eq. (158) to s2tf
                          s2tf(n,tj,m2)=s2tf(n,tj,m2)
     1                     -sum(dtg(ng1,2:3,ir,jr))*blk(3,k3)
     2                     *sigt(ng,iso)*sigf(ng1,iso2)
c Eq. (158) to sigf-sigt
c s2ft comprises 2 contributions: Eq158ft+Eq177ft
c here is the 1st contribution resulted from Eq. (158) to s2ft                   
                          s2ft(n,tj,m2)=s2ft(n,tj,m2)
     1                     -sum(dtg(ng1,2:3,ir,jr))*blk(3,k3)
     2                     *sigf(ng,iso)*sigt(ng1,iso2)
c Eq. (158) for sigf-sigf
c s2ff comprises 4 contributions: Eq158ff+Eq160ff+Eq177ff+Eq179ff
c here is the 1st contribution resulted from Eq. (158) to s2ff
                          s2ff(n,tj,m2)=s2ff(n,tj,m2)
     1                     -sum(dtg(ng1,2:3,ir,jr))*blk(3,k3)
     2                     *sigf(ng,iso)*sigf(ng1,iso2)
c Eq. (160) for sigt-sigf                   
c s2tf comprises 2 contributions: Eq158tf+Eq160tf
c here is the 2nd contribution resulted from Eq. (160) to s2tf.
c thus, the computation for s2ff is completed. 
                          do lg1=1,neg
                            if(sigf(ng1,iso2).ne.0.0)then
                              st=sum(dfg(2:3,lg1,ng1,ir,jr))
     1                         *chivec(lg1,i2)*blk(3,k3)
     2                         *nusigf(ng1,iso2)/sigf(ng1,iso2)
                              s2tf(n,tj,m2)=s2tf(n,tj,m2)
     1                         +st*sigt(ng,iso)*sigf(ng1,iso2)
                            end if
                          end do ! lg1
c Eq. (160) for sigt-nu
c where m2=1,...,Jnu                   
c for s2tnu, Eq. (160) is the only contribution, therefore, this is the
c final result
                          do lg1=1,neg
                            if(sigf(ng1,iso2).ne.0.0)then
                              st=sum(dfg(2:3,lg1,ng1,ir,jr))
     1                         *chivec(lg1,i2)*blk(3,k3)
     2                         *nusigf(ng1,iso2)*sigt(ng,iso)
                              s2tnu(n,tj,m2)=s2tnu(n,tj,m2)
     1                         +st
                            end if
                          end do ! lg1
c Eq. (160) to sigf-sigf                   
c s2ff comprises 4 contributions: Eq158ff+Eq160ff+Eq177ff+Eq179ff
c here is the 2nd contribution resulted from Eq. (160) to s2ff
                          do lg1=1,neg
                            if(sigf(ng1,iso2).ne.0.0)then
                              st=sum(dfg(2:3,lg1,ng1,ir,jr))
     1                         *chivec(lg1,i2)*blk(3,k3)
     2                         *nusigf(ng1,iso2)/sigf(ng1,iso2)
                              s2ff(n,tj,m2)=s2ff(n,tj,m2)
     1                         +st*sigf(ng,iso)*sigf(ng1,iso2)
                            end if
                          end do ! lg1
c Eq. (160) for sigf-nu
c s2fnu comprises 2 contributions: Eq160fnu+Eq179fnu
c here is the 1st contribution resulted from Eq. (160) to s2fnu                   
                          do lg1=1,neg
                            if(sigf(ng1,iso2).ne.0.0)then
                              st=sum(dfg(2:3,lg1,ng1,ir,jr))
     1                         *chivec(lg1,i2)*blk(3,k3)
     2                         *sigf(ng1,iso2)
                              s2fnu(n,tj,m2)=s2fnu(n,tj,m2)
     1                         +st*sigf(ng,iso)
     2                         *nusigf(ng1,iso2)/sigf(ng1,iso2)
                            end if
                          end do ! lg1
c Eq. (159) + Eq. (158) for sigt-sigs
c where m2=1,...,Jsigs
c s2ts comprises 2 contributions: Eq158ts+Eq159ts
c both contributions are accounted for as follows.
c thus, the computation for s2ts is completed.
                          do lg1=1,neg 
                            do l=0,isct
c Eq. 159ts
                              m2=noff5+neg*neg*l+neg*(ng1-1)+lg1 
                              st=dsg(2,ng1,lg1,l,ir,jr)*blk(3,k3)
     1                         +dsg(3,ng1,lg1,l,ir,jr)*blk(3,k3)
                               if(l.eq.0)then
c Eq. (158) for the contribution of sigs(l=0) to sigt
                                 st=st-sum(dtg(lg1,2:3,ir,jr))
     1                            *blk(3,k3)
                               end if     
                              s2ts(n,tj,m2)=
     1                         s2ts(n,tj,m2)+st
     2                         *sigt(ng,iso)*sigs(ng1,lg1,l,iso2)
                            end do ! l
                          end do ! lg1
c Eq. (159) + Eq. (158) for sigf-sigs
c s2fs comprises 4 contributions: Eq158fs+Eq159fs+Eq178fs+Eq177fs
c here is the 1st and 2nd contributions resulted from Eq. (158) and Eq. (159) 
                          do lg1=1,neg 
                            do l=0,isct
c Eq. 159fs
                              m2=noff5+neg*neg*l+neg*(ng1-1)+lg1
                              st=dsg(2,ng1,lg1,l,ir,jr)*blk(3,k3)
     1                         +dsg(3,ng1,lg1,l,ir,jr)*blk(3,k3)
                              if(l.eq.0)then
c Eq. (158) for the contributions of both sigf and sigs(l=0) to sigt
                                 st=st-sum(dtg(lg1,2:3,ir,jr))
     1                            *blk(3,k3)
                              end if     
                              s2fs(n,tj,m2)=
     1                         s2fs(n,tj,m2)+st*sigf(ng,iso)
     2                         *sigs(ng1,lg1,l,iso2)
                            end do ! l
                          end do ! lg1
c
                        end do ! ng1
                      end do ! ir
                    end do ! jr
                  end do ! j2
                end do ! i2
              end if ! okf.or.oka
            end do ! ng 
c
c nu
c
c loop over all the parameters in vector nu=[nu(1,1),nu(2,1),...]
c indexed by tj,where tj=1,...,Jnu.
            do ng=1,neg
              tj=neg*(k2-1)+ng 
c
c reading 2nd-order adjoint functions u(2,g;1,j) and u(2,g;2,j)
c computed from Eqs. (183) and (184), respectively. Use "nf" to find the
c correct files to read from.  
              nf=1*noff2+noff3+ng 
              write(afile,'("fo2/",i5.5,"_rmflux")')nf
              inquire(file=afile, exist=okf)
              if(okf)then 
                call rddantm2o(1,nf,n)
                if(iangflux.eq.1)call rddanta2o(1,nf,n)
              endif
              write(afile,'("ad2/",i5.5,"_amflux")')nf
              inquire(file=afile, exist=oka)
              if(oka)then
                call rddantm2o(2,nf,n)
                if(iangflux.eq.1)call rddanta2o(2,nf,n)
              endif
c the subroutine "integrals_2nd" is used to compute the integrals of 
c dtg, dfg, and dsg, for each parameter nu_tj
              if(okf.or.oka)then
                call integrals_2nd(0,rlip,fmom,amom,
     1           afreg,afadj,fmom2,amom2,afreg2,afadj2,neg0,ndir0,
     2           nitm0,njtm0,dtg,dsg,dfg)
c to compute the mixed second derivatives: keep nu_tj unchanged,
c loop over nu_m2 
                k3=0
                do i2=1,nm
                  do j2=1,ncb(i2)
                    k3=k3+1
                    iso2=nm+k3
                    noff5=neg*neg*nmom*(k3-1)
                    do ir=1,nr
                      do jr=1,nz
                        if(mat(ir,jr).ne.i2)cycle
                        do ng1=1,neg
                          m2=neg*(k3-1)+ng1                    
c Eq. (177) to nu-sigt
c for s2nut, Eq. (177) is the only contribution, therefore, this is the
c final result                   
                          if(sigf(ng,iso).ne.0.0)then
                            s2nut(n,tj,m2)=s2nut(n,tj,m2)
     1                       -sum(dtg(ng1,2:3,ir,jr))*blk(3,k3)
     2                       *nusigf(ng,iso)/sigf(ng,iso)
     3                       *sigt(ng1,iso2)
                          end if
c Eq. (177) to nu-sigf
c s2nuf comprises 2 contributions: Eq177nuf+Eq179nuf
c here is the 1st contribution resulted from Eq. (177) to s2nuf                    
                          if(sigf(ng,iso).ne.0.0)then
                            s2nuf(n,tj,m2)=s2nuf(n,tj,m2)
     1                       -sum(dtg(ng1,2:3,ir,jr))*blk(3,k3)
     2                       *nusigf(ng,iso)/sigf(ng,iso)
     3                       *sigf(ng1,iso2)
                          end if
c Eq. (178) + Eq. (177) to nu-sigs
c s2nus comprises 2 contributions: Eq177nus+Eq178nus
c both contributions are accounted for as follows.
c thus, the computation for s2nus is completed.
                          do lg1=1,neg 
                            do l=0,isct
                              m2=noff5+neg*neg*l+neg*(ng1-1)+lg1
c Eq. (178)
                              st=dsg(2,ng1,lg1,l,ir,jr)*blk(3,k3)
     1                         +dsg(3,ng1,lg1,l,ir,jr)*blk(3,k3)
c Eq. (177) for the contribution of sigs(l=0) to sigt
                               if(l.eq.0)then
                                 st=st-sum(dtg(lg1,2:3,ir,jr))
     1                            *blk(3,k3)
                               end if     
                             if(sigf(ng,iso).ne.0.0)then
                               s2nus(n,tj,m2)=
     1                          s2nus(n,tj,m2)+st*nusigf(ng,iso)
     2                          /sigf(ng,iso)*sigs(ng1,lg1,l,iso2)
                              end if
                            end do ! l
                          end do ! lg1
c Eq. (179) to nu-nu                   
c for s2nunu, Eq. (179) is the only contribution, therefore, this is the
c final result 
                          m2=neg*(k3-1)+ng1
                          do lg1=1,neg
                            if(sigf(ng,iso).ne.0.0)then
                              st=sum(dfg(2:3,lg1,ng1,ir,jr))
     1                         *chivec(lg1,i2)*blk(3,k3)
     2                         *nusigf(ng1,iso2)
     3                         *nusigf(ng,iso)/sigf(ng,iso)
                              s2nunu(n,tj,m2)=
     1                         s2nunu(n,tj,m2)+st
                            end if
                          end do ! lg1
c Eq. (179) to nu-sigf
c s2nuf comprises 2 contributions: Eq177nuf+Eq179nuf
c here is the 2nd contribution resulted from Eq. (179) to s2nuf,
c thus, the computation for s2nuf is completed.                    
                          do lg1=1,neg
                            if(sigf(ng1,iso2).ne.0.0)then
                              st=sum(dfg(2:3,lg1,ng1,ir,jr))
     1                         *chivec(lg1,i2)*blk(3,k3)
     2                         *nusigf(ng1,iso2)/sigf(ng1,iso2)
                              if(iso2.eq.iso.and.ng1.eq.ng)then
                                st=st+dfg(1,lg1,ng1,ir,jr)
     1                           *chivec(lg1,i2)*blk(3,k3)
                              end if
                            end if        
                            if(sigf(ng,iso).ne.0.0)then
                              s2nuf(n,tj,m2)=s2nuf(n,tj,m2)
     1                         +st*nusigf(ng,iso)/sigf(ng,iso)
     2                         *sigf(ng1,iso2)
                            end if
                          end do ! lg1
ccc
                        end do ! ng1
                      end do ! ir
                    end do ! jr
                  end do ! j2
                end do ! i2
              end if ! okf.or.oka
            end do ! ng 
c
c sigf
c
c loop over all the parameters in vector sigf=[sigf(1,1),sigf(2,1),...]
c indexed by tj, where tj=1,...,Jsigf.
            do ng=1,neg
              tj=neg*(k2-1)+ng 
c
c reading 2nd-order adjoint functions u(2,g;1,j) and u(2,g;2,j)
c computed from Eqs. (183) and (184), respectively. Use "nf" to find the
c correct files to read from.  
              nf=2*noff2+noff3+ng 
              write(afile,'("fo2/",i5.5,"_rmflux")')nf
              inquire(file=afile, exist=okf)
              if(okf)then 
                call rddantm2o(1,nf,n)
                if(iangflux.eq.1)call rddanta2o(1,nf,n)
              endif
              write(afile,'("ad2/",i5.5,"_amflux")')nf
              inquire(file=afile, exist=oka)
              if(oka)then
                call rddantm2o(2,nf,n)
                if(iangflux.eq.1)call rddanta2o(2,nf,n)
              endif
c the subroutine "integrals_2nd" is used to compute the integrals of 
c dtg, dfg, and dsg, for each parameter tj
              if(okf.or.oka)then
                call integrals_2nd(0,rlip,fmom,amom,
     1           afreg,afadj,fmom2,amom2,afreg2,afadj2,neg0,ndir0,
     2           nitm0,njtm0,dtg,dsg,dfg)
c to compute the mixed second derivatives: keep sigf_tj unchanged,
c loop over m2 
                k3=0
                do i2=1,nm
                  do j2=1,ncb(i2)
                    k3=k3+1
                    iso2=nm+k3
                    noff4=neg*neg*(k3-1)
                    noff5=neg*neg*nmom*(k3-1)
                    do ir=1,nr
                      do jr=1,nz
                        if(mat(ir,jr).ne.i2)cycle
                        do ng1=1,neg
                          m2=neg*(k3-1)+ng1                    
c Eq. (177) for sigf-sigt
c s2ft comprises two contributions: Eq158ft+Eq177ft
c here is the 2nd contribution resulted from Eq. (177) to s2ft
c thus, the computation for s2ft is completed.                     
                          if(sigf(ng,iso).ne.0.0)then
                            s2ft(n,tj,m2)=s2ft(n,tj,m2)
     1                       -sum(dtg(ng1,2:3,ir,jr))*blk(3,k3)
     2                       *sigf(ng,iso)*sigt(ng1,iso2)
                          end if
c Eq. (177) for sigf-sigf
c s2ff comprises 4 contributions: Eq158ff+Eq160ff+Eq177ff+Eq179ff
c here is the 3rd contribution resulted from Eq. (177) to s2ff
                          if(sigf(ng,iso).ne.0.0)then
                            s2ff(n,tj,m2)=s2ff(n,tj,m2)
     1                       -sum(dtg(ng1,2:3,ir,jr))*blk(3,k3)
     2                       *sigf(ng,iso)*sigf(ng1,iso2)
                          end if
c Eq. (178) + Eq. (177) for sigf-sigs
c s2fs comprises 4 contributions: Eq158fs+Eq159fs+Eq178fs+Eq177fs
c here are the 3rd and 4th contributions resulted from Eq. (177) and Eq. (178)
c thus, the computation for s2fs is completed. 
                          do lg1=1,neg 
                            do l=0,isct
                              m2=noff5+neg*neg*l+neg*(ng1-1)+lg1
c Eq. (178) 
                              st=dsg(2,ng1,lg1,l,ir,jr)*blk(3,k3)
     1                         +dsg(3,ng1,lg1,l,ir,jr)*blk(3,k3)
c Eq. (177) to account for the contribution of sigs(l=0) to sigt
                               if(l.eq.0)then
                                 st=st-sum(dtg(lg1,2:3,ir,jr))
     1                            *blk(3,k3)
                               end if     
                               s2fs(n,tj,m2)=
     1                          s2fs(n,tj,m2)+st*sigf(ng,iso)
     2                          *sigs(ng1,lg1,l,iso2)
                            end do ! l
                          end do ! lg1
c Eq. (179) for sigf-nu
c s2fnu comprises 2 contributions: Eq160fnu+Eq179fnu
c here is the 2nd contribution resulted from Eq. (179) to s2fnu.
c thus, the computation for s2fnu is completed.                    
                          m2=neg*(k3-1)+ng1
                          do lg1=1,neg
                            if(sigf(ng1,iso2).ne.0.0)then
                              st=sum(dfg(2:3,lg1,ng1,ir,jr))
     1                         *chivec(lg1,i2)*blk(3,k3)
     2                         *sigf(ng1,iso2)
                              if(iso2.eq.iso.and.ng1.eq.ng)then
                                st=st+dfg(1,lg1,ng1,ir,jr)
     1                           *chivec(lg1,i2)*blk(3,k3)
                              end if  
                              s2fnu(n,tj,m2)=s2fnu(n,tj,m2)
     1                         +st*sigf(ng,iso)
     2                         *nusigf(ng1,iso2)/sigf(ng1,iso2)
                            end if
                          end do ! lg1
c Eq. (179) for sigf-sigf
c s2ff comprises 4 contributions: Eq158ff+Eq160ff+Eq177ff+Eq179ff
c here is the 4th contribution resulted from Eq. (179) to s2ff,
c thus, the computation for s2ff is completed.                   
                          do lg1=1,neg
                            if(sigf(ng1,iso2).ne.0.0)then
                              st=sum(dfg(2:3,lg1,ng1,ir,jr))
     1                         *chivec(lg1,i2)*blk(3,k3)
     2                         *nusigf(ng1,iso2)/sigf(ng1,iso2)
                              s2ff(n,tj,m2)=s2ff(n,tj,m2)
     1                         +st*sigf(ng,iso)*sigf(ng1,iso2)
                            end if
                          end do ! lg1
ccc
                        end do ! ng1
                      end do ! ir
                    end do ! jr
                  end do ! j2
                end do ! i2
              end if ! okf.or.oka
            end do ! ng 
c
c sigs(l=0) to sigt
c
            do l=0,0
              do lg=1,neg 
c                
c reading 2nd-order adjoint functions psi(2,g;1,j) and psi(2,g;2,j)
c with respect to each parameter sigs(ng,lg,l=0,iso). Use "nf" to find the
c correct files to read from. Note that for l=0, the adjoint functions
c computed using Eqs. (164) and (165) with respect to sigs(ng,lg,l=0,iso)
c are the same as the adjoint functions computed using Eqs. (164) and
c (165) with respect to the parameter sigt(lg,iso), since the sources for
c the 2nd-LASS are the same.
c
c For the sources used to solve for psi(2,g;1,j) and psi(2,g;2,j) with
c respect to sigt, refer to Eqs. (37) and (38) in Ref. 2.
c
c For the sources used to solve for psi(2,g;1,j) and psi(2,g;2,j) with
c respect to sigs(l=0), refer to Eqs. (36) and (37) in Ref. 3.
c 
c Therefore, no need to run forward/adjoint PARTISN to solve 2nd-LASS for
c the adjoint functions psi(2,g;1,j) and psi(2,g;2,j) for sigs(ng,lg,l=0,tj).
c Instead, the 2nd-order adjoint functions psi(2,g;1,j) and psi(2,g;2,j) 
c obtained for the total cross sections can be direclty used. 
c 
c The correspondence between the parameter sigs(ng,lg,l=0,iso) and the
c parameter sigt(lg,iso) that having the same adjoint functions are as
c follows:
c  sigs(ng=1,lg=1,l=0,iso=3)  corresponds to sigt(lg=1,iso=3);
c  sigs(ng=1,lg=2,l=0,iso=3)  corresponds to sigt(lg=2,iso=3); 
c  sigs(ng=1,lg=30,l=0,iso=3) corresponds to sigt(lg=30,iso=3); 
c  sigs(ng=2,lg=1,l=0,iso=3)  corresponds to sigt(lg=1,iso=3); 
c  sigs(ng=2,lg=2,l=0,iso=3)  corresponds to sigt(lg=2,iso=3); 
c  sigs(ng=2,lg=30,l=0,iso=3) corresponds to sigt(lg=30,iso=3);
c and so on so forth. Note that the correspondence between the two
c parameters depends on lg, regardless of ng.   
c
                nf=neg*(k2-1)+lg
                write(afile,'("fo2/",i5.5,"_rmflux")')nf
                inquire(file=afile, exist=okf)
                if(okf) then 
                  call rddantm2o(1,nf,n)
                  if(iangflux.eq.1.and.l.eq.0)call rddanta2o(1,nf,n)
                endif
                write(afile,'("ad2/",i5.5,"_amflux")')nf
                inquire(file=afile, exist=oka)
                if (oka) then
                  call rddantm2o(2,nf,n)
                  if(iangflux.eq.1.and.l.eq.0)call rddanta2o(2,nf,n)
                endif
                if (okf.or.oka)then  
                  call integrals_2nd(0,rlip,fmom,amom,
     1              afreg,afadj,fmom2,amom2,afreg2,afadj2,neg0,ndir0,
     2              nitm0,njtm0,dtg,dsg,dfg)
c loop over each 0th order scattering cross section parameter
c sigs(ng,lg,l=0,iso), indexed by tj.
                  do ng=1,neg 
                    if(sigs(ng,lg,l,iso).ne.0.d0)then
                      noff5=neg*neg*nmom*(k2-1)
                      tj=noff5+neg*neg*l+neg*(ng-1)+lg 
c to compute the mixed second derivatives, for each sigs(ng,lg,l=0,iso)
c loop over m2 
                      k3=0
                      do i2=1,nm
                        do j2=1,ncb(i2)
c for testing
!                      do i2=1,1
!                        do j2=1,1
                          k3=k3+1
                          iso2=nm+k3
                          noff6=neg*neg*nmom*(k3-1)
                          noff7=neg*neg*(k3-1)
                          do ir=1,nr
                            do jr=1,nz
                              if(mat(ir,jr).ne.i2)cycle
                              do ng1=1,neg
                                tj=noff5+neg*neg*l+neg*(ng-1)+lg
                                m2=neg*(k3-1)+ng1 
c Eq. (158) for sigs-sigt                  
c s2st comprises 2 contributions: Eq158st+Eq167st
c here is the 1st contribution resulted from Eq. (158) to s2st, 
c due to sigs(l=0) is part of sigt, where m2=1,...,Jsigf 
                                s2st(n,tj,m2)=
     1                           s2st(n,tj,m2)
     2                           -sum(dtg(ng1,2:3,ir,jr))*blk(3,k3)
     3                           *sigs(ng,lg,l,iso)*sigt(ng1,iso2)
c Eq. (158) for sigs-sigf 
c s2sf comprises 4 contributions: Eq158sf+Eq160sf+Eq167sf+Eq1s9sf
c here is the 1st contribution resulted from Eq. (158) to s2sf,
c due to sigs(l=0) is part of sigt, where m2=1,...,Jsigf 
                                s2sf(n,tj,m2)=
     1                           s2sf(n,tj,m2)
     2                            -sum(dtg(ng1,2:3,ir,jr))*blk(3,k3)
     3                            *sigs(ng,lg,l,iso)*sigf(ng1,iso2)
c Eq. (160) for sigs-sigf
c s2sf comprises 4 contributions: Eq158sf+Eq160sf+Eq167sf+Eq1s9sf
c here is the 2nd contribution resulted from Eq. (160) to s2sf,
c due to both sigs(l=0) and sigf are parts of sigt                    
                                do lg1=1,neg
                                  if(sigf(ng1,iso2).ne.0.0)then
                                    st=sum(dfg(2:3,lg1,ng1,ir,jr))
     1                              *chivec(lg1,i2)*blk(3,k3)
     2                               *nusigf(ng1,iso2)/sigf(ng1,iso2)
                                    s2sf(n,tj,m2)=
     1                               s2sf(n,tj,m2)+st
     2                               *sigs(ng,lg,l,iso)*sigf(ng1,iso2)
                                  end if
                                end do ! lg1
c Eq. (160) for sigs-nu 
c s2snu comprises 2 contributions: Eq160snu+Eq169snu
c here is the 1st contribution resulted from Eq. (160) to s2snu, 
c due to sigs(l=0) is part of sigt
                                do lg1=1,neg
                                  if(sigf(ng1,iso2).ne.0.0)then
                                    st=sum(dfg(2:3,lg1,ng1,ir,jr))
     1                               *chivec(lg1,i2)*blk(3,k3)
     2                               *sigf(ng1,iso2)
                                    s2snu(n,tj,m2)=
     1                               s2snu(n,tj,m2)
     2                               +st*sigs(ng,lg,l,iso)
     3                               *nusigf(ng1,iso2)/sigf(ng1,iso2)
                                  end if
                                end do ! lg1
c Eq. (158) for sigs-sigs
c s2ss comprises 4 contributions: Eq158ss+Eq159ss+Eq167ss+Eq168ss
c here is the 1st contribution resulted from Eq. (158) to s2ss,
c due to the contribution of both sigs_j(l=0) and sigs_m2(l=0) to sigt
                                do lg1=1,neg
                                  do l1=0,0
c due to s2ss(0:nrrr,neg*neg,neg*neg*nmom*nel,0:isct), tj has to be
c within the size of neg*neg  
                                    tj=neg*(ng-1)+lg 
                                    m2=noff6+neg*neg*l1+neg*(ng1-1)+lg1
                                    s2ss(n,tj,m2,l)=
     1                               s2ss(n,tj,m2,l)
     2                               -sum(dtg(lg1,2:3,ir,jr))*blk(3,k3)
     3                               *sigs(ng,lg,l,iso)
     4                               *sigs(ng1,lg1,l1,iso2)
                                  end do ! l1
                                end do ! lg1                                     
c Eq. (159) for sigs_j-sigs_m2
c s2ss comprises 4 contributions: Eq158ss+Eq159ss+Eq167ss+Eq168ss
c here is the 2nd contribution resulted from Eq. (158) to s2ss,
c due to the contribution from sigs_j(l=0) to sigt
                                do lg1=1,neg 
                                  do l1=0,isct
                                    tj=neg*(ng-1)+lg 
                                    m2=noff6+neg*neg*l1+neg*(ng1-1)+lg1
                                    st=dsg(2,ng1,lg1,l1,ir,jr)*blk(3,k3)
     1                               +dsg(3,ng1,lg1,l1,ir,jr)*blk(3,k3)
                                    s2ss(n,tj,m2,l)=
     1                               s2ss(n,tj,m2,l)+st
     2                               *sigs(ng,lg,l,iso)
     3                               *sigs(ng1,lg1,l1,iso2)
                                  end do ! l1
                                end do ! lg1
ccc
                              end do ! ng1
                            end do ! jr
                          end do ! ir
                        end do ! j2
                      end do ! i2
                    end if ! sigs.ne.0.d0
                  end do ! ng 
                end if ! okf.or.oka
              end do ! lg                
            end do ! l
c
c sigs
c
c reading 2nd-order adjoint functions theta(2,g;1,j) and theta(2,g;2,j)
c computed from Eqs. (173) and (174), respectively. Use "nf" to find the
c correct files to read from. 
            do l=0,isct
c for testing
!            do l=0,0
              do ng=1,neg 
                do lg=1,neg 
c loop for each parameter in vector sigs=[...sigs(ng,lg,l,i)...], indexed by
c tj.
                  noff5=neg*neg*nmom*(k2-1)
                  tj=noff5+neg*neg*l+neg*(ng-1)+lg                
c                
                  if(sigs(ng,lg,l,iso).ne.0.d0)then
                    nf=3*noff2+tj
                    write(afile,'("fo2/",i5.5,"_rmflux")')nf
                    inquire(file=afile, exist=okf)
                    if(okf) then 
                      call rddantm2o(1,nf,n)
                      if(iangflux.eq.1)call rddanta2o(1,nf,n)
                    endif
                    write(afile,'("ad2/",i5.5,"_amflux")')nf
                    inquire(file=afile, exist=oka)
                    if (oka) then
                      call rddantm2o(2,nf,n)
                      if(iangflux.eq.1)call rddanta2o(2,nf,n)
                    endif
                    if (okf.or.oka)then  
                      call integrals_2nd(0,rlip,fmom,amom,
     1                 afreg,afadj,fmom2,amom2,afreg2,afadj2,neg0,ndir0,
     2                 nitm0,njtm0,dtg,dsg,dfg)
c to compute the mixed second derivatives: keep s_tj unchanged,
c loop over m2 
                      k3=0
                      do i2=1,nm
                        do j2=1,ncb(i2)
c for testing
!                      do i2=1,1
!                        do j2=1,1
                          k3=k3+1
                          iso2=nm+k3
                          noff6=neg*neg*nmom*(k3-1)
                          noff7=neg*neg*(k3-1)
                          do ir=1,nr
                            do jr=1,nz
                              if(mat(ir,jr).ne.i2)cycle
                              do ng1=1,neg
                                tj=noff5+neg*neg*l+neg*(ng-1)+lg
                                m2=neg*(k3-1)+ng1 
c Eq. (167) for sigs-sigt
c s2st comprises 2 contributions: Eq158st+Eq167st
c here is the 2nd contribution resulted from Eq. (167) to s2st,
c thus, the computation for s2st is completed. 
                                s2st(n,tj,m2)=
     1                           s2st(n,tj,m2)
     2                           -sum(dtg(ng1,2:3,ir,jr))*blk(3,k3)
     3                           *sigs(ng,lg,l,iso)*sigt(ng1,iso2)
c Eq. (167) for sigs-sigf
c s2sf comprises 4 contributions: Eq158sf+Eq160sf+Eq167sf+Eq1s9sf
c here is the 3rd contribution resulted from Eq. (167) to s2sf                  
                                s2sf(n,tj,m2)=
     1                           s2sf(n,tj,m2)
     2                           -sum(dtg(ng1,2:3,ir,jr))*blk(3,k3)
     3                           *sigs(ng,lg,l,iso)*sigf(ng1,iso2)
c Eq. (169) for sigs-sigf                   
c s2sf comprises 4 contributions: Eq158sf+Eq160sf+Eq167sf+Eq1s9sf
c here is the 4th contribution resulted from Eq. (169) to s2sf,
c thus, the computation for s2sf is completed.
                                do lg1=1,neg
                                  if(sigf(ng1,iso2).ne.0.0)then
                                    st=sum(dfg(2:3,lg1,ng1,ir,jr))
     1                               *chivec(lg1,i2)*blk(3,k3)
     2                               *nusigf(ng1,iso2)/sigf(ng1,iso2)
                                    s2sf(n,tj,m2)=
     1                               s2sf(n,tj,m2)+st
     2                               *sigs(ng,lg,l,iso)*sigf(ng1,iso2)
                                  end if
                                end do ! lg1
c Eq. (169) for sigs-nu
c s2snu comprises 2 contributions: Eq160snu+Eq169snu
c here is the 2nd contribution resulted from Eq. (169) to s2snu,
c thus, the computation for s2snu is completed.                  
                                do lg1=1,neg
                                  if(sigf(ng1,iso2).ne.0.0)then
                                    st=sum(dfg(2:3,lg1,ng1,ir,jr))
     1                               *chivec(lg1,i2)*blk(3,k3)
     2                               *sigf(ng1,iso2)
                                    s2snu(n,tj,m2)=
     1                               s2snu(n,tj,m2)
     2                               +st*sigs(ng,lg,l,iso)
     3                               *nusigf(ng1,iso2)/sigf(ng1,iso2)
                                  end if
                                end do ! lg1
c Eq. (167) for sigs-sigs
c s2ss comprises 4 contributions: Eq158ss+Eq159ss+Eq167ss+Eq168ss
c here is the 3rd contribution resulted from Eq. (167) to s2ss,
c due to the contribution of sigs_m2(l=0) to sigt
                                do lg1=1,neg
                                  do l1=0,0
                                    tj=neg*(ng-1)+lg
                                    m2=noff6+neg*neg*l1+neg*(ng1-1)+lg1
                                    s2ss(n,tj,m2,l)=
     1                               s2ss(n,tj,m2,l)
     2                               -sum(dtg(lg1,2:3,ir,jr))*blk(3,k3)
     3                               *sigs(ng,lg,l,iso)
     4                               *sigs(ng1,lg1,l1,iso2)
                                  end do ! l1
                                end do ! lg1 
c Eq. (168) for sigs_j-sigs_m2
c s2ss comprises 4 contributions: Eq158ss+Eq159ss+Eq167ss+Eq168ss
c here is the 4th contribution resulted from Eq. (168) to s2ss,
c thus, the computation for s2ss is completed.
                                do lg1=1,neg 
                                  do l1=0,isct
                                    tj=neg*(ng-1)+lg
                                    m2=noff6+neg*neg*l1+neg*(ng1-1)+lg1
                                    st=dsg(2,ng1,lg1,l1,ir,jr)*blk(3,k3)
     1                               +dsg(3,ng1,lg1,l1,ir,jr)*blk(3,k3)
                                    s2ss(n,tj,m2,l)=
     1                               s2ss(n,tj,m2,l)+st
     2                               *sigs(ng,lg,l,iso)
     3                               *sigs(ng1,lg1,l1,iso2)
                                  end do ! l1
                                end do ! lg1
c
                              end do ! ng1
                            end do ! jr
                          end do ! ir
                        end do ! j2
                      end do ! i2
                    end if ! okf.or.oka
                  end if ! sigs.ne.0.d0
                end do ! lg                
              end do ! ng 
            end do ! l
c due to the extremely large size of the 2nd-order sensitivity matrix 
c for s2ss, which is Jsigs x Jsigs, here only a portion of the full matrix
c is output.   
c Output 2nd-order mixed sensitivities with respect to the same
c scattering cross sections into a matrix named s2ssdiag, which is
c comprised by the diagonal elements of the full matrix s2ss 
            do l=0,isct
              noff8=neg*neg*nmom*(k2-1)+neg*neg*l
              do ng=1,neg
                do lg=1,neg
                  tj=neg*(ng-1)+lg
                  s2ssdiag(n,ng,lg,l,k2)=s2ss(n,tj,noff8+tj,l)  
                end do ! lg 
              end do ! ng
            end do ! l 
c re-initialize the matrix for next isotope
            s2ss(0:nrrr,1:neg*neg,1:neg*neg*nmom*nel,0:isct)=0.d0
c
c end do for blocks of sigt,nu,sigf,sigs
          end do ! j1
        end do ! i1

c for leakage sensitivities, divide by leakage
        if(n.eq.0)then
          fissd=1.d0/lkg
        end if    
        if(n.eq.0)then
          s2tt(n,1:neg*nel,1:neg*nel)=
     1     s2tt(n,1:neg*nel,1:neg*nel)*fissd
          s2ts(n,1:neg*nel,1:neg*neg*nmom*nel)=
     1     s2ts(n,1:neg*nel,1:neg*neg*nmom*nel)*fissd
          s2ft(n,1:neg*nel,1:neg*nel)=
     1     s2ft(n,1:neg*nel,1:neg*nel)*fissd
          s2nut(n,1:neg*nel,1:neg*nel)=
     1     s2nut(n,1:neg*nel,1:neg*nel)*fissd
          s2tf(n,1:neg*nel,1:neg*nel)=
     1     s2tf(n,1:neg*nel,1:neg*nel)*fissd
          s2tnu(n,1:neg*nel,1:neg*nel)=
     1     s2tnu(n,1:neg*nel,1:neg*nel)*fissd
          s2st(n,1:neg*neg*nmom*nel,1:neg*nel)=
     1     s2st(n,1:neg*neg*nmom*nel,1:neg*nel)*fissd
          s2ssdiag(n,1:neg,1:neg,0:isct,1:nel2)=
     1     s2ssdiag(n,1:neg,1:neg,0:isct,1:nel2)*fissd
          s2nus(n,1:neg*nel,1:neg*neg*nmom*nel)=
     1     s2nus(n,1:neg*nel,1:neg*neg*nmom*nel)*fissd
          s2fs(n,1:neg*nel,1:neg*neg*nmom*nel)=
     1     s2fs(n,1:neg*nel,1:neg*neg*nmom*nel)*fissd
          s2sf(n,1:neg*neg*nmom*nel,1:neg*nel)=
     1     s2sf(n,1:neg*neg*nmom*nel,1:neg*nel)*fissd
          s2snu(n,1:neg*neg*nmom*nel,1:neg*nel)=
     1     s2snu(n,1:neg*neg*nmom*nel,1:neg*nel)*fissd
          s2nunu(n,1:neg*nel,1:neg*nel)=
     1     s2nunu(n,1:neg*nel,1:neg*nel)*fissd
          s2ff(n,1:neg*nel,1:neg*nel)=
     1     s2ff(n,1:neg*nel,1:neg*nel)*fissd
          s2nuf(n,1:neg*nel,1:neg*nel)=
     1     s2nuf(n,1:neg*nel,1:neg*nel)*fissd
          s2fnu(n,1:neg*nel,1:neg*nel)=
     1     s2fnu(n,1:neg*nel,1:neg*nel)*fissd
        end if
c
c output
        write(*,'(/,"writing result arrays...")')
        if(n.eq.0)then
          if(icalc.eq.0)then
            write(ius2tt,'(/," leakage sensitivities",10x,"leakage=",
     1       1pe15.7)')lkg
            write(ius2tf,'(/," leakage sensitivities",10x,"leakage=",
     1       1pe15.7)')lkg
            write(ius2tnu,'(/," leakage sensitivities",10x,"leakage=",
     1       1pe15.7)')lkg
            write(ius2ts,'(/," leakage sensitivities",10x,"leakage=",
     1       1pe15.7)')lkg
            write(ius2nut,'(/," leakage sensitivities",10x,"leakage=",
     1       1pe15.7)')lkg
            write(ius2nunu,'(/," leakage sensitivities",10x,"leakage=",
     1       1pe15.7)')lkg
            write(ius2nuf,'(/," leakage sensitivities",10x,"leakage=",
     1       1pe15.7)')lkg
            write(ius2nus,'(/," leakage sensitivities",10x,"leakage=",
     1       1pe15.7)')lkg
            write(ius2ft,'(/," leakage sensitivities",10x,"leakage=",
     1       1pe15.7)')lkg
            write(ius2fs,'(/," leakage sensitivities",10x,"leakage=",
     1       1pe15.7)')lkg
            write(ius2fnu,'(/," leakage sensitivities",10x,"leakage=",
     1       1pe15.7)')lkg
            write(ius2ff,'(/," leakage sensitivities",10x,"leakage=",
     1       1pe15.7)')lkg
            write(ius2st,'(/," leakage sensitivities",10x,"leakage=",
     1       1pe15.7)')lkg
            write(ius2sf,'(/," leakage sensitivities",10x,"leakage=",
     1       1pe15.7)')lkg
            write(ius2snu,'(/," leakage sensitivities",10x,"leakage=",
     1       1pe15.7)')lkg
            write(ius2ss,'(/," leakage sensitivities",10x,"leakage=",
     1       1pe15.7)')lkg
          end if
        end if
c k2 and m are used to keep track of materials for writing the header.
        k2=0
        m=1
        do i=1,nel
          k2=k2+1
          write(ius2tt,'(/,1x,a24," isotope",i4," in material",i4)')
     1     zaidfull(i),k2,m  
          write(ius2tf,'(/,1x,a24," isotope",i4," in material",i4)')
     1     zaidfull(i),k2,m 
          write(ius2tnu,'(/,1x,a24," isotope",i4," in material",i4)')
     1     zaidfull(i),k2,m 
          write(ius2ts,'(/,1x,a24," isotope",i4," in material",i4)')
     1     zaidfull(i),k2,m 
          write(ius2nut,'(/,1x,a24," isotope",i4," in material",i4)')
     1     zaidfull(i),k2,m 
          write(ius2nunu,'(/,1x,a24," isotope",i4," in material",i4)')
     1     zaidfull(i),k2,m 
          write(ius2nuf,'(/,1x,a24," isotope",i4," in material",i4)')
     1     zaidfull(i),k2,m 
          write(ius2nus,'(/,1x,a24," isotope",i4," in material",i4)')
     1     zaidfull(i),k2,m 
          write(ius2ft,'(/,1x,a24," isotope",i4," in material",i4)')
     1     zaidfull(i),k2,m
          write(ius2fs,'(/,1x,a24," isotope",i4," in material",i4)')
     1     zaidfull(i),k2,m 
          write(ius2fnu,'(/,1x,a24," isotope",i4," in material",i4)')
     1     zaidfull(i),k2,m 
          write(ius2ff,'(/,1x,a24," isotope",i4," in material",i4)')
     1     zaidfull(i),k2,m 
c write mixed 2nd-order sensitivities matrix for sigt-sigt
          write(ius2tt,'("2nd-order mixed sensitivities w.r.to "
     1     "sigt_j of isotope",i2," in material",i2," and sigt_m2 for"
     2     "  all isotopes (rows are for sigt_j, columns are for "
     3     "sigt_m2)")')k2,m
          write(ius2tt,'(9999i14)')(ng,ng=1,neg*nel)
          do ng=1,neg
            tj=neg*(i-1)+ng
            write(ius2tt,'(i6,1p9999e14.6)')tj,
     1        (s2tt(n,tj,m2),m2=1,neg*nel)
          end do ! ng
c write mixed 2nd-order sensitivities matrix for sigt-sigf
          write(ius2tf,'("2nd-order mixed sensitivities w.r.to "
     1     "sigt_j of isotope",i2," in material",i2," and sigf_m2 for"
     2     " all isotopes (rows are for sigt_j, columns are for "
     3     "sigf_m2)")')k2,m
          write(ius2tf,'(9999i14)')(ng,ng=1,neg*nel)
          do ng=1,neg
            tj=neg*(i-1)+ng
            write(ius2tf,'(i6,1p9999e14.6)')tj,
     1        (s2tf(n,tj,m2),m2=1,neg*nel)
          end do ! ng          
c write mixed 2nd-order sensitivities matrix for sigt-nu
          write(ius2tnu,'("2nd-order mixed sensitivities w.r.to "
     1     "sigt_j of isotope",i2," in material",i2," and nu_m2 for all"
     2     " isotopes (rows are for sigt_j, columns are for "
     3     "nu_m2)")')k2,m
          write(ius2tnu,'(9999i14)')(ng,ng=1,neg*nel)
          do ng=1,neg
            tj=neg*(i-1)+ng
            write(ius2tnu,'(i6,1p9999e14.6)')tj,
     1        (s2tnu(n,tj,m2),m2=1,neg*nel)
          end do ! ng
c write mixed 2nd-order sensitivities matrix for sigt-sigs
          k4=0
          m4=1
          do j=1,nel
            k4=k4+1
            do l=0,isct
              noff2=neg*neg*nmom*(j-1)+neg*neg*l
              write(ius2ts,'("2nd-order mixed sensitivities w.r.to "
     1         "sigt_j of isotope",i2," in material",i2," and sigs_m2"
     2         " of isotope",i2," in material",i2,", order=",i2," (rows"
     3         " are sigt_j,columns are sigs_m2)")')k2,m,k4,m4,l
              write(ius2ts,'(9999i14)')(lg,lg=noff2+1,noff2+neg*neg)
              do ng=1,neg
                tj=neg*(i-1)+ng
                write(ius2ts,'(i6,1p9999e14.6)')tj,
     1          (s2ts(n,tj,m2),m2=noff2+1,noff2+neg*neg)
              end do ! ng
            end do ! l
            if(k4.eq.ncb(m4))then
              k4=0
              m4=m4+1
            end if
          end do ! j
c write mixed 2nd-order sensitivities matrix for nu-sigt
          write(ius2nut,'("2nd-order mixed sensitivities w.r.to "
     1     "nu_j of isotope",i2," in material",i2," and sigt_m2 for"
     2     " all isotopes (rows are for nu_j, columns are for "
     3     "sigt_m2)")')k2,m
          write(ius2nut,'(9999i14)')(ng,ng=1,neg*nel)
          do ng=1,neg
            tj=neg*(i-1)+ng
            write(ius2nut,'(i6,1p9999e14.6)')tj,
     1        (s2nut(n,tj,m2),m2=1,neg*nel)
          end do ! ng 
c write mixed 2nd-order sensitivities matrix for nu-nu
          write(ius2nunu,'("2nd-order mixed sensitivities w.r.to "
     1     "nu_j of isotope",i2," in material",i2," and nu_m2 for all"
     2     " isotopes (rows are for nu_j, columns are for "
     3     "nu_m2)")')k2,m
          write(ius2nunu,'(9999i14)')(ng,ng=1,neg*nel)
          do ng=1,neg
            tj=neg*(i-1)+ng
            write(ius2nunu,'(i6,1p9999e14.6)')tj,
     1        (s2nunu(n,tj,m2),m2=1,neg*nel)
          end do ! ng         
c write mixed 2nd-order sensitivities matrix for nu-sigf
          write(ius2nuf,'("2nd-order mixed sensitivities w.r.to "
     1     "nu_j of isotope",i2," in material",i2," and sigf_m2 for all"
     2     " isotopes (rows are for nu_j, columns are for "
     3     "sigf_m2)")')k2,m
          write(ius2nuf,'(9999i14)')(ng,ng=1,neg*nel)
          do ng=1,neg
            tj=neg*(i-1)+ng
            write(ius2nuf,'(i6,1p9999e14.6)')tj,
     1        (s2nuf(n,tj,m2),m2=1,neg*nel)
          end do ! ng 
c write mixed 2nd-order sensitivities matrix for nu-sigs
          k4=0
          m4=1
          do j=1,nel
            k4=k4+1
            do l=0,isct
              noff2=neg*neg*nmom*(j-1)+neg*neg*l
              write(ius2nus,'("2nd-order mixed sensitivities w.r.to "
     1         "nu_j of isotope",i2," in material",i2," and sigs_m2 of "
     2         "isotope",i2," in material",i2,", order=",i2," (rows "
     3         "are for nu_j,columns are for sigs_m2)")')k2,m,k4,m4,l
              write(ius2nus,'(9999i14)')(lg,lg=noff2+1,noff2+neg*neg)
              do ng=1,neg
                tj=neg*(i-1)+ng
                write(ius2nus,'(i6,1p9999e14.6)')tj,
     1          (s2nus(n,tj,m2),m2=noff2+1,noff2+neg*neg)
              end do ! ng
            end do ! l
            if(k4.eq.ncb(m4))then
              k4=0
              m4=m4+1
            end if
          end do ! j
c write mixed 2nd-order sensitivities matrix for sigf-sigt
          write(ius2ft,'("2nd-order mixed sensitivities w.r.to "
     1     "sigf_j of isotope",i2," in material",i2," and sigt_m2 for"
     2     " all isotopes (rows are for sigf_j, columns are for "
     3     "sigt_m2)")')k2,m
          write(ius2ft,'(9999i14)')(ng,ng=1,neg*nel)
          do ng=1,neg
            tj=neg*(i-1)+ng
            write(ius2ft,'(i6,1p9999e14.6)')tj,
     1        (s2ft(n,tj,m2),m2=1,neg*nel)
          end do ! ng 
c write mixed 2nd-order sensitivities matrix for sigf-sigs
          k4=0
          m4=1
          do j=1,nel
            k4=k4+1
            do l=0,isct
              noff2=neg*neg*nmom*(j-1)+neg*neg*l
              write(ius2fs,'("2nd-order mixed sensitivities w.r.to "
     1         "sigf_j of isotope",i2," in material",i2," and sigs_m2"
     2         " of isotope",i2," in material",i2,", order=",i2," (rows"
     3         " are for sigf_j,columns are for sigs_m2)")')k2,m,k4,m4,l
              write(ius2fs,'(9999i14)')(lg,lg=noff2+1,noff2+neg*neg)
              do ng=1,neg
                tj=neg*(i-1)+ng
                write(ius2fs,'(i6,1p9999e14.6)')tj,
     1          (s2fs(n,tj,m2),m2=noff2+1,noff2+neg*neg)
              end do ! ng
            end do ! l
            if(k4.eq.ncb(m4))then
              k4=0
              m4=m4+1
            end if
          end do ! j
c write mixed 2nd-order sensitivities matrix for sigf-nu
          write(ius2fnu,'("2nd-order mixed sensitivities w.r.to "
     1     "sigf_j of isotope",i2," in material",i2," and nu_m2 for all"
     2     " isotopes (rows are for sigf_j, columns are for "
     3     "nu_m2)")')k2,m
          write(ius2fnu,'(9999i14)')(ng,ng=1,neg*nel)
          do ng=1,neg
            tj=neg*(i-1)+ng
            write(ius2fnu,'(i6,1p9999e14.6)')tj,
     1        (s2fnu(n,tj,m2),m2=1,neg*nel)
          end do ! ng 
c write mixed 2nd-order sensitivities matrix for sigf-sigf
          write(ius2ff,'("2nd-order mixed sensitivities w.r.to "
     1     "sigf_j of isotope",i2," in material",i2," and sigf_m2 for"
     2     " all isotopes (rows are for sigf_j, columns are for "
     3     "sigf_m2)")')k2,m
          write(ius2ff,'(9999i14)')(ng,ng=1,neg*nel)
          do ng=1,neg
            tj=neg*(i-1)+ng
            write(ius2ff,'(i6,1p9999e14.6)')tj,
     1        (s2ff(n,tj,m2),m2=1,neg*nel)
          end do ! ng 
c write mixed 2nd-order sensitivities matrix for sigs-sigt         
          do l=0,isct
            noff2=neg*neg*nmom*(i-1)+neg*neg*l
            write(ius2st,'(/,1x,a24," isotope",i2," in material",i2)')
     1       zaidfull(i),k2,m
            write(ius2st,'("2nd-order mixed sensitivities w.r.to "
     1       "sigs_j of isotope",i2," in material",i2," and sigt_m2 for"
     2       " all isotopes, order=",i2," of sigs (rows are for sigs_j,"
     3       " columns are for sigt_m2)")')k2,m,l
            write(ius2st,'(9999i14)')(lg,lg=1,neg*nel)
            do ng=1,neg
              do lg=1,neg
                tj=noff2+neg*(ng-1)+lg
                write(ius2st,'(i6,1p9999e14.6)')tj,
     1           (s2st(n,tj,m2),m2=1,neg*nel)
              end do ! ng
            end do ! ng
          end do ! l
c write mixed 2nd-order sensitivities matrix for sigs-sigf         
          do l=0,isct
            noff2=neg*neg*nmom*(i-1)+neg*neg*l
            write(ius2sf,'(/,1x,a24," isotope",i2," in material",i2)')
     1       zaidfull(i),k2,m
            write(ius2sf,'("2nd-order mixed sensitivities w.r.to "
     1       "sigs_j of isotope",i2," in material",i2," and sigf_m2 for"
     2       " all isotopes, order=",i2," of sigs_j (rows are for"
     3       " sigs_j, columns are for sigf_m2)")')k2,m,l
            write(ius2sf,'(9999i14)')(lg,lg=1,neg*nel)
            do ng=1,neg
              do lg=1,neg
                tj=noff2+neg*(ng-1)+lg
                write(ius2sf,'(i6,1p9999e14.6)')tj,
     1           (s2sf(n,tj,m2),m2=1,neg*nel)
              end do ! ng
            end do ! ng
          end do ! l
c write mixed 2nd-order sensitivities matrix for sigs-nu         
          do l=0,isct
            noff2=neg*neg*nmom*(i-1)+neg*neg*l
            write(ius2snu,'(/,1x,a24," isotope",i2," in material",i2)')
     1       zaidfull(i),k2,m
            write(ius2snu,'("2nd-order mixed sensitivities w.r.to "
     1       "sigs_j of isotope",i2," in material",i2," and nu_m2 for"
     2       " all isotopes, order=",i2," of sigs_j (rows are for "
     3       "sigs_j, columns are for nu_m2)")')k2,m,l
            write(ius2snu,'(9999i14)')(lg,lg=1,neg*nel)
            do ng=1,neg
              do lg=1,neg
                tj=noff2+neg*(ng-1)+lg
                write(ius2snu,'(i6,1p9999e14.6)')tj,
     1           (s2snu(n,tj,m2),m2=1,neg*nel)
              end do ! ng
            end do ! ng
          end do ! l
c write unmixed 2nd-order sensitivities matrix with respect to the same
c scattering cross sections, which are the diagonal elements of the full
c matrix s2ss
          do l=0,isct
            write(ius2ss,'(/,1x,a24," isotope",i2," in material",i2)')
     1       zaidfull(i),k2,m
            write(ius2ss,'("2nd-order unmixed sensitivities w.r.to the"
     1       " same scattering cross sections of isotope",i2," in "
     2       "material",i2,", order=",i2," (rows are for ng, columns "
     3       "are for lg)")')k2,m,l
            write(ius2ss,'(9999i14)')(lg,lg=1,neg)
            do ng=1,neg
              write(ius2ss,'(i6,1p9999e14.6)')ng,
     1         (s2ssdiag(n,ng,lg,l,i),lg=1,neg)
            end do ! ng
          end do ! l
ccc
          if(k2.eq.ncb(m))then
            k2=0
            m=m+1
          end if
        end do ! i
c
      end do ! n
c
      if(nrrr.eq.0)then
        write(*,'("no sens_rr file because there are no reaction ",
     1   "rates specified.")')
        write(iuo,'("no sens_rr file because there are no reaction ",
     1   "rates specified.")')
      end if
c
      return
      end subroutine
