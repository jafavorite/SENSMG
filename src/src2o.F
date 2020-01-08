      subroutine src2o(ncalcs,niso0)
      use COMS
c construct sources for 2nd-level adjoint sensitivity system (2nd-LASS)
c for the lakage response with respect to various model parameters
c
c the sources were based on equations in
c Ref. 1: D. G. Cacuci, "Application of the Second-order Comprehensive
c Adjoint Sensitivity Analysis Methodology to Compute 1st- and 2nd-Order
c Sensitivities of Flux Functionals in a Multiplying System With
c Source," Nucl. Sci. Eng., 193, 555-600 (2019); 
chttps://doi.org/10.1080/00295639.2018.1553910.
c
c expressions for the source terms derived from the 2nd-LASS equations
c in Ref. 1 can be found from the following references:
c
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
c This code is for computing 2nd-ASAM sources with respect to sigt, nu,
c sigf, and sigs. (2nd-ASAM sources with respect to isotopic number 
c density N(i) are also available, will be integrated into this 
c subroutine as part of the delivery).
c 
c The parameters in this subroutine are organized as follows:   
c       sigt
c       nu
c       sigf
c       sigs
c       N(i) 
c
c The definitions of parameters can be found in Ref. 2.
c
c for sigt - the vector sigt=[t_1,...,t_Jsigt] is defined by Eq. (13)
c            of Ref. 2, where Jsigt=neg*nel. Specifically, the parameters
c            sigt(ng,i), ng=1,...,neg;i=1,...,nel, are organized as:
c
c            sigt=[sigt(ng=1,  i=1),   !
c                  sigt(ng=2,  i=1),   !for ng=1,...,neg; i=1 
c                  ...                 !
c                  sigt(ng=neg,i=1),   !
c                  sigt(ng=1,  i=2),   !
c                  sigt(ng=2,  i=2),   !for ng=1,...,neg; i=2
c                  ...                 !
c                  sigt(ng=neg,i=2),   !
c                  ...                 
c                  sigt(ng,    i),     
c                  ...                 
c                  sigt(ng=1,  i=nel), !
c                  sigt(ng=2,  i=nel), !for ng=1,...,neg; i=nel
c                  ...                 !
c                  sigt(ng=neg,i=nel)] !
c
c for nu   - the vector nu=[nu_1,...,nu_Jnu] is defined by Eq. (18) of
c            Ref. 2, where Jnu=neg*nel. Note: the vector nu includes all  
c            the parameters nu(ng,i), for ng=1,...,neg, and for all isotopes
c            i=1,...,nel. The values of nu  will be zero for those 
c            non-fission isotopes. This is necessary to reserve at most 
c            spaces, so that it works even if all isotopes are fissionable.     
c            
c            are organized as:
c              nu=[nu(ng=1,  i=1),   !
c                  nu(ng=2,  i=1),   !for ng=1,...,neg; i=1 
c                  ...               !
c                  nu(ng=neg,i=1),   !
c                  nu(ng=1,  i=2),   !
c                  nu(ng=2,  i=2),   !for ng=1,...,neg; i=2
c                  ...               !
c                  nu(ng=neg,i=2),   !
c                  ...               
c                  nu(ng,    i),     
c                  ...               
c                  nu(ng=1,  i=nel), !
c                  nu(ng=2,  i=nel), !for ng=1,...,neg; i=nel
c                  ...               !
c                  nu(ng=neg,i=nel)] !
c
c for sigf - the vector sigf=[f_1,...,f_Jsigf] is defined by Eq. (17)
c            of Ref. 2, where Jsigf=neg*nel. Specifically, the parameters
c            sigf(ng,i), ng=1,...,neg;i=1,...,nel, are organized as:
c
c            sigf=[sigf(ng=1,  i=1),   !
c                  sigf(ng=2,  i=1),   !for ng=1,...,neg; i=1 
c                  ...                 !
c                  sigf(ng=neg,i=1),   !
c                  sigf(ng=1,  i=2),   !
c                  sigf(ng=2,  i=2),   !for ng=1,...,neg; i=2
c                  ...                 !
c                  sigf(ng=neg,i=2),   !
c                  ...                 
c                  sigf(ng,    i),     
c                  ...                 
c                  sigf(ng=1,  i=nel), !
c                  sigf(ng=2,  i=nel), !for ng=1,...,neg; i=nel
c                  ...                 !
c                  sigf(ng=neg,i=nel)] !
c
c for sigs - the vector sigs=[s_1,...,s_Jsigs] is defined by Eq. (15)
c            in Ref. 2, where Js=neg*neg*nmom*nel. The parameters
c            sigs(ng,lg,l,i) for lg=1,neg;ng=1,neg;l=0,isct;i=1,nel,
c            are organized as:
c
c           sigs=[sigs(ng=1,lg=1,l=0,i=1),     !                 
c                 sigs(ng=1,lg=2,l=0,i=1),     !lg=1,neg;ng=1;l=0;i=1 
c                 ...                          !
c                 sigs(ng=1,lg=neg,l=0,i=1),   !
c                 sigs(ng=2,lg=1,l=0,i=1),     ! 
c                 sigs(ng=2,lg=2,l=0,i=1),     !lg=1,neg;ng=2;l=0;i=1
c                 ...                          !
c                 sigs(ng=2,lg=neg,l=0,i=1),   !
c                 ...
c                 sigs(ng=neg,lg=1,l=0,i=1),   ! 
c                 sigs(ng=neg,lg=2,l=0,i=1),   !lg=1,neg;ng=neg;l=0;i=1
c                 ...                          !
c                 sigs(ng=1,lg=1,l=1,i=1),     !
c                 sigs(ng=1,lg=2,l=1,i=1),     !lg=1,neg;ng=1;l=1;i=1
c                 ...                          !
c                 sigs(ng=1,lg=neg,l=1,i=1),   !
c                 sigs(ng=2,lg=1,l=1,i=1),     !
c                 sigs(ng=2,lg=2,l=1,i=1),     !lg=1,neg;ng=2;l=1;i=1
c                 ...                          !
c                 sigs(ng=2,lg=neg,l=1,i=1),   !
c                 ...
c                 sigs(ng,lg,l,i),      
c                 ...
c                 sigs(ng=neg,lg=1,l=isct,i=nel),  !
c                 sigs(ng=eng,lg=2,l=isct,i=nel),  !lg=1,neg;ng=neg;
c                 ...                              ! l=isct;i=nel
c                 sigs(ng=eng,lg=neg,l=isct,i=bel) !
c                 ]    
c
c globally, the assembled parameter table are ordered as follows:
c                 [[sigt],[nu],[sigf],[sigs],[N]]
c
c in this subroutine, the parameters are indexed by "k3". 
c
      implicit none
      integer ncalcs ! =nrrr+1
      integer niso0  ! number of isotopes in snxedt
      real*8 neusf(3,neg,0:nmom-1,it,jt),neusa(3,neg,0:nmom-1,it,jt),
     1 neusf_scat(neg,0:nmom-1,it,jt),neusa_scat(neg,0:nmom-1,it,jt), 
     2 st,rlip(0:isct)
      integer iso,ir,jr,i1,i,j1,j,l,ll,k,k2,k3,n,ng,lg,
     1 jf(neg*neg*nmom*nel+3*neg*nel),ja(neg*neg*nmom*nel+3*neg*nel),
     2 noff1,noff2,noff3,noff4,noff5

c this section is to compute partisn sources.
      do l=0,isct
        rlip(l)=real(2*l+1)
      end do ! l
c
c compute moments expansion of flux (or use moments from file)
c construct source using moments expansion of flux
c
c noff1 = 3*neg*nel+neg*neg*nmom*nel
c       defined the "at most" partisn forward or adjoint inputs
c       files (note that input files for N(i) are not counted
c       yet). so far, it consists the following:
c       neg*nel          - for sigt
c       neg*nel          - for nu
c       neg*nel          - for sigf
c       neg*neg*nmom*nel - for sigs  
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
      noff1=3*neg*nel+neg*neg*nmom*nel
      noff2=neg*nel
      jf(1:noff1)=0
      ja(1:noff1)=0
      neusf_scat(1:neg,0:nmom-1,1:it,1:jt)=0.d0
      neusa_scat(1:neg,0:nmom-1,1:it,1:jt)=0.d0
      write(*,'(/,"writing partisn inputs for second order. ",
     1 "there are at most",i7," forward files and",i7," adjoint files.",
     2 " nel=",i3)')noff1,noff1,nel
      write(iuo,'(/,"writing partisn inputs for second order. ",
     1 "there are at most",i7," forward files and",i7," adjoint files.",
     2 " nel=",i3)')noff1,noff1,nel
      if(noff1.gt.99999)then
        write(*,'("error. max number of inputs is 99999. increase ",
     1   "i5.5 format statements in wrdant2o.")')
        write(iuo,'("error. max number of inputs is 99999. increase ",
     1   "i5.5 format statements in wrdant2o.")')
        call stoponerror
      end if
      do n=0,ncalcs-1
        k2=0
        do i1=1,nm
          do j1=1,ncb(i1)
c for testing
c        do i1=1,1
c          do j1=1,1
            k2=k2+1
            iso=nm+k2 ! nm is the offset to the isotope cross sections
c offset between isotopes
            noff3=neg*(k2-1)
            noff4=neg*neg*nmom*(k2-1)
c the ng loop is placed here for the consideration of writing partisn
c inputs. Note that each energy group corresponds to a parameter, thus
c a forward and an adjoint partisn input are needed.   
            do ng=1,neg
c arrays need to be initialized inside the ng loop
              neusf(1:3,1:neg,0:nmom-1,1:it,1:jt)=0.d0
              neusa(1:3,1:neg,0:nmom-1,1:it,1:jt)=0.d0
              do ir=1,nr
                do jr=1,nz
                  if(mat(ir,jr).ne.i1)cycle
c the OMP PARALLEL is removed in consideration of fewer array operations
c and smaller arrays. not sure the performance comparing with the parallel.
c please consider to add them back if necessary.      
                  do i=iindex(ir-1)+1,iindex(ir)
                    do j=jindex(jr-1)+1,jindex(jr)
                      ll=-1
                      do l=0,isct
                        do k=0,l
                          ll=ll+1
c source for sigt of each parameter
                          k3=noff3+ng
c   Eq. (164) in Ref. 1
                          neusf(1,ng,ll,i,j)=-fmom(ng,ll,i,j)*blk(3,k2)
     1                     *rlip(l)
                          if(jf(k3).eq.0)then
                            if(neusf(1,ng,ll,i,j).ne.0.d0)jf(k3)=1
                          end if
c   Eq. (165) in Ref. 1
                          neusa(1,ng,ll,i,j)=-amom(ng,ll,i,j)*blk(3,k2)
     1                     *rlip(l)     
                          if(ja(k3).eq.0)then
                            if(neusa(1,ng,ll,i,j).ne.0.d0)ja(k3)=1
                          end if
c source for nu of each parameter
                          k3=noff2+noff3+ng
c   Eq. (183) in Ref. 1 for nu
                          if(ll.eq.0)then
                            do lg=1,neg
                              st=blk(3,k2)*sigf(ng,iso)
                              neusf(2,lg,ll,i,j)=
     1                         +fmom(ng,ll,i,j)*chivec(lg,i1)*st
                            end do ! lg
                          end if
                          if(jf(k3).eq.0)then
                            if(neusf(2,ng,ll,i,j).ne.0.d0)jf(k3)=1
                          end if
c   Eq. (184) in Ref. 1 for nu
                          if(ll.eq.0)then
                            st=blk(3,k2)*sigf(ng,iso)
                            do lg=1,neg
                             neusa(2,ng,ll,i,j)=neusa(2,ng,ll,i,j)
     1                        +amom(lg,ll,i,j)*chivec(lg,i1)*st
                            end do ! lg
                          end if
                          if(ja(k3).eq.0)then
                            if(neusa(2,ng,ll,i,j).ne.0.d0)ja(k3)=1
                          end if
c source for sigf of each parameter
                          k3=2*noff2+noff3+ng
c   Eq. (183) in Ref. 1 for sigf
                          if(ll.eq.0)then
                            do lg=1,neg
                              if(sigf(ng,iso).ne.0.d0)then 
                               st=blk(3,k2)*nusigf(ng,iso)/sigf(ng,iso)
                               neusf(3,lg,ll,i,j)=
     1                          +fmom(ng,ll,i,j)*chivec(lg,i1)*st
                              end if
                            end do ! lg
                          end if
                          if(jf(k3).eq.0)then
                            if(neusf(3,ng,ll,i,j).ne.0.d0)jf(k3)=1
                          end if
c   Eq. (184)  in Ref. 1 for sigf
                          if(ll.eq.0)then
                            if(sigf(ng,iso).ne.0.d0)then
                              st=blk(3,k2)*nusigf(ng,iso)/sigf(ng,iso)
                              do lg=1,neg
                                neusa(3,ng,ll,i,j)=neusa(3,ng,ll,i,j)
     1                           +amom(lg,ll,i,j)*chivec(lg,i1)*st
                              end do ! lg
                            end if
                          end if
                          if(ja(k3).eq.0)then
                            if(neusa(3,ng,ll,i,j).ne.0.d0)ja(k3)=1
                          end if
c for one-d, exit k loop (still in l loop).
                          if(calc2d.eq.-1.or.calc2d.eq.0)exit
                        end do ! k
                      end do ! l
                    end do ! j
                  end do ! i
                end do ! jr
              end do ! ir
c write partisn forward/adjoint inputs for each parameter of sigt, nu, sigf. 
c has to be here because the parameters are defined for each energy group.  
              k3=noff3+ng
              if(jf(k3).eq.1)then
                call wrdant2o(0,k3,niso0,
     1           neusf(1,1:neg,0:nmom-1,1:it,1:jt))
              end if
              if(ja(k3).eq.1)then
                call wrdant2o(1,k3,niso0,
     1           neusa(1,1:neg,0:nmom-1,1:it,1:jt))
              end if
              k3=noff2+noff3+ng
              if(jf(k3).eq.1)then
                call wrdant2o(0,k3,niso0,
     1           neusf(2,1:neg,0:nmom-1,1:it,1:jt))
              end if
              if(ja(k3).eq.1)then
                call wrdant2o(1,k3,niso0,
     1           neusa(2,1:neg,0:nmom-1,1:it,1:jt))
              end if
              k3=2*noff2+noff3+ng
              if(jf(k3).eq.1)then
                call wrdant2o(0,k3,niso0,
     1           neusf(3,1:neg,0:nmom-1,1:it,1:jt))
              end if
              if(ja(k3).eq.1)then
                call wrdant2o(1,k3,niso0,
     1           neusa(3,1:neg,0:nmom-1,1:it,1:jt))
              end if
            end do ! ng
c
c compute sources for sigs
c
c the l loop is placed outside other loops (i.e., ng, lg, ir, jr, i, j),
c it is in the consideration of writing the partisn inputs. for the same
c reason, the ng loop and lg loop are placed outside the ir, jr, i, and j
c loops.   
            do l=0,isct
c for testing
c            do l=3,3
              noff5=neg*neg*l
              do ng=1,neg
                do lg=1,neg
                  do ir=1,nr
                    do jr=1,nz
                      if(mat(ir,jr).ne.i1)cycle
                        do i=iindex(ir-1)+1,iindex(ir)
                          do j=jindex(jr-1)+1,jindex(jr)
c
c source for sigs(ng,lg,l,i)
c 
c k3 is used to index each parameters of sigs(ng,lg,l,i).
c it is first offset by "3*noff2" so it points to the first  
c parameter of sigs; then offest by "noff4" to point to the correct 
c isotope; then offset by "noff5" to point to the correct l-th order 
c of sigs; then offset by "neg*(ng-1)" to point to the correct energy
c group inside the l-th order sub-block; finally, "lg" accounts for 
c the incremental in each ng loop. Thus, in explicit expresion,  
c k3=3*neg*nel+neg*neg*nmom*(k2-1)+neg*neg*l+neg*(ng-1)+lg 
c   
                            k3=3*noff2+noff4+noff5+neg*(ng-1)+lg
c Eq. (173) in Ref. 1 for sigs
                            if(sigs(ng,lg,l,iso).ne.0.d0)then
                              neusf_scat(ng,l,i,j)=
     1                         fmom(lg,l,i,j)*blk(3,k2)*rlip(l)
c partisn requires neusf_scat(ng,0,i,j) must have a non-zero value  
                              if(l.gt.0)neusf_scat(ng,0,i,j)=
     1                           neusf_scat(ng,l,i,j)*1.d-24
                              jf(k3)=1
                            end if
c Eq. (174) in Ref. 1 for sigs
                            if(sigs(ng,lg,l,iso).ne.0.d0)then
                              neusa_scat(lg,l,i,j)=
     1                         amom(ng,l,i,j)*blk(3,k2)*rlip(l)     
c partisn requires neusa_scat(ng,0,i,j) must have a non-zero value
                              if(l.gt.0)neusa_scat(lg,0,i,j)=
     1                           neusa_scat(lg,l,i,j)*1.d-24
                              ja(k3)=1
                            end if
c for 2-D cylinders, those sources need a two-dimensional treatment.
                          end do ! j
                        end do ! i
                      end do ! jr
                    end do ! ir
c write partisn forward/adjoint inputs for each parameter of sigs. 
c has to be here because the parameters are defined at the lg energy-group.
c 
c if only want to write partisn files for self-scatter cross sections,use:
c                  if(k3.eq.(3*noff2+noff4+noff5+neg*(ng-1)+ng))then    
c write partisn input files for all scttering cross section parameters 
                    if(jf(k3).eq.1)then
                      call wrdant2o(0,k3,niso0,
     1                 neusf_scat(1:neg,0:nmom-1,1:it,1:jt))
                    end if
                    if(ja(k3).eq.1)then
                      call wrdant2o(1,k3,niso0,
     1                 neusa_scat(1:neg,0:nmom-1,1:it,1:jt))
                    end if
c                  end if ! for self-scatter cross sections only
c reset only if they are not zero  
                  if(jf(k3).eq.1)
     1              neusf_scat(1:neg,0:nmom-1,1:it,1:jt)=0.d0
                  if(ja(k3).eq.1)
     1              neusa_scat(1:neg,0:nmom-1,1:it,1:jt)=0.d0
                end do ! lg
              end do ! ng
            end do ! l
c
          end do ! j1
        end do ! i1
      end do ! n
c
c write partisn input
      write(*,'(/,"written partisn inputs for second order. ",
     1 "there were ",i7," forward files and ",i7," adjoint files.",
     2 " nel=",i5)')sum(jf(1:noff1)),sum(ja(1:noff1)),nel
      write(iuo,'(/,"written partisn inputs for second order. ",
     1 "there were ",i7," forward files and ",i7," adjoint files.",
     2 " nel=",i3)')sum(jf(1:noff1)),sum(ja(1:noff1)),nel
c
      return
      end subroutine
