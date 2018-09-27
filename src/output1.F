      subroutine output1(sden,stot,sabs,scap,schi,snu,sfis,
     1 ssct,sscti,sscts,sscto,sins,ssctt)
c output indirect sensitivities and totals.
      use COMS
      implicit none
      real*8 sden(neg),stot(neg),sabs(neg),scap(neg),
     1 sfis(neg),snu(neg,neg),schi(neg,neg),
     2 ssct(neg,neg,0:isct),sscti(neg,0:isct),
     3 sscto(neg,0:isct),sscts(neg,0:isct),sins(neg),
     4 ssctt(neg)
      integer l,lg,ng,n1,n2
c head is the heading " group  density...." length must be greater than
c 6+14*9+14*3*isct = 132+42*isct. 510 is ok for isct.le.9.
      character head*510,cchinorm(0:2)*21,ccna(0:2)*5
c there are several instances of 9999e14.6 and 9999i14. the number
c of energy groups neg must be .le. 9999 and the scattering order
c isct must satisfy 13+3*isct.le.9999 or isct.lt.3300.
      data cchinorm/"no normalization     ","full normalization   ",
     1 "partial normalization"/
      data ccna/"_(nn)","_(fn)","_(pn)"/
c
c summary table first, just density and reaction totals
      write(ius,'(" summary (sums over groups)")')
      write(ius,'(" density",4x,1pe16.6,/," total",6x,e16.6,/,
     1 " abs",8x,e16.6,/," (n,g)",6x,e16.6,/," chi",a5,3x,e16.6,/,
     2 " nu",9x,e16.6,/," fiss",7x,e16.6)')sum(sden(1:neg)),
     3 sum(stot(1:neg)),sum(sabs(1:neg)),
     4 sum(scap(1:neg)),ccna(ichinorm),
     5 sum(schi(1:neg,1:neg)),sum(snu(1:neg,1:neg)),
     6 sum(sfis(1:neg))
      do l=0,isct
        write(ius,'(" in-scat-",i1,2x,1pe16.6,/," self-scat-",i1,
     1   e16.6,/," out-scat-",i1,1x,e16.6)')
     2   l,sum(sscti(1:neg,l)),l,sum(sscts(1:neg,l)),
     3   l,sum(sscto(1:neg,l))
      end do ! l
      if(icalc.eq.0.or.icalc.eq.3.or.icalc.eq.4)then
        write(ius,'(" source",5x,1pe16.6)')sum(sins(1:neg))
      end if 
      write(ius,'(" ssctt",6x,1pe16.6)')sum(ssctt(1:neg))
c now the details
      head=" "
      write(head,'(1x,"group",2x,"e_lower",7x,"e_upper",7x,
     1 "e_average",5x,"density",7x,"total",9x,"abs",11x,"(n,g)",
     1 9x,"chi",a5,6x,"nu",12x,"fiss")')ccna(ichinorm)
      n1=6+10*14+3 ! i6,10e14.6, 3 spaces
      do l=0,isct
        n2=n1+14*3
        write(head(n1:n2),'("in-scat-",i1,5x,"self-scat-",i1,
     1   3x,"out-scat-",i1,2x)')l,l,l
        n1=n2
      end do ! l
      if(icalc.eq.0.or.icalc.eq.3.or.icalc.eq.4)then
        n2=n1+14
        write(head(n1:n2),'("source",8x)')
        n1=n2
      end if 
      n2=n1+14
      write(head(n1:n2),'("ssctt ",8x)')
      n1=n2
      write(ius,'(a)')trim(head)
      do ng=1,neg
        if(icalc.eq.0.or.icalc.eq.3.or.icalc.eq.4)then
          write(ius,'(i6,1p9999e14.6)')ng,ebins(1:3,ng),
     1     sden(ng),stot(ng),sabs(ng),scap(ng),
     2     sum(schi(ng,1:neg)),
     3     sum(snu(1:neg,ng)),sfis(ng),
     4     (sscti(ng,l),sscts(ng,l),sscto(ng,l),
     5     l=0,isct),sins(ng),ssctt(ng)
        else
          write(ius,'(i6,1p9999e14.6)')ng,ebins(1:3,ng),
     1     sden(ng),stot(ng),sabs(ng),scap(ng),
     2     sum(schi(ng,1:neg)),
     3     sum(snu(1:neg,ng)),sfis(ng),
     4     (sscti(ng,l),sscts(ng,l),sscto(ng,l),
     5     l=0,isct),ssctt(ng)
        end if 
      end do ! ng
c 42x is 14*3 for the energy groups.
      if(icalc.eq.0.or.icalc.eq.3.or.icalc.eq.4)then
        write(ius,'(" total",42x,1p9999e14.6)')sum(sden(1:neg)),
     1   sum(stot(1:neg)),sum(sabs(1:neg)),
     2   sum(scap(1:neg)),
     3   sum(schi(1:neg,1:neg)),sum(snu(1:neg,1:neg)),
     4   sum(sfis(1:neg)),
     5   (sum(sscti(1:neg,l)),sum(sscts(1:neg,l)),
     6   sum(sscto(1:neg,l)),l=0,isct),sum(sins(1:neg)),
     7   sum(ssctt(1:neg))
      else
        write(ius,'(" total",42x,1p9999e14.6)')sum(sden(1:neg)),
     1   sum(stot(1:neg)),sum(sabs(1:neg)),
     2   sum(scap(1:neg)),
     3   sum(schi(1:neg,1:neg)),sum(snu(1:neg,1:neg)),
     4   sum(sfis(1:neg)),
     5   (sum(sscti(1:neg,l)),sum(sscts(1:neg,l)),
     6   sum(sscto(1:neg,l)),l=0,isct),sum(ssctt(1:neg))
      end if
      if(sum(abs(schi(1:neg,1:neg))).ne.0.d0)then
        write(ius,'("fission chi sensitivity matrix",
     1   " (columns are from, rows are to), ",a)')
     2   trim(cchinorm(ichinorm))
        write(ius,'(9999i14)')(ng,ng=1,neg)
        do ng=1,neg
          write(ius,'(i6,1p9999e14.6)')ng,(schi(ng,lg),lg=1,neg)
        end do ! ng
      end if
      do l=0,isct
        write(ius,'("scattering sensitivity matrix, order=",i2,
     1   " (columns are outscattering, rows are inscattering)")')l
        write(ius,'(9999i14)')(ng,ng=1,neg)
        do ng=1,neg
          write(ius,'(i6,1p9999e14.6)')ng,(ssct(ng,lg,l),lg=1,neg)
        end do ! ng
      end do ! l
c
      return
      end subroutine
c
      subroutine outputh(iu,n,irrx1,irrx2,icalc,rxnratio)
      implicit none
      real*8 rxnratio
      integer iu,n,irrx1(2),irrx2(2),icalc
c
      if(irrx2(1).gt.0)then
        write(iu,'(/," ratio",i3,": numer zaid=",i6," rxn=",i2,
     1   "; denom zaid=",i6," rxn=",i2,10x,"ratio=",1pe14.6)')n,
     2   irrx1(1:2),irrx2(1:2),rxnratio
      else if(icalc.eq.1.or.icalc.eq.2)then
        write(iu,'(/," ratio",i3,": numer zaid=",i6," rxn=",i2,
     1   "; denom = flux",22x,"ratio=",1pe14.6)')n,
     1   irrx1(1:2),rxnratio
      else if(icalc.eq.0.or.icalc.eq.3.or.icalc.eq.4)then
        write(iu,'(/," ratio",i3,": numer zaid=",i6," rxn=",i2,
     1   "; denom = none",22x,"ratio=",1pe14.6)')n,
     1   irrx1(1:2),rxnratio
      end if
      return
      end subroutine
