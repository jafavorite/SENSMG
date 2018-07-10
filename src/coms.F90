module F77KINDS
  !
  ! Define various compiler-dependent F77 variable types
  !
  implicit none

  COMPLEX*16 F77KINDS_c16
  COMPLEX*8 F77KINDS_c8

  REAL*8 F77KINDS_r8

  REAL*4 F77KINDS_r4
  INTEGER*4 F77KINDS_i4
  LOGICAL*4 F77KINDS_l4

  INTEGER*2 F77KINDS_i2
  LOGICAL*2 F77KINDS_l2

  INTEGER*1 F77KINDS_i1
  LOGICAL*1 F77KINDS_l1

end module F77KINDS
!--------------------------------------------------------------------------------
module F90KINDS
  !
  ! Dynamically define the compiler-dependent F90 kind parameters
  !
  use  F77KINDS

  implicit none

  INTEGER, PARAMETER :: &
       I1KIND  = KIND(F77KINDS_i1), & ! INTEGER*1
       I2KIND  = KIND(F77KINDS_i2), & ! INTEGER*2
       I4KIND  = KIND(F77KINDS_i4), & ! INTEGER*4
       L1KIND  = KIND(F77KINDS_l1), & ! LOGICAL*1
       L2KIND  = KIND(F77KINDS_l2), & ! LOGICAL*2
       L4KIND  = KIND(F77KINDS_l4), & ! LOGICAL*4
       R4KIND  = KIND(F77KINDS_r4), & ! REAL*4
       R8KIND  = KIND(F77KINDS_r8), & ! REAL*8
       CX8KIND = KIND(F77KINDS_c8), & ! COMPLEX*8
       CX16KIND= KIND(F77KINDS_c16)   ! COMPLEX*16

end module F90KINDS
!--------------------------------------------------------------------------------
module MISC_CONSTANTS
  !
  ! Define miscellaneous parameters
  !
  use F90KINDS

  real(R8KIND) :: avg,pi,pi43
! MCNP value for Avogadro's number
! parameter (avg=0.60220434469282d0)
! PARTISN value for Avogadro's number
  parameter (avg=0.602214129d0)
  parameter (pi=3.1415926535898d+0,pi43=4.d0*pi/3.d0)
  character(27), parameter :: abc = "abcdefghijklmnopqrstuvwxyz "
  character(62), parameter :: abc123 = &
     "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"
  integer, parameter :: iui =  2,  &  ! input file, rdmdl
                        iuc =  7,  &  ! control file, rdctrl
                        iuk =  8,  &  ! forward output for keff, rddantk
                        iug =  9,  &  ! gendir, rdgendir
                        iug2= 11,  &  ! energy group data file or detector efficiency data file
                        iun = 12,  &  ! sources or misc output files
                        iuo = 56,  &  ! sensmg.log
                        ius = 57,  &  ! output sensitivities, calcsens
                        iur = 58,  &  ! output sensitivities, calcsens_r
                        iup = 59,  &  ! partisn input file, wrdantga,
                                      ! wrdantnm, wrdantxs, wrodninp.
                                      ! different file is open each time
                        iue = 60,  &  ! snxedt file, rdsnxedth, rdsnxedt
                        ium = 61,  &  ! rmflux or amflux file, rddantm
                        iux = 62,  &  ! macrsx file, rdmacrxs, wrmacrxs
                        iux2= 63,  &  ! macrsx file to write, wrmacrxs
                        iuga= 64,  &  ! gen. adjoint moments, rdgmom, wrgmom
                                      ! gen. adjoint ang. fluxes, rdgang, wrgang
                        iua = 65,  &  ! raflxm or aaflxm file, rddanta
                        iuf = 66,  &  ! fixsrc, wrfixsrc
                        iub = 67,  &  ! senslx, binary version of sens_l_x; senssm
                        iut0= 89,  &  ! input to makemg
                        iut = 90      ! messages to control script, open and close

end module MISC_CONSTANTS
!--------------------------------------------------------------------------------
module ARRAY_DIM_CONSTANTS
  !
  ! Define hard-wired array dimension parameters
  !
  ! maxedits is the number of edits available from partisn, =20 in partisn 8
  !
  use F90KINDS

  integer, parameter :: &
       maxedits=20

end module ARRAY_DIM_CONSTANTS
!--------------------------------------------------------------------------------
module GEOM
  !
  ! geometry and materials
  !
  use F90KINDS

  real(R8KIND) :: suvo
  integer isn,it,jt,nel,nxs,niso,nr,nz,nm,lastr,ilkg
  integer isct,ndir,nmom
  integer calc2d ! -1/0/1 slab/sphere/cylinder
  integer imisc ! 0/1 sources/misc used to compute neutron source
  integer ialphan ! 0/1 0/1 do not/do include (alpha,n) sources
  integer irrcmr,irrcmz ! indices of coarse mesh to use for reaction rate calcs
  integer nrrr ! number of reaction rate ratios to compute
  integer nrrx ! number of unique reaction rates
  integer nflux ! 0/1 flux is not/is one of the calculations
  integer nag ! number of alpha-particle groups in sources4c calculation
  integer iplane,jplane ! r and z surfaces to write bs files (temporary)
  integer cellsol ! 0/1 partisn's cellsol

end module GEOM
!--------------------------------------------------------------------------------
module ENERG
  !
  use F90KINDS

  ! neg is the number of groups in a neutron problem.
  integer neg,lng

end module ENERG
!--------------------------------------------------------------------------------
module XSECS
  !
  ! sources and cross section files
  !
  use F90KINDS
  use ARRAY_DIM_CONSTANTS

  ! fissdata = 0/1/2 fission transfer matrix/chi matrix/chi vector
  integer maxup,fissdata
  character(220) :: ndi,datapath,sens_data
  character(8)   :: libname
  character(8)   :: detname
  character(8)   :: cedits(0:maxedits)

end module XSECS
!--------------------------------------------------------------------------------
module FLUXES
  !
  ! fluxes and detected quantities
  !
  use F90KINDS

  real(R8KIND) :: keff,alpha,lkg,yasym,sm2

end module FLUXES
!--------------------------------------------------------------------------------
module VAR
  !
  ! ictrl is the control character read from the control file
  ! iter is the number of iterations of the generalized adjoint solver
  ! icalc = 0/1/2 fixed-source/keff/alpha
  ! nofxup is the flag to tell whether to write "nofxup=1"
  ! aflxfrm = 0/1 no/yes Use the angular flux formulation (ievt=2 & dsasrch=2 only).
  ! idbgw is the flag to tell whether to write the debugging files
  ! epsi is the convergence criteria epsi to use in all dant calculations
  ! epsig is the convergence criteria for the successive approximations of the
  !  generalized adjoint functions.
  ! idbgw is the flag to tell whether to write verification files.
  ! iangflux = 0/1 use only moments/use angular fluxes for sigt term
  ! ichinorm = fission spectrum sensitivity normalization;
  !  0/1/2 unnormalized/full norm/partial norm
  ! isrcacc_no = 0/1/2/3 turns off source acceleration for
  !  neither/forward/adjoint/both.
  ! iaflux = 0/1 don't/do use the afluxx/afluxy partisn feature (8_27 and later).
  ! num_threads is the number of threads to use for inner products
  use F90KINDS

  integer ictrl,iter,icalc,nofxup,aflxfrm,idbgw,iver,iangflux,ichinorm,isrcacc_no, &
   iaflux,num_threads
  real(R8KIND) :: epsi,epsig

end module VAR
!--------------------------------------------------------------------------------
module CHARAC
  !
  ! characters
  !
  use F90KINDS
  character(120) :: id
  character(120) :: ifile
  character(120) :: partisn

end module CHARAC
!--------------------------------------------------------------------------------
module COMS
  !
  use F90KINDS
  use MISC_CONSTANTS
  use ARRAY_DIM_CONSTANTS
  use GEOM
  use ENERG
  use XSECS
  use FLUXES
  use CHARAC
  use VAR

  implicit none

! arrays allocated in allocate_arrays_1
  real(R8KIND), allocatable, dimension(:)       :: r, z, dr, dz, rho, rxnratio
  real(R8KIND), allocatable, dimension(:,:)     :: mass, vol, blk
  integer, allocatable, dimension(:)            :: iints, jints, ncb, ismat, isan
  integer, allocatable, dimension(:,:)          :: mat, ifcel, irrr, irri, irrx

! arrays allocated in allocate_arrays
  real(R8KIND), allocatable, dimension(:)       :: rfm, zfm, vel, dir, wgt, &
    eta, atwt, nsrc, deteff
  real(R8KIND), allocatable, dimension(:,:)     :: dv, sar, sigt, nusigf, siga, sigf, sigc, &
    rrxs, ebins, nsrcf, sfiso, saniso, chisrc
  real(R8KIND), allocatable, dimension(:,:,:)   :: scalar, chi
  real(R8KIND), allocatable, dimension(:,:,:,:) :: fmom, amom, gmom, afreg, afadj, &
    afgad, sigs
  integer, allocatable, dimension(:)            :: iindex, jindex
  integer, allocatable, dimension(:,:,:,:)      :: scgr
  character, allocatable, dimension(:)          :: zaidfull*24

! arrays allocated in allocate_arrays3
  real(R8KIND), allocatable, dimension(:,:,:) :: omia
  real(R8KIND), allocatable, dimension(:,:,:,:) :: afregas, afadjas, forsa, adjsa

  logical :: force_alloc = .FALSE.
  integer :: nitm
  integer :: njtm

  CONTAINS

  subroutine allocate_arrays_1()
    !
    ! This subroutine is used to allocate dynamic arrays using nr, nz,
    ! nm, nel, nrrr
    !
    integer :: ierr

    ierr = 0

    if (force_alloc) then
       deallocate(r, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("ERROR.  cannot deallocate array: ",a,".")') &
               'r(0:nr)'
          call stoponerror
       end if
    end if
    if (.NOT. allocated(r)) &
         allocate(r(0:nr), STAT=ierr)
    if(ierr /= 0)then
       write(*,'("ERROR.  cannot allocate array: ",a,".")') &
            'r(0:nr)'
       call stoponerror
    end if

    if (force_alloc) then
       deallocate(z, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("ERROR.  cannot deallocate array: ",a,".")') &
               'z(0:nz)'
          call stoponerror
       end if
    end if
    if (.NOT. allocated(z)) &
         allocate(z(0:nz), STAT=ierr)
    if(ierr /= 0)then
       write(*,'("ERROR.  cannot allocate array: ",a,".")') &
            'z(0:nz)'
       call stoponerror
    end if

    if (force_alloc) then
       deallocate(dr, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("ERROR.  cannot deallocate array: ",a,".")') &
               'dr(nr)'
          call stoponerror
       end if
    end if
    if (.NOT. allocated(dr)) &
         allocate(dr(nr), STAT=ierr)
    if(ierr /= 0)then
       write(*,'("ERROR.  cannot allocate array: ",a,".")') &
            'dr(nr)'
       call stoponerror
    end if

    if (force_alloc) then
       deallocate(dz, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("ERROR.  cannot deallocate array: ",a,".")') &
               'dz(nz)'
          call stoponerror
       end if
    end if
    if (.NOT. allocated(dz)) &
         allocate(dz(nz), STAT=ierr)
    if(ierr /= 0)then
       write(*,'("ERROR.  cannot allocate array: ",a,".")') &
            'dz(nz)'
       call stoponerror
    end if

    if (force_alloc) then
       deallocate(rho, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("ERROR.  cannot deallocate array: ",a,".")') &
               'rho(0:nm)'
          call stoponerror
       end if
    end if
    if (.NOT. allocated(rho)) &
         allocate(rho(0:nm), STAT=ierr)
    if(ierr /= 0)then
       write(*,'("ERROR.  cannot allocate array: ",a,".")') &
            'rho(0:nm)'
       call stoponerror
    end if

    if (force_alloc) then
       deallocate(rxnratio, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("ERROR.  cannot deallocate array: ",a,".")') &
               'rxnratio(nrrr)'
          call stoponerror
       end if
    end if
    if (.NOT. allocated(rxnratio)) &
         allocate(rxnratio(nrrr), STAT=ierr)
    if(ierr /= 0)then
       write(*,'("ERROR.  cannot allocate array: ",a,".")') &
            'rxnratio(nrrr)'
       call stoponerror
    end if

    if (force_alloc) then
       deallocate(mass, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("ERROR.  cannot deallocate array: ",a,".")') &
               'mass(nr,nz)'
          call stoponerror
       end if
    end if
    if (.NOT. allocated(mass)) &
         allocate(mass(nr,nz), STAT=ierr)
    if(ierr /= 0)then
       write(*,'("ERROR.  cannot allocate array: ",a,".")') &
            'mass(nr,nz)'
       call stoponerror
    end if

    if (force_alloc) then
       deallocate(vol, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("ERROR.  cannot deallocate array: ",a,".")') &
               'vol(nr,nz)'
          call stoponerror
       end if
    end if
    if (.NOT. allocated(vol)) &
         allocate(vol(nr,nz), STAT=ierr)
    if(ierr /= 0)then
       write(*,'("ERROR.  cannot allocate array: ",a,".")') &
            'vol(nr,nz)'
       call stoponerror
    end if

    if (force_alloc) then
       deallocate(blk, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("ERROR.  cannot deallocate array: ",a,".")') &
               'blk(3,nel)'
          call stoponerror
       end if
    end if
    if (.NOT. allocated(blk)) &
         allocate(blk(3,nel), STAT=ierr)
    if(ierr /= 0)then
       write(*,'("ERROR.  cannot allocate array: ",a,".")') &
            'blk(3,nel)'
       call stoponerror
    end if

    if (force_alloc) then
       deallocate(iints, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("ERROR.  cannot deallocate array: ",a,".")') &
               'iints(nr)'
          call stoponerror
       end if
    end if
    if (.NOT. allocated(iints)) &
         allocate(iints(nr), STAT=ierr)
    if(ierr /= 0)then
       write(*,'("ERROR.  cannot allocate array: ",a,".")') &
            'iints(nr)'
       call stoponerror
    end if

    if (force_alloc) then
       deallocate(jints, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("ERROR.  cannot deallocate array: ",a,".")') &
               'jints(nz)'
          call stoponerror
       end if
    end if
    if (.NOT. allocated(jints)) &
         allocate(jints(nz), STAT=ierr)
    if(ierr /= 0)then
       write(*,'("ERROR.  cannot allocate array: ",a,".")') &
            'jints(nz)'
       call stoponerror
    end if

    if (force_alloc) then
       deallocate(ncb, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("ERROR.  cannot deallocate array: ",a,".")') &
               'ncb(nm)'
          call stoponerror
       end if
    end if
    if (.NOT. allocated(ncb)) &
         allocate(ncb(nm), STAT=ierr)
    if(ierr /= 0)then
       write(*,'("ERROR.  cannot allocate array: ",a,".")') &
            'ncb(nm)'
       call stoponerror
    end if

    if (force_alloc) then
       deallocate(ismat, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("ERROR.  cannot deallocate array: ",a,".")') &
               'ismat(0:nm)'
          call stoponerror
       end if
    end if
    if (.NOT. allocated(ismat)) &
         allocate(ismat(0:nm), STAT=ierr)
    if(ierr /= 0)then
       write(*,'("ERROR.  cannot allocate array: ",a,".")') &
            'ismat(0:nm)'
       call stoponerror
    end if

    if (force_alloc) then
       deallocate(isan, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("ERROR.  cannot deallocate array: ",a,".")') &
               'isan(0:nm)'
          call stoponerror
       end if
    end if
    if (.NOT. allocated(isan)) &
         allocate(isan(0:nm), STAT=ierr)
    if(ierr /= 0)then
       write(*,'("ERROR.  cannot allocate array: ",a,".")') &
            'isan(0:nm)'
       call stoponerror
    end if

    if (force_alloc) then
       deallocate(mat, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("ERROR.  cannot deallocate array: ",a,".")') &
               'mat(0:nr+1,0:nz+1)'
          call stoponerror
       end if
    end if
    if (.NOT. allocated(mat)) &
         allocate(mat(0:nr+1,0:nz+1), STAT=ierr)
    if(ierr /= 0)then
       write(*,'("ERROR.  cannot allocate array: ",a,".")') &
            'mat(0:nr+1,0:nz+1)'
       call stoponerror
    end if

    if (force_alloc) then
       deallocate(ifcel, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("ERROR.  cannot deallocate array: ",a,".")') &
               'ifcel(0:nr+1,0:nz+1)'
          call stoponerror
       end if
    end if
    if (.NOT. allocated(ifcel)) &
         allocate(ifcel(0:nr+1,0:nz+1), STAT=ierr)
    if(ierr /= 0)then
       write(*,'("ERROR.  cannot allocate array: ",a,".")') &
            'ifcel(0:nr+1,0:nz+1)'
       call stoponerror
    end if

    if (force_alloc) then
       deallocate(irrr, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("ERROR.  cannot deallocate array: ",a,".")') &
               'irrr(4,nrrr)'
          call stoponerror
       end if
    end if
    if (.NOT. allocated(irrr)) &
         allocate(irrr(4,nrrr), STAT=ierr)
    if(ierr /= 0)then
       write(*,'("ERROR.  cannot allocate array: ",a,".")') &
            'irrr(4,nrrr)'
       call stoponerror
    end if

    if (force_alloc) then
       deallocate(irri, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("ERROR.  cannot deallocate array: ",a,".")') &
               'irri(2,nrrr)'
          call stoponerror
       end if
    end if
    if (.NOT. allocated(irri)) &
         allocate(irri(2,nrrr), STAT=ierr)
    if(ierr /= 0)then
       write(*,'("ERROR.  cannot allocate array: ",a,".")') &
            'irri(2,nrrr)'
       call stoponerror
    end if

    if (force_alloc) then
       deallocate(irrx, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("ERROR.  cannot deallocate array: ",a,".")') &
               'irrx(2,2*nrrr)'
          call stoponerror
       end if
    end if
    if (.NOT. allocated(irrx)) &
         allocate(irrx(2,2*nrrr), STAT=ierr)
    if(ierr /= 0)then
       write(*,'("ERROR.  cannot allocate array: ",a,".")') &
            'irrx(2,2*nrrr)'
       call stoponerror
    end if

    return
  end subroutine allocate_arrays_1


  subroutine allocate_arrays()
    !
    ! This subroutine is used to allocate several dynamic arrays
    !
    integer :: ierr

    ierr = 0

    if (force_alloc) then
       deallocate(ebins, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("ERROR.  cannot deallocate array: ",a,".")') &
               'ebins(3,neg)'
          call stoponerror
       end if
    end if
    if (.NOT. allocated(ebins)) &
         allocate(ebins(3,neg), STAT=ierr)
    if(ierr /= 0)then
       write(*,'("ERROR.  cannot allocate array: ",a,".")') &
            'ebins(3,neg)'
       call stoponerror
    end if

    if (force_alloc) then
       deallocate(rrxs, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("ERROR.  cannot deallocate array: ",a,".")') &
               'rrxs(neg,0:nrrx+nflux)'
          call stoponerror
       end if
    end if
    if (.NOT. allocated(rrxs)) &
         allocate(rrxs(neg,0:nrrx+nflux), STAT=ierr)
    if(ierr /= 0)then
       write(*,'("ERROR.  cannot allocate array: ",a,".")') &
            'rrxs(neg,0:nrrx+nflux)'
       call stoponerror
    end if

    if (force_alloc) then
       deallocate(atwt, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("ERROR.  cannot deallocate array: ",a,".")') &
               'atwt(niso)'
          call stoponerror
       end if
    end if
    if (.NOT. allocated(atwt)) &
         allocate(atwt(niso), STAT=ierr)
    if(ierr /= 0)then
       write(*,'("ERROR.  cannot allocate array: ",a,".")') &
            'atwt(niso)'
       call stoponerror
    end if

    if (force_alloc) then
       deallocate(sigt, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("ERROR.  cannot deallocate array: ",a,".")') &
               'sigt(neg,0:nxs)'
          call stoponerror
       end if
    end if
    if (.NOT. allocated(sigt)) &
         allocate(sigt(neg,0:nxs), STAT=ierr)
    if(ierr /= 0)then
       write(*,'("ERROR.  cannot allocate array: ",a,".")') &
            'sigt(neg,0:nxs)'
       call stoponerror
    end if

    if (force_alloc) then
       deallocate(sigs, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("ERROR.  cannot deallocate array: ",a,".")') &
               'sigs(neg,neg,0:isct,0:nxs)'
          call stoponerror
       end if
    end if
    if (.NOT. allocated(sigs)) &
         allocate(sigs(neg,neg,0:isct,0:nxs), STAT=ierr)
    if(ierr /= 0)then
       write(*,'("ERROR.  cannot allocate array: ",a,".")') &
            'sigs(neg,neg,0:isct,0:nxs)'
       call stoponerror
    end if

    if (force_alloc) then
       deallocate(nusigf, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("ERROR.  cannot deallocate array: ",a,".")') &
               'nusigf(neg,0:nxs)'
          call stoponerror
       end if
    end if
    if (.NOT. allocated(nusigf)) &
         allocate(nusigf(neg,0:nxs), STAT=ierr)
    if(ierr /= 0)then
       write(*,'("ERROR.  cannot allocate array: ",a,".")') &
            'nusigf(neg,0:nxs)'
       call stoponerror
    end if

    if (force_alloc) then
       deallocate(chi, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("ERROR.  cannot deallocate array: ",a,".")') &
               'chi(neg,neg,0:nxs)'
          call stoponerror
       end if
    end if
    if (.NOT. allocated(chi)) &
         allocate(chi(neg,neg,0:nxs), STAT=ierr)
    if(ierr /= 0)then
       write(*,'("ERROR.  cannot allocate array: ",a,".")') &
            'chi(neg,neg,0:nxs)'
       call stoponerror
    end if

    if (force_alloc) then
       deallocate(vel, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("ERROR.  cannot deallocate array: ",a,".")') &
               'vel(neg)'
          call stoponerror
       end if
    end if
    if (.NOT. allocated(vel)) &
         allocate(vel(neg), STAT=ierr)
    if(ierr /= 0)then
       write(*,'("ERROR.  cannot allocate array: ",a,".")') &
            'vel(neg)'
       call stoponerror
    end if

    if (force_alloc) then
       deallocate(siga, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("ERROR.  cannot deallocate array: ",a,".")') &
               'siga(neg,0:nxs)'
          call stoponerror
       end if
    end if
    if (.NOT. allocated(siga)) &
         allocate(siga(neg,0:nxs), STAT=ierr)
    if(ierr /= 0)then
       write(*,'("ERROR.  cannot allocate array: ",a,".")') &
            'siga(neg,0:nxs)'
       call stoponerror
    end if

    if (force_alloc) then
       deallocate(sigf, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("ERROR.  cannot deallocate array: ",a,".")') &
               'sigf(neg,0:nxs)'
          call stoponerror
       end if
    end if
    if (.NOT. allocated(sigf)) &
         allocate(sigf(neg,0:nxs), STAT=ierr)
    if(ierr /= 0)then
       write(*,'("ERROR.  cannot allocate array: ",a,".")') &
            'sigf(neg,0:nxs)'
       call stoponerror
    end if

    if (force_alloc) then
       deallocate(sigc, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("ERROR.  cannot deallocate array: ",a,".")') &
               'sigc(neg,0:nxs)'
          call stoponerror
       end if
    end if
    if (.NOT. allocated(sigc)) &
         allocate(sigc(neg,0:nxs), STAT=ierr)
    if(ierr /= 0)then
       write(*,'("ERROR.  cannot allocate array: ",a,".")') &
            'sigc(neg,0:nxs)'
       call stoponerror
    end if

    if (force_alloc) then
       deallocate(scgr, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("ERROR.  cannot deallocate array: ",a,".")') &
               'scgr(neg,0:isct,0:nxs,2)'
          call stoponerror
       end if
    end if
    if (.NOT. allocated(scgr)) &
         allocate(scgr(neg,0:isct,0:nxs,2), STAT=ierr)
    if(ierr /= 0)then
       write(*,'("ERROR.  cannot allocate array: ",a,".")') &
            'scgr(neg,0:isct,0:nxs,2)'
       call stoponerror
    end if

    if (force_alloc) then
       deallocate(rfm, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("ERROR.  cannot deallocate array: ",a,".")') &
               'rfm(0:nitm)'
          call stoponerror
       end if
    end if
    if (.NOT. allocated(rfm)) &
         allocate(rfm(0:nitm), STAT=ierr)
    if(ierr /= 0)then
       write(*,'("ERROR.  cannot allocate array: ",a,".")') &
            'rfm(0:nitm)'
       call stoponerror
    end if

    if (force_alloc) then
       deallocate(dv, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("ERROR.  cannot deallocate array: ",a,".")') &
               'dv(nitm,njtm)'
          call stoponerror
       end if
    end if
    if (.NOT. allocated(dv)) &
         allocate(dv(nitm,njtm), STAT=ierr)
    if(ierr /= 0)then
       write(*,'("ERROR.  cannot allocate array: ",a,".")') &
            'dv(nitm,njtm)'
       call stoponerror
    end if
    
    if (force_alloc) then
       deallocate(zfm, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("ERROR.  cannot deallocate array: ",a,".")') &
               'zfm(0:njtm)'
          call stoponerror
       end if
    end if
    if (.NOT. allocated(zfm)) &
         allocate(zfm(0:njtm), STAT=ierr)
    if(ierr /= 0)then
       write(*,'("ERROR.  cannot allocate array: ",a,".")') &
            'zfm(0:njtm)'
       call stoponerror
    end if

!   if (force_alloc) then
!      deallocate(sar, STAT=ierr)
!      if(ierr /= 0)then
!         write(*,'("ERROR.  cannot deallocate array: ",a,".")') &
!              'sar(nrmx,njtm)'
!         call stoponerror
!      end if
!   end if
!   if (.NOT. allocated(sar)) &
!        allocate(sar(nrmx,njtm), STAT=ierr)
!   if(ierr /= 0)then
!      write(*,'("ERROR.  cannot allocate array: ",a,".")') &
!           'sar(nrmx,njtm)'
!      call stoponerror
!   end if

!   if (force_alloc) then
!      deallocate(sah, STAT=ierr)
!      if(ierr /= 0)then
!         write(*,'("ERROR.  cannot deallocate array: ",a,".")') &
!              'sah(nr,nz)'
!         call stoponerror
!      end if
!   end if
!   if (.NOT. allocated(sah)) &
!        allocate(sah(nr,nz), STAT=ierr)
!   if(ierr /= 0)then
!      write(*,'("ERROR.  cannot allocate array: ",a,".")') &
!           'sah(nr,nz)'
!      call stoponerror
!   end if

    if (force_alloc) then
       deallocate(iindex, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("ERROR.  cannot deallocate array: ",a,".")') &
               'iindex(0:nr)'
          call stoponerror
       end if
    end if
    if (.NOT. allocated(iindex)) &
         allocate(iindex(0:nr), STAT=ierr)
    if(ierr /= 0)then
       write(*,'("ERROR.  cannot allocate array: ",a,".")') &
            'iindex(0:nr)'
       call stoponerror
    end if

    if (force_alloc) then
       deallocate(jindex, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("ERROR.  cannot deallocate array: ",a,".")') &
               'jindex(0:nz)'
          call stoponerror
       end if
    end if
    if (.NOT. allocated(jindex)) &
         allocate(jindex(0:nz), STAT=ierr)
    if(ierr /= 0)then
       write(*,'("ERROR.  cannot allocate array: ",a,".")') &
            'jindex(0:nz)'
       call stoponerror
    end if

    if (force_alloc) then
       deallocate(nsrc, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("ERROR.  cannot deallocate array: ",a,".")') &
               'nsrc(0:nm)'
          call stoponerror
       end if
    end if
    if (.NOT. allocated(nsrc)) &
         allocate(nsrc(0:nm), STAT=ierr)
    if(ierr /= 0)then
       write(*,'("ERROR.  cannot allocate array: ",a,".")') &
            'nsrc(0:nm)'
       call stoponerror
    end if

    if (force_alloc) then
       deallocate(deteff, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("ERROR.  cannot deallocate array: ",a,".")') &
               'deteff(1:neg)'
          call stoponerror
       end if
    end if
    if (.NOT. allocated(deteff)) &
         allocate(deteff(1:neg), STAT=ierr)
    if(ierr /= 0)then
       write(*,'("ERROR.  cannot allocate array: ",a,".")') &
            'deteff(1:neg)'
       call stoponerror
    end if

    if (force_alloc) then
       deallocate(nsrcf, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("ERROR.  cannot deallocate array: ",a,".")') &
               'nsrcf(neg,0:nm)'
          call stoponerror
       end if
    end if
    if (.NOT. allocated(nsrcf)) &
         allocate(nsrcf(neg,0:nm), STAT=ierr)
    if(ierr /= 0)then
       write(*,'("ERROR.  cannot allocate array: ",a,".")') &
            'nsrcf(neg,0:nm)'
       call stoponerror
    end if

    if (force_alloc) then
       deallocate(sfiso, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("ERROR.  cannot deallocate array: ",a,".")') &
               'sfiso(neg,nel)'
          call stoponerror
       end if
    end if
    if (.NOT. allocated(sfiso)) &
         allocate(sfiso(neg,nel), STAT=ierr)
    if(ierr /= 0)then
       write(*,'("ERROR.  cannot allocate array: ",a,".")') &
            'sfiso(neg,nel)'
       call stoponerror
    end if

    if (force_alloc) then
       deallocate(saniso, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("ERROR.  cannot deallocate array: ",a,".")') &
               'saniso(neg,nel)'
          call stoponerror
       end if
    end if
    if (.NOT. allocated(saniso)) &
         allocate(saniso(neg,nel), STAT=ierr)
    if(ierr /= 0)then
       write(*,'("ERROR.  cannot allocate array: ",a,".")') &
            'saniso(neg,nel)'
       call stoponerror
    end if

    if (force_alloc) then
       deallocate(chisrc, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("ERROR.  cannot deallocate array: ",a,".")') &
               'chisrc(neg,nm)'
          call stoponerror
       end if
    end if
    if (.NOT. allocated(chisrc)) &
         allocate(chisrc(neg,nm), STAT=ierr)
    if(ierr /= 0)then
       write(*,'("ERROR.  cannot allocate array: ",a,".")') &
            'chisrc(neg,nm)'
       call stoponerror
    end if

    if (force_alloc) then
       deallocate(fmom, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("ERROR.  cannot deallocate array: ",a,".")') &
               'fmom(neg,0:nmom-1,nitm,njtm)'
          call stoponerror
       end if
    end if
    if (.NOT. allocated(fmom)) &
         allocate(fmom(neg,0:nmom-1,nitm,njtm), STAT=ierr)
    if(ierr /= 0)then
       write(*,'("ERROR.  cannot allocate array: ",a,".")') &
            'fmom(neg,0:nmom-1,nitm,njtm)'
       call stoponerror
    end if

    if (force_alloc) then
       deallocate(amom, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("ERROR.  cannot deallocate array: ",a,".")') &
               'amom(neg,0:nmom-1,nitm,njtm)'
          call stoponerror
       end if
    end if
    if (.NOT. allocated(amom)) &
         allocate(amom(neg,0:nmom-1,nitm,njtm), STAT=ierr)
    if(ierr /= 0)then
       write(*,'("ERROR.  cannot allocate array: ",a,".")') &
            'amom(neg,0:nmom-1,nitm,njtm)'
       call stoponerror
    end if

    if (nrrx.gt.0)then
      if (force_alloc) then
         deallocate(gmom, STAT=ierr)
         if(ierr /= 0)then
            write(*,'("ERROR.  cannot deallocate array: ",a,".")') &
                 'gmom(neg,0:nmom-1,nitm,njtm)'
            call stoponerror
         end if
      end if
      if (.NOT. allocated(gmom)) &
           allocate(gmom(neg,0:nmom-1,nitm,njtm), STAT=ierr)
      if(ierr /= 0)then
         write(*,'("ERROR.  cannot allocate array: ",a,".")') &
              'gmom(neg,0:nmom-1,nitm,njtm)'
         call stoponerror
      end if
    end if

! these are only needed if there are angular fluxes.
  if (iangflux == 1) then
    if (force_alloc) then
       deallocate(afreg, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("ERROR.  cannot deallocate array: ",a,".")') &
               'afreg(neg,ndir,nitm,njtm)'
          call stoponerror
       end if
    end if
    if (.NOT. allocated(afreg)) &
         allocate(afreg(neg,ndir,nitm,njtm), STAT=ierr)
    if(ierr /= 0)then
       write(*,'("ERROR.  cannot allocate array: ",a,".")') &
            'afreg(neg,ndir,nitm,njtm)'
       call stoponerror
    end if

    if (force_alloc) then
       deallocate(afadj, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("ERROR.  cannot deallocate array: ",a,".")') &
               'afadj(neg,ndir,nitm,njtm)'
          call stoponerror
       end if
    end if
    if (.NOT. allocated(afadj)) &
         allocate(afadj(neg,ndir,nitm,njtm), STAT=ierr)
    if(ierr /= 0)then
       write(*,'("ERROR.  cannot allocate array: ",a,".")') &
            'afadj(neg,ndir,nitm,njtm)'
       call stoponerror
    end if

    if (nrrx.gt.0)then
      if (force_alloc) then
         deallocate(afgad, STAT=ierr)
         if(ierr /= 0)then
            write(*,'("ERROR.  cannot deallocate array: ",a,".")') &
                 'afgad(neg,ndir,nitm,njtm)'
            call stoponerror
         end if
      end if
      if (.NOT. allocated(afgad)) &
           allocate(afgad(neg,ndir,nitm,njtm), STAT=ierr)
      if(ierr /= 0)then
         write(*,'("ERROR.  cannot allocate array: ",a,".")') &
              'afgad(neg,ndir,nitm,njtm)'
         call stoponerror
      end if
    end if

  else if(iangflux==0)then
    if (force_alloc) then
       deallocate(afreg, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("ERROR.  cannot deallocate array: ",a,".")') &
               'afreg(1,1,1,1)'
          call stoponerror
       end if
    end if
    if (.NOT. allocated(afreg)) &
         allocate(afreg(1,1,1,1), STAT=ierr)
    if(ierr /= 0)then
       write(*,'("ERROR.  cannot allocate array: ",a,".")') &
            'afreg(1,1,1,1)'
       call stoponerror
    end if

    if (force_alloc) then
       deallocate(afadj, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("ERROR.  cannot deallocate array: ",a,".")') &
               'afadj(1,1,1,1)'
          call stoponerror
       end if
    end if
    if (.NOT. allocated(afadj)) &
         allocate(afadj(1,1,1,1), STAT=ierr)
    if(ierr /= 0)then
       write(*,'("ERROR.  cannot allocate array: ",a,".")') &
            'afadj(1,1,1,1)'
       call stoponerror
    end if

    if (force_alloc) then
       deallocate(afgad, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("ERROR.  cannot deallocate array: ",a,".")') &
               'afgad(1,1,1,1)'
          call stoponerror
       end if
    end if
    if (.NOT. allocated(afgad)) &
         allocate(afgad(1,1,1,1), STAT=ierr)
    if(ierr /= 0)then
       write(*,'("ERROR.  cannot allocate array: ",a,".")') &
            'afgad(1,1,1,1)'
       call stoponerror
    end if

  end if ! iangflux

    if (force_alloc) then
       deallocate(dir, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("ERROR.  cannot deallocate array: ",a,".")') &
               'dir(ndir)'
          call stoponerror
       end if
    end if
    if (.NOT. allocated(dir)) &
         allocate(dir(ndir), STAT=ierr)
    if(ierr /= 0)then
       write(*,'("ERROR.  cannot allocate array: ",a,".")') &
            'dir(ndir)'
       call stoponerror
    end if
 
    if (force_alloc) then
       deallocate(wgt, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("ERROR.  cannot deallocate array: ",a,".")') &
               'wgt(ndir)'
          call stoponerror
       end if
    end if
    if (.NOT. allocated(wgt)) &
         allocate(wgt(ndir), STAT=ierr)
    if(ierr /= 0)then
       write(*,'("ERROR.  cannot allocate array: ",a,".")') &
            'wgt(ndir)'
       call stoponerror
    end if
 
    if (force_alloc) then
       deallocate(eta, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("ERROR.  cannot deallocate array: ",a,".")') &
               'eta(ndir)'
          call stoponerror
       end if
    end if
    if (.NOT. allocated(eta)) &
         allocate(eta(ndir), STAT=ierr)
    if(ierr /= 0)then
       write(*,'("ERROR.  cannot allocate array: ",a,".")') &
            'eta(ndir)'
       call stoponerror
    end if
 
!   if (force_alloc) then
!      deallocate(nadjs, STAT=ierr)
!      if(ierr /= 0)then
!         write(*,'("ERROR.  cannot deallocate array: ",a,".")') &
!              'nadjs(neg,ndir)'
!         call stoponerror
!      end if
!   end if
!   if (.NOT. allocated(nadjs)) &
!        allocate(nadjs(neg,ndir), STAT=ierr)
!   if(ierr /= 0)then
!      write(*,'("ERROR.  cannot allocate array: ",a,".")') &
!           'nadjs(neg,ndir)'
!      call stoponerror
!   end if
 
    if (force_alloc) then
       deallocate(zaidfull, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("ERROR.  cannot deallocate array: ",a,".")') &
               'zaidfull(niso)'
          call stoponerror
       end if
    end if
    if (.NOT. allocated(zaidfull)) &
         allocate(zaidfull(niso), STAT=ierr)
    if(ierr /= 0)then
       write(*,'("ERROR.  cannot allocate array: ",a,".")') &
            'zaidfull(niso)'
       call stoponerror
    end if

    return
  end subroutine allocate_arrays

  subroutine allocate_arrays3()
    !
    ! This subroutine is used to allocate several dynamic arrays
    ! Experimental, called from integrals_r_2d_as to get more memory (not
    ! successful)
    !
    integer :: ierr

    ierr = 0

    if (allocated(fmom)) &
       deallocate(fmom, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("ERROR.  cannot deallocate array: ",a,".")') &
               'fmom(neg,0:nmom-1,nitm,njtm)'
          call stoponerror
       else
          write(*,'("deallocated fmom")')
       end if

    if (allocated(amom)) &
       deallocate(amom, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("ERROR.  cannot deallocate array: ",a,".")') &
               'amom(neg,0:nmom-1,nitm,njtm)'
          call stoponerror
       else
          write(*,'("deallocated amom")')
       end if

    if (allocated(gmom)) &
       deallocate(gmom, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("ERROR.  cannot deallocate array: ",a,".")') &
               'gmom(neg,0:nmom-1,nitm,njtm)'
          call stoponerror
       else
          write(*,'("deallocated gmom")')
       end if

    if (force_alloc) then
       deallocate(omia, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("ERROR.  cannot deallocate array: ",a,".")') &
               'omia(neg,0:max(nr,nz),max(nitm,njtm))'
          call stoponerror
       end if
    end if
    if (.NOT. allocated(omia)) &
         allocate(omia(neg,0:max(nr,nz),max(nitm,njtm)), STAT=ierr)
    if(ierr /= 0)then
       write(*,'("ERROR.  cannot allocate array: ",a,".")') &
            'omia(neg,0:max(nr,nz),max(nitm,njtm))'
       write(*,'("ierror=",i3)')ierr
       call stoponerror
    end if

    if (force_alloc) then
       deallocate(afregas, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("ERROR.  cannot deallocate array: ",a,".")') &
               'afregas(neg,ndir,0:max(nr,nz),max(nitm,njtm))'
          call stoponerror
       end if
    end if
    if (.NOT. allocated(afregas)) &
         allocate(afregas(neg,ndir,0:max(nr,nz),max(nitm,njtm)), STAT=ierr)
    if(ierr /= 0)then
       write(*,'("ERROR.  cannot allocate array: ",a,".")') &
            'afregas(neg,ndir,0:max(nr,nz),max(nitm,njtm))'
       write(*,'("ierror=",i3)')ierr
       call stoponerror
    end if

    if (force_alloc) then
       deallocate(afadjas, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("ERROR.  cannot deallocate array: ",a,".")') &
               'afadjas(neg,ndir,0:max(nr,nz),max(nitm,njtm))'
          call stoponerror
       end if
    end if
    if (.NOT. allocated(afadjas)) &
         allocate(afadjas(neg,ndir,0:max(nr,nz),max(nitm,njtm)), STAT=ierr)
    if(ierr /= 0)then
       write(*,'("ERROR.  cannot allocate array: ",a,".")') &
            'afadjas(neg,ndir,0:max(nr,nz),max(nitm,njtm))'
       write(*,'("ierror=",i3)')ierr
       call stoponerror
    end if

    if (force_alloc) then
       deallocate(forsa, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("ERROR.  cannot deallocate array: ",a,".")') &
               'forsa(neg,nmom,0:max(nr,nz),max(nitm,njtm))'
          call stoponerror
       end if
    end if
    if (.NOT. allocated(forsa)) &
         allocate(forsa(neg,nmom,0:max(nr,nz),max(nitm,njtm)), STAT=ierr)
    if(ierr /= 0)then
       write(*,'("ERROR.  cannot allocate array: ",a,".")') &
            'forsa(neg,nmom,0:max(nr,nz),max(nitm,njtm))'
       write(*,'("ierror=",i3)')ierr
       call stoponerror
    end if

    if (force_alloc) then
       deallocate(adjsa, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("ERROR.  cannot deallocate array: ",a,".")') &
               'adjsa(neg,nmom,0:max(nr,nz),max(nitm,njtm))'
          call stoponerror
       end if
    end if
    if (.NOT. allocated(adjsa)) &
         allocate(adjsa(neg,nmom,0:max(nr,nz),max(nitm,njtm)), STAT=ierr)
    if(ierr /= 0)then
       write(*,'("ERROR.  cannot allocate array: ",a,".")') &
            'adjsa(neg,nmom,0:max(nr,nz),max(nitm,njtm))'
       write(*,'("ierror=",i3)')ierr
       call stoponerror
    end if

    return
  end subroutine allocate_arrays3
end module COMS
