!===============================================================================
!   Name    : RESEVOIRS 
!   Function: A (multi-)reservoir analysis of EMC3-EIRENE simulation results	
!   Contains: RESERVOIRS
!   Subroutines  : 	multi_reservoir_analysis(ZONE_TYPE, Lc)
!			setup_reservoirs(isetup, Lcmax)
!			run_analysis()
!  
!   Authors: Heinke Frerichs --hfrerichs@wisc.edu 
!            Ian Waters --iwaters@wisc.edu
!   GitHub Page: https://github.com/aeschylus314/RESERVOIRS 
! 
!   TODO: Finish Zone Type Definition of Reservoirs
!===============================================================================

module RESERVOIRS
  use PHYSICAL_CELL
  implicit none
  private

  ! definition of reservoirs
  integer, parameter :: &
     ZONE_TYPE         = 1, &
     CONNECTION_LENGTH = 2

  ! map plasma cell ic to reservoir IRES(ic)
  integer, dimension(:), allocatable :: IRES

  integer :: resn  ! total number of reservoirs
  integer, parameter :: iu = 82   ! Dummy variable for reading in datafile

  real*8, dimension(:), allocatable :: SRES, NRES ! SRES(i) = source in the ith reservoir
                                                  ! NRES(i) = # of particles in the ith reservoir
  real*8, dimension(:), allocatable :: Lc_cell    ! Lc_cell(i) = Connection length of ith cell
  public :: setup_reservoirs, run_analysis        ! Subroutines
  contains

  !---------------------------------------------------------------------
  ! Sets up definition of reservoirs: 1=Confined Region, 2=SOL
  !---------------------------------------------------------------------
  subroutine setup_reservoirs(isetup, Lcmax)
  IMPLICIT NONE
  integer, intent(in) :: isetup
  real*8, intent(in) :: Lcmax 

  integer :: ic, iz

  if (allocated(IRES)) deallocate(IRES)
  allocate (IRES(NC_PL));  IRES = -1

  select case (isetup)
   case (ZONE_TYPE)
     resn = 2
     do ic=1,NC_PL
        iz = IZCELL(ic)
        
        !Currently Set Up for a diverted tokamak with 1 toroidal zone. Come back to. 
        if (iz == 0) then
            IRES(ic) = 1
        else
            IRES(ic) = 2
        end if
     end do

   case (CONNECTION_LENGTH)
     resn = 2
     
     ! Loads the connection length for each cell from the connection length file
     allocate (Lc_cell(NC_PL))
     open  (iu, file='CONNECTION_LENGTH')
     read  (iu, *) Lc_cell
     close (iu)

     do ic=1,NC_PL
        if (Lc_cell(ic) >= Lcmax) then
            IRES(ic) = 1
        else
            IRES(ic) = 2
        end if
     end do

   case default
     write (6, *) 'error in subroutine setup_reservoires, invalid isetup = ', isetup
     stop
  end select

  end subroutine setup_reservoirs
  !---------------------------------------------------------------------
  ! Adds up the source and total particle number in each reservoir
  !---------------------------------------------------------------------
  subroutine run_analysis()
  use PLASMA_PARAM
  use BOUNDARY_COND
  use SOURCE_V_PL
  IMPLICIT NONE

  real*8  :: Sp       ! Normalized source density in a cell
  integer :: ic, i    ! ic = cell #, i = reservoir #

  if (allocated(SRES)) deallocate(SRES)
  allocate (SRES(0:resn));  SRES = 0.d0
  if (allocated(NRES)) deallocate(NRES)
  allocate (NRES(0:resn));  NRES = 0.d0

  ! core domain
  SRES(0) = PFLUX_B_SF(1,0)  +  SP_MAIN_OLD * (PFLUX_TOTAL(1) - PFLUX_B_SF(1,0))

  ! edge domain
  do ic=1, NC_PL
     i = IRES(ic)
     if (i < 1  .or.  i > resn) then
         cycle
     else
         Sp      = VSOUP0(ic,1) * (PFLUX_TOTAL(1) - PFLUX_B_SF(1,0))
         SRES(i) = SRES(i) + Sp          * VOLCEL(ic)
         NRES(i) = NRES(i) + DENS0(ic,1) * VOLCEL(ic)
     end if
  end do

  ! display results
  write (6, *) "reservoir      source strength [A]    # particles"
  write (6, *) 0, SRES(0)
  do i=1,resn
     write (6, *) i, SRES(i), NRES(i)
  end do

  end subroutine run_analysis
  !---------------------------------------------------------------------

end module RESERVOIRS
!===============================================================================

!===============================================================================
subroutine multi_reservoir_analysis(ZONE_TYPE, Lc)
  use RESERVOIRS
  implicit none
  integer, intent(in) :: ZONE_TYPE ! 1 for zone separation, 2 for separation by 
                                   ! connection length
  real*8, intent(in) :: Lc         ! Max Connection Length

  call setup_reservoirs(ZONE_TYPE, Lc)
  call run_analysis()

end subroutine multi_reservoir_analysis
!===============================================================================
