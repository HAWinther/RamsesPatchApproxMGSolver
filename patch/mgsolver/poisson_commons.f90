module poisson_commons 
  use amr_commons
  use poisson_parameters

#if defined(MGSOLVER)

  ! Standard gravity arrays
  real(dp),allocatable,dimension(:),target  ::phi_array, phi_old_array    ! Potential
  real(dp),allocatable,dimension(:,:),target::f_array                     ! 3-force
  real(dp),allocatable,dimension(:)  ::rho                                ! Density

  ! Modified gravity arrays
  real(dp),allocatable,dimension(:),target  ::phi_mg_array, phi_mg_old_array ! Potential
  real(dp),allocatable,dimension(:,:),target::f_mg_array                     ! 3-force

  ! Pointers for allowing to use the same methods for both equations
  real(dp),pointer :: phi(:),phi_old(:),f(:,:)
  real(dp),pointer :: phi_mg(:),phi_mg_old(:),f_mg(:,:)

#else

  real(dp),allocatable,dimension(:)  ::phi,phi_old       ! Potential
  real(dp),allocatable,dimension(:)  ::rho               ! Density
  real(dp),allocatable,dimension(:,:)::f                 ! 3-force

#endif

  real(dp),allocatable,dimension(:)  ::rho_top   ! Density at last CIC level                                 

  ! Multigrid lookup table for amr -> mg index mapping
  integer, allocatable, dimension(:) :: lookup_mg   ! Lookup table

  ! Communicator arrays for multigrid levels
  type(communicator), allocatable, dimension(:,:) :: active_mg
  type(communicator), allocatable, dimension(:,:) :: emission_mg

  ! Minimum MG level
  integer :: levelmin_mg

  ! Multigrid safety switch
  logical, allocatable, dimension(:) :: safe_mode

  ! Multipole coefficients
  real(dp),dimension(1:ndim+1)::multipole

end module poisson_commons
