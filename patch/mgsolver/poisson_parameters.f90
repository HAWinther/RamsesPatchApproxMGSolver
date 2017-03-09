module poisson_parameters
  use amr_parameters

#if defined(MGSOLVER)

  ! Turn on solver. If off then just standard RAMSES
  logical :: mg_solver_active = .false.

  ! Flag for using the Poisson solver for the MG scalar field equation
  logical :: use_poisson_solver_for_mg = .false.

  ! Use Winther & Ferreira 2014 screening method
  logical :: include_screening = .false.

  ! Convergence criterion for modified gravity Poisson solver
  real(dp) :: epsilon_mg = 1.0D-6
  
  ! Parameters for your model
  real(dp), dimension(1:5) :: myparam = 0.0d0

#if defined(SYMMETRON)
  
  ! Parameters for symmetron (Range in Mpc/h, coupling, scale-factor for symmetry breaking)
  real(dp) :: lambda_0 = 1.0
  real(dp) :: beta_0   = 1.0
  real(dp) :: assb     = 0.5

#elif defined(FOFR)

  ! Parameters for Hu-Sawicky f(R) = R + R0 * fR0/n * (R0/R)^(n)
  real(dp) :: n_fofr  = 1.0
  real(dp) :: fofr0   = 1.0d-5

#endif

  ! Flag to keep track of if we have started solving
  ! For high a the range / coupling is often so small that solving is pointless
  ! We start solving when |Geff(kmax,a)/G - 1| > eps_to_start_solving
  logical :: mg_solver_has_started = .false.
  real(dp) :: eps_to_start_solving = 0.001
#endif

  ! Convergence criterion for Poisson solvers
  real(dp)::epsilon=1.0D-4

  ! Type of force computation
  integer ::gravity_type=0

  ! Gravity parameters
  real(dp),dimension(1:10)::gravity_params=0.0

  ! Maximum level for CIC dark matter interpolation
  integer :: cic_levelmax=0

  ! Min level for CG solver
  ! level < cg_levelmin uses fine multigrid
  ! level >=cg_levelmin uses conjugate gradient
  integer :: cg_levelmin=999

  ! Gauss-Seidel smoothing sweeps for fine multigrid
  integer, parameter :: ngs_fine   = 2
  integer, parameter :: ngs_coarse = 2

  ! Number of multigrid cycles for coarse levels *in safe mode*
  !   1 is the fastest,
  !   2 is slower but can give much better convergence in some cases
  integer, parameter :: ncycles_coarse_safe = 1

end module poisson_parameters
