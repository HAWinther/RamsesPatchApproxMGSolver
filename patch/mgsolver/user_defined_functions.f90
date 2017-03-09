!########################################################################
!                                                                       # 
! Fast approximate method to simulate structure formation in            #
! modified gravity models with screening where the screening depends    #
! on the value of the newtonian potential (chamleon models).            #
! Such models can be most easily described in the m(a), beta(a)         #
! formatlism of Brax et al.                                             #
!                                                                       # 
! See Winther & Ferreira for details about the method                   #
!                                                                       # 
! Implement the three functions scalar_coupling, scalar_mass_term, and  #
! critical_potential below or use the define FOFR or SYMMETRON for the  #
! Hu-Sawicky f(R) model / Symmetron model                               #
!                                                                       # 
! Hans A. Winther (2015) University of Oxford                           #
!                                                                       # 
!########################################################################

!########################################################################
! The function a^2 m^2(a) B0^2                                          #
!########################################################################

function calc_scalar_mass_term(anow)
  use poisson_parameters
  use poisson_commons
  implicit none
  real(dp) :: calc_scalar_mass_term
  real(dp) :: mass, anow, RR, RR0

  ! The mass m(a) (or m(a)/H0 to be accurate)
#if defined(SYMMETRON)

  real(dp) :: scalar_mass_term_symmetron 
  mass = scalar_mass_term_symmetron(anow)

#elif defined(FOFR)

  real(dp) :: scalar_mass_term_fofr
  mass = scalar_mass_term_fofr(anow)

#else

  real(dp) :: scalar_mass_term_mymodel
  mass = scalar_mass_term_mymodel(anow)

#endif

  ! Go from (m/H0) ==> (m B0 * a)^2
  calc_scalar_mass_term = (anow * mass * boxlen_ini/2998.0)**2

end function

!########################################################################
! The function 2beta(a)^2                                               #
!########################################################################

function calc_scalar_coupling(anow)
  use poisson_parameters
  use poisson_commons
  implicit none
  real(dp) :: calc_scalar_coupling
  real(dp) :: anow, beta

  ! The coupling beta(a)
#if defined(SYMMETRON)

  real(dp) :: scalar_coupling_symmetron
  beta = scalar_coupling_symmetron(anow)

#elif defined(FOFR)

  real(dp) :: scalar_coupling_fofr
  beta = scalar_coupling_fofr(anow)

#else

  real(dp) :: scalar_coupling_mymodel
  beta = scalar_coupling_mymodel(anow)

#endif

  ! Transform to 2beta(a)^2
  calc_scalar_coupling = 2.0 * beta * beta

end function

!########################################################################
! The screening function dR/R(Phi)                                      #
!########################################################################

function calc_screening(phi_newton_now, rho_now, anow)
  use poisson_parameters
  use poisson_commons
  implicit none
  integer  :: icell_amr
  real(dp) :: calc_screening, anow, rho_now
  real(dp) :: phi_newton, phi_newton_now, phi_crit
 
  ! The critical potential for screening
#if defined(SYMMETRON)

  real(dp) :: phi_crit_symmetron
  phi_crit = phi_crit_symmetron(anow)

#elif defined(FOFR)

  real(dp) :: phi_crit_fofr
  phi_crit = phi_crit_fofr(anow)

#else

  real(dp) :: phi_crit_mymodel
  phi_crit = phi_crit_mymodel(anow)

#endif

  ! Newtonian potential
  phi_newton = phi_newton_now * (boxlen_ini/2998.0/anow)**2

  ! No screening in voids (PhiNewton > 0)
  if(phi_newton >= 0.0d0) then
    calc_screening = 1.0d0
    return
  endif

  ! Calculate screening condition
  phi_crit = abs(phi_crit/phi_newton)
  if(phi_crit < 1.0d0) then
    calc_screening = phi_crit
  else
    calc_screening = 1.0d0
  endif

end function

!###################################################################
! Functions for Symmetron gravity                                  #
!###################################################################

#if defined(SYMMETRON)

function phi_crit_symmetron(anow)
  use amr_commons
  use poisson_parameters
  use poisson_commons
  implicit none
  real(dp) :: anow, phi_crit_symmetron

  ! Phi_critical for the symmetron
  phi_crit_symmetron = 3.0d0*omega_m/assb**3 * (anow * lambda_0)**2

  return
end function

function scalar_mass_term_symmetron(anow)
  use poisson_parameters
  use poisson_commons
  implicit none
  real(dp) :: scalar_mass_term_symmetron, anow

  ! m(a)/H0 for the symmetron
  if(anow>assb)then
    scalar_mass_term_symmetron = 2998.0/(lambda_0*boxlen_ini) * sqrt(1.0-(assb/anow)**3) 
  else
    scalar_mass_term_symmetron = 0.0d0
  endif

  return
end function

function scalar_coupling_symmetron(anow)
  use poisson_parameters
  use poisson_commons
  implicit none
  real(dp) :: scalar_coupling_symmetron
  real(dp) :: beta, anow

  ! beta(a) for the symmetron
  if(anow>assb)then
    scalar_coupling_symmetron = beta_0 * sqrt(1.0-(assb/anow)**3)
  else
    scalar_coupling_symmetron = 0.0d0
  endif

  return
end function

#elif defined(FOFR)

!###################################################################
! Functions for f(R) gravity                                       #
!###################################################################

function phi_crit_fofr(anow)
  use amr_commons
  use poisson_parameters
  use poisson_commons
  implicit none
  real(dp) :: anow, R0overR, phi_crit_fofr

  ! Phi_critical for Hu-Sawicky f(R)
  R0overR  = (omega_m + 4.0 * omega_l)/(omega_m/anow**3 + 4.0 * omega_l)
  phi_crit_fofr = 1.5d0 * fofr0 * (R0overR)**(n_fofr+1.0)

  return
end function

function scalar_mass_term_fofr(anow)
  use poisson_parameters
  use poisson_commons
  implicit none
  real(dp) :: scalar_mass_term_fofr
  real(dp) :: RR, RR0, anow

  ! m(a)/H0 for Hu-Sawicky f(R)
  RR   = 3.0 * ( omega_m/anow**3 + 4.0*omega_l )
  RR0  = 3.0 * ( omega_m + 4.0 * omega_l ) 
  scalar_mass_term_fofr = RR0/(3.0*(n_fofr+1.0)*fofr0) * (RR/RR0)**(n_fofr + 2.0)
  scalar_mass_term_fofr = sqrt(scalar_mass_term_fofr)

  return
end function

function scalar_coupling_fofr(anow)
  use poisson_parameters
  use poisson_commons
  implicit none
  real(dp) :: scalar_coupling_fofr, anow

  ! beta(a) for Hu-Sawicky f(R)
  scalar_coupling_fofr = 1.0/sqrt(6.0) 

  return
end function

#else

!###################################################################
! Functions for your model                                         #
!###################################################################

function phi_crit_mymodel(anow)
  use amr_commons
  use poisson_parameters
  use poisson_commons
  implicit none
  real(dp) :: anow, R0overR, phi_crit_mymodel

  ! Phi_critical for your model
  phi_crit_mymodel = 1.0

  return
end function

function scalar_mass_term_mymodel(anow)
  use poisson_parameters
  use poisson_commons
  implicit none
  real(dp) :: scalar_mass_term_mymodel
  real(dp) :: RR, RR0, anow

  ! m(a)/H0 for your model
  scalar_mass_term_mymodel = 0.0

  return
end function

function scalar_coupling_mymodel(anow)
  use poisson_parameters
  use poisson_commons
  implicit none
  real(dp) :: scalar_coupling_mymodel, anow

  ! beta(a) for your model
  scalar_coupling_mymodel = 1.0

  return
end function

#endif

