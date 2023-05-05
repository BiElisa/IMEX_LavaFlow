!********************************************************************************
!> \brief Constitutive equations
!********************************************************************************
MODULE constitutive_2d

  USE parameters_2d, ONLY : wp, sp ,tolh
  USE parameters_2d, ONLY : n_eqns , n_vars , n_solid
  USE parameters_2d, ONLY : rheology_flag , rheology_model , energy_flag ,      &
       liquid_flag , gas_flag , alpha_flag , slope_correction_flag ,            &
       curvature_term_flag
  USE parameters_2d, ONLY : temperature_profile_flag ! EB : add       

  IMPLICIT none

  !> flag used for size of implicit non linear-system
  LOGICAL, ALLOCATABLE :: implicit_flag(:)

  !> map from implicit variables to original variables 
  INTEGER, ALLOCATABLE :: implicit_map(:)
  
  !> gravitational acceleration 
  REAL(wp) :: grav
  REAL(wp) :: inv_grav

  !> drag coefficients (Voellmy-Salm model)
  REAL(wp) :: mu
  REAL(wp) :: xi

  !> drag coefficients (B&W model)
  REAL(wp) :: friction_factor

  !> drag coefficients (plastic model)
  REAL(wp) :: tau
  
  !> Velocity profile parameter  
  REAL*8 :: beta_vel   ! EB : added

  !--- Temperature parameters
  
  !> Temperature profile parameters
  REAL*8 :: zeta , eta , theta , beta_T ! EB : added

  !> evironment temperature [K]
  REAL(wp) :: T_env

  !> radiative coefficient
  REAL(wp) :: rad_coeff

  !> friction coefficient
  REAL(wp) :: frict_coeff

  !> reference temperature [K]
  REAL(wp) :: T_ref

  !> reference kinematic viscosity [m2/s]
  REAL(wp) :: nu_ref

  !> viscosity parameter [K-1] (b in Table 1 Costa & Macedonio, 2005)
  REAL(wp) :: visc_par

  !> yield strength for lava rheology [kg m-1 s-2] (Eq.4 Kelfoun & Varga, 2015)
  REAL(wp) :: tau0

  !> fluid specific heat [J kg-1 K-1]
  REAL(wp) :: c_p

  !> coefficient for the convective term in the temperature equation for lava
  REAL(wp) :: convective_term_coeff
  
  !> atmospheric heat trasnfer coefficient [W m-2 K-1] (lambda in C&M, 2005)
  REAL(wp) :: atm_heat_transf_coeff

  !> fractional area of the exposed inner core (f in C&M, 2005)
  REAL(wp) :: exp_area_fract

  !> coefficient for the radiative term in the temperature equation for lava
  REAL(wp) :: radiative_term_coeff

  !> coefficient for the conductive term in the temperature equation for lava
  REAL(wp) :: conductive_term_coeff

  !> Stephan-Boltzmann constant [W m-2 K-4]
  REAL(wp), PARAMETER :: SBconst = 5.67E-8_wp

  !> emissivity (eps in Costa & Macedonio, 2005)
  REAL(wp) :: emissivity
  
  !> temperature of lava-ground interface (T_gr in Costa & Macedonio, 2005)
  REAL(wp) :: T_ground

  !> thermal boundary layer fraction of total thickness
  REAL(wp) :: enne
  
  ! EB : added
  ! Articolo 1 : T_max.   (for Paper [B,dMV,DB. Modiffied SW model. AMM, 2021.] )
  ! Articolo 2 : thermal_conductivity_fluid , thermal_conductivity_soil , emme , 
  !              rho_soil , c_p_soil , T_soil , zeta, eta , theta, beta_T.
  
  !> Maximum temperature (if piecewise linear temperature profile and
  !> none thermal exchange)
  REAL*8 :: T_max
  
  !> thermal conductivity of fluid [W m-1 K-1] (k in Costa & Macedonio, 2005)
  REAL*8 :: thermal_conductivity_fluid
  
  !> thermal conductivity of soil [W m-1 K-1] 
  REAL*8 :: thermal_conductivity_soil
  
  !> thermal boundary layer fraction of the soil, thickness: h*emme
  REAL(wp) :: emme ! (EB forse rinominare come emme_thermBLS)

  !> density of soil
  REAL*8 :: rho_soil
  
  !> specific heat of soil
  REAL*8 :: c_p_soil
  
  !> temperature of deep soil, does not change
  REAL*8 :: T_soil

  ! EB : end added

  !> Specific heat of carrier phase (gas or liquid)
  REAL(wp) :: sp_heat_c  ! ( initialized from input)   

  !> Density of carrier phase in substrate ( units: kg m-3 )
  REAL(wp) :: rho_c_sub
  
  !> Ambient density of air ( units: kg m-3 )
  REAL(wp) :: rho_a_amb

  !> Specific heat of air (units: J K-1 kg-1)
  REAL(wp) :: sp_heat_a

  !> Specific gas constant of air (units: J kg-1 K-1)
  REAL(wp) :: sp_gas_const_a

  !> Kinematic viscosity of air (units: m2 s-1)
  REAL(wp) :: kin_visc_a

  !> Kinematic viscosity of liquid (units: m2 s-1)
  REAL(wp) :: kin_visc_l

  !> Kinematic viscosity of carrier phase (units: m2 s-1)
  REAL(wp) :: kin_visc_c

  !> Temperature of ambient air (units: K)
  REAL(wp) :: T_ambient

  !> Density of sediments ( units: kg m-3 )
  REAL(wp), ALLOCATABLE :: rho_s(:)

  !> Reciprocal of density of sediments ( units: kg m-3 )
  REAL(wp), ALLOCATABLE :: inv_rho_s(:)

  !> Reciprocal of density of sediments ( units: kg m-3 )
  COMPLEX(wp), ALLOCATABLE :: c_inv_rho_s(:)

  !> Diameter of sediments ( units: m )
  REAL(wp), ALLOCATABLE :: diam_s(:)

  !> Specific heat of solids (units: J K-1 kg-1)
  REAL(wp), ALLOCATABLE :: sp_heat_s(:)

  !> ambient pressure (units: Pa)
  REAL(wp) :: pres

  !> reciprocal of ambient pressure (units: Pa)
  REAL(wp) :: inv_pres

  !> liquid density (units: kg m-3)
  REAL(wp) :: rho_l

  !> reciprocal of liquid density (units: kg m-3)
  REAL(wp) :: inv_rho_l

  !> Sepcific heat of liquid (units: J K-1 kg-1)
  REAL(wp) :: sp_heat_l

CONTAINS

  !******************************************************************************
  !> \brief Initialization of relaxation flags
  !
  !> This subroutine set the number and the flags of the non-hyperbolic
  !> terms.
  !> \date 07/09/2012 EB : forse cambiare intro: numero termini non iperbolici da trattare implicitamente
  !******************************************************************************

  SUBROUTINE init_problem_param

    USE parameters_2d, ONLY : n_nh
    IMPLICIT NONE

    integer :: i,j

    ALLOCATE( implicit_flag(n_eqns) )

    implicit_flag(1:n_eqns) = .FALSE.
    implicit_flag(2) = .TRUE.
    implicit_flag(3) = .TRUE.

    ! Temperature
    IF ( rheology_model .EQ. 3 ) THEN

       implicit_flag(4) = .TRUE.

    END IF
       
    ! Solid volume fraction
    implicit_flag(5:4+n_solid) = .FALSE.

    n_nh = COUNT( implicit_flag )

    ALLOCATE( implicit_map(n_nh) )

    j=0
    DO i=1,n_eqns

       IF ( implicit_flag(i) ) THEN

          j=j+1
          implicit_map(j) = i

       END IF

    END DO

    WRITE(*,*) 'Implicit equations =',n_nh
    
    RETURN
    
  END SUBROUTINE init_problem_param

  !******************************************************************************
  !> \brief Physical variables
  !
  !> This subroutine evaluates from the conservative local variables qj
  !> the local physical variables  (\f$h,u,v,\alpha_s,\rho_m,T,\alpha_l \f$).
  !> \param[in]    r_qj        real conservative variables 
  !> \param[out]   r_h         real-value flow thickness 
  !> \param[out]   r_u         real-value flow x-velocity 
  !> \param[out]   r_v         real-value flow y-velocity
  !> \param[out]   r_alphas    real-value solid volume fractions
  !> \param[out]   r_rho_m     real-value flow density
  !> \param[out]   r_T         real-value flow temperature 
  !> \param[out]   r_alphal    real-value liquid volume fraction
  !> \param[out]   r_red_grav  real-value reduced gravity
  !
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !> \date 2019/12/13
  !******************************************************************************

  SUBROUTINE r_phys_var(r_qj , r_h , r_u , r_v , r_alphas , r_rho_m , r_T ,     &
       r_alphal , r_red_grav )

    USE parameters_2d, ONLY : eps_sing , eps_sing4
    IMPLICIT none

    REAL(wp), INTENT(IN) :: r_qj(n_vars)       !< real-value conservative var
    REAL(wp), INTENT(OUT) :: r_h               !< real-value flow thickness
    REAL(wp), INTENT(OUT) :: r_u               !< real-value x-velocity
    REAL(wp), INTENT(OUT) :: r_v               !< real-value y-velocity
    REAL(wp), INTENT(OUT) :: r_alphas(n_solid) !< real-value solid volume fracts
    REAL(wp), INTENT(OUT) :: r_rho_m           !< real-value mixture density
    REAL(wp), INTENT(OUT) :: r_T               !< real-value temperature
    REAL(wp), INTENT(OUT) :: r_alphal          !< real-value liquid volume fract
    REAL(wp), INTENT(OUT) :: r_red_grav        !< real-value reduced gravity

    REAL(wp) :: r_inv_rhom
    REAL(wp) :: r_xs(n_solid)     !< real-value solid mass fraction
    REAL(wp) :: r_xs_tot

    REAL(wp) :: r_Ri            !< real-value Richardson number
    REAL(wp) :: r_xl            !< real-value liquid mass fraction
    REAL(wp) :: r_xc            !< real-value carrier phase mass fraction
    REAL(wp) :: r_alphac        !< real-value carrier phase volume fraction
    REAL(wp) :: r_sp_heat_mix   !< Specific heat of mixture
    REAL(wp) :: r_rho_c         !< real-value carrier phase density [kg/m3]
    REAL(wp) :: r_inv_rho_c

    REAL(wp) :: inv_qj1
 
    ! compute solid mass fractions
    ! IF ( r_qj(1) .GT. eps_sing ) THEN
    IF ( r_qj(1) .GT. 0.0_wp ) THEN

       inv_qj1 = 1.0_wp / r_qj(1)

       r_xs(1:n_solid) = r_qj(5:4+n_solid) * inv_qj1

       IF ( SUM( r_qj(5:4+n_solid) ) .EQ. r_qj(1) ) THEN

          r_xs(1:n_solid) = r_xs(1:n_solid) / SUM( r_xs(1:n_solid) )

       END IF
       
    ELSE

       r_xs(1:n_solid) = 0.0_wp

    END IF

    r_xs_tot = SUM(r_xs)

    IF ( gas_flag .AND. liquid_flag ) THEN

       ! compute liquid mass fraction
       IF ( r_qj(1) .GT. eps_sing ) THEN

          r_xl = r_qj(n_vars) * inv_qj1

       ELSE

          r_xl = 0.0_wp

       END IF

       ! compute carrier phase (gas) mass fraction
       r_xc =  1.0_wp - r_xs_tot - r_xl

       ! specific heat of the mixutre: mass average of sp. heat pf phases
       r_sp_heat_mix = DOT_PRODUCT( r_xs(1:n_solid) , sp_heat_s(1:n_solid) )    &
            + r_xl * sp_heat_l + r_xc * sp_heat_c

    ELSE

       ! compute carrier phase (gas or liquid) mass fraction
       r_xc = 1.0_wp - r_xs_tot

       ! specific heaf of the mixutre: mass average of sp. heat pf phases
       r_sp_heat_mix = DOT_PRODUCT( r_xs(1:n_solid) , sp_heat_s(1:n_solid) )    &
            + r_xc * sp_heat_c

    END IF

    ! compute temperature from energy
    IF ( r_qj(1) .GT. eps_sing ) THEN

       IF ( energy_flag ) THEN

          r_T = ( r_qj(4) - 0.5_wp * ( r_qj(2)**2 + r_qj(3)**2 ) * inv_qj1 ) /  &
               ( r_qj(1) * r_sp_heat_mix ) 

       ELSE

          r_T = r_qj(4) / ( r_qj(1) * r_sp_heat_mix ) 

       END IF

       IF ( r_T .LE. 0.0_wp ) r_T = T_ambient

    ELSE

       r_T = T_ambient

    END IF

    IF ( gas_flag ) THEN

       ! carrier phase is gas
       r_rho_c =  pres / ( sp_gas_const_a * r_T )
       r_inv_rho_c = sp_gas_const_a * r_T * inv_pres
       sp_heat_c = sp_heat_a

    ELSE

       r_rho_c = rho_l
       r_inv_rho_c = inv_rho_l
       sp_heat_c = sp_heat_l

    END IF

    r_inv_rhom = DOT_PRODUCT( r_xs(1:n_solid) , inv_rho_s(1:n_solid) )       &
         + r_xc * r_inv_rho_c

    IF ( gas_flag .AND. liquid_flag ) THEN

       r_inv_rhom = r_inv_rhom + r_xl * inv_rho_l
       r_alphal = r_xl * r_rho_m * inv_rho_l

    END IF

    ! mixture density
    r_rho_m = 1.0_wp / r_inv_rhom

    ! convert from mass fraction to volume fraction
    r_alphas(1:n_solid) = r_xs(1:n_solid) * r_rho_m * inv_rho_s(1:n_solid)

!!$    IF ( ( r_qj(1) .GT. 0.0_wp ) .AND. ( ABS( r_alphas(1) - 0.1_wp ) .GT. 1.0E-5_wp ) ) THEN
!!$
!!$       WRITE(*,*) 'qui',r_alphas(1) - 0.1_wp,r_qj(1),r_qj(5)
!!$       READ(*,*)
!!$
!!$    END IF
    
    ! convert from mass fraction to volume fraction
    r_alphac = r_xc * r_rho_m * r_inv_rho_c

    r_h = r_qj(1) * r_inv_rhom

    ! reduced gravity
    r_red_grav = ( r_rho_m - rho_a_amb ) * r_inv_rhom * grav

    ! velocity components
    IF ( r_qj(1) .GT. eps_sing ) THEN

       r_u = r_qj(2) * inv_qj1
       r_v = r_qj(3) * inv_qj1

    ELSE

       r_u = SQRT(2.0_wp) * r_qj(1) * r_qj(2) / SQRT( r_qj(1)**4 + eps_sing4 )
       r_v = SQRT(2.0_wp) * r_qj(1) * r_qj(3) / SQRT( r_qj(1)**4 + eps_sing4 )

    END IF

    ! Richardson number
    IF ( ( r_u**2 + r_v**2 ) .GT. 0.0_wp ) THEN

       r_Ri = r_red_grav * r_h / ( r_u**2 + r_v**2 )

    ELSE

       r_Ri = 0.0_wp

    END IF

    RETURN

  END SUBROUTINE r_phys_var

  !******************************************************************************
  !> \brief Physical variables
  !
  !> This subroutine evaluates from the conservative local variables qj
  !> the local physical variables  (\f$h,u,v,T,\rho_m,red grav,\alpha_s \f$).
  !> \param[in]    c_qj      complex conservative variables 
  !> \param[out]   h         complex-value flow thickness 
  !> \param[out]   u         complex-value flow x-velocity 
  !> \param[out]   v         complex-value flow y-velocity
  !> \param[out]   T         complex-value flow temperature 
  !> \param[out]   rho_m     complex-value flow density
  !> \param[out]   alphas    complex-value solid volume fractions
  !> \param[out]   inv_rhom  complex-value mixture density reciprocal
  !
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !> \date 2019/12/13
  !******************************************************************************

  SUBROUTINE c_phys_var( c_qj , h , u , v , T , rho_m , alphas , inv_rhom )

    USE COMPLEXIFY
    USE parameters_2d, ONLY : eps_sing , eps_sing4
    IMPLICIT none

    COMPLEX(wp), INTENT(IN) :: c_qj(n_vars)
    COMPLEX(wp), INTENT(OUT) :: h               !< height [m]
    COMPLEX(wp), INTENT(OUT) :: u               !< velocity (x direction) [m s-1]
    COMPLEX(wp), INTENT(OUT) :: v               !< velocity (y direction) [m s-1]
    COMPLEX(wp), INTENT(OUT) :: T               !< temperature [K]
    COMPLEX(wp), INTENT(OUT) :: rho_m           !< mixture density [kg m-3]
    COMPLEX(wp), INTENT(OUT) :: alphas(n_solid) !< sediment volume fractions
    COMPLEX(wp), INTENT(OUT) :: inv_rhom       !< 1/mixture density [kg-1 m3]

    COMPLEX(wp) :: xs(n_solid)             !< sediment mass fractions
    COMPLEX(wp) :: xs_tot                  !< sum of solid mass fraction
    COMPLEX(wp) :: xl                      !< liquid mass fraction
    COMPLEX(wp) :: xc                      !< carrier phase mass fraction
    COMPLEX(wp) :: sp_heat_mix             !< Specific heat of mixture
    COMPLEX(wp) :: inv_cqj1                !< reciprocal of 1st cons. variable
    COMPLEX(wp) :: inv_rho_c               !< carrier phase reciprocal

    ! compute solid mass fractions
    ! IF ( REAL(c_qj(1)) .GT. eps_sing ) THEN
    IF ( REAL(c_qj(1)) .GT. 0.0_wp ) THEN

       inv_cqj1 = 1.0_wp / c_qj(1)

       xs(1:n_solid) = c_qj(5:4+n_solid) * inv_cqj1
       
    ELSE

       inv_cqj1 = CMPLX(0.0_wp,0.0_wp,wp)
       xs(1:n_solid) = CMPLX(0.0_wp,0.0_wp,wp)

    END IF

    xs_tot = SUM(xs)

    ! compute carrier phase (gas or liquid) mass fraction
    xc = CMPLX(1.0_wp,0.0_wp,wp) - xs_tot

    ! specific heaf of the mixutre: mass average of sp. heat pf phases
    sp_heat_mix = DOT_PRODUCT( xs(1:n_solid) , sp_heat_s(1:n_solid) )           &
         + xc * sp_heat_c

    IF ( gas_flag .AND. liquid_flag ) THEN

       ! compute liquid mass fraction
       xl = c_qj(n_vars) * inv_cqj1

       ! compute carrier phase (gas) mass fraction
       xc =  xc - xl

       ! specific heaf of the mixutre: mass average of sp. heat pf phases
       sp_heat_mix = sp_heat_mix + xl * sp_heat_l

    END IF

    ! compute temperature from energy
    IF ( REAL(c_qj(1)) .GT. eps_sing ) THEN

       IF ( energy_flag ) THEN

          T = ( c_qj(4) - 0.5_wp * ( c_qj(2)**2 + c_qj(3)**2 ) * inv_cqj1 ) /   &
               ( c_qj(1) * sp_heat_mix ) 

       ELSE

          T = c_qj(4) / ( c_qj(1) * sp_heat_mix ) 
          
       END IF

       IF ( REAL(T) .LE. 0.0_wp ) T = CMPLX(T_ambient,0.0_wp,wp)

    ELSE

       T = CMPLX(T_ambient,0.0_wp,wp)

    END IF

    IF ( gas_flag ) THEN

       ! carrier phase is gas
       inv_rho_c = sp_gas_const_a * T * inv_pres

    ELSE

       inv_rho_c = CMPLX(inv_rho_l,0.0_wp,wp)
    
    END IF

    inv_rhom = DOT_PRODUCT( xs(1:n_solid) , c_inv_rho_s(1:n_solid) )            &
         + xc * inv_rho_c

    IF ( gas_flag .AND. liquid_flag ) inv_rhom = inv_rhom + xl * inv_rho_l

    rho_m = 1.0_wp / inv_rhom

    ! convert from mass fraction to volume fraction
    alphas(1:n_solid) = rho_m * xs(1:n_solid) * c_inv_rho_s(1:n_solid)
    
    h = c_qj(1) * inv_rhom

    ! velocity components
    IF ( REAL( c_qj(1) ) .GT. eps_sing ) THEN

       u = c_qj(2) * inv_cqj1
       v = c_qj(3) * inv_cqj1

    ELSE

       u = SQRT(2.0_wp) * c_qj(1) * c_qj(2) / SQRT( c_qj(1)**4 + eps_sing4 )
       v = SQRT(2.0_wp) * c_qj(1) * c_qj(3) / SQRT( c_qj(1)**4 + eps_sing4 )

    END IF

    RETURN

  END SUBROUTINE c_phys_var


  !******************************************************************************
  !> \brief Physical variables
  !
  !> This subroutine evaluates from the physical real-value local variables qpj, 
  !> all the (real-valued ) variables that define the physical state and that are
  !> needed to compute the explicit equations terms.
  !> \param[in]    qpj          real-valued physical variables 
  !> \param[out]   r_Ri         real-valued Richardson number 
  !> \param[out]   r_rho_m      real-valued mixture density 
  !> \param[out]   r_rho_c      real-valued carrier phase density 
  !> \param[out]   r_red_grav   real-valued reduced gravity
  !
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !> \date 10/10/2019
  !******************************************************************************

  SUBROUTINE mixt_var(qpj,r_Ri,r_rho_m,r_rho_c,r_red_grav)

    IMPLICIT none

    REAL(wp), INTENT(IN) :: qpj(n_vars+2) !< real-value physical variables
    REAL(wp), INTENT(OUT) :: r_Ri         !< real-value Richardson number
    REAL(wp), INTENT(OUT) :: r_rho_m      !< real-value mixture density [kg/m3]
    REAL(wp), INTENT(OUT) :: r_rho_c !< real-value carrier phase density [kg/m3]
    REAL(wp), INTENT(OUT) :: r_red_grav   !< real-value reduced gravity
    REAL(wp) :: r_u                       !< real-value x-velocity
    REAL(wp) :: r_v                       !< real-value y-velocity
    REAL(wp) :: r_h                       !< real-value flow thickness
    REAL(wp) :: r_alphas(n_solid)         !< real-value solid volume fractions
    REAL(wp) :: r_T                       !< real-value temperature [K]
    REAL(wp) :: r_alphal                  !< real-value liquid volume fraction

    REAL(wp) :: alphas_tot                !< total solid fraction

    r_h = qpj(1)

    IF ( qpj(1) .LE. 0.0_wp ) THEN

       r_u = 0.0_wp
       r_v = 0.0_wp
       r_T = T_ambient
       r_alphas(1:n_solid) = 0.0_wp
       r_red_grav = 0.0_wp
       r_Ri = 0.0_wp
       r_rho_m = rho_a_amb
       IF ( gas_flag .AND. liquid_flag ) r_alphal = 0.0_wp

       RETURN

    END IF

    r_u = qpj(n_vars+1)
    r_v = qpj(n_vars+2)
    r_T = qpj(4)

    IF ( alpha_flag ) THEN

       r_alphas(1:n_solid) = qpj(5:4+n_solid)

    ELSE

       r_alphas(1:n_solid) = qpj(5:4+n_solid) / qpj(1)

    END IF

    alphas_tot = SUM(r_alphas)

    IF ( gas_flag ) THEN

       ! continuous phase is air
       r_rho_c =  pres / ( sp_gas_const_a * r_T )

    ELSE

       ! continuous phase is liquid
       r_rho_c = rho_l

    END IF

    IF ( gas_flag .AND. liquid_flag ) THEN

       IF ( alpha_flag ) THEN

          r_alphal = qpj(n_vars)

       ELSE

          r_alphal = qpj(n_vars) / qpj(1)

       END IF

       ! density of mixture of carrier (gas), liquid and solids
       r_rho_m = ( 1.0_wp - alphas_tot - r_alphal ) * r_rho_c                   &
            + DOT_PRODUCT( r_alphas , rho_s ) + r_alphal * rho_l

    ELSE

       ! density of mixture of carrier phase and solids
       r_rho_m = ( 1.0_wp - alphas_tot ) * r_rho_c + DOT_PRODUCT( r_alphas ,    &
            rho_s ) 

    END IF

    ! reduced gravity
    r_red_grav = ( r_rho_m - rho_a_amb ) / r_rho_m * grav

    

    ! Richardson number
    IF ( ( r_u**2 + r_v**2 ) .GT. 0.0_wp ) THEN

       r_Ri = MIN(1.D15,r_red_grav * r_h / ( r_u**2 + r_v**2 ))

    ELSE

       r_Ri = 1.D10

    END IF

    RETURN

  END SUBROUTINE mixt_var

  !******************************************************************************
  !> \brief Conservative to physical variables
  !
  !> This subroutine evaluates from the conservative variables qc the 
  !> array of physical variables qp:\n
  !> - qp(1) = \f$ h \f$
  !> - qp(2) = \f$ hu \f$
  !> - qp(3) = \f$ hv \f$
  !> - qp(4) = \f$ T \f$
  !> - qp(5:4+n_solid) = \f$ alphas(1:n_solid) \f$
  !> - qp(n_vars) = \f$ alphal \f$
  !> - qp(n_vars+1) = \f$ u \f$
  !> - qp(n_vars+2) = \f$ v \f$
  !> .
  !> The physical variables are those used for the linear reconstruction at the
  !> cell interfaces. Limiters are applied to the reconstructed slopes.
  !> \param[in]     qc     local conservative variables 
  !> \param[out]    qp     local physical variables  
  !> \param[out]    p_dyn  local dynamic pressure
  !
  !> \date 2019/11/11
  !
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************

  SUBROUTINE qc_to_qp(qc,qp,p_dyn)

    IMPLICIT none

    REAL(wp), INTENT(IN) :: qc(n_vars)
    REAL(wp), INTENT(OUT) :: qp(n_vars+2)
    REAL(wp), INTENT(OUT) :: p_dyn

    REAL(wp) :: r_h               !< real-value flow thickness
    REAL(wp) :: r_u               !< real-value x-velocity
    REAL(wp) :: r_v               !< real-value y-velocity
    REAL(wp) :: r_alphas(n_solid) !< real-value solid volume fractions
    REAL(wp) :: r_rho_m           !< real-value mixture density [kg/m3]
    REAL(wp) :: r_T               !< real-value temperature [K]
    REAL(wp) :: r_alphal          !< real-value liquid volume fraction
    REAL(wp) :: r_red_grav
    
    CALL r_phys_var(qc,r_h,r_u,r_v,r_alphas,r_rho_m,r_T,r_alphal,r_red_grav)

    qp(1) = r_h
    
    qp(2) = r_h*r_u
    qp(3) = r_h*r_v
    
    qp(4) = r_T

    IF ( alpha_flag ) THEN

       qp(5:4+n_solid) = r_alphas(1:n_solid)
       IF ( gas_flag .AND. liquid_flag ) qp(n_vars) = r_alphal
     
    ELSE

       qp(5:4+n_solid) = r_alphas(1:n_solid) * r_h
       IF ( gas_flag .AND. liquid_flag ) qp(n_vars) = r_alphal * r_h

    END IF
    
    qp(n_vars+1) = r_u
    qp(n_vars+2) = r_v

    p_dyn = 0.5_wp * r_rho_m * ( r_u**2 + r_v**2 )
    
    RETURN

  END SUBROUTINE qc_to_qp

  !******************************************************************************
  !> \brief Physical to conservative variables
  !
  !> This subroutine evaluates the conservative real_value variables qc from the 
  !> array of real_valued physical variables qp:\n
  !> - qp(1) = \f$ h \f$
  !> - qp(2) = \f$ h*u \f$
  !> - qp(3) = \f$ h*v \f$
  !> - qp(4) = \f$ T \f$
  !> - qp(5:4+n_s) = \f$ alphas(1:n_s) \f$
  !> - qp(n_vars) = \f$ alphal \f$
  !> - qp(n_vars+1) = \f$ u \f$
  !> - qp(n_vars+2) = \f$ v \f$
  !> .
  !> \param[in]    qp      physical variables  
  !> \param[in]    B       local topography
  !> \param[out]   qc      conservative variables
  !
  !> \date 2019/11/18
  !
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************

  SUBROUTINE qp_to_qc(qp,qc)

    IMPLICIT none

    REAL(wp), INTENT(IN) :: qp(n_vars+2)
    REAL(wp), INTENT(OUT) :: qc(n_vars)

    REAL(wp) :: r_sp_heat_mix
    REAL(wp) :: sum_sl

    REAL(wp) :: r_u               !< real-value x-velocity
    REAL(wp) :: r_v               !< real-value y-velocity
    REAL(wp) :: r_h               !< real-value flow thickness
    REAL(wp) :: r_hu              !< real-value volumetric x-flow
    REAL(wp) :: r_hv              !< real-value volumetric y-flow
    REAL(wp) :: r_alphas(n_solid) !< real-value solid volume fractions
    REAL(wp) :: r_xl              !< real-value liquid mass fraction
    REAL(wp) :: r_xc              !< real-value carrier phase mass fraction
    REAL(wp) :: r_T               !< real-value temperature [K]
    REAL(wp) :: r_alphal          !< real-value liquid volume fraction
    REAL(wp) :: r_alphac          !< real-value carrier phase volume fraction
    REAL(wp) :: r_rho_m           !< real-value mixture density [kg/m3]
    REAL(wp) :: r_rho_c           !< real-value carrier phase density [kg/m3]
    REAL(wp) :: r_xs(n_solid)     !< real-value solid mass fraction
 
    REAL(wp) :: r_alphas_rhos(n_solid)
    REAL(wp) :: alphas_tot

    REAL(wp) :: r_inv_rhom

    r_h = qp(1)

    IF ( r_h .GT. 0.0_wp ) THEN

       r_hu = qp(2)
       r_hv = qp(3)

       r_u = qp(n_vars+1)
       r_v = qp(n_vars+2)

    ELSE

       r_hu = 0.0_wp
       r_hv = 0.0_wp

       r_u = 0.0_wp
       r_v = 0.0_wp

       qc(1:n_vars) = 0.0_wp
       RETURN

    END IF

    r_T  = qp(4)

    IF ( alpha_flag ) THEN

       r_alphas(1:n_solid) = qp(5:4+n_solid)
   
    ELSE

       r_alphas(1:n_solid) = qp(5:4+n_solid) / qp(1)

    END IF

    alphas_tot = SUM(r_alphas)

    r_alphas_rhos(1:n_solid) = r_alphas(1:n_solid) * rho_s(1:n_solid)
 
    IF ( gas_flag ) THEN

       ! carrier phase is gas
       r_rho_c = pres / ( sp_gas_const_a * r_T )
       sp_heat_c = sp_heat_a

    ELSE

       ! carrier phase is liquid
       r_rho_c = rho_l
       sp_heat_c = sp_heat_l

    END IF

    IF ( gas_flag .AND. liquid_flag ) THEN

       ! mixture of gas, liquid and solid
       IF ( alpha_flag ) THEN
       
          r_alphal = qp(n_vars)

       ELSE

          r_alphal = qp(n_vars) / qp(1)

       END IF

       ! check and correction on dispersed phases volume fractions
       IF ( ( alphas_tot + r_alphal ) .GT. 1.0_wp ) THEN

          sum_sl = alphas_tot + r_alphal
          r_alphas(1:n_solid) = r_alphas(1:n_solid) / sum_sl
          r_alphal = r_alphal / sum_sl

       ELSEIF ( ( alphas_tot + r_alphal ) .LT. 0.0_wp ) THEN

          r_alphas(1:n_solid) = 0.0_wp
          r_alphal = 0.0_wp

       END IF

       ! carrier phase volume fraction
       r_alphac = 1.0_wp - alphas_tot - r_alphal

       ! volume averaged mixture density: carrier (gas) + solids + liquid
       r_rho_m = r_alphac * r_rho_c + SUM( r_alphas_rhos(1:n_solid) )           &
            + r_alphal * rho_l

       r_inv_rhom = 1.0_wp / r_rho_m

       ! liquid mass fraction
       r_xl = r_alphal * rho_l * r_inv_rhom

       ! solid mass fractions
       r_xs(1:n_solid) = r_alphas_rhos(1:n_solid) * r_inv_rhom

       ! carrier (gas) mass fraction
       r_xc = r_alphac * r_rho_c * r_inv_rhom

       ! mass averaged mixture specific heat
       r_sp_heat_mix =  DOT_PRODUCT( r_xs , sp_heat_s ) + r_xl * sp_heat_l      &
            + r_xc * sp_heat_c

    ELSE

       ! mixture of carrier phase ( gas or liquid ) and solid

       ! check and corrections on dispersed phases
       IF ( alphas_tot .GT. 1.0_wp ) THEN

          r_alphas(1:n_solid) = r_alphas(1:n_solid) / alphas_tot

       ELSEIF ( alphas_tot .LT. 0.0_wp ) THEN

          r_alphas(1:n_solid) = 0.0_wp

       END IF

       ! carrier (gas or liquid) volume fraction
       r_alphac = 1.0_wp - alphas_tot 

       ! volume averaged mixture density: carrier (gas or liquid) + solids
       r_rho_m = r_alphac * r_rho_c + SUM( r_alphas_rhos(1:n_solid) )

       r_inv_rhom = 1.0_wp / r_rho_m

       ! solid mass fractions
       r_xs(1:n_solid) = r_alphas_rhos(1:n_solid) * r_inv_rhom

       ! carrier (gas or liquid) mass fraction
       r_xc = r_alphac * r_rho_c * r_inv_rhom

       ! mass averaged mixture specific heat
       r_sp_heat_mix =  DOT_PRODUCT( r_xs , sp_heat_s ) + r_xc * sp_heat_c

    END IF

    qc(1) = r_rho_m * r_h 
    qc(2) = r_rho_m * r_hu
    qc(3) = r_rho_m * r_hv

    IF ( energy_flag ) THEN

       IF ( r_h .GT. 0.0_wp ) THEN

          ! total energy (internal and kinetic)
          qc(4) = r_h * r_rho_m * ( r_sp_heat_mix * r_T                         &
               + 0.5_wp * ( r_u**2 + r_v**2 ) )

       ELSE

          qc(4) = 0.0_wp

       END IF

    ELSE

       ! internal energy
       qc(4) = r_h * r_rho_m * r_sp_heat_mix * r_T 

    END IF

    qc(5:4+n_solid) = r_h * r_alphas(1:n_solid) * rho_s(1:n_solid)
    qc(5:4+n_solid) = r_xs * qc(1)

    IF ( gas_flag .AND. liquid_flag ) qc(n_vars) = r_h * r_alphal * rho_l

    RETURN

  END SUBROUTINE qp_to_qc

  !******************************************************************************
  !> \brief Additional Physical variables
  !
  !> This subroutine evaluates from the physical local variables qpj, the two
  !> additional local variables qp2j = (h+B,u,v). 
  !> \param[in]    qpj    real-valued physical variables 
  !> \param[in]    Bj     real-valued local topography 
  !> \param[out]   qp2j   real-valued physical variables 
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \date 10/10/2019
  !******************************************************************************

  SUBROUTINE qp_to_qp2(qpj,Bj,qp2j)

    IMPLICIT none

    REAL(wp), INTENT(IN) :: qpj(n_vars+2)
    REAL(wp), INTENT(IN) :: Bj
    REAL(wp), INTENT(OUT) :: qp2j(3)

    qp2j(1) = qpj(1) + Bj

    IF ( qpj(1) .LE. 0.0_wp ) THEN

       qp2j(2) = 0.0_wp
       qp2j(3) = 0.0_wp

    ELSE

       qp2j(2) = qpj(2)/qpj(1)
       qp2j(3) = qpj(3)/qpj(1)

    END IF

    RETURN

  END SUBROUTINE qp_to_qp2

  !******************************************************************************
  !> \brief Local Characteristic speeds x direction
  !
  !> This subroutine computes from the physical variable evaluates the largest 
  !> positive and negative characteristic speed in the x-direction. 
  !> \param[in]     qpj           array of local physical variables
  !> \param[out]    vel_min       minimum x-velocity
  !> \param[out]    vel_max       maximum x-velocity
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \date 05/12/2017
  !******************************************************************************

  SUBROUTINE eval_local_speeds_x(qpj,grav_coeff,vel_min,vel_max)

    IMPLICIT none

    REAL(wp), INTENT(IN) :: qpj(n_vars+2)
    REAL(wp), INTENT(IN) :: grav_coeff

    REAL(wp), INTENT(OUT) :: vel_min(n_vars) , vel_max(n_vars)

    REAL(wp) :: r_h          !< real-value flow thickness [m]
    REAL(wp) :: r_u          !< real-value x-velocity [m s-1]
    REAL(wp) :: r_v          !< real-value y-velocity [m s-1]
    REAL(wp) :: r_Ri         !< real-value Richardson number
    REAL(wp) :: r_rho_m      !< real-value mixture density [kg m-3]
    REAL(wp) :: r_rho_c      !< real-value carrier phase density [kg m-3]
    REAL(wp) :: r_red_grav   !< real-value reduced gravity [m s-2]
    REAL(wp) :: r_celerity

    CALL mixt_var(qpj,r_Ri,r_rho_m,r_rho_c,r_red_grav)

    r_h = qpj(1)
    r_u = qpj(n_vars+1)
    r_v = qpj(n_vars+2)

    ! EB : modified with beta_vel

    IF ( r_red_grav * r_h .LT. 0.0_wp ) THEN

       vel_min(1:n_eqns) = beta_vel * r_u
       vel_max(1:n_eqns) = beta_vel * r_u

    ELSE

       r_celerity = SQRT( r_red_grav * r_h * grav_coeff )
       vel_min(1:n_eqns) = beta_vel * r_u - ( beta_vel * (beta_vel - 1.0_wp) * r_u**2 + r_celerity )
       vel_max(1:n_eqns) = beta_vel * r_u + ( beta_vel * (beta_vel - 1.0_wp) * r_u**2 + r_celerity )

    END IF
    
    ! EB : end modified with beta_vel

    RETURN

  END SUBROUTINE eval_local_speeds_x

  !******************************************************************************
  !> \brief Local Characteristic speeds y direction
  !
  !> This subroutine computes from the physical variable evaluates the largest 
  !> positive and negative characteristic speed in the y-direction. 
  !> \param[in]     qpj           array of local physical variables
  !> \param[out]    vel_min       minimum y-velocity
  !> \param[out]    vel_max       maximum y-velocity
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \date 05/12/2017
  !******************************************************************************

  SUBROUTINE eval_local_speeds_y(qpj,grav_coeff,vel_min,vel_max)

    IMPLICIT none

    REAL(wp), INTENT(IN)  :: qpj(n_vars+2)
    REAL(wp), INTENT(IN) :: grav_coeff
    REAL(wp), INTENT(OUT) :: vel_min(n_vars) , vel_max(n_vars)

    REAL(wp) :: r_h          !< real-value flow thickness
    REAL(wp) :: r_u          !< real-value x-velocity
    REAL(wp) :: r_v          !< real-value y-velocity
    REAL(wp) :: r_Ri         !< real-value Richardson number
    REAL(wp) :: r_rho_m      !< real-value mixture density [kg/m3]
    REAL(wp) :: r_rho_c      !< real-value carrier phase density [kg/m3]
    REAL(wp) :: r_red_grav   !< real-value reduced gravity
    REAL(wp) :: r_celerity
    
    CALL mixt_var(qpj,r_Ri,r_rho_m,r_rho_c,r_red_grav)

    r_h = qpj(1)
    r_u = qpj(n_vars+1)
    r_v = qpj(n_vars+2)
    
    ! EB : modified with beta_vel

    IF ( r_red_grav * r_h .LT. 0.0_wp ) THEN

       vel_min(1:n_eqns) = beta_vel * r_v
       vel_max(1:n_eqns) = beta_vel * r_v

    ELSE

       r_celerity = SQRT( r_red_grav * r_h * grav_coeff )
       vel_min(1:n_eqns) = beta_vel * r_v - ( beta_vel * (beta_vel - 1.0_wp) * r_v**2 + r_celerity )
       vel_max(1:n_eqns) = beta_vel * r_v + ( beta_vel * (beta_vel - 1.0_wp) * r_v**2 + r_celerity )

    END IF
    
    ! EB : end modified with beta_vel

    RETURN

  END SUBROUTINE eval_local_speeds_y

  !******************************************************************************
  !> \brief Hyperbolic Fluxes
  !
  !> This subroutine evaluates the numerical fluxes given the conservative
  !> variables qcj and physical variables qpj.
  !> \date 01/06/2012
  !> \param[in]     qcj      real local conservative variables 
  !> \param[in]     qpj      real local physical variables 
  !> \param[in]     dir      direction of the flux (1=x,2=y)
  !> \param[out]    flux     real  fluxes    
  !
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************

  SUBROUTINE eval_fluxes (qcj,qpj,grav_coeff,dir,flux,source_jk) ! EB : added 'source_jk'

    USE parameters_2d, ONLY : eps_sing

    IMPLICIT none

    REAL(wp), INTENT(IN) :: qcj(n_vars)
    REAL(wp), INTENT(IN) :: qpj(n_vars+2)
    REAL(wp), INTENT(IN) :: grav_coeff
    INTEGER, INTENT(IN) :: dir
    REAL(wp), INTENT(IN) :: source_jk ! EB : added

    REAL(wp), INTENT(OUT) :: flux(n_eqns)

    REAL(wp) :: r_h          !< real-value flow thickness [m]
    REAL(wp) :: r_u          !< real-value x-velocity [m s-1]
    REAL(wp) :: r_v          !< real-value y-velocity [m s-1]
    REAL(wp) :: r_Ri         !< real-value Richardson number
    REAL(wp) :: r_rho_m      !< real-value mixture density [kg m-3]
    REAL(wp) :: r_rho_c      !< real-value carrier phase density [kg m-3]
    REAL(wp) :: r_red_grav   !< real-value reduced gravity [m s-2]

    pos_thick:IF ( qpj(1) .GT. eps_sing ) THEN

       r_h = qpj(1)
       r_u = qpj(n_vars+1)
       r_v = qpj(n_vars+2)

       CALL mixt_var(qpj,r_Ri,r_rho_m,r_rho_c,r_red_grav)

       IF ( dir .EQ. 1 ) THEN

          ! Mass flux in x-direction: u * ( rhom * h )
          flux(1) = r_u * qcj(1)
          
          ! EB : modified adding beta_vel
          
          ! x-momentum flux in x-direction + hydrostatic pressure term
          flux(2) = beta_vel * r_u * qcj(2) + 0.5_wp * r_rho_m * grav_coeff * r_red_grav   &
               * r_h**2

          ! y-momentum flux in x-direction: u * ( rho * h * v )
          flux(3) = beta_vel * r_u * qcj(3)
          
          ! EB : end (modified adding beta_vel)

          IF ( energy_flag ) THEN

             ! ENERGY flux in x-direction
             flux(4) = r_u * ( qcj(4) + 0.5_wp * r_rho_m * grav_coeff           &
                  * r_red_grav * r_h**2 )

          ELSE
          
             ! EB : modified. If the cell is on the vent, no temperature profile
             !      is present. Otherwise, beta_T is added far from the vent

             IF ( source_jk .GT. 0.0_wp ) THEN
             
                ! Temperature flux in x-direction: u * ( rhom * Cp * h * T )
                flux(4) = r_u * qcj(4) 
             
             ELSE
             
                ! Temperature flux in x-direction: beta_T * u * ( rhom * Cp * h * T )
                flux(4) = beta_T * r_u * qcj(4) + ( 1 - beta_T ) * T_soil * qcj(2)
             
             END IF

             ! EB : end (modified with beta_T)

          END IF

          ! Mass flux of solid in x-direction: u * ( h * alphas * rhos )
          flux(5:4+n_solid) = r_u * qcj(5:4+n_solid)

          ! Solid flux can't be larger than total flux
          IF ( ( flux(1) .GT. 0.0_wp ) .AND. ( SUM(flux(5:4+n_solid)) / flux(1) &
               .GT. 1.0_wp ) ) THEN

             flux(5:4+n_solid) = &
                  flux(5:4+n_solid) / SUM(flux(5:4+n_solid)) * flux(1)

          END IF

          ! Mass flux of liquid in x-direction: u * ( h * alphal * rhol )
          IF ( gas_flag .AND. liquid_flag ) flux(n_vars) = r_u * qcj(n_vars)

       ELSEIF ( dir .EQ. 2 ) THEN

          ! flux G (derivated wrt y in the equations)
          flux(1) = r_v * qcj(1)
          
          ! EB : modified adding beta_vel

          flux(2) = beta_vel * r_v * qcj(2)

          flux(3) = beta_vel * r_v * qcj(3) + 0.5_wp * r_rho_m * grav_coeff * r_red_grav   &
               * r_h**2
          
          ! EB : end (modified adding beta_vel)

          IF ( energy_flag ) THEN

             ! ENERGY flux in x-direction
             flux(4) = r_v * ( qcj(4) + 0.5_wp * r_rho_m * grav_coeff           &
                  * r_red_grav * r_h**2 )

          ELSE
          
             ! EB : modified. If the cell is on the vent, no temperature profile
             !      is present. Otherwise, beta_T is added far from the vent

             IF ( source_jk .GT. 0.0_wp ) THEN
             
                ! Temperature flux in x-direction: u * ( rhom * Cp * h * T )
                flux(4) = r_v * qcj(4) 
             
             ELSE
             
                ! Temperature flux in x-direction: beta_T * u * ( rhom * Cp * h * T )
                flux(4) = beta_T * r_v * qcj(4) + ( 1 - beta_T) * T_soil * qcj(3)
             
             END IF

             ! EB : end (modified with beta_T)

          END IF

          ! Mass flux of solid in y-direction: v * ( h * alphas * rhos )
          flux(5:4+n_solid) = r_v * qcj(5:4+n_solid)

          ! Solid flux can't be larger than total flux
          IF ( ( flux(1) .GT. 0.0_wp ) .AND. ( SUM(flux(5:4+n_solid)) / flux(1) &
               .GT. 1.0_wp ) ) THEN

             flux(5:4+n_solid) = &
                  flux(5:4+n_solid) / SUM(flux(5:4+n_solid)) * flux(1)

          END IF

          ! Mass flux of liquid in x-direction: u * ( h * alphal * rhol )
          IF ( gas_flag .AND. liquid_flag ) flux(n_vars) = r_v * qcj(n_vars)

       END IF

    ELSE

       flux(1:n_eqns) = 0.0_wp

    ENDIF pos_thick

    RETURN

  END SUBROUTINE eval_fluxes

  !******************************************************************************
  !> \brief Explicit source term
  !
  !> This subroutine evaluates the non-hyperbolic terms to be treated explicitely
  !> in the DIRK numerical scheme (e.g. gravity,source of mass). The sign of the
  !> terms is taken with the terms on the right-hand side of the equations.
  !> \date 2019/12/13
  !> \param[in]     B_primej_x         local x-slope
  !> \param[in]     B_primej_y         local y_slope
  !> \param[in]     B_secondj_xx       local 2nd derivative in x-direction
  !> \param[in]     B_secondj_xy       local 2nd derivative in xy-direction
  !> \param[in]     B_secondj_yy       local 2nd derivative in y-direction
  !> \param[in]     grav_coeff         correction factor for topography slope
  !> \param[in]     d_grav_coeff_dx    x-derivative of grav_coeff
  !> \param[in]     d_grav_coeff_dy    y-derivative of grav_coeff
  !> \param[in]     source_xy          local source of mass
  !> \param[in]     qpj                physical variables
  !> \param[in]     time               simlation time (needed for source)
  !> \param[in]     cell_fract_jk      fraction of cell contributing to source
  !> \param[out]    expl_term          explicit term
  !
  !> @author 
  !> Mattia de' Michieli Vitturi
  ! EB : credo che 'source_xy' non sia più inutilizzato
  !******************************************************************************

  SUBROUTINE eval_expl_terms( Bprimej_x, Bprimej_y, Bsecondj_xx , Bsecondj_xy , &
       Bsecondj_yy, grav_coeff, d_grav_coeff_dx , d_grav_coeff_dy , source_xy , &
       qpj, expl_term, time, cell_fract_jk )

    USE parameters_2d, ONLY : vel_source , T_source ,            &
         time_param , bottom_radial_source_flag

    
    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: Bprimej_x
    REAL(wp), INTENT(IN) :: Bprimej_y
    REAL(wp), INTENT(IN) :: Bsecondj_xx
    REAL(wp), INTENT(IN) :: Bsecondj_xy
    REAL(wp), INTENT(IN) :: Bsecondj_yy
    REAL(wp), INTENT(IN) :: grav_coeff
    REAL(wp), INTENT(IN) :: d_grav_coeff_dx
    REAL(wp), INTENT(IN) :: d_grav_coeff_dy
    
    REAL(wp), INTENT(IN) :: source_xy

    REAL(wp), INTENT(IN) :: qpj(n_vars+2)      !< local physical variables 
    REAL(wp), INTENT(OUT) :: expl_term(n_eqns) !< local explicit forces 

    REAL(wp), INTENT(IN) :: time
    REAL(wp), INTENT(IN) :: cell_fract_jk
    
    REAL(wp) :: r_h          !< real-value flow thickness
    REAL(wp) :: r_u          !< real-value x-velocity
    REAL(wp) :: r_v          !< real-value y-velocity
    REAL(wp) :: r_Ri         !< real-value Richardson number
    REAL(wp) :: r_rho_m      !< real-value mixture density [kg/m3]
    REAL(wp) :: r_rho_c      !< real-value carrier phase density [kg/m3]
    REAL(wp) :: r_red_grav   !< real-value reduced gravity
    REAL(wp) :: r_alphas(n_solid)  !< real-value solid volume fractions
    REAL(wp) :: r_xs(n_solid)      !< real-value solid mass fractions 
    REAL(wp) :: r_alphal     !< real-value liquid volume fraction      
    REAL(wp) :: r_xl         !< real-value liquid mass fraction
    REAL(wp) :: r_alphac     !< real-value carrier phase volume fraction
    REAL(wp) :: r_xc         !< real-values carrier phase mass fraction
    
    REAL(wp) :: alphas_tot   !< total volume fraction of solid
    REAL(wp) :: sum_sl       !< sum of liquid and solid volume fractions
    REAL(wp) :: r_sp_heat_mix !< real_value mixture specific heat
    
    REAL(wp) :: t_rem
    REAL(wp) :: t_coeff
    REAL(wp) :: h_dot

    REAL(wp) :: r_tilde_grav
    REAL(wp) :: centr_force_term
    
    expl_term(1:n_eqns) = 0.0_wp

    IF ( ( qpj(1) .LE. 0.0_wp ) .AND. ( cell_fract_jk .EQ. 0.0_wp ) ) RETURN

    r_h = qpj(1)
    r_u = qpj(n_vars+1)
    r_v = qpj(n_vars+2)

    CALL mixt_var(qpj,r_Ri,r_rho_m,r_rho_c,r_red_grav)

    IF ( curvature_term_flag ) THEN

       centr_force_term = Bsecondj_xx * r_u**2 + Bsecondj_xy * r_u * r_v +      &
            Bsecondj_yy * r_v**2

    ELSE

       centr_force_term = 0.0_wp

    END IF
    
    r_tilde_grav = r_red_grav + centr_force_term

    ! units of dqc(2)/dt [kg m-1 s-2]
    expl_term(2) = - grav_coeff * r_rho_m * r_tilde_grav * r_h * Bprimej_x      &
         + grav_coeff * r_red_grav * r_rho_m * 0.5_wp * r_h**2 * d_grav_coeff_dx 
    
    ! units of dqc(3)/dt [kg m-1 s-2]
    expl_term(3) = - grav_coeff * r_rho_m * r_tilde_grav * r_h * Bprimej_y      &
         + grav_coeff * r_red_grav * r_rho_m * 0.5_wp * r_h**2 * d_grav_coeff_dx

    IF ( energy_flag ) THEN

       expl_term(4) = expl_term(2) * r_u + expl_term(3) * r_v

    ELSE

       expl_term(4) = 0.0_wp

    END IF
    
    ! ----------- ADDITIONAL EXPLICIT TERMS FOR BOTTOM RADIAL SOURCE ------------ 

    IF ( .NOT.bottom_radial_source_flag ) THEN

       RETURN

    END IF

    t_rem = MOD( time + time_param(4) , time_param(1) )

    IF ( time_param(3) .EQ. 0.0_wp ) THEN

       IF ( t_rem .LE. time_param(2) ) THEN

          t_coeff = 1.0_wp

       ELSE

          t_coeff = 0.0_wp

       END IF
          
    ELSE

       IF ( t_rem .LE. time_param(3) ) THEN

          t_coeff = ( t_rem / time_param(3) ) 

       ELSEIF ( t_rem .LE. time_param(2) - time_param(3) ) THEN

          t_coeff = 1.0_wp
          
       ELSEIF ( t_rem .LE. time_param(2) ) THEN
          
          t_coeff = 1.0_wp - ( t_rem - time_param(2) + time_param(3) ) /        &
               time_param(3)
          
       ELSE
          
          t_coeff = 0.0_wp
          
       END IF

    END IF

    h_dot = cell_fract_jk * vel_source
        
    IF ( gas_flag ) THEN

       ! carrier phase is gas
       r_rho_c = pres / ( sp_gas_const_a * t_source )
       sp_heat_c = sp_heat_a

    ELSE

       ! carrier phase is liquid
       r_rho_c = rho_l
       sp_heat_c = sp_heat_l

    END IF

       ! mixture of carrier phase ( gas or liquid ) and solid

       ! check and corrections on dispersed phases
       IF ( alphas_tot .GT. 1.0_wp ) THEN

          r_alphas(1:n_solid) = r_alphas(1:n_solid) / alphas_tot

       ELSEIF ( alphas_tot .LT. 0.0_wp ) THEN

          r_alphas(1:n_solid) = 0.0_wp

       END IF

    ! carrier (gas or liquid) volume fraction
    r_alphac = 1.0_wp - alphas_tot 

    ! volume averaged mixture density: carrier (gas or liquid) + solids
    r_rho_m = r_alphac * r_rho_c + DOT_PRODUCT( r_alphas , rho_s ) 

    ! solid mass fractions
    r_xs(1:n_solid) = r_alphas(1:n_solid) * rho_s(1:n_solid) / r_rho_m

    ! carrier (gas or liquid) mass fraction
    r_xc = r_alphac * r_rho_c / r_rho_m

    ! mass averaged mixture specific heat
    r_sp_heat_mix =  DOT_PRODUCT( r_xs , sp_heat_s ) + r_xc * sp_heat_c

    
    expl_term(1) = expl_term(1) + t_coeff * h_dot * r_rho_m
    expl_term(2) = expl_term(2) + 0.0_wp
    expl_term(3) = expl_term(3) + 0.0_wp

    IF ( energy_flag ) THEN

       expl_term(4) = expl_term(4) + t_coeff * h_dot * r_rho_m * r_sp_heat_mix  &
            * t_source

    ELSE ! EB : non ho capito questo 'ELSE'

       expl_term(4) = expl_term(4) + t_coeff * h_dot * r_rho_m * r_sp_heat_mix  &
            * t_source

    END IF
  
    RETURN

  END SUBROUTINE eval_expl_terms
  
  !******************************************************************************
  !> \brief Implicit source terms
  !
  !> This subroutine evaluates the source terms  of the system of equations,
  !> both for real or complex inputs, that are treated implicitely in the DIRK
  !> numerical scheme.
  !> \date 01/06/2012
  !> \param[in]     Bprimej_x       topography slope in x-direction
  !> \param[in]     Bprimej_y       topography slope in y-direction
  !> \param[in]     c_qj            complex conservative variables 
  !> \param[in]     r_qj            real conservative variables 
  !> \param[out]    c_nh_term_impl  complex non-hyperbolic terms     
  !> \param[out]    r_nh_term_impl  real non-hyperbolic terms
  !
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************

  SUBROUTINE eval_implicit_terms( Bprimej_x, Bprimej_y, source_jk, c_qj, c_nh_term_impl,   &
       r_qj , r_nh_term_impl ) ! EB : added "source_jk"

    USE COMPLEXIFY 

    USE parameters_2d, ONLY : four_thirds , neg_four_thirds

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: Bprimej_x
    REAL(wp), INTENT(IN) :: Bprimej_y
    COMPLEX(wp), INTENT(IN), OPTIONAL :: c_qj(n_vars)
    COMPLEX(wp), INTENT(OUT), OPTIONAL :: c_nh_term_impl(n_eqns)
    REAL(wp), INTENT(IN), OPTIONAL :: r_qj(n_vars)
    REAL(wp), INTENT(OUT), OPTIONAL :: r_nh_term_impl(n_eqns)
    REAL(wp), INTENT(IN) :: source_jk ! EB : added

    COMPLEX(wp) :: h                       !< height [m]
    COMPLEX(wp) :: inv_h                   !< 1/height [m-1]
    COMPLEX(wp) :: u                       !< velocity (x direction) [m/s]
    COMPLEX(wp) :: v                       !< velocity (y direction) [m/s]
    COMPLEX(wp) :: w                       !< velocity (z direction) [m/s]
    COMPLEX(wp) :: T                       !< temperature [K]
    COMPLEX(wp) :: rho_m                   !< mixture density [kg/m3]
    COMPLEX(wp) :: alphas(n_solid)         !< sediment volume fractions
    COMPLEX(wp) :: inv_rho_m               !< 1/mixture density [kg-1 m3]
 
    COMPLEX(wp) :: qj(n_vars)
    COMPLEX(wp) :: nh_term(n_eqns)
    COMPLEX(wp) :: source_term(n_eqns)

    COMPLEX(wp) :: mod_vel
    COMPLEX(wp) :: mod_vel2
    COMPLEX(wp) :: gamma
    REAL(wp) :: h_threshold

    INTEGER :: i

    !--- Lahars rheology model variables

    !> Temperature in C
    COMPLEX(wp) :: Tc

    COMPLEX(wp) :: expA , expB

    !> 1st param for fluid viscosity empirical relationship (O'Brian et al, 1993)
    COMPLEX(wp) :: alpha1    ! (units: kg m-1 s-1 )
    
    !> Fluid dynamic viscosity (units: kg m-1 s-1 )
    COMPLEX(wp) :: fluid_visc

    !> Total friction slope (dimensionless): s_f = s_v+s_td+s_y
    COMPLEX(wp) :: s_f

    !> Viscous slope component of total Friction (dimensionless)
    COMPLEX(wp) :: s_v

    !> Turbulent dispersive slope component of total friction (dimensionless)
    COMPLEX(wp) :: s_td

    COMPLEX(wp) :: temp_term

    IF ( present(c_qj) .AND. present(c_nh_term_impl) ) THEN

       qj = c_qj

    ELSEIF ( present(r_qj) .AND. present(r_nh_term_impl) ) THEN

       DO i = 1,n_vars

          qj(i) = CMPLX( r_qj(i),0.0_wp,wp )

       END DO

    ELSE

       WRITE(*,*) 'Constitutive, eval_fluxes: problem with arguments'
       STOP

    END IF

    ! initialize the source terms
    source_term(1:n_eqns) = CMPLX(0.0_wp,0.0_wp,wp)

    IF (rheology_flag) THEN

       CALL c_phys_var(qj,h,u,v,T,rho_m,alphas,inv_rho_m)

       IF ( slope_correction_flag ) THEN

          w = u * Bprimej_x + v * Bprimej_y

       ELSE

          w = CMPLX( 0.0_wp , 0.0_wp , wp )

       END IF

       mod_vel2 = u**2 + v**2 + w**2
       mod_vel = SQRT( mod_vel2 )

       IF ( rheology_model .EQ. 1 ) THEN
          ! Voellmy Salm rheology

          IF ( REAL(mod_vel) .NE. 0.0_wp ) THEN 

             ! IMPORTANT: grav3_surf is always negative 
             source_term(2) = source_term(2) - rho_m * ( u / mod_vel ) *        &
                  ( grav / xi ) * mod_vel2

             source_term(3) = source_term(3) - rho_m * ( v / mod_vel ) *        &
                  ( grav / xi ) * mod_vel2

          ENDIF

       ELSEIF ( rheology_model .EQ. 2 ) THEN

          ! Plastic rheology
          IF ( REAL(mod_vel) .NE. 0.0_wp ) THEN 

             source_term(2) = source_term(2) - rho_m * tau * ( u / mod_vel )

             source_term(3) = source_term(3) - rho_m * tau * ( v / mod_vel )

          ENDIF

       ELSEIF ( rheology_model .EQ. 3 ) THEN

          h_threshold = 1.0E-10_wp

          ! Temperature dependent rheology
          IF ( REAL(h) .GT. h_threshold ) THEN

             ! Equation 6 from Costa & Macedonio, 2005
             gamma = 3.0_wp * nu_ref / h * EXP( - visc_par * ( T - T_ref ) )

          ELSE

             ! Equation 6 from Costa & Macedonio, 2005
             gamma = 3.0_wp * nu_ref / h_threshold * EXP( - visc_par            &
                  * ( T - T_ref ) )

          END IF

          IF ( REAL(mod_vel) .NE. 0.0_wp ) THEN 

             ! Last R.H.S. term in equation 2 from Costa & Macedonio, 2005
             source_term(2) = source_term(2) - rho_m * gamma * u

             ! Last R.H.S. term in equation 3 from Costa & Macedonio, 2005
             source_term(3) = source_term(3) - rho_m * gamma * v

          ENDIF
          
       ELSEIF ( rheology_model .EQ. 5 ) THEN

          tau = 1.0E-3_wp / ( 1.0_wp + 10.0_wp * h ) * mod_vel

          IF ( REAL(mod_vel) .NE. 0.0_wp ) THEN

             source_term(2) = source_term(2) - rho_m * tau * ( u / mod_vel )
             source_term(3) = source_term(3) - rho_m * tau * ( v / mod_vel )

          END IF


       ELSEIF ( rheology_model .EQ. 6 ) THEN

          IF ( REAL(mod_vel) .NE. 0.0_wp ) THEN 

             source_term(2) = source_term(2) - rho_m * u * friction_factor *    &
                  mod_vel

             source_term(3) = source_term(3) - rho_m * v * friction_factor *    &
                  mod_vel

          ENDIF

       ENDIF
       
       ! EB : added the source terms for the thermal heat exchanges
       
       h_threshold = 1.0E-10_wp
       
       IF ( h .GT. h_threshold ) THEN

          conductive_term_coeff = thermal_conductivity_fluid * enne / h
          
       ELSE
       
          conductive_term_coeff = 0.0_wp

       END IF          
       
       IF ( temperature_profile_flag ) THEN
       
          ! At the vent: no temperature profile ( T = T_ground = T_surf )
          ! and no conduction
          IF ( source_jk .GT. 0.0_wp ) THEN 
          
             source_term(4) = source_term(4)                                   &
               - radiative_term_coeff * ( T**4 - T_env**4 )                    &
               - convective_term_coeff * ( T - T_env )          
          
          ELSE
          
             source_term(4) = source_term(4)                                   &
               - radiative_term_coeff *                                        &
                 ( ( zeta * T + ( 1 - zeta ) * T_soil)**4 - T_env**4 )         &
               - convective_term_coeff *                                       &
                 ( zeta * T + ( 1 - zeta ) * T_soil - T_env )                  &
               - conductive_term_coeff * ( zeta - eta ) * ( T - T_soil )
          
          END IF
               
       ELSE
       
          ! At the vent: no temperature profile ( T = T_ground = T_surf )
          ! and no conduction
          IF ( source_jk .GT. 0.0_wp ) THEN 
          
             source_term(4) = source_term(4)                                   &
               - radiative_term_coeff * ( T**4 - T_env**4 )                    &
               - convective_term_coeff * ( T - T_env )                         
          
          ELSE
          
             source_term(4) = source_term(4)                                   &
               - radiative_term_coeff * ( T**4 - T_env**4 )                    &
               - convective_term_coeff * ( T - T_env )                         &
               - conductive_term_coeff * ( T - T_ground )
          
          END IF
    
       END IF
               
       ! EB : end (added the source terms for the thermal heat exchanges)

    ENDIF

    nh_term = source_term

    IF ( present(c_qj) .AND. present(c_nh_term_impl) ) THEN

       c_nh_term_impl = nh_term

    ELSEIF ( present(r_qj) .AND. present(r_nh_term_impl) ) THEN

       r_nh_term_impl = REAL( nh_term )

    END IF

    RETURN

  END SUBROUTINE eval_implicit_terms

  !******************************************************************************
  !> \brief Non-Hyperbolic semi-implicit terms
  !
  !> This subroutine evaluates the non-hyperbolic terms that are solved
  !> semi-implicitely by the solver. For example, any discontinuous term that
  !> appears in the friction terms.
  !> \date 20/01/2018
  !> \param[in]     grav3_surf         gravity correction 
  !> \param[in]     qcj                real conservative variables 
  !> \param[out]    nh_semi_impl_term  real non-hyperbolic terms
  !
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************

  SUBROUTINE eval_nh_semi_impl_terms( Bprimej_x , Bprimej_y , Bsecondj_xx ,     &
       Bsecondj_xy , Bsecondj_yy , grav_coeff , qcj , nh_semi_impl_term )

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: Bprimej_x
    REAL(wp), INTENT(IN) :: Bprimej_y
    REAL(wp), INTENT(IN) :: Bsecondj_xx
    REAL(wp), INTENT(IN) :: Bsecondj_xy
    REAL(wp), INTENT(IN) :: Bsecondj_yy
    REAL(wp), INTENT(IN) :: grav_coeff

    REAL(wp), INTENT(IN) :: qcj(n_vars)
    REAL(wp), INTENT(OUT) :: nh_semi_impl_term(n_eqns)

    REAL(wp) :: source_term(n_eqns)

    REAL(wp) :: mod_vel

    REAL(wp) :: h_threshold

    !--- Lahars rheology model variables

    !> Yield strenght (units: kg m-1 s-2)
    REAL(wp) :: tau_y

    !> Yield slope component of total friction (dimensionless)
    REAL(wp) :: s_y

    REAL(wp) :: r_h               !< real-value flow thickness
    REAL(wp) :: r_u               !< real-value x-velocity
    REAL(wp) :: r_v               !< real-value y-velocity
    REAL(wp) :: r_w               !< real_value z-velocity
    REAL(wp) :: r_alphas(n_solid) !< real-value solid volume fractions
    REAL(wp) :: r_rho_m           !< real-value mixture density [kg/m3]
    REAL(wp) :: r_T               !< real-value temperature [K]
    REAL(wp) :: r_alphal          !< real-value liquid volume fraction
    REAL(wp) :: r_red_grav
    
    REAL(wp) :: temp_term
    REAL(wp) :: centr_force_term

    ! initialize and evaluate the forces terms
    source_term(1:n_eqns) = 0.0_wp

    IF (rheology_flag) THEN

       CALL r_phys_var(qcj,r_h,r_u,r_v,r_alphas,r_rho_m,r_T,r_alphal,r_red_grav)
    
       ! Voellmy Salm rheology
       IF ( rheology_model .EQ. 1 ) THEN
          
          IF ( slope_correction_flag ) THEN

             r_w = r_u * Bprimej_x + r_v * Bprimej_y

          ELSE

             r_w = 0.0_wp

          END IF

          mod_vel = SQRT( r_u**2 + r_v**2 + r_w**2 )
          
          IF ( mod_vel .GT. 0.0_wp ) THEN

             IF ( curvature_term_flag ) THEN

                ! centrifugal force term: (u,v)^T*Hessian*(u,v)
                centr_force_term = Bsecondj_xx * r_u**2 + Bsecondj_xy * r_u *   &
                     r_v + Bsecondj_yy * r_v**2

             ELSE

                centr_force_term = 0.0_wp 

             END IF

             temp_term = r_rho_m *  mu * r_h * grav_coeff * ( r_red_grav +      &
                  centr_force_term ) / mod_vel
             
             ! units of dqc(2)/dt=d(rho h v)/dt (kg m-1 s-2)
             source_term(2) = source_term(2) - temp_term * r_u

             ! units of dqc(3)/dt=d(rho h v)/dt (kg m-1 s-2)
             source_term(3) = source_term(3) - temp_term * r_v

          END IF

          ! Plastic rheology
       ELSEIF ( rheology_model .EQ. 2 ) THEN


          ! Temperature dependent rheology
       ELSEIF ( rheology_model .EQ. 3 ) THEN

          mod_vel = SQRT( r_u**2 + r_v**2 )
          
          IF ( mod_vel .GT. 0.0_wp ) THEN

             ! units of dqc(2)/dt [kg m-1 s-2]
             source_term(2) = source_term(2) - tau0 * r_u / mod_vel

             ! units of dqc(3)/dt [kg m-1 s-2]
             source_term(3) = source_term(3) - tau0 * r_v / mod_vel

          END IF
   
       ELSEIF ( rheology_model .EQ. 5 ) THEN

       ENDIF

    ENDIF

    nh_semi_impl_term = source_term

    RETURN

  END SUBROUTINE eval_nh_semi_impl_terms

END MODULE constitutive_2d


