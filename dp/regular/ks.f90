module KUST_STI
! Description:
!    Contains the subroutines necessary for the Kustaanheimo-Stiefel
!    formulation.
! 
! References:
! [1] Stiefel, E. L. and Scheifele, G. "Linear and Regular Celestial Mechanics",
!     Springer-Verlag, 1971.
! 
! Author:
!    Davide Amato
!    Space Dynamics Group - Technical University of Madrid
!    d.amato@upm.es
!
! ==============================================================================

use KINDS, only: dk
implicit none


contains


! ==============================================================================
! 01. INTEGRATION PROCEDURES
! ==============================================================================

subroutine KS_RHS(neq,s,u,udot)
! Description:
!    Computes the value of the right-hand side of the equations of motion of the
!    Kustaanheimo-Stiefel formulation. The equation for the perturbed harmonic
!    oscillator integrated here is written as:
!    
!    u'' + (h/2) * u = -(u/2) * V + (u^2/2) L^T * F
!    
!    where F = P - dV/dx is the total perturbation. This is equivalent to
!    Eq. (9, 53) of Ref. [1].
!
! ==============================================================================

! MODULES
use PHYS_CONST,    only: GE_nd,RE_nd,ERR_constant_nd
use SETTINGS,      only: eqs,insgrav,isun,imoon,idrag,iSRP
use PERTURBATIONS, only: PPOTENTIAL,PACC_ICRF,PACC_MMEIAUE
use AUXILIARIES,   only: coordSyst

! VARIABLES
implicit none
! Arguments
integer,intent(in)       ::  neq                  ! Number of equations.
real(dk),intent(in)      ::  s                    ! Value of fictitious time.
real(dk),intent(in)      ::  u(1:neq)             ! KS state vector.
real(dk),intent(out)     ::  udot(1:neq)          ! RHS of EoM's, ND.
! Local variables
! -- State variables
real(dk)    ::  x(1:4),xdot(1:4)                  ! Radius and velocity in R^4, ND
real(dk)    ::  rmag                              ! Radius magnitude, ND
real(dk)    ::  t                                 ! Physical time, ND
! -- Perturbations
real(dk)    ::  F(1:4)               ! Total perturbation
real(dk)    ::  P(1:4)               ! Non-potential part
real(dk)    ::  Vpot                 ! Perturbing potential (dimensionless)
real(dk)    ::  mdVdr(1:4),dVdt      ! Partial derivatives of the potential
real(dk)    ::  gradV_sph(1:3)       ! Gradient of V in spherical coordinates
!real(dk)    ::  rVpot               ! Perturbing potential term
!real(dk)    ::  drVdu(1:4)          ! Derivatives of the perturbing potential term
! -- Misc
real(dk)    ::  L(1:4,1:4)           ! KS matrix

! ==============================================================================

! STATE VECTOR DICTIONARY
! u(1:4)        u1,...,u4; KS-position, in R^4
! u(5:8)        u1',...,u4'; KS-velocity, in R^4
! u(9)          h; (-total energy) = (-Keplerian energy) + (-potential)
! u(10)         t; non-dimensional physical time

! ==============================================================================
! 01. COMPUTE CARTESIAN COORDINATES
! ==============================================================================

call KS2CART(u,x,xdot)
rmag = sqrt(dot_product(x,x))
t = u(10)

! ==============================================================================
! 02. COMPUTE POTENTIAL PERTURBATIONS
! ==============================================================================

Vpot = 0._dk; mdVdr = 0._dk; dVdt = 0._dk
Vpot = PPOTENTIAL(insgrav,GE_nd,RE_nd,x(1:3),rmag,t)
mdVdr(1:3) = PACC_ICRF(insgrav,0,0,0,0,x(1:3),xdot(1:3),rmag,t,gradV_sph)
dVdt = gradV_sph(3) * ERR_constant_nd

! ==============================================================================
! 03. COMPUTE NON-POTENTIAL PERTURBATIONS
! ==============================================================================

P = 0._dk
select case(trim(coordSyst))
  case('ICRF')
    P(1:3) = PACC_ICRF(0,isun,imoon,idrag,iSRP,x(1:3),xdot(1:3),rmag,t)
  
  case('MMEIAUE')
    P(1:3) = PACC_MMEIAUE(0,isun,imoon,iSRP,x(1:3),xdot(1:3),rmag,t)

end select
F = mdVdr + P

! ==============================================================================
! 04. EVALUATE RIGHT-HAND SIDE
! ==============================================================================

udot = 0._dk; L = 0._dk

! KS matrix
L = KSMAT(u(1:4))

! Velocities
udot(1:4) = u(5:8)

! Accelerations
udot(5:8) = -0.5_dk * ( u(1:4) * (u(9) + Vpot) - rmag * matmul(transpose(L),F) )

! Total energy
udot(9) = -rmag * dVdt - 2._dk * dot_product(u(5:8),matmul(transpose(L),P))

! Time
udot(10) = rmag

end subroutine KS_RHS


function KSMAT(u)
! Description:
!    Computes the KS-matrix from the K-S state vector.
! 
! ==============================================================================

! VARIABLES
implicit none
! Arguments
real(dk),intent(in)  ::  u(:)                ! K-S state vector
real(dk)             ::  KSMAT(1:4,1:4)      ! K-S-matrix

! ==============================================================================

KSMAT = RESHAPE([ u(1),  u(2),  u(3),  u(4)&
                &,-u(2),  u(1),  u(4), -u(3)&
                &,-u(3), -u(4),  u(1),  u(2)&
                &, u(4), -u(3),  u(2), -u(1)],[4,4])

end function KSMAT 


subroutine KS2CART(u,x,xdot)
! Description:
!    Computes the Cartesian state (dimensionless) from the KS state vector.
! 
! ==============================================================================

! VARIABLES
implicit none
! Arguments
real(dk),intent(in)   ::  u(:)                  ! KS (extended) state vector
real(dk),intent(out)  ::  x(:),xdot(:)          ! Position and velocity, ND
! Locals
real(dk)              ::  r                     ! Radius magnitude, ND

! ==============================================================================

! Position
x(1) = u(1)**2 - u(2)**2 - u(3)**2 + u(4)**2
x(2) = 2._dk*(u(1)*u(2) - u(3)*u(4))
x(3) = 2._dk*(u(1)*u(3) + u(2)*u(4))
r    = u(1)**2 + u(2)**2 + u(3)**2 + u(4)**2

! Velocity
xdot(1) = 2._dk*(u(1)*u(5) - u(2)*u(6) - u(3)*u(7) + u(4)*u(8))/r
xdot(2) = 2._dk*(u(2)*u(5) + u(1)*u(6) - u(4)*u(7) - u(3)*u(8))/r
xdot(3) = 2._dk*(u(3)*u(5) + u(4)*u(6) + u(1)*u(7) + u(2)*u(8))/r

! If x, xdot are in R^4 their fourth component is null.
if (size(x,1) > 3)    x(4:) = 0._dk
if (size(xdot,1) > 3) xdot(4:) = 0._dk

end subroutine KS2CART


subroutine CART2KS(R0_d,V0_d,t0_d,mu,DU,TU,u0,Vpot_d)
! Description:
!     Transforms Cartesian position and velocity into the KS state vector. 
!
! ==============================================================================

! VARIABLES
implicit none
! Arguments
real(dk),intent(in)     ::  R0_d(1:3),V0_d(1:3)  ! Dimensional position and velocity [km,km/s]
real(dk),intent(in)     ::  t0_d                 ! Dimensional time [s]
real(dk),intent(in)     ::  mu                   ! Grav parameter of the main body [km^3/s^2]
real(dk),intent(in)     ::  DU,TU                ! Ref. quantities for ND [km,1/s]
real(dk),intent(in)     ::  Vpot_d               ! Dimensional potential [km^2/s^2]
real(dk),intent(out)    ::  u0(:)                ! Extended KS state vector

! Locals
! -- State vectors and associated quantities
real(dk)            ::  x0(1:4)              ! Radius vector in R^4, ND
real(dk)            ::  xdot0(1:4)           ! Velocity vector in R^4, ND
real(dk)            ::  t0                   ! Initial time, ND
real(dk)            ::  r                    ! Radius magnitude, ND
! -- Energies and potentials
real(dk)            ::  h                    ! Total energy, ND
real(dk)            ::  Ksq                  ! Grav parameter of the main body, ND
real(dk)            ::  Kin                  ! Kinetic energy, ND
real(dk)            ::  Vpot                 ! Perturbing potential, ND

! ==============================================================================

! STATE VECTOR DICTIONARY
! u(1:4)        u1,...,u4; KS-position, in R^4
! u(5:8)        u1',...,u4'; KS-velocity, in R^4
! u(9)          h; (-total energy) = (-Keplerian energy) + (-potential)
! u(10)         t; non-dimensional physical time

! ==============================================================================
! 01. INITIALIZE POSITION AND VELOCITY IN R^4, NON-DIMENSIONALIZATIONS
! ==============================================================================

x0(1:3) = R0_d/DU
x0(4)   = 0._dk
r       = sqrt(dot_product(x0,x0))
t0      = t0_d*TU

xdot0(1:3) = V0_d/(DU*TU)
xdot0(4)   = 0._dk

! Also, non-dimensionalize potential and grav parameter
Vpot = Vpot_d/(DU*TU)**2
Ksq  = mu/(DU**3*TU**2)

! ==============================================================================
! 02. INITIALIZE KS-POSITION AND KS-VELOCITY (state vector u)
! ==============================================================================

! KS-POSITION
if ( x0(1) >= 0.) then
    u0(1) = 0._dk
    u0(4) = sqrt(.5_dk*(r + x0(1)) - u0(1)**2)
    u0(2) = (x0(2)*u0(1) + x0(3)*u0(4))/(r + x0(1))
    u0(3) = (x0(3)*u0(1) - x0(2)*u0(4))/(r + x0(1))
else
    u0(2) = 0._dk
    u0(3) = sqrt(.5_dk*(r - x0(1)) - u0(2)**2)
    u0(1) = (x0(2)*u0(2) + x0(3)*u0(3))/(r - x0(1))
    u0(4) = (x0(3)*u0(2) - x0(2)*u0(3))/(r - x0(1))
end if

! KS-VELOCITY
u0(5) = .5_dk*dot_product(u0(1:3),xdot0(1:3))
u0(6) = .5_dk*dot_product([-u0(2),u0(1),u0(4)],xdot0(1:3))
u0(7) = .5_dk*dot_product([-u0(3),-u0(4),u0(1)],xdot0(1:3))
u0(8) = .5_dk*dot_product([u0(4),-u0(3),u0(2)],xdot0(1:3))

! Total energy
Kin   = .5_dk*dot_product(xdot0,xdot0)
h     = Ksq/r - Kin - Vpot
u0(9) = h

! Initial time
u0(10) = t0

end subroutine CART2KS


subroutine KS_EVT(neq,phi,u,ng,roots)

! MODULES
use AUXILIARIES, only: MJD0,MJDnext,MJDf,DU,TU
use PHYS_CONST,  only: secsPerDay,RE,reentry_radius_nd
use SETTINGS,    only: eqs

! VARIABLES
implicit none
! Arguments IN
integer,intent(in)    ::  neq
integer,intent(in)    ::  ng
real(dk),intent(in)   ::  phi
real(dk),intent(in)   ::  u(1:neq)
! Arguments OUT
real(dk),intent(out)  ::  roots(1:ng)

! Locals
real(dk)  ::  t         ! Current time [-]
real(dk)  ::  rmag

! ==============================================================================

roots = 1._dk

! Get time
!flag_time = eqs - 2
t = u(10)

! ==============================================================================
! 01. Next timestep
! ==============================================================================

roots(1) = t - (MJDnext - MJD0)*secsPerDay*TU

! ==============================================================================
! 02. Stop integration
! ==============================================================================

roots(2) = t - (MJDf - MJD0)*secsPerDay*TU

! ==============================================================================
! 03. Re-entry
! ==============================================================================

rmag  = dot_product(u(1:4),u(1:4))

roots(3) = rmag - reentry_radius_nd

end subroutine KS_EVT


function KS_TE2TIME(u,s,flag_time)
! Description:
!    Gets the value of physical time from the KS state vector "u" and
!    fictitious time "s".
!
! ==============================================================================

! VARIABLES
implicit none
! Arguments IN
real(dk),intent(in)     ::  u(1:10),s
integer,intent(in)      ::  flag_time
! Function definition
real(dk)                ::  KS_TE2TIME

! ==============================================================================

if ( flag_time == 0 ) then
    ! Physical time
    KS_TE2TIME = u(10)

else if  ( flag_time == 1 ) then
    ! Constant time element (TBD)

end if

end function KS_TE2TIME


end module KUST_STI