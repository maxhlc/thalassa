module NSGRAV
! Description:
!    Contains procedures for the calculation of perturbations due to the non-
!    sphericity of the main body. INITIALIZE_NSGRAV reads main body data from a
!    data file, and initializes coefficient matrices.
!    The calculation of the perturbing potential, perturbing acceleration, and
!    the time derivative of the potential in the body-fixed frame (the latter is
!    needed in regularized formulations) takes place in PINES_NSG.
!    NORMFACT gives the normalization factor for the gravitational coefficients.
! 
! Author:
!    Davide Amato
!    The University of Arizona
!    davideamato@email.arizona.edu
! 
! Revisions:
!    180806: Overhaul and implementation of Pines method.
!    181204: Use GMST (IERS 2006 conventions) rather than ERA. Consider Earth
!            Rotation Rate derived from IERS 2006 conventions.
! 
! ==============================================================================

! MODULES
use KINDS,      only: dk
use SETTINGS,   only: gdeg,gord
use IO,         only: id_earth, id_eop
use PHYS_CONST, only: qk, GE, RE, flatt, omegaE, secsPerDay, secsPerSidDay, twopi, ERR_constant, DAS2R, DMAS2R

implicit none

! VARIABLES
! Spherical harmonics (unnormalized) 
integer               ::  maxDeg, maxOrd
real(dk),allocatable  ::  Cnm(:,:),Snm(:,:)

! Pines algorithm arrays
real(dk),allocatable  ::  Anm(:,:),Dnm(:,:),Enm(:,:),Fnm(:,:),Gnm(:,:)
real(dk),allocatable  ::  Rm(:),Im(:),Pn(:)
real(dk),allocatable  ::  Aux1(:),Aux2(:),Aux3(:),Aux4(:,:)

! earth orientation parameters array
real(dk),allocatable  ::  eop(:,:)

! Pre-computed rotation matrices
logical               :: usePrecomputed = .FALSE.
integer               :: nrotation
real(dk), allocatable :: precomputedMJD(:)
real(dk), allocatable :: precomputedBPN_T(:, :, :)
real(dk), allocatable :: precomputedW_T(:, :, :)

contains

subroutine INITIALIZE_EOP(eopfile)

! formal parameters
character(len=*), intent(in)  ::  eopfile

! local parameters
integer        :: io, nlines, i
character(200) :: line
character(200) :: format
real(dk)       :: mjd, pmx, pmy, UT1_UTC, LOD, dx00, dy00

! set number of lines to zero
nlines = 0

! count lines
open(unit=id_eop,file=trim(eopfile),status='old',action='read')
do
  read(id_eop,*,iostat=io)
  if (io/=0) exit
  nlines = nlines + 1
end do
close(id_eop)

! rewind 
rewind(id_eop)

! allocate
if (allocated(eop)) deallocate(eop)
allocate(eop(1:nlines, 1:7))

! read in eop data
open(unit=id_eop,file=trim(eopfile),status='old',action='read')
do i = 1,nlines

  ! read line
  read(id_eop, '(7X,F8.2,3X,F9.6,10X,F9.6,12X,F10.7,11X,F7.4,11X,F9.3,10X,F9.3,67X)')mjd,pmx,pmy,UT1_UTC,LOD,dx00,dy00

  ! add to EOP array
  eop(i,1) = mjd
  eop(i,2) = pmx
  eop(i,3) = pmy
  eop(i,4) = UT1_UTC
  eop(i,5) = LOD
  eop(i,6) = dx00
  eop(i,7) = dy00

enddo

close(id_eop)

end subroutine INITIALIZE_EOP

subroutine DEINITIALIZE_EOP()

if (allocated(eop)) deallocate(eop)

end subroutine DEINITIALIZE_EOP

subroutine INITIALIZE_ROTATION(MJDstart, MJDend, period)
  implicit none

  ! Parameters
  real(dk), intent(in) :: MJDstart
  real(dk), intent(in) :: MJDend
  real(dk), intent(in) :: period

  ! Locals
  integer  :: irotation
  real(dk) :: tempMJD, MJD_TT, UTCjd1, UTCjd2, TTjd1, TTjd2
  integer  :: eopIdx
  real(dk) :: tempBPN_T(3,3), tempW_T(3,3)

  ! Calculate number of periods
  ! NOTE: adds an additional point after the final time
  nrotation = int((MJDend - MJDstart) / period) + 2

  ! Ensure arrays are not initialised
  call DEINITIALIZE_ROTATION()

  ! Allocate arrays
  allocate(precomputedMJD(1:nrotation))
  allocate(precomputedBPN_T(1:3, 1:3, 1:nrotation))
  allocate(precomputedW_T(1:3, 1:3, 1:nrotation))

  ! Iterate through periods
  do irotation = 1, nrotation
    ! Calculate MJD
    tempMJD = MJDstart + (irotation - 1) * period

    ! Preprocess dates and EOP index
    call UTCtoPARAMS(tempMJD, MJD_TT, UTCjd1, UTCjd2, TTjd1, TTjd2, eopIdx)

    ! Calculate rotation matrices
    call GCRStoCIRS_MATRIX(TTjd1, TTjd2, eopIdx, tempBPN_T)
    call TIRStoITRS_MATRIX(TTjd1, TTjd2, eopIdx, tempW_T)

    ! Store results
    precomputedMJD(irotation) = tempMJD
    precomputedBPN_T(1:3, 1:3, irotation) =  tempBPN_T
    precomputedW_T(1:3, 1:3, irotation) = tempW_T
  end do

  ! Enable precomputed rotation matrices
  usePrecomputed = .TRUE.
end subroutine INITIALIZE_ROTATION

subroutine DEINITIALIZE_ROTATION()
  ! Disable precomputed rotation matrices
  usePrecomputed = .FALSE.

  ! Deallocate rotation dates
  if (allocated(precomputedMJD)) deallocate(precomputedMJD)

  ! Deallocate rotation matrices
  if (allocated(precomputedBPN_T)) deallocate(precomputedBPN_T)
  if (allocated(precomputedW_T)) deallocate(precomputedW_T)
end subroutine DEINITIALIZE_ROTATION

subroutine GET_ROTATION(MJD_UTC, BPN_T, W_T)
  implicit none

  ! Formals
  real(dk), intent(in)  :: MJD_UTC
  real(dk), intent(out) :: BPN_T(3,3), W_T(3,3)

  ! Locals
  integer  :: irotation
  real(dk) :: MJD1, MJD2
  real(dk) :: fraction

  ! Extract index
  ! TODO: replace brute force with better search algorithm?
  do irotation = 1, nrotation - 1
    ! Extract dates
    MJD1 = precomputedMJD(irotation)
    MJD2 = precomputedMJD(irotation + 1)

    ! Check current date is within interval
    if ((MJD_UTC .ge. MJD1) .and. (MJD_UTC .lt. MJD2)) then
      ! Break do loop
      exit
    end if
  end do

  ! Interpolate matrices
  fraction = (MJD_UTC - MJD1) / (MJD2 - MJD1)
  BPN_T = (1.0 - fraction) * precomputedBPN_T(:, :, irotation) + fraction * precomputedBPN_T(:, :, irotation + 1)
  W_T = (1.0 - fraction) * precomputedW_T(:, :, irotation) + fraction * precomputedW_T(:, :, irotation + 1)
end subroutine GET_ROTATION

subroutine INITIALIZE_NSGRAV(earthFile)
! Description:
!    Reads Earth gravity data from a text file. Initializes the gravitational
!    parameter, equatorial radius, flattening, rotational velocity, spherical
!    harmonics coefficients, and auxiliary matrices for Pines' algorithm. The
!    latter computes the Earth potential in Cartesian coordinates, avoiding
!    singularities due to the spherical harmonics.
!    Part of this subroutine is due to Hodei Urrutxua (Universidad Rey Juan
!    Carlos, Madrid, Spain) and Claudio Bombardelli (Universidad Politécnica de
!    Madrid, Madrid, Spain).
! 
! Author:
!    Davide Amato
!    The University of Arizona
!    davideamato@email.arizona.edu
! 
! Revisions:
!    180806: Subroutine created from part of READ_PHYS().
! 
! ==============================================================================

! Arguments
character(len=*),intent(in)  ::  earthFile

! Locals
character(len=72)  ::  dummy
real(dk)           ::  invFlatt
integer            ::  i,j,l,m,n


! ==============================================================================

open(unit=id_earth,file=trim(earthFile),status='old',action='read')
read(id_earth,'(a)') (dummy, i=1,4)
read(id_earth,'(a37,i3)') dummy, maxDeg
read(id_earth,'(a37,i3)') dummy, maxOrd
read(id_earth,'(a36,e22.15)') dummy, GE
read(id_earth,'(a36,e22.15)') dummy, RE
read(id_earth,'(a36,e22.15)') dummy, invFlatt
read(id_earth,'(a36,e22.15)') dummy, omegaE

flatt    = 1._dk/invFlatt

! Read and de-normalize spherical harmonics coefficients and auxiliary matrices
! for Pines' algorithm.
! Initialize
l = 1; m = 0;
allocate(Cnm(1:maxDeg,0:maxDeg)); Cnm = 0._dk
allocate(Snm(1:maxDeg,0:maxDeg)); Snm = 0._dk
allocate(Anm(0:maxDeg+2, 0:maxDeg+2))
allocate(Dnm(1:maxDeg,   0:maxDeg))
allocate(Gnm(1:maxDeg,   0:maxDeg))
allocate(Enm(1:maxDeg,   1:maxDeg))
allocate(Fnm(1:maxDeg,   1:maxDeg))
allocate(Rm(0:maxDeg)    )
allocate(Im(0:maxDeg)    )
allocate(Pn(0:maxDeg + 1))
allocate(Aux1(1:maxDeg+1))
allocate(Aux2(1:maxDeg+1))
allocate(Aux3(1:maxDeg+1))
allocate(Aux4(1:maxDeg+1, 0:maxDeg+1))

read(id_earth,'(a)') (dummy, i=1,2)
do i=1,maxDeg
  do j=0,minval([i,maxOrd])
    read(id_earth,'(2(1x,i2),2(1x,e24.17))') l,m,Cnm(i,j),Snm(i,j)
    Cnm(i,j) = Cnm(i,j)/NORMFACT(i,j)
    Snm(i,j) = Snm(i,j)/NORMFACT(i,j)
  end do
end do

close(id_earth)

! Fill coefficient arrays for Pines algorithm
Anm(:,:) = 0._dk
Dnm(:,:) = 0._dk
Gnm(:,:) = 0._dk
Enm(:,:) = 0._dk
Fnm(:,:) = 0._dk
Anm(0,0) = 1._dk
Anm(1,1) = 1._dk
do n = 1,maxDeg + 1
  Aux1(n) = 2._dk*n + 1._dk
  Aux2(n) = Aux1(n) / (n+1._dk)
  Aux3(n) = n / (n+1._dk)
  do m = 0, n-1
    Aux4(n,m) = (n+m+1._dk)
  end do
end do

! Earth Rotation Rate (revolutions per tropical day)
secsPerSidDay = twopi/omegaE
ERR_constant = secsPerDay/secsPerSidDay

end subroutine INITIALIZE_NSGRAV

subroutine DEINITIALIZE_NSGRAV()
! Description:
!    Deallocate memory used for Earth gravity data
! 
! Author:
!    Max Hallgarten La Casta
!    Imperial College London
!    m.hallgarten-la-casta21@imperial.ac.uk
! 
! ==============================================================================

! Deallocate memory
if (allocated(Cnm)) deallocate(Cnm)
if (allocated(Snm)) deallocate(Snm)
if (allocated(Anm)) deallocate(Anm)
if (allocated(Dnm)) deallocate(Dnm)
if (allocated(Gnm)) deallocate(Gnm)
if (allocated(Enm)) deallocate(Enm)
if (allocated(Fnm)) deallocate(Fnm)
if (allocated(Rm)) deallocate(Rm)
if (allocated(Im)) deallocate(Im)
if (allocated(Pn)) deallocate(Pn)
if (allocated(Aux1)) deallocate(Aux1)
if (allocated(Aux2)) deallocate(Aux2)
if (allocated(Aux3)) deallocate(Aux3)
if (allocated(Aux4)) deallocate(Aux4)
end subroutine


function NORMFACT(l,m)
! Description:
!    Normalization factor for the spherical harmonics coefficients and
!    associated Legendre functions:
! 
!    sqrt( (l + m)! / ( (2 - delta_{0,m}) * (2n  + 1) * (n - m)! ) )
! 
! Reference:
!    [1] Montenbruck, O., Gill, E., "Satellite Orbits", p. 58, Springer, 2000.
! 
! Author:
!    Davide Amato
!    The University of Arizona
!    davideamato@email.arizona.edu
! 
! ==============================================================================

! Arguments and function definition
integer,intent(in)  ::  l,m
real(qk)            ::  NORMFACT
! Locals
real(qk)            ::  lr,mr
real(qk)            ::  kron
real(qk)            ::  numer,denom

! ==============================================================================
lr = real(l,qk)
mr = real(m,qk)

numer = gamma(lr + mr + 1._qk)

if (m == 0) then
  kron = 1._qk
else
  kron = 0._qk
end if

denom = (2._qk - kron) * (2._qk*lr + 1._qk) * gamma(lr - mr + 1._qk)

NORMFACT = sqrt(numer/denom)

end function NORMFACT

subroutine GCRStoCIRS_MATRIX(TTjd1, TTjd2, eopIdx, BPN_T)
  ! No implicit variables
  implicit none
  
  ! Formals
  real(dk), intent(in)  :: TTjd1, TTjd2
  integer,  intent(in)  :: eopIdx
  real(dk), intent(out) :: BPN_T(3,3)
  
  ! Locals
  real(dk)  :: dx00, dy00
  real(dk)  :: X, Y, S
  
  ! CIP offsets wrt IAU 20000
  dx00 = eop(eopIdx, 5) * DMAS2R
  dy00 = eop(eopIdx, 6) * DMAS2R
  
  ! CIP and CIO, IAU 2000A (Represents the PN matrix in Vallado)
  call iau_XYS00B(TTjd1, TTjd2, X, Y, S)
  
  ! Add CIP corrections
  X = X + dx00;
  Y = Y + dy00;
  
  ! GCRS to CIRS matrix (celestial to intermediate matrix [BPN]' = [N]'[P]'[B]')
  call iau_C2IXYS(X, Y, S, BPN_T)
end subroutine GCRStoCIRS_MATRIX

subroutine ERA_MATRIX(UTCjd1, UTCjd2, eopIdx, R_T)
  ! No implicit variables
  implicit none

  ! External functions
  real(dk), external :: iau_ERA00

  ! Formals
  real(dk), intent(in)  :: UTCjd1, UTCjd2
  integer,  intent(in)  :: eopIdx
  real(dk), intent(out) :: R_T(3,3)

  ! Locals
  real(dk) :: era
  real(dk) :: UT1_UTC, UT1jd1, UT1Jd2
  integer  :: J

  ! Convert UTC to UT1
  UT1_UTC = eop(eopIdx, 4) 
  call iau_UTCUT1(UTCjd1, UTCjd2, UT1_UTC, UT1jd1, UT1Jd2, J)

  ! Earth Rotation Angle  
  era = iau_ERA00(UT1jd1, UT1jd2)

  ! Rotation matrix about pole
  R_T(1,1) = 1
  R_T(1,2) = 0
  R_T(1,3) = 0
  R_T(2,1) = 0
  R_T(2,2) = 1
  R_T(2,3) = 0
  R_T(3,1) = 0
  R_T(3,2) = 0
  R_T(3,3) = 1
  call iau_RZ(era, R_T)
end subroutine ERA_MATRIX

subroutine TIRStoITRS_MATRIX(TTjd1, TTjd2, eopIdx, W_T)
  ! No implicit variables
  implicit none
  
  ! External functions
  real(dk), external :: iau_SP00

  ! Formals
  real(dk), intent(in)  :: TTjd1, TTjd2
  integer,  intent(in)  :: eopIdx
  real(dk), intent(out) :: W_T(3,3)
  
  ! Locals
  real(dk)  :: Xp, Yp

  ! Polar motion
  Xp = eop(eopIdx, 2) * DAS2R
  Yp = eop(eopIdx, 3) * DAS2R
  
  ! Polar motion matrix (TIRS->ITRS, IERS 2003) (W_T matrix)
  call iau_POM00(Xp, Yp, iau_SP00(TTjd1, TTjd2), W_T)
end subroutine TIRStoITRS_MATRIX

subroutine UTCtoPARAMS(MJD_UTC, MJD_TT, UTCjd1, UTCjd2, TTjd1, TTjd2, eopIdx)
  ! External subroutines
  use PHYS_CONST, only: UTC2TT

  ! Formals
  real(dk), intent(in)  :: MJD_UTC
  real(dk), intent(out) :: UTCjd1, UTCjd2
  real(dk), intent(out) :: MJD_TT
  real(dk), intent(out) :: TTjd1, TTjd2
  integer,  intent(out) :: eopIdx

  ! Locals
  integer :: imjdCurr, imjd0, imjdDiff

  ! Get UTC and TT
  MJD_TT  = UTC2TT(MJD_UTC)
  UTCjd1  = 2400000.5_dk
  UTCjd2  = MJD_UTC
  TTjd1   = 2400000.5_dk
  TTjd2   = MJD_TT

  ! Get EOP data index
  imjdCurr = int(MJD_UTC)
  imjd0 = int(eop(1,1))

  ! Check for unprovided date
  if (imjdCurr < imjd0) then
    write(*,*)imjdCurr,imjd0,'Requested date is before start of Earth Orientation Data'
  endif

  ! Index of requested MJD
  imjdDiff = imjdCurr - imjd0
  eopIdx = 1 + imjdDIff
end subroutine UTCtoPARAMS

subroutine GCRFtoITRF_MATRIX_(W_T, R_T, BPN_T, gcrf_itrf)
  ! No implict variables
  implicit none

  ! Formals
  real(dk), intent(in)  :: W_T(3,3), R_T(3,3), BPN_T(3,3)
  real(dk), intent(out) :: gcrf_itrf(3,3)

  ! Locals
  real(dk) :: temp(3,3)

  ! Multiply [R]',[BPN]',and [W]'
  call iau_RXR(R_T, BPN_T, temp)
  call iau_RXR(W_T, temp, gcrf_itrf)
end subroutine GCRFtoITRF_MATRIX_

subroutine GCRFtoITRF_MATRIX(MJD_UTC, gcrf_itrf)
  ! No implicit variables
  implicit none

  ! Formals
  real(dk), intent(in)  :: MJD_UTC
  real(dk), intent(out) :: gcrf_itrf(3,3)

  ! Locals
  integer   :: eopIdx
  real(dk)  :: UTCjd1, UTCjd2, TTjd1, TTjd2, MJD_TT
  real(dk)  :: BPN_T(3,3), R_T(3,3), W_T(3,3)

  ! Preprocess dates and EOP index
  call UTCtoPARAMS(MJD_UTC, MJD_TT, UTCjd1, UTCjd2, TTjd1, TTjd2, eopIdx)

  ! GCRS to CIRS matrix (celestial to intermediate matrix [BPN]' = [N]'[P]'[B]')
  call GCRStoCIRS_MATRIX(TTjd1, TTjd2, eopIdx, BPN_T)

  ! ERA matrix
  call ERA_MATRIX(UTCjd1, UTCjd2, eopIdx, R_T)

  ! Polar motion matrix (TIRS->ITRS, IERS 2003) (W_T matrix)
  call TIRStoITRS_MATRIX(TTjd1, TTjd2, eopIdx, W_T)

  ! Multiply [R]',[BPN]',and [W]'
  call GCRFtoITRF_MATRIX_(W_T, R_T, BPN_T, gcrf_itrf)
end subroutine GCRFtoITRF_MATRIX

subroutine GCRFtoITRF(MJD_UTC, Rgcrf, Ritrf, gcrf_itrf)
  implicit none

  ! Parameters
  real(dk), intent(in)  :: MJD_UTC
  real(dk), intent(in)  :: Rgcrf(3)
  real(dk), intent(out) :: Ritrf(3)
  real(dk), intent(out) :: gcrf_itrf(3,3)

  ! Locals
  integer   :: eopIdx
  real(dk)  :: UTCjd1, UTCjd2, TTjd1, TTjd2, MJD_TT
  real(dk)  :: BPN_T(3,3), R_T(3,3), W_T(3,3)

  ! Get rotation matrix
  if (usePrecomputed) then
    ! Preprocess dates and EOP index
    call UTCtoPARAMS(MJD_UTC, MJD_TT, UTCjd1, UTCjd2, TTjd1, TTjd2, eopIdx)
    
    ! Calculate ERA
    call ERA_MATRIX(UTCjd1, UTCjd2, eopIdx, R_T)

    ! Interpolate matrices
    call GET_ROTATION(MJD_UTC, BPN_T, W_T)

    ! Calculate rotation matrix
    call GCRFtoITRF_MATRIX_(W_T, R_T, BPN_T, gcrf_itrf)
  else
    ! Calculate rotation matrix directly
    call GCRFtoITRF_MATRIX(MJD_UTC, gcrf_itrf)
  endif

  ! Rotate position vector
  call iau_RXP(gcrf_itrf, Rgcrf, Ritrf)
end subroutine GCRFtoITRF


subroutine PINES_NSG(GM,RE,rIn,tau,FIn,pot,dPot)
! Description:
!    Compute the perturbing acceleration, perturbing potential (optional), and
!    time derivative of the perturbing potential in the body-fixed frame
!    (optional), given the position vector wrt the non-spherical body, its
!    gravitational parameter and its potential coefficients.
!    Uses the method described in the Reference to perform the calculation in
!    Cartesian coordinates, thereby avoiding associated Legendre functions and
!    polynomials. The formulation is non-singular everywhere for r > RE, where
!    RE is the radius of the non-spherical body.
! 
!    Adapted from code developed by Hodei Urrutxua (Universidad Rey Juan Carlos,
!    Madrid, Spain) and Claudio Bombardelli (Universidad Politécnica de Madrid,
!    Madrid, Spain).
! 
! Author:
!    Davide Amato
!    The University of Arizona
!    davideamato@email.arizona.edu
! 
! Reference:
!    S. Pines, "Uniform Representation of the Gravitational Potential and its
!    derivatives," AIAA J. 11 (11), pp. 1508-1511, 1973.
! 
! Revisions:
!    180806: First working version of subroutine.
!    181204: Use GMST (IERS 2006 conventions) rather than ERA. Consider Earth
!            Rotation Rate derived from IERS 2006 conventions.
!
! ==============================================================================

! Use associations
use AUXILIARIES, only: MJD0, TU, T2MJD
use PHYS_CONST,  only: delta_JD_MJD
use PHYS_CONST,  only: ERR_IAU06, UTC2TT

! Arguments
real(dk),intent(in)   ::  GM             ! Gravitational parameter
real(dk),intent(in)   ::  RE             ! Equatorial radius
real(dk),intent(in)   ::  rIn(1:3)       ! Position in the inertial frame
real(dk),intent(in)   ::  tau            ! Physical time
real(dk),intent(out),optional  ::  FIn(1:3)  ! Perturbing acceleration in the inertial frame
real(dk),intent(out),optional  ::  pot   ! Perturbing potential
real(dk),intent(out),optional  ::  dPot  ! Time derivative of the potential in body-fixed frame
! Locals
real(dk)  ::  rNorm,rNormSq  ! Norm of position vector and square
real(dk)  ::  s,t,u  ! Direction cosines in the body-fixed frame
real(dk)  ::  rho    ! = equatorial radius / r 
real(dk)  ::  F(1:3),a1,a2,a3,a4  ! Acceleration and its components 
integer   ::  n,m    ! Harmonic indices
logical   ::  skip_EFG
! GMST-related quantities
real(dk)  ::  MJD_UTC, MJD_TT      ! UTC and TT dates
real(dk)  ::  GMST,cosGMST,sinGMST ! GMST and its trig functions
real(dk)  ::  ERR, ERR_nd          ! Earth Rotation Rate [rad/s, -]

! SOFA routines
real(dk) :: iau_GMST06

! gcrf -> itrf conversion
real(dk) :: Rgcrf(3), Ritrf(3)
real(dk) :: Vgcrf(3), Vitrf(3)
real(dk) :: Agcrf(3), Aitrf(3)
real(dk) :: dummy
real(dk) :: gcrf_itrf(3,3)
real(dk) :: itrf_gcrf(3,3)

! ==============================================================================

rNormSq = dot_product(rIn,rIn)
rNorm = sqrt(rNormSq)
rho = RE/rNorm

! ==============================================================================
! 01. Transform from inertial to body-fixed frame
! ==============================================================================

! get UTC and TT
MJD_UTC = T2MJD(tau)
MJD_TT  = UTC2TT(MJD_UTC)

! get rotation matrix from gcrf to itrf
call GCRFtoITRF(MJD_UTC, rIn, Ritrf, gcrf_itrf)

! transpose to get rotation from itrf to gcrf
call iau_TR(gcrf_itrf, itrf_gcrf)

! non-dimensionalized satellite body-fixed coordinates
u = Ritrf(3)/rNorm
s = Ritrf(1)/rNorm
t = Ritrf(2)/rNorm

! ==============================================================================
! 02. Fill in coefficient matrices and auxiliary vectors
! ==============================================================================

! Fill A Matrix
Anm(0,0) = 1._dk
Anm(1,1) = 1._dk
Anm(1,0) = u
do n = 1, gdeg + 1
  Anm(n+1,n+1) = Aux1(n) * Anm(n,n) ! Fill the diagonal
  Anm(n+1,0) = Aux2(n) * u * Anm(n,0) - Aux3(n) * Anm(n-1,0) ! Fill the 1st column
  Anm(n+1,n) = u * Anm(n+1,n+1)    ! Fill the subdiagonal

end do
do n = 2, gdeg + 1 ! Fill remaining elements
  do m = 0, n - 2
    Anm(n+1,m+1) = Aux4(n,m) * Anm(n,m) + u * Anm(n,m+1)
  end do

end do

! Fill R, I, and P vectors
Rm(0) = 1._dk
Im(0) = 0._dk
Pn(0) = GM / rNorm
Pn(1) = rho * Pn(0)
do n = 1, gdeg
  Rm(n)   = s  * Rm(n-1) - t * Im(n-1)
  Im(n)   = s  * Im(n-1) + t * Rm(n-1)
  Pn(n+1) = rho * Pn(n)

end do

! Fill D, E, and F matrices
skip_EFG = present(pot) .and. .not.(present(FIn)) .and. .not.(present(dPot))
do m = 1, gord
  do n = m, gdeg
    Dnm(n,m) = Cnm(n,m)*Rm(m)   + Snm(n,m)*Im(m)
    if (.not.(skip_EFG)) then
      Enm(n,m) = Cnm(n,m)*Rm(m-1) + Snm(n,m)*Im(m-1)
      Fnm(n,m) = Snm(n,m)*Rm(m-1) - Cnm(n,m)*Im(m-1)

    end if
  end do

end do
do n = 1, gdeg
  Dnm(n,0) = Cnm(n,0)*Rm(0)  !+ S(n,0)*I(0) = 0
end do

! ==============================================================================
! 03. Perturbing potential
! ==============================================================================

if (present(pot)) then
  pot = 0._dk
  do m = 1, gord
    do n = m, gdeg
      pot = pot + Pn(n) * Anm(n,m) * Dnm(n,m)
		
    end do
	
  end do
  do n = 1, gdeg
  pot = pot + Pn(n) * Anm(n,0) * Dnm(n,0)
	
  end do
  pot = -pot  ! Change the sign to get the potential
end if

! ==============================================================================
! 04. Perturbing acceleration
! ==============================================================================

if (present(FIn)) then
  a1 = 0._dk; a2 = 0._dk; a3 = 0._dk; a4 = 0._dk
  do m = 1, gord
    do n = m, gdeg
      a1 = a1 + Pn(n+1) * Anm(n,m) * m * Enm(n,m)
      a2 = a2 + Pn(n+1) * Anm(n,m) * m * Fnm(n,m)
      a3 = a3 + Pn(n+1) * Anm(n,m+1)   * Dnm(n,m)
      a4 = a4 - Pn(n+1) * Anm(n+1,m+1) * Dnm(n,m)
    end do
  end do
  do n = 1, gdeg
    a3 = a3 + Pn(n+1) * Anm(n,1)     * Dnm(n,0)
    a4 = a4 - Pn(n+1) * Anm(n+1,1)   * Dnm(n,0)

  end do
  F = [a1, a2, a3] + [s, t, u] * a4
  F = F / RE

  ! Transform to inertial frame
  call iau_RXP(itrf_gcrf, F, FIn)

end if

! ==============================================================================
! 05. Time derivative of potential in body-fixed frame
! ==============================================================================

if(present(dPot)) then
  dPot = 0._dk
  ERR = ERR_IAU06(0._dk, MJD_TT)
  ERR_nd = ERR / TU
  do m = 1, gord
    do n = m, gdeg
      Gnm(n,m) = m * ( t * Enm(n,m) - s * Fnm(n,m) )
      Gnm(n,m) = ERR_nd * Gnm(n,m)
      dPot = dPot + Pn(n) * Anm(n,m) * Gnm(n,m)
      
    end do
    
  end do
  dPot = -dPot
end if

end subroutine PINES_NSG




end module NSGRAV
