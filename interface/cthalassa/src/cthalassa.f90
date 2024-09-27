module CTHALASSA
    ! Description:
    !    Interface to call THALASSA propagations via C
    ! 
    ! Author:
    !    Max Hallgarten La Casta
    !    Imperial College London
    !    m.hallgarten-la-casta21@imperial.ac.uk

    ! Load modules
    use, intrinsic :: iso_c_binding
    use            :: CTHALASSA_TYPES

    implicit none

    ! Strings for filepaths
    character(len=:), allocatable :: phys_path
    character(len=:), allocatable :: earth_path
    character(len=:), allocatable :: kernel_path
    character(len=:), allocatable :: eop_path

    contains

        subroutine THALASSA_OPEN(model, paths) BIND(C)
            ! Description:
            !    Opens the THALASSA interface
            ! 
            ! Author:
            !    Max Hallgarten La Casta
            !    Imperial College London
            !    m.hallgarten-la-casta21@imperial.ac.uk

            ! Import THALASSA modules
            use PHYS_CONST, only: READ_PHYS
            use NSGRAV,     only: INITIALIZE_NSGRAV, INITIALIZE_EOP
            use SUN_MOON,   only: INITIALIZE_LEGENDRE, GslSun, GslMoon
            use SETTINGS,   only: isun, imoon, iephem

            ! Subroutine arguments
            type(THALASSA_PHYSICALMODEL_STRUCT), intent(in) :: model
            type(THALASSA_PATHS_STRUCT),         intent(in) :: paths

            ! Dump log to NULL
            call OPEN_NULL_LOG()

            ! Load physical model
            call LOAD_PHYSICALMODEL(model)

            ! Allocate memory for file paths
            allocate(character(paths%phys_path_len) :: phys_path)
            allocate(character(paths%earth_path_len) :: earth_path)
            allocate(character(paths%kernel_path_len) :: kernel_path)
            allocate(character(paths%eop_path_len) :: eop_path)

            ! Parse strings
            call PTR_TO_STR(paths%phys_path,   paths%phys_path_len,   phys_path)
            call PTR_TO_STR(paths%earth_path,  paths%earth_path_len,  earth_path)
            call PTR_TO_STR(paths%kernel_path, paths%kernel_path_len, kernel_path)
            call PTR_TO_STR(paths%eop_path,    paths%eop_path_len,    eop_path)

            ! Load physical model data
            call READ_PHYS(phys_path)

            ! Load Earth model data
            call INITIALIZE_NSGRAV(earth_path)

            ! Load EOP data
            call INITIALIZE_EOP(eop_path)

            ! Initialize Legendre coefficients, if needed
            if (isun > 1) then
                call INITIALIZE_LEGENDRE(isun, GslSun)
            end if
            if (imoon > 1) then
                call INITIALIZE_LEGENDRE(imoon, GslMoon)
            end if

            ! Load SPICE kernels
            if (iephem == 1) then
                call FURNSH(kernel_path)
            end if
        end subroutine THALASSA_OPEN

        subroutine THALASSA_CLOSE() BIND(C)
            ! Description:
            !    Closes the THALASSA interface
            ! 
            ! Author:
            !    Max Hallgarten La Casta
            !    Imperial College London
            !    m.hallgarten-la-casta21@imperial.ac.uk

            ! Import THALASSA modules
            use SETTINGS, only: isun, imoon, iephem
            use NSGRAV,   only: DEINITIALIZE_NSGRAV, DEINITIALIZE_EOP
            use SUN_MOON, only: DEINITIALIZE_LEGENDRE, GslSun, GslMoon

            ! Deallocate memory for file paths
            if (allocated(phys_path)) deallocate(phys_path)
            if (allocated(earth_path)) deallocate(earth_path)
            if (allocated(kernel_path)) deallocate(kernel_path)
            if (allocated(eop_path)) deallocate(eop_path)

            ! Deallocate memory for Earth model data
            call DEINITIALIZE_NSGRAV()

            ! Deallocate memory for EOP
            call DEINITIALIZE_EOP()

            ! Deallocate memory for Legendre coefficients
            if (isun > 1) then
                call DEINITIALIZE_LEGENDRE(GslSun)
            end if
            if (imoon > 1) then
                call DEINITIALIZE_LEGENDRE(GslMoon)
            end if
            
            ! Unload SPICE kernels
            if (iephem == 1) then
                call UNLOAD(kernel_path)
            end if

            ! Close log
            call CLOSE_NULL_LOG()
        end subroutine THALASSA_CLOSE

        subroutine THALASSA_RUN(ntime, times, inputstate, outputstates, object, propagator) bind(C)
            ! Description:
            !    Runs a THALASSA propagation with verbose output
            ! 
            ! Author:
            !    Max Hallgarten La Casta
            !    Imperial College London
            !    m.hallgarten-la-casta21@imperial.ac.uk

            ! Load Fortran modules
            use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan

            ! Load THALASSA modules
            use KINDS,       only: dk
            use AUXILIARIES, only: MJD0, MJDvector, useMJDVector
            use PROPAGATE,   only: DPROP_REGULAR

            ! Subroutine parameters
            integer(c_size_t), intent(in)                   :: ntime
            real(c_double), intent(in), dimension(ntime)    :: times
            real(c_double), intent(in)                      :: inputstate(1:6)
            real(c_double), intent(out), dimension(6*ntime) :: outputstates
            type(THALASSA_OBJECT_STRUCT), intent(in)        :: object
            type(THALASSA_PROPAGATOR_STRUCT), intent(in)    :: propagator

            ! Locals
            ! Integration span and dt [solar days]
            real(dk) :: tspan, tstep
            ! Initial conditions (EMEJ2000)
            real(dk) :: R0(1:3), V0(1:3)
            ! Measurement of CPU time, diagnostics
            integer  :: exitcode
            ! Trajectory and number of propagated states
            real(dk), allocatable :: cart(:, :)
            integer(c_size_t)     :: ncart
            ! Function calls and integration steps
            integer :: int_steps, tot_calls

            ! Fill output with NaNs
            outputstates = ieee_value(0., ieee_quiet_nan)

            ! Load propagator settings (ignore tspan and tstep from propagator structure)
            call LOAD_PROPAGATOR(propagator)

            ! Load object properties
            call LOAD_OBJECT(object)

            ! Force use of time vector 
            useMJDVector = 1

            ! Load times
            MJD0 = times(1)
            MJDvector = times
            tspan = times(size(times)) - times(1)
            tstep = times(2) - times(1) ! use size of first step to avoid allocation issues

            ! Load initial conditions
            R0 = inputstate(1:3)
            V0 = inputstate(4:6)

            ! Propagate orbit
            call DPROP_REGULAR(R0, V0, tspan, tstep, cart, int_steps, tot_calls, exitcode)

            ! Extract number of times propagated
            ncart = size(cart, 1)

            ! Check for early termination
            if ( ncart .lt. ntime ) then
                ! Remove last state vector as early termination results
                ! in inconsistency with the time vector
                ncart = ncart - 1
            end if

            ! Copy output to the output matrix
            outputstates(1:6*ncart) = reshape(transpose(cart(:, 2:7)), (/ 6*ncart /))
        end subroutine

        subroutine PTR_TO_STR(ptr, ptr_len, str)
            ! Description:
            !    Converts a C char array pointer to a Fortran string
            ! 
            ! Author:
            !    Max Hallgarten La Casta
            !    Imperial College London
            !    m.hallgarten-la-casta21@imperial.ac.uk

            ! Subroutine parameters
            character(kind=c_char), intent(in)  :: ptr(*)
            integer(c_size_t),      intent(in)  :: ptr_len
            character(len=ptr_len), intent(out) :: str

            ! Locals
            integer(c_size_t) :: i

            ! Iterate through characters
            do i = 1, ptr_len
                str(i:i) = ptr(i)
            end do
        end subroutine PTR_TO_STR

        subroutine OPEN_NULL_LOG()
            ! Description:
            !    Point the log file id to the null
            ! 
            ! Author:
            !    Max Hallgarten La Casta
            !    Imperial College London
            !    m.hallgarten-la-casta21@imperial.ac.uk

            ! Load THALASSA modules
            use IO, only: id_log

            ! Locals
            character(*), parameter :: nulunix='/dev/null', nulwin='NUL'
            integer                 :: ios

            ! Open null file
            ! Modified from: https://scivision.github.io/fortran2018-examples/sourcefile/devnull.f90.html (accessed 29/11/2022)
            open(unit=id_log, file=nulunix, status='old', iostat=ios, action='write')
            if (ios /= 0) open(unit=id_log, file=nulwin, status='old', iostat=ios, action='write')
            if (ios /= 0) error stop 'could not open a NULL file handle'
        end subroutine OPEN_NULL_LOG

        subroutine CLOSE_NULL_LOG()
            ! Description:
            !    Closes the log file
            ! 
            ! Author:
            !    Max Hallgarten La Casta
            !    Imperial College London
            !    m.hallgarten-la-casta21@imperial.ac.uk

            ! Load THALASSA modules
            use IO, only: id_log

            ! Close log
            close(id_log)
        end subroutine CLOSE_NULL_LOG

        subroutine LOAD_PHYSICALMODEL(model)
            ! Description:
            !    Load physical model parameters
            ! 
            ! Author:
            !    Max Hallgarten La Casta
            !    Imperial College London
            !    m.hallgarten-la-casta21@imperial.ac.uk

            ! Import THALASSA modules
            use SETTINGS, only: insgrav, isun, imoon, idrag, iF107, iSRP, iephem, gdeg, gord, Mord, Sord

            ! Subroutine parameters
            type(THALASSA_PHYSICALMODEL_STRUCT), intent(in) :: model

            ! Set physical model parameters
            insgrav = model%insgrav
            isun    = model%isun
            imoon   = model%imoon
            idrag   = model%idrag
            iF107   = model%iF107
            iSRP    = model%iSRP
            iephem  = model%iephem
            gdeg    = model%gdeg
            gord    = model%gord

            ! Set Legendre expansion orders
            Sord = isun
            Mord = imoon
        end subroutine LOAD_PHYSICALMODEL

        subroutine LOAD_PROPAGATOR(propagator, tspan, tstep)
            ! Description:
            !    Loads propagator parameters
            ! 
            ! Author:
            !    Max Hallgarten La Casta
            !    Imperial College London
            !    m.hallgarten-la-casta21@imperial.ac.uk

            ! Import THALASSA modules
            use KINDS,    only: dk
            use SETTINGS, only: mxstep, tol, imcoll, eqs

            ! Subroutine parameters
            type(THALASSA_PROPAGATOR_STRUCT), intent(in)            :: propagator
            real(dk),                         intent(out), optional :: tspan, tstep

            ! Set integrator settings
            mxstep = propagator%mxstep
            tol    = propagator%tol
            if(present(tspan)) tspan  = propagator%tspan
            if(present(tstep)) tstep  = propagator%tstep
            imcoll = propagator%imcoll

            ! Set equations of motion
            eqs = propagator%eqs
        end subroutine LOAD_PROPAGATOR

        subroutine LOAD_OBJECT(object)
            ! Description:
            !    Loads object parameters
            ! 
            ! Author:
            !    Max Hallgarten La Casta
            !    Imperial College London
            !    m.hallgarten-la-casta21@imperial.ac.uk

            ! Load THALASSA modules
            use PHYS_CONST, only: SCMass, ADrag, ASRP, CD, CR, A2M_Drag, A2M_SRP

            ! Subroutine parameters
            type(THALASSA_OBJECT_STRUCT), intent(in) :: object

            ! Set object properties
            SCmass = object%mass
            ADrag  = object%area_drag
            ASRP   = object%area_srp
            CD     = object%cd
            CR     = object%cr
            
            ! Calculate area-to-mass ratios
            A2M_Drag = ADrag/SCMass
            A2M_SRP  = ASRP/SCmass
        end subroutine LOAD_OBJECT

        subroutine LOAD_STATE(state, mjd, RV)
            ! Description:
            !    Loads an object state from a struct
            ! 
            ! Author:
            !    Max Hallgarten La Casta
            !    Imperial College London
            !    m.hallgarten-la-casta21@imperial.ac.uk

            ! Load THALASSA modules
            use KINDS, only: dk

            ! Subroutine parameters
            type(THALASSA_STATE_STRUCT), intent(in)  :: state
            real(dk),                    intent(out) :: mjd, RV(1:6)

            ! Load epoch and state
            mjd = state%mjd
            RV = state%RV
        end subroutine LOAD_STATE

        subroutine SAVE_STATE(state, mjd, RV)
            ! Description:
            !    Saves an object state into a struct
            ! 
            ! Author:
            !    Max Hallgarten La Casta
            !    Imperial College London
            !    m.hallgarten-la-casta21@imperial.ac.uk

            ! Load THALASSA modules
            use KINDS, only: dk

            ! Subroutine parameters
            type(THALASSA_STATE_STRUCT), intent(out) :: state
            real(dk),                    intent(in)  :: mjd, RV(1:6)

            ! Load epoch and state
            state%mjd = mjd
            state%RV = RV
        end subroutine SAVE_STATE

end module CTHALASSA
