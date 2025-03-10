! Description:
!> @file
!!   Example calling direct model hydrometeor scattering simulation using
!!   pre-defined particle optical properties for a MW sensor.
!
!> @brief
!!   Example calling direct model hydrometeor scattering simulation using
!!   pre-defined particle optical properties for a MW sensor.
!!
!! @details
!!   This example is most easily run using the run_example_hydro_mw_file_fwd.sh
!!   script (see the user guide).
!!
!!   This program requires the following files:
!!     the file containing input profiles (e.g. prof.dat)
!!     the file containing input hydrometeor data
!!     the RTTOV optical depth coefficient file
!!     the corresponding hydrotable file
!!
!!   The output is written to a file called example_hydro_mw_file_fwd_output.dat
!!
!!   This program is very similar to example_fwd.F90 with the addition
!!   of the hydrometeor scattering. See also example_hydro_visir_file_fwd.F90 for
!!   an example for a visible/IR sensor.
!!
!!   You may wish to base your own code on this example in which case you
!!   should edit in particular the sections of code marked as follows:
!!       !================================
!!       !======Read =====start===========
!!            code to be modified
!!       !======Read ===== end ===========
!!       !================================
!!
!
! Copyright:
!    This software was developed within the context of
!    the EUMETSAT Satellite Application Facility on
!    Numerical Weather Prediction (NWP SAF), under the
!    Cooperation Agreement dated 7 September 2021, between
!    EUMETSAT and the Met Office, UK, by one or more partners
!    within the NWP SAF. The partners in the NWP SAF are
!    the Met Office, ECMWF, DWD and MeteoFrance.
!
!    Copyright 2024, EUMETSAT, All Rights Reserved.
!
PROGRAM example_hydro_mw_file_fwd

  ! rttov_const contains useful RTTOV constants
  USE rttov_const, ONLY :            &
         errorstatus_success,        &
         errorstatus_fatal,          &
         alloc_flag,                 &
         dealloc_flag,               &
         platform_name,              &
         inst_name,                  &
         interp_rochon_wfn,          &
         hydro_index_rain,           &
         hydro_index_snow,           &
         hydro_index_graupel,        &
         hydro_index_clw,            &
         hydro_index_ciw,            &
         thermal_solver_delta_edd,   &
         cloud_overlap_2col_weighted

  ! rttov_types contains definitions of all RTTOV data types
  USE rttov_types, ONLY :     &
         rttov_chanprof,      &
         rttov_options,       &
         rttov_coefs,         &
         rttov_profile,       &
         rttov_emis_refl,     &
         rttov_transmission,  &
         rttov_radiance

  ! jpim, jprv and jplm are the RTTOV integer, real and logical KINDs
  USE rttov_kinds, ONLY : jpim, jprv, jplm

  USE rttov_unix_env, ONLY : rttov_exit

  IMPLICIT NONE

#include "rttov_direct.interface"
#include "rttov_parallel_direct.interface"
#include "rttov_read_coefs.interface"
#include "rttov_dealloc_coefs.interface"
#include "rttov_alloc_direct.interface"
#include "rttov_user_check_options.interface"
#include "rttov_print_opts.interface"
#include "rttov_print_profile.interface"
#include "rttov_skipcommentline.interface"

  !--------------------------
  !
  INTEGER(KIND=jpim), PARAMETER :: iup   = 20   ! unit for input profile file
  INTEGER(KIND=jpim), PARAMETER :: ioout = 21   ! unit for output

  ! RTTOV variables/structures
  !====================
  TYPE(rttov_options)              :: opts                     ! Options structure
  TYPE(rttov_coefs)                :: coefs                    ! Coefficients structure
  TYPE(rttov_chanprof),    POINTER :: chanprof(:)    => NULL() ! Input channel/profile list
  TYPE(rttov_profile),     POINTER :: profiles(:)    => NULL() ! Input profiles
  TYPE(rttov_emis_refl),   POINTER :: emis_refl(:)   => NULL() ! Input/output emissivities/reflectances
  TYPE(rttov_transmission)         :: transmission             ! Output transmittances
  TYPE(rttov_radiance)             :: radiance                 ! Output radiances

  INTEGER(KIND=jpim)               :: errorstatus              ! Return error status of RTTOV subroutine calls

  INTEGER(KIND=jpim) :: alloc_status
  CHARACTER(LEN=25)  :: NameOfRoutine = 'example_hydro_mw_file_fwd'
  INTEGER(KIND=jpim) :: nsurfaces = 1

  ! variables for input
  !====================
  CHARACTER(LEN=256) :: coef_filename
  CHARACTER(LEN=256) :: hydrotable_filename
  CHARACTER(LEN=256) :: prof_filename
  CHARACTER(LEN=256) :: hydro_prof_filename
  INTEGER(KIND=jpim) :: nthreads
  INTEGER(KIND=jpim) :: nlevels, nlayers
  INTEGER(KIND=jpim) :: nprof
  INTEGER(KIND=jpim) :: nchannels
  INTEGER(KIND=jpim) :: nchanprof
  INTEGER(KIND=jpim), ALLOCATABLE :: channel_list(:)

  INTEGER(KIND=jpim) :: isurf, iprof, ichan, j, joff, ilay
  INTEGER            :: ios

  !------------------------------------------------! Vito
  REAL(KIND=jprv)              :: input_emissivity ! Vito
  REAL(KIND=jprv), ALLOCATABLE :: emissivity_val(:)! Vito
  REAL(KIND=jprv)              :: ifac             ! Vito
  REAL(KIND=jprv)              :: zenith           ! Vito
  !------------------------------------------------! Vito
  !-------------------------------------------------------------------------! V
  Character (len=100)  :: fmt,fout1,path ! Variables for storing the output ! V
  REAL(Kind=jprv)      :: dummy    ! For data type convertion               ! V
  !-------------------------------------------------------------------------! V

  !- End of header --------------------------------------------------------

  ! The usual steps to take when running RTTOV are as follows:
  !   1. Specify required RTTOV options
  !   2. Read coefficients
  !   3. Allocate RTTOV input and output structures
  !   4. Set up the chanprof array with the channels/profiles to simulate
  !   5. Read input profile(s)
  !   6. Set up surface emissivity and/or reflectance
  !   7. Call rttov_direct and store results
  !   8. Deallocate all structures and arrays

  ! If nthreads is greater than 1 the parallel RTTOV interface is used.
  ! To take advantage of multi-threaded execution you must have compiled
  ! RTTOV with openmp enabled. See the user guide and the compiler flags.

  errorstatus = errorstatus_success

  !=====================================================
  !========== Interactive inputs == start ==============

  WRITE(0,*) 'enter path of coefficient file'
  READ(*,*) coef_filename
  WRITE(0,*) 'enter path of hydrotable file'
  READ(*,*) hydrotable_filename
  WRITE(0,*) 'enter path of file containing (clear-sky) profile data'
  READ(*,*) prof_filename
  WRITE(0,*) 'enter path of file containing hydrometeor profile data'
  READ(*,*) hydro_prof_filename
  !WRITE(0,*) 'enter number of profiles'
  !READ(*,*) nprof
  WRITE(0,*) 'enter number of profile pressure half-levels'
  READ(*,*) nlevels
  WRITE(0,*) 'enter number of channels to simulate per profile'
  READ(*,*) nchannels
  ALLOCATE(channel_list(nchannels))
  WRITE(0,*) 'enter space-separated channel list'
  READ(*,*,iostat=ios) channel_list(:)
  WRITE(0,*) 'enter number of threads to use'
  READ(*,*) nthreads

  nlayers = nlevels - 1
  WRITE(*,*) 'RUN/    nlayers=nlevel-1: ', nlayers

  !--------------------------------------------------------------! Vito
  READ(*,*) input_emissivity                                     ! Vito
  ALLOCATE (emissivity_val(nchannels))                           ! Vito
  emissivity_val = input_emissivity                              ! Vito
  READ(*,*) ifac                                                 ! Vito
  WRITE(*,*) 'CHECK SOLAR: ', opts % rt_all % solar
  READ(*,*) zenith                                               ! Vito
  !-------------------------------------------------------------------------! V


  ! --------------------------------------------------------------------------
  ! 1. Initialise RTTOV options structure
  ! --------------------------------------------------------------------------

  ! Set the interpolation method
  opts % interpolation % interp_mode = interp_rochon_wfn

  ! Include atmospheric refraction
  opts % rt_all % refraction = .FALSE. ! vito changed for consistency with rttov-scatt v13

  ! Enable hydrometeors:
  opts % scatt % hydrometeors = .TRUE.

  ! Select scattering solver for thermal radiation:
  opts % scatt % thermal_solver = thermal_solver_delta_edd

  ! Select cloud overlap parameterisation, and
  ! input a single hydrometeor fraction profile
  opts % cloud_overlap % overlap_param  = cloud_overlap_2col_weighted
  opts % cloud_overlap % per_hydro_frac = .FALSE.
  opts % cloud_overlap % col_threshold = 0.0010  ! Vito

  ! Set the relevant flag to .TRUE. when supplying a profile of the given trace
  ! gas (ensure the rtcoef file supports the gas):
  opts % rt_all % o3_data = .FALSE.


  !opts % rt_all % co2_data         = .FALSE.   ! Vito
  !opts % rt_all % n2o_data         = .FALSE.   ! Vito
  !opts % rt_all % ch4_data         = .FALSE.   ! Vito
  !opts % rt_all % co_data          = .FALSE.   ! Vito
  !opts % rt_all % so2_data         = .FALSE.   ! Vito
  !opts % clw_absorption % clw_data = .FALSE.   ! Vito

  ! Enable printing of warnings
  opts % config % verbose = .TRUE.

  !========== Interactive inputs == end ==============
  !===================================================


  ! --------------------------------------------------------------------------
  ! 2. Read coefficients
  ! --------------------------------------------------------------------------

  ! Read optical depth and hydrotable files together
  CALL rttov_read_coefs(errorstatus, coefs, opts,  &
                        file_coef = coef_filename, &
                        file_hydrotable = hydrotable_filename)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'fatal error reading coefficients'
    CALL rttov_exit(errorstatus)
  ENDIF
  print*, 'hydrotable_filename'

  ! Ensure input number of channels is not higher than number stored in coefficient file
  IF (nchannels > coefs % coef % fmv_nchn) THEN
    nchannels = coefs % coef % fmv_nchn
  ENDIF


  ! move opening file hre for nprof:
  OPEN(iup, file=TRIM(prof_filename), status='old', iostat=ios)
  WRITE(*,*) 'OPENING PROFILE CLEAR SKY'
  IF (ios /= 0) THEN
    WRITE(*,*) 'error opening profile file ios= ', ios
    CALL rttov_exit(errorstatus_fatal)
  ENDIF
  CALL rttov_skipcommentline(iup, errorstatus)
  READ(iup,*) nprof       

  print*, 'nprof: ', nprof
  ! --------------------------------------------------------------------------
  ! 3. Allocate RTTOV input and output structures
  ! --------------------------------------------------------------------------

  ! Determine the total number of radiances to simulate (nchanprof).
  ! In this example we simulate all specified channels for each profile, but
  ! in general one can simulate a different number of channels for each profile.

  nchanprof = nchannels * nprof

  ! Allocate structures for rttov_direct
  CALL rttov_alloc_direct( &
        errorstatus,       &
        alloc_flag,        &
        nchanprof,         &
        nprof,             &
        nlevels,           &
        nsurfaces,         &
        opts,              &
        coefs,             &
        chanprof,          &
        profiles,          &
        emis_refl,         &
        transmission,      &
        radiance,          &
        init = .TRUE._jplm)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'allocation error for rttov_direct structures'
    CALL rttov_exit(errorstatus)
  ENDIF

  WRITE(*,*) 'OK HERE? ALLOCATED DIRECT IN TEST'
  ! --------------------------------------------------------------------------
  ! 4. Build the list of profile/channel indices in chanprof
  ! --------------------------------------------------------------------------

  j = 0
  DO iprof = 1, nprof
    DO ichan = 1, nchannels
      j = j + 1
      chanprof(j)%prof = iprof
      chanprof(j)%chan = channel_list(ichan)
    ENDDO
  ENDDO


  ! --------------------------------------------------------------------------
  ! Optional: Ensure the options and coefficients are consistent
  ! --------------------------------------------------------------------------

  CALL rttov_user_check_options(errorstatus, opts, coefs, chanprof)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'error in rttov options'
    CALL rttov_exit(errorstatus)
  ENDIF


  ! --------------------------------------------------------------------------
  ! 5. Read profile data
  ! --------------------------------------------------------------------------

  !===============================================
  !========== Read profiles == start =============

  isurf = 1  ! Surface index

  !OPEN(iup, file=TRIM(prof_filename), status='old', iostat=ios)
  !WRITE(*,*) 'OPENING PROFILE CLEAR SKY'
  !IF (ios /= 0) THEN
  !  WRITE(*,*) 'error opening profile file ios= ', ios
  !  CALL rttov_exit(errorstatus_fatal)
  !ENDIF
  !CALL rttov_skipcommentline(iup, errorstatus)
  !READ(iup,*) nprof                                                           !

  ! Read gas units for profiles
  !READ(iup,*) profiles(1) % gas_units
  !profiles(:) % gas_units = profiles(1) % gas_units
  CALL rttov_skipcommentline(iup, errorstatus)

  ! Loop over all profiles and read data for each one
  DO iprof = 1, nprof

    print*, 'iprof loop: ', iprof

    READ(iup,*) profiles(iprof) % gas_units
    CALL rttov_skipcommentline(iup, errorstatus)

    ! Read pressure (hPa), temp (K), WV (gas units ppmv or kg/kg - as read above)
    READ(iup,*) profiles(iprof) % p_half(:)
    CALL rttov_skipcommentline(iup, errorstatus)
    READ(iup,*) profiles(iprof) % t(:)
    CALL rttov_skipcommentline(iup, errorstatus)
    READ(iup,*) profiles(iprof) % q(:)
    CALL rttov_skipcommentline(iup, errorstatus)

    ! Near-surface air variables
    READ(iup,*) profiles(iprof) % near_surface(isurf) % t2m, &
                profiles(iprof) % near_surface(isurf) % q2m, &
                profiles(iprof) % near_surface(isurf) % wind_u10m, &
                profiles(iprof) % near_surface(isurf) % wind_v10m
    CALL rttov_skipcommentline(iup, errorstatus)

    ! Skin variables
    READ(iup,*) profiles(iprof) % skin(isurf) % t,        &
    !            profiles(iprof) % skin(isurf) % salinity, &
                profiles(iprof) % skin(isurf) % fastem
    CALL rttov_skipcommentline(iup, errorstatus)

    ! Surface type and water type
    !READ(iup,*) profiles(iprof) % skin(isurf) % surftype, &
    !            profiles(iprof) % skin(isurf) % watertype
    !CALL rttov_skipcommentline(iup, errorstatus)

    profiles(iprof) % skin(isurf) % surftype = 1
    profiles(iprof) % skin(isurf) % watertype= 1

    ! Elevation, latitude and longitude
    READ(iup,*) profiles(iprof) % elevation!, &
    !            profiles(iprof) % latitude,  &
    !            profiles(iprof) % longitude
    CALL rttov_skipcommentline(iup, errorstatus)

    profiles(iprof) % zenangle = zenith 
    ! Satellite angles
    !READ(iup,*) profiles(iprof) % zenangle, &
    !            profiles(iprof) % azangle
    !CALL rttov_skipcommentline(iup, errorstatus)

  ENDDO
  CLOSE(iup)


  ! ---------------------------------------
  ! Specify the hydrometeor input profiles
  ! ---------------------------------------

  ! The MW hydrotable files contain optical parameters for 5 pre-defined particle types:
  ! rain, snow, graupel, cloud liquid water, cloud ice water.
  !
  ! The profiles(:)%hydro(:,:) array is of size (ntypes,nlayers) where ntypes is 5 for the
  ! default MW hydrotable files supplied with RTTOV. Each column corresponds to the particle
  ! types listed above, in that order.
  ! Hydrometeor concentrations are assigned to the column for the relevant particle. You can
  ! specify non-zero concentrations for any combination of available particle types in each
  ! layer.
  !
  ! Finally the vertical profile of total layer cloud fraction (for all hydrometeors in the
  ! layer) is specified in profiles(:)%hydro_frac(1,1:nlayers).

  OPEN(iup, file=TRIM(hydro_prof_filename), status='old', iostat=ios)
  WRITE(*,*) 'OPENING HYDRO PROFILE'
  IF (ios /= 0) THEN
    WRITE(*,*) 'error opening profile file ios= ', ios
    CALL rttov_exit(errorstatus_fatal)
  ENDIF
  CALL rttov_skipcommentline(iup, errorstatus)

  WRITE(*,*) 'PRINT OK'

  ! Read flag to indicate hydrometeor units (T => kg/kg; F => g/m^3)
  ! default is True (kg/kg)                                                      ! Vito
  !READ(iup,*) profiles(1)%mmr_hydro                                             ! Vito
  profiles(:)%mmr_hydro = .True.   !profiles(1)%mmr_hydro                        ! Vito
  !CALL rttov_skipcommentline(iup, errorstatus)

  WRITE(*,*) 'PRINT OK'
  ! Loop over all profiles and read data for each one
  DO iprof = 1, nprof

    ! Read cloud fraction and hydrometeor concentration profiles.
    DO ilay = 1, nlayers
      READ(iup,*) profiles(iprof) % hydro_frac(1,ilay),           &  ! cloud cover (0-1)
                  profiles(iprof) % hydro(hydro_index_clw,ilay),  &  ! cloud liquid water (kg/kg)
                  profiles(iprof) % hydro(hydro_index_ciw,ilay),  &  ! cloud ice water (kg/kg)
                  profiles(iprof) % hydro(hydro_index_rain,ilay), &  ! rain (kg/kg)
                  profiles(iprof) % hydro(hydro_index_snow,ilay), &  ! snow (kg/kg)
                  profiles(iprof) % hydro(hydro_index_graupel,ilay)  ! grau (kg/kg)
                  WRITE(*,*) 'ilay read:', ilay
      ! graupel (kg/kg) (not in data file)
      !profiles(iprof) % hydro(hydro_index_graupel,ilay) = 0._jprv
    ENDDO
    CALL rttov_skipcommentline(iup, errorstatus)

  ENDDO
  CLOSE(iup)

  !========== Read profiles == end =============
  !=============================================
  opts % scatt % hydrometeors = .TRUE.
  !opts % scatt % hydrometeors = .FALSE.

  ! --------------------------------------------------------------------------
  ! 6. Specify surface emissivity and reflectance
  ! --------------------------------------------------------------------------

  ! The emis_refl structure was initialised in the allocation call above: all
  ! numerical members are set to zero, and the calc_* members are set to true.

  ! In this example we leave RTTOV to provide all surface emissivities and
  ! reflectances.

  emis_refl(isurf) % emis_in(:)   = input_emissivity     ! Vito
  emis_refl(isurf) % calc_emis(:) = .FALSE.              ! Vito

  WRITE(*,*) 'PRINT OK'
  ! --------------------------------------------------------------------------
  ! 7. Call RTTOV forward model
  ! --------------------------------------------------------------------------
  IF (nthreads <= 1) THEN
    CALL rttov_direct(           &
            errorstatus,         &! out   error flag
            opts,                &! in    options structure
            coefs,               &! in    coefficients structure
            chanprof,            &! in    channel and profile index structure
            profiles,            &! in    profile array
            emis_refl,           &! inout emissivities/reflectances
            transmission,        &! inout computed transmittances
            radiance)             ! inout computed radiances
  ELSE
    CALL rttov_parallel_direct(  &
            errorstatus,         &! out   error flag
            opts,                &! in    options structure
            coefs,               &! in    coefficients structure
            chanprof,            &! in    channel and profile index structure
            profiles,            &! in    profile array
            emis_refl,           &! inout emissivities/reflectances
            transmission,        &! inout computed transmittances
            radiance,            &! inout computed radiances
            nthreads = nthreads)  ! in    number of threads to use
  ENDIF

  IF (errorstatus /= errorstatus_success) THEN
    WRITE (*,*) 'rttov_direct error'
    CALL rttov_exit(errorstatus)
  ENDIF

  !-----------------------------------  ! V
  ! Define the path and a string to assign to the name of the output                 ! V
  fmt =  '(I5)'                                                                      ! V
  ! Convert the angle to a string and put it to the filenames of the output          ! V
  WRITE(fout1,fmt) INT(zenith)                                                       ! V
  ! Trim `fout1` to remove leading and trailing whitespace
  fout1 = adjustl(trim(fout1))


  !=====================================================
  !============== Output results == start ==============

  ! Open output file where results are written
  !OPEN(ioout, file='output_'//NameOfRoutine//'.dat', status='unknown', form='formatted', iostat=ios)
  OPEN(ioout, file='output_'//trim(fout1)//'_'//trim(prof_filename), status='unknown', form='formatted', iostat=ios)
  IF (ios /= 0) THEN
    WRITE(*,*) 'error opening the output file ios= ', ios
    CALL rttov_exit(errorstatus_fatal)
  ENDIF

  !WRITE(ioout,*)' -----------------'
  !WRITE(ioout,*)' Instrument ', inst_name(coefs % coef % id_inst)
  !WRITE(ioout,*)' -----------------'
  !WRITE(ioout,*)' '
  CALL rttov_print_opts(opts)

  DO iprof = 1, nprof

    joff = (iprof - 1) * nchannels

    !WRITE(ioout,*)' '
    !WRITE(ioout,*)' Profile ', iprof

    !CALL rttov_print_profile(profiles(iprof), lu=ioout)

    !WRITE(ioout,777)'CHANNELS PROCESSED FOR SAT ', platform_name(coefs % coef % id_platform), coefs % coef % id_sat
    !WRITE(ioout,111) (chanprof(j) % chan, j = 1+joff, nchannels+joff)
    !WRITE(ioout,*)' '
    !WRITE(ioout,*)'CALCULATED BRIGHTNESS TEMPERATURES (K):'
    !WRITE(ioout,222) (radiance % bt(j), j = 1+joff, nchannels+joff)
    WRITE(ioout,222) (radiance % bt(j), j = 1+joff, nchannels+joff)
    !WRITE(ioout,*)' '
    !WRITE(ioout,*)'CALCULATED RADIANCES (mW/m2/sr/cm-1):'
    !WRITE(ioout,444) (radiance % total(j), j = 1+joff, nchannels+joff)
    !WRITE(ioout,*)' '
    !WRITE(ioout,*)'CALCULATED SURFACE EMISSIVITIES:'
    !WRITE(ioout,444) (emis_refl(isurf) % emis_out(j), j = 1+joff, nchannels+joff)
  ENDDO

  ! Close output file
  CLOSE(ioout, iostat=ios)
  IF (ios /= 0) THEN
    WRITE(*,*) 'error closing the output file ios= ', ios
    CALL rttov_exit(errorstatus_fatal)
  ENDIF


  !--------------------------------------------------------------------------------------------! Vito
  ! The following should be the same as Surface2SpaceTransmittance_CS.txt   == tau_total       ! V
  OPEN(21,   file='surface2space_transm_'//trim(fout1)//'_'//trim(prof_filename), status='unknown', form='formatted', iostat=ios)      ! V
  DO iprof = 1, nprof                                                                          ! V
    joff = (iprof-1_jpim) * nchannels                                                          ! V
    WRITE(21,4444) (transmission % tau_total(j), j = 1+joff, nchannels+joff)                   ! V
  ENDDO                                                                                        ! V
  CLOSE(21)                                                                      ! Vito 


  !============== Output results == end ==============
  !=====================================================


  ! --------------------------------------------------------------------------
  ! 8. Deallocate all RTTOV arrays and structures
  ! --------------------------------------------------------------------------
  DEALLOCATE (channel_list, stat=alloc_status)
  IF (alloc_status /= 0) THEN
    WRITE(*,*) 'mem dellocation error'
  ENDIF

  ! Deallocate structures for rttov_direct
  CALL rttov_alloc_direct( &
        errorstatus,       &
        dealloc_flag,      &
        nchanprof,         &
        nprof,             &
        nlevels,           &
        nsurfaces,         &
        opts,              &
        coefs,             &
        chanprof,          &
        profiles,          &
        emis_refl,         &
        transmission,      &
        radiance)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'deallocation error for rttov_direct structures'
    CALL rttov_exit(errorstatus)
  ENDIF

  CALL rttov_dealloc_coefs(errorstatus, coefs)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'coefs deallocation error'
  ENDIF


! Format definitions for output
111  FORMAT(1X,24I8)
222  FORMAT(1X,24F8.2)
444  FORMAT(1X,24F8.3)
777  FORMAT(/,A,A9,I3)
4446 FORMAT(1X,I3,24F8.4)      ! Vito
4444 FORMAT(1X,24F8.4)         ! Vito

END PROGRAM example_hydro_mw_file_fwd
