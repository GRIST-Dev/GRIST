!===================================================================================================
!
!  Created by LiXiaohan on 19/08/10, adopted from CAM5
!
!  Provide constituent distributions and properties to the radiation and 
!  cloud microphysics routines.
!  
!  The logic to control which constituents are used in the climate calculations
!  and which are used in diagnostic radiation calculations is contained in this module.
!
!===================================================================================================

 module grist_phys_prop

    use grist_handle_error,                 only: endrun
    use grist_constants,                    only: r8, i4
    use radconstants,                       only: nrh, nlwbands, nswbands, idx_sw_diag
    use grist_mpi
    use grist_wrap_nf

    implicit none
    private
    save

    public :: physprop_accum_unique_files ,  &! Make a list of the unique set of files that contain properties
              physprop_get_id             ,  &! Return ID used to access the property data from the input files
              physprop_init               ,  &! Initialization -- read the input datasets
              physprop_get                    ! Return data for specified ID

    integer, parameter, public :: ot_length = 32

    ! Data from one input dataset is stored in a structure of type(physprop_type).
    type :: physprop_type
       character(len=256) :: sourcefile ! Absolute pathname of data file.
       character(len=ot_length)  :: opticsmethod ! one of {hygro,nonhygro}
    
       ! for hygroscopic species of externally mixed aerosols
       real(r8), allocatable :: sw_hygro_ext(:,:)
       real(r8), allocatable :: sw_hygro_ssa(:,:)
       real(r8), allocatable :: sw_hygro_asm(:,:)
       real(r8), allocatable :: lw_hygro_abs(:,:)
    
       ! for nonhygroscopic species of externally mixed aerosols
       real(r8), allocatable :: sw_nonhygro_ext(:)
       real(r8), allocatable :: sw_nonhygro_ssa(:)
       real(r8), allocatable :: sw_nonhygro_asm(:)
       real(r8), allocatable :: sw_nonhygro_scat(:)
       real(r8), allocatable :: sw_nonhygro_ascat(:)
       real(r8), allocatable :: lw_abs(:)
    
       ! complex refractive index
       complex(r8), allocatable :: refindex_aer_sw(:)
       complex(r8), allocatable :: refindex_aer_lw(:)
    
       ! for radius-dependent mass-specific quantities
       real(r8), allocatable :: r_sw_ext(:,:)
       real(r8), allocatable :: r_sw_scat(:,:)
       real(r8), allocatable :: r_sw_ascat(:,:)
       real(r8), allocatable :: r_lw_abs(:,:)
       real(r8), allocatable :: mu(:)
    
       ! for modal optics
       real(r8), allocatable :: extpsw(:,:,:,:) ! specific extinction
       real(r8), allocatable :: abspsw(:,:,:,:) ! specific absorption
       real(r8), allocatable :: asmpsw(:,:,:,:) ! asymmetry factor
       real(r8), allocatable :: absplw(:,:,:,:) ! specific absorption
       real(r8), allocatable :: refrtabsw(:,:)  ! table of real refractive indices for aerosols visible
       real(r8), allocatable :: refitabsw(:,:)  ! table of imag refractive indices for aerosols visible
       real(r8), allocatable :: refrtablw(:,:)  ! table of real refractive indices for aerosols infrared
       real(r8), allocatable :: refitablw(:,:)  ! table of imag refractive indices for aerosols infrared
    
       ! microphysics parameters.
       character(len=32) :: aername ! for output of number concentration
       real(r8) :: density_aer      ! density of aerosol (kg/m3)
       real(r8) :: hygro_aer        ! hygroscopicity of aerosol
       real(r8) :: dryrad_aer       ! number mode radius (m) of aerosol size distribution
       real(r8) :: dispersion_aer   ! geometric standard deviation of aerosol size distribution
       real(r8) :: num_to_mass_aer  ! ratio of number concentration to mass concentration (#/kg)
                                    ! *** Is this actually (kg/#) ???
       ! mode parameters
       integer :: ncoef       ! number of Chebyshev coefficients
       integer :: prefr       ! dimension in table of real refractive indices
       integer :: prefi       ! dimension in table of imag refractive indices
       real(r8) :: sigmag     ! geometric standard deviation of the number distribution for aerosol mode
       real(r8) :: dgnum      ! geometric dry mean diameter of the number distribution for aerosol mode
       real(r8) :: dgnumlo    ! lower limit of dgnum
       real(r8) :: dgnumhi    ! upper limit of dgnum
       real(r8) :: rhcrystal  ! crystalization relative humidity for mode
       real(r8) :: rhdeliques ! deliquescence relative humidity for mode
    
       integer :: n_mu_samples ! LiXH adds

    endtype physprop_type

    ! This module stores data in an array of physprop_type structures.  The way this data
    ! is accessed outside the module is via a physprop ID, which is an index into the array.
    integer :: numphysprops = 0 ! an incremental total across ALL clim and diag constituents

    type (physprop_type), allocatable :: physprop(:)

    integer, parameter :: maxuniquefiles = 50
    character(len=256) :: uniquefilenames(maxuniquefiles)

 contains

    subroutine physprop_accum_unique_files(radname, type)

    ! Count number of aerosols in input radname array.  Aerosols are identified
    ! as strings with a ".nc" suffix.
    ! Construct a cumulative list of unique filenames containing physical property data.

    character(len=*), intent(in)  :: radname(:)
    character(len=1), intent(in)  :: type(:)

    integer :: ncnst, i
    character(len=*), parameter :: subname = 'physprop_accum_unique_files'
    !------------------------------------------------------------------------------------
  
    ncnst = ubound(radname, 1)

    do i = 1, ncnst

       ! check if radname is either a bulk aerosol or a mode
       if (type(i) == 'A' .or. type(i) == 'M') then

          ! check if this filename has been used by another aerosol.  If not
          ! then add it to the list of unique names.
          if (physprop_get_id(radname(i)) < 0) then
             numphysprops = numphysprops + 1
             if (numphysprops > maxuniquefiles) then
                print*, 'rank = ', mpi_rank() 
                print*, subname//': request for more than ',maxuniquefiles, ' values'
                call endrun(subname//': need to increase maxuniquefiles value')
             end if
             uniquefilenames(numphysprops) = trim(radname(i))
          endif

       endif
    enddo

    end subroutine physprop_accum_unique_files


    subroutine physprop_init()
    use grist_physics_iofile,       only: getfile
    ! Read properties from the aerosol data files.

    ! ***N.B.*** The calls to physprop_accum_unique_files must be made before calling
    !            this init routine.  physprop_accum_unique_files is responsible for building
    !            the list of files to be read here.

    ! Local variables
    integer            :: fileindex
    character(len=256) :: locfn ! path to actual file used
    character(len=32)  :: aername_str ! string read from netCDF file -- may contain trailing
                                      ! nulls which aren't dealt with by trim()
    
    integer :: ierr ! error codes from mpi

    !------------------------------------------------------------------------------------

    allocate(physprop(numphysprops))

    do fileindex = 1, numphysprops

       call getfile(uniquefilenames(fileindex), locfn, 0)
       physprop(fileindex)%sourcefile = locfn
 
       !--------------------LiXH modified the io manuscript--------------------> 
       !P.s. Only ~.nc can be used, open/close is merged in aerosol_optics_init 

       ! Open the physprop file
       !call cam_pio_openfile(nc_id, locfn, PIO_NOWRITE)
 
       call aerosol_optics_init(physprop(fileindex))
 
       ! Close the physprop file
       !call pio_closefile(nc_id)
       !<-------------------LiXH modified the io manuscript--------------------- 
 
    end do

    if(mpi_rank()==0)then
        print*, 'get all radiation constituent files successfully'
        print*, ''
    end if

    end subroutine physprop_init


    integer function physprop_get_id(filename)

    ! Look for filename in the global list of unique filenames (module data uniquefilenames).
    ! If found, return it's index in the list.  Otherwise return -1.

    character(len=*), intent(in) :: filename
    integer iphysprop

    physprop_get_id = -1
    do iphysprop = 1, numphysprops
      if(trim(uniquefilenames(iphysprop)) == trim(filename) ) then
        physprop_get_id = iphysprop
        return
      endif
    enddo

    end function physprop_get_id


! Purpose : Return requested properties for specified ID.
    subroutine physprop_get(id, sourcefile, opticstype,             &
    sw_hygro_ext, sw_hygro_ssa, sw_hygro_asm, lw_hygro_abs,         &
    sw_nonhygro_ext, sw_nonhygro_ssa, sw_nonhygro_asm,              &
    sw_nonhygro_scat, sw_nonhygro_ascat, lw_abs,                    &
    refindex_aer_sw, refindex_aer_lw,                               &
    r_sw_ext, r_sw_scat, r_sw_ascat, r_lw_abs, mu,                  &
    extpsw, abspsw, asmpsw, absplw, refrtabsw,                      &
    refitabsw, refrtablw, refitablw,                                &
    aername, density_aer, hygro_aer, dryrad_aer, dispersion_aer,    &
    num_to_mass_aer, ncoef, prefr, prefi, sigmag,                   &
    dgnum, dgnumlo, dgnumhi, rhcrystal, rhdeliques,                 &
    n_mu_samples)

    ! io
    integer,                            intent(in)  :: id
    character(len=256),       optional, intent(out) :: sourcefile ! Absolute pathname of data file.
    character(len=ot_length), optional, intent(out) :: opticstype
    real(r8),          optional, intent(out) :: sw_hygro_ext(:,:)
    real(r8),          optional, intent(out) :: sw_hygro_ssa(:,:) 
    real(r8),          optional, intent(out) :: sw_hygro_asm(:,:) 
    real(r8),          optional, intent(out) :: lw_hygro_abs(:,:)         
    real(r8),          optional, intent(out) :: sw_nonhygro_ext(:)
    real(r8),          optional, intent(out) :: sw_nonhygro_ssa(:)
    real(r8),          optional, intent(out) :: sw_nonhygro_asm(:)
    real(r8),          optional, intent(out) :: sw_nonhygro_scat(:)
    real(r8),          optional, intent(out) :: sw_nonhygro_ascat(:)
    real(r8),          optional, intent(out) :: lw_abs(:)         
    complex(r8),       optional, intent(out) :: refindex_aer_sw(:)
    complex(r8),       optional, intent(out) :: refindex_aer_lw(:)
    real(r8),          optional, intent(out) :: r_sw_ext(:,:)
    real(r8),          optional, intent(out) :: r_sw_scat(:,:)
    real(r8),          optional, intent(out) :: r_sw_ascat(:,:)
    real(r8),          optional, intent(out) :: r_lw_abs(:,:)
    real(r8),          optional, intent(out) :: mu(:)
    real(r8),          optional, intent(out) :: extpsw(:,:,:,:)
    real(r8),          optional, intent(out) :: abspsw(:,:,:,:)
    real(r8),          optional, intent(out) :: asmpsw(:,:,:,:)
    real(r8),          optional, intent(out) :: absplw(:,:,:,:)
    real(r8),          optional, intent(out) :: refrtabsw(:,:)
    real(r8),          optional, intent(out) :: refitabsw(:,:)
    real(r8),          optional, intent(out) :: refrtablw(:,:)
    real(r8),          optional, intent(out) :: refitablw(:,:)
    character(len=20), optional, intent(out) :: aername           
    real(r8),          optional, intent(out) :: density_aer       
    real(r8),          optional, intent(out) :: hygro_aer         
    real(r8),          optional, intent(out) :: dryrad_aer        
    real(r8),          optional, intent(out) :: dispersion_aer
    real(r8),          optional, intent(out) :: num_to_mass_aer
    integer,           optional, intent(out) :: ncoef
    integer,           optional, intent(out) :: prefr
    integer,           optional, intent(out) :: prefi
    real(r8),          optional, intent(out) :: sigmag
    real(r8),          optional, intent(out) :: dgnum
    real(r8),          optional, intent(out) :: dgnumlo
    real(r8),          optional, intent(out) :: dgnumhi
    real(r8),          optional, intent(out) :: rhcrystal
    real(r8),          optional, intent(out) :: rhdeliques
    ! LiXH adds for GRIST
    integer,           optional, intent(out) :: n_mu_samples

    ! local
    character(len=*), parameter :: subname = 'physprop_get'

    if (id <= 0 .or. id > numphysprops) then
       print*, subname//': illegal ID value: ', id,numphysprops,' My Rank:',mpi_rank()
       call endrun('physprop_get: ID out of range')
    end if

    if (present(sourcefile))        sourcefile        =  physprop(id)%sourcefile
    if (present(opticstype))        opticstype        =  physprop(id)%opticsmethod
    if (present(sw_hygro_ext))      sw_hygro_ext      =  physprop(id)%sw_hygro_ext
    if (present(sw_hygro_ssa))      sw_hygro_ssa      =  physprop(id)%sw_hygro_ssa
    if (present(sw_hygro_asm))      sw_hygro_asm      =  physprop(id)%sw_hygro_asm
    if (present(lw_hygro_abs))      lw_hygro_abs      =  physprop(id)%lw_hygro_abs
    if (present(sw_nonhygro_ext))   sw_nonhygro_ext   =  physprop(id)%sw_nonhygro_ext
    if (present(sw_nonhygro_ssa))   sw_nonhygro_ssa   =  physprop(id)%sw_nonhygro_ssa
    if (present(sw_nonhygro_asm))   sw_nonhygro_asm   =  physprop(id)%sw_nonhygro_asm
    if (present(sw_nonhygro_scat))  sw_nonhygro_scat  =  physprop(id)%sw_nonhygro_scat
    if (present(sw_nonhygro_ascat)) sw_nonhygro_ascat =  physprop(id)%sw_nonhygro_ascat
    if (present(lw_abs))            lw_abs            =  physprop(id)%lw_abs

    if (present(refindex_aer_sw))   refindex_aer_sw   =  physprop(id)%refindex_aer_sw
    if (present(refindex_aer_lw))   refindex_aer_lw   =  physprop(id)%refindex_aer_lw

    if (present(r_sw_ext))          r_sw_ext      =  physprop(id)%r_sw_ext
    if (present(r_sw_scat))         r_sw_scat     =  physprop(id)%r_sw_scat
    if (present(r_sw_ascat))        r_sw_ascat    =  physprop(id)%r_sw_ascat
    if (present(r_lw_abs))          r_lw_abs      =  physprop(id)%r_lw_abs
    if (present(mu))                mu            =  physprop(id)%mu

    if (present(extpsw))            extpsw        =  physprop(id)%extpsw
    if (present(abspsw))            abspsw        =  physprop(id)%abspsw
    if (present(asmpsw))            asmpsw        =  physprop(id)%asmpsw
    if (present(absplw))            absplw        =  physprop(id)%absplw
    if (present(refrtabsw))         refrtabsw     =  physprop(id)%refrtabsw
    if (present(refitabsw))         refitabsw     =  physprop(id)%refitabsw
    if (present(refrtablw))         refrtablw     =  physprop(id)%refrtablw
    if (present(refitablw))         refitablw     =  physprop(id)%refitablw

    if (present(aername))         aername         =  physprop(id)%aername
    if (present(density_aer))     density_aer     =  physprop(id)%density_aer
    if (present(hygro_aer))       hygro_aer       =  physprop(id)%hygro_aer
    if (present(dryrad_aer))      dryrad_aer      =  physprop(id)%dryrad_aer
    if (present(dispersion_aer))  dispersion_aer  =  physprop(id)%dispersion_aer
    if (present(num_to_mass_aer)) num_to_mass_aer =  physprop(id)%num_to_mass_aer

    if (present(ncoef))           ncoef           =  physprop(id)%ncoef
    if (present(prefr))           prefr           =  physprop(id)%prefr
    if (present(prefi))           prefi           =  physprop(id)%prefi
    if (present(sigmag))          sigmag          =  physprop(id)%sigmag
    if (present(dgnum))           dgnum           =  physprop(id)%dgnum
    if (present(dgnumlo))         dgnumlo         =  physprop(id)%dgnumlo
    if (present(dgnumhi))         dgnumhi         =  physprop(id)%dgnumhi
    if (present(rhcrystal))       rhcrystal       =  physprop(id)%rhcrystal
    if (present(rhdeliques))      rhdeliques      =  physprop(id)%rhdeliques

    ! LiXH adds for GRIST
    if (present(n_mu_samples))    n_mu_samples    =  physprop(id)%n_mu_samples
    end subroutine physprop_get


! Purpose: Determine the opticstype, then call the appropriate routine to read the data.
    subroutine aerosol_optics_init(phys_prop)
! io
    type(physprop_type), intent(inout) :: phys_prop  ! data after interp onto cam rh mesh
! local
    integer, parameter                 :: omode = 0  ! nf_nowrite   
    integer :: ierr
    integer :: ncid, varid, status
    integer :: opticslength_id, opticslength
    integer :: op_type_id

    call wrap_open(trim(phys_prop%sourcefile), omode, ncid)
    call wrap_inq_dimid(ncid, 'opticsmethod_len', opticslength_id)
    call wrap_inq_dimlen(ncid, opticslength_id, opticslength)
    if ( opticslength .gt. ot_length ) then
       if(mpi_rank()==0)print*,' optics type length in'//phys_prop%sourcefile//' excedes maximum length of 32' 
       call endrun('In aerosol_optics_init')
    endif
    call wrap_inq_varid(ncid, 'opticsmethod', op_type_id)
    call wrap_get_var_text(ncid, op_type_id,phys_prop%opticsmethod)

    select case (phys_prop%opticsmethod)
    case ('zero')
       call zero_optics_init(phys_prop, ncid)

    case ('hygro')
       call hygro_optics_init(phys_prop, ncid)

    case ('hygroscopic')
       call hygroscopic_optics_init(phys_prop, ncid) 

    case ('nonhygro')
       call nonhygro_optics_init(phys_prop, ncid)

    case ('insoluble')
       call insoluble_optics_init(phys_prop, ncid) 

    case ('volcanic_radius')
       call volcanic_radius_optics_init(phys_prop, ncid)

    case ('volcanic')
       call volcanic_optics_init(phys_prop, ncid)

    case ('modal')
       call modal_optics_init(phys_prop, ncid)

    case default
       if(mpi_rank()==0)print*,'aerosol_optics_init: unsupported optics type'//&
            trim(phys_prop%opticsmethod)//' in file '//phys_prop%sourcefile
       call endrun('In aerosol_optics_init')
    end select

    call wrap_close(ncid)

    end subroutine aerosol_optics_init


! Purpose : Read optics data of type 'hygro' and interpolate it to CAM's rh mesh.
    subroutine hygro_optics_init(phys_prop, nc_id)
    ! io
    type (physprop_type), intent(inout) :: phys_prop  ! data after interp onto cam rh mesh
    integer,              intent(in)    :: nc_id      ! indentifier for netcdf file
    ! local
    integer :: rh_idx_id, lw_band_id, sw_band_id
    integer :: kbnd, krh
    integer :: rh_id, sw_ext_id, sw_ssa_id, sw_asm_id, lw_ext_id
    integer :: nbnd, swbands
    ! temp data from hygroscopic file before interpolation onto cam-rh-mesh
    integer  :: nfilerh ! number of rh values in file
    real(r8), allocatable, dimension(:) :: frh
    real(r8), allocatable, dimension(:,:)  :: fsw_ext
    real(r8), allocatable, dimension(:,:)  :: fsw_ssa
    real(r8), allocatable, dimension(:,:)  :: fsw_asm

    real(r8) :: rh ! real rh value on cam rh mesh (indexvalue)
 
    allocate(phys_prop%sw_hygro_ext(nrh,nswbands))
    allocate(phys_prop%sw_hygro_ssa(nrh,nswbands))
    allocate(phys_prop%sw_hygro_asm(nrh,nswbands))
    allocate(phys_prop%lw_abs(nlwbands))

    call wrap_inq_dimid(nc_id, 'rh_idx', rh_idx_id)
    call wrap_inq_dimlen(nc_id, rh_idx_id, nfilerh)
    call wrap_inq_dimid(nc_id, 'lw_band', lw_band_id)
    call wrap_inq_dimid(nc_id, 'sw_band', sw_band_id)
    call wrap_inq_dimlen(nc_id, lw_band_id, nbnd)

    if (nbnd .ne. nlwbands) call endrun(phys_prop%sourcefile// &
         ' has the wrong number of lwbands')

    call wrap_inq_dimlen(nc_id, sw_band_id, swbands) 

    if(swbands .ne. nswbands) call endrun(phys_prop%sourcefile// &
          ' has the wrong number of sw bands')

    call wrap_inq_varid(nc_id, 'rh'    , rh_id)
    call wrap_inq_varid(nc_id, 'ext_sw', sw_ext_id)
    call wrap_inq_varid(nc_id, 'ssa_sw', sw_ssa_id)
    call wrap_inq_varid(nc_id, 'asm_sw', sw_asm_id)
    call wrap_inq_varid(nc_id, 'abs_lw', lw_ext_id)

    ! specific optical properties on file's rh mesh
    allocate(fsw_ext(nfilerh,nswbands))
    allocate(fsw_asm(nfilerh,nswbands))
    allocate(fsw_ssa(nfilerh,nswbands))
    allocate(frh(nfilerh))

    call wrap_get_var_realx(nc_id, rh_id, frh)
    call wrap_get_var_realx(nc_id, sw_ext_id, fsw_ext)
    call wrap_get_var_realx(nc_id, sw_ssa_id, fsw_ssa)
    call wrap_get_var_realx(nc_id, sw_asm_id, fsw_asm)
    call wrap_get_var_realx(nc_id, lw_ext_id, phys_prop%lw_abs)

    ! interpolate onto cam's rh mesh
    do kbnd = 1,nswbands
       do krh = 1, nrh
          rh = 1.0_r8 / nrh * (krh - 1)
          phys_prop%sw_hygro_ext(krh,kbnd) = &
             exp_interpol( frh, fsw_ext(:,kbnd) / fsw_ext(1,kbnd), rh ) &
             * fsw_ext(1, kbnd)
          phys_prop%sw_hygro_ssa(krh,kbnd) = &
             lin_interpol( frh, fsw_ssa(:,kbnd) / fsw_ssa(1,kbnd), rh ) &
             * fsw_ssa(1, kbnd)
          phys_prop%sw_hygro_asm(krh,kbnd) = &
             lin_interpol( frh, fsw_asm(:,kbnd) / fsw_asm(1,kbnd), rh ) &
             * fsw_asm(1, kbnd)
       enddo
    enddo

    deallocate (fsw_ext, fsw_asm, fsw_ssa, frh)

    ! read refractive index data if available
    call refindex_aer_init(phys_prop, nc_id)

    ! read bulk aero props
    call bulk_props_init(phys_prop, nc_id)

    end subroutine hygro_optics_init


    subroutine zero_optics_init(phys_prop, nc_id)
    ! io
    type (physprop_type), intent(inout) :: phys_prop  ! storage for file data
    integer,              intent(in) :: nc_id         ! indentifier for netcdf file
    ! local
    integer :: lw_band_id, sw_band_id
    integer :: sw_ext_id, sw_ssa_id, sw_asm_id, lw_ext_id
    integer :: swbands, nbnd
    integer :: ierr ! error flag
    !------------------------------------------------------------------------------------

    ! perhaps this doesn't even need allocated.
    allocate (phys_prop%sw_nonhygro_ext(nswbands))
    allocate (phys_prop%sw_nonhygro_ssa(nswbands))
    allocate (phys_prop%sw_nonhygro_asm(nswbands))
    allocate (phys_prop%lw_abs(nlwbands))

    phys_prop%sw_nonhygro_ext = 0._r8
    phys_prop%sw_nonhygro_ssa = 0._r8
    phys_prop%sw_nonhygro_asm = 0._r8
    phys_prop%lw_abs = 0._r8

    end subroutine zero_optics_init


! Purpose : Read optics data of type 'nonhygro'
    subroutine insoluble_optics_init(phys_prop, nc_id)
    ! io
    type (physprop_type), intent(inout) :: phys_prop  ! storage for file data
    integer,              intent(in)    :: nc_id      ! indentifier for netcdf file
    ! local
    integer :: lw_band_id, sw_band_id
    integer :: sw_ext_id, sw_ssa_id, sw_asm_id, lw_ext_id
    integer :: swbands, nbnd
    integer :: ierr ! error flag
    integer :: start(2), count(2)

    allocate (phys_prop%sw_nonhygro_ext(nswbands))
    allocate (phys_prop%sw_nonhygro_ssa(nswbands))
    allocate (phys_prop%sw_nonhygro_asm(nswbands))
    allocate (phys_prop%lw_abs(nlwbands))

    call wrap_inq_dimid(nc_id, 'lw_band', lw_band_id)
    call wrap_inq_dimid(nc_id, 'sw_band', sw_band_id)
    call wrap_inq_dimlen(nc_id, lw_band_id, nbnd)

    if (nbnd .ne. nlwbands) call endrun(phys_prop%sourcefile// &
         ' has the wrong number of lwbands')

    call wrap_inq_dimlen(nc_id, sw_band_id, swbands)

    if(swbands .ne. nswbands) call endrun(phys_prop%sourcefile// &
          ' has the wrong number of sw bands')

    call wrap_inq_varid(nc_id, 'ext_sw', sw_ext_id)
    call wrap_inq_varid(nc_id, 'ssa_sw', sw_ssa_id)
    call wrap_inq_varid(nc_id, 'asm_sw', sw_asm_id)
    call wrap_inq_varid(nc_id, 'abs_lw', lw_ext_id)

    start = 1
    count=(/1,swbands/)

    call wrap_get_vara_realx (nc_id, sw_ext_id, start, count, phys_prop%sw_nonhygro_ext)
    call wrap_get_vara_realx (nc_id, sw_ssa_id, start, count, phys_prop%sw_nonhygro_ssa)
    call wrap_get_vara_realx (nc_id, sw_asm_id, start, count, phys_prop%sw_nonhygro_asm)
    count = (/1,nbnd/)
    call wrap_get_vara_realx (nc_id, lw_ext_id, start, count, phys_prop%lw_abs)

    ! read refractive index data if available
    call refindex_aer_init(phys_prop, nc_id)

    ! read bulk aero props
    call bulk_props_init(phys_prop, nc_id)

    end subroutine insoluble_optics_init


! Purpose : Read optics data of type 'volcanic_radius'
    subroutine volcanic_radius_optics_init(phys_prop, nc_id)
    ! io
    type (physprop_type), intent(inout) :: phys_prop  ! storage for file data
    integer,              intent(in)    :: nc_id      ! indentifier for netcdf file
    ! local
    integer :: lw_band_id, sw_band_id, mu_id, mu_did
    integer :: sw_ext_id, sw_scat_id, sw_ascat_id, lw_abs_id
    integer :: swbands, nbnd, n_mu_samples
    integer :: ierr ! error flag

    call wrap_inq_dimid(nc_id, 'mu_samples', mu_did)
    call wrap_inq_dimlen(nc_id, mu_did, n_mu_samples)

    !------------LiXH adds for GRIST-------------->
    phys_prop%n_mu_samples = n_mu_samples
    !<-----------LiXH adds for GRIST---------------

    allocate (phys_prop%r_sw_ext(nswbands,n_mu_samples))
    allocate (phys_prop%r_sw_scat(nswbands,n_mu_samples))
    allocate (phys_prop%r_sw_ascat(nswbands,n_mu_samples))
    allocate (phys_prop%r_lw_abs(nlwbands,n_mu_samples))
    allocate (phys_prop%mu(n_mu_samples))

    call wrap_inq_dimid(nc_id, 'lw_band', lw_band_id)
    call wrap_inq_dimid(nc_id, 'sw_band', sw_band_id)
    call wrap_inq_dimlen(nc_id, lw_band_id, nbnd)

    if (nbnd .ne. nlwbands) call endrun(phys_prop%sourcefile// &
         ' has the wrong number of lwbands')

    call wrap_inq_dimlen(nc_id, sw_band_id, swbands) 

    if (swbands .ne. nswbands) call endrun(phys_prop%sourcefile// &
         ' has the wrong number of sw bands')

    call wrap_inq_varid(nc_id, 'bext_sw', sw_ext_id)
    call wrap_inq_varid(nc_id, 'bsca_sw', sw_scat_id)
    call wrap_inq_varid(nc_id, 'basc_sw', sw_ascat_id)
    call wrap_inq_varid(nc_id, 'babs_lw', lw_abs_id)
    call wrap_inq_varid(nc_id, 'mu_samples', mu_id)

    call wrap_get_var_realx(nc_id, sw_ext_id, phys_prop%r_sw_ext)
    call wrap_get_var_realx(nc_id, sw_scat_id, phys_prop%r_sw_scat)
    call wrap_get_var_realx(nc_id, sw_ascat_id, phys_prop%r_sw_ascat)
    call wrap_get_var_realx(nc_id, lw_abs_id, phys_prop%r_lw_abs)
    call wrap_get_var_realx(nc_id, mu_id, phys_prop%mu)

    ! read bulk aero props
    call bulk_props_init(phys_prop, nc_id)

    end subroutine volcanic_radius_optics_init


! Purpose : Read optics data of type 'volcanic'
    subroutine volcanic_optics_init(phys_prop, nc_id)
    ! io
    type (physprop_type), intent(inout) :: phys_prop  ! storage for file data
    integer,              intent(in)    :: nc_id      ! indentifier for netcdf file

    ! local
    integer :: lw_band_id, sw_band_id
    integer :: sw_ext_id, sw_scat_id, sw_ascat_id, lw_abs_id
    integer :: swbands, nbnd
    integer :: ierr ! error flag

    allocate (phys_prop%sw_nonhygro_ext(nswbands))
    allocate (phys_prop%sw_nonhygro_scat(nswbands))
    allocate (phys_prop%sw_nonhygro_ascat(nswbands))
    allocate (phys_prop%lw_abs(nlwbands))

    call wrap_inq_dimid(nc_id, 'lw_band', lw_band_id)
    call wrap_inq_dimid(nc_id, 'sw_band', sw_band_id)
    call wrap_inq_dimlen(nc_id, lw_band_id, nbnd)

    if (nbnd .ne. nlwbands) call endrun(phys_prop%sourcefile// &
         ' has the wrong number of lwbands')

    call wrap_inq_dimlen(nc_id, sw_band_id, swbands) 

    if (swbands .ne. nswbands) call endrun(phys_prop%sourcefile// &
         ' has the wrong number of sw bands')

    call wrap_inq_varid(nc_id, 'bext_sw', sw_ext_id)
    call wrap_inq_varid(nc_id, 'bsca_sw', sw_scat_id)
    call wrap_inq_varid(nc_id, 'basc_sw', sw_ascat_id)
    call wrap_inq_varid(nc_id, 'babs_lw', lw_abs_id)
 
    call wrap_get_var_realx(nc_id, sw_ext_id, phys_prop%sw_nonhygro_ext)
    call wrap_get_var_realx(nc_id, sw_scat_id, phys_prop%sw_nonhygro_scat)
    call wrap_get_var_realx(nc_id, sw_ascat_id, phys_prop%sw_nonhygro_ascat)
    call wrap_get_var_realx(nc_id, lw_abs_id, phys_prop%lw_abs)
 
    ! read bulk aero props
    call bulk_props_init(phys_prop, nc_id)

    end subroutine volcanic_optics_init


! Purpose : Read optics data of type 'hygroscopic' and interpolate it to CAM's rh mesh.
    subroutine hygroscopic_optics_init(phys_prop, nc_id)
    ! io
    type (physprop_type), intent(inout) :: phys_prop  ! data after interp onto cam rh mesh
    integer,              intent(in)    :: nc_id      ! indentifier for netcdf file
    ! local
    integer :: ierr ! error flag
    integer :: rh_idx_id, lw_band_id, sw_band_id
    integer :: kbnd, krh
    integer :: rh_id, sw_ext_id, sw_ssa_id, sw_asm_id, lw_ext_id
    integer :: nbnd, swbands
    ! temp data from hygroscopic file before interpolation onto cam-rh-mesh
    integer  :: nfilerh ! number of rh values in file
    real(r8), allocatable, dimension(:)    :: frh
    real(r8), allocatable, dimension(:,:)  :: fsw_ext
    real(r8), allocatable, dimension(:,:)  :: fsw_ssa
    real(r8), allocatable, dimension(:,:)  :: fsw_asm
    real(r8), allocatable, dimension(:,:)  :: flw_abs

    real(r8) :: rh ! real rh value on cam rh mesh (indexvalue)
    character(len=*), parameter :: sub = 'hygroscopic_optics_init'

    allocate(phys_prop%sw_hygro_ext(nrh,nswbands))
    allocate(phys_prop%sw_hygro_ssa(nrh,nswbands))
    allocate(phys_prop%sw_hygro_asm(nrh,nswbands))
    allocate(phys_prop%lw_hygro_abs(nrh,nlwbands))

    call wrap_inq_dimid(nc_id, 'rh_idx', rh_idx_id)
    call wrap_inq_dimlen(nc_id, rh_idx_id, nfilerh)
    call wrap_inq_dimid(nc_id, 'lw_band', lw_band_id)
    call wrap_inq_dimlen(nc_id, lw_band_id, nbnd)

    if (nbnd .ne. nlwbands) call endrun(phys_prop%sourcefile// &
         ' has the wrong number of lwbands')

    call wrap_inq_dimid(nc_id, 'sw_band', sw_band_id)
    call wrap_inq_dimlen(nc_id, sw_band_id, swbands) 

    if(swbands .ne. nswbands) call endrun(phys_prop%sourcefile// &
          ' has the wrong number of sw bands')

    call wrap_inq_varid(nc_id, 'rh'    , rh_id)
    call wrap_inq_varid(nc_id, 'ext_sw', sw_ext_id)
    call wrap_inq_varid(nc_id, 'ssa_sw', sw_ssa_id)
    call wrap_inq_varid(nc_id, 'asm_sw', sw_asm_id)
    call wrap_inq_varid(nc_id, 'abs_lw', lw_ext_id)

    ! specific optical properties on file's rh mesh
    allocate(fsw_ext(nfilerh,nswbands))
    allocate(fsw_asm(nfilerh,nswbands))
    allocate(fsw_ssa(nfilerh,nswbands))
    allocate(flw_abs(nfilerh,nlwbands))
    allocate(frh(nfilerh))

    call wrap_get_var_realx(nc_id, rh_id, frh)
    call wrap_get_var_realx(nc_id, sw_ext_id, fsw_ext)
    call wrap_get_var_realx(nc_id, sw_ssa_id, fsw_ssa)
    call wrap_get_var_realx(nc_id, sw_asm_id, fsw_asm)
    call wrap_get_var_realx(nc_id, lw_ext_id, flw_abs)

    ! interpolate onto cam's rh mesh
    do kbnd = 1,nswbands
       do krh = 1, nrh
          rh = 1.0_r8 / nrh * (krh - 1)
          phys_prop%sw_hygro_ext(krh,kbnd) = &
             exp_interpol( frh, fsw_ext(:,kbnd) / fsw_ext(1,kbnd), rh ) &
             * fsw_ext(1, kbnd)
          phys_prop%sw_hygro_ssa(krh,kbnd) = &
             lin_interpol( frh, fsw_ssa(:,kbnd) / fsw_ssa(1,kbnd), rh ) &
             * fsw_ssa(1, kbnd)
          phys_prop%sw_hygro_asm(krh,kbnd) = &
             lin_interpol( frh, fsw_asm(:,kbnd) / fsw_asm(1,kbnd), rh ) &
             * fsw_asm(1, kbnd)
       enddo
    enddo
    do kbnd = 1,nlwbands
       do krh = 1, nrh
          rh = 1.0_r8 / nrh * (krh - 1)
          phys_prop%lw_hygro_abs(krh,kbnd) = &
             exp_interpol( frh, flw_abs(:,kbnd) / flw_abs(1,kbnd), rh ) &
             * flw_abs(1, kbnd)
       enddo
    enddo

    deallocate (fsw_ext, fsw_asm, fsw_ssa, flw_abs, frh)

    ! read refractive index data if available
    call refindex_aer_init(phys_prop, nc_id)

    ! read bulk aero props
    call bulk_props_init(phys_prop, nc_id)

    end subroutine hygroscopic_optics_init


! Purpose : Read optics data of type 'nonhygro'
    subroutine nonhygro_optics_init(phys_prop, nc_id)
    ! io
    type (physprop_type), intent(inout) :: phys_prop  ! storage for file data
    integer,              intent(in)    :: nc_id      ! indentifier for netcdf file
    ! local
    integer :: lw_band_id, sw_band_id
    integer :: sw_ext_id, sw_ssa_id, sw_asm_id, lw_ext_id
    integer :: swbands, nbnd
    integer :: ierr ! error flag

    allocate (phys_prop%sw_nonhygro_ext(nswbands))
    allocate (phys_prop%sw_nonhygro_ssa(nswbands))
    allocate (phys_prop%sw_nonhygro_asm(nswbands))
    allocate (phys_prop%lw_abs(nlwbands))

    call wrap_inq_dimid(nc_id, 'lw_band', lw_band_id)
    call wrap_inq_dimid(nc_id, 'sw_band', sw_band_id)
    call wrap_inq_dimlen(nc_id, lw_band_id, nbnd)

    if (nbnd .ne. nlwbands) call endrun(phys_prop%sourcefile// &
         ' has the wrong number of lwbands')

    call wrap_inq_dimlen(nc_id, sw_band_id, swbands) 

    if (swbands .ne. nswbands) call endrun(phys_prop%sourcefile// &
         ' has the wrong number of sw bands')

    call wrap_inq_varid(nc_id, 'ext_sw', sw_ext_id)
    call wrap_inq_varid(nc_id, 'ssa_sw', sw_ssa_id)
    call wrap_inq_varid(nc_id, 'asm_sw', sw_asm_id)
    call wrap_inq_varid(nc_id, 'abs_lw', lw_ext_id)

    call wrap_get_var_realx(nc_id, sw_ext_id, phys_prop%sw_nonhygro_ext)
    call wrap_get_var_realx(nc_id, sw_ssa_id, phys_prop%sw_nonhygro_ssa)
    call wrap_get_var_realx(nc_id, sw_asm_id, phys_prop%sw_nonhygro_asm)
    call wrap_get_var_realx(nc_id, lw_ext_id, phys_prop%lw_abs)

    ! read refractive index data if available
    call refindex_aer_init(phys_prop, nc_id)

    ! read bulk aero props
    call bulk_props_init(phys_prop, nc_id)

    end subroutine nonhygro_optics_init


! Purpose: read refractive indices of aerosol
    subroutine refindex_aer_init(phys_prop, nc_id)

    include 'netcdf.inc'

    ! io
    type (physprop_type), intent(inout) :: phys_prop  ! storage for file data
    integer,              intent(in)    :: nc_id      ! indentifier for netcdf file
    ! local
    integer :: i
    integer :: istat1, istat2, istat3                 ! status flags
    integer :: vid_real, vid_im                       ! variable ids
    real(r8), pointer :: ref_real(:), ref_im(:)       ! tmp storage for components of complex index

    ! assume that the dimensions lw_band and sw_band have already been checked
    ! by the calling subroutine

    ! Check that the variables are present before allocating storage and reading.
    ! Since we're setting complex data values, both the real and imaginary parts must
    ! be present or neither will be read.

    istat1 = nf_inq_varid (nc_id, 'refindex_real_aer_sw', vid_real)
    istat2 = nf_inq_varid (nc_id, 'refindex_im_aer_sw',   vid_im)

    if(istat1 .eq. NF_NOERR .and. istat2 .eq. NF_NOERR)then
        allocate(ref_real(nswbands), ref_im(nswbands))

        call wrap_get_var_realx(nc_id, vid_real, ref_real)

        call wrap_get_var_realx(nc_id, vid_im, ref_im)

        ! successfully read refindex data -- set complex values in physprop object
        allocate(phys_prop%refindex_aer_sw(nswbands))
        do i = 1, nswbands
           phys_prop%refindex_aer_sw(i) = cmplx(ref_real(i), abs(ref_im(i)))
        end do

        deallocate(ref_real, ref_im)

    end if

    istat1 = nf_inq_varid (nc_id, 'refindex_real_aer_lw', vid_real)
    istat2 = nf_inq_varid (nc_id, 'refindex_im_aer_lw',   vid_im)

    if(istat1 .eq. NF_NOERR .and. istat2 .eq. NF_NOERR)then
        allocate(ref_real(nlwbands), ref_im(nlwbands))

        call wrap_get_var_realx(nc_id, vid_real, ref_real)

        call wrap_get_var_realx(nc_id, vid_im, ref_im)

        ! successfully read refindex data -- set complex values in physprop object
        allocate(phys_prop%refindex_aer_lw(nlwbands))
        do i = 1, nlwbands
           phys_prop%refindex_aer_lw(i) = cmplx(ref_real(i), abs(ref_im(i)))
        end do

        deallocate(ref_real, ref_im)

    end if

    end subroutine refindex_aer_init


! Purpose: read optics data for modal aerosols
    subroutine modal_optics_init(props, nc_id)
    ! io
    type (physprop_type), intent(inout) :: props      ! storage for file data
    integer,              intent(in)    :: nc_id      ! indentifier for netcdf file

    ! local
    integer :: ierr
    integer :: did
    integer :: ival
    integer :: vid
    real(r8):: tmp(1)
    real(r8), pointer :: rval(:,:,:,:,:) ! temp array used to eliminate a singleton dimension
    character(len=*), parameter :: subname = 'modal_optics_init'

    ! Check dimensions for number of lw and sw bands
    call wrap_inq_dimid(nc_id, 'lw_band', did)
    call wrap_inq_dimlen(nc_id, did, ival)
    if (ival .ne. nlwbands) call endrun(subname//':'//props%sourcefile// &
         ' has the wrong number of lwbands')

    call wrap_inq_dimid(nc_id, 'sw_band', did)
    call wrap_inq_dimlen(nc_id, did, ival) 
    if (ival .ne. nswbands) call endrun(subname//':'//props%sourcefile// &
         ' has the wrong number of sw bands')

    call wrap_inq_dimid(nc_id, 'coef_number', did)
    call wrap_inq_dimlen(nc_id, did, props%ncoef) 

    call wrap_inq_dimid(nc_id, 'refindex_real', did)
    call wrap_inq_dimlen(nc_id, did, props%prefr)

    call wrap_inq_dimid(nc_id, 'refindex_im', did)
    call wrap_inq_dimlen(nc_id, did, props%prefi)

    ! Allocate arrays
    allocate( props%extpsw(props%ncoef,props%prefr,props%prefi,nswbands), &
              props%abspsw(props%ncoef,props%prefr,props%prefi,nswbands), &
              props%asmpsw(props%ncoef,props%prefr,props%prefi,nswbands), &
              props%absplw(props%ncoef,props%prefr,props%prefi,nlwbands), &
              props%refrtabsw(props%prefr,nswbands), &
              props%refitabsw(props%prefi,nswbands), &
              props%refrtablw(props%prefr,nlwbands), &
              props%refitablw(props%prefi,nlwbands)  )

    ! allocate temp to remove the mode dimension from the sw variables
    allocate(rval(props%ncoef,props%prefr,props%prefi,1,nswbands))

    call wrap_inq_varid(nc_id, 'extpsw', vid)
    call wrap_get_var_realx(nc_id, vid, rval)
    props%extpsw = rval(:,:,:,1,:)
 
    call wrap_inq_varid(nc_id, 'abspsw', vid)
    call wrap_get_var_realx(nc_id, vid, rval)
    props%abspsw = rval(:,:,:,1,:)
    
    call wrap_inq_varid(nc_id, 'asmpsw', vid)
    call wrap_get_var_realx(nc_id, vid, rval)
    props%asmpsw = rval(:,:,:,1,:)
 
    deallocate(rval)

    ! allocate temp to remove the mode dimension from the lw variables
    allocate(rval(props%ncoef,props%prefr,props%prefi,1,nlwbands))

    call wrap_inq_varid(nc_id, 'absplw', vid)
    call wrap_get_var_realx(nc_id, vid, rval)
    props%absplw = rval(:,:,:,1,:)
 
    deallocate(rval)

    call wrap_inq_varid(nc_id, 'refindex_real_sw', vid)
    call wrap_get_var_realx(nc_id, vid, props%refrtabsw)
 
    call wrap_inq_varid(nc_id, 'refindex_im_sw', vid)
    call wrap_get_var_realx(nc_id, vid, props%refitabsw)
 
    call wrap_inq_varid(nc_id, 'refindex_real_lw', vid)
    call wrap_get_var_realx(nc_id, vid, props%refrtablw)
 
    call wrap_inq_varid(nc_id, 'refindex_im_lw', vid)
    call wrap_get_var_realx(nc_id, vid, props%refitablw)
 
    call wrap_inq_varid(nc_id, 'sigmag', vid)
    call wrap_get_var_realx(nc_id, vid, tmp)
    props%sigmag = tmp(1)
 
    call wrap_inq_varid(nc_id, 'dgnum', vid)
    call wrap_get_var_realx(nc_id, vid, tmp)
    props%dgnum = tmp(1)

    call wrap_inq_varid(nc_id, 'dgnumlo', vid)
    call wrap_get_var_realx(nc_id, vid, tmp)
    props%dgnumlo = tmp(1)

    call wrap_inq_varid(nc_id, 'dgnumhi', vid)
    call wrap_get_var_realx(nc_id, vid, tmp)
    props%dgnumhi = tmp(1)

    call wrap_inq_varid(nc_id, 'rhcrystal', vid)
    call wrap_get_var_realx(nc_id, vid, tmp)
    props%rhcrystal = tmp(1)

    call wrap_inq_varid(nc_id, 'rhdeliques', vid)
    call wrap_get_var_realx(nc_id, vid, tmp)
    props%rhdeliques = tmp(1)
 
    end subroutine modal_optics_init


! Purpose: read props for bulk aerosols
    subroutine bulk_props_init(physprop, nc_id)
    ! io
    type (physprop_type), intent(inout) :: physprop ! storage for file data
    integer,              intent(in)    :: nc_id    ! indentifier for netcdf file
    ! local
    integer :: ierr
    integer :: vid
    real(r8):: tmp(1)
    logical :: debug = .true.
    character(len=*), parameter :: subname = 'bulk_props_init'

    call wrap_inq_varid(nc_id, 'name', vid)
    call wrap_get_var_text(nc_id, vid, physprop%aername)

    ! use GLC function to remove trailing nulls and blanks.
    ! physprop%aername = aername_str(:GLC(aername_str))
    call wrap_inq_varid(nc_id, 'density', vid)
    call wrap_get_var_realx(nc_id, vid, tmp)
    physprop%density_aer = tmp(1)

    call wrap_inq_varid(nc_id, 'sigma_logr', vid)
    call wrap_get_var_realx(nc_id, vid, tmp)
    physprop%dispersion_aer = tmp(1)

    call wrap_inq_varid(nc_id, 'dryrad', vid)
    call wrap_get_var_realx(nc_id, vid, tmp)
    physprop%dryrad_aer = tmp(1)

    call wrap_inq_varid(nc_id, 'hygroscopicity', vid)
    call wrap_get_var_realx(nc_id, vid, tmp)
    physprop%hygro_aer = tmp(1)

    call wrap_inq_varid(nc_id, 'num_to_mass_ratio', vid)
    call wrap_get_var_realx(nc_id, vid, tmp)
    physprop%num_to_mass_aer = tmp(1)

    ! Output select data to log file
    if (debug .and. mpi_rank()==0) then
       if (trim(physprop%aername) == 'SULFATE') then
          print*, '_______ hygroscopic growth in visible band _______'
          call aer_optics_log_rh('SO4', physprop%sw_hygro_ext(:,idx_sw_diag), &
             physprop%sw_hygro_ssa(:,idx_sw_diag), physprop%sw_hygro_asm(:,idx_sw_diag))
       end if
       print*, subname//': finished for ', trim(physprop%aername)
    end if

    end subroutine bulk_props_init


! Purpose: interpolates f(x) to point y
!          assuming f(x) = f(x0) exp a(x - x0)
!          where a = ( ln f(x1) - ln f(x0) ) / (x1 - x0)
!          x0 <= x <= x1
!          assumes x is monotonically increasing
! Author: D. Fillmore
    function exp_interpol(x, f, y) result(g)

    implicit none

    real(r8), intent(in), dimension(:) :: x  ! grid points
    real(r8), intent(in), dimension(:) :: f  ! grid function values
    real(r8), intent(in) :: y                ! interpolation point
    real(r8) :: g                            ! interpolated function value

    integer :: k  ! interpolation point index
    integer :: n  ! length of x
    real(r8) :: a

    n = size(x)

    ! find k such that x(k) < y =< x(k+1)
    ! set k = 1 if y <= x(1)  and  k = n-1 if y > x(n)

    if (y <= x(1)) then
      k = 1
    else if (y >= x(n)) then
      k = n - 1
    else
      k = 1
      do while (y > x(k+1) .and. k < n)
        k = k + 1
      end do
    end if

    ! interpolate
    a = (  log( f(k+1) / f(k) )  ) / ( x(k+1) - x(k) )
    g = f(k) * exp( a * (y - x(k)) )
    return
    end function exp_interpol


! Purpose: interpolates f(x) to point y
!          assuming f(x) = f(x0) + a * (x - x0)
!          where a = ( f(x1) - f(x0) ) / (x1 - x0)
!          x0 <= x <= x1
!          assumes x is monotonically increasing
! Author: D. Fillmore
    function lin_interpol(x, f, y) result(g)

    implicit none

    real(r8), intent(in), dimension(:) :: x  ! grid points
    real(r8), intent(in), dimension(:) :: f  ! grid function values
    real(r8), intent(in) :: y                ! interpolation point
    real(r8) :: g                            ! interpolated function value

    integer :: k  ! interpolation point index
    integer :: n  ! length of x
    real(r8) :: a

    n = size(x)

    ! find k such that x(k) < y =< x(k+1)
    ! set k = 1 if y <= x(1)  and  k = n-1 if y > x(n)

    if (y <= x(1)) then
      k = 1
    else if (y >= x(n)) then
      k = n - 1
    else
      k = 1
      do while (y > x(k+1) .and. k < n)
        k = k + 1
      end do
    end if

    ! interpolate
    a = (  f(k+1) - f(k) ) / ( x(k+1) - x(k) )
    g = f(k) + a * (y - x(k))
    return
    end function lin_interpol


! Purpose : write out aerosol optical properties for a set of test rh values
!           to test hygroscopic growth interpolation
! Author: D. Fillmore
    subroutine aer_optics_log_rh(name, ext, ssa, asm)
    ! io
    character(len=*), intent(in) :: name
    real(r8), intent(in) :: ext(nrh)
    real(r8), intent(in) :: ssa(nrh)
    real(r8), intent(in) :: asm(nrh)
    ! local
    integer :: krh_test
    integer, parameter :: nrh_test = 36
    integer :: krh
    real(r8) :: rh
    real(r8) :: rh_test(nrh_test)
    real(r8) :: exti
    real(r8) :: ssai
    real(r8) :: asmi
    real(r8) :: wrh

    do krh_test = 1, nrh_test
       rh_test(krh_test) = sqrt(sqrt(sqrt(sqrt(((krh_test - 1.0_r8) / (nrh_test - 1))))))
    enddo
    write(*,'(2x, a)') name
    write(*,'(2x, a, 4x, a, 4x, a, 4x, a)') '   rh', 'ext (m^2 kg^-1)', '  ssa', '  asm'

    ! loop through test rh values
    do krh_test = 1, nrh_test
       ! find corresponding rh index
       rh = rh_test(krh_test)
       krh = min(floor( (rh) * nrh ) + 1, nrh - 1)
       wrh = (rh) *nrh - krh
       exti = ext(krh + 1) * (wrh + 1) - ext(krh) * wrh
       ssai = ssa(krh + 1) * (wrh + 1) - ssa(krh) * wrh
       asmi = asm(krh + 1) * (wrh + 1) - asm(krh) * wrh
       write(*, '(2x, f5.3, 4x, f13.3, 4x, f5.3, 4x, f5.3)') rh_test(krh_test), exti, ssai, asmi
    end do

    end subroutine aer_optics_log_rh


 end module grist_phys_prop
