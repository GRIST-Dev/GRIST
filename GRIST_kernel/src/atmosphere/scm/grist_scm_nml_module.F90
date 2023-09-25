!======================================================
!
!  Created by LiXiaohan on 19/4/10.
!  NAMELIST module for GRIST_SCM
!======================================================

 module grist_scm_nml_module
    use grist_mpi
    use grist_constants,            only: i4, r8
    use grist_scm_comm_module,      only: start_ymd, start_tod,         &
                                          scm_lat, scm_lon,             &
                                          scm_file_name,                &
                                          scam_file_name,               &
                                          scm_test_name,                &
                                          scm_relaxation
    
    implicit none
    private
    public  :: grist_scm_read_nml

    contains

    subroutine grist_scm_read_nml

    character(len=300) :: filename
    integer (i4)       :: fileunit
    namelist /time_para/ start_ymd,              &
                         start_tod
    namelist /scm_para/  scm_lat,                &
                         scm_lon,                &
                         scm_file_name,          &
                         scam_file_name,         &
                         scm_test_name,          &
                         scm_relaxation

!================================================
!       set namelist defaults
!================================================
    start_ymd       = 0
    start_tod       = 0
    scm_lat         = 0._r8
    scm_lon         = 0._r8
    scm_file_name   = ''
    scam_file_name  = ''
    scm_test_name   = ''
    scm_relaxation  = .false.

!================================================
!       read namelist file
!================================================
    filename = "grist_scm.nml"
    fileunit = 11
    open  (fileunit, status='old',file=filename)
    read  (fileunit, nml=time_para)
    read  (fileunit, nml=scm_para)
    close (fileunit)
    
    if(mpi_rank() .eq. 0)then
        print*,"================================================"
        print*,"                                                "
        print*,"    Parameters of GRIST Single Column Model     "
        print*,"                                                "
        print*,"       Latitude:"         , scm_lat
        print*,"       Longitude:"        , scm_lon
        print*,"       Do Relaxation:"    , scm_relaxation
        print*,"================================================"
     end if
     return

     end subroutine grist_scm_read_nml

 end module grist_scm_nml_module
