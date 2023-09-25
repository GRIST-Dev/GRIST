 module grist_mps_nml_module

   use grist_constants,   only: i4, r8, day2sec
   use grist_time_manager,only: get_curr_calday
   use grist_mpi
   use grist_handle_error,only: endrun

   implicit none

   private

   public  :: set_mps_vars,   &

              CU_ml_model,  & ! CU machine learning
              mlFilePath

   character(1024)    :: CU_ml_model! CU machine learning
   character(1024)    :: mlFilePath

  contains

  subroutine set_mps_vars()

!================================================
! global vars have been defined in the header
!================================================

! local
  character(len=300) :: filename
  character(len=300) :: buffer
  integer (i4)       :: fileunit
  integer (i4)       :: leap_year_count, yr
  real(r8)           :: juld2, juld1
 
  namelist /mps_para/ CU_ml_model,  & ! CU machine learning
                      mlFilePath

    filename = "grist_mps.nml"

    fileunit = 1

    open  (fileunit, status='old',file=filename)
    read  (fileunit, nml=mps_para)
    close (fileunit)

  !  step_ra = 10800/idint(model_timestep)

    if(mpi_rank() .eq. 0)then
       print*,"**********************************************************"
       print*,"     The MPS is used following: ", trim(filename)
       print*,"**********************************************************"
    end if

    return
  end subroutine set_mps_vars

  end module grist_mps_nml_module
