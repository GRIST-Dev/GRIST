
!----------------------------------------------------------------------------
! Created on Aug 11 2023
! Author: Siyuan Chen
! Version 1.0
! Description: 
!----------------------------------------------------------------------------
module grist_tracer_module_vars_control_mixed

    use grist_constants,            only: r4 => ns
    use grist_domain_types,         only: global_domain
    use grist_nml_module,           only: nlev, nlevp
    use grist_data_utils,           only: scalar_r8_to_r4,  &
                                          scalar_r4_to_r8
! data
    use grist_hpe_constants,        only: eta_face_a_r8 => eta_face_a,  &
                                          eta_face_b_r8 => eta_face_b,  &
                                          eta_full_a_r8 => eta_full_a,  &
                                          eta_full_b_r8 => eta_full_b

    use grist_dycore_vars_module,           only: dycoreVarSurface, dycoreVarCellFace, dycoreVarCellFull
    use grist_tracer_transport_vars_module, only: tracerVarEdgeFull, tracerVarCellFace, tracerVarCellFull

    use grist_lib
    
    implicit none

    real(r4), allocatable   ::  eta_face_a(:)
    real(r4), allocatable   ::  eta_face_b(:)
    real(r4), allocatable   ::  eta_full_a(:)
    real(r4), allocatable   ::  eta_full_b(:)

    real(r4), allocatable   ::  edt_leng(:)
    real(r4), allocatable   ::  edp_leng(:)
    real(r4), allocatable   ::  edp_nr(:,:)
    real(r4), allocatable   ::  plg_nr(:,:)
    real(r4), allocatable   ::  plg_areag(:)
    real(r4), allocatable   ::  vtx_p(:,:)
    real(r4), allocatable   ::  edt_v_weight_prime_cell(:,:,:)

contains

    subroutine grist_tracer_module_vars_mixed_init(mesh)

        implicit none

        type(global_domain), intent(inout)  ::  mesh

        if(.not.allocated(eta_face_a_r8)) stop "not allocated eta_face_a_r8"
        if(.not.allocated(eta_face_b_r8)) stop "not allocated eta_face_b_r8"
        if(.not.allocated(eta_full_a_r8)) stop "not allocated eta_full_a_r8"
        if(.not.allocated(eta_full_b_r8)) stop "not allocated eta_full_b_r8"

        allocate(eta_face_a(nlevp)); eta_face_a = eta_face_a_r8
        allocate(eta_face_b(nlevp)); eta_face_b = eta_face_b_r8
        allocate(eta_full_a(nlev )); eta_full_a = eta_full_a_r8
        allocate(eta_full_b(nlev )); eta_full_b = eta_full_b_r8

        allocate(edt_leng(mesh%ne_full))
        allocate(edp_leng(mesh%ne_full))
        allocate(edp_nr(3, mesh%ne_full))
        allocate(plg_nr(mesh%maxvnb,mesh%nv_full))
        allocate(plg_areag(mesh%nv_full))
        allocate(vtx_p(3,mesh%nv_full))
        allocate(edt_v_weight_prime_cell(mesh%maxvnb,2,mesh%ne_full))

        edt_leng                = mesh%edt_leng
        edp_leng                = mesh%edp_leng
        edp_nr                  = mesh%edp_nr
        plg_nr                  = mesh%plg_nr
        plg_areag               = mesh%plg_areag
        vtx_p                   = mesh%vtx_p
        edt_v_weight_prime_cell = mesh%edt_v_weight_prime_cell
        
    end subroutine grist_tracer_module_vars_mixed_init

    subroutine grist_tracer_transport_vars_r8_to_r4

        implicit none

        call scalar_r8_to_r4(tracerVarEdgeFull%scalar_normal_velocity_avg_adv )
        call scalar_r8_to_r4(tracerVarEdgeFull%scalar_normal_mass_flux_avg_adv)
        call scalar_r8_to_r4(tracerVarCellFull%scalar_tracer_mass_n           )
        call scalar_r8_to_r4(tracerVarCellFull%scalar_tracer_mxrt_n           )

    end subroutine grist_tracer_transport_vars_r8_to_r4

    subroutine grist_tracer_transport_vars_r4_to_r8

        implicit none

        call scalar_r4_to_r8(tracerVarEdgeFull%scalar_normal_velocity_avg_adv )
        call scalar_r4_to_r8(tracerVarEdgeFull%scalar_normal_mass_flux_avg_adv)
        call scalar_r4_to_r8(tracerVarCellFull%scalar_tracer_mass_n           )
        call scalar_r4_to_r8(tracerVarCellFull%scalar_tracer_mxrt_n           )

    end subroutine grist_tracer_transport_vars_r4_to_r8

end module grist_tracer_module_vars_control_mixed