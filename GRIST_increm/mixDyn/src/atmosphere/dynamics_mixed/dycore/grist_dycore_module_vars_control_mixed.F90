module grist_dycore_module_vars_control_mixed

    use grist_constants,        only: r4 => ns
    use grist_domain_types,     only: global_domain
    use grist_nml_module,       only: nlev, nlevp, www_damping_coef_r8 => www_damping_coef, zd_r8 => zd
! data
    use grist_hpe_constants,    only: eta_full_r8   => eta_full,    &
                                      eta_face_r8   => eta_face,    &
                                      eta_full_a_r8 => eta_full_a,  &
                                      eta_full_b_r8 => eta_full_b,  &
                                      eta_face_a_r8 => eta_face_a,  &
                                      eta_face_b_r8 => eta_face_b,  &
                                      deta_full_r8  => deta_full,   &
                                      deta_face_r8  => deta_face

    use grist_dycore_diffusion_module,  only: perot_weight_at_pc_r8 => perot_weight_at_pc
                                      
    implicit none

    real(r4)                    ::  www_damping_coef
    real(r4)                    ::  zd

    real(r4), allocatable       ::  eta_full  (:)
    real(r4), allocatable       ::  eta_face  (:)
    real(r4), allocatable       ::  eta_full_a(:)
    real(r4), allocatable       ::  eta_full_b(:)
    real(r4), allocatable       ::  eta_face_a(:)
    real(r4), allocatable       ::  eta_face_b(:)
    real(r4), allocatable       ::  deta_full (:)
    real(r4), allocatable       ::  deta_face (:)

    real(r4), allocatable       ::  edt_leng(:)
    real(r4), allocatable       ::  edp_leng(:)
    real(r4), allocatable       ::  edp_nr(:,:)
    real(r4), allocatable       ::  vtx_p(:,:)
    real(r4), allocatable       ::  plg_nr(:)
    real(r4), allocatable       ::  plg_areag(:)
    real(r4), allocatable       ::  tri_areag(:)
    real(r4), allocatable       ::  edp_trsk_edpl_on_edge(:,:)
    !real(r4), allocatable       ::  edp_trsk_on_edge(:,:)
    !real(r4), allocatable       ::  edp_edpl_on_edge(:,:)
    real(r4), allocatable       ::  edt_v_weight_prime_cell(:,:,:)
    real(r4), allocatable       ::  tri_c_lat(:)
    real(r4), allocatable       ::  tri_kite_area(:,:)

    real(r4), allocatable       ::  perot_weight_at_pc(:,:,:)

contains

    subroutine grist_dycore_module_vars_mixed_init(mesh, phase)

        implicit none

        type(global_domain), intent(inout)  ::  mesh
        integer,             intent(in)     ::  phase

        if(phase.eq.1)then

            www_damping_coef = www_damping_coef_r8
            zd               = zd_r8

            if(.not.allocated(eta_full_r8  )) stop "not allocated eta_full_r8  "
            if(.not.allocated(eta_face_r8  )) stop "not allocated eta_face_r8  "
            if(.not.allocated(eta_full_a_r8)) stop "not allocated eta_full_a_r8"
            if(.not.allocated(eta_full_b_r8)) stop "not allocated eta_full_b_r8"
            if(.not.allocated(eta_face_a_r8)) stop "not allocated eta_face_a_r8"
            if(.not.allocated(eta_face_b_r8)) stop "not allocated eta_face_b_r8"
            if(.not.allocated(deta_full_r8 )) stop "not allocated deta_full_r8 "
            if(.not.allocated(deta_face_r8 )) stop "not allocated deta_face_r8 "

            allocate(eta_full  (nlev) ); eta_full   = eta_full_r8
            allocate(eta_face  (nlevp)); eta_face   = eta_face_r8
            allocate(eta_full_a(nlev) ); eta_full_a = eta_full_a_r8
            allocate(eta_full_b(nlev) ); eta_full_b = eta_full_b_r8
            allocate(eta_face_a(nlevp)); eta_face_a = eta_face_a_r8
            allocate(eta_face_b(nlevp)); eta_face_b = eta_face_b_r8
            allocate(deta_full(nlev)  ); deta_full  = deta_full_r8
            allocate(deta_face(nlevp) ); deta_face  = deta_face_r8

            allocate(edt_leng(mesh%ne_full))
            allocate(edp_leng(mesh%ne_full))
            allocate(edp_nr(3, mesh%ne_full))
            allocate(vtx_p (3, mesh%nv_full))
            allocate(plg_areag(mesh%nv_full))
            allocate(tri_areag(mesh%nt_full))
            allocate(edp_trsk_edpl_on_edge(16,mesh%ne_full))
            !allocate(edp_trsk_on_edge(16,mesh%ne_full))
            !allocate(edp_edpl_on_edge(16,mesh%ne_full))
            allocate(edt_v_weight_prime_cell(mesh%maxvnb,2,mesh%ne_full))
            allocate(tri_c_lat(mesh%nt_full))
            allocate(tri_kite_area(3,mesh%nt_full))
            
            edt_leng                = mesh%edt_leng
            edp_leng                = mesh%edp_leng
            edp_nr                  = mesh%edp_nr
            vtx_p                   = mesh%vtx_p
            plg_areag               = mesh%plg_areag
            tri_areag               = mesh%tri_areag
            edp_trsk_edpl_on_edge   = mesh%edp_trsk_edpl_on_edge
            !edp_trsk_on_edge        = mesh%edp_trsk_on_edge
            !edp_edpl_on_edge        = mesh%edp_edpl_on_edge
            edt_v_weight_prime_cell = mesh%edt_v_weight_prime_cell
            tri_c_lat               = mesh%tri_c_lat
            tri_kite_area           = mesh%tri_kite_area

        else

            ! allocate(perot_weight_at_pc(3,8,mesh%nv_halo(2)))

            ! perot_weight_at_pc      = perot_weight_at_pc_r8

        end if
        
    end subroutine grist_dycore_module_vars_mixed_init

end module grist_dycore_module_vars_control_mixed
