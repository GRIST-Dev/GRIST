!======================================================
!
!  Created by LiXiaohan on 19/5/20.
!  Adopted from CAM_EUL
!  
!  Purpose:
!  to handle Semi-Lagrangian transport of tracer 
!  in the context of Eulerian Spectral dynamics.
!======================================================

 module grist_scm_scanslt

    use grist_constants,        only: i4, r8
    use grist_nml_module,       only: nlev, nlevp, ntracer
    use grist_handle_error,     only: endrun
    use grist_mpi

    implicit none

    private
 
    public :: scanslt_run

    integer, parameter  :: pmap = 20000             ! max dimension of evenly spaced vert. 
                                                    ! grid used by SLT code to map the departure pts into true
                                                    ! model levels.
    logical             :: plimdr = .true.

    real(r8), allocatable :: lbasdz(:,:,:)          ! vert (full levels) deriv wghts
    real(r8), allocatable :: lbassd(:,:,:)          ! vert (half levels) deriv wghts
    real(r8), allocatable :: detai(:)               ! intervals between vert half levs.
    integer  :: kdpmpf(pmap)                        ! artificial full vert grid indices
    integer  :: kdpmph(pmap)                        ! artificial half vert grid indices

 contains
   

! Purpose: Driving routine for semi-lagrangian transport.
    subroutine scanslt_run(nstep, qminus, dt,     etadot, &
                           detam, etamid, etaint, qfcst )
! io
    integer,  intent(in) :: nstep
    real(r8), intent(in) :: dt
    real(r8), intent(in) :: etadot(nlevp)
    real(r8), intent(in) :: etamid(nlev)
    real(r8), intent(in) :: etaint(nlevp)
    real(r8), intent(in) :: qminus(ntracer, nlev)
    real(r8), intent(inout) :: detam(nlev)     ! delta eta at levels
    real(r8), intent(out)   :: qfcst(ntracer, nlev)
! local
    integer :: k, iter
    real(r8) :: sigmp(nlev)

    integer, parameter :: itermx=4
    integer, parameter :: itermn=1

    allocate(lbasdz(4,2,nlev));  lbasdz = 0._r8
    allocate(lbassd(4,2,nlevp)); lbassd = 0._r8
    allocate(detai(nlevp));      detai  = 0._r8

    call grdini(pmap    ,etamid  ,etaint  ,             &
               lbasdz   ,lbassd  ,detam   ,detai    ,   &
               kdpmpf   ,kdpmph  )

    do k=1,nlev
        sigmp(k) = etamid(k)
    end do

    if (nstep .le. 1) then
       iter = itermx
    else
       iter = itermn
    end if
   
    call sltb1(pmap    ,dt      ,iter   ,etadot , &
               etamid  ,etaint  ,detam  ,detai  , &
               lbasdz  ,lbassd  ,kdpmpf ,kdpmph , &
               sigmp   ,qfcst   ,qminus )

    deallocate(lbasdz)
    deallocate(lbassd)
    deallocate(detai)

    end subroutine scanslt_run


! Purpose: 
! Initialize model and extended grid parameters
! Initialize weights for Lagrange cubic derivative estimates
! Initialize weights for Lagrange cubic interpolant
    subroutine grdini(pmap    ,etamid  ,etaint  ,          &
                      lbasdz  ,lbassd  ,detam   ,detai   , &
                      kdpmpf  ,kdpmph     )

! io
    integer,  intent(in) :: pmap              ! dimension of artificial vert. grid
    real(r8), intent(in) :: etamid(nlev)      ! full-level model vertical grid
    real(r8), intent(in) :: etaint(nlevp)     ! half-level model vertical grid

    real(r8), intent(out) :: lbasdz(4,2,nlev) ! vertical (full levels) deriv weights
    real(r8), intent(out) :: lbassd(4,2,nlevp)! vertical (half levels) deriv weights
    real(r8), intent(out) :: detam (nlev)     ! intervals between vertical full levs.
    real(r8), intent(out) :: detai (nlevp)    ! intervals between vertical half levs.

    integer,  intent(out) :: kdpmpf(pmap)     ! artificial full vertical grid indices
    integer,  intent(out) :: kdpmph(pmap)     ! artificial half vertical grid indices
! local
    integer k                 ! index
    real(r8) etamln(nlev)         ! log(etamid)
    real(r8) etailn(nlevp)        ! log(etaint)
    real(r8) detamln(nlev)        ! dlog(etamid)
    real(r8) detailn(nlevp)       ! dlog(etaint)

    detai(:)=0._r8
    kdpmpf(:)=0._r8
    kdpmph(:)=0._r8
    call basdz(nlev    ,etamid  ,lbasdz  )
    call basdz(nlevp   ,etaint  ,lbassd  )

! Compute interval lengths in vertical grids.
    do k = 1,nlev
       etamln(k) = log(etamid(k))
    end do
    do k = 1,nlevp
       etailn(k) = log(etaint(k))
    end do
    do k = 1,nlev-1
       detam  (k) = etamid(k+1) - etamid(k)
       detamln(k) = etamln(k+1) - etamln(k)
    end do
    do k = 1,nlev
       detai  (k) = etaint(k+1) - etaint(k)
       detailn(k) = etailn(k+1) - etailn(k)
    end do
!
! Build artificial evenly spaced vertical grid for use in determining
! vertical position of departure point.
! Build one grid for full model levels and one for half levels.
!
    call vrtmap(nlev    ,pmap    ,etamln  ,detamln ,kdpmpf  )
    call vrtmap(nlevp   ,pmap    ,etailn  ,detailn ,kdpmph  )

    return
    end subroutine grdini


! Purpose: 
! Map indices of an artificial evenly spaced (in log) vertical grid to
! the indices of the log of the model vertical grid.  The resultant
! array of mapped indices will be used by "kdpfnd" to find the vertical
! location of any departure point relative to the model grid.
    subroutine vrtmap (pkdim   ,pmap    ,sigln   ,dsigln  ,kdpmap  )
! io
    integer, intent(in) :: pkdim             ! dimension of "sigln" and "dsigln"
    integer, intent(in) :: pmap              ! dimension of "kdpmap"
    real(r8), intent(in) :: sigln (pkdim)    ! model levels (log(eta))
    real(r8), intent(in) :: dsigln(pkdim)    ! intervals between model levels (log)
    integer, intent(out) :: kdpmap(pmap)     ! array of mapped indices
! local
    integer imin              ! |
    integer k                 ! |-- indices
    integer kk                ! |
    integer newmap            ! estimated value of "pmap"

    real(r8) del              ! artificial grid interval
    real(r8) dp               ! artificial departure point
    real(r8) eps              ! epsilon factor

    eps = 1.e-05_r8
    del = ( sigln(pkdim) - sigln(1) )/real(pmap,r8)
    imin = ismin( pkdim-1,dsigln(:pkdim-1), 1 )
    if (del + eps  >=  dsigln(imin)) then
       newmap = ( sigln(pkdim) - sigln(1) )/dsigln(imin) + 1
       print*,' VRTMAP:  Not enough artificial grid intervals.'
       print*,' Currently, "pmap" is set to ',pmap
       print*,' Reset parameter "pmap" to at least ',newmap
       print*,'rank=',mpi_rank()
       call endrun("in vrtmap of grist_scm_scanslt.f90")
    end if

    kdpmap(1) = 1
    do kk = 2,pmap
       dp = sigln(1) + real(kk-1,r8)*del
       do k = 1,pkdim-1
          if(dp > sigln(k) + eps) then
             kdpmap(kk) = k
          end if
       end do
    end do

    return
    end subroutine vrtmap




! Purpose: 
! Compute weights for the calculation of derivative estimates at two
! center points of the four point stencil for each interval in the
! unequally spaced vertical grid (as defined by the array sig).
! Estimates are from differentiating a Lagrange cubic polynomial
! through the four point stencil.
    subroutine basdz(pkdim   ,sig     ,lbasdz  )
! io
    integer , intent(in) :: pkdim      ! vertical dimension
    real(r8), intent(in) :: sig(pkdim) ! sigma levels (actually a generic vert. coord)
    real(r8), intent(out):: lbasdz(4,2,pkdim) ! vertical interpolation weights
! local
    integer kk                ! index

    do kk = 2,pkdim-2
       call lcdbas( sig(kk-1:kk+2), lbasdz(1:4,1,kk), lbasdz(1:4,2,kk) )
    end do

    end subroutine basdz


! Purpose: 
! Calculate weights used to evaluate derivative estimates at the
! inner grid points of a four point stencil based on Lagrange
! cubic polynomial through four unequally spaced points.
    subroutine lcdbas(grd     ,dbas2   ,dbas3   )
! io
    real(r8), intent(in) :: grd(4)    ! grid stencil
    real(r8), intent(out):: dbas2(4)  ! derivatives at grid point 2.
    real(r8), intent(out):: dbas3(4)  ! derivatives at grid point 3.

!  grd    Coordinate values of four points in stencil.
!  dbas2  Derivatives of the four basis functions at grid point 2.
!  dbas3  Derivatives of the four basis functions at grid point 3.

! local
    real(r8) x1                   !  |
    real(r8) x2                   !  |- grid values
    real(r8) x3                   !  |
    real(r8) x4                   !  |
    real(r8) x1mx2                !  |
    real(r8) x1mx3                !  |
    real(r8) x1mx4                !  |- differences of grid values
    real(r8) x2mx3                !  |
    real(r8) x2mx4                !  |
    real(r8) x3mx4                !  |

    x1 = grd(1)
    x2 = grd(2)
    x3 = grd(3)
    x4 = grd(4)
    x1mx2 = x1 - x2
    x1mx3 = x1 - x3
    x1mx4 = x1 - x4
    x2mx3 = x2 - x3
    x2mx4 = x2 - x4
    x3mx4 = x3 - x4

    dbas2(1) =   x2mx3 * x2mx4 / ( x1mx2 * x1mx3 * x1mx4 )
    dbas2(2) =   -1._r8/x1mx2 + 1._r8/x2mx3 + 1._r8/x2mx4
    dbas2(3) = - x1mx2 * x2mx4 / ( x1mx3 * x2mx3 * x3mx4 )
    dbas2(4) =   x1mx2 * x2mx3 / ( x1mx4 * x2mx4 * x3mx4 )

    dbas3(1) = - x2mx3 * x3mx4 / ( x1mx2 * x1mx3 * x1mx4 )
    dbas3(2) =   x1mx3 * x3mx4 / ( x1mx2 * x2mx3 * x2mx4 )
    dbas3(3) =   -1._r8/x1mx3 - 1._r8/x2mx3 + 1._r8/x3mx4
    dbas3(4) = - x1mx3 * x2mx3 / ( x1mx4 * x2mx4 * x3mx4 )

    return
    end subroutine lcdbas



    subroutine sltb1(pmap    ,dt      ,iterdp  ,wb      , &
                     sig     ,sigh    ,dsig    ,dsigh   , &
                     lbasdz  ,lbassd  ,kdpmpf  ,kdpmph  , &
                     sigmp   ,fbout   ,qminus  )

! io
    integer , intent(in) :: pmap                ! artificial vert grid dim.
    real(r8), intent(in) :: dt                  ! time step (seconds)
    integer , intent(in) :: iterdp              ! iteration count
    real(r8), intent(in) :: wb(nlevp)           ! eta-dot
    real(r8), intent(in) :: sig  (nlev)         ! vertical full levels
    real(r8), intent(in) :: sigh (nlevp)        ! vertical half levels
    real(r8), intent(in) :: dsig (nlev)         ! inc. between full levs
    real(r8), intent(in) :: dsigh(nlevp)        ! inc. between half levs
    real(r8), intent(in) :: lbasdz(4,2,nlev)    ! vert full level deriv wts
    real(r8), intent(in) :: lbassd(4,2,nlevp)   ! vert half level deriv wts
    integer , intent(in) :: kdpmpf(pmap)        ! artificial vert grid index
    integer , intent(in) :: kdpmph(pmap)        ! artificial vert grid index
    real(r8), intent(in) :: qminus(ntracer,nlev)! constituents
    real(r8), intent(inout) :: sigmp(nlev)      ! vert coord of mid-point
    real(r8), intent(out) :: fbout(ntracer,nlev)! advected constituents
! local
    integer m, k
    integer kdp(nlev)               ! vertical dep point index
    real(r8) fhr(ntracer,nlev)      ! horizontal interpolants
    real(r8) sigdp(nlev)            ! vertical   departure pt. coord.
    real(r8) fhst(ntracer,nlev)     ! derivative at top of interval
    real(r8) fhsb(ntracer,nlev)     ! derivative at bot of interval
    real(r8) wst(nlevp)             ! w derivative at top of interval
    real(r8) wsb(nlevp)             ! w derivative at bot of interval

    do m = 1, ntracer
       do k = 1, nlev
          fhr(m,k) = qminus(m,k)
       end do
    end do

    ! Compute vertical derivatives of vertical wind
    call cubzdr(nlevp, wb, lbassd, wst, wsb )

    ! Compute departure points and corresponding indices.
    call vrtdep(pmap    ,dt      ,iterdp  ,wb      ,wst     , &
                wsb     ,sig     ,sigh    ,dsigh   ,kdpmpf  , &
                kdpmph  ,sigmp   ,sigdp   ,kdp )

    ! Vertical derivatives of scalar fields.
    ! Loop over constituents.
    do m = 1, ntracer
        call cubzdr(nlev   ,fhr(m,:)  ,lbasdz  ,fhst(m,:)  , fhsb(m,:) )
    end do

    if( plimdr )then
        call limdz(fhr   ,dsig    ,fhst    ,fhsb  )
    end if


    ! Vertical interpolation of scalar fields.
    call herzin(nlev    ,ntracer   ,fhr     ,fhst    ,fhsb    , &
                sig     ,dsig      ,sigdp   ,kdp     ,fbout   )

    end subroutine sltb1


! Purpose: 
! Vertical derivative estimates for a vertical slice using Lagrangian
! cubic formulas.  
!
! Method: 
! Derivatives are set to zero at the top and bottom.
! At the "inner nodes" of the top and bottom intervals, a "one sided"
! estimate is used.
    subroutine cubzdr(pkdim   ,f    ,lbasdz  ,dfz1   , dfz2)
! io
    integer,  intent(in) :: pkdim             ! vertical     dimension
    real(r8), intent(in) :: f(pkdim)          ! constituent field
    real(r8), intent(in) :: lbasdz(4,2,pkdim) ! vertical interpolation weights

    real(r8), intent(out) :: dfz1(pkdim)      ! derivative at top of interval
    real(r8), intent(out) :: dfz2(pkdim)      ! derivative at bot of interval
! local

    integer k               ! indices

    do k=2,pkdim-2
! Lagrangian derivative estimates (cubic) for the two center nodes in a
! four node stencil.

       dfz1(k) = lbasdz(1,1,k)*f(k-1) +  &
          lbasdz(2,1,k)*f(k)   +  &
          lbasdz(3,1,k)*f(k+1) +  &
          lbasdz(4,1,k)*f(k+2)

       dfz2(k) = lbasdz(1,2,k)*f(k-1) +  &
          lbasdz(2,2,k)*f(k)   +  &
          lbasdz(3,2,k)*f(k+1) +  &
          lbasdz(4,2,k)*f(k+2)
    end do

! Constrain derivatives to zero at top and bottom of vertical grid.
! At the interior nodes of the intervals at the top and bottom of the
! vertical grid, use the derivative estimate at that same node for the
! adjacent interval.  (This is a "one-sided" estimate for that node.)
!
    dfz1(1) = 0.0_r8
    dfz2(1) = dfz1(2)
    dfz1(pkdim-1) = dfz2(pkdim-2)
    dfz2(pkdim-1) = 0.0_r8
! LiXH add:
    dfz1(pkdim) = 0._r8
    dfz2(pkdim) = 0._r8
    return

    end subroutine cubzdr



! Purpose: Compute vertical departure point and departure point index.
    subroutine vrtdep(pmap    ,dt      ,iterdp  ,wb      ,wst     , &
                      wsb     ,sig     ,sigh    ,dsigh   ,kdpmpf  , & 
                      kdpmph  ,sigmp   ,sigdp   ,kdp              )
! io
    integer , intent(in) :: pmap           ! dimension of artificial vert grid
    real(r8), intent(in) :: dt             ! time step (seconds)
    integer , intent(in) :: iterdp         ! number of iterations
    real(r8), intent(in) :: wb (nlevp)     ! vertical velocity
    real(r8), intent(in) :: wst(nlevp)     ! z-derivative of wb at top of interval
    real(r8), intent(in) :: wsb(nlevp)     ! z-derivative of wb at bot of interval
    real(r8), intent(in) :: sig  (nlev )   ! sigma values of model full levels
    real(r8), intent(in) :: sigh (nlevp)   ! sigma values of model half levels
    real(r8), intent(in) :: dsigh(nlevp)   ! increment between half levels
    integer , intent(in) :: kdpmpf(pmap)   ! artificial grid indices
    integer , intent(in) :: kdpmph(pmap)   ! artificial grid indices
    real(r8), intent(inout) :: sigmp(nlev) ! vert coords of traj mid-points
    real(r8), intent(out) :: sigdp(nlev)   ! vert coords of traj departure points
    integer , intent(out) :: kdp(nlev)     ! vertical departure point indices
! local
    integer iter             
    integer k                 
    real(r8) wmp(nlev)   ! vert vel. at midpoint

    do iter = 1,iterdp
! Compute midpoint indices in half-index sigma-level arrays (use kdp
! as temporary storage).

        call kdpfnd(nlevp   ,pmap    ,sigh    ,sigmp   ,kdpmph  ,  kdp)

! Interpolate sigma dot field to trajectory midpoints using Hermite
! cubic interpolant.
        call herzin(nlevp   ,1       ,wb      ,wst     ,wsb     , &
                    sigh    ,dsigh   ,sigmp   ,kdp     ,wmp     )

! Update estimate of trajectory midpoint.
        do k = 1,nlev
            sigmp(k) = sig(k) - .5_r8*dt*wmp(k)
        end do

! Restrict vertical midpoints to be between the top and bottom half-
! index sigma levels.
        call vdplim(nlevp   ,sigh    ,sigmp  )
    end do

! Compute trajectory endpoints.
    do k = 1,nlev
        sigdp(k) = sig(k) - dt*wmp(k)
    end do

! Restrict vertical departure points to be between the top and bottom
! full-index sigma levels.

    call vdplim(nlev    ,sig     ,sigdp )

! Vertical indices for trajectory endpoints that point into full-index
! sigma level arrays.
!
    call kdpfnd(nlev    ,pmap    ,sig     ,sigdp   ,kdpmpf  , kdp )

    return
    end subroutine vrtdep


! Purpose: 
! Restrict vertical departure points to be between the top and bottom
! sigma levels of the "full-" or "half-" level grid
    subroutine vdplim(pkdim   ,sig     ,sigdp    )
! io
    integer , intent(in)    :: pkdim              ! vertical dimension
    real(r8), intent(in)    :: sig(pkdim)         ! vertical coordinate of model grid
    real(r8), intent(inout) :: sigdp(nlev)        ! vertical coords. of departure points.
! local
    integer k                 ! index

    do k=1,nlev
       if (sigdp(k) < sig(1)) then
          sigdp(k) = sig(1)
       end if
       if (sigdp(k) >= sig(pkdim)) then
          sigdp(k) = sig(pkdim)*(1._r8 - 10._r8*epsilon(sigdp))
       end if
    end do

    return
    end subroutine vdplim


! Purpose: 
! Determine vertical departure point indices that point into a grid
! containing the full or half sigma levels.  Use an artificial evenly 
! spaced vertical grid to map into the true model levels.
! 
! Method: 
! Indices are computed assuming the the sigdp values have
! been constrained so that sig(1) .le. sigdp(i,j) .lt. sig(pkdim).
    subroutine kdpfnd(pkdim   ,pmap    ,sig     ,sigdp   ,kdpmap  ,kdp  )
! io
    integer , intent(in) :: pkdim             ! dimension of "sig"
    integer , intent(in) :: pmap              ! dimension of "kdpmap"
    real(r8), intent(in) :: sig  (pkdim)      ! vertical grid coordinates
    integer , intent(in) :: kdpmap(pmap)      ! array of model grid indices which
    real(r8), intent(in) :: sigdp(nlev)       ! vertical coords. of departure points
    integer , intent(out):: kdp(nlev)         ! vertical index for each dep. pt.
! local
    integer k,ii            ! indices
    real(r8) rdel             ! recip. of interval in artificial grid
    real(r8) sig1ln           ! ln (sig(1))

    rdel   = real(pmap,r8)/( log(sig(pkdim)) - log(sig(1)) )
    sig1ln = log( sig(1) )

    do k=1,nlev
! First guess of the departure point's location in the model grid
        ii = max0(1,min0(pmap,int((log(sigdp(k))-sig1ln)*rdel+1._r8)))
        kdp(k) = kdpmap(ii)

! Determine if location is in next interval
        if(sigdp(k) >= sig(kdp(k)+1)) then
           kdp(k) = kdp(k) + 1
        end if
    end do

    return
    end subroutine kdpfnd


! Purpose: 
! Interpolate field on vertical slice to vertical departure point using
! Hermite cubic interpolation.
    subroutine herzin(pkdim   ,pf      ,f       ,fst     ,fsb     , &
                      sig     ,dsig    ,sigdp   ,kdp     ,fdp     )
! io
    integer, intent(in)  :: pkdim           ! vertical dimension
    integer, intent(in)  :: pf              ! dimension (number of fields)
    real(r8), intent(in) :: f    (pf,pkdim) ! fields
    real(r8), intent(in) :: fst  (pf,pkdim) ! z-derivatives at top edge of interval
    real(r8), intent(in) :: fsb  (pf,pkdim) ! z-derivatives at bot edge of interval
    real(r8), intent(in) :: sig  (pkdim)    ! vertical grid coordinates
    real(r8), intent(in) :: dsig (pkdim)    ! intervals between vertical grid pts.
    real(r8), intent(in) :: sigdp(nlev)     ! vertical coord. of departure point
    integer, intent(in)  :: kdp  (nlev)     ! vertical index  of departure point

    real(r8), intent(out) :: fdp(pf,nlev)   ! z-interpolants
! local
    integer k,m             ! indices
    real(r8) dzk              ! vert interval containing the dep. pt.
    real(r8) zt               
    real(r8) zb               
    real(r8) ht (1)
    real(r8) hb (1)            
    real(r8) dht(1)            
    real(r8) dhb(1)            

    do k=1,nlev
       dzk = dsig(kdp(k))
       zt = ( sig(kdp(k)+1) - sigdp(k) )/dzk
       zb = 1._r8 - zt
       ht (1) = ( 3.0_r8 - 2.0_r8*zt )*zt**2
       hb (1) = ( 3.0_r8 - 2.0_r8*zb )*zb**2
       dht(1) = -dzk*( zt - 1._r8 )*zt**2
       dhb(1) =  dzk*( zb - 1._r8 )*zb**2

       do m=1,pf
           fdp(m,k) = f(m,kdp(k))* ht(1) +  &
           fst(m,kdp(k))*dht(1) +  &
           f(m,kdp(k)+1)* hb(1) +  &
           fsb(m,kdp(k))*dhb(1)
       end do
    end do

    return
    end subroutine herzin


! Purpose: 
! Apply SCMO limiter to vertical derivative estimates on a vertical slice.
    subroutine limdz(f       ,dsig    ,fst     ,fsb     )

! io
    real(r8), intent(in) :: f(ntracer,nlev)   ! input field
    real(r8), intent(in) :: dsig(nlev)        ! size of vertical interval

    real(r8), intent(inout) :: fst(ntracer,nlev) ! z-derivative at top of interval
    real(r8), intent(inout) :: fsb(ntracer,nlev) ! z-derivative at bot of interval
! local
    integer k                 ! vertical    index
    integer m                 ! constituent index
    integer plevm1
    real(r8) rdsig                ! 1./dsig
    real(r8) deli           ! simple linear derivative
    real(r8) fac,tmp1,tmp2
   
    plevm1 = nlev - 1
    fac = 3._r8*(1._r8 - 10._r8*epsilon(fac))

    do m = 1, ntracer
       do k = 1, nlev-1
          rdsig = 1.0_r8/dsig(k)
          deli  = ( f(m,k+1) - f(m,k) )*rdsig
          tmp1 = fac*deli
          tmp2 = abs( tmp1 )
          if( deli*fst(m,k)   <= 0.0_r8  ) fst(m,k) = 0._r8
          if( deli*fsb(m,k)   <= 0.0_r8  ) fsb(m,k) = 0._r8
          if( abs(    fst(m,k) ) >  tmp2 ) fst(m,k) = tmp1
          if( abs(    fsb(m,k) ) >  tmp2 ) fsb(m,k) = tmp1
       end do
    end do

    return
    end subroutine limdz


   integer function ismin(n,sx,incx)
! io
   real(r8), intent(in) :: sx(n)  ! array to be searched           

   integer n                      ! number of values to be searched
   integer incx                   ! increment through array        
! local
   integer imin                   ! tmp index of min value         
   integer i                      ! loop index                     

   real(r8) min                   ! minimum value                  
   
   min = sx(1)
   imin = 1
   do 100 i=1+incx,n,incx
      if (min .gt. sx(i)) then
         min = sx(i)
         imin = i
      endif
100   continue
      ismin = imin
      return
   end function ismin
 

 end module grist_scm_scanslt
