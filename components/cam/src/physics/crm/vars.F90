module vars
  use grid
  use params, only: crm_rknd
#ifdef CRM
#ifdef MODAL_AERO
  use modal_aero_data,   only: ntot_amode
#endif
#endif

  implicit none
  !--------------------------------------------------------------------
  ! prognostic variables:

  real(crm_rknd), allocatable :: u   (:,:,:,:) ! x-wind
  real(crm_rknd), allocatable :: v   (:,:,:,:) ! y-wind
  real(crm_rknd), allocatable :: w   (:,:,:,:) ! z-wind
  real(crm_rknd), allocatable :: t   (:,:,:,:) ! liquid/ice water static energy

  !--------------------------------------------------------------------
  ! diagnostic variables:

  real(crm_rknd), allocatable :: p       (:,:,:,:)     ! perturbation pressure (from Poison eq)
  real(crm_rknd), allocatable :: tabs    (:,:,:,:)                 ! temperature
  real(crm_rknd), allocatable :: qv      (:,:,:,:)                ! water vapor
  real(crm_rknd), allocatable :: qcl     (:,:,:,:)                ! liquid water  (condensate)
  real(crm_rknd), allocatable :: qpl     (:,:,:,:)                ! liquid water  (precipitation)

  real(crm_rknd), allocatable :: qci     (:,:,:,:)                ! ice water  (condensate)
  real(crm_rknd), allocatable :: qpi     (:,:,:,:)                ! ice water  (precipitation)
  real(crm_rknd), allocatable :: tke2    (:,:,:,:)   ! SGS TKE
  real(crm_rknd), allocatable :: tk2     (:,:,:,:) ! SGS eddyviscosity

  !--------------------------------------------------------------------
  ! time-tendencies for prognostic variables

  real(crm_rknd), allocatable :: dudt   (:,:,:,:,:)
  real(crm_rknd), allocatable :: dvdt   (:,:,:,:,:)
  real(crm_rknd), allocatable :: dwdt   (:,:,:,:,:)

  !----------------------------------------------------------------
  ! Temporary storage array:

  real(crm_rknd), allocatable :: misc(:,:,:,:)
  !------------------------------------------------------------------
  ! fluxes at the top and bottom of the domain:

  real(crm_rknd), allocatable :: fluxbu  (:,:,:)
  real(crm_rknd), allocatable :: fluxbv  (:,:,:)
  real(crm_rknd), allocatable :: fluxbt  (:,:,:)
  real(crm_rknd), allocatable :: fluxbq  (:,:,:)
  real(crm_rknd), allocatable :: fluxtu  (:,:,:)
  real(crm_rknd), allocatable :: fluxtv  (:,:,:)
  real(crm_rknd), allocatable :: fluxtt  (:,:,:)
  real(crm_rknd), allocatable :: fluxtq  (:,:,:)
  real(crm_rknd), allocatable :: fzero   (:,:,:)
  real(crm_rknd), allocatable :: precsfc (:,:,:) ! surface precip. rate
  real(crm_rknd), allocatable :: precssfc(:,:,:) ! surface ice precip. rate

  !-----------------------------------------------------------------
  ! profiles

  real(crm_rknd), allocatable :: t0   (:,:)
  real(crm_rknd), allocatable :: q0   (:,:)
  real(crm_rknd), allocatable :: qv0  (:,:)
  real(crm_rknd), allocatable :: tabs0(:,:)
  real(crm_rknd), allocatable :: tv0  (:,:)
  real(crm_rknd), allocatable :: u0   (:,:)
  real(crm_rknd), allocatable :: v0   (:,:)
  real(crm_rknd), allocatable :: tg0  (:,:)
  real(crm_rknd), allocatable :: qg0  (:,:)
  real(crm_rknd), allocatable :: ug0  (:,:)
  real(crm_rknd), allocatable :: vg0  (:,:)
  real(crm_rknd), allocatable :: p0   (:,:)
  real(crm_rknd), allocatable :: tke0 (:,:)
  real(crm_rknd), allocatable :: t01  (:,:)
  real(crm_rknd), allocatable :: q01  (:,:)
  real(crm_rknd), allocatable :: qp0  (:,:)
  real(crm_rknd), allocatable :: qn0  (:,:)
  !-----------------------------------------------------------------
  ! reference vertical profiles:
  real(crm_rknd), allocatable :: prespot(:,:)  ! (1000./pres)**R/cp
  real(crm_rknd), allocatable :: rho    (:,:)   ! air density at pressure levels,kg/m3
  real(crm_rknd), allocatable :: rhow   (:,:)   ! air density at vertical velocity levels,kg/m3
  real(crm_rknd), allocatable :: bet    (:,:)   ! = ggr/tv0
  real(crm_rknd), allocatable :: gamaz  (:,:) ! ggr/cp*z
  real(crm_rknd), allocatable :: wsub   (:,:)   ! Large-scale subsidence velocity,m/s
  real(crm_rknd), allocatable :: qtend  (:,:) ! Large-scale tendency for total water
  real(crm_rknd), allocatable :: ttend  (:,:) ! Large-scale tendency for temp.
  real(crm_rknd), allocatable :: utend  (:,:) ! Large-scale tendency for u
  real(crm_rknd), allocatable :: vtend  (:,:) ! Large-scale tendency for v

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !---------------------------------------------------------------------
  !  Horizontally varying stuff (as a function of xy)
  !
  real(crm_rknd), allocatable :: sstxy    (:,:,:) !  surface temperature xy-distribution
  real(crm_rknd), allocatable :: fcory    (:,:)      !  Coriolis parameter xy-distribution
  real(crm_rknd), allocatable :: fcorzy   (:,:)      !  z-Coriolis parameter xy-distribution
  real(crm_rknd), allocatable :: latitude (:,:,:)      ! latitude (degrees,:)
  real(crm_rknd), allocatable :: longitude(:,:,:)      ! longitude(degrees,:)
  real(crm_rknd), allocatable :: prec_xy  (:,:,:) ! mean precip. rate for outout
  real(crm_rknd), allocatable :: pw_xy    (:,:,:)   ! precipitable water
  real(crm_rknd), allocatable :: cw_xy    (:,:,:)   ! cloud water path
  real(crm_rknd), allocatable :: iw_xy    (:,:,:)   ! ice water path
  real(crm_rknd), allocatable :: cld_xy   (:,:,:)   ! cloud frequency
  real(crm_rknd), allocatable :: u200_xy  (:,:,:) ! u-wind at 200 mb
  real(crm_rknd), allocatable :: usfc_xy  (:,:,:) ! u-wind at at the surface
  real(crm_rknd), allocatable :: v200_xy  (:,:,:) ! v-wind at 200 mb
  real(crm_rknd), allocatable :: vsfc_xy  (:,:,:) ! v-wind at the surface
  real(crm_rknd), allocatable :: w500_xy  (:,:,:) ! w at 500 mb
  !----------------------------------------------------------------------
  ! Vertical profiles of quantities sampled for statitistics purposes:

  real(crm_rknd), allocatable :: w_max(:)
  real(crm_rknd), allocatable :: u_max(:)

  real(crm_rknd), allocatable :: twsb(:,:)
  real(crm_rknd), allocatable :: precflux(:,:)
  real(crm_rknd), allocatable :: uwle(:,:)
  real(crm_rknd), allocatable :: uwsb(:,:)
  real(crm_rknd), allocatable :: vwle(:,:)
  real(crm_rknd), allocatable :: vwsb(:,:)
  real(crm_rknd), allocatable :: tkelediss(:,:)
  real(crm_rknd), allocatable :: tdiff(:,:)
  real(crm_rknd), allocatable :: tlat(:,:)
  real(crm_rknd), allocatable :: tlatqi(:,:)
  real(crm_rknd), allocatable :: qifall(:,:)
  real(crm_rknd), allocatable :: qpfall(:,:)

  ! energy conservation diagnostics:
  real(8), allocatable :: total_water_evap(:)
  real(8), allocatable :: total_water_prec(:)

  real(crm_rknd), allocatable :: CF3D(:,:,:,:)  ! Cloud fraction
  ! =1.0 when there is no fractional cloudiness scheme
  ! = cloud fraction produced by fractioal cloudiness scheme when avaiable

  ! 850 mbar horizontal winds
  real(crm_rknd), allocatable :: u850_xy(:,:,:) ! zonal velocity at 850 mb
  real(crm_rknd), allocatable :: v850_xy(:,:,:) ! meridional velocity at 850 mb

  ! Surface pressure
  real(crm_rknd), allocatable :: psfc_xy(:,:,:) ! pressure (in millibar) at lowest grid point
  ! Saturated water vapor path, useful for computing column relative humidity
  real(crm_rknd), allocatable :: swvp_xy(:,:,:)  ! saturated water vapor path (wrt water)
  ! Cloud and echo top heights, and cloud top temperature (instantaneous)
  real(crm_rknd), allocatable :: cloudtopheight(:,:,:)
  real(crm_rknd), allocatable :: echotopheight (:,:,:)
  real(crm_rknd), allocatable :: cloudtoptemp  (:,:,:)

  ! END UW ADDITIONS
  !===========================================================================
#if (defined CRM && defined MODAL_AERO)
  real(crm_rknd), allocatable :: naer (:,:,:)     ! Aerosol number concentration [/m3]
  real(crm_rknd), allocatable :: vaer (:,:,:)     ! aerosol volume concentration [m3/m3]
  real(crm_rknd), allocatable :: hgaer(:,:,:)    ! hygroscopicity of aerosol mode
#endif


contains

  subroutine allocate_vars(nx_gl_in,ny_gl_in,nz_gl_in,ncrms)
    use openacc_utils
    implicit none
    integer, intent(in) :: nx_gl_in,ny_gl_in,nz_gl_in,ncrms
    integer :: test_out
    real(crm_rknd) :: zero
    YES3D = 0
    if (ny_gl_in.gt.1) YES3D = 1

    allocate( u(ncrms,dimx1_u:dimx2_u,dimy1_u:dimy2_u,nz_gl_in-1)  )
    allocate( v(ncrms,dimx1_v:dimx2_v,dimy1_v:dimy2_v,nz_gl_in-1)  )
    allocate( w(ncrms,dimx1_w:dimx2_w,dimy1_w:dimy2_w,nz_gl_in )  )
    allocate( t(ncrms,dimx1_s:dimx2_s,dimy1_s:dimy2_s,nz_gl_in-1)  )
    allocate( p       (ncrms,0:nx_gl_in, (1-YES3D):ny_gl_in, nz_gl_in-1)      )
    allocate( tabs(ncrms,nx_gl_in, ny_gl_in, nz_gl_in-1)                  )
    allocate( qv(ncrms,  nx_gl_in, ny_gl_in, nz_gl_in-1)                 )
    allocate( qcl(ncrms, nx_gl_in, ny_gl_in, nz_gl_in-1)                 )
    allocate( qpl(ncrms, nx_gl_in, ny_gl_in, nz_gl_in-1)                 )
    allocate( qci(ncrms, nx_gl_in, ny_gl_in, nz_gl_in-1)                 )
    allocate( qpi(ncrms, nx_gl_in, ny_gl_in, nz_gl_in-1)                 )
    allocate( tke2(ncrms,dimx1_s:dimx2_s, dimy1_s:dimy2_s, nz_gl_in-1)    )
    allocate( tk2  (ncrms,0:nx_gl_in+1, (1-YES3D):nyp1, nz_gl_in-1)  )
    allocate( dudt(ncrms,nx_gl_in+1, ny_gl_in, nz_gl_in-1, 3) )
    allocate( dvdt(ncrms,nx_gl_in, nyp1, nz_gl_in-1, 3) )
    allocate( dwdt(ncrms,nx_gl_in, ny_gl_in  , nz_gl_in,  3) )
    allocate( misc(ncrms,nx_gl_in, ny_gl_in, nz_gl_in) )
    allocate( fluxbu(ncrms,nx_gl_in,ny_gl_in) )
    allocate( fluxbv(ncrms,nx_gl_in,ny_gl_in) )
    allocate( fluxbt(ncrms,nx_gl_in,ny_gl_in) )
    allocate( fluxbq(ncrms,nx_gl_in,ny_gl_in) )
    allocate( fluxtu(ncrms,nx_gl_in,ny_gl_in) )
    allocate( fluxtv(ncrms,nx_gl_in,ny_gl_in) )
    allocate( fluxtt(ncrms,nx_gl_in,ny_gl_in) )
    allocate( fluxtq(ncrms,nx_gl_in,ny_gl_in) )
    allocate( fzero(ncrms,nx_gl_in,ny_gl_in) )
    allocate( precsfc(ncrms,nx_gl_in,ny_gl_in)  )
    allocate( precssfc(ncrms,nx_gl_in,ny_gl_in)  )
    allocate( t0(ncrms,nz_gl_in-1) )
    allocate( q0(ncrms,nz_gl_in-1) )
    allocate( qv0(ncrms,nz_gl_in-1) )
    allocate( tabs0(ncrms,nz_gl_in-1) )
    allocate( tv0(ncrms,nz_gl_in-1) )
    allocate( u0(ncrms,nz_gl_in-1) )
    allocate( v0(ncrms,nz_gl_in-1) )
    allocate( tg0(ncrms,nz_gl_in-1) )
    allocate( qg0(ncrms,nz_gl_in-1) )
    allocate( ug0(ncrms,nz_gl_in-1) )
    allocate( vg0(ncrms,nz_gl_in-1) )
    allocate( p0(ncrms,nz_gl_in-1) )
    allocate( tke0(ncrms,nz_gl_in-1) )
    allocate( t01(ncrms,nz_gl_in-1) )
    allocate( q01(ncrms,nz_gl_in-1) )
    allocate( qp0(ncrms,nz_gl_in-1) )
    allocate( qn0(ncrms,nz_gl_in-1) )
   
    allocate( prespot(ncrms,nz_gl_in-1)   )
    allocate( rho(ncrms,nz_gl_in-1)     )
    allocate( rhow(ncrms,nz_gl_in )    )
    allocate( bet(ncrms,nz_gl_in-1)     )
    allocate( gamaz(ncrms,nz_gl_in-1)  )
    allocate( wsub(ncrms,nz_gl_in )    )
    allocate( qtend(ncrms,nz_gl_in-1)  )
    allocate( ttend(ncrms,nz_gl_in-1)  )
    allocate( utend(ncrms,nz_gl_in-1)  )
    allocate( vtend(ncrms,nz_gl_in-1)  )
    allocate( sstxy    (ncrms,0:nx_gl_in,(1-YES3D):ny_gl_in)   )
    allocate( fcory(ncrms,0:ny_gl_in)       )
    allocate( fcorzy(ncrms,ny_gl_in)       )
    allocate( latitude(ncrms,nx_gl_in,ny_gl_in)        )
    allocate( longitude(ncrms,nx_gl_in,ny_gl_in)        )
    allocate( prec_xy(ncrms,nx_gl_in,ny_gl_in)  )
    allocate( pw_xy(ncrms,nx_gl_in,ny_gl_in)    )
    allocate( cw_xy(ncrms,nx_gl_in,ny_gl_in)    )
    allocate( iw_xy(ncrms,nx_gl_in,ny_gl_in)    )
    allocate( cld_xy(ncrms,nx_gl_in,ny_gl_in)    )
    allocate( u200_xy(ncrms,nx_gl_in,ny_gl_in)  )
    allocate( usfc_xy(ncrms,nx_gl_in,ny_gl_in)  )
    allocate( v200_xy(ncrms,nx_gl_in,ny_gl_in)  )
    allocate( vsfc_xy(ncrms,nx_gl_in,ny_gl_in)  )
    allocate( w500_xy(ncrms,nx_gl_in,ny_gl_in)  )
    allocate( twsb(ncrms,nz_gl_in) )
    allocate( precflux(ncrms,nz_gl_in) )
    allocate( uwle(ncrms,nz_gl_in) )
    allocate( uwsb(ncrms,nz_gl_in) )
    allocate( vwle(ncrms,nz_gl_in) )
    allocate( vwsb(ncrms,nz_gl_in) )
    allocate( tkelediss(ncrms,nz_gl_in) )
    allocate( tdiff    (ncrms,nz_gl_in) )
    allocate( tlat(ncrms,nz_gl_in) )
    allocate( tlatqi(ncrms,nz_gl_in) )
    allocate( qifall(ncrms,nz_gl_in) )
    allocate( qpfall(ncrms,nz_gl_in) )
    allocate( cf3d(ncrms,1:nx_gl_in, 1:ny_gl_in, 1:nz_gl_in-1)   )
    allocate( u850_xy(ncrms,nx_gl_in,ny_gl_in)  )
    allocate( v850_xy(ncrms,nx_gl_in,ny_gl_in)  )
    allocate( psfc_xy(ncrms,nx_gl_in,ny_gl_in)  )
    allocate( swvp_xy(ncrms,nx_gl_in,ny_gl_in)   )
    allocate( cloudtopheight(ncrms,nx_gl_in,ny_gl_in) )
    allocate( echotopheight(ncrms,nx_gl_in,ny_gl_in) )
    allocate( cloudtoptemp(ncrms,nx_gl_in,ny_gl_in) )
    allocate( u_max(ncrms) )
    allocate( w_max(ncrms) )
    allocate( total_water_evap(ncrms) )
    allocate( total_water_prec(ncrms) )
#if (defined CRM && defined MODAL_AERO)
    allocate( naer(ncrms,nz_gl_in-1, ntot_amode) )
    allocate( vaer(ncrms,nz_gl_in-1, ntot_amode) )
    allocate( hgaer(ncrms,nz_gl_in-1, ntot_amode) )
#endif

    call prefetch( u )
    call prefetch( v )
    call prefetch( w )
    call prefetch( t )
    call prefetch( p )
    call prefetch( tabs )
    call prefetch( qv )
    call prefetch( qcl )
    call prefetch( qpl )
    call prefetch( qci )
    call prefetch( qpi )
    call prefetch( tke2 )
    call prefetch( tk2 )
    call prefetch( dudt )
    call prefetch( dvdt )
    call prefetch( dwdt )
    call prefetch( misc )
    call prefetch( fluxbu )
    call prefetch( fluxbv )
    call prefetch( fluxbt )
    call prefetch( fluxbq )
    call prefetch( fluxtu )
    call prefetch( fluxtv )
    call prefetch( fluxtt )
    call prefetch( fluxtq )
    call prefetch( fzero )
    call prefetch( precsfc )
    call prefetch( precssfc )
    call prefetch( t0 )
    call prefetch( q0 )
    call prefetch( qv0 )
    call prefetch( tabs0 )
    call prefetch( tv0 )
    call prefetch( u0 )
    call prefetch( v0 )
    call prefetch( tg0 )
    call prefetch( qg0 )
    call prefetch( ug0 )
    call prefetch( vg0 )
    call prefetch( p0 )
    call prefetch( tke0 )
    call prefetch( t01 )
    call prefetch( q01 )
    call prefetch( qp0 )
    call prefetch( qn0 )
    call prefetch( prespot )
    call prefetch( rho )
    call prefetch( rhow )
    call prefetch( bet )
    call prefetch( gamaz )
    call prefetch( wsub )
    call prefetch( qtend )
    call prefetch( ttend )
    call prefetch( utend )
    call prefetch( vtend )
    call prefetch( sstxy )
    call prefetch( fcory )
    call prefetch( fcorzy )
    call prefetch( latitude )
    call prefetch( longitude )
    call prefetch( prec_xy )
    call prefetch( pw_xy )
    call prefetch( cw_xy )
    call prefetch( iw_xy )
    call prefetch( cld_xy )
    call prefetch( u200_xy )
    call prefetch( usfc_xy )
    call prefetch( v200_xy )
    call prefetch( vsfc_xy )
    call prefetch( w500_xy )
    call prefetch( twsb )
    call prefetch( precflux )
    call prefetch( uwle )
    call prefetch( uwsb )
    call prefetch( vwle )
    call prefetch( vwsb )
    call prefetch( tkelediss )
    call prefetch( tdiff )
    call prefetch( tlat )
    call prefetch( tlatqi )
    call prefetch( qifall )
    call prefetch( qpfall )
    call prefetch( CF3D )
    call prefetch( u850_xy )
    call prefetch( v850_xy )
    call prefetch( psfc_xy )
    call prefetch( swvp_xy )
    call prefetch( cloudtopheight )
    call prefetch( echotopheight )
    call prefetch( cloudtoptemp )
    call prefetch( u_max )
    call prefetch( w_max )
    call prefetch( total_water_evap )
    call prefetch( total_water_prec )
#if (defined CRM && defined MODAL_AERO)
    call prefetch( naer )
    call prefetch( vaer )
    call prefetch( hgaer  )
#endif

    zero = 0

    u = zero
    v = zero
    w = zero
    t = zero
    p = zero
    tabs = zero
    qv = zero
    qcl = zero
    qpl = zero
    qci = zero
    qpi = zero
    tke2 = zero
    tk2 = zero
    dudt = zero
    dvdt = zero
    dwdt = zero
    misc = zero
    fluxbu = zero
    fluxbv = zero
    fluxbt = zero
    fluxbq = zero
    fluxtu = zero
    fluxtv = zero
    fluxtt = zero
    fluxtq = zero
    fzero = zero
    precsfc = zero
    precssfc = zero
    t0 = zero
    q0 = zero
    qv0 = zero
    tabs0 = zero
    tv0 = zero
    u0 = zero
    v0 = zero
    tg0 = zero
    qg0 = zero
    ug0 = zero
    vg0 = zero
    p0 = zero
    tke0 = zero
    t01 = zero
    q01 = zero
    qp0 = zero
    qn0 = zero
    prespot = zero
    rho = zero
    rhow = zero
    bet = zero
    gamaz = zero
    wsub = zero
    qtend = zero
    ttend = zero
    utend = zero
    vtend = zero
    sstxy = zero
    fcory = zero
    fcorzy = zero
    latitude = zero
    longitude = zero
    prec_xy = zero
    pw_xy = zero
    cw_xy = zero
    iw_xy = zero
    cld_xy = zero
    u200_xy = zero
    usfc_xy = zero
    v200_xy = zero
    vsfc_xy = zero
    w500_xy = zero
    twsb = zero
    precflux = zero
    uwle = zero
    uwsb = zero
    vwle = zero
    vwsb = zero
    tkelediss = zero
    tdiff = zero
    tlat = zero
    tlatqi = zero
    qifall = zero
    qpfall = zero
    CF3D = 1.
    u850_xy = zero
    v850_xy = zero
    psfc_xy = zero
    swvp_xy = zero
    cloudtopheight = zero
    echotopheight = zero
    cloudtoptemp = zero
    u_max = zero
    w_max = zero
    total_water_evap = zero
    total_water_prec = zero
#if (defined CRM && defined MODAL_AERO)
    naer = zero
    vaer = zero
    hgaer = zero
#endif
  end subroutine allocate_vars


  subroutine deallocate_vars()
    implicit none

    deallocate(u)
    deallocate(v)
    deallocate(w)
    deallocate(t)
    deallocate(p)
    deallocate(tabs)
    deallocate(qv)
    deallocate(qcl)
    deallocate(qpl)
    deallocate(qci)
    deallocate(qpi)
    deallocate(dudt)
    deallocate(dvdt)
    deallocate(dwdt)
    deallocate(misc)

    deallocate( fluxbu )
    deallocate( fluxbv )
    deallocate( fluxbt )
    deallocate( fluxbq )
    deallocate( fluxtu )
    deallocate( fluxtv )
    deallocate( fluxtt )
    deallocate( fluxtq )
    deallocate( fzero )
    deallocate( precsfc )
    deallocate( precssfc )
    deallocate( t0 )
    deallocate( q0 )
    deallocate( qv0 )
    deallocate( tabs0 )
    deallocate( tv0 )
    deallocate( u0 )
    deallocate( v0 )
    deallocate( tg0 )
    deallocate( qg0 )
    deallocate( ug0 )
    deallocate( vg0 )
    deallocate( p0 )
    deallocate( tke0 )
    deallocate( t01 )
    deallocate( q01 )
    deallocate( qp0 )
    deallocate( qn0 )
    deallocate( prespot )
    deallocate( rho )
    deallocate( rhow )
    deallocate( bet )
    deallocate( gamaz )
    deallocate( wsub )
    deallocate( qtend )
    deallocate( ttend )
    deallocate( utend )
    deallocate( vtend )
    deallocate( sstxy )
    deallocate( fcory )
    deallocate( fcorzy )
    deallocate( latitude )
    deallocate( longitude )
    deallocate( prec_xy )
    deallocate( pw_xy )
    deallocate( cw_xy )
    deallocate( iw_xy )
    deallocate( cld_xy )
    deallocate( u200_xy )
    deallocate( usfc_xy )
    deallocate( v200_xy )
    deallocate( vsfc_xy )
    deallocate( w500_xy )
    deallocate( twsb )
    deallocate( precflux )
    deallocate( uwle )
    deallocate( uwsb )
    deallocate( vwle )
    deallocate( vwsb )
    deallocate( tkelediss )
    deallocate( tdiff )
    deallocate( tlat )
    deallocate( tlatqi )
    deallocate( qifall )
    deallocate( qpfall )
    deallocate( CF3D )
    deallocate( u850_xy )
    deallocate( v850_xy )
    deallocate( psfc_xy )
    deallocate( swvp_xy )
    deallocate( cloudtopheight )
    deallocate( echotopheight )
    deallocate( cloudtoptemp )
    deallocate( u_max )
    deallocate( w_max )
    deallocate( total_water_evap )
    deallocate( total_water_prec )
#if (defined CRM && defined MODAL_AERO)
    deallocate( naer )
    deallocate( vaer )
    deallocate( hgaer  )
#endif
end subroutine deallocate_vars

end module vars
