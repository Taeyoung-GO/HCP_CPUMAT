!=============================================================================
! HCP_CPUMAT.f90  -  Consolidated HCP Crystal Plasticity UMAT for ABAQUS
!=============================================================================
! HCP phase only | Sinh slip law | Kocks-Mecking hardening | GND optional
!
! ---- PROPS (31 entries) ---------------------------------------------------
!  [1-3]   phi1, Phi, phi2        Bunge Euler angles (degrees)
!  [4]     nslip                  # HCP slip systems: 3/6/12/18/24/30
!  [5]     caratio                c/a ratio
!  [6]     temperature            K  (-1 = use ABAQUS field)
!  [7-11]  C11,C12,C13,C33,C44   Elastic constants (MPa)
!  [12]    alpha0                 Sinh pre-factor  (0 = parametric)
!  [13]    beta0                  Sinh exponent    (0 = parametric)
!  [14]    rhom                   Mobile disloc. density (1/um^2)
!  [15]    DeltaF                 Activation energy (J)
!  [16]    nu0                    Attempt frequency (1/s)
!  [17]    AV0                    Activation volume (0 = parametric, b^3)
!  [18]    gamma0                 Geometric factor for AV
!  [19]    k1                     KM storage coeff (1/um)
!  [20]    k2_const               KM annihilation  (0 = parametric)
!  [21]    X_km                   Annihilation factor
!  [22]    g_km                   Boltzmann activation factor
!  [23]    D_km                   Drag stress parameter (J)
!  [24]    pdot0                  Reference strain rate (1/s)
!  [25]    tauc0_basal            Initial CRSS, Basal        (MPa)
!  [26]    tauc0_prism            Initial CRSS, Prismatic    (MPa)
!  [27]    tauc0_pyr1             Initial CRSS, Pyramidal-1  (MPa)
!  [28]    tauc0_pyr2             Initial CRSS, Pyramidal-2  (MPa)
!  [29]    rho0                   Initial SSD density (1/um^2)
!  [30]    gf                     Geometric/Taylor hardening factor
!  [31]    gndmodel               0=off  1=element-average curl(Fp)
!
! ---- STATEV (NSTATV = 4*nslip + 25) --------------------------------------
!  [1:ns]           gammasum   Cumulative slip per system
!  [ns+1:2ns]       tauc       Current CRSS (MPa)
!  [2ns+1:3ns]      ssd        SSD density (1/um^2)
!  [3ns+1:4ns]      gnd        GND density (1/um^2)
!  [4ns+1:4ns+9]    Fp         Plastic deformation gradient (row-major 3x3)
!  [4ns+10:4ns+18]  gmat       Orientation matrix (row-major 3x3)
!  [4ns+19:4ns+24]  sig_t      Cauchy stress at previous step (6-vec, MPa)
!  [4ns+25]         evmp       Equivalent plastic strain
!  (ns = nslip)
!=============================================================================

!=============================================================================
! MODULE: hcp_constants  -  Physical constants and numerical settings
!=============================================================================
      module hcp_constants
      implicit none
      real(8), parameter :: PI     = 3.14159
      real(8), parameter :: KB     = 1.38064852d-23   ! Boltzmann (J/K)
      real(8), parameter :: SMALL  = 1.0d-15
      real(8), parameter :: TOL_NR = 1.0d-8           ! NR convergence (MPa)
      integer, parameter :: MAXITER = 200              ! Max NR iterations
      integer, parameter :: MAXNS  = 30               ! Max slip systems
      real(8), parameter :: CUTBACK = 0.25d0
      end module hcp_constants

!=============================================================================
! MODULE: hcp_slip_data  -  HCP slip system vectors (Miller-Bravais 4-index)
!   Ordering:  1-3   Basal      {0001}<11-20>
!              4-6   Prismatic  {10-10}<11-20>
!              7-12  Pyramidal-1 {10-11}<11-20>
!             13-24  Pyramidal-2 {11-22}<11-23>  (first 12)
!             25-30  Pyramidal-2 {11-22}<11-23>  (last 6)
!=============================================================================
      module hcp_slip_data
      implicit none
      ! Slip directions [u v t w]  (30 systems)
      real(8), parameter :: DIR3H(30,4) = reshape( [&
        2d0,-1d0,-1d0, 0d0,  -1d0, 2d0,-1d0, 0d0,  -1d0,-1d0, 2d0, 0d0, & ! Basal 1-3
        2d0,-1d0,-1d0, 0d0,  -1d0, 2d0,-1d0, 0d0,  -1d0,-1d0, 2d0, 0d0, & ! Prism 4-6
       -1d0, 2d0,-1d0, 0d0,  -2d0, 1d0, 1d0, 0d0,  -1d0,-1d0, 2d0, 0d0, & ! Pyr-1 7-9
        1d0,-2d0, 1d0, 0d0,   2d0,-1d0,-1d0, 0d0,   1d0, 1d0,-2d0, 0d0, & ! Pyr-1 10-12
       -2d0, 1d0, 1d0, 3d0,  -1d0,-1d0, 2d0, 3d0,  -1d0,-1d0, 2d0, 3d0, & ! Pyr-2 13-15
        1d0,-2d0, 1d0, 3d0,   1d0,-2d0, 1d0, 3d0,   2d0,-1d0,-1d0, 3d0, & ! Pyr-2 16-18
        2d0,-1d0,-1d0, 3d0,   1d0, 1d0,-2d0, 3d0,   1d0, 1d0,-2d0, 3d0, & ! Pyr-2 19-21
       -1d0, 2d0,-1d0, 3d0,  -1d0, 2d0,-1d0, 3d0,  -2d0, 1d0, 1d0, 3d0, & ! Pyr-2 22-24
       -1d0,-1d0, 2d0, 3d0,   1d0,-2d0, 1d0, 3d0,   2d0,-1d0,-1d0, 3d0, & ! Pyr-2 25-27
        1d0, 1d0,-2d0, 3d0,  -1d0, 2d0,-1d0, 3d0,  -2d0, 1d0, 1d0, 3d0  & ! Pyr-2 28-30
      ], [30,4])
      ! Slip plane normals [h k i l]  (30 systems)
      real(8), parameter :: NOR3H(30,4) = reshape( [&
        0d0, 0d0, 0d0, 1d0,   0d0, 0d0, 0d0, 1d0,   0d0, 0d0, 0d0, 1d0, & ! Basal 1-3
        0d0, 1d0,-1d0, 0d0,  -1d0, 0d0, 1d0, 0d0,   1d0,-1d0, 0d0, 0d0, & ! Prism 4-6
        1d0, 0d0,-1d0, 1d0,   0d0, 1d0,-1d0, 1d0,  -1d0, 1d0, 0d0, 1d0, & ! Pyr-1 7-9
       -1d0, 0d0, 1d0, 1d0,   0d0,-1d0, 1d0, 1d0,   1d0,-1d0, 0d0, 1d0, & ! Pyr-1 10-12
        1d0, 0d0,-1d0, 1d0,   1d0, 0d0,-1d0, 1d0,   0d0, 1d0,-1d0, 1d0, & ! Pyr-2 13-15
        0d0, 1d0,-1d0, 1d0,  -1d0, 1d0, 0d0, 1d0,  -1d0, 1d0, 0d0, 1d0, & ! Pyr-2 16-18
       -1d0, 0d0, 1d0, 1d0,  -1d0, 0d0, 1d0, 1d0,   0d0,-1d0, 1d0, 1d0, & ! Pyr-2 19-21
        0d0,-1d0, 1d0, 1d0,   1d0,-1d0, 0d0, 1d0,   1d0,-1d0, 0d0, 1d0, & ! Pyr-2 22-24
        1d0, 1d0,-2d0, 2d0,  -1d0, 2d0,-1d0, 2d0,  -2d0, 1d0, 1d0, 2d0, & ! Pyr-2 25-27
        1d0,-2d0, 1d0, 2d0,  -1d0,-1d0, 2d0, 2d0,   2d0,-1d0,-1d0, 2d0  & ! Pyr-2 28-30
      ], [30,4])
      end module hcp_slip_data

!=============================================================================
! MODULE: hcp_globals  -  Global state arrays (shared UMAT <-> UEXTERNALDB)
!=============================================================================
      module hcp_globals
      implicit none
      integer, save :: numel = 0, numpt = 0
      integer, save :: g_nslip = 0
      logical, save :: g_initialized = .false.
      ! Stored per (element, ip): Fp, gnd, slip directions/normals, burgers
      real(8), allocatable, save :: glob_Fp(:,:,:,:)      ! (nel,npt,3,3)
      real(8), allocatable, save :: glob_gnd(:,:,:)       ! (nel,npt,ns)
      real(8), allocatable, save :: glob_dirc(:,:,:,:)    ! (nel,npt,ns,3)
      real(8), allocatable, save :: glob_norc(:,:,:,:)    ! (nel,npt,ns,3)
      real(8), allocatable, save :: glob_burgerv(:,:,:)   ! (nel,npt,ns)
      real(8), allocatable, save :: glob_coords(:,:,:)    ! (nel,npt,3)
      integer, allocatable, save :: glob_ipinit(:,:)      ! (nel,npt)
      end module hcp_globals

!=============================================================================
! SUBROUTINE UMAT  -  ABAQUS user material interface
!=============================================================================
      subroutine umat(stress,statev,ddsdde,sse,spd,scd,rpl,ddsddt,
     +   drplde,drpldt,stran,dstran,time,dtime,temp,dtemp,predef,
     +   dpred,cmname,ndi,nshr,ntens,nstatv,props,nprops,coords,
     +   drot,pnewdt,celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,
     +   jstep,kinc)

      use hcp_constants
      use hcp_globals, only: numel, numpt, g_initialized, g_nslip,
     +   glob_Fp, glob_gnd, glob_dirc, glob_norc, glob_burgerv,
     +   glob_coords, glob_ipinit
      implicit none

      character(len=80), intent(in) :: cmname
      integer, intent(in) :: ndi,nshr,ntens,nstatv,nprops
      integer, intent(in) :: noel,npt,layer,kspt
      integer, intent(in) :: jstep(4),kinc
      real(8), intent(in) :: stran(ntens),dstran(ntens),time(2),dtime
      real(8), intent(in) :: temp,dtemp,predef(*),dpred(*)
      real(8), intent(in) :: props(nprops),coords(3),drot(3,3)
      real(8), intent(in) :: dfgrd0(3,3),dfgrd1(3,3),celent
      real(8), intent(inout) :: stress(ntens),statev(nstatv)
      real(8), intent(inout) :: pnewdt
      real(8), intent(out) :: ddsdde(ntens,ntens)
      real(8), intent(out) :: sse,spd,scd,rpl
      real(8), intent(out) :: ddsddt(ntens),drplde(ntens),drpldt

      !--- Local variables
      integer :: ns, is, i, j, cpconv
      real(8) :: phi1, PHI, phi2, caratio, mattemp
      real(8) :: C11,C12,C13,C33,C44,C66
      real(8) :: Cc(6,6), Cs(6,6), rot4(6,6)
      real(8) :: gmatinv_0(3,3), gmatinv_t(3,3), gmatinv(3,3)
      real(8) :: dirc(MAXNS,3), norc(MAXNS,3)
      real(8) :: dirs_t(MAXNS,3), nors_t(MAXNS,3)
      real(8) :: Schmid(MAXNS,3,3), Schmidvec(MAXNS,6)
      real(8) :: SxS(MAXNS,6,6)
      real(8) :: burgerv(MAXNS), tauc0(MAXNS)
      real(8) :: rho0, gf
      real(8) :: gammasum_t(MAXNS), tauc_t(MAXNS)
      real(8) :: ssd_t(MAXNS), gnd_t(MAXNS)
      real(8) :: Fp_t(3,3), Fp(3,3)
      real(8) :: sigma_t(6), sigma(6), jacobi(6,6)
      real(8) :: evmp_t, evmp
      real(8) :: gammasum(MAXNS), tauc(MAXNS)
      real(8) :: ssd(MAXNS), gnd(MAXNS)
      real(8) :: slipparam(10), hparam(10)
      integer :: sv0, gndmodel

      !--- Read PROPS
      phi1    = props(1);  PHI     = props(2);  phi2    = props(3)
      ns      = int(props(4))
      caratio = props(5)
      if (props(6) < 0d0) then;  mattemp = temp
      else;                       mattemp = props(6);  end if
      C11 = props(7);  C12 = props(8);  C13 = props(9)
      C33 = props(10); C44 = props(11); C66 = (C11-C12)/2d0
      slipparam(1)  = props(12)  ! alpha0
      slipparam(2)  = props(13)  ! beta0
      slipparam(3)  = props(14)  ! rhom
      slipparam(4)  = props(15)  ! DeltaF
      slipparam(5)  = props(16)  ! nu0
      slipparam(6)  = props(17)  ! AV0
      slipparam(7)  = props(18)  ! gamma0
      hparam(1)     = props(19)  ! k1
      hparam(2)     = props(20)  ! k2_const
      hparam(3)     = props(21)  ! X_km
      hparam(4)     = props(22)  ! g_km
      hparam(5)     = props(23)  ! D_km
      hparam(6)     = props(24)  ! pdot0
      rho0          = props(29)
      gf            = props(30)
      gndmodel      = int(props(31))

      ! Burgers vectors and initial CRSS
      burgerv = 0d0;  tauc0 = 0d0
      call hcp_set_burgers_crss(ns, caratio, props(25), props(26),
     +   props(27), props(28), burgerv, tauc0)

      !--- Track global element/IP count for UEXTERNALDB
      if (noel > numel) numel = noel
      if (npt  > numpt) numpt = npt
      g_nslip = ns

      !--- Elasticity tensor (HCP crystal frame, Voigt 6x6)
      call hcp_elastic_tensor(C11,C12,C13,C33,C44,C66,Cc)

      !--- Orientation matrix from Euler angles (crystal->sample transpose)
      call euler2gmat(phi1, PHI, phi2, gmatinv_0)

      !--- Recover STATEV
      sv0 = 0
      gammasum_t(1:ns) = statev(sv0+1:sv0+ns);    sv0=sv0+ns
      tauc_t(1:ns)     = statev(sv0+1:sv0+ns);    sv0=sv0+ns
      ssd_t(1:ns)      = statev(sv0+1:sv0+ns);    sv0=sv0+ns
      gnd_t(1:ns)      = statev(sv0+1:sv0+ns);    sv0=sv0+ns
      Fp_t(1,1)=statev(sv0+1); Fp_t(1,2)=statev(sv0+2)
      Fp_t(1,3)=statev(sv0+3); Fp_t(2,1)=statev(sv0+4)
      Fp_t(2,2)=statev(sv0+5); Fp_t(2,3)=statev(sv0+6)
      Fp_t(3,1)=statev(sv0+7); Fp_t(3,2)=statev(sv0+8)
      Fp_t(3,3)=statev(sv0+9); sv0=sv0+9
      gmatinv_t(1,1)=statev(sv0+1); gmatinv_t(1,2)=statev(sv0+2)
      gmatinv_t(1,3)=statev(sv0+3); gmatinv_t(2,1)=statev(sv0+4)
      gmatinv_t(2,2)=statev(sv0+5); gmatinv_t(2,3)=statev(sv0+6)
      gmatinv_t(3,1)=statev(sv0+7); gmatinv_t(3,2)=statev(sv0+8)
      gmatinv_t(3,3)=statev(sv0+9); sv0=sv0+9
      sigma_t(1:6) = statev(sv0+1:sv0+6);         sv0=sv0+6
      evmp_t       = statev(sv0+1)

      !--- First-increment initialization (guard: glob_ipinit allocated by UEXTERNALDB LOP=2/KINC=1)
      if (.not. allocated(glob_ipinit) .or. glob_ipinit(noel,npt) == 0) then
        gmatinv_t = gmatinv_0
        Fp_t = reshape([1d0,0d0,0d0, 0d0,1d0,0d0, 0d0,0d0,1d0],[3,3])
        tauc_t(1:ns) = tauc0(1:ns)
        ssd_t(1:ns)  = rho0
        gnd_t(1:ns)  = 0d0
        gammasum_t   = 0d0;  sigma_t = 0d0;  evmp_t = 0d0
        glob_ipinit(noel,npt) = 1
      end if

      !--- Read GND from global array (updated by UEXTERNALDB after previous step)
      if (gndmodel == 1 .and. allocated(glob_gnd)) then
        if (noel <= size(glob_gnd,1) .and. npt <= size(glob_gnd,2)) then
          gnd_t(1:ns) = glob_gnd(noel,npt,1:ns)
        end if
      end if

      !--- Build slip system vectors in crystal frame
      call hcp_init_slipvectors(ns, caratio, dirc, norc)

      !--- Rotate slip systems to sample frame (using gmatinv_t)
      call rotate_slip_systems(ns, gmatinv_t, dirc, norc, dirs_t, nors_t)

      !--- Build Schmid tensors
      call build_schmid(ns, dirs_t, nors_t, Schmid, Schmidvec, SxS)

      !--- Store slip vectors in global array for UEXTERNALDB
      if (allocated(glob_dirc)) then
        glob_dirc(noel,npt,1:ns,:) = dirs_t(1:ns,:)
        glob_norc(noel,npt,1:ns,:) = nors_t(1:ns,:)
        glob_burgerv(noel,npt,1:ns) = burgerv(1:ns)
        glob_coords(noel,npt,:) = coords
      end if

      !--- Rotate elasticity to sample frame
      call rotord4sig(gmatinv_t, rot4)
      Cs = matmul(rot4, matmul(Cc, transpose(rot4)))

      !--- CP solver
      call hcp_solve(noel, npt, dfgrd0, dfgrd1, dtime,
     +   ns, caratio, mattemp, Cs, gf, burgerv,
     +   tauc0, slipparam, hparam,
     +   Fp_t, gmatinv_t, sigma_t, evmp_t,
     +   gammasum_t, tauc_t, ssd_t, gnd_t,
     +   Schmid, Schmidvec, SxS,
     +   sigma, jacobi, Fp, gmatinv,
     +   gammasum, tauc, ssd, evmp,
     +   cpconv)

      !--- Cutback on non-convergence
      if (cpconv == 0) then
        pnewdt = CUTBACK
        return
      end if

      !--- Update global Fp for GND calculation
      if (allocated(glob_Fp)) glob_Fp(noel,npt,:,:) = Fp

      !--- Store STATEV
      sv0 = 0
      statev(sv0+1:sv0+ns) = gammasum(1:ns); sv0=sv0+ns
      statev(sv0+1:sv0+ns) = tauc(1:ns);     sv0=sv0+ns
      statev(sv0+1:sv0+ns) = ssd(1:ns);      sv0=sv0+ns
      statev(sv0+1:sv0+ns) = gnd_t(1:ns);    sv0=sv0+ns  ! gnd updated by UEXTERNALDB
      statev(sv0+1)=Fp(1,1); statev(sv0+2)=Fp(1,2); statev(sv0+3)=Fp(1,3)
      statev(sv0+4)=Fp(2,1); statev(sv0+5)=Fp(2,2); statev(sv0+6)=Fp(2,3)
      statev(sv0+7)=Fp(3,1); statev(sv0+8)=Fp(3,2); statev(sv0+9)=Fp(3,3)
      sv0=sv0+9
      statev(sv0+1)=gmatinv(1,1); statev(sv0+2)=gmatinv(1,2)
      statev(sv0+3)=gmatinv(1,3); statev(sv0+4)=gmatinv(2,1)
      statev(sv0+5)=gmatinv(2,2); statev(sv0+6)=gmatinv(2,3)
      statev(sv0+7)=gmatinv(3,1); statev(sv0+8)=gmatinv(3,2)
      statev(sv0+9)=gmatinv(3,3); sv0=sv0+9
      statev(sv0+1:sv0+6) = sigma(1:6); sv0=sv0+6
      statev(sv0+1) = evmp

      !--- Return stress and consistent tangent
      stress(1:6) = sigma(1:6)
      ddsdde(1:6,1:6) = jacobi(1:6,1:6)

      return
      end subroutine umat

!=============================================================================
! SUBROUTINE UEXTERNALDB  -  Initialization and GND update
!=============================================================================
      subroutine uexternaldb(lop,lrestart,time,dtime,kstep,kinc)

      use hcp_constants
      use hcp_globals
      implicit none
      integer, intent(in) :: lop,lrestart,kstep,kinc
      real(8), intent(in) :: time(2),dtime

      integer :: iel, ipt, is, ns
      integer :: maxnel_est, maxnpt_est
      real(8) :: Fp_avg(3,3), dFp(3,3), L_char
      real(8) :: coords_avg(3), dcrd(3), spread_sq
      real(8) :: sn(3), nn(3), b
      real(8) :: proj

      !=========================================================
      ! LOP=0: Called once at analysis start - allocate arrays
      !=========================================================
      if (lop == 0) then
        g_initialized = .false.
        numel = 0;  numpt = 0;  g_nslip = 0
        return
      end if

      !=========================================================
      ! LOP=1: Start of step (before UMAT calls) - no-op here
      !=========================================================
      if (lop == 1) then
        return
      end if

      !=========================================================
      ! LOP=2: Called after each increment.
      !   KINC=1: allocate global arrays (UMAT has now set numel/numpt)
      !   KINC>1: compute GND from element-average Fp variation
      !=========================================================
      if (lop == 2) then

        !--- Allocate on first completed increment
        if (.not. g_initialized .and. numel > 0 .and. numpt > 0) then
          maxnel_est = numel + 100
          maxnpt_est = numpt
          allocate(glob_Fp(maxnel_est, maxnpt_est, 3, 3))
          allocate(glob_gnd(maxnel_est, maxnpt_est, MAXNS))
          allocate(glob_dirc(maxnel_est, maxnpt_est, MAXNS, 3))
          allocate(glob_norc(maxnel_est, maxnpt_est, MAXNS, 3))
          allocate(glob_burgerv(maxnel_est, maxnpt_est, MAXNS))
          allocate(glob_coords(maxnel_est, maxnpt_est, 3))
          allocate(glob_ipinit(maxnel_est, maxnpt_est))
          glob_Fp = 0d0;  glob_gnd = 0d0;  glob_dirc = 0d0
          glob_norc = 0d0;  glob_burgerv = 0d0;  glob_coords = 0d0
          glob_ipinit = 0
          do iel = 1, maxnel_est
            do ipt = 1, maxnpt_est
              glob_Fp(iel,ipt,1,1) = 1d0
              glob_Fp(iel,ipt,2,2) = 1d0
              glob_Fp(iel,ipt,3,3) = 1d0
            end do
          end do
          g_initialized = .true.
          write(*,*) '[HCP_CPUMAT] Allocated: numel=',numel,' numpt=',numpt
          return
        end if

        if (.not. g_initialized) return
        if (g_nslip == 0) return

        ns = g_nslip

        ! Loop over all elements
        do iel = 1, numel

          ! Compute element-average Fp
          Fp_avg = 0d0
          coords_avg = 0d0
          do ipt = 1, numpt
            Fp_avg     = Fp_avg     + glob_Fp(iel,ipt,:,:)
            coords_avg = coords_avg + glob_coords(iel,ipt,:)
          end do
          Fp_avg     = Fp_avg     / dble(numpt)
          coords_avg = coords_avg / dble(numpt)

          ! Characteristic length: RMS spread of IP coordinates
          spread_sq = 0d0
          do ipt = 1, numpt
            dcrd = glob_coords(iel,ipt,:) - coords_avg
            spread_sq = spread_sq + dot_product(dcrd,dcrd)
          end do
          L_char = sqrt(spread_sq / dble(numpt))
          if (L_char < 1d-12) L_char = 1d0   ! fallback 1 um

          ! GND at each IP from Fp deviation
          do ipt = 1, numpt
            dFp = glob_Fp(iel,ipt,:,:) - Fp_avg

            do is = 1, ns
              b  = glob_burgerv(iel,ipt,is)
              sn = glob_dirc(iel,ipt,is,:)
              nn = glob_norc(iel,ipt,is,:)
              if (b < 1d-20) cycle

              ! Project dFp onto slip system Schmid tensor: tr(dFp . (m x n))
              proj = sum(dFp * spread(sn,2,3) * spread(nn,1,3))
              glob_gnd(iel,ipt,is) = abs(proj) / (b * L_char)
            end do
          end do

        end do  ! element loop
        return
      end if

      !=========================================================
      ! LOP=3: End of analysis - deallocate
      !=========================================================
      if (lop == 3) then
        if (allocated(glob_Fp))      deallocate(glob_Fp)
        if (allocated(glob_gnd))     deallocate(glob_gnd)
        if (allocated(glob_dirc))    deallocate(glob_dirc)
        if (allocated(glob_norc))    deallocate(glob_norc)
        if (allocated(glob_burgerv)) deallocate(glob_burgerv)
        if (allocated(glob_coords))  deallocate(glob_coords)
        if (allocated(glob_ipinit))  deallocate(glob_ipinit)
        g_initialized = .false.
        return
      end if

      return
      end subroutine uexternaldb

!=============================================================================
! SUBROUTINE hcp_solve  -  Main CP solver (predictor + NR loop)
!=============================================================================
      subroutine hcp_solve(noel, npt, dfgrd0, dfgrd1, dt,
     +   ns, caratio, mattemp, Cs, gf, burgerv,
     +   tauc0, slipparam, hparam,
     +   Fp_t, gmatinv_t, sigma_t, evmp_t,
     +   gammasum_t, tauc_t, ssd_t, gnd_t,
     +   Schmid, Schmidvec, SxS,
     +   sigma, jacobi, Fp, gmatinv,
     +   gammasum, tauc, ssd, evmp,
     +   cpconv)

      use hcp_constants
      use hcp_globals, only: g_nslip
      implicit none
      integer, intent(in) :: noel, npt, ns
      real(8), intent(in) :: dfgrd0(3,3), dfgrd1(3,3), dt
      real(8), intent(in) :: caratio, mattemp, gf
      real(8), intent(in) :: Cs(6,6), burgerv(MAXNS), tauc0(MAXNS)
      real(8), intent(in) :: slipparam(10), hparam(10)
      real(8), intent(in) :: Fp_t(3,3), gmatinv_t(3,3), sigma_t(6)
      real(8), intent(in) :: evmp_t
      real(8), intent(in) :: gammasum_t(MAXNS), tauc_t(MAXNS)
      real(8), intent(in) :: ssd_t(MAXNS), gnd_t(MAXNS)
      real(8), intent(in) :: Schmid(MAXNS,3,3), Schmidvec(MAXNS,6)
      real(8), intent(in) :: SxS(MAXNS,6,6)
      real(8), intent(out) :: sigma(6), jacobi(6,6)
      real(8), intent(out) :: Fp(3,3), gmatinv(3,3)
      real(8), intent(out) :: gammasum(MAXNS), tauc(MAXNS)
      real(8), intent(out) :: ssd(MAXNS), evmp
      integer, intent(out) :: cpconv

      ! Local
      real(8) :: F(3,3), F_t(3,3), Finv(3,3), detF
      real(8) :: Fdot(3,3), L(3,3), W(3,3)
      real(8) :: dstran(6), dstran33(3,3)
      real(8) :: sigmatr(6), sigma0(6)
      real(8) :: tau(MAXNS), tautr(MAXNS)
      real(8) :: tauceff(MAXNS), rhofor(MAXNS)
      real(8) :: rhotot(MAXNS)
      real(8) :: Lp(3,3), Dp(3,3), dstranp33(3,3)
      real(8) :: invJ(6,6)
      real(8) :: gammadot(MAXNS), pdot
      real(8) :: dRot(3,3), We(3,3), Le(3,3)
      real(8) :: Fe(3,3), invFp(3,3), detFp
      real(8) :: dtauc(MAXNS), dssd(MAXNS)
      real(8) :: sigma33(3,3), dummy6(6)
      integer :: is, iter, oiter
      real(8) :: maxx, MAXXCR = 1.0d-4

      g_nslip = ns
      cpconv  = 1

      !--- Deformation gradient (remove thermal/residual - simplified here)
      F   = dfgrd1;  F_t = dfgrd0

      !--- Velocity gradient
      Fdot = (F - F_t) / dt
      call inv3x3(F, Finv, detF)
      if (detF < SMALL) then;  cpconv=0;  return;  end if
      L = matmul(Fdot, Finv)

      !--- Strain increment
      dstran33 = (L + transpose(L)) * 0.5d0 * dt
      call matvec6(dstran33, dstran)
      dstran(4:6) = 2d0 * dstran(4:6)

      !--- Spin
      W = (L - transpose(L)) * 0.5d0

      !--- Trial stress
      sigmatr = sigma_t + matmul(Cs, dstran)

      !--- RSS on trial state
      do is = 1, ns
        tautr(is) = dot_product(Schmidvec(is,1:6), sigmatr)
      end do

      !--- Effective CRSS (Taylor law: tauc0 + G*b*gf*sqrt(rhofor))
      rhotot(1:ns) = ssd_t(1:ns) + gnd_t(1:ns)
      do is = 1, ns
        rhofor(is) = max(rhotot(is), SMALL)
      end do
      call hcp_crss(ns, gf, burgerv, tauc_t, rhofor, tauceff)

      !--- Check if plasticity is needed
      maxx = maxval(abs(tautr(1:ns)) / tauceff(1:ns))
      if (maxx <= MAXXCR) then
        sigma  = sigmatr
        jacobi = Cs
        Fp     = Fp_t
        call euler_spin_update(gmatinv_t, W, dt, gmatinv)
        gammasum = gammasum_t;  tauc = tauc_t;  ssd = ssd_t;  evmp = evmp_t
        return
      end if

      !--- Predictor: weighted average of elastic and former step
      sigma0 = 0.5d0*sigma_t + 0.5d0*sigmatr

      !--- Outer loop: update state variables (semi-implicit)
      sigma = sigma0
      do is = 1, ns
        tau(is) = dot_product(Schmidvec(is,1:6), sigma)
      end do
      gammasum = gammasum_t;  tauc = tauc_t;  ssd = ssd_t

      oiter = 0
      do while (oiter < MAXITER .and. cpconv == 1)

        !--- Inner NR loop
        call hcp_innerloop(ns, caratio, mattemp, Cs,
     +     burgerv, slipparam, dt, sigmatr,
     +     tauceff, Schmid, Schmidvec, SxS,
     +     sigma, tau, cpconv,
     +     gammadot, Lp, Dp, dstranp33, invJ, iter)

        if (iter >= MAXITER) then;  cpconv=0;  return;  end if
        if (any(gammadot/=gammadot)) then;  cpconv=0;  return;  end if

        !--- Equivalent plastic strain rate
        pdot = sqrt(2d0/3d0 * sum(Dp*Dp))

        !--- Update cumulative slip
        do is = 1, ns
          gammasum(is) = gammasum_t(is) + gammadot(is)*dt
        end do
        evmp = evmp_t + pdot*dt

        !--- Hardening (Kocks-Mecking)
        call hcp_km_hardening(ns, mattemp, dt, burgerv, hparam,
     +     gammadot, pdot, rhofor, ssd_t,
     +     dtauc, dssd)

        tauc(1:ns) = tauc_t(1:ns) + dtauc(1:ns)
        ssd(1:ns)  = ssd_t(1:ns)  + dssd(1:ns)
        ssd(1:ns)  = max(ssd(1:ns), 0d0)

        !--- Recompute forest density and CRSS
        rhotot(1:ns) = ssd(1:ns) + gnd_t(1:ns)
        do is = 1, ns
          rhofor(is) = max(rhotot(is), SMALL)
        end do
        call hcp_crss(ns, gf, burgerv, tauc, rhofor, tauceff)

        if (any(tauc < 0d0)) then;  cpconv=0;  return;  end if
        oiter = oiter + 1
        exit  ! Single outer iteration (explicit state update); change for semi-implicit
      end do

      !--- Consistent tangent
      jacobi = matmul(invJ, Cs)
      if (any(jacobi/=jacobi)) then;  cpconv=0;  return;  end if

      !--- Update plastic deformation gradient
      call inv3x3(F, Finv, detF)
      Fp = Fp_t + matmul(matmul(Finv, Lp), matmul(F, Fp_t))*dt
      call deter3x3(Fp, detFp)
      if (detFp <= SMALL) then;  cpconv=0;  return;  end if
      Fp = Fp / detFp**(1d0/3d0)

      !--- Orientation update via elastic spin
      call inv3x3(Fp, invFp, detFp)
      Fe = matmul(F, invFp)
      Le = L - Lp
      We = (Le - transpose(Le)) * 0.5d0
      dRot = reshape([1d0,0d0,0d0,0d0,1d0,0d0,0d0,0d0,1d0],[3,3]) + We*dt
      gmatinv = matmul(dRot, gmatinv_t)

      return
      end subroutine hcp_solve

!=============================================================================
! SUBROUTINE hcp_innerloop  -  Newton-Raphson stress update
!=============================================================================
      subroutine hcp_innerloop(ns, caratio, mattemp, Cs,
     +   burgerv, slipparam, dt, sigmatr,
     +   tauceff, Schmid, Schmidvec, SxS,
     +   sigma, tau, cpconv,
     +   gammadot, Lp, Dp, dstranp33, invJ, iter)

      use hcp_constants
      implicit none
      integer, intent(in)    :: ns
      real(8), intent(in)    :: caratio, mattemp, dt
      real(8), intent(in)    :: Cs(6,6), burgerv(MAXNS), slipparam(10)
      real(8), intent(in)    :: tauceff(MAXNS), sigmatr(6)
      real(8), intent(in)    :: Schmid(MAXNS,3,3), Schmidvec(MAXNS,6)
      real(8), intent(in)    :: SxS(MAXNS,6,6)
      real(8), intent(inout) :: sigma(6), tau(MAXNS)
      integer, intent(inout) :: cpconv
      real(8), intent(out)   :: gammadot(MAXNS), Lp(3,3), Dp(3,3)
      real(8), intent(out)   :: dstranp33(3,3), invJ(6,6)
      integer, intent(out)   :: iter

      real(8) :: Pmat(6,6), dstranp(6)
      real(8) :: J(6,6), psi(6), psinorm, dsigma(6)
      real(8) :: dgdt(MAXNS), dgdtc(MAXNS)
      integer :: is, err
      real(8), parameter :: I6(6,6) = reshape([&
        1d0,0d0,0d0,0d0,0d0,0d0, 0d0,1d0,0d0,0d0,0d0,0d0,&
        0d0,0d0,1d0,0d0,0d0,0d0, 0d0,0d0,0d0,1d0,0d0,0d0,&
        0d0,0d0,0d0,0d0,1d0,0d0, 0d0,0d0,0d0,0d0,0d0,1d0],[6,6])

      psinorm = 1d0;  iter = 0

      do while (psinorm >= TOL_NR .and. iter < MAXITER .and. cpconv==1)

        !--- Sinh slip law
        call hcp_sinhslip(ns, caratio, mattemp, burgerv, slipparam,
     +     tau, tauceff, Schmid, SxS, dt,
     +     gammadot, dgdt, dgdtc, Lp, Dp, Pmat)

        !--- Plastic strain increment
        dstranp33 = Dp * dt
        call matvec6(dstranp33, dstranp)
        dstranp(4:6) = 2d0 * dstranp(4:6)

        !--- NR Jacobian: J = I + Cs * Pmat
        J = I6 + matmul(Cs, Pmat)

        !--- Invert J
        call nolapinverse(J, invJ, 6, err)
        if (err /= 0) then;  cpconv=0;  return;  end if

        !--- Residual: psi = sigmatr - sigma - Cs*dstranp
        psi = sigmatr - sigma - matmul(Cs, dstranp)
        psinorm = sqrt(sum(psi*psi))

        !--- Stress correction
        dsigma = matmul(invJ, psi)
        sigma  = sigma + dsigma

        !--- Update RSS
        do is = 1, ns
          tau(is) = dot_product(Schmidvec(is,1:6), sigma)
        end do

        iter = iter + 1
      end do

      return
      end subroutine hcp_innerloop

!=============================================================================
! SUBROUTINE hcp_sinhslip  -  Sinh slip law for HCP
!   gamma_dot = alpha * sinh(beta * (|tau| - tauc)) * sign(tau)
!   if |tau| >= tauc, else 0
!=============================================================================
      subroutine hcp_sinhslip(ns, caratio, T, burgerv, sparam,
     +   tau, tauceff, Schmid, SxS, dt,
     +   gammadot, dgdt, dgdtc, Lp, Dp, Pmat)

      use hcp_constants
      implicit none
      integer, intent(in) :: ns
      real(8), intent(in) :: caratio, T, dt
      real(8), intent(in) :: burgerv(MAXNS), sparam(10)
      real(8), intent(in) :: tau(MAXNS), tauceff(MAXNS)
      real(8), intent(in) :: Schmid(MAXNS,3,3), SxS(MAXNS,6,6)
      real(8), intent(out) :: gammadot(MAXNS), dgdt(MAXNS), dgdtc(MAXNS)
      real(8), intent(out) :: Lp(3,3), Dp(3,3), Pmat(6,6)

      real(8) :: alpha0, beta0, rhom, DeltaF, nu0, AV0, gamma0
      real(8) :: alpha, beta, AV, lambda
      real(8) :: abstau, sgntau, drv
      integer :: is

      alpha0  = sparam(1);  beta0  = sparam(2)
      rhom    = sparam(3);  DeltaF = sparam(4);  nu0 = sparam(5)
      AV0     = sparam(6);  gamma0 = sparam(7)

      gammadot=0d0; dgdt=0d0; dgdtc=0d0
      Lp=0d0; Dp=0d0; Pmat=0d0

      do is = 1, ns
        abstau = abs(tau(is))
        sgntau = sign(1d0, tau(is))

        !--- Pre-factor alpha
        if (alpha0 == 0d0) then
          alpha = rhom * burgerv(is)**2 * nu0 *
     +            exp(-DeltaF / (KB * T))
        else
          alpha = alpha0
        end if

        !--- Activation volume and beta
        if (AV0 == 0d0) then
          if (tauceff(is) > SMALL) then
            lambda = 1d0 / sqrt(max(tauceff(is)/1d0, SMALL))  ! rough estimate
          else
            lambda = 1d0
          end if
          AV = gamma0 * lambda * burgerv(is)**2
        else
          AV = AV0 * burgerv(is)**3
        end if
        if (beta0 == 0d0) then
          beta = AV / (KB * T)
        else
          beta = beta0
        end if

        !--- c/a correction for Pyramidal-2 systems (is > 12)
        if (is > 12) then
          alpha = caratio*caratio * alpha
          beta  = caratio*caratio * beta
        end if

        !--- Threshold condition
        if (abstau >= tauceff(is)) then
          drv = beta * (abstau - tauceff(is))
          gammadot(is) = alpha * sinh(drv) * sgntau
          dgdt(is)     = alpha * beta * cosh(drv)
          dgdtc(is)    = -alpha * beta * cosh(drv) * sgntau
          Pmat = Pmat + dt * dgdt(is) * SxS(is,:,:)
          Lp   = Lp   + gammadot(is) * Schmid(is,:,:)
        end if
      end do

      Dp   = (Lp + transpose(Lp)) * 0.5d0
      Pmat = (Pmat + transpose(Pmat)) * 0.5d0

      return
      end subroutine hcp_sinhslip

!=============================================================================
! SUBROUTINE hcp_crss  -  Taylor dislocation hardening
!   tauc_eff = tauc0 + gf * b * sqrt(rhofor)
!
!   NOTE: 'gf' here is the combined factor gf_dimensioned = alpha * G  [MPa/um]
!         so the user sets PROPS(30) = alpha * G12  (e.g., 0.3 * 32000 = 9600 MPa/um)
!         and burgerv in [um], rhofor in [1/um^2].
!=============================================================================
      subroutine hcp_crss(ns, gf, burgerv, tauc0, rhofor, tauceff)
      use hcp_constants
      implicit none
      integer, intent(in) :: ns
      real(8), intent(in) :: gf, burgerv(MAXNS), tauc0(MAXNS)
      real(8), intent(in) :: rhofor(MAXNS)
      real(8), intent(out) :: tauceff(MAXNS)
      integer :: is
      do is = 1, ns
        tauceff(is) = tauc0(is) + gf * burgerv(is) * sqrt(max(rhofor(is),0d0))
      end do
      return
      end subroutine hcp_crss

!=============================================================================
! SUBROUTINE hcp_km_hardening  -  Kocks-Mecking SSD evolution
!   drho = (k1/b * sqrt(rhofor) - k2 * rho) * |gammadot| * dt
!=============================================================================
      subroutine hcp_km_hardening(ns, T, dt, burgerv, hparam,
     +   gammadot, pdot, rhofor, ssd_t, dtauc, dssd)

      use hcp_constants
      implicit none
      integer, intent(in) :: ns
      real(8), intent(in) :: T, dt, pdot
      real(8), intent(in) :: burgerv(MAXNS), hparam(10)
      real(8), intent(in) :: gammadot(MAXNS), rhofor(MAXNS), ssd_t(MAXNS)
      real(8), intent(out) :: dtauc(MAXNS), dssd(MAXNS)

      real(8) :: k1, k2_const, X_km, g_km, D_km, pdot0
      real(8) :: k2, b, KBT
      integer :: is

      k1       = hparam(1);  k2_const = hparam(2);  X_km = hparam(3)
      g_km     = hparam(4);  D_km     = hparam(5);  pdot0 = hparam(6)

      KBT = KB * T * 1d12  ! J -> MPa*um^3 (1J = 10^12 MPa*um^3)

      dtauc = 0d0;  dssd = 0d0

      do is = 1, ns
        b = burgerv(is)
        if (b < SMALL) cycle

        if (k2_const == 0d0) then
          if (pdot > SMALL .and. pdot0 > SMALL .and. D_km > SMALL) then
            k2 = k1*X_km*b/g_km * (1d0 - KBT/(D_km*b**3)*log(pdot/pdot0))
          else
            k2 = k1*X_km*b/g_km
          end if
          k2 = max(k2, 0d0)
        else
          k2 = k2_const
        end if

        dssd(is) = (k1/b * sqrt(max(rhofor(is),0d0)) -
     +              k2 * ssd_t(is)) * abs(gammadot(is)) * dt
      end do
      ! dtauc updated implicitly via ssd -> rhofor -> tauceff
      dtauc = 0d0

      return
      end subroutine hcp_km_hardening

!=============================================================================
! SUBROUTINE hcp_init_slipvectors  -  HCP slip directions/normals (Cartesian)
!   Miller-Bravais [uvtw] -> Cartesian using c/a ratio
!=============================================================================
      subroutine hcp_init_slipvectors(ns, caratio, dirc, norc)
      use hcp_slip_data
      use hcp_constants
      implicit none
      integer, intent(in)  :: ns
      real(8), intent(in)  :: caratio
      real(8), intent(out) :: dirc(MAXNS,3), norc(MAXNS,3)

      integer :: is
      real(8) :: d(3), n(3), dmag, nmag

      dirc=0d0;  norc=0d0

      do is = 1, min(ns, 30)
        ! [uvtw] -> orthogonal: x=3u/2, y=(u+2v)*sqrt(3)/2, z=w*(c/a)
        d(1) = 3d0 * DIR3H(is,1) / 2d0
        d(2) = (DIR3H(is,1) + 2d0*DIR3H(is,2)) * sqrt(3d0) / 2d0
        d(3) = DIR3H(is,4) * caratio

        ! (hkil) -> orthogonal: x=h, y=(h+2k)/sqrt(3), z=l/(c/a)
        n(1) = NOR3H(is,1)
        n(2) = (NOR3H(is,1) + 2d0*NOR3H(is,2)) / sqrt(3d0)
        n(3) = NOR3H(is,4) / caratio

        dmag = sqrt(sum(d*d));  if (dmag < SMALL) dmag=1d0
        nmag = sqrt(sum(n*n));  if (nmag < SMALL) nmag=1d0

        dirc(is,:) = d / dmag
        norc(is,:) = n / nmag
      end do

      return
      end subroutine hcp_init_slipvectors

!=============================================================================
! SUBROUTINE rotate_slip_systems  -  Rotate crystal-frame slip vectors to sample
!=============================================================================
      subroutine rotate_slip_systems(ns, gmatinv, dirc, norc, dirs, nors)
      use hcp_constants
      implicit none
      integer, intent(in)  :: ns
      real(8), intent(in)  :: gmatinv(3,3), dirc(MAXNS,3), norc(MAXNS,3)
      real(8), intent(out) :: dirs(MAXNS,3), nors(MAXNS,3)
      integer :: is
      real(8) :: td(3), tn(3), dm, nm

      dirs=0d0;  nors=0d0
      do is = 1, ns
        td = matmul(gmatinv, dirc(is,:))
        tn = matmul(gmatinv, norc(is,:))
        dm = sqrt(sum(td*td));  if (dm < SMALL) dm=1d0
        nm = sqrt(sum(tn*tn));  if (nm < SMALL) nm=1d0
        dirs(is,:) = td / dm
        nors(is,:) = tn / nm
      end do
      return
      end subroutine rotate_slip_systems

!=============================================================================
! SUBROUTINE build_schmid  -  Schmid tensors from slip vectors
!=============================================================================
      subroutine build_schmid(ns, dirs, nors, Schmid, Schmidvec, SxS)
      use hcp_constants
      implicit none
      integer, intent(in)  :: ns
      real(8), intent(in)  :: dirs(MAXNS,3), nors(MAXNS,3)
      real(8), intent(out) :: Schmid(MAXNS,3,3), Schmidvec(MAXNS,6), SxS(MAXNS,6,6)

      integer :: is, i, j
      real(8) :: sn(3), nn(3), SNij(3,3), NSij(3,3), sni(6), nsi(6)

      Schmid=0d0; Schmidvec=0d0; SxS=0d0

      do is = 1, ns
        sn = dirs(is,:);  nn = nors(is,:)
        do i=1,3
          do j=1,3
            SNij(i,j) = sn(i)*nn(j)
            NSij(j,i) = nn(j)*sn(i)
            Schmid(is,i,j) = SNij(i,j)
          end do
        end do
        call gmatvec6(SNij, sni)
        call gmatvec6(NSij, nsi)
        Schmidvec(is,:) = sni
        do i=1,6;  do j=1,6
          SxS(is,i,j) = sni(i)*nsi(j)
        end do;  end do
      end do
      return
      end subroutine build_schmid

!=============================================================================
! SUBROUTINE hcp_set_burgers_crss  -  Assign burgers vector and initial CRSS
!=============================================================================
      subroutine hcp_set_burgers_crss(ns, caratio,
     +   tauc_basal, tauc_prism, tauc_pyr1, tauc_pyr2,
     +   burgerv, tauc0)
      use hcp_constants
      implicit none
      integer, intent(in)  :: ns
      real(8), intent(in)  :: caratio
      real(8), intent(in)  :: tauc_basal, tauc_prism, tauc_pyr1, tauc_pyr2
      real(8), intent(out) :: burgerv(MAXNS), tauc0(MAXNS)

      real(8) :: b_a, b_ca
      ! Burgers magnitude: <a>-type = 1 (normalized), <c+a>-type = sqrt(1+3*(c/a)^2)
      ! In practice user supplies in physical units (um); here we use ratio.
      b_a  = 1d0                        ! normalized basal/prism/pyr1 Burgers
      b_ca = sqrt(1d0 + 3d0*caratio**2) ! normalized pyr2 Burgers

      burgerv = 0d0;  tauc0 = 0d0

      ! Basal 1-3
      burgerv(1:min(ns,3))  = b_a;   tauc0(1:min(ns,3))  = tauc_basal
      ! Prismatic 4-6
      if (ns >= 6) then
        burgerv(4:6)  = b_a;         tauc0(4:6)  = tauc_prism
      end if
      ! Pyramidal-1 7-12
      if (ns >= 12) then
        burgerv(7:12) = b_a;         tauc0(7:12) = tauc_pyr1
      end if
      ! Pyramidal-2 13+
      if (ns >= 13) then
        burgerv(13:ns) = b_ca;       tauc0(13:ns) = tauc_pyr2
      end if

      return
      end subroutine hcp_set_burgers_crss

!=============================================================================
! SUBROUTINE euler2gmat  -  Bunge Euler angles -> gmatinv (sample->crystal)
!=============================================================================
      subroutine euler2gmat(phi1_deg, Phi_deg, phi2_deg, gmatinv)
      use hcp_constants
      implicit none
      real(8), intent(in)  :: phi1_deg, Phi_deg, phi2_deg
      real(8), intent(out) :: gmatinv(3,3)
      real(8) :: p1, P, p2, g(3,3)

      p1 = phi1_deg * PI / 180d0
      P  = Phi_deg  * PI / 180d0
      p2 = phi2_deg * PI / 180d0

      g(1,1) =  cos(p1)*cos(p2) - sin(p1)*sin(p2)*cos(P)
      g(2,1) = -cos(p1)*sin(p2) - sin(p1)*cos(p2)*cos(P)
      g(3,1) =  sin(p1)*sin(P)
      g(1,2) =  sin(p1)*cos(p2) + cos(p1)*sin(p2)*cos(P)
      g(2,2) = -sin(p1)*sin(p2) + cos(p1)*cos(p2)*cos(P)
      g(3,2) = -cos(p1)*sin(P)
      g(1,3) =  sin(p2)*sin(P)
      g(2,3) =  cos(p2)*sin(P)
      g(3,3) =  cos(P)

      gmatinv = transpose(g)  ! sample->crystal = g^T
      return
      end subroutine euler2gmat

!=============================================================================
! SUBROUTINE euler_spin_update  -  Orientation update for elastic case
!=============================================================================
      subroutine euler_spin_update(gmatinv_t, W, dt, gmatinv)
      implicit none
      real(8), intent(in)  :: gmatinv_t(3,3), W(3,3), dt
      real(8), intent(out) :: gmatinv(3,3)
      real(8) :: domega(3), dW33(3,3)

      domega(1) = (W(1,2) - W(2,1)) * dt / 2d0
      domega(2) = (W(3,1) - W(1,3)) * dt / 2d0
      domega(3) = (W(2,3) - W(3,2)) * dt / 2d0

      dW33 = 0d0
      dW33(1,2) =  domega(1);  dW33(2,1) = -domega(1)
      dW33(1,3) = -domega(2);  dW33(3,1) =  domega(2)
      dW33(2,3) =  domega(3);  dW33(3,2) = -domega(3)

      gmatinv = gmatinv_t + matmul(dW33, gmatinv_t)
      return
      end subroutine euler_spin_update

!=============================================================================
! SUBROUTINE hcp_elastic_tensor  -  HCP stiffness matrix (Voigt 6x6, crystal frame)
!=============================================================================
      subroutine hcp_elastic_tensor(C11,C12,C13,C33,C44,C66,Cc)
      implicit none
      real(8), intent(in)  :: C11,C12,C13,C33,C44,C66
      real(8), intent(out) :: Cc(6,6)
      Cc = 0d0
      Cc(1,1)=C11; Cc(1,2)=C12; Cc(1,3)=C13
      Cc(2,1)=C12; Cc(2,2)=C11; Cc(2,3)=C13
      Cc(3,1)=C13; Cc(3,2)=C13; Cc(3,3)=C33
      Cc(4,4)=C44; Cc(5,5)=C44; Cc(6,6)=C66
      return
      end subroutine hcp_elastic_tensor

!=============================================================================
! SUBROUTINE rotord4sig  -  Bond transformation matrix for Voigt 6x6 rotation
!   Transforms stiffness: Cs = rot4 * Cc * rot4^T
!   Voigt index pairs: 1=(11),2=(22),3=(33),4=(12),5=(13),6=(23)
!=============================================================================
      subroutine rotord4sig(Q, rot4)
      implicit none
      real(8), intent(in)  :: Q(3,3)
      real(8), intent(out) :: rot4(6,6)
      integer :: i, j
      ! Voigt pair indices
      integer :: p(6,2)

      p(1,:) = [1,1];  p(2,:) = [2,2];  p(3,:) = [3,3]
      p(4,:) = [1,2];  p(5,:) = [1,3];  p(6,:) = [2,3]

      rot4 = 0d0
      do i = 1, 6
        do j = 1, 6
          ! Normal-normal terms
          if (i <= 3 .and. j <= 3) then
            rot4(i,j) = Q(p(i,1),p(j,1)) * Q(p(i,2),p(j,2))
          ! Normal-shear terms
          else if (i <= 3 .and. j > 3) then
            rot4(i,j) = Q(p(i,1),p(j,1))*Q(p(i,2),p(j,2)) +
     +                  Q(p(i,1),p(j,2))*Q(p(i,2),p(j,1))
          ! Shear-normal terms
          else if (i > 3 .and. j <= 3) then
            rot4(i,j) = Q(p(i,1),p(j,1)) * Q(p(i,2),p(j,2))
          ! Shear-shear terms
          else
            rot4(i,j) = Q(p(i,1),p(j,1))*Q(p(i,2),p(j,2)) +
     +                  Q(p(i,1),p(j,2))*Q(p(i,2),p(j,1))
          end if
        end do
      end do
      return
      end subroutine rotord4sig

!=============================================================================
! SUBROUTINE matvec6  -  Convert symmetric 3x3 matrix to 6-vector (Voigt)
!=============================================================================
      subroutine matvec6(A, v)
      implicit none
      real(8), intent(in)  :: A(3,3)
      real(8), intent(out) :: v(6)
      v(1)=A(1,1); v(2)=A(2,2); v(3)=A(3,3)
      v(4)=A(1,2); v(5)=A(1,3); v(6)=A(2,3)
      return
      end subroutine matvec6

!=============================================================================
! SUBROUTINE gmatvec6  -  Same as matvec6 (alias kept for clarity)
!=============================================================================
      subroutine gmatvec6(A, v)
      implicit none
      real(8), intent(in)  :: A(3,3)
      real(8), intent(out) :: v(6)
      call matvec6(A, v)
      return
      end subroutine gmatvec6

!=============================================================================
! SUBROUTINE inv3x3  -  Analytic 3x3 matrix inverse
!=============================================================================
      subroutine inv3x3(A, Ainv, det)
      implicit none
      real(8), intent(in)  :: A(3,3)
      real(8), intent(out) :: Ainv(3,3), det

      det = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2)) &
          - A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)) &
          + A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))

      if (abs(det) < 1d-15) then;  Ainv=0d0;  return;  end if

      Ainv(1,1) = (A(2,2)*A(3,3)-A(2,3)*A(3,2))/det
      Ainv(1,2) = (A(1,3)*A(3,2)-A(1,2)*A(3,3))/det
      Ainv(1,3) = (A(1,2)*A(2,3)-A(1,3)*A(2,2))/det
      Ainv(2,1) = (A(2,3)*A(3,1)-A(2,1)*A(3,3))/det
      Ainv(2,2) = (A(1,1)*A(3,3)-A(1,3)*A(3,1))/det
      Ainv(2,3) = (A(1,3)*A(2,1)-A(1,1)*A(2,3))/det
      Ainv(3,1) = (A(2,1)*A(3,2)-A(2,2)*A(3,1))/det
      Ainv(3,2) = (A(1,2)*A(3,1)-A(1,1)*A(3,2))/det
      Ainv(3,3) = (A(1,1)*A(2,2)-A(1,2)*A(2,1))/det
      return
      end subroutine inv3x3

!=============================================================================
! SUBROUTINE deter3x3  -  Determinant of 3x3 matrix
!=============================================================================
      subroutine deter3x3(A, det)
      implicit none
      real(8), intent(in)  :: A(3,3)
      real(8), intent(out) :: det
      det = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2)) &
          - A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)) &
          + A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))
      return
      end subroutine deter3x3

!=============================================================================
! SUBROUTINE nolapinverse  -  LU-based inversion of n x n matrix
!=============================================================================
      subroutine nolapinverse(A, Ainv, n, err)
      implicit none
      integer, intent(in)  :: n
      real(8), intent(in)  :: A(n,n)
      real(8), intent(out) :: Ainv(n,n)
      integer, intent(out) :: err

      real(8) :: L(n,n), U(n,n), x(n), b(n)
      real(8) :: y(n), s
      integer :: i, j, k

      err = 0
      L = 0d0; U = 0d0

      ! LU decomposition (Doolittle, no pivoting)
      do i = 1, n
        do j = i, n
          s = A(i,j)
          do k = 1, i-1;  s = s - L(i,k)*U(k,j);  end do
          U(i,j) = s
        end do
        if (abs(U(i,i)) < 1d-30) then;  err=1;  Ainv=0d0;  return;  end if
        do j = i+1, n
          s = A(j,i)
          do k = 1, i-1;  s = s - L(j,k)*U(k,i);  end do
          L(j,i) = s / U(i,i)
        end do
        L(i,i) = 1d0
      end do

      ! Solve for each column of the inverse
      do j = 1, n
        b = 0d0;  b(j) = 1d0
        ! Forward substitution Ly=b
        do i = 1, n
          y(i) = b(i)
          do k = 1, i-1;  y(i) = y(i) - L(i,k)*y(k);  end do
        end do
        ! Back substitution Ux=y
        do i = n, 1, -1
          x(i) = y(i)
          do k = i+1, n;  x(i) = x(i) - U(i,k)*x(k);  end do
          x(i) = x(i) / U(i,i)
        end do
        Ainv(:,j) = x
      end do

      return
      end subroutine nolapinverse
