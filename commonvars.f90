!-----------------------------------------------------------------------
!  Copyright 2017 Mikhail Osintcev
!  This file is part of the EMtool developed at NCSU
!-----------------------------------------------------------------------
! This module contains global variables and objects

module commonvars
  
  implicit none;

  ! How many auxiliary tasks may be used
  integer, parameter :: auxmaxcount = 20;

  ! For convergence checking global multiplier of all dimensions in points
  integer :: conv;
  
  !--------------CONSTANTS--------------
    real*8, parameter :: Pi = 3.141592653589793238462643d0;
    real*8, parameter  :: cc = 1.0d0;     !speed of light

  !----------------DOMAIN---------------
    real*8  :: ftime;
    real*8  :: xsize, ysize, zsize; ! domain parameters
    real*8  :: cntx, cnty, cntz;    ! center of the universe

    integer :: maxa,maxb,maxc
    
  !-------------DIPOLE SOURCE---------------
    real*8 :: a, v, p, q
    real*8 :: alphaaa;
    real*8 :: om;
    real*8 :: Qx, Qy, Qz;
    real*8, allocatable :: cDNu(:);
    real*8, allocatable :: cNu(:);

    real*8, allocatable :: MEJx(:,:,:), MEJy(:,:,:), MEJz(:,:,:), MHx(:,:,:), MHy(:,:,:), MHz(:,:,:);
    real*8, allocatable :: DMEJx(:,:,:,:), DMEJy(:,:,:,:), DMEJz(:,:,:,:), DMHx(:,:,:,:), DMHy(:,:,:,:), DMHz(:,:,:,:);
    integer :: source_mainalloc = 0;

    real*8, allocatable :: auxMEJx(:,:,:), auxMEJy(:,:,:), auxMEJz(:,:,:), auxMHx(:,:,:), auxMHy(:,:,:), auxMHz(:,:,:);
    real*8, allocatable :: auxDMEJx(:,:,:,:), auxDMEJy(:,:,:,:), auxDMEJz(:,:,:,:), auxDMHx(:,:,:,:), auxDMHy(:,:,:,:), auxDMHz(:,:,:,:);
    integer :: source_auxalloc = 0;
    
  !-------------FATAL ERRORS--------------
    integer :: fatalerr=0;
    ! 0 - No errors

    ! 100 - Error in the initialization of the problem, zero points
    ! 101 - Error in the initialization of the problem, zero differentials
    
    ! 201 - Trying to update the E point out of the mesh
    ! 202 - Trying to update the H point out of the mesh
    ! 203 - Incorrect use of Component in updateEpoint function
    ! 204 - Incorrect use of Component in updateHpoint function
    ! 205 - Illegal boundary condition for E in fillEboundaryPoint function
    ! 206 - Illegal boundary condition for H in fillEboundaryPoint function
    ! 207 - Incorrect use of Component in fillEboundaryPoint function
    ! 208 - Incorrect use of Component in fillHboundaryPoint function
    ! 209 - Incorrect use of Component in getCurrentEPoint function
    ! 210 - Incorrect Current Type is used in getCurrentEPoint function
    ! 211 - Incorrect use of Component in getCurrentEPoint function
    ! 212 - Incorrect Current Type is used in getCurrentEPoint function
    ! 213 - Incorrect use of Component in JEbuildPoint function
    ! 214 - Incorrect use of Component in JHbuildPoint function
    ! 215 - Incorrect Problem Type is used while problem initiation
    ! 216 - Incorrect Effective Currents Type in JEffbuild function

    ! 301 - Trying to update the E point out of the mesh in PML
    ! 302 - Trying to update the H point out of the mesh in PML
    ! 303 - Incorrect use of Component in updateEpoint function
    ! 304 - Incorrect use of Component in updateHpoint function
    ! 305 - Illegal boundary condition for E in fillEboundaryPointPML function
    ! 306 - Illegal boundary condition for H in fillHboundaryPointPML function
    ! 307 - Incorrect use of Component in fillEboundaryPointPML function
    ! 308 - Incorrect use of Component in fillHboundaryPointPML function
    ! 309 - Incorrect use of Component in getCurrentEPointPML function
    ! 310 - Incorrect Current Type is used in getCurrentEPointPML function
    ! 311 - Incorrect use of Component in getCurrentHPointPML function
    ! 312 - Incorrect Current Type is used in getCurrentHPointPML function
    ! 313 - PML can not get the currents from buffer

    ! 400 - Error in the initialization of the solution, zero differentials
    ! 401 - Incorrect solution type during the solution initialization
    ! 402 - Auxiliary problems size correction should be > 0
    ! 403 - There is no free timers more!

    ! 1001 - Trying to write to the file that isn't open

  !------------BUILD OPTIONS--------------
    integer :: checkupdated = 0; ! Check that all the points are updated once
    integer :: checkbounds = 0;    ! Chech that all points are in the array bounds
   
  !------------GLOBAL TASKS--------------
    integer :: globaltask = 0;
    ! 0 - Calculate fields propagation
    ! 1 - Check convergence of the scheme     

  !-----------PARALLELIZATION------------
    integer :: doparallel = 0;
      ! 0 - No parallel computing
      ! 1 - Parallel computing (don't forget -openmp key while compiling)
    integer*4 :: nprocs;
    integer*4, allocatable :: balance(:);

    !-------------TIMERS-------------------
    real*8 :: tstarts(0:35);
    real*8 :: tends(0:35);
    real*8 :: times(0:35);
    character(len=10) :: timestitles(0:35);

  !----------------OUTPUT----------------
    !Console output type
    integer :: propconsoleoutput = 2;
      ! 0 - Totally silent
      ! 1 - Print only start and finish of the modelling
      ! 2 - Print start, finish and the results of each time step (steps, error, error point, max value);
      ! 3 - Print start, finish and the results of each time step (steps, error, remaining time);  
    integer :: screen_erroroutputsteps = 1;
    integer :: file_erroroutputsteps = 1;
    integer :: file_pmlerrors = 1;
    integer :: file_writelocations = 1;
    integer :: file_auxproblemsoutput = 0;           ! number of aux problems
    integer :: file_auxproblemsoutputstep = 1;       ! step in aux problems output
    integer :: file_auxproblemswritelocations = 1;   ! write locations or not in aux problems
    integer :: file_compositedivergence = 1;         ! write divergence
    integer :: file_staticfieldmax = 1;              ! write static field
    integer :: file_mainfields = 0;                  ! write main problem field values
    integer :: file_mainpmlfields = 0;               ! write pml of main problem field values
    
    !File output type
    integer :: fileoutput = 1;
      ! 0 - No file output
      ! 1 - Print Error (err_value)

    type, public :: tProblemStarter
       integer :: ptype;
       integer :: new_id;
       integer :: schemetype;
         ! 0 - Classic Yee 2-2
       integer :: currentstype;
         ! 0 - Source
         ! 1 - Cutted Source
         ! 2 - Buffer
         ! 3 - Cutted Buffer
         ! 4 - Cutted Buffer with Poisson
       integer :: pmlcurrentstype;
       integer :: effcurrentstype;
         ! 0 - common
         ! 1 - with poisson
       integer :: sourcetype;
         ! 0 - Dipole
         ! 1 - SDipole
       integer :: Eboundarytype;
       integer :: Hboundarytype;
         ! -1 - no boundaries
         ! 0 - zeros
         ! 1 - analytical solution 
         ! 2 - scheme
         ! 3 - pml
         ! 4 - buffer
         ! 5 - sommerfeld abc
         ! 6 - higdon abc
         ! 7 - betz-mittra abc
         ! 8 - mur abc
         ! 9 - givoli-neta abc Ey_x=0 2th order in time, exact elsewhere
         ! 10 - givoli-neta abc Ey_x=0, unsplit pml elsewhere
         ! 11 - givoli-neta abc Ey_x=0 1th order in time, exact elsewhere
         ! 12 - givoli-neta abc Ey_x=0 2th order in time, big aux task elsewhere
         ! 13 - Hagstrom-Warburton abc Ey_x=0, with exact elsewhere  
       integer :: PML_thickness;
       real*8  :: PML_param;
       ! Auxiliary offset
       integer :: bufoffsetx, bufoffsety, bufoffsetz;
       integer :: maindx, maindy, maindz;
       ! Effective currents params
         integer :: Namu;
         integer :: Ndmu;
         real*8 :: auxmaxtime;
         real*8 :: sgm_s;
       !-----------MODELLING VALUES----------
       integer :: Nx;        ! number of point for x
       integer :: Ny;        ! number of point for y
       integer :: Nz;        ! number of point for z
       integer :: Nt;        ! number of point for t
       real*8  :: hhx, hhy, hhz, ht;   ! differentials
       ! Poisson parameters
       integer :: npx;
       integer :: npy;
       integer :: npz;
       ! ABC Parameters
       real*8 :: hig1, hig2;
       integer :: dtmode, dxmode, dxauxmode, dtauxmode;
       integer :: gnNum;
       integer :: hwNum;
       integer :: hwEvaE;
       integer :: hwEvaP;
       integer :: hwEvaEa;
       integer :: hwEvaPa;
    end type tProblemStarter

    type, public :: tSolutionStarter
       integer :: sid
       !-----------Main problem size------------
       integer :: Nx;  ! number of point for x
       integer :: Ny;  ! number of point for y
       integer :: Nz;  ! number of point for z
       !-----------Aux problems size------------
       integer :: dNxaux;  ! number of point for x
       integer :: dNyaux;  ! number of point for y
       integer :: dNzaux;  ! number of point for z
       ! Solution type
         integer :: soltype;
         ! 0 - single solution with analytical bounds
         ! 1 - single solution with PML bounds
         ! 2 - one auxiliary problem that gives boundary points with analytical bounds
         ! 3 - one auxiliary problem with PML bounds
         ! 4 - one aux problem, PML bounds and Effective currents
         ! 5 - lacunaes without Poisson
         ! 6 - lacunaes with Poisson
         ! 7 - pml without deviding
         ! 8 - sommerfeld abc
         ! 9 - higdon abc
         ! 10 - betz-mittra abc
         ! 11 - mur abc
       ! Scheme type
         integer :: main_schemetype;
         integer :: aux_schemetype;
         ! 0 - Yee 2-2 scheme
       ! Source type
         integer :: sourcetype;
         ! 0 - Dipole
         ! 1 - SDipole
       ! PML Parameters  
         integer :: PML_thickness;
         real*8  :: PML_param;
       ! Poisson parameters
         integer :: npx;
         integer :: npy;
         integer :: npz;
       ! Effective currents params
         integer :: Namu;
         integer :: Ndmu;
         real*8 :: auxmaxtime;
         real*8 :: sgm_s;
         integer :: auxgaptime;
       ! ABC Parameters
         real*8 :: hig1, hig2;
         integer :: dtmode, dxmode, dxauxmode, dtauxmode;
         integer :: gnNum;
         integer :: hwNum;
         integer :: hwEvaE;
         integer :: hwEvaP;
         integer :: hwEvaEa;
         integer :: hwEvaPa;
    end type tSolutionStarter

    type, public :: tFieldError
       real*8 :: val, absval, fval;
       integer :: x, y, z;
    end type tFieldError

end module commonvars






