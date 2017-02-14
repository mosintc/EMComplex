!-----------------------------------------------------------------------
!  Copyright 2017 Mikhail Osintcev
!  This file is part of the EMtool developed at NCSU
!-----------------------------------------------------------------------
! This module contains part of implementation of tSolution class

module solutionclass
! Source class that describes analytical solution and sources  
  use commonvars
  use problemclass
  use meshclass
  use bufferclass
  use timersclass
  implicit none;
  
  type, public :: tSolution
     integer :: sid;
     integer :: soltype;
     ! Currents type
     integer :: main_currentstype;
     integer :: main_pmlcurrentstype;
     integer :: aux_currentstype;
     integer :: aux_pmlcurrentstype;
     ! Effective Currents type
     integer :: aux_effcurrentstype;
     ! Scheme type
     integer :: main_schemetype;
     integer :: aux_schemetype;
     ! Source type
     integer :: main_sourcetype;
     integer :: aux_sourcetype;
     ! Main problem boundaries
     integer :: mainEboundarytype;
     integer :: mainHboundarytype;
         ! You can choose one of following ABC for E field and (-1, 0, 1, 2) for H field!
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
         ! 9 - givoli-neta abc Ey_x=0, Exact elsewhere
         ! 10 - givoli-neta abc Ey_x=0, Unsplit pml elsewhere
         ! 12 - givoli-neta abc Ey_x=0, big aux task elsewhere
         ! 13 - hagstrow-warburton abc Ey_x=0, Exact elsewhere
         ! 14 - hagstrow-warburton abc Ey_x=0, Unsplit PML elsewhere
         ! 15 - HWeva abc Ey_x=0, Exact elsewhere
         ! 16 - HWEva ABC Ey_x=Nx, Unpslit PML elsewhere     
     ! Mesh description
     integer :: Nx, Ny, Nz, Nt;      ! discretization values
     real*8  :: hhx, hhy, hhz, ht;   ! differentials
     ! Problems
     type(TProblem) :: mainproblem;
     ! Poisson parameters
     integer :: npx, npy, npz;     
     ! Aux Problems
       integer :: auxnum;
       class(TProblem), allocatable :: auxproblems(:);
       type(TTimers) :: timers;
       integer :: dNxaux, dNyaux, dNzaux;      ! discretization values
       integer :: auxEboundarytype;
       integer :: auxHboundarytype;
     ! Effective currents params
       real*8 :: auxmaxtime;
       real*8 :: sgm_s;
       integer :: Namu;
       integer :: Ndmu;
       integer :: auxlifetime;
       integer :: auxgaptime;
     ! Buffer
       class(TBuffer), pointer :: buffer;
     ! PMLs Parameters  
       integer :: main_PML_thickness;
       real*8  :: main_PML_param;
       integer :: aux_PML_thickness;
       real*8  :: aux_PML_param;
       ! ABC Parameters
       real*8 :: hig1, hig2;
       integer :: dtmode, dxmode, dxauxmode, dtauxmode;
       integer :: gnNum;
       integer :: hwNum;
       integer :: hwEvaE;
       integer :: hwEvaP;
       integer :: hwEvaEa;
       integer :: hwEvaPa;
    contains
      procedure solution_init;
      procedure solution_destroy;
      procedure solution_DestroyAuxProblems
      procedure solution_build;
      procedure solution_DoStep;
      ! Propagation procedures
      procedure propagateSingleProblem;
      procedure propagateOneAuxProblem;
      procedure PropagateOneAuxProblemBufCurrents;
      procedure PropagateQuasiLacunaes;
      procedure PropagateLacunaesPoisson;
      procedure PropagateLacunaesPoissonTiming;
      procedure PropagatePurePML;
      procedure PropagatePurePMLwithABC;
      procedure PropagateLacunasNoEffCurrentsNoPoisson;
      procedure PropagateLacunasNoEffCurrentsNoPoissonOneSource;
      ! Output procedures
      procedure solution_report
      procedure solution_printstep
      procedure solution_writeoutput
      procedure solution_writecheckinfo
      procedure solution_WriteTimes
      ! Auxiliary problems
      procedure addAuxProblem
      procedure dropAuxProblem
      procedure manageAuxProblems
  end type tSolution

contains

!Solution_Init------------------------------------------------------------------   
    subroutine solution_Init(this, starter, pbuffer)
    !Initialize the solution
      class(tSolution) :: this
      class(tSolutionStarter) :: starter
      type(tProblemStarter) :: pst
      class(TBuffer), pointer :: pbuffer;
      integer :: i;
      this%sid = starter%sid;
      this%soltype = starter%soltype;
      if ((starter%Nx>0).AND.(starter%Ny>0).AND.(starter%Nz>0)) then
         this%Nx = starter%Nx-1*conv;
         this%Ny = starter%Ny-1*conv;
         this%Nz = starter%Nz-1*conv;
         this%hhx = xsize/this%Nx;
         this%hhy = ysize/this%Ny;
         this%hhz = zsize/this%Nz;
         this%ht = this%hhx/3.0d0/cc ! 0.05d0/dfloat(conv);   !1/cc/sqrt(1/this%hhx**2+1/this%hhy**2+1/this%hhz**2);
         this%Nt = nint(ftime / this%ht);
      else
         fatalerr = 400;
         write(*,*) 'FATAL ERROR, Solution is initiated with zero points in one dimensions';
         write(*,*) 'Nx', starter%Nx, 'Ny', starter%Ny, 'Nz', starter%Nz;
      endif

      this%npx = starter%npx;
      this%npy = starter%npy;
      this%npz = starter%npz;

      if ((starter%dNxaux>=0).AND.(starter%dNyaux>=0).AND.(starter%dNzaux>=0)) then
         this%dNxaux = starter%dNxaux;
         this%dNyaux = starter%dNyaux;
         this%dNzaux = starter%dNzaux;
      else
         fatalerr = 402;
         write(*,*) 'FATAL ERROR, Auxiliary problems size correction should be > 0!';   
      endif

      this%main_schemetype = starter%main_schemetype;
      this%main_sourcetype = starter%sourcetype;
      this%aux_schemetype = starter%aux_schemetype;
      this%aux_sourcetype = starter%sourcetype;
      
      this%main_PML_thickness = starter%PML_thickness;
      this%main_PML_param = starter%PML_param;
      this%aux_PML_thickness = starter%PML_thickness;
      this%aux_PML_param = starter%PML_param;
      this%auxmaxtime = starter%auxmaxtime;
      this%sgm_s=starter%sgm_s;
      this%auxgaptime = starter%auxgaptime;     
      
      this%auxlifetime = nint((1.0*sqrt(xsize**2+ysize**2+zsize**2)/cc+2*this%auxmaxtime)/this%ht)+this%auxgaptime;

      this%hig1 = starter%hig1;
      this%hig2 = starter%hig2;
      this%gnNum = starter%gnNum;
      this%hwNum = starter%hwNum;
      this%hwEvaE = starter%hwEvaE;
      this%hwEvaP = starter%hwEvaP;
      this%hwEvaEa = starter%hwEvaEa;
      this%hwEvaPa = starter%hwEvaPa;

      this%dtmode = starter%dtmode;
      this%dxmode = starter%dxmode;
      this%dxauxmode = starter%dxauxmode;
      this%dtauxmode = starter%dtauxmode;

      allocate(cDNu(0:23));
         cDNu(0) = -51480.0d0*(q**15);
         cDNu(1) = 360360.0d0*(q**14);
         cDNu(2) = -1081080.0d0*(q**13);
         cDNu(3) = 1801800.0d0*(q**12);
         cDNu(4) = -1801800.0d0*(q**11);
         cDNu(5) = 1081080.0d0*(q**10);
         cDNu(6) = -360360.0d0*(q**9);
         cDNu(7) = 51480.0d0*(q**8);
         
         cDNu(8) = -720720.0d0*(q**15)
         cDNu(9) = 4684680.0d0*(q**14)
         cDNu(10) = -12972960.0d0*(q**13)
         cDNu(11) = 19819800.0d0*(q**12)
         cDNu(12) = -18018000.0d0*(q**11)
         cDNu(13) = 9729720.0d0*(q**10)
         cDNu(14) = -2882880.0d0*(q**9)
         cDNu(15) = 360360.0d0*(q**8)
         
         cDNu(16) = -9369360.0d0*(q**15)
         cDNu(17) = 56216160.0d0*(q**14)
         cDNu(18) = -142702560.0d0*(q**13)
         cDNu(19) = 198198000.0d0*(q**12)
         cDNu(20) = -162162000.0d0*(q**11)
         cDNu(21) = 77837760.0d0*(q**10)
         cDNu(22) = -20180160.0d0*(q**9)
         cDNu(23) = 2162160.0d0*(q**8)
         
         allocate(cNu(0:7));
         cNu(0) = -3432.0d0*(q**15);
         cNu(1) = 25740.0d0*(q**14);
         cNu(2) = -83160.0d0*(q**13);
         cNu(3) = 150150.0d0*(q**12);
         cNu(4) = -163800.0d0*(q**11);
         cNu(5) = 108108.0d0*(q**10);
         cNu(6) = -40040.0d0*(q**9);
         cNu(7) = 6435.0d0*(q**8);

         timestitles(0) = 'MainSource';
         timestitles(1) = 'Main E int';
         timestitles(2) = 'Main H int';
         timestitles(3) = 'Main E bnd';
         timestitles(4) = 'Main H bnd';
         timestitles(5) = 'Poisson';
         timestitles(6) = 'Manage Aux';
         timestitles(7) = 'Output';
         timestitles(8) = 'Analytical';
         timestitles(9) = 'GetError';
         timestitles(10) = 'AuxSources';
         timestitles(11) = 'Aux E int';
         timestitles(12) = 'Aux H int';
         timestitles(13) = 'Aux E bnd';
         timestitles(14) = 'Aux H bnd';
         timestitles(15) = 'PML E int';
         timestitles(16) = 'PML H int';
         timestitles(17) = 'BufClear';
         timestitles(18) = 'SaveFields';
         timestitles(19) = 'AuxSaveFie';
      ! Initialize problems
      select case (this%soltype)     
      case (0) ! Single solution - Analytical E boundaries
         this%main_currentstype = 0; ! Currents from Source
         pst%ptype = 0;
         pst%new_id = 0;
         pst%Nx = this%Nx;        
         pst%Ny = this%Ny;        
         pst%Nz = this%Nz;       
         pst%Nt = this%Nt;      
         pst%hhx = this%hhx;
         pst%hhy = this%hhy;
         pst%hhz = this%hhz;
         pst%ht = this%ht;
         pst%npx = this%npx;        
         pst%npy = this%npy;       
         pst%npz = this%npz;
         pst%schemetype = this%main_schemetype;
         pst%sourcetype = this%main_sourcetype;
         pst%currentstype = this%main_currentstype;     
         pst%effcurrentstype = 0;  ! Common effective currents
         pst%Eboundarytype = 1;    ! E - Analytical solution  
         pst%Hboundarytype = 2;    ! H - Scheme
         pst%Namu = starter%Namu;
         pst%Ndmu = starter%Ndmu;
         pst%auxmaxtime = this%auxmaxtime;
         pst%sgm_s = this%sgm_s;
         this%buffer => pbuffer;
         call this%buffer%buffer_init(this%Nx,this%Ny,this%Nz);
         call this%mainproblem%problem_init(pst, this%buffer)
      case (1) ! Single solution - PML E boundaries
         this%main_currentstype = 0; ! Currents from Source
         this%main_pmlcurrentstype = 0; ! Currents from Source
         pst%ptype = 0;
         pst%new_id = 0;
         pst%Nx = this%Nx;       
         pst%Ny = this%Ny;       
         pst%Nz = this%Nz;       
         pst%Nt = this%Nt;       
         pst%hhx = this%hhx;
         pst%hhy = this%hhy;
         pst%hhz = this%hhz;
         pst%ht = this%ht;
         pst%npx = this%npx;        
         pst%npy = this%npy;       
         pst%npz = this%npz;
         pst%schemetype = this%main_schemetype;
         pst%sourcetype = this%main_sourcetype;
         pst%currentstype = this%main_currentstype;     
         pst%effcurrentstype = 0;  ! Common effective currents
         pst%Eboundarytype = 3;    ! E - Unsplit PML
         pst%Hboundarytype = 2;    ! H - Scheme
         pst%PML_thickness = this%main_PML_thickness;
         pst%PML_param = this%main_PML_param;
         pst%Namu = starter%Namu;
         pst%Ndmu = starter%Ndmu;
         pst%auxmaxtime = this%auxmaxtime;
         pst%pmlcurrentstype = this%main_pmlcurrentstype;
         pst%sgm_s = this%sgm_s;
         this%buffer => pbuffer;
         call this%buffer%buffer_init(this%Nx,this%Ny,this%Nz);
         call this%mainproblem%problem_init(pst, this%buffer);
      case (100) ! Unsplit PML, Lacunaes, the main source is cut, NO Poisson
         !Buffer
         this%buffer => pbuffer;
         call this%buffer%buffer_init(this%Nx + 2*this%dNxaux,this%Ny + 2*this%dNyaux,this%Nz + 2*this%dNzaux);
         !Main problem
         this%main_currentstype = 0; ! Currents from Source
         pst%ptype = 0;
         pst%new_id = 0;
         pst%Nx = this%Nx;       
         pst%Ny = this%Ny;       
         pst%Nz = this%Nz;       
         pst%Nt = this%Nt;       
         pst%hhx = this%hhx;
         pst%hhy = this%hhy;
         pst%hhz = this%hhz;
         pst%npx = this%npx;        
         pst%npy = this%npy;       
         pst%npz = this%npz; 
         pst%ht = this%ht;
         pst%schemetype = this%main_schemetype;
         pst%sourcetype = this%main_sourcetype;
         pst%currentstype = this%main_currentstype;         
         pst%Eboundarytype = 4;    ! E - From aux problems
         pst%Hboundarytype = 2;    ! H - Scheme
         pst%PML_thickness = this%main_PML_thickness;
         pst%PML_param = this%main_PML_param;
         pst%Namu = starter%Namu;
         pst%Ndmu = starter%Ndmu;
         pst%auxmaxtime = this%auxmaxtime;
         pst%sgm_s = this%sgm_s;   
         pst%bufoffsetx=this%dNxaux;
         pst%bufoffsety=this%dNyaux;
         pst%bufoffsetz=this%dNzaux;
         pst%maindx=0;
         pst%maindy=0;
         pst%maindz=0;
         call this%mainproblem%problem_init(pst, this%buffer)
         call this%mainproblem%UploadMainMuToBuffer;            
         ! Aux problems
         this%auxEboundarytype = 3;    ! E - PML
         this%auxHboundarytype = 2;    ! H - Scheme
         this%aux_currentstype = 3;    ! Currents from Source cutted by Theta
         this%aux_pmlcurrentstype = 3; ! Currents from Source cutted by Theta
         this%auxnum = 0;
         allocate(this%auxproblems(0:auxmaxcount-1));
         call this%timers%timers_init(auxmaxcount);
         call this%AddAuxProblem;
      case (101) ! Pure Lacunas, the main source is cut, NO Poisson
          this%auxmaxtime = 0.5*0.268*sqrt(xsize**2+ysize**2+zsize**2)/cc;
         write(*,*) 'AUXMAXTIME', this%auxmaxtime;
         this%dNxaux = nint(1.268*this%Nx*sqrt(3.0d0)/2);        
         this%dNyaux = nint(1.268*this%Ny*sqrt(3.0d0)/2);       
         this%dNzaux = nint(1.268*this%Nz*sqrt(3.0d0)/2);
         this%auxlifetime = nint(1.268*sqrt(3.0d0)*this%Nx*this%hhx/this%ht);
         !Buffer
         this%buffer => pbuffer;
         call this%buffer%buffer_init(this%Nx + 2*this%dNxaux,this%Ny + 2*this%dNyaux,this%Nz + 2*this%dNzaux);
         !Main problem
         this%main_currentstype = 0; ! Currents from Source
         pst%ptype = 0;
         pst%new_id = 0;
         pst%Nx = this%Nx;       
         pst%Ny = this%Ny;       
         pst%Nz = this%Nz;       
         pst%Nt = this%Nt;       
         pst%hhx = this%hhx;
         pst%hhy = this%hhy;
         pst%hhz = this%hhz;
         pst%npx = this%npx;        
         pst%npy = this%npy;       
         pst%npz = this%npz;
         pst%ht = this%ht;
         pst%schemetype = this%main_schemetype;
         pst%sourcetype = this%main_sourcetype;
         pst%currentstype = this%main_currentstype;         
         pst%Eboundarytype = 4;    ! E - From aux problems
         pst%Hboundarytype = 2;    ! H - Scheme
         pst%PML_thickness = this%main_PML_thickness;
         pst%PML_param = this%main_PML_param;
         pst%Namu = starter%Namu;
         pst%Ndmu = starter%Ndmu;
         pst%auxmaxtime = this%auxmaxtime;
         pst%sgm_s = this%sgm_s;   
         pst%bufoffsetx=this%dNxaux;
         pst%bufoffsety=this%dNyaux;
         pst%bufoffsetz=this%dNzaux;
         pst%maindx=0;
         pst%maindy=0;
         pst%maindz=0;
         call this%mainproblem%problem_init(pst, this%buffer)
         call this%mainproblem%UploadMainMuToBuffer;            
         ! Aux problems
         this%auxEboundarytype = 0;    ! E - No Boundary
         this%auxHboundarytype = 2;    ! H - Scheme
         this%aux_currentstype = 1;    ! Currents from Source cutted by Theta
         this%aux_pmlcurrentstype = 1; ! Currents from Source cutted by Theta
         this%auxnum = 0;
         allocate(this%auxproblems(0:auxmaxcount-1));
         call this%timers%timers_init(auxmaxcount);
         call this%AddAuxProblem;
      case (102) ! Pure Lacunas, the main source is cut, NO Poisson, MAIN SOURCE FROM BUFFER FOR ALL AUX PROBLEMS
         this%auxmaxtime = 0.5*0.268*sqrt(xsize**2+ysize**2+zsize**2)/cc;
         write(*,*) 'AUXMAXTIME', this%auxmaxtime;
         this%dNxaux = nint(1.268*this%Nx*sqrt(3.0d0)/2);        
         this%dNyaux = nint(1.268*this%Ny*sqrt(3.0d0)/2);       
         this%dNzaux = nint(1.268*this%Nz*sqrt(3.0d0)/2);
         this%auxlifetime = nint(1.268*sqrt(3.0d0)*this%Nx*this%hhx/this%ht);
         !Buffer
         this%buffer => pbuffer;
         call this%buffer%buffer_init(this%Nx + 2*this%dNxaux,this%Ny + 2*this%dNyaux,this%Nz + 2*this%dNzaux);
         !Main problem
         this%main_currentstype = 0; ! Currents from Source
         pst%ptype = 0;
         pst%new_id = 0;
         pst%Nx = this%Nx;       
         pst%Ny = this%Ny;       
         pst%Nz = this%Nz;       
         pst%Nt = this%Nt;       
         pst%hhx = this%hhx;
         pst%hhy = this%hhy;
         pst%hhz = this%hhz;
         pst%npx = this%npx;        
         pst%npy = this%npy;       
         pst%npz = this%npz;
         pst%ht = this%ht;
         pst%schemetype = this%main_schemetype;
         pst%sourcetype = this%main_sourcetype;
         pst%currentstype = this%main_currentstype;         
         pst%Eboundarytype = 4;    ! E - From aux problems
         pst%Hboundarytype = 2;    ! H - Scheme
         pst%PML_thickness = this%main_PML_thickness;
         pst%PML_param = this%main_PML_param;
         pst%Namu = starter%Namu;
         pst%Ndmu = starter%Ndmu;
         pst%auxmaxtime = this%auxmaxtime;
         pst%sgm_s = this%sgm_s;   
         pst%bufoffsetx=this%dNxaux;
         pst%bufoffsety=this%dNyaux;
         pst%bufoffsetz=this%dNzaux;
         pst%maindx=0;
         pst%maindy=0;
         pst%maindz=0;
         call this%mainproblem%problem_init(pst, this%buffer)
         call this%mainproblem%UploadMainMuToBuffer;            
         ! Aux problems
         this%auxEboundarytype = -1;    ! E - No Boundary
         this%auxHboundarytype = -1;    ! H - No Boundary
         this%aux_currentstype = 3;    ! Currents from buffer cutted by Theta
         this%aux_pmlcurrentstype = 3; ! Currents from buffer cutted by Theta
         this%auxnum = 0;
         allocate(this%auxproblems(0:auxmaxcount-1));
         call this%timers%timers_init(auxmaxcount);
         call this%AddAuxProblem;     
      case (2) ! Main problem and one auxiliary problem - Analytical E boundaries
         !Buffer
         this%buffer => pbuffer;
         call this%buffer%buffer_init(this%Nx + 2*this%dNxaux,this%Ny + 2*this%dNyaux,this%Nz + 2*this%dNzaux);
         !Main problem
         this%main_currentstype = 0; ! Currents from Source
         this%main_pmlcurrentstype = 0; ! Currents from Source
         pst%ptype = 0;
         pst%new_id = 0;
         pst%Nx = this%Nx;        
         pst%Ny = this%Ny;        
         pst%Nz = this%Nz;      
         pst%Nt = this%Nt;        
         pst%hhx = this%hhx;
         pst%hhy = this%hhy;
         pst%hhz = this%hhz;
         pst%ht = this%ht;
         pst%npx = this%npx;        
         pst%npy = this%npy;       
         pst%npz = this%npz;
         pst%schemetype = this%main_schemetype;
         pst%sourcetype = this%main_sourcetype;
         pst%currentstype = this%main_currentstype;
         pst%effcurrentstype = 0;  ! Common effective currents
         pst%Eboundarytype = 4;    ! E - From Buffer
         pst%Hboundarytype = 2;    ! H - Scheme
         pst%PML_thickness = this%main_PML_thickness;
         pst%PML_param = this%main_PML_param;
         pst%Namu = starter%Namu;
         pst%Ndmu = starter%Ndmu;
         pst%auxmaxtime = this%auxmaxtime;
         pst%pmlcurrentstype = this%main_pmlcurrentstype;
         pst%sgm_s = this%sgm_s;
         pst%bufoffsetx=this%dNxaux;
         pst%bufoffsety=this%dNyaux;
         pst%bufoffsetz=this%dNzaux;      
         call this%mainproblem%problem_init(pst, this%buffer)
         !Aux problem
         this%auxEboundarytype = 1;   ! E - Analytical values
         this%auxHboundarytype = 2;   ! H - Scheme
         this%aux_currentstype = 0;   ! Currents from Source
         this%aux_effcurrentstype = 0 ! Common effective currents
         allocate(this%auxproblems(0:0));
         call this%timers%timers_init(auxmaxcount);
         this%auxnum = 0;
         this%auxlifetime = this%Nt+1;
         call this%AddAuxProblem;
      case (3) ! Main problem and one auxiliary problem - PML E boundaries
         !Buffer
         this%buffer => pbuffer;
         call this%buffer%buffer_init(this%Nx + 2*this%dNxaux,this%Ny + 2*this%dNyaux,this%Nz + 2*this%dNzaux);
         !Main problem
         this%main_currentstype = 0; ! Currents from Source
         pst%ptype = 0;
         pst%new_id = 0;
         pst%Nx = this%Nx;       
         pst%Ny = this%Ny;        
         pst%Nz = this%Nz;       
         pst%Nt = this%Nt;       
         pst%hhx = this%hhx;
         pst%hhy = this%hhy;
         pst%hhz = this%hhz;
         pst%ht = this%ht;
         pst%npx = this%npx;        
         pst%npy = this%npy;       
         pst%npz = this%npz;
         pst%schemetype = this%main_schemetype;
         pst%sourcetype = this%main_sourcetype;
         pst%currentstype = this%main_currentstype;
         pst%effcurrentstype = 0;  ! Common effective currents
         pst%Eboundarytype = 4;    ! E - From aux problems
         pst%Hboundarytype = 2;    ! H - Scheme
         pst%PML_thickness = this%main_PML_thickness;
         pst%PML_param = this%main_PML_param;
         pst%Namu = starter%Namu;
         pst%Ndmu = starter%Ndmu;
         pst%auxmaxtime = this%auxmaxtime;
         pst%sgm_s = this%sgm_s;      
         pst%bufoffsetx=this%dNxaux;
         pst%bufoffsety=this%dNyaux;
         pst%bufoffsetz=this%dNzaux;      
         call this%mainproblem%problem_init(pst, this%buffer)
         !Aux problem
         this%auxEboundarytype = 3;   ! E - Unsplit PML
         this%auxHboundarytype = 2;   ! H - Scheme
         this%aux_currentstype = 0;   ! Currents from Source
         this%aux_effcurrentstype = 0 ! Common effective currents
         allocate(this%auxproblems(0:0));
         call this%timers%timers_init(auxmaxcount);
         this%auxnum = 0;
         this%auxlifetime = this%Nt+1;
         call this%AddAuxProblem;
      case (4) ! Main problem and one auxiliary problem (Effective currents, PML E boundaries)
         !Buffer
         this%buffer => pbuffer;
         call this%buffer%buffer_init(this%Nx + 2*this%dNxaux,this%Ny + 2*this%dNyaux,this%Nz + 2*this%dNzaux);
         !Main problem
         this%main_currentstype = 0; ! Currents from Source
         pst%ptype = 0;
         pst%new_id = 0;
         pst%Nx = this%Nx;
         pst%Ny = this%Ny;
         pst%Nz = this%Nz;
         pst%Nt = this%Nt;
         pst%hhx = this%hhx;
         pst%hhy = this%hhy;
         pst%hhz = this%hhz;
         pst%ht = this%ht;
         pst%npx = this%npx;        
         pst%npy = this%npy;       
         pst%npz = this%npz;
         pst%schemetype = this%main_schemetype;
         pst%sourcetype = this%main_sourcetype;
         pst%currentstype = this%main_currentstype;
         pst%effcurrentstype = 0;  ! Common effective currents
         pst%Eboundarytype = 4;    ! E - From aux problems
         pst%Hboundarytype = 2;    ! H - Scheme
         pst%PML_thickness = this%main_PML_thickness;
         pst%PML_param = this%main_PML_param;
         pst%Namu = starter%Namu;
         pst%Ndmu = starter%Ndmu;
         pst%auxmaxtime = this%auxmaxtime;
         pst%sgm_s = this%sgm_s;   
         pst%bufoffsetx=this%dNxaux;
         pst%bufoffsety=this%dNyaux;
         pst%bufoffsetz=this%dNzaux;      
         call this%mainproblem%problem_init(pst, this%buffer)
         !Aux problems
         this%auxEboundarytype = 3;   ! E - Unsplit PML
         this%auxHboundarytype = 2;   ! H - Scheme
         this%aux_currentstype = 2;   ! Effective currents
         this%aux_effcurrentstype = 0 ! Common effective currents
         allocate(this%auxproblems(0:0));
         call this%timers%timers_init(auxmaxcount);
         this%auxnum = 0;
         this%auxlifetime = this%Nt+1;
         call this%AddAuxProblem;
      case (5) ! Main problem and one auxiliary problems cutted by Theta function (Effective currents, PML E boundaries)   
         !Buffer
         this%buffer => pbuffer;
         call this%buffer%buffer_init(this%Nx + 2*this%dNxaux,this%Ny + 2*this%dNyaux,this%Nz + 2*this%dNzaux);
         !Main problem
         this%main_currentstype = 0; ! Currents from Source
         pst%ptype = 0;
         pst%new_id = 0;
         pst%Nx = this%Nx;        
         pst%Ny = this%Ny;       
         pst%Nz = this%Nz;       
         pst%Nt = this%Nt;      
         pst%hhx = this%hhx;
         pst%hhy = this%hhy;
         pst%hhz = this%hhz;
         pst%ht = this%ht;
         pst%npx = this%npx;        
         pst%npy = this%npy;       
         pst%npz = this%npz;      
         pst%schemetype = this%main_schemetype;
         pst%sourcetype = this%main_sourcetype;
         pst%currentstype = this%main_currentstype;
         pst%effcurrentstype = 0;  ! Common effective currents
         pst%Eboundarytype = 4;    ! E - From aux problems and static field
         pst%Hboundarytype = 2;    ! H - Scheme
         pst%PML_thickness = this%main_PML_thickness;
         pst%PML_param = this%main_PML_param;
         pst%Namu = starter%Namu;
         pst%Ndmu = starter%Ndmu;
         pst%auxmaxtime = this%auxmaxtime;
         pst%sgm_s = this%sgm_s;   
         pst%bufoffsetx=this%dNxaux;
         pst%bufoffsety=this%dNyaux;
         pst%bufoffsetz=this%dNzaux;      
         call this%mainproblem%problem_init(pst, this%buffer)
         !Aux problems
         this%auxEboundarytype = 3;   ! E - Buffer cutted with Theta
         this%auxHboundarytype = 2;   ! H - Scheme
         this%aux_currentstype = 3;   ! Effective currents cutted by Theta
         this%aux_effcurrentstype = 0 ! Common effective currents
         this%auxnum = 0;
         allocate(this%auxproblems(0:auxmaxcount-1));
         call this%timers%timers_init(auxmaxcount);
         call this%AddAuxProblem;
      case (6) ! With Poisson Correction! Main problem and one auxiliary problems cutted by Theta function (Effective currents, PML E boundaries)    
         !Buffer
         this%buffer => pbuffer;
         call this%buffer%buffer_init(this%Nx + 2*this%dNxaux,this%Ny + 2*this%dNyaux,this%Nz + 2*this%dNzaux);
         !Main problem
         this%main_currentstype = 0; ! Currents from Source
         pst%ptype = 0;
         pst%new_id = 0;
         pst%Nx = this%Nx;       
         pst%Ny = this%Ny;       
         pst%Nz = this%Nz;       
         pst%Nt = this%Nt;       
         pst%hhx = this%hhx;
         pst%hhy = this%hhy;
         pst%hhz = this%hhz;
         pst%npx = this%npx;        
         pst%npy = this%npy;       
         pst%npz = this%npz; 
         pst%ht = this%ht;
         pst%schemetype = this%main_schemetype;
         pst%sourcetype = this%main_sourcetype;
         pst%currentstype = this%main_currentstype;
         pst%Eboundarytype = 4;    ! E - From aux problems
         pst%effcurrentstype = 1;  ! Effective currents with Poisson
         pst%Hboundarytype = 2;    ! H - Scheme
         pst%PML_thickness = this%main_PML_thickness;
         pst%PML_param = this%main_PML_param;
         pst%Namu = starter%Namu;
         pst%Ndmu = starter%Ndmu;
         pst%auxmaxtime = this%auxmaxtime;
         pst%sgm_s = this%sgm_s;   
         pst%bufoffsetx=this%dNxaux;
         pst%bufoffsety=this%dNyaux;
         pst%bufoffsetz=this%dNzaux;
         pst%maindx=0;
         pst%maindy=0;
         pst%maindz=0;
         call this%mainproblem%problem_init(pst, this%buffer)
         call this%mainproblem%UploadMainMuToBuffer;            
         ! Aux problems
         this%auxEboundarytype = 3;   ! E - PML
         this%auxHboundarytype = 2;   ! H - Scheme
         this%aux_currentstype = 4;   ! Effective currents cutted by Theta with Poisson correction
         this%aux_effcurrentstype = 0 ! Common effective currents
         this%auxnum = 0;
         allocate(this%auxproblems(0:auxmaxcount-1));
         call this%timers%timers_init(auxmaxcount);
         call this%AddAuxProblem;
      case (7) ! Pure PML (non devided PML)
         this%main_currentstype = 0; ! Currents from Source
         pst%ptype = 0;
         pst%new_id = 0;
         pst%Nx = this%Nx;       
         pst%Ny = this%Ny;       
         pst%Nz = this%Nz;       
         pst%Nt = this%Nt;       
         pst%hhx = this%hhx;
         pst%hhy = this%hhy;
         pst%hhz = this%hhz;
         pst%ht = this%ht;
         pst%npx = this%npx;        
         pst%npy = this%npy;       
         pst%npz = this%npz;
         pst%schemetype = this%main_schemetype;
         pst%sourcetype = this%main_sourcetype;
         pst%currentstype = this%main_currentstype;     
         pst%effcurrentstype = 0;  ! Common effective currents
         pst%Eboundarytype = 3;    ! E - Unsplit PML
         pst%Hboundarytype = 2;    ! H - Scheme
         pst%PML_thickness = this%main_PML_thickness;
         pst%PML_param = this%main_PML_param;
         pst%Namu = starter%Namu;
         pst%Ndmu = starter%Ndmu;
         pst%auxmaxtime = this%auxmaxtime;
         pst%sgm_s = this%sgm_s;
         this%buffer => pbuffer;
         call this%buffer%buffer_init(this%Nx,this%Ny,this%Nz);
         call this%mainproblem%problem_init(pst, this%buffer)
      case (8) ! Single solution - Sommerfeld ABC 
         this%main_currentstype = 0; ! Currents from Source
         pst%ptype = 0;
         pst%new_id = 0;
         pst%Nx = this%Nx;        
         pst%Ny = this%Ny;        
         pst%Nz = this%Nz;       
         pst%Nt = this%Nt;      
         pst%hhx = this%hhx;
         pst%hhy = this%hhy;
         pst%hhz = this%hhz;
         pst%ht = this%ht;
         pst%npx = this%npx;        
         pst%npy = this%npy;       
         pst%npz = this%npz;
         pst%schemetype = this%main_schemetype;
         pst%sourcetype = this%main_sourcetype;
         pst%currentstype = this%main_currentstype;     
         pst%effcurrentstype = 0;  ! Common effective currents
         pst%Eboundarytype = 5;    ! E - Analytical solution  
         pst%Hboundarytype = 2;    ! H - Scheme
         pst%Namu = starter%Namu;
         pst%Ndmu = starter%Ndmu;
         pst%auxmaxtime = this%auxmaxtime;
         pst%sgm_s = this%sgm_s;
         this%buffer => pbuffer;
         call this%buffer%buffer_init(this%Nx,this%Ny,this%Nz);
         call this%mainproblem%problem_init(pst, this%buffer)
      case (9) ! Single solution - Higdon ABC 
         this%main_currentstype = 0; ! Currents from Source
         pst%ptype = 0;
         pst%new_id = 0;
         pst%Nx = this%Nx;        
         pst%Ny = this%Ny;        
         pst%Nz = this%Nz;       
         pst%Nt = this%Nt;      
         pst%hhx = this%hhx;
         pst%hhy = this%hhy;
         pst%hhz = this%hhz;
         pst%ht = this%ht;
         pst%npx = this%npx;        
         pst%npy = this%npy;       
         pst%npz = this%npz;
         pst%schemetype = this%main_schemetype;
         pst%sourcetype = this%main_sourcetype;
         pst%currentstype = this%main_currentstype;     
         pst%effcurrentstype = 0;  ! Common effective currents
         pst%Eboundarytype = 6;    ! E - Analytical solution  
         pst%Hboundarytype = 2;    ! H - Scheme
         pst%Namu = starter%Namu;
         pst%Ndmu = starter%Ndmu;
         pst%auxmaxtime = this%auxmaxtime;
         pst%sgm_s = this%sgm_s;
         pst%hig1 = this%hig1;
         pst%hig2 = this%hig2;
         pst%gnNum = this%gnNum;
         this%buffer => pbuffer;
         call this%buffer%buffer_init(this%Nx,this%Ny,this%Nz);
         call this%mainproblem%problem_init(pst, this%buffer)
      case (10) ! Single solution - Betz-Mittra ABC 
         this%main_currentstype = 0; ! Currents from Source
         pst%ptype = 0;
         pst%new_id = 0;
         pst%Nx = this%Nx;        
         pst%Ny = this%Ny;        
         pst%Nz = this%Nz;       
         pst%Nt = this%Nt;      
         pst%hhx = this%hhx;
         pst%hhy = this%hhy;
         pst%hhz = this%hhz;
         pst%ht = this%ht;
         pst%npx = this%npx;        
         pst%npy = this%npy;       
         pst%npz = this%npz;
         pst%schemetype = this%main_schemetype;
         pst%sourcetype = this%main_sourcetype;
         pst%currentstype = this%main_currentstype;     
         pst%effcurrentstype = 0;  ! Common effective currents
         pst%Eboundarytype = 7;    ! E - Analytical solution  
         pst%Hboundarytype = 2;    ! H - Scheme
         pst%Namu = starter%Namu;
         pst%Ndmu = starter%Ndmu;
         pst%auxmaxtime = this%auxmaxtime;
         pst%sgm_s = this%sgm_s;
         pst%hig1 = this%hig1;
         pst%hig2 = this%hig2;
         pst%gnNum = this%gnNum;
         this%buffer => pbuffer;
         call this%buffer%buffer_init(this%Nx,this%Ny,this%Nz);
         call this%mainproblem%problem_init(pst, this%buffer)
      case (11) ! Single solution - Mur ABC 
         this%main_currentstype = 0; ! Currents from Source
         pst%ptype = 0;
         pst%new_id = 0;
         pst%Nx = this%Nx;        
         pst%Ny = this%Ny;        
         pst%Nz = this%Nz;       
         pst%Nt = this%Nt;      
         pst%hhx = this%hhx;
         pst%hhy = this%hhy;
         pst%hhz = this%hhz;
         pst%ht = this%ht;
         pst%npx = this%npx;        
         pst%npy = this%npy;       
         pst%npz = this%npz;
         pst%schemetype = this%main_schemetype;
         pst%sourcetype = this%main_sourcetype;
         pst%currentstype = this%main_currentstype;     
         pst%effcurrentstype = 0;  ! Common effective currents
         pst%Eboundarytype = 8;    ! E - Analytical solution  
         pst%Hboundarytype = 2;    ! H - Scheme
         pst%Namu = starter%Namu;
         pst%Ndmu = starter%Ndmu;
         pst%auxmaxtime = this%auxmaxtime;
         pst%sgm_s = this%sgm_s;
         this%buffer => pbuffer;
         call this%buffer%buffer_init(this%Nx,this%Ny,this%Nz);
         call this%mainproblem%problem_init(pst, this%buffer)
      case (12) ! Single solution - Givoli-Neta High-order ABC on X=0 and exact solution elsewhere
         this%main_currentstype = 0; ! Currents from Source
         pst%ptype = 0;
         pst%new_id = 0;
         pst%Nx = this%Nx;        
         pst%Ny = this%Ny;        
         pst%Nz = this%Nz;       
         pst%Nt = this%Nt;      
         pst%hhx = this%hhx;
         pst%hhy = this%hhy;
         pst%hhz = this%hhz;
         pst%ht = this%ht;
         pst%npx = this%npx;        
         pst%npy = this%npy;       
         pst%npz = this%npz;
         pst%schemetype = this%main_schemetype;
         pst%sourcetype = this%main_sourcetype;
         pst%currentstype = this%main_currentstype;     
         pst%effcurrentstype = 0;  ! Common effective currents
         pst%Eboundarytype = 9;    ! E - Givoli-Neta with exact elsewhere  
         pst%Hboundarytype = 2;    ! H - Scheme
         pst%Namu = starter%Namu;
         pst%Ndmu = starter%Ndmu;
         pst%auxmaxtime = this%auxmaxtime;
         pst%sgm_s = this%sgm_s;
         pst%gnNum = this%gnNum;
         pst%dtmode = this%dtmode;
         pst%dxmode = this%dxmode;
         pst%dxauxmode = this%dxauxmode;
         pst%dtauxmode = this%dtauxmode;
         this%buffer => pbuffer;
         call this%buffer%buffer_init(this%Nx,this%Ny,this%Nz);
         call this%mainproblem%problem_init(pst, this%buffer);
      case (13) ! Single solution - Givoli-Neta High-order ABC on X=0 and Unsplit PML solution elsewhere  
         this%main_currentstype = 0; ! Currents from Source
         this%main_pmlcurrentstype = 0; ! Currents from Source
         pst%ptype = 0;
         pst%new_id = 0;
         pst%Nx = this%Nx;       
         pst%Ny = this%Ny;       
         pst%Nz = this%Nz;       
         pst%Nt = this%Nt;       
         pst%hhx = this%hhx;
         pst%hhy = this%hhy;
         pst%hhz = this%hhz;
         pst%ht = this%ht;
         pst%npx = this%npx;        
         pst%npy = this%npy;       
         pst%npz = this%npz;
         pst%schemetype = this%main_schemetype;
         pst%sourcetype = this%main_sourcetype;
         pst%currentstype = this%main_currentstype;     
         pst%effcurrentstype = 0;  ! Common effective currents
         pst%Eboundarytype = 10;    ! E - Givoli-Neta with Unsplit PML elsewhere  
         pst%Hboundarytype = 2;    ! H - Scheme
         pst%PML_thickness = this%main_PML_thickness;
         pst%PML_param = this%main_PML_param;
         pst%Namu = starter%Namu;
         pst%Ndmu = starter%Ndmu;
         pst%auxmaxtime = this%auxmaxtime;
         pst%pmlcurrentstype = this%main_pmlcurrentstype;
         pst%sgm_s = this%sgm_s;
         pst%gnNum = this%gnNum;
         pst%dtmode = this%dtmode;
         pst%dxmode = this%dxmode;
         pst%dxauxmode = this%dxauxmode;
         pst%dtauxmode = this%dtauxmode;
         this%buffer => pbuffer;
         call this%buffer%buffer_init(this%Nx,this%Ny,this%Nz);
         call this%mainproblem%problem_init(pst, this%buffer)
      case (14) ! Main problem with Givoli-Neta at EY_x=0 and one big auxiliary problem without boundaries
         !Buffer
         this%buffer => pbuffer;
         call this%buffer%buffer_init(this%Nx + 2*this%dNxaux,this%Ny + 2*this%dNyaux,this%Nz + 2*this%dNzaux);
         !Main problem
         this%main_currentstype = 0; ! Currents from Source
         this%main_pmlcurrentstype = 0; ! Currents from Source
         pst%ptype = 0;
         pst%new_id = 0;
         pst%Nx = this%Nx;       
         pst%Ny = this%Ny;        
         pst%Nz = this%Nz;       
         pst%Nt = this%Nt;       
         pst%hhx = this%hhx;
         pst%hhy = this%hhy;
         pst%hhz = this%hhz;
         pst%ht = this%ht;
         pst%npx = this%npx;        
         pst%npy = this%npy;       
         pst%npz = this%npz;
         pst%schemetype = this%main_schemetype;
         pst%sourcetype = this%main_sourcetype;
         pst%currentstype = this%main_currentstype;
         pst%effcurrentstype = 0;  ! Common effective currents
         pst%Eboundarytype = 11;   ! E - From aux problems and Givoli-Neta
         pst%Hboundarytype = 2;    ! H - Scheme
         pst%PML_thickness = this%main_PML_thickness;
         pst%PML_param = this%main_PML_param;
         pst%Namu = starter%Namu;
         pst%Ndmu = starter%Ndmu;
         pst%auxmaxtime = this%auxmaxtime;
         pst%sgm_s = this%sgm_s;
         pst%gnNum = this%gnNum;
         pst%dtmode = this%dtmode;
         pst%dxmode = this%dxmode;
         pst%dxauxmode = this%dxauxmode;
         pst%dtauxmode = this%dtauxmode;
         pst%bufoffsetx=this%dNxaux;
         pst%bufoffsety=this%dNyaux;
         pst%bufoffsetz=this%dNzaux;
         call this%mainproblem%problem_init(pst, this%buffer)
         !Aux problem
         this%auxEboundarytype = -1;   ! E - no boundaries at all
         this%auxHboundarytype = 2;    ! H - Scheme
         this%aux_currentstype = 0;    ! Currents from Source
         this%aux_effcurrentstype = 0  ! Common effective currents
         allocate(this%auxproblems(0:0));
         call this%timers%timers_init(auxmaxcount);
         this%auxnum = 0;
         this%auxlifetime = this%Nt+1;
         call this%AddAuxProblem;
      case (15) ! Lacunaes with Givoli-Neta X=0 + PML elsewhere, the main source is cut by Theta, NO Poisson
         !Buffer
         this%buffer => pbuffer;
         call this%buffer%buffer_init(this%Nx + 2*this%dNxaux,this%Ny + 2*this%dNyaux,this%Nz + 2*this%dNzaux);
         !Main problem
         this%main_currentstype = 0; ! Currents from Source
         pst%ptype = 0;
         pst%new_id = 0;
         pst%Nx = this%Nx;       
         pst%Ny = this%Ny;       
         pst%Nz = this%Nz;       
         pst%Nt = this%Nt;       
         pst%hhx = this%hhx;
         pst%hhy = this%hhy;
         pst%hhz = this%hhz;
         pst%npx = this%npx;        
         pst%npy = this%npy;       
         pst%npz = this%npz; 
         pst%ht = this%ht;
         pst%schemetype = this%main_schemetype;
         pst%sourcetype = this%main_sourcetype;
         pst%currentstype = this%main_currentstype;         
         pst%Eboundarytype = 4;    ! E - From aux problems
         pst%Hboundarytype = 2;    ! H - Scheme
         pst%PML_thickness = this%main_PML_thickness;
         pst%PML_param = this%main_PML_param;
         pst%Namu = starter%Namu;
         pst%Ndmu = starter%Ndmu;
         pst%auxmaxtime = this%auxmaxtime;
         pst%sgm_s = this%sgm_s;   
         pst%bufoffsetx=this%dNxaux;
         pst%bufoffsety=this%dNyaux;
         pst%bufoffsetz=this%dNzaux;
         pst%maindx=0;
         pst%maindy=0;
         pst%maindz=0;
         call this%mainproblem%problem_init(pst, this%buffer)
         call this%mainproblem%UploadMainMuToBuffer;            
         ! Aux problems
         this%auxEboundarytype = 10;   ! E - Givoli-Neta + PML
         this%auxHboundarytype = 2;    ! H - Scheme
         this%aux_currentstype = 1;    ! Currents from Source cutted by Theta
         this%aux_pmlcurrentstype = 1; ! Currents from Source cutted by Theta
         this%auxnum = 0;
         allocate(this%auxproblems(0:auxmaxcount-1));
         call this%timers%timers_init(auxmaxcount);
         call this%AddAuxProblem;   
      case (16) ! Main problem and one auxiliary problem without boundaries
         !Buffer
         this%buffer => pbuffer;
         call this%buffer%buffer_init(this%Nx + 2*this%dNxaux,this%Ny + 2*this%dNyaux,this%Nz + 2*this%dNzaux);
         !Main problem
         this%main_currentstype = 0; ! Currents from Source
         pst%ptype = 0;
         pst%new_id = 0;
         pst%Nx = this%Nx;       
         pst%Ny = this%Ny;       
         pst%Nz = this%Nz;       
         pst%Nt = this%Nt;       
         pst%hhx = this%hhx;
         pst%hhy = this%hhy;
         pst%hhz = this%hhz;
         pst%npx = this%npx;        
         pst%npy = this%npy;       
         pst%npz = this%npz; 
         pst%ht = this%ht;
         pst%schemetype = this%main_schemetype;
         pst%sourcetype = this%main_sourcetype;
         pst%currentstype = this%main_currentstype;
         pst%Eboundarytype = 4;    ! E - From aux problems
         pst%effcurrentstype = 0;  ! Effective currents no Poisson
         pst%Hboundarytype = 2;    ! H - Scheme
         pst%PML_thickness = this%main_PML_thickness;
         pst%PML_param = this%main_PML_param;
         pst%Namu = starter%Namu;
         pst%Ndmu = starter%Ndmu;
         pst%auxmaxtime = this%auxmaxtime;
         pst%sgm_s = this%sgm_s;   
         pst%bufoffsetx=this%dNxaux;
         pst%bufoffsety=this%dNyaux;
         pst%bufoffsetz=this%dNzaux;
         pst%maindx=0;
         pst%maindy=0;
         pst%maindz=0;
         call this%mainproblem%problem_init(pst, this%buffer)
         call this%mainproblem%UploadMainMuToBuffer;            
         ! Aux problems
         this%auxEboundarytype = 0;   ! E - No Boundary
         this%auxHboundarytype = 2;   ! H - Scheme
         this%aux_currentstype = 1;   ! Original source cutted by Theta
         this%aux_effcurrentstype = 0 ! Common effective currents
         this%auxnum = 0;
         allocate(this%auxproblems(0:auxmaxcount-1));
         call this%timers%timers_init(auxmaxcount);
         call this%AddAuxProblem;            
      case (17) ! Single solution - Hagstrom-Warburton ABC on X=0 and exact solution elsewhere
         this%main_currentstype = 0; ! Currents from Source
         pst%ptype = 0;
         pst%new_id = 0;
         pst%Nx = this%Nx;        
         pst%Ny = this%Ny;        
         pst%Nz = this%Nz;       
         pst%Nt = this%Nt;      
         pst%hhx = this%hhx;
         pst%hhy = this%hhy;
         pst%hhz = this%hhz;
         pst%ht = this%ht;
         pst%npx = this%npx;        
         pst%npy = this%npy;       
         pst%npz = this%npz;
         pst%schemetype = this%main_schemetype;
         pst%sourcetype = this%main_sourcetype;
         pst%currentstype = this%main_currentstype;     
         pst%effcurrentstype = 0;  ! Common effective currents
         pst%Eboundarytype = 13;   ! E - Hagstrom-Warburton with exact elsewhere  
         pst%Hboundarytype = 2;    ! H - Scheme
         pst%Namu = starter%Namu;
         pst%Ndmu = starter%Ndmu;
         pst%auxmaxtime = this%auxmaxtime;
         pst%sgm_s = this%sgm_s;
         pst%hwNum = this%hwNum;
         pst%dtmode = this%dtmode;
         pst%dxmode = this%dxmode;
         pst%dxauxmode = this%dxauxmode;
         pst%dtauxmode = this%dtauxmode;
         this%buffer => pbuffer;
         call this%buffer%buffer_init(this%Nx,this%Ny,this%Nz);
         call this%mainproblem%problem_init(pst, this%buffer);
      case (18) ! Single solution - Hagstrom-Warburton High-order ABC on X=0 and Unsplit PML solution elsewhere  
         this%main_currentstype = 0; ! Currents from Source
         this%main_pmlcurrentstype = 0; ! Currents from Source
         pst%ptype = 0;
         pst%new_id = 0;
         pst%Nx = this%Nx;       
         pst%Ny = this%Ny;       
         pst%Nz = this%Nz;       
         pst%Nt = this%Nt;       
         pst%hhx = this%hhx;
         pst%hhy = this%hhy;
         pst%hhz = this%hhz;
         pst%ht = this%ht;
         pst%npx = this%npx;        
         pst%npy = this%npy;       
         pst%npz = this%npz;
         pst%schemetype = this%main_schemetype;
         pst%sourcetype = this%main_sourcetype;
         pst%currentstype = this%main_currentstype;     
         pst%effcurrentstype = 0;  ! Common effective currents
         pst%Eboundarytype = 14;   ! E - Hagstrom-Warburton with Unsplit PML elsewhere  
         pst%Hboundarytype = 2;    ! H - Scheme
         pst%PML_thickness = this%main_PML_thickness;
         pst%PML_param = this%main_PML_param;
         pst%Namu = starter%Namu;
         pst%Ndmu = starter%Ndmu;
         pst%auxmaxtime = this%auxmaxtime;
         pst%pmlcurrentstype = this%main_pmlcurrentstype;
         pst%sgm_s = this%sgm_s;
         pst%hwNum = this%hwNum;
         pst%dtmode = this%dtmode;
         pst%dxmode = this%dxmode;
         pst%dxauxmode = this%dxauxmode;
         pst%dtauxmode = this%dtauxmode;
         this%buffer => pbuffer;
         call this%buffer%buffer_init(this%Nx,this%Ny,this%Nz);
         call this%mainproblem%problem_init(pst, this%buffer)
      case (19) ! Lacunaes with Hagstrom-Warburton X=0 + PML elsewhere, the main source is cut by Theta, NO Poisson
         !Buffer
         this%buffer => pbuffer;
         call this%buffer%buffer_init(this%Nx + 2*this%dNxaux,this%Ny + 2*this%dNyaux,this%Nz + 2*this%dNzaux);
         !Main problem
         this%main_currentstype = 0; ! Currents from Source
         pst%ptype = 0;
         pst%new_id = 0;
         pst%Nx = this%Nx;       
         pst%Ny = this%Ny;       
         pst%Nz = this%Nz;       
         pst%Nt = this%Nt;       
         pst%hhx = this%hhx;
         pst%hhy = this%hhy;
         pst%hhz = this%hhz;
         pst%npx = this%npx;        
         pst%npy = this%npy;       
         pst%npz = this%npz; 
         pst%ht = this%ht;
         pst%schemetype = this%main_schemetype;
         pst%sourcetype = this%main_sourcetype;
         pst%currentstype = this%main_currentstype;         
         pst%Eboundarytype = 4;    ! E - From aux problems
         pst%Hboundarytype = 2;    ! H - Scheme
         pst%PML_thickness = this%main_PML_thickness;
         pst%PML_param = this%main_PML_param;
         pst%Namu = starter%Namu;
         pst%Ndmu = starter%Ndmu;
         pst%auxmaxtime = this%auxmaxtime;
         pst%sgm_s = this%sgm_s;   
         pst%bufoffsetx=this%dNxaux;
         pst%bufoffsety=this%dNyaux;
         pst%bufoffsetz=this%dNzaux;
         pst%maindx=0;
         pst%maindy=0;
         pst%maindz=0;
         call this%mainproblem%problem_init(pst, this%buffer)
         call this%mainproblem%UploadMainMuToBuffer;            
         ! Aux problems
         this%auxEboundarytype = 14;   ! E - Hagstrom-Warburton + PML
         this%auxHboundarytype = 2;    ! H - Scheme
         this%aux_currentstype = 1;    ! Currents from Source cutted by Theta
         this%aux_pmlcurrentstype = 1; ! Currents from Source cutted by Theta
         this%auxnum = 0;
         allocate(this%auxproblems(0:auxmaxcount-1));
         call this%timers%timers_init(auxmaxcount);
         call this%AddAuxProblem;
      case (20) ! Single solution - HWeva ABC on Ey_X=Nx, Exact solution elsewhere
         this%main_currentstype = 0; ! Currents from Source
         pst%ptype = 0;
         pst%new_id = 0;
         pst%Nx = this%Nx;        
         pst%Ny = this%Ny;        
         pst%Nz = this%Nz;       
         pst%Nt = this%Nt;      
         pst%hhx = this%hhx;
         pst%hhy = this%hhy;
         pst%hhz = this%hhz;
         pst%ht = this%ht;
         pst%npx = this%npx;        
         pst%npy = this%npy;       
         pst%npz = this%npz;
         pst%schemetype = this%main_schemetype;
         pst%sourcetype = this%main_sourcetype;
         pst%currentstype = this%main_currentstype;     
         pst%effcurrentstype = 0;  ! Common effective currents
         pst%Eboundarytype = 15;   ! E - Hagstrom-Warburton with exact elsewhere  
         pst%Hboundarytype = 2;    ! H - Scheme
         pst%Namu = starter%Namu;
         pst%Ndmu = starter%Ndmu;
         pst%auxmaxtime = this%auxmaxtime;
         pst%hig1 = this%hig1;
         pst%hig2 = this%hig2;
         pst%sgm_s = this%sgm_s;
         pst%hwNum = this%hwNum;
         pst%hwEvaE = this%hwEvaE;
         pst%hwEvaP = this%hwEvaP;
         pst%hwEvaEa = this%hwEvaEa;
         pst%hwEvaPa = this%hwEvaPa;
         pst%dtmode = this%dtmode;
         pst%dxmode = this%dxmode;
         pst%dxauxmode = this%dxauxmode;
         pst%dtauxmode = this%dtauxmode;
         this%buffer => pbuffer;
         call this%buffer%buffer_init(this%Nx,this%Ny,this%Nz);
         call this%mainproblem%problem_init(pst, this%buffer);
      case (21) ! Single solution - HWeva ABC on Ey_X=Nx, Unsplit PML elsewhere  
         this%main_currentstype = 0; ! Currents from Source
         this%main_pmlcurrentstype = 0; ! Currents from Source
         pst%ptype = 0;
         pst%new_id = 0;
         pst%Nx = this%Nx;       
         pst%Ny = this%Ny;       
         pst%Nz = this%Nz;       
         pst%Nt = this%Nt;       
         pst%hhx = this%hhx;
         pst%hhy = this%hhy;
         pst%hhz = this%hhz;
         pst%ht = this%ht;
         pst%npx = this%npx;        
         pst%npy = this%npy;       
         pst%npz = this%npz;
         pst%schemetype = this%main_schemetype;
         pst%sourcetype = this%main_sourcetype;
         pst%currentstype = this%main_currentstype;     
         pst%effcurrentstype = 0;  ! Common effective currents
         pst%Eboundarytype = 16;   ! E - HWEva with Unsplit PML elsewhere  
         pst%Hboundarytype = 2;    ! H - Scheme
         pst%PML_thickness = this%main_PML_thickness;
         pst%PML_param = this%main_PML_param;
         pst%Namu = starter%Namu;
         pst%Ndmu = starter%Ndmu;
         pst%auxmaxtime = this%auxmaxtime;
         pst%pmlcurrentstype = this%main_pmlcurrentstype;
         pst%sgm_s = this%sgm_s;
         pst%hig1 = this%hig1;
         pst%hig2 = this%hig2;
         pst%hwNum = this%hwNum;
         pst%hwEvaE = this%hwEvaE;
         pst%hwEvaP = this%hwEvaP;
         pst%hwEvaEa = this%hwEvaEa;
         pst%hwEvaPa = this%hwEvaPa;
         pst%dtmode = this%dtmode;
         pst%dxmode = this%dxmode;
         pst%dxauxmode = this%dxauxmode;
         pst%dtauxmode = this%dtauxmode;
         this%buffer => pbuffer;
         call this%buffer%buffer_init(this%Nx,this%Ny,this%Nz);
         call this%mainproblem%problem_init(pst, this%buffer)
      case (22) ! HWeva ABC on Ey_X=Nx, Unsplit PML elsewhere, lacunaes, no poisson, no effective currents, Source 0 only
         !Buffer
         this%buffer => pbuffer;
         call this%buffer%buffer_init(this%Nx + 2*this%dNxaux,this%Ny + 2*this%dNyaux,this%Nz + 2*this%dNzaux);
         !Main problem
         this%main_currentstype = 0; ! Currents from Source
         pst%ptype = 0;
         pst%new_id = 0;
         pst%Nx = this%Nx;       
         pst%Ny = this%Ny;       
         pst%Nz = this%Nz;       
         pst%Nt = this%Nt;       
         pst%hhx = this%hhx;
         pst%hhy = this%hhy;
         pst%hhz = this%hhz;
         pst%npx = this%npx;        
         pst%npy = this%npy;       
         pst%npz = this%npz; 
         pst%ht = this%ht;
         pst%schemetype = this%main_schemetype;
         pst%sourcetype = this%main_sourcetype;
         pst%currentstype = this%main_currentstype;         
         pst%Eboundarytype = 4;    ! E - From aux problems
         pst%Hboundarytype = 2;    ! H - Scheme
         pst%PML_thickness = this%main_PML_thickness;
         pst%PML_param = this%main_PML_param;
         pst%Namu = starter%Namu;
         pst%Ndmu = starter%Ndmu;
         pst%auxmaxtime = this%auxmaxtime;
         pst%sgm_s = this%sgm_s;   
         pst%bufoffsetx=this%dNxaux;
         pst%bufoffsety=this%dNyaux;
         pst%bufoffsetz=this%dNzaux;
         pst%maindx=0;
         pst%maindy=0;
         pst%maindz=0;
         call this%mainproblem%problem_init(pst, this%buffer)
         call this%mainproblem%UploadMainMuToBuffer;            
         ! Aux problems
         this%auxEboundarytype = 16;   ! E - HWEva + PML
         this%auxHboundarytype = 2;    ! H - Scheme
         this%aux_currentstype = 1;    ! Currents from Source cutted by Theta
         this%aux_pmlcurrentstype = 1; ! Currents from Source cutted by Theta
         this%auxnum = 0;
         allocate(this%auxproblems(0:auxmaxcount-1));
         call this%timers%timers_init(auxmaxcount);
         call this%AddAuxProblem;
      case (23) ! HWeva ABC on Ey_X=Nx, Unsplit PML elsewhere, lacunaes, Poisson
         !Buffer
         this%buffer => pbuffer;
         call this%buffer%buffer_init(this%Nx + 2*this%dNxaux,this%Ny + 2*this%dNyaux,this%Nz + 2*this%dNzaux);
         !Main problem
         this%main_currentstype = 0; ! Currents from Source
         pst%ptype = 0;
         pst%new_id = 0;
         pst%Nx = this%Nx;       
         pst%Ny = this%Ny;       
         pst%Nz = this%Nz;       
         pst%Nt = this%Nt;       
         pst%hhx = this%hhx;
         pst%hhy = this%hhy;
         pst%hhz = this%hhz;
         pst%npx = this%npx;        
         pst%npy = this%npy;       
         pst%npz = this%npz; 
         pst%ht = this%ht;
         pst%schemetype = this%main_schemetype;
         pst%sourcetype = this%main_sourcetype;
         pst%currentstype = this%main_currentstype;         
         pst%Eboundarytype = 4;    ! E - From aux problems
         pst%effcurrentstype = 1;  ! Effective currents with Poisson
         pst%Hboundarytype = 2;    ! H - Scheme
         pst%PML_thickness = this%main_PML_thickness;
         pst%PML_param = this%main_PML_param;
         pst%Namu = starter%Namu;
         pst%Ndmu = starter%Ndmu;
         pst%auxmaxtime = this%auxmaxtime;
         pst%sgm_s = this%sgm_s;   
         pst%bufoffsetx=this%dNxaux;
         pst%bufoffsety=this%dNyaux;
         pst%bufoffsetz=this%dNzaux;
         pst%maindx=0;
         pst%maindy=0;
         pst%maindz=0;
         call this%mainproblem%problem_init(pst, this%buffer)
         call this%mainproblem%UploadMainMuToBuffer;            
         ! Aux problems
         this%auxEboundarytype = 16;   ! E - HWEva + PML
         this%auxHboundarytype = 2;    ! H - Scheme
         this%aux_currentstype = 4;    ! Effective currents cutted by Theta with Poisson correction
         this%aux_effcurrentstype = 0 ! Common effective currents
         this%aux_pmlcurrentstype = 4; ! Effective currents cutted by Theta with Poisson correction
         this%auxnum = 0;
         allocate(this%auxproblems(0:auxmaxcount-1));
         call this%timers%timers_init(auxmaxcount);
         call this%AddAuxProblem;
      case (25) ! Sommerfeld ABCe, lacunaes, Poisson    
         !Buffer
         this%buffer => pbuffer;
         call this%buffer%buffer_init(this%Nx + 2*this%dNxaux,this%Ny + 2*this%dNyaux,this%Nz + 2*this%dNzaux);
         !Main problem
         this%main_currentstype = 0; ! Currents from Source
         pst%ptype = 0;
         pst%new_id = 0;
         pst%Nx = this%Nx;       
         pst%Ny = this%Ny;       
         pst%Nz = this%Nz;       
         pst%Nt = this%Nt;       
         pst%hhx = this%hhx;
         pst%hhy = this%hhy;
         pst%hhz = this%hhz;
         pst%npx = this%npx;        
         pst%npy = this%npy;       
         pst%npz = this%npz; 
         pst%ht = this%ht;
         pst%schemetype = this%main_schemetype;
         pst%sourcetype = this%main_sourcetype;
         pst%currentstype = this%main_currentstype;
         pst%Eboundarytype = 4;    ! E - From aux problems
         pst%effcurrentstype = 1;  ! Effective currents with Poisson
         pst%Hboundarytype = 2;    ! H - Scheme
         pst%PML_thickness = this%main_PML_thickness;
         pst%PML_param = this%main_PML_param;
         pst%Namu = starter%Namu;
         pst%Ndmu = starter%Ndmu;
         pst%auxmaxtime = this%auxmaxtime;
         pst%sgm_s = this%sgm_s;   
         pst%bufoffsetx=this%dNxaux;
         pst%bufoffsety=this%dNyaux;
         pst%bufoffsetz=this%dNzaux;
         pst%maindx=0;
         pst%maindy=0;
         pst%maindz=0;
         call this%mainproblem%problem_init(pst, this%buffer)
         call this%mainproblem%UploadMainMuToBuffer;            
         ! Aux problems
         this%auxEboundarytype = 0;   ! E - Low order ABC
         this%auxHboundarytype = 2;   ! H - Scheme
         this%aux_currentstype = 4;   ! Effective currents cutted by Theta with Poisson correction
         this%aux_effcurrentstype = 0 ! Common effective currents
         this%auxnum = 0;
         allocate(this%auxproblems(0:auxmaxcount-1));
         call this%timers%timers_init(auxmaxcount);
         call this%AddAuxProblem;   
      case default
         fatalerr = 401;
         write(*,*) 'Incorrect solution type during the solution initialization', this%soltype;
      end select
      call this%solution_report;
    end subroutine solution_Init;

    
!solution_report----------------------------------------------------------    
    subroutine solution_report(this)
    ! Print on the screen details of the solution
      class(tSolution) :: this;
      integer :: i;
      write(*,*) '------SOLUTION REPORT-----'
      write(*,*) 'Solution ID:', this%sid;
      write(*,*) 'Type:', this%soltype;
      select case (this%soltype)
      case(0)
         write(*,*) 'Single solution with analytical bounds';
      case(1)
         write(*,*) 'Single solution with unsplit PML bounds';
      case(100)
         write(*,*) 'Unsplit PML, Lacunaes, the main source is cut, NO Poisson';
      case(101)
         write(*,*) 'Pure Lacunas, the main source is cut, NO Poisson';
      case(102)
         write(*,*) 'Pure Lacunas, the main source is cut, NO Poisson, One Source for Aux Problems';   
      case(2)
         write(*,*) 'Main problem and one aux problem with analytical bounds';
      case(3)
         write(*,*) 'Main problem and one aux problem with unsplit PML bounds';
      case(4)
         write(*,*) 'Main problem, one aux problem with unsplit PML bounds and Effective currents';
      case(5)
         write(*,*) 'Quasi-lacunaes, no Poisson, aux problems with unsplit PML';
      case(6)
         write(*,*) 'Lacunaes with Poisson correction, aux problems with unsplit PML';
      case(7)
         write(*,*) 'Pure PML';
      case(8)
         write(*,*) 'Single solution with Sommerfeld ABC';
      case(9)
         write(*,*) 'Single solution with Higdon ABC';
      case(10)
         write(*,*) 'Single solution with Betz-Mittra ABC';
      case(11)
         write(*,*) 'Single solution with Mur ABC';
      case(12)
         write(*,*) 'Givoli-Neta ABC Single solution and exact elsewhere';
      case(13)
         write(*,*) 'Givoli-Neta ABC Single solution and Unsplit PML elsewhere';   
      case(14)
         write(*,*) 'Givoli-Neta in main problem and one big aux problem without boundaries';
      case(15)
         write(*,*) 'Givoli-Neta X=0 and PML elsewhere, Lacunaes, the main source is cut, NO Poisson';   
      case(16)
         write(*,*) 'Main problem and big aux problems without boundaries';
      case(17)
         write(*,*) 'Hagstrom-Warburton ABC Single solution and exact solution elsewhere';
      case(18)
         write(*,*) 'Hagstrom-Warburton ABC Single solution and Unsplit PML elsewhere';
      case(19)
         write(*,*) 'Hagstrom-Warburton X=0 and PML elsewhere, Lacunaes, the main source is cut, NO Poisson';
      case(20)
         write(*,*) 'HWEva Ey_X=0 ABC Single solution and exact solution elsewhere';
      case(21)
         write(*,*) 'HWEva Ey_X=0 ABC Single solution and Unsplit PML elsewhere';
      case(22)
         write(*,*) 'HWEva Ey_X=0 ABC and PML elsewhere, Lacunaes, the main source is cut, NO Poisson';
      case(23)
         write(*,*) 'HWeva ABC on Ey_X=Nx, Unsplit PML elsewhere, lacunaes, Poisson';
      case(25)
         write(*,*) 'Sommerfeld ABCe, lacunaes, Poisson';    
      case default
         write(*,*) 'UNKNOWN type of solution'
      end select
      write(*,*) 'Nx:', this%Nx;
      write(*,*) 'Ny:', this%Ny;
      write(*,*) 'Nz:', this%Nz;
      write(*,*) 'Nt:', this%Nt;
      write(*,*) 'hhx:', this%hhx;
      write(*,*) 'hhy:', this%hhy;
      write(*,*) 'hhz:', this%hhz;
      write(*,*) 'ht:', this%ht;
      write(*,*) 'NXaux:', this%dNXaux;
      write(*,*) 'NYaux:', this%dNYaux;
      write(*,*) 'NZaux:', this%dNZaux;
      select case (this%main_schemetype)
      case(0)
         write(*,*) 'Scheme: Classic Yee 2-2 scheme';
      case default
         write(*,*) 'Scheme: UNKNOWN scheme type';
      end select
      select case (this%main_sourcetype)
      case(0)
         write(*,*) 'Source: Dipole';
      case(1)
         write(*,*) 'Source: SDipole';
      case default
         write(*,*) 'Source: UNKNOWN source type';   
      end select
      if (this%soltype==5.OR.this%soltype==6.OR.this%soltype==14.OR.this%soltype==15.OR.this%soltype==16.OR.this%soltype==19.OR.this%soltype==22.OR.this%soltype==23.OR.this%soltype==101.OR.this%soltype==100.OR.this%soltype==25) then
         write(*,*) 'Domain diameter:', sqrt(xsize**2+ysize**2+zsize**2)/cc
         write(*,*) 'Aux problems standrad life time:', this%auxlifetime;
         write(*,*) 'Each aux will last more for (gap):', this%auxgaptime;
         write(*,*) 'Start and drop times:'
         write(*,*)  1, 'Aux task time:', 0,':',this%auxlifetime-nint(this%auxmaxtime/this%ht);
         do i=1,3
            write(*,*) i+1, 'Aux task time:', nint((this%sgm_s*this%auxmaxtime+(1+this%sgm_s)*this%auxmaxtime*(i-1))/this%ht),':',nint((this%sgm_s*this%auxmaxtime+(1+this%sgm_s)*this%auxmaxtime*(i-1))/this%ht) + this%auxlifetime;
         enddo
      endif
    end subroutine solution_report

    
!Solution_Destroy----------------------------------------------------------    
    subroutine solution_Destroy(this)
    ! Destroy solution object
      class(tSolution) :: this;
      select case (this%soltype)
      case (0) ! Analytical E boundaries
         call this%mainproblem%problem_destroy; 
      case (1) ! PML E boundaries
         call this%mainproblem%problem_destroy;
      case (100) ! Unsplit PML, Lacunaes, the main source is cut, NO Poisson
         call this%mainproblem%problem_destroy;
         call this%solution_DestroyAuxProblems;
         call this%timers%timers_destroy;
      case (101) ! Pure Lacunas, the main source is cut, NO Poisson
         call this%mainproblem%problem_destroy;
         call this%solution_DestroyAuxProblems;
         call this%timers%timers_destroy;
      case (102) ! Pure Lacunas, the main source is cut, NO Poisson, One Source for all aux problems
         call this%mainproblem%problem_destroy;
         call this%solution_DestroyAuxProblems;
         call this%timers%timers_destroy;         
      case (2) ! Auxiliary problem and analytical E bounds
         call this%mainproblem%problem_destroy;
         call this%solution_DestroyAuxProblems;
         call this%timers%timers_destroy;
      case (3) ! Auxiliary problem and PML E bounds
         call this%mainproblem%problem_destroy;
         call this%solution_DestroyAuxProblems;
         call this%timers%timers_destroy; 
      case (4) ! Auxiliary problem, Effective currents and PML E bounds
         call this%mainproblem%problem_destroy;
         call this%solution_DestroyAuxProblems;
         call this%timers%timers_destroy;
      case (5) ! Auxiliary problems, Effective currents and PML E bounds
         call this%mainproblem%problem_destroy;
         call this%solution_DestroyAuxProblems;
         call this%timers%timers_destroy;
      case (6) ! Auxiliary problems, Effective currents and PML E bounds
         call this%mainproblem%problem_destroy;
         call this%solution_DestroyAuxProblems;
         call this%timers%timers_destroy;
      case (7) ! Pure PML
         call this%mainproblem%problem_destroy;
      case (8) ! Sommerfeld ABC
         call this%mainproblem%problem_destroy;
      case (9) ! Higdon ABC
         call this%mainproblem%problem_destroy;
      case (10) ! Betz-Mittra ABC
         call this%mainproblem%problem_destroy;
      case (11) ! Mur ABC
         call this%mainproblem%problem_destroy;
      case (12) ! Givoli-Neta ABC X=0 with exact elsewhere
         call this%mainproblem%problem_destroy;
      case (13) ! Givoli-Neta ABC X=0 with UnsplitPML elsewhere 
         call this%mainproblem%problem_destroy;   
      case (14) ! Givoli-Neta in main problem and one big aux problem without boundaries
         call this%mainproblem%problem_destroy;
         call this%solution_DestroyAuxProblems;
         call this%timers%timers_destroy;
      case (15) ! Givoli-Neta X=0 + PML elsewhere with Lacunaes, the main source is cut by Theta, NO Poisson
         call this%mainproblem%problem_destroy;
         call this%solution_DestroyAuxProblems;
         call this%timers%timers_destroy;   
      case (16) ! Main problem and one big aux problem without boundaries
         call this%mainproblem%problem_destroy;
         call this%solution_DestroyAuxProblems;
         call this%timers%timers_destroy;
      case (17) ! Single solution - Hagstrom-Warburton ABC on X=0 and exact solution elsewhere
         call this%mainproblem%problem_destroy;
      case (18) ! Hagstrom-Warburton ABC X=0 with UnsplitPML elsewhere 
         call this%mainproblem%problem_destroy;
      case (19) ! Hagstrom-Warburton X=0 + PML elsewhere with Lacunaes, the main source is cut by Theta, NO Poisson
         call this%mainproblem%problem_destroy;
         call this%solution_DestroyAuxProblems;
         call this%timers%timers_destroy;
      case (20) ! HWEva ABC X=0 with Exact solution elsewhere 
         call this%mainproblem%problem_destroy;
      case (21) ! HWEva ABC X=0 with UnsplitPML elsewhere 
         call this%mainproblem%problem_destroy;
      case (22) ! HWEva Ey_X=0 ABC and PML elsewhere, Lacunaes, the main source is cut, NO Poisson
         call this%mainproblem%problem_destroy;
         call this%solution_DestroyAuxProblems;
         call this%timers%timers_destroy;
      case (23) ! HWeva ABC on Ey_X=Nx, Unsplit PML elsewhere, lacunaes, Poisson
         call this%mainproblem%problem_destroy;
         call this%solution_DestroyAuxProblems;
         call this%timers%timers_destroy;
      case (25) ! Sommerfeld ABCe, lacunaes, Poisson
         call this%mainproblem%problem_destroy;
         call this%solution_DestroyAuxProblems;
         call this%timers%timers_destroy;   
      end select 
    end subroutine solution_Destroy;

    
!Solution_DestroyAuxProblems----------------------------------------------------------
    subroutine solution_DestroyAuxProblems(this)
    ! Destroy solution object
      class(tSolution) :: this;
      integer :: i;
      do i=0,auxmaxcount-1
         if (this%timers%isused(i)==1) then
            call this%auxproblems(i)%problem_destroy;
         endif    
      enddo
      deallocate(this%auxproblems);
    end subroutine solution_DestroyAuxProblems;
    
    
!Solution_Build----------------------------------------------------------    
    subroutine solution_Build(this)
    ! Build the solution
      class(tSolution) :: this;
      integer :: t, i;
      real*8 :: ttime, tstart, tend;
      ttime = 0.0d0;
      
      do t=1,this%Nt
          call cpu_time(tstart)
          !$ tstart = OMP_get_wtime();
          
          call cpu_time(tstarts(17))
          !$ tstarts(17) = OMP_get_wtime(); 
               call this%buffer%clearFields;
          call cpu_time(tends(17))
          !$ tends(17) = OMP_get_wtime(); 

          if ((mod(t,file_erroroutputsteps)==0)) then
              call cpu_time(tstarts(8))
              !$ tstarts(8) = OMP_get_wtime();   
                 call this%mainproblem%getAnalyticSolution(t)
              call cpu_time(tends(8))
              !$ tends(8) = OMP_get_wtime(); 
          endif
          
          call this%solution_DoStep(t);
          
          if ((mod(t,screen_erroroutputsteps)==0).OR.(mod(t,file_erroroutputsteps)==0)) then
              call cpu_time(tstarts(9))
              !$ tstarts(9) = OMP_get_wtime();   
                  call this%mainproblem%getCurrentError(t);
              call cpu_time(tends(9))
              !$ tends(9) = OMP_get_wtime(); 
           endif
           
          call cpu_time(tstarts(7))
          !$ tstarts(7) = OMP_get_wtime();    
             call this%solution_WriteOutput(t)
          call cpu_time(tends(7))
          !$ tends(7) = OMP_get_wtime();
          
          call cpu_time(tend)
          !$ tend = OMP_get_wtime();  
          ttime = ttime + (tend-tstart);

          call cpu_time(tstarts(7))
          !$ tstarts(7) = OMP_get_wtime();   
          call this%solution_printstep(t, ttime);
          call cpu_time(tends(7))
          !$ tends(7) = OMP_get_wtime();
      enddo
    end subroutine solution_Build;


!Solution_DoStep----------------------------------------------------------    
    subroutine solution_DoStep(this,t)
    ! Perform one step of the solution
      class(tSolution) :: this;
      integer, intent(in) :: t;
      select case (this%soltype)
      case (0) ! Analytical E boundaries
         call this%propagateSingleProblem(t); 
      case (1) ! PML E boundaries
         call this%propagateSingleProblem(t);
      case (100) ! Unsplit PML, Lacunaes, the main source is cut, NO Poisson
         call this%PropagateLacunasNoEffCurrentsNoPoissonOneSource(t);
      case (101) ! Pure Lacunas, the main source is cut, NO Poisson
         call this%PropagateLacunasNoEffCurrentsNoPoisson(t);
      case (102) ! Pure Lacunas, the main source is cut, NO Poisson, One source for all aux problems
         call this%PropagateLacunasNoEffCurrentsNoPoissonOneSource(t);
      case (2) ! Auxiliary problem and analytical E bounds
         call this%propagateOneAuxProblem(t);
      case (3) ! Auxiliary problem and PML E bounds
         call this%PropagateOneAuxProblem(t);
      case (4) ! Auxiliary problem, Effective currents and PML E bounds
         call this%PropagateOneAuxProblemBufCurrents(t);
      case (5) ! Auxiliary problems, Effective currents and PML E bounds
         call this%PropagateQuasiLacunaes(t);
      case (6) ! Auxiliary problems, Effective currents and PML E bounds, Poisson correction
         call this%PropagateLacunaesPoisson(t);
      case (7) ! Pure PML
         call this%propagatePurePML(t);
      case (8) ! Sommerfeld ABC
         call this%propagateSingleProblem(t);
      case (9) ! Higdon ABC
         call this%propagateSingleProblem(t);
      case (10) ! Betz-Mittra ABC
         call this%propagateSingleProblem(t);
      case (11) ! Mur ABC
         call this%propagateSingleProblem(t);
      case (12) ! Givoli-Neta ABC X=0 with exact elsewhere
         call this%propagateSingleProblem(t);
      case (13) ! Givoli-Neta ABC X=0 with UnsplitPML elsewhere
         call this%propagatePurePMLwithABC(t);  
      case (14) ! Givoli-Neta in main problem and one big aux problem without boundaries
         call this%PropagateOneAuxProblem(t);
      case (15) ! Givoli-Neta X=0 + PML elsewhere with Lacunaes, the main source is cut by Theta, NO Poisson
         call this%PropagateLacunasNoEffCurrentsNoPoisson(t);   
      case (16) ! Main problem and one big aux problem without boundaries
         call this%PropagateLacunasNoEffCurrentsNoPoisson(t);
      case (17) ! Hagstrom-Warburton X=0 with exact elsewhere
         call this%propagateSingleProblem(t);
      case (18) ! Hagstrom-Warburton ABC X=0 with UnsplitPML elsewhere
         call this%propagatePurePMLwithABC(t);
      case (19) ! Hagstrom-Warburton ABC X=0 + PML elsewhere with Lacunaes, the main source is cut by Theta, NO Poisson
         call this%PropagateLacunasNoEffCurrentsNoPoisson(t);
      case (20) ! HWEva ABC Ey_x=Nx with exact elsewhere
         call this%propagateSingleProblem(t);
      case (21) ! HWEva ABC Ey_x=Nx with UnsplitPML elsewhere
         call this%propagatePurePMLwithABC(t);
      case (22) ! HWEva ABC Ey_x=Nx + PML elsewhere with Lacunaes, the main source is cut by Theta, NO Poisson
         call this%PropagateLacunasNoEffCurrentsNoPoisson(t);
      case (23) ! HWeva ABC on Ey_X=Nx, Unsplit PML elsewhere, lacunaes, Poisson
         call this%PropagateLacunaesPoisson(t);
      case (25) ! Sommerfeld ABCe, lacunaes, Poisson
         call this%PropagateLacunaesPoisson(t);   
      end select
    end subroutine solution_DoStep;

!PropagateSingleProblem----------------------------------------------------------    
    subroutine propagateSingleProblem(this,t)
    ! Propagate single main problem
      class(tSolution) :: this;
      integer, intent(in) :: t;
      integer :: i, j, k;
      call this%mainproblem%problem_DoIndependentStep(t)
    end subroutine PropagateSingleProblem;   

    
!PropagateWithAuxProblems----------------------------------------------------------    
    subroutine propagateOneAuxProblem(this,t)
    ! Propagate problem with auxiliary problems
      class(tSolution) :: this;
      integer :: i;
      integer, intent(in) :: t;
      call this%buffer%clearFields;    
      call this%auxproblems(0)%problem_DoIndependentStep(t)
      call this%auxproblems(0)%addEtobuffer;       
      call this%mainproblem%problem_DoIndependentStep(t)
    end subroutine PropagateOneAuxProblem;


!PropagateOneAuxProblemBufCurrents----------------------------------------------------------    
    subroutine PropagateOneAuxProblemBufCurrents(this,t)
    ! Propagate problem with effective currents for auxiliary problems
      class(tSolution) :: this;
      integer :: i;
      integer, intent(in) :: t;
      call this%buffer%clearFields;
      call this%mainproblem%saveFields(0);
      ! Main problem part
      call this%mainproblem%getSourceCurrents(t);          ! Fill J arrays with actual values of currents       
      call this%mainproblem%updateEInterior;               ! Calculate E field everywhere except outer boundary
      call this%mainproblem%JEeffBuild(t);
      call this%mainproblem%updateHInterior;               ! Calculate H field everywhere except outer boundary
      call this%mainproblem%JHeffBuild(t);
      call this%mainproblem%saveEffCurrentsToBuffer;
      ! Aux problem part
      call this%auxproblems(0)%problem_DoIndependentStep(t);
      call this%auxproblems(0)%addEtobuffer;
      ! Main problem finish     
      call this%mainproblem%fillEBoundary(t);              ! Fill E boundary with current Eboundary condition
      call this%mainproblem%fillHBoundary(t);              ! Fill H boundary with current Hboundary condition
      call this%solution_WriteCheckInfo(t);
    end subroutine PropagateOneAuxProblemBufCurrents;


!PropagateQuasiLacunaes----------------------------------------------------------    
    subroutine PropagateQuasiLacunaes(this,t)
    ! Propagate problem with effective currents for auxiliary problems
      class(tSolution) :: this;
      integer :: i, k;
      integer, intent(in) :: t;
      if (checkupdated==1) then
         call this%mainproblem%Ef%mesh_clearChA;           ! Clear checking arrays
         call this%mainproblem%Hf%mesh_clearChA;           ! Clear checking arrays
      endif
      call this%buffer%clearFields;       
      ! Main problem part
      call this%mainproblem%saveFields(0);
      call this%mainproblem%getSourceCurrents(t);          ! Fill J arrays with actual values of currents       
      call this%mainproblem%updateEInterior;               ! Calculate E field everywhere except outer boundary
      call this%mainproblem%JEeffBuild(t);
      call this%mainproblem%updateHInterior;               ! Calculate H field everywhere except outer boundary
      call this%mainproblem%JHeffBuild(t);
      call this%mainproblem%saveEffCurrentsToBuffer;
      ! Aux problem part
      do k=0,auxmaxcount-1
         if (this%timers%isused(k)==1) then
            call this%auxproblems(k)%problem_DoIndependentStep(t);
            call this%auxproblems(k)%addEtobuffer;
         endif
         !  The following steps allow to accumulate E static field, but it doesn't work
         
         !   if (this%timers%timer(k)>0.AND.k<this%auxnum) then
         !      write(*,*) 'Problem doing step:', k;
         !      call this%auxproblems(k)%problem_DoIndependentStep(t);
         !       write(*,*) 'Problem adding to buffer:', k;
         !      call this%auxproblems(k)%addEtobuffer;
         !   endif
         !   if (k<this%auxnum) then
         !      write(*,*) 'Problem doing step:', k;
         !      call this%auxproblems(k)%problem_DoIndependentStep(t);
         !   endif   
      enddo;   
      ! Main problem finish     
      call this%mainproblem%fillEBoundary(t);              ! Fill E boundary with current Eboundary condition
      call this%mainproblem%fillHBoundary(t);              ! Fill H boundary with current Hboundary condition
      if (checkupdated==1) then
         call this%mainproblem%Ef%mesh_checkarrays;        ! Check all the points are updated once
         call this%mainproblem%Hf%mesh_checkarrays;        ! Check all the points are updated once
      endif
      call this%solution_WriteCheckInfo(t);
      call this%ManageAuxProblems(t);
    end subroutine PropagateQuasiLacunaes;      


!PropagateLacunaesPoisson----------------------------------------------------------    
    subroutine PropagateLacunaesPoisson(this,t)
    ! Propagate problem with effective currents for auxiliary problems
      class(tSolution) :: this;
      integer :: i, j, k;
      integer, intent(in) :: t;
      real*8 :: mex;

      ! Main problem part
      call cpu_time(tstarts(5))
      !$ tstarts(5) = OMP_get_wtime(); 
      call this%mainproblem%poisson%Ecmpst1_to_Ecmpst0;
   call cpu_time(tends(5))
      !$ tends(5) = OMP_get_wtime();

   call cpu_time(tstarts(18))
      !$ tstarts(18) = OMP_get_wtime();  
      call this%mainproblem%saveFields(0);
   call cpu_time(tends(18))
      !$ tends(18) = OMP_get_wtime(); 

   call cpu_time(tstarts(0))
      !$ tstarts(0) = OMP_get_wtime();    
      call this%mainproblem%getSourceCurrents(t);          ! Fill J arrays with actual values of currents
   call cpu_time(tends(0))
      !$ tends(0) = OMP_get_wtime(); 

   call cpu_time(tstarts(1))
      !$ tstarts(1) = OMP_get_wtime();
      call this%mainproblem%updateEInterior;               ! Calculate E field everywhere except outer boundary
   call cpu_time(tends(1))
     !$ tends(1) = OMP_get_wtime();
   
      !Poisson correction!
   call cpu_time(tstarts(5))
      !$ tstarts(5) = OMP_get_wtime();   
      call this%mainproblem%UploadFieldsToPoisson;
      call this%mainproblem%poisson%Esolver(t);
      call this%mainproblem%poisson%gradientEcalc(t);
      call this%mainproblem%poisson%Ecpmst1_formation(t);
      call this%mainproblem%addEcompositesToBuffer;
      call this%mainproblem%JEeffBuild(t);
   call cpu_time(tends(5))
   !$ tends(5) = OMP_get_wtime();

   call cpu_time(tstarts(2))
      !$ tstarts(2) = OMP_get_wtime();
      call this%mainproblem%updateHInterior;               ! Calculate H field everywhere except outer boundary
   call cpu_time(tends(2))
      !$ tends(2) = OMP_get_wtime();

   call cpu_time(tstarts(5))
      !$ tstarts(5) = OMP_get_wtime();
      call this%mainproblem%JHeffBuild(t);
      call this%mainproblem%saveEffCurrentsToBuffer;
   call cpu_time(tends(5))
   !$ tends(5) = OMP_get_wtime();
   
      ! Aux problem part
      do k=0,auxmaxcount-1
         if (this%timers%isused(k)==1) then
            call this%auxproblems(k)%problem_DoIndependentStep(t);
            call cpu_time(tstarts(3))
            !$ tstarts(3) = OMP_get_wtime();
                 call this%auxproblems(k)%addEtobuffer;
            call cpu_time(tends(3))
            !$ tends(3) = OMP_get_wtime(); 
         endif
      enddo;
      ! Main problem finish
    call cpu_time(tstarts(3))
      !$ tstarts(3) = OMP_get_wtime();   
       call this%mainproblem%fillEBoundary(t);              ! Fill E boundary with current Eboundary condition
    call cpu_time(tends(3))
    !$ tends(3) = OMP_get_wtime();

    call cpu_time(tstarts(4))
      !$ tstarts(4) = OMP_get_wtime();
      call this%mainproblem%fillHBoundary(t);              ! Fill H boundary with current Hboundary condition
    call cpu_time(tends(4))
      !$ tends(4) = OMP_get_wtime();

    call cpu_time(tstarts(7))
         !$ tstarts(7) = OMP_get_wtime();
         call this%solution_WriteCheckInfo(t);
    call cpu_time(tends(7))
    !$ tends(7) = OMP_get_wtime();

    call cpu_time(tstarts(6))
      !$ tstarts(6) = OMP_get_wtime();
      call this%ManageAuxProblems(t);
     call cpu_time(tends(6))
     !$ tends(6) = OMP_get_wtime();
     do k=0,31
        times(k) = times(k) + (tends(k)-tstarts(k));
        tstarts(k) = 0.0d0;
        tends(k) = 0.0d0;
      enddo 
    end subroutine PropagateLacunaesPoisson;


!PropagateLacunasNoEffCurrentsNoPoisson----------------------------------------------------------    
    subroutine PropagateLacunasNoEffCurrentsNoPoisson(this,t)
    ! Propagate problem with effective currents for auxiliary problems
      class(tSolution) :: this;
      integer, intent(in) :: t;
      integer :: i, j, k;
      real*8 :: mex;
      call this%buffer%clearFields;
      ! Main problem part
      call this%mainproblem%saveFields(0);


      call this%mainproblem%getSourceCurrents(t);          ! Fill J arrays with actual values of currents
  
    call cpu_time(tstarts(1))
      !$ tstarts(1) = OMP_get_wtime();
      call this%mainproblem%updateEInterior;               ! Calculate E field everywhere except outer boundary        
    call cpu_time(tends(1))
      !$ tends(1) = OMP_get_wtime();

    call cpu_time(tstarts(2))
      !$ tstarts(2) = OMP_get_wtime();
      call this%mainproblem%updateHInterior;               ! Calculate H field everywhere except outer boundary
    call cpu_time(tends(2))
      !$ tends(2) = OMP_get_wtime(); 
    
      ! Aux problem part
      do k=0,auxmaxcount-1
         if (this%timers%isused(k)==1) then
            call this%auxproblems(k)%problem_DoIndependentStep(t);
            call this%auxproblems(k)%addEtobuffer;
         endif
      enddo;
      ! Main problem finish
      
    call cpu_time(tstarts(3))
      !$ tstarts(3) = OMP_get_wtime();        
      call this%mainproblem%fillEBoundary(t);              ! Fill E boundary with current Eboundary condition
    call cpu_time(tends(3))
      !$ tends(3) = OMP_get_wtime(); 

    call cpu_time(tstarts(4))
      !$ tstarts(4) = OMP_get_wtime();
      call this%mainproblem%fillHBoundary(t);              ! Fill H boundary with current Hboundary condition
    call cpu_time(tends(4))
           !$ tends(4) = OMP_get_wtime();
    call this%solution_WriteCheckInfo(t)

    call cpu_time(tstarts(6))
      !$ tstarts(6) = OMP_get_wtime();
      call this%ManageAuxProblems(t);
    call cpu_time(tends(6))
           !$ tends(6) = OMP_get_wtime();
    
      do k=0,31
         times(k) = times(k) + (tends(k)-tstarts(k));
         tstarts(k) = 0.0d0;
         tends(k) = 0.0d0;
      enddo  
    end subroutine PropagateLacunasNoEffCurrentsNoPoisson;  


!PropagateLacunasNoEffCurrentsNoPoissonOneSource----------------------------------------------------------    
    subroutine PropagateLacunasNoEffCurrentsNoPoissonOneSource(this,t)
    ! Propagate problem with effective currents for auxiliary problems
      class(tSolution) :: this;
      integer, intent(in) :: t;
      integer :: i, j, k;
      real*8 :: mex;
      call cpu_time(tstarts(18))
      !$ tstarts(18) = OMP_get_wtime(); 
          call this%mainproblem%saveFields(0);
      call cpu_time(tends(18))
      !$ tends(18) = OMP_get_wtime(); 
      
    call cpu_time(tstarts(0))
      !$ tstarts(0) = OMP_get_wtime();  
      call this%mainproblem%getSourceCurrents(t);             ! Fill J arrays with actual values of currents
    call cpu_time(tends(0))
      !$ tends(0) = OMP_get_wtime();

    call cpu_time(tstarts(10))
      !$ tstarts(10) = OMP_get_wtime();  
      call this%mainproblem%saveOriginalCurrentsToBuffer(); !!! This is the main difference
    call cpu_time(tends(10))
      !$ tends(10) = OMP_get_wtime();

    call cpu_time(tstarts(1))
      !$ tstarts(1) = OMP_get_wtime();
      call this%mainproblem%updateEInterior;                  ! Calculate E field everywhere except outer boundary        
    call cpu_time(tends(1))
      !$ tends(1) = OMP_get_wtime();

    call cpu_time(tstarts(2))
      !$ tstarts(2) = OMP_get_wtime();
      call this%mainproblem%updateHInterior;                  ! Calculate H field everywhere except outer boundary
    call cpu_time(tends(2))
      !$ tends(2) = OMP_get_wtime(); 

      ! Aux problem part
      do k=0,auxmaxcount-1
         if (this%timers%isused(k)==1) then
            call this%auxproblems(k)%problem_DoIndependentStep(t);
            call cpu_time(tstarts(3))
            !$ tstarts(3) = OMP_get_wtime();
                call this%auxproblems(k)%addEtobuffer;
            call cpu_time(tends(3))
            !$ tends(3) = OMP_get_wtime(); 
         endif
      enddo;
      ! Main problem finish     

    call cpu_time(tstarts(3))
      !$ tstarts(3) = OMP_get_wtime();      
      call this%mainproblem%fillEBoundary(t);              ! Fill E boundary with current Eboundary condition
    call cpu_time(tends(3))
      !$ tends(3) = OMP_get_wtime(); 

    call cpu_time(tstarts(4))
      !$ tstarts(4) = OMP_get_wtime();
      call this%mainproblem%fillHBoundary(t);              ! Fill H boundary with current Hboundary condition
    call cpu_time(tends(4))
           !$ tends(4) = OMP_get_wtime();

    call cpu_time(tstarts(7))
         !$ tstarts(7) = OMP_get_wtime();
         call this%solution_WriteCheckInfo(t);
    call cpu_time(tends(7))
    !$ tends(7) = OMP_get_wtime();

    call cpu_time(tstarts(6))
      !$ tstarts(6) = OMP_get_wtime();
       call this%ManageAuxProblems(t);
    call cpu_time(tends(6))
           !$ tends(6) = OMP_get_wtime();

    do k=0,31
       times(k) = times(k) + (tends(k)-tstarts(k));
       tstarts(k) = 0.0d0;
       tends(k) = 0.0d0;
    enddo  
    end subroutine PropagateLacunasNoEffCurrentsNoPoissonOneSource;      

!PropagateEffCurrentsMANYAuxProblemsPoisson----------------------------------------------------------    
    subroutine PropagateLacunaesPoissonTiming(this,t)
    ! Propagate problem with effective currents for auxiliary problems
      class(tSolution) :: this;
      integer :: i, j, k;
      integer, intent(in) :: t;
      real*8 :: mex;
         if (checkupdated==1) then
           call this%mainproblem%Ef%mesh_clearChA;           ! Clear checking arrays
           call this%mainproblem%Hf%mesh_clearChA;           ! Clear checking arrays
         endif

       call cpu_time(tstarts(0))
       !$ tstarts(0) = OMP_get_wtime();
          call this%buffer%clearFields;
       call cpu_time(tends(0))
       !$ tends(0) = OMP_get_wtime();  

       ! Main problem part
       
       call cpu_time(tstarts(1))
       !$ tstarts(1) = OMP_get_wtime();
          call this%mainproblem%poisson%Ecmpst1_to_Ecmpst0;      
       call cpu_time(tends(1))
       !$ tends(1) = OMP_get_wtime();  

       call cpu_time(tstarts(2))
       !$ tstarts(2) = OMP_get_wtime();
          call this%mainproblem%saveFields(0);
       call cpu_time(tends(2))
       !$ tends(2) = OMP_get_wtime(); 
       
       if ((mod(t,file_erroroutputsteps)==0)) then
          call cpu_time(tstarts(3))
          !$ tstarts(3) = OMP_get_wtime();
             call this%mainproblem%getAnalyticSolution(t)
          call cpu_time(tends(3))
          !$ tends(3) = OMP_get_wtime(); 
       endif

       call cpu_time(tstarts(4))
       !$ tstarts(4) = OMP_get_wtime();
          call this%mainproblem%getSourceCurrents(t);          ! Fill J arrays with actual values of currents
       call cpu_time(tends(4))
       !$ tends(4) = OMP_get_wtime();

       call cpu_time(tstarts(5))
       !$ tstarts(5) = OMP_get_wtime();
          call this%mainproblem%updateEInterior;               ! Calculate E field everywhere except outer boundary
       call cpu_time(tends(5))
       !$ tends(5) = OMP_get_wtime();
   
             
       !Poisson correction!
          call cpu_time(tstarts(6))
          !$ tstarts(6) = OMP_get_wtime();
             call this%mainproblem%UploadFieldsToPoisson;
          call cpu_time(tends(6))
          !$ tends(6) = OMP_get_wtime();

          call cpu_time(tstarts(7))
          !$ tstarts(7) = OMP_get_wtime();
             call this%mainproblem%poisson%Esolver(t);
          call cpu_time(tends(7))
          !$ tends(7) = OMP_get_wtime();

          call cpu_time(tstarts(8))
          !$ tstarts(8) = OMP_get_wtime();
             call this%mainproblem%poisson%gradientEcalc(t);
          call cpu_time(tends(8))
          !$ tends(8) = OMP_get_wtime();

          call cpu_time(tstarts(9))
          !$ tstarts(9) = OMP_get_wtime();
             call this%mainproblem%poisson%Ecpmst1_formation(t);
          call cpu_time(tends(9))
          !$ tends(9) = OMP_get_wtime();

          call cpu_time(tstarts(10))
          !$ tstarts(10) = OMP_get_wtime();
             call this%mainproblem%addEcompositesToBuffer;
          call cpu_time(tends(10))
          !$ tends(10) = OMP_get_wtime();

          call cpu_time(tstarts(11))
          !$ tstarts(11) = OMP_get_wtime();
             call this%mainproblem%JEeffBuild(t);
          call cpu_time(tends(11))
          !$ tends(11) = OMP_get_wtime();

          call cpu_time(tstarts(12))
          !$ tstarts(12) = OMP_get_wtime();
             call this%mainproblem%updateHInterior;               ! Calculate H field everywhere except outer boundary
          call cpu_time(tends(12))
          !$ tends(12) = OMP_get_wtime();

          call cpu_time(tstarts(13))
          !$ tstarts(13) = OMP_get_wtime();
             call this%mainproblem%JHeffBuild(t);
          call cpu_time(tends(13))
          !$ tends(13) = OMP_get_wtime();

          call cpu_time(tstarts(14))
          !$ tstarts(14) = OMP_get_wtime();
             call this%mainproblem%saveEffCurrentsToBuffer;
          call cpu_time(tends(14))
          !$ tends(14) = OMP_get_wtime();

        ! Aux problem part
       
        do k=0,auxmaxcount-1
           if (this%timers%isused(k)==1) then
              call this%auxproblems(k)%problem_DoIndependentStepTiming(t);
        
              call cpu_time(tstarts(16))
              !$ tstarts(16) = OMP_get_wtime();
                 call this%auxproblems(k)%addEtobuffer;
              call cpu_time(tends(16))
              !$ tends(16) = OMP_get_wtime();
           endif
        enddo;
      
        ! Main problem finish

        call cpu_time(tstarts(17))
        !$ tstarts(17) = OMP_get_wtime();
           call this%mainproblem%fillEBoundary(t);              ! Fill E boundary with current Eboundary condition
        call cpu_time(tends(17))
        !$ tends(17) = OMP_get_wtime();

        call cpu_time(tstarts(18))
        !$ tstarts(18) = OMP_get_wtime();
           call this%mainproblem%fillHBoundary(t);              ! Fill H boundary with current Hboundary condition
        call cpu_time(tends(18))
        !$ tends(18) = OMP_get_wtime();

         if ((mod(t,screen_erroroutputsteps)==0).OR.(mod(t,file_erroroutputsteps)==0)) then
            call cpu_time(tstarts(19))
            !$ tstarts(19) = OMP_get_wtime();
               call this%mainproblem%getCurrentError(t);
            call cpu_time(tends(19))
            !$ tends(19) = OMP_get_wtime();
         endif

         call cpu_time(tstarts(20))
         !$ tstarts(20) = OMP_get_wtime();
            call this%solution_WriteOutput(t)
         call cpu_time(tends(20))
         !$ tends(20) = OMP_get_wtime();

       !  call cpu_time(tstarts(21))
       !  !$ tstarts(21) = OMP_get_wtime();
       !     call this%solution_WriteCheckInfo(t)
       !  call cpu_time(tends(21))
       !  !$ tends(21) = OMP_get_wtime();
         
       !  if (checkupdated==1) then
       !     call this%mainproblem%Ef%mesh_checkarrays;        ! Check all the points are updated once
       !     call this%mainproblem%Hf%mesh_checkarrays;        ! Check all the points are updated once
       !  endif

         call cpu_time(tstarts(23))
         !$ tstarts(23) = OMP_get_wtime();
            call this%ManageAuxProblems(t);
         call cpu_time(tends(23))
         !$ tends(23) = OMP_get_wtime();
         
         call cpu_time(tstarts(22))
         !$ tstarts(22) = OMP_get_wtime();
          !  call this%solution_printstep(t);
         call cpu_time(tends(22))
         !$ tends(22) = OMP_get_wtime();

         do k=0,31
            times(k) = times(k) + (tends(k)-tstarts(k));
         enddo         
      call this%solution_WriteTimes;
    end subroutine PropagateLacunaesPoissonTiming;


!PropagatePurePML----------------------------------------------------------    
    subroutine propagatePurePML(this,t)
    ! Propagate single main problem
      class(tSolution) :: this;
      integer :: i, j, k;
      integer, intent(in) :: t;
      call this%mainproblem%getSourceCurrentsPML(t); 
      call this%mainproblem%PML%updateEInteriorPMLEverywhere;               
      call this%mainproblem%PML%updateHInteriorPMLEverywhere;
    end subroutine propagatePurePML;


!PropagatePurePMLwithABC----------------------------------------------------    
    subroutine propagatePurePMLwithABC(this,t)
    ! Propagate single main problem
      class(tSolution) :: this;
      integer :: i, j, k;
      integer, intent(in) :: t;
      call this%mainproblem%getSourceCurrentsPML(t); 
      call this%mainproblem%PML%updateEInteriorPMLEverywhere;               
      call this%mainproblem%PML%updateHInteriorPMLEverywhere;

      call this%mainproblem%getSourceCurrents(t);          ! Fill J arrays with actual values of currents
      call this%mainproblem%saveFields(1);
      call this%mainproblem%saveFields(0);
         
      if (this%mainproblem%EboundaryType == 9.OR.this%mainproblem%EboundaryType == 10.OR.this%mainproblem%EboundaryType == 11.OR.this%mainproblem%EboundaryType == 12) then
         call this%mainproblem%saveGivoliNetaAuxVariables
      endif;
      if (this%mainproblem%EboundaryType == 14) then
         call this%mainproblem%saveHWAuxVariables
      endif
      if (this%mainproblem%EboundaryType == 16) then
         call this%mainproblem%saveHWEvaAuxVariables
      endif
         
      call this%mainproblem%updateEInterior;               ! Calculate E field everywhere except outer boundary
      call this%mainproblem%updateHInterior;               ! Calculate H field everywhere except outer boundary
      select case (this%mainproblem%EboundaryType)
      case (10)   
         call this%mainproblem%fillEBoundaryABCGivoliNetaPML(t);
      case (14)
         call this%mainproblem%fillEBoundaryABCHWPML(t);
      case (16)
         call this%mainproblem%fillEBoundaryABCHWEvaPML(t);      
      end select   
      call this%mainproblem%fillHBoundary(t);
    end subroutine propagatePurePMLwithABC;    
    
    
!solution_WriteTimes----------------------------------------------------------    
    subroutine solution_WriteTimes(this)
    ! Write how much time each function performs
      class(tSolution) :: this;
      real*8 :: alltimes;
      integer :: k;
      alltimes = 0;
      do k=0,31
         alltimes = alltimes+times(k);
      enddo
      write(*,fmt=6455) 'buffer%clearFields' , times(0), times(0)/alltimes;
      write(*,fmt=6455) 'mainproblem%poisson%Ecmpst1_to_Ecmpst0', times(1), times(1)/alltimes;
      write(*,fmt=6455) 'mainproblem%saveFields' , times(2), times(2)/alltimes;
      write(*,fmt=6455) 'mainproblem%getAnalyticSolution(t)' , times(3), times(3)/alltimes;
      write(*,fmt=6455) 'mainproblem%getSourceCurrents(t)' , times(4), times(4)/alltimes;
      write(*,fmt=6455) 'mainproblem%updateEInterior' , times(5), times(5)/alltimes;
      write(*,fmt=6455) 'mainproblem%UploadFieldsToPoisson' , times(6), times(6)/alltimes;
      write(*,fmt=6455) 'mainproblem%poisson%Esolver(t)' , times(7), times(7)/alltimes;
      write(*,fmt=6455) 'mainproblem%poisson%gradientEcalc(t)' , times(8), times(8)/alltimes;
      write(*,fmt=6455) 'mainproblem%poisson%Ecpmst1_formation(t)' , times(9), times(9)/alltimes;
      write(*,fmt=6455) 'mainproblem%addEcompositesToBuffer' , times(10), times(10)/alltimes;
      write(*,fmt=6455) 'mainproblem%JEeffBuild(t)' , times(11), times(11)/alltimes;
      write(*,fmt=6455) 'mainproblem%updateHInterior' , times(12), times(12)/alltimes;
      write(*,fmt=6455) 'mainproblem%JHeffBuild(t)' , times(13), times(13)/alltimes;
      write(*,fmt=6455) 'mainproblem%saveEffCurrentsToBuffer' , times(14), times(14)/alltimes;

      write(*,fmt=6455) 'auxproblems(k)%getSourceCurrents' , times(24), times(24)/alltimes;
      write(*,fmt=6455) 'auxproblems(k)%saveFields' , times(25), times(25)/alltimes;
      write(*,fmt=6455) 'auxproblems(k)%updateEInterior' , times(26), times(26)/alltimes;
      write(*,fmt=6455) 'auxproblems(k)%updateHInterior' , times(27), times(27)/alltimes;
      write(*,fmt=6455) 'auxproblems(k)%uploadFieldsToPML' , times(28), times(28)/alltimes;
      write(*,fmt=6455) 'auxproblems(k)%PML%updateEinterior' , times(15), times(15)/alltimes;
      write(*,fmt=6455) 'auxproblems(k)%PML%updateHinterior' , times(29), times(29)/alltimes;

      write(*,fmt=6455) 'auxproblems(k)%fillEBoundary(t)' , times(30), times(30)/alltimes;
      write(*,fmt=6455) 'auxproblems(k)%fillHBoundary(t)' , times(31), times(31)/alltimes;

      write(*,fmt=6455) 'auxproblems(k)%addEtobuffer' , times(16), times(16)/alltimes;
      write(*,fmt=6455) 'mainproblem%fillEBoundary(t)' , times(17), times(17)/alltimes;
      write(*,fmt=6455) 'mainproblem%fillHBoundary(t)' , times(18), times(18)/alltimes;
      write(*,fmt=6455) 'mainproblem%getCurrentError' , times(19), times(19)/alltimes;
      write(*,fmt=6455) 'solution_WriteOutput(t)' , times(20), times(20)/alltimes;
      write(*,fmt=6455) 'solution_WriteCheckInfo(t)' , times(21), times(21)/alltimes;
      write(*,fmt=6455) 'solution_printstep(t)' , times(22), times(22)/alltimes;
      write(*,fmt=6455) 'ManageAuxProblems(t)' , times(23), times(23)/alltimes;   
      6455 format(A50,',',E22.15,',',E22.15);
    end subroutine solution_WriteTimes
    

    
!solution_WriteCheckInfo----------------------------------------------------------    
    subroutine solution_WriteCheckInfo(this, t)
    ! Build the solution
      class(tSolution) :: this;
      integer, intent(in) :: t;
      integer :: i,k;
      if ((mod(t,file_auxproblemsoutputstep)==0)) then
         if (file_auxproblemsoutput>0) then
            do i=0,9
               if (this%timers%isused(i)==1.AND.this%auxproblems(i)%problem_id<file_auxproblemsoutput) then
                  k = this%auxproblems(i)%problem_id;
                  if (file_auxproblemswritelocations==1) then
                     write(unit=1000+k*100,   fmt=10) maxval(abs(this%auxproblems(i)%Ef%X)), maxloc(abs(this%auxproblems(i)%Ef%X));
                     write(unit=1000+k*100+1, fmt=10) maxval(abs(this%auxproblems(i)%Ef%Y)), maxloc(abs(this%auxproblems(i)%Ef%Y));
                     write(unit=1000+k*100+2, fmt=10) maxval(abs(this%auxproblems(i)%Ef%Z)), maxloc(abs(this%auxproblems(i)%Ef%Z));
                     write(unit=1000+k*100+3, fmt=10) maxval(abs(this%auxproblems(i)%Hf%X)), maxloc(abs(this%auxproblems(i)%Hf%X));
                     write(unit=1000+k*100+4, fmt=10) maxval(abs(this%auxproblems(i)%Hf%Y)), maxloc(abs(this%auxproblems(i)%Hf%Y));
                     write(unit=1000+k*100+5, fmt=10) maxval(abs(this%auxproblems(i)%Hf%Z)), maxloc(abs(this%auxproblems(i)%Hf%Z));
                     write(unit=1000+k*100+6, fmt=10) maxval(abs(this%auxproblems(i)%JE%X)), maxloc(abs(this%auxproblems(i)%JE%X));
                     write(unit=1000+k*100+7, fmt=10) maxval(abs(this%auxproblems(i)%JE%Y)), maxloc(abs(this%auxproblems(i)%JE%Y));
                     write(unit=1000+k*100+8, fmt=10) maxval(abs(this%auxproblems(i)%JE%Z)), maxloc(abs(this%auxproblems(i)%JE%Z));
                     write(unit=1000+k*100+9, fmt=10) maxval(abs(this%auxproblems(i)%JH%X)), maxloc(abs(this%auxproblems(i)%JH%X));
                     write(unit=1000+k*100+10, fmt=10) maxval(abs(this%auxproblems(i)%JH%Y)), maxloc(abs(this%auxproblems(i)%JH%Y));
                     write(unit=1000+k*100+11, fmt=10) maxval(abs(this%auxproblems(i)%JH%Z)), maxloc(abs(this%auxproblems(i)%JH%Z));
                     write(unit=1000+k*100+12, fmt=20) this%auxproblems(i)%thetabig(this%auxproblems(i)%Ti05(t)-(1.0+this%auxproblems(i)%sgm_s)*this%auxproblems(i)%auxmaxtime*this%auxproblems(i)%problem_id)
                     write(unit=2000+k*100,   fmt=10) maxval(abs(this%auxproblems(i)%PML%Ef%X)), maxloc(abs(this%auxproblems(i)%PML%Ef%X));
                     write(unit=2000+k*100+1, fmt=10) maxval(abs(this%auxproblems(i)%PML%Ef%Y)), maxloc(abs(this%auxproblems(i)%PML%Ef%Y));
                     write(unit=2000+k*100+2, fmt=10) maxval(abs(this%auxproblems(i)%PML%Ef%Z)), maxloc(abs(this%auxproblems(i)%PML%Ef%Z));
                     write(unit=2000+k*100+3, fmt=10) maxval(abs(this%auxproblems(i)%PML%Hf%X)), maxloc(abs(this%auxproblems(i)%PML%Hf%X));
                     write(unit=2000+k*100+4, fmt=10) maxval(abs(this%auxproblems(i)%PML%Hf%Y)), maxloc(abs(this%auxproblems(i)%PML%Hf%Y));
                     write(unit=2000+k*100+5, fmt=10) maxval(abs(this%auxproblems(i)%PML%Hf%Z)), maxloc(abs(this%auxproblems(i)%PML%Hf%Z));
                     write(unit=3000+k*100,   fmt=10) maxval(abs(this%auxproblems(i)%PML%JfE%X)), maxloc(abs(this%auxproblems(i)%PML%JfE%X));
                     write(unit=3000+k*100+1, fmt=10) maxval(abs(this%auxproblems(i)%PML%JfE%Y)), maxloc(abs(this%auxproblems(i)%PML%JfE%Y));
                     write(unit=3000+k*100+2, fmt=10) maxval(abs(this%auxproblems(i)%PML%JfE%Z)), maxloc(abs(this%auxproblems(i)%PML%JfE%Z));
                     write(unit=3000+k*100+3, fmt=10) maxval(abs(this%auxproblems(i)%PML%JfH%X)), maxloc(abs(this%auxproblems(i)%PML%JfH%X));
                     write(unit=3000+k*100+4, fmt=10) maxval(abs(this%auxproblems(i)%PML%JfH%Y)), maxloc(abs(this%auxproblems(i)%PML%JfH%Y));
                     write(unit=3000+k*100+5, fmt=10) maxval(abs(this%auxproblems(i)%PML%JfH%Z)), maxloc(abs(this%auxproblems(i)%PML%JfH%Z));
                  else
                     write(unit=1000+k*100,   fmt=20) maxval(abs(this%auxproblems(i)%Ef%X));
                     write(unit=1000+k*100+1, fmt=20) maxval(abs(this%auxproblems(i)%Ef%Y));
                     write(unit=1000+k*100+2, fmt=20) maxval(abs(this%auxproblems(i)%Ef%Z));
                     write(unit=1000+k*100+3, fmt=20) maxval(abs(this%auxproblems(i)%Hf%X));
                     write(unit=1000+k*100+4, fmt=20) maxval(abs(this%auxproblems(i)%Hf%Y));
                     write(unit=1000+k*100+5, fmt=20) maxval(abs(this%auxproblems(i)%Hf%Z));
                     write(unit=1000+k*100+6, fmt=20) maxval(abs(this%auxproblems(i)%JE%X));
                     write(unit=1000+k*100+7, fmt=20) maxval(abs(this%auxproblems(i)%JE%Y));
                     write(unit=1000+k*100+8, fmt=20) maxval(abs(this%auxproblems(i)%JE%Z));
                     write(unit=1000+k*100+9, fmt=20) maxval(abs(this%auxproblems(i)%JH%X));
                     write(unit=1000+k*100+10, fmt=20) maxval(abs(this%auxproblems(i)%JH%Y));
                     write(unit=1000+k*100+11, fmt=20) maxval(abs(this%auxproblems(i)%JH%Z));
                     write(unit=1000+k*100+12, fmt=20) this%auxproblems(i)%thetabig(this%auxproblems(i)%Ti05(t)-(1.0+this%auxproblems(i)%sgm_s)*this%auxproblems(i)%auxmaxtime*this%auxproblems(i)%problem_id)
                     write(unit=2000+k*100,   fmt=20) maxval(abs(this%auxproblems(i)%PML%Ef%X));
                     write(unit=2000+k*100+1, fmt=20) maxval(abs(this%auxproblems(i)%PML%Ef%Y));
                     write(unit=2000+k*100+2, fmt=20) maxval(abs(this%auxproblems(i)%PML%Ef%Z));
                     write(unit=2000+k*100+3, fmt=20) maxval(abs(this%auxproblems(i)%PML%Hf%X));
                     write(unit=2000+k*100+4, fmt=20) maxval(abs(this%auxproblems(i)%PML%Hf%Y));
                     write(unit=2000+k*100+5, fmt=20) maxval(abs(this%auxproblems(i)%PML%Hf%Z));
                     write(unit=3000+k*100,   fmt=20) maxval(abs(this%auxproblems(i)%PML%JfE%X));
                     write(unit=3000+k*100+1, fmt=20) maxval(abs(this%auxproblems(i)%PML%JfE%Y));
                     write(unit=3000+k*100+2, fmt=20) maxval(abs(this%auxproblems(i)%PML%JfE%Z));
                     write(unit=3000+k*100+3, fmt=20) maxval(abs(this%auxproblems(i)%PML%JfH%X));
                     write(unit=3000+k*100+4, fmt=20) maxval(abs(this%auxproblems(i)%PML%JfH%Y));
                     write(unit=3000+k*100+5, fmt=20) maxval(abs(this%auxproblems(i)%PML%JfH%Z));
                  endif
               endif
            enddo
         endif
      endif   
      10 format(E22.15,',',I3,',',I3,',',I3);
      20 format(E22.15); 
    end subroutine solution_WriteCheckInfo;    

!solution_PrintStep----------------------------------------------------------    
    subroutine solution_PrintStep(this, t, ttime)
    ! Print step information
      class(tSolution) :: this;
      integer, intent(in) :: t;
      integer :: k;
      real*8 :: ttime
      integer :: ds, hs, mn, sc, ms;
      if (ttime<60.0d0) then
         write(*,100) t, this%Nt, INT(ttime);
      elseif (ttime>=60.0d0.and.ttime<3600.0d0) then
         mn = AINT(ttime/60.0d0);
         sc = AINT(ttime - mn*60.0d0);
         write(*,101) t, this%Nt, mn, sc
      elseif (ttime>=3600.0d0.and.ttime<86400.0d0) then
         hs = AINT(ttime/3600.0d0);
         mn = AINT( (ttime - hs*3600.0d0)/60.0d0 );
         sc = AINT( (ttime - hs*3600.0d0 - mn*60.0d0) );
         write(*,102) t, this%Nt, hs, mn, sc
      elseif (ttime>=86400.0d0) then
         ds = AINT(ttime/86400.0d0);
         hs = AINT((ttime - ds*86400.0d0)/3600.0d0);
         mn = AINT( (ttime - ds*86400.0d0 - hs*3600.0d0)/60.0d0 );
         sc = AINT(ttime - ds*86400.0d0 - hs*3600.0d0 - mn*60.0d0);
         write(*,103) t, this%Nt, ds, hs, mn, sc
      endif
      100 format ('Step',I6,1x,'|',1x,I6,1x,'|',1x,I2,'s');
      101 format ('Step',I6,1x,'|',1x,I6,1x,'|',1x,I2,'m',1x,I2,'s');
      102 format ('Step',I6,1x,'|',1x,I6,1x,'|',1x,I2,'h',1x,I2,'m',1x,I2,'s');
      103 format ('Step',I6,1x,'|',1x,I6,1x,'|',1x,I2,'d',1x,I2,'h',1x,I2,'m',1x,I2,'s'); 
      select case (this%soltype)
      case(0) ! Analytical E boundaries
         call this%mainproblem%problem_PrintCurrentErrors(t);
         !call this%mainproblem%problem_PrintMaximumFields(t);
         !call this%mainproblem%problem_PrintMaximumCurrents(t);
      case(1) ! PML E boundaries
         call this%mainproblem%problem_PrintCurrentErrors(t);
      case (100) ! Unsplit PML, Lacunaes, the main source is cut, NO Poisson
         call this%mainproblem%problem_PrintCurrentErrors(t);
         !do k=0,auxmaxcount-1
         !   if (this%timers%isused(k)==1) then
         !      write(*,fmt=50) k, this%timers%isused(k), this%auxproblems(k)%problem_id, this%timers%timer(k)
         !   endif 
         !enddo        
      case(101) ! Pure Lacunas, the main source is cut, NO Poisson
         call this%mainproblem%problem_PrintCurrentErrors(t);
      case(102) ! Pure Lacunas, the main source is cut, NO Poisson, Original Sources
        call this%mainproblem%problem_PrintCurrentErrors(t);      
      case(2) ! One auxiliary problem and analytical E bounds
        call this%mainproblem%problem_PrintCurrentErrors(t);
      case(3) ! One auxiliary problem and PML E bounds
        call this%mainproblem%problem_PrintCurrentErrors(t);
      case(4) ! One auxiliary problem, Effective currents and PML E bounds
        call this%mainproblem%problem_PrintCurrentErrors(t);
      case(5) ! Quasi-lacunaes
         call this%mainproblem%problem_PrintCurrentErrors(t);
         do k=0,auxmaxcount-1
            if (this%timers%isused(k)==1) then
               write(*,fmt=50) k, this%timers%isused(k), this%auxproblems(k)%problem_id, this%timers%timer(k)
            !else
            !   write(*,fmt=51) k, this%timers%isused(k), this%timers%timer(k)
            endif
         enddo
      case(6) ! Lacunaes with Poisson
         call this%mainproblem%problem_PrintCurrentErrors(t);
         do k=0,auxmaxcount-1
            if (this%timers%isused(k)==1) then
               write(*,fmt=50) k, this%timers%isused(k), this%auxproblems(k)%problem_id, this%timers%timer(k)
            !else
            !   write(*,fmt=51) k, this%timers%isused(k), this%timers%timer(k)
            endif
         enddo
      case(7) ! Pure PML
         call this%mainproblem%problem_PrintPMLCurrentErrors(t);
      case(8) ! Sommerfeld ABC
         call this%mainproblem%problem_PrintCurrentErrors(t);
      case(9) ! Higdon ABC
         call this%mainproblem%problem_PrintCurrentErrors(t);
      case(10) ! Betz-Mittra ABC
         call this%mainproblem%problem_PrintCurrentErrors(t);
      case(11) ! Mur ABC
         call this%mainproblem%problem_PrintCurrentErrors(t);
      case(12) ! Givoli-Neta ABC X=0 with exact elsewhere
         call this%mainproblem%problem_PrintCurrentErrors(t);
      case(13) ! Givoli-Neta ABC X=0 with UnsplitPML elsewhere
         call this%mainproblem%problem_PrintCurrentErrors(t);   
      case(14) ! Givoli-Neta and one big auxiliary problem without boundaries
         call this%mainproblem%problem_PrintCurrentErrors(t);
      case(15) ! Givoli-Neta ABC + Lacunaes with Poisson
         call this%mainproblem%problem_PrintCurrentErrors(t);
         do k=0,auxmaxcount-1
            if (this%timers%isused(k)==1) then
               write(*,fmt=50) k, this%timers%isused(k), this%auxproblems(k)%problem_id, this%timers%timer(k)
            !else
            !   write(*,fmt=51) k, this%timers%isused(k), this%timers%timer(k)
            endif
         enddo   
      case(16) ! One big auxiliary problem without boundaries
         call this%mainproblem%problem_PrintCurrentErrors(t);
         do k=0,auxmaxcount-1
            if (this%timers%isused(k)==1) then
               write(*,fmt=50) k, this%timers%isused(k), this%auxproblems(k)%problem_id, this%timers%timer(k)
            !else
            !   write(*,fmt=51) k, this%timers%isused(k), this%timers%timer(k)
            endif
         enddo    
      case(17) ! Hagstrom-Warburton ABC with exact elsewhere
         call this%mainproblem%problem_PrintCurrentErrors(t);
      case(18) ! Hagstrom-Warburton ABC with UnsplitPML elsewhere
         call this%mainproblem%problem_PrintCurrentErrors(t);
      case(19) ! Hagstrom-Warburton ABC + Lacunaes with Poisson
         call this%mainproblem%problem_PrintCurrentErrors(t);
         do k=0,auxmaxcount-1
            if (this%timers%isused(k)==1) then
               write(*,fmt=50) k, this%timers%isused(k), this%auxproblems(k)%problem_id, this%timers%timer(k)
            !else
            !   write(*,fmt=51) k, this%timers%isused(k), this%timers%timer(k)
            endif
         enddo
      case(20) ! HWEva ABC Ey_X=Nx with exact elsewhere
         call this%mainproblem%problem_PrintCurrentErrors(t);
      case(21) ! HWEva ABC Ey_X=Nx with Unsplit PML elsewhere
         call this%mainproblem%problem_PrintCurrentErrors(t);
      case(22) ! HWEva ABC Ey_x=Nx + PML elsewhere with Lacunaes, the main source is cut by Theta, NO Poisson
         call this%mainproblem%problem_PrintCurrentErrors(t);
         do k=0,auxmaxcount-1
            if (this%timers%isused(k)==1) then
               write(*,fmt=50) k, this%timers%isused(k), this%auxproblems(k)%problem_id, this%timers%timer(k)
            !else
            !   write(*,fmt=51) k, this%timers%isused(k), this%timers%timer(k)
            endif
         enddo
      case(23) ! HWEva ABC Ey_x=Nx + PML elsewhere with Lacunaes, the main source is cut by Theta, WITH Poisson
         call this%mainproblem%problem_PrintCurrentErrors(t);
         do k=0,auxmaxcount-1
            if (this%timers%isused(k)==1) then
               write(*,fmt=50) k, this%timers%isused(k), this%auxproblems(k)%problem_id, this%timers%timer(k)
            !else
            !   write(*,fmt=51) k, this%timers%isused(k), this%timers%timer(k)
            endif
         enddo
      case(24) ! One aux problem, HWeva ABC on Ey_X=Nx, Unsplit PML elsewhere, lacunaes, Poisson
         call this%mainproblem%problem_PrintCurrentErrors(t);
      case(25) ! Sommerfeld ABC, lacunaes with Poisson
         call this%mainproblem%problem_PrintCurrentErrors(t);
         do k=0,auxmaxcount-1
            if (this%timers%isused(k)==1) then
               write(*,fmt=50) k, this%timers%isused(k), this%auxproblems(k)%problem_id, this%timers%timer(k)
            !else
            !   write(*,fmt=51) k, this%timers%isused(k), this%timers%timer(k)
            endif
         enddo   
      end select
      50 format('Slot: ',I2,';    Is Used: ',I1,';    Aux#: ',I2,';    Timer: ',I4);
      51 format('Slot: ',I2,';    Is Used: ',I1,';    Aux#:  -;    Timer: ',I4);
    end subroutine solution_PrintStep;

!WriteOutput----------------------------------------------
  subroutine solution_WriteOutput(this,t)
    !Print the beginning of simulation
    class(tSolution) :: this;
    integer :: t;
    logical :: isopened;
    if ((mod(t,file_erroroutputsteps)==0)) then
       if (file_writelocations==1) then
          write(unit=100, fmt=10) this%mainproblem%Ex_err%val, this%mainproblem%Ex_err%x, this%mainproblem%Ex_err%y, this%mainproblem%Ex_err%z;   
          write(unit=101, fmt=10) this%mainproblem%Ey_err%val, this%mainproblem%Ey_err%x, this%mainproblem%Ey_err%y, this%mainproblem%Ey_err%z;   
          write(unit=102, fmt=10) this%mainproblem%Ez_err%val, this%mainproblem%Ez_err%x, this%mainproblem%Ez_err%y, this%mainproblem%Ez_err%z;   
          write(unit=103, fmt=10) this%mainproblem%Hx_err%val, this%mainproblem%Hx_err%x, this%mainproblem%Hx_err%y, this%mainproblem%Hx_err%z;   
          write(unit=104, fmt=10) this%mainproblem%Hy_err%val, this%mainproblem%Hy_err%x, this%mainproblem%Hy_err%y, this%mainproblem%Hy_err%z;   
          write(unit=105, fmt=10) this%mainproblem%Hz_err%val, this%mainproblem%Hz_err%x, this%mainproblem%Hz_err%y, this%mainproblem%Hz_err%z;
       else
          write(unit=100, fmt=20) this%mainproblem%Ex_err%val;   
          write(unit=101, fmt=20) this%mainproblem%Ey_err%val;   
          write(unit=102, fmt=20) this%mainproblem%Ez_err%val;   
          write(unit=103, fmt=20) this%mainproblem%Hx_err%val;   
          write(unit=104, fmt=20) this%mainproblem%Hy_err%val;   
          write(unit=105, fmt=20) this%mainproblem%Hz_err%val;
       endif

       if (file_pmlerrors==1) then
          if (file_writelocations==1) then
             write(unit=200, fmt=10) this%mainproblem%PML%Ex_err%val, this%mainproblem%PML%Ex_err%x, this%mainproblem%PML%Ex_err%y, this%mainproblem%PML%Ex_err%z;
             write(unit=201, fmt=10) this%mainproblem%PML%Ey_err%val, this%mainproblem%PML%Ey_err%x, this%mainproblem%PML%Ey_err%y, this%mainproblem%PML%Ey_err%z;   
             write(unit=202, fmt=10) this%mainproblem%PML%Ez_err%val, this%mainproblem%PML%Ez_err%x, this%mainproblem%PML%Ez_err%y, this%mainproblem%PML%Ez_err%z;
             write(unit=203, fmt=10) this%mainproblem%PML%Hx_err%val, this%mainproblem%PML%Hx_err%x, this%mainproblem%PML%Hx_err%y, this%mainproblem%PML%Hx_err%z;
             write(unit=204, fmt=10) this%mainproblem%PML%Hy_err%val, this%mainproblem%PML%Hy_err%x, this%mainproblem%PML%Hy_err%y, this%mainproblem%PML%Hy_err%z;
             write(unit=205, fmt=10) this%mainproblem%PML%Hz_err%val, this%mainproblem%PML%Hz_err%x, this%mainproblem%PML%Hz_err%y, this%mainproblem%PML%Hz_err%z;
          else
             write(unit=200, fmt=20) this%mainproblem%PML%Ex_err%val;
             write(unit=201, fmt=20) this%mainproblem%PML%Ey_err%val;   
             write(unit=202, fmt=20) this%mainproblem%PML%Ez_err%val;
             write(unit=203, fmt=20) this%mainproblem%PML%Hx_err%val;
             write(unit=204, fmt=20) this%mainproblem%PML%Hy_err%val;
             write(unit=205, fmt=20) this%mainproblem%PML%Hz_err%val; 
          endif
       endif

       if (file_mainfields==1) then
          if (file_writelocations==1) then
             write(unit=106, fmt=10) maxval(abs(this%mainproblem%Ef%X)), maxloc(abs(this%mainproblem%Ef%X));
             write(unit=107, fmt=10) maxval(abs(this%mainproblem%Ef%Y)), maxloc(abs(this%mainproblem%Ef%Y));
             write(unit=108, fmt=10) maxval(abs(this%mainproblem%Ef%Z)), maxloc(abs(this%mainproblem%Ef%Z));
             write(unit=109, fmt=10) maxval(abs(this%mainproblem%Hf%X)), maxloc(abs(this%mainproblem%Hf%X));
             write(unit=110, fmt=10) maxval(abs(this%mainproblem%Hf%Y)), maxloc(abs(this%mainproblem%Hf%Y));
             write(unit=111, fmt=10) maxval(abs(this%mainproblem%Hf%Z)), maxloc(abs(this%mainproblem%Hf%Z));
             write(unit=112, fmt=10) maxval(abs(this%mainproblem%Je%X)), maxloc(abs(this%mainproblem%Je%X));
             write(unit=113, fmt=10) maxval(abs(this%mainproblem%Je%Y)), maxloc(abs(this%mainproblem%Je%Y));
             write(unit=114, fmt=10) maxval(abs(this%mainproblem%Je%Z)), maxloc(abs(this%mainproblem%Je%Z));
             write(unit=115, fmt=10) maxval(abs(this%mainproblem%Jh%X)), maxloc(abs(this%mainproblem%Jh%X));
             write(unit=116, fmt=10) maxval(abs(this%mainproblem%Jh%Y)), maxloc(abs(this%mainproblem%Jh%Y));
             write(unit=117, fmt=10) maxval(abs(this%mainproblem%Jh%Z)), maxloc(abs(this%mainproblem%Jh%Z));
          else
             write(unit=106, fmt=10) maxval(abs(this%mainproblem%Ef%X));
             write(unit=107, fmt=10) maxval(abs(this%mainproblem%Ef%Y));
             write(unit=108, fmt=10) maxval(abs(this%mainproblem%Ef%Z));
             write(unit=109, fmt=10) maxval(abs(this%mainproblem%Hf%X));
             write(unit=110, fmt=10) maxval(abs(this%mainproblem%Hf%Y));
             write(unit=111, fmt=10) maxval(abs(this%mainproblem%Hf%Z));
             write(unit=112, fmt=10) maxval(abs(this%mainproblem%Je%X));
             write(unit=113, fmt=10) maxval(abs(this%mainproblem%Je%Y));
             write(unit=114, fmt=10) maxval(abs(this%mainproblem%Je%Z));
             write(unit=115, fmt=10) maxval(abs(this%mainproblem%Jh%X));
             write(unit=116, fmt=10) maxval(abs(this%mainproblem%Jh%Y));
             write(unit=117, fmt=10) maxval(abs(this%mainproblem%Jh%Z));
          endif
       endif

       if (file_mainpmlfields==1) then
          if (file_writelocations==1) then
             write(unit=206, fmt=10) maxval(abs(this%mainproblem%PML%Ef%X)), maxloc(abs(this%mainproblem%PML%Ef%X));
             write(unit=207, fmt=10) maxval(abs(this%mainproblem%PML%Ef%Y)), maxloc(abs(this%mainproblem%PML%Ef%Y));
             write(unit=208, fmt=10) maxval(abs(this%mainproblem%PML%Ef%Z)), maxloc(abs(this%mainproblem%PML%Ef%Z));
             write(unit=209, fmt=10) maxval(abs(this%mainproblem%PML%Hf%X)), maxloc(abs(this%mainproblem%PML%Hf%X));
             write(unit=210, fmt=10) maxval(abs(this%mainproblem%PML%Hf%Y)), maxloc(abs(this%mainproblem%PML%Hf%Y));
             write(unit=211, fmt=10) maxval(abs(this%mainproblem%PML%Hf%Z)), maxloc(abs(this%mainproblem%PML%Hf%Z));
             write(unit=212, fmt=10) maxval(abs(this%mainproblem%PML%JfE%X)), maxloc(abs(this%mainproblem%PML%JfE%X));
             write(unit=213, fmt=10) maxval(abs(this%mainproblem%PML%JfE%Y)), maxloc(abs(this%mainproblem%PML%JfE%Y));
             write(unit=214, fmt=10) maxval(abs(this%mainproblem%PML%JfE%Z)), maxloc(abs(this%mainproblem%PML%JfE%Z));
             write(unit=215, fmt=10) maxval(abs(this%mainproblem%PML%JfH%X)), maxloc(abs(this%mainproblem%PML%JfH%X));
             write(unit=216, fmt=10) maxval(abs(this%mainproblem%PML%JfH%Y)), maxloc(abs(this%mainproblem%PML%JfH%Y));
             write(unit=217, fmt=10) maxval(abs(this%mainproblem%PML%JfH%Z)), maxloc(abs(this%mainproblem%PML%JfH%Z));
          else
             write(unit=206, fmt=10) maxval(abs(this%mainproblem%PML%Ef%X));
             write(unit=207, fmt=10) maxval(abs(this%mainproblem%PML%Ef%Y));
             write(unit=208, fmt=10) maxval(abs(this%mainproblem%PML%Ef%Z));
             write(unit=209, fmt=10) maxval(abs(this%mainproblem%PML%Hf%X));
             write(unit=210, fmt=10) maxval(abs(this%mainproblem%PML%Hf%Y));
             write(unit=211, fmt=10) maxval(abs(this%mainproblem%PML%Hf%Z));
             write(unit=212, fmt=10) maxval(abs(this%mainproblem%PML%JfE%X));
             write(unit=213, fmt=10) maxval(abs(this%mainproblem%PML%JfE%Y));
             write(unit=214, fmt=10) maxval(abs(this%mainproblem%PML%JfE%Z));
             write(unit=215, fmt=10) maxval(abs(this%mainproblem%PML%JfH%X));
             write(unit=216, fmt=10) maxval(abs(this%mainproblem%PML%JfH%Y));
             write(unit=217, fmt=10) maxval(abs(this%mainproblem%PML%JfH%Z));
          endif
       endif
       
    endif
    if (file_staticfieldmax==1) then
       if (file_writelocations==1) then
          write(unit=400, fmt=10) maxval(abs(this%buffer%Estatic%X)), maxloc(abs(this%buffer%Estatic%X));
          write(unit=401, fmt=10) maxval(abs(this%buffer%Estatic%Y)), maxloc(abs(this%buffer%Estatic%Y));
          write(unit=402, fmt=10) maxval(abs(this%buffer%Estatic%Z)), maxloc(abs(this%buffer%Estatic%Z));
          write(unit=403, fmt=10) maxval(abs(this%buffer%Hstatic%X)), maxloc(abs(this%buffer%Hstatic%X));
          write(unit=404, fmt=10) maxval(abs(this%buffer%Hstatic%Y)), maxloc(abs(this%buffer%Hstatic%Y));
          write(unit=405, fmt=10) maxval(abs(this%buffer%Hstatic%Z)), maxloc(abs(this%buffer%Hstatic%Z));
       else
          write(unit=400, fmt=20) maxval(abs(this%buffer%Estatic%X));
          write(unit=401, fmt=20) maxval(abs(this%buffer%Estatic%Y));
          write(unit=402, fmt=20) maxval(abs(this%buffer%Estatic%Z));
          write(unit=403, fmt=20) maxval(abs(this%buffer%Hstatic%X));
          write(unit=404, fmt=20) maxval(abs(this%buffer%Hstatic%Y));
          write(unit=405, fmt=20) maxval(abs(this%buffer%Hstatic%Z));          
       endif
    endif
    10 format(E22.15,',',I3,',',I3,',',I3);
    20 format(E22.15);
  end subroutine solution_WriteOutput


!AddAuxProblem----------------------------------------------------------    
    subroutine AddAuxProblem(this)
    ! Add new auxiliary problem
      class(tSolution) :: this;
      type(tProblemStarter) :: pst;
      integer :: pos;
      pst%ptype = 1;
      pst%new_id = this%auxnum;
      pst%Nx = this%Nx + 2*this%dNxaux;        ! number of points for x
      pst%Ny = this%Ny + 2*this%dNyaux;        ! number of points for y
      pst%Nz = this%Nz + 2*this%dNzaux;        ! number of points for z
      pst%Nt = this%Nt;        ! number of points for t
      pst%hhx = this%hhx;
      pst%hhy = this%hhy;
      pst%hhz = this%hhz;
      pst%ht = this%ht;
      pst%schemetype = this%aux_schemetype;
      pst%sourcetype = this%aux_sourcetype;
      pst%PML_thickness = this%aux_PML_thickness;
      pst%PML_param = this%aux_PML_param;
      pst%bufoffsetx=0;
      pst%bufoffsety=0;
      pst%bufoffsetz=0;
      pst%maindx=this%dNxaux;
      pst%maindy=this%dNyaux;
      pst%maindz=this%dNzaux;
      pst%Namu = this%Namu;
      pst%Ndmu = this%Ndmu;
      pst%auxmaxtime = this%auxmaxtime;
      pst%hig1 = this%hig1;
      pst%hig2 = this%hig2;
      pst%gnNum = this%gnNum;
      pst%hwNum = this%hwNum;
      pst%hwEvaE = this%hwEvaE;
      pst%hwEvaP = this%hwEvaP;
      pst%hwEvaEa = this%hwEvaEa;
      pst%hwEvaPa = this%hwEvaPa;
      pst%dtmode = this%dtmode;
      pst%dxmode = this%dxmode;
      pst%dxauxmode = this%dxauxmode;
      pst%dtauxmode = this%dtauxmode;
      pst%sgm_s = this%sgm_s;
      pst%currentstype = this%aux_currentstype;
      pst%pmlcurrentstype = this%aux_pmlcurrentstype;
      pst%effcurrentstype = this%aux_effcurrentstype;
      pst%Eboundarytype = this%auxEboundarytype;
      pst%Hboundarytype = this%auxHboundarytype;
      pos = this%timers%getfreetimer();
      !pos = this%auxnum;
      if (pos==-1) then
         fatalerr = 403;
         write(*,*) 'FATAL ERROR: There is no free timers more, enlarge auxmaxcount!', this%timers%howmanyisused()
      endif
      call this%auxproblems(pos)%problem_init(pst, this%buffer);
      this%timers%isused(pos) = 1;
      if (pst%new_id==0) then
         this%timers%timer(pos) = this%auxlifetime-nint(this%auxmaxtime/this%ht);
      else
         this%timers%timer(pos) = this%auxlifetime+1;
      endif
      this%auxnum = this%auxnum+1;
    end subroutine AddAuxProblem;


!DropAuxProblem----------------------------------------------------------    
    subroutine DropAuxProblem(this, Num)
    ! Add new auxiliary problem
      class(tSolution) :: this;
      integer, intent(in) :: Num;
      call this%auxproblems(Num)%problem_destroy;
      call this%timers%closetimer(Num);
    end subroutine DropAuxProblem;
    

!ManageAuxProblems----------------------------------------------------------    
    subroutine ManageAuxProblems(this, t)
    ! Add and destroy auxiliary problems if necessary
      class(tSolution) :: this;
      integer, intent(in) :: t;
      integer :: k;
      ! Add new problems
      if (t/=0.AND.mod(t,nint((this%sgm_s*this%auxmaxtime+(1+this%sgm_s)*this%auxmaxtime*(this%auxnum-1))/this%ht))==0) then
         call this%AddAuxProblem;
         write(*,*) 'Creating new aux problem', this%auxnum-1;
         write(*,*) 'Now we have aux problems', this%timers%howmanyisused();
      endif
      ! Refresh timers and drop unnecessary aux tasks
      do k=0,auxmaxcount-1
         if (this%timers%isused(k)==1) then
            this%timers%timer(k) = this%timers%timer(k)-1;
            if (this%timers%timer(k)==0) then
               if (this%soltype==5) then
                  write(*,*) 'Add E to Estatic from the problem', k;
                  call this%auxproblems(k)%addEstaticToBuffer;
               endif
               write(*,*) 'Dropping aux problem', k;
               call this%DropAuxProblem(k);
            endif
         endif
      enddo;
   end subroutine ManageAuxProblems;     
  
end module solutionclass
