module problemclass
  use commonvars
  use meshclass
  use auxmeshclass
  use edgeauxmeshclass
  use sourceclass
  use pmlclass
  use parallel
  use bufferclass
  use poissonclass
  implicit none;
  
  type, public :: tProblem
    integer problem_type;
      ! 0 - main problem
      ! 1 - auxiliary problem
    integer :: problem_id;
    integer :: schemetype;
      ! 0 - Yee 2-2 scheme
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
    integer :: Eboundarytype;  
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
      ! 14 - Hagstrom-Warburton abc Ey_x=0, Unsplit PML elsewhere
      ! 15 - HWEva ABC Ey_x=Nx, exact elsewhere
      ! 16 - HWEva ABC Ey_x=Nx, Unpslit PML elsewhere
    integer :: Hboundarytype;
    integer :: Nx, Ny, Nz, Nt;      ! discretization values
    real*8  :: hhx, hhy, hhz, ht;   ! differentials
    type(tMesh) Ef, Hf;             ! actual fields of the interior problem
    type(tMesh) Ean, Han, Je, Jh;   ! fields and currents from analytical solution
    type(tMesh) Efold, Hfold;       ! previous values of the actual fields
    type(tMesh) Efold2, Hfold2;     ! previous previous values of the actual fields
    type(tMesh) JEeff, JHeff;       ! Effective currents of the interior problem
    type(tSource) Source;           ! Source object that provides analytical solution and currents
    type(tPML) PML;                 ! PML object
    ! Poisson parameters
    integer :: npx, npy, npz; 
    type(tPoisson) Poisson;         ! Poisson object
    ! Physical points and arrays
    real*8, allocatable :: xi(:), yi(:), zi(:), ti(:);
    real*8, allocatable :: xi05(:), yi05(:), zi05(:), ti05(:);
    real*8, allocatable :: mup(:,:,:);
    ! Errors variables
    type(tFieldError) :: Ex_err, Ey_err, Ez_err, Hx_err, Hy_err, Hz_err;
    type(tFieldError) :: Ex_maxerr, Ey_maxerr, Ez_maxerr, Hx_maxerr, Hy_maxerr, Hz_maxerr;
    character(len=2) :: err_f;
    real*8 :: err, abserr;
    ! ABC Parameters
    integer :: dtmode, dxmode, dxauxmode, dtauxmode;
    real*8 :: hig1, hig2;             ! Angles of Higdon ABC
    ! Givoli-Neta Aux Variables
    integer :: gnNum;                                                                                                             ! Number of Givoli-Neta angles
    type(tAuxMesh), allocatable :: GNx(:), GNy(:), GNz(:), GNxold(:), GNyold(:), GNzold(:), GNxold2(:), GNyold2(:), GNzold2(:);   ! Givoli-Neta Aux Variables
    real*8, allocatable :: gnA(:), GNs(:), GNalp(:), GNbet(:);                                                                    ! Givoli-Neta coefficients
    ! Hagstrom-Warburton Aux Variables
    integer :: hwNum;                                                                                                             ! Number of HW angles
    type(tAuxMesh), allocatable :: hwX(:), hwY(:), hwZ(:), hwXold(:), hwYold(:), hwZold(:), hwXold2(:), hwYold2(:), hwZold2(:);   ! HW Aux Variables
    real*8, allocatable :: hwAng(:), hwAc(:), hwL(:,:), hwA(:,:), hwB(:,:), hwC(:,:), hwD(:,:);                                   ! HW coefficients and matrices
    real*8, allocatable :: hwFnew(:), hwF(:), hwFold(:), hwGF(:);                                                                 ! HW additional arrays
    ! HWeva additional Aux Variables
    integer :: hwEvaE;
    integer :: hwEvaP;                                                                                                            ! Number of HW evanescent variables
    real*8, allocatable :: hwSigmas(:), hwM(:,:), hwN(:,:), hwS(:,:), hwV(:,:)                                                    ! HW coefficients and matrices
    ! HWEva Edge Aux Veriables
    integer :: hwEvaEa;
    integer :: hwEvaPa;                                                                                                           ! Number of HW evanescent variables for edges
    type(tEdgeAuxMesh), allocatable :: hwYEdgeNx(:,:), hwYEdgeNxOld(:,:), hwYEdgeNxOld2(:,:);                                     ! HW Edge Aux Variables
    real*8, allocatable :: hwAngA(:), hwAcA(:), hwLa(:,:), hwAa(:,:), hwBa(:,:), hwCa(:,:), hwDa(:,:);                            ! HW Edge coefficients and matrices
    real*8, allocatable :: hwFAnew(:), hwFA(:), hwFAold(:), hwGFA(:);                                                             ! HW additional arrays
    real*8, allocatable :: hwSigmasA(:), hwMa(:,:), hwNa(:,:), hwSa(:,:), hwVa(:,:);
    ! Buffer
    class(TBuffer), pointer :: buffer;
    ! Aux paramenters
    integer :: bufoffsetx, bufoffsety, bufoffsetz;
    integer :: maindx, maindy, maindz; 
    real*8 :: auxmaxtime;
    real*8 :: sgm_s;
    ! Effective currents params
    integer :: Namu;
ц    integer :: Ndmu;
    integer :: Nlc;
    real*8 :: Dlc;
    contains
      procedure problem_Init
      procedure problem_Destroy
      procedure initPhysPoints
      procedure initMikeSource;
      procedure updateEInterior
      procedure updateHInterior
      procedure getSourceCurrents
      procedure getAnalyticSolution
      procedure problem_DoIndependentStep
      procedure problem_DoIndependentStepTiming
      procedure fillEBoundary
      procedure fillEBoundaryABCSommerfeld
      procedure fillEBoundaryABCHigdon
      procedure fillEBoundaryABCBetzMittra
      procedure fillEBoundaryABCMur
      procedure fillEBoundaryABCGivoliNeta
      procedure fillEBoundaryABCGivoliNetaPML
      procedure fillEBoundaryABCGivoliNetaPure
      procedure fillEBoundaryABCHW
      procedure fillEBoundaryABCHWPML
      procedure fillEBoundaryABCHWEva
      procedure fillEBoundaryABCHWEvaPML
      procedure fillHBoundary
      procedure savefields
      ! point functions
      procedure updateEpoint
      procedure updateHpoint
      procedure fillEBoundaryPoint
      procedure fillHBoundaryPoint
      procedure getCurrentEPoint
      procedure getCurrentHPoint
      procedure checkEPointUpdated
      procedure JEeffbuildPoint
      procedure JHeffbuildPoint
      procedure JEeffbuildPointPoisson
      ! error functions
      procedure getCurrentError
      procedure getCurrentErrorPML
      procedure GetDifference
      procedure SaveToMaxError
      ! pml functions
      procedure getSourceCurrentsPML
      procedure getAnalyticSolutionPML
      procedure uploadFieldsToPML
      procedure getCurrentEPointPML
      procedure getCurrentHPointPML
      ! buffer functions
      procedure addEToBuffer
      procedure addHToBuffer
      procedure addEstaticToBuffer
      procedure addEcompositesToBuffer
      procedure saveEffCurrentsToBuffer
      procedure saveOriginalCurrentsToBuffer
      ! effective currents formation
      procedure JEeffBuild
      procedure JHeffBuild
      ! Mu and Hat functions
      procedure P4
      procedure P4inv
      procedure thetabig
      procedure Muu
      procedure InitMuArray
      ! Poisson procedures
      procedure UploadFieldsToPoisson;
      ! output procedures
      procedure problem_PrintCurrentErrors
      procedure problem_PrintMaximumAnalyticalFields
      procedure problem_PrintMaximumFields
      procedure problem_PrintPMLCurrentErrors
      procedure problem_PrintMaximumCurrents
      procedure UploadMainMuToBuffer;
      procedure problem_report;
      procedure WriteFieldToFile;
      ! Givoli-Neta ABC functions
      procedure initGivoliNetaAuxVariables;
      procedure saveGivoliNetaAuxVariables;
      procedure destroyGivoliNetaAuxVariables;
      procedure GivoliNetaEYmain;
      procedure GivoliNetaEYaux;
      ! Hagstrom-Warburton ABC functions
      procedure initHWAuxVariables;
      procedure saveHWAuxVariables;
      procedure destroyHWAuxVariables;
      procedure HWEyUpdate;
      ! HWEva ABC functions
      procedure initHWEvaAuxVariables;
      procedure saveHWEvaAuxVariables;
      procedure destroyHWEvaAuxVariables;
      procedure initHWEvaEdgeAuxVariables;
      procedure saveHWEvaEdgeAuxVariables;
      procedure destroyHWEvaEdgeAuxVariables;
      procedure HWEvaEyUpdate;
      procedure HWEvaEySomUpdate;
      procedure HWEvaEyHigdonUpdate;
      procedure HWEvaEySuperUpdate;
      procedure HWEvaEyBetzMittraUpdate;
      procedure HWEvaEyHigdon3DUpdate;
      ! Timer Functions
      procedure starttimer;
      procedure endtimer;
  end type tProblem

contains

!-- Sommerfeld ABC, Higdon ABC, Betz-Mittra ABC, Mur ABC are in the following file:
  include 'simpleABCs.f90'

!-- Givoli-Neta ABC is in the following file:
  include 'GivoliNetaABC.f90'

!-- Hagstrom-Warburton Basic ABC is in the following file:
  include 'HagstromWarburtonBasicABC.f90'  

  
!problem_init------------------------------------------------------------------   
    subroutine problem_Init(this, pstarter, pbuffer)
    !Initialize the problem
      class(tProblem) :: this
      class(tProblemStarter), intent(in) :: pstarter
      class(TBuffer), pointer :: pbuffer;
      integer :: i,j,k;
      if (pstarter%ptype < 2) then
         this%problem_type = pstarter%ptype;
      else
         fatalerr = 215;
         write(*,*) 'FATAL ERROR, Incorrect Problem Type is used while problem initiation!', pstarter%ptype;
      endif
      this%problem_id = pstarter%new_id;
      this%Eboundarytype = pstarter%Eboundarytype;
      this%Hboundarytype = pstarter%Hboundarytype;
      this%currentstype = pstarter%currentstype;
      this%pmlcurrentstype = pstarter%pmlcurrentstype;
      this%effcurrentstype = pstarter%effcurrentstype;
      this%schemetype = pstarter%schemetype;
      this%Source%sourcetype = pstarter%sourcetype;
      if ((pstarter%Nx>0).AND.(pstarter%Nx>0).AND.(pstarter%Nx>0).AND.(pstarter%Nt>0)) then
         this%Nx = pstarter%Nx;
         this%Ny = pstarter%Ny;
         this%Nz = pstarter%Nz;
         this%Nt = pstarter%Nt;        
      else
         fatalerr = 100;
         write(*,*) 'FATAL ERROR, Problem is initiated with 0 points in one dimension';
         write(*,*) 'Nx', pstarter%Nx, 'Ny', pstarter%Ny, 'Nz', pstarter%Nz, 'Nt', pstarter%Nt;
      endif
      if ((pstarter%hhx>0).AND.(pstarter%hhy>0).AND.(pstarter%hhz>0).AND.(pstarter%ht>0)) then
         this%hhx = pstarter%hhx;
         this%hhy = pstarter%hhy;
         this%hhz = pstarter%hhz;
         this%ht = pstarter%ht; 
      else
         fatalerr = 101;
         write(*,*) 'FATAL ERROR, Problem is initiated with zero differentials';
         write(*,*) 'Nx', pstarter%hhx, 'Ny', pstarter%hhy, 'Nz', pstarter%hhz, 'Nt', pstarter%ht;
      endif

      this%npx = pstarter%npx;
      this%npy = pstarter%npy;
      this%npz = pstarter%npz;  

      this%bufoffsetx = pstarter%bufoffsetx;
      this%bufoffsety = pstarter%bufoffsety;
      this%bufoffsetz = pstarter%bufoffsetz;

      this%maindx = pstarter%maindx;
      this%maindy = pstarter%maindy;
      this%maindz = pstarter%maindz;

      this%hig1 = pstarter%hig1;
      this%hig2 = pstarter%hig2;
      this%gnNum = pstarter%gnNum;
      this%hwNum = pstarter%hwNum;
      this%hwEvaE = pstarter%hwEvaE;
      this%hwEvaP = pstarter%hwEvaP;
      this%hwEvaEa = pstarter%hwEvaEa;
      this%hwEvaPa = pstarter%hwEvaPa;

      this%dtmode = pstarter%dtmode;
      this%dxmode = pstarter%dxmode;
      this%dxauxmode = pstarter%dxauxmode;
      this%dtauxmode = pstarter%dtauxmode;
      
      this%buffer => pbuffer;
     
      this%auxmaxtime = pstarter%auxmaxtime;
      this%sgm_s=pstarter%sgm_s;

      this%Ex_maxerr%absval = 0.0d0;
      this%Ey_maxerr%absval = 0.0d0;
      this%Ez_maxerr%absval = 0.0d0;
      this%Hx_maxerr%absval = 0.0d0;
      this%Hy_maxerr%absval = 0.0d0;
      this%Hz_maxerr%absval = 0.0d0;
      this%Ex_maxerr%fval = 0.0d0;
      this%Ey_maxerr%fval = 0.0d0;
      this%Ez_maxerr%fval = 0.0d0;
      this%Hx_maxerr%fval = 0.0d0;
      this%Hy_maxerr%fval = 0.0d0;
      this%Hz_maxerr%fval = 0.0d0;
      this%Ex_maxerr%val = 0.0d0;
      this%Ey_maxerr%val = 0.0d0;
      this%Ez_maxerr%val = 0.0d0;
      this%Hx_maxerr%val = 0.0d0;
      this%Hy_maxerr%val = 0.0d0;
      this%Hz_maxerr%val = 0.0d0;
     
      call this%Ef%mesh_init(this%Nx,this%Ny,this%Nz,0,this%problem_id);
      call this%Hf%mesh_init(this%Nx,this%Ny,this%Nz,1,this%problem_id);
      call this%Ean%mesh_init(this%Nx,this%Ny,this%Nz,0,this%problem_id);
      call this%Han%mesh_init(this%Nx,this%Ny,this%Nz,1,this%problem_id);
      call this%Je%mesh_init(this%Nx,this%Ny,this%Nz,0,this%problem_id);
      call this%Jh%mesh_init(this%Nx,this%Ny,this%Nz,1,this%problem_id);      
      call this%Efold%mesh_init(this%Nx,this%Ny,this%Nz,0,this%problem_id);
      call this%Hfold%mesh_init(this%Nx,this%Ny,this%Nz,1,this%problem_id);
      call this%Efold2%mesh_init(this%Nx,this%Ny,this%Nz,0,this%problem_id);
      call this%Hfold2%mesh_init(this%Nx,this%Ny,this%Nz,1,this%problem_id);
      call this%JEeff%mesh_init(this%Nx,this%Ny,this%Nz,0,this%problem_id);
      call this%JHeff%mesh_init(this%Nx,this%Ny,this%Nz,1,this%problem_id);
      
      allocate(this%xi(0:this%Nx));
      allocate(this%yi(0:this%Ny));
      allocate(this%zi(0:this%Nz));
      allocate(this%ti(0:this%Nt));
      allocate(this%xi05(0:this%Nx+1));
      allocate(this%yi05(0:this%Ny+1));
      allocate(this%zi05(0:this%Nz+1));
      allocate(this%ti05(0:(this%Nt+1)));
      allocate(this%mup(0:2*(this%Nx+1),0:2*(this%Ny+1),0:2*(this%Nz+1)));
      
      call this%initPhysPoints;
      
      if (this%Source%sourcetype == 0) then
         if (this%problem_type==0) then
            if (source_mainalloc == 0) then
               call this%initMikeSource(0);
               source_mainalloc = 1;
            endif
         else if (this%problem_type==1) then
            if (source_auxalloc == 0) then
               call this%initMikeSource(1);
               source_auxalloc = 1;
            endif
         endif   
      endif
      this%Source%problem_type = this%problem_type;
      
      ! Poisson initialization
      if (this%problem_type==0.AND.this%effcurrentstype==1) then
         call this%Poisson%poisson_init(this%Nx, this%Ny, this%Nz, this%hhx, this%hhy, this%hhz, this%ht, this%problem_id, this%npx, this%npy, this%npz);
      endif
      
      this%Namu = pstarter%Namu;
      this%Ndmu = pstarter%Ndmu;
      this%Nlc = this%Namu + this%Ndmu;
      this%Dlc=this%xi(this%Nx-this%Namu)-this%xi(this%Nx-this%Nlc);
      
      ! Mu array initialization
      if (this%problem_type==0) then
         call this%InitMuArray;
      endif
      
      !PML Initialization
      if ((this%Eboundarytype==3).OR.(this%Eboundarytype==10).OR.(this%Eboundarytype==14).OR.(this%Eboundarytype==16)) then
         call this%PML%pml_init(this%Nx, this%Ny, this%Nz, this%hhx, this%hhy, this%hhz, this%ht, pstarter%PML_thickness, pstarter%PML_param, this%problem_id*100+1);
      endif
      
      ! Givoli-Neta coefficients initialization
      if ((this%Eboundarytype==9).OR.(this%Eboundarytype==10).OR.(this%Eboundarytype==11).OR.(this%Eboundarytype==12)) then
         call this%initGivoliNetaAuxVariables;
      endif
      
      ! Hagstrom-Warburton coefficients initialization
      if (this%Eboundarytype==13.OR.this%Eboundarytype==14) then
         call this%initHWAuxVariables;
      endif

      ! HWeva coefficients initialization
      if (this%Eboundarytype==15.OR.this%Eboundarytype==16) then
         call this%initHWEvaAuxVariables;
      endif
      
      call this%problem_report;
    end subroutine problem_Init;

!initMikeSource------------------------------------------------------------------   
    subroutine initMikeSource(this, mode)
    !Initialize the problem
      class(tProblem) :: this
      integer, intent(in) :: mode;
      integer :: i, j, k;
      if (mode==0) then
         allocate(MEJx(1:this%Nx, 0:this%Ny, 0:this%Nz));
         allocate(MEJy(0:this%Nx, 1:this%Ny, 0:this%Nz));
         allocate(MEJz(0:this%Nx, 0:this%Ny, 1:this%Nz));
         allocate(MHx(0:this%Nx, 1:this%Ny, 1:this%Nz));
         allocate(MHy(1:this%Nx, 0:this%Ny, 1:this%Nz));
         allocate(MHz(1:this%Nx, 1:this%Ny, 0:this%Nz));
         allocate(DMEJx(0:14,1:this%Nx, 0:this%Ny, 0:this%Nz));
         allocate(DMEJy(0:14,0:this%Nx, 1:this%Ny, 0:this%Nz));
         allocate(DMEJz(0:14,0:this%Nx, 0:this%Ny, 1:this%Nz));
         allocate(DMHx(0:14,0:this%Nx, 1:this%Ny, 1:this%Nz));
         allocate(DMHy(0:14,1:this%Nx, 0:this%Ny, 1:this%Nz));
         allocate(DMHz(0:14,1:this%Nx, 1:this%Ny, 0:this%Nz));
         do k=0,this%Nz
            do j=0,this%Ny
               do i=1,this%Nx
                  MEJx(i,j,k)=getM(this%xi05(i), this%yi(j), this%zi(k));
                  DMEJx(0,i,j,k)=getDM(this%xi05(i), this%yi(j), this%zi(k), 1);
                  DMEJx(1,i,j,k)=getDM(this%xi05(i), this%yi(j), this%zi(k), 2);
                  DMEJx(2,i,j,k)=getDM(this%xi05(i), this%yi(j), this%zi(k), 3);
                  DMEJx(3,i,j,k)=getDM(this%xi05(i), this%yi(j), this%zi(k), 31);
                  DMEJx(4,i,j,k)=getDM(this%xi05(i), this%yi(j), this%zi(k), 32);
                  DMEJx(5,i,j,k)=getDM(this%xi05(i), this%yi(j), this%zi(k), 21);
                  DMEJx(6,i,j,k)=getDM(this%xi05(i), this%yi(j), this%zi(k), 11);
                  DMEJx(7,i,j,k)=getDM(this%xi05(i), this%yi(j), this%zi(k), 22);
                  DMEJx(8,i,j,k)=getDM(this%xi05(i), this%yi(j), this%zi(k), 33);
                  DMEJx(9,i,j,k)=getDM(this%xi05(i), this%yi(j), this%zi(k), 222);
                  DMEJx(10,i,j,k)=getDM(this%xi05(i), this%yi(j), this%zi(k), 333);
                  DMEJx(11,i,j,k)=getDM(this%xi05(i), this%yi(j), this%zi(k), 322);
                  DMEJx(12,i,j,k)=getDM(this%xi05(i), this%yi(j), this%zi(k), 311);
                  DMEJx(13,i,j,k)=getDM(this%xi05(i), this%yi(j), this%zi(k), 211);
                  DMEJx(14,i,j,k)=getDM(this%xi05(i), this%yi(j), this%zi(k), 332);
               enddo
            enddo
         enddo
         do k=0,this%Nz
            do j=1,this%Ny
               do i=0,this%Nx
                  MEJy(i,j,k)=getM(this%xi(i), this%yi05(j), this%zi(k));
                  DMEJy(0,i,j,k)=getDM(this%xi(i), this%yi05(j), this%zi(k), 1);
                  DMEJy(1,i,j,k)=getDM(this%xi(i), this%yi05(j), this%zi(k), 2);
                  DMEJy(2,i,j,k)=getDM(this%xi(i), this%yi05(j), this%zi(k), 3);
                  DMEJy(3,i,j,k)=getDM(this%xi(i), this%yi05(j), this%zi(k), 31);
                  DMEJy(4,i,j,k)=getDM(this%xi(i), this%yi05(j), this%zi(k), 32);
                  DMEJy(5,i,j,k)=getDM(this%xi(i), this%yi05(j), this%zi(k), 21);
                  DMEJy(6,i,j,k)=getDM(this%xi(i), this%yi05(j), this%zi(k), 11);
                  DMEJy(7,i,j,k)=getDM(this%xi(i), this%yi05(j), this%zi(k), 22);
                  DMEJy(8,i,j,k)=getDM(this%xi(i), this%yi05(j), this%zi(k), 33);
                  DMEJy(9,i,j,k)=getDM(this%xi(i), this%yi05(j), this%zi(k), 222);
                  DMEJy(10,i,j,k)=getDM(this%xi(i), this%yi05(j), this%zi(k), 333);
                  DMEJy(11,i,j,k)=getDM(this%xi(i), this%yi05(j), this%zi(k), 322);
                  DMEJy(12,i,j,k)=getDM(this%xi(i), this%yi05(j), this%zi(k), 311);
                  DMEJy(13,i,j,k)=getDM(this%xi(i), this%yi05(j), this%zi(k), 211);
                  DMEJy(14,i,j,k)=getDM(this%xi(i), this%yi05(j), this%zi(k), 332);
               enddo
            enddo
         enddo
         do k=1,this%Nz
            do j=0,this%Ny
               do i=0,this%Nx
                  MEJz(i,j,k)=getM(this%xi(i), this%yi(j), this%zi05(k));
                  DMEJz(0,i,j,k)=getDM(this%xi(i), this%yi(j), this%zi05(k), 1);
                  DMEJz(1,i,j,k)=getDM(this%xi(i), this%yi(j), this%zi05(k), 2);
                  DMEJz(2,i,j,k)=getDM(this%xi(i), this%yi(j), this%zi05(k), 3);
                  DMEJz(3,i,j,k)=getDM(this%xi(i), this%yi(j), this%zi05(k), 31);
                  DMEJz(4,i,j,k)=getDM(this%xi(i), this%yi(j), this%zi05(k), 32);
                  DMEJz(5,i,j,k)=getDM(this%xi(i), this%yi(j), this%zi05(k), 21);
                  DMEJz(6,i,j,k)=getDM(this%xi(i), this%yi(j), this%zi05(k), 11);
                  DMEJz(7,i,j,k)=getDM(this%xi(i), this%yi(j), this%zi05(k), 22);
                  DMEJz(8,i,j,k)=getDM(this%xi(i), this%yi(j), this%zi05(k), 33);
                  DMEJz(9,i,j,k)=getDM(this%xi(i), this%yi(j), this%zi05(k), 222);
                  DMEJz(10,i,j,k)=getDM(this%xi(i), this%yi(j), this%zi05(k), 333);
                  DMEJz(11,i,j,k)=getDM(this%xi(i), this%yi(j), this%zi05(k), 322);
                  DMEJz(12,i,j,k)=getDM(this%xi(i), this%yi(j), this%zi05(k), 311);
                  DMEJz(13,i,j,k)=getDM(this%xi(i), this%yi(j), this%zi05(k), 211);
                  DMEJz(14,i,j,k)=getDM(this%xi(i), this%yi(j), this%zi05(k), 332);
               enddo
            enddo
         enddo
         
!!! H !!!
         do k=1,this%Nz
            do j=1,this%Ny
               do i=0,this%Nx
                  MHx(i,j,k)=getM(this%xi(i), this%yi05(j), this%zi05(k));
                  DMHx(0,i,j,k)=getDM(this%xi(i), this%yi05(j), this%zi05(k), 1);
                  DMHx(1,i,j,k)=getDM(this%xi(i), this%yi05(j), this%zi05(k), 2);
                  DMHx(2,i,j,k)=getDM(this%xi(i), this%yi05(j), this%zi05(k), 3);
                  DMHx(3,i,j,k)=getDM(this%xi(i), this%yi05(j), this%zi05(k), 31);
                  DMHx(4,i,j,k)=getDM(this%xi(i), this%yi05(j), this%zi05(k), 32);
                  DMHx(5,i,j,k)=getDM(this%xi(i), this%yi05(j), this%zi05(k), 21);
                  DMHx(6,i,j,k)=getDM(this%xi(i), this%yi05(j), this%zi05(k), 11);
                  DMHx(7,i,j,k)=getDM(this%xi(i), this%yi05(j), this%zi05(k), 22);
                  DMHx(8,i,j,k)=getDM(this%xi(i), this%yi05(j), this%zi05(k), 33);
                  DMHx(9,i,j,k)=getDM(this%xi(i), this%yi05(j), this%zi05(k), 222);
                  DMHx(10,i,j,k)=getDM(this%xi(i), this%yi05(j), this%zi05(k), 333);
                  DMHx(11,i,j,k)=getDM(this%xi(i), this%yi05(j), this%zi05(k), 322);
                  DMHx(12,i,j,k)=getDM(this%xi(i), this%yi05(j), this%zi05(k), 311);
                  DMHx(13,i,j,k)=getDM(this%xi(i), this%yi05(j), this%zi05(k), 211);
                  DMHx(14,i,j,k)=getDM(this%xi(i), this%yi05(j), this%zi05(k), 332);
               enddo
            enddo
         enddo
         do k=1,this%Nz
            do j=0,this%Ny
               do i=1,this%Nx
                  MHy(i,j,k)=getM(this%xi05(i), this%yi(j), this%zi05(k));
                  DMHy(0,i,j,k)=getDM(this%xi05(i), this%yi(j), this%zi05(k), 1);
                  DMHy(1,i,j,k)=getDM(this%xi05(i), this%yi(j), this%zi05(k), 2);
                  DMHy(2,i,j,k)=getDM(this%xi05(i), this%yi(j), this%zi05(k), 3);
                  DMHy(3,i,j,k)=getDM(this%xi05(i), this%yi(j), this%zi05(k), 31);
                  DMHy(4,i,j,k)=getDM(this%xi05(i), this%yi(j), this%zi05(k), 32);
                  DMHy(5,i,j,k)=getDM(this%xi05(i), this%yi(j), this%zi05(k), 21);
                  DMHy(6,i,j,k)=getDM(this%xi05(i), this%yi(j), this%zi05(k), 11);
                  DMHy(7,i,j,k)=getDM(this%xi05(i), this%yi(j), this%zi05(k), 22);
                  DMHy(8,i,j,k)=getDM(this%xi05(i), this%yi(j), this%zi05(k), 33);
                  DMHy(9,i,j,k)=getDM(this%xi05(i), this%yi(j), this%zi05(k), 222);
                  DMHy(10,i,j,k)=getDM(this%xi05(i), this%yi(j), this%zi05(k), 333);
                  DMHy(11,i,j,k)=getDM(this%xi05(i), this%yi(j), this%zi05(k), 322);
                  DMHy(12,i,j,k)=getDM(this%xi05(i), this%yi(j), this%zi05(k), 311);
                  DMHy(13,i,j,k)=getDM(this%xi05(i), this%yi(j), this%zi05(k), 211);
                  DMHy(14,i,j,k)=getDM(this%xi05(i), this%yi(j), this%zi05(k), 332);
               enddo
            enddo
         enddo
         do k=0,this%Nz
            do j=1,this%Ny
               do i=1,this%Nx
                  MHz(i,j,k)=getM(this%xi05(i), this%yi05(j), this%zi(k));
                  DMHz(0,i,j,k)=getDM(this%xi05(i), this%yi05(j), this%zi(k), 1);
                  DMHz(1,i,j,k)=getDM(this%xi05(i), this%yi05(j), this%zi(k), 2);
                  DMHz(2,i,j,k)=getDM(this%xi05(i), this%yi05(j), this%zi(k), 3);
                  DMHz(3,i,j,k)=getDM(this%xi05(i), this%yi05(j), this%zi(k), 31);
                  DMHz(4,i,j,k)=getDM(this%xi05(i), this%yi05(j), this%zi(k), 32);
                  DMHz(5,i,j,k)=getDM(this%xi05(i), this%yi05(j), this%zi(k), 21);
                  DMHz(6,i,j,k)=getDM(this%xi05(i), this%yi05(j), this%zi(k), 11);
                  DMHz(7,i,j,k)=getDM(this%xi05(i), this%yi05(j), this%zi(k), 22);
                  DMHz(8,i,j,k)=getDM(this%xi05(i), this%yi05(j), this%zi(k), 33);
                  DMHz(9,i,j,k)=getDM(this%xi05(i), this%yi05(j), this%zi(k), 222);
                  DMHz(10,i,j,k)=getDM(this%xi05(i), this%yi05(j), this%zi(k), 333);
                  DMHz(11,i,j,k)=getDM(this%xi05(i), this%yi05(j), this%zi(k), 322);
                  DMHz(12,i,j,k)=getDM(this%xi05(i), this%yi05(j), this%zi(k), 311);
                  DMHz(13,i,j,k)=getDM(this%xi05(i), this%yi05(j), this%zi(k), 211);
                  DMHz(14,i,j,k)=getDM(this%xi05(i), this%yi05(j), this%zi(k), 332);
               enddo
            enddo
         enddo
         write(*,*) 'Main mike source matrices inited!';
      else if (mode==1) then
         allocate(auxMEJx(1:this%Nx, 0:this%Ny, 0:this%Nz));
         allocate(auxMEJy(0:this%Nx, 1:this%Ny, 0:this%Nz));
         allocate(auxMEJz(0:this%Nx, 0:this%Ny, 1:this%Nz));
         allocate(auxMHx(0:this%Nx, 1:this%Ny, 1:this%Nz));
         allocate(auxMHy(1:this%Nx, 0:this%Ny, 1:this%Nz));
         allocate(auxMHz(1:this%Nx, 1:this%Ny, 0:this%Nz));
         allocate(auxDMEJx(0:14,1:this%Nx, 0:this%Ny, 0:this%Nz));
         allocate(auxDMEJy(0:14,0:this%Nx, 1:this%Ny, 0:this%Nz));
         allocate(auxDMEJz(0:14,0:this%Nx, 0:this%Ny, 1:this%Nz));
         allocate(auxDMHx(0:14,0:this%Nx, 1:this%Ny, 1:this%Nz));
         allocate(auxDMHy(0:14,1:this%Nx, 0:this%Ny, 1:this%Nz));
         allocate(auxDMHz(0:14,1:this%Nx, 1:this%Ny, 0:this%Nz));
         do k=0,this%Nz
            do j=0,this%Ny
               do i=1,this%Nx
                  auxMEJx(i,j,k)=getM(this%xi05(i), this%yi(j), this%zi(k));
                  auxDMEJx(0,i,j,k)=getDM(this%xi05(i), this%yi(j), this%zi(k), 1);
                  auxDMEJx(1,i,j,k)=getDM(this%xi05(i), this%yi(j), this%zi(k), 2);
                  auxDMEJx(2,i,j,k)=getDM(this%xi05(i), this%yi(j), this%zi(k), 3);
                  auxDMEJx(3,i,j,k)=getDM(this%xi05(i), this%yi(j), this%zi(k), 31);
                  auxDMEJx(4,i,j,k)=getDM(this%xi05(i), this%yi(j), this%zi(k), 32);
                  auxDMEJx(5,i,j,k)=getDM(this%xi05(i), this%yi(j), this%zi(k), 21);
                  auxDMEJx(6,i,j,k)=getDM(this%xi05(i), this%yi(j), this%zi(k), 11);
                  auxDMEJx(7,i,j,k)=getDM(this%xi05(i), this%yi(j), this%zi(k), 22);
                  auxDMEJx(8,i,j,k)=getDM(this%xi05(i), this%yi(j), this%zi(k), 33);
                  auxDMEJx(9,i,j,k)=getDM(this%xi05(i), this%yi(j), this%zi(k), 222);
                  auxDMEJx(10,i,j,k)=getDM(this%xi05(i), this%yi(j), this%zi(k), 333);
                  auxDMEJx(11,i,j,k)=getDM(this%xi05(i), this%yi(j), this%zi(k), 322);
                  auxDMEJx(12,i,j,k)=getDM(this%xi05(i), this%yi(j), this%zi(k), 311);
                  auxDMEJx(13,i,j,k)=getDM(this%xi05(i), this%yi(j), this%zi(k), 211);
                  auxDMEJx(14,i,j,k)=getDM(this%xi05(i), this%yi(j), this%zi(k), 332);
               enddo
            enddo
         enddo
         do k=0,this%Nz
            do j=1,this%Ny
               do i=0,this%Nx
                  auxMEJy(i,j,k)=getM(this%xi(i), this%yi05(j), this%zi(k));
                  auxDMEJy(0,i,j,k)=getDM(this%xi(i), this%yi05(j), this%zi(k), 1);
                  auxDMEJy(1,i,j,k)=getDM(this%xi(i), this%yi05(j), this%zi(k), 2);
                  auxDMEJy(2,i,j,k)=getDM(this%xi(i), this%yi05(j), this%zi(k), 3);
                  auxDMEJy(3,i,j,k)=getDM(this%xi(i), this%yi05(j), this%zi(k), 31);
                  auxDMEJy(4,i,j,k)=getDM(this%xi(i), this%yi05(j), this%zi(k), 32);
                  auxDMEJy(5,i,j,k)=getDM(this%xi(i), this%yi05(j), this%zi(k), 21);
                  auxDMEJy(6,i,j,k)=getDM(this%xi(i), this%yi05(j), this%zi(k), 11);
                  auxDMEJy(7,i,j,k)=getDM(this%xi(i), this%yi05(j), this%zi(k), 22);
                  auxDMEJy(8,i,j,k)=getDM(this%xi(i), this%yi05(j), this%zi(k), 33);
                  auxDMEJy(9,i,j,k)=getDM(this%xi(i), this%yi05(j), this%zi(k), 222);
                  auxDMEJy(10,i,j,k)=getDM(this%xi(i), this%yi05(j), this%zi(k), 333);
                  auxDMEJy(11,i,j,k)=getDM(this%xi(i), this%yi05(j), this%zi(k), 322);
                  auxDMEJy(12,i,j,k)=getDM(this%xi(i), this%yi05(j), this%zi(k), 311);
                  auxDMEJy(13,i,j,k)=getDM(this%xi(i), this%yi05(j), this%zi(k), 211);
                  auxDMEJy(14,i,j,k)=getDM(this%xi(i), this%yi05(j), this%zi(k), 332);
               enddo
            enddo
         enddo
         do k=1,this%Nz
            do j=0,this%Ny
               do i=0,this%Nx
                  auxMEJz(i,j,k)=getM(this%xi(i), this%yi(j), this%zi05(k));
                  auxDMEJz(0,i,j,k)=getDM(this%xi(i), this%yi(j), this%zi05(k), 1);
                  auxDMEJz(1,i,j,k)=getDM(this%xi(i), this%yi(j), this%zi05(k), 2);
                  auxDMEJz(2,i,j,k)=getDM(this%xi(i), this%yi(j), this%zi05(k), 3);
                  auxDMEJz(3,i,j,k)=getDM(this%xi(i), this%yi(j), this%zi05(k), 31);
                  auxDMEJz(4,i,j,k)=getDM(this%xi(i), this%yi(j), this%zi05(k), 32);
                  auxDMEJz(5,i,j,k)=getDM(this%xi(i), this%yi(j), this%zi05(k), 21);
                  auxDMEJz(6,i,j,k)=getDM(this%xi(i), this%yi(j), this%zi05(k), 11);
                  auxDMEJz(7,i,j,k)=getDM(this%xi(i), this%yi(j), this%zi05(k), 22);
                  auxDMEJz(8,i,j,k)=getDM(this%xi(i), this%yi(j), this%zi05(k), 33);
                  auxDMEJz(9,i,j,k)=getDM(this%xi(i), this%yi(j), this%zi05(k), 222);
                  auxDMEJz(10,i,j,k)=getDM(this%xi(i), this%yi(j), this%zi05(k), 333);
                  auxDMEJz(11,i,j,k)=getDM(this%xi(i), this%yi(j), this%zi05(k), 322);
                  auxDMEJz(12,i,j,k)=getDM(this%xi(i), this%yi(j), this%zi05(k), 311);
                  auxDMEJz(13,i,j,k)=getDM(this%xi(i), this%yi(j), this%zi05(k), 211);
                  auxDMEJz(14,i,j,k)=getDM(this%xi(i), this%yi(j), this%zi05(k), 332);
               enddo
            enddo
         enddo
         
!!! H !!!
         do k=1,this%Nz
            do j=1,this%Ny
               do i=0,this%Nx
                  auxMHx(i,j,k)=getM(this%xi(i), this%yi05(j), this%zi05(k));
                  auxDMHx(0,i,j,k)=getDM(this%xi(i), this%yi05(j), this%zi05(k), 1);
                  auxDMHx(1,i,j,k)=getDM(this%xi(i), this%yi05(j), this%zi05(k), 2);
                  auxDMHx(2,i,j,k)=getDM(this%xi(i), this%yi05(j), this%zi05(k), 3);
                  auxDMHx(3,i,j,k)=getDM(this%xi(i), this%yi05(j), this%zi05(k), 31);
                  auxDMHx(4,i,j,k)=getDM(this%xi(i), this%yi05(j), this%zi05(k), 32);
                  auxDMHx(5,i,j,k)=getDM(this%xi(i), this%yi05(j), this%zi05(k), 21);
                  auxDMHx(6,i,j,k)=getDM(this%xi(i), this%yi05(j), this%zi05(k), 11);
                  auxDMHx(7,i,j,k)=getDM(this%xi(i), this%yi05(j), this%zi05(k), 22);
                  auxDMHx(8,i,j,k)=getDM(this%xi(i), this%yi05(j), this%zi05(k), 33);
                  auxDMHx(9,i,j,k)=getDM(this%xi(i), this%yi05(j), this%zi05(k), 222);
                  auxDMHx(10,i,j,k)=getDM(this%xi(i), this%yi05(j), this%zi05(k), 333);
                  auxDMHx(11,i,j,k)=getDM(this%xi(i), this%yi05(j), this%zi05(k), 322);
                  auxDMHx(12,i,j,k)=getDM(this%xi(i), this%yi05(j), this%zi05(k), 311);
                  auxDMHx(13,i,j,k)=getDM(this%xi(i), this%yi05(j), this%zi05(k), 211);
                  auxDMHx(14,i,j,k)=getDM(this%xi(i), this%yi05(j), this%zi05(k), 332);
               enddo
            enddo
         enddo
         do k=1,this%Nz
            do j=0,this%Ny
               do i=1,this%Nx
                  auxMHy(i,j,k)=getM(this%xi05(i), this%yi(j), this%zi05(k));
                  auxDMHy(0,i,j,k)=getDM(this%xi05(i), this%yi(j), this%zi05(k), 1);
                  auxDMHy(1,i,j,k)=getDM(this%xi05(i), this%yi(j), this%zi05(k), 2);
                  auxDMHy(2,i,j,k)=getDM(this%xi05(i), this%yi(j), this%zi05(k), 3);
                  auxDMHy(3,i,j,k)=getDM(this%xi05(i), this%yi(j), this%zi05(k), 31);
                  auxDMHy(4,i,j,k)=getDM(this%xi05(i), this%yi(j), this%zi05(k), 32);
                  auxDMHy(5,i,j,k)=getDM(this%xi05(i), this%yi(j), this%zi05(k), 21);
                  auxDMHy(6,i,j,k)=getDM(this%xi05(i), this%yi(j), this%zi05(k), 11);
                  auxDMHy(7,i,j,k)=getDM(this%xi05(i), this%yi(j), this%zi05(k), 22);
                  auxDMHy(8,i,j,k)=getDM(this%xi05(i), this%yi(j), this%zi05(k), 33);
                  auxDMHy(9,i,j,k)=getDM(this%xi05(i), this%yi(j), this%zi05(k), 222);
                  auxDMHy(10,i,j,k)=getDM(this%xi05(i), this%yi(j), this%zi05(k), 333);
                  auxDMHy(11,i,j,k)=getDM(this%xi05(i), this%yi(j), this%zi05(k), 322);
                  auxDMHy(12,i,j,k)=getDM(this%xi05(i), this%yi(j), this%zi05(k), 311);
                  auxDMHy(13,i,j,k)=getDM(this%xi05(i), this%yi(j), this%zi05(k), 211);
                  auxDMHy(14,i,j,k)=getDM(this%xi05(i), this%yi(j), this%zi05(k), 332);
               enddo
            enddo
         enddo
         do k=0,this%Nz
            do j=1,this%Ny
               do i=1,this%Nx
                  auxMHz(i,j,k)=getM(this%xi05(i), this%yi05(j), this%zi(k));
                  auxDMHz(0,i,j,k)=getDM(this%xi05(i), this%yi05(j), this%zi(k), 1);
                  auxDMHz(1,i,j,k)=getDM(this%xi05(i), this%yi05(j), this%zi(k), 2);
                  auxDMHz(2,i,j,k)=getDM(this%xi05(i), this%yi05(j), this%zi(k), 3);
                  auxDMHz(3,i,j,k)=getDM(this%xi05(i), this%yi05(j), this%zi(k), 31);
                  auxDMHz(4,i,j,k)=getDM(this%xi05(i), this%yi05(j), this%zi(k), 32);
                  auxDMHz(5,i,j,k)=getDM(this%xi05(i), this%yi05(j), this%zi(k), 21);
                  auxDMHz(6,i,j,k)=getDM(this%xi05(i), this%yi05(j), this%zi(k), 11);
                  auxDMHz(7,i,j,k)=getDM(this%xi05(i), this%yi05(j), this%zi(k), 22);
                  auxDMHz(8,i,j,k)=getDM(this%xi05(i), this%yi05(j), this%zi(k), 33);
                  auxDMHz(9,i,j,k)=getDM(this%xi05(i), this%yi05(j), this%zi(k), 222);
                  auxDMHz(10,i,j,k)=getDM(this%xi05(i), this%yi05(j), this%zi(k), 333);
                  auxDMHz(11,i,j,k)=getDM(this%xi05(i), this%yi05(j), this%zi(k), 322);
                  auxDMHz(12,i,j,k)=getDM(this%xi05(i), this%yi05(j), this%zi(k), 311);
                  auxDMHz(13,i,j,k)=getDM(this%xi05(i), this%yi05(j), this%zi(k), 211);
                  auxDMHz(14,i,j,k)=getDM(this%xi05(i), this%yi05(j), this%zi(k), 332);
               enddo
            enddo
         enddo
         write(*,*) 'Aux mike source matrices inited!';
      endif
    end subroutine initMikeSource;
    
!problem_report------------------------------------------------------------------        
    subroutine problem_report(this)
    !Write to the screen all the parameters of the problem
      class(tProblem) :: this
      integer :: k;
      write(*,*) '------PROBLEM REPORT-----'
      write(*,*) 'Problem ID:', this%problem_id;
      write(*,*) 'Problem Type:', this%problem_type;
      write(*,*) 'Scheme Type:', this%schemetype;
      write(*,*) 'Currents Type:', this%currentstype;
      write(*,*) 'Effective Currents Type:', this%effcurrentstype;
      write(*,*) 'E boundary:', this%Eboundarytype;
      write(*,*) 'H boundary:', this%Hboundarytype;
      write(*,*) 'Nx:', this%Nx;
      write(*,*) 'Ny:', this%Ny;
      write(*,*) 'Nz:', this%Nz;
      write(*,*) 'Nt:', this%Nt;
      write(*,*) 'hhx:', this%hhx;
      write(*,*) 'hhy:', this%hhy;
      write(*,*) 'hhz:', this%hhz;
      write(*,*) 'ht:', this%ht;
      write(*,*) 'Xi(0)', this%xi(0);
      write(*,*) 'Xi(Nx)', this%xi(this%Nx);
      write(*,*) 'Yi(0)', this%yi(0);
      write(*,*) 'Yi(Nz)', this%yi(this%Ny);
      write(*,*) 'Zi(0)', this%zi(0);
      write(*,*) 'Zi(Nz)', this%zi(this%Nz);
      write(*,*) 'Ti(0)', this%ti(0);
      write(*,*) 'Ti(Nt)', this%ti(this%Nt);

      !write(*,*) '====Xi===='
      !do k=0,this%Nx
      !   write(*,*) this%Xi(k);   
      !enddo
      !write(*,*) '====Xi05===='
      !do k=0,this%Nx+1
      !   write(*,*) this%Xi05(k);   
      !enddo
      !write(*,*) '====Yi===='
      !do k=0,this%Ny
      !   write(*,*) this%Yi(k);   
      !enddo
      !write(*,*) '====Yi05===='
      !do k=0,this%Ny+1
      !   write(*,*) this%Yi05(k);   
      !enddo
      !write(*,*) '====Zi===='
      !do k=0,this%Nz
      !   write(*,*) this%Zi(k);   
      !enddo
      !write(*,*) '====Zi05===='
      !do k=0,this%Nz+1
      !   write(*,*) this%Zi05(k);   
      !enddo
      !write(*,*) '====Ti===='
      !do k=0,this%Nt
      !   write(*,*) this%Ti(k);   
      !enddo
      !write(*,*) '====Ti05===='
      !do k=0,this%Nt
      !   write(*,*) this%Ti05(k);   
      !enddo
      
      write(*,*) 'Buffer Offset X:', this%bufoffsetx;
      write(*,*) 'Buffer Offset Y:', this%bufoffsety;
      write(*,*) 'Buffer Offset Z:', this%bufoffsetz;

      write(*,*) 'Main dX:', this%maindx;
      write(*,*) 'Main dY:', this%maindy;
      write(*,*) 'Main dZ:', this%maindz;

      write(*,*) 'AuxMaxTime:', this%auxmaxtime;
      write(*,*) 'SigmaS:', this%sgm_s;


      if (this%problem_type==0) then

         write(*,*) 'Namu', this%Namu
         write(*,*) 'Ndmu', this%Ndmu
         write(*,*) 'Nlc', this%Nlc
         write(*,*) 'Dlc', this%Dlc

         if (this%Eboundarytype==4)then
            write(*,*) 'Poisson'
            write(*,*) '  Poisson npx:', this%poisson%npx;
            write(*,*) '  Poisson npy:', this%poisson%npy;
            write(*,*) '  Poisson npz:', this%poisson%npz;
            write(*,*) '  Poisson Nin:', this%poisson%Nin;
            write(*,403) this%poisson%Nx, this%poisson%Ny, this%poisson%Nz;
403         format ('   Poisson Nx, Ny, Nz: ',I4,', ',I4,', ',I4)
            write(*,404) this%poisson%hhx, this%poisson%hhy, this%poisson%hhz, this%poisson%ht;
404         format ('   Poisson hhx, hhy, hhz, ht: ',E14.7,', ',E14.7,', ',E14.7,', ',E14.7)
         endif
      endif

      if ((this%Eboundarytype==3).OR.(this%Eboundarytype==10).OR.(this%Eboundarytype==14).OR.(this%Eboundarytype==16)) then
         write(*,*) 'PML'
         write(*,*) '  PML ID:', this%PML%pml_id;
         write(*,*) '  PML E boundary:', this%PML%EboundarytypePML;
         write(*,*) '  PML H boundary:', this%PML%HboundarytypePML;
         write(*,*) '  PML NumPml:', this%PML%NumPml;
         write(*,*) '  PML Apml:', this%PML%Apml;
         write(*,405) this%PML%Nx, this%PML%Ny, this%PML%Nz;
         405 format ('   PML Nx, Ny, Nz, Nt: ',I4,', ',I4,', ',I4)
         write(*,406) this%PML%hhx, this%PML%hhy, this%PML%hhz, this%PML%ht;
         406      format ('   PML hhx, hhy, hhz, ht: ',E10.3,', ',E10.3,', ',E10.3,', ',E10.3)
         write(*,*) '  PML Currents type', this%pmlcurrentstype;
         write(*,*) '  PMLXi(0)', this%PML%xi(0);
         write(*,*) '  PMLXi(Nx)', this%PML%xi(this%PML%Nx);
         write(*,*) '  PMLYi(0)', this%PML%yi(0);
         write(*,*) '  PMLYi(Nz)', this%PML%yi(this%PML%Ny);
         write(*,*) '  PMLZi(0)', this%PML%zi(0);
         write(*,*) '  PMLZi(Nz)', this%PML%zi(this%PML%Nz);

         !write(*,*) '====PML Xi===='
         !do k=0,this%PML%Nx
         !  write(*,*) this%PML%Xi(k);   
         !enddo
         !write(*,*) '====PML Xi05===='
         !do k=0,this%PML%Nx+1
         !  write(*,*) this%PML%Xi05(k);   
         !enddo
         
         write(*,*) '  PMLsgmxi(0)', this%PML%sgmxi(0);
         write(*,*) '  PMLsgmxi(Nx)', this%PML%sgmxi(this%PML%Nx);
         write(*,*) '  PMLsgmyi(0)', this%PML%sgmyi(0);
         write(*,*) '  PMLsgmyi(Nz)', this%PML%sgmyi(this%PML%Ny);
         write(*,*) '  PMLsgmzi(0)', this%PML%sgmzi(0);
         write(*,*) '  PMLsgmzi(Nz)', this%PML%sgmzi(this%PML%Nz);

         !write(*,*) '====PML Sigma Xi===='
         !do k=0,this%PML%Nx
         !  write(*,*) this%PML%sgmXi(k);   
         !enddo
         !write(*,*) '====PML Sigma Xi05===='
         !do k=0,this%PML%Nx
         !  write(*,*) this%PML%sgmXi05(k);   
         !enddo

         !write(*,*) '====PML Sigma Yi===='
         !do k=0,this%PML%Ny
         !  write(*,*) this%PML%sgmYi(k);   
         !enddo
         !write(*,*) '====PML Sigma Yi05===='
         !do k=0,this%PML%Ny
         !  write(*,*) this%PML%sgmYi05(k);   
         !enddo

         !write(*,*) '====PML Sigma Zi===='
         !do k=0,this%PML%Nz
         !  write(*,*) this%PML%sgmZi(k);   
         !enddo
         !write(*,*) '====PML Sigma Zi05===='
         !do k=0,this%PML%Nz
         !  write(*,*) this%PML%sgmZi05(k);   
         !enddo
      endif

      ! Write Higdon angles
      if ((this%Eboundarytype==6).OR.(this%Eboundarytype==7).OR.(this%Eboundarytype==15).OR.(this%Eboundarytype==16)) then
         write(*,*) 'Higdon angle 1: ', this%hig1;
         write(*,*) 'Higdon angle 2: ', this%hig2;
      endif

      ! Write Givoli-Neta Angles
      if ((this%Eboundarytype==9).OR.(this%Eboundarytype==10).OR.(this%Eboundarytype==11)) then
         write(*,*) 'Givoli-Neta Aux Variables: ', this%GNnum;
         do k=1,this%GNnum
            write(*,501) k, this%gnA(k);
         enddo   
      endif
501   format (' Givoli-Neta Angle(',I2,'): ', E14.7);

      ! Write Hagstrom-Warburton Angles
      if (this%Eboundarytype==13.OR.this%Eboundarytype==14) then
         write(*,*) 'Hagstrom-Warburton Aux Variables: ', this%hwNum;
         do k=0,this%hwNum-1
            write(*,502) k, this%hwAng(k);
         enddo
      endif
502   format (' Hagstrom-Warburton Angle(',I2,'): ', E14.7);

      ! Write HWEva Angles and Sigmas
      if (this%Eboundarytype==15.OR.(this%Eboundarytype==16)) then
         write(*,*) 'HWEva Aux Angles: ', this%hwEvaP;
         write(*,*) 'HWEva Aux Sigmas: ', this%hwEvaE;
         do k=0,this%hwEvaP
            write(*,503) k, this%hwAng(k);
         enddo
         do k=1,this%hwEvaE
            write(*,504) k, this%hwSigmas(k);
         enddo


         !write(*,*) 'HWEva Aux Aux Angles: ', this%hwEvaPa;
         !write(*,*) 'HWEva Aux Aux Sigmas: ', this%hwEvaEa;
         !do k=0,this%hwEvaPa
         !   write(*,505) k, this%hwAngA(k);
         !enddo
         !do k=1,this%hwEvaEa
         !   write(*,506) k, this%hwSigmasA(k);
         !enddo
      endif

503   format (' HWEva Angle(',I2,'): ', E14.7);
504   format (' HWEva Sigma(',I2,'): ', E14.7);
505   format (' HWEva Aux Angle(',I2,'): ', E14.7);
506   format (' HWEva Aux Sigma(',I2,'): ', E14.7);

      ! Write difference aproximations modes for ABCs
      if ((this%Eboundarytype==9).OR.(this%Eboundarytype==10).OR.(this%Eboundarytype==11).OR.(this%Eboundarytype==13).OR.(this%Eboundarytype==14).OR.(this%Eboundarytype==15).OR.(this%Eboundarytype==16)) then
         write(*,*) 'DT approximation order: ', this%dtmode;
         write(*,*) 'DX approximation order: ', this%dxmode;
         write(*,*) 'DT approximation order in aux functions: ', this%dtauxmode;
         write(*,*) 'DX approximation order in aux functions: ', this%dxauxmode;
      endif
      
      write(*,*) '------End Problem Report-----'
    end subroutine problem_report
    
    

!problem_destroy------------------------------------------------------------------      
    subroutine problem_Destroy(this)
    !Destroy the problem structures and free the memory  
      class(tProblem) :: this
      call this%Ef%mesh_destroy;
      call this%Hf%mesh_destroy;
      call this%Ean%mesh_destroy;
      call this%Han%mesh_destroy;
      call this%Je%mesh_destroy;
      call this%Jh%mesh_destroy;
      call this%Efold%mesh_destroy;
      call this%Hfold%mesh_destroy;
      call this%Efold2%mesh_destroy;
      call this%Hfold2%mesh_destroy;
      call this%JEeff%mesh_destroy;
      call this%JHeff%mesh_destroy;
      if ((this%Eboundarytype==3).OR.(this%Eboundarytype==10).OR.(this%Eboundarytype==14).OR.(this%Eboundarytype==16)) then
         call this%PML%pml_destroy;
      endif;
      if (this%problem_type==0.AND.this%effcurrentstype==1) then
         call this%poisson%poisson_destroy;
      endif
      if ((this%Eboundarytype==9).OR.(this%Eboundarytype==10).OR.(this%Eboundarytype==11)) then
         call this%destroyGivoliNetaAuxVariables;
      endif;
      if (this%Eboundarytype==13.OR.this%Eboundarytype==14) then
         call this%destroyHWAuxVariables;
      endif;
      if (this%Eboundarytype==15.OR.this%Eboundarytype==16) then
         call this%destroyHWEvaAuxVariables;
      endif;
      deallocate(this%xi,this%yi,this%zi, this%ti, this%xi05, this%yi05, this%zi05, this%ti05, this%mup);

    end subroutine problem_Destroy;

!initPhysPoints------------------------------------------------------------------      
    subroutine initPhysPoints(this)
    !Initialize arrays with physical points
      class(tProblem) :: this
      integer k;      
      do k=0,this%Nt
         this%Ti(k)=this%ht*k;    
      enddo
      do k=0,this%Nt+1
         this%Ti05(k)=this%ht*(k-0.5d0);
      enddo     
      do k=0,this%Nx
         this%Xi(k)=-this%Nx*this%hhx*cntx+this%hhx*(k);
         this%Xi05(k)=-this%Nx*this%hhx*cntx+this%hhx*(k-0.5d0);    
      enddo
      this%Xi05(this%Nx+1)=-this%Nx*this%hhx*cntx+this%hhx*(this%Nx+1.0d0-0.5d0);  
      do k=0,this%Ny
         this%Yi(k)=-this%Ny*this%hhy*cnty+this%hhy*(k);
         this%Yi05(k)=-this%Ny*this%hhy*cnty+this%hhy*(k-0.5d0);
      enddo
      this%Yi05(this%Ny+1)=-this%Ny*this%hhy*cnty+this%hhy*(this%Ny+1.0d0-0.5d0);  
      do k=0,this%Nz
         this%Zi(k)=-this%Nz*this%hhz*cntz+this%hhz*(k);
         this%Zi05(k)=-this%Nz*this%hhz*cntz+this%hhz*(k-0.5d0);
      enddo
      this%Zi05(this%Nz+1)=-this%Nz*this%hhz*cntz+this%hhz*(this%Nz+1.0d0-0.5d0);
    end subroutine initPhysPoints;
    

!Problem_DoStep------------------------------------------------------------------      
    subroutine problem_DoIndependentStep(this, t)
    !Do one time step of the propagation
      class(tProblem) :: this
      integer, intent(in) :: t;
      integer :: k;
     
      if (this%EboundaryType == 10.OR.this%EboundaryType == 14.OR.this%EboundaryType == 16) then
         call this%getSourceCurrentsPML(t);
         call this%PML%updateEInteriorPMLEverywhere;
         call this%PML%updateHInteriorPMLEverywhere;
      endif

           call this%starttimer(0, 10);
      call this%getSourceCurrents(t);          ! Fill J arrays with actual values of currents
           call this%endtimer(0, 10);
      

            call this%starttimer(18, 19);
      if (this%EboundaryType == 6.OR.this%EboundaryType == 7.OR.this%EboundaryType == 8.OR.this%EboundaryType == 9&
           &.OR.this%EboundaryType == 10.OR.this%EboundaryType == 11.OR.this%EboundaryType == 13.OR.this%EboundaryType == 14.OR.this%EboundaryType == 15.OR.this%EboundaryType == 16) then
         call this%saveFields(1);
      endif
      call this%saveFields(0);
            call this%endtimer(18, 19);
      
      if (this%EboundaryType == 9.OR.this%EboundaryType == 10.OR.this%EboundaryType == 11) then
         call this%saveGivoliNetaAuxVariables
      endif
      if (this%EboundaryType == 13.OR.this%EboundaryType == 14) then
         call this%saveHWAuxVariables
      endif
      if (this%EboundaryType == 15.OR.this%EboundaryType == 16) then
         call this%saveHWEvaAuxVariables
      endif
    
           call this%starttimer(1, 11);
      call this%updateEInterior;               ! Calculate E field everywhere except outer boundary
           call this%endtimer(1, 11);
      
           call this%starttimer(2, 12);
      call this%updateHInterior;               ! Calculate H field everywhere except outer boundary
           call this%endtimer(2, 12);
      
           call this%starttimer(3, 13);                        
      if ((this%Eboundarytype==3).OR.(this%Hboundarytype==3)) then
         call this%uploadFieldsToPML;          ! Transfer field values from the problem to PML
         call this%PML%pml_DoStep(t);          ! PML single step
      endif;
      select case (this%EboundaryType)         ! Fill E boundary with current Eboundary condition
      case (5)
         call this%fillEBoundaryABCSommerfeld(t);              
      case (6)
         call this%fillEBoundaryABCHigdon(t);             
      case (7)
         call this%fillEBoundaryABCBetzMittra(t);
      case (8)
         call this%fillEBoundaryABCMur(t);
      case (9)
         call this%fillEBoundaryABCGivoliNeta(t);
      case (10)
         call this%fillEBoundaryABCGivoliNetaPML(t);
      case (11)
         call this%fillEBoundaryABCGivoliNetaPure(t);
      case (13)
         call this%fillEBoundaryABCHW(t);
      case (14)
         call this%fillEBoundaryABCHWPML(t);
      case (15)
         call this%fillEBoundaryABCHWEva(t);
      case (16)
         call this%fillEBoundaryABCHWEvaPML(t); 
      case default
         call this%fillEBoundary(t);             
      end select
           call this%endtimer(3, 13);    

      call this%starttimer(4, 14);  
           call this%fillHBoundary(t);              ! Fill H boundary with current Hboundary condition
      call this%endtimer(4, 14);  

      do k=0,31
         times(k) = times(k) + (tends(k)-tstarts(k));
         tstarts(k) = 0.0d0;
         tends(k) = 0.0d0;
      enddo
       
    end subroutine problem_DoIndependentStep;

!Problem_StartTimer------------------------------------------------------------------      
    subroutine starttimer(this, t1, t2)
      class(tProblem) :: this
      integer :: t1, t2;
      if (this%problem_type==0) then
         call cpu_time(tstarts(t1))
         !$ tstarts(t1) = OMP_get_wtime();
      else
         call cpu_time(tstarts(t2))
         !$ tstarts(t2) = OMP_get_wtime();
      endif
    end subroutine starttimer;

!Problem_EndTimer------------------------------------------------------------------      
    subroutine endtimer(this, t1, t2)
      class(tProblem) :: this
      integer :: t1, t2;
      if (this%problem_type==0) then
         call cpu_time(tends(t1))
         !$ tends(t1) = OMP_get_wtime();  
      else
         call cpu_time(tends(t2))
         !$ tends(t2) = OMP_get_wtime();  
      endif
    end subroutine endtimer;     

!Problem_DoIndependentStepTiming------------------------------------------------------------------      
    subroutine problem_DoIndependentStepTiming(this, t)
    !Do one time step of the propagation
      class(tProblem) :: this
      integer, intent(in) :: t;
      integer :: k;
      if (checkupdated==1) then
         call this%Ef%mesh_clearChA;           ! Clear checking arrays
         call this%Hf%mesh_clearChA;           ! Clear checking arrays
         if ((this%Eboundarytype==3).OR.(this%Hboundarytype==3)) then
            call this%PML%Ef%mesh_clearChA;       ! Clear checking arrays
            call this%PML%Hf%mesh_clearChA;       ! Clear checking arrays
         endif
      endif

      call cpu_time(tstarts(24))
      !$ tstarts(24) = OMP_get_wtime();
         call this%getSourceCurrents(t);          ! Fill J arrays with actual values of currents
         call cpu_time(tends(24))
      !$ tends(24) = OMP_get_wtime();
      
      call cpu_time(tstarts(25))
      !$ tstarts(25) = OMP_get_wtime();       
         call this%saveFields(0);
      call cpu_time(tends(25))
      !$ tends(25) = OMP_get_wtime();   
      
      call cpu_time(tstarts(26))
      !$ tstarts(26) = OMP_get_wtime();      
         call this%updateEInterior;               ! Calculate E field everywhere except outer boundary
      call cpu_time(tends(26))
      !$ tends(26) = OMP_get_wtime();   

      call cpu_time(tstarts(27))
      !$ tstarts(27) = OMP_get_wtime();      
         call this%updateHInterior;               ! Calculate H field everywhere except outer boundary
      call cpu_time(tends(27))
      !$ tends(27) = OMP_get_wtime();   

      if ((this%Eboundarytype==3).OR.(this%Hboundarytype==3)) then
         call cpu_time(tstarts(28))
         !$ tstarts(28) = OMP_get_wtime();      
            call this%uploadFieldsToPML;          ! Transfer field values from the problem to PML
         call cpu_time(tends(28))
         !$ tends(28) = OMP_get_wtime();   

         call cpu_time(tstarts(29))
         !$ tstarts(29) = OMP_get_wtime();      
            call this%PML%pml_DoStep(t);          ! PML single step
         call cpu_time(tends(29))
         !$ tends(29) = OMP_get_wtime();   
      endif; 

      call cpu_time(tstarts(30))
      !$ tstarts(30) = OMP_get_wtime();      
         call this%fillEBoundary(t);              ! Fill E boundary with current Eboundary condition
      call cpu_time(tends(30))
      !$ tends(30) = OMP_get_wtime();    

      call cpu_time(tstarts(31))
      !$ tstarts(31) = OMP_get_wtime();  
         call this%fillHBoundary(t);              ! Fill H boundary with current Hboundary condition
      call cpu_time(tends(31))
      !$ tends(31) = OMP_get_wtime();

      if (checkupdated==1) then
         if ((this%Eboundarytype==3).OR.(this%Hboundarytype==3)) then
            call this%PML%Ef%mesh_checkarrays;    ! Check all the points are updated once
            call this%PML%Hf%mesh_checkarrays;    ! Check all the points are updated once
         endif
         call this%Ef%mesh_checkarrays;        ! Check all the points are updated once
         call this%Hf%mesh_checkarrays;        ! Check all the points are updated once
      endif
    end subroutine Problem_DoIndependentStepTiming;
    

!updateEInterior------------------------------------------------------------------
    subroutine updateEInterior(this)
    ! Update all the E fields in the interior
      class(tProblem) :: this
      integer :: i, j ,k
      !$OMP PARALLEL DO SCHEDULE(GUIDED)
      do k=1,this%Nz-1
         do j=1,this%Ny-1
            do i=1,this%Nx
               call this%updateEpoint(i, j, k, 1); ! Ex update
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
      !$OMP PARALLEL DO SCHEDULE(GUIDED)
      do k=1,this%Nz-1
         do j=1,this%Ny
            do i=1,this%Nx-1
               call this%updateEpoint(i, j, k, 2); ! Ey update
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
      !$OMP PARALLEL DO SCHEDULE(GUIDED)
      do k=1,this%Nz
         do j=1,this%Ny-1
            do i=1,this%Nx-1
               call this%updateEpoint(i, j, k, 3); ! Ez update
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
    end subroutine updateEInterior;


!updateHInterior------------------------------------------------------------------
    subroutine updateHInterior(this)
    ! Update all the H fields in the interior
      class(tProblem) :: this
      integer :: i, j ,k
      !$OMP PARALLEL DO SCHEDULE(GUIDED)
      do k=2,this%Nz-1
         do j=2,this%Ny-1
            do i=1,this%Nx-1
               call this%updateHpoint(i, j, k, 1); ! Hx update
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
      !$OMP PARALLEL DO SCHEDULE(GUIDED)
      do k=2,this%Nz-1
         do j=1,this%Ny-1
            do i=2,this%Nx-1
               call this%updateHpoint(i, j, k, 2); ! Hy update
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
      !$OMP PARALLEL DO SCHEDULE(GUIDED)
      do k=1,this%Nz-1
         do j=2,this%Ny-1
            do i=2,this%Nx-1
               call this%updateHpoint(i, j, k, 3); ! Hz update
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
    end subroutine updateHInterior;  

    
!updateEpoint------------------------------------------------------------------
    subroutine updateEpoint(this, Px, Py, Pz, component)
    ! Update E point with scheme
      class(tProblem) :: this
      integer, intent(in) :: Px, Py, Pz, component;
      integer :: inc;
      inc = 0;
      select case(component)
      case(1) ! EX update
     !    if (checkbounds==1) then
     !       ! Check all the necessary points are correct
     !       inc = inc + this%Ef%mesh_checkpoint(Px, Py, Pz, component);
     !       inc = inc + this%Hf%mesh_checkpoint(Px, Py+1, Pz, 3);
     !       inc = inc + this%Hf%mesh_checkpoint(Px, Py, Pz, 3);
     !       inc = inc + this%Hf%mesh_checkpoint(Px, Py, Pz+1, 2);
     !       inc = inc + this%Hf%mesh_checkpoint(Px, Py, Pz, 2);
     !    endif
     !    if (inc==0) then
            this%Ef%X(Px,Py,Pz)=this%Ef%X(Px,Py,Pz) + cc*this%ht*( (this%Hf%Z(Px,Py+1,Pz)-this%Hf%Z(Px,Py,Pz))/this%hhy - (this%Hf%Y(Px,Py,Pz+1)-this%Hf%Y(Px,Py,Pz))/this%hhz )-this%ht*this%Je%X(Px,Py,Pz)*4.0d0*Pi/cc;
     !       if (checkupdated==1) then
     !          this%Ef%chX(Px,Py,Pz) = this%Ef%chX(Px,Py,Pz)+1;   ! this string is important for update checking
     !       endif     
     !    else
     !       fatalerr = 201;
     !       write(*,*) 'Alert! Ex point update bounds error.'
     !    endif 
      case(2) ! EY update
     !    if (checkbounds==1) then
     !      ! Check all the necessary points are correct
     !       inc = inc + this%Ef%mesh_checkpoint(Px, Py, Pz, component);
     !       inc = inc + this%Hf%mesh_checkpoint(Px, Py, Pz+1, 1);
     !       inc = inc + this%Hf%mesh_checkpoint(Px, Py, Pz, 1);
     !       inc = inc + this%Hf%mesh_checkpoint(Px+1, Py, Pz, 3);
     !       inc = inc + this%Hf%mesh_checkpoint(Px, Py, Pz, 3);
     !    endif     
     !    if (inc==0) then
            this%Ef%Y(Px,Py,Pz)=this%Ef%Y(Px,Py,Pz) + cc*this%ht*( (this%Hf%X(Px,Py,Pz+1)-this%Hf%X(Px,Py,Pz))/this%hhz - (this%Hf%Z(Px+1,Py,Pz)-this%Hf%Z(Px,Py,Pz))/this%hhx )-this%ht*this%Je%Y(Px,Py,Pz)*4.0d0*Pi/cc;
     !       if (checkupdated==1) then
     !          this%Ef%chY(Px,Py,Pz) = this%Ef%chY(Px,Py,Pz)+1;   ! this string is important for update checking
     !       endif
     !    else
     !       fatalerr = 201;
     !       write(*,*) 'Alert! Ey point update bounds error.'
     !    endif 
      case(3) ! EZ update
     !    if (checkbounds==1) then
            ! Check all the necessary points are correct
     !       inc = inc + this%Ef%mesh_checkpoint(Px, Py, Pz, component);
     !       inc = inc + this%Hf%mesh_checkpoint(Px+1, Py, Pz, 2);
     !       inc = inc + this%Hf%mesh_checkpoint(Px, Py, Pz, 2);
     !       inc = inc + this%Hf%mesh_checkpoint(Px, Py+1, Pz, 1);
     !       inc = inc + this%Hf%mesh_checkpoint(Px, Py, Pz, 1);
     !    endif
     !    if (inc==0) then
            this%Ef%Z(Px,Py,Pz)=this%Ef%Z(Px,Py,Pz) + cc*this%ht*( (this%Hf%Y(Px+1,Py,Pz)-this%Hf%Y(Px,Py,Pz))/this%hhx - (this%Hf%X(Px,Py+1,Pz)-this%Hf%X(Px,Py,Pz))/this%hhy )-this%ht*this%Je%Z(Px,Py,Pz)*4.0d0*Pi/cc;
     !       if (checkupdated==1) then
     !          this%Ef%chZ(Px,Py,Pz) = this%Ef%chZ(Px,Py,Pz)+1;   ! this string is important for update checking
     !       endif   
     !    else
     !       fatalerr = 201;
     !       write(*,*) 'Alert! Ez point update bounds error.'
     !    endif
      case default ! Incorrect
         fatalerr = 203;
         write(*,*) 'Alert! Incorrect use of Component in updateEpoint function!', component
      end select
    end subroutine updateEpoint;   


!CheckEPointUpdated------------------------------------------------------------------ 
   integer function checkEpointUpdated(this, Px, Py, Pz, component) 
   !Check that point in E fields is updated only one time
     class(tProblem) :: this
     integer :: Px, Py, Pz, component, res
     res = 0;
     select case (component)
     case (1)
        if (this%Ef%chX(Px,Py,Pz)/=1) then
           res=1;
        endif
     case (2)
        if (this%Ef%chY(Px,Py,Pz)/=1) then
           res=1;
        endif
     case (3)
        if (this%Ef%chZ(Px,Py,Pz)/=1) then
           res=1;
        endif
     end select
     checkEpointUpdated = res;
   end function checkEpointUpdated;   

    
!updateHpoint------------------------------------------------------------------
    subroutine updateHpoint(this, Px, Py, Pz, component)
    !Initialize arrays with physical points
      class(tProblem) :: this
      integer, intent(in) :: Px, Py, Pz, component;
      integer :: inc;
      inc = 0;
      select case(component)
      case(1) ! Hx update
        ! if (checkbounds==1) then
        !    ! Check all the necessary points are correct
        !    inc = inc + this%Hf%mesh_checkpoint(Px, Py, Pz, component);
        !    inc = inc + this%Ef%mesh_checkpoint(Px, Py, Pz, 3);
        !    inc = inc + this%Ef%mesh_checkpoint(Px, Py-1, Pz, 3);
        !    inc = inc + this%Ef%mesh_checkpoint(Px, Py, Pz, 2);
        !    inc = inc + this%Ef%mesh_checkpoint(Px, Py, Pz-1, 2);
        ! endif 
        ! if (inc==0) then
            this%Hf%X(Px,Py,Pz)=this%Hf%X(Px,Py,Pz) - cc*this%ht*( (this%Ef%Z(Px,Py,Pz)-this%Ef%Z(Px,Py-1,Pz))/this%hhy - (this%Ef%Y(Px,Py,Pz)-this%Ef%Y(Px,Py,Pz-1))/this%hhz ) -this%ht*this%Jh%X(Px,Py,Pz)*4.0d0*Pi/cc;
        !    if (checkupdated==1) then
        !       this%Hf%chX(Px,Py,Pz) = this%Hf%chX(Px,Py,Pz)+1; ! this string is important for update checking
        !    endif   
        ! else
        !    fatalerr = 202;
        !    write(*,*) 'Alert! Hx point update bounds error.'
        ! endif   
      case(2) ! Hy update
        ! if (checkbounds==1) then
        !    ! Check all the necessary points are correct
        !    inc = inc + this%Hf%mesh_checkpoint(Px, Py, Pz, component);
        !    inc = inc + this%Ef%mesh_checkpoint(Px, Py, Pz, 1);
        !    inc = inc + this%Ef%mesh_checkpoint(Px, Py, Pz-1, 1);
        !    inc = inc + this%Ef%mesh_checkpoint(Px, Py, Pz, 3);
        !    inc = inc + this%Ef%mesh_checkpoint(Px-1, Py, Pz, 3);
        ! endif 
        ! if (inc==0) then
            this%Hf%Y(Px,Py,Pz)=this%Hf%Y(Px,Py,Pz) - cc*this%ht*( (this%Ef%X(Px,Py,Pz)-this%Ef%X(Px,Py,Pz-1))/this%hhz - (this%Ef%Z(Px,Py,Pz)-this%Ef%Z(Px-1,Py,Pz))/this%hhx ) - this%ht*this%Jh%Y(Px,Py,Pz)*4.0d0*Pi/cc;
        !    if (checkupdated==1) then
        !       this%Hf%chY(Px,Py,Pz) = this%Hf%chY(Px,Py,Pz)+1; ! this string is important for update checking
        !    endif   
        ! else
        !    fatalerr = 202;
        !    write(*,*) 'Alert! Hy point update bounds error.'
        ! endif   
      case(3) ! Hz update
        ! if (checkbounds==1) then
        !    ! Check all the necessary points are correct
        !    inc = inc + this%Hf%mesh_checkpoint(Px, Py, Pz, component);
        !    inc = inc + this%Ef%mesh_checkpoint(Px, Py, Pz, 2);
        !    inc = inc + this%Ef%mesh_checkpoint(Px-1, Py, Pz, 2);
        !    inc = inc + this%Ef%mesh_checkpoint(Px, Py, Pz, 1);
        !    inc = inc + this%Ef%mesh_checkpoint(Px, Py-1, Pz, 1);
        ! endif 
        ! if (inc==0) then
            this%Hf%Z(Px,Py,Pz)=this%Hf%Z(Px,Py,Pz) - cc*this%ht*( (this%Ef%Y(Px,Py,Pz)-this%Ef%Y(Px-1,Py,Pz))/this%hhx - (this%Ef%X(Px,Py,Pz)-this%Ef%X(Px,Py-1,Pz))/this%hhy ) - this%ht*this%Jh%Z(Px,Py,Pz)*4.0d0*Pi/cc;
        !    if (checkupdated==1) then
        !       this%Hf%chZ(Px,Py,Pz) = this%Hf%chZ(Px,Py,Pz)+1; ! this string is important for update checking
        !    endif  
        ! else
        !    fatalerr = 202;
        !    write(*,*) 'Alert! Hz point update bounds error.'
        ! endif
      case default ! Incorrect
         fatalerr = 204;
         write(*,*) 'Alert! Incorrect use of Component in updateHpoint function!', component   
      end select
    end subroutine updateHpoint;
    
    
!getCurrentEPoint------------------------------------------------------------------
    real*8 function getCurrentEPoint(this, Px, Py, Pz, t, component)
    !Fill Je sources with results from analytical solution
      class(tProblem) :: this
      integer, intent(in) :: Px, Py, Pz;
      integer, intent(in) :: t;
      integer, intent(in) :: component;
      real*8 :: res, t0, t1;
      real*8 :: argthetabig;
      real*8 :: muEx, muEy, muEz, ejE;
      res = 0;
      select case (this%currentstype)
      case (0) ! Analytic source  
         select case (component)
         case(1) ! Je X point
            res = this%Source%getJpoint(this%xi05(Px), this%yi(Py), this%zi(Pz), this%Ti05(t), Px, Py, Pz, 1);
         case(2) ! Je Y point
            res = this%Source%getJpoint(this%xi(Px), this%yi05(Py), this%zi(Pz), this%Ti05(t), Px, Py, Pz, 2);
         case(3) ! Je Z point
            res = this%Source%getJpoint(this%xi(Px), this%yi(Py), this%zi05(Pz), this%Ti05(t), Px, Py, Pz, 3);
         case default
            fatalerr = 209;
            write(*,*) 'Alert! Incorrect use of Component in getCurrentEPoint function!', component   
         end select
      case (1) ! Analytic source cutted with Theta function
         select case (component)
         case(1) ! Je X point
            res = this%Source%getJpoint(this%xi05(Px), this%yi(Py), this%zi(Pz), this%Ti05(t), Px, Py, Pz, 1);
         case(2) ! Je Y point
            res = this%Source%getJpoint(this%xi(Px), this%yi05(Py), this%zi(Pz), this%Ti05(t), Px, Py, Pz, 2);
         case(3) ! Je Z point
            res = this%Source%getJpoint(this%xi(Px), this%yi(Py), this%zi05(Pz), this%Ti05(t), Px, Py, Pz, 3);
         case default
            fatalerr = 209;
            write(*,*) 'Alert! Incorrect use of Component in getCurrentEPoint function!', component   
         end select
         argthetabig = this%Ti05(t)-(1.0+this%sgm_s)*this%auxmaxtime*this%problem_id;
         res = res * this%thetabig(argthetabig);
      case (2) ! Sources from buffer  
         select case (component)
         case(1) ! Je X point
            res = this%buffer%Je%X(Px,Py,Pz);
         case(2) ! Je Y point
            res = this%buffer%Je%Y(Px,Py,Pz);
         case(3) ! Je Z point
            res = this%buffer%Je%Z(Px,Py,Pz);
         case default
            fatalerr = 209;
            write(*,*) 'Alert! Incorrect use of Component in getCurrentEPoint function!', component   
         end select
      case (3) ! Sources from buffer cutted with Theta function 
         select case (component)
         case(1) ! Je X point
            res = this%buffer%Je%X(Px,Py,Pz);
         case(2) ! Je Y point
            res = this%buffer%Je%Y(Px,Py,Pz);
         case(3) ! Je Z point
            res = this%buffer%Je%Z(Px,Py,Pz);
         case default
            fatalerr = 209;
            write(*,*) 'Alert! Incorrect use of Component in getCurrentEPoint function!', component   
         end select
         argthetabig = this%Ti05(t)-(1.0d0+this%sgm_s)*this%auxmaxtime*this%problem_id;
         res = res * this%thetabig(argthetabig);
      case (4) ! Sources from buffer cutted with Theta function and Poisson correction
         argthetabig = this%Ti05(t)-(1.0d0+this%sgm_s)*this%auxmaxtime*this%problem_id;
         t1=this%thetabig(argthetabig+this%ht/2.0d0);
         t0=this%thetabig(argthetabig-this%ht/2.0d0);
         select case (component)
         case(1) ! Je X point
            if (Px>=2+this%maindx.AND.PX<=this%NX-2-this%maindx.AND.Py>=2+this%maindy.AND.PY<=this%Ny-2-this%maindy.AND.Pz>=2+this%maindz.AND.Pz<=this%Nz-2-this%maindz) then
               muEx=this%buffer%mainMu%X(Px, Py, Pz);
               ejE=(1.0d0-muEx)*(t1*this%buffer%Ec1%X(Px,Py,Pz)-t0*this%buffer%Ec0%X(Px,Py,Pz))/this%ht;
               res = this%buffer%Je%X(Px,Py,Pz)* this%thetabig(argthetabig) + ejE*cc/4.0d0/Pi;
            else
              res = 0.0d0;  
            endif
         case(2) ! Je Y point
            if (Px>=2+this%maindx.AND.PX<=this%NX-2-this%maindx.AND.Py>=2+this%maindy.AND.PY<=this%Ny-2-this%maindy.AND.Pz>=2+this%maindz.AND.Pz<=this%Nz-2-this%maindz) then
               muEy=this%buffer%mainMu%Y(Px, Py, Pz);
               ejE=(1.0d0-muEy)*(t1*this%buffer%Ec1%Y(Px,Py,Pz)-t0*this%buffer%Ec0%Y(Px,Py,Pz))/this%ht;
               res = this%buffer%Je%Y(Px,Py,Pz)* this%thetabig(argthetabig) + ejE*cc/4.0d0/Pi;
            else
              res = 0.0d0;  
            endif   
         case(3) ! Je Z point
            if (Px>=2+this%maindx.AND.PX<=this%NX-2-this%maindx.AND.Py>=2+this%maindy.AND.PY<=this%Ny-2-this%maindy.AND.Pz>=2+this%maindz.AND.Pz<=this%Nz-2-this%maindz) then
               muEz=this%buffer%mainMu%Z(Px, Py, Pz);
               ejE=(1.0-muEz)*(t1*this%buffer%Ec1%Z(Px,Py,Pz)-t0*this%buffer%Ec0%Z(Px,Py,Pz))/this%ht;
               res = this%buffer%Je%Z(Px,Py,Pz)* this%thetabig(argthetabig) + ejE*cc/4.0d0/Pi;     
            else
              res = 0.0d0;  
            endif
         case default
            fatalerr = 209;
            write(*,*) 'Alert! Incorrect use of Component in getCurrentEPoint function!', component   
         end select
      case default
         fatalerr = 210;
         write(*,*) 'Alert! Incorrect Current Type is used in getCurrentEPoint function!', this%currentstype
      end select   
      getCurrentEPoint = res;
    end function getCurrentEPoint;


!getCurrentHPoint------------------------------------------------------------------
    real*8 function getCurrentHPoint(this, Px, Py, Pz, t, component)
    !Fill Je sources with results from analytical solution
      class(tProblem) :: this
      integer, intent(in) :: Px, Py, Pz;
      integer, intent(in) :: t;
      integer, intent(in) :: component;
      real*8 :: res;
      real*8 :: argthetabig;
      res = 0;
      select case (this%currentstype)
      case (0) ! Analytic source is zero
         select case (component)
         case(1) ! Jh X point
            res = 0.0d0;
         case(2) ! Jh Y point
            res = 0.0d0;
         case(3) ! Jh Z point
            res = 0.0d0;
         case default
            fatalerr = 211;
            write(*,*) 'Alert! Incorrect use of Component in getCurrentHPoint function!', component   
         end select
      case (1) ! Analytic source is zero, but cutted with P4 function
         select case (component)
         case(1) ! Jh X point
            res = 0.0d0;
         case(2) ! Jh Y point
            res = 0.0d0; 
         case(3) ! Jh Z point
            res = 0.0d0;
         case default
            fatalerr = 211;
            write(*,*) 'Alert! Incorrect use of Component in getCurrentHPoint function!', component   
         end select
         argthetabig = this%Ti05(t)-(1.0+this%sgm_s)*this%auxmaxtime*this%problem_id;
         res = res * this%thetabig(argthetabig);
      case (2) ! Sources from buffer  
         select case (component)
         case(1) ! Jh X point
            res = this%buffer%Jh%X(Px,Py,Pz);
         case(2) ! Jh Y point
            res = this%buffer%Jh%Y(Px,Py,Pz);
         case(3) ! Jh Z point
            res = this%buffer%Jh%Z(Px,Py,Pz);
         case default
            fatalerr = 211;
            write(*,*) 'Alert! Incorrect use of Component in getCurrentHPoint function!', component   
         end select
      case (3) ! Sources from buffer cutted with Theta function 
         select case (component)
         case(1) ! Jh X point
            res = this%buffer%Jh%X(Px,Py,Pz);
         case(2) ! Jh Y point
            res = this%buffer%Jh%Y(Px,Py,Pz);
         case(3) ! Jh Z point
            res = this%buffer%Jh%Z(Px,Py,Pz);
         case default
            fatalerr = 211;
            write(*,*) 'Alert! Incorrect use of Component in getCurrentHPoint function!', component   
         end select
         argthetabig = this%Ti05(t)-(1.0d0+this%sgm_s)*this%auxmaxtime*this%problem_id;
         res = res * this%thetabig(argthetabig);
      case (4) ! Sources from buffer cutted with Theta function with Poisson correction (no difference as wihtout Poisson)
         argthetabig = this%Ti05(t)-(1.0+this%sgm_s)*this%auxmaxtime*this%problem_id;
         select case (component)
         case(1) ! Jh X point
            if (Px>=2+this%maindx.AND.PX<=this%NX-2-this%maindx.AND.Py>=2+this%maindy.AND.PY<=this%Ny-2-this%maindy.AND.Pz>=2+this%maindz.AND.Pz<=this%Nz-2-this%maindz) then
               res = this%buffer%Jh%X(Px,Py,Pz)*this%thetabig(argthetabig);
            else
               res = 0.0d0;
            endif
         case(2) ! Jh Y point
            if (Px>=2+this%maindx.AND.PX<=this%NX-2-this%maindx.AND.Py>=2+this%maindy.AND.PY<=this%Ny-2-this%maindy.AND.Pz>=2+this%maindz.AND.Pz<=this%Nz-2-this%maindz) then
               res = this%buffer%Jh%Y(Px,Py,Pz)*this%thetabig(argthetabig);
            else
               res = 0.0d0;
            endif   
         case(3) ! Jh Z point
            if (Px>=2+this%maindx.AND.PX<=this%NX-2-this%maindx.AND.Py>=2+this%maindy.AND.PY<=this%Ny-2-this%maindy.AND.Pz>=2+this%maindz.AND.Pz<=this%Nz-2-this%maindz) then
               res = this%buffer%Jh%Z(Px,Py,Pz)*this%thetabig(argthetabig);
            else
               res = 0.0d0;
            endif   
         case default
            fatalerr = 211;
            write(*,*) 'Alert! Incorrect use of Component in getCurrentHPoint function!', component   
         end select
      case default
         fatalerr = 212;
         write(*,*) 'Alert! Incorrect Current Type is used in getCurrentHPoint function!', this%currentstype
      end select   
      getCurrentHPoint = res;
    end function getCurrentHPoint;
    
    
!getSourceCurrents------------------------------------------------------------------
    subroutine getSourceCurrents(this, t)
    !Fill Jan sources with results from analytical solution
      class(tProblem) :: this
      integer :: i,j,k
      integer, intent(in) :: t;
      ! JE Currents
      !$OMP PARALLEL DO SCHEDULE(GUIDED)
      do k=0,this%Nz
         do j=0,this%Ny
            do i=1,this%Nx
               this%JE%X(i,j,k)=this%getCurrentEPoint(i,j,k,t,1); ! Je X
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
      !$OMP PARALLEL DO SCHEDULE(GUIDED)
      do k=0,this%Nz
         do j=1,this%Ny
            do i=0,this%Nx
               this%JE%Y(i,j,k)=this%getCurrentEPoint(i,j,k,t,2); ! Je Y
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
      !$OMP PARALLEL DO SCHEDULE(GUIDED)
      do k=1,this%Nz
         do j=0,this%Ny
            do i=0,this%Nx
               this%JE%Z(i,j,k)=this%getCurrentEPoint(i,j,k,t,3); ! Je Z
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
      ! JH Currents
      !$OMP PARALLEL DO SCHEDULE(GUIDED)
      do k=1,this%Nz
         do j=1,this%Ny
            do i=0,this%Nx
               this%JH%X(i,j,k)=this%getCurrentHPoint(i,j,k,t,1); ! Jh X
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
      !$OMP PARALLEL DO SCHEDULE(GUIDED)
      do k=1,this%Nz
         do j=0,this%Ny
            do i=1,this%Nx
               this%JH%Y(i,j,k)=this%getCurrentHPoint(i,j,k,t,2); ! Jh Y
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
      !$OMP PARALLEL DO SCHEDULE(GUIDED)
      do k=0,this%Nz
         do j=1,this%Ny
            do i=1,this%Nx
               this%JH%Z(i,j,k)=this%getCurrentHPoint(i,j,k,t,3); ! Jh Z
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
    end subroutine getSourceCurrents;

    
!getCurrentEPointPML------------------------------------------------------------------
    real*8 function getCurrentEPointPML(this, Sx, Sy, Sz, t, component)
    !Fill Jan sources with results from analytical solution
      class(tProblem) :: this
      integer, intent(in) :: Sx, Sy, Sz;
      integer :: Px, Py, Pz;
      integer, intent(in) :: t;
      integer, intent(in) :: component;
      real*8 :: res, t0, t1;
      real*8 :: argthetabig;
      real*8 :: muEx, muEy, muEz, ejE;
      res = 0;
      Px = Sx;
      Py = Sy;
      Pz = Sz;
      select case (this%pmlcurrentstype)
      case (0) ! Original source  
         select case (component)
         case(1) ! Je X point
            res = this%Source%getJpoint(this%PML%xi05(Px), this%PML%yi(Py), this%PML%zi(Pz), this%Ti05(t), Px, Py, Pz, 1);
         case(2) ! Je Y point
            res = this%Source%getJpoint(this%PML%xi(Px), this%PML%yi05(Py), this%PML%zi(Pz), this%Ti05(t), Px, Py, Pz, 2);
         case(3) ! Je Z point
            res = this%Source%getJpoint(this%PML%xi(Px), this%PML%yi(Py), this%PML%zi05(Pz), this%Ti05(t), Px, Py, Pz, 3);
         case default
            fatalerr = 309;
            write(*,*) 'Alert! Incorrect use of Component in getCurrentEPointPML function!', component   
         end select
      case (1) ! Original source cutted by Mu function
         select case (component)
         case(1) ! Je X point
            res = this%Source%getJpoint(this%PML%xi05(Px), this%PML%yi(Py), this%PML%zi(Pz), this%Ti05(t), Px, Py, Pz, 1);
         case(2) ! Je Y point
            res = this%Source%getJpoint(this%PML%xi(Px), this%PML%yi05(Py), this%PML%zi(Pz), this%Ti05(t), Px, Py, Pz, 2);
         case(3) ! Je Z point
            res = this%Source%getJpoint(this%PML%xi(Px), this%PML%yi(Py), this%PML%zi05(Pz), this%Ti05(t), Px, Py, Pz, 3);
         case default
            fatalerr = 309;
            write(*,*) 'Alert! Incorrect use of Component in getCurrentEPointPML function!', component   
         end select
         argthetabig = this%Ti05(t)-(1.0+this%sgm_s)*this%auxmaxtime*this%problem_id;
         res = res * this%thetabig(argthetabig); 
      case (2) ! Sources from buffer
         Px = Px - this%PML%NumPML;
         Py = Py - this%PML%NumPML;
         Pz = Pz - this%PML%NumPML;
         select case (component)
         case(1) ! Je X point
            res = this%buffer%Je%X(Px,Py,Pz);
         case(2) ! Je Y point
            res = this%buffer%Je%Y(Px,Py,Pz);
         case(3) ! Je Z point
            res = this%buffer%Je%Z(Px,Py,Pz);
         case default
            fatalerr = 209;
            write(*,*) 'Alert! Incorrect use of Component in getCurrentEPoint function!', component   
         end select
      case (3) ! Sources from buffer cutted with Theta function
         Px = Px - this%PML%NumPML;
         Py = Py - this%PML%NumPML;
         Pz = Pz - this%PML%NumPML;
         select case (component)
         case(1) ! Je X point
            res = this%buffer%Je%X(Px,Py,Pz);
         case(2) ! Je Y point
            res = this%buffer%Je%Y(Px,Py,Pz);
         case(3) ! Je Z point
            res = this%buffer%Je%Z(Px,Py,Pz);
         case default
            fatalerr = 209;
            write(*,*) 'Alert! Incorrect use of Component in getCurrentEPoint function!', component   
         end select
         argthetabig = this%Ti05(t)-(1.0d0+this%sgm_s)*this%auxmaxtime*this%problem_id;
         res = res * this%thetabig(argthetabig);
      case (4) ! Sources from buffer cutted with Theta function and Poisson correction
         Px = Px - this%PML%NumPML;
         Py = Py - this%PML%NumPML;
         Pz = Pz - this%PML%NumPML;
         argthetabig = this%Ti05(t)-(1.0d0+this%sgm_s)*this%auxmaxtime*this%problem_id;
         t1=this%thetabig(argthetabig+this%ht/2.0d0);
         t0=this%thetabig(argthetabig-this%ht/2.0d0);
         select case (component)
         case(1) ! Je X point
            if (Px>=2+this%maindx.AND.PX<=this%NX-2-this%maindx.AND.Py>=2+this%maindy.AND.PY<=this%Ny-2-this%maindy.AND.Pz>=2+this%maindz.AND.Pz<=this%Nz-2-this%maindz) then
               muEx=this%buffer%mainMu%X(Px, Py, Pz);
               ejE=(1.0d0-muEx)*(t1*this%buffer%Ec1%X(Px,Py,Pz)-t0*this%buffer%Ec0%X(Px,Py,Pz))/this%ht;
               res = this%buffer%Je%X(Px,Py,Pz)* this%thetabig(argthetabig) + ejE*cc/4.0d0/Pi;
            else
              res = 0.0d0;  
            endif
         case(2) ! Je Y point
            if (Px>=2+this%maindx.AND.PX<=this%NX-2-this%maindx.AND.Py>=2+this%maindy.AND.PY<=this%Ny-2-this%maindy.AND.Pz>=2+this%maindz.AND.Pz<=this%Nz-2-this%maindz) then
               muEy=this%buffer%mainMu%Y(Px, Py, Pz);
               ejE=(1.0d0-muEy)*(t1*this%buffer%Ec1%Y(Px,Py,Pz)-t0*this%buffer%Ec0%Y(Px,Py,Pz))/this%ht;
               res = this%buffer%Je%Y(Px,Py,Pz)* this%thetabig(argthetabig) + ejE*cc/4.0d0/Pi;
            else
              res = 0.0d0;  
            endif   
         case(3) ! Je Z point
            if (Px>=2+this%maindx.AND.PX<=this%NX-2-this%maindx.AND.Py>=2+this%maindy.AND.PY<=this%Ny-2-this%maindy.AND.Pz>=2+this%maindz.AND.Pz<=this%Nz-2-this%maindz) then
               muEz=this%buffer%mainMu%Z(Px, Py, Pz);
               ejE=(1.0-muEz)*(t1*this%buffer%Ec1%Z(Px,Py,Pz)-t0*this%buffer%Ec0%Z(Px,Py,Pz))/this%ht;
               res = this%buffer%Je%Z(Px,Py,Pz)* this%thetabig(argthetabig) + ejE*cc/4.0d0/Pi;     
            else
              res = 0.0d0;  
            endif
         case default
            fatalerr = 209;
            write(*,*) 'Alert! Incorrect use of Component in getCurrentEPoint function!', component   
         end select
      case default
         fatalerr = 310;
         write(*,*) 'Alert! Incorrect Current Type is used in getCurrentEPointPML function!', this%currentstype
      end select
      getCurrentEPointPML = res;
    end function getCurrentEPointPML;


!getCurrentHPointPML------------------------------------------------------------------
    real*8 function getCurrentHPointPML(this, Sx, Sy, Sz, t, component)
    !Fill Jan sources with results from analytical solution
      class(tProblem) :: this
      integer, intent(in) :: Sx, Sy, Sz;
      integer :: Px, Py, Pz;
      integer, intent(in) :: t;
      integer, intent(in) :: component;
      real*8 :: res;
      real*8 :: argthetabig;
      res = 0;
      Px = Sx;
      Py = Sy;
      Pz = Sz;
      select case (this%pmlcurrentstype)
      case (0) ! Original sources don't have magnetic sources  
         select case (component)
         case(1) ! Jh X point
            res = 0;
         case(2) ! Jh Y point
            res = 0;
         case(3) ! Jh Z point
            res = 0;
         case default
            fatalerr = 311;
            write(*,*) 'Alert! Incorrect use of Component in getCurrentHPointPML function!', component   
         end select
      case (1) ! Original sources don't have magnetic sources    
         select case (component)
         case(1) ! Jh X point
            res = 0;
         case(2) ! Jh Y point
            res = 0;
         case(3) ! Jh Z point
            res = 0;
         case default
            fatalerr = 311;
            write(*,*) 'Alert! Incorrect use of Component in getCurrentHPointPML function!', component   
         end select
         argthetabig = this%Ti05(t)-(1.0+this%sgm_s)*this%auxmaxtime*this%problem_id;
         res = res * this%thetabig(argthetabig);
      case (2) ! Sources from buffer
         Px = Px - this%PML%NumPML;
         Py = Py - this%PML%NumPML;
         Pz = Pz - this%PML%NumPML;
         select case (component)
         case(1) ! Jh X point
            res = this%buffer%Jh%X(Px,Py,Pz);
         case(2) ! Jh Y point
            res = this%buffer%Jh%Y(Px,Py,Pz);
         case(3) ! Jh Z point
            res = this%buffer%Jh%Z(Px,Py,Pz);
         case default
            fatalerr = 211;
            write(*,*) 'Alert! Incorrect use of Component in getCurrentHPoint function!', component   
         end select
      case (3) ! Sources from buffer cutted with Theta function
         Px = Px - this%PML%NumPML;
         Py = Py - this%PML%NumPML;
         Pz = Pz - this%PML%NumPML;
         select case (component)
         case(1) ! Jh X point
            res = this%buffer%Jh%X(Px,Py,Pz);
         case(2) ! Jh Y point
            res = this%buffer%Jh%Y(Px,Py,Pz);
         case(3) ! Jh Z point
            res = this%buffer%Jh%Z(Px,Py,Pz);
         case default
            fatalerr = 211;
            write(*,*) 'Alert! Incorrect use of Component in getCurrentHPoint function!', component   
         end select
         argthetabig = this%Ti05(t)-(1.0d0+this%sgm_s)*this%auxmaxtime*this%problem_id;
         res = res * this%thetabig(argthetabig);
      case (4) ! Sources from buffer cutted with Theta function with Poisson correction (no difference as wihtout Poisson)
         Px = Px - this%PML%NumPML;
         Py = Py - this%PML%NumPML;
         Pz = Pz - this%PML%NumPML;
         argthetabig = this%Ti05(t)-(1.0+this%sgm_s)*this%auxmaxtime*this%problem_id;
         select case (component)
         case(1) ! Jh X point
            if (Px>=2+this%maindx.AND.PX<=this%NX-2-this%maindx.AND.Py>=2+this%maindy.AND.PY<=this%Ny-2-this%maindy.AND.Pz>=2+this%maindz.AND.Pz<=this%Nz-2-this%maindz) then
               res = this%buffer%Jh%X(Px,Py,Pz)*this%thetabig(argthetabig);
            else
               res = 0.0d0;
            endif
         case(2) ! Jh Y point
            if (Px>=2+this%maindx.AND.PX<=this%NX-2-this%maindx.AND.Py>=2+this%maindy.AND.PY<=this%Ny-2-this%maindy.AND.Pz>=2+this%maindz.AND.Pz<=this%Nz-2-this%maindz) then
               res = this%buffer%Jh%Y(Px,Py,Pz)*this%thetabig(argthetabig);
            else
               res = 0.0d0;
            endif   
         case(3) ! Jh Z point
            if (Px>=2+this%maindx.AND.PX<=this%NX-2-this%maindx.AND.Py>=2+this%maindy.AND.PY<=this%Ny-2-this%maindy.AND.Pz>=2+this%maindz.AND.Pz<=this%Nz-2-this%maindz) then
               res = this%buffer%Jh%Z(Px,Py,Pz)*this%thetabig(argthetabig);
            else
               res = 0.0d0;
            endif   
         case default
            fatalerr = 211;
            write(*,*) 'Alert! Incorrect use of Component in getCurrentHPoint function!', component   
         end select   
      case default
         fatalerr = 312;
         write(*,*) 'Alert! Incorrect Current Type is used in getCurrentHPointPML function!', this%currentstype
      end select   
      getCurrentHPointPML = res;
    end function getCurrentHPointPML;    
    

!getSourceCurrentsPML------------------------------------------------------------------
    subroutine getSourceCurrentsPML(this, t)
    !Upload J sources to PML
      class(tProblem) :: this
      integer :: i,j,k
      integer, intent(in) :: t;
      !JE currents
      do k=0,this%PML%Nz
         do j=0,this%PML%Ny
            do i=1,this%PML%Nx
               if ((i>=this%PML%NumPML+1).AND.(i<=this%PML%Nx-this%PML%NumPML).AND.(j>=this%PML%NumPML+1).AND.(j<=this%PML%Ny-this%PML%NumPML-1).AND.(k>=this%PML%NumPML+1).AND.(k<=this%PML%Nz-this%PML%NumPML-1)) then
                  this%PML%JfE%X(i,j,k)=this%getCurrentEPointPML(i,j,k,t,1); ! Jx
               endif   
            enddo
         enddo
      enddo
      do k=0,this%PML%Nz
         do j=1,this%PML%Ny
            do i=0,this%PML%Nx
               if ((i>=this%PML%NumPML+1).AND.(i<=this%PML%Nx-this%PML%NumPML-1).AND.(j>=this%PML%NumPML+1).AND.(j<=this%PML%Ny-this%PML%NumPML).AND.(k>=this%PML%NumPML+1).AND.(k<=this%PML%Nz-this%PML%NumPML-1)) then
                  this%PML%JfE%Y(i,j,k)=this%getCurrentEPointPML(i,j,k,t,2); ! Jy
               endif     
            enddo
         enddo
      enddo
      do k=1,this%PML%Nz
         do j=0,this%PML%Ny
            do i=0,this%PML%Nx
               if ((i>=this%PML%NumPML+1).AND.(i<=this%PML%Nx-this%PML%NumPML-1).AND.(j>=this%PML%NumPML+1).AND.(j<=this%PML%Ny-this%PML%NumPML-1).AND.(k>=this%PML%NumPML+1).AND.(k<=this%PML%Nz-this%PML%NumPML)) then
                  this%PML%JfE%Z(i,j,k)=this%getCurrentEPointPML(i,j,k,t,3); ! Jz
               endif   
            enddo
         enddo
      enddo
      !JH currents
      do k=1,this%PML%Nz
         do j=1,this%PML%Ny
            do i=0,this%PML%Nx
               if ((i>=this%PML%NumPML+1).AND.(i<=this%PML%Nx-this%PML%NumPML-1).AND.(j>=this%PML%NumPML+2).AND.(j<=this%PML%Ny-this%PML%NumPML-1).AND.(k>=this%PML%NumPML+2).AND.(k<=this%PML%Nz-this%PML%NumPML-1)) then
                  this%PML%JfH%X(i,j,k)=this%getCurrentHPointPML(i,j,k,t,1); ! Jx
               endif   
            enddo
         enddo
      enddo
      do k=1,this%PML%Nz
         do j=0,this%PML%Ny
            do i=1,this%PML%Nx
               if ((i>=this%PML%NumPML+2).AND.(i<=this%PML%Nx-this%PML%NumPML-1).AND.(j>=this%PML%NumPML+1).AND.(j<=this%PML%Ny-this%PML%NumPML-1).AND.(k>=this%PML%NumPML+2).AND.(k<=this%PML%Nz-this%PML%NumPML-1)) then
                  this%PML%JfH%Y(i,j,k)=this%getCurrentHPointPML(i,j,k,t,2); ! Jy
               endif     
            enddo
         enddo
      enddo
      do k=0,this%PML%Nz
         do j=1,this%PML%Ny
            do i=1,this%PML%Nx
               if ((i>=this%PML%NumPML+2).AND.(i<=this%PML%Nx-this%PML%NumPML-1).AND.(j>=this%PML%NumPML+2).AND.(j<=this%PML%Ny-this%PML%NumPML-1).AND.(k>=this%PML%NumPML+1).AND.(k<=this%PML%Nz-this%PML%NumPML-1)) then  
                  this%PML%JfH%Z(i,j,k)=this%getCurrentHPointPML(i,j,k,t,3); ! Jz
               endif   
            enddo
         enddo
      enddo
    end subroutine getSourceCurrentsPML;    

    
!getAnalyticSolution------------------------------------------------------------------
    subroutine getAnalyticSolution(this, t)
    ! Fill the analytical solution field components in the time point
      class(tProblem) :: this
      integer :: i,j,k
      integer, intent(in) :: t;
      ! E field filling
      !$OMP PARALLEL DO SCHEDULE(GUIDED)
      do k=0,this%Nz
         do j=0,this%Ny
            do i=1,this%Nx
               this%Ean%X(i,j,k)=this%Source%getEpoint(this%xi05(i),this%yi(j),this%zi(k), this%Ti(t), i, j, k, 1) ! Ex
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
      !$OMP PARALLEL DO SCHEDULE(GUIDED)
      do k=0,this%Nz
         do j=1,this%Ny
            do i=0,this%Nx
               this%Ean%Y(i,j,k)=this%Source%getEpoint(this%xi(i),this%yi05(j),this%zi(k), this%Ti(t), i, j, k, 2) ! Ey
             enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
      !$OMP PARALLEL DO SCHEDULE(GUIDED)
      do k=1,this%Nz
         do j=0,this%Ny
            do i=0,this%Nx
               this%Ean%Z(i,j,k)=this%Source%getEpoint(this%xi(i),this%yi(j),this%zi05(k), this%Ti(t), i, j, k, 3) ! Ey
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
      ! H field filling
      !$OMP PARALLEL DO SCHEDULE(GUIDED)
      do k=1,this%Nz
         do j=1,this%Ny
            do i=0,this%Nx
               this%Han%X(i,j,k)=this%Source%getHpoint(this%xi(i),this%yi05(j),this%zi05(k), this%Ti05(t+1), i, j, k, 1) ! Hx
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
      !$OMP PARALLEL DO SCHEDULE(GUIDED)
      do k=1,this%Nz
         do j=0,this%Ny
            do i=1,this%Nx
               this%Han%Y(i,j,k)=this%Source%getHpoint(this%xi05(i),this%yi(j),this%zi05(k), this%Ti05(t+1), i, j, k, 2) ! Hy
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
      !$OMP PARALLEL DO SCHEDULE(GUIDED)
      do k=0,this%Nz
         do j=1,this%Ny
            do i=1,this%Nx
               this%Han%Z(i,j,k)=this%Source%getHpoint(this%xi05(i),this%yi05(j),this%zi(k), this%Ti05(t+1), i, j, k, 3) ! Hz
             enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
    end subroutine getAnalyticSolution


!saveFields------------------------------------------------------------------
    subroutine saveFields(this, mode)
    ! Save actual fields to the _old fields
      class(tProblem) :: this
      integer, intent(in) :: mode 
      integer :: i,j,k
      select case (mode)
      case(0)
         ! E field filling
         do k=0,this%Nz
            do j=0,this%Ny
               do i=1,this%Nx
                  this%Efold%X(i,j,k)=this%Ef%X(i,j,k); ! Ex
               enddo
            enddo
         enddo
         do k=0,this%Nz
            do j=1,this%Ny
               do i=0,this%Nx
                  this%Efold%Y(i,j,k)=this%Ef%Y(i,j,k); ! Ey
               enddo
            enddo
         enddo
         do k=1,this%Nz
            do j=0,this%Ny
               do i=0,this%Nx
                  this%Efold%Z(i,j,k)=this%Ef%Z(i,j,k); ! Ey
               enddo
            enddo
         enddo
         ! H field filling
         do k=1,this%Nz
            do j=1,this%Ny
               do i=0,this%Nx
                  this%Hfold%X(i,j,k)=this%Hf%X(i,j,k); ! Hx
               enddo
            enddo
         enddo
         do k=1,this%Nz
            do j=0,this%Ny
               do i=1,this%Nx
                  this%Hfold%Y(i,j,k)=this%Hf%Y(i,j,k); ! Hy
               enddo
            enddo
         enddo
         do k=0,this%Nz
            do j=1,this%Ny
               do i=1,this%Nx
                  this%Hfold%Z(i,j,k)=this%Hf%Z(i,j,k); ! Hz
               enddo
            enddo
         enddo
      case(1)
         ! E field filling
         do k=0,this%Nz
            do j=0,this%Ny
               do i=1,this%Nx
                  this%Efold2%X(i,j,k)=this%Efold%X(i,j,k); ! Ex
               enddo
            enddo
         enddo
         do k=0,this%Nz
            do j=1,this%Ny
               do i=0,this%Nx
                  this%Efold2%Y(i,j,k)=this%Efold%Y(i,j,k); ! Ey
               enddo
            enddo
         enddo
         do k=1,this%Nz
            do j=0,this%Ny
               do i=0,this%Nx
                  this%Efold2%Z(i,j,k)=this%Efold%Z(i,j,k); ! Ey
               enddo
            enddo
         enddo
         ! H field filling
         do k=1,this%Nz
            do j=1,this%Ny
               do i=0,this%Nx
                  this%Hfold2%X(i,j,k)=this%Hfold%X(i,j,k); ! Hx
               enddo
            enddo
         enddo
         do k=1,this%Nz
            do j=0,this%Ny
               do i=1,this%Nx
                  this%Hfold2%Y(i,j,k)=this%Hfold%Y(i,j,k); ! Hy
               enddo
            enddo
         enddo
         do k=0,this%Nz
            do j=1,this%Ny
               do i=1,this%Nx
                  this%Hfold2%Z(i,j,k)=this%Hfold%Z(i,j,k); ! Hz
               enddo
            enddo
         enddo
      end select
    end subroutine saveFields


!getAnalyticSolutionPML---------------------------------------------------------------
    subroutine getAnalyticSolutionPML(this, t)
    ! Fill the analytical solution field components in the time point
      class(tProblem) :: this
      integer :: i,j,k
      integer, intent(in) :: t;
      ! E field filling
      do k=0,this%PML%Nz
         do j=0,this%PML%Ny
            do i=1,this%PML%Nx
               this%PML%Ean%X(i,j,k)=this%Source%getEpoint(this%PML%xi05(i),this%PML%yi(j),this%PML%zi(k), this%Ti(t), i, j, k, 1) ! Ex
            enddo
         enddo
      enddo
      do k=0,this%PML%Nz
         do j=1,this%PML%Ny
            do i=0,this%PML%Nx
               this%PML%Ean%Y(i,j,k)=this%Source%getEpoint(this%PML%xi(i),this%PML%yi05(j),this%PML%zi(k), this%Ti(t), i, j, k, 2) ! Ey
            enddo
         enddo
      enddo
      do k=1,this%PML%Nz
         do j=0,this%PML%Ny
            do i=0,this%PML%Nx
               this%PML%Ean%Z(i,j,k)=this%Source%getEpoint(this%PML%xi(i),this%PML%yi(j),this%PML%zi05(k), this%Ti(t), i, j, k, 3) ! Ey
            enddo
         enddo
      enddo

      ! H field filling
      do k=1,this%PML%Nz
         do j=1,this%PML%Ny
            do i=0,this%PML%Nx
               this%PML%Han%X(i,j,k)=this%Source%getHpoint(this%PML%xi(i),this%PML%yi05(j),this%PML%zi05(k), this%Ti05(t+1), i, j, k, 1) ! Hx
            enddo
         enddo
      enddo
      do k=1,this%PML%Nz
         do j=0,this%PML%Ny
            do i=1,this%PML%Nx
               this%PML%Han%Y(i,j,k)=this%Source%getHpoint(this%PML%xi05(i),this%PML%yi(j),this%PML%zi05(k), this%Ti05(t+1), i, j, k, 2) ! Hy
            enddo
         enddo
      enddo
      do k=0,this%PML%Nz
         do j=1,this%PML%Ny
            do i=1,this%PML%Nx
               this%PML%Han%Z(i,j,k)=this%Source%getHpoint(this%PML%xi05(i),this%PML%yi05(j),this%PML%zi(k), this%Ti05(t+1), i, j, k, 3) ! Hz
            enddo
         enddo
      enddo
    end subroutine getAnalyticSolutionPML    


!fillEBoundaryPoint------------------------------------------------------------------
  subroutine  fillEBoundaryPoint(this, t, Px, Py, Pz, component)
    ! Fill E point with some boundary condition
    class(tProblem) :: this
    integer, intent(in) :: t;
    integer, intent(in) :: Px, Py, Pz, component;
    select case (this%Eboundarytype)
    case (0) ! Zeros
       select case (component)
       case(1)
          this%Ef%X(Px, Py, Pz) = 0;
       case(2)
          this%Ef%Y(Px, Py, Pz) = 0;
       case(3)
          this%Ef%Z(Px, Py, Pz) = 0;
       case default ! illegal type
          fatalerr = 207;
          write(*,*) 'Alert! Incorrect use of Component in fillEboundaryPoint function!', component     
       end select
    case (1) ! Analytical solution
       select case (component)
       case(1)
          this%Ef%X(Px, Py, Pz) = this%Source%getEpoint(this%xi05(Px), this%yi(Py),  this%zi(Pz), this%Ti(t), Px, Py, Pz, component);
       case(2)
          this%Ef%Y(Px, Py, Pz) = this%Source%getEpoint(this%xi(Px), this%yi05(Py),  this%zi(Pz), this%Ti(t), Px, Py, Pz, component);
       case(3)
          this%Ef%Z(Px, Py, Pz) = this%Source%getEpoint(this%xi(Px), this%yi(Py),  this%zi05(Pz), this%Ti(t), Px, Py, Pz, component);
       case default ! illegal type
          fatalerr = 207;
          write(*,*) 'Alert! Incorrect use of Component in fillEboundaryPoint function!', component    
       end select
    case (2) ! Scheme
       select case (component)
       case(1)
          call this%updateEpoint(Px, Py, Pz, 1);
       case(2)
          call this%updateEpoint(Px, Py, Pz, 2);
       case(3)
          call this%updateEpoint(Px, Py, Pz, 3);
       case default ! illegal type
          fatalerr = 207;
          write(*,*) 'Alert! Incorrect use of Component in fillEboundaryPoint function!', component    
       end select
    case (3) ! PML
       select case (component)
       case(1)
          this%Ef%X(Px, Py, Pz)=this%PML%Ef%X(Px+this%PML%NumPml, Py+this%PML%NumPml, Pz+this%PML%NumPml);
       case(2)
          this%Ef%Y(Px, Py, Pz)=this%PML%Ef%Y(Px+this%PML%NumPml, Py+this%PML%NumPml, Pz+this%PML%NumPml);
       case(3)
          this%Ef%Z(Px, Py, Pz)=this%PML%Ef%Z(Px+this%PML%NumPml, Py+this%PML%NumPml, Pz+this%PML%NumPml);
       case default ! illegal type
          fatalerr = 207;
          write(*,*) 'Alert! Incorrect use of Component in fillEboundaryPoint function!', component    
       end select
    case (4) ! Buffer
       select case (component)
       case(1)
          this%Ef%X(Px, Py, Pz)=this%buffer%Ef%X(Px+this%bufoffsetx, Py+this%bufoffsety, Pz+this%bufoffsetz);
       case(2)
          this%Ef%Y(Px, Py, Pz)=this%buffer%Ef%Y(Px+this%bufoffsetx, Py+this%bufoffsety, Pz+this%bufoffsetz);
       case(3)
          this%Ef%Z(Px, Py, Pz)=this%buffer%Ef%Z(Px+this%bufoffsetx, Py+this%bufoffsety, Pz+this%bufoffsetz);
       case default ! illegal type
          fatalerr = 207;
          write(*,*) 'Alert! Incorrect use of Component in fillEboundaryPoint function!', component    
       end select
    case (5) ! Buffer and Static
       select case (component)
       case(1)
          this%Ef%X(Px, Py, Pz)=this%buffer%Ef%X(Px+this%bufoffsetx, Py+this%bufoffsety, Pz+this%bufoffsetz)+this%buffer%Estatic%X(Px+this%bufoffsetx, Py+this%bufoffsety, Pz+this%bufoffsetz);
       case(2)
          this%Ef%Y(Px, Py, Pz)=this%buffer%Ef%Y(Px+this%bufoffsetx, Py+this%bufoffsety, Pz+this%bufoffsetz)+this%buffer%Estatic%Y(Px+this%bufoffsetx, Py+this%bufoffsety, Pz+this%bufoffsetz);
       case(3)
          this%Ef%Z(Px, Py, Pz)=this%buffer%Ef%Z(Px+this%bufoffsetx, Py+this%bufoffsety, Pz+this%bufoffsetz)+this%buffer%Estatic%Z(Px+this%bufoffsetx, Py+this%bufoffsety, Pz+this%bufoffsetz);
       case default ! illegal type
          fatalerr = 207;
          write(*,*) 'Alert! Incorrect use of Component in fillEboundaryPoint function!', component    
       end select    
    case default ! illegal type
       fatalerr = 205;
       write(*,*) 'Alert! Illegal boundary condition for E in fillEboundaryPoint function!', this%Eboundarytype   
    end select
    !if (checkupdated==1) then
    !   if (this%Eboundarytype /= 2) then
    !      select case (component)
    !      case(1)
    !         this%Ef%chX(Px, Py, Pz) = this%Ef%chX(Px, Py, Pz)+1;
    !      case(2)
    !         this%Ef%chY(Px, Py, Pz) = this%Ef%chY(Px, Py, Pz)+1;
    !      case(3)
    !         this%Ef%chZ(Px, Py, Pz) = this%Ef%chZ(Px, Py, Pz)+1;
    !      end select
    !   endif
    !endif  
  end subroutine fillEBoundaryPoint


!fillHBoundaryPoint------------------------------------------------------------------
  subroutine  fillHBoundaryPoint(this, t, Px, Py, Pz, component)
    ! Fill H point with some boundary condition
    class(tProblem) :: this
    integer, intent(in) :: t;
    integer, intent(in) :: Px, Py, Pz, component;
    select case (this%Hboundarytype)
    case (0) ! Zeros
       select case (component)
       case(1)
          this%Hf%X(Px, Py, Pz) = 0;
       case(2)
          this%Hf%Y(Px, Py, Pz) = 0;
       case(3)
          this%Hf%Z(Px, Py, Pz) = 0;
       case default ! illegal type
          fatalerr = 208;
          write(*,*) 'Alert! Incorrect use of Component in fillHboundaryPoint function!', component    
       end select
    case (1) ! Analytical solution
       select case (component)
       case(1)
          this%Hf%X(Px, Py, Pz) = this%Source%getHpoint(this%xi(Px), this%yi05(Py),  this%zi05(Pz), this%Ti05(t+1), Px, Py, Pz, component);      
       case(2)
          this%Hf%Y(Px, Py, Pz) = this%Source%getHpoint(this%xi05(Px), this%yi(Py),  this%zi05(Pz), this%Ti05(t+1), Px, Py, Pz, component);
       case(3)
          this%Hf%Z(Px, Py, Pz) = this%Source%getHpoint(this%xi05(Px), this%yi05(Py),  this%zi(Pz), this%Ti05(t+1), Px, Py, Pz, component);
       case default ! illegal type
          fatalerr = 208;
          write(*,*) 'Alert! Incorrect use of Component in fillHboundaryPoint function!', component       
       end select
    case (2) ! Scheme
       select case (component)
       case(1)
          call this%updateHpoint(Px, Py, Pz, 1);
       case(2)
          call this%updateHpoint(Px, Py, Pz, 2);
       case(3)
          call this%updateHpoint(Px, Py, Pz, 3);
       case default ! illegal type
          fatalerr = 208;
          write(*,*) 'Alert! Incorrect use of Component in fillHboundaryPoint function!', component       
       end select
    case (3) ! PML
       select case (component)
       case(1)
          this%Hf%X(Px, Py, Pz)=this%PML%Hf%X(Px+this%PML%NumPml, Py+this%PML%NumPml, Pz+this%PML%NumPml);
       case(2)
          this%Hf%Y(Px, Py, Pz)=this%PML%Hf%Y(Px+this%PML%NumPml, Py+this%PML%NumPml, Pz+this%PML%NumPml);
       case(3)
          this%Hf%Z(Px, Py, Pz)=this%PML%Hf%Z(Px+this%PML%NumPml, Py+this%PML%NumPml, Pz+this%PML%NumPml);
       case default ! illegal type
          fatalerr = 208;
          write(*,*) 'Alert! Incorrect use of Component in fillHboundaryPoint function!', component    
       end select
    case (4) ! Buffer
       select case (component)
       case(1)
          this%Hf%X(Px, Py, Pz)=this%buffer%Hf%X(Px+this%bufoffsetx, Py+this%bufoffsety, Pz+this%bufoffsetz);
       case(2)
          this%Hf%Y(Px, Py, Pz)=this%buffer%Hf%Y(Px+this%bufoffsetx, Py+this%bufoffsety, Pz+this%bufoffsetz);
       case(3)
          this%Hf%Z(Px, Py, Pz)=this%buffer%Hf%Z(Px+this%bufoffsetx, Py+this%bufoffsety, Pz+this%bufoffsetz);
       case default ! illegal type
          fatalerr = 208;
          write(*,*) 'Alert! Incorrect use of Component in fillHboundaryPoint function!', component    
       end select    
    case default ! illegal type
       fatalerr = 206;
       write(*,*) 'Alert! Illegal boundary condition for H in fillEboundaryPoint function!', this%Hboundarytype
    end select
    !if (checkupdated==1) then
    !   if (this%Hboundarytype /= 2) then  
    !      select case (component)
    !      case(1)
    !         this%Hf%chX(Px, Py, Pz) = this%Hf%chX(Px, Py, Pz)+1;
    !      case(2)
    !         this%Hf%chY(Px, Py, Pz) = this%Hf%chY(Px, Py, Pz)+1;
    !      case(3)
    !         this%Hf%chZ(Px, Py, Pz) = this%Hf%chZ(Px, Py, Pz)+1;
    !      end select
    !   endif
    !endif
  end subroutine fillHBoundaryPoint

  
!fillEBoundary------------------------------------------------------------------
  subroutine  fillEBoundary(this,t)
    !Fill E boundary with some boundary condition
    class(tProblem) :: this
    integer, intent(in) :: t
    integer :: i ,j ,k
    if (this%EBoundaryType>-1) then
       ! Fill EX Boundaries
       !$OMP PARALLEL DO SCHEDULE(GUIDED)
       do k=0,this%Nz
          do i=1,this%Nx
             call this%fillEBoundaryPoint(t, i, 0, k, 1);
             call this%fillEBoundaryPoint(t, i, this%Ny, k, 1);
          enddo
       enddo
       !$OMP END PARALLEL DO
       !$OMP PARALLEL DO SCHEDULE(GUIDED)
       do j=1,this%Ny-1
          do i=1,this%Nx
             call this%fillEBoundaryPoint(t, i, j, 0, 1); 
             call this%fillEBoundaryPoint(t, i, j, this%Nz, 1);
          enddo
       enddo
      !$OMP END PARALLEL DO
       ! Fill EY Boundaries
       !$OMP PARALLEL DO SCHEDULE(GUIDED)
       do j=1,this%Ny
          do i=0,this%Nx
             call this%fillEBoundaryPoint(t, i, j, 0, 2);
             call this%fillEBoundaryPoint(t, i, j, this%Nz, 2);
          enddo
       enddo
       !$OMP END PARALLEL DO
       !$OMP PARALLEL DO SCHEDULE(GUIDED)
       do k=1,this%Nz-1
          do j=1,this%Ny
             call this%fillEBoundaryPoint(t, 0, j, k, 2);
             call this%fillEBoundaryPoint(t, this%Nx, j, k, 2);
          enddo
       enddo
       !$OMP END PARALLEL DO
       ! Fill EZ Boundaries
       !$OMP PARALLEL DO SCHEDULE(GUIDED)
       do k=1,this%Nz
          do j=0,this%Ny
             call this%fillEBoundaryPoint(t, 0, j, k, 3);
             call this%fillEBoundaryPoint(t, this%Nx, j, k, 3);
          enddo
       enddo
       !$OMP END PARALLEL DO
       !$OMP PARALLEL DO SCHEDULE(GUIDED)
       do k=1,this%Nz
          do i=1,this%Nx-1
             call this%fillEBoundaryPoint(t, i, 0, k, 3);
             call this%fillEBoundaryPoint(t, i, this%Ny, k, 3);
          enddo
       enddo
       !$OMP END PARALLEL DO
    endif
  end subroutine fillEBoundary
 

!fillHBoundary------------------------------------------------------------------
  subroutine  fillHBoundary(this, t)
    ! Fill H boundary with some boundary condition
    class(tProblem) :: this
    integer, intent(in) :: t
    integer :: i ,j ,k
    if (this%HBoundaryType>-1) then
       ! HX Boundaries Update
       !$OMP PARALLEL DO SCHEDULE(GUIDED)
       do j=1,this%Ny
          do i=0,this%Nx
             call this%fillHBoundaryPoint(t, i, j, 1, 1);
             call this%fillHBoundaryPoint(t, i, j, this%Nz, 1);
          enddo
       enddo
       !$OMP END PARALLEL DO
       !$OMP PARALLEL DO SCHEDULE(GUIDED)
       do k=2,this%Nz-1
          do i=0,this%Nx
             call this%fillHBoundaryPoint(t, i, 1, k, 1);
             call this%fillHBoundaryPoint(t, i, this%Ny, k, 1);
          enddo
       enddo
       !$OMP END PARALLEL DO
       !$OMP PARALLEL DO SCHEDULE(GUIDED)
       do k=2,this%Nz-1
          do j=2,this%Ny-1
             call this%fillHBoundaryPoint(t, 0, j, k, 1);
             call this%fillHBoundaryPoint(t, this%Nx, j, k, 1);
          enddo
       enddo
       !$OMP END PARALLEL DO
       !$OMP PARALLEL DO SCHEDULE(GUIDED)
       ! HY Boundaries Update
       do k=1,this%Nz
          do j=0,this%Ny
             call this%fillHBoundaryPoint(t, 1, j, k, 2);
             call this%fillHBoundaryPoint(t, this%Nx, j, k, 2);
          enddo
       enddo
       !$OMP END PARALLEL DO
       !$OMP PARALLEL DO SCHEDULE(GUIDED)
       do j=0,this%Ny
          do i=2,this%Nx-1
             call this%fillHBoundaryPoint(t, i, j, 1, 2);
             call this%fillHBoundaryPoint(t, i, j, this%Nz, 2);
          enddo
       enddo
       !$OMP END PARALLEL DO
       !$OMP PARALLEL DO SCHEDULE(GUIDED)
       do k=2,this%Nz-1
          do i=2,this%Nx-1
             call this%fillHBoundaryPoint(t, i, 0, k, 2);
             call this%fillHBoundaryPoint(t, i, this%Ny, k, 2);
          enddo
       enddo
       !$OMP END PARALLEL DO
       ! HZ Boundaries Update
       !$OMP PARALLEL DO SCHEDULE(GUIDED)
       do k=0,this%Nz
          do i=1,this%Nx
             call this%fillHBoundaryPoint(t, i, 1, k, 3);
             call this%fillHBoundaryPoint(t, i, this%Ny, k, 3);
          enddo
       enddo
       !$OMP END PARALLEL DO
       !$OMP PARALLEL DO SCHEDULE(GUIDED)
       do k=0,this%Nz
          do j=2,this%Ny-1
             call this%fillHBoundaryPoint(t, 1, j, k, 3);
             call this%fillHBoundaryPoint(t, this%Nx, j, k, 3);
          enddo
       enddo
       !$OMP END PARALLEL DO
       !$OMP PARALLEL DO SCHEDULE(GUIDED)
       do j=2,this%Ny-1
          do i=2,this%Nx-1
             call this%fillHBoundaryPoint(t, i, j, 0, 3);
             call this%fillHBoundaryPoint(t, i, j, this%Nz, 3);
          enddo
       enddo
       !$OMP END PARALLEL DO
    endif
  end subroutine fillHBoundary


!uploadFieldsToPML------------------------------------------------------------------
    subroutine uploadFieldsToPML(this)
    ! Transfer all the fields without boundaries to PML
      class(tProblem) :: this
      integer :: i,j,k
      ! E field transfering
      !$OMP PARALLEL DO SCHEDULE(GUIDED)
      do k=1,this%Nz-1
         do j=1,this%Ny-1
            do i=1,this%Nx
               this%PML%Ef%X(i+this%PML%NumPml,j+this%PML%NumPml,k+this%PML%NumPml) = this%Ef%X(i,j,k); ! Ex
               !if (checkupdated==1) then
               !   this%PML%Ef%chX(i+this%PML%NumPml,j+this%PML%NumPml,k+this%PML%NumPml) = this%PML%Ef%chX(i+this%PML%NumPml,j+this%PML%NumPml,k+this%PML%NumPml)+1;   ! this string is important for update checking
               !endif
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
      !$OMP PARALLEL DO SCHEDULE(GUIDED)
      do k=1,this%Nz-1
         do j=1,this%Ny
            do i=1,this%Nx-1
               this%PML%Ef%Y(i+this%PML%NumPml,j+this%PML%NumPml,k+this%PML%NumPml) = this%Ef%Y(i,j,k); ! Ey
               !if (checkupdated==1) then
               !   this%PML%Ef%chY(i+this%PML%NumPml,j+this%PML%NumPml,k+this%PML%NumPml) = this%PML%Ef%chY(i+this%PML%NumPml,j+this%PML%NumPml,k+this%PML%NumPml)+1;   ! this string is important for update checking
               !endif
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
      !$OMP PARALLEL DO SCHEDULE(GUIDED)
      do k=1,this%Nz
         do j=1,this%Ny-1
            do i=1,this%Nx-1
               this%PML%Ef%Z(i+this%PML%NumPml,j+this%PML%NumPml,k+this%PML%NumPml) = this%Ef%Z(i,j,k); ! Ez
               !if (checkupdated==1) then
               !   this%PML%Ef%chZ(i+this%PML%NumPml,j+this%PML%NumPml,k+this%PML%NumPml) = this%PML%Ef%chZ(i+this%PML%NumPml,j+this%PML%NumPml,k+this%PML%NumPml)+1;   ! this string is important for update checking
               !endif
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
      
      ! H field transfering
      !$OMP PARALLEL DO SCHEDULE(GUIDED)
      do k=2,this%Nz-1
         do j=2,this%Ny-1
            do i=1,this%Nx-1
               this%PML%Hf%X(i+this%PML%NumPml,j+this%PML%NumPml,k+this%PML%NumPml) = this%Hf%X(i,j,k); ! Hx
               !if (checkupdated==1) then
               !   this%PML%Hf%chX(i+this%PML%NumPml,j+this%PML%NumPml,k+this%PML%NumPml) = this%PML%Hf%chX(i+this%PML%NumPml,j+this%PML%NumPml,k+this%PML%NumPml)+1;   ! this string is important for update checking
               !endif
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
      !$OMP PARALLEL DO SCHEDULE(GUIDED)
      do k=2,this%Nz-1
         do j=1,this%Ny-1
            do i=2,this%Nx-1
               this%PML%Hf%Y(i+this%PML%NumPml,j+this%PML%NumPml,k+this%PML%NumPml) = this%Hf%Y(i,j,k); ! Hy
               !if (checkupdated==1) then
               !   this%PML%Hf%chY(i+this%PML%NumPml,j+this%PML%NumPml,k+this%PML%NumPml) = this%PML%Hf%chY(i+this%PML%NumPml,j+this%PML%NumPml,k+this%PML%NumPml)+1;   ! this string is important for update checking
               !endif
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
      !$OMP PARALLEL DO SCHEDULE(GUIDED)
      do k=1,this%Nz-1
         do j=2,this%Ny-1
            do i=2,this%Nx-1
               this%PML%Hf%Z(i+this%PML%NumPml,j+this%PML%NumPml,k+this%PML%NumPml) = this%Hf%Z(i,j,k); ! Hz
               !if (checkupdated==1) then
               !   this%PML%Hf%chZ(i+this%PML%NumPml,j+this%PML%NumPml,k+this%PML%NumPml) = this%PML%Hf%chZ(i+this%PML%NumPml,j+this%PML%NumPml,k+this%PML%NumPml)+1;   ! this string is important for update checking
               !endif
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
    end subroutine uploadFieldsToPML 


!getCurrentError------------------------------------------------------------------
  subroutine getCurrentError(this, t)
  ! Calculate the maximum relative and asolute errors in the current time-point by each field component 
    use commonvars
    class(tProblem) :: this
    integer, intent(in) :: t;
    integer i, j, k
    real*8 :: error;
    error = 0
    this%err = 0
    this%Ex_err%absval = 0.0d0
    this%Ey_err%absval = 0.0d0
    this%Ez_err%absval = 0.0d0
    this%Hx_err%absval = 0.0d0
    this%Hy_err%absval = 0.0d0
    this%Hz_err%absval = 0.0d0
    this%Ex_err%fval = 0.0d0
    this%Ey_err%fval = 0.0d0
    this%Ez_err%fval = 0.0d0
    this%Hx_err%fval = 0.0d0
    this%Hy_err%fval = 0.0d0
    this%Hz_err%fval = 0.0d0
    ! H field 
    do k=1,this%Nz
       do j=1,this%Ny
          do i=0,this%Nx
             if (abs(this%Hf%X(i,j,k)-this%Han%X(i,j,k)) > this%Hx_err%absval) then
                this%Hx_err%absval = abs(this%Hf%X(i,j,k)-this%Han%X(i,j,k));
                this%Hx_err%x=i
                this%Hx_err%y=j
                this%Hx_err%z=k
             endif
             if (abs(this%Hf%X(i,j,k)) > this%Hx_err%fval) then
                this%Hx_err%fval = abs(this%Hf%X(i,j,k)); 
             endif
          enddo
       enddo
    enddo
    do k=1,this%Nz
       do j=0,this%Ny
          do i=1,this%Nx
             if (abs(this%Hf%Y(i,j,k)-this%Han%Y(i,j,k)) > this%Hy_err%absval) then
                this%Hy_err%absval = abs(this%Hf%Y(i,j,k)-this%Han%Y(i,j,k));
                this%Hy_err%x=i
                this%Hy_err%y=j
                this%Hy_err%z=k
             endif
             if (abs(this%Hf%Y(i,j,k)) > this%Hy_err%fval) then
                this%Hy_err%fval = abs(this%Hf%Y(i,j,k)); 
             endif
          enddo
       enddo
    enddo
    do k=0,this%Nz
       do j=1,this%Ny
          do i=1,this%Nx
             if (abs(this%Hf%Z(i,j,k)-this%Han%Z(i,j,k)) > this%Hz_err%absval) then
                this%Hz_err%absval = abs(this%Hf%Z(i,j,k)-this%Han%Z(i,j,k));
                this%Hz_err%x=i
                this%Hz_err%y=j
                this%Hz_err%z=k
             endif
             if (abs(this%Hf%Z(i,j,k)) > this%Hz_err%fval) then
                this%Hz_err%fval = abs(this%Hf%Z(i,j,k)); 
             endif
          enddo
       enddo
    enddo

    ! E field 
    do k=0,this%Nz
       do j=0,this%Ny
          do i=1,this%Nx
             if (abs(this%Ef%X(i,j,k)-this%Ean%X(i,j,k)) > this%Ex_err%absval) then
                this%Ex_err%absval = abs(this%Ef%X(i,j,k)-this%Ean%X(i,j,k));
                this%Ex_err%x=i
                this%Ex_err%y=j
                this%Ex_err%z=k
             endif
             if (abs(this%Ef%X(i,j,k)) > this%Ex_err%fval) then
                this%Ex_err%fval = abs(this%Ef%X(i,j,k)); 
             endif
          enddo
       enddo
    enddo
    do k=0,this%Nz
       do j=1,this%Ny
          do i=0,this%Nx
             if (abs(this%Ef%Y(i,j,k)-this%Ean%Y(i,j,k)) > this%Ey_err%absval) then
                this%Ey_err%absval = abs(this%Ef%Y(i,j,k)-this%Ean%Y(i,j,k));
                this%Ey_err%x=i
                this%Ey_err%y=j
                this%Ey_err%z=k
             endif
             if (abs(this%Ef%Y(i,j,k)) > this%Ey_err%fval) then
                this%Ey_err%fval = abs(this%Ef%Y(i,j,k)); 
             endif
          enddo
       enddo
    enddo
    do k=1,this%Nz
       do j=0,this%Ny
          do i=0,this%Nx
             if (abs(this%Ef%Z(i,j,k)-this%Ean%Z(i,j,k)) > this%Ez_err%absval) then
                this%Ez_err%absval = abs(this%Ef%Z(i,j,k)-this%Ean%Z(i,j,k));
                this%Ez_err%x=i
                this%Ez_err%y=j
                this%Ez_err%z=k
             endif
             if (abs(this%Ef%Z(i,j,k)) > this%Ez_err%fval) then
                this%Ez_err%fval = abs(this%Ef%Z(i,j,k)); 
             endif
          enddo
       enddo
    enddo
    if (this%Hx_err%fval /= 0.0d0.AND.this%Hx_err%absval/=this%Hx_err%fval) then
       this%Hx_err%val = this%Hx_err%absval / this%Hx_err%fval;
    else
       this%Hx_err%val = -1;
    endif
    if (this%Hy_err%fval /= 0.0d0.AND.this%Hy_err%absval/=this%Hy_err%fval) then
       this%Hy_err%val = this%Hy_err%absval / this%Hy_err%fval;
    else
       this%Hy_err%val = -1;
    endif
    if (this%Hz_err%fval /= 0.0d0.AND.this%Hz_err%absval/=this%Hz_err%fval) then
       this%Hz_err%val = this%Hz_err%absval / this%Hz_err%fval;
    else
       this%Hz_err%val = -1;
    endif
    if (this%Ex_err%fval /= 0.0d0.AND.this%Ex_err%absval/=this%Ex_err%fval) then
       this%Ex_err%val = this%Ex_err%absval / this%Ex_err%fval;
    else
       this%Ex_err%val = -1;
    endif
    if (this%Ey_err%fval /= 0.0d0.AND.this%Ey_err%absval/=this%Ey_err%fval) then
       this%Ey_err%val = this%Ey_err%absval / this%Ey_err%fval;
    else
       this%Ey_err%val = -1;
    endif
    if (this%Ez_err%fval /= 0.0d0.AND.this%Ez_err%absval/=this%Ez_err%fval) then
       this%Ez_err%val = this%Ez_err%absval / this%Ez_err%fval;
    else
       this%Ez_err%val = -1;
    endif     
    this%err = max(this%Ex_err%val, this%Ey_err%val, this%Ez_err%val, this%Hx_err%val, this%Hy_err%val, this%Hz_err%val);
    this%abserr = max(this%Ex_err%absval, this%Ey_err%absval, this%Ez_err%absval, this%Hx_err%absval, this%Hy_err%absval, this%Hz_err%absval);

    if (t>this%Nt*0.1) then
       if (this%Ex_err%val > this%Ex_maxerr%val) then
          call this%SaveToMaxError(0);
       endif
       if (this%Ey_err%val > this%Ey_maxerr%val) then
          call this%SaveToMaxError(1);
       endif
       if (this%Ez_err%val > this%Ez_maxerr%val) then
          call this%SaveToMaxError(2);
       endif
       if (this%Hx_err%val > this%Hx_maxerr%val) then
          call this%SaveToMaxError(3);
       endif
       if (this%Hy_err%val > this%Hy_maxerr%val) then
          call this%SaveToMaxError(4);
       endif
       if (this%Hz_err%val > this%Hz_maxerr%val) then
          call this%SaveToMaxError(5);
       endif
    endif   
  end subroutine getCurrentError


!SaveToMaxError
  subroutine SaveToMaxError(this, component)
    class(tProblem) :: this
    integer, intent(in) :: component;
    select case (component)
    case (0)
       this%Ex_maxerr%val = this%Ex_err%val;
       this%Ex_maxerr%absval = this%Ex_err%absval;
       this%Ex_maxerr%fval = this%Ex_err%fval;
       this%Ex_maxerr%X = this%Ex_err%X;
       this%Ex_maxerr%Y = this%Ex_err%Y;
       this%Ex_maxerr%Z = this%Ex_err%Z;
    case (1)
       this%Ey_maxerr%val = this%Ey_err%val;
       this%Ey_maxerr%absval = this%Ey_err%absval;
       this%Ey_maxerr%fval = this%Ey_err%fval;
       this%Ey_maxerr%X = this%Ey_err%X;
       this%Ey_maxerr%Y = this%Ey_err%Y;
       this%Ey_maxerr%Z = this%Ey_err%Z;
    case (2)
       this%Ez_maxerr%val = this%Ez_err%val;
       this%Ez_maxerr%absval = this%Ez_err%absval;
       this%Ez_maxerr%fval = this%Ez_err%fval;
       this%Ez_maxerr%X = this%Ez_err%X;
       this%Ez_maxerr%Y = this%Ez_err%Y;
       this%Ez_maxerr%Z = this%Ez_err%Z;
    case (3)
       this%Hx_maxerr%val = this%Hx_err%val;
       this%Hx_maxerr%absval = this%Hx_err%absval;
       this%Hx_maxerr%fval = this%Hx_err%fval;
       this%Hx_maxerr%X = this%Hx_err%X;
       this%Hx_maxerr%Y = this%Hx_err%Y;
       this%Hx_maxerr%Z = this%Hx_err%Z;
    case (4)
       this%Hy_maxerr%val = this%Hy_err%val;
       this%Hy_maxerr%absval = this%Hy_err%absval;
       this%Hy_maxerr%fval = this%Hy_err%fval;
       this%Hy_maxerr%X = this%Hy_err%X;
       this%Hy_maxerr%Y = this%Hy_err%Y;
       this%Hy_maxerr%Z = this%Hy_err%Z;
    case (5)
       this%Hz_maxerr%val = this%Hz_err%val;
       this%Hz_maxerr%absval = this%Hz_err%absval;
       this%Hz_maxerr%fval = this%Hz_err%fval;
       this%Hz_maxerr%X = this%Hz_err%X;
       this%Hz_maxerr%Y = this%Hz_err%Y;
       this%Hz_maxerr%Z = this%Hz_err%Z;   
    end select
  end subroutine SaveToMaxError  

!getCurrentErrorPML------------------------------------------------------------------
   subroutine getCurrentErrorPML(this, t)
   ! Calculate the maximum Error of the PML in the current time-point by each field component
    use commonvars
    class(tProblem) :: this
    integer, intent(in) :: t;
    integer i, j, k
    real*8 :: error;
    error = 0.0d0;
    this%PML%err = 0.0d0;
    this%PML%Ex_err%absval = 0.0d0;
    this%PML%Ey_err%absval = 0.0d0;
    this%PML%Ez_err%absval = 0.0d0;
    this%PML%Hx_err%absval = 0.0d0;
    this%PML%Hy_err%absval = 0.0d0;
    this%PML%Hz_err%absval = 0.0d0;
    this%PML%Ex_err%fval = 0.0d0;
    this%PML%Ey_err%fval = 0.0d0;
    this%PML%Ez_err%fval = 0.0d0;
    this%PML%Hx_err%fval = 0.0d0;
    this%PML%Hy_err%fval = 0.0d0;
    this%PML%Hz_err%fval = 0.0d0;
    ! H field filling
    do k=1+(this%PML%NumPml),this%PML%Nz-(this%PML%NumPml)
       do j=1+(this%PML%NumPml),this%PML%Ny-(this%PML%NumPml)
          do i=0+(this%PML%NumPml),this%PML%Nx-(this%PML%NumPml)
             if (abs(this%PML%Hf%X(i,j,k)-this%PML%Han%X(i,j,k)) > this%PML%Hx_err%absval) then
                this%PML%Hx_err%absval = abs(this%PML%Hf%X(i,j,k)-this%PML%Han%X(i,j,k));
                this%PML%Hx_err%x=i
                this%PML%Hx_err%y=j
                this%PML%Hx_err%z=k
             endif
             if (abs(this%PML%Hf%X(i,j,k)) > this%PML%Hx_err%fval) then
                this%PML%Hx_err%fval = abs(this%PML%Hf%X(i,j,k)); 
             endif
          enddo
       enddo
    enddo
    do k=1+(this%PML%NumPml),this%PML%Nz-(this%PML%NumPml)
       do j=0+(this%PML%NumPml),this%PML%Ny-(this%PML%NumPml)
          do i=1+(this%PML%NumPml),this%PML%Nx-(this%PML%NumPml)
             if (abs(this%PML%Hf%Y(i,j,k)-this%PML%Han%Y(i,j,k)) > this%PML%Hy_err%absval) then
                this%PML%Hy_err%absval = abs(this%PML%Hf%Y(i,j,k)-this%PML%Han%Y(i,j,k));
                this%PML%Hy_err%x=i
                this%PML%Hy_err%y=j
                this%PML%Hy_err%z=k
             endif
             if (abs(this%PML%Hf%Y(i,j,k)) > this%PML%Hy_err%fval) then
                this%PML%Hy_err%fval = abs(this%PML%Hf%Y(i,j,k)); 
             endif
          enddo
       enddo
    enddo
    do k=0+(this%PML%NumPml),this%PML%Nz-(this%PML%NumPml)
       do j=1+(this%PML%NumPml),this%PML%Ny-(this%PML%NumPml)
          do i=1+(this%PML%NumPml),this%PML%Nx-(this%PML%NumPml)
             if (abs(this%PML%Hf%Z(i,j,k)-this%PML%Han%Z(i,j,k)) > this%PML%Hz_err%absval) then
                this%PML%Hz_err%absval = abs(this%PML%Hf%Z(i,j,k)-this%PML%Han%Z(i,j,k));
                this%PML%Hz_err%x=i
                this%PML%Hz_err%y=j
                this%PML%Hz_err%z=k
             endif
             if (abs(this%PML%Hf%Z(i,j,k)) > this%PML%Hz_err%fval) then
                this%PML%Hz_err%fval = abs(this%PML%Hf%Z(i,j,k)); 
             endif
          enddo
       enddo
    enddo

    ! E field filling
    do k=0+(this%PML%NumPml),this%PML%Nz-(this%PML%NumPml)
       do j=0+(this%PML%NumPml),this%PML%Ny-(this%PML%NumPml)
          do i=1+(this%PML%NumPml),this%PML%Nx-(this%PML%NumPml)
             if (abs(this%PML%Ef%X(i,j,k)-this%PML%Ean%X(i,j,k)) > this%PML%Ex_err%absval) then
                this%PML%Ex_err%absval = abs(this%PML%Ef%X(i,j,k)-this%PML%Ean%X(i,j,k));
                this%PML%Ex_err%x=i
                this%PML%Ex_err%y=j
                this%PML%Ex_err%z=k
             endif
             if (abs(this%PML%Ef%X(i,j,k)) > this%PML%Ex_err%fval) then
                this%PML%Ex_err%fval = abs(this%PML%Ef%X(i,j,k)); 
             endif
          enddo
       enddo
    enddo
    do k=0+(this%PML%NumPml),this%PML%Nz-(this%PML%NumPml)
       do j=1+(this%PML%NumPml),this%PML%Ny-(this%PML%NumPml)
          do i=0+(this%PML%NumPml),this%PML%Nx-(this%PML%NumPml)
             if (abs(this%PML%Ef%Y(i,j,k)-this%PML%Ean%Y(i,j,k)) > this%PML%Ey_err%absval) then
                this%PML%Ey_err%absval = abs(this%PML%Ef%Y(i,j,k)-this%PML%Ean%Y(i,j,k));
                this%PML%Ey_err%x=i
                this%PML%Ey_err%y=j
                this%PML%Ey_err%z=k
             endif
             if (abs(this%PML%Ef%Y(i,j,k)) > this%PML%Ey_err%fval) then
                this%PML%Ey_err%fval = abs(this%PML%Ef%Y(i,j,k)); 
             endif
          enddo
       enddo
    enddo
    do k=1+(this%PML%NumPml),this%PML%Nz-(this%PML%NumPml)
       do j=0+(this%PML%NumPml),this%PML%Ny-(this%PML%NumPml)
          do i=0+(this%PML%NumPml),this%PML%Nx-(this%PML%NumPml)
             if (abs(this%PML%Ef%Z(i,j,k)-this%PML%Ean%Z(i,j,k)) > this%PML%Ez_err%absval) then
                this%PML%Ez_err%absval = abs(this%PML%Ef%Z(i,j,k)-this%PML%Ean%Z(i,j,k));
                this%PML%Ez_err%x=i
                this%PML%Ez_err%y=j
                this%PML%Ez_err%z=k
             endif
             if (abs(this%PML%Ef%Z(i,j,k)) > this%PML%Ez_err%fval) then
                this%PML%Ez_err%fval = abs(this%PML%Ef%Z(i,j,k)); 
             endif
          enddo
       enddo
    enddo
    if (this%PML%Hx_err%fval /= 0) then
       this%PML%Hx_err%val = this%PML%Hx_err%absval / this%PML%Hx_err%fval;
    else
       this%PML%Hx_err%val = -1;
    endif
    if (this%PML%Hy_err%fval /= 0) then
       this%PML%Hy_err%val = this%PML%Hy_err%absval / this%PML%Hy_err%fval;
    else
       this%PML%Hy_err%val = -1;
    endif
    if (this%PML%Hz_err%fval /= 0) then
       this%PML%Hz_err%val = this%PML%Hz_err%absval / this%PML%Hz_err%fval;
    else
       this%PML%Hz_err%val = -1;
    endif
    if (this%PML%Ex_err%fval /= 0) then
       this%PML%Ex_err%val = this%PML%Ex_err%absval / this%PML%Ex_err%fval;
    else
       this%PML%Ex_err%val = -1;
    endif
    if (this%PML%Ey_err%fval /= 0) then
       this%PML%Ey_err%val = this%PML%Ey_err%absval / this%PML%Ey_err%fval;
    else
       this%PML%Ey_err%val = -1;
    endif
    if (this%PML%Ez_err%fval /= 0) then
       this%PML%Ez_err%val = this%PML%Ez_err%absval / this%PML%Ez_err%fval;
    else
       this%PML%Ez_err%val = -1;
    endif     
    this%PML%err = max(this%PML%Ex_err%val, this%PML%Ey_err%val, this%PML%Ez_err%val, this%PML%Hx_err%val, this%PML%Hy_err%val, this%PML%Hz_err%val);
    this%PML%abserr = max(this%PML%Ex_err%absval, this%PML%Ey_err%absval, this%PML%Ez_err%absval, this%PML%Hx_err%absval, this%PML%Hy_err%absval, this%PML%Hz_err%absval);
  end subroutine getCurrentErrorPML


!P4------------------------------------------------------------------------------------------  
  real*8 function P4(this, x)
    ! Transition function for hat function and Mu functions
    class(tProblem) :: this
    real*8 :: x;	
    if (x <= 0d0) then
       P4=0.0d0
    elseif ((x>=0.0d0).AND.(x<=1.0d0)) then
       P4=1.0d0-(1.0d0-35.0d0*x**4+84.0d0*x**5-70.0d0*x**6+20.0d0*x**7)
    else
       P4=1.0d0
    endif
  end function P4

  
!P4inv----------------------------------------------------------------------------------------
  real*8 function P4inv(this, x)
    ! Inverse transition function for hat function and Mu functions
    class(tProblem) :: this
    real*8 :: x;	
    if (x<=0d0) then
       P4inv=1.0d0
    elseif ((x>=0.0d0).AND.(x<=1.0d0)) then
       P4inv=(1.0d0-35.0d0*x**4+84.0d0*x**5-70.0d0*x**6+20.0d0*x**7)
    else
       P4inv=0.0d0
    endif
  end function P4inv
  
  
!ThetaBig------------------------------------------------------------------------------------
  real*8 function thetabig(this, x)
    ! Hat function for auxuliary problems
    class(tProblem) :: this
    real*8 :: x;
    if (abs(x)<=this%sgm_s*this%auxmaxtime) then
	thetabig=1.0d0
    elseif (abs(x)>=this%auxmaxtime) then
	thetabig=0.0d0
    else
	thetabig=this%P4inv((abs(x)-this%sgm_s*this%auxmaxtime)/(this%auxmaxtime*(1.0d0-this%sgm_s)))
    endif
  end function thetabig

  
!Muu-------------------------------------------------------------------------------------------
  real*8 function muu(this, xlc, ylc, zlc)
    !Mu function for generation of effective currents
    class(tProblem) :: this
    real*8, intent(in) :: xlc, ylc, zlc;
    real*8 :: rxy,rxyz;
    !Planes---------
      !X plane
    	if (xlc>=this%Xi(this%Nx-this%Nlc).AND.ylc>this%Yi(this%Nlc).AND.ylc<this%Yi(this%Ny-this%Nlc).AND.zlc>this%Zi(this%Nlc).AND.zlc<this%Zi(this%Nz-this%Nlc)) then
	   muu=this%P4((xlc-this%Xi(this%Nx-this%Nlc))/this%Dlc);
        elseif (xlc<=this%Xi(this%Nlc).AND.ylc>this%Yi(this%Nlc).AND.ylc<this%Yi(this%Ny-this%Nlc).AND.zlc>this%Zi(this%Nlc).AND.zlc<this%Zi(this%Nz-this%Nlc)) then
           muu=this%P4((this%Xi(this%Nlc)-xlc)/this%Dlc)
      !Y plane
	elseif (ylc>=this%Yi(this%Ny-this%Nlc).AND.xlc>this%Xi(this%Nlc).AND.xlc<this%Xi(this%Nx-this%Nlc).AND.zlc>this%Zi(this%Nlc).AND.zlc<this%Zi(this%Nz-this%Nlc)) then
	   muu=this%P4((ylc-this%Yi(this%Ny-this%Nlc))/this%Dlc)
        elseif (ylc<=this%Yi(this%Nlc).AND.xlc>this%Xi(this%Nlc).AND.xlc<this%Xi(this%Nx-this%Nlc).AND.zlc>this%Zi(this%Nlc).AND.zlc<this%Zi(this%Nz-this%Nlc)) then
	   muu=this%P4((this%Yi(this%Nlc)-ylc)/this%Dlc)
      !Z plane
	elseif (zlc>=this%Zi(this%Nz-this%Nlc).AND.xlc>this%Xi(this%Nlc).AND.xlc<this%Xi(this%Nx-this%Nlc).AND.ylc>this%Yi(this%Nlc).AND.ylc<this%Yi(this%Ny-this%Nlc)) then
           muu=this%P4((zlc-this%Zi(this%Nz-this%Nlc))/this%Dlc)
        elseif (zlc<=this%Zi(this%Nlc).AND.xlc>this%Xi(this%Nlc).AND.xlc<this%Xi(this%Nx-this%Nlc).AND.ylc>this%Yi(this%Nlc).AND.ylc<this%Yi(this%Ny-this%Nlc)) then
           muu=this%P4((this%Zi(this%Nlc)-zlc)/this%Dlc)
    !Edges---------        
      !Vertical edges
        elseif (xlc>=this%Xi(this%Nx-this%Nlc).AND.ylc>=this%Yi(this%Ny-this%Nlc).AND.zlc>this%Zi(this%Nlc).AND.zlc<this%Zi(this%Nz-this%Nlc)) then
           rxy=sqrt((xlc-this%Xi(this%Nx-this%Nlc))**2+(ylc-this%Yi(this%Ny-this%Nlc))**2);
           muu=this%P4(rxy/this%Dlc);
        elseif (xlc<=this%Xi(this%Nlc).AND.ylc>=this%Yi(this%Ny-this%Nlc).AND.zlc>this%Zi(this%Nlc).AND.zlc<this%Zi(this%Nz-this%Nlc)) then
           rxy=sqrt((xlc-this%Xi(this%Nlc))**2+(ylc-this%Yi(this%Ny-this%Nlc))**2);
           muu=this%P4(rxy/this%Dlc);
        elseif (xlc>=this%Xi(this%Nx-this%Nlc).AND.ylc<=this%Yi(this%Nlc).AND.zlc>this%Zi(this%Nlc).AND.zlc<this%Zi(this%Nz-this%Nlc)) then
           rxy=sqrt((xlc-this%Xi(this%Nx-this%Nlc))**2+(ylc-this%Yi(this%Nlc))**2); 
           muu=this%P4(rxy/this%Dlc);
        elseif (xlc<=this%Xi(this%Nlc).AND.ylc<=this%Yi(this%Nlc).AND.zlc>this%Zi(this%Nlc).AND.zlc<this%Zi(this%Nz-this%Nlc)) then
           rxy=sqrt((xlc-this%Xi(this%Nlc))**2+(ylc-this%Yi(this%Nlc))**2);
           muu=this%P4(rxy/this%Dlc);
      !Horizontal edges, upper plane      
	elseif (xlc>=this%Xi(this%Nx-this%Nlc).AND.ylc>this%Yi(this%Nlc).AND.ylc<this%Yi(this%Ny-this%Nlc).AND.zlc>=this%Zi(this%Nz-this%Nlc)) then
           rxy=sqrt((xlc-this%Xi(this%Nx-this%Nlc))**2+(zlc-this%Zi(this%Nz-this%Nlc))**2) 
           muu=this%P4(rxy/this%Dlc);
        elseif (xlc<=this%Xi(this%Nlc).AND.ylc>this%Yi(this%Nlc).AND.ylc<this%Yi(this%Ny-this%Nlc).AND.zlc>=this%Zi(this%Nz-this%Nlc)) then
           rxy=sqrt((xlc-this%Xi(this%Nlc))**2+(zlc-this%Zi(this%Nz-this%Nlc))**2); 
           muu=this%P4(rxy/this%Dlc);
	elseif (xlc>this%Xi(this%Nlc).AND.xlc<this%Xi(this%Nx-this%Nlc).AND.ylc>=this%Yi(this%Ny-this%Nlc).AND.zlc>=this%Zi(this%Nz-this%Nlc)) then
           rxy=sqrt((ylc-this%Yi(this%Ny-this%Nlc))**2+(zlc-this%Zi(this%Nz-this%Nlc))**2); 
           muu=this%P4(rxy/this%Dlc);
	elseif (xlc>this%Xi(this%Nlc).AND.xlc<this%Xi(this%Nx-this%Nlc).AND.ylc<=this%Yi(this%nlc).AND.zlc>=this%Zi(this%Nz-this%Nlc)) then
           rxy=sqrt((ylc-this%Yi(this%Nlc))**2+(zlc-this%Zi(this%Nz-this%Nlc))**2);
           muu=this%P4(rxy/this%Dlc);
      !Horizontal edges, lower plane      
	elseif (xlc>=this%Xi(this%Nx-this%Nlc).AND.ylc>this%Yi(this%Nlc).AND.ylc<this%Yi(this%Ny-this%Nlc).AND.zlc<=this%Zi(this%Nlc)) then
           rxy=sqrt((xlc-this%Xi(this%Nx-this%Nlc))**2+(zlc-this%Zi(this%Nlc))**2);
           muu=this%P4(rxy/this%Dlc)
        elseif (xlc<=this%Xi(this%Nlc).AND.ylc>this%Yi(this%Nlc).AND.ylc<this%Yi(this%Ny-this%Nlc).AND.zlc<=this%Zi(this%Nlc)) then
           rxy=sqrt((xlc-this%Xi(this%Nlc))**2+(zlc-this%Zi(this%Nlc))**2); 
           muu=this%P4(rxy/this%Dlc);
	elseif (ylc>=this%yi(this%Ny-this%Nlc).AND.xlc>this%Xi(this%Nlc).AND.xlc<this%Xi(this%Nx-this%Nlc).AND.zlc<=this%Zi(this%Nz-this%Nlc)) then
           rxy=sqrt((ylc-this%Yi(this%Ny-this%Nlc))**2+(zlc-this%Zi(this%Nlc))**2) 
           muu=this%P4(rxy/this%Dlc);
	elseif (ylc<=this%Yi(this%Nlc).AND.xlc>this%Xi(this%Nlc).AND.xlc<this%Xi(this%Nx-this%Nlc).AND.zlc<=this%Zi(this%Nlc)) then
           rxy=sqrt((ylc-this%Yi(this%Nlc))**2+(zlc-this%Zi(this%Nlc))**2) 
           muu=this%P4(rxy/this%Dlc);
    !Corners--------
      !Upper plane 
	elseif (xlc>=this%Xi(this%Nx-this%Nlc).AND.ylc>=this%Yi(this%Ny-this%Nlc).AND.zlc>=this%Zi(this%Nz-this%Nlc)) then
           rxyz=sqrt((xlc-this%Xi(this%Nx-this%Nlc))**2+(ylc-this%Yi(this%Ny-this%Nlc))**2+(zlc-this%Zi(this%Nz-this%Nlc))**2);
           muu=this%P4(rxyz/this%Dlc);
        elseif (xlc<=this%Xi(this%Nlc).AND.ylc>=this%Yi(this%Ny-this%Nlc).AND.zlc>=this%Zi(this%Nz-this%Nlc)) then
           rxyz=sqrt((xlc-this%Xi(this%Nlc))**2+(ylc-this%Yi(this%Ny-this%Nlc))**2+(zlc-this%Zi(this%Nz-this%Nlc))**2);
           muu=this%P4(rxyz/this%Dlc);
        elseif (xlc>=this%Xi(this%Nx-this%Nlc).AND.ylc<=this%Yi(this%Nlc).AND.zlc>=this%Zi(this%Nz-this%Nlc)) then
	   rxyz=sqrt((xlc-this%Xi(this%Nx-this%Nlc))**2+(ylc-this%Yi(this%Nlc))**2+(zlc-this%Zi(this%Nz-this%Nlc))**2);
           muu=this%P4(rxyz/this%Dlc);
        elseif (xlc<=this%Xi(this%Nlc).AND.ylc<=this%Yi(this%Nlc).AND.zlc>=this%Zi(this%Nz-this%Nlc)) then
           rxyz=sqrt((xlc-this%Xi(this%Nlc))**2+(ylc-this%Yi(this%Nlc))**2+(zlc-this%Zi(this%Nz-this%Nlc))**2);
           muu=this%P4(rxyz/this%Dlc);
      !Lower plane
 	elseif (xlc>=this%Xi(this%Nx-this%Nlc).AND.ylc>=this%Yi(this%Ny-this%Nlc).AND.zlc<=this%Zi(this%Nlc)) then
           rxyz=sqrt((xlc-this%Xi(this%Nx-this%Nlc))**2+(ylc-this%Yi(this%Ny-this%Nlc))**2+(zlc-this%Zi(this%Nlc))**2);
           muu=this%P4(rxyz/this%Dlc);
        elseif (xlc<=this%Xi(this%Nlc).AND.ylc>=this%Yi(this%Ny-this%Nlc).AND.zlc<=this%Zi(this%Nlc)) then
	   rxyz=sqrt((xlc-this%Xi(this%Nlc))**2+(ylc-this%Yi(this%Ny-this%Nlc))**2+(zlc-this%Zi(this%Nlc))**2);
           muu=this%P4(rxyz/this%Dlc);
        elseif (xlc>=this%Xi(this%Nx-this%Nlc).AND.ylc<=this%Yi(this%Nlc).AND.zlc<=this%Zi(this%Nlc)) then
           rxyz=sqrt((xlc-this%Xi(this%Nx-this%Nlc))**2+(ylc-this%Yi(this%Nlc))**2+(zlc-this%Zi(this%Nlc))**2);
           muu=this%P4(rxyz/this%Dlc);
        elseif (xlc<=this%Xi(this%Nlc).AND.ylc<=this%Yi(this%Nlc).AND.zlc<=this%Zi(this%Nlc)) then
           rxyz=sqrt((xlc-this%Xi(this%Nlc))**2+(ylc-this%Yi(this%Nlc))**2+(zlc-this%Zi(this%Nlc))**2);
           muu=this%P4(rxyz/this%Dlc);
	else
   	   muu=0.0
    endif
   ! muu=1;
    end function muu

!InitMuArray-------------------------------------------------------------------------------------------
    subroutine InitMuArray(this)
      !Initialize mup array to store all the values of mu function for effective currents
      class(tProblem) :: this
      integer :: i,j,k
      do k=0,this%Nz+1
         do j=0,this%Ny+1
            do i=0,this%Nx+1
               this%mup(2*i,2*j,2*k)=this%muu(this%xi05(i), this%yi05(j), this%zi05(k));
               if (i<this%NX+1) then
                  this%mup(2*i+1,2*j,2*k)=this%muu(this%xi(i), this%yi05(j), this%zi05(k));
               endif
               if (j<this%NY+1) then
                  this%mup(2*i,2*j+1,2*k)=this%muu(this%xi05(i), this%yi(j), this%zi05(k));
               endif
               if (k<this%NZ+1) then
                  this%mup(2*i,2*j,2*k+1)=this%muu(this%xi05(i), this%yi05(j), this%zi(k));
               endif

               if (i<this%NX+1.AND.j<this%NY+1) then
                  this%mup(2*i+1,2*j+1,2*k)=this%muu(this%xi(i), this%yi(j), this%zi05(k));
               endif
               if (i<this%NX+1.AND.k<this%NZ+1) then
                  this%mup(2*i+1,2*j,2*k+1)=this%muu(this%xi(i), this%yi05(j), this%zi(k));
               endif
               if (j<this%NY+1.AND.k<this%NZ+1) then
                  this%mup(2*i,2*j+1,2*k+1)=this%muu(this%xi05(i), this%yi(j), this%zi(k));
               endif
               if (i<this%NX+1.AND.j<this%NY+1.AND.k<this%NZ+1) then
                  this%mup(2*i+1,2*j+1,2*k+1)=this%muu(this%xi(i), this%yi(j), this%zi(k));
               endif
            enddo
         enddo
      enddo  
    end subroutine InitMuArray  
  
!Problem_PrintCurrentErrors------------------------------------------------------------------      
    subroutine problem_PrintCurrentErrors(this, t)
    !Print current errors to the screen
      class(tProblem) :: this
      integer, intent(in) :: t;
      if ((mod(t,screen_erroroutputsteps)==0)) then
         write (*,*) 'errorHx=', this%Hx_err%absval, this%Hx_err%val
         write (*,*) 'errorHy=', this%Hy_err%absval, this%Hy_err%val
         write (*,*) 'errorHz=', this%Hz_err%absval, this%Hz_err%val
         write (*,*) 'errorEx=', this%Ex_err%absval, this%Ex_err%val
         write (*,*) 'errorEy=', this%Ey_err%absval, this%Ey_err%val
         write (*,*) 'errorEz=', this%Ez_err%absval, this%Ez_err%val
      endif
    end subroutine problem_PrintCurrentErrors;


!Problem_PrintMaximumAnalyticalFields------------------------------------------------------------------      
    subroutine problem_PrintMaximumAnalyticalFields(this, t)
    !Print current errors to the screen
      class(tProblem) :: this
      integer, intent(in) :: t;
      integer :: i,j,k;
      real*8 :: maxEx, maxEy, maxEz, maxHx, maxHy, maxHz;
      maxEx = 0.0d0;
      maxEy = 0.0d0;
      maxEz = 0.0d0;
      maxHx = 0.0d0;
      maxHy = 0.0d0;
      maxHz = 0.0d0;
      
      ! E field 
    do k=0,this%Nz
       do j=0,this%Ny
          do i=1,this%Nx
             if (this%Ean%X(i,j,k)>maxEx) then
                maxEx = this%Ean%X(i,j,k);
             endif
          enddo
       enddo
    enddo
    do k=0,this%Nz
       do j=1,this%Ny
          do i=0,this%Nx
             if (this%Ean%Y(i,j,k)>maxEy) then
                maxEy = this%Ean%Y(i,j,k);
             endif
          enddo
       enddo
    enddo
    do k=1,this%Nz
       do j=0,this%Ny
          do i=0,this%Nx
             if (this%Ean%Z(i,j,k)>maxEz) then
                maxEz = this%Ean%Z(i,j,k);
             endif
          enddo
       enddo
    enddo
    
      ! H field 
    do k=1,this%Nz
       do j=1,this%Ny
          do i=0,this%Nx
             if (this%Han%X(i,j,k)>maxHx) then
                maxHx = this%Han%X(i,j,k);
             endif
          enddo
       enddo
    enddo
    do k=1,this%Nz
       do j=0,this%Ny
          do i=1,this%Nx
             if (this%Han%Y(i,j,k)>maxHy) then
                maxHy = this%Han%Y(i,j,k);
             endif
          enddo
       enddo
    enddo
    do k=0,this%Nz
       do j=1,this%Ny
          do i=1,this%Nx
             if (this%Han%Z(i,j,k)>maxHz) then
                maxHz = this%Han%Z(i,j,k);
             endif
          enddo
       enddo
    enddo

    !write (*,*) 'x: ', this%xi(14);
    !write (*,*) 'y: ', this%yi(13);
    !write (*,*) 'z: ', this%zi05(12);
    !write (*,*) 't: ', this%ti(120);
    !write (*,*) 'Ez(14,13,12,120)=', this%Source%getEpoint(this%xi(14),this%yi(13),this%zi05(12), this%Ti(120), 3) ! Ey
    
    write (*,*) 'maxEx: ', maxEx;
    write (*,*) 'maxEy: ', maxEy;
    write (*,*) 'maxEz: ', maxEz;
    write (*,*) 'maxHx: ', maxHx;
    write (*,*) 'maxHy: ', maxHy;
    write (*,*) 'maxHz: ', maxHz;
  end subroutine problem_PrintMaximumAnalyticalFields;


!Problem_PrintMaximumCurrents------------------------------------------------------------------      
    subroutine problem_PrintMaximumCurrents(this, t)
    !Print current errors to the screen
      class(tProblem) :: this
      integer, intent(in) :: t;
      integer :: i,j,k;
      real*8 :: maxJx, maxJy, maxJz;
      maxJx = 0.0d0;
      maxJy = 0.0d0;
      maxJz = 0.0d0;
      
      ! E field 
    do k=0,this%Nz
       do j=0,this%Ny
          do i=1,this%Nx
             if (this%Je%X(i,j,k)>maxJx) then
                maxJx = this%Je%X(i,j,k);
             endif
          enddo
       enddo
    enddo
    do k=0,this%Nz
       do j=1,this%Ny
          do i=0,this%Nx
            if (this%Je%Y(i,j,k)>maxJy) then
                maxJy = this%Je%Y(i,j,k);
             endif
          enddo
       enddo
    enddo
    do k=1,this%Nz
       do j=0,this%Ny
          do i=0,this%Nx
             if (this%Je%Z(i,j,k)>maxJz) then
                maxJz = this%Je%Z(i,j,k);
             endif
          enddo
       enddo
    enddo

    !write (*,*) 'x: ', this%xi(14);
    !write (*,*) 'y: ', this%yi(13);
    !write (*,*) 'z: ', this%zi05(12);
    !write (*,*) 't: ', this%ti(120);
    !write (*,*) 'Ez(14,13,12,120)=', this%Source%getEpoint(this%xi(14),this%yi(13),this%zi05(12), this%Ti(120), 3) ! Ey
    
    write (*,*) 'maxJx: ', maxJx;
    write (*,*) 'maxJy: ', maxJy;
    write (*,*) 'maxJz: ', maxJz;
  end subroutine problem_PrintMaximumCurrents;  


!Problem_PrintMaximumFields------------------------------------------------------------------      
    subroutine problem_PrintMaximumFields(this, t)
    !Print current errors to the screen
      class(tProblem) :: this
      integer, intent(in) :: t;
      integer :: i,j,k;
      real*8 :: maxEx, maxEy, maxEz, maxHx, maxHy, maxHz;
      maxEx = 0.0d0;
      maxEy = 0.0d0;
      maxEz = 0.0d0;
      maxHx = 0.0d0;
      maxHy = 0.0d0;
      maxHz = 0.0d0;
      
    ! E field 
    do k=0,this%Nz
       do j=0,this%Ny
          do i=1,this%Nx
             if (this%Ef%X(i,j,k)>maxEx) then
                maxEx = this%Ef%X(i,j,k);
             endif
          enddo
       enddo
    enddo
    do k=0,this%Nz
       do j=1,this%Ny
          do i=0,this%Nx
             if (this%Ef%Y(i,j,k)>maxEy) then
                maxEy = this%Ef%Y(i,j,k);
             endif
          enddo
       enddo
    enddo
    do k=1,this%Nz
       do j=0,this%Ny
          do i=0,this%Nx
             if (this%Ef%Z(i,j,k)>maxEz) then
                maxEz = this%Ef%Z(i,j,k);
             endif
          enddo
       enddo
    enddo
    
    ! H field 
    do k=1,this%Nz
       do j=1,this%Ny
          do i=0,this%Nx
             if (this%Hf%X(i,j,k)>maxHx) then
                maxHx = this%Hf%X(i,j,k);
             endif
          enddo
       enddo
    enddo
    do k=1,this%Nz
       do j=0,this%Ny
          do i=1,this%Nx
             if (this%Hf%Y(i,j,k)>maxHy) then
                maxHy = this%Hf%Y(i,j,k);
             endif
          enddo
       enddo
    enddo
    do k=0,this%Nz
       do j=1,this%Ny
          do i=1,this%Nx
             if (this%Hf%Z(i,j,k)>maxHz) then
                maxHz = this%Hf%Z(i,j,k);
             endif
          enddo
       enddo
    enddo

    !write (*,*) 'x: ', this%xi(14);
    !write (*,*) 'y: ', this%yi(13);
    !write (*,*) 'z: ', this%zi05(12);
    !write (*,*) 't: ', this%ti(120);
    !write (*,*) 'Ez(14,13,12,120)=', this%Source%getEpoint(this%xi(14),this%yi(13),this%zi05(12), this%Ti(120), 3) ! Ey
    
    write (*,*) 'maxEx: ', maxEx;
    write (*,*) 'maxEy: ', maxEy;
    write (*,*) 'maxEz: ', maxEz;
    write (*,*) 'maxHx: ', maxHx;
    write (*,*) 'maxHy: ', maxHy;
    write (*,*) 'maxHz: ', maxHz;
    end subroutine problem_PrintMaximumFields;  


!Problem_PrintPMLCurrentErrors------------------------------------------------------------------      
    subroutine problem_PrintPMLCurrentErrors(this, t)
    !Print current errors to the screen
      class(tProblem) :: this
      integer, intent(in) :: t;
      if ((mod(t,screen_erroroutputsteps)==0)) then
         write (*,*) 'PMLerrorHx=', this%PML%Hx_err%absval, this%PML%Hx_err%val
         write (*,*) 'PMLerrorHy=', this%PML%Hy_err%absval, this%PML%Hy_err%val
         write (*,*) 'PMLerrorHz=', this%PML%Hz_err%absval, this%PML%Hz_err%val
         write (*,*) 'PMLerrorEx=', this%PML%Ex_err%absval, this%PML%Ex_err%val
         write (*,*) 'PMLerrorEy=', this%PML%Ey_err%absval, this%PML%Ey_err%val
         write (*,*) 'PMLerrorEz=', this%PML%Ez_err%absval, this%PML%Ez_err%val
      endif
    end subroutine problem_PrintPMLCurrentErrors;    


!addEtobuffer------------------------------------------------------------------
    subroutine addEtobuffer(this)
    ! Add E values to buffer
      class(tProblem) :: this
      integer :: i,j,k
      do k=0,this%Nz
         do j=0,this%Ny
            do i=1,this%Nx
               this%buffer%Ef%X(i,j,k)=this%buffer%Ef%X(i,j,k)+this%Ef%X(i,j,k); ! Ex
            enddo
         enddo
      enddo
      do k=0,this%Nz
         do j=1,this%Ny
            do i=0,this%Nx
               this%buffer%Ef%Y(i,j,k)=this%buffer%Ef%Y(i,j,k)+this%Ef%Y(i,j,k); ! Ey
            enddo
         enddo
      enddo
      do k=1,this%Nz
         do j=0,this%Ny
            do i=0,this%Nx
               this%buffer%Ef%Z(i,j,k)=this%buffer%Ef%Z(i,j,k)+this%Ef%Z(i,j,k); ! Ez
            enddo
         enddo
      enddo
    end subroutine addEtobuffer 


!addHtobuffer------------------------------------------------------------------
    subroutine addHtobuffer(this)
      class(tProblem) :: this
      integer :: i,j,k
      do k=1,this%Nz
         do j=1,this%Ny
            do i=0,this%Nx
               this%buffer%Hf%X(i,j,k)=this%buffer%Hf%X(i,j,k)+this%Hf%X(i,j,k); ! Hx
            enddo
         enddo
      enddo
      do k=1,this%Nz
         do j=0,this%Ny
            do i=1,this%Nx
               this%buffer%Hf%Y(i,j,k)=this%buffer%Hf%Y(i,j,k)+this%Hf%Y(i,j,k); ! Hy
            enddo
         enddo
      enddo
      do k=0,this%Nz
         do j=1,this%Ny
            do i=1,this%Nx
               this%buffer%Hf%Z(i,j,k)=this%buffer%Hf%Z(i,j,k)+this%Hf%Z(i,j,k); ! Hz
             enddo
         enddo
      enddo
    end subroutine addHtobuffer    


!addEstaticToBuffer------------------------------------------------------------------
    subroutine addEstaticToBuffer(this)
    ! Add E values to buffer
      class(tProblem) :: this
      integer :: i,j,k
      do k=0,this%Nz
         do j=0,this%Ny
            do i=1,this%Nx
               this%buffer%Estatic%X(i,j,k)=this%buffer%Estatic%X(i,j,k)+this%Ef%X(i,j,k); ! Ex
             enddo
         enddo
      enddo
      do k=0,this%Nz
         do j=1,this%Ny
            do i=0,this%Nx
               this%buffer%Estatic%Y(i,j,k)=this%buffer%Estatic%Y(i,j,k)+this%Ef%Y(i,j,k); ! Ey
             enddo
         enddo
      enddo
      do k=1,this%Nz
         do j=0,this%Ny
            do i=0,this%Nx
               this%buffer%Estatic%Z(i,j,k)=this%buffer%Estatic%Z(i,j,k)+this%Ef%Z(i,j,k); ! Ez
            enddo
         enddo
      enddo
    end subroutine addEstaticToBuffer
    
    
!JEeffbuildPoint------------------------------------------------------------------
  subroutine  JEeffbuildPoint(this, t, Px, Py, Pz, component)
    ! Build JE effective current at the point
    class(tProblem) :: this
    integer, intent(in) :: t;
    integer, intent(in) :: Px, Py, Pz, component;
    real*8 :: muEx, muEy, muEz;
    real*8 :: muHxr, muHxl, muHyr, muHyl, muHzr, muHzl;
    real*8 :: DH;
    select case (component)
    case(1) !JEeff_X
      muEx=  this%mup(2*Px,2*Py+1,2*Pz+1)     !this%muu(this%xi05(Px), this%yi(Py), this%zi(Pz));
      muHyr= this%mup(2*Px,2*Py+1,2*(Pz+1))   !this%muu(this%xi05(Px), this%yi(Py), this%zi05(Pz+1));
      muHyl= this%mup(2*Px,2*Py+1,2*Pz)       !this%muu(this%xi05(Px), this%yi(Py), this%zi05(Pz));
      muHzr= this%mup(2*Px,2*(Py+1),2*Pz+1)   !this%muu(this%xi05(Px), this%yi05(Py+1), this%zi(Pz));
      muHzl= this%mup(2*Px,2*Py,2*Pz+1)       !this%muu(this%xi05(Px), this%yi05(Py), this%zi(Pz));
      DH = (muHzr*this%Hf%Z(Px,Py+1,Pz)-muHzl*this%Hf%Z(Px,Py,Pz))/this%hhy-(muHyr*this%Hf%Y(Px,Py,Pz+1)-muHyl*this%Hf%Y(Px,Py,Pz))/this%hhz;
      this%JEeff%X(Px,Py,Pz) = ( DH - muEx*(this%Ef%X(Px,Py,Pz) - this%Efold%X(Px,Py,Pz))/cc/this%ht )*cc/4.0d0/Pi;
    case(2) !JEeff_Y
      muEy=  this%mup(2*Px+1,2*Py,2*Pz+1)     !this%muu(this%xi(Px),this%yi05(Py),this%zi(Pz));
      muHxr= this%mup(2*Px+1,2*Py,2*(Pz+1))   !this%muu(this%xi(Px),this%yi05(Py),this%zi05(Pz+1));
      muHxl= this%mup(2*Px+1,2*Py,2*Pz)       !this%muu(this%xi(Px),this%yi05(Py),this%zi05(Pz));
      muHzr= this%mup(2*(Px+1),2*Py,2*Pz+1)   !this%muu(this%xi05(Px+1),this%yi05(Py),this%zi(Pz));
      muHzl= this%mup(2*Px,2*Py,2*Pz+1)       !this%muu(this%xi05(Px),this%yi05(Py),this%zi(Pz));
      DH = (muHxr*this%Hf%X(Px,Py,Pz+1)-muHxl*this%Hf%X(Px,Py,Pz))/this%hhz-(muHzr*this%Hf%Z(Px+1,Py,Pz)-muHzl*this%Hf%Z(Px,Py,Pz))/this%hhx;
      this%JEeff%Y(Px,Py,Pz) = ( DH - muEy*(this%Ef%Y(Px,Py,Pz) - this%Efold%Y(Px,Py,Pz))/cc/this%ht )*cc/4.0d0/Pi;
    case(3) !JEeff_Z
      muEz= this%mup(2*Px+1,2*Py+1,2*Pz)      !this%muu(this%xi(Px),this%yi(Py),this%zi05(Pz));
      muHyr=this%mup(2*(Px+1),2*Py+1,2*Pz)    !this%muu(this%xi05(Px+1),this%yi(Py),this%zi05(Pz));
      muHyl=this%mup(2*Px,2*Py+1,2*Pz)        !this%muu(this%xi05(Px),this%yi(Py),this%zi05(Pz));
      muHxr=this%mup(2*Px+1,2*(Py+1),2*Pz)    !this%muu(this%xi(Px),this%yi05(Py+1),this%zi05(Pz));
      muHxl=this%mup(2*Px+1,2*Py,2*Pz)        !this%muu(this%xi(Px),this%yi05(Py),this%zi05(Pz));
      DH = (muHyr*this%Hf%Y(Px+1,Py,Pz)-muHyl*this%Hf%Y(Px,Py,Pz))/this%hhx-(muHxr*this%Hf%X(Px,Py+1,Pz)-muHxl*this%Hf%X(Px,Py,Pz))/this%hhy;
      this%JEeff%Z(Px,Py,Pz) = ( DH - muEz*(this%Ef%Z(Px,Py,Pz) - this%Efold%Z(Px,Py,Pz))/cc/this%ht )*cc/4.0d0/Pi;
    case default ! illegal type
       fatalerr = 213;
       write(*,*) 'Alert! Incorrect use of Component in JEbuildPoint function!', component     
    end select
  end subroutine JEeffbuildPoint


!JHeffbuildPoint------------------------------------------------------------------
  subroutine  JHeffbuildPoint(this, t, Px, Py, Pz, component)
    ! Build JE effective current at the point
    class(tProblem) :: this
    integer, intent(in) :: t;
    integer, intent(in) :: Px, Py, Pz, component;
    real*8 :: muExr, muExl, muEyr, muEyl, muEzr, muEzl;
    real*8 :: muHx, muHy, muHz;
    real*8 :: DE;
    select case (component)
    case(1) !JHeff_X
      muHx=  this%mup(2*Px+1,2*Py,2*Pz)       !this%muu(this%xi(Px), this%yi05(Py), this%zi05(Pz));
      muEzr= this%mup(2*Px+1,2*Py+1,2*Pz)     !this%muu(this%xi(Px), this%yi(Py), this%zi05(Pz));
      muEzl= this%mup(2*Px+1,2*(Py-1)+1,2*Pz) !this%muu(this%xi(Px), this%yi(Py-1), this%zi05(Pz));
      muEyr= this%mup(2*Px+1,2*Py,2*Pz+1)     !this%muu(this%xi(Px), this%yi05(Py), this%zi(Pz));
      muEyl= this%mup(2*Px+1,2*Py,2*(Pz-1)+1) !this%muu(this%xi(Px), this%yi05(Py), this%zi(Pz-1));
      DE = (muEzr*this%Ef%Z(Px,Py,Pz)-muEzl*this%Ef%Z(Px,Py-1,Pz))/this%hhy-(muEyr*this%Ef%Y(Px,Py,Pz)-muEyl*this%Ef%Y(Px,Py,Pz-1))/this%hhz;
      this%JHeff%X(Px,Py,Pz) = -( DE + muHx*(this%Hf%X(Px,Py,Pz) - this%Hfold%X(Px,Py,Pz))/cc/this%ht )*cc/4.0d0/Pi;
    case(2) !JHeff_Y
      muHy=  this%mup(2*Px,2*Py+1,2*Pz)        !this%muu(this%xi05(Px), this%yi(Py), this%zi05(Pz));
      muEzr= this%mup(2*Px+1,2*Py+1,2*Pz)      !this%muu(this%xi(Px), this%yi(Py), this%zi05(Pz));
      muEzl= this%mup(2*(Px-1)+1,2*Py+1,2*Pz)  !this%muu(this%xi(Px-1), this%yi(Py), this%zi05(Pz));
      muExr= this%mup(2*Px,2*Py+1,2*Pz+1)      !this%muu(this%xi05(Px), this%yi(Py), this%zi(Pz));
      muExl= this%mup(2*Px,2*Py+1,2*(Pz-1)+1)  !this%muu(this%xi05(Px), this%yi(Py), this%zi(Pz-1));
      DE = (muExr*this%Ef%X(Px,Py,Pz)-muExl*this%Ef%X(Px,Py,Pz-1))/this%hhz-(muEzr*this%Ef%Z(Px,Py,Pz)-muEzl*this%Ef%Z(Px-1,Py,Pz))/this%hhx;
      this%JHeff%Y(Px,Py,Pz) = -( DE + muHy*(this%Hf%Y(Px,Py,Pz) - this%Hfold%Y(Px,Py,Pz))/cc/this%ht )*cc/4.0d0/Pi;
    case(3) !JHeff_Z
      muHz=  this%mup(2*Px,2*Py,2*Pz+1)         !this%muu(this%xi05(Px), this%yi05(Py), this%zi(Pz));
      muEyr= this%mup(2*Px+1,2*Py,2*Pz+1)       !this%muu(this%xi(Px), this%yi05(Py), this%zi(Pz));
      muEyl= this%mup(2*(Px-1)+1,2*Py,2*Pz+1)   !this%muu(this%xi(Px-1), this%yi05(Py), this%zi(Pz));
      muExr= this%mup(2*Px,2*Py+1,2*Pz+1)       !this%muu(this%xi05(Px), this%yi(Py), this%zi(Pz));
      muExl= this%mup(2*Px,2*(Py-1)+1,2*Pz+1)   !this%muu(this%xi05(Px), this%yi(Py-1), this%zi(Pz));
      DE = (muEyr*this%Ef%Y(Px,Py,Pz)-muEyl*this%Ef%Y(Px-1,Py,Pz))/this%hhx-(muExr*this%Ef%X(Px,Py,Pz)-muExl*this%Ef%X(Px,Py-1,Pz))/this%hhy;
      this%JHeff%Z(Px,Py,Pz) = -( DE + muHz*(this%Hf%Z(Px,Py,Pz) - this%Hfold%Z(Px,Py,Pz))/cc/this%ht )*cc/4.0d0/Pi;
    case default ! illegal type
       fatalerr = 214;
       write(*,*) 'Alert! Incorrect use of Component in JHbuildPoint function!', component     
    end select
  end subroutine JHeffbuildPoint  


!JEeffbuildPointPoisson------------------------------------------------------------------
  subroutine  JEeffbuildPointPoisson(this, t, Px, Py, Pz, component)
    ! Build JE effective current at the point
    class(tProblem) :: this
    integer, intent(in) :: t;
    integer, intent(in) :: Px, Py, Pz, component;
    real*8 :: muEx, muEy, muEz;
    real*8 :: muHxr, muHxl, muHyr, muHyl, muHzr, muHzl;
    real*8 :: DH;
    select case (component)
    case(1) !JEeff_X
      muEx=  this%mup(2*Px,2*Py+1,2*Pz+1)     !this%muu(this%xi05(Px), this%yi(Py), this%zi(Pz));
      muHyr= this%mup(2*Px,2*Py+1,2*(Pz+1))   !this%muu(this%xi05(Px), this%yi(Py), this%zi05(Pz+1));
      muHyl= this%mup(2*Px,2*Py+1,2*Pz)       !this%muu(this%xi05(Px), this%yi(Py), this%zi05(Pz));
      muHzr= this%mup(2*Px,2*(Py+1),2*Pz+1)   !this%muu(this%xi05(Px), this%yi05(Py+1), this%zi(Pz));
      muHzl= this%mup(2*Px,2*Py,2*Pz+1)       !this%muu(this%xi05(Px), this%yi05(Py), this%zi(Pz));    
      DH = (muHzr*this%Hf%Z(Px,Py+1,Pz)-muHzl*this%Hf%Z(Px,Py,Pz))/this%hhy-(muHyr*this%Hf%Y(Px,Py,Pz+1)-muHyl*this%Hf%Y(Px,Py,Pz))/this%hhz;
      this%JEeff%X(Px,Py,Pz) = ( DH - (this%poisson%Ec1%X(Px,Py,Pz) - this%poisson%Ec0%X(Px,Py,Pz))/cc/this%ht)*cc/4.0d0/Pi;
   case(2) !JEeff_Y
      muEy=  this%mup(2*Px+1,2*Py,2*Pz+1)     !this%muu(this%xi(Px),this%yi05(Py),this%zi(Pz));
      muHxr= this%mup(2*Px+1,2*Py,2*(Pz+1))   !this%muu(this%xi(Px),this%yi05(Py),this%zi05(Pz+1));
      muHxl= this%mup(2*Px+1,2*Py,2*Pz)       !this%muu(this%xi(Px),this%yi05(Py),this%zi05(Pz));
      muHzr= this%mup(2*(Px+1),2*Py,2*Pz+1)   !this%muu(this%xi05(Px+1),this%yi05(Py),this%zi(Pz));
      muHzl= this%mup(2*Px,2*Py,2*Pz+1)       !this%muu(this%xi05(Px),this%yi05(Py),this%zi(Pz));
      DH = (muHxr*this%Hf%X(Px,Py,Pz+1)-muHxl*this%Hf%X(Px,Py,Pz))/this%hhz-(muHzr*this%Hf%Z(Px+1,Py,Pz)-muHzl*this%Hf%Z(Px,Py,Pz))/this%hhx;
      this%JEeff%Y(Px,Py,Pz) = ( DH - (this%poisson%Ec1%Y(Px,Py,Pz) - this%poisson%Ec0%Y(Px,Py,Pz))/cc/this%ht)*cc/4.0d0/Pi;
   case(3) !JEeff_Z
      muEz=  this%mup(2*Px+1,2*Py+1,2*Pz)      !this%muu(this%xi(Px),this%yi(Py),this%zi05(Pz));
      muHyr= this%mup(2*(Px+1),2*Py+1,2*Pz)    !this%muu(this%xi05(Px+1),this%yi(Py),this%zi05(Pz));
      muHyl= this%mup(2*Px,2*Py+1,2*Pz)        !this%muu(this%xi05(Px),this%yi(Py),this%zi05(Pz));
      muHxr= this%mup(2*Px+1,2*(Py+1),2*Pz)    !this%muu(this%xi(Px),this%yi05(Py+1),this%zi05(Pz));
      muHxl= this%mup(2*Px+1,2*Py,2*Pz)        !this%muu(this%xi(Px),this%yi05(Py),this%zi05(Pz));   
      DH = (muHyr*this%Hf%Y(Px+1,Py,Pz)-muHyl*this%Hf%Y(Px,Py,Pz))/this%hhx-(muHxr*this%Hf%X(Px,Py+1,Pz)-muHxl*this%Hf%X(Px,Py,Pz))/this%hhy;
      this%JEeff%Z(Px,Py,Pz) = ( DH - (this%poisson%Ec1%Z(Px,Py,Pz) - this%poisson%Ec0%Z(Px,Py,Pz))/cc/this%ht)*cc/4.0d0/Pi;
    case default ! illegal type
       fatalerr = 213;
       write(*,*) 'Alert! Incorrect use of Component in JEbuildPoint function!', component     
    end select
    1210 format(E14.7) 
  end subroutine JEeffbuildPointPoisson
  

!JEeffbuild----------------------------------------------------------------------
  subroutine  JEeffbuild(this, t)
    ! Build JE effective currents
    class(tProblem) :: this
    integer, intent(in) :: t;
    integer :: i, j, k;
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do k=1,this%Nz-1
       do j=1,this%Ny-1
          do i=1,this%Nx
             select case (this%effcurrentstype)  ! JEx
             case (0)
                call this%JEeffbuildPoint(t,i,j,k,1); 
             case (1)
                call this%JEeffbuildPointPoisson(t,i,j,k,1);
             case default
                fatalerr = 216;
                write(*,*) 'Alert! Incorrect Effective Currents Type in JEffbuild function!', this%effcurrentstype    
             end select   
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do k=1,this%Nz-1
       do j=1,this%Ny
          do i=1,this%Nx-1
             select case (this%effcurrentstype)  ! JEy
             case (0)
                call this%JEeffbuildPoint(t,i,j,k,2); 
             case (1)
                call this%JEeffbuildPointPoisson(t,i,j,k,2);
             case default
                fatalerr = 216;
                write(*,*) 'Alert! Incorrect Effective Currents Type in JEffbuild function!', this%effcurrentstype    
             end select 
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do k=1,this%Nz
       do j=1,this%Ny-1
          do i=1,this%Nx-1
             select case (this%effcurrentstype)  ! JEz
             case (0)
                call this%JEeffbuildPoint(t,i,j,k,3); 
             case (1)
                call this%JEeffbuildPointPoisson(t,i,j,k,3);
             case default
                fatalerr = 216;
                write(*,*) 'Alert! Incorrect Effective Currents Type in JEffbuild function!', this%effcurrentstype    
             end select 
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
    this%JEeff%X(this%Nx-1:this%Nx,:,:)=0.0d0;
    this%JEeff%X(1:2,:,:)=0.0d0;
    this%JEeff%X(:,this%Ny-1:this%Ny,:)=0.0d0;
    this%JEeff%X(:,0:1,:)=0.0d0;
    this%JEeff%X(:,:,this%Nz-1:this%Nz)=0.0d0;
    this%JEeff%X(:,:,0:1)=0.0d0;

    this%JEeff%Y(this%Nx-1:this%Nx,:,:)=0.0d0;
    this%JEeff%Y(0:1,:,:)=0.0d0;
    this%JEeff%Y(:,this%Ny-1:this%Ny,:)=0.0d0;
    this%JEeff%Y(:,1:2,:)=0.0d0;
    this%JEeff%Y(:,:,this%Nz-1:this%Nz)=0.0d0;
    this%JEeff%Y(:,:,0:1)=0.0d0;
      
    this%JEeff%Z(this%Nx-1:this%Nx,:,:)=0.0d0;
    this%JEeff%Z(0:1,:,:)=0.0d0;
    this%JEeff%Z(:,this%Ny-1:this%Ny,:)=0.0d0;
    this%JEeff%Z(:,0:1,:)=0.0d0;
    this%JEeff%Z(:,:,this%Nz-1:this%Nz)=0.0d0;
    this%JEeff%Z(:,:,1:2)=0.0d0;
  end subroutine JEeffbuild


!JHeffbuild-------------------------------------------------------------------------
  subroutine  JHeffbuild(this, t)
    ! Build JE effective current at the point
    class(tProblem) :: this
    integer, intent(in) :: t;
    integer :: i, j, k;
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do k=2,this%Nz-1
       do j=2,this%Ny-1
          do i=1,this%Nx-1
             call this%JHeffbuildPoint(t,i,j,k,1); ! JHx
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do k=2,this%Nz-1
       do j=1,this%Ny-1
          do i=2,this%Nx-1
             call this%JHeffbuildPoint(t,i,j,k,2); ! JHy
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do k=1,this%Nz-1
       do j=2,this%Ny-1
          do i=2,this%Nx-1
             call this%JHeffbuildPoint(t,i,j,k,3); ! JHz
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
  end subroutine JHeffbuild


!saveEffCurrentsToBuffer-------------------------------------------------------------------------
  subroutine saveEffCurrentsToBuffer(this)
    ! Save all effective currents to buffer
    class(tProblem) :: this
    integer :: i, j, k;
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do k=0,this%Nz
       do j=0,this%Ny
          do i=1,this%Nx
             this%buffer%Je%X(i+this%bufoffsetX,j+this%bufoffsetY,k+this%bufoffsetZ)=this%JEeff%X(i,j,k);!this%Je%X(i,j,k);   ! ! JEx
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do k=0,this%Nz
       do j=1,this%Ny
          do i=0,this%Nx
             this%buffer%Je%Y(i+this%bufoffsetX,j+this%bufoffsetY,k+this%bufoffsetZ)=this%JEeff%Y(i,j,k);!this%Je%Y(i,j,k);     ! ! JEy
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do k=1,this%Nz
       do j=0,this%Ny
          do i=0,this%Nx
             this%buffer%Je%Z(i+this%bufoffsetX,j+this%bufoffsetY,k+this%bufoffsetZ)=this%JEeff%Z(i,j,k);!this%Je%Z(i,j,k);       ! ! JEz
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do k=1,this%Nz
       do j=1,this%Ny
          do i=0,this%Nx
             this%buffer%Jh%X(i+this%bufoffsetX,j+this%bufoffsetY,k+this%bufoffsetZ)=this%JHeff%X(i,j,k);!0;       ! ! JHx
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do k=1,this%Nz
       do j=0,this%Ny
          do i=1,this%Nx
             this%buffer%Jh%Y(i+this%bufoffsetX,j+this%bufoffsetY,k+this%bufoffsetZ)=this%JHeff%Y(i,j,k);!0;         ! ! JHy
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do k=0,this%Nz
       do j=1,this%Ny
          do i=1,this%Nx
             this%buffer%Jh%Z(i+this%bufoffsetX,j+this%bufoffsetY,k+this%bufoffsetZ)=this%JHeff%Z(i,j,k);!0;         ! ! JHz
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
  end subroutine saveEffCurrentsToBuffer


!saveOriginalCurrentsToBuffer-------------------------------------------------------------------------
  subroutine saveOriginalCurrentsToBuffer(this)
    ! Save all original currents to buffer
    class(tProblem) :: this
    integer :: i, j, k;
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do k=0,this%Nz
       do j=0,this%Ny
          do i=1,this%Nx
             this%buffer%Je%X(i+this%bufoffsetX,j+this%bufoffsetY,k+this%bufoffsetZ)=this%Je%X(i,j,k);!this%Je%X(i,j,k);   ! ! JEx
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do k=0,this%Nz
       do j=1,this%Ny
          do i=0,this%Nx
             this%buffer%Je%Y(i+this%bufoffsetX,j+this%bufoffsetY,k+this%bufoffsetZ)=this%Je%Y(i,j,k);!this%Je%Y(i,j,k);     ! ! JEy
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do k=1,this%Nz
       do j=0,this%Ny
          do i=0,this%Nx
             this%buffer%Je%Z(i+this%bufoffsetX,j+this%bufoffsetY,k+this%bufoffsetZ)=this%Je%Z(i,j,k);!this%Je%Z(i,j,k);       ! ! JEz
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do k=1,this%Nz
       do j=1,this%Ny
          do i=0,this%Nx
             this%buffer%Jh%X(i+this%bufoffsetX,j+this%bufoffsetY,k+this%bufoffsetZ)=0; !0;       ! ! JHx
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do k=1,this%Nz
       do j=0,this%Ny
          do i=1,this%Nx
             this%buffer%Jh%Y(i+this%bufoffsetX,j+this%bufoffsetY,k+this%bufoffsetZ)=0; !0;         ! ! JHy
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do k=0,this%Nz
       do j=1,this%Ny
          do i=1,this%Nx
             this%buffer%Jh%Z(i+this%bufoffsetX,j+this%bufoffsetY,k+this%bufoffsetZ)=0; !0;         ! ! JHz
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
  end subroutine saveOriginalCurrentsToBuffer  


!UploadFieldsToPoisson------------------------------------------------------------------------------------
  subroutine UploadFieldsToPoisson(this)
  ! Upload current fields to Poisson problem
      class(tProblem) :: this
      integer :: i,j,k
      do k=0,this%Nz
         do j=0,this%Ny
            do i=1,this%Nx
               this%poisson%Ef%X(i,j,k)=this%Ef%X(i,j,k); ! Ex
            enddo
         enddo
      enddo
      do k=0,this%Nz
         do j=1,this%Ny
            do i=0,this%Nx
               this%poisson%Ef%Y(i,j,k)=this%Ef%Y(i,j,k); ! Ey
            enddo
         enddo
      enddo
      do k=1,this%Nz
         do j=0,this%Ny
            do i=0,this%Nx
               this%poisson%Ef%Z(i,j,k)=this%Ef%Z(i,j,k); ! Ez
            enddo
         enddo
      enddo
  end subroutine UploadFieldsToPoisson


!addEcompositesToBuffer-------------------------------------------------------------------------
  subroutine addEcompositesToBuffer(this)
    ! Save all effective currents to buffer
    class(tProblem) :: this
    integer :: i, j, k;
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do k=0,this%Nz
       do j=0,this%Ny
          do i=1,this%Nx
             this%buffer%Ec0%X(i+this%bufoffsetX,j+this%bufoffsetY,k+this%bufoffsetZ)=this%poisson%Ec0%X(i,j,k)!Ec0
             this%buffer%Ec1%X(i+this%bufoffsetX,j+this%bufoffsetY,k+this%bufoffsetZ)=this%poisson%Ec1%X(i,j,k)!Ec1
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do k=0,this%Nz
       do j=1,this%Ny
          do i=0,this%Nx
             this%buffer%Ec0%Y(i+this%bufoffsetX,j+this%bufoffsetY,k+this%bufoffsetZ)=this%poisson%Ec0%Y(i,j,k)!Ec0
             this%buffer%Ec1%Y(i+this%bufoffsetX,j+this%bufoffsetY,k+this%bufoffsetZ)=this%poisson%Ec1%Y(i,j,k)!Ec1
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do k=1,this%Nz
       do j=0,this%Ny
          do i=0,this%Nx
             this%buffer%Ec0%Z(i+this%bufoffsetX,j+this%bufoffsetY,k+this%bufoffsetZ)=this%poisson%Ec0%Z(i,j,k)!Ec0
             this%buffer%Ec1%Z(i+this%bufoffsetX,j+this%bufoffsetY,k+this%bufoffsetZ)=this%poisson%Ec1%Z(i,j,k)!Ec1
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
  end subroutine addEcompositesToBuffer

    
!UploadMainMuToBuffer-------------------------------------------------------------------------
  subroutine UploadMainMuToBuffer(this)
    ! Save all effective currents to buffer
    class(tProblem) :: this
    integer :: i, j, k;
      do k=0,this%Nz
         do j=0,this%Ny
            do i=1,this%Nx
               this%buffer%mainMu%X(i+this%bufoffsetX,j+this%bufoffsetY,k+this%bufoffsetZ)=this%muu(this%xi05(i),this%yi(j),this%zi(k));
             enddo
         enddo
      enddo
      do k=0,this%Nz
         do j=1,this%Ny
            do i=0,this%Nx
               this%buffer%mainMu%Y(i+this%bufoffsetX,j+this%bufoffsetY,k+this%bufoffsetZ)=this%muu(this%xi(i),this%yi05(j),this%zi(k));
             enddo
         enddo
      enddo
      do k=1,this%Nz
         do j=0,this%Ny
            do i=0,this%Nx
               this%buffer%mainMu%Z(i+this%bufoffsetX,j+this%bufoffsetY,k+this%bufoffsetZ)=this%muu(this%xi(i),this%yi(j),this%zi05(k));
            enddo
         enddo
      enddo      
    end subroutine UploadMainMuToBuffer


!GetDifference-------------------------------------------------------------------------
  real*8 function GetDifference(this, component)
    ! Save all effective currents to buffer
    class(tProblem) :: this
    integer :: component;
    integer :: i, j, k;
    real*8 :: mx;
    mx = 0.0d0;
    select case (component)
    case(0)
       do k=0,this%Nz
         do j=0,this%Ny
            do i=1,this%Nx
               if ((this%Ef%X(i,j,k)-this%Efold%X(i,j,k))>mx) then
                  mx = this%Ef%X(i,j,k)-this%Efold%X(i,j,k);
               endif
             enddo
         enddo
      enddo
    case(1)
       do k=0,this%Nz
         do j=1,this%Ny
            do i=0,this%Nx
               if ((this%Ef%Y(i,j,k)-this%Efold%Y(i,j,k))>mx) then
                  mx = this%Ef%Y(i,j,k)-this%Efold%Y(i,j,k);
               endif
             enddo
         enddo
      enddo
    case(2)
       do k=1,this%Nz
         do j=0,this%Ny
            do i=0,this%Nx
               if ((this%Ef%Z(i,j,k)-this%Efold%Z(i,j,k))>mx) then
                  mx = this%Ef%Z(i,j,k)-this%Efold%Z(i,j,k);
               endif
             enddo
         enddo
      enddo
   end select
   GetDifference = mx;
    end function GetDifference   


!WriteFieldToFile-------------------------------------------------------------------------
    subroutine WriteFieldToFile(this, filenum, component)
    ! Write all field components to file
      class(tProblem) :: this
      integer :: filenum;
      integer :: component;
      integer :: i,j,k;
      select case(component)
      case(0)
         do k=0,this%Nz
            do j=0,this%Ny
               do i=1,this%Nx
                  if (this%Ef%X(i,j,k)<1.0E-16) then
                        write(unit=filenum, fmt=3400) 0.0d0, i, j, k
                  else
                        write(unit=filenum, fmt=3400) this%Ef%X(i,j,k), i, j, k
                  endif
               enddo
            enddo
         enddo
      case(1)
         do k=0,this%Nz
            do j=1,this%Ny
               do i=0,this%Nx
                  if (this%Ef%Y(i,j,k)<1.0E-16) then
                        write(unit=filenum, fmt=3400) 0.0d0, i, j, k
                  else
                        write(unit=filenum, fmt=3400) this%Ef%Y(i,j,k), i, j, k
                  endif
               enddo
            enddo
         enddo
      case(2)
         do k=1,this%Nz
            do j=0,this%Ny
               do i=0,this%Nx
                  if (this%Ef%Z(i,j,k)<1.0E-16) then
                        write(unit=filenum, fmt=3400) 0.0d0, i, j, k
                  else
                        write(unit=filenum, fmt=3400) this%Ef%Z(i,j,k), i, j, k
                  endif
               enddo
            enddo
         enddo
      case(3)
         do k=1,this%Nz
            do j=1,this%Ny
               do i=0,this%Nx
                  if (this%Hf%X(i,j,k)<1.0E-16) then
                        write(unit=filenum, fmt=3400) 0.0d0, i, j, k
                  else
                        write(unit=filenum, fmt=3400) this%Hf%X(i,j,k), i, j, k
                  endif
               enddo
            enddo
         enddo
      case(4)
         do k=1,this%Nz
            do j=0,this%Ny
               do i=1,this%Nx
                  if (this%Hf%Y(i,j,k)<1.0E-16) then
                        write(unit=filenum, fmt=3400) 0.0d0, i, j, k
                  else
                        write(unit=filenum, fmt=3400) this%Hf%Y(i,j,k), i, j, k
                  endif
               enddo
            enddo
         enddo
      case(5)
         do k=0,this%Nz
            do j=1,this%Ny
               do i=1,this%Nx
                  if (this%Hf%Z(i,j,k)<1.0E-16) then
                        write(unit=filenum, fmt=3400) 0.0d0, i, j, k
                  else
                        write(unit=filenum, fmt=3400) this%Hf%Z(i,j,k), i, j, k
                  endif
               enddo
            enddo
         enddo   
      end select
      3400 format(E22.15,',',I2,',',I2,',',I2);  
    end subroutine WriteFieldToFile


!initHWEvaAuxVariables---------------------------
    subroutine initHWEvaAuxVariables(this)
    ! Allocate arrays for HWeva ABC aux variables
      class(tProblem) :: this;
      integer :: i, j, info;
      real*8, allocatable :: ipiv(:), Ematr(:,:);

      if (this%hwEvaE>0) then
         allocate(this%hwX(0:this%hwEvaP+this%hwEvaE+1));
         allocate(this%hwY(0:this%hwEvaP+this%hwEvaE+1));
         allocate(this%hwZ(0:this%hwEvaP+this%hwEvaE+1));
         allocate(this%hwXold(0:this%hwEvaP+this%hwEvaE+1));
         allocate(this%hwYold(0:this%hwEvaP+this%hwEvaE+1));
         allocate(this%hwZold(0:this%hwEvaP+this%hwEvaE+1));
         allocate(this%hwXold2(0:this%hwEvaP+this%hwEvaE+1));
         allocate(this%hwYold2(0:this%hwEvaP+this%hwEvaE+1));
         allocate(this%hwZold2(0:this%hwEvaP+this%hwEvaE+1));
         do i=0,this%hwEvaP+this%hwEvaE+1
            call this%hwX(i)%auxmesh_init(this%Nx, this%Ny, this%Nz, 0, this%problem_id);
            call this%hwY(i)%auxmesh_init(this%Nx, this%Ny, this%Nz, 1, this%problem_id);
            call this%hwZ(i)%auxmesh_init(this%Nx, this%Ny, this%Nz, 2, this%problem_id);
            call this%hwXold(i)%auxmesh_init(this%Nx, this%Ny, this%Nz, 0, this%problem_id);
            call this%hwYold(i)%auxmesh_init(this%Nx, this%Ny, this%Nz, 1, this%problem_id);
            call this%hwZold(i)%auxmesh_init(this%Nx, this%Ny, this%Nz, 2, this%problem_id);
            call this%hwXold2(i)%auxmesh_init(this%Nx, this%Ny, this%Nz, 0, this%problem_id);
            call this%hwYold2(i)%auxmesh_init(this%Nx, this%Ny, this%Nz, 1, this%problem_id);
            call this%hwZold2(i)%auxmesh_init(this%Nx, this%Ny, this%Nz, 2, this%problem_id);
         end do
      else
         allocate(this%hwX(0:this%hwEvaP));
         allocate(this%hwY(0:this%hwEvaP));
         allocate(this%hwZ(0:this%hwEvaP));
         allocate(this%hwXold(0:this%hwEvaP));
         allocate(this%hwYold(0:this%hwEvaP));
         allocate(this%hwZold(0:this%hwEvaP));
         allocate(this%hwXold2(0:this%hwEvaP));
         allocate(this%hwYold2(0:this%hwEvaP));
         allocate(this%hwZold2(0:this%hwEvaP));
         do i=0,this%hwEvaP
            call this%hwX(i)%auxmesh_init(this%Nx, this%Ny, this%Nz, 0, this%problem_id);
            call this%hwY(i)%auxmesh_init(this%Nx, this%Ny, this%Nz, 1, this%problem_id);
            call this%hwZ(i)%auxmesh_init(this%Nx, this%Ny, this%Nz, 2, this%problem_id);
            call this%hwXold(i)%auxmesh_init(this%Nx, this%Ny, this%Nz, 0, this%problem_id);
            call this%hwYold(i)%auxmesh_init(this%Nx, this%Ny, this%Nz, 1, this%problem_id);
            call this%hwZold(i)%auxmesh_init(this%Nx, this%Ny, this%Nz, 2, this%problem_id);
            call this%hwXold2(i)%auxmesh_init(this%Nx, this%Ny, this%Nz, 0, this%problem_id);
            call this%hwYold2(i)%auxmesh_init(this%Nx, this%Ny, this%Nz, 1, this%problem_id);
            call this%hwZold2(i)%auxmesh_init(this%Nx, this%Ny, this%Nz, 2, this%problem_id);
         end do
      endif 

      ! Set angles for HWEva
      allocate(this%hwAng(0:this%hwEvaP));
      this%hwAng(0) = 0.0d0;
      if (this%hwEvaP>0) then
         do i=1,this%hwEvaP
            this%hwAng(i) = 0.0d0 !i*90.0d0/(this%hwEvaP+1.0d0);
         enddo
      endif

      ! Set basic coefficients for HW
      allocate(this%hwAC(0:this%hwEvaP));
      do i=0,this%hwEvaP
         this%hwAC(i) = cos(this%hwAng(i)/180.0d0*PI);
      enddo

      ! Set sigmas for HWEva
      if (this%hwEvaE>0) then  
         allocate(this%hwSigmas(1:this%hwEvaE));
         do i=1,this%hwEvaE
            this%hwSigmas(i) = 1.0d0;
         enddo
      endif   

      !Fill coefficients for HW

      if (this%hwEvaE>0) then ! With evanescent waves
         allocate(this%hwL(0:this%hwEvaP+1, 0:this%hwEvaP+1));
         allocate(this%hwM(0:this%hwEvaP+1, 0:this%hwEvaP+1));
         allocate(this%hwA(0:this%hwEvaP+this%hwEvaE+1, 0:this%hwEvaP+this%hwEvaE+1));
         allocate(this%hwB(0:this%hwEvaP+this%hwEvaE+1, 0:this%hwEvaP+this%hwEvaE+1));
         allocate(this%hwC(0:this%hwEvaP+this%hwEvaE+1, 0:this%hwEvaP+this%hwEvaE+1));
         allocate(this%hwD(0:this%hwEvaP+this%hwEvaE+1, 0:this%hwEvaP+this%hwEvaE+1));
         allocate(Ematr(0:this%hwEvaP+this%hwEvaE+1, 0:this%hwEvaP+this%hwEvaE+1));
         allocate(this%hwN(1:this%hwEvaE, 1:this%hwEvaE));
         allocate(this%hwS(1:this%hwEvaE, 1:this%hwEvaE));
         allocate(this%hwV(1:this%hwEvaE, 1:this%hwEvaE));
         do i=0,this%hwEvaP+1
            do j=0,this%hwEvaP+1
               this%hwL(i,j)=0.0d0;
               this%hwM(i,j)=0.0d0;
            enddo
         enddo
         do i=0,this%hwEvaP+this%hwEvaE+1
            do j=0,this%hwEvaP+this%hwEvaE+1
               this%hwA(i,j)=0.0d0;
               this%hwB(i,j)=0.0d0;
               this%hwC(i,j)=0.0d0;
               this%hwD(i,j)=0.0d0;
               Ematr(i,j)=0.0d0;
               if (i==j) then
                  Ematr(i,j)=1.0d0;
               endif
            enddo
         enddo
         do i=1,this%hwEvaE
            do j=1,this%hwEvaE
               this%hwN(i,j)=0.0d0;
               this%hwS(i,j)=0.0d0;
               this%hwV(i,j)=0.0d0;
            enddo
         enddo
      else                   ! No evanescent waves
         allocate(this%hwL(0:this%hwEvaP, 0:this%hwEvaP));
         allocate(this%hwM(0:this%hwEvaP, 0:this%hwEvaP));
         allocate(this%hwA(0:this%hwEvaP, 0:this%hwEvaP));
         allocate(this%hwB(0:this%hwEvaP, 0:this%hwEvaP));
         allocate(this%hwC(0:this%hwEvaP, 0:this%hwEvaP));
         allocate(this%hwD(0:this%hwEvaP, 0:this%hwEvaP));
         allocate(Ematr(0:this%hwEvaP, 0:this%hwEvaP));
         do i=0,this%hwEvaP
            do j=0,this%hwEvaP
               this%hwL(i,j)=0.0d0;
               this%hwM(i,j)=0.0d0;
               this%hwA(i,j)=0.0d0;
               this%hwB(i,j)=0.0d0;
               this%hwC(i,j)=0.0d0;
               this%hwD(i,j)=0.0d0;
               Ematr(i,j)=0.0d0;
               if (i==j) then
                  Ematr(i,j)=1.0d0;
               endif
            enddo
         enddo
      endif
      

      ! Filling L and M arrays
      if (this%hwEvaP>0) then
         this%hwL(1,0) = 2.0d0*this%hwAC(1)*(1.0d0-this%hwAC(0)**2);
         this%hwL(1,1) = this%hwAC(0)*(1.0d0+2.0d0*this%hwAC(0)*this%hwAC(1)+this%hwAC(1)**2);
         this%hwL(1,2) = this%hwAC(0)*(1.0d0-this%hwAC(1)**2);
         this%hwM(1,0) = 2.0d0*this%hwAC(1);
         this%hwM(1,1) = this%hwAC(0);
         this%hwM(1,2) = this%hwAC(0);      
         if (this%hwEvaP>1) then
            do i=2,this%hwEvaP
               this%hwL(i,i-1) = this%hwAC(i)*(1.0d0-this%hwAC(i-1)**2);
               this%hwL(i,i) = this%hwAC(i)*(1.0d0+this%hwAC(i-1)**2)+this%hwAC(i-1)*(1.0d0+this%hwAC(i)**2);
               this%hwM(i,i-1) = this%hwAC(i);
               this%hwM(i,i) = this%hwAC(i)+this%hwAC(i-1);
               if (this%hwEvaE>0) then
                  this%hwM(i,i+1) = this%hwAC(i-1);
                  this%hwL(i,i+1) = this%hwAC(i-1)*(1.0d0-this%hwAC(i)**2);
               else
                  if (i/=this%hwEvaP) then
                     this%hwM(i,i+1) = this%hwAC(i-1);
                     this%hwL(i,i+1) = this%hwAC(i-1)*(1.0d0-this%hwAC(i)**2);
                  endif
               endif
            enddo
         endif
      endif
      

      ! Filling N, S, V arrays
      if (this%hwEvaE>0) then
         this%hwL(this%hwEvaP+1,this%hwEvaP) = (1.0d0-this%hwAC(this%hwEvaP)**2)/this%hwAC(this%hwEvaP);
         this%hwL(this%hwEvaP+1,this%hwEvaP+1) = (1.0d0+this%hwAC(this%hwEvaP)**2)/this%hwAC(this%hwEvaP);
         this%hwM(this%hwEvaP+1,this%hwEvaP) = 1.0d0/this%hwAC(this%hwEvaP);
         this%hwM(this%hwEvaP+1,this%hwEvaP+1) = 1.0d0/this%hwAC(this%hwEvaP);
         this%hwN(1,1) = 1.0d0/(this%hwAC(this%hwEvaP)*this%hwSigmas(1));
         this%hwN(1,2) = 1.0d0/(this%hwAC(this%hwEvaP)*this%hwSigmas(1));       
         this%hwV(1,1) = 1.0d0/(this%hwAC(this%hwEvaP)*this%hwSigmas(1));
         this%hwV(1,2) = 1.0d0/(this%hwAC(this%hwEvaP)*this%hwSigmas(1));
         this%hwS(1,1) = -this%hwSigmas(1)/this%hwAC(this%hwEvaP);
         this%hwS(1,2) = this%hwSigmas(1)/this%hwAC(this%hwEvaP);
         if (this%hwEvaE>1) then
            do i=2,this%hwEvaE
               this%hwN(i,i-1) = 1.0d0/this%hwSigmas(i-1);
               this%hwN(i,i) = 1.0d0/this%hwSigmas(i-1) + 1.0d0/this%hwSigmas(i);
               this%hwV(i,i-1) = 1.0d0/this%hwSigmas(i-1);
               this%hwV(i,i) = 1.0d0/this%hwSigmas(i-1) + 1.0d0/this%hwSigmas(i);
               this%hwS(i,i-1) = this%hwSigmas(i-1);
               this%hwS(i,i) = -(this%hwSigmas(i-1) + this%hwSigmas(i));
               if (i/=this%hwEvaE) then
                  this%hwN(i,i+1) = 1.0d0/this%hwSigmas(i);
                  this%hwV(i,i+1) = 1.0d0/this%hwSigmas(i);
                  this%hwS(i,i+1) = this%hwSigmas(i);
               endif
            enddo
         endif
      endif

      ! Fill A, B, C, D with L and M
      if (this%hwEvaE>0) then
         do i=0,this%hwEvaP+1
            do j=0,this%hwEvaP+1
               this%hwA(i,j) = this%hwL(i,j);
               this%hwB(i,j) = -2.0d0*this%hwL(i,j);
               this%hwC(i,j) = this%hwL(i,j);
               this%hwD(i,j) = (cc**2)*(this%ht**2)*this%hwM(i,j);
            enddo
         enddo
      else
         do i=0,this%hwEvaP
            do j=0,this%hwEvaP
               this%hwA(i,j) = this%hwL(i,j);
               this%hwB(i,j) = -2.0d0*this%hwL(i,j);
               this%hwC(i,j) = this%hwL(i,j);
               this%hwD(i,j) = (cc**2)*(this%ht**2)*this%hwM(i,j);
            enddo
         enddo 
      endif
         
      ! Fill A, B, C, D with N, V, S
      if (this%hwEvaE>0) then
         do i=1,this%hwEvaE
            do j=1,this%hwEvaE
               this%hwA(i+this%hwEvaP+1,j+this%hwEvaP) = this%hwN(i,j);
               this%hwB(i+this%hwEvaP+1,j+this%hwEvaP) = -2.0d0*this%hwN(i,j)-(cc**2)*(this%ht**2)*this%hwS(i,j);
               this%hwC(i+this%hwEvaP+1,j+this%hwEvaP) = this%hwN(i,j);
               this%hwD(i+this%hwEvaP+1,j+this%hwEvaP) = (cc**2)*(this%ht**2)*this%hwV(i,j);
            enddo
         enddo   
      endif

      ! More koefficients for sean variable
!      this%hwA(0, 0) = 3.0d0*this%hwAC(0);
!      this%hwA(0, 1) = -3.0d0*this%hwAC(0);
!      this%hwB(0, 0) = -4.0d0*this%hwAC(0);
!      this%hwB(0, 1) = 4.0d0*this%hwAC(0);
!      this%hwC(0, 0) = this%hwAC(0);
!      this%hwC(0, 1) = -this%hwAC(0);

      select case(this%dtmode)
      case(1)
          this%hwA(0,0)=1.0d0*this%hwAC(0);
          this%hwB(0,0)=-1.0d0*this%hwAC(0);
          this%hwC(0,0)=0.0d0;
          if (this%hwEvaP>0) then
             this%hwA(0,1)=-1.0d0;
             this%hwB(0,1)=1.0d0;
             this%hwC(0,1)=0.0d0;
          endif
      case(2)
          this%hwA(0, 0) = 3.0d0*this%hwAC(0);
          this%hwB(0, 0) = -4.0d0*this%hwAC(0);
          this%hwC(0, 0) = this%hwAC(0);
          if (this%hwEvaP>0) then
             this%hwA(0, 1) = -3.0d0*this%hwAC(0);
             this%hwB(0, 1) = 4.0d0*this%hwAC(0);
             this%hwC(0, 1) = -this%hwAC(0);
          endif
      end select

      if (this%hwEvaE>0) then
         this%hwA(this%hwEvaP+1, this%hwEvaP+this%hwEvaE+1) = this%hwAC(this%hwEvaP);
         this%hwA(this%hwEvaP+2, this%hwEvaP+this%hwEvaE+1) = -3.0d0/2.0d0*cc*this%ht;
         this%hwB(this%hwEvaP+1, this%hwEvaP+this%hwEvaE+1) = -2.0d0*this%hwAC(this%hwEvaP);
         this%hwB(this%hwEvaP+2, this%hwEvaP+this%hwEvaE+1) = 2.0d0*cc*this%ht;
         this%hwC(this%hwEvaP+1, this%hwEvaP+this%hwEvaE+1) = this%hwAC(this%hwEvaP);
         this%hwC(this%hwEvaP+2, this%hwEvaP+this%hwEvaE+1) = -1.0d0/2.0d0*cc*this%ht;
      endif   

      write(600,*) 'Angles';
      do i=0,this%hwEvaP
         WRITE(600,520) (this%hwAng(i))
      enddo

      write(600,*) 'A koeffs';
      do i=0,this%hwEvaP
         WRITE(600,520) (this%hwAC(i))
      enddo

      if (this%hwEvaE>0) then
         write(600,*) 'L array';
         do i=0,this%hwEvaP+1
            do j=0,this%hwEvaP+1
               WRITE(600, 520, advance="no") this%hwL(i,j)
            enddo
            write(600,*) '';
         enddo
         write(600,*) 'M array';
         do i=0,this%hwEvaP+1
            do j=0,this%hwEvaP+1
               WRITE(600, 520, advance="no") this%hwM(i,j)
            enddo
            write(600,*) '';
         enddo
         write(600,*) 'N array';
         do i=1,this%hwEvaE
            do j=1,this%hwEvaE
               WRITE(600, 520, advance="no") this%hwN(i,j)
            enddo
            write(600,*) '';
         enddo
         write(600,*) 'V array';
         do i=1,this%hwEvaE
            do j=1,this%hwEvaE
               WRITE(600, 520, advance="no") this%hwV(i,j)
            enddo
            write(600,*) '';
         enddo
         write(600,*) 'S array';
         do i=1,this%hwEvaE
            do j=1,this%hwEvaE
               WRITE(600, 520, advance="no") this%hwS(i,j)
            enddo
            write(600,*) '';
         enddo
         write(600,*) 'A array';
         do i=0,this%hwEvaP+this%hwEvaE+1
            do j=0,this%hwEvaP+this%hwEvaE+1
               WRITE(600, 520, advance="no") this%hwA(i,j)
            enddo
            write(600,*) '';
         enddo
         write(600,*) 'B array';
         do i=0,this%hwEvaP+this%hwEvaE+1
            do j=0,this%hwEvaP+this%hwEvaE+1
               WRITE(600, 520, advance="no") this%hwB(i,j)
            enddo
            write(600,*) '';
         enddo
         write(600,*) 'C array';
         do i=0,this%hwEvaP+this%hwEvaE+1
            do j=0,this%hwEvaP+this%hwEvaE+1
               WRITE(600, 520, advance="no") this%hwC(i,j)
            enddo
            write(600,*) '';
         enddo
         write(600,*) 'D array';
         do i=0,this%hwEvaP+this%hwEvaE+1
            do j=0,this%hwEvaP+this%hwEvaE+1
               WRITE(600, 520, advance="no") this%hwD(i,j)
            enddo
            write(600,*) '';
         enddo
      else
         write(600,*) 'L array';
         do i=0,this%hwEvaP
            do j=0,this%hwEvaP
               WRITE(600, 520, advance="no") this%hwL(i,j)
            enddo
            write(600,*) '';
         enddo
         write(600,*) 'M array';
         do i=0,this%hwEvaP
            do j=0,this%hwEvaP
               WRITE(600, 520, advance="no") this%hwM(i,j)
            enddo
            write(600,*) '';
         enddo
         write(600,*) 'A array';
         do i=0,this%hwEvaP
            do j=0,this%hwEvaP
               WRITE(600, 520, advance="no") this%hwA(i,j)
            enddo
            write(600,*) '';
         enddo
         write(600,*) 'B array';
         do i=0,this%hwEvaP
            do j=0,this%hwEvaP
               WRITE(600, 520, advance="no") this%hwB(i,j)
            enddo
            write(600,*) '';
         enddo
         write(600,*) 'C array';
         do i=0,this%hwEvaP
            do j=0,this%hwEvaP
               WRITE(600, 520, advance="no") this%hwC(i,j)
            enddo
            write(600,*) '';
         enddo
         write(600,*) 'D array';
         do i=0,this%hwEvaP
            do j=0,this%hwEvaP
               WRITE(600, 520, advance="no") this%hwD(i,j)
            enddo
            write(600,*) '';
         enddo
      endif
      
520   format (E11.4, ' ');      

      if (this%hwEvaE>0) then
         allocate(ipiv(0:this%hwEvaP+this%hwEvaE+1));
         call dgesv(this%hwEvaP+this%hwEvaE+2,this%hwEvaP+this%hwEvaE+2,this%hwA,this%hwEvaP+this%hwEvaE+2,ipiv,Ematr,this%hwEvaP+this%hwEvaE+2,info);
         do i=0,this%hwEvaP+this%hwEvaE+1 
            do j=0,this%hwEvaP+this%hwEvaE+1
               this%hwA(i,j)=Ematr(i,j);
            enddo
         enddo
         write(600,*) 'A^-1 array';
         do i=0,this%hwEvaP+this%hwEvaE+1
            do j=0,this%hwEvaP+this%hwEvaE+1
               WRITE(600, 520, advance="no") this%hwA(i,j)
            enddo
            write(600,*) '';
         enddo
         allocate(this%hwFnew(0:this%hwEvaP+this%hwEvaE+1 ), this%hwF(0:this%hwEvaP+this%hwEvaE+1 ), this%hwFold(0:this%hwEvaP+this%hwEvaE+1 ), this%hwGF(0:this%hwEvaP+this%hwEvaE+1 ));
         do i=0,this%hwEvaP+this%hwEvaE+1
            this%hwFnew(i)=0.0d0;
            this%hwF(i)=0.0d0;
            this%hwFold(i)=0.0d0;
            this%hwGF(i)=0.0d0;
         enddo
      else
         allocate(ipiv(0:this%hwEvaP));
         call dgesv(this%hwEvaP+1,this%hwEvaP+1,this%hwA,this%hwEvaP+1,ipiv,Ematr,this%hwEvaP+1,info);
         do i=0,this%hwEvaP
            do j=0,this%hwEvaP
               this%hwA(i,j)=Ematr(i,j);
            enddo
         enddo
         write(600,*) 'A^-1 array';
         do i=0,this%hwEvaP
            do j=0,this%hwEvaP
               WRITE(600, 520, advance="no") this%hwA(i,j)
            enddo
            write(600,*) '';
         enddo
         allocate(this%hwFnew(0:this%hwEvaP), this%hwF(0:this%hwEvaP), this%hwFold(0:this%hwEvaP), this%hwGF(0:this%hwEvaP));
         do i=0,this%hwEvaP
            this%hwFnew(i)=0.0d0;
            this%hwF(i)=0.0d0;
            this%hwFold(i)=0.0d0;
            this%hwGF(i)=0.0d0;
         enddo
      endif
      deallocate(ipiv);
      deallocate(Ematr);
      !call this%initHWEvaEdgeAuxVariables(); 
    end subroutine initHWEvaAuxVariables


!initHWEvaEdgeAuxVariables---------------------------
    subroutine initHWEvaEdgeAuxVariables(this)
    ! Allocate arrays for HWeva ABC Edge Aux Variables
      class(tProblem) :: this;
      integer :: i, j, info;
      real*8, allocatable :: ipiv(:), Ematr(:,:);
      integer :: im, ia;
      if (this%hwEvaE>0) then
         im = this%hwEvaP+this%hwEvaE+1;
      else
         im = this%hwEvaP;
      endif

      if (this%hwEvaEa>0) then
         ia = this%hwEvaPa+this%hwEvaEa+1;
      else
         ia = this%hwEvaPa;
      endif 

      allocate(this%hwYEdgeNX(0:im,0:ia));
      allocate(this%hwYEdgeNXold(0:im,0:ia));
      allocate(this%hwYEdgeNXold2(0:im,0:ia));
      do i=0,im
         do j=0,ia
            call this%hwYEdgeNX(i,j)%EdgeAuxMesh_init(this%Ny, this%Nz, 1, this%problem_id);
            call this%hwYEdgeNXold(i,j)%EdgeAuxMesh_init(this%Ny, this%Nz, 1, this%problem_id);
            call this%hwYEdgeNXold2(i,j)%EdgeAuxMesh_init(this%Ny, this%Nz, 1, this%problem_id);
         enddo
      enddo

      ! Set angles for HWEva on Edges
      allocate(this%hwAngA(0:this%hwEvaPa));
      this%hwAngA(0) = 0.0d0;
      if (this%hwEvaPa>0) then
         do i=1,this%hwEvaPa
            this%hwAngA(i) = 0.0d0 !i*90.0d0/(this%hwEvaP+1.0d0);
         enddo
      endif

      ! Set basic coefficients for HW on Edges
      allocate(this%hwACa(0:this%hwEvaPa));
      do i=0,this%hwEvaPa
         this%hwACa(i) = cos(this%hwAngA(i)/180.0d0*PI);
      enddo

      ! Set sigmas for HWEva on Edges
      if (this%hwEvaEa>0) then  
         allocate(this%hwSigmasA(1:this%hwEvaEa));
         do i=1,this%hwEvaEa
            this%hwSigmasA(i) = 1.0d0;
         enddo
      endif   

      !Fill coefficients for HW

      if (this%hwEvaEa>0) then ! With evanescent waves
         allocate(this%hwLa(0:this%hwEvaPa+1, 0:this%hwEvaPa+1));
         allocate(this%hwMa(0:this%hwEvaPa+1, 0:this%hwEvaPa+1));
         allocate(this%hwAa(0:this%hwEvaPa+this%hwEvaEa+1, 0:this%hwEvaPa+this%hwEvaEa+1));
         allocate(this%hwBa(0:this%hwEvaPa+this%hwEvaEa+1, 0:this%hwEvaPa+this%hwEvaEa+1));
         allocate(this%hwCa(0:this%hwEvaPa+this%hwEvaEa+1, 0:this%hwEvaPa+this%hwEvaEa+1));
         allocate(this%hwDa(0:this%hwEvaPa+this%hwEvaEa+1, 0:this%hwEvaPa+this%hwEvaEa+1));
         allocate(Ematr(0:this%hwEvaPa+this%hwEvaEa+1, 0:this%hwEvaPa+this%hwEvaEa+1));
         allocate(this%hwNa(1:this%hwEvaEa, 1:this%hwEvaEa));
         allocate(this%hwSa(1:this%hwEvaEa, 1:this%hwEvaEa));
         allocate(this%hwVa(1:this%hwEvaEa, 1:this%hwEvaEa));
         do i=0,this%hwEvaPa+1
            do j=0,this%hwEvaPa+1
               this%hwLa(i,j)=0.0d0;
               this%hwMa(i,j)=0.0d0;
            enddo
         enddo
         do i=0,this%hwEvaPa+this%hwEvaEa+1
            do j=0,this%hwEvaPa+this%hwEvaEa+1
               this%hwAa(i,j)=0.0d0;
               this%hwBa(i,j)=0.0d0;
               this%hwCa(i,j)=0.0d0;
               this%hwDa(i,j)=0.0d0;
               Ematr(i,j)=0.0d0;
               if (i==j) then
                  Ematr(i,j)=1.0d0;
               endif
            enddo
         enddo
         do i=1,this%hwEvaEa
            do j=1,this%hwEvaEa
               this%hwNa(i,j)=0.0d0;
               this%hwSa(i,j)=0.0d0;
               this%hwVa(i,j)=0.0d0;
            enddo
         enddo
      else                   ! No evanescent waves
         allocate(this%hwLa(0:this%hwEvaPa, 0:this%hwEvaPa));
         allocate(this%hwMa(0:this%hwEvaPa, 0:this%hwEvaPa));
         allocate(this%hwAa(0:this%hwEvaPa, 0:this%hwEvaPa));
         allocate(this%hwBa(0:this%hwEvaPa, 0:this%hwEvaPa));
         allocate(this%hwCa(0:this%hwEvaPa, 0:this%hwEvaPa));
         allocate(this%hwDa(0:this%hwEvaPa, 0:this%hwEvaPa));
         allocate(Ematr(0:this%hwEvaPa, 0:this%hwEvaPa));
         do i=0,this%hwEvaPa
            do j=0,this%hwEvaPa
               this%hwLa(i,j)=0.0d0;
               this%hwMa(i,j)=0.0d0;
               this%hwAa(i,j)=0.0d0;
               this%hwBa(i,j)=0.0d0;
               this%hwCa(i,j)=0.0d0;
               this%hwDa(i,j)=0.0d0;
               Ematr(i,j)=0.0d0;
               if (i==j) then
                  Ematr(i,j)=1.0d0;
               endif
            enddo
         enddo
      endif
      

      ! Filling L and M arrays
      if (this%hwEvaPa>0) then
         this%hwLa(1,0) = 2.0d0*this%hwACa(1)*(1.0d0-this%hwACa(0)**2);
         this%hwLa(1,1) = this%hwACa(0)*(1.0d0+2.0d0*this%hwACa(0)*this%hwACa(1)+this%hwACa(1)**2);
         this%hwLa(1,2) = this%hwACa(0)*(1.0d0-this%hwACa(1)**2);
         this%hwMa(1,0) = 2.0d0*this%hwACa(1);
         this%hwMa(1,1) = this%hwACa(0);
         this%hwMa(1,2) = this%hwACa(0);      
         if (this%hwEvaPa>1) then
            do i=2,this%hwEvaPa
               this%hwLa(i,i-1) = this%hwACa(i)*(1.0d0-this%hwACa(i-1)**2);
               this%hwLa(i,i) = this%hwACa(i)*(1.0d0+this%hwACa(i-1)**2)+this%hwACa(i-1)*(1.0d0+this%hwACa(i)**2);
               this%hwMa(i,i-1) = this%hwACa(i);
               this%hwMa(i,i) = this%hwACa(i)+this%hwACa(i-1);
               if (this%hwEvaEa>0) then
                  this%hwMa(i,i+1) = this%hwACa(i-1);
                  this%hwLa(i,i+1) = this%hwACa(i-1)*(1.0d0-this%hwACa(i)**2);
               else
                  if (i/=this%hwEvaPa) then
                     this%hwMa(i,i+1) = this%hwACa(i-1);
                     this%hwLa(i,i+1) = this%hwACa(i-1)*(1.0d0-this%hwACa(i)**2);
                  endif
               endif
            enddo
         endif
      endif
      

      ! Filling N, S, V arrays
      if (this%hwEvaEa>0) then
         this%hwLa(this%hwEvaPa+1,this%hwEvaPa) = (1.0d0-this%hwACa(this%hwEvaPa)**2)/this%hwACa(this%hwEvaPa);
         this%hwLa(this%hwEvaPa+1,this%hwEvaPa+1) = (1.0d0+this%hwACa(this%hwEvaPa)**2)/this%hwACa(this%hwEvaPa);
         this%hwMa(this%hwEvaPa+1,this%hwEvaPa) = 1.0d0/this%hwACa(this%hwEvaPa);
         this%hwMa(this%hwEvaPa+1,this%hwEvaPa+1) = 1.0d0/this%hwACa(this%hwEvaPa);
         this%hwNa(1,1) = 1.0d0/(this%hwACa(this%hwEvaPa)*this%hwSigmasA(1));
         this%hwNa(1,2) = 1.0d0/(this%hwACa(this%hwEvaPa)*this%hwSigmasA(1));       
         this%hwVa(1,1) = 1.0d0/(this%hwACa(this%hwEvaPa)*this%hwSigmasA(1));
         this%hwVa(1,2) = 1.0d0/(this%hwACa(this%hwEvaPa)*this%hwSigmasA(1));
         this%hwSa(1,1) = -this%hwSigmasA(1)/this%hwACa(this%hwEvaPa);
         this%hwSa(1,2) = this%hwSigmasA(1)/this%hwACa(this%hwEvaPa);
         if (this%hwEvaEa>1) then
            do i=2,this%hwEvaEa
               this%hwNa(i,i-1) = 1.0d0/this%hwSigmasA(i-1);
               this%hwNa(i,i) = 1.0d0/this%hwSigmasA(i-1) + 1.0d0/this%hwSigmasA(i);
               this%hwVa(i,i-1) = 1.0d0/this%hwSigmasA(i-1);
               this%hwVa(i,i) = 1.0d0/this%hwSigmasA(i-1) + 1.0d0/this%hwSigmasA(i);
               this%hwSa(i,i-1) = this%hwSigmasA(i-1);
               this%hwSa(i,i) = -(this%hwSigmasA(i-1) + this%hwSigmasA(i));
               if (i/=this%hwEvaEa) then
                  this%hwNa(i,i+1) = 1.0d0/this%hwSigmasA(i);
                  this%hwVa(i,i+1) = 1.0d0/this%hwSigmasA(i);
                  this%hwSa(i,i+1) = this%hwSigmasA(i);
               endif
            enddo
         endif
      endif

      ! Fill A, B, C, D with L and M
      if (this%hwEvaEa>0) then
         do i=0,this%hwEvaPa+1
            do j=0,this%hwEvaPa+1
               this%hwAa(i,j) = this%hwLa(i,j);
               this%hwBa(i,j) = -2.0d0*this%hwLa(i,j);
               this%hwCa(i,j) = this%hwLa(i,j);
               this%hwDa(i,j) = (cc**2)*(this%ht**2)*this%hwMa(i,j);
            enddo
         enddo
      else
         do i=0,this%hwEvaPa
            do j=0,this%hwEvaPa
               this%hwAa(i,j) = this%hwLa(i,j);
               this%hwBa(i,j) = -2.0d0*this%hwLa(i,j);
               this%hwCa(i,j) = this%hwLa(i,j);
               this%hwDa(i,j) = (cc**2)*(this%ht**2)*this%hwMa(i,j);
            enddo
         enddo 
      endif
         
      ! Fill A, B, C, D with N, V, S
      if (this%hwEvaEa>0) then
         do i=1,this%hwEvaEa
            do j=1,this%hwEvaEa
               this%hwAa(i+this%hwEvaPa+1,j+this%hwEvaPa) = this%hwNa(i,j);
               this%hwBa(i+this%hwEvaPa+1,j+this%hwEvaPa) = -2.0d0*this%hwNa(i,j)-(cc**2)*(this%ht**2)*this%hwSa(i,j);
               this%hwCa(i+this%hwEvaPa+1,j+this%hwEvaPa) = this%hwNa(i,j);
               this%hwDa(i+this%hwEvaPa+1,j+this%hwEvaPa) = (cc**2)*(this%ht**2)*this%hwVa(i,j);
            enddo
         enddo   
      endif

      ! More koefficients for sean variable
!      this%hwA(0, 0) = 3.0d0*this%hwAC(0);
!      this%hwA(0, 1) = -3.0d0*this%hwAC(0);
!      this%hwB(0, 0) = -4.0d0*this%hwAC(0);
!      this%hwB(0, 1) = 4.0d0*this%hwAC(0);
!      this%hwC(0, 0) = this%hwAC(0);
!      this%hwC(0, 1) = -this%hwAC(0);

      select case(this%dtmode)
      case(1)
          this%hwAa(0,0)=1.0d0*this%hwACa(0);
          this%hwBa(0,0)=-1.0d0*this%hwACa(0);
          this%hwCa(0,0)=0.0d0;
          if (this%hwEvaPa>0) then
             this%hwAa(0,1)=-1.0d0;
             this%hwBa(0,1)=1.0d0;
             this%hwCa(0,1)=0.0d0;
          endif
      case(2)
          this%hwAa(0, 0) = 3.0d0*this%hwACa(0);
          this%hwBa(0, 0) = -4.0d0*this%hwACa(0);
          this%hwCa(0, 0) = this%hwACa(0);
          if (this%hwEvaPa>0) then
             this%hwAa(0, 1) = -3.0d0*this%hwACa(0);
             this%hwBa(0, 1) = 4.0d0*this%hwACa(0);
             this%hwCa(0, 1) = -this%hwACa(0);
          endif
      end select

      if (this%hwEvaEa>0) then
         this%hwAa(this%hwEvaPa+1, this%hwEvaPa+this%hwEvaEa+1) = this%hwACa(this%hwEvaPa);
         this%hwAa(this%hwEvaPa+2, this%hwEvaPa+this%hwEvaEa+1) = -3.0d0/2.0d0*cc*this%ht;
         this%hwBa(this%hwEvaPa+1, this%hwEvaPa+this%hwEvaEa+1) = -2.0d0*this%hwACa(this%hwEvaPa);
         this%hwBa(this%hwEvaPa+2, this%hwEvaPa+this%hwEvaEa+1) = 2.0d0*cc*this%ht;
         this%hwCa(this%hwEvaPa+1, this%hwEvaPa+this%hwEvaEa+1) = this%hwACa(this%hwEvaPa);
         this%hwCa(this%hwEvaPa+2, this%hwEvaPa+this%hwEvaEa+1) = -1.0d0/2.0d0*cc*this%ht;
      endif   

      write(601,*) 'Angles';
      do i=0,this%hwEvaPa
         WRITE(601,520) (this%hwAngA(i))
      enddo

      write(601,*) 'A koeffs';
      do i=0,this%hwEvaPa
         WRITE(601,520) (this%hwACa(i))
      enddo

      if (this%hwEvaEa>0) then
         write(601,*) 'La array';
         do i=0,this%hwEvaPa+1
            do j=0,this%hwEvaPa+1
               WRITE(601, 520, advance="no") this%hwLa(i,j)
            enddo
            write(601,*) '';
         enddo
         write(601,*) 'Ma array';
         do i=0,this%hwEvaPa+1
            do j=0,this%hwEvaPa+1
               WRITE(601, 520, advance="no") this%hwMa(i,j)
            enddo
            write(601,*) '';
         enddo
         write(601,*) 'Na array';
         do i=1,this%hwEvaEa
            do j=1,this%hwEvaEa
               WRITE(601, 520, advance="no") this%hwNa(i,j)
            enddo
            write(601,*) '';
         enddo
         write(601,*) 'Va array';
         do i=1,this%hwEvaEa
            do j=1,this%hwEvaEa
               WRITE(601, 520, advance="no") this%hwVa(i,j)
            enddo
            write(601,*) '';
         enddo
         write(601,*) 'Sa array';
         do i=1,this%hwEvaEa
            do j=1,this%hwEvaEa
               WRITE(601, 520, advance="no") this%hwSa(i,j)
            enddo
            write(601,*) '';
         enddo
         write(601,*) 'Aa array';
         do i=0,this%hwEvaPa+this%hwEvaEa+1
            do j=0,this%hwEvaPa+this%hwEvaEa+1
               WRITE(601, 520, advance="no") this%hwAa(i,j)
            enddo
            write(601,*) '';
         enddo
         write(601,*) 'Ba array';
         do i=0,this%hwEvaPa+this%hwEvaEa+1
            do j=0,this%hwEvaPa+this%hwEvaEa+1
               WRITE(601, 520, advance="no") this%hwBa(i,j)
            enddo
            write(601,*) '';
         enddo
         write(601,*) 'Ca array';
         do i=0,this%hwEvaPa+this%hwEvaEa+1
            do j=0,this%hwEvaPa+this%hwEvaEa+1
               WRITE(601, 520, advance="no") this%hwCa(i,j)
            enddo
            write(601,*) '';
         enddo
         write(601,*) 'Da array';
         do i=0,this%hwEvaPa+this%hwEvaEa+1
            do j=0,this%hwEvaPa+this%hwEvaEa+1
               WRITE(601, 520, advance="no") this%hwDa(i,j)
            enddo
            write(601,*) '';
         enddo
      else
         write(601,*) 'La array';
         do i=0,this%hwEvaPa
            do j=0,this%hwEvaPa
               WRITE(601, 520, advance="no") this%hwLa(i,j)
            enddo
            write(601,*) '';
         enddo
         write(601,*) 'Ma array';
         do i=0,this%hwEvaPa
            do j=0,this%hwEvaPa
               WRITE(601, 520, advance="no") this%hwMa(i,j)
            enddo
            write(601,*) '';
         enddo
         write(601,*) 'Aa array';
         do i=0,this%hwEvaPa
            do j=0,this%hwEvaPa
               WRITE(601, 520, advance="no") this%hwAa(i,j)
            enddo
            write(601,*) '';
         enddo
         write(601,*) 'Ba array';
         do i=0,this%hwEvaPa
            do j=0,this%hwEvaPa
               WRITE(601, 520, advance="no") this%hwBa(i,j)
            enddo
            write(601,*) '';
         enddo
         write(601,*) 'Ca array';
         do i=0,this%hwEvaPa
            do j=0,this%hwEvaPa
               WRITE(601, 520, advance="no") this%hwCa(i,j)
            enddo
            write(601,*) '';
         enddo
         write(601,*) 'Da array';
         do i=0,this%hwEvaPa
            do j=0,this%hwEvaPa
               WRITE(601, 520, advance="no") this%hwDa(i,j)
            enddo
            write(601,*) '';
         enddo
      endif
      
520   format (E11.4, ' ');      

      if (this%hwEvaEa>0) then
         allocate(ipiv(0:this%hwEvaPa+this%hwEvaEa+1));
         call dgesv(this%hwEvaPa+this%hwEvaEa+2,this%hwEvaPa+this%hwEvaEa+2,this%hwAa,this%hwEvaPa+this%hwEvaEa+2,ipiv,Ematr,this%hwEvaPa+this%hwEvaEa+2,info);
         do i=0,this%hwEvaPa+this%hwEvaEa+1 
            do j=0,this%hwEvaPa+this%hwEvaEa+1
               this%hwAa(i,j)=Ematr(i,j);
            enddo
         enddo
         write(601,*) 'Aa^-1 array';
         do i=0,this%hwEvaPa+this%hwEvaEa+1
            do j=0,this%hwEvaPa+this%hwEvaEa+1
               WRITE(601, 520, advance="no") this%hwAa(i,j)
            enddo
            write(601,*) '';
         enddo
         allocate(this%hwFAnew(0:this%hwEvaPa+this%hwEvaEa+1 ), this%hwFA(0:this%hwEvaPa+this%hwEvaEa+1 ), this%hwFAold(0:this%hwEvaPa+this%hwEvaEa+1 ), this%hwGFA(0:this%hwEvaPa+this%hwEvaEa+1 ));
         do i=0,this%hwEvaPa+this%hwEvaEa+1
            this%hwFAnew(i)=0.0d0;
            this%hwFA(i)=0.0d0;
            this%hwFAold(i)=0.0d0;
            this%hwGFA(i)=0.0d0;
         enddo
      else
         allocate(ipiv(0:this%hwEvaPa));
         call dgesv(this%hwEvaPa+1,this%hwEvaPa+1,this%hwAa,this%hwEvaPa+1,ipiv,Ematr,this%hwEvaPa+1,info);
         do i=0,this%hwEvaPa
            do j=0,this%hwEvaPa
               this%hwAa(i,j)=Ematr(i,j);
            enddo
         enddo
         write(601,*) 'Aa^-1 array';
         do i=0,this%hwEvaPa
            do j=0,this%hwEvaPa
               WRITE(601, 520, advance="no") this%hwAa(i,j)
            enddo
            write(601,*) '';
         enddo
         allocate(this%hwFAnew(0:this%hwEvaPa), this%hwFA(0:this%hwEvaPa), this%hwFAold(0:this%hwEvaPa), this%hwGFA(0:this%hwEvaPa));
         do i=0,this%hwEvaPa
            this%hwFAnew(i)=0.0d0;
            this%hwFA(i)=0.0d0;
            this%hwFAold(i)=0.0d0;
            this%hwGFA(i)=0.0d0;
         enddo
      endif
      deallocate(ipiv);
      deallocate(Ematr);
    end subroutine initHWEvaEdgeAuxVariables   


!destroyHWEvaAuxVariables---------------------------
    subroutine destroyHWEvaAuxVariables(this)
    ! Destroy arrays for Givoli-Neta ABC aux variables
      class(tProblem) :: this;
      integer :: i;
      if (this%hwEvaE>0) then
         do i=0,this%hwEvaP+this%hwEvaE+1
            call this%hwX(i)%auxmesh_destroy;
            call this%hwY(i)%auxmesh_destroy;
            call this%hwZ(i)%auxmesh_destroy;
            call this%hwXOld(i)%auxmesh_destroy;
            call this%hwYOld(i)%auxmesh_destroy;
            call this%hwZOld(i)%auxmesh_destroy;
            call this%hwXOld2(i)%auxmesh_destroy;
            call this%hwYOld2(i)%auxmesh_destroy;
            call this%hwZOld2(i)%auxmesh_destroy;
         end do
      else
         do i=0,this%hwEvaP
            call this%hwX(i)%auxmesh_destroy;
            call this%hwY(i)%auxmesh_destroy;
            call this%hwZ(i)%auxmesh_destroy;
            call this%hwXOld(i)%auxmesh_destroy;
            call this%hwYOld(i)%auxmesh_destroy;
            call this%hwZOld(i)%auxmesh_destroy;
            call this%hwXOld2(i)%auxmesh_destroy;
            call this%hwYOld2(i)%auxmesh_destroy;
            call this%hwZOld2(i)%auxmesh_destroy;
         end do
      endif
      deallocate(this%hwX, this%hwY, this%hwZ, this%hwXOld, this%hwYOld, this%hwZOld, this%hwXOld2, this%hwYOld2, this%hwZOld2);
      deallocate(this%hwAng);
      deallocate(this%hwAC);
      deallocate(this%hwL);
      deallocate(this%hwM);
      deallocate(this%hwA);
      deallocate(this%hwB);
      deallocate(this%hwC);
      deallocate(this%hwD);
      deallocate(this%hwFnew, this%hwF, this%hwFold, this%hwGF);
      if (this%hwEvaE>0) then
         deallocate(this%hwSigmas);
         deallocate(this%hwN);
         deallocate(this%hwS);
         deallocate(this%hwV);
      endif
      
      !call this%destroyHWEvaEdgeAuxVariables();
    end subroutine destroyHWEvaAuxVariables


!destroyHWEvaEdgeAuxVariables---------------------------
    subroutine destroyHWEvaEdgeAuxVariables(this)
    ! Destroy arrays for Givoli-Neta ABC Edge aux variables
      class(tProblem) :: this;
      integer :: i,j;
      integer :: im, ia;
      if (this%hwEvaE>0) then
         im = this%hwEvaP+this%hwEvaE+1;
      else
         im = this%hwEvaP;
      endif

      if (this%hwEvaEa>0) then
         ia = this%hwEvaPa+this%hwEvaEa+1;
      else
         ia = this%hwEvaPa;
      endif 
      do i=0,im
         do j=0,ia
            call this%hwYEdgeNX(i,j)%EdgeAuxMesh_destroy;
            call this%hwYEdgeNXold(i,j)%EdgeAuxMesh_destroy;
            call this%hwYEdgeNXold2(i,j)%EdgeAuxMesh_destroy;
         enddo
      enddo
      deallocate(this%hwYEdgeNx, this%hwYEdgeNxOld, this%hwYEdgeNxOld2);
      deallocate(this%hwAngA);
      deallocate(this%hwSigmasA);
      deallocate(this%hwACa);
      deallocate(this%hwLa);
      deallocate(this%hwAa);
      deallocate(this%hwBa);
      deallocate(this%hwCa);
      deallocate(this%hwDa);
      deallocate(this%hwNa);
      deallocate(this%hwSa);
      deallocate(this%hwVa);
      deallocate(this%hwMa);
      deallocate(this%hwFAnew, this%hwFA, this%hwFAold, this%hwGFA);
    end subroutine destroyHWEvaEdgeAuxVariables    


!saveHWEvaAuxVariables---------------------------
    subroutine saveHWEvaAuxVariables(this)
    ! Save arrays for HW ABC aux variables
      class(tProblem) :: this;
      integer :: i, j, k, v, im;
      if (this%hwEvaE>0) then
         im=this%hwEvaP+this%hwEvaE+1
      else
         im=this%hwEvaP
      endif   
      do v=0,im
         ! HWx field filling
         !$OMP PARALLEL DO SCHEDULE(GUIDED)
         do k=0,this%Nz
            do j=0,this%Ny
               do i=1,this%Nx
                  this%HWxOld2(v)%Y(i,j,k)=this%HWxOld(v)%Y(i,j,k);
                  this%HWxOld(v)%Y(i,j,k)=this%HWx(v)%Y(i,j,k);
                  this%HWxOld2(v)%Z(i,j,k)=this%HWxOld(v)%Z(i,j,k);
                  this%HWxOld(v)%Z(i,j,k)=this%HWx(v)%Z(i,j,k);
               enddo
            enddo
         enddo
         !$OMP END PARALLEL DO
         ! HWy field filling
         !$OMP PARALLEL DO SCHEDULE(GUIDED)
         do k=0,this%Nz
            do j=1,this%Ny
               do i=0,this%Nx
                  this%HWyOld2(v)%X(i,j,k)=this%HWyOld(v)%X(i,j,k);
                  this%HWyOld(v)%X(i,j,k)=this%HWy(v)%X(i,j,k);
                  this%HWyOld2(v)%Z(i,j,k)=this%HWyOld(v)%Z(i,j,k);
                  this%HWyOld(v)%Z(i,j,k)=this%HWy(v)%Z(i,j,k);
               enddo
            enddo
         enddo
         !$OMP END PARALLEL DO
         ! HWz field filling
         !$OMP PARALLEL DO SCHEDULE(GUIDED)
         do k=1,this%Nz
            do j=0,this%Ny
               do i=0,this%Nx
                  this%HWzOld2(v)%X(i,j,k)=this%HWzOld(v)%X(i,j,k);
                  this%HWzOld(v)%X(i,j,k)=this%HWz(v)%X(i,j,k);
                  this%HWzOld2(v)%Y(i,j,k)=this%HWzOld(v)%Y(i,j,k);
                  this%HWzOld(v)%Y(i,j,k)=this%HWz(v)%Y(i,j,k);
               enddo
            enddo
         enddo
         !$OMP END PARALLEL DO
      enddo
      !call this%saveHWEvaEdgeAuxVariables();
    end subroutine saveHWEvaAuxVariables


!saveHWEvaEdgeAuxVariables---------------------------
    subroutine saveHWEvaEdgeAuxVariables(this)
    ! Save arrays for HW ABC Edge aux variables
      class(tProblem) :: this;
      integer :: i, j, k, v, va, im, ia;
      if (this%hwEvaE>0) then
         im = this%hwEvaP+this%hwEvaE+1;
      else
         im = this%hwEvaP;
      endif

      if (this%hwEvaEa>0) then
         ia = this%hwEvaPa+this%hwEvaEa+1;
      else
         ia = this%hwEvaPa;
      endif    
      do v=0,im
         do va=0,ia
            ! HWx EdgeNX field filling
            !$OMP PARALLEL DO SCHEDULE(GUIDED)
            do j=1,this%Ny
               this%hwYEdgeNXold2(v,va)%fbottom(j)=this%hwYEdgeNXold(v,va)%fbottom(j);
               this%hwYEdgeNXold(v,va)%fbottom(j)=this%hwYEdgeNX(v,va)%fbottom(j);
               this%hwYEdgeNXold2(v,va)%fup(j)=this%hwYEdgeNXold(v,va)%fup(j);
               this%hwYEdgeNXold(v,va)%fup(j)=this%hwYEdgeNX(v,va)%fup(j);
            enddo
            !$OMP END PARALLEL DO
            !$OMP PARALLEL DO SCHEDULE(GUIDED)
            do k=0,this%Nz
               this%hwYEdgeNXold2(v,va)%fleft(k)=this%hwYEdgeNXold(v,va)%fleft(k);
               this%hwYEdgeNXold(v,va)%fleft(k)=this%hwYEdgeNX(v,va)%fleft(k);
               this%hwYEdgeNXold2(v,va)%fright(k)=this%hwYEdgeNXold(v,va)%fright(k);
               this%hwYEdgeNXold(v,va)%fright(k)=this%hwYEdgeNX(v,va)%fright(k);
            enddo
            !$OMP END PARALLEL DO        
      enddo
      enddo
    end subroutine saveHWEvaEdgeAuxVariables    


!fillEBoundaryABCHWEva------------------------------------------------------------------
  subroutine  fillEBoundaryABCHWEva(this,t)
    !Fill E boundary with ABC Hagstrom-Warburton Ey_x=0, exact elsewhere
    class(tProblem) :: this
    integer, intent(in) :: t
    integer :: i ,j ,k, v, it  
    !================ EX =================
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do k=0,this%Nz
       do i=1,this%Nx
          this%Ef%X(i, 0, k) = this%Source%getEpoint(this%xi05(i), this%yi(0),  this%zi(k), this%Ti(t), i, 0, k, 1);
          this%Ef%X(i, this%Ny, k) = this%Source%getEpoint(this%xi05(i), this%yi(this%Ny),  this%zi(k), this%Ti(t), i, this%Ny, k, 1);
       enddo
    enddo
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do j=1,this%Ny-1
       do i=1,this%Nx
          this%Ef%X(i, j, 0) = this%Source%getEpoint(this%xi05(i), this%yi(j),  this%zi(0), this%Ti(t), i, j, 0, 1);
          this%Ef%X(i, j, this%Nz) = this%Source%getEpoint(this%xi05(i), this%yi(j),  this%zi(this%Nz), this%Ti(t), i, j, this%Nz, 1);
       enddo
    enddo
    !$OMP END PARALLEL DO
    !================ EY =================
    ! Fill EY Boundaries
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do j=1,this%Ny
       do i=0,this%Nx-1
          this%Ef%Y(i, j, 0) = this%Source%getEpoint(this%xi(i), this%yi05(j),  this%zi(0), this%Ti(t), i, j, 0, 2);
          this%Ef%Y(i, j, this%Nz) = this%Source%getEpoint(this%xi(i), this%yi05(j),  this%zi(this%Nz), this%Ti(t), i, j, this%Nz, 2);
       enddo
    enddo
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do k=0,this%Nz
       do j=1,this%Ny
          this%Ef%Y(0, j, k) = this%Source%getEpoint(this%xi(0), this%yi05(j),  this%zi(k), this%Ti(t), 0, j, k, 2);
!          this%Ef%Y(this%Nx, j, k) = this%Source%getEpoint(this%xi(this%Nx), this%yi05(j),  this%zi(k), this%Ti(t), 2);
       enddo
    enddo
    !$OMP END PARALLEL DO
    !!!*****************Ey_x=0 update*****************
    call this%HWEvaEyHigdonUpdate(t);
    !================ EZ =================
    !$OMP PARALLEL DO SCHEDULE(GUIDED)   
    do k=1,this%Nz
       do j=0,this%Ny
          this%Ef%Z(0, j, k) = this%Source%getEpoint(this%xi(0), this%yi(j),  this%zi05(k), this%Ti(t), 0, j, k, 3);
          this%Ef%Z(this%Nx, j, k) = this%Source%getEpoint(this%xi(this%Nx), this%yi(j),  this%zi05(k), this%Ti(t), this%Nx, j, k, 3);
       enddo
    enddo
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do k=1,this%Nz
       do i=1,this%Nx-1
          this%Ef%Z(i, 0, k) = this%Source%getEpoint(this%xi(i), this%yi(0),  this%zi05(k), this%Ti(t), i, 0, k, 3);
          this%Ef%Z(i, this%Ny, k) = this%Source%getEpoint(this%xi(i), this%yi(this%Ny),  this%zi05(k), this%Ti(t), i, this%Ny, k, 3);
       enddo
    enddo
    !$OMP END PARALLEL DO  
  end subroutine fillEBoundaryABCHWEva

  
!HWEvaEyUpdate------------------------------------------------------------------
  subroutine HWEvaEyUpdate(this)
    !Fill Ey_x=0 boundary with ABC Hagstrom-Warburton
    class(tProblem) :: this
    integer :: j ,k, v, w, im
    real*8 :: de, df1, df2, tp;
    character*8 :: trans;
    real*8, allocatable :: temp(:);
    
    if (this%hwEvaE>0) then
       im=this%hwEvaP+this%hwEvaE+1
    else
       im=this%hwEvaP
    endif
    allocate(temp(0:im));
    !=========================== EY == X BOUNDARIES ============================
    !EYFace X=0 -----
    do k=0,this%Nz
       do j=1,this%Ny
          select case(this%dxmode)
          case(1)
             de = (this%EfOld%Y(this%Nx, j, k) - this%EfOld%Y(this%Nx-1, j, k))/(this%hhx);
          case(2)
             de = (3.0d0*this%EfOld%Y(this%Nx, j, k) - 4.0d0*this%EfOld%Y(this%Nx-1, j, k) + this%EfOld%Y(this%Nx-2, j, k) )/(2.0d0*this%hhx);
          end select
          do v=0,im                  ! Vectors formation
             this%hwF(v) = this%hwYOld(v)%X(this%Nx, j, k);
             this%hwFold(v) = this%hwYold2(v)%X(this%Nx, j, k);
             select case(this%dxauxmode)
             case(1) !===1st order===
               if (j==1.AND.k==0) then                             ! j=1, k=0
                  df1 = (this%hwYOld(v)%X(this%Nx,1,0)-2.0d0*this%hwYOld(v)%X(this%Nx,2,0)+this%hwYOld(v)%X(this%Nx,3,0))/(this%hhy**2);   ! forward one-sided in Y
                  df2 = (this%hwYOld(v)%X(this%Nx,1,0)-2.0d0*this%hwYOld(v)%X(this%Nx,1,1)+this%hwYOld(v)%X(this%Nx,1,2))/(this%hhz**2);   ! forward one-sided in Z
               elseif (j==this%Ny.AND.k==0) then                   ! j=Ny, k=0
                  df1 = (this%hwYOld(v)%X(this%Nx,this%Ny,0)-2.0d0*this%hwYOld(v)%X(this%Nx,this%Ny-1,0)+this%hwYOld(v)%X(this%Nx,this%Ny-2,0))/(this%hhy**2);     ! backward one-sided in Y
                  df2 = (this%hwYOld(v)%X(this%Nx,this%Ny,0)-2.0d0*this%hwYOld(v)%X(this%Nx,this%Ny,1)+this%hwYOld(v)%X(this%Nx,this%Ny,2))/(this%hhz**2);         ! forward one-sided in Z 
               elseif (j==1.AND.k==this%Nz) then                   ! j=1, k=Nz
                  df1 = (this%hwYOld(v)%X(this%Nx,1,this%Nz)-2.0d0*this%hwYOld(v)%X(this%Nx,2,this%Nz)+this%hwYOld(v)%X(this%Nx,3,this%Nz))/(this%hhy**2);         ! forward one-sided in Y
                  df2 = (this%hwYOld(v)%X(this%Nx,1,this%Nz)-2.0d0*this%hwYOld(v)%X(this%Nx,1,this%Nz-1)+this%hwYOld(v)%X(this%Nx,1,this%Nz-2))/(this%hhz**2);     ! backward one-sided in Z 
               elseif (j==this%Ny.AND.k==this%Nz) then             ! j=Ny, k=Nz
                  df1 = (this%hwYOld(v)%X(this%Nx,this%Ny,this%Nz)-2.0d0*this%hwYOld(v)%X(this%Nx,this%Ny-1,this%Nz)+this%hwYOld(v)%X(this%Nx,this%Ny-2,this%Nz))/(this%hhy**2);  ! backward one-sided in Y
                  df2 = (this%hwYOld(v)%X(this%Nx,this%Ny,this%Nz)-2.0d0*this%hwYOld(v)%X(this%Nx,this%Ny,this%Nz-1)+this%hwYOld(v)%X(this%Nx,this%Ny,this%Nz-2))/(this%hhz**2);  ! backward one-sided in Z 
               elseif (j>1.AND.j<this%Ny.AND.k==0) then            ! j>1, j<Ny, k=0
                  df1 = (this%hwYOld(v)%X(this%Nx,j+1,0)-2.0d0*this%hwYOld(v)%X(this%Nx,j,0)+this%hwYOld(v)%X(this%Nx,j-1,0))/(this%hhy**2);                                      ! central in Y
                  df2 = (this%hwYOld(v)%X(this%Nx,j,0)-2.0d0*this%hwYOld(v)%X(this%Nx,j,1)+this%hwYOld(v)%X(this%Nx,j,2))/(this%hhz**2);                                          ! forward one-sided in Z      
               elseif (j>1.AND.j<this%Ny.AND.k==this%Nz) then      ! j>1, j<Ny, k=Nz
                  df1 = (this%hwYOld(v)%X(this%Nx,j+1,this%Nz)-2.0d0*this%hwYOld(v)%X(this%Nx,j,this%Nz)+this%hwYOld(v)%X(this%Nx,j-1,this%Nz))/(this%hhy**2);                    ! central in Y
                  df2 = (this%hwYOld(v)%X(this%Nx,j,this%Nz)-2.0d0*this%hwYOld(v)%X(this%Nx,j,this%Nz-1)+this%hwYOld(v)%X(this%Nx,j,this%Nz-2))/(this%hhz**2);                    ! backward one-sided Z 
               elseif (j==1.AND.k>0.AND.k<this%Nz) then            ! j=1, k>0, k<Nz
                  df1 = (this%hwYOld(v)%X(this%Nx,1,k)-2.0d0*this%hwYOld(v)%X(this%Nx,2,k)+this%hwYOld(v)%X(this%Nx,3,k))/(this%hhy**2);                                           ! forward one-sided in Y
                  df2 = (this%hwYOld(v)%X(this%Nx,1,k+1)-2.0d0*this%hwYOld(v)%X(this%Nx,1,k)+this%hwYOld(v)%X(this%Nx,1,k-1))/(this%hhz**2);                                       ! central in Z 
               elseif (j==this%Ny.AND.k>0.AND.k<this%Nz) then      ! j=Ny, k>0, k<Nz
                  df1 = (this%hwYOld(v)%X(this%Nx,this%Ny,k)-2.0d0*this%hwYOld(v)%X(this%Nx,this%Ny-1,k)+this%hwYOld(v)%X(this%Nx,this%Ny-2,k))/(this%hhy**2);    ! backward one-sided in Y
                  df2 = (this%hwYOld(v)%X(this%Nx,this%Ny,k+1)-2.0d0*this%hwYOld(v)%X(this%Nx,this%Ny,k)+this%hwYOld(v)%X(this%Nx,this%Ny,k-1))/(this%hhz**2);    ! central in Z 
               else
                  df1 = (this%hwYOld(v)%X(this%Nx,j+1,k)-2.0d0*this%hwYOld(v)%X(this%Nx,j,k)+this%hwYOld(v)%X(this%Nx,j-1,k))/(this%hhy**2);    ! central in Y
                  df2 = (this%hwYOld(v)%X(this%Nx,j,k+1)-2.0d0*this%hwYOld(v)%X(this%Nx,j,k)+this%hwYOld(v)%X(this%Nx,j,k-1))/(this%hhz**2);    ! central in Z 
               endif
             case(2) !===2nd order===
                if (j==1.AND.k==0) then                             ! j=1, k=0
                   df1 = (2.0d0*this%hwYOld(v)%X(this%Nx,1,0)-5.0d0*this%hwYOld(v)%X(this%Nx,2,0)+4.0d0*this%hwYOld(v)%X(this%Nx,3,0)-this%hwYOld(v)%X(this%Nx,4,0))/(this%hhy**2);   ! forward one-sided in Y
                   df2 = (2.0d0*this%hwYOld(v)%X(this%Nx,1,0)-5.0d0*this%hwYOld(v)%X(this%Nx,1,1)+4.0d0*this%hwYOld(v)%X(this%Nx,1,2)-this%hwYOld(v)%X(this%Nx,1,3))/(this%hhz**2);   ! forward one-sided in Z
                elseif (j==this%Ny.AND.k==0) then                   ! j=Ny, k=0
                   df1 = (2.0d0*this%hwYOld(v)%X(this%Nx,this%Ny,0)-5.0d0*this%hwYOld(v)%X(this%Nx,this%Ny-1,0)+4.0d0*this%hwYOld(v)%X(this%Nx,this%Ny-2,0)-this%hwYOld(v)%X(this%Nx,this%Ny-3,0))/(this%hhy**2);    ! backward one-sided in Y
                   df2 = (2.0d0*this%hwYOld(v)%X(this%Nx,this%Ny,0)-5.0d0*this%hwYOld(v)%X(this%Nx,this%Ny,1)+4.0d0*this%hwYOld(v)%X(this%Nx,this%Ny,2)-this%hwYOld(v)%X(this%Nx,this%Ny,3))/(this%hhz**2);          ! forward one-sided in Z 
                elseif (j==1.AND.k==this%Nz) then                   ! j=1, k=Nz
                   df1 = (2.0d0*this%hwYOld(v)%X(this%Nx,1,this%Nz)-5.0d0*this%hwYOld(v)%X(this%Nx,2,this%Nz)+4.0d0*this%hwYOld(v)%X(this%Nx,3,this%Nz)-this%hwYOld(v)%X(this%Nx,4,this%Nz))/(this%hhy**2);          ! forward one-sided in Y
                   df2 = (2.0d0*this%hwYOld(v)%X(this%Nx,1,this%Nz)-5.0d0*this%hwYOld(v)%X(this%Nx,1,this%Nz-1)+4.0d0*this%hwYOld(v)%X(this%Nx,1,this%Nz-2)-this%hwYOld(v)%X(this%Nx,1,this%Nz-3))/(this%hhz**2);    ! backward one-sided in Z 
                elseif (j==this%Ny.AND.k==this%Nz) then             ! j=Ny, k=Nz
                   df1 = (2.0d0*this%hwYOld(v)%X(this%Nx,this%Ny,this%Nz)-5.0d0*this%hwYOld(v)%X(this%Nx,this%Ny-1,this%Nz)+4.0d0*this%hwYOld(v)%X(this%Nx,this%Ny-2,this%Nz)-this%hwYOld(v)%X(this%Nx,this%Ny-3,this%Nz))/(this%hhy**2);  ! backward one-sided in Y
                   df2 = (2.0d0*this%hwYOld(v)%X(this%Nx,this%Ny,this%Nz)-5.0d0*this%hwYOld(v)%X(this%Nx,this%Ny,this%Nz-1)+4.0d0*this%hwYOld(v)%X(this%Nx,this%Ny,this%Nz-2)-this%hwYOld(v)%X(this%Nx,this%Ny,this%Nz-3))/(this%hhz**2);  ! backward one-sided in Z 
                elseif (j>1.AND.j<this%Ny.AND.k==0) then            ! j>1, j<Ny, k=0
                   df1 = (this%hwYOld(v)%X(this%Nx,j+1,0)-2.0d0*this%hwYOld(v)%X(this%Nx,j,0)+this%hwYOld(v)%X(this%Nx,j-1,0))/(this%hhy**2);                                                                  ! central in Y
                   df2 = (2.0d0*this%hwYOld(v)%X(this%Nx,j,0)-5.0d0*this%hwYOld(v)%X(this%Nx,j,1)+4.0d0*this%hwYOld(v)%X(this%Nx,j,2)-this%hwYOld(v)%X(this%Nx,j,3))/(this%hhz**2);                                  ! forward one-sided in Z      
                elseif (j>1.AND.j<this%Ny.AND.k==this%Nz) then      ! j>1, j<Ny, k=Nz
                   df1 = (this%hwYOld(v)%X(this%Nx,j+1,this%Nz)-2.0d0*this%hwYOld(v)%X(this%Nx,j,this%Nz)+this%hwYOld(v)%X(this%Nx,j-1,this%Nz))/(this%hhy**2);                                                ! central in Y
                   df2 = (2.0d0*this%hwYOld(v)%X(this%Nx,j,this%Nz)-5.0d0*this%hwYOld(v)%X(this%Nx,j,this%Nz-1)+4.0d0*this%hwYOld(v)%X(this%Nx,j,this%Nz-2)-this%hwYOld(v)%X(this%Nx,j,this%Nz-3))/(this%hhz**2);    ! backward one-sided Z 
                elseif (j==1.AND.k>0.AND.k<this%Nz) then            ! j=1, k>0, k<Nz
                   df1 = (2.0d0*this%hwYOld(v)%X(this%Nx,1,k)-5.0d0*this%hwYOld(v)%X(this%Nx,2,k)+4.0d0*this%hwYOld(v)%X(this%Nx,3,k)-this%hwYOld(v)%X(this%Nx,4,k))/(this%hhy**2);                                  ! forward one-sided in Y
                   df2 = (this%hwYOld(v)%X(this%Nx,1,k+1)-2.0d0*this%hwYOld(v)%X(this%Nx,1,k)+this%hwYOld(v)%X(this%Nx,1,k-1))/(this%hhz**2);                                                                  ! central in Z 
                elseif (j==this%Ny.AND.k>0.AND.k<this%Nz) then      ! j=Ny, k>0, k<Nz
                   df1 = (2.0d0*this%hwYOld(v)%X(this%Nx,this%Ny,k)-5.0d0*this%hwYOld(v)%X(this%Nx,this%Ny-1,k)+4.0d0*this%hwYOld(v)%X(this%Nx,this%Ny-2,k)-this%hwYOld(v)%X(this%Nx,this%Ny-3,k))/(this%hhy**2);    ! backward one-sided in Y
                   df2 = (this%hwYOld(v)%X(this%Nx,this%Ny,k+1)-2.0d0*this%hwYOld(v)%X(this%Nx,this%Ny,k)+this%hwYOld(v)%X(this%Nx,this%Ny,k-1))/(this%hhz**2);                                                ! central in Z 
                else
                   df1 = (this%hwYOld(v)%X(this%Nx,j+1,k)-2.0d0*this%hwYOld(v)%X(this%Nx,j,k)+this%hwYOld(v)%X(this%Nx,j-1,k))/(this%hhy**2);    ! central in Y
                   df2 = (this%hwYOld(v)%X(this%Nx,j,k+1)-2.0d0*this%hwYOld(v)%X(this%Nx,j,k)+this%hwYOld(v)%X(this%Nx,j,k-1))/(this%hhz**2);    ! central in Z 
                endif
             end select
             this%hwGF(v) = df1+df2;
          enddo

          do v=0,im
             temp(v) = 0.0d0;
          enddo
          ! B*Y^(n)
          do v=0,im
             tp = 0.0d0;
             do w=0,im
                tp = tp - this%hwB(v,w)*this%hwF(w);
             enddo
             temp(v) = temp(v)+tp;
          enddo
          ! C*Y^(n-1)
          do v=0,im
             tp = 0.0d0;
             do w=0,im
                tp = tp - this%hwC(v,w)*this%hwFold(w);
             enddo
             temp(v) = temp(v)+tp;
          enddo
          ! D*GY^(n)
          do v=0,im
             tp = 0.0d0;
             do w=0,im
               tp = tp + this%hwD(v,w)*this%hwGF(w);
             enddo
             temp(v) = temp(v)+tp;
          enddo
          !add dx_u
          select case(this%dtmode)
          case(1)
             temp(0) = temp(0)-1.0d0*cc*this%ht*de;
          case(2)
             temp(0) = temp(0)-2.0d0*cc*this%ht*de;
          end select
          ! A^(-1)*()
          do v=0,im
             tp = 0.0d0;
             do w=0,im
                tp=tp + this%hwA(v,w)*temp(w);
             enddo
             this%hwY(v)%X(this%Nx, j, k) = tp;
          enddo
          this%Ef%Y(this%Nx, j, k) = this%hwY(0)%X(this%Nx, j, k);
       enddo
    enddo
    deallocate(temp);
  end subroutine HWEvaEyUpdate
  

!HWEvaEySomUpdate------------------------------------------------------------------
  subroutine HWEvaEySomUpdate(this)
    !Fill Ey_x=0 boundary with ABC Hagstrom-Warburton
    class(tProblem) :: this
    integer :: j ,k, v, w, im
    real*8 :: de, df1, df2, tp, gx, gy, gz;
    character*8 :: trans;
    real*8, allocatable :: temp(:);
    gx = (1.0d0-cc*this%ht/this%hhx)/(1.0d0+cc*this%ht/this%hhx);
    gy = (1.0d0-cc*this%ht/this%hhy)/(1.0d0+cc*this%ht/this%hhy);
    gz = (1.0d0-cc*this%ht/this%hhz)/(1.0d0+cc*this%ht/this%hhz);
    if (this%hwEvaE>0) then
       im=this%hwEvaP+this%hwEvaE+1
    else
       im=this%hwEvaP
    endif
    allocate(temp(0:im));
    !=========================== EY == X BOUNDARIES ============================
    !EYFace X=0 -----
    do k=1,this%Nz-1
       do j=2,this%Ny-1
          select case(this%dxmode)
          case(1)
             de = (this%EfOld%Y(this%Nx, j, k) - this%EfOld%Y(this%Nx-1, j, k))/(this%hhx);
          case(2)
             de = (3.0d0*this%EfOld%Y(this%Nx, j, k) - 4.0d0*this%EfOld%Y(this%Nx-1, j, k) + this%EfOld%Y(this%Nx-2, j, k) )/(2.0d0*this%hhx);
          end select
          do v=0,im                  ! Vectors formation
             this%hwF(v) = this%hwYOld(v)%X(this%Nx, j, k);
             this%hwFold(v) = this%hwYold2(v)%X(this%Nx, j, k);
             select case(this%dxauxmode)
             case(1) !===1st order===
               if (j==1.AND.k==0) then                             ! j=1, k=0
                  df1 = (this%hwYOld(v)%X(this%Nx,1,0)-2.0d0*this%hwYOld(v)%X(this%Nx,2,0)+this%hwYOld(v)%X(this%Nx,3,0))/(this%hhy**2);   ! forward one-sided in Y
                  df2 = (this%hwYOld(v)%X(this%Nx,1,0)-2.0d0*this%hwYOld(v)%X(this%Nx,1,1)+this%hwYOld(v)%X(this%Nx,1,2))/(this%hhz**2);   ! forward one-sided in Z
               elseif (j==this%Ny.AND.k==0) then                   ! j=Ny, k=0
                  df1 = (this%hwYOld(v)%X(this%Nx,this%Ny,0)-2.0d0*this%hwYOld(v)%X(this%Nx,this%Ny-1,0)+this%hwYOld(v)%X(this%Nx,this%Ny-2,0))/(this%hhy**2);     ! backward one-sided in Y
                  df2 = (this%hwYOld(v)%X(this%Nx,this%Ny,0)-2.0d0*this%hwYOld(v)%X(this%Nx,this%Ny,1)+this%hwYOld(v)%X(this%Nx,this%Ny,2))/(this%hhz**2);         ! forward one-sided in Z 
               elseif (j==1.AND.k==this%Nz) then                   ! j=1, k=Nz
                  df1 = (this%hwYOld(v)%X(this%Nx,1,this%Nz)-2.0d0*this%hwYOld(v)%X(this%Nx,2,this%Nz)+this%hwYOld(v)%X(this%Nx,3,this%Nz))/(this%hhy**2);         ! forward one-sided in Y
                  df2 = (this%hwYOld(v)%X(this%Nx,1,this%Nz)-2.0d0*this%hwYOld(v)%X(this%Nx,1,this%Nz-1)+this%hwYOld(v)%X(this%Nx,1,this%Nz-2))/(this%hhz**2);     ! backward one-sided in Z 
               elseif (j==this%Ny.AND.k==this%Nz) then             ! j=Ny, k=Nz
                  df1 = (this%hwYOld(v)%X(this%Nx,this%Ny,this%Nz)-2.0d0*this%hwYOld(v)%X(this%Nx,this%Ny-1,this%Nz)+this%hwYOld(v)%X(this%Nx,this%Ny-2,this%Nz))/(this%hhy**2);  ! backward one-sided in Y
                  df2 = (this%hwYOld(v)%X(this%Nx,this%Ny,this%Nz)-2.0d0*this%hwYOld(v)%X(this%Nx,this%Ny,this%Nz-1)+this%hwYOld(v)%X(this%Nx,this%Ny,this%Nz-2))/(this%hhz**2);  ! backward one-sided in Z 
               elseif (j>1.AND.j<this%Ny.AND.k==0) then            ! j>1, j<Ny, k=0
                  df1 = (this%hwYOld(v)%X(this%Nx,j+1,0)-2.0d0*this%hwYOld(v)%X(this%Nx,j,0)+this%hwYOld(v)%X(this%Nx,j-1,0))/(this%hhy**2);                                      ! central in Y
                  df2 = (this%hwYOld(v)%X(this%Nx,j,0)-2.0d0*this%hwYOld(v)%X(this%Nx,j,1)+this%hwYOld(v)%X(this%Nx,j,2))/(this%hhz**2);                                          ! forward one-sided in Z      
               elseif (j>1.AND.j<this%Ny.AND.k==this%Nz) then      ! j>1, j<Ny, k=Nz
                  df1 = (this%hwYOld(v)%X(this%Nx,j+1,this%Nz)-2.0d0*this%hwYOld(v)%X(this%Nx,j,this%Nz)+this%hwYOld(v)%X(this%Nx,j-1,this%Nz))/(this%hhy**2);                    ! central in Y
                  df2 = (this%hwYOld(v)%X(this%Nx,j,this%Nz)-2.0d0*this%hwYOld(v)%X(this%Nx,j,this%Nz-1)+this%hwYOld(v)%X(this%Nx,j,this%Nz-2))/(this%hhz**2);                    ! backward one-sided Z 
               elseif (j==1.AND.k>0.AND.k<this%Nz) then            ! j=1, k>0, k<Nz
                  df1 = (this%hwYOld(v)%X(this%Nx,1,k)-2.0d0*this%hwYOld(v)%X(this%Nx,2,k)+this%hwYOld(v)%X(this%Nx,3,k))/(this%hhy**2);                                           ! forward one-sided in Y
                  df2 = (this%hwYOld(v)%X(this%Nx,1,k+1)-2.0d0*this%hwYOld(v)%X(this%Nx,1,k)+this%hwYOld(v)%X(this%Nx,1,k-1))/(this%hhz**2);                                       ! central in Z 
               elseif (j==this%Ny.AND.k>0.AND.k<this%Nz) then      ! j=Ny, k>0, k<Nz
                  df1 = (this%hwYOld(v)%X(this%Nx,this%Ny,k)-2.0d0*this%hwYOld(v)%X(this%Nx,this%Ny-1,k)+this%hwYOld(v)%X(this%Nx,this%Ny-2,k))/(this%hhy**2);    ! backward one-sided in Y
                  df2 = (this%hwYOld(v)%X(this%Nx,this%Ny,k+1)-2.0d0*this%hwYOld(v)%X(this%Nx,this%Ny,k)+this%hwYOld(v)%X(this%Nx,this%Ny,k-1))/(this%hhz**2);    ! central in Z 
               else
                  df1 = (this%hwYOld(v)%X(this%Nx,j+1,k)-2.0d0*this%hwYOld(v)%X(this%Nx,j,k)+this%hwYOld(v)%X(this%Nx,j-1,k))/(this%hhy**2);    ! central in Y
                  df2 = (this%hwYOld(v)%X(this%Nx,j,k+1)-2.0d0*this%hwYOld(v)%X(this%Nx,j,k)+this%hwYOld(v)%X(this%Nx,j,k-1))/(this%hhz**2);    ! central in Z 
               endif
             case(2) !===2nd order===
                if (j==1.AND.k==0) then                             ! j=1, k=0
                   df1 = (2.0d0*this%hwYOld(v)%X(this%Nx,1,0)-5.0d0*this%hwYOld(v)%X(this%Nx,2,0)+4.0d0*this%hwYOld(v)%X(this%Nx,3,0)-this%hwYOld(v)%X(this%Nx,4,0))/(this%hhy**2);   ! forward one-sided in Y
                   df2 = (2.0d0*this%hwYOld(v)%X(this%Nx,1,0)-5.0d0*this%hwYOld(v)%X(this%Nx,1,1)+4.0d0*this%hwYOld(v)%X(this%Nx,1,2)-this%hwYOld(v)%X(this%Nx,1,3))/(this%hhz**2);   ! forward one-sided in Z
                elseif (j==this%Ny.AND.k==0) then                   ! j=Ny, k=0
                   df1 = (2.0d0*this%hwYOld(v)%X(this%Nx,this%Ny,0)-5.0d0*this%hwYOld(v)%X(this%Nx,this%Ny-1,0)+4.0d0*this%hwYOld(v)%X(this%Nx,this%Ny-2,0)-this%hwYOld(v)%X(this%Nx,this%Ny-3,0))/(this%hhy**2);    ! backward one-sided in Y
                   df2 = (2.0d0*this%hwYOld(v)%X(this%Nx,this%Ny,0)-5.0d0*this%hwYOld(v)%X(this%Nx,this%Ny,1)+4.0d0*this%hwYOld(v)%X(this%Nx,this%Ny,2)-this%hwYOld(v)%X(this%Nx,this%Ny,3))/(this%hhz**2);          ! forward one-sided in Z 
                elseif (j==1.AND.k==this%Nz) then                   ! j=1, k=Nz
                   df1 = (2.0d0*this%hwYOld(v)%X(this%Nx,1,this%Nz)-5.0d0*this%hwYOld(v)%X(this%Nx,2,this%Nz)+4.0d0*this%hwYOld(v)%X(this%Nx,3,this%Nz)-this%hwYOld(v)%X(this%Nx,4,this%Nz))/(this%hhy**2);          ! forward one-sided in Y
                   df2 = (2.0d0*this%hwYOld(v)%X(this%Nx,1,this%Nz)-5.0d0*this%hwYOld(v)%X(this%Nx,1,this%Nz-1)+4.0d0*this%hwYOld(v)%X(this%Nx,1,this%Nz-2)-this%hwYOld(v)%X(this%Nx,1,this%Nz-3))/(this%hhz**2);    ! backward one-sided in Z 
                elseif (j==this%Ny.AND.k==this%Nz) then             ! j=Ny, k=Nz
                   df1 = (2.0d0*this%hwYOld(v)%X(this%Nx,this%Ny,this%Nz)-5.0d0*this%hwYOld(v)%X(this%Nx,this%Ny-1,this%Nz)+4.0d0*this%hwYOld(v)%X(this%Nx,this%Ny-2,this%Nz)-this%hwYOld(v)%X(this%Nx,this%Ny-3,this%Nz))/(this%hhy**2);  ! backward one-sided in Y
                   df2 = (2.0d0*this%hwYOld(v)%X(this%Nx,this%Ny,this%Nz)-5.0d0*this%hwYOld(v)%X(this%Nx,this%Ny,this%Nz-1)+4.0d0*this%hwYOld(v)%X(this%Nx,this%Ny,this%Nz-2)-this%hwYOld(v)%X(this%Nx,this%Ny,this%Nz-3))/(this%hhz**2);  ! backward one-sided in Z 
                elseif (j>1.AND.j<this%Ny.AND.k==0) then            ! j>1, j<Ny, k=0
                   df1 = (this%hwYOld(v)%X(this%Nx,j+1,0)-2.0d0*this%hwYOld(v)%X(this%Nx,j,0)+this%hwYOld(v)%X(this%Nx,j-1,0))/(this%hhy**2);                                                                  ! central in Y
                   df2 = (2.0d0*this%hwYOld(v)%X(this%Nx,j,0)-5.0d0*this%hwYOld(v)%X(this%Nx,j,1)+4.0d0*this%hwYOld(v)%X(this%Nx,j,2)-this%hwYOld(v)%X(this%Nx,j,3))/(this%hhz**2);                                  ! forward one-sided in Z      
                elseif (j>1.AND.j<this%Ny.AND.k==this%Nz) then      ! j>1, j<Ny, k=Nz
                   df1 = (this%hwYOld(v)%X(this%Nx,j+1,this%Nz)-2.0d0*this%hwYOld(v)%X(this%Nx,j,this%Nz)+this%hwYOld(v)%X(this%Nx,j-1,this%Nz))/(this%hhy**2);                                                ! central in Y
                   df2 = (2.0d0*this%hwYOld(v)%X(this%Nx,j,this%Nz)-5.0d0*this%hwYOld(v)%X(this%Nx,j,this%Nz-1)+4.0d0*this%hwYOld(v)%X(this%Nx,j,this%Nz-2)-this%hwYOld(v)%X(this%Nx,j,this%Nz-3))/(this%hhz**2);    ! backward one-sided Z 
                elseif (j==1.AND.k>0.AND.k<this%Nz) then            ! j=1, k>0, k<Nz
                   df1 = (2.0d0*this%hwYOld(v)%X(this%Nx,1,k)-5.0d0*this%hwYOld(v)%X(this%Nx,2,k)+4.0d0*this%hwYOld(v)%X(this%Nx,3,k)-this%hwYOld(v)%X(this%Nx,4,k))/(this%hhy**2);                                  ! forward one-sided in Y
                   df2 = (this%hwYOld(v)%X(this%Nx,1,k+1)-2.0d0*this%hwYOld(v)%X(this%Nx,1,k)+this%hwYOld(v)%X(this%Nx,1,k-1))/(this%hhz**2);                                                                  ! central in Z 
                elseif (j==this%Ny.AND.k>0.AND.k<this%Nz) then      ! j=Ny, k>0, k<Nz
                   df1 = (2.0d0*this%hwYOld(v)%X(this%Nx,this%Ny,k)-5.0d0*this%hwYOld(v)%X(this%Nx,this%Ny-1,k)+4.0d0*this%hwYOld(v)%X(this%Nx,this%Ny-2,k)-this%hwYOld(v)%X(this%Nx,this%Ny-3,k))/(this%hhy**2);    ! backward one-sided in Y
                   df2 = (this%hwYOld(v)%X(this%Nx,this%Ny,k+1)-2.0d0*this%hwYOld(v)%X(this%Nx,this%Ny,k)+this%hwYOld(v)%X(this%Nx,this%Ny,k-1))/(this%hhz**2);                                                ! central in Z 
                else
                   df1 = (this%hwYOld(v)%X(this%Nx,j+1,k)-2.0d0*this%hwYOld(v)%X(this%Nx,j,k)+this%hwYOld(v)%X(this%Nx,j-1,k))/(this%hhy**2);    ! central in Y
                   df2 = (this%hwYOld(v)%X(this%Nx,j,k+1)-2.0d0*this%hwYOld(v)%X(this%Nx,j,k)+this%hwYOld(v)%X(this%Nx,j,k-1))/(this%hhz**2);    ! central in Z 
                endif
             end select
             this%hwGF(v) = df1+df2;
          enddo

          do v=0,im
             temp(v) = 0.0d0;
          enddo
          ! B*Y^(n)
          do v=0,im
             tp = 0.0d0;
             do w=0,im
                tp = tp - this%hwB(v,w)*this%hwF(w);
             enddo
             temp(v) = temp(v)+tp;
          enddo
          ! C*Y^(n-1)
          do v=0,im
             tp = 0.0d0;
             do w=0,im
                tp = tp - this%hwC(v,w)*this%hwFold(w);
             enddo
             temp(v) = temp(v)+tp;
          enddo
          ! D*GY^(n)
          do v=0,im
             tp = 0.0d0;
             do w=0,im
               tp = tp + this%hwD(v,w)*this%hwGF(w);
             enddo
             temp(v) = temp(v)+tp;
          enddo
          !add dx_u
          select case(this%dtmode)
          case(1)
             temp(0) = temp(0)-1.0d0*cc*this%ht*de;
          case(2)
             temp(0) = temp(0)-2.0d0*cc*this%ht*de;
          end select
          ! A^(-1)*()
          do v=0,im
             tp = 0.0d0;
             do w=0,im
                tp=tp + this%hwA(v,w)*temp(w);
             enddo
             this%hwY(v)%X(this%Nx, j, k) = tp;
          enddo
          this%Ef%Y(this%Nx, j, k) = this%hwY(0)%X(this%Nx, j, k);
       enddo
    enddo

    do j=2,this%Ny-1
       do v=0,im  ! Vectors formation
          !k=0
          this%hwY(v)%X(this%Nx, j, 0) = gz*this%hwY(v)%X(this%Nx, j, 0) - gz*this%hwY(v)%X(this%Nx, j, 1) + this%hwY(v)%X(this%Nx, j, 1);
          !k=Nz
          this%hwY(v)%X(this%Nx, j, this%Nz) = gz*this%hwY(v)%X(this%Nx, j, this%Nz) - gz*this%hwY(v)%X(this%Nx, j, this%Nz-1) + this%hwY(v)%X(this%Nx, j, this%Nz-1);
       enddo
       this%Ef%Y(this%Nx, j, 0) = this%hwY(0)%X(this%Nx, j, 0);
       this%Ef%Y(this%Nx, j, this%Nz) = this%hwY(0)%X(this%Nx, j, this%Nz);
    enddo   

    do k=1,this%Nz-1
       do v=0,im  ! Vectors formation
          !j=1
          this%hwY(v)%X(this%Nx, 1, k) = gy*this%hwY(v)%X(this%Nx, 1, k) - gy*this%hwY(v)%X(this%Nx, 2, k) + this%hwY(v)%X(this%Nx, 2, k);
          !j=Ny
          this%hwY(v)%X(this%Nx, this%Ny, k) = gy*this%hwY(v)%X(this%Nx, this%Ny, k) - gy*this%hwY(v)%X(this%Nx, this%Ny-1, k) + this%hwY(v)%X(this%Nx, this%Ny-1, k);
       enddo
       this%Ef%Y(this%Nx, 1, k) = this%hwY(0)%X(this%Nx, 1, k);
       this%Ef%Y(this%Nx, this%Ny, k) = this%hwY(0)%X(this%Nx, this%Ny, k);
    enddo

    do v=0,im  ! Vectors formation
       !j=1 k=0
          tp = gy*this%hwY(v)%X(this%Nx, 1, 0) - gy*this%hwY(v)%X(this%Nx, 2, 0) + this%hwY(v)%X(this%Nx, 2, 0);
          this%hwY(v)%X(this%Nx, 1, 0) = (tp + gz*this%hwY(v)%X(this%Nx, 1, 0) - gz*this%hwY(v)%X(this%Nx, 1, 1) + this%hwY(v)%X(this%Nx, 1, 1))/2.0d0;     
       !j=1 k=Nz
          tp = gy*this%hwY(v)%X(this%Nx, 1, this%Nz) - gy*this%hwY(v)%X(this%Nx, 2, this%Nz) + this%hwY(v)%X(this%Nx, 2, this%Nz);
          this%hwY(v)%X(this%Nx, 1, this%Nz) = (tp + gz*this%hwY(v)%X(this%Nx, 1, this%Nz) - gz*this%hwY(v)%X(this%Nx, 1, this%Nz-1) + this%hwY(v)%X(this%Nx, 1, this%Nz-1))/2.0d0;
       !j=Ny k=0
          tp = gy*this%hwY(v)%X(this%Nx, this%Ny, 0) - gy*this%hwY(v)%X(this%Nx, this%Ny-1, 0) + this%hwY(v)%X(this%Nx, this%Ny-1, 0);
          this%hwY(v)%X(this%Nx, this%Ny, 0) = (tp + gz*this%hwY(v)%X(this%Nx, this%Ny, 0) - gz*this%hwY(v)%X(this%Nx, this%Ny, 1) + this%hwY(v)%X(this%Nx, this%Ny, 1))/2.0d0;
       !j=Ny k=Nz
          tp = gy*this%hwY(v)%X(this%Nx, this%Ny, this%Nz) - gy*this%hwY(v)%X(this%Nx, this%Ny-1, this%Nz) + this%hwY(v)%X(this%Nx, this%Ny-1, this%Nz);
          this%hwY(v)%X(this%Nx, this%Ny, this%Nz) = (tp + gz*this%hwY(v)%X(this%Nx, this%Ny, this%Nz) - gz*this%hwY(v)%X(this%Nx, this%Ny, this%Nz-1) + this%hwY(v)%X(this%Nx, this%Ny, this%Nz-1))/2.0d0;
    enddo

    this%Ef%Y(this%Nx, 1, 0) = this%hwY(0)%X(this%Nx, 1, 0);
    this%Ef%Y(this%Nx, 1, this%Nz) = this%hwY(0)%X(this%Nx, 1, this%Nz);
    this%Ef%Y(this%Nx, this%Ny, 0) = this%hwY(0)%X(this%Nx, this%Ny, 0);
    this%Ef%Y(this%Nx, this%Ny, this%Nz) = this%hwY(0)%X(this%Nx, this%Ny, this%Nz);   
   
    deallocate(temp);
  end subroutine HWEvaEySomUpdate


!HWEvaEyHigdonUpdate------------------------------------------------------------------
  subroutine HWEvaEyHigdonUpdate(this, t)
    !Fill Ey_x=0 boundary with ABC Hagstrom-Warburton, Higdon for boundaries of aux functions
    class(tProblem) :: this
    integer :: t;
    integer :: j ,k, v, w, im
    real*8 :: de, df1, df2, tp;
    character*8 :: trans;
    real*8, allocatable :: temp(:);
    real*8 :: gx1, gx2, gy1, gy2, gz1, gz2;
    write(*,*) 'Do Higdon';
    gx1 = (1.0d0-cc*this%ht/(cos(this%hig1/180.0d0*PI)*this%hhx))/(1.0d0+cc*this%ht/(cos(this%hig1/180.0d0*PI)*this%hhx));
    gx2 = (1.0d0-cc*this%ht/(cos(this%hig2/180.0d0*PI)*this%hhx))/(1.0d0+cc*this%ht/(cos(this%hig2/180.0d0*PI)*this%hhx));
    gy1 = (1.0d0-cc*this%ht/(cos(this%hig1/180.0d0*PI)*this%hhy))/(1.0d0+cc*this%ht/(cos(this%hig1/180.0d0*PI)*this%hhy));
    gy2 = (1.0d0-cc*this%ht/(cos(this%hig2/180.0d0*PI)*this%hhy))/(1.0d0+cc*this%ht/(cos(this%hig2/180.0d0*PI)*this%hhy));
    gz1 = (1.0d0-cc*this%ht/(cos(this%hig1/180.0d0*PI)*this%hhx))/(1.0d0+cc*this%ht/(cos(this%hig1/180.0d0*PI)*this%hhz));
    gz2 = (1.0d0-cc*this%ht/(cos(this%hig2/180.0d0*PI)*this%hhx))/(1.0d0+cc*this%ht/(cos(this%hig2/180.0d0*PI)*this%hhz));
    if (this%hwEvaE>0) then
       im=this%hwEvaP+this%hwEvaE+1
    else
       im=this%hwEvaP
    endif
    allocate(temp(0:im));
    !=========================== EY == X BOUNDARIES ============================
    !EYFace X=0 -----
    do k=1,this%Nz-1
       do j=2,this%Ny-1
          select case(this%dxmode)
          case(1)
             de = (this%EfOld%Y(this%Nx, j, k) - this%EfOld%Y(this%Nx-1, j, k))/(this%hhx);
          case(2)
             de = (3.0d0*this%EfOld%Y(this%Nx, j, k) - 4.0d0*this%EfOld%Y(this%Nx-1, j, k) + this%EfOld%Y(this%Nx-2, j, k) )/(2.0d0*this%hhx);
          end select
          do v=0,im                  ! Vectors formation
             this%hwF(v) = this%hwYOld(v)%X(this%Nx, j, k);
             this%hwFold(v) = this%hwYold2(v)%X(this%Nx, j, k);
             select case(this%dxauxmode)
             case(1) !===1st order===
               if (j==1.AND.k==0) then                             ! j=1, k=0
                  df1 = (this%hwYOld(v)%X(this%Nx,1,0)-2.0d0*this%hwYOld(v)%X(this%Nx,2,0)+this%hwYOld(v)%X(this%Nx,3,0))/(this%hhy**2);   ! forward one-sided in Y
                  df2 = (this%hwYOld(v)%X(this%Nx,1,0)-2.0d0*this%hwYOld(v)%X(this%Nx,1,1)+this%hwYOld(v)%X(this%Nx,1,2))/(this%hhz**2);   ! forward one-sided in Z
               elseif (j==this%Ny.AND.k==0) then                   ! j=Ny, k=0
                  df1 = (this%hwYOld(v)%X(this%Nx,this%Ny,0)-2.0d0*this%hwYOld(v)%X(this%Nx,this%Ny-1,0)+this%hwYOld(v)%X(this%Nx,this%Ny-2,0))/(this%hhy**2);     ! backward one-sided in Y
                  df2 = (this%hwYOld(v)%X(this%Nx,this%Ny,0)-2.0d0*this%hwYOld(v)%X(this%Nx,this%Ny,1)+this%hwYOld(v)%X(this%Nx,this%Ny,2))/(this%hhz**2);         ! forward one-sided in Z 
               elseif (j==1.AND.k==this%Nz) then                   ! j=1, k=Nz
                  df1 = (this%hwYOld(v)%X(this%Nx,1,this%Nz)-2.0d0*this%hwYOld(v)%X(this%Nx,2,this%Nz)+this%hwYOld(v)%X(this%Nx,3,this%Nz))/(this%hhy**2);         ! forward one-sided in Y
                  df2 = (this%hwYOld(v)%X(this%Nx,1,this%Nz)-2.0d0*this%hwYOld(v)%X(this%Nx,1,this%Nz-1)+this%hwYOld(v)%X(this%Nx,1,this%Nz-2))/(this%hhz**2);     ! backward one-sided in Z 
               elseif (j==this%Ny.AND.k==this%Nz) then             ! j=Ny, k=Nz
                  df1 = (this%hwYOld(v)%X(this%Nx,this%Ny,this%Nz)-2.0d0*this%hwYOld(v)%X(this%Nx,this%Ny-1,this%Nz)+this%hwYOld(v)%X(this%Nx,this%Ny-2,this%Nz))/(this%hhy**2);  ! backward one-sided in Y
                  df2 = (this%hwYOld(v)%X(this%Nx,this%Ny,this%Nz)-2.0d0*this%hwYOld(v)%X(this%Nx,this%Ny,this%Nz-1)+this%hwYOld(v)%X(this%Nx,this%Ny,this%Nz-2))/(this%hhz**2);  ! backward one-sided in Z 
               elseif (j>1.AND.j<this%Ny.AND.k==0) then            ! j>1, j<Ny, k=0
                  df1 = (this%hwYOld(v)%X(this%Nx,j+1,0)-2.0d0*this%hwYOld(v)%X(this%Nx,j,0)+this%hwYOld(v)%X(this%Nx,j-1,0))/(this%hhy**2);                                      ! central in Y
                  df2 = (this%hwYOld(v)%X(this%Nx,j,0)-2.0d0*this%hwYOld(v)%X(this%Nx,j,1)+this%hwYOld(v)%X(this%Nx,j,2))/(this%hhz**2);                                          ! forward one-sided in Z      
               elseif (j>1.AND.j<this%Ny.AND.k==this%Nz) then      ! j>1, j<Ny, k=Nz
                  df1 = (this%hwYOld(v)%X(this%Nx,j+1,this%Nz)-2.0d0*this%hwYOld(v)%X(this%Nx,j,this%Nz)+this%hwYOld(v)%X(this%Nx,j-1,this%Nz))/(this%hhy**2);                    ! central in Y
                  df2 = (this%hwYOld(v)%X(this%Nx,j,this%Nz)-2.0d0*this%hwYOld(v)%X(this%Nx,j,this%Nz-1)+this%hwYOld(v)%X(this%Nx,j,this%Nz-2))/(this%hhz**2);                    ! backward one-sided Z 
               elseif (j==1.AND.k>0.AND.k<this%Nz) then            ! j=1, k>0, k<Nz
                  df1 = (this%hwYOld(v)%X(this%Nx,1,k)-2.0d0*this%hwYOld(v)%X(this%Nx,2,k)+this%hwYOld(v)%X(this%Nx,3,k))/(this%hhy**2);                                           ! forward one-sided in Y
                  df2 = (this%hwYOld(v)%X(this%Nx,1,k+1)-2.0d0*this%hwYOld(v)%X(this%Nx,1,k)+this%hwYOld(v)%X(this%Nx,1,k-1))/(this%hhz**2);                                       ! central in Z 
               elseif (j==this%Ny.AND.k>0.AND.k<this%Nz) then      ! j=Ny, k>0, k<Nz
                  df1 = (this%hwYOld(v)%X(this%Nx,this%Ny,k)-2.0d0*this%hwYOld(v)%X(this%Nx,this%Ny-1,k)+this%hwYOld(v)%X(this%Nx,this%Ny-2,k))/(this%hhy**2);    ! backward one-sided in Y
                  df2 = (this%hwYOld(v)%X(this%Nx,this%Ny,k+1)-2.0d0*this%hwYOld(v)%X(this%Nx,this%Ny,k)+this%hwYOld(v)%X(this%Nx,this%Ny,k-1))/(this%hhz**2);    ! central in Z 
               else
                  df1 = (this%hwYOld(v)%X(this%Nx,j+1,k)-2.0d0*this%hwYOld(v)%X(this%Nx,j,k)+this%hwYOld(v)%X(this%Nx,j-1,k))/(this%hhy**2);    ! central in Y
                  df2 = (this%hwYOld(v)%X(this%Nx,j,k+1)-2.0d0*this%hwYOld(v)%X(this%Nx,j,k)+this%hwYOld(v)%X(this%Nx,j,k-1))/(this%hhz**2);    ! central in Z 
               endif
             case(2) !===2nd order===
                if (j==1.AND.k==0) then                             ! j=1, k=0
                   df1 = (2.0d0*this%hwYOld(v)%X(this%Nx,1,0)-5.0d0*this%hwYOld(v)%X(this%Nx,2,0)+4.0d0*this%hwYOld(v)%X(this%Nx,3,0)-this%hwYOld(v)%X(this%Nx,4,0))/(this%hhy**2);   ! forward one-sided in Y
                   df2 = (2.0d0*this%hwYOld(v)%X(this%Nx,1,0)-5.0d0*this%hwYOld(v)%X(this%Nx,1,1)+4.0d0*this%hwYOld(v)%X(this%Nx,1,2)-this%hwYOld(v)%X(this%Nx,1,3))/(this%hhz**2);   ! forward one-sided in Z
                elseif (j==this%Ny.AND.k==0) then                   ! j=Ny, k=0
                   df1 = (2.0d0*this%hwYOld(v)%X(this%Nx,this%Ny,0)-5.0d0*this%hwYOld(v)%X(this%Nx,this%Ny-1,0)+4.0d0*this%hwYOld(v)%X(this%Nx,this%Ny-2,0)-this%hwYOld(v)%X(this%Nx,this%Ny-3,0))/(this%hhy**2);    ! backward one-sided in Y
                   df2 = (2.0d0*this%hwYOld(v)%X(this%Nx,this%Ny,0)-5.0d0*this%hwYOld(v)%X(this%Nx,this%Ny,1)+4.0d0*this%hwYOld(v)%X(this%Nx,this%Ny,2)-this%hwYOld(v)%X(this%Nx,this%Ny,3))/(this%hhz**2);          ! forward one-sided in Z 
                elseif (j==1.AND.k==this%Nz) then                   ! j=1, k=Nz
                   df1 = (2.0d0*this%hwYOld(v)%X(this%Nx,1,this%Nz)-5.0d0*this%hwYOld(v)%X(this%Nx,2,this%Nz)+4.0d0*this%hwYOld(v)%X(this%Nx,3,this%Nz)-this%hwYOld(v)%X(this%Nx,4,this%Nz))/(this%hhy**2);          ! forward one-sided in Y
                   df2 = (2.0d0*this%hwYOld(v)%X(this%Nx,1,this%Nz)-5.0d0*this%hwYOld(v)%X(this%Nx,1,this%Nz-1)+4.0d0*this%hwYOld(v)%X(this%Nx,1,this%Nz-2)-this%hwYOld(v)%X(this%Nx,1,this%Nz-3))/(this%hhz**2);    ! backward one-sided in Z 
                elseif (j==this%Ny.AND.k==this%Nz) then             ! j=Ny, k=Nz
                   df1 = (2.0d0*this%hwYOld(v)%X(this%Nx,this%Ny,this%Nz)-5.0d0*this%hwYOld(v)%X(this%Nx,this%Ny-1,this%Nz)+4.0d0*this%hwYOld(v)%X(this%Nx,this%Ny-2,this%Nz)-this%hwYOld(v)%X(this%Nx,this%Ny-3,this%Nz))/(this%hhy**2);  ! backward one-sided in Y
                   df2 = (2.0d0*this%hwYOld(v)%X(this%Nx,this%Ny,this%Nz)-5.0d0*this%hwYOld(v)%X(this%Nx,this%Ny,this%Nz-1)+4.0d0*this%hwYOld(v)%X(this%Nx,this%Ny,this%Nz-2)-this%hwYOld(v)%X(this%Nx,this%Ny,this%Nz-3))/(this%hhz**2);  ! backward one-sided in Z 
                elseif (j>1.AND.j<this%Ny.AND.k==0) then            ! j>1, j<Ny, k=0
                   df1 = (this%hwYOld(v)%X(this%Nx,j+1,0)-2.0d0*this%hwYOld(v)%X(this%Nx,j,0)+this%hwYOld(v)%X(this%Nx,j-1,0))/(this%hhy**2);                                                                  ! central in Y
                   df2 = (2.0d0*this%hwYOld(v)%X(this%Nx,j,0)-5.0d0*this%hwYOld(v)%X(this%Nx,j,1)+4.0d0*this%hwYOld(v)%X(this%Nx,j,2)-this%hwYOld(v)%X(this%Nx,j,3))/(this%hhz**2);                                  ! forward one-sided in Z      
                elseif (j>1.AND.j<this%Ny.AND.k==this%Nz) then      ! j>1, j<Ny, k=Nz
                   df1 = (this%hwYOld(v)%X(this%Nx,j+1,this%Nz)-2.0d0*this%hwYOld(v)%X(this%Nx,j,this%Nz)+this%hwYOld(v)%X(this%Nx,j-1,this%Nz))/(this%hhy**2);                                                ! central in Y
                   df2 = (2.0d0*this%hwYOld(v)%X(this%Nx,j,this%Nz)-5.0d0*this%hwYOld(v)%X(this%Nx,j,this%Nz-1)+4.0d0*this%hwYOld(v)%X(this%Nx,j,this%Nz-2)-this%hwYOld(v)%X(this%Nx,j,this%Nz-3))/(this%hhz**2);    ! backward one-sided Z 
                elseif (j==1.AND.k>0.AND.k<this%Nz) then            ! j=1, k>0, k<Nz
                   df1 = (2.0d0*this%hwYOld(v)%X(this%Nx,1,k)-5.0d0*this%hwYOld(v)%X(this%Nx,2,k)+4.0d0*this%hwYOld(v)%X(this%Nx,3,k)-this%hwYOld(v)%X(this%Nx,4,k))/(this%hhy**2);                                  ! forward one-sided in Y
                   df2 = (this%hwYOld(v)%X(this%Nx,1,k+1)-2.0d0*this%hwYOld(v)%X(this%Nx,1,k)+this%hwYOld(v)%X(this%Nx,1,k-1))/(this%hhz**2);                                                                  ! central in Z 
                elseif (j==this%Ny.AND.k>0.AND.k<this%Nz) then      ! j=Ny, k>0, k<Nz
                   df1 = (2.0d0*this%hwYOld(v)%X(this%Nx,this%Ny,k)-5.0d0*this%hwYOld(v)%X(this%Nx,this%Ny-1,k)+4.0d0*this%hwYOld(v)%X(this%Nx,this%Ny-2,k)-this%hwYOld(v)%X(this%Nx,this%Ny-3,k))/(this%hhy**2);    ! backward one-sided in Y
                   df2 = (this%hwYOld(v)%X(this%Nx,this%Ny,k+1)-2.0d0*this%hwYOld(v)%X(this%Nx,this%Ny,k)+this%hwYOld(v)%X(this%Nx,this%Ny,k-1))/(this%hhz**2);                                                ! central in Z 
                else
                   df1 = (this%hwYOld(v)%X(this%Nx,j+1,k)-2.0d0*this%hwYOld(v)%X(this%Nx,j,k)+this%hwYOld(v)%X(this%Nx,j-1,k))/(this%hhy**2);    ! central in Y
                   df2 = (this%hwYOld(v)%X(this%Nx,j,k+1)-2.0d0*this%hwYOld(v)%X(this%Nx,j,k)+this%hwYOld(v)%X(this%Nx,j,k-1))/(this%hhz**2);    ! central in Z 
                endif
             end select
             this%hwGF(v) = df1+df2;
          enddo

          do v=0,im
             temp(v) = 0.0d0;
          enddo
          ! B*Y^(n)
          do v=0,im
             tp = 0.0d0;
             do w=0,im
                tp = tp - this%hwB(v,w)*this%hwF(w);
             enddo
             temp(v) = temp(v)+tp;
          enddo
          ! C*Y^(n-1)
          do v=0,im
             tp = 0.0d0;
             do w=0,im
                tp = tp - this%hwC(v,w)*this%hwFold(w);
             enddo
             temp(v) = temp(v)+tp;
          enddo
          ! D*GY^(n)
          do v=0,im
             tp = 0.0d0;
             do w=0,im
               tp = tp + this%hwD(v,w)*this%hwGF(w);
             enddo
             temp(v) = temp(v)+tp;
          enddo
          !add dx_u
          select case(this%dtmode)
          case(1)
             temp(0) = temp(0)-1.0d0*cc*this%ht*de;
          case(2)
             temp(0) = temp(0)-2.0d0*cc*this%ht*de;
          end select
          ! A^(-1)*()
          do v=0,im
             tp = 0.0d0;
             do w=0,im
                tp=tp + this%hwA(v,w)*temp(w);
             enddo
             this%hwY(v)%X(this%Nx, j, k) = tp;
          enddo
          this%Ef%Y(this%Nx, j, k) = this%hwY(0)%X(this%Nx, j, k);
       enddo
    enddo
   
    do j=2,this%Ny-1
       do v=0,im  ! Vectors formation
          !k=0
          this%hwY(v)%X(this%Nx, j, 0) = -(gz1+gz2)*this%hwY(v)%X(this%Nx, j, 1)-gz1*gz2*this%hwY(v)%X(this%Nx, j, 2)  &
          +(gz1+gz2)*this%hwYOld(v)%X(this%Nx, j, 0)+(2.0d0+2.0d0*gz1*gz2)*this%hwYOld(v)%X(this%Nx, j, 1)+(gz1+gz2)*this%hwYOld(v)%X(this%Nx, j, 2)  &
          -gz1*gz2*this%hwYOld2(v)%X(this%Nx, j, 0)-(gz1+gz2)*this%hwYOld2(v)%X(this%Nx, j, 1)-this%hwYOld2(v)%X(this%Nx, j, 2);
          !k=Nz         
          this%hwY(v)%X(this%Nx, j, this%Nz) = -(gz1+gz2)*this%hwY(v)%X(this%Nx, j, this%Nz-1)-gz1*gz2*this%hwY(v)%X(this%Nx, j, this%Nz-2)  &
          +(gz1+gz2)*this%hwYOld(v)%X(this%Nx, j, this%Nz)+(2.0d0+2.0d0*gz1*gz2)*this%hwYOld(v)%X(this%Nx, j, this%Nz-1)+(gz1+gz2)*this%hwYOld(v)%X(this%Nx, j, this%Nz-2)  &
          -gz1*gz2*this%hwYOld2(v)%X(this%Nx, j, this%Nz)-(gz1+gz2)*this%hwYOld2(v)%X(this%Nx, j, this%Nz-1)-this%hwYOld2(v)%X(this%Nx, j, this%Nz-2);
       enddo
       this%Ef%Y(this%Nx, j, 0) = this%hwY(0)%X(this%Nx, j, 0);
       this%Ef%Y(this%Nx, j, this%Nz) = this%hwY(0)%X(this%Nx, j, this%Nz);
    enddo   

    do k=1,this%Nz-1
       do v=0,im  ! Vectors formation
          !j=1
          this%hwY(v)%X(this%Nx, 1, k) = -(gy1+gy2)*this%hwY(v)%X(this%Nx, 2, k)-gy1*gy2*this%hwY(v)%X(this%Nx, 3, k)  &
          +(gy1+gy2)*this%hwYOld(v)%X(this%Nx, 1, k)+(2.0d0+2.0d0*gy1*gy2)*this%hwYOld(v)%X(this%Nx, 2, k)+(gy1+gy2)*this%hwYOld(v)%X(this%Nx, 3, k)  &
          -gy1*gy2*this%hwYOld2(v)%X(this%Nx, 1, k)-(gy1+gy2)*this%hwYOld2(v)%X(this%Nx, 2, k)-this%hwYOld2(v)%X(this%Nx, 3, k);
         !j=Ny
          this%hwY(v)%X(this%Nx, this%Ny, k) = -(gy1+gy2)*this%hwY(v)%X(this%Nx, this%Ny-1, k)-gy1*gy2*this%hwY(v)%X(this%Nx, this%Ny-2, k)  &
          +(gy1+gy2)*this%hwYOld(v)%X(this%Nx, this%Ny, k)+(2.0d0+2.0d0*gy1*gy2)*this%hwYOld(v)%X(this%Nx, this%Ny-1, k)+(gy1+gy2)*this%hwYOld(v)%X(this%Nx, this%Ny-2, k)  &
          -gy1*gy2*this%hwYOld2(v)%X(this%Nx, this%Ny, k)-(gy1+gy2)*this%hwYOld2(v)%X(this%Nx, this%Ny-1, k)-this%hwYOld2(v)%X(this%Nx, this%Ny-2, k);
       enddo
       this%Ef%Y(this%Nx, 1, k) = this%hwY(0)%X(this%Nx, 1, k);
       this%Ef%Y(this%Nx, this%Ny, k) = this%hwY(0)%X(this%Nx, this%Ny, k);
    enddo

    do v=0,im  ! Vectors formation
       !j=1 k=0
          tp = -(gy1+gy2)*this%hwY(v)%X(this%Nx, 2, 0)-gy1*gy2*this%hwY(v)%X(this%Nx, 3, 0)  &
          +(gy1+gy2)*this%hwYOld(v)%X(this%Nx, 1, 0)+(2.0d0+2.0d0*gy1*gy2)*this%hwYOld(v)%X(this%Nx, 2, 0)+(gy1+gy2)*this%hwYOld(v)%X(this%Nx, 3, 0)  &
          -gy1*gy2*this%hwYOld2(v)%X(this%Nx, 1, 0)-(gy1+gy2)*this%hwYOld2(v)%X(this%Nx, 2, 0)-this%hwYOld2(v)%X(this%Nx, 3, 0);
          this%hwY(v)%X(this%Nx, 1, 0) = (tp -(gz1+gz2)*this%hwY(v)%X(this%Nx, 1, 1)-gz1*gz2*this%hwY(v)%X(this%Nx, 1, 2)  &
          +(gz1+gz2)*this%hwYOld(v)%X(this%Nx, 1, 0)+(2.0d0+2.0d0*gz1*gz2)*this%hwYOld(v)%X(this%Nx, 1, 1)+(gz1+gz2)*this%hwYOld(v)%X(this%Nx, 1, 2)  &
          -gz1*gz2*this%hwYOld2(v)%X(this%Nx, 1, 0)-(gz1+gz2)*this%hwYOld2(v)%X(this%Nx, 1, 1)-this%hwYOld2(v)%X(this%Nx, 1, 2))/2.0d0;
       !j=1 k=Nz
          tp = -(gy1+gy2)*this%hwY(v)%X(this%Nx, 2, this%Nz)-gy1*gy2*this%hwY(v)%X(this%Nx, 3, this%Nz)  &
          +(gy1+gy2)*this%hwYOld(v)%X(this%Nx, 1, this%Nz)+(2.0d0+2.0d0*gy1*gy2)*this%hwYOld(v)%X(this%Nx, 2, this%Nz)+(gy1+gy2)*this%hwYOld(v)%X(this%Nx, 3, this%Nz)  &
          -gy1*gy2*this%hwYOld2(v)%X(this%Nx, 1, this%Nz)-(gy1+gy2)*this%hwYOld2(v)%X(this%Nx, 2, this%Nz)-this%hwYOld2(v)%X(this%Nx, 3, this%Nz);
          this%hwY(v)%X(this%Nx, 1, this%Nz) = (tp -(gz1+gz2)*this%hwY(v)%X(this%Nx, 1, this%Nz-1)-gz1*gz2*this%hwY(v)%X(this%Nx, 1, this%Nz-2)  &
          +(gz1+gz2)*this%hwYOld(v)%X(this%Nx, 1, this%Nz)+(2.0d0+2.0d0*gz1*gz2)*this%hwYOld(v)%X(this%Nx, 1, this%Nz-1)+(gz1+gz2)*this%hwYOld(v)%X(this%Nx, 1, this%Nz-2)  &
          -gz1*gz2*this%hwYOld2(v)%X(this%Nx, 1, this%Nz)-(gz1+gz2)*this%hwYOld2(v)%X(this%Nx, 1, this%Nz-1)-this%hwYOld2(v)%X(this%Nx, 1, this%Nz-2))/2.0d0;
       !j=Ny k=0
          tp = -(gy1+gy2)*this%hwY(v)%X(this%Nx, this%Ny-1, 0)-gy1*gy2*this%hwY(v)%X(this%Nx, this%Ny-2, 0)  &
          +(gy1+gy2)*this%hwYOld(v)%X(this%Nx, this%Ny, 0)+(2.0d0+2.0d0*gy1*gy2)*this%hwYOld(v)%X(this%Nx, this%Ny-1, 0)+(gy1+gy2)*this%hwYOld(v)%X(this%Nx, this%Ny-2, 0)  &
          -gy1*gy2*this%hwYOld2(v)%X(this%Nx, this%Ny, 0)-(gy1+gy2)*this%hwYOld2(v)%X(this%Nx, this%Ny-1, 0)-this%hwYOld2(v)%X(this%Nx, this%Ny-2, 0);
          this%hwY(v)%X(this%Nx, this%Ny, 0) = (tp -(gz1+gz2)*this%hwY(v)%X(this%Nx, this%Ny, 1)-gz1*gz2*this%hwY(v)%X(this%Nx, this%Ny, 2)  &
          +(gz1+gz2)*this%hwYOld(v)%X(this%Nx, this%Ny, 0)+(2.0d0+2.0d0*gz1*gz2)*this%hwYOld(v)%X(this%Nx, this%Ny, 1)+(gz1+gz2)*this%hwYOld(v)%X(this%Nx, this%Ny, 2)  &
          -gz1*gz2*this%hwYOld2(v)%X(this%Nx, this%Ny, 0)-(gz1+gz2)*this%hwYOld2(v)%X(this%Nx, this%Ny, 1)-this%hwYOld2(v)%X(this%Nx, this%Ny, 2))/2.0d0;
       !j=Ny k=Nz
          tp = -(gy1+gy2)*this%hwY(v)%X(this%Nx, this%Ny-1, this%Nz)-gy1*gy2*this%hwY(v)%X(this%Nx, this%Ny-2, this%Nz)  &
          +(gy1+gy2)*this%hwYOld(v)%X(this%Nx, this%Ny, this%Nz)+(2.0d0+2.0d0*gy1*gy2)*this%hwYOld(v)%X(this%Nx, this%Ny-1, this%Nz)+(gy1+gy2)*this%hwYOld(v)%X(this%Nx, this%Ny-2, this%Nz)  &
          -gy1*gy2*this%hwYOld2(v)%X(this%Nx, this%Ny, this%Nz)-(gy1+gy2)*this%hwYOld2(v)%X(this%Nx, this%Ny-1, this%Nz)-this%hwYOld2(v)%X(this%Nx, this%Ny-2, this%Nz);
          this%hwY(v)%X(this%Nx, this%Ny, this%Nz) = (tp - (gz1+gz2)*this%hwY(v)%X(this%Nx, this%Ny, this%Nz-1)-gz1*gz2*this%hwY(v)%X(this%Nx, this%Ny, this%Nz-2)  &
          +(gz1+gz2)*this%hwYOld(v)%X(this%Nx, this%Ny, this%Nz)+(2.0d0+2.0d0*gz1*gz2)*this%hwYOld(v)%X(this%Nx, this%Ny, this%Nz-1)+(gz1+gz2)*this%hwYOld(v)%X(this%Nx, this%Ny, this%Nz-2)  &
          -gz1*gz2*this%hwYOld2(v)%X(this%Nx, this%Ny, this%Nz)-(gz1+gz2)*this%hwYOld2(v)%X(this%Nx, this%Ny, this%Nz-1)-this%hwYOld2(v)%X(this%Nx, this%Ny, this%Nz-2) )/2.0d0;
       enddo
       
    this%Ef%Y(this%Nx, 1, 0) = this%hwY(0)%X(this%Nx, 1, 0);
    this%Ef%Y(this%Nx, 1, this%Nz) = this%hwY(0)%X(this%Nx, 1, this%Nz);
    this%Ef%Y(this%Nx, this%Ny, 0) = this%hwY(0)%X(this%Nx, this%Ny, 0);
    this%Ef%Y(this%Nx, this%Ny, this%Nz) = this%hwY(0)%X(this%Nx, this%Ny, this%Nz);   
   
    deallocate(temp);
  end subroutine HWEvaEyHigdonUpdate


!HWEvaEyBetzMittraUpdate------------------------------------------------------------------
  subroutine HWEvaEyBetzMittraUpdate(this, t)
    !Fill Ey_x=0 boundary with ABC Hagstrom-Warburton, Betz-Mittra for boundaries of aux functions
    class(tProblem) :: this
    integer :: t;
    integer :: j ,k, v, w, im
    real*8 :: de, df1, df2, tp, tp1, tp2;
    character*8 :: trans;
    real*8, allocatable :: temp(:);
    real*8 :: alp, ro1x, ro2x, betax, gx1, gx2, ro1y, ro2y, betay, gy1, gy2, ro1z, ro2z, betaz, gz1, gz2;
    write(*,*) 'Do Betz-Mittra';
    alp = 0.1d0/this%hhx;
    ro1x = cc*this%ht/(cos(this%hig1/180.0d0*PI)*this%hhx);
    ro2x = cc*this%ht/(cos(this%hig2/180.0d0*PI)*this%hhx);
    betax = (1.0d0+ro1x)/(1.0d0+ro1x*(1.0d0+alp*this%hhx));
    gx1 = (1.0d0-ro1x)/(1.0d0+ro1x*(1.0d0+alp*this%hhx));
    gx2 = (1.0d0-ro2x)/(1.0d0+ro2x*(1.0d0+alp*this%hhx));

    ro1y = cc*this%ht/(cos(this%hig1/180.0d0*PI)*this%hhy);
    ro2y = cc*this%ht/(cos(this%hig2/180.0d0*PI)*this%hhy);
    betay = (1.0d0+ro1y)/(1.0d0+ro1y*(1.0d0+alp*this%hhy));
    gy1 = (1.0d0-ro1y)/(1.0d0+ro1y*(1.0d0+alp*this%hhy));
    gy2 = (1.0d0-ro2y)/(1.0d0+ro2y*(1.0d0+alp*this%hhy));

    ro1z = cc*this%ht/(cos(this%hig1/180.0d0*PI)*this%hhz);
    ro2z = cc*this%ht/(cos(this%hig2/180.0d0*PI)*this%hhz);
    betaz = (1.0d0+ro1z)/(1.0d0+ro1z*(1.0d0+alp*this%hhz));
    gz1 = (1.0d0-ro1z)/(1.0d0+ro1z*(1.0d0+alp*this%hhz));
    gz2 = (1.0d0-ro2z)/(1.0d0+ro2z*(1.0d0+alp*this%hhz));
    
    if (this%hwEvaE>0) then
       im=this%hwEvaP+this%hwEvaE+1
    else
       im=this%hwEvaP
    endif
    allocate(temp(0:im));
    !=========================== EY == X BOUNDARIES ============================
    !EYFace X=0 -----
    do k=1,this%Nz-1
       do j=2,this%Ny-1
          select case(this%dxmode)
          case(1)
             de = (this%EfOld%Y(this%Nx, j, k) - this%EfOld%Y(this%Nx-1, j, k))/(this%hhx);
          case(2)
             de = (3.0d0*this%EfOld%Y(this%Nx, j, k) - 4.0d0*this%EfOld%Y(this%Nx-1, j, k) + this%EfOld%Y(this%Nx-2, j, k) )/(2.0d0*this%hhx);
          end select
          do v=0,im                  ! Vectors formation
             this%hwF(v) = this%hwYOld(v)%X(this%Nx, j, k);
             this%hwFold(v) = this%hwYold2(v)%X(this%Nx, j, k);
             df1 = (this%hwYOld(v)%X(this%Nx,j+1,k)-2.0d0*this%hwYOld(v)%X(this%Nx,j,k)+this%hwYOld(v)%X(this%Nx,j-1,k))/(this%hhy**2);    ! central in Y
             df2 = (this%hwYOld(v)%X(this%Nx,j,k+1)-2.0d0*this%hwYOld(v)%X(this%Nx,j,k)+this%hwYOld(v)%X(this%Nx,j,k-1))/(this%hhz**2);    ! central in Z 
             this%hwGF(v) = df1+df2;
          enddo

          do v=0,im
             temp(v) = 0.0d0;
          enddo
          ! B*Y^(n)
          do v=0,im
             tp = 0.0d0;
             do w=0,im
                tp = tp - this%hwB(v,w)*this%hwF(w);
             enddo
             temp(v) = temp(v)+tp;
          enddo
          ! C*Y^(n-1)
          do v=0,im
             tp = 0.0d0;
             do w=0,im
                tp = tp - this%hwC(v,w)*this%hwFold(w);
             enddo
             temp(v) = temp(v)+tp;
          enddo
          ! D*GY^(n)
          do v=0,im
             tp = 0.0d0;
             do w=0,im
               tp = tp + this%hwD(v,w)*this%hwGF(w);
             enddo
             temp(v) = temp(v)+tp;
          enddo
          !add dx_u
          select case(this%dtmode)
          case(1)
             temp(0) = temp(0)-1.0d0*cc*this%ht*de;
          case(2)
             temp(0) = temp(0)-2.0d0*cc*this%ht*de;
          end select
          ! A^(-1)*()
          do v=0,im
             tp = 0.0d0;
             do w=0,im
                tp=tp + this%hwA(v,w)*temp(w);
             enddo
             this%hwY(v)%X(this%Nx, j, k) = tp;
          enddo
          this%Ef%Y(this%Nx, j, k) = this%hwY(0)%X(this%Nx, j, k);
       enddo
    enddo

    WRITE(*,*) 'im', im
    
    do j=2,this%Ny-1
       do v=0,im  ! Vectors formation
          !k=0
          this%hwY(v)%X(this%Nx, j, 0) = -(gz1+gz2)*this%hwY(v)%X(this%Nx, j, 1)-gz1*gz2*this%hwY(v)%X(this%Nx, j, 2)  &
          +(gz1+gz2)*this%hwYOld(v)%X(this%Nx, j, 0)+((betaz+1.0d0)+2.0d0*gz1*gz2)*this%hwYOld(v)%X(this%Nx, j, 1)+(gz1+gz2*betaz)*this%hwYOld(v)%X(this%Nx, j, 2)  &
          -gz1*gz2*this%hwYOld2(v)%X(this%Nx, j, 0)-(gz1+gz2*betaz)*this%hwYOld2(v)%X(this%Nx, j, 1)-betaz*this%hwYOld2(v)%X(this%Nx, j, 2);
          !k=Nz         
          this%hwY(v)%X(this%Nx, j, this%Nz) = -(gz1+gz2)*this%hwY(v)%X(this%Nx, j, this%Nz-1)-gz1*gz2*this%hwY(v)%X(this%Nx, j, this%Nz-2)  &
          +(gz1+gz2)*this%hwYOld(v)%X(this%Nx, j, this%Nz)+((betaz+1.0d0)+2.0d0*gz1*gz2)*this%hwYOld(v)%X(this%Nx, j, this%Nz-1)+(gz1+gz2*betaz)*this%hwYOld(v)%X(this%Nx, j, this%Nz-2)  &
          -gz1*gz2*this%hwYOld2(v)%X(this%Nx, j, this%Nz)-(gz1+gz2*betaz)*this%hwYOld2(v)%X(this%Nx, j, this%Nz-1)-betaz*this%hwYOld2(v)%X(this%Nx, j, this%Nz-2);         
       enddo
       this%Ef%Y(this%Nx, j, 0) = this%hwY(0)%X(this%Nx, j, 0);
       this%Ef%Y(this%Nx, j, this%Nz) = this%hwY(0)%X(this%Nx, j, this%Nz);
    enddo   

    do k=1,this%Nz-1
       do v=0,im  ! Vectors formation
          !j=1
          this%hwY(v)%X(this%Nx, 1, k) = -(gy1+gy2)*this%hwY(v)%X(this%Nx, 2, k)-gy1*gy2*this%hwY(v)%X(this%Nx, 3, k)  &
          +(gy1+gy2)*this%hwYOld(v)%X(this%Nx, 1, k)+((betay+1.0d0)+2.0d0*gy1*gy2)*this%hwYOld(v)%X(this%Nx, 2, k)+(gy1+gy2*betay)*this%hwYOld(v)%X(this%Nx, 3, k)  &
          -gy1*gy2*this%hwYOld2(v)%X(this%Nx, 1, k)-(gy1+gy2*betay)*this%hwYOld2(v)%X(this%Nx, 2, k)-betay*this%hwYOld2(v)%X(this%Nx, 3, k);
         !j=Ny
          this%hwY(v)%X(this%Nx, this%Ny, k) = -(gy1+gy2)*this%hwY(v)%X(this%Nx, this%Ny-1, k)-gy1*gy2*this%hwY(v)%X(this%Nx, this%Ny-2, k)  &
          +(gy1+gy2)*this%hwYOld(v)%X(this%Nx, this%Ny, k)+((betay+1.0d0)+2.0d0*gy1*gy2)*this%hwYOld(v)%X(this%Nx, this%Ny-1, k)+(gy1+gy2*betay)*this%hwYOld(v)%X(this%Nx, this%Ny-2, k)  &
          -gy1*gy2*this%hwYOld2(v)%X(this%Nx, this%Ny, k)-(gy1+gy2*betay)*this%hwYOld2(v)%X(this%Nx, this%Ny-1, k)-betay*this%hwYOld2(v)%X(this%Nx, this%Ny-2, k);
       enddo
       this%Ef%Y(this%Nx, 1, k) = this%hwY(0)%X(this%Nx, 1, k);
       this%Ef%Y(this%Nx, this%Ny, k) = this%hwY(0)%X(this%Nx, this%Ny, k);
    enddo

    do v=0,im  ! Vectors formation
       !j=1 k=0
          tp1 = -(gy1+gy2)*this%hwY(v)%X(this%Nx, 2, 0)-gy1*gy2*this%hwY(v)%X(this%Nx, 3, 0)  &
          +(gy1+gy2)*this%hwYOld(v)%X(this%Nx, 1, 0)+((betay+1.0d0)+2.0d0*gy1*gy2)*this%hwYOld(v)%X(this%Nx, 2, 0)+(gy1+gy2*betay)*this%hwYOld(v)%X(this%Nx, 3, 0)  &
          -gy1*gy2*this%hwYOld2(v)%X(this%Nx, 1, 0)-(gy1+gy2*betay)*this%hwYOld2(v)%X(this%Nx, 2, 0)-betay*this%hwYOld2(v)%X(this%Nx, 3, 0);
          tp2 = -(gz1+gz2)*this%hwY(v)%X(this%Nx, 1, 1)-gz1*gz2*this%hwY(v)%X(this%Nx, 1, 2)  &
          +(gz1+gz2)*this%hwYOld(v)%X(this%Nx, 1, 0)+((betaz+1.0d0)+2.0d0*gz1*gz2)*this%hwYOld(v)%X(this%Nx, 1, 1)+(gz1+gz2*betaz)*this%hwYOld(v)%X(this%Nx, 1, 2)  &
          -gz1*gz2*this%hwYOld2(v)%X(this%Nx, 1, 0)-(gz1+gz2*betaz)*this%hwYOld2(v)%X(this%Nx, 1, 1)-betaz*this%hwYOld2(v)%X(this%Nx, 1, 2);
          this%hwY(v)%X(this%Nx, 1, 0) = (tp1+tp2)/2.0d0;
       !j=1 k=Nz         
          tp1 = -(gy1+gy2)*this%hwY(v)%X(this%Nx, 2, this%Nz)-gy1*gy2*this%hwY(v)%X(this%Nx, 3, this%Nz)  &
          +(gy1+gy2)*this%hwYOld(v)%X(this%Nx, 1, this%Nz)+((betay+1.0d0)+2.0d0*gy1*gy2)*this%hwYOld(v)%X(this%Nx, 2, this%Nz)+(gy1+gy2*betay)*this%hwYOld(v)%X(this%Nx, 3, this%Nz)  &
          -gy1*gy2*this%hwYOld2(v)%X(this%Nx, 1, this%Nz)-(gy1+gy2*betay)*this%hwYOld2(v)%X(this%Nx, 2, this%Nz)-betay*this%hwYOld2(v)%X(this%Nx, 3, this%Nz);
          tp2 = -(gz1+gz2)*this%hwY(v)%X(this%Nx, 1, this%Nz-1)-gz1*gz2*this%hwY(v)%X(this%Nx, 1, this%Nz-2)  &
          +(gz1+gz2)*this%hwYOld(v)%X(this%Nx, 1, this%Nz)+((betaz+1.0d0)+2.0d0*gz1*gz2)*this%hwYOld(v)%X(this%Nx, 1, this%Nz-1)+(gz1+gz2*betaz)*this%hwYOld(v)%X(this%Nx, 1, this%Nz-2)  &
          -gz1*gz2*this%hwYOld2(v)%X(this%Nx, 1, this%Nz)-(gz1+gz2*betaz)*this%hwYOld2(v)%X(this%Nx, 1, this%Nz-1)-betaz*this%hwYOld2(v)%X(this%Nx, 1, this%Nz-2);         
          this%hwY(v)%X(this%Nx, 1, this%Nz) = (tp1+tp2)/2.0d0;
       !j=Ny k=0
          tp1 = -(gy1+gy2)*this%hwY(v)%X(this%Nx, this%Ny-1, 0)-gy1*gy2*this%hwY(v)%X(this%Nx, this%Ny-2, 0)  &
          +(gy1+gy2)*this%hwYOld(v)%X(this%Nx, this%Ny, 0)+((betay+1.0d0)+2.0d0*gy1*gy2)*this%hwYOld(v)%X(this%Nx, this%Ny-1, 0)+(gy1+gy2*betay)*this%hwYOld(v)%X(this%Nx, this%Ny-2, 0)  &
          -gy1*gy2*this%hwYOld2(v)%X(this%Nx, this%Ny, 0)-(gy1+gy2*betay)*this%hwYOld2(v)%X(this%Nx, this%Ny-1, 0)-betay*this%hwYOld2(v)%X(this%Nx, this%Ny-2, 0);
          tp2 = -(gz1+gz2)*this%hwY(v)%X(this%Nx, this%Ny, 1)-gz1*gz2*this%hwY(v)%X(this%Nx, this%Ny, 2)  &
          +(gz1+gz2)*this%hwYOld(v)%X(this%Nx, this%Ny, 0)+((betaz+1.0d0)+2.0d0*gz1*gz2)*this%hwYOld(v)%X(this%Nx, this%Ny, 1)+(gz1+gz2*betaz)*this%hwYOld(v)%X(this%Nx, this%Ny, 2)  &
          -gz1*gz2*this%hwYOld2(v)%X(this%Nx, this%Ny, 0)-(gz1+gz2*betaz)*this%hwYOld2(v)%X(this%Nx, this%Ny, 1)-betaz*this%hwYOld2(v)%X(this%Nx, this%Ny, 2);
          this%hwY(v)%X(this%Nx, this%Ny, 0) = (tp1+tp2)/2.0d0;
       !j=Ny k=Nz
          tp1 = -(gy1+gy2)*this%hwY(v)%X(this%Nx, this%Ny-1, this%Nz)-gy1*gy2*this%hwY(v)%X(this%Nx, this%Ny-2, this%Nz)  &
          +(gy1+gy2)*this%hwYOld(v)%X(this%Nx, this%Ny, this%Nz)+((betay+1.0d0)+2.0d0*gy1*gy2)*this%hwYOld(v)%X(this%Nx, this%Ny-1, this%Nz)+(gy1+gy2*betay)*this%hwYOld(v)%X(this%Nx, this%Ny-2, this%Nz)  &
          -gy1*gy2*this%hwYOld2(v)%X(this%Nx, this%Ny, this%Nz)-(gy1+gy2*betay)*this%hwYOld2(v)%X(this%Nx, this%Ny-1, this%Nz)-betay*this%hwYOld2(v)%X(this%Nx, this%Ny-2, this%Nz);
          tp2 = -(gz1+gz2)*this%hwY(v)%X(this%Nx, this%Ny, this%Nz-1)-gz1*gz2*this%hwY(v)%X(this%Nx, this%Ny, this%Nz-2)  &
          +(gz1+gz2)*this%hwYOld(v)%X(this%Nx, this%Ny, this%Nz)+((betaz+1.0d0)+2.0d0*gz1*gz2)*this%hwYOld(v)%X(this%Nx, this%Ny, this%Nz-1)+(gz1+gz2*betaz)*this%hwYOld(v)%X(this%Nx, this%Ny, this%Nz-2)  &
          -gz1*gz2*this%hwYOld2(v)%X(this%Nx, this%Ny, this%Nz)-(gz1+gz2*betaz)*this%hwYOld2(v)%X(this%Nx, this%Ny, this%Nz-1)-betaz*this%hwYOld2(v)%X(this%Nx, this%Ny, this%Nz-2);
          this%hwY(v)%X(this%Nx, this%Ny, this%Nz) = (tp1+tp2)/2.0d0;
       enddo       
    this%Ef%Y(this%Nx, 1, 0) = this%hwY(0)%X(this%Nx, 1, 0);
    this%Ef%Y(this%Nx, 1, this%Nz) = this%hwY(0)%X(this%Nx, 1, this%Nz);
    this%Ef%Y(this%Nx, this%Ny, 0) = this%hwY(0)%X(this%Nx, this%Ny, 0);
    this%Ef%Y(this%Nx, this%Ny, this%Nz) = this%hwY(0)%X(this%Nx, this%Ny, this%Nz);     
    deallocate(temp);
  end subroutine HWEvaEyBetzMittraUpdate


!HWEvaEySuperUpdate------------------------------------------------------------------
  subroutine HWEvaEySuperUpdate(this,t)
    !Fill Ey_x=0 boundary with ABC Hagstrom-Warburton, HW-ABC for boundaries of aux functions
    class(tProblem) :: this
    integer :: t;
    integer :: j ,k, v, va, w, wa, im, ia
    real*8 :: de, df1, df2, tp;
    character*8 :: trans;
    real*8, allocatable :: temp(:), tempa(:);
    real*8 :: gx, gy, gz, gx1, gx2, gy1, gy2, gz1, gz2;
    write(*,*) 'Do Super 2 edges';
    gx1 = (1.0d0-cc*this%ht/(cos(this%hig1/180.0d0*PI)*this%hhx))/(1.0d0+cc*this%ht/(cos(this%hig1/180.0d0*PI)*this%hhx));
    gx2 = (1.0d0-cc*this%ht/(cos(this%hig2/180.0d0*PI)*this%hhx))/(1.0d0+cc*this%ht/(cos(this%hig2/180.0d0*PI)*this%hhx));
    gy1 = (1.0d0-cc*this%ht/(cos(this%hig1/180.0d0*PI)*this%hhy))/(1.0d0+cc*this%ht/(cos(this%hig1/180.0d0*PI)*this%hhy));
    gy2 = (1.0d0-cc*this%ht/(cos(this%hig2/180.0d0*PI)*this%hhy))/(1.0d0+cc*this%ht/(cos(this%hig2/180.0d0*PI)*this%hhy));
    gz1 = (1.0d0-cc*this%ht/(cos(this%hig1/180.0d0*PI)*this%hhx))/(1.0d0+cc*this%ht/(cos(this%hig1/180.0d0*PI)*this%hhz));
    gz2 = (1.0d0-cc*this%ht/(cos(this%hig2/180.0d0*PI)*this%hhx))/(1.0d0+cc*this%ht/(cos(this%hig2/180.0d0*PI)*this%hhz));

    gx = (1.0d0-cc*this%ht/this%hhx)/(1.0d0+cc*this%ht/this%hhx);
    gy = (1.0d0-cc*this%ht/this%hhy)/(1.0d0+cc*this%ht/this%hhy);
    gz = (1.0d0-cc*this%ht/this%hhz)/(1.0d0+cc*this%ht/this%hhz);
    if (this%hwEvaE>0) then
       im=this%hwEvaP+this%hwEvaE+1
    else
       im=this%hwEvaP
    endif
    allocate(temp(0:im));

    if (this%hwEvaEa>0) then
       ia=this%hwEvaPa+this%hwEvaEa+1
    else
       ia=this%hwEvaPa
    endif
    allocate(tempa(0:ia));
    
    !=========================== EY == X BOUNDARIES ============================
    !EYFace X=0 -----

    ! FILL FACE
    do k=1,this%Nz-1
       do j=2,this%Ny-1
          select case(this%dxmode)
          case(1)
             de = (this%EfOld%Y(this%Nx, j, k) - this%EfOld%Y(this%Nx-1, j, k))/(this%hhx);
          case(2)
             de = (3.0d0*this%EfOld%Y(this%Nx, j, k) - 4.0d0*this%EfOld%Y(this%Nx-1, j, k) + this%EfOld%Y(this%Nx-2, j, k) )/(2.0d0*this%hhx);
          end select
          do v=0,im                  ! Vectors formation
             this%hwF(v) = this%hwYOld(v)%X(this%Nx, j, k);
             this%hwFold(v) = this%hwYold2(v)%X(this%Nx, j, k);
             df1 = (this%hwYOld(v)%X(this%Nx,j+1,k)-2.0d0*this%hwYOld(v)%X(this%Nx,j,k)+this%hwYOld(v)%X(this%Nx,j-1,k))/(this%hhy**2);    ! central in Y
             df2 = (this%hwYOld(v)%X(this%Nx,j,k+1)-2.0d0*this%hwYOld(v)%X(this%Nx,j,k)+this%hwYOld(v)%X(this%Nx,j,k-1))/(this%hhz**2);    ! central in Z 
             this%hwGF(v) = df1+df2;                   
          enddo
          do v=0,im
             temp(v) = 0.0d0;
          enddo
          ! B*Y^(n)
          do v=0,im
             tp = 0.0d0;
             do w=0,im
                tp = tp - this%hwB(v,w)*this%hwF(w);
             enddo
             temp(v) = temp(v)+tp;
          enddo
          ! C*Y^(n-1)
          do v=0,im
             tp = 0.0d0;
             do w=0,im
                tp = tp - this%hwC(v,w)*this%hwFold(w);
             enddo
             temp(v) = temp(v)+tp;
          enddo
          ! D*GY^(n)
          do v=0,im
             tp = 0.0d0;
             do w=0,im
               tp = tp + this%hwD(v,w)*this%hwGF(w);
             enddo
             temp(v) = temp(v)+tp;
          enddo
          !add dx_u
          select case(this%dtmode)
          case(1)
             temp(0) = temp(0)-1.0d0*cc*this%ht*de;
          case(2)
             temp(0) = temp(0)-2.0d0*cc*this%ht*de;
          end select
          ! A^(-1)*()
          do v=0,im
             tp = 0.0d0;
             do w=0,im
                tp=tp + this%hwA(v,w)*temp(w);
             enddo
             this%hwY(v)%X(this%Nx, j, k) = tp;
          enddo
          this%Ef%Y(this%Nx, j, k) = this%hwY(0)%X(this%Nx, j, k);
       enddo
    enddo

    ! FILL X=NX, Z=0 Z=NZ EDGES
    do j=2,this%Ny-1

       do v=1,im ! For every Aux Finction except the first one
          ! Z=0 EDGE
          select case(this%dxmode)
          case(1)
             de = (this%hwYOld(v)%X(this%Nx, j, 0) - this%hwYOld(v)%X(this%Nx, j, 1))/(this%hhz);
          case(2)
             de = (3.0d0*this%hwYOld(v)%X(this%Nx, j, 0) - 4.0d0*this%hwYOld(v)%X(this%Nx, j, 1) + this%hwYOld(v)%X(this%Nx, j, 2) )/(2.0d0*this%hhz);
          end select
          do va=0,ia  ! Vectors formation
             this%hwFA(va) = this%hwYEdgeNxOld(v,va)%fbottom(j);
             this%hwFAold(va) = this%hwYEdgeNxOld2(v,va)%fbottom(j);
             this%hwGFA(va) = (this%hwYEdgeNxOld(v,va)%fbottom(j+1)-2.0d0*this%hwYEdgeNxOld(v,va)%fbottom(j)+this%hwYEdgeNxOld(v,va)%fbottom(j-1))/(this%hhy**2);    ! central in Y
          enddo
          do va=0,ia
             tempa(va) = 0.0d0;
          enddo
          ! Ba*Y^(n)
          do va=0,ia
             tp = 0.0d0;
             do wa=0,ia
                tp = tp - this%hwBa(va,wa)*this%hwFA(wa);
             enddo
             tempa(va) = tempa(va)+tp;
          enddo
          ! Ca*Y^(n-1)
          do va=0,ia
             tp = 0.0d0;
             do wa=0,ia
                tp = tp - this%hwCa(va,wa)*this%hwFAold(wa);
             enddo
             tempa(va) = tempa(va)+tp;
          enddo
          ! Da*GY^(n)
          do va=0,ia
             tp = 0.0d0;
             do wa=0,ia
               tp = tp + this%hwDa(va,wa)*this%hwGFA(wa);
             enddo
             tempa(va) = tempa(va)+tp;
          enddo
          !add dx_u
          select case(this%dtmode)
          case(1)
             tempa(0) = tempa(0)-1.0d0*cc*this%ht*de;
          case(2)
             tempa(0) = tempa(0)-2.0d0*cc*this%ht*de;
          end select
          ! Aa^(-1)*()
          do va=0,ia
             tp = 0.0d0;
             do wa=0,ia
                tp=tp + this%hwAa(va,wa)*tempa(wa);
             enddo
             this%hwYEdgeNx(v,va)%fbottom(j) = tp;
          enddo     
          this%hwY(v)%X(this%Nx, j, 0) = this%hwYEdgeNx(v,0)%fbottom(j);
          
          ! Z=NZ EDGE
          select case(this%dxmode)
          case(1)
             de = (this%hwYOld(v)%X(this%Nx, j, this%Nz) - this%hwYOld(v)%X(this%Nx, j, this%Nz-1))/(this%hhz);
          case(2)
             de = (3.0d0*this%hwYOld(v)%X(this%Nx, j, this%Nz) - 4.0d0*this%hwYOld(v)%X(this%Nx, j, this%Nz-1) + this%hwYOld(v)%X(this%Nx, j, this%Nz-2) )/(2.0d0*this%hhz);
          end select
          do va=0,ia  ! Vectors formation
             this%hwFA(va) = this%hwYEdgeNxOld(v,va)%fup(j);
             this%hwFAold(va) = this%hwYEdgeNxOld2(v,va)%fup(j);
             this%hwGFA(va) = (this%hwYEdgeNxOld(v,va)%fup(j+1)-2.0d0*this%hwYEdgeNxOld(v,va)%fup(j)+this%hwYEdgeNxOld(v,va)%fup(j-1))/(this%hhy**2);    ! central in Y
          enddo
          do va=0,ia
             tempa(va) = 0.0d0;
          enddo
          ! Ba*Y^(n)
          do va=0,ia
             tp = 0.0d0;
             do wa=0,ia
                tp = tp - this%hwBa(va,wa)*this%hwFA(wa);
             enddo
             tempa(va) = tempa(va)+tp;
          enddo
          ! Ca*Y^(n-1)
          do va=0,ia
             tp = 0.0d0;
             do wa=0,ia
                tp = tp - this%hwCa(va,wa)*this%hwFAold(wa);
             enddo
             tempa(va) = tempa(va)+tp;
          enddo
          ! Da*GY^(n)
          do va=0,ia
             tp = 0.0d0;
             do wa=0,ia
               tp = tp + this%hwDa(va,wa)*this%hwGFA(wa);
             enddo
             tempa(va) = tempa(va)+tp;
          enddo
          !add dx_u
          select case(this%dtmode)
          case(1)
             tempa(0) = tempa(0)-1.0d0*cc*this%ht*de;
          case(2)
             tempa(0) = tempa(0)-2.0d0*cc*this%ht*de;
          end select
          ! Aa^(-1)*()
          do va=0,ia
             tp = 0.0d0;
             do wa=0,ia
                tp=tp + this%hwAa(va,wa)*tempa(wa);
             enddo
             this%hwYEdgeNx(v,va)%fup(j) = tp;
          enddo     
          this%hwY(v)%X(this%Nx, j, this%Nz) = this%hwYEdgeNx(v,0)%fup(j);
  
       enddo
    enddo

   ! MAIN FILL X=NX EDGE
    do j=2,this%Ny-1
       !Z=0 EDGE
       this%hwY(0)%X(this%Nx, j, 0) = this%Source%getEpoint(this%xi(this%Nx), this%yi05(j),  this%zi(0), this%Ti(t), this%Nx, j, 0, 2);
       this%Ef%Y(this%Nx, j, 0) = this%hwY(0)%X(this%Nx, j, 0);
       !Z=Nz EDGE
       this%hwY(0)%X(this%Nx, j, this%Nz) = this%Source%getEpoint(this%xi(this%Nx), this%yi05(j),  this%zi(this%Nz), this%Ti(t), this%Nx, j, this%Nz, 2);
       this%Ef%Y(this%Nx, j, this%Nz) = this%hwY(0)%X(this%Nx, j, this%Nz);
    enddo   

    ! FILL ENDS of THE EDGE
    do v=1,im     ! For every Aux Function except the field
       do va=0,ia ! For every AuxA Function fill with Sommerfeld
          ! z=0 j=1   
          this%hwYEdgeNx(v,va)%fbottom(1) = gy*this%hwYEdgeNx(v,va)%fbottom(1) - gy*this%hwYEdgeNx(v,va)%fbottom(2) + this%hwYEdgeNxOld(v,va)%fbottom(2);
          ! z=0 j=Ny
          this%hwYEdgeNx(v,va)%fbottom(this%Ny) = gy*this%hwYEdgeNx(v,va)%fbottom(this%Ny) - gy*this%hwYEdgeNx(v,va)%fbottom(this%Ny-1) + this%hwYEdgeNxOld(v,va)%fbottom(this%Ny-1);
          ! z=Nz j=1   
          this%hwYEdgeNx(v,va)%fup(1) = gy*this%hwYEdgeNx(v,va)%fup(1) - gy*this%hwYEdgeNx(v,va)%fup(2) + this%hwYEdgeNxOld(v,va)%fup(2);
          ! z=Nz j=Ny
          this%hwYEdgeNx(v,va)%fup(this%Ny) = gy*this%hwYEdgeNx(v,va)%fup(this%Ny) - gy*this%hwYEdgeNx(v,va)%fup(this%Ny-1) + this%hwYEdgeNxOld(v,va)%fup(this%Ny-1);
       enddo
       this%hwY(v)%X(this%Nx, 1, 0) = this%hwYEdgeNx(v,0)%fbottom(1);
       this%hwY(v)%X(this%Nx, this%Ny, 0) = this%hwYEdgeNx(v,0)%fbottom(this%Ny);
       this%hwY(v)%X(this%Nx, 1, this%Nz) = this%hwYEdgeNx(v,0)%fup(1);
       this%hwY(v)%X(this%Nx, this%Ny, this%Nz) = this%hwYEdgeNx(v,0)%fup(this%Ny);
    enddo
   
    this%hwY(0)%X(this%Nx, 1, 0) = this%Source%getEpoint(this%xi(this%Nx), this%yi05(1),  this%zi(0), this%Ti(t), this%Nx, 1, 0, 2);
    this%hwY(0)%X(this%Nx, this%Ny, 0) = this%Source%getEpoint(this%xi(this%Nx), this%yi05(this%Ny),  this%zi(0), this%Ti(t), this%Nx, this%Ny, 0, 2);
    this%Ef%Y(this%Nx, 1, 0) = this%hwY(0)%X(this%Nx, 1, 0);
    this%Ef%Y(this%Nx, this%Ny, 0) = this%hwY(0)%X(this%Nx, this%Ny, 0);

    this%hwY(0)%X(this%Nx, 1, this%Nz) = this%Source%getEpoint(this%xi(this%Nx), this%yi05(1),  this%zi(this%Nz), this%Ti(t), this%Nx, 1, this%Nz, 2);
    this%hwY(0)%X(this%Nx, this%Ny, this%Nz) = this%Source%getEpoint(this%xi(this%Nx), this%yi05(this%Ny),  this%zi(this%Nz), this%Ti(t), this%Nx, this%Ny, this%Nz, 2);
    this%Ef%Y(this%Nx, 1, this%Nz) = this%hwY(0)%X(this%Nx, 1, this%Nz);
    this%Ef%Y(this%Nx, this%Ny, this%Nz) = this%hwY(0)%X(this%Nx, this%Ny, this%Nz);


    
!    do j=2,this%Ny-1
!       do v=0,im  ! Vectors formation
          !k=0
!          this%hwY(v)%X(this%Nx, j, 0) = -(gz1+gz2)*this%hwY(v)%X(this%Nx, j, 1)-gz1*gz2*this%hwY(v)%X(this%Nx, j, 2)  &
!          +(gz1+gz2)*this%hwYOld(v)%X(this%Nx, j, 0)+(2.0d0+2.0d0*gz1*gz2)*this%hwYOld(v)%X(this%Nx, j, 1)+(gz1+gz2)*this%hwYOld(v)%X(this%Nx, j, 2)  &
!          -gz1*gz2*this%hwYOld2(v)%X(this%Nx, j, 0)-(gz1+gz2)*this%hwYOld2(v)%X(this%Nx, j, 1)-this%hwYOld2(v)%X(this%Nx, j, 2);
          !k=Nz         
!          this%hwY(v)%X(this%Nx, j, this%Nz) = -(gz1+gz2)*this%hwY(v)%X(this%Nx, j, this%Nz-1)-gz1*gz2*this%hwY(v)%X(this%Nx, j, this%Nz-2)  &
!          +(gz1+gz2)*this%hwYOld(v)%X(this%Nx, j, this%Nz)+(2.0d0+2.0d0*gz1*gz2)*this%hwYOld(v)%X(this%Nx, j, this%Nz-1)+(gz1+gz2)*this%hwYOld(v)%X(this%Nx, j, this%Nz-2)  &
!          -gz1*gz2*this%hwYOld2(v)%X(this%Nx, j, this%Nz)-(gz1+gz2)*this%hwYOld2(v)%X(this%Nx, j, this%Nz-1)-this%hwYOld2(v)%X(this%Nx, j, this%Nz-2);
!       enddo
!       this%Ef%Y(this%Nx, j, 0) = this%hwY(0)%X(this%Nx, j, 0);
!       this%Ef%Y(this%Nx, j, this%Nz) = this%hwY(0)%X(this%Nx, j, this%Nz);
!    enddo   

    do k=1,this%Nz-1
       do v=0,im  ! Vectors formation
          !j=1
          this%hwY(v)%X(this%Nx, 1, k) = -(gy1+gy2)*this%hwY(v)%X(this%Nx, 2, k)-gy1*gy2*this%hwY(v)%X(this%Nx, 3, k)  &
          +(gy1+gy2)*this%hwYOld(v)%X(this%Nx, 1, k)+(2.0d0+2.0d0*gy1*gy2)*this%hwYOld(v)%X(this%Nx, 2, k)+(gy1+gy2)*this%hwYOld(v)%X(this%Nx, 3, k)  &
          -gy1*gy2*this%hwYOld2(v)%X(this%Nx, 1, k)-(gy1+gy2)*this%hwYOld2(v)%X(this%Nx, 2, k)-this%hwYOld2(v)%X(this%Nx, 3, k);
         !j=Ny
          this%hwY(v)%X(this%Nx, this%Ny, k) = -(gy1+gy2)*this%hwY(v)%X(this%Nx, this%Ny-1, k)-gy1*gy2*this%hwY(v)%X(this%Nx, this%Ny-2, k)  &
          +(gy1+gy2)*this%hwYOld(v)%X(this%Nx, this%Ny, k)+(2.0d0+2.0d0*gy1*gy2)*this%hwYOld(v)%X(this%Nx, this%Ny-1, k)+(gy1+gy2)*this%hwYOld(v)%X(this%Nx, this%Ny-2, k)  &
          -gy1*gy2*this%hwYOld2(v)%X(this%Nx, this%Ny, k)-(gy1+gy2)*this%hwYOld2(v)%X(this%Nx, this%Ny-1, k)-this%hwYOld2(v)%X(this%Nx, this%Ny-2, k);
       enddo
       this%Ef%Y(this%Nx, 1, k) = this%hwY(0)%X(this%Nx, 1, k);
       this%Ef%Y(this%Nx, this%Ny, k) = this%hwY(0)%X(this%Nx, this%Ny, k);
    enddo

    do v=0,im  ! Vectors formation
       if (v==0) then
       !j=1 k=0
          tp = -(gy1+gy2)*this%hwY(v)%X(this%Nx, 2, 0)-gy1*gy2*this%hwY(v)%X(this%Nx, 3, 0)  &
          +(gy1+gy2)*this%hwYOld(v)%X(this%Nx, 1, 0)+(2.0d0+2.0d0*gy1*gy2)*this%hwYOld(v)%X(this%Nx, 2, 0)+(gy1+gy2)*this%hwYOld(v)%X(this%Nx, 3, 0)  &
          -gy1*gy2*this%hwYOld2(v)%X(this%Nx, 1, 0)-(gy1+gy2)*this%hwYOld2(v)%X(this%Nx, 2, 0)-this%hwYOld2(v)%X(this%Nx, 3, 0);
          this%hwY(v)%X(this%Nx, 1, 0) = (tp -(gz1+gz2)*this%hwY(v)%X(this%Nx, 1, 1)-gz1*gz2*this%hwY(v)%X(this%Nx, 1, 2)  &
          +(gz1+gz2)*this%hwYOld(v)%X(this%Nx, 1, 0)+(2.0d0+2.0d0*gz1*gz2)*this%hwYOld(v)%X(this%Nx, 1, 1)+(gz1+gz2)*this%hwYOld(v)%X(this%Nx, 1, 2)  &
          -gz1*gz2*this%hwYOld2(v)%X(this%Nx, 1, 0)-(gz1+gz2)*this%hwYOld2(v)%X(this%Nx, 1, 1)-this%hwYOld2(v)%X(this%Nx, 1, 2))/2.0d0;
       !j=1 k=Nz
          tp = -(gy1+gy2)*this%hwY(v)%X(this%Nx, 2, this%Nz)-gy1*gy2*this%hwY(v)%X(this%Nx, 3, this%Nz)  &
          +(gy1+gy2)*this%hwYOld(v)%X(this%Nx, 1, this%Nz)+(2.0d0+2.0d0*gy1*gy2)*this%hwYOld(v)%X(this%Nx, 2, this%Nz)+(gy1+gy2)*this%hwYOld(v)%X(this%Nx, 3, this%Nz)  &
          -gy1*gy2*this%hwYOld2(v)%X(this%Nx, 1, this%Nz)-(gy1+gy2)*this%hwYOld2(v)%X(this%Nx, 2, this%Nz)-this%hwYOld2(v)%X(this%Nx, 3, this%Nz);
          this%hwY(v)%X(this%Nx, 1, this%Nz) = (tp -(gz1+gz2)*this%hwY(v)%X(this%Nx, 1, this%Nz-1)-gz1*gz2*this%hwY(v)%X(this%Nx, 1, this%Nz-2)  &
          +(gz1+gz2)*this%hwYOld(v)%X(this%Nx, 1, this%Nz)+(2.0d0+2.0d0*gz1*gz2)*this%hwYOld(v)%X(this%Nx, 1, this%Nz-1)+(gz1+gz2)*this%hwYOld(v)%X(this%Nx, 1, this%Nz-2)  &
          -gz1*gz2*this%hwYOld2(v)%X(this%Nx, 1, this%Nz)-(gz1+gz2)*this%hwYOld2(v)%X(this%Nx, 1, this%Nz-1)-this%hwYOld2(v)%X(this%Nx, 1, this%Nz-2))/2.0d0;
       !j=Ny k=0
          tp = -(gy1+gy2)*this%hwY(v)%X(this%Nx, this%Ny-1, 0)-gy1*gy2*this%hwY(v)%X(this%Nx, this%Ny-2, 0)  &
          +(gy1+gy2)*this%hwYOld(v)%X(this%Nx, this%Ny, 0)+(2.0d0+2.0d0*gy1*gy2)*this%hwYOld(v)%X(this%Nx, this%Ny-1, 0)+(gy1+gy2)*this%hwYOld(v)%X(this%Nx, this%Ny-2, 0)  &
          -gy1*gy2*this%hwYOld2(v)%X(this%Nx, this%Ny, 0)-(gy1+gy2)*this%hwYOld2(v)%X(this%Nx, this%Ny-1, 0)-this%hwYOld2(v)%X(this%Nx, this%Ny-2, 0);
          this%hwY(v)%X(this%Nx, this%Ny, 0) = (tp -(gz1+gz2)*this%hwY(v)%X(this%Nx, this%Ny, 1)-gz1*gz2*this%hwY(v)%X(this%Nx, this%Ny, 2)  &
          +(gz1+gz2)*this%hwYOld(v)%X(this%Nx, this%Ny, 0)+(2.0d0+2.0d0*gz1*gz2)*this%hwYOld(v)%X(this%Nx, this%Ny, 1)+(gz1+gz2)*this%hwYOld(v)%X(this%Nx, this%Ny, 2)  &
          -gz1*gz2*this%hwYOld2(v)%X(this%Nx, this%Ny, 0)-(gz1+gz2)*this%hwYOld2(v)%X(this%Nx, this%Ny, 1)-this%hwYOld2(v)%X(this%Nx, this%Ny, 2))/2.0d0;
       !j=Ny k=Nz
          tp = -(gy1+gy2)*this%hwY(v)%X(this%Nx, this%Ny-1, this%Nz)-gy1*gy2*this%hwY(v)%X(this%Nx, this%Ny-2, this%Nz)  &
          +(gy1+gy2)*this%hwYOld(v)%X(this%Nx, this%Ny, this%Nz)+(2.0d0+2.0d0*gy1*gy2)*this%hwYOld(v)%X(this%Nx, this%Ny-1, this%Nz)+(gy1+gy2)*this%hwYOld(v)%X(this%Nx, this%Ny-2, this%Nz)  &
          -gy1*gy2*this%hwYOld2(v)%X(this%Nx, this%Ny, this%Nz)-(gy1+gy2)*this%hwYOld2(v)%X(this%Nx, this%Ny-1, this%Nz)-this%hwYOld2(v)%X(this%Nx, this%Ny-2, this%Nz);
          this%hwY(v)%X(this%Nx, this%Ny, this%Nz) = (tp - (gz1+gz2)*this%hwY(v)%X(this%Nx, this%Ny, this%Nz-1)-gz1*gz2*this%hwY(v)%X(this%Nx, this%Ny, this%Nz-2)  &
          +(gz1+gz2)*this%hwYOld(v)%X(this%Nx, this%Ny, this%Nz)+(2.0d0+2.0d0*gz1*gz2)*this%hwYOld(v)%X(this%Nx, this%Ny, this%Nz-1)+(gz1+gz2)*this%hwYOld(v)%X(this%Nx, this%Ny, this%Nz-2)  &
          -gz1*gz2*this%hwYOld2(v)%X(this%Nx, this%Ny, this%Nz)-(gz1+gz2)*this%hwYOld2(v)%X(this%Nx, this%Ny, this%Nz-1)-this%hwYOld2(v)%X(this%Nx, this%Ny, this%Nz-2) )/2.0d0;
       endif
    enddo

    this%Ef%Y(this%Nx, 1, 0) = this%hwY(0)%X(this%Nx, 1, 0);
    this%Ef%Y(this%Nx, 1, this%Nz) = this%hwY(0)%X(this%Nx, 1, this%Nz);
    this%Ef%Y(this%Nx, this%Ny, 0) = this%hwY(0)%X(this%Nx, this%Ny, 0);
    this%Ef%Y(this%Nx, this%Ny, this%Nz) = this%hwY(0)%X(this%Nx, this%Ny, this%Nz);   
   
    deallocate(temp);
    deallocate(tempa);
  end subroutine HWEvaEySuperUpdate


!HWEvaEyHigdon3DUpdate------------------------------------------------------------------
  subroutine HWEvaEyHigdon3DUpdate(this, t)
    !Fill Ey_x=0 boundary with ABC Hagstrom-Warburton, Higdon for boundaries of aux functions
    class(tProblem) :: this
    integer :: t;
    integer :: i, j ,k, v, w, im
    real*8 :: de, df1, df2, tp;
    character*8 :: trans;
    real*8, allocatable :: temp(:);
    real*8 :: gx1, gx2, gy1, gy2, gz1, gz2;
    write(*,*) 'Do 3D Higdon';
    gx1 = (1.0d0-cc*this%ht/(cos(this%hig1/180.0d0*PI)*this%hhx))/(1.0d0+cc*this%ht/(cos(this%hig1/180.0d0*PI)*this%hhx));
    gx2 = (1.0d0-cc*this%ht/(cos(this%hig2/180.0d0*PI)*this%hhx))/(1.0d0+cc*this%ht/(cos(this%hig2/180.0d0*PI)*this%hhx));
    gy1 = (1.0d0-cc*this%ht/(cos(this%hig1/180.0d0*PI)*this%hhy))/(1.0d0+cc*this%ht/(cos(this%hig1/180.0d0*PI)*this%hhy));
    gy2 = (1.0d0-cc*this%ht/(cos(this%hig2/180.0d0*PI)*this%hhy))/(1.0d0+cc*this%ht/(cos(this%hig2/180.0d0*PI)*this%hhy));
    gz1 = (1.0d0-cc*this%ht/(cos(this%hig1/180.0d0*PI)*this%hhx))/(1.0d0+cc*this%ht/(cos(this%hig1/180.0d0*PI)*this%hhz));
    gz2 = (1.0d0-cc*this%ht/(cos(this%hig2/180.0d0*PI)*this%hhx))/(1.0d0+cc*this%ht/(cos(this%hig2/180.0d0*PI)*this%hhz));
    if (this%hwEvaE>0) then
       im=this%hwEvaP+this%hwEvaE+1
    else
       im=this%hwEvaP
    endif
    allocate(temp(0:im));
    !=========================== EY == X=Nx BOUNDARIES ============================
    !EYFace X=Nx -----
    do k=1,this%Nz-1
       do j=2,this%Ny-1
          select case(this%dxmode)
          case(1)
             de = (this%EfOld%Y(this%Nx, j, k) - this%EfOld%Y(this%Nx-1, j, k))/(this%hhx);
          case(2)
             de = (3.0d0*this%EfOld%Y(this%Nx, j, k) - 4.0d0*this%EfOld%Y(this%Nx-1, j, k) + this%EfOld%Y(this%Nx-2, j, k) )/(2.0d0*this%hhx);
          end select
          do v=0,im                  ! Vectors formation
             this%hwF(v) = this%hwYOld(v)%X(this%Nx, j, k);
             this%hwFold(v) = this%hwYold2(v)%X(this%Nx, j, k);
             df1 = (this%hwYOld(v)%X(this%Nx,j+1,k)-2.0d0*this%hwYOld(v)%X(this%Nx,j,k)+this%hwYOld(v)%X(this%Nx,j-1,k))/(this%hhy**2);    ! central in Y
             df2 = (this%hwYOld(v)%X(this%Nx,j,k+1)-2.0d0*this%hwYOld(v)%X(this%Nx,j,k)+this%hwYOld(v)%X(this%Nx,j,k-1))/(this%hhz**2);    ! central in Z 
             this%hwGF(v) = df1+df2;
          enddo
          do v=0,im
             temp(v) = 0.0d0;
          enddo
          ! B*Y^(n)
          do v=0,im
             tp = 0.0d0;
             do w=0,im
                tp = tp - this%hwB(v,w)*this%hwF(w);
             enddo
             temp(v) = temp(v)+tp;
          enddo
          ! C*Y^(n-1)
          do v=0,im
             tp = 0.0d0;
             do w=0,im
                tp = tp - this%hwC(v,w)*this%hwFold(w);
             enddo
             temp(v) = temp(v)+tp;
          enddo
          ! D*GY^(n)
          do v=0,im
             tp = 0.0d0;
             do w=0,im
               tp = tp + this%hwD(v,w)*this%hwGF(w);
             enddo
             temp(v) = temp(v)+tp;
          enddo
          !add dx_u
          select case(this%dtmode)
          case(1)
             temp(0) = temp(0)-1.0d0*cc*this%ht*de;
          case(2)
             temp(0) = temp(0)-2.0d0*cc*this%ht*de;
          end select
          ! A^(-1)*()
          do v=0,im
             tp = 0.0d0;
             do w=0,im
                tp=tp + this%hwA(v,w)*temp(w);
             enddo
             this%hwY(v)%X(this%Nx, j, k) = tp;
         enddo
          this%Ef%Y(this%Nx, j, k) = this%hwY(0)%X(this%Nx, j, k);
       enddo
    enddo

    do k=1,this%Nz-1
      do v=0,im  ! Vectors formation
         !j=1
         this%hwY(v)%X(this%Nx, 1, k) = -(gy1+gy2)*this%hwY(v)%X(this%Nx, 2, k)-gy1*gy2*this%hwY(v)%X(this%Nx, 3, k)  &
         +(gy1+gy2)*this%hwYOld(v)%X(this%Nx, 1, k)+(2.0d0+2.0d0*gy1*gy2)*this%hwYOld(v)%X(this%Nx, 2, k)+(gy1+gy2)*this%hwYOld(v)%X(this%Nx, 3, k)  &
         -gy1*gy2*this%hwYOld2(v)%X(this%Nx, 1, k)-(gy1+gy2)*this%hwYOld2(v)%X(this%Nx, 2, k)-this%hwYOld2(v)%X(this%Nx, 3, k);
         !j=Ny
         this%hwY(v)%X(this%Nx, this%Ny, k) = -(gy1+gy2)*this%hwY(v)%X(this%Nx, this%Ny-1, k)-gy1*gy2*this%hwY(v)%X(this%Nx, this%Ny-2, k)  &
         +(gy1+gy2)*this%hwYOld(v)%X(this%Nx, this%Ny, k)+(2.0d0+2.0d0*gy1*gy2)*this%hwYOld(v)%X(this%Nx, this%Ny-1, k)+(gy1+gy2)*this%hwYOld(v)%X(this%Nx, this%Ny-2, k)  &
         -gy1*gy2*this%hwYOld2(v)%X(this%Nx, this%Ny, k)-(gy1+gy2)*this%hwYOld2(v)%X(this%Nx, this%Ny-1, k)-this%hwYOld2(v)%X(this%Nx, this%Ny-2, k);
      enddo
      this%Ef%Y(this%Nx, 1, k) = this%hwY(0)%X(this%Nx, 1, k);
      this%Ef%Y(this%Nx, this%Ny, k) = this%hwY(0)%X(this%Nx, this%Ny, k);
    enddo
    
    do j=2,this%Ny-1
       do v=0,im  ! Vectors formation
          !k=0
          this%hwY(v)%X(this%Nx, j, 0) = -(gz1+gz2)*this%hwY(v)%X(this%Nx, j, 1)-gz1*gz2*this%hwY(v)%X(this%Nx, j, 2)  &
          +(gz1+gz2)*this%hwYOld(v)%X(this%Nx, j, 0)+(2.0d0+2.0d0*gz1*gz2)*this%hwYOld(v)%X(this%Nx, j, 1)+(gz1+gz2)*this%hwYOld(v)%X(this%Nx, j, 2)  &
          -gz1*gz2*this%hwYOld2(v)%X(this%Nx, j, 0)-(gz1+gz2)*this%hwYOld2(v)%X(this%Nx, j, 1)-this%hwYOld2(v)%X(this%Nx, j, 2);
          !k=Nz         
!          this%hwY(v)%X(this%Nx, j, this%Nz) = -(gz1+gz2)*this%hwY(v)%X(this%Nx, j, this%Nz-1)-gz1*gz2*this%hwY(v)%X(this%Nx, j, this%Nz-2)  &
!          +(gz1+gz2)*this%hwYOld(v)%X(this%Nx, j, this%Nz)+(2.0d0+2.0d0*gz1*gz2)*this%hwYOld(v)%X(this%Nx, j, this%Nz-1)+(gz1+gz2)*this%hwYOld(v)%X(this%Nx, j, this%Nz-2)  &
!         -gz1*gz2*this%hwYOld2(v)%X(this%Nx, j, this%Nz)-(gz1+gz2)*this%hwYOld2(v)%X(this%Nx, j, this%Nz-1)-this%hwYOld2(v)%X(this%Nx, j, this%Nz-2);
       enddo
       this%Ef%Y(this%Nx, j, 0) = this%hwY(0)%X(this%Nx, j, 0);
!       this%Ef%Y(this%Nx, j, this%Nz) = this%hwY(0)%X(this%Nx, j, this%Nz);
    enddo   

    do v=0,im  ! Vectors formation
       !j=1 k=0
          tp = -(gy1+gy2)*this%hwY(v)%X(this%Nx, 2, 0)-gy1*gy2*this%hwY(v)%X(this%Nx, 3, 0)  &
          +(gy1+gy2)*this%hwYOld(v)%X(this%Nx, 1, 0)+(2.0d0+2.0d0*gy1*gy2)*this%hwYOld(v)%X(this%Nx, 2, 0)+(gy1+gy2)*this%hwYOld(v)%X(this%Nx, 3, 0)  &
          -gy1*gy2*this%hwYOld2(v)%X(this%Nx, 1, 0)-(gy1+gy2)*this%hwYOld2(v)%X(this%Nx, 2, 0)-this%hwYOld2(v)%X(this%Nx, 3, 0);
          this%hwY(v)%X(this%Nx, 1, 0) = (tp -(gz1+gz2)*this%hwY(v)%X(this%Nx, 1, 1)-gz1*gz2*this%hwY(v)%X(this%Nx, 1, 2)  &
          +(gz1+gz2)*this%hwYOld(v)%X(this%Nx, 1, 0)+(2.0d0+2.0d0*gz1*gz2)*this%hwYOld(v)%X(this%Nx, 1, 1)+(gz1+gz2)*this%hwYOld(v)%X(this%Nx, 1, 2)  &
          -gz1*gz2*this%hwYOld2(v)%X(this%Nx, 1, 0)-(gz1+gz2)*this%hwYOld2(v)%X(this%Nx, 1, 1)-this%hwYOld2(v)%X(this%Nx, 1, 2))/2.0d0;
       !j=1 k=Nz
!          tp = -(gy1+gy2)*this%hwY(v)%X(this%Nx, 2, this%Nz)-gy1*gy2*this%hwY(v)%X(this%Nx, 3, this%Nz)  &
!          +(gy1+gy2)*this%hwYOld(v)%X(this%Nx, 1, this%Nz)+(2.0d0+2.0d0*gy1*gy2)*this%hwYOld(v)%X(this%Nx, 2, this%Nz)+(gy1+gy2)*this%hwYOld(v)%X(this%Nx, 3, this%Nz)  &
!          -gy1*gy2*this%hwYOld2(v)%X(this%Nx, 1, this%Nz)-(gy1+gy2)*this%hwYOld2(v)%X(this%Nx, 2, this%Nz)-this%hwYOld2(v)%X(this%Nx, 3, this%Nz);
!          this%hwY(v)%X(this%Nx, 1, this%Nz) = (tp -(gz1+gz2)*this%hwY(v)%X(this%Nx, 1, this%Nz-1)-gz1*gz2*this%hwY(v)%X(this%Nx, 1, this%Nz-2)  &
!          +(gz1+gz2)*this%hwYOld(v)%X(this%Nx, 1, this%Nz)+(2.0d0+2.0d0*gz1*gz2)*this%hwYOld(v)%X(this%Nx, 1, this%Nz-1)+(gz1+gz2)*this%hwYOld(v)%X(this%Nx, 1, this%Nz-2)  &
!          -gz1*gz2*this%hwYOld2(v)%X(this%Nx, 1, this%Nz)-(gz1+gz2)*this%hwYOld2(v)%X(this%Nx, 1, this%Nz-1)-this%hwYOld2(v)%X(this%Nx, 1, this%Nz-2))/2.0d0;
       !j=Ny k=0
          tp = -(gy1+gy2)*this%hwY(v)%X(this%Nx, this%Ny-1, 0)-gy1*gy2*this%hwY(v)%X(this%Nx, this%Ny-2, 0)  &
          +(gy1+gy2)*this%hwYOld(v)%X(this%Nx, this%Ny, 0)+(2.0d0+2.0d0*gy1*gy2)*this%hwYOld(v)%X(this%Nx, this%Ny-1, 0)+(gy1+gy2)*this%hwYOld(v)%X(this%Nx, this%Ny-2, 0)  &
          -gy1*gy2*this%hwYOld2(v)%X(this%Nx, this%Ny, 0)-(gy1+gy2)*this%hwYOld2(v)%X(this%Nx, this%Ny-1, 0)-this%hwYOld2(v)%X(this%Nx, this%Ny-2, 0);
          this%hwY(v)%X(this%Nx, this%Ny, 0) = (tp -(gz1+gz2)*this%hwY(v)%X(this%Nx, this%Ny, 1)-gz1*gz2*this%hwY(v)%X(this%Nx, this%Ny, 2)  &
          +(gz1+gz2)*this%hwYOld(v)%X(this%Nx, this%Ny, 0)+(2.0d0+2.0d0*gz1*gz2)*this%hwYOld(v)%X(this%Nx, this%Ny, 1)+(gz1+gz2)*this%hwYOld(v)%X(this%Nx, this%Ny, 2)  &
          -gz1*gz2*this%hwYOld2(v)%X(this%Nx, this%Ny, 0)-(gz1+gz2)*this%hwYOld2(v)%X(this%Nx, this%Ny, 1)-this%hwYOld2(v)%X(this%Nx, this%Ny, 2))/2.0d0;
       !j=Ny k=Nz
!          tp = -(gy1+gy2)*this%hwY(v)%X(this%Nx, this%Ny-1, this%Nz)-gy1*gy2*this%hwY(v)%X(this%Nx, this%Ny-2, this%Nz)  &
!          +(gy1+gy2)*this%hwYOld(v)%X(this%Nx, this%Ny, this%Nz)+(2.0d0+2.0d0*gy1*gy2)*this%hwYOld(v)%X(this%Nx, this%Ny-1, this%Nz)+(gy1+gy2)*this%hwYOld(v)%X(this%Nx, this%Ny-2, this%Nz)  &
!          -gy1*gy2*this%hwYOld2(v)%X(this%Nx, this%Ny, this%Nz)-(gy1+gy2)*this%hwYOld2(v)%X(this%Nx, this%Ny-1, this%Nz)-this%hwYOld2(v)%X(this%Nx, this%Ny-2, this%Nz);
!          this%hwY(v)%X(this%Nx, this%Ny, this%Nz) = (tp - (gz1+gz2)*this%hwY(v)%X(this%Nx, this%Ny, this%Nz-1)-gz1*gz2*this%hwY(v)%X(this%Nx, this%Ny, this%Nz-2)  &
!          +(gz1+gz2)*this%hwYOld(v)%X(this%Nx, this%Ny, this%Nz)+(2.0d0+2.0d0*gz1*gz2)*this%hwYOld(v)%X(this%Nx, this%Ny, this%Nz-1)+(gz1+gz2)*this%hwYOld(v)%X(this%Nx, this%Ny, this%Nz-2)  &
!          -gz1*gz2*this%hwYOld2(v)%X(this%Nx, this%Ny, this%Nz)-(gz1+gz2)*this%hwYOld2(v)%X(this%Nx, this%Ny, this%Nz-1)-this%hwYOld2(v)%X(this%Nx, this%Ny, this%Nz-2) )/2.0d0;
       enddo
       
    this%Ef%Y(this%Nx, 1, 0) = this%hwY(0)%X(this%Nx, 1, 0);
!    this%Ef%Y(this%Nx, 1, this%Nz) = this%hwY(0)%X(this%Nx, 1, this%Nz);
    this%Ef%Y(this%Nx, this%Ny, 0) = this%hwY(0)%X(this%Nx, this%Ny, 0);
!    this%Ef%Y(this%Nx, this%Ny, this%Nz) = this%hwY(0)%X(this%Nx, this%Ny, this%Nz);
    

    !=========================== EY == Z=Nz BOUNDARIES ============================
    !EYFace Z=Nx -----
    do i=1,this%Nx-1
       do j=2,this%Ny-1
          select case(this%dxmode)
          case(1)
             de = (this%EfOld%Y(i, j, this%Nz) - this%EfOld%Y(i, j, this%Nz-1))/(this%hhz);
          case(2)
             de = (3.0d0*this%EfOld%Y(i, j, this%Nz) - 4.0d0*this%EfOld%Y(i, j, this%Nz-1) + this%EfOld%Y(i, j, this%Nz-2) )/(2.0d0*this%hhz);
          end select
          do v=0,im                  ! Vectors formation
             this%hwF(v) = this%hwYOld(v)%Z(i, j, this%Nz);
             this%hwFold(v) = this%hwYold2(v)%Z(i, j, this%Nz);
             df1 = (this%hwYOld(v)%Z(i,j+1,this%Nz)-2.0d0*this%hwYOld(v)%Z(i,j,this%Nz)+this%hwYOld(v)%Z(i,j-1,this%Nz))/(this%hhy**2);    ! central in Y
             df2 = (this%hwYOld(v)%Z(i+1,j,this%Nz)-2.0d0*this%hwYOld(v)%Z(i,j,this%Nz)+this%hwYOld(v)%Z(i-1,j,this%Nz))/(this%hhx**2);    ! central in X 
             this%hwGF(v) = df1+df2;
          enddo

          do v=0,im
             temp(v) = 0.0d0;
          enddo
          ! B*Y^(n)
          do v=0,im
             tp = 0.0d0;
             do w=0,im
                tp = tp - this%hwB(v,w)*this%hwF(w);
             enddo
             temp(v) = temp(v)+tp;
          enddo
          ! C*Y^(n-1)
          do v=0,im
             tp = 0.0d0;
             do w=0,im
                tp = tp - this%hwC(v,w)*this%hwFold(w);
             enddo
             temp(v) = temp(v)+tp;
          enddo
          ! D*GY^(n)
          do v=0,im
             tp = 0.0d0;
             do w=0,im
               tp = tp + this%hwD(v,w)*this%hwGF(w);
             enddo
             temp(v) = temp(v)+tp;
          enddo
          !add dx_u
          select case(this%dtmode)
          case(1)
             temp(0) = temp(0)-1.0d0*cc*this%ht*de;
          case(2)
             temp(0) = temp(0)-2.0d0*cc*this%ht*de;
          end select
          ! A^(-1)*()
          do v=0,im
             tp = 0.0d0;
             do w=0,im
                tp=tp + this%hwA(v,w)*temp(w);
             enddo
             this%hwY(v)%Z(i, j, this%Nz) = tp;
          enddo
          this%Ef%Y(i, j, this%Nz) = this%hwY(0)%Z(i, j, this%Nz);
       enddo
    enddo

    do i=1,this%Nx-1
       do v=0,im  
          !j=1
          this%hwY(v)%Z(i, 1, this%Nz) = -(gy1+gy2)*this%hwY(v)%Z(i, 2, this%Nz)-gy1*gy2*this%hwY(v)%Z(i, 3, this%Nz)  &
          +(gy1+gy2)*this%hwYOld(v)%Z(i, 1, this%Nz)+(2.0d0+2.0d0*gy1*gy2)*this%hwYOld(v)%Z(i, 2, this%Nz)+(gy1+gy2)*this%hwYOld(v)%Z(i, 3, this%Nz)  &
          -gy1*gy2*this%hwYOld2(v)%Z(i, 1, this%Nz)-(gy1+gy2)*this%hwYOld2(v)%Z(i, 2, this%Nz)-this%hwYOld2(v)%Z(i, 3, this%Nz);
         !j=Ny
          this%hwY(v)%Z(i, this%Ny, this%Nz) = -(gy1+gy2)*this%hwY(v)%Z(i, this%Ny-1, this%Nz)-gy1*gy2*this%hwY(v)%Z(i, this%Ny-2, this%Nz)  &
          +(gy1+gy2)*this%hwYOld(v)%Z(i, this%Ny, this%Nz)+(2.0d0+2.0d0*gy1*gy2)*this%hwYOld(v)%Z(i, this%Ny-1, this%Nz)+(gy1+gy2)*this%hwYOld(v)%Z(i, this%Ny-2, this%Nz)  &
          -gy1*gy2*this%hwYOld2(v)%Z(i, this%Ny, this%Nz)-(gy1+gy2)*this%hwYOld2(v)%Z(i, this%Ny-1, this%Nz)-this%hwYOld2(v)%Z(i, this%Ny-2, this%Nz);
       enddo
       this%Ef%Y(i, 1, this%Nz) = this%hwY(0)%Z(i, 1, this%Nz);
       this%Ef%Y(i, this%Ny, this%Nz) = this%hwY(0)%Z(i, this%Ny, this%Nz);
    enddo
    
   
    do j=2,this%Ny-1
       do v=0,im  
          !i=0
          this%hwY(v)%Z(0, j, this%Nz) = -(gx1+gx2)*this%hwY(v)%Z(1, j, this%Nz)-gx1*gx2*this%hwY(v)%Z(2, j, this%Nz)  &
          +(gx1+gx2)*this%hwYOld(v)%Z(0, j, this%Nz)+(2.0d0+2.0d0*gx1*gx2)*this%hwYOld(v)%Z(1, j, this%Nz)+(gx1+gx2)*this%hwYOld(v)%Z(2, j, this%Nz)  &
          -gx1*gx2*this%hwYOld2(v)%Z(0, j, this%Nz)-(gx1+gx2)*this%hwYOld2(v)%Z(1, j, this%Nz)-this%hwYOld2(v)%Z(2, j, this%Nz);
          !i=Nx
!         this%hwY(v)%Z(this%Nx, j, this%Nz) = -(gx1+gx2)*this%hwY(v)%Z(this%Nx-1, j, this%Nz)-gx1*gx2*this%hwY(v)%Z(this%Nx-2, j, this%Nz)  &
!         +(gx1+gx2)*this%hwYOld(v)%Z(this%Nx, j, this%Nz)+(2.0d0+2.0d0*gx1*gx2)*this%hwYOld(v)%Z(this%Nx-1, j, this%Nz)+(gx1+gx2)*this%hwYOld(v)%Z(this%Nx-2, j, this%Nz)  &
!         -gx1*gx2*this%hwYOld2(v)%Z(this%Nx, j, this%Nz)-(gx1+gx2)*this%hwYOld2(v)%Z(this%Nx-1, j, this%Nz)-this%hwYOld2(v)%Z(this%Nx-2, j, this%Nz);
       enddo
       this%Ef%Y(0, j, this%Nz) = this%hwY(0)%Z(0, j, this%Nz);
!      this%Ef%Y(this%Nx, j, this%Nz) = this%hwY(0)%Z(this%Nx, j, this%Nz);
    enddo   

    do v=0,im  ! Vectors formation
       !i=0 j=1
          tp = -(gy1+gy2)*this%hwY(v)%Z(0, 2, this%Nz)-gy1*gy2*this%hwY(v)%Z(0, 3, this%Nz)  &
          +(gy1+gy2)*this%hwYOld(v)%Z(0, 1, this%Nz)+(2.0d0+2.0d0*gy1*gy2)*this%hwYOld(v)%Z(0, 2, this%Nz)+(gy1+gy2)*this%hwYOld(v)%Z(0, 3, this%Nz)  &
          -gy1*gy2*this%hwYOld2(v)%Z(0, 1, this%Nz)-(gy1+gy2)*this%hwYOld2(v)%Z(0, 2, this%Nz)-this%hwYOld2(v)%Z(0, 3, this%Nz); 
          this%hwY(v)%Z(0, 1, this%Nz) = (tp -(gx1+gx2)*this%hwY(v)%Z(1, 1, this%Nz)-gx1*gx2*this%hwY(v)%Z(2, 1, this%Nz)  &
          +(gx1+gx2)*this%hwYOld(v)%Z(0, 1, this%Nz)+(2.0d0+2.0d0*gz1*gz2)*this%hwYOld(v)%Z(1, 1, this%Nz)+(gx1+gx2)*this%hwYOld(v)%Z(2, 1, this%Nz)  &
          -gx1*gx2*this%hwYOld2(v)%Z(0, 1, this%Nz)-(gx1+gx2)*this%hwYOld2(v)%Z(1, 1, this%Nz)-this%hwYOld2(v)%Z(2, 1, this%Nz))/2.0d0;
       !i=Nx j=1
!          tp = -(gy1+gy2)*this%hwY(v)%Z(this%Nx, 2, this%Nz)-gy1*gy2*this%hwY(v)%Z(this%Nx, 3, this%Nz)  &
!          +(gy1+gy2)*this%hwYOld(v)%Z(this%Nx, 1, this%Nz)+(2.0d0+2.0d0*gy1*gy2)*this%hwYOld(v)%Z(this%Nx, 2, this%Nz)+(gy1+gy2)*this%hwYOld(v)%Z(this%Nx, 3, this%Nz)  &
!          -gy1*gy2*this%hwYOld2(v)%Z(this%Nx, 1, this%Nz)-(gy1+gy2)*this%hwYOld2(v)%Z(this%Nx, 2, this%Nz)-this%hwYOld2(v)%Z(this%Nx, 3, this%Nz);        
!          this%hwY(v)%Z(this%Nx, 1, this%Nz) = (tp -(gx1+gx2)*this%hwY(v)%Z(this%Nx-1, 1, this%Nz)-gx1*gx2*this%hwY(v)%Z(this%Nx-2, 1, this%Nz)  &
!          +(gx1+gx2)*this%hwYOld(v)%Z(this%Nx, 1, this%Nz)+(2.0d0+2.0d0*gx1*gx2)*this%hwYOld(v)%Z(this%Nx-1, 1, this%Nz)+(gx1+gx2)*this%hwYOld(v)%Z(this%Nx-2, 1, this%Nz)  &
!          -gx1*gx2*this%hwYOld2(v)%Z(this%Nx, 1, this%Nz)-(gx1+gx2)*this%hwYOld2(v)%Z(this%Nx-1, 1, this%Nz)-this%hwYOld2(v)%Z(this%Nx-2, 1, this%Nz))/2.0d0;
       !i=0 j=Ny
          tp = -(gy1+gy2)*this%hwY(v)%Z(0, this%Ny-1, this%Nz)-gy1*gy2*this%hwY(v)%Z(0, this%Ny-2, this%Nz)  &
          +(gy1+gy2)*this%hwYOld(v)%Z(0, this%Ny, this%Nz)+(2.0d0+2.0d0*gy1*gy2)*this%hwYOld(v)%Z(0, this%Ny-1, this%Nz)+(gy1+gy2)*this%hwYOld(v)%Z(0, this%Ny-2, this%Nz)  &
          -gy1*gy2*this%hwYOld2(v)%Z(0, this%Ny, this%Nz)-(gy1+gy2)*this%hwYOld2(v)%Z(0, this%Ny-1, this%Nz)-this%hwYOld2(v)%Z(0, this%Ny-2, this%Nz);
          this%hwY(v)%Z(0, this%Ny, this%Nz) = (tp -(gx1+gx2)*this%hwY(v)%Z(1, this%Ny, this%Nz)-gx1*gx2*this%hwY(v)%Z(2, this%Ny, this%Nz)  &
          +(gx1+gx2)*this%hwYOld(v)%Z(0, this%Ny, this%Nz)+(2.0d0+2.0d0*gx1*gx2)*this%hwYOld(v)%Z(1, this%Ny, this%Nz)+(gx1+gx2)*this%hwYOld(v)%Z(2, this%Ny, this%Nz)  &
          -gx1*gx2*this%hwYOld2(v)%Z(0, this%Ny, this%Nz)-(gx1+gx2)*this%hwYOld2(v)%Z(1, this%Ny, this%Nz)-this%hwYOld2(v)%Z(2, this%Ny, this%Nz))/2.0d0;
       !i=Nx j=Ny
!          tp = -(gy1+gy2)*this%hwY(v)%Z(this%Nx, this%Ny-1, this%Nz)-gy1*gy2*this%hwY(v)%Z(this%Nx, this%Ny-2, this%Nz)  &
!          +(gy1+gy2)*this%hwYOld(v)%Z(this%Nx, this%Ny, this%Nz)+(2.0d0+2.0d0*gy1*gy2)*this%hwYOld(v)%Z(this%Nx, this%Ny-1, this%Nz)+(gy1+gy2)*this%hwYOld(v)%Z(this%Nx, this%Ny-2, this%Nz)  &
!          -gy1*gy2*this%hwYOld2(v)%Z(this%Nx, this%Ny, this%Nz)-(gy1+gy2)*this%hwYOld2(v)%Z(this%Nx, this%Ny-1, this%Nz)-this%hwYOld2(v)%Z(this%Nx, this%Ny-2, this%Nz);  
!          this%hwY(v)%Z(this%Nx, this%Ny, this%Nz) = (tp - (gx1+gx2)*this%hwY(v)%Z(this%Nx-1, this%Ny, this%Nz)-gx1*gx2*this%hwY(v)%Z(this%Nx-2, this%Ny, this%Nz)  &
!          +(gx1+gx2)*this%hwYOld(v)%Z(this%Nx, this%Ny, this%Nz)+(2.0d0+2.0d0*gx1*gx2)*this%hwYOld(v)%Z(this%Nx-1, this%Ny, this%Nz)+(gx1+gx2)*this%hwYOld(v)%Z(this%Nx-2, this%Ny, this%Nz)  &
!          -gx1*gx2*this%hwYOld2(v)%Z(this%Nx, this%Ny, this%Nz)-(gx1+gx2)*this%hwYOld2(v)%Z(this%Nx-1, this%Ny, this%Nz)-this%hwYOld2(v)%Z(this%Nx-2, this%Ny, this%Nz) )/2.0d0;
       enddo
       
    this%Ef%Y(0, 1, this%Nz) = this%hwY(0)%Z(0, 1, this%Nz);
!    this%Ef%Y(this%Nx, 1, this%Nz) = this%hwY(0)%Z(this%Nx, 1, this%Nz);
    this%Ef%Y(0, this%Ny, this%Nz) = this%hwY(0)%Z(0, this%Ny, this%Nz);
!    this%Ef%Y(this%Nx, this%Ny, this%Nz) = this%hwY(0)%Z(this%Nx, this%Ny, this%Nz);


    ! CROSS EDGE

    do j=2,this%Ny-1
       do v=0,im  
          tp = -(gx1+gx2)*this%hwY(v)%Z(this%Nx-1, j, this%Nz)-gx1*gx2*this%hwY(v)%Z(this%Nx-2, j, this%Nz)  &
          +(gx1+gx2)*this%hwYOld(v)%Z(this%Nx, j, this%Nz)+(2.0d0+2.0d0*gx1*gx2)*this%hwYOld(v)%Z(this%Nx-1, j, this%Nz)+(gx1+gx2)*this%hwYOld(v)%Z(this%Nx-2, j, this%Nz)  &
          -gx1*gx2*this%hwYOld2(v)%Z(this%Nx, j, this%Nz)-(gx1+gx2)*this%hwYOld2(v)%Z(this%Nx-1, j, this%Nz)-this%hwYOld2(v)%Z(this%Nx-2, j, this%Nz);                   

          tp = (tp -(gz1+gz2)*this%hwY(v)%X(this%Nx, j, this%Nz-1)-gz1*gz2*this%hwY(v)%X(this%Nx, j, this%Nz-2)  &
          +(gz1+gz2)*this%hwYOld(v)%X(this%Nx, j, this%Nz)+(2.0d0+2.0d0*gz1*gz2)*this%hwYOld(v)%X(this%Nx, j, this%Nz-1)+(gz1+gz2)*this%hwYOld(v)%X(this%Nx, j, this%Nz-2)  &
          -gz1*gz2*this%hwYOld2(v)%X(this%Nx, j, this%Nz)-(gz1+gz2)*this%hwYOld2(v)%X(this%Nx, j, this%Nz-1)-this%hwYOld2(v)%X(this%Nx, j, this%Nz-2))/2.0d0;
          this%hwY(v)%X(this%Nx, j, this%Nz) = tp;
          this%hwY(v)%Z(this%Nx, j, this%Nz) = tp;
       enddo
       this%Ef%Y(this%Nx, j, this%Nz) = this%hwY(0)%X(this%Nx, j, this%Nz);
    enddo

    ! CROSS CORNERS

    do v=0,im  ! Vectors formation
       ! x=Nx face
         tp = -(gy1+gy2)*this%hwY(v)%X(this%Nx, 2, this%Nz)-gy1*gy2*this%hwY(v)%X(this%Nx, 3, this%Nz)  &
         +(gy1+gy2)*this%hwYOld(v)%X(this%Nx, 1, this%Nz)+(2.0d0+2.0d0*gy1*gy2)*this%hwYOld(v)%X(this%Nx, 2, this%Nz)+(gy1+gy2)*this%hwYOld(v)%X(this%Nx, 3, this%Nz)  &
         -gy1*gy2*this%hwYOld2(v)%X(this%Nx, 1, this%Nz)-(gy1+gy2)*this%hwYOld2(v)%X(this%Nx, 2, this%Nz)-this%hwYOld2(v)%X(this%Nx, 3, this%Nz);
         tp = tp -(gz1+gz2)*this%hwY(v)%X(this%Nx, 1, this%Nz-1)-gz1*gz2*this%hwY(v)%X(this%Nx, 1, this%Nz-2)  &
         +(gz1+gz2)*this%hwYOld(v)%X(this%Nx, 1, this%Nz)+(2.0d0+2.0d0*gz1*gz2)*this%hwYOld(v)%X(this%Nx, 1, this%Nz-1)+(gz1+gz2)*this%hwYOld(v)%X(this%Nx, 1, this%Nz-2)  &
         -gz1*gz2*this%hwYOld2(v)%X(this%Nx, 1, this%Nz)-(gz1+gz2)*this%hwYOld2(v)%X(this%Nx, 1, this%Nz-1)-this%hwYOld2(v)%X(this%Nx, 1, this%Nz-2);
       
       ! z=Nz face
         tp = tp -(gy1+gy2)*this%hwY(v)%Z(this%Nx, 2, this%Nz)-gy1*gy2*this%hwY(v)%Z(this%Nx, 3, this%Nz)  &
         +(gy1+gy2)*this%hwYOld(v)%Z(this%Nx, 1, this%Nz)+(2.0d0+2.0d0*gy1*gy2)*this%hwYOld(v)%Z(this%Nx, 2, this%Nz)+(gy1+gy2)*this%hwYOld(v)%Z(this%Nx, 3, this%Nz)  &
         -gy1*gy2*this%hwYOld2(v)%Z(this%Nx, 1, this%Nz)-(gy1+gy2)*this%hwYOld2(v)%Z(this%Nx, 2, this%Nz)-this%hwYOld2(v)%Z(this%Nx, 3, this%Nz);        
         tp = tp -(gx1+gx2)*this%hwY(v)%Z(this%Nx-1, 1, this%Nz)-gx1*gx2*this%hwY(v)%Z(this%Nx-2, 1, this%Nz)  &
         +(gx1+gx2)*this%hwYOld(v)%Z(this%Nx, 1, this%Nz)+(2.0d0+2.0d0*gx1*gx2)*this%hwYOld(v)%Z(this%Nx-1, 1, this%Nz)+(gx1+gx2)*this%hwYOld(v)%Z(this%Nx-2, 1, this%Nz)  &
         -gx1*gx2*this%hwYOld2(v)%Z(this%Nx, 1, this%Nz)-(gx1+gx2)*this%hwYOld2(v)%Z(this%Nx-1, 1, this%Nz)-this%hwYOld2(v)%Z(this%Nx-2, 1, this%Nz);

         tp = tp / 4.0d0;

         this%hwY(v)%X(this%Nx, 1, this%Nz) = tp;
         this%hwY(v)%Z(this%Nx, 1, this%Nz) = tp;
       
       ! x=Nx face
          tp = -(gy1+gy2)*this%hwY(v)%X(this%Nx, this%Ny-1, this%Nz)-gy1*gy2*this%hwY(v)%X(this%Nx, this%Ny-2, this%Nz)  &
          +(gy1+gy2)*this%hwYOld(v)%X(this%Nx, this%Ny, this%Nz)+(2.0d0+2.0d0*gy1*gy2)*this%hwYOld(v)%X(this%Nx, this%Ny-1, this%Nz)+(gy1+gy2)*this%hwYOld(v)%X(this%Nx, this%Ny-2, this%Nz)  &
          -gy1*gy2*this%hwYOld2(v)%X(this%Nx, this%Ny, this%Nz)-(gy1+gy2)*this%hwYOld2(v)%X(this%Nx, this%Ny-1, this%Nz)-this%hwYOld2(v)%X(this%Nx, this%Ny-2, this%Nz);
          tp = tp - (gz1+gz2)*this%hwY(v)%X(this%Nx, this%Ny, this%Nz-1)-gz1*gz2*this%hwY(v)%X(this%Nx, this%Ny, this%Nz-2)  &
          +(gz1+gz2)*this%hwYOld(v)%X(this%Nx, this%Ny, this%Nz)+(2.0d0+2.0d0*gz1*gz2)*this%hwYOld(v)%X(this%Nx, this%Ny, this%Nz-1)+(gz1+gz2)*this%hwYOld(v)%X(this%Nx, this%Ny, this%Nz-2)  &
          -gz1*gz2*this%hwYOld2(v)%X(this%Nx, this%Ny, this%Nz)-(gz1+gz2)*this%hwYOld2(v)%X(this%Nx, this%Ny, this%Nz-1)-this%hwYOld2(v)%X(this%Nx, this%Ny, this%Nz-2);

       ! z=Nz face
          
          tp = tp -(gy1+gy2)*this%hwY(v)%Z(this%Nx, this%Ny-1, this%Nz)-gy1*gy2*this%hwY(v)%Z(this%Nx, this%Ny-2, this%Nz)  &
          +(gy1+gy2)*this%hwYOld(v)%Z(this%Nx, this%Ny, this%Nz)+(2.0d0+2.0d0*gy1*gy2)*this%hwYOld(v)%Z(this%Nx, this%Ny-1, this%Nz)+(gy1+gy2)*this%hwYOld(v)%Z(this%Nx, this%Ny-2, this%Nz)  &
          -gy1*gy2*this%hwYOld2(v)%Z(this%Nx, this%Ny, this%Nz)-(gy1+gy2)*this%hwYOld2(v)%Z(this%Nx, this%Ny-1, this%Nz)-this%hwYOld2(v)%Z(this%Nx, this%Ny-2, this%Nz);  
          tp = tp - (gx1+gx2)*this%hwY(v)%Z(this%Nx-1, this%Ny, this%Nz)-gx1*gx2*this%hwY(v)%Z(this%Nx-2, this%Ny, this%Nz)  &
          +(gx1+gx2)*this%hwYOld(v)%Z(this%Nx, this%Ny, this%Nz)+(2.0d0+2.0d0*gx1*gx2)*this%hwYOld(v)%Z(this%Nx-1, this%Ny, this%Nz)+(gx1+gx2)*this%hwYOld(v)%Z(this%Nx-2, this%Ny, this%Nz)  &
          -gx1*gx2*this%hwYOld2(v)%Z(this%Nx, this%Ny, this%Nz)-(gx1+gx2)*this%hwYOld2(v)%Z(this%Nx-1, this%Ny, this%Nz)-this%hwYOld2(v)%Z(this%Nx-2, this%Ny, this%Nz) ;

          tp = tp / 4.0d0;
        
          this%hwY(v)%X(this%Nx, this%Ny, this%Nz) = tp;
          this%hwY(v)%Z(this%Nx, this%Ny, this%Nz) = tp;
    enddo

    this%Ef%Y(this%Nx, 1, this%Nz) = this%hwY(0)%X(this%Nx, 1, this%Nz);
    this%Ef%Y(this%Nx, this%Ny, this%Nz) = this%hwY(0)%X(this%Nx, this%Ny, this%Nz);

    
    deallocate(temp);
  end subroutine HWEvaEyHigdon3DUpdate


!fillEBoundaryABCHWEvaPML------------------------------------------------------------------
  subroutine  fillEBoundaryABCHWEvaPML(this,t)
    !Fill E boundary with ABC Hagstrom-Warburton High-order 2th order on Ey_x=Nx, Unsplit PML elsewhere
    class(tProblem) :: this
    integer, intent(in) :: t
    integer :: i ,j ,k, v, it  
    !================ EX =================
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do k=0,this%Nz
       do i=1,this%Nx
          this%Ef%X(i, 0, k) = this%PML%Ef%X(i+this%PML%NumPml, 0+this%PML%NumPml, k+this%PML%NumPml);
          this%Ef%X(i, this%Ny, k) = this%PML%Ef%X(i+this%PML%NumPml, this%Ny+this%PML%NumPml, k+this%PML%NumPml);
       enddo
    enddo
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do j=1,this%Ny-1
       do i=1,this%Nx
          this%Ef%X(i, j, 0) = this%PML%Ef%X(i+this%PML%NumPml, j+this%PML%NumPml, 0+this%PML%NumPml);
          this%Ef%X(i, j, this%Nz) = this%PML%Ef%X(i+this%PML%NumPml, j+this%PML%NumPml, this%Nz+this%PML%NumPml);
       enddo
    enddo
    !$OMP END PARALLEL DO
    !================ EY =================
    ! Fill EY Boundaries
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do j=1,this%Ny
       do i=0,this%Nx-1
          this%Ef%Y(i, j, 0) = this%PML%Ef%Y(i+this%PML%NumPml, j+this%PML%NumPml, 0+this%PML%NumPml);
          this%Ef%Y(i, j, this%Nz) = this%PML%Ef%Y(i+this%PML%NumPml, j+this%PML%NumPml, this%Nz+this%PML%NumPml);
       enddo
    enddo
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do k=0,this%Nz
       do j=1,this%Ny
          this%Ef%Y(0, j, k) = this%PML%Ef%Y(0+this%PML%NumPml, j+this%PML%NumPml, k+this%PML%NumPml);
!          this%Ef%Y(this%Nx, j, k) = this%PML%Ef%Y(this%Nx+this%PML%NumPml, j+this%PML%NumPml, k+this%PML%NumPml);
       enddo
    enddo
    !$OMP END PARALLEL DO
    ! Update Ey_x=Nx with HWEva
    call this%HWEvaEyHigdonUpdate(t);
    !================ EZ =================
    !$OMP PARALLEL DO SCHEDULE(GUIDED)   
    do k=1,this%Nz
       do j=0,this%Ny
          this%Ef%Z(0, j, k) = this%PML%Ef%Z(0+this%PML%NumPml, j+this%PML%NumPml, k+this%PML%NumPml);
          this%Ef%Z(this%Nx, j, k) = this%PML%Ef%Z(this%Nx+this%PML%NumPml, j+this%PML%NumPml, k+this%PML%NumPml);
       enddo
    enddo
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do k=1,this%Nz
       do i=1,this%Nx-1
          this%Ef%Z(i, 0, k) = this%PML%Ef%Z(i+this%PML%NumPml, 0+this%PML%NumPml, k+this%PML%NumPml);
         this%Ef%Z(i, this%Ny, k) = this%PML%Ef%Z(i+this%PML%NumPml, this%Ny+this%PML%NumPml, k+this%PML%NumPml);
       enddo
    enddo
    !$OMP END PARALLEL DO      
  end subroutine fillEBoundaryABCHWEvaPML  
    
end module problemclass
