!-----------------------------------------------------------------------
!  Copyright 2017 Mikhail Osintcev
!  This file is part of the EMtool developed at NCSU
!-----------------------------------------------------------------------
! This is the main module of the program
  
! I insert all the modules here to avoid compilation by module  
  include 'commonvars.f90'
  include 'problemclass.f90'
  include 'output.f90'
  include 'meshclass.f90'
  include 'sourceclass.f90'
  include 'pmlclass.f90'
  include 'newdipole.f90'
  include 'sdipole.f90'
  include 'parallel.f90'
  include 'solutionclass.f90'
  include 'bufferclass.f90'
  include 'timersclass.f90'
  include 'poissonclass.f90'
  include 'auxmeshclass.f90'
  include 'edgeauxmeshclass.f90'

! The program starts here  
PROGRAM main
  use commonvars
  use solutionclass
  use output
  use meshclass
  use bufferclass
  use parallel
    
  implicit none;
  
  call DoSimpleFieldModelling(); ! This function allows to build one solution with one grid
  !call DoCheckConvergence();    ! This function allows to build one solution with 1st, 2nd and 4th grids and calculate convergence in Cauchy meaning

END PROGRAM main

!DoSimpleFieldModelling------------------------------------------------------
subroutine DoSimpleFieldModelling()
! This function allows to build one solution with one grid with estimation of the error 
  use commonvars
  use solutionclass
  use output
  use meshclass
  use bufferclass
  use parallel
  class(TBuffer), pointer :: pbuffer;
  class(TBuffer), allocatable, target :: buffers(:);
  real*8 :: start_timer, end_timer, st_start, st_stop;
  type(TSolutionStarter) :: pst;
  type(TSolution) :: sol1;
  integer :: t, i, k;  
  ! Allocate memory for buffers
    allocate(buffers(0:10));
    pbuffer => buffers(0);
    
  ! Divergence free source parameters, good test: size 2, points 24, ftime 50
    a = 0.01;       ! amplification multiplier of dipole
    v = 1;          ! amplification frequency of the dipole
    p = 0.5;        ! smoothing factor of dipole field
    q = 0.1;        ! smoothing factor of amplification sin

  ! Non divergence free source parameters, good test: size 5.55, points 38, ftime 100
    alphaaa=22.0d0;
    om=1.0d0;

  ! Parallelization
    doparallel = 1;  ! 0 - off, 1 - on

  ! Convergence  
    conv = 1;  ! The number of the grid. Place here 2 or 4 to use the second or the third grid.

  ! Output options
    ! Main parameters
    screen_erroroutputsteps = 10*conv;        ! Sreen output frequency of the errors in time steps.
    file_erroroutputsteps = 10*conv;          ! File output frequency of the errors in time steps.
    ! Optional parameters
    file_pmlerrors = 1;                       ! File output frequency of the errors in pml in time steps.
    file_writelocations = 1;                  ! Write to files position of the maximum of the errors.
    file_auxproblemsoutput = 9;               ! amount of auxproblems for output
    file_auxproblemsoutputstep = 1*conv;      ! step in aux problems output
    file_auxproblemswritelocations = 0;       ! write locations or not in aux problems 
    file_compositedivergence = 0;             ! write divergence to file or not
    file_staticfieldmax = 0;                  ! write static field maximums to files
    file_mainfields = 1;                      ! write main problem field values
    file_mainpmlfields = 0;                   ! write pml of main problem field values

  ! Solution initialization
  ! Here you must set up all the parameters properly!  
    pst%sid = 0;
    pst%soltype = 1;
      ! 0 - single solution with analytical bounds
      ! 1 - single solution with PML bounds
      ! 2 - AuxTest 1 - one auxiliary problem with analytical bounds
      ! 3 - AuxTest 2 - one auxiliary problem with PML bounds
      ! 4 - AuxTest 3 - one aux problem, PML bounds and Effective currents
      ! 5 - lacunaes without Poisson and Unsplit PML
      ! 6 - lacunaes with Poisson and Unsplit PML
      ! 7 - pml without separation
      ! 8 - sommerfeld abc
      ! 9 - higdon abc
      ! 10 - betz-mittra abc
      ! 11 - mur abc
      ! 12 - givoli-neta abc ey_X=0, exact elsewhere
      ! 13 - givoli-neta abc ey_X=0, Unsplit PML elsewhere
      ! 14 - givoli-neta abc ey_X=0, one large auxiliary problem without boundaries
      ! 15 - givoli-neta abc ey_X=0, Unsplit PML elsewhere, lacunaes, no poisson, no effective currents
      ! 16 - one big auxiliary problem without boundaries
      ! 17 - hagstrom-warburton ABC on ey_X=0 and exact solution elsewhere
      ! 18 - hagstrom-warburton ABC on ey_X=0, Unsplit PML elsewhere
      ! 19 - hagstrom-warburton ABC on ey_X=0, lacunaes, no poisson, no effective currents
      ! 20 - HWeva ABC on Ey_X=Nx, Exact solution elsewhere
      ! 21 - HWeva ABC on Ey_X=Nx, Unsplit PML elsewhere
      ! 22 - HWeva ABC on Ey_X=Nx, Unsplit PML elsewhere, lacunaes, no poisson, no effective currents, Source 0 only
      ! 23 - HWeva ABC on Ey_X=Nx, Unsplit PML elsewhere, lacunaes, Poisson
      ! 25 - Sommerfeld ABCe, lacunaes, Poisson
      ! 100 - Unsplit PML, lacunaes, no poisson, no effective currents, the main source is cut Source 0 only
      ! 101 - Large Auxiliary Problems + Lacunas - Lacunaes, no poisson, no effective currents, Source 0 only
      
  ! Higdon ABC Parameters
    pst%hig1 = 0.0d0;
    pst%hig2 = 0.0d0;

  ! Givoli-Neta ABC Parameters  
    pst%gnNum = 5;      ! Number of Givoli-Neta angles

  ! Hagstrom-Warburton ABC Parameters  
    pst%hwNum = 2;      ! Number of Hagstrom-Warburton angles
  ! Hagstrom-Warburton EVA ABC Parameters
    pst%hwEvaP = 4;      ! Number of aux variables in HWEva
    pst%hwEvaE = 0;      ! Number of evanescent koeffs in HWEva
    pst%hwEvaPa = 4;     ! Number of edge aux variables in HWEva
    pst%hwEvaEa = 0;     ! Number of edge evanescent koeffs in HWEva 

  ! ABC orders, these parameters is better not to change:
    pst%dtmode = 2;
    pst%dxmode = 2;
    pst%dxauxmode = 2;
    pst%dtauxmode = 2;

  ! PML Parameters
    pst%PML_thickness = 10*conv;
    pst%PML_param = 1.0d0*10.0*1.0*1.0; ! 0.5  

  ! Domain options
    xsize = 6.5d0;       ! x size of computation domain
    ysize = 6.5d0;       ! y size of computation domain
    zsize = 6.5d0;       ! z size of computation domain
    cntx = 0.5d0;        ! x center point of the solution
    cnty = 0.5d0;        ! y center point of the solution
    cntz = 0.5d0;        ! z center point of the solution
    ftime = 100;         ! total simulation time in 'seconds'    
    
  ! Number of points in the domain
    pst%Nx = 66*conv;    !38 ! 61
    pst%Ny = 66*conv;    
    pst%Nz = 66*conv;

  ! Scheme and Source
    pst%main_schemetype = 0;
      ! 0 - Yee 2-2 scheme
    pst%sourcetype = 0;
      ! 0 - Dipole
      ! 1 - SDipole 

  ! Lacunas and auxiliary problems  
    ! Number of points for Poisson
      pst%npx = 16*conv-1;        
      pst%npy = 16*conv-1;
      pst%npz = 16*conv-1;
    
    ! Difference between main problem and aux problems
      pst%dNxaux = 0*conv;
      pst%dNyaux = 0*conv;
      pst%dNzaux = 0*conv; 

    ! Effective currents options
      pst%Namu = 3*conv;
      pst%Ndmu = 4*conv;
 
    ! Auxiliary tasks options
      pst%aux_schemetype = 0;
      pst%auxmaxtime = 40.0d0;  ! The total simulation time of auxiliary problems
      pst%sgm_s = 0.8d0;        ! Sigma parameter for auxiliary problems 
      pst%auxgaptime = 10*conv; ! The additional time for each auxiliary problem (in time steps)

  ! Initialize timers for time estimation
  do k=0,31
     times(k) = 0;
  enddo
    
  call InitializeOutput;          ! Open files and prepare for output
  call InitParallelParameters;    ! Initiate OpenMP

  call sol1%solution_init(pst, pbuffer);   ! Solution initialization

  pause ! Pause before the start of simulation

  call cpu_time(start_timer)
  !$ start_timer = OMP_get_wtime();
  call PrintStart(sol1);                   

  call sol1%solution_build;                ! THE MAIN ALGORITHM

  call cpu_time(end_timer)
  !$ end_timer = OMP_get_wtime();
  call PrintFinish(end_timer-start_timer, sol1)
  
  call sol1%solution_destroy;    ! Destoy the object and free memory

  call FinalizeOutput;           ! Final output
  deallocate(buffers);
  
end subroutine DoSimpleFieldModelling




!DoCheckConvergence------------------------------------------------------
subroutine DoCheckConvergence()
! This function allows to build one solution with 1st, 2nd and 4th grids and calculate convergence in Cauchy meaning
  use commonvars
  use solutionclass
  use output
  use meshclass
  use bufferclass
  use parallel
  ! Check convergence on the grids
  class(TBuffer), pointer :: pbuffer, pbuffer2, pbuffer4;
  class(TBuffer), allocatable, target :: buffers(:);
  real*8 :: start_timer, end_timer, st_start, st_stop;
  type(TSolutionStarter) :: pst;
  type(TSolution) :: sol1, sol3, sol5;
  integer :: t, t2, t4, i;
  integer :: params(0:11);
  ! Buffer
    allocate(buffers(0:10));
    pbuffer => buffers(0);
    pbuffer2 => buffers(1);
    pbuffer4 => buffers(2);
    
  ! Divergence free source parameters, good test: size 2, points 24, ftime 50
    a = 0.01;       ! amplification multiplier of dipole
    v = 1;          ! amplification frequency of the dipole
    p = 0.5;        ! smoothing factor of dipole field
    q = 0.1;        ! smoothing factor of amplification sin

  ! Non divergence free source parameters, good test: size 5.55, points 38, ftime 100
    alphaaa=22.0d0;
    om=1.0d0;

  ! Parallelization
    doparallel = 1;  ! 0 - off, 1 - on

  ! Convergence  
    conv = 1;  ! The number of the grid. Place here 2 or 4 to use the second or the third grid.

  ! Output options
    ! Main parameters
    screen_erroroutputsteps = 10*conv;        ! Sreen output frequency of the errors in time steps.
    file_erroroutputsteps = 10*conv;          ! File output frequency of the errors in time steps.
    ! Optional parameters
    file_pmlerrors = 1;                       ! File output frequency of the errors in pml in time steps.
    file_writelocations = 1;                  ! Write to files position of the maximum of the errors.
    file_auxproblemsoutput = 9;               ! amount of auxproblems for output
    file_auxproblemsoutputstep = 1*conv;      ! step in aux problems output
    file_auxproblemswritelocations = 0;       ! write locations or not in aux problems 
    file_compositedivergence = 0;             ! write divergence to file or not
    file_staticfieldmax = 0;                  ! write static field maximums to files
    file_mainfields = 1;                      ! write main problem field values
    file_mainpmlfields = 0;                   ! write pml of main problem field values

  ! Solution initialization
  ! Here you must set up all the parameters properly!  
    pst%sid = 0;
    pst%soltype = 1;
      ! 0 - single solution with analytical bounds
      ! 1 - single solution with PML bounds
      ! 2 - AuxTest 1 - one auxiliary problem with analytical bounds
      ! 3 - AuxTest 2 - one auxiliary problem with PML bounds
      ! 4 - AuxTest 3 - one aux problem, PML bounds and Effective currents
      ! 5 - lacunaes without Poisson and Unsplit PML
      ! 6 - lacunaes with Poisson and Unsplit PML
      ! 7 - pml without separation
      ! 8 - sommerfeld abc
      ! 9 - higdon abc
      ! 10 - betz-mittra abc
      ! 11 - mur abc
      ! 12 - givoli-neta abc ey_X=0, exact elsewhere
      ! 13 - givoli-neta abc ey_X=0, Unsplit PML elsewhere
      ! 14 - givoli-neta abc ey_X=0, one large auxiliary problem without boundaries
      ! 15 - givoli-neta abc ey_X=0, Unsplit PML elsewhere, lacunaes, no poisson, no effective currents
      ! 16 - one big auxiliary problem without boundaries
      ! 17 - hagstrom-warburton ABC on ey_X=0 and exact solution elsewhere
      ! 18 - hagstrom-warburton ABC on ey_X=0, Unsplit PML elsewhere
      ! 19 - hagstrom-warburton ABC on ey_X=0, lacunaes, no poisson, no effective currents
      ! 20 - HWeva ABC on Ey_X=Nx, Exact solution elsewhere
      ! 21 - HWeva ABC on Ey_X=Nx, Unsplit PML elsewhere
      ! 22 - HWeva ABC on Ey_X=Nx, Unsplit PML elsewhere, lacunaes, no poisson, no effective currents, Source 0 only
      ! 23 - HWeva ABC on Ey_X=Nx, Unsplit PML elsewhere, lacunaes, Poisson
      ! 25 - Sommerfeld ABCe, lacunaes, Poisson
      ! 100 - Unsplit PML, lacunaes, no poisson, no effective currents, the main source is cut Source 0 only
      ! 101 - Large Auxiliary Problems + Lacunas - Lacunaes, no poisson, no effective currents, Source 0 only
      
  ! Higdon ABC Parameters
    pst%hig1 = 0.0d0;
    pst%hig2 = 0.0d0;

  ! Givoli-Neta ABC Parameters  
    pst%gnNum = 5;      ! Number of Givoli-Neta angles

  ! Hagstrom-Warburton ABC Parameters  
    pst%hwNum = 2;      ! Number of Hagstrom-Warburton angles
  ! Hagstrom-Warburton EVA ABC Parameters
    pst%hwEvaP = 4;      ! Number of aux variables in HWEva
    pst%hwEvaE = 0;      ! Number of evanescent koeffs in HWEva
    pst%hwEvaPa = 4;     ! Number of edge aux variables in HWEva
    pst%hwEvaEa = 0;     ! Number of edge evanescent koeffs in HWEva 

  ! ABC orders, these parameters is better not to change:
    pst%dtmode = 2;
    pst%dxmode = 2;
    pst%dxauxmode = 2;
    pst%dtauxmode = 2;

  ! PML Parameters
    pst%PML_thickness = 10*conv;
    pst%PML_param = 1.0d0*10.0*1.0*1.0; ! 0.5  

  ! Domain options
    xsize = 6.5d0;      ! x size of computation domain
    ysize = 6.5d0;      ! y size of computation domain
    zsize = 6.5d0;      ! z size of computation domain
    cntx = 0.5d0;       ! x center point of the solution
    cnty = 0.5d0;       ! y center point of the solution
    cntz = 0.5d0;       ! z center point of the solution
    ftime = 20;         ! total simulation time in 'seconds'    
    
  ! Number of points in the domain
    pst%Nx = 66*conv;    !38 ! 61
    pst%Ny = 66*conv;    
    pst%Nz = 66*conv;

  ! Scheme and Source
    pst%main_schemetype = 0;
      ! 0 - Yee 2-2 scheme
    pst%sourcetype = 1;
      ! 0 - Dipole
      ! 1 - SDipole 

  ! Lacunas and auxiliary problems  
    ! Number of points for Poisson
      pst%npx = 16*conv-1;        
      pst%npy = 16*conv-1;
      pst%npz = 16*conv-1;
    
    ! Difference between main problem and aux problems
      pst%dNxaux = 0*conv;
      pst%dNyaux = 0*conv;
      pst%dNzaux = 0*conv; 

    ! Effective currents options
      pst%Namu = 3*conv;
      pst%Ndmu = 4*conv;
 
    ! Auxiliary tasks options
      pst%aux_schemetype = 0;
      pst%auxmaxtime = 40.0d0;  ! The total simulation time of auxiliary problems
      pst%sgm_s = 0.8d0;        ! Sigma parameter for auxiliary problems 
      pst%auxgaptime = 10*conv; ! The additional time for each auxiliary problem (in time steps)
  
  call InitializeOutput;
  call InitParallelParameters;

  params(0) = pst%Nx;
  params(1) = pst%Ny;
  params(2) = pst%Nz;
  params(3) = pst%npx;
  params(4) = pst%npy;
  params(5) = pst%npz;
  params(6) = pst%dNxaux;
  params(7) = pst%dNyaux;
  params(8) = pst%dNzaux;
  params(9) = pst%PML_thickness;
  params(10) = pst%Namu;
  params(11) = pst%Ndmu;

  call sol1%solution_init(pst, pbuffer);

  pst%sid = 3;
  pst%Nx = (params(0)-1)*3+1;
  pst%Ny = (params(1)-1)*3+1;
  pst%Nz = (params(2)-1)*3+1;
  pst%npx = (params(3)-1)*3+1;
  pst%npy = (params(4)-1)*3+1;
  pst%npz = (params(5)-1)*3+1;
  pst%dNxaux = params(6)*3;
  pst%dNyaux = params(7)*3;
  pst%dNzaux = params(8)*3;
  pst%PML_thickness = params(9)*3;
  pst%Namu = params(10)*3;
  pst%Ndmu = params(11)*3;
  
  call sol3%solution_init(pst, pbuffer2);

  pst%sid = 5;
  pst%Nx = (params(0)-1)*5+1;
  pst%Ny = (params(1)-1)*5+1;
  pst%Nz = (params(2)-1)*5+1;
  pst%npx = (params(3)-1)*5+1;
  pst%npy = (params(4)-1)*5+1;
  pst%npz = (params(5)-1)*5+1;
  pst%dNxaux = params(6)*5;
  pst%dNyaux = params(7)*5;
  pst%dNzaux = params(8)*5;
  pst%PML_thickness = params(9)*5;
  pst%Namu = params(10)*5;
  pst%Ndmu = params(11)*5;
  
  call sol5%solution_init(pst, pbuffer4);
  pause

  write(*,*) 'SOL 1  ====Xi===='
  do k=0,sol1%mainproblem%Nx
     write(*,*) sol1%mainproblem%Xi(k);   
  enddo
  write(*,*) 'SOL 3  ====Xi===='
  do k=0,sol3%mainproblem%Nx
     if (MOD(k,3)==0) then
        write(*,*) sol3%mainproblem%Xi(k);
     endif
  enddo
  write(*,*) 'SOL 5  ====Xi===='
  do k=0,sol5%mainproblem%Nx
     if (MOD(k,5)==0) then
        write(*,*) sol5%mainproblem%Xi(k);
     endif
  enddo

  write(*,*) 'SOL 1  ====Xi05===='
  do k=1,sol1%mainproblem%Nx
     write(*,*) sol1%mainproblem%Xi05(k);   
  enddo
  write(*,*) 'SOL 3  ====Xi05===='
  do k=2,sol3%mainproblem%Nx
     if (MOD(k-2,3)==0) then
        write(*,*) sol3%mainproblem%Xi05(k);
     endif
  enddo
  write(*,*) 'SOL 5  ====Xi05===='
  do k=3,sol5%mainproblem%Nx
     if (MOD(k-3,5)==0) then
        write(*,*) sol5%mainproblem%Xi05(k);
     endif
  enddo
  
  call cpu_time(start_timer)
  !$ start_timer = OMP_get_wtime();
  call PrintStart(sol1);

  do t=1,sol1%Nt
   
     call sol1%solution_DoStep(t);
     write(*,*) 'Sol 1 step';
     do i=1,3
        t2=(t-1)*3+i;
        call sol3%solution_DoStep(t2);
        write(*,*) 'Sol 3 step';
     enddo
     do i=1,5
        t4=(t-1)*5+i;
        call sol5%solution_DoStep(t4);
        write(*,*) 'Sol 5 step';
     enddo
     call CheckConvCauchy(sol1,sol3,sol5,t)
  enddo
 
  call cpu_time(end_timer)
  !$ end_timer = OMP_get_wtime();
  call PrintFinish(end_timer-start_timer, sol1)
  
  call sol1%solution_destroy;
  call sol3%solution_destroy;
  call sol5%solution_destroy;

  call FinalizeOutput;
  deallocate(buffers);
end subroutine DoCheckConvergence


!CheckConvCauchy------------------------------------------
  subroutine CheckConvCauchy(sol,sol3,sol5,t)
  ! Calculate the maximum Error in the current time-point based on three grids
    use commonvars
    use solutionclass
    use output
    use meshclass
    use bufferclass
    use parallel
    type(tSolution) :: sol, sol3, sol5
    real*8 :: a1, a2, a3, cur, maxEy, d1, d2;
    integer i, j, k, i3, j3, k3, i5, j5, k5,Is,Js,Ks;
    integer, intent(in) :: t;
    real*8, dimension(6,2) :: localdiffs, diffs;
    real*8, dimension(6) :: speeds;
    is=0;
    js=0;
    ks=0;
    do k=0,sol%Nz
       do j=1,sol%Ny
          do i=0,sol%Nx
             i3 = i*3;
             j3 = j*3-1;
             k3 = k*3;
             i5 = i*5;
             j5 = j*5-2;
             k5 = k*5;             
             a1 = sol%mainproblem%Ef%Y(i,j,k);
             if (abs(a1)>d1) then
                d1 = abs(a1);
                Is=i;
                Js=j;
                Ks=k;
             endif
          enddo
       enddo
    enddo
    if (sol3%mainproblem%Ef%Y(is*3,js*3-1,ks*3)-sol5%mainproblem%Ef%Y(is*5,js*5-2,ks*5)/=0) then
       write(*,*) (sol%mainproblem%Ef%Y(is,js,ks)-sol3%mainproblem%Ef%Y(is*3,js*3-1,ks*3))/(sol3%mainproblem%Ef%Y(is*3,js*3-1,ks*3)-sol5%mainproblem%Ef%Y(is*5,js*5-2,ks*5));
    endif    
    640 format('Ey Couchy Minimum Convergence:', E14.7)
    641 format('d1:', E14.7)
    642 format('d2:', E14.7)
end subroutine CheckConvCauchy
      
      

