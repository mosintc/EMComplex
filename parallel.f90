module Parallel
  use commonvars;
  use omp_lib;
contains
  
!AddParBalance---------------------------------------------------------
  subroutine AddParBalance()
  !Add information about balance in calculations
    integer :: proc;
    !$   proc = OMP_get_thread_num()+1;
    !$   balance(proc)=balance(proc)+1;
  end subroutine AddParBalance

!InitParallelParameters-------------------------------------------------
  subroutine InitParallelParameters();
 
    !$  nprocs = OMP_get_max_threads();
    !$  if ((doparallel == 1).AND.(nprocs>0)) then
    !$     allocate(balance(1:nprocs));
    !$     do i=1,nprocs
    !$        balance(i)=0;
    !$     enddo
    !$  else
    !$    call OMP_set_dynamic(.FALSE.)
    !$    call OMP_set_num_threads(1)
    !$    nprocs = 1;
    !$    allocate(balance(1:1));
    !$    balance(1)=0;
    !$  endif
    
  end subroutine InitParallelParameters;
  
end module Parallel
