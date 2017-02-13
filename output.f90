module output
! Module for output results to the screen of file  
  use problemclass
  use solutionclass
  implicit none;
  integer :: cfile;
  integer :: files(0:499);
 contains

   subroutine openfile(unt, fl)
     integer :: unt;
     character (len=*) :: fl;
     open(unit=unt, file=fl, status='replace', action='write');
     cfile = cfile+1;
     files(cfile) = unt;
   end subroutine openfile
   
  
!InitializeOutput
   subroutine InitializeOutput
     character (len=100) :: fname;
     integer :: i, k;
     do i=0,499
        files(i) = -1;
     enddo
     cfile = -1;
     
     ! Main problem output 100
     call openfile(100, 'Ex_error.txt');
     call openfile(101, 'Ey_error.txt');
     call openfile(102, 'Ez_error.txt');
     call openfile(103, 'Hx_error.txt');
     call openfile(104, 'Hy_error.txt');
     call openfile(105, 'Hz_error.txt');
     call openfile(120, 'Timing.txt');
     if (file_mainfields==1) then
        call openfile(106, './output/Ex_main.txt');
        call openfile(107, './output/Ey_main.txt');
        call openfile(108, './output/Ez_main.txt');
        call openfile(109, './output/Hx_main.txt');
        call openfile(110, './output/Hy_main.txt');
        call openfile(111, './output/Hz_main.txt');
        call openfile(112, './output/JEx_main.txt');
        call openfile(113, './output/JEy_main.txt');
        call openfile(114, './output/JEz_main.txt');
        call openfile(115, './output/JHx_main.txt');
        call openfile(116, './output/JHy_main.txt');
        call openfile(117, './output/JHz_main.txt');
     endif

     ! PML output 200
     if (file_pmlerrors==1) then
        call openfile(200, './output/PML_Ex_error.txt');
        call openfile(201, './output/PML_Ey_error.txt');
        call openfile(202, './output/PML_Ez_error.txt');
        call openfile(203, './output/PML_Hx_error.txt');
        call openfile(204, './output/PML_Hy_error.txt');
        call openfile(205, './output/PML_Hz_error.txt');   
     endif
     if (file_mainpmlfields==1) then
        call openfile(206, './output/Ex_mainPML.txt');
        call openfile(207, './output/Ey_mainPML.txt');
        call openfile(208, './output/Ez_mainPML.txt');
        call openfile(209, './output/Hx_mainPML.txt');
        call openfile(210, './output/Hy_mainPML.txt');
        call openfile(211, './output/Hz_mainPML.txt');
        call openfile(212, './output/JEx_mainPML.txt');
        call openfile(213, './output/JEy_mainPML.txt');
        call openfile(214, './output/JEz_mainPML.txt');
        call openfile(215, './output/JHx_mainPML.txt');
        call openfile(216, './output/JHy_mainPML.txt');
        call openfile(217, './output/JHz_mainPML.txt');
     endif

     !Auxiliary problems output - 1000
     if (file_auxproblemsoutput>0) then
        do i=0,file_auxproblemsoutput-1
           write(fname, "(A15,I1,A4)") './output/Ex_aux', i, '.txt';
           call openfile(1000+i*100, fname);
           fname = '';
           write(fname, "(A15,I1,A4)") './output/Ey_aux', i, '.txt';
           call openfile(1000+i*100+1, fname);
           fname = '';
           write(fname, "(A15,I1,A4)") './output/Ez_aux', i, '.txt';
           call openfile(1000+i*100+2, fname);
           fname = '';
           write(fname, "(A15,I1,A4)") './output/Hx_aux', i, '.txt';
           call openfile(1000+i*100+3, fname);
           fname = '';
           write(fname, "(A15,I1,A4)") './output/Hy_aux', i, '.txt';
           call openfile(1000+i*100+4, fname);
           fname = '';
           write(fname, "(A15,I1,A4)") './output/Hz_aux', i, '.txt';
           call openfile(1000+i*100+5, fname);
           fname = '';
           write(fname, "(A16,I1,A4)") './output/JEx_aux', i, '.txt';
           call openfile(1000+i*100+6, fname);
           fname = '';
           write(fname, "(A16,I1,A4)") './output/JEy_aux', i, '.txt';
           call openfile(1000+i*100+7, fname);
           fname = '';
           write(fname, "(A16,I1,A4)") './output/JEz_aux', i, '.txt';
           call openfile(1000+i*100+8, fname);
           fname = '';
           write(fname, "(A16,I1,A4)") './output/JHx_aux', i, '.txt';
           call openfile(1000+i*100+9, fname);
           fname = '';
           write(fname, "(A16,I1,A4)") './output/JHy_aux', i, '.txt';
           call openfile(1000+i*100+10, fname);
           fname = '';
           write(fname, "(A16,I1,A4)") './output/JHz_aux', i, '.txt';
           call openfile(1000+i*100+11, fname);
           fname = '';
           write(fname, "(A14,I1,A4)") './output/Theta', i, '.txt';
           call openfile(1000+i*100+12, fname);
           fname = '';
        enddo
           do i=0,file_auxproblemsoutput-1
              write(fname, "(A19,I1,A4)") './output/Ex_pml_aux', i, '.txt';
              call openfile(2000+i*100, fname);
              fname = '';
              write(fname, "(A19,I1,A4)") './output/Ey_pml_aux', i, '.txt';
              call openfile(2000+i*100+1, fname);
              fname = '';
              write(fname, "(A19,I1,A4)") './output/Ez_pml_aux', i, '.txt';
              call openfile(2000+i*100+2, fname);
              fname = '';
              write(fname, "(A19,I1,A4)") './output/Hx_pml_aux', i, '.txt';
              call openfile(2000+i*100+3, fname);
              fname = '';
              write(fname, "(A19,I1,A4)") './output/Hy_pml_aux', i, '.txt';
              call openfile(2000+i*100+4, fname);
              fname = '';
              write(fname, "(A19,I1,A4)") './output/Hz_pml_aux', i, '.txt';
              call openfile(2000+i*100+5, fname);
              fname = '';

              write(fname, "(A19,I1,A4)") './output/JEx_pml_aux', i, '.txt';
              call openfile(3000+i*100, fname);
              fname = '';
              write(fname, "(A19,I1,A4)") './output/JEy_pml_aux', i, '.txt';
              call openfile(3000+i*100+1, fname);
              fname = '';
              write(fname, "(A19,I1,A4)") './output/JEz_pml_aux', i, '.txt';
              call openfile(3000+i*100+2, fname);
              fname = '';
              write(fname, "(A19,I1,A4)") './output/JHx_pml_aux', i, '.txt';
              call openfile(3000+i*100+3, fname);
              fname = '';
              write(fname, "(A19,I1,A4)") './output/JHy_pml_aux', i, '.txt';
              call openfile(3000+i*100+4, fname);
              fname = '';
              write(fname, "(A19,I1,A4)") './output/JHz_pml_aux', i, '.txt';
              call openfile(3000+i*100+5, fname);
              fname = '';
           enddo
     endif

     ! Divergence 300
     if (file_compositedivergence==1) then
        call openfile(300, './output/DivEComposite.txt');
     endif   
     
     ! Static Field 400
     if (file_staticfieldmax==1) then
        call openfile(400, './output/static_Ex_max.txt');
        call openfile(401, './output/static_Ey_max.txt');
        call openfile(402, './output/static_Ez_max.txt');
        call openfile(403, './output/static_Hx_max.txt');
        call openfile(404, './output/static_Hy_max.txt');
        call openfile(405, './output/static_Hz_max.txt');
     endif
  
     ! Free output 500
     call openfile(500, './output/temp.txt');

     ! ABCs 600
     call openfile(600, './output/hw.txt');
     call openfile(601, './output/hw_edge.txt');
     do k=0,31
       timestitles(k) = 'NA';
     enddo
  end subroutine InitializeOutput;

!FinalizeOutput
  subroutine FinalizeOutput
    ! Close all the files
    logical :: isopened;
    integer :: i;
    if (cfile>-1) then
       do i=0,cfile-1
          inquire(files(i),opened=isopened)
          if (isopened==.TRUE.) then
             close(unit=files(i))
          endif   
       enddo
    endif
  end subroutine FinalizeOutput
 

!PrintStart-----------------------------------------------
  subroutine PrintStart(sol)
    !Print the beginning of simulation
    type(tSolution) :: sol;
    integer :: k;
    write(*,*) '--------------COMPUTATION PARAMETERS--------------'
    call PrintParameters(sol%mainproblem)
    write(*,*) '------ Simulation started'
  end subroutine PrintStart

  
!PrintParameters------------------------------------------    
  subroutine PrintParameters(prob)
    !Print modelling parameters
    type(tProblem) :: prob;
    integer :: i,j, pn;
    select case (doparallel)
    case (0)
      write(*,*) "No parallelization"
    case (1)
      if (nprocs > 0) then
         write(*,31) nprocs;
         31 format ('Prallelization on threads: ',I2)
      else
         write(*,*) 'Parallelization in unavaible. Please use -openmp directive while compiling.';
      endif
    end select
    write(*,22) xsize, ysize, ysize;
    22 format ('Size of Domain: (',F7.3,',',F7.3,',',F7.3,' )');
    write(*,23) ftime;
    23 format ('Simulation Time: ',F7.3);
  end subroutine PrintParameters 
 
  
!PrintFinish------------------------------------------
  subroutine PrintFinish(simtime, sol)
    !Print information at the end of simulation
    integer :: i, k;
    real*8, intent(in) :: simtime;
    type(tSolution) :: sol;
    write(*,*) '------ Simulation finished'
    call PrintParameters(sol%mainproblem);
    write(*,*) 'Maximum errors:'
    write(*,89) sol%mainproblem%Hx_maxerr%val;
    write(*,90) sol%mainproblem%Hy_maxerr%val;
    write(*,91) sol%mainproblem%Hz_maxerr%val;
    write(*,86) sol%mainproblem%Ex_maxerr%val;
    write(*,87) sol%mainproblem%Ey_maxerr%val;
    write(*,88) sol%mainproblem%Ez_maxerr%val;
    86  format ('errormaxEX:',1x,E22.15)
    87  format ('errormaxEY:',1x,E22.15)
    88  format ('errormaxEZ:',1x,E22.15)
    89  format ('errormaxHX:',1x,E22.15)
    90  format ('errormaxHY:',1x,E22.15)
    91  format ('errormaxHZ:',1x,E22.15)
    write(*,92) simtime
    92  format ('Time elapsed:',1x,E14.7)

    write(unit=120, fmt=94) simtime;
    do k=0,31
       write(*,93) k, timestitles(k), times(k);
       write(unit=120, fmt=94) times(k);
    enddo
    write(unit=120, fmt=95) 'TOTAL TIME';
    do k=0,31
       write(unit=120, fmt=95) timestitles(k);
    enddo
    93  format ('Timer',I2,' - ',A10,':',1x,E14.7)
    94  format (E14.7)
    95  format (A10)   
  end subroutine PrintFinish
  
end module output
