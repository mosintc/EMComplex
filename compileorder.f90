chmod +x ./compileall

ifort -c -qopenmp commonvars.f90
ifort -c -qopenmp dipole.f90
ifort -c -qopenmp sdipole.f90
ifort -c -qopenmp newdipole.f90
ifort -c -qopenmp timersclass.f90
ifort -c -qopenmp meshclass.f90
ifort -c -qopenmp auxmeshclass.f90
ifort -c -qopenmp edgeauxmeshclass.f90
ifort -c -qopenmp sourceclass.f90
ifort -c -qopenmp bufferclass.f90
ifort -c -qopenmp parallel.f90
ifort -c -qopenmp poissonclass.f90
ifort -c -qopenmp sourceclass.f90
ifort -c -qopenmp pmlclass.f90
ifort -c -qopenmp problemclass.f90
ifort -c -qopenmp solutionclass.f90
ifort -o main main.f90 fftsg.f fftsg3d.f -qopenmp -mkl
ifort -o main main.f90 fftsg.f fftsg3d.f -qopenmp -mkl


ifort -c -qopenmp commonvars.f90 -mmic
ifort -c -qopenmp dipole.f90 -mmic
ifort -c -qopenmp sdipole.f90 -mmic
ifort -c -qopenmp timersclass.f90 -mmic
ifort -c -qopenmp meshclass.f90 -mmic
ifort -c -qopenmp auxmeshclass.f90 -mmic
ifort -c -qopenmp edgeauxmeshclass.f90 -mmic
ifort -c -qopenmp sourceclass.f90 -mmic
ifort -c -qopenmp bufferclass.f90 -mmic
ifort -c -qopenmp parallel.f90 -mmic
ifort -c -qopenmp poissonclass.f90 -mmic
ifort -c -qopenmp pmlclass.f90 -mmic
ifort -c -qopenmp problemclass.f90 -mmic
ifort -c -qopenmp solutionclass.f90 -mmic
ifort -o main main.f90 fftsg.f fftsg3d.f -qopenmp -mkl -mmic
ifort -o main main.f90 fftsg.f fftsg3d.f -qopenmp -mkl -mmic

ifort -c commonvars.f90
ifort -c dipole.f90
ifort -c newdipole.f90
ifort -c sdipole.f90
ifort -c timersclass.f90
ifort -c meshclass.f90
ifort -c auxmeshclass.f90
ifort -c edgeauxmeshclass.f90
ifort -c sourceclass.f90
ifort -c bufferclass.f90
ifort -c parallel.f90
ifort -c poissonclass.f90
ifort -c pmlclass.f90
ifort -c problemclass.f90
ifort -c solutionclass.f90
ifort -o main main.f90 fftsg.f fftsg3d.f -mkl
ifort -o main main.f90 fftsg.f fftsg3d.f -mkl
