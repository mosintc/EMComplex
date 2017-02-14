!-----------------------------------------------------------------------
!  Copyright 2017 Mikhail Osintcev
!  This file is part of the EMtool developed at NCSU
!-----------------------------------------------------------------------
! This module contains implementation of non divergence free source for T. Hagstrom

module testsolution
  ! Computate the analytic solution of the dipole
  implicit none
  real*8, parameter :: Pi = 3.141592653589793238462643d0;
contains

! ==== ELECTRIC COMPONENTS ====
real*8 function Ex_exact(x,y,z,t,alpha,om,c)
implicit none
  real*8, intent(in) :: x,y,z,t,alpha,om,c
  real*8 ro1,ro2,cosine,sine,wp,p_p,p_pp,p_ppp,p_pppp,psi,E_r,E_theta,Eperp
  ro1=sqrt(x**2+y**2+z**2)
  ro2=sqrt(x**2+y**2)
  cosine=z/ro1;sine=ro2/ro1
  if (t-ro1>=0.0d0) then
     call steady_oscillate(t-ro1,om,wp,p_p,p_pp,p_ppp,p_pppp,0)
     psi=1.0d0-exp(-alpha*ro1**5)
     E_r=2.0d0*(psi/ro1**3)*cosine*(wp+p_p*ro1/c)-5.0d0*alpha*ro1**2*exp(-alpha*ro1**5)*cosine*(wp+p_p*ro1/c+p_pp*ro1**2/(c**2))
     E_theta=(psi/ro1**3)*sine*(wp+p_p*ro1/c+p_pp*ro1**2/(c**2))
     Eperp=E_r*sine+E_theta*cosine
     Ex_exact=Eperp*x/ro2 
  else
     Ex_exact=0.0d0
  end if	
end function Ex_exact


real*8 function Ey_exact(x,y,z,t,alpha,om,c)
  implicit none
  real*8, intent(in) :: x,y,z,t,alpha,om,c
  real*8 ro1,ro2,cosine,sine,wp,p_p,p_pp,p_ppp,p_pppp,psi,E_r,E_theta,Eperp
  ro1=sqrt(x**2+y**2+z**2)
  ro2=sqrt(x**2+y**2)
  cosine=z/ro1;sine=ro2/ro1
  if (t-ro1>=0.0d0) then
     call steady_oscillate(t-ro1,om,wp,p_p,p_pp,p_ppp,p_pppp,0)
     psi=1.0d0-exp(-alpha*ro1**5)
     E_r=2.0d0*(psi/ro1**3)*cosine*(wp+p_p*ro1/c)-5.0d0*alpha*ro1**2*exp(-alpha*ro1**5)*cosine*(wp+p_p*ro1/c+p_pp*ro1**2/(c**2))
     E_theta=(psi/ro1**3)*sine*(wp+p_p*ro1/c+p_pp*ro1**2/(c**2))
     Eperp=E_r*sine+E_theta*cosine
     Ey_exact=Eperp*y/ro2 
  else
     Ey_exact=0.0d0
  end if	
end function Ey_exact


real*8 function Ez_exact(x,y,z,t,alpha,om,c)
  implicit none
  real*8, intent(in) :: x,y,z,t,alpha,om,c
  real*8 ro1,ro2,cosine,sine,wp,p_p,p_pp,p_ppp,p_pppp,psi,E_r,E_theta,Eperp
  ro1=sqrt(x**2+y**2+z**2)
  ro2=sqrt(x**2+y**2)
  cosine=z/ro1;sine=ro2/ro1
  if (t-ro1>=0.0d0) then
     call steady_oscillate(t-ro1,om,wp,p_p,p_pp,p_ppp,p_pppp,0)
     psi=1.0d0-exp(-alpha*ro1**5)
     E_r=2.0d0*(psi/ro1**3)*cosine*(wp+p_p*ro1/c)-5.0d0*alpha*ro1**2*exp(-alpha*ro1**5)*cosine*(wp+p_p*ro1/c+p_pp*ro1**2/(c**2))
     E_theta=(psi/ro1**3)*sine*(wp+p_p*ro1/c+p_pp*ro1**2/(c**2))
     Ez_exact=E_r*cosine-E_theta*sine
  else
     Ez_exact=0.0d0
  end if	
end function Ez_exact

! ==== MAGNETIC COMPONENTS ====

real*8 function Hx_exact(x,y,z,t,alpha,om,c)
  implicit none
  real*8, intent(in) :: x,y,z,t,alpha,om,c
  real*8 ro1,ro2,cosine,sine,wp,p_p,p_pp,p_ppp,p_pppp,psi,E_r,E_theta,Eperp,sinp
  ro1=sqrt(x**2+y**2+z**2)
  ro2=sqrt(x**2+y**2)
  cosine=z/ro1;sine=ro2/ro1
  sinp=y/ro2
  if (t-ro1>=0.0d0) then
     call steady_oscillate(t-ro1,om,wp,p_p,p_pp,p_ppp,p_pppp,0)
     psi=1.0d0-exp(-alpha*ro1**5);
     Hx_exact=-(psi/(ro1**2)/c)*sine*(p_p+p_pp*ro1/c)*sinp;
  else
     Hx_exact=0.0d0
  end if
end function Hx_exact


real*8 function Hy_exact(x,y,z,t,alpha,om,c)
  implicit none
  real*8, intent(in) :: x,y,z,t,alpha,om,c
  real*8 ro1,ro2,cosine,sine,wp,p_p,p_pp,p_ppp,p_pppp,psi,E_r,E_theta,Eperp,cosp
  ro1=sqrt(x**2+y**2+z**2)
  ro2=sqrt(x**2+y**2)
  cosine=z/ro1;sine=ro2/ro1
  cosp=x/ro2
  if (t-ro1>=0.0d0) then
     call steady_oscillate(t-ro1,om,wp,p_p,p_pp,p_ppp,p_pppp,0)
     psi=1.0d0-exp(-alpha*ro1**5);
     Hy_exact=(psi/(ro1**2)/c)*sine*(p_p+p_pp*ro1/c)*cosp;
  else
     Hy_exact=0.0d0
  end if
end function Hy_exact


real*8 function Hz_exact(x,y,z,t,alpha,om,c)
  implicit none
  real*8, intent(in) :: x,y,z,t,alpha,om,c
  Hz_exact=0.0d0;
end function Hz_exact

! ==== CURRENTS ====
   
real*8 function Jx_current(x,y,z,t,alpha,om,c)
  implicit none
  real*8, intent(in) :: x,y,z,t,alpha,om,c
  real*8 ro1,ro2,cosine,sine,wp,p_p,p_pp,p_ppp,p_pppp,psi,E_r,E_theta,Eperp,sinp,psir
  real*8 ttt,j_r,j_theta,Jphi,jperp
  integer wr;
  ro1=dsqrt(x**2+y**2+z**2)
  ro2=dsqrt(x**2+y**2)
  cosine=z/ro1;sine=ro2/ro1 
  if (t-ro1>=0.0d0) then 
     call steady_oscillate(t-ro1,om,wp,p_p,p_pp,p_ppp,p_pppp,0)
     psir=5.0d0*alpha*(ro1**4)*exp(-alpha*(ro1**5))
     j_r=psir*cosine*(p_p+p_pp*ro1/c+p_ppp*ro1**2/(c**2))/(ro1**2) 
     j_theta=-psir*sine*(p_p+p_pp*ro1/c)/(ro1**2)/c
     jperp=j_r*sine+j_theta*cosine
     Jx_current=(jperp*x/ro2)*c/4.0d0/Pi 
  else
     Jx_current=0.0d0
  end if
end function Jx_current


real*8 function Jy_current(x,y,z,t,alpha,om,c)
  implicit none
  real*8, intent(in) :: x,y,z,t,alpha,om,c
  real*8 ro1,ro2,cosine,sine,wp,p_p,p_pp,p_ppp,p_pppp,psi,E_r,E_theta,Eperp,sinp,psir
  real*8 ttt,j_r,j_theta,Jphi,jperp
  ro1=sqrt(x**2+y**2+z**2)
  ro2=sqrt(x**2+y**2)
  cosine=z/ro1;sine=ro2/ro1
  if (t-ro1>=0.0d0) then
     call steady_oscillate(t-ro1,om,wp,p_p,p_pp,p_ppp,p_pppp,0)
     psir=5.0d0*alpha*(ro1**4)*exp(-alpha*(ro1**5))
     j_r=psir*cosine*(p_p+p_pp*ro1/c+p_ppp*ro1**2/(c**2))/(ro1**2) 
     j_theta=-psir*sine*(p_p+p_pp*ro1/c)/(ro1**2)/c
     jperp=j_r*sine+j_theta*cosine
     Jy_current=jperp*y/ro2*c/4.0d0/Pi 
  else
     Jy_current=0.0d0
  end if
end function Jy_current


real*8 function Jz_current(x,y,z,t,alpha,om,c)
  implicit none
  real*8, intent(in) :: x,y,z,t,alpha,om,c
  real*8 ro1,ro2,cosine,sine,wp,p_p,p_pp,p_ppp,p_pppp,psi,E_r,E_theta,Eperp,sinp,psir
  real*8 ttt,j_r,j_theta,Jphi,jperp
  ro1=sqrt(x**2+y**2+z**2)
  ro2=sqrt(x**2+y**2)
  cosine=z/ro1;sine=ro2/ro1
  if (t-ro1>=0.0d0) then
     call steady_oscillate(t-ro1,om,wp,p_p,p_pp,p_ppp,p_pppp,0)
     psir=5.0d0*alpha*(ro1**4)*exp(-alpha*(ro1**5))
     j_r=psir*cosine*(p_p+p_pp*ro1/c+p_ppp*ro1**2/(c**2))/(ro1**2) 
     j_theta=-psir*sine*(p_p+p_pp*ro1/c)/(ro1**2)/c
     Jz_current=(j_r*cosine-j_theta*sine)*c/4.0d0/Pi
  else
     Jz_current=0.0d0
  end if
end function Jz_current

! ==== AUXILIARY FUNCTION ====

subroutine steady_oscillate(t,om,wp,p_p,p_pp,p_ppp,p_pppp,Qtst)
  implicit none
  real*8, intent(in) :: t,om;
  real*8, intent(out) :: wp,p_p,p_pp,p_ppp,p_pppp;
  real*8 om1,om2
  integer :: Qtst;
  om1=om/2.0d0;om2=om1/dsqrt(2.0d0);
  wp = dsin(om1 * t) ** 4 + dsin(om2 * t) ** 4;
  p_p = 0.4D1 * dsin(om1 * t) ** 3 * dcos(om1 * t) * om1 + 0.4D1 * dsin(om2 * t) ** 3 * dcos(om2 * t) * om2;
  p_pp = 0.12D2 * dsin(om1 * t) ** 2 * dcos(om1 * t) ** 2 * om1 ** 2 - 0.4D1 * dsin(om1 * t) ** 4 * om1 ** 2 + 0.12D2 * dsin(om2 * t) ** 2 * dcos(om2 * t) ** 2 * om2 ** 2 - 0.4D1 * dsin(om2 * t) ** 4 * om2 ** 2;
  p_ppp = 0.24D2 * dsin(om1 * t) * dcos(om1 * t) ** 3 * om1 ** 3 - 0.40D2 * dsin(om1 * t) ** 3 * dcos(om1 * t) * om1 ** 3 +0.24D2 * dsin(om2 * t) * dcos(om2 * t) ** 3 * om2 ** 3 - 0.40D2 * dsin(om2 * t) ** 3* dcos(om2 * t) * om2 ** 3;
  p_pppp = 0.24D2 * dcos(om1 * t) ** 4 * om1 ** 4 - 0.192D3 * dsin(om1 * t) ** 2 * dcos(om1 * t) ** 2 * om1 ** 4 + 0.40D2 * dsin(om1 * t) ** 4*om1 ** 4+0.24D2*dcos(om2 * t) ** 4 * om2 ** 4 - 0.192D3*dsin(om2 * t) ** 2 * dcos(om2 * t) ** 2 * om2 ** 4 +0.40D2 * dsin(om2*t) ** 4 * om2 ** 4;
end subroutine steady_oscillate 
      
end module testsolution
