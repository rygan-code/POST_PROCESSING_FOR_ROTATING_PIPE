
!-----------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------
subroutine   allocate_flow_data
  use flow_data
  implicit none
!-----------------------------------------------------------------------

! allocate flow data space
  allocate(f(nx,ny,nz,5),fn(nx,ny,nz,5), &
    du(nx,ny,nz,5),Amu(nx,ny,nz),Amu_t(nx,ny,nz))

  allocate(d(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP), &
    u(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP), &
    v(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP), &
    w(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP), &
    vor_x(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP), &
    vor_y(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP), &
    vor_z(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP), &
    Q_criterion(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP), &
    T(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP)  )

  d=1.d0; u=1.d0; v=1.d0; w=1.d0; T=1.d0       ! initial as 1.0
  Amu=0.d0                                     ! initial as 0
  du=0.d0                       ! initial as 0
  Amu_t=0.d0

!-----Coordinate and Jacobian coefficients   (Akx1=Akx/Ajac) -----------
  allocate( Axx(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP), &
    Ayy(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP), &
    Azz(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP), &
    Akx(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP), &
    Aky(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP), &
    Akz(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP), &
    Aix(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP), &
    Aiy(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP), &
    Aiz(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP), &
    Asx(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP), &
    Asy(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP), &
    Asz(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP), &
    Ajac(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP), &
    Akx1(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP), &
    Aky1(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP), &
    Akz1(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP), &
    Aix1(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP), &
    Aiy1(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP), &
    Aiz1(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP), &
    Asx1(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP), &
    Asy1(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP), &
    Asz1(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP))

  if(Scheme%Scheme_Invis == OCFD_Scheme_Hybrid) then
    allocate(Rhybrid(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP))
    Rhybrid=0.d0
  endif

!  ------initial as 1.0 --------
  Axx=1.d0;  Ayy=1.d0;  Azz=1.d0
  Akx=1.d0;  Aky=1.d0;  Akz=1.d0
  Aix=1.d0;  Aiy=1.d0;  Aiz=1.d0
  Asx=1.d0;  Asy=1.d0;  Asz=1.d0
  Ajac=1.d0
  Akx1=1.d0; Aky1=1.d0; Akz1=1.d0
  Aix1=1.d0; Aiy1=1.d0; Aiz1=1.d0
  Asx1=1.d0; Asy1=1.d0; Asz1=1.d0
end

subroutine   deallocate_flow_data
  use flow_data
  implicit none
  deallocate(f,fn,du,Amu,Amu_t,d,u,v,T,vor_x,vor_y,vor_z,Axx,Ayy,Akx,Aky,Aix,Aiy,Ajac,Akx1,Aky1,Aix1,Aiy1)
  if(Scheme%Scheme_Invis == OCFD_Scheme_Hybrid) then
    deallocate(Rhybrid)
  endif
end


subroutine   allocate_fourier_data
  use fourier_data
  implicit none
!-----------------------------------------------------------------------

  allocate( d_f(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP), &
  u_f(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP), &
  v_f(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP), &
  w_f(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP), &
  T_f(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP), &
  vor_xf(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP), &
  vor_yf(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP), &
  vor_zf(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP), &
  Q_criterion_f(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP))
end

subroutine   deallocate_fourier_data
  use fourier_data
  implicit none
  deallocate(d_f,u_f,v_f,w_f,T_f,vor_xf,vor_yf,vor_zf,Q_criterion_f)
end