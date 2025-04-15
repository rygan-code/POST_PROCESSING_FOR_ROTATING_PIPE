

!-------------Viscous term -----------------------------------------------------------------
subroutine cal_viscous_term(sigma_x,sigma_y,sigma_z)
  Use flow_data
  implicit none

  real(kind=OCFD_REAL_KIND),allocatable,dimension(:,:,:,:)::  Ev1,Ev2,Ev3
  real(kind=OCFD_REAL_KIND),allocatable,dimension(:,:,:)::  uk,ui,us,vk,vi,vs,wk,wi,ws,Tk,Ti,Ts
  real(kind=OCFD_REAL_KIND)::  div, ux,uy,uz,vx,vy,vz,wx,wy,wz,Tx,Ty,Tz,Amu1,Amuk, &
    s11,s12,s13,s22,s23,s33,E1,E2,E3
  real(kind=OCFD_REAL_KIND):: s11d, s12d, s13d, s22d, s23d, s33d, Iso, up, down,delta
  real(kind=OCFD_REAL_KIND),parameter:: Prt=0.9d0, tmp2_3=2.d0/3.d0, cw = 0.325d0
  real*8,intent(out),dimension(nx,ny,nz) :: sigma_x,sigma_y,sigma_z
  integer:: i,j,k,m
  character(len=50) :: filename
!------------------------------------------------------------------------

  allocate(Ev1(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP,3), &
    Ev2(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP,3), &
    Ev3(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP,3) )
  allocate(uk(nx,ny,nz),ui(nx,ny,nz),us(nx,ny,nz),vk(nx,ny,nz),vi(nx,ny,nz),vs(nx,ny,nz), &
    wk(nx,ny,nz),wi(nx,ny,nz),ws(nx,ny,nz),Tk(nx,ny,nz),Ti(nx,ny,nz),Ts(nx,ny,nz))
  call comput_Amu         ! comput viscous coefficient (by using Sutherland eq.)

  call OCFD_dx0(u,uk,Scheme%Scheme_Vis)
  call OCFD_dx0(v,vk,Scheme%Scheme_Vis)
  call OCFD_dx0(w,wk,Scheme%Scheme_Vis)
  call OCFD_dx0(T,Tk,Scheme%Scheme_Vis)
  call OCFD_dy0(u,ui,Scheme%Scheme_Vis)
  call OCFD_dy0(v,vi,Scheme%Scheme_Vis)
  call OCFD_dy0(w,wi,Scheme%Scheme_Vis)
  call OCFD_dy0(T,Ti,Scheme%Scheme_Vis)
  call OCFD_dz0(u,us,Scheme%Scheme_Vis)
  call OCFD_dz0(v,vs,Scheme%Scheme_Vis)
  call OCFD_dz0(w,ws,Scheme%Scheme_Vis)
  call OCFD_dz0(T,Ts,Scheme%Scheme_Vis)

!-------------------------------------------------------------

  do k=1,nz
    do j=1,ny
      do i=1,nx
        ux=uk(i,j,k)*Akx(i,j,k)+ui(i,j,k)*Aix(i,j,k)+us(i,j,k)*Asx(i,j,k)
        vx=vk(i,j,k)*Akx(i,j,k)+vi(i,j,k)*Aix(i,j,k)+vs(i,j,k)*Asx(i,j,k)
        wx=wk(i,j,k)*Akx(i,j,k)+wi(i,j,k)*Aix(i,j,k)+ws(i,j,k)*Asx(i,j,k)
        Tx=Tk(i,j,k)*Akx(i,j,k)+Ti(i,j,k)*Aix(i,j,k)+Ts(i,j,k)*Asx(i,j,k)

        uy=uk(i,j,k)*Aky(i,j,k)+ui(i,j,k)*Aiy(i,j,k)+us(i,j,k)*Asy(i,j,k)
        vy=vk(i,j,k)*Aky(i,j,k)+vi(i,j,k)*Aiy(i,j,k)+vs(i,j,k)*Asy(i,j,k)
        wy=wk(i,j,k)*Aky(i,j,k)+wi(i,j,k)*Aiy(i,j,k)+ws(i,j,k)*Asy(i,j,k)
        Ty=Tk(i,j,k)*Aky(i,j,k)+Ti(i,j,k)*Aiy(i,j,k)+Ts(i,j,k)*Asy(i,j,k)

        uz=uk(i,j,k)*Akz(i,j,k)+ui(i,j,k)*Aiz(i,j,k)+us(i,j,k)*Asz(i,j,k)
        vz=vk(i,j,k)*Akz(i,j,k)+vi(i,j,k)*Aiz(i,j,k)+vs(i,j,k)*Asz(i,j,k)
        wz=wk(i,j,k)*Akz(i,j,k)+wi(i,j,k)*Aiz(i,j,k)+ws(i,j,k)*Asz(i,j,k)
        Tz=Tk(i,j,k)*Akz(i,j,k)+Ti(i,j,k)*Aiz(i,j,k)+Ts(i,j,k)*Asz(i,j,k)
        div=ux+vy+wz

        Amu1=Amu(i,j,k)            !  In this version, turbulence model is not supported !
        Amuk=Cp*Amu(i,j,k)/Pr

        s11=(2.d0*ux-tmp2_3*div)* Amu1          ! tmp2_3=2.d0/3.d0
        s22=(2.d0*vy-tmp2_3*div)* Amu1
        s33=(2.d0*wz-tmp2_3*div)* Amu1

        s12=(uy+vx)*Amu1
        s13=(uz+wx)*Amu1
        s23=(vz+wy)*Amu1

        Ev1(i,j,k,1)=Akx1(i,j,k)*s11+Aky1(i,j,k)*s12+Akz1(i,j,k)*s13
        Ev1(i,j,k,2)=Akx1(i,j,k)*s12+Aky1(i,j,k)*s22+Akz1(i,j,k)*S23
        Ev1(i,j,k,3)=Akx1(i,j,k)*s13+Aky1(i,j,k)*s23+Akz1(i,j,k)*s33

        Ev2(i,j,k,1)=Aix1(i,j,k)*s11+Aiy1(i,j,k)*s12+Aiz1(i,j,k)*s13
        Ev2(i,j,k,2)=Aix1(i,j,k)*s12+Aiy1(i,j,k)*s22+Aiz1(i,j,k)*S23
        Ev2(i,j,k,3)=Aix1(i,j,k)*s13+Aiy1(i,j,k)*s23+Aiz1(i,j,k)*s33

        Ev3(i,j,k,1)=Asx1(i,j,k)*s11+Asy1(i,j,k)*s12+Asz1(i,j,k)*s13
        Ev3(i,j,k,2)=Asx1(i,j,k)*s12+Asy1(i,j,k)*s22+Asz1(i,j,k)*S23
        Ev3(i,j,k,3)=Asx1(i,j,k)*s13+Asy1(i,j,k)*s23+Asz1(i,j,k)*s33

      enddo
    enddo
  enddo


!------  x,y,z direction
  do m=1,3
    call OCFD_dx0(Ev1(1-LAP,1-LAP,1-LAP,m),uk,Scheme%Scheme_Vis)
    call OCFD_dy0(Ev2(1-LAP,1-LAP,1-LAP,m),ui,Scheme%Scheme_Vis)
    call OCFD_dz0(Ev3(1-LAP,1-LAP,1-LAP,m),us,Scheme%Scheme_Vis)
    if(m.eq.0)then
      sigma_x = (uk(i,j,k)+ui(i,j,k)+us(i,j,k))*Ajac(i,j,k)
    else if(m.eq.1)then
      sigma_y = (uk(i,j,k)+ui(i,j,k)+us(i,j,k))*Ajac(i,j,k)
    else if(m.eq.2)then
      sigma_z = (uk(i,j,k)+ui(i,j,k)+us(i,j,k))*Ajac(i,j,k)
    end if
  enddo

  deallocate( Ev1,Ev2,Ev3,uk,ui,us,vk,vi,vs,wk,wi,ws,Tk,Ti,Ts)
end

!----------------------------------------------------------
! viscous coefficient: Sutherland Eq.
subroutine comput_Amu
  use flow_data
  implicit none
  real(kind=OCFD_REAL_KIND):: Tsb,tmpR
  integer:: i,j,k
  Tsb=110.4d0/Ref_Amu_T0
  TmpR=1.d0/Re
  do k=1,nz
    do j=1,ny
      do i=1,nx
        Amu(i,j,k)=TmpR*(1.d0+Tsb)*sqrt(T(i,j,k)**3)  /(Tsb+T(i,j,k))
      enddo
    enddo
  enddo
end