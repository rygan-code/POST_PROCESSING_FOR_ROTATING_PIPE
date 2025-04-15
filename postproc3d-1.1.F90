program post
  use flow_para
  use flow_data
  use fourier_data
  
  implicit none
  real*8,allocatable,dimension(:,:,:,:):: u_wavelet,v_wavelet,w_wavelet,vor_x_wavelet,vor_y_wavelet,vor_z_wavelet
  real*8,allocatable,dimension(:,:,:):: temp,wxx,wxy,wxz,wyx,wyy,wyz,wzx,wzy,wzz,ux_induced,uy_induced,uz_induced
  real*8,allocatable,dimension(:,:)::d_av_x, u_av_x, v_av_x, w_av_x, T_av_x,energy,helicity
  real*8,allocatable,dimension(:)::d_av_xz, u_av_xz, v_av_xz, w_av_xz, T_av_xz
  real*8:: p00, Tw, R, nu, muw, u_tau, u_b, Re_b, Re_tau,get_mu, lambda, scale, position
  integer:: num,i,j,k,m,n
  
  call mpi_init(ierr)                       ! initial of MPI
  call mpi_comm_rank(MPI_COMM_WORLD,my_id,ierr)   ! get my_id
  if(id.eq.0)print*, "Post-processing code for OpenCFD2,  output tecplot-type data"
  ! read(*,*) nx,ny,nz
  open(66,file='grid2d.in')
  read(66,*)
  read(66,*)
  read(66,*) nx, nz, ny
  close(66)

  print*, " nx= ",nx," ny= ",ny," nz= ",nz
  LAP = 0
  call   allocate_flow_data        ! f, fn, d,u,v,w,T, Axx,Ayy,....,Ajac

  ! print*, "Please input Tw"
  ! read(*,*) Tw
  Tw=1

  call set_default_parameters
  call set_parameters
  p00=1.d0/(gamma*Ma*Ma)


  call read_mesh(nx,ny,nz,Axx,Ayy,Azz)
  call comput_Jacobian3d

  open(77,file="opencfd.dat",form="unformatted")
  read(77) Istep
  read(77) tt
  if(id.eq.0) print*, "Istep, tt=", Istep, tt

  call read3d(77,nx,ny,nz,d)
  call read3d(77,nx,ny,nz,u)
  call read3d(77,nx,ny,nz,v)
  call read3d(77,nx,ny,nz,w)
  call read3d(77,nx,ny,nz,T)
  close(77)

  call cal_vor
  ! allocate(wxx(nx,ny,nz),& 
  !          wxy(nx,ny,nz),&
  !          wxz(nx,ny,nz),&
  !          wyx(nx,ny,nz),&
  !          wyy(nx,ny,nz),&
  !          wyz(nx,ny,nz),&
  !          wzx(nx,ny,nz),&
  !          wzy(nx,ny,nz),&
  !          wzz(nx,ny,nz))
  ! allocate(ux_induced(nx,ny,nz),uy_induced(nx,ny,nz),uz_induced(nx,ny,nz))
  ! call cal_grad(vor_x,wxx,wxy,wxz)
  ! call cal_grad(vor_y,wyx,wyy,wyz)
  ! call cal_grad(vor_z,wzx,wzy,wzz)
  ! call poisson_solver(wzy-wyz,ux_induced)
  ! call poisson_solver(wxz-wzx,uy_induced)
  ! call poisson_solver(wyx-wxy,uz_induced)
  ! call write_3D(nx,ny,nz,Axx,Ayy,Azz,d,ux_induced,uy_induced,uz_induced,T,P00,Ma,Q_criterion)
  ! deallocate(wxx,wxy,wxz,wyx,wyy,wyz,wzx,wzy,wzz,ux_induced,uy_induced,uz_induced)

  ! average field in x-z plane
  if(id.eq.0)then
    allocate(d_av_x(ny,nz),u_av_x(ny,nz),v_av_x(ny,nz),w_av_x(ny,nz),T_av_x(ny,nz))
    allocate(d_av_xz(ny),u_av_xz(ny),v_av_xz(ny),  w_av_xz(ny),T_av_xz(ny))
    allocate(temp(nx,ny,nz))
    d_av_x = sum(d,dim=1)/nx
    d_av_xz = sum(d_av_x, dim=2)/nz

    u_av_x = sum(u,dim=1)/nx
    u_av_xz = sum(u_av_x, dim=2)/nz

    v_av_x = sum(v,dim=1)/nx
    v_av_xz = sum(v_av_x, dim=2)/nz

    w_av_x = sum(w,dim=1)/nx
    w_av_xz = sum(w_av_x, dim=2)/nz

    T_av_x = sum(T,dim=1)/nx
    T_av_xz = sum(T_av_x, dim=2)/nz

    ! Radius
    R = sqrt((Ayy(1,ny,1)-Ayy(1,1,1))**2 + (Azz(1,ny,1)-Azz(1,1,1))**2)
    ! Viscousity at T_ref or T0, ref density is 1.0
    nu = 1.d0/Re
    ! Viscousity at Tw
    muw = get_mu(nu,Tw,Ref_Amu_T0)
    ! wall shear velocity
    u_tau = sqrt(muw*(u_av_xz(ny-1)-u_av_xz(ny))/ &
            sqrt((Ayy(1,ny,1)-Ayy(1,ny-1,1))**2 + (Azz(1,ny,1)-Azz(1,ny-1,1))**2)/d_av_xz(ny))
    ! Bulk velocity
    u_b = sum(u_av_xz)/ny
    ! Bulk Re
    Re_b = 2*R*u_b/nu
    ! wall shear Re
    Re_tau = R*u_tau*d_av_xz(ny)/muw

    ! lambda
    lambda = 8.d0*muw*(u_av_xz(ny-1)-u_av_xz(ny))/ &
    sqrt((Ayy(1,ny,1)-Ayy(1,ny-1,1))**2 + (Azz(1,ny,1)-Azz(1,ny-1,1))**2)/(d_av_xz(ny)*u_b**2)

    print*, "Re_b = ", Re_b, "Re_tau = ", Re_tau, "lambda = ", lambda, "ub = ", u_b, "u_tau = ", u_tau
    print*, muw, (u_av_xz(ny-1)-u_av_xz(ny)), d_av_xz(ny), sqrt((Ayy(1,ny,1)-Ayy(1,ny-1,1))**2 + (Azz(1,ny,1)-Azz(1,ny-1,1))**2)

    open(99, file="yplus.dat")
    do i = ny,1,-1
      write(99,*)sqrt((Ayy(1,i,1)-Ayy(1,ny,1))**2+(Azz(1,i,1)-Azz(1,ny,1))**2)*&
      u_tau*d_av_xz(i)/get_mu(nu,T_av_xz(i),Ref_Amu_T0),u_av_xz(i)/u_tau
    end do
    close(99)

    open(100, file="velocity_profile.dat")
    do i = 1,ny/2
      write(100,*)sqrt((Ayy(1,i,1)-Ayy(1,1,1))**2+(Azz(1,i,1)-Azz(1,1,1))**2),u_av_xz(i)/u_tau
    end do
    close(100)
end if
call mpi_barrier(MPI_COMM_WORLD)


!-------------------------- Fourier analysis----------------
  call allocate_fourier_data
  allocate(u_wavelet(16,ny,ny/2,16))
  allocate(v_wavelet(16,ny,ny/2,16))
  allocate(w_wavelet(16,ny,ny/2,16))
  allocate(vor_x_wavelet(16,ny,ny/2,16))
  allocate(vor_y_wavelet(16,ny,ny/2,16))
  allocate(vor_z_wavelet(16,ny,ny/2,16))
  allocate(energy(ny,ny/2))
  allocate(helicity(ny,ny/2))
  ! call fourier_transform(nx,ny,nz,Axx,Ayy,Azz,Q_criterion,Q_criterion_f)
  !$omp parallel do
  do m=1, nx/16
    do n=1, nz/16
      call wavelet_transform(ny,Axx((m-1)*16+1:m*16,:,(n-1)*16+1:n*16),Ayy((m-1)*16+1:m*16,:,(n-1)*16+1:n*16),&
                            Azz((m-1)*16+1:m*16,:,(n-1)*16+1:n*16),u((m-1)*16+1:m*16,:,(n-1)*16+1:n*16),u_wavelet)
      if(isnan(sum(u_wavelet))) then
        print*, "u_wavelet is nan"
        print*, m,n
        stop
      end if
      call wavelet_transform(ny,Axx((m-1)*16+1:m*16,:,(n-1)*16+1:n*16),Ayy((m-1)*16+1:m*16,:,(n-1)*16+1:n*16),&
                            Azz((m-1)*16+1:m*16,:,(n-1)*16+1:n*16),v((m-1)*16+1:m*16,:,(n-1)*16+1:n*16),v_wavelet)
      call wavelet_transform(ny,Axx((m-1)*16+1:m*16,:,(n-1)*16+1:n*16),Ayy((m-1)*16+1:m*16,:,(n-1)*16+1:n*16),&
                            Azz((m-1)*16+1:m*16,:,(n-1)*16+1:n*16),w((m-1)*16+1:m*16,:,(n-1)*16+1:n*16),w_wavelet)
      call wavelet_transform(ny,Axx((m-1)*16+1:m*16,:,(n-1)*16+1:n*16),Ayy((m-1)*16+1:m*16,:,(n-1)*16+1:n*16),&
                            Azz((m-1)*16+1:m*16,:,(n-1)*16+1:n*16),vor_x((m-1)*16+1:m*16,:,(n-1)*16+1:n*16),vor_x_wavelet)
      call wavelet_transform(ny,Axx((m-1)*16+1:m*16,:,(n-1)*16+1:n*16),Ayy((m-1)*16+1:m*16,:,(n-1)*16+1:n*16),&
                            Azz((m-1)*16+1:m*16,:,(n-1)*16+1:n*16),vor_y((m-1)*16+1:m*16,:,(n-1)*16+1:n*16),vor_y_wavelet)
      call wavelet_transform(ny,Axx((m-1)*16+1:m*16,:,(n-1)*16+1:n*16),Ayy((m-1)*16+1:m*16,:,(n-1)*16+1:n*16),&
                            Azz((m-1)*16+1:m*16,:,(n-1)*16+1:n*16),vor_z((m-1)*16+1:m*16,:,(n-1)*16+1:n*16),vor_z_wavelet)
      print*, m,n
      do i = 1 , 16
        do j = 1, 16
          energy = energy + sum((u_wavelet*u_wavelet+v_wavelet*v_wavelet+&
                            w_wavelet*w_wavelet)*sum(Ajac((m-1)*16+i,:,(n-1)*16+j)))
          helicity = helicity + sum((u_wavelet*vor_x_wavelet+v_wavelet*vor_y_wavelet+&
                                w_wavelet*vor_z_wavelet)*sum(Ajac((m-1)*16+i,:,(n-1)*16+j)))
        end do
      end do
    end do
  end do

  open(44,file='invarians.dat')
  write(44,*) 'variables= sacle,position,energy,helicity'
  write(44,*) 'zone j=',ny, 'k=',ny
  do k=1,ny/2 !position
    do j=1,ny !scale
      scale=sqrt((Ayy(1,j,1)-Ayy(1,1,1))**2 + (Azz(1,j,1)-Azz(1,1,1))**2)
      position=sqrt(Ayy(1,j,1)**2 + Azz(1,j,1)**2)
      write(44,'(10f16.8)')scale,position,energy(j,k),helicity(j,k)
    end do
  end do
  close(44)
  
  ! call deallocate_fourier_data
  num=1
  do while( num .ne. 0)
    print*, "Plot 2D-plane/3D:  0  quit,  1 i-section,  2 j-section,   3 k-section , 4 3D   "
    read(*,*) num
    if(num .eq. 1) then
      call write_i(nx,ny,nz,Axx,Ayy,Azz,d,u,v,w,T,P00)
    else if (num .eq. 2) then
      call write_j(nx,ny,nz,Axx,Ayy,Azz,d,u,v,w,T,P00)
    else if (num .eq. 3) then
      call write_k(nx,ny,nz,Axx,Ayy,Azz,d,u,v,w,T,P00)
    else if (num .eq. 4) then
      call write_3D(nx,ny,nz,Axx,Ayy,Azz,d,u,v,w,T,P00,Ma,Q_criterion)
    endif
  enddo
  deallocate(d_av_x, u_av_x, v_av_x, w_av_x, T_av_x)
  deallocate(d_av_xz,u_av_xz,v_av_xz, w_av_xz,T_av_xz)
  deallocate(temp,u_wavelet,v_wavelet,w_wavelet,vor_x_wavelet,vor_y_wavelet,vor_z_wavelet,energy,helicity)

  call deallocate_flow_data
  call mpi_finalize(ierr)
end program post

!==================================
subroutine read3d(no,nx,ny,nz,u)
  implicit none
  integer:: no,nx,ny,nz,k
  real*8:: u(nx,ny,nz)
  print*, "read 3d data ..."
  do k=1,nz
    read (no)  u(:,:,k)
  enddo
end

!========================================
subroutine write_i(nx,ny,nz,x,y,z,d,u,v,w,T,P00)
  implicit  none
  integer:: nx,ny,nz,i,j,k
  real*8:: p00
  real*8,dimension(nx,ny,nz)::x,y,z,d,u,v,w,T
  print*, 'please input i ...'
  read(*,*) i

  open(101,file='flow2d_i.dat')
  write(101,*) 'variables= r,y,z,d,u,v,w,p,T'
  write(101,*) 'zone j=',ny, 'k=',nz
  do k=1,nz
    do j=1,ny
      write(101,'(10f16.8)') x(i,j,k),y(i,j,k),z(i,j,k),   &
        d(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k),p00*d(i,j,k)*T(i,j,k),  &
        T(i,j,k)
    enddo
  enddo
  close(101)
end

!c---------------------------------------------------
subroutine write_j(nx,ny,nz,x,y,z,d,u,v,w,T,P00)
  implicit  none
  integer:: nx,ny,nz,i,j,k
  real*8:: p00
  real*8,dimension(nx,ny,nz)::x,y,z,d,u,v,w,T
  print*, 'please input j ...'
  read(*,*) j

  open(44,file='flow2d_j.dat')
  write(44,*) 'variables= x,y,z,d,u,v,w,p,T'
  write(44,*) 'zone i=',nx, 'k=',nz
  do k=1,nz
    do i=1,nx
      write(44,'(10f16.8)') x(i,j,k),y(i,j,k),z(i,j,k),   &
        d(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k),p00*d(i,j,k)*T(i,j,k),  &
        T(i,j,k)
    enddo
  enddo
  close(44)
end
!c---------------------------------------------------
subroutine write_k(nx,ny,nz,x,y,z,d,u,v,w,T,P00)
  implicit  none
  integer:: nx,ny,nz,i,j,k
  real*8:: p00
  real*8,dimension(nx,ny,nz)::x,y,z,d,u,v,w,T
  print*, 'please input k ...'
  read(*,*) k

  open(44,file='flow2d_k.dat')
  write(44,*) 'variables= x,y,z,d,u,v,w,p,T'
  write(44,*) 'zone i=',nx, 'j=',ny
  do j=1,ny
    do i=1,nx
      write(44,'(10f16.8)') x(i,j,k),y(i,j,k),z(i,j,k),   &
        d(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k),p00*d(i,j,k)*T(i,j,k),  &
        T(i,j,k)
    enddo
  enddo
  close(44)
end

subroutine write_3D(nx,ny,nz,x,y,z,d,u,v,w,T,P00,Ma,Q_criterion)
  implicit  none
  integer:: ib, ie, jb,je,kb,ke, istep, jstep,kstep
  integer:: nx,ny,nz,i,j,k,Iflag_Q
  real*8:: p00,Ma
  real*8,dimension(nx,ny,nz)::x,y,z,d,u,v,w,T,Q_criterion

  print*, "please input ib, ie, jb, je, kb, ke, istep, jstep, kstep"
  ib = 1; ie = 512; jb = 1; je = 256; kb = 1; ke = 512; istep = 4; jstep = 4; kstep = 4
  ! read(*,*) ib, ie, jb, je, kb, ke, istep, jstep, kstep
  open(44,file='flow3d.dat')
  write(44,*) 'variables= x,y,z,d,u,v,w,p,T,Mach,Q_criterion'
  write(44,*) 'zone i=', (ie-ib)/istep +1 , 'j=', (je-jb)/jstep +1, 'k=', (ke-kb)/kstep+2
  do k=kb,ke,kstep
    do j=jb,je,jstep
      do i=ib,ie,istep
        write(44,'(11f16.8)') x(i,j,k),y(i,j,k),z(i,j,k),   &
          d(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k),p00*d(i,j,k)*T(i,j,k),  &
          T(i,j,k),Ma*U(i,j,k)/sqrt(T(i,j,k)),Q_criterion(i,j,k)
      enddo
    enddo
  enddo
  k = kb
  do j=jb,je,jstep
    do i=ib,ie,istep
      write(44,'(11f16.8)') x(i,j,k),y(i,j,k),z(i,j,k),   &
      d(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k),p00*d(i,j,k)*T(i,j,k),  &
      T(i,j,k),Ma*U(i,j,k)/sqrt(T(i,j,k)),Q_criterion(i,j,k)
    end do
  end do
  close(44)
end
!--------------------------------------------
subroutine read_mesh(nx,ny,nz,x,y,z)
  implicit none
  integer:: nx,ny,nz,mesh_type,i,j,k
  integer,parameter:: GRID1D=10, GRID2D_PLANE=20, GRID2D_AXIAL_SYMM=21, GRID3D=30, GRID_AND_JACOBIAN3D=31
  real*8:: x(nx,ny,nz),y(nx,ny,nz),z(nx,ny,nz)
  real*8,allocatable:: x1d(:),y1d(:),z1d(:),x2d(:,:),y2d(:,:)
  allocate(x1d(nx),y1d(ny),z1d(nz),x2d(nx,ny),y2d(nx,ny))
  print*, "read  mesh ......"
  print*, "pleas input mesh_type,   10: 1D mesh;  20: 2D-plane mesh;  21: 2D-AxialSymmetry mesh;  30: 3D mesh"
  ! read(*,*)  mesh_type
  mesh_type=30

  open(56,file='OCFD-grid.dat',form='unformatted')
  if(mesh_type ==GRID1D) then
    read(56) x1d
    read(56) y1d
    read(56) z1d
    do k=1,nz
      do j=1,ny
        do i=1,nx
          x(i,j,k)=x1d(i)
          y(i,j,k)=y1d(j)
          z(i,j,k)=z1d(k)
        enddo
      enddo
    enddo
  else if( mesh_type ==GRID2D_PLANE) then
    read(56) x2d
    read(56) y2d
    read(56) z1d
    do k=1,nz
      do j=1,ny
        do i=1,nx
          x(i,j,k)=x2d(i,j)
          y(i,j,k)=y2d(i,j)
          z(i,j,k)=z1d(k)
        enddo
      enddo
    enddo
  else if( mesh_type ==GRID2D_AXIAL_SYMM) then
    read(56) x2d
    read(56) y2d
    read(56) z1d
    do k=1,nz
      do j=1,ny
        do i=1,nx
          x(i,j,k)=x2d(i,j)
          y(i,j,k)=y2d(i,j)*cos(z1d(k))
          z(i,j,k)=y2d(i,j)*sin(z1d(k))
        enddo
      enddo
    enddo
  else
    call read3d(56,nx,ny,nz,x)
    call read3d(56,nx,ny,nz,y)
    call read3d(56,nx,ny,nz,z)
  endif
  close(56)
  deallocate(x1d,y1d,z1d,x2d,y2d)
end

function get_mu(mu_ref, T, T0) result(mu)
  implicit none
  real(8), intent(in) :: mu_ref, T, T0
  real(8) :: mu

  mu = mu_ref * sqrt(T**3)*(1.d0+110.4d0/T0)/(110.4d0/T0+T)
end function get_mu

subroutine cal_vor
  use flow_para
  use flow_data
  implicit none
  ! integer,intent(in)::nx,ny,nz
  real*8,dimension(nx,ny,nz)::ui,vi,wi,uk,vk,wk,us,ws,vs,ux,uy,uz,vx,vy,vz,wx,wy,wz
  integer::i,j,k

  call OCFD_dx0(u,uk,Scheme%Scheme_Vis)
  call OCFD_dx0(v,vk,Scheme%Scheme_Vis)
  call OCFD_dx0(w,wk,Scheme%Scheme_Vis)
  call OCFD_dy0(u,ui,Scheme%Scheme_Vis)
  call OCFD_dy0(v,vi,Scheme%Scheme_Vis)
  call OCFD_dy0(w,wi,Scheme%Scheme_Vis)
  call OCFD_dz0(u,us,Scheme%Scheme_Vis)
  call OCFD_dz0(v,vs,Scheme%Scheme_Vis)
  call OCFD_dz0(w,ws,Scheme%Scheme_Vis)

  ux=uk*Akx+ui*Aix+us*Asx
  vx=vk*Akx+vi*Aix+vs*Asx
  wx=wk*Akx+wi*Aix+ws*Asx
  uy=uk*Aky+ui*Aiy+us*Asy
  vy=vk*Aky+vi*Aiy+vs*Asy
  wy=wk*Aky+wi*Aiy+ws*Asy
  uz=uk*Akz+ui*Aiz+us*Asz
  vz=vk*Akz+vi*Aiz+vs*Asz
  wz=wk*Akz+wi*Aiz+ws*Asz
  Q_criterion=ux*vy+ux*wz+vy*wz -uy*vx-uz*wx-vz*wy  !! TK=Q=II(UX)
  vor_x=wy-vz
  vor_y=uz-wx
  vor_z=vx-uy
end subroutine cal_vor

subroutine cal_grad(fin,foutx,fouty,foutz)
  use flow_para
  use flow_data
  implicit none
  real*8,intent(in),dimension(nx,ny,nz)::fin
  real*8,dimension(nx,ny,nz)::fink,fini,fins
  real*8,dimension(nx,ny,nz),intent(out)::foutx,fouty,foutz

  integer::i,j,k

  call OCFD_dx0(fin,fink,Scheme%Scheme_Vis)
  call OCFD_dy0(fin,fini,Scheme%Scheme_Vis)
  call OCFD_dz0(fin,fins,Scheme%Scheme_Vis)

  foutx=fink*Akx+fini*Aix+fins*Asx
  fouty=fink*Aky+fini*Aiy+fins*Asy
  foutz=fink*Akz+fini*Aiz+fins*Asz
    
end subroutine cal_grad

! spectral analysis
subroutine fourier_transform(nx,ny,nz,x,y,z,fin,fout)
  integer:: nx,ny,nz,i,j,kmax
  real*8,dimension(nx,ny,nz)::fin,fout,x,y,z,r,dr
  real*8::rmin,rmax,dr_max,pi=3.14159265358979323846

  r = sqrt(y**2 + z**2); dr = r(1,1,2)-r(1,1,1)
  do i = 1 , ny-1
    dr(:,i,:) = r(:,i+1,:)-r(:,i,:)
  end do
  print*,r(1,1,1),r(1,ny,1)
  dr(:,ny,:) = dr(:,ny-1,:)
  dr_max = maxval(dr(1,:,1))
  rmin = r(1,1,1); rmax = r(1,ny,1); L = rmax-rmin
  kmax = min(int(pi/dr_max),ny)
  print*,kmax
  do i = 1 , kmax
    do j = 1, ny
      fout(:,i,:) = fout(:,i,:) + 2/L*fin(:,j,:)*sin(i*pi*(r(i,j,1)-rmin)/L)*dr(1,j,1)
    end do
  end do
end subroutine

subroutine wavelet_transform(ny,x,y,z,fin,fout)
  integer:: nx,ny,nz,i,j,kmax
  real*8,dimension(16,ny,16)::fin,x,y,z,r,dr
  real*8,dimension(16,ny,ny/2,16)::fout
  real*8,dimension(ny,ny,ny/2)::psi
  real*8::rmin,rmax,dr_max,pi=3.14159265358979323846,position,scale,mexican_hat_wavelet

  r = sqrt(y**2 + z**2)
  do i = 1 , ny-1   ! scale
    scale = (r(1,i+1,1)-r(1,1,1))
    do j = 1, ny/2    ! position
      position = r(1,j*2,1)
      do k = 1, ny
        psi(k,i,j) = mexican_hat_wavelet((r(1,k,1)-position)/scale)
        fout(:,i,j,:) = fout(:,i,j,:) + scale**(-0.5)*psi(k,i,j)*fin(:,k,:)
      end do
    end do
  end do
end subroutine

function mexican_hat_wavelet(x) result(psi)
  implicit none
  real(8), intent(in) :: x
  real(8) :: psi

  psi = (1-x**2)*exp(-0.5*x**2)
end function mexican_hat_wavelet

! Finite difference for viscous terms  (Centeral Schemes)

subroutine OCFD_dx0(f,fx, Num_Scheme)
  use flow_para
  implicit none
  integer Num_Scheme,i,j,k
  real(kind=OCFD_REAL_KIND)::   f(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP),fx(nx,ny,nz)
  real(kind=OCFD_REAL_KIND)::   a1,a2,a3,c1,c2,c3,c4
  a1=1.d0/(60.d0*hx)             ! 6th centeral scheme
  a2=-3.d0/(20.d0*hx)
  a3=3.d0/(4.d0*hx)

  c1=0.8d0/hx                    ! 8th centreral scheme
  c2=-0.2d0/hx
  c3=3.80952380952380952d-2/hx
  c4=-3.571428571428571428d-3/hx


!------Scheme for inner points  ( CD6, CD8)
  if(Num_Scheme == OCFD_Scheme_CD6 ) then
    do k=1,nz
      do j=1,ny
        do i=1,nx
          fx(i,j,k)=a1*(f(i+3,j,k)-f(i-3,j,k)) +a2*(f(i+2,j,k)-f(i-2,j,k)) +a3*(f(i+1,j,k)-f(i-1,j,k))
        enddo
      enddo
    enddo

  else if (Num_Scheme == OCFD_Scheme_CD8) then
    do k=1,nz
      do j=1,ny
        do i=1,nx
          fx(i,j,k)=c1*(f(i+1,j,k)-f(i-1,j,k)) +c2*(f(i+2,j,k)-f(i-2,j,k))  &
            +c3*(f(i+3,j,k)-f(i-3,j,k)) +c4*(f(i+4,j,k)-f(i-4,j,k))
        enddo
      enddo
    enddo
    ! print*,fx
  else
    print*, 'This Numerical Scheme is not supported in viscous terms !'
    print*, 'Only CD6 or CD8 can be used in viscous terms'
    stop
  endif

!---------Boundary Scheme ------- (low-order scheme)----
!---------- i- boundary  ---------------------------
  if(npx .eq. 0 .and. Para%Iperiodic_X .eq. 0) then
    do k=1,nz
      do j=1,ny
          fx(1,j,k)=a1*(f(4,j,k)-f(nx-2,j,k)) +a2*(f(3,j,k)-f(nx-1,j,k)) +a3*(f(2,j,k)-f(nx,j,k))
          fx(2,j,k)=a1*(f(5,j,k)-f(nx-1,j,k)) +a2*(f(4,j,k)-f(nx,j,k)) +a3*(f(3,j,k)-f(1,j,k))
          fx(3,j,k)=a1*(f(6,j,k)-f(nx,j,k)) +a2*(f(5,j,k)-f(1,j,k)) +a3*(f(4,j,k)-f(2,j,k))
      enddo
    enddo
    if(Num_Scheme == OCFD_Scheme_CD8) then
      do k=1,nz
        do j=1,ny
          fx(1,j,k)=c1*(f(2,j,k)-f(nx,j,k)) +c2*(f(3,j,k)-f(nx-1,j,k))  &
            +c3*(f(4,j,k)-f(nx-2,j,k)) +c4*(f(5,j,k)-f(nx-3,j,k))
          fx(2,j,k)=c1*(f(3,j,k)-f(1,j,k)) +c2*(f(4,j,k)-f(nx,j,k))  &
            +c3*(f(5,j,k)-f(nx-1,j,k)) +c4*(f(6,j,k)-f(nx-2,j,k))
          fx(3,j,k)=c1*(f(4,j,k)-f(2,j,k)) +c2*(f(5,j,k)-f(1,j,k))  &
            +c3*(f(6,j,k)-f(nx,j,k)) +c4*(f(7,j,k)-f(nx-1,j,k))
          fx(4,j,k)=c1*(f(5,j,k)-f(3,j,k)) +c2*(f(6,j,k)-f(2,j,k))  &
            +c3*(f(7,j,k)-f(1,j,k)) +c4*(f(8,j,k)-f(nx,j,k))
        enddo
      enddo
    endif
  endif
!--------- i+ boundary ------------------------
  if(npx .eq. npx0-1 .and. Para%Iperiodic_X .eq. 0) then
    do k=1,nz
      do j=1,ny
        fx(nx,j,k)=(f(nx-2,j,k)-4.d0*f(nx-1,j,k)  +3.d0*f(nx,j,k))/(2.d0*hx)  ! 2nd one-side scheme
        fx(nx-1,j,k)=(f(nx,j,k)-f(nx-2,j,k))/(2.d0*hx)                             ! 2nd centeral scheme
        fx(nx-2,j,k)=(8.d0*(f(nx-1,j,k)-f(nx-3,j,k)) - (f(nx,j,k)-f(nx-4,j,k)))/(12.d0*hx)   ! 4th central scheme
      enddo
    enddo
    if(Num_Scheme == OCFD_Scheme_CD8) then
      do k=1,nz
        do j=1,ny
          fx(nx-3,j,k)=a1*(f(nx,j,k)-f(nx-6,j,k)) +a2*(f(nx-1,j,k)-f(nx-5,j,k))  +a3*(f(nx-2,j,k)-f(nx-4,j,k))	  ! 6th centeral scheme
        enddo
      enddo
    endif
  endif


end


!c----------------------------------------------------------

subroutine OCFD_dy0(f,fy,Num_Scheme )
  use flow_para
  implicit none
  integer Num_Scheme,i,j,k
  real(kind=OCFD_REAL_KIND):: f(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP), fy(nx,ny,nz)
  real(kind=OCFD_REAL_KIND):: a1,a2,a3,c1,c2,c3,c4
  a1=1.d0/(60.d0*hy)
  a2=-3.d0/(20.d0*hy)
  a3=3.d0/(4.d0*hy)
  c1=0.8d0/hy
  c2=-0.2d0/hy
  c3=3.80952380952380952d-2/hy
  c4=-3.571428571428571428d-3/hy

!------Scheme for inner points  ( CD6, CD8)
  if(Num_Scheme == OCFD_Scheme_CD6 ) then
    do k=1,nz
      do j=1,ny
        do i=1,nx
          fy(i,j,k)=a1*(f(i,j+3,k)-f(i,j-3,k))  +a2*(f(i,j+2,k)-f(i,j-2,k)) +a3*(f(i,j+1,k)-f(i,j-1,k))
        enddo
      enddo
    enddo

  else if (Num_Scheme == OCFD_Scheme_CD8) then
    do k=1,nz
      do j=1,ny
        do i=1,nx
          fy(i,j,k)=c1*(f(i,j+1,k)-f(i,j-1,k)) +c2*(f(i,j+2,k)-f(i,j-2,k))  &
            +c3*(f(i,j+3,k)-f(i,j-3,k)) +c4*(f(i,j+4,k)-f(i,j-4,k))
        enddo
      enddo
    enddo
  else

    print*, 'This Numerical Scheme is not supported in viscous terms !'
    print*, 'Only CD6 or CD8 can be used in viscous terms'
    stop
  endif

!---------Boundary Scheme ------- (low-order scheme)----
!---------- j- boundary  ---------------------------
  if(npy .eq. 0 .and. Para%Iperiodic_Y .eq. 0) then
    do k=1,nz
      do i=1,nx
        fy(i,1,k)=(-3.d0*f(i,1,k)+4.d0*f(i,2,k)-f(i,3,k))/(2.d0*hy)           ! 2nd one-side scheme
        fy(i,2,k)=(f(i,3,k)-f(i,1,k))/(2.d0*hy)                             ! 2nd centeral scheme
        fy(i,3,k)=(8.d0*(f(i,4,k)-f(i,2,k)) - (f(i,5,k)-f(i,1,k)))/(12.d0*hy)   ! 4th central scheme
      enddo
    enddo
    if(Num_Scheme == OCFD_Scheme_CD8) then
      do k=1,nz
        do i=1,nx
          fy(i,4,k)=a1*(f(i,7,k)-f(i,1,k)) +a2*(f(i,6,k)-f(i,2,k))  +a3*(f(i,5,k)-f(i,3,k))	  ! 6th centeral scheme
        enddo
      enddo
    endif
  endif
!--------- j+ boundary ------------------------
  if(npy .eq. npy0-1 .and. Para%Iperiodic_Y .eq. 0) then
    do k=1,nz
      do i=1,nx
        fy(i,ny,k)=(f(i,ny-2,k)-4.d0*f(i,ny-1,k)  +3.d0*f(i,ny,k))/(2.d0*hy)  ! 2nd one-side scheme
        fy(i,ny-1,k)=(f(i,ny,k)-f(i,ny-2,k))/(2.d0*hy)                             ! 2nd centeral scheme
        fy(i,ny-2,k)=(8.d0*(f(i,ny-1,k)-f(i,ny-3,k)) - (f(i,ny,k)-f(i,ny-4,k)))/(12.d0*hy)   ! 4th central scheme
      enddo
    enddo
    if(Num_Scheme == OCFD_Scheme_CD8) then
      do k=1,nz
        do i=1,nx
          fy(i,ny-3,k)=a1*(f(i,ny,k)-f(i,ny-6,k)) +a2*(f(i,ny-1,k)-f(i,ny-5,k))  +a3*(f(i,ny-2,k)-f(i,ny-4,k))	  ! 6th centeral scheme
        enddo
      enddo
    endif
  endif

end


!c----------------------------------------------------------

subroutine OCFD_dz0(f,fz,Num_Scheme )
  use flow_para
  implicit none
  integer Num_Scheme,i,j,k
  real(kind=OCFD_REAL_KIND):: f(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP), fz(nx,ny,nz)
  real(kind=OCFD_REAL_KIND):: a1,a2,a3,c1,c2,c3,c4
  a1=1.d0/(60.d0*hz)
  a2=-3.d0/(20.d0*hz)
  a3=3.d0/(4.d0*hz)
  c1=0.8d0/hz
  c2=-0.2d0/hz
  c3=3.80952380952380952d-2/hz
  c4=-3.571428571428571428d-3/hz

!------Scheme for inner points  ( CD6, CD8)
  if(Num_Scheme == OCFD_Scheme_CD6 ) then
    do k=1,nz
      do j=1,ny
        do i=1,nx
          fz(i,j,k)=a1*(f(i,j,k+3)-f(i,j,k-3))  +a2*(f(i,j,k+2)-f(i,j,k-2)) +a3*(f(i,j,k+1)-f(i,j,k-1))
        enddo
      enddo
    enddo

  else if (Num_Scheme == OCFD_Scheme_CD8) then
    do k=1,nz
      do j=1,ny
        do i=1,nx
          fz(i,j,k)=c1*(f(i,j,k+1)-f(i,j,k-1)) +c2*(f(i,j,k+2)-f(i,j,k-2))  &
            +c3*(f(i,j,k+3)-f(i,j,k-3)) +c4*(f(i,j,k+4)-f(i,j,k-4))
        enddo
      enddo
    enddo
  else

    print*, 'This Numerical Scheme is not supported in viscous terms !'
    print*, 'Only CD6 or CD8 can be used in viscous terms'
    stop
  endif

!---------Boundary Scheme ------- (low-order scheme)----
!---------- k- boundary  ---------------------------
  if(npz .eq. 0 .and. Para%Iperiodic_Z .eq. 0) then
    do j=1,ny
      do i=1,nx
        fz(i,j,1)=(-3.d0*f(i,j,1)+4.d0*f(i,j,2)-f(i,j,3))/(2.d0*hz)           ! 2nd one-side scheme
        fz(i,j,2)=(f(i,j,3)-f(i,j,1))/(2.d0*hz)                             ! 2nd centeral scheme
        fz(i,j,3)=(8.d0*(f(i,j,4)-f(i,j,2)) - (f(i,j,5)-f(i,j,1)))/(12.d0*hz)   ! 4th central scheme
      enddo
    enddo
    if(Num_Scheme == OCFD_Scheme_CD8) then
      do j=1,ny
        do i=1,nx
          fz(i,j,4)=a1*(f(i,j,7)-f(i,j,1)) +a2*(f(i,j,6)-f(i,j,2))  +a3*(f(i,j,5)-f(i,j,3))	  ! 6th centeral scheme
        enddo
      enddo
    endif
  endif
!--------- k+ boundary ------------------------
  if(npz .eq. npz0-1 .and. Para%Iperiodic_Z .eq. 0) then
    do j=1,ny
      do i=1,nx
        fz(i,j,nz)=(f(i,j,nz-2)-4.d0*f(i,j,nz-1)  +3.d0*f(i,j,nz))/(2.d0*hz)  ! 2nd one-side scheme
        fz(i,j,nz-1)=(f(i,j,nz)-f(i,j,nz-2))/(2.d0*hz)                             ! 2nd centeral scheme
        fz(i,j,nz-2)=(8.d0*(f(i,j,nz-1)-f(i,j,nz-3)) - (f(i,j,nz)-f(i,j,nz-4)))/(12.d0*hz)   ! 4th central scheme
      enddo
    enddo
    if(Num_Scheme == OCFD_Scheme_CD8) then
      do j=1,ny
        do i=1,nx
          fz(i,j,nz-3)=a1*(f(i,j,nz)-f(i,j,nz-6)) +a2*(f(i,j,nz-1)-f(i,j,nz-5))  +a3*(f(i,j,nz-2)-f(i,j,nz-4))	  ! 6th centeral scheme
        enddo
      enddo
    endif
  endif

end