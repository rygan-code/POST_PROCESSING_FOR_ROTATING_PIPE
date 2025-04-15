subroutine cal_NL_uiuj(uiuj,NL_uiuj)
    use flow_para
    use flow_data
    implicit none
    real*8, dimension(nx,ny,nz,3,3),intent(in) :: uiuj
    real*8, dimension(nx,ny,nz,3,3),intent(out) :: NL_uiuj
    real*8, dimension(nx,ny,nz) :: uiuj_x, uiuj_y, uiuj_z

    call cal_grad(uiuj(:,:,:,1,1),uiuj_x,uiuj_y,uiuj_z)
    NL_uiuj(:,:,:,1,1) = u*uiuj_x+v*uiuj_y+w*uiuj_z
    call cal_grad(uiuj(:,:,:,1,2),uiuj_x,uiuj_y,uiuj_z)
    NL_uiuj(:,:,:,1,2) = u*uiuj_x+v*uiuj_y+w*uiuj_z
    call cal_grad(uiuj(:,:,:,1,3),uiuj_x,uiuj_y,uiuj_z)
    NL_uiuj(:,:,:,1,3) = u*uiuj_x+v*uiuj_y+w*uiuj_z
    call cal_grad(uiuj(:,:,:,2,2),uiuj_x,uiuj_y,uiuj_z)
    NL_uiuj(:,:,:,2,2) = u*uiuj_x+v*uiuj_y+w*uiuj_z
    call cal_grad(uiuj(:,:,:,2,3),uiuj_x,uiuj_y,uiuj_z)
    NL_uiuj(:,:,:,2,3) = u*uiuj_x+v*uiuj_y+w*uiuj_z
    call cal_grad(uiuj(:,:,:,3,3),uiuj_x,uiuj_y,uiuj_z)
    NL_uiuj(:,:,:,3,3) = u*uiuj_x+v*uiuj_y+w*uiuj_z
    NL_uiuj(:,:,:,2,1) = NL_uiuj(:,:,:,1,2)
    NL_uiuj(:,:,:,3,1) = NL_uiuj(:,:,:,1,3)
    NL_uiuj(:,:,:,3,2) = NL_uiuj(:,:,:,2,3)
end subroutine

subroutine cal_Ps_uiuj(Ps_uiuj)
    use flow_para
    use flow_data
    implicit none
    real*8, dimension(nx,ny,nz) :: pressure, pressure_x, pressure_y, pressure_z
    real*8, dimension(nx,ny,nz,3,3),intent(out) :: Ps_uiuj

    pressure = 1.d0/(gamma*Ma*Ma)*d*T
    call cal_grad(pressure,pressure_x,pressure_y,pressure_z)
    Ps_uiuj(:,:,:,1,1) = pressure_x*u+pressure_x*u
    Ps_uiuj(:,:,:,1,2) = pressure_x*v+pressure_y*u
    Ps_uiuj(:,:,:,1,3) = pressure_x*w+pressure_z*u
    Ps_uiuj(:,:,:,2,1) = pressure_y*u+pressure_x*v
    Ps_uiuj(:,:,:,2,2) = pressure_y*v+pressure_y*v
    Ps_uiuj(:,:,:,2,3) = pressure_y*w+pressure_z*v
    Ps_uiuj(:,:,:,3,1) = pressure_z*u+pressure_x*w
    Ps_uiuj(:,:,:,3,2) = pressure_z*v+pressure_y*w
    Ps_uiuj(:,:,:,3,3) = pressure_z*w+pressure_z*w
end subroutine

subroutine cal_Vs_uiuj(Vs_uiuj)
    use flow_para
    use flow_data
    implicit none
    real*8, dimension(nx,ny,nz) :: sigma_x, sigma_y, sigma_z
    real*8, dimension(nx,ny,nz,3,3),intent(out) :: Vs_uiuj

    call comput_Amu
    call cal_viscous_term(sigma_x,sigma_y,sigma_z)
    Vs_uiuj(:,:,:,1,1) = sigma_x*u+sigma_x*u
    Vs_uiuj(:,:,:,1,2) = sigma_x*v+sigma_y*u
    Vs_uiuj(:,:,:,1,3) = sigma_x*w+sigma_z*u
    Vs_uiuj(:,:,:,2,1) = sigma_y*u+sigma_x*v
    Vs_uiuj(:,:,:,2,2) = sigma_y*v+sigma_y*v
    Vs_uiuj(:,:,:,2,3) = sigma_y*w+sigma_z*v
    Vs_uiuj(:,:,:,3,1) = sigma_z*u+sigma_x*w
    Vs_uiuj(:,:,:,3,2) = sigma_z*v+sigma_y*w
    Vs_uiuj(:,:,:,3,3) = sigma_z*w+sigma_z*w
    
end subroutine

subroutine cal_Tv_uiuj(Tv_uiuj)
    use flow_para
    use flow_data
    implicit none
    real*8, dimension(nx,ny,nz) :: dux, duy, duz
    real*8, dimension(nx,ny,nz) :: ux, uy, uz, vx, vy, vz, wx, wy, wz
    real*8, dimension(nx,ny,nz) :: sigma_x, sigma_y, sigma_z
    real*8, dimension(nx,ny,nz) :: pressure, pressure_x, pressure_y, pressure_z
    real*8, dimension(nx,ny,nz,3,3),intent(out) :: Tv_uiuj

    pressure = 1.d0/(gamma*Ma*Ma)*d*T
    call cal_grad(pressure,pressure_x,pressure_y,pressure_z)
    call comput_Amu
    call cal_viscous_term(sigma_x,sigma_y,sigma_z)
    call cal_grad(u,ux,uy,uz)
    call cal_grad(v,vx,vy,vz)
    call cal_grad(w,wx,wy,wz)

    dux =  -(ux*u+uy*v+uz*w) - pressure_x + sigma_x
    duy =  -(vx*u+vy*v+vz*w) - pressure_y + sigma_y
    duz =  -(wx*u+wy*v+wz*w) - pressure_z + sigma_z

    Tv_uiuj(:,:,:,1,1) = dux*u+dux*u
    Tv_uiuj(:,:,:,1,2) = dux*v+duy*u
    Tv_uiuj(:,:,:,1,3) = dux*w+duz*u
    Tv_uiuj(:,:,:,2,1) = duy*u+dux*v
    Tv_uiuj(:,:,:,2,2) = duy*v+duy*v
    Tv_uiuj(:,:,:,2,3) = duy*w+duz*v
    Tv_uiuj(:,:,:,3,1) = duz*u+dux*w
    Tv_uiuj(:,:,:,3,2) = duz*v+duy*w
    Tv_uiuj(:,:,:,3,3) = duz*w+duz*w
end subroutine

subroutine fourier_transform_x(nx,ny,nz,x,fin,fout)
    implicit none
    integer,intent(in) :: nx,ny,nz
    real*8, dimension(nx,ny,nz),intent(in) :: fin,x
    complex*16 :: coefficient
    complex*16, dimension(nx,ny,nz),intent(out) :: fout
    integer :: i,j,k
    real*8 :: pi = 3.14159265358979323846
    real*8 :: deltax,L,k_freq
    
    L = x(nx,1,1)-x(1,1,1)
    deltax = L/nx
    fout = (0.d0,0.d0)
    do i = 1, nx
        if (i <= nx/2+1) then
            k_freq = 2.d0*pi*dble(i-1)/L
        else
            k_freq = 2.d0*pi*dble(i-1-nx)/L
        end if
        do j = 1, nx
            coefficient = exp(-(0.d0,1.d0)*k_freq*real(j-1,kind=8))
            fout(i,:,:) = fout(i,:,:) + coefficient*fin(j,:,:)/nx
        end do
    end do

end subroutine

subroutine fourier_inverse_transform_x(nx,ny,nz,x,fin,fout)
    implicit none
    integer,intent(in) :: nx,ny,nz
    real*8, dimension(nx,ny,nz),intent(in) :: x
    complex*16 :: coefficient
    complex*16, dimension(nx,ny,nz),intent(in) :: fin
    complex*16, dimension(nx,ny,nz) :: temp
    real*8, dimension(nx,ny,nz),intent(out) :: fout
    integer :: i,j,k
    real*8 :: pi = 3.14159265358979323846
    real*8 :: deltax,L,k_freq
    
    L = x(nx,1,1)-x(1,1,1)
    deltax = L/nx
    fout = (0.d0,0.d0)
    do i = 1, nx
        if (i <= nx/2+1) then
            k_freq = 2.d0*pi*dble(i-1)/L
        else
            k_freq = 2.d0*pi*dble(i-1-nx)/L
        end if
        do j = 1, nx
            coefficient = exp((0.d0,1.d0)*k_freq*real(j-1,kind=8))
            fout(j,:,:) = fout(j,:,:) + coefficient*fin(i,:,:)
        end do
    end do
    fout = real(temp)

end subroutine

subroutine gaussian_filter(nx,ny,nz,fin,fout,x,y,z,filter_length)
    implicit none
    integer,intent(in) :: nx,ny,nz
    real*8, dimension(nx,ny,nz),intent(in) :: x,y,z
    real*8, dimension(nx,ny,nz) :: fin,fout
    real*8, dimension(ny,nz) :: filter_core
    complex*16, dimension(nx,ny,nz) :: fin_fourier,fout_fourier
    integer :: i,j,k,ii,jj,kk
    real*8 :: sum,weight,filter_length
    real*8 :: pi = 3.14159265358979323846,temp

    call fourier_transform_x(nx,ny,nz,x,fin,fin_fourier)
    temp = (6.0d0/filter_length**2/pi)**0.5d0
    do j = 1, ny
        do k = 1, nz
            do i = 1, nx
                filter_core = temp**2*exp(-((y(1,:,:)-y(1,j,k))**2+(z(1,:,:)-z(1,j,k))**2)/filter_length**2)*&
                exp(-i**2*filter_length**2/24)
                fout_fourier(i,j,k) = sum(fin_fourier(i,:,:)*filter_core)
            end do
        end do
    end do
    call fourier_inverse_transform_x(nx,ny,nz,x,fout_fourier,fout)
    
end subroutine

subroutine test_uiuj!(Tv_uiuj,Vs_uiuj,Ps_uiuj,NL_uiuj)
    use flow_para
    use flow_data
    implicit none
    real*8, dimension(nx,ny,nz) :: Tv_uiuj,Vs_uiuj,Ps_uiuj,NL_uiuj,uiuj,test
    real*8, dimension(nx,ny,nz) :: uu, uy, uz
    complex*16, dimension(nx,ny,nz) :: ux
    integer :: i,j,k
    real*8 :: pi = 3.14159265358979323846

    do i = 1 ,nx
        uu(i,:,:) = sin(2*pi*real(i)/nx)
    end do
    call fourier_transform_x(nx,ny,nz,Axx,uu,ux)
    call fourier_inverse_transform_x(nx,ny,nz,Axx,ux,uy)
    print *, uu(:,1,1)
    print *, ux(:,1,1)
    print *, uy(:,1,1)


end subroutine