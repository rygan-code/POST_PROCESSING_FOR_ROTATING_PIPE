subroutine biot_savart(nx,ny,nz,x,y,z,wx,wy,wz,Jac1,ux,uy,uz)
implicit none

integer, intent(in) :: nx,ny,nz
real*8 , intent(in),dimension(nx,ny,nz):: x,y,z,wx,wy,wz,Jac1
real*8, intent(out),dimension(nx,ny,nz):: ux,uy,uz
real*8, dimension(nx,ny,nz):: rx,ry,rz,r
integer :: i,j,k
  do i = 1, nx  ! vorticity in (i,j,k)
    do j = 1, ny
      do k = 1, nz
        rx = mod(x-x(i,j,k)+x(nx,1,1)/2,x(nx,1,1))-x(nx,1,1)/2
        ry = y-y(i,j,k)
        rz = z-z(i,j,k)
        r = sqrt(rx**2+ry**2+rz**2)
        ux = ux + (wy*rz-wz*ry)/r**3*Jac1(i,j,k)
        uy = uy + (wz*rx-wx*rz)/r**3*Jac1(i,j,k)
        uz = uz + (wx*ry-wy*rx)/r**3*Jac1(i,j,k)

        ! far field
        ux = ux + 2*1.202*(wy*rz-wz*ry)/r**3*Jac1(i,j,k)/x(nx,1,1)**3
        uy = uy + 2*1.202*(wz*.0-wx*rz)/r**3*Jac1(i,j,k)/x(nx,1,1)**3
        uz = uz + 2*1.202*(wx*ry-wy*.0)/r**3*Jac1(i,j,k)/x(nx,1,1)**3
      end do
    end do
  end do
end subroutine

subroutine poisson_solver_diffusion(source,phi)
use flow_para
use flow_data
implicit none
  real*8 ,intent(in),dimension(nx,ny,nz) :: source
  real*8, dimension(nx,ny,nz) :: phix,phiy,phiz,phixx,phiyy,phizz,residual
  real*8, dimension(nx,ny,nz) :: phixk,phixi,phixs,phiyk,phiyi,phiys,phizk,phizi,phizs
  real*8, dimension(nx,ny,nz),intent(out) :: phi
  real*8 :: delta_t=0.001
  integer :: i,max_iter=10000

  do i = 1, max_iter
    call cal_grad(phi,phix,phiy,phiz)
    call OCFD_dx0(phix,phixk,Scheme%Scheme_Vis)
    call OCFD_dy0(phix,phixi,Scheme%Scheme_Vis)
    call OCFD_dz0(phix,phixs,Scheme%Scheme_Vis)
    call OCFD_dx0(phiy,phiyk,Scheme%Scheme_Vis)
    call OCFD_dy0(phiy,phiyi,Scheme%Scheme_Vis)
    call OCFD_dz0(phiy,phiys,Scheme%Scheme_Vis)
    call OCFD_dx0(phiz,phizk,Scheme%Scheme_Vis)
    call OCFD_dy0(phiz,phizi,Scheme%Scheme_Vis)
    call OCFD_dz0(phiz,phizs,Scheme%Scheme_Vis)
    phixx=phixk*Akx+phixi*Aix+phixs*Asx
    phiyy=phiyk*Aky+phiyi*Aiy+phiys*Asy
    phizz=phizk*Akz+phizi*Aiz+phizs*Asz

    phi(1:nx,2:ny-1,2:nz-1) = phi(1:nx,2:ny-1,2:nz-1) + &
        (source(1:nx,2:ny-1,2:nz-1)-phixx(1:nx,2:ny-1,2:nz-1)-phiyy(1:nx,2:ny-1,2:nz-1)-phizz(1:nx,2:ny-1,2:nz-1))*delta_t
    residual = source-phixx-phiyy-phizz
    print*,residual
    if (maxval(abs(residual)) < 1.d-6)then
      exit
      print*, 'iteration:',i,'residual:',maxval(abs(residual))
    end if
  end do


end subroutine

subroutine poisson_solver_gaussian_seidel(source,phi,initial_value)
use flow_para
use flow_data
implicit none
  real*8 ,intent(in),dimension(nx,ny,nz) :: source
  real*8, dimension(3,3) :: matrix,matrix_inv
  real*8, dimension(nx,ny,nz) :: lame_x,lame_y,lame_z,initial_value,phi_old
  real*8, dimension(nx,ny,nz),intent(out) :: phi
  real*8 :: tolerance,residual,ax_p,ax_m,ay_p,ay_m,az_p,az_m,det,ay
  integer :: iter,max_iter=10000,i,j,k
  
  do i=1, nx
    do j=1, ny
      do k=1, nz
        matrix = reshape([Akx(i,j,k), Aky(i,j,k), Akz(i,j,k),&
                          Aix(i,j,k), Aiy(i,j,k), Aiz(i,j,k),&
                          Asx(i,j,k), Asy(i,j,k), Asz(i,j,k)], [3, 3])
        det = matrix(1,1)*(matrix(2,2)*matrix(3,3) - matrix(2,3)*matrix(3,2)) &
            - matrix(1,2)*(matrix(2,1)*matrix(3,3) - matrix(2,3)*matrix(3,1)) &
            + matrix(1,3)*(matrix(2,1)*matrix(3,2) - matrix(2,2)*matrix(3,1)) 

        matrix_inv(1,1) =  (matrix(2,2)*matrix(3,3) - matrix(2,3)*matrix(3,2))
        matrix_inv(1,2) = -(matrix(1,2)*matrix(3,3) - matrix(1,3)*matrix(3,2))
        matrix_inv(1,3) =  (matrix(1,2)*matrix(2,3) - matrix(1,3)*matrix(2,2))
        
        matrix_inv(2,1) = -(matrix(2,1)*matrix(3,3) - matrix(2,3)*matrix(3,1))
        matrix_inv(2,2) =  (matrix(1,1)*matrix(3,3) - matrix(1,3)*matrix(3,1))
        matrix_inv(2,3) = -(matrix(1,1)*matrix(2,3) - matrix(1,3)*matrix(2,1))
        
        matrix_inv(3,1) =  (matrix(2,1)*matrix(3,2) - matrix(2,2)*matrix(3,1))
        matrix_inv(3,2) = -(matrix(1,1)*matrix(3,2) - matrix(1,2)*matrix(3,1))
        matrix_inv(3,3) =  (matrix(1,1)*matrix(2,2) - matrix(1,2)*matrix(2,1))
        matrix_inv = matrix_inv / det
        lame_x(i,j,k) = sqrt((matrix_inv(1,1))**2+(matrix_inv(1,2))**2+(matrix_inv(1,3))**2)
        lame_y(i,j,k) = sqrt((matrix_inv(2,1))**2+(matrix_inv(2,2))**2+(matrix_inv(2,3))**2)
        lame_z(i,j,k) = sqrt((matrix_inv(3,1))**2+(matrix_inv(3,2))**2+(matrix_inv(3,3))**2)
      end do
    end do
  end do
  print*,111
  ! if(isnan(sum(lame_x)+sum(lame_y)+sum(lame_z)))print*,1
  phi=initial_value
  tolerance = 1.d-6
  do iter=1,max_iter
  residual = 0

    phi_old = phi
    do i=1,nx
      do j=2,ny-1
        do k=1,nz
          if(j.gt.2.and.j.lt.ny-1)then
            ax_p=lame_y(mod(i+1-1,nx)+1,j,k)*lame_z(mod(i+1-1,nx)+1,j,k)/lame_x(mod(i+1-1,nx)+1,j,k)
            ax_m=lame_y(mod(i-1-1,nx)+1,j,k)*lame_z(mod(i-1-1,nx)+1,j,k)/lame_x(mod(i-1-1,nx)+1,j,k)
            ay_p=lame_x(i,j+1,k)*lame_z(i,j+1,k)/lame_y(i,j+1,k)
            ay_m=lame_x(i,j-1,k)*lame_z(i,j-1,k)/lame_y(i,j-1,k)
            az_p=lame_x(i,j,mod(k+1-1,nz)+1)*lame_y(i,j,mod(k+1-1,nz)+1)/lame_z(i,j,mod(k+1-1,nz)+1)
            az_m=lame_x(i,j,mod(k-1-1,nz)+1)*lame_y(i,j,mod(k-1-1,nz)+1)/lame_z(i,j,mod(k-1-1,nz)+1)

            phi(i,j,k) = (ax_p*phi_old(mod(i+2-1,nx)+1,j,k)+ax_m*phi_old(mod(i-2-1,nx)+1,j,k)+&         
                          ay_p*phi_old(i,j+2,k)+ay_m*phi_old(i,j-2,k)+&       
                          az_p*phi_old(i,j,mod(k+2-1,nz)+1)+az_m*phi_old(i,j,mod(k-2-1,nz)+1)-&
                          4*lame_x(i,j,k)*lame_y(i,j,k)*lame_z(i,j,k)*source(i,j,k)) &
                          /((ax_p+ax_m)+(ay_p+ay_m)+(az_p+az_m))
          else if(j.eq.2)then            
            ax_p=lame_y(mod(i+1-1,nx)+1,j,k)*lame_z(mod(i+1-1,nx)+1,j,k)/lame_x(mod(i+1-1,nx)+1,j,k)
            ax_m=lame_y(mod(i-1-1,nx)+1,j,k)*lame_z(mod(i-1-1,nx)+1,j,k)/lame_x(mod(i-1-1,nx)+1,j,k)
            ay_p=lame_x(i,j+2,k)*lame_z(i,j+2,k)/lame_y(i,j+2,k)
            ay  =lame_x(i,j+1,k)*lame_z(i,j+1,k)/lame_y(i,j+1,k)
            ay_m=lame_x(i,j,k)*lame_z(i,j,k)/lame_y(i,j,k)
            az_p=lame_x(i,j,mod(k+1-1,nz)+1)*lame_y(i,j,mod(k+1-1,nz)+1)/lame_z(i,j,mod(k+1-1,nz)+1)
            az_m=lame_x(i,j,mod(k-1-1,nz)+1)*lame_y(i,j,mod(k-1-1,nz)+1)/lame_z(i,j,mod(k-1-1,nz)+1)

            phi(i,j,k) = (ax_p*phi_old(mod(i+2-1,nx)+1,j,k)+ax_m*phi_old(mod(i-2-1,nx)+1,j,k)+&         
                          -ay_p*phi_old(i,j+3,k)+4*ay*phi_old(i,j+2,k)+(ay_p-3*ay_m)*phi_old(i,j+1,k)+3*ay_m*phi_old(i,j-1,k)+&       
                          az_p*phi_old(i,j,mod(k+2-1,nz)+1)+az_m*phi_old(i,j,mod(k-2-1,nz)+1)-&
                          4*lame_x(i,j,k)*lame_y(i,j,k)*lame_z(i,j,k)*source(i,j,k)) &
                          /((ax_p+ax_m)+4*ay+(az_p+az_m))
          else if(j.eq.ny-1)then            
            ax_p=lame_y(mod(i+1-1,nx)+1,j,k)*lame_z(mod(i+1-1,nx)+1,j,k)/lame_x(mod(i+1-1,nx)+1,j,k)
            ax_m=lame_y(mod(i-1-1,nx)+1,j,k)*lame_z(mod(i-1-1,nx)+1,j,k)/lame_x(mod(i-1-1,nx)+1,j,k)
            ay_p=lame_x(i,j,k)*lame_z(i,j,k)/lame_y(i,j,k)
            ay  =lame_x(i,j-1,k)*lame_z(i,j-1,k)/lame_y(i,j-1,k)
            ay_m=lame_x(i,j-2,k)*lame_z(i,j-2,k)/lame_y(i,j-2,k)
            az_p=lame_x(i,j,mod(k+1-1,nz)+1)*lame_y(i,j,mod(k+1-1,nz)+1)/lame_z(i,j,mod(k+1-1,nz)+1)
            az_m=lame_x(i,j,mod(k-1-1,nz)+1)*lame_y(i,j,mod(k-1-1,nz)+1)/lame_z(i,j,mod(k-1-1,nz)+1)

            phi(i,j,k) = (ax_p*phi_old(mod(i+2-1,nx)+1,j,k)+ax_m*phi_old(mod(i-2-1,nx)+1,j,k)+&         
                          -ay_m*phi_old(i,j-3,k)+4*ay*phi_old(i,j-2,k)+(ay_m-3*ay_p)*phi_old(i,j-1,k)+3*ay_p*phi_old(i,j+1,k)+&       
                          az_p*phi_old(i,j,mod(k+2-1,nz)+1)+az_m*phi_old(i,j,mod(k-2-1,nz)+1)-&
                          4*lame_x(i,j,k)*lame_y(i,j,k)*lame_z(i,j,k)*source(i,j,k)) &
                          /((ax_p+ax_m)+4*ay+(az_p+az_m))
          end if
        end do
      end do
    end do
    residual = sum(abs(phi-phi_old))
    print *,"interations",iter,"residual=",residual
    if(residual.lt.tolerance) then
      print*, "Converged in ", iter, " iterations"
      exit
    end if
  end do
end subroutine