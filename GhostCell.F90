subroutine Jac_Ghost_boundary
   use flow_data
   implicit none
   integer:: i,j,k
   call Jac_Ghost_Extent_2nd    ! 2nd order exterpolation , in case Para%Ghost_Cell(1) .eq. 1
end

subroutine Flow_Ghost_boundary
   use flow_data
   implicit none
   call Flow_Ghost_Extent_2nd  ! 2nd exterpolation, in case Scheme%Scheme_boundary(:) .eq. 1 or 2
end

!-----------------------------------------------
! 2nd Exter-ploation Ghost Cell , in case Para%Ghost_Cell(:) ==1 or 2
subroutine Flow_Ghost_Extent_2nd
   use flow_data
   implicit none
   integer i,j,k,k1,i1,j1
   Real(kind=OCFD_REAL_KIND),parameter:: s1=0.8d0, s2=1.2d0       ! limter

   if (npy==0) then
      do k = 1,nz
         do j = 1-LAP,0
            do i = 1, nx
               j1 = 2-j

               u(i,j,k) = 2.d0*u(i,1,k)-u(i,j1,k)
               v(i,j,k) = 2.d0*v(i,1,k)-v(i,j1,k)
               w(i,j,k) = 2.d0*w(i,1,k)-w(i,j1,k)
               d(i,j,k) = 2.d0*d(i,1,k)-d(i,j1,k)
               T(i,j,k) = 2.d0*T(i,1,k)-T(i,j1,k)
            end do
         end do
      end do
   endif

   if (npy==npy0-1) then
      do k=1,nz
         do j=ny+1, ny+LAP
            do i=1, nx
               j1=2*ny-j

               u(i,j,k) = 2.d0*u(i,ny,k)-u(i,j1,k)
               v(i,j,k) = 2.d0*v(i,ny,k)-v(i,j1,k)
               w(i,j,k) = 2.d0*w(i,ny,k)-w(i,j1,k)
               d(i,j,k) = 2.d0*d(i,ny,k)-d(i,j1,k)
               T(i,j,k) = 2.d0*T(i,ny,k)-T(i,j1,k)
            enddo
         enddo
      enddo
   end if

   if (Para%Iperiodic_Z == 1 .and. Para%Z_mpi == 0) then
      if(npz==0) then   ! k-

         do k=1-LAP,0
            do j=1,ny
               do i=1,nx
                  k1=nz+k

                  u(i,j,k)=u(i,j,k1)
                  v(i,j,k)=v(i,j,k1)
                  w(i,j,k)=w(i,j,k1)

                  d(i,j,k)=d(i,j,k1)
                  T(i,j,k)=T(i,j,k1)

               enddo
            enddo
         enddo
      end if
      if (npz==npz0-1) then ! k+
         do k=nz+1,nz+LAP
            do j=1,ny
               do i=1,nx
                  k1=k-nz

                  u(i,j,k)=u(i,j,k1)
                  v(i,j,k)=v(i,j,k1)
                  w(i,j,k)=w(i,j,k1)

                  d(i,j,k)=d(i,j,k1)
                  T(i,j,k)=T(i,j,k1)

               enddo
            enddo
         enddo
      endif
   end if
end

!--------------------------------------------------
! Ghost Cell for Jacobian coefficients ; 2nd order exterpolation
subroutine Jac_Ghost_Extent_2nd
   use flow_data
   implicit none
   integer:: i,j,k, i1,j1,k1
   if(npy .eq. 0 .and. Para%Ghost_Cell(3) .eq. 2) then  !i-
      ! do k=1,nz
      !    do j=1,ny
      !       do i=1-LAP, 0
      !          i1= 2-i
      !          Axx(i,j,k)=2.d0*Axx(1,j,k)-Axx(i1,j,k)
      !          Ayy(i,j,k)=2.d0*Ayy(1,j,k)-Ayy(i1,j,k)
      !          Azz(i,j,k)=2.d0*Azz(1,j,k)-Azz(i1,j,k)
      !       enddo
      !    enddo
      ! enddo

      call comput_Jacobian3d_Ghost(3)           ! Comput Jocabian coefficient of the Ghost Cells

   endif

   if(npy .eq. npy0-1 .and.  Para%Ghost_Cell(4) .eq. 2) then  ! i+
      ! do k=1,nz
      !    do j=1,ny
      !       do i=nx+1, nx+LAP
      !          i1=2*nx-i
      !          Axx(i,j,k)=2.d0*Axx(nx,j,k)-Axx(i1,j,k)
      !          Ayy(i,j,k)=2.d0*Ayy(nx,j,k)-Ayy(i1,j,k)
      !          Azz(i,j,k)=2.d0*Azz(nx,j,k)-Azz(i1,j,k)
      !       enddo
      !    enddo
      ! enddo

      call comput_Jacobian3d_Ghost(4)           ! Comput Jocabian coefficient of the Ghost Cells
   endif

   If (Para%Iperiodic_Z == 1 .and. Para%Z_mpi == 0) then
      if(npz==0) then   ! k-
         do k=1-LAP,0
            do j=1,ny
               do i=1,nx
                  k1=nz+k
                  Axx(i,j,k)=Axx(i,j,k1)
                  Ayy(i,j,k)=Ayy(i,j,k1)
                  Azz(i,j,k)=Azz(i,j,k1)

                  Ajac(i,j,k)=Ajac(i,j,k1)
                  Akx(i,j,k)=Akx(i,j,k1)
                  Aky(i,j,k)=Aky(i,j,k1)
                  Akz(i,j,k)=Akz(i,j,k1)
                  Aix(i,j,k)=Aix(i,j,k1)
                  Aiy(i,j,k)=Aiy(i,j,k1)
                  Aiz(i,j,k)=Aiz(i,j,k1)
                  Asx(i,j,k)=Asx(i,j,k1)
                  Asy(i,j,k)=Asy(i,j,k1)
                  Asz(i,j,k)=Asz(i,j,k1)

               enddo
            enddo
         enddo
      end if
      if (npz==npz0-1) then ! k+
         do k=nz+1,nz+LAP
            do j=1,ny
               do i=1,nx
                  k1=k-nz
                  Axx(i,j,k)=Axx(i,j,k1)
                  Ayy(i,j,k)=Ayy(i,j,k1)
                  Azz(i,j,k)=Azz(i,j,k1)

                  Ajac(i,j,k)=Ajac(i,j,k1)
                  Akx(i,j,k)=Akx(i,j,k1)
                  Aky(i,j,k)=Aky(i,j,k1)
                  Akz(i,j,k)=Akz(i,j,k1)
                  Aix(i,j,k)=Aix(i,j,k1)
                  Aiy(i,j,k)=Aiy(i,j,k1)
                  Aiz(i,j,k)=Aiz(i,j,k1)
                  Asx(i,j,k)=Asx(i,j,k1)
                  Asy(i,j,k)=Asy(i,j,k1)
                  Asz(i,j,k)=Asz(i,j,k1)
               enddo
            enddo
         enddo
      endif
   end if

end

!--------Comput Jocabian coefficients at Ghost Cells
subroutine comput_Jacobian3d_Ghost(nb)
   Use flow_data
   implicit none
   integer:: nb, ib,ie,jb,je,kb,ke         ! Jocabian data range
   integer::  ib1, ie1, jb1,je1,kb1,ke1               ! coordinate data range
   real(kind=OCFD_REAL_KIND):: xi1,xj1,xk1, yi1, yj1, yk1 , zi1, zj1, zk1, Jac1
   integer:: i,j,k

   if(nb .eq. 1) then            ! i-
      ib=1-LAP;   ie=0 ;  jb=1;  je=ny ;  kb=1;  ke=nz
      ib1=ib;   ie1=1 ;  jb1=jb;  je1=je ;  kb1=kb;  ke1=ke
   else if (nb .eq. 2) then      ! i+
      ib=nx+1;   ie=nx+LAP ;   jb=1;  je=ny ;  kb=1;  ke=nz
      ib1=nx;    ie1=ie ;     jb1=jb;  je1=je ;  kb1=kb;  ke1=ke
   else if(nb .eq. 3) then
      ib=1;      ie=nx ;      jb=1-LAP;  je=0 ;   kb=1;  ke=nz
      ib1=ib;    ie1=ie ;     jb1=jb;    je1=1 ;  kb1=kb;  ke1=ke
   else if(nb .eq. 4) then
      ib=1;      ie=nx ;      jb=ny+1;   je=ny+LAP ;   kb=1;  ke=nz
      ib1=ib;    ie1=ie ;     jb1=ny;    je1=je ;  kb1=kb;  ke1=ke
   else if(nb .eq. 5) then
      ib=1;      ie=nx ;      jb=1;     je=ny ;   kb=1-LAP;  ke=0
      ib1=ib;    ie1=ie ;     jb1=jb;   je1=je ;  kb1=kb;    ke1=1
   else if(nb .eq. 6) then
      ib=1;      ie=nx ;      jb=1;    je=ny ;     kb=nz+1;  ke=nz+LAP
      ib1=ib;    ie1=ie ;     jb1=jb;    je1=je ;  kb1=nz;  ke1=ke
   endif


   do k=kb,ke
      do j=jb,je
         do i=ib,ie

            if(i .eq. ib1 ) then
               xi1=(-3.d0*Axx(i,j,k)+4.d0*Axx(i+1,j,k)-Axx(i+2,j,k))/(2.d0*hx)     ! 2nd one-side scheme
               yi1=(-3.d0*Ayy(i,j,k)+4.d0*Ayy(i+1,j,k)-Ayy(i+2,j,k))/(2.d0*hx)
               zi1=(-3.d0*Azz(i,j,k)+4.d0*Azz(i+1,j,k)-Azz(i+2,j,k))/(2.d0*hx)
            else if (i .eq. ie1) then
               xi1=(Axx(i-2,j,k)-4.d0*Axx(i-1,j,k)  +3.d0*Axx(i,j,k))/(2.d0*hx)  ! 2nd one-side scheme
               yi1=(Ayy(i-2,j,k)-4.d0*Ayy(i-1,j,k)  +3.d0*Ayy(i,j,k))/(2.d0*hx)  ! 2nd one-side scheme
               zi1=(Azz(i-2,j,k)-4.d0*Azz(i-1,j,k)  +3.d0*Azz(i,j,k))/(2.d0*hx)  ! 2nd one-side scheme
            else if (i .eq. ib1+1  .or. i .eq. ie1-1) then
               xi1=(Axx(i+1,j,k)-Axx(i-1,j,k))/(2.d0*hx)                             ! 2nd centeral scheme
               yi1=(Ayy(i+1,j,k)-Ayy(i-1,j,k))/(2.d0*hx)
               zi1=(Azz(i+1,j,k)-Azz(i-1,j,k))/(2.d0*hx)
            else
               xi1=(8.d0*(Axx(i+1,j,k)-Axx(i-1,j,k)) - (Axx(i+2,j,k)-Axx(i-2,j,k)))/(12.d0*hx)   ! 4th central scheme
               yi1=(8.d0*(Ayy(i+1,j,k)-Ayy(i-1,j,k)) - (Ayy(i+2,j,k)-Ayy(i-2,j,k)))/(12.d0*hx)
               zi1=(8.d0*(Azz(i+1,j,k)-Azz(i-1,j,k)) - (Azz(i+2,j,k)-Azz(i-2,j,k)))/(12.d0*hx)
            endif

            if(j .eq. jb1 ) then
               xj1=(-3.d0*Axx(i,j,k)+4.d0*Axx(i,j+1,k)-Axx(i,j+2,k))/(2.d0*hy)     ! 2nd one-side scheme
               yj1=(-3.d0*Ayy(i,j,k)+4.d0*Ayy(i,j+1,k)-Ayy(i,j+2,k))/(2.d0*hy)
               zj1=(-3.d0*Azz(i,j,k)+4.d0*Azz(i,j+1,k)-Azz(i,j+2,k))/(2.d0*hy)

            else if (j .eq. je1)	 then
               xj1=(Axx(i,j-2,k)-4.d0*Axx(i,j-1,k)  +3.d0*Axx(i,j,k))/(2.d0*hy)  ! 2nd one-side scheme
               yj1=(Ayy(i,j-2,k)-4.d0*Ayy(i,j-1,k)  +3.d0*Ayy(i,j,k))/(2.d0*hy)
               zj1=(Azz(i,j-2,k)-4.d0*Azz(i,j-1,k)  +3.d0*Azz(i,j,k))/(2.d0*hy)
            else if (j .eq. jb1+1  .or. j .eq. je1-1)  then
               xj1=(Axx(i,j+1,k)-Axx(i,j-1,k))/(2.d0*hy)     		 ! 2nd centeral scheme
               yj1=(Ayy(i,j+1,k)-Ayy(i,j-1,k))/(2.d0*hy)
               zj1=(Azz(i,j+1,k)-Azz(i,j-1,k))/(2.d0*hy)
            else
               xj1=(8.d0*(Axx(i,j+1,k)-Axx(i,j-1,k)) - (Axx(i,j+2,k)-Axx(i,j-2,k)))/(12.d0*hy)   ! 4th central scheme
               yj1=(8.d0*(Ayy(i,j+1,k)-Ayy(i,j-1,k)) - (Ayy(i,j+2,k)-Ayy(i,j-2,k)))/(12.d0*hy)
               zj1=(8.d0*(Azz(i,j+1,k)-Azz(i,j-1,k)) - (Azz(i,j+2,k)-Azz(i,j-2,k)))/(12.d0*hy)
            endif

            if(k .eq. kb1 ) then
               xk1=(-3.d0*Axx(i,j,k)+4.d0*Axx(i,j,k+1)-Axx(i,j,k+2))/(2.d0*hz)     ! 2nd one-side scheme
               yk1=(-3.d0*Ayy(i,j,k)+4.d0*Ayy(i,j,k+1)-Ayy(i,j,k+2))/(2.d0*hz)
               zk1=(-3.d0*Azz(i,j,k)+4.d0*Azz(i,j,k+1)-Azz(i,j,k+2))/(2.d0*hz)
            else if (k .eq. ke1)	 then
               xk1=(Axx(i,j,k-2)-4.d0*Axx(i,j,k-1)  +3.d0*Axx(i,j,k))/(2.d0*hz)  ! 2nd one-side scheme
               yk1=(Ayy(i,j,k-2)-4.d0*Ayy(i,j,k-1)  +3.d0*Ayy(i,j,k))/(2.d0*hz)
               zk1=(Azz(i,j,k-2)-4.d0*Azz(i,j,k-1)  +3.d0*Azz(i,j,k))/(2.d0*hz)
            else if (k .eq. kb1+1  .or. k .eq. ke1-1)  then
               xk1=(Axx(i,j,k+1)-Axx(i,j,k-1))/(2.d0*hz)                             ! 2nd centeral scheme
               yk1=(Ayy(i,j,k+1)-Ayy(i,j,k-1))/(2.d0*hz)
               zk1=(Azz(i,j,k+1)-Azz(i,j,k-1))/(2.d0*hz)
            else
               xk1=(8.d0*(Axx(i,j,k+1)-Axx(i,j,k-1)) - (Axx(i,j,k+2)-Axx(i,j,k-2)))/(12.d0*hz)   ! 4th central scheme
               yk1=(8.d0*(Ayy(i,j,k+1)-Ayy(i,j,k-1)) - (Ayy(i,j,k+2)-Ayy(i,j,k-2)))/(12.d0*hz)
               zk1=(8.d0*(Azz(i,j,k+1)-Azz(i,j,k-1)) - (Azz(i,j,k+2)-Azz(i,j,k-2)))/(12.d0*hz)
            endif

            Jac1=1.d0/(xi1*yj1*zk1+yi1*zj1*xk1+zi1*xj1*yk1-zi1*yj1*xk1-yi1*xj1*zk1-xi1*zj1*yk1)   ! 1./Jocabian = d(x,y,z)/d(i,j,k)
            Ajac(i,j,k)=Jac1
            Akx(i,j,k)=Jac1*(yj1*zk1-zj1*yk1)
            Aky(i,j,k)=Jac1*(zj1*xk1-xj1*zk1)
            Akz(i,j,k)=Jac1*(xj1*yk1-yj1*xk1)
            Aix(i,j,k)=Jac1*(yk1*zi1-zk1*yi1)
            Aiy(i,j,k)=Jac1*(zk1*xi1-xk1*zi1)
            Aiz(i,j,k)=Jac1*(xk1*yi1-yk1*xi1)
            Asx(i,j,k)=Jac1*(yi1*zj1-zi1*yj1)
            Asy(i,j,k)=Jac1*(zi1*xj1-xi1*zj1)
            Asz(i,j,k)=Jac1*(xi1*yj1-yi1*xj1)
            if(Jac1 .lt. 0) then
               print*, "      comput_Jacobian3d_Ghost, Jocabian < 0 !!! , Jac=", Jac1
               print*, "      i,j,k=", i_offset(npx)+i-1, j_offset(npy)+j-1, k_offset(npz)+k-1
            endif
         enddo
      enddo
   enddo

end


