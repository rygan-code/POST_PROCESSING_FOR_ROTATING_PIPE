!-----OpenCFD-SC version 2 --------------------------------------------------------------------
! Copyright by Li Xinliang, LHD, Institute of Mechanics, CAS, lixl@imech.ac.cn
! Codey by Li Xinliang, 2021-2
!----------------------------------------------------------------------------------------------

! Double precsion (real*8)  or  Single Precision (real*4)       (default: double precision)
module OCFD_precision
  use mpi
  implicit none
!------For Doubleprecision  (real*8)--------------------------------------------------------------------------
  integer,parameter::OCFD_REAL_KIND=8,  OCFD_DATA_TYPE=MPI_DOUBLE_PRECISION   ! double precison computing
! ------For Single precision (real*4)-------------------------------------------------------------------------
!     integer,parameter::OCFD_REAL_KIND=4,  OCFD_DATA_TYPE=MPI_REAL             !  single precision computing
end module OCFD_precision


 !------parameters used in OpenCFD---------------------------------------
module OCFD_constants
  Use OCFD_precision
  implicit none
  integer,parameter:: Nvars=5           ! 5 conservative variables, 5 Equations
  integer,parameter:: OCFD_Turb_None=0                                              ! Turbulence model

  integer,parameter::  OCFD_Scheme_WENO5=1,OCFD_Scheme_WENO7=2,  &                  ! Schemes
    OCFD_Scheme_OMP6=3, OCFD_Scheme_UD7L=4,   &
    OCFD_Scheme_CD6=5,OCFD_Scheme_CD8=6,   &
    OCFD_Scheme_Hybrid=10,  &
    OCFD_Scheme_USER=99

  integer,parameter:: OCFD_Split_SW=1, OCFD_Split_LLF=2

  integer,parameter:: BC_None=0,BC_Blunt2d=1,  BC_BoundaryLayer=2, BC_SweptCorner=3, BC_User_Def=99

  integer,parameter:: GRID1D=10, GRID2D_PLANE=20, GRID2D_AXIAL_SYMM=21, GRID3D=30, GRID_AND_JACOBIAN3D=31

  integer,parameter:: DEL_LIFT=1, DEL_RIGHT=2, DEL_NONE=0             !   WENO5-type boundary Scheme
  integer,parameter:: OCFD_ANA_USER=99, OCFD_ANA_time_average=100, &
    OCFD_ANA_Q=101, OCFD_ANA_BOX=102, OCFD_ANA_SAVEDATA=103, OCFD_ANA_Corner=104
end

!----mpi parameters------------------------------------------------
module Para_mpi
  implicit none
  integer:: my_id, npx,npy,npz, npx0,npy0,npz0, nprocs, &         ! npx0, npy0 : zone numbers in x- and y- direction
    MPI_COMM_X,MPI_COMM_Y,MPI_COMM_Z,   MPI_COMM_XY,MPI_COMM_YZ,  MPI_COMM_XZ, &
    ID_XM1,ID_XP1,ID_YM1,ID_YP1, ID_ZM1, ID_ZP1                  ! ID_XM1:  MPI ID for zone-1 (lift one) in x-direction
  integer:: LAP                                                    ! Overlap length for Block-Block
  integer,allocatable,dimension(:):: i_offset,j_offset, k_offset, i_nn,j_nn, k_nn
  integer:: TYPE_LAPX1,TYPE_LAPY1,TYPE_LAPZ1,TYPE_LAPX2,TYPE_LAPY2,TYPE_LAPZ2           ! user defined MPI DATA TYPE (for Send & Recv)
end


!-----------parameters for numerical schemes------------------------
module Scheme_Para
  use OCFD_precision
  implicit none
  TYPE TYPE_Scheme
    integer:: Scheme_Invis, Scheme_Vis
    integer:: Bound_index(2,3), Scheme_boundary(6)    !Bound_index(:,:)  if using boundary scheme(0/1)
    !Scheme_boundary(:) : boundary scheme type ( i-, i+, j-, j+, k-, k+:  0 : default, WENO5 Del-substencil type; -1 full-WENO5 type)
    integer:: Ka1,Kb1,Ka2,Kb2                      ! [Ka1,Kb1]  stencil for flux+ ;   [Ka2,Kb2] stencil for flux-
    integer:: Ka1_H1,Kb1_H1,Ka2_H1,Kb2_H1,Ka1_H2,Kb1_H2,Ka2_H2,Kb2_H2,Ka1_H3,Kb1_H3,Ka2_H3,Kb2_H3       ! [Ka,Kb] for Hybrid scheme (scheme1, scheme2, scheme3)
    Real(kind=OCFD_REAL_KIND):: UD7L_Diss, UD7L(8) , Hybrid_para(100)
    ! UD7L_Diss:   dissipation coefficient (0-1;  0 : CD8,   1: UD7)
    ! UD7L coefficients for 7th low-dissipation upwind difference sheme
  end TYPE
  TYPE (TYPE_Scheme):: Scheme

end


!-----------flow parameters ----------------------------------------------
module flow_para
  use OCFD_constants
  use Para_mpi
  use Scheme_para
  implicit none


  TYPE TYPE_Para         ! input parameters
    integer::   Iperiodic_X,Iperiodic_Y,Iperiodic_Z,Z_mpi          ! 1/0 : If periodic or not
    integer::   IF_Scheme_Character , IF_Viscous ,Istep_Show, Istep_Save, Flux_Splitting,&
      Iflag_Gridtype ,  &                          ! 0: Grid only, 1: Grid and Jocabian
      IBC , IF_Mass_Force, & ! boundary conditon ;  IF_Mase_Force :  0 (default): no mass force;   1 with mass force
      ANA_Number, &          ! number of analysis process
      Nfiltering, &          ! Number of filtering
      Ghost_cell(6)          ! 0 : non-GhostCell (default);  1 Ghost Cell (Auto, exterpolation);  2 Ghost cell (user defined)
    Real(kind=OCFD_REAL_KIND):: Periodic_ISpan(3), Periodic_JSpan(3) , Periodic_KSpan(3)		    ! Span in peridoic direction
    Real(kind=OCFD_REAL_KIND):: omega, Ma, gamma
    integer :: KRK
  end TYPE

  TYPE (TYPE_Para):: Para

  integer:: nx_global, ny_global,nz_global,  nx, ny, nz, Istep, iistep
  Real(kind=OCFD_REAL_KIND):: Cp, Cv, tt, hx, hy, hz
end

!---------flow data---------------------------------------------------------
module flow_data
  use flow_para
  implicit none

  real(kind=OCFD_REAL_KIND),allocatable,dimension(:,:,:):: Axx,Ayy,Azz,  &   ! Coordinates
    Akx,Aky,Akz,Aix,Aiy,Aiz,Asx,Asy,Asz, Ajac , &       !  Jacobian coefficients
    Akx1,Aky1,Akz1,Aix1,Aiy1,Aiz1,Asx1,Asy1,Asz1          ! Akx1=Akx/Ajac
  real(kind=OCFD_REAL_KIND),allocatable:: f(:,:,:,:),fn(:,:,:,:), du(:,:,:,:), Amu(:,:,:),Amu_t(:,:,:)
  real(kind=OCFD_REAL_KIND),allocatable,dimension(:,:,:):: d,u,v,w,T
  real(kind=OCFD_REAL_KIND),allocatable,dimension(:,:,:):: vor_x,vor_y,vor_z,Q_criterion
  real(kind=OCFD_REAL_KIND),allocatable,dimension(:,:,:):: Rhybrid         ! Index of hybrid scheme
  Real(kind=OCFD_REAL_KIND):: Re, Ma, gamma, Pr , Ref_Amu_T0 ,dt, End_time
end

module fourier_data
  use flow_para
  implicit none

  real(kind=OCFD_REAL_KIND),allocatable,dimension(:,:,:):: d_f,u_f,v_f,w_f,T_f
  real(kind=OCFD_REAL_KIND),allocatable,dimension(:,:,:):: vor_xf,vor_yf,vor_zf,Q_criterion_f
end

!------------Set parameters --------------------------------
subroutine set_parameters
  use flow_para
  implicit none

  hx=1.d0/(nx_global-1.d0)
  hy=1.d0/(ny_global-1.d0)
  hz=1.d0/(nz_global-1.d0)

  Cv=1.d0/(Para%gamma*(Para%gamma-1.d0)*Para%Ma*Para%Ma)
  Cp=Cv*Para%gamma

!---------Set Scheme%bound_index       (指示6个边界 是否采用边界格式）
! Scheme%Scheme_boundary(:)==-1    ! Ghost-Cell type boundary   (Do not use boundary scheme)
  Scheme%bound_index(:,:)=0             ! default :  not use boundary scheme

  if(npx .eq. 0 .and. Para%Iperiodic_X .eq. 0      )  Scheme%bound_index(1,1)= 1    !  i-
  if(npx .eq. npx0-1 .and. Para%Iperiodic_X .eq.0  )  Scheme%bound_index(2,1)= 1	 !  i+

  if(npy .eq. 0 .and. Para%Iperiodic_Y .eq. 0      )  Scheme%bound_index(1,2)= 1      !  j-
  if(npy .eq. npy0-1 .and. Para%Iperiodic_Y .eq. 0 )  Scheme%bound_index(2,2)= 1	     !  j+

  if(npz .eq. 0 .and. Para%Iperiodic_Z .eq. 0      )  Scheme%bound_index(1,3)= 1      !  k-
  if(npz .eq. npz0-1 .and. Para%Iperiodic_Z .eq. 0 )  Scheme%bound_index(2,3)= 1	     !  k+





! ----[Ka1,Kb1]:  Stencil of positive flux F(i+1/2) : [i+Ka1, i+Kb1] ;
!     [Ka2,Kb2]:  Stencil of negative flux ;

  select case (Scheme%Scheme_Invis )
   case (OCFD_Scheme_WENO5)
    Scheme%Ka1=-2 ;  Scheme%Kb1=2      ! Stencil for F(i+1/2) of WENO5+ : [i-2, ..., i+2]
    Scheme%Ka2=-1 ;  Scheme%Kb2=3      ! Stencil for F(i+1/2) of WENO5- : [i-1, ..., i+3]
   case( OCFD_Scheme_WENO7)
    Scheme%Ka1=-3; Scheme%Kb1=3
    Scheme%Ka2=-2; Scheme%Kb2=4
   case( OCFD_Scheme_OMP6)
    Scheme%Ka1=-3; Scheme%Kb1=4
    Scheme%Ka2=-3; Scheme%Kb2=4
   case( OCFD_Scheme_UD7L) 		   ! low-dissipative Upwind Difference scheme
    Scheme%Ka1=-3; Scheme%Kb1=4
    Scheme%Ka2=-3; Scheme%Kb2=4
   case( OCFD_Scheme_Hybrid) 		   !  Hybrid scheme
    Scheme%Ka1=-3; Scheme%Kb1=4
    Scheme%Ka2=-3; Scheme%Kb2=4

   case default
    print*, "The Inviscous Scheme is not supported"
    stop
  end select

!       [Ka,Kb] for hybrid scheme
  Scheme%Ka1_H1=-3; Scheme%Kb1_H1=4; Scheme%Ka2_H1=-3; Scheme%Kb2_H1=4      ! UD7L
  Scheme%Ka1_H2=-3; Scheme%Kb1_H2=3; Scheme%Ka2_H2=-2; Scheme%Kb2_H2=4      ! WENO7
  Scheme%Ka1_H3=-2; Scheme%Kb1_H3=2; Scheme%Ka2_H3=-1; Scheme%Kb2_H3=3      ! WENO5

  call set_scheme_para

end

subroutine set_scheme_para
  use flow_para
  implicit none
  real(kind=OCFD_REAL_KIND):: UD7(8),CD8(8)

  ! Use 6th order central if UD7L_diss = -1
  if (Scheme%UD7L_Diss == -1.d0) then
    Scheme%UD7L = (/ 0.d0, 1.d0, -8.d0, 37.d0, 37.d0, -8.d0, 1.d0, 0.d0 /)/60.d0
  else
    UD7=(/ -3.d0, 25.d0, -101.d0, 319.d0,  214.d0, -38.d0, 4.d0, 0.d0  /)       ! coefficients for UD7
    CD8=(/ -3.d0,  29.d0, -139.d0,   533.d0,   533.d0, -139.d0,  29.d0, -3.d0 /)        ! coefficients for CD8
    Scheme%UD7L=Scheme%UD7L_Diss*UD7/420.d0 + (1.d0-Scheme%UD7L_Diss)*CD8/840.d0
  end if
end

subroutine comput_Jacobian3d
   Use flow_data
   implicit none
   real(kind=OCFD_REAL_KIND), allocatable,dimension(:,:,:):: xi,xj,xk, yi, yj, yk, zi, zj, zk
   real(kind=OCFD_REAL_KIND):: xi1,xj1,xk1, yi1, yj1, yk1 , zi1, zj1, zk1, Jac1
   integer:: i,j,k

   allocate (xi(nx,ny,nz),xj(nx,ny,nz),xk(nx,ny,nz),   &
      yi(nx,ny,nz),yj(nx,ny,nz),yk(nx,ny,nz),   &
      zi(nx,ny,nz),zj(nx,ny,nz),zk(nx,ny,nz))

   call OCFD_dx0(Axx(1:nx,1:ny,1:nz),xi,Scheme%Scheme_Vis)
   call OCFD_dx0(Ayy(1:nx,1:ny,1:nz),yi,Scheme%Scheme_Vis)
   call OCFD_dx0(Azz(1:nx,1:ny,1:nz),zi,Scheme%Scheme_Vis)
   call OCFD_dy0(Axx(1:nx,1:ny,1:nz),xj,Scheme%Scheme_Vis)
   call OCFD_dy0(Ayy(1:nx,1:ny,1:nz),yj,Scheme%Scheme_Vis)
   call OCFD_dy0(Azz(1:nx,1:ny,1:nz),zj,Scheme%Scheme_Vis)
   call OCFD_dz0(Axx(1:nx,1:ny,1:nz),xk,Scheme%Scheme_Vis)
   call OCFD_dz0(Ayy(1:nx,1:ny,1:nz),yk,Scheme%Scheme_Vis)
   call OCFD_dz0(Azz(1:nx,1:ny,1:nz),zk,Scheme%Scheme_Vis)


   do k=1,nz
      do j=1,ny
         do i=1,nx
            xi1=xi(i,j,k); xj1=xj(i,j,k); xk1=xk(i,j,k)
            yi1=yi(i,j,k); yj1=yj(i,j,k); yk1=yk(i,j,k)
            zi1=zi(i,j,k); zj1=zj(i,j,k); zk1=zk(i,j,k)
            Jac1=1.d0/(xi1*yj1*zk1+yi1*zj1*xk1+zi1*xj1*yk1-zi1*yj1*xk1-yi1*xj1*zk1-xi1*zj1*yk1)   ! 1./Jocabian = d(x,y,z)/d(i,j,k)
            Ajac(i,j,k)=Jac1
            ! unit vector in xyz direction
            Akx(i,j,k)=Jac1*(yj1*zk1-zj1*yk1)   
            Aky(i,j,k)=Jac1*(zj1*xk1-xj1*zk1)
            Akz(i,j,k)=Jac1*(xj1*yk1-yj1*xk1)
            Aix(i,j,k)=Jac1*(yk1*zi1-zk1*yi1)
            Aiy(i,j,k)=Jac1*(zk1*xi1-xk1*zi1)
            Aiz(i,j,k)=Jac1*(xk1*yi1-yk1*xi1)
            Asx(i,j,k)=Jac1*(yi1*zj1-zi1*yj1)
            Asy(i,j,k)=Jac1*(zi1*xj1-xi1*zj1)
            Asz(i,j,k)=Jac1*(xi1*yj1-yi1*xj1)
            ! print*,i,j,k,xi1,xj1,xk1,yi1,yj1,yk1,zi1,zj1,zk1
            ! if(Jac1 .lt. 0) then
            !    print*, " Jocabian < 0 !!! , Jac=", Jac1
            !    print*, "i,j,k=", i_offset(npx)+i-1, j_offset(npy)+j-1, k_offset(npz)+k-1
            ! endif
         enddo
      enddo
   enddo
   deallocate ( xi,xj,xk, yi, yj, yk, zi, zj, zk)
end

subroutine set_default_parameters
    use flow_para
    use flow_data
  !--------Set defalut parameter ----------------
    gamma=1.4d0; Para%gamma=gamma
    Re=1000.d0
    Ma=0.3d0; Para%Ma=Ma
    Pr=0.7d0
    npx0=1
    npy0=1
    npz0=1
    LAP=0
    Scheme%Scheme_Invis = OCFD_Scheme_WENO5
    Scheme%Scheme_Vis = OCFD_Scheme_CD8
    Iperiodic_X=0
    Iperiodic_Y=0
    Iperiodic_Z=0
    IF_Scheme_Character=0
    Ref_Amu_T0=288.15d0
    dt=1.d-4
    End_time=100.d0
    Iflag_Gridtype=GRID3D
    IF_Viscous=1

  end subroutine