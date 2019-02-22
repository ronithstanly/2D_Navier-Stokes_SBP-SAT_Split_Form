  !**************************************************************************
  !**************************************************************************
  !!!!**************** SBP - SAT - Split Fluxes *************************!!!!
  !**************************************************************************
  !**************************************************************************
  ! This code solves the 2-D Compressible Navier-Stokes equations using 
  ! Summation by Parts (SBP) operators, Simultaneous Approximation Terms (SAT) 
  ! boundary conditions and Split-form fluxes. Second, fourth and sixth-order 
  ! accurate SBP operators are avialable. The default test case selected is the 
  ! Lid Driven cavity at Re=100
  ! Author: Ronith Stanly
  ! Lecturer: Prof.Steven Frankel
  ! Last Updated: 22 February, 2019
  ! CFD Lab, Technion, Israel
  !**************************************************************************
  !**************************************************************************

  ! This module defines the KIND types of all the variables used in the code: 
  ! I4B, I2B and I1B for integer variables, SP and DP for real variables (and
  ! SPC and DPC for corresponding complex cases), and LGT for the default 
  ! logical type. This follows the convention used the Numerical Recipes for 
  ! Fortran 90 types module 'nrtype', pp. 1361
  MODULE types_vars
    ! Symbolic names for kind types of 4-, 2- and 1-byte integers:   
    INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
    INTEGER, PARAMETER :: I2B = SELECTED_INT_KIND(4)
    INTEGER, PARAMETER :: I1B = SELECTED_INT_KIND(2)
    ! Symbolic names for kind types of single- and double-precison reals
    INTEGER, PARAMETER :: SP = KIND(1.0)
    INTEGER, PARAMETER :: DP = KIND(1.0D0)
    ! Symbolic names for kind types of single- and double-precison complex
    INTEGER, PARAMETER :: SPC = KIND((1.0,1.0))
    INTEGER, PARAMETER :: DPC = KIND((1.0D0,1.0D0))
    ! Symbolic name for kind type of default logical
    INTEGER, PARAMETER :: LOGIC = KIND(.true.)
    ! Frequently used mathematical constants (with precision to spare)
    REAL(DP), PARAMETER :: zero  = 0.0_dp
    REAL(DP), PARAMETER :: half  = 0.5_dp
    REAL(DP), PARAMETER :: one   = 1.0_dp
    REAL(DP), PARAMETER :: two   = 2.0_dp
    REAL(DP), PARAMETER :: three = 3.0_dp
    REAL(DP), PARAMETER :: four  = 4.0_dp
    REAL(DP), PARAMETER :: pi    = 3.141592653589793238462643383279502884197_dp
    REAL(DP), PARAMETER :: pio2  = 1.57079632679489661923132169163975144209858_dp
    REAL(DP), PARAMETER :: twopi = 6.283185307179586476925286766559005768394_dp
  END MODULE types_vars

  !**************************************************************************
  !**************************************************************************
  !**************************************************************************
  !**************************************************************************
  ! This module defines and allocates the variables needed in the current simulation
  ! If needed, add new variables at the beginning of the module, then allocate 
  ! them in the subroutine memalloc
  MODULE variables
    USE types_vars
    ! Add new variables here
    INTEGER :: nel, ntimes, nptsx, nptsy, nptst, iwave, tick, nvar, ic_user, porder
    REAL(DP) :: a, b, Dx, Dy, t, cfl, cfl_ip, u_rk3, u_rk4, l2_rk3
    REAL(DP) :: c, d, Dt, gam, rgas
    REAL(DP) :: ay, by, alpha_split, beta_split

    !1D arrays
    REAL(DP), ALLOCATABLE, DIMENSION(:) :: x, y, time, h, l2_rk4
    REAL(DP), ALLOCATABLE, DIMENSION(:) :: e_x_0, e_x_n, e_y_0, e_y_n

    !2D arrays
    REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: Hmatinv, Qmat, Mmat, Smat, Bmat, D1mat, D2mat, Dx1, Q_tran, Hinvx, Hinvy
    REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: Dy1, g_x_0, g_x_n, ke, temp, sos, mach, rho, uvel, vvel, pre, e, et
    REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: g_y_0, g_y_n, correction

    REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: rho_u, f2, df2dx, f6, df6dx, f7, df7dx
    REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: rho_v, g2, dg2dy, g6, dg6dy, g7, dg7dy

    !3D arrays
    REAL(DP), ALLOCATABLE, DIMENSION(:,:,:) :: f, u, uold, dfdx, k1, k2, k3, k4, u_ic, u_dd, rhs, u_e, dfdx_e
    REAL(DP), ALLOCATABLE, DIMENSION(:,:,:) :: g, dgdy, temp1d, SAT_x0, SAT_xn, dgdy_e, SAT_y0, SAT_yn 
    REAL(DP), ALLOCATABLE, DIMENSION(:,:,:) :: f_visc, g_visc, dfdx_visc, dgdy_visc

    REAL(DP), ALLOCATABLE, DIMENSION(:,:,:) :: phi, f1, df1dx, split_x1, split_x2, split_x3, split_x4
    REAL(DP), ALLOCATABLE, DIMENSION(:,:,:) :: split_x5, split_x6, split_x7, f3, df3dx, df4dx, f5, df5dx
    REAL(DP), ALLOCATABLE, DIMENSION(:,:,:) :: pre_mat_x, dpdx,f4
    REAL(DP), ALLOCATABLE, DIMENSION(:,:,:) :: g1, dg1dy, split_y1, split_y2, split_y3, split_y4
    REAL(DP), ALLOCATABLE, DIMENSION(:,:,:) :: split_y5, split_y6, split_y7, g3, dg3dy, dg4dy, g5, dg5dy
    REAL(DP), ALLOCATABLE, DIMENSION(:,:,:) :: pre_mat_y, dpdy, g4
    CONTAINS


!************* Allocate memory*******************
    SUBROUTINE memalloc
      ! Allocate memory for grid, solution, and flux function

      !1D arrays
      ALLOCATE(h(1:2),l2_rk4(1:2))
      ALLOCATE(x(0:nptsx),e_x_0(0:nptsx),e_x_n(0:nptsx),e_y_0(0:nptsy),e_y_n(0:nptsy))

      ALLOCATE(y(0:nptsy))

      !2D arrays
      ALLOCATE(Dx1(0:nptsx,0:nptsx),Hinvx(0:nptsx,0:nptsx),Hinvy(0:nptsy,0:nptsy))

      ALLOCATE(rho(0:nptsx,0:nptsy),uvel(0:nptsx,0:nptsy),vvel(0:nptsx,0:nptsy),pre(0:nptsx,0:nptsy))
      ALLOCATE(e(0:nptsx,0:nptsy),et(0:nptsx,0:nptsy),Dy1(0:nptsy,0:nptsy)) 
      ALLOCATE(g_x_0(1:nvar,0:nptsy),g_x_n(1:nvar,0:nptsy))
      ALLOCATE(ke(0:nptsx,0:nptsy),temp(0:nptsx,0:nptsy),sos(0:nptsx,0:nptsy),mach(0:nptsx,0:nptsy))
      ALLOCATE(g_y_0(1:nvar,0:nptsx),g_y_n(1:nvar,0:nptsx))

      ALLOCATE(rho_u(0:nptsx,0:nptsy),f2(0:nptsx,0:nptsy))
      ALLOCATE(df2dx(0:nptsx,0:nptsy),f6(0:nptsx,0:nptsy),df6dx(0:nptsx,0:nptsy),f7(0:nptsx,0:nptsy))
      ALLOCATE(df7dx(0:nptsx,0:nptsy))
      ALLOCATE(rho_v(0:nptsx,0:nptsy),g2(0:nptsx,0:nptsy),dg2dy(0:nptsx,0:nptsy),g6(0:nptsx,0:nptsy))
      ALLOCATE(dg6dy(0:nptsx,0:nptsy),g7(0:nptsx,0:nptsy),dg7dy(0:nptsx,0:nptsy))

      !3D arrays
      ALLOCATE(f(1:nvar,0:nptsx,0:nptsy),u(1:nvar,0:nptsx,0:nptsy),u_ic(1:nvar,0:nptsx,0:nptsy))
      ALLOCATE(u_dd(1:nvar,0:nptsx,0:nptsy),uold(1:nvar,0:nptsx,0:nptsy),dfdx(1:nvar,0:nptsx,0:nptsy))
      ALLOCATE(k1(1:nvar,0:nptsx,0:nptsy),k2(1:nvar,0:nptsx,0:nptsy),k3(1:nvar,0:nptsx,0:nptsy))
      ALLOCATE(k4(1:nvar,0:nptsx,0:nptsy),g(1:nvar,0:nptsx,0:nptsy),rhs(1:nvar,0:nptsx,0:nptsy))
      ALLOCATE(u_e(1:nvar,0:nptsx,0:nptsy),dfdx_e(1:nvar,0:nptsx,0:nptsy),dgdy(1:nvar,0:nptsx,0:nptsy))
      ALLOCATE(SAT_x0(1:nvar,0:nptsx,0:nptsy),SAT_xn(1:nvar,0:nptsx,0:nptsy),dgdy_e(1:nvar,0:nptsx,0:nptsy))
      ALLOCATE(SAT_y0(1:nvar,0:nptsx,0:nptsy),SAT_yn(1:nvar,0:nptsx,0:nptsy))
      ALLOCATE(f_visc(1:nvar,0:nptsx,0:nptsy),g_visc(1:nvar,0:nptsx,0:nptsy))
      ALLOCATE(dfdx_visc(1:nvar,0:nptsx,0:nptsy),dgdy_visc(1:nvar,0:nptsx,0:nptsy))

      ALLOCATE(phi(1:nvar,0:nptsx,0:nptsy),f1(1:nvar,0:nptsx,0:nptsy),df1dx(1:nvar,0:nptsx,0:nptsy))
      ALLOCATE(split_x1(1:nvar,0:nptsx,0:nptsy),split_x2(1:nvar,0:nptsx,0:nptsy),split_x3(1:nvar,0:nptsx,0:nptsy))
      ALLOCATE(split_x4(1:nvar,0:nptsx,0:nptsy),split_x5(1:nvar,0:nptsx,0:nptsy),split_x6(1:nvar,0:nptsx,0:nptsy))
      ALLOCATE(split_x7(1:nvar,0:nptsx,0:nptsy),f3(1:nvar,0:nptsx,0:nptsy),df3dx(1:nvar,0:nptsx,0:nptsy))  
      ALLOCATE(f4(1:nvar,0:nptsx,0:nptsy),df4dx(1:nvar,0:nptsx,0:nptsy),f5(1:nvar,0:nptsx,0:nptsy)) 
      ALLOCATE(df5dx(1:nvar,0:nptsx,0:nptsy),pre_mat_x(1:nvar,0:nptsx,0:nptsy),dpdx(1:nvar,0:nptsx,0:nptsy))
      ALLOCATE(g1(1:nvar,0:nptsx,0:nptsy),dg1dy(1:nvar,0:nptsx,0:nptsy),split_y1(1:nvar,0:nptsx,0:nptsy))
      ALLOCATE(split_y2(1:nvar,0:nptsx,0:nptsy),split_y3(1:nvar,0:nptsx,0:nptsy),split_y4(1:nvar,0:nptsx,0:nptsy))
      ALLOCATE(split_y5(1:nvar,0:nptsx,0:nptsy),split_y6(1:nvar,0:nptsx,0:nptsy),split_y7(1:nvar,0:nptsx,0:nptsy))
      ALLOCATE(g3(1:nvar,0:nptsx,0:nptsy),dg3dy(1:nvar,0:nptsx,0:nptsy),g4(1:nvar,0:nptsx,0:nptsy))
      ALLOCATE(dg4dy(1:nvar,0:nptsx,0:nptsy),g5(1:nvar,0:nptsx,0:nptsy),dg5dy(1:nvar,0:nptsx,0:nptsy))
      ALLOCATE(pre_mat_y(1:nvar,0:nptsx,0:nptsy),dpdy(1:nvar,0:nptsx,0:nptsy))

    END SUBROUTINE memalloc


!************* Deallocate memory (end of the program)*******************
    SUBROUTINE dealloc
      ! Deallocate memory for grid, solution, and flux function
      DEALLOCATE(x,u,f,uold,dfdx,k1,k2,k3,k4,u_dd,u_e)
    END SUBROUTINE dealloc

  END MODULE variables

  !**************************************************************************
  !**************************************************************************
  MODULE subroutines
    USE types_vars
    USE variables
    CONTAINS

!************* Get the inputs *******************
   SUBROUTINE inputs
    IMPLICIT NONE  ! Forces explicit type declaration to avoid errors

    ! Assume 2D domain is x[-5,5],y[-5,5] and time
    a = 0.0d0; b = 1.0d0; ay=0.0d0; by=1.0d0; c = 0.0d0; d = 6.0d0
    ! Kennedy & Gruber Split parameters
    alpha_split=0.25d0; beta_split=0.25d0

    ! Number of variables (or no. of eqs. in the system), here 2 (for 2D Navier-Stokes it is 4)
    nvar=4
    gam=1.4d0 ! gamma
  !  rgas=8.314d0 ! Gas constant from ideal gas law
    rgas=8.314598/28.97d0 

    ! Read from screen input information for program
    WRITE(*,*) 'Assuming equal number of cewlls in both directions'
    WRITE(*,*) 'Please input the number of cells/volumes/elements:'
    READ(*,*) nel
    nptsx = nel
    nptsy = nel

    WRITE(*,*) 'Please input the desired CFL number'
    READ(*,*) cfl
    cfl_ip=cfl
    Dx = (b-a)/FLOAT(nptsx)
    Dy= (by-ay)/FLOAT(nptsy)
    !Dt = (cfl_ip*Dx)
    !Dt = (cfl_ip * Dx)/a
    !nptst = ABS(d-c)/Dt

    WRITE(*,*) 'Enter preferred order of accuracy 2,4 or 6?'
    READ(*,*) porder

    WRITE(*,*) 'Time-step size=', Dt
    WRITE(*,*) 'Domain will be discretized with ', nptsx**2 
    WRITE(*,*) 'for the given CFL number=', cfl_ip

    ! Echo print your input to make sure it is correct
    WRITE(*,*) 'Your 2D domain is from, x=', a, ' to ', b, 'y=',ay,' to ',by
    WRITE(*,*) 'and time is from ',c,'s to ', d,'s'

  END SUBROUTINE inputs

!*********** Generate a 2D grid ***************
    SUBROUTINE grid2d
    IMPLICIT NONE
    INTEGER :: i,j

    ! Generate grid
    h(1)=Dx
    h(2)=Dy
    DO i = 0, nptsx
      x(i) = a + i*Dx
    END DO 

    DO j=0, nptsy
      y(j) = ay + j*Dy
    END DO

    END SUBROUTINE grid2d


!************** Provide intial condition *************************
!***************** Lid-Driven Cavity *****************************
    SUBROUTINE init2d
      IMPLICIT NONE
      INTEGER :: i, j
      REAL, DIMENSION(:,:) :: rad(0:nptsx,0:nptsy) 

      DO i=0, nptsx
        DO j=0, nptsy
          rho(i,j)=1.0d0

          uvel(i,j)=1.0d0
          vvel(i,j)=0.0d0
          
          pre(i,j)=1.0d0

          sos=SQRT(gam*pre(i,j)/rho(i,j))

        END DO
      END DO

    Dt=cfl_ip*MIN(Dx/MAXVAL(uvel+sos),Dy/MAXVAL(vvel+sos))
    nptst = ABS(d-c)/Dt

    OPEN(1, file = 'initial.dat', status = 'replace')   
    WRITE(1,'(A)') 'VARIABLES = "X","Y","rho","uvel","vvel","pre"'
    WRITE(1,*) 'ZONE I=',nptsx+1,', J=',nptsy+1,', ZONETYPE=ORDERED,'
    WRITE(1,*) 'DATAPACKING=POINT, SOLUTIONTIME=0.0'
    DO i = 0, nptsx 
      DO j=0, nptsy      
        WRITE(1, *) x(i), y(j), rho(i,j), uvel(i,j), vvel(i,j), pre(i,j)
      END DO
    END DO
    CLOSE(1)

    END SUBROUTINE init2d


!************** Construct solution vector *************************
    SUBROUTINE solvec
      IMPLICIT NONE
      INTEGER :: i,j

      DO i = 0, nptsx
        DO j=0, nptsy
          e(i,j)=pre(i,j)/((gam-1.0d0))
          et(i,j)=0.5d0*(uvel(i,j)**2.0d0+vvel(i,j)**2.0d0) + e(i,j)
          u(1,i,j)=rho(i,j)
          u(2,i,j)=rho(i,j)*uvel(i,j)
          u(3,i,j)=rho(i,j)*vvel(i,j)
          u(4,i,j)=rho(i,j)*et(i,j)
        END DO
      END DO  

      OPEN(20, file = 'initial_solvec.dat', status = 'replace')   
      WRITE(20,'(A)') 'VARIABLES = "X","Y","u1","u2","u3"'
      WRITE(20,*) 'ZONE I=',nptsx+1,', J=',nptsy+1,', ZONETYPE=ORDERED,'
      WRITE(20,*) 'DATAPACKING=POINT, SOLUTIONTIME=0.0'
      DO i = 0, nptsx       
        DO j=0, nptsy
          WRITE(20, *) x(i), y(j), u(1,i,j), u(2,i,j), u(3,i,j)
        END DO
      END DO
      CLOSE(20) 
   
      u_ic=u ! To compute l2 with initial condition

    END SUBROUTINE solvec

!************** Compute viscous flux vectors f_visc and g_visc *************************
    SUBROUTINE flux_visc
      IMPLICIT NONE
      INTEGER :: i, j
      REAL :: C1, S, Cp, Pr  
      ! 1-D arrays
      REAL, DIMENSION(:) :: termx(0:nptsx), termx1(0:nptsx)
      REAL, DIMENSION(:) :: termy(0:nptsy), termy1(0:nptsy)
      REAL, DIMENSION(:) :: mu(0:nptsx,0:nptsy),kappa(0:nptsx,0:nptsy),lambda(0:nptsx,0:nptsy)
      ! 2-D arrays
      REAL, DIMENSION(:,:) :: tau_xx(0:nptsx,0:nptsy),tau_xy(0:nptsx,0:nptsy)
      REAL, DIMENSION(:,:) :: tau_yx(0:nptsx,0:nptsy),tau_yy(0:nptsx,0:nptsy)
      REAL, DIMENSION(:,:) :: dudx(0:nptsx,0:nptsy),dvdy(0:nptsx,0:nptsy)
      REAL, DIMENSION(:,:) :: dvdx(0:nptsx,0:nptsy),dudy(0:nptsx,0:nptsy)
      REAL, DIMENSION(:,:) :: dtdx(0:nptsx,0:nptsy),dtdy(0:nptsx,0:nptsy)

      ! Coefficients to compute dynamic viscousity
      C1=1.458d0*10.d0**(-6.0d0)
      S=110.4d0
      Cp=1.0049d0 ! For air at 300K
      Pr=0.707d0  ! For air at 300K

      DO i=0, nptsx
        DO j=0, nptsy
          ! Molecular viscosity
           mu(i,j)=0.01d0 ! Re=100
          ! Bulk viscosity
          lambda(i,j)=-(2.0d0/3.0d0)*mu(i,j)
          ! Thermal conductivity
          kappa(i,j)=(mu(i,j)*Cp)/Pr
        END DO
      END DO

      ! dudx
      DO j=0, nptsy
        termx(0:nptsx)=uvel(0:nptsx,j)
        termx1(0:nptsx)=MATMUL(Dx1(0:nptsx,0:nptsx),termx(0:nptsx))
        dudx(0:nptsx,j)=termx1(0:nptsx)
      END DO

      ! dvdx
      DO j=0, nptsy
        termx(0:nptsx)=vvel(0:nptsx,j)
        termx1(0:nptsx)=MATMUL(Dx1(0:nptsx,0:nptsx),termx(0:nptsx))
        dvdx(0:nptsx,j)=termx1(0:nptsx)
      END DO

      ! dvdy
      DO i=0, nptsx
        termy(0:nptsy)=vvel(i,0:nptsy)
        termy1(0:nptsy)=MATMUL(Dy1(0:nptsy,0:nptsy),termy(0:nptsy))
        dvdy(i,0:nptsy)=termy1(0:nptsy)
      END DO

      ! dudy
      DO i=0, nptsx
        termy(0:nptsy)=uvel(i,0:nptsy)
        termy1(0:nptsy)=MATMUL(Dy1(0:nptsy,0:nptsy),termy(0:nptsy))
        dudy(i,0:nptsy)=termy1(0:nptsy)
      END DO

      ! Stresses
      DO i=0, nptsx
        DO j=0, nptsy
          ! Normal Stresses
          tau_xx(i,j)=lambda(i,j)*(dudx(i,j)+dvdy(i,j))+2.0d0*mu(i,j)*dudx(i,j)
          tau_yy(i,j)=lambda(i,j)*(dudx(i,j)+dvdy(i,j))+2.0d0*mu(i,j)*dvdy(i,j)
          ! Shear stresses
          tau_xy(i,j)=mu(i,j)*(dvdx(i,j)+dudy(i,j))
          tau_yx(i,j)=tau_xy(i,j)
        END DO
      END DO

      ! Thermal gradients
      ! dtdx
      DO j=0, nptsy
        termx(0:nptsx)=temp(0:nptsx,j)
        termx1(0:nptsx)=MATMUL(Dx1(0:nptsx,0:nptsx),termx(0:nptsx))
        dtdx(0:nptsx,j)=termx1(0:nptsx)
      END DO
    
      ! dtdy
      DO i=0, nptsx
        termy(0:nptsy)=temp(i,0:nptsy)
        termy1(0:nptsy)=MATMUL(Dy1(0:nptsy,0:nptsy),termy(0:nptsy))
        dtdy(i,0:nptsy)=termy1(0:nptsy)
      END DO

      ! viscous fluxes
      DO i=0, nptsx
        DO j=0, nptsy
        f_visc(1,i,j)=0.0d0
        f_visc(2,i,j)=tau_xx(i,j) 
        f_visc(3,i,j)=tau_xy(i,j)
        f_visc(4,i,j)=kappa(i,j)*dtdx(i,j)+uvel(i,j)*tau_xx(i,j)+vvel(i,j)*tau_xy(i,j)

        g_visc(1,i,j)=0.0d0
        g_visc(2,i,j)=tau_yx(i,j)
        g_visc(3,i,j)=tau_yy(i,j)
        g_visc(4,i,j)=kappa(i,j)*dtdy(i,j)+uvel(i,j)*tau_yx(i,j)+vvel(i,j)*tau_yy(i,j)
        END DO
      END DO

    END SUBROUTINE flux_visc


!------------------------------------------------------

  SUBROUTINE matrix(idir,npts)
  IMPLICIT NONE
  INTEGER :: i,j
  INTEGER :: idir, npts
  REAL(DP) :: x1

  allocate (Hmatinv(0:npts,0:npts),Qmat(0:npts,0:npts))
  allocate (Mmat(0:npts,0:npts),Smat(0:npts,0:npts))
  allocate (Bmat(0:npts,0:npts))
  allocate (D1mat(0:npts,0:npts),D2mat(0:npts,0:npts))

! Construct scheme matrices

  select case (porder)
!****** 2nd order accurate *****
  case(2)
! Hmatrixinv
  Hmatinv = 0.0
  Hmatinv(0,0) = 2.0/h(idir)
  do i = 1,npts-1
    Hmatinv(i,i) = 1.0/h(idir)
  end do
  Hmatinv(npts,npts) = 2.0/h(idir)

! Qmatrix
  Qmat = 0.0
  Qmat(0,0) = -0.5
  Qmat(0,1) = 0.5
  do i = 1, npts-1
    Qmat(i,i-1) = -0.5
    Qmat(i,i) = 0.0
    Qmat(i,i+1) = 0.5
  end do
  Qmat(npts,npts-1) = -0.5
  Qmat(npts,npts) = 0.5

! Mmatrix
  Mmat = 0.0
  Mmat(0,0) = 1.0
  Mmat(0,1) = -1.0
  Mmat(1,0) = -1.0
  do i = 1, npts-1
    Mmat(i,i-1) = -1.0
    Mmat(i,i) = 2.0
    Mmat(i,i+1) = -1.0
  end do
  Mmat(npts,npts-1) = -1.0
  Mmat(npts,npts) = 1.0
  Mmat(:,:) = Mmat(:,:)/h(idir)

! Bmatrix
  Bmat = 0.0
  Bmat(0,0) = -1.0
  Bmat(npts,npts) = 1.0

! Smatrix
  Smat = 0.0
  Smat(0,0) = 3./2.
  Smat(0,1) = -2.0
  Smat(0,2) = 0.5
  do i = 1,npts-1
    Smat(i,i) = 1.0
  end do
  Smat(npts,npts) = 3./2.
  Smat(npts,npts-1) = -2.0
  Smat(npts,npts-2) = 0.5
  Smat(:,:) = Smat(:,:)/h(idir)

! D1matrix
  D1mat = matmul(Hmatinv,Qmat) 
! D2matrix
  D2mat = matmul(Hmatinv,(-Mmat+matmul(Bmat,Smat)))


!****** 4th order accurate *****
  case(4)
! Hmatinv
  Hmatinv = 0.0
  Hmatinv(0,0) = 1./(17./48.)
  Hmatinv(1,1) = 1./(59./48.)
  Hmatinv(2,2) = 1./(43./48.)
  Hmatinv(3,3) = 1./(49./48.)
  do i = 4, npts-4
    Hmatinv(i,i) = 1.0
  end do
  do i = 0, 3
    Hmatinv(npts-i,npts-i) = Hmatinv(i,i)
  end do
  Hmatinv = Hmatinv/h(idir)

! Qmatrix
  Qmat = 0.0
  Qmat(0,0) = -0.5
  Qmat(0,1) = 59./96.
  Qmat(0,2) = -1./12.
  Qmat(0,3) = -1./32.
  Qmat(1,0) = -59./96.
  Qmat(2,0) = 1./12.
  Qmat(3,0) = 1./32.
  Qmat(1,1) = 0.0
  Qmat(1,2) = 59./96.
  Qmat(2,1) = -59./96.
  Qmat(2,2) = 0.0
  Qmat(2,3) = 59./96.
  Qmat(2,4) = -1./12.
  Qmat(3,1) = 0.0
  Qmat(3,2) = -59./96.
  Qmat(3,3) = 0.0
  Qmat(3,4) = 2./3.
  Qmat(3,5) = -1./12.
  do i = 4, npts-4
    Qmat(i,i-2) = 1./12.
    Qmat(i,i-1) = -2./3.
    Qmat(i,i) = 0.0
    Qmat(i,i+1) = 2./3.
    Qmat(i,i+2) = -1./12.
  end do
  Qmat(npts,npts) = 0.5
  Qmat(npts,npts-1) = -59./96.
  Qmat(npts,npts-2) = 1./12.
  Qmat(npts,npts-3) = 1./32.
  Qmat(npts-1,npts) = 59./96.
  Qmat(npts-1,npts-2) = -59./96.
  Qmat(npts-2,npts) = -1./12.
  Qmat(npts-2,npts-1) = 59./96.
  Qmat(npts-2,npts-3) = -59./96.
  Qmat(npts-2,npts-4) = 1./12.
  Qmat(npts-3,npts) = -1./32.
  Qmat(npts-3,npts-2) = 59./96.
  Qmat(npts-3,npts-4) = -2./3.
  Qmat(npts-3,npts-5) = 1./12.

! D1matrix
  D1mat = matmul(Hmatinv,Qmat) 

! Dmat2
  D2mat = 0.0
  D2mat(0,0) = 2.0
  D2mat(0,1) = -5.0
  D2mat(0,2) = 4.0
  D2mat(0,3) = -1.0
  D2mat(1,0) = 1.0
  D2mat(2,0) = -4./43.
  D2mat(3,0) = -1./49.
  D2mat(1,1) = -2.0
  D2mat(1,2) = 1.0
  D2mat(1,3) = 0.0
  D2mat(2,1) = 59./43.
  D2mat(2,2) = -110./43.
  D2mat(2,3) = 59./43.
  D2mat(2,4) = -4./43.
  D2mat(3,1) = 0.0
  D2mat(3,2) = 59./49.
  D2mat(3,3) = -118./49.
  D2mat(3,4) = 64./49.
  D2mat(3,5) = -4./49.
  do i = 4,npts
    D2mat(i,i) = 4./3.
    D2mat(i,i-1) = -1./12.
  end do
  do i = 4,npts-3
    D2mat(i,i+1) = -5./2.
    D2mat(i,i+2) = 4./3.
    D2mat(i,i+3) = 1./12.
  end do
  D2mat(npts-1,npts) = -5./2.
  D2mat(npts-2,npts) = 4./3.
  D2mat(npts-2,npts-1) = -5./2.

  D2mat = D2mat/h(idir)**2

! Smat
  Smat = 0.0
  Smat(0,0) = 11./6.
  Smat(0,1) = -3.
  Smat(0,2) = 3./2.
  Smat(0,3) = -1./3.
  Smat(1,1) = 1.0
  Smat(npts,npts) = 11./6.
  Smat(npts,npts-1) = -3.
  Smat(npts,npts-2) = 3./2.
  Smat(npts,npts-3) = -1./3.
  Smat(npts-1,npts-1) = 1.0
  Smat = Smat/h(idir)

!****** 6th order accurate *****
  case(6)
! Hmatinv
  Hmatinv = 0.0
  Hmatinv(0,0) = 1./(13649./43200.)
  Hmatinv(1,1) = 1./(12013./8640.)
  Hmatinv(2,2) = 1./(2711./4320.)
  Hmatinv(3,3) = 1./(5359./4320.)
  Hmatinv(4,4) = 1./(7877./8640.)
  Hmatinv(5,5) = 1./(43801./43200.)
  do i = 6, npts-6
    Hmatinv(i,i) = 1.0
  end do
  do i = 0, 5
    Hmatinv(npts-i,npts-i) = Hmatinv(i,i)
  end do
  Hmatinv = Hmatinv/h(idir)

! Qmatrix
  x1 = 342523./518400.
  Qmat = 0.0
  Qmat(0,0) = -0.5
  Qmat(0,1) = x1 - 953./16200.
  Qmat(1,0) = -Qmat(0,1)
  Qmat(0,2) = -4.0*x1 + 715489./259200.
  Qmat(2,0) = -Qmat(0,2)
  Qmat(0,3) = 6.0*x1 - 62639./14400.
  Qmat(3,0) = -Qmat(0,3)
  Qmat(0,4) = -4.0*x1 + 147127./51840.
  Qmat(4,0) = -Qmat(0,4)
  Qmat(0,5) = x1 - 89387./129600.
  Qmat(5,0) = -Qmat(0,5)
  Qmat(1,2) = 10.*x1 - 57139./8640.
  Qmat(2,1) = -Qmat(1,2)
  Qmat(1,3) = -20.*x1 + 745733./51840.
  Qmat(3,1) = -Qmat(1,3)
  Qmat(1,4) = 15.*x1 - 18343./1728.
  Qmat(4,1) = -Qmat(1,4)
  Qmat(1,5) = -4.*x1 + 240569./86400.
  Qmat(5,1) = -Qmat(1,5)
  Qmat(2,3) = 20.*x1 - 176839./12960.
  Qmat(3,2) = -Qmat(2,3)
  Qmat(2,4) = -20.*x1 + 242111./17280.
  Qmat(4,2) = -Qmat(2,4)
  Qmat(2,5) = 6.*x1 - 182261./43200.
  Qmat(5,2) = -Qmat(2,5)
  Qmat(3,4) = 10.*x1 - 165041./25920.
  Qmat(4,3) = -Qmat(3,4)
  Qmat(3,5) = -4.*x1 + 710473./259200.
  Qmat(5,3) = -Qmat(3,5)
  Qmat(3,6) = 1./60.
  Qmat(6,3) = -Qmat(3,6)
  Qmat(4,5) = x1
  Qmat(5,4) = -Qmat(4,5)
  Qmat(4,6) = -3./20.
  Qmat(4,7) = 1./60.
  Qmat(5,6) = 3./4.
  Qmat(5,7) = -3./20.
  Qmat(5,8) = 1./60.

  do i = 6, npts-6
    Qmat(i,i-3) = -1./60.
    Qmat(i,i-2) = 3./20.
    Qmat(i,i-1) = -3./4.
    Qmat(i,i) = 0.0
    Qmat(i,i+1) = 3./4.
    Qmat(i,i+2) = -3./20.
    Qmat(i,i+3) = 1./60.
  end do
  do i = 0,5
    do j = 0,5
      Qmat(npts-i,npts-j) = -Qmat(i,j)
    end do
  end do
  Qmat(npts-3,npts-6) = -Qmat(3,6)
  Qmat(npts-4,npts-6) = -Qmat(4,6)
  Qmat(npts-4,npts-7) = -Qmat(4,7)
  Qmat(npts-5,npts-6) = -Qmat(5,6)
  Qmat(npts-5,npts-7) = -Qmat(5,7)
  Qmat(npts-5,npts-8) = -Qmat(5,8)

! D1matrix
  D1mat = matmul(Hmatinv,Qmat) 

  end select

  END SUBROUTINE matrix


!************** Compute flux vector in x *************************
    SUBROUTINE flux
    IMPLICIT NONE
      INTEGER :: i, j, ivar
      REAL, DIMENSION(:) :: termx(0:nptsx), termx1(0:nptsx), termx2(0:nptsx), termx3(0:nptsx)
      REAL, DIMENSION(:) :: termx4(0:nptsx), termx5(0:nptsx), termx6(0:nptsx), termx7(0:nptsx)
      REAL, DIMENSION(:) :: termx8(0:nptsx), termx9(0:nptsx), termx10(0:nptsx), termx11(0:nptsx)
      REAL, DIMENSION(:) :: termx12(0:nptsx), termx13(0:nptsx)

       phi(1,0:nptsx,0:nptsy)=1.0d0
       phi(2,0:nptsx,0:nptsy)=uold(2,0:nptsx,0:nptsy)/uold(1,0:nptsx,0:nptsy) ! u
       phi(3,0:nptsx,0:nptsy)=uold(3,0:nptsx,0:nptsy)/uold(1,0:nptsx,0:nptsy) ! v
       phi(4,0:nptsx,0:nptsy)=0.5d0*((uold(2,0:nptsx,0:nptsy)/uold(1,0:nptsx,0:nptsy))**2.0d0+ & 
                                     (uold(3,0:nptsx,0:nptsy)/uold(1,0:nptsx,0:nptsy))**2.0d0)

       rho_u(0:nptsx,0:nptsy)=uold(1,0:nptsx,0:nptsy)*uold(2,0:nptsx,0:nptsy)/uold(1,0:nptsx,0:nptsy)

!********** Compute 1st split flux term *********
       DO ivar=1, nvar
         f1(ivar,0:nptsx,0:nptsy)=rho_u(0:nptsx,0:nptsy)*phi(ivar,0:nptsx,0:nptsy)
       END DO

       DO ivar=1, nvar
         DO j=0, nptsy
          termx(0:nptsx)=f1(ivar,0:nptsx,j)
          termx1(0:nptsx)=MATMUL(Dx1(0:nptsx,0:nptsx),termx(0:nptsx))
          df1dx(ivar,0:nptsx,j)=termx1(0:nptsx)
         END DO
       END DO

       split_x1(1:nvar,0:nptsx,0:nptsy)=df1dx(1:nvar,0:nptsx,0:nptsy)

!********** Compute 2nd split flux term *********
       f2(0:nptsx,0:nptsy)=rho_u(0:nptsx,0:nptsy)

       DO j=0, nptsy
         termx2(0:nptsx)=f2(0:nptsx,j)
         termx3(0:nptsx)=MATMUL(Dx1(0:nptsx,0:nptsx),termx2(0:nptsx))
         df2dx(0:nptsx,j)=termx3(0:nptsx)
       END DO

       DO ivar=1, nvar
         split_x2(ivar,0:nptsx,0:nptsy)=phi(ivar,0:nptsx,0:nptsy)*df2dx(0:nptsx,0:nptsy)
       END DO

!********** Compute 3rd split flux term *********
       DO ivar=1,nvar
         f3(ivar,0:nptsx,0:nptsy)=uold(1,0:nptsx,0:nptsy)*phi(ivar,0:nptsx,0:nptsy) ! rho*phi
       END DO

       DO ivar=1, nvar
         DO j=0, nptsy
          termx4(0:nptsx)=f3(ivar,0:nptsx,j)
          termx5(0:nptsx)=MATMUL(Dx1(0:nptsx,0:nptsx),termx4(0:nptsx))
          df3dx(ivar,0:nptsx,j)=termx5(0:nptsx)
         END DO
       END DO

       DO ivar=1, nvar
         split_x3(ivar,0:nptsx,0:nptsy)=uold(2,0:nptsx,0:nptsy)/uold(1,0:nptsx,0:nptsy)*df3dx(ivar,0:nptsx,0:nptsy) ! u*df3dx
       END DO

!********** Compute 4th split flux term *********
       DO ivar=1,nvar
         f4(ivar,0:nptsx,0:nptsy)=uold(2,0:nptsx,0:nptsy)/uold(1,0:nptsx,0:nptsy)*phi(ivar,0:nptsx,0:nptsy) ! u*phi
       END DO

       DO ivar=1, nvar
         DO j=0, nptsy
          termx6(0:nptsx)=f4(ivar,0:nptsx,j)
          termx7(0:nptsx)=MATMUL(Dx1(0:nptsx,0:nptsx),termx6(0:nptsx))
          df4dx(ivar,0:nptsx,j)=termx7(0:nptsx)
         END DO
       END DO

       DO ivar=1, nvar
         split_x4(ivar,0:nptsx,0:nptsy)=uold(1,0:nptsx,0:nptsy)*df4dx(ivar,0:nptsx,0:nptsy) ! rho*df4dx
       END DO

!********** Compute 5th split flux term *********
       DO ivar=1,nvar
         f5(ivar,0:nptsx,0:nptsy)=phi(ivar,0:nptsx,0:nptsy) ! phi
       END DO

       DO ivar=1, nvar
         DO j=0, nptsy
          termx8(0:nptsx)=f5(ivar,0:nptsx,j)
          termx9(0:nptsx)=MATMUL(Dx1(0:nptsx,0:nptsx),termx8(0:nptsx))
          df5dx(ivar,0:nptsx,j)=termx9(0:nptsx)
         END DO
       END DO

       DO ivar=1, nvar
         split_x5(ivar,0:nptsx,0:nptsy)=rho_u(0:nptsx,0:nptsy)*df5dx(ivar,0:nptsx,0:nptsy) ! rho_u*df5dx
       END DO

!********** Compute 6th split flux term *********
       f6(0:nptsx,0:nptsy)=uold(2,0:nptsx,0:nptsy)/uold(1,0:nptsx,0:nptsy) ! u

       DO j=0, nptsy
          termx10(0:nptsx)=f6(0:nptsx,j)
          termx11(0:nptsx)=MATMUL(Dx1(0:nptsx,0:nptsx),termx10(0:nptsx))
          df6dx(0:nptsx,j)=termx11(0:nptsx)
       END DO

       DO ivar=1, nvar
         split_x6(ivar,0:nptsx,0:nptsy)=uold(1,0:nptsx,0:nptsy)*phi(ivar,0:nptsx,0:nptsy)*df6dx(0:nptsx,0:nptsy) ! rho*phi*df6dx
       END DO

!********** Compute 7th split flux term *********
        f7(0:nptsx,0:nptsy)=uold(1,0:nptsx,0:nptsy) ! rho

       DO j=0, nptsy
          termx12(0:nptsx)=f7(0:nptsx,j)
          termx13(0:nptsx)=MATMUL(Dx1(0:nptsx,0:nptsx),termx12(0:nptsx))
          df7dx(0:nptsx,j)=termx13(0:nptsx)
       END DO

       DO ivar=1, nvar
         split_x7(ivar,0:nptsx,0:nptsy)=uold(2,0:nptsx,0:nptsy)/uold(1,0:nptsx,0:nptsy)*phi(ivar,0:nptsx,0:nptsy) &
                                                                             *df7dx(0:nptsx,0:nptsy) ! u*phi*df7dx
       END DO

       dfdx= alpha_split*split_x1 + beta_split*(split_x2+split_x3+split_x4) + & 
             (1-alpha_split-2*beta_split)*(split_x5+split_x6+split_x7)

    END SUBROUTINE flux


!************** Derivative of pressure terms with x *************
     SUBROUTINE dpre_dx
      INTEGER :: i, j, ivar, k
      REAL, DIMENSION(:) :: termx14(0:nptsx), termx15(0:nptsx)
    
     pre_mat_x(1,0:nptsx,0:nptsy)=0.0d0
     pre_mat_x(2,0:nptsx,0:nptsy)=pre(0:nptsx,0:nptsy) 
     pre_mat_x(3,0:nptsx,0:nptsy)=0.0d0
     pre_mat_x(4,0:nptsx,0:nptsy)=(uold(2,0:nptsx,0:nptsy)/uold(1,0:nptsx,0:nptsy)*pre(0:nptsx,0:nptsy))*(1.0d0/(gam-1.0d0)+1.0d0) 

      DO ivar=1, nvar
       DO j=0, nptsy
          termx14(0:nptsx)=pre_mat_x(ivar,0:nptsx,j)
          termx15(0:nptsx)=MATMUL(Dx1(0:nptsx,0:nptsx),termx14(0:nptsx))
          dpdx(ivar,0:nptsx,j)=termx15(0:nptsx)
       END DO
      END DO

    END SUBROUTINE dpre_dx


!************** Compute flux vector in y *************************
    SUBROUTINE glux
    IMPLICIT NONE
      INTEGER :: i, j, ivar
      REAL, DIMENSION(:) :: termy(0:nptsy), termy1(0:nptsy), termy2(0:nptsy), termy3(0:nptsy)
      REAL, DIMENSION(:) :: termy4(0:nptsy), termy5(0:nptsy), termy6(0:nptsy), termy7(0:nptsy)
      REAL, DIMENSION(:) :: termy8(0:nptsy), termy9(0:nptsy), termy10(0:nptsy), termy11(0:nptsy)
      REAL, DIMENSION(:) :: termy12(0:nptsy), termy13(0:nptsy)

       phi(1,0:nptsx,0:nptsy)=1.0d0
       phi(2,0:nptsx,0:nptsy)=uold(2,0:nptsx,0:nptsy)/uold(1,0:nptsx,0:nptsy) ! u
       phi(3,0:nptsx,0:nptsy)=uold(3,0:nptsx,0:nptsy)/uold(1,0:nptsx,0:nptsy) ! v
       phi(4,0:nptsx,0:nptsy)=0.5d0*((uold(2,0:nptsx,0:nptsy)/uold(1,0:nptsx,0:nptsy))**2.0d0+ & 
                                     (uold(3,0:nptsx,0:nptsy)/uold(1,0:nptsx,0:nptsy))**2.0d0)

       rho_v(0:nptsx,0:nptsy)=uold(1,0:nptsx,0:nptsy)*uold(3,0:nptsx,0:nptsy)/uold(1,0:nptsx,0:nptsy)

!********** Compute 1st split flux term *********
       DO ivar=1, nvar
         g1(ivar,0:nptsx,0:nptsy)=rho_v(0:nptsx,0:nptsy)*phi(ivar,0:nptsx,0:nptsy)
       END DO

       DO ivar=1, nvar
         DO i=0, nptsx
          termy(0:nptsy)=g1(ivar,i,0:nptsy)
          termy1(0:nptsy)=MATMUL(Dy1(0:nptsy,0:nptsy),termy(0:nptsy))
          dg1dy(ivar,i,0:nptsy)=termy1(0:nptsy)
         END DO
       END DO

       split_y1(1:nvar,0:nptsx,0:nptsy)=dg1dy(1:nvar,0:nptsx,0:nptsy)

!********** Compute 2nd split flux term *********
       g2(0:nptsx,0:nptsy)=rho_v(0:nptsx,0:nptsy)

       DO i=0, nptsx
         termy2(0:nptsy)=g2(i,0:nptsy)
         termy3(0:nptsy)=MATMUL(Dy1(0:nptsy,0:nptsy),termy2(0:nptsy))
         dg2dy(i,0:nptsy)=termy3(0:nptsy)
       END DO

       DO ivar=1, nvar
         split_y2(ivar,0:nptsx,0:nptsy)=phi(ivar,0:nptsx,0:nptsy)*dg2dy(0:nptsx,0:nptsy)
       END DO

!********** Compute 3rd split flux term *********
       DO ivar=1,nvar
         g3(ivar,0:nptsx,0:nptsy)=uold(1,0:nptsx,0:nptsy)*phi(ivar,0:nptsx,0:nptsy) ! rho*phi
       END DO

       DO ivar=1, nvar
         DO i=0, nptsx
          termy4(0:nptsy)=g3(ivar,i,0:nptsy)
          termy5(0:nptsy)=MATMUL(Dy1(0:nptsy,0:nptsy),termy4(0:nptsy))
          dg3dy(ivar,i,0:nptsy)=termy5(0:nptsy)
         END DO
       END DO

       DO ivar=1, nvar
         split_y3(ivar,0:nptsx,0:nptsy)=uold(3,0:nptsx,0:nptsy)/uold(1,0:nptsx,0:nptsy)*dg3dy(ivar,0:nptsx,0:nptsy) ! v*dg3dy
       END DO

!********** Compute 4th split flux term *********
       DO ivar=1,nvar
         g4(ivar,0:nptsx,0:nptsy)=uold(3,0:nptsx,0:nptsy)/uold(1,0:nptsx,0:nptsy)*phi(ivar,0:nptsx,0:nptsy) ! v*phi
       END DO

       DO ivar=1, nvar
         DO i=0, nptsx
          termy6(0:nptsy)=g4(ivar,i,0:nptsy)
          termy7(0:nptsy)=MATMUL(Dy1(0:nptsy,0:nptsy),termy6(0:nptsy))
          dg4dy(ivar,i,0:nptsy)=termy7(0:nptsy)
         END DO
       END DO

       DO ivar=1, nvar
         split_y4(ivar,0:nptsx,0:nptsy)=uold(1,0:nptsx,0:nptsy)*dg4dy(ivar,0:nptsx,0:nptsy) ! rho*dg4dy
       END DO

!********** Compute 5th split flux term *********
       DO ivar=1,nvar
         g5(ivar,0:nptsx,0:nptsy)=phi(ivar,0:nptsx,0:nptsy) ! phi
       END DO

       DO ivar=1, nvar
         DO i=0, nptsx
          termy8(0:nptsy)=g5(ivar,i,0:nptsy)
          termy9(0:nptsy)=MATMUL(Dy1(0:nptsy,0:nptsy),termy8(0:nptsy))
          dg5dy(ivar,i,0:nptsy)=termy9(0:nptsy)
         END DO
       END DO

       DO ivar=1, nvar
         split_y5(ivar,0:nptsx,0:nptsy)=rho_v(0:nptsx,0:nptsy)*dg5dy(ivar,0:nptsx,0:nptsy) ! rho_v*dg5dy
       END DO

!********** Compute 6th split flux term *********
       g6(0:nptsx,0:nptsy)=uold(3,0:nptsx,0:nptsy)/uold(1,0:nptsx,0:nptsy) ! v

       DO i=0, nptsx
          termy10(0:nptsy)=g6(i,0:nptsy)
          termy11(0:nptsy)=MATMUL(Dy1(0:nptsy,0:nptsy),termy10(0:nptsy))
          dg6dy(i,0:nptsy)=termy11(0:nptsy)
       END DO

       DO ivar=1, nvar
         split_y6(ivar,0:nptsx,0:nptsy)=uold(1,0:nptsx,0:nptsy)*phi(ivar,0:nptsx,0:nptsy)*dg6dy(0:nptsx,0:nptsy) ! rho*phi*df6dx
       END DO

!********** Compute 7th split flux term *********
        g7(0:nptsx,0:nptsy)=uold(1,0:nptsx,0:nptsy) ! rho

       DO i=0, nptsx
          termy12(0:nptsy)=g7(i,0:nptsy)
          termy13(0:nptsy)=MATMUL(Dy1(0:nptsy,0:nptsy),termy12(0:nptsy))
          dg7dy(i,0:nptsy)=termy13(0:nptsy)
       END DO

       DO ivar=1, nvar
         split_y7(ivar,0:nptsx,0:nptsy)=uold(3,0:nptsx,0:nptsy)/uold(1,0:nptsx,0:nptsy)*phi(ivar,0:nptsx,0:nptsy) &
                                                                             *dg7dy(0:nptsx,0:nptsy) ! v*phi*dg7dy
       END DO

       dgdy= alpha_split*split_y1 + beta_split*(split_y2+split_y3+split_y4) + & 
             (1-alpha_split-2*beta_split)*(split_y5+split_y6+split_y7)

    END SUBROUTINE glux


!************** Derivative of pressure terms with y *************
     SUBROUTINE dpre_dy
      INTEGER :: i, j, ivar, k
      REAL, DIMENSION(:) :: termy14(0:nptsy), termy15(0:nptsy)
    
     pre_mat_y(1,0:nptsx,0:nptsy)=0.0d0
     pre_mat_y(2,0:nptsx,0:nptsy)=0.0d0
     pre_mat_y(3,0:nptsx,0:nptsy)=pre(0:nptsx,0:nptsy) 
     pre_mat_y(4,0:nptsx,0:nptsy)=(uold(3,0:nptsx,0:nptsy)/uold(1,0:nptsx,0:nptsy)*pre(0:nptsx,0:nptsy))*(1.0d0/(gam-1.0d0)+1.0d0) 

      DO ivar=1, nvar
       DO i=0, nptsx
          termy14(0:nptsy)=pre_mat_y(ivar,i,0:nptsy)
          termy15(0:nptsy)=MATMUL(Dy1(0:nptsy,0:nptsy),termy14(0:nptsy))
          dpdy(ivar,i,0:nptsy)=termy15(0:nptsy)
       END DO
      END DO

    END SUBROUTINE dpre_dy


!************** Derivative of viscous flux wrt x *************
     SUBROUTINE dflux_dx_visc
      INTEGER :: i, j, ivar, k
      REAL, DIMENSION(:) :: termx(0:nptsx), termx1(0:nptsx)

      DO ivar=1, nvar
       DO j=0, nptsy
          termx(0:nptsx)=f_visc(ivar,0:nptsx,j)
          termx1(0:nptsx)=MATMUL(Dx1(0:nptsx,0:nptsx),termx(0:nptsx))
          dfdx_visc(ivar,0:nptsx,j)=termx1(0:nptsx)
       END DO
      END DO

    END SUBROUTINE dflux_dx_visc


!************** Derivative of viscous flux wrt y *************
     SUBROUTINE dglux_dy_visc
      INTEGER :: i, j, ivar, k, l
      REAL, DIMENSION(:) :: termy(0:nptsy), termy1(0:nptsy)

      DO ivar= 1, nvar
        DO i = 0, nptsx
          termy(0:nptsy) = g_visc(ivar, i, 0:nptsy)
          termy1(0:nptsy) = MATMUL(Dy1(0:nptsy, 0:nptsy), termy(0:nptsy))
          dgdy_visc(ivar, i, 0:nptsy) = termy1(0:nptsy)
        END DO
      END DO

    END SUBROUTINE dglux_dy_visc


!************** SAT Left (left on x axis) *************
    SUBROUTINE SAT_x_left
      INTEGER :: i, j, ivar, k
      REAL(DP) :: sigma
      ALLOCATE(temp1d(1:nvar,0:nptsx,0:nptsy),correction(1:nvar,0:nptsy))

       sigma=1.7d0 ! Energy stable for sigma>=0.5

      ! Target values at left boundary
      ! Perioidc BC; so left boundary=right, u_0=u_n
        g_x_0(1,0:nptsy)=1.0d0
        g_x_0(2,0:nptsy)=0.0d0
        g_x_0(3,0:nptsy)=0.0d0
        g_x_0(4,0:nptsy)=2.5d0

      hinv_x0=Hinvx(0,0)

      correction(1:nvar,0:nptsy)=sigma*hinv_x0*(uold(1:nvar,0,0:nptsy)-g_x_0(1:nvar,0:nptsy))

      ! Intializing SAT_x0 as a (nx X ny) matrix of zeros
      SAT_x0(1:nvar,0:nptsx,0:nptsy)=0.0d0
      ! Replacing the left-hand side of the matrix (i.e., x=0 of all rows) with SAT_x0 values; the remaining values are zeros

      DO j=0, nptsy
        SAT_x0(1:nvar,0,j)=correction(1:nvar,j)
      END DO


     DEALLOCATE(temp1d,correction)

    END SUBROUTINE SAT_x_left


!************** SAT Right (right on x axis) *************
    SUBROUTINE SAT_x_right
      INTEGER :: i, j, ivar, k
      REAL(DP) :: sigma
      ALLOCATE(temp1d(1:nvar,0:nptsx,0:nptsy),correction(1:nvar,0:nptsy))

       sigma=1.7d0  ! Energy stable for sigma>=0.5

       ! Target values at right boundary
      ! Perioidc BC; so right boundary=left, u_0=u_n
        g_x_n(1,0:nptsy)=1.0d0
        g_x_n(2,0:nptsy)=0.0d0
        g_x_n(3,0:nptsy)=0.0d0
        g_x_n(4,0:nptsy)=2.5d0

     hinv_xn=Hinvx(nptsx,nptsx)

      correction(1:nvar,0:nptsy)=sigma*hinv_xn*(uold(1:nvar,nptsx,0:nptsy)-g_x_n(1:nvar,0:nptsy))

      ! Intializing SAT_xn as a (nx X ny) matrix of zeros
      SAT_xn(1:nvar,0:nptsx,0:nptsy)=0.0d0
      ! Replacing the right-hand side of the matrix (i.e., x=nptsx of all rows) with SAT_xn values; the remaining values are zeros
      DO j=0, nptsy
        SAT_xn(1:nvar,nptsx,j)=correction(1:nvar,j)
      END DO

    DEALLOCATE(temp1d,correction)

    END SUBROUTINE SAT_x_right


!************** SAT Bottom (on Y axis) *************
    SUBROUTINE SAT_y_bottom
      INTEGER :: i, j, k
      REAL(DP) :: sigma
      REAL, DIMENSION(:) :: termy(0:nptsy), termy1(0:nptsy)
      ALLOCATE(temp1d(1:nvar,0:nptsx,0:nptsy),correction(1:nvar,0:nptsx))

       sigma=1.7d0 ! Energy stable for sigma>=0.5

      ! Target value at bottom boundary
        g_y_0(1,0:nptsx)=1.0d0
        g_y_0(2,0:nptsx)=0.0d0
        g_y_0(3,0:nptsx)=0.0d0
        g_y_0(4,0:nptsx)=2.5d0

      hinv_y0=Hinvy(0,0)

      correction(1:nvar,0:nptsx)=sigma*hinv_y0*(uold(1:nvar,0:nptsx,0)-g_y_0(1:nvar,0:nptsx))

      ! Intializing SAT_y0 as a (nx X ny) matrix of zeros
      SAT_y0(1:nvar,0:nptsx,0:nptsy)=0.0d0
      ! Replacing the top side of the matrix (i.e., y=0 of all rows) with SAT_y0 values; the remaining values are zeros
      DO i=0, nptsx
        SAT_y0(1:nvar,i,0)=correction(1:nvar,i)
      END DO

      DEALLOCATE(temp1d,correction)

    END SUBROUTINE SAT_y_bottom


!************** SAT Top(on Y axis) *************
    SUBROUTINE SAT_y_top
      INTEGER :: i, j, k
      REAL(DP) :: sigma
      REAL, DIMENSION(:) :: termy(0:nptsy), termy1(0:nptsy)
      ALLOCATE(temp1d(1:nvar,0:nptsx,0:nptsy),correction(1:nvar,0:nptsx))

	  sigma=0.5d0  ! Energy stable for sigma>=0.5


      ! Target value at top boundary
        g_y_n(1,0:nptsx)=1.0d0
        g_y_n(2,0:nptsx)=1.0d0
        g_y_n(3,0:nptsx)=0.0d0
        g_y_n(4,0:nptsx)=3.0d0

      hinv_yn=Hinvy(nptsy,nptsy)

      correction(1:nvar,0:nptsx)=sigma*hinv_yn*(uold(1:nvar,0:nptsx,nptsy)-g_y_n(1:nvar,0:nptsx))

      ! Intializing SAT_y0 as a (nx X ny) matrix of zeros
      SAT_yn(1:nvar,0:nptsx,0:nptsy)=0.0d0
      ! Replacing the bottom side of the matrix (i.e., y=nptsy of all rows) with SAT_yn values; the remaining values are zeros
      DO i=0, nptsx
        SAT_yn(1:nvar,i,nptsy)=correction(1:nvar,i)
      END DO

      DEALLOCATE(temp1d,correction)

    END SUBROUTINE SAT_y_top


!************** RHS *************
    SUBROUTINE right_hand_side

     ! Inviscid
     CALL flux
     CALL dpre_dx
     CALL glux
     CALL dpre_dy

     ! Viscous
     CALL flux_visc
     CALL dflux_dx_visc 
     CALL dglux_dy_visc
     ! SATs
     CALL SAT_x_left
     CALL SAT_x_right
     CALL SAT_y_bottom
     CALL SAT_y_top

     rhs(1:nvar,0:nptsx,0:nptsy)=-dfdx(1:nvar,0:nptsx,0:nptsy)-dpdx(1:nvar,0:nptsx,0:nptsy) +&
                                   dfdx_visc(1:nvar,0:nptsx,0:nptsy)- &
                                   dgdy(1:nvar,0:nptsx,0:nptsy)-dpdy(1:nvar,0:nptsx,0:nptsy) + &
                                   dgdy_visc(1:nvar,0:nptsx,0:nptsy)- &
                                   SAT_x0(1:nvar,0:nptsx,0:nptsy)-SAT_xn(1:nvar,0:nptsx,0:nptsy)- &
                                   SAT_y0(1:nvar,0:nptsx,0:nptsy)-SAT_yn(1:nvar,0:nptsx,0:nptsy)

    END SUBROUTINE right_hand_side


!************** Decompose the solution vector back to primitive variables after each dt *************
! This also computes other essential quantities for visualization. 
! This should be called inside RK-4 (after computing uold at each dt)since this decomposes pressure 
! from the solution vector which will then be needed to compute the flux in the next dt.
    SUBROUTINE decomp(ut)
      INTEGER :: i,j
      REAL(DP), DIMENSION(:, 0:, 0:) :: ut

      DO i=0, nptsx
        DO j=0, nptsy
          rho(i,j)=ut(1,i,j)
          uvel(i,j)=ut(2,i,j)/ut(1,i,j)
          vvel(i,j)=ut(3,i,j)/ut(1,i,j)
          et(i,j)=ut(4,i,j)/ut(1,i,j)
          ke(i,j)=0.5d0*(uvel(i,j)**2.0d0+vvel(i,j)**2.0d0)
          temp(i,j)=(gam-1.0d0)*(et(i,j)-ke(i,j))/rgas
          pre(i,j)=rho(i,j)*rgas*temp(i,j)
          sos(i,j)=SQRT(gam*pre(i,j)/rho(i,j)) ! Speed of sound
          mach(i,j)=uvel(i,j)/sos(i,j)
        END DO
      END DO


    OPEN(unit = 2, file = 'sol.dat', status = 'replace')
    WRITE(2,'(A)') 'VARIABLES = "X", "Y", "rho", "uvel", "vvel", "pre"'
      WRITE(2,*) 'ZONE I=',nptsx+1,', J=',nptsy+1,', ZONETYPE=ORDERED,'
    WRITE(2,*) 'DATAPACKING=POINT, SOLUTIONTIME=0.0'
    DO i = 0, nptsx
      DO j=0, nptsy
     ! WRITE(2, *) x(k), u(1,k), u(2,k)
       WRITE(2, *) x(i), x(j), rho(i,j), uvel(i,j), vvel(i,j), pre(i,j)
      END DO
    END DO
    CLOSE(2)

    END SUBROUTINE decomp


!************** Time integration *************
  ! RK-4
  SUBROUTINE rk4
    USE types_vars
    USE variables
    
    IMPLICIT NONE
    INTEGER :: i,j, k

    uold=u   
      OPEN(unit = 44, file = 'convergence.dat', status = 'replace')
    DO j=0, nptst

!******************************************** STEP 1 **********************************************
       CALL right_hand_side
       k1=Dt*rhs
       uold=u+(k1/2.0d0)
       CALL decomp(uold)

!******************************************** STEP 2 **********************************************
       CALL right_hand_side
       k2= Dt*rhs
       uold=u+(k2/2.0d0)
       CALL decomp(uold) 
     
!******************************************** STEP 3 **********************************************
       CALL right_hand_side
       k3=Dt*rhs
       uold=u+k3
       CALL decomp(uold)
       
!******************************************** STEP 4 **********************************************
       CALL right_hand_side
       k4=Dt*rhs
       uold=u+((1.0d0/6.0d0)*(k1+(2.0d0*k2)+(2.0d0*k3)+k4))
       CALL decomp(uold) 
       u=uold


       t=t+Dt
       WRITE(*,*) 't=',t, ', Max of u1=', MAXVAL(u(1,:,:)), ', Max of u2=', MAXVAL(u(2,:,:)), ', Max of u3=', MAXVAL(u(3,:,:)), &
                                                                                              ', Max of u4=', MAXVAL(u(4,:,:))
       WRITE(*,*) 't=',t, ', Min of u1=', MINVAL(u(1,:,:)), ', Min of u2=', MINVAL(u(2,:,:)), ', Min of u3=', MINVAL(u(3,:,:)), &
                                                                                              ', Min of u4=', MINVAL(u(4,:,:))
      WRITE(44,*) SUM(u(2,:,:))


    END DO
      CLOSE(44)
  END SUBROUTINE rk4

  END MODULE subroutines



  !**************************************************************************
  !**************************************************************************
  !**************************************************************************
  !**************************************************************************
  ! MAIN PROGRAM
  ! Program template 
  PROGRAM CNSE2d_Split_form_lid_cavity
    USE types_vars
    USE variables
    USE subroutines
    INTEGER :: t1, rate, t2

    CALL system_clock(t1, rate)  

    CALL inputs
    CALL memalloc

    CALL grid2d
    CALL init2d

    CALL solvec
    CALL decomp(u)

    CALL matrix(1,nptsx)
    Dx1=D1mat
    Hinvx=Hmatinv
 !   DEALLOCATE(Hmatinv,Qmat,Q_tran,D1mat,D2mat)
    DEALLOCATE(Hmatinv,Qmat,D1mat,D2mat,Mmat,Smat,Bmat)

    CALL matrix(2,nptsy)
    Dy1=D1mat
    Hinvy=Hmatinv
 !   DEALLOCATE(Hmatinv,Qmat,Q_tran,D1mat,D2mat)
    DEALLOCATE(Hmatinv,Qmat,D1mat,D2mat,Mmat,Smat,Bmat)

    CALL rk4

    WRITE(*,*) 'Nx=',nptsx,', Ny=',nptsy, ', CFL=',cfl_ip 
    WRITE(*,*) 'Dx=', Dx, ', Dy=', Dy, ', Dt=', Dt

    CALL dealloc

    CALL system_clock(t2)
    WRITE(*,*) 'Elapsed time=', REAL(t2 - t1)/REAL(rate), 's'

  END PROGRAM CNSE2d_Split_form_lid_cavity
