!****************************************************************
!
!	Improved conservative phase-field LBM solver
!	for two-phase flows in a periodic domain (2D)
!
!----------------------------------------------------------------
!	Based on the following paper:
!
!	A. Fakhari, T. Mitchell, C. Leonardi, and D. Bolster,
!	"Improved locality of the phase-field lattice-Boltzmann model
!	for immiscible fluids at high density ratios",
!	Physical Review E 96, 053301 (2017)
!----------------------------------------------------------------
!
!	written by Abbas Fakhari 10/30/2016
!	03/23/2017:	updated					 
!	10/18/2018:	minor updates
!
!****************************************************************

MODULE SHARE
	IMPLICIT NONE
	INTEGER :: t, X, Y

	INTEGER,PARAMETER :: L0 = 128

	INTEGER,PARAMETER :: tf   = 10000
	INTEGER,PARAMETER :: step = tf/10

	INTEGER,PARAMETER :: Nx = L0
	INTEGER,PARAMETER :: Ny = L0

	INTEGER,PARAMETER :: X0 = L0/2 + 1
	INTEGER,PARAMETER :: Y0 = L0/2 + 1

	INTEGER,PARAMETER :: ex(0:8) = [0, 1, 0,-1, 0, 1,-1,-1, 1]
	INTEGER,PARAMETER :: ey(0:8) = [0, 0, 1, 0,-1, 1, 1,-1,-1]
	REAL(8),PARAMETER :: Wa(0:8) = [16,4, 4, 4, 4, 1, 1, 1, 1] / 36.d0

	REAL(8),PARAMETER :: R = L0/8.d0

	REAL(8),PARAMETER :: Rhol  = 0.001d0
	REAL(8),PARAMETER :: Rhoh  = 1
	REAL(8),PARAMETER :: dRho3 = (Rhoh - Rhol)/3

	REAL(8),PARAMETER :: Sigma = 0.01d0

	REAL(8),PARAMETER :: tau = 0.3d0 + 0.5d0
	REAL(8),PARAMETER :: s8 = 1.d0/tau

	REAL(8),PARAMETER :: W    = 4
	REAL(8),PARAMETER :: Beta = 12.d0 * Sigma/W
	REAL(8),PARAMETER :: k    = 1.5d0 * Sigma*W

	REAL(8),PARAMETER :: M   = 0.02d0
	REAL(8),PARAMETER :: w_c = 1.d0/(0.5d0 + 3.d0*M)

	REAL(8) :: h(0:8,0:Nx+1,0:Ny+1), g(0:8,0:Nx+1,0:Ny+1)
	REAL(8) :: Gamma(0:8), Ga_Wa(0:8), heq(0:8), geq(0:8), hlp(0:8), eF(0:8)
	REAL(8) :: C(0:Nx+1,0:Ny+1), P(Nx,Ny), mu(Nx,Ny), DcDx(Nx,Ny), DcDy(Nx,Ny)
	REAL(8) :: Rho(Nx,Ny), Ux(Nx,Ny), Uy(Nx,Ny), ni(Nx,Ny), nj(Nx,Ny)
	REAL(8) :: tmp, Ri, Fx, Fy

END MODULE SHARE

!**********************************************************************

PROGRAM Conservative_PhaseField_Periodic_2D
	USE SHARE
	IMPLICIT NONE

	CALL Initialize_distributions

	OPEN (1, file = 'LBM.out')
	WRITE(1,*) 'Variables = X, Y, Ux, Uy, C, P'

	!=========================================================================================
	PRINT '(/A/)', '   tf    Sigma     W      M      R      tau    s8    Rhol   Rhoh     L0'
	PRINT '(I6,F8.4,F7.1,2F8.3,4F7.3,I7)', tf, Sigma, W, M, R, tau, s8, Rhol, Rhoh, L0

	PRINT '(/A6,5A12,A12/)', 't', 'C_min', 'C_max', 'Ux_max', 'Uy_max', '|U_max|', 'Mass_C'
	!=========================================================================================

	DO t = 0, tf

		IF( ISNAN(C(2,2)) )THEN
			PRINT '(/A,I5/)', '!!! THE PROGRAM DIVERGED AT t =', t
			STOP
		ELSEIF( MOD(t,step)==0 )THEN
			CALL RESULTS_Output
			PRINT '(I7,2F12.6,3E12.3,F11.2)', t, MINVAL(C), MAXVAL(C), MAXVAL(ABS(Ux)), MAXVAL(ABS(Uy)), DSQRT( MAXVAL(Ux**2+Uy**2) ), SUM(C(1:Nx,1:Ny))
		END IF

		CALL Improved_PhaseField_h_g

	END DO

	CALL CPU_TIME( tmp )
	PRINT '(/A)',' *****************************************'
	PRINT '(A,F12.1)', ' Time Elapsed:', tmp
	PRINT '( A)',' *****************************************'

END

!**********************************************************************

SUBROUTINE Initialize_distributions
	USE SHARE
	IMPLICIT NONE

	P  = 0
	Ux = 0
	Uy = 0

	DO Y = 0, Ny+1	!!1, Ny
	DO X = 0, Nx+1	!!1, Nx

		Ri = DSQRT( (X-(X0-0.5d0))**2.d0 + (Y-(Y0-0.5d0))**2.d0 )

	!	C(X,Y) = 0.5d0 + 0.5d0 * TANH(2*(R-Ri)/W)	!drop
		C(X,Y) = 0.5d0 - 0.5d0 * TANH(2*(R-Ri)/W)	!bubble

	END DO
	END DO

!!	CALL Boundary_Conditions_C( C )

	CALL Chemical_Potential

	CALL Isotropic_Gradient( C, DcDx, DcDy )

	CALL normal_FD

	DO Y = 1, Ny
	DO X = 1, Nx

		Rho(X,Y) = Rhol + C(X,Y) * (Rhoh - Rhol)

	!	P(X,Y) = P(X,Y) + C(X,Y) * Sigma/R /(Rho(X,Y)/3)	!in 2D (drop)
		P(X,Y) = P(X,Y) - C(X,Y) * Sigma/R /(Rho(X,Y)/3)	!in 2D (bubble)

		CALL Equilibrium_new( Ux(X,Y), Uy(X,Y) )

		Gamma(:) = Ga_Wa(:) + Wa(:)

		!*******************	heq

		eF(:)  = ( 1.d0 - 4.d0*(C(X,Y)-0.5d0)**2.d0 )/W * ( ex(:)*ni(X,Y) + ey(:)*nj(X,Y) )

		hlp(:) = Wa(:) * eF(:)

		h(:, X,Y) = C(X,Y)*Gamma(:) - 0.5d0 * hlp(:)

		!*******************	geq

		g(:, X,Y) = P(X,Y) * Wa(:) + Ga_Wa(:)

	END DO
	END DO

END

!**********************************************************************

SUBROUTINE Equilibrium_new( U, V )
	USE SHARE, ONLY: ex, ey, Wa, Ga_Wa
	IMPLICIT NONE
	REAL(8), INTENT(IN) :: U, V

	REAL(8) :: U2, eU(0:8)

	U2 = U*U + V*V

	eU(:) = ex(:) * U  + ey(:) * V

	Ga_Wa(:) = Wa(:) * ( eU(:)*(3.d0 + 4.5d0*eU(:)) - 1.5d0*U2 )

END

!**********************************************************************

SUBROUTINE Improved_PhaseField_h_g
	USE SHARE, ONLY: h, g
	IMPLICIT NONE

	CALL Collision_h_g

	CALL Boundary_Conditions_f( h )
	CALL Boundary_Conditions_f( g )

	CALL Streaming( h )
	CALL Streaming( g )

	CALL Macroscopic_Properties_h
	CALL Macroscopic_Properties_g

END

!**********************************************************************

SUBROUTINE Collision_h_g
	USE SHARE
	IMPLICIT NONE

	REAL(8) :: FpX, FpY, Fmx, FmY

!!	CALL Boundary_Conditions_C( C )				!no need (already called)
!!	CALL Isotropic_Gradient( C, DcDx, DcDy )	!no need (already called)

	CALL normal_FD

	DO Y = 1, Ny
	DO X = 1, Nx

		CALL Equilibrium_new( Ux(X,Y), Uy(X,Y) )

		Gamma(:) = Ga_Wa(:) + Wa(:)

	    !*******************	COLLISION (h)

		eF(:)  = ( 1.d0 - 4.d0*(C(X,Y)-0.5d0)**2.d0 )/W * ( ex(:)*ni(X,Y) + ey(:)*nj(X,Y) )

		hlp(:) = Wa(:) * eF(:)

		heq(:) = C(X,Y)*Gamma(:) - 0.5d0 * hlp(:)

		h(:, X,Y) = h(:, X,Y) * (1.d0-w_c) + heq(:) * w_c + hlp(:)

		!*******************	COLLISION (g)	******************
		!*******************	calculate the forcing terms

		FpX = - P(X,Y) * dRho3 * DcDx(X,Y)
		FpY = - P(X,Y) * dRho3 * DcDy(X,Y)

		geq(:) = P(X,Y) * Wa(:) + Ga_Wa(:)
		CALL Calculate_Viscous_Force( tau, DcDx(X,Y), DcDy(X,Y), g(:,X,Y)-geq(:), FmX, FmY )

		Fx = mu(X,Y) * DcDx(X,Y) + FpX + FmX
		Fy = mu(X,Y) * DcDy(X,Y) + FpY + FmY

		eF(:) = ex(:) * Fx + ey(:) * Fy

		hlp(:) = 3.d0 * Wa(:) * eF(:) / Rho(X,Y)

		geq(:) = P(X,Y) * Wa(:) + Ga_Wa(:) - 0.5d0 * hlp(:)

		g(:, X,Y) = g(:, X,Y) * (1-s8) + geq(:) * s8 + hlp(:)

	END DO
	END DO

END

!**********************************************************************

SUBROUTINE Boundary_Conditions_f( f )
	USE SHARE, ONLY: Nx, Ny
	IMPLICIT NONE
	REAL(8), INTENT(INOUT) :: f(0:8, 0:Nx+1,0:Ny+1)

!-- left and right boundaries
	f(:,  0  ,:) = f(:, Nx,:)	!periodic
	f(:, Nx+1,:) = f(:, 1 ,:)	!periodic

!-- bottom and top boundaries
	f(:, :, 0  ) = f(:, :,Ny)	!periodic
	f(:, :,Ny+1) = f(:, :,1 )	!periodic

END

!**********************************************************************

SUBROUTINE Boundary_Conditions_C( A )
	USE SHARE, ONLY: Nx, Ny
	IMPLICIT NONE
	REAL(8),INTENT(INOUT) :: A(0:Nx+1,0:Ny+1)

	CALL PeriodicXY_C( A )

END

!**********************************************************************

SUBROUTINE PeriodicXY_C( A )
	USE SHARE, ONLY: Nx, Ny
	IMPLICIT NONE
	REAL(8),INTENT(INOUT):: A(0:Nx+1,0:Ny+1)

!-- bottom and top boundaries
	A(:,  0 ) = A(:,Ny  )	! periodic
	A(:,Ny+1) = A(:, 1  )	! periodic

!-- left and right boundaries
	A(  0 ,:) = A(Nx  ,:)	! periodic
	A(Nx+1,:) = A( 1  ,:)	! periodic

END

!**********************************************************************

SUBROUTINE Streaming( f )
	USE SHARE, ONLY: X, Y, Nx, Ny, ex, ey
	IMPLICIT NONE
	REAL(8), INTENT(INOUT) :: f(0:8, 0:Nx+1,0:Ny+1)
	
	INTEGER :: I
	REAL(8) :: fnew(8,Nx,Ny)

	DO Y=1,Ny
	DO X=1,Nx
		DO I = 1, 8
			fnew(I,X,Y) = f(I,X-ex(I),Y-ey(I))
		END DO
	END DO
	END DO

	f(1:8,1:Nx,1:Ny) = fnew

END

!**********************************************************************

SUBROUTINE Macroscopic_Properties_h
	USE SHARE
	IMPLICIT NONE

	C(1:Nx,1:Ny) = SUM( h(:, 1:Nx,1:Ny), 1 )

	Rho = Rhol + C(1:Nx,1:Ny) * (Rhoh - Rhol)

END

!**********************************************************************

SUBROUTINE Macroscopic_Properties_g
	USE SHARE
	IMPLICIT NONE

	REAL(8) :: FpX, FpY, Fmx, FmY

	CALL Boundary_Conditions_C( C )

	CALL Chemical_Potential

	CALL Isotropic_Gradient( C, DcDx, DcDy )

	DO Y = 1, Ny
	DO X = 1, Nx

		P(X,Y) = SUM( g(:, X,Y) )

		FpX = - P(X,Y) * dRho3 * DcDx(X,Y)
		FpY = - P(X,Y) * dRho3 * DcDy(X,Y)

		CALL Equilibrium_new( Ux(X,Y), Uy(X,Y) )
		geq(:) = P(X,Y) * Wa(:) + Ga_Wa(:)
		CALL Calculate_Viscous_Force( tau, DcDx(X,Y), DcDy(X,Y), g(:,X,Y)-geq(:), FmX, FmY )

		Fx = mu(X,Y) * DcDx(X,Y) + FpX + FmX
		Fy = mu(X,Y) * DcDy(X,Y) + FpY + FmY

		Ux(X,Y) = g(1,X,Y)-g(3,X,Y)+g(5,X,Y)-g(6,X,Y)-g(7,X,Y)+g(8,X,Y) + 0.5d0*Fx/Rho(X,Y)
		Uy(X,Y) = g(2,X,Y)-g(4,X,Y)+g(5,X,Y)+g(6,X,Y)-g(7,X,Y)-g(8,X,Y) + 0.5d0*Fy/Rho(X,Y)

	END DO
	END DO

END

!**********************************************************************

SUBROUTINE Chemical_Potential
	USE SHARE, ONLY: X, Y, Nx, Ny, Beta, k, C, mu
	IMPLICIT NONE

	REAL(8) :: D2C

	DO Y = 1, Ny
	DO X = 1, Nx

		D2C = ( C(X-1,Y-1)+C(X+1,Y-1)+C(X-1,Y+1)+C(X+1,Y+1) &
			+4*(C(X  ,Y-1)+C(X-1,Y  )+C(X+1,Y  )+C(X  ,Y+1)) - 20*C(X,Y) )/6

		mu(X,Y) = 4*Beta * C(X,Y) * (C(X,Y)-1.d0) * (C(X,Y)-0.5d0) - k*D2C

	END DO
	END DO

END

!**********************************************************************

SUBROUTINE normal_FD
	USE SHARE, ONLY: X, Y, Nx, Ny, C, DcDx, DcDy, tmp, ni, nj
	IMPLICIT NONE

	DO Y = 1, Ny
	DO X = 1, Nx

		tmp = DSQRT( DcDx(X,Y)**2 + DcDy(X,Y)**2 + 1.d-32 )

		ni(X,Y) = DcDx(X,Y) / tmp
		nj(X,Y) = DcDy(X,Y) / tmp

	END DO
	END DO

END

!**********************************************************************

SUBROUTINE Isotropic_Gradient( C, DcDx, DcDy )
	USE SHARE, ONLY: X, Y, Nx, Ny
	IMPLICIT NONE
	REAL(8),INTENT(IN) :: C(0:Nx+1,0:Ny+1)
	REAL(8),INTENT(OUT):: DcDx(Nx,Ny), DcDy(Nx,Ny)

	DO Y = 1, Ny
	DO X = 1, Nx

		DcDx(X,Y) = (C(X+1,Y  ) - C(X-1,Y  ))/3 + ( C(X+1,Y-1) + C(X+1,Y+1) - C(X-1,Y-1) - C(X-1,Y+1))/12
		DcDy(X,Y) = (C(X  ,Y+1) - C(X  ,Y-1))/3 + ( C(X-1,Y+1) + C(X+1,Y+1) - C(X-1,Y-1) - C(X+1,Y-1))/12

	END DO
	END DO

END

!**********************************************************************

SUBROUTINE Calculate_Viscous_Force( tau, DcDx, DcDy, gneq, FmX, FmY )
	IMPLICIT NONE
	REAL(8),INTENT(IN)  :: tau, DcDx, DcDy, gneq(0:8)
	REAL(8),INTENT(OUT) :: FmX, FmY

	CALL Calculate_Viscous_Force_BGK( tau, DcDx, DcDy, gneq, FmX, FmY )

END

!**********************************************************************

SUBROUTINE Calculate_Viscous_Force_BGK( tau, DcDx, DcDy, gneq, FmX, FmY )
	USE SHARE, ONLY: Rhoh, Rhol
	IMPLICIT NONE
	REAL(8),INTENT(IN)  :: tau, DcDx, DcDy, gneq(0:8)
	REAL(8),INTENT(OUT) :: FmX, FmY

	REAL(8) :: sxx, sxy, syy

	CALL Calculate_Stress_Tensor_BGK( gneq(1:), sxx, sxy, syy )

	FmX = (0.5d0-tau)/tau * (sxx*DcDx+sxy*DcDy) * (Rhoh-Rhol)
	FmY = (0.5d0-tau)/tau * (sxy*DcDx+syy*DcDy) * (Rhoh-Rhol)

END

!**********************************************************************

SUBROUTINE Calculate_Stress_Tensor_BGK( gneq, sxx, sxy, syy )
	USE SHARE, ONLY: ex, ey
	IMPLICIT NONE
	REAL(8),INTENT(IN) :: gneq(1:8)
	REAL(8),INTENT(OUT):: sxx, sxy, syy

	sxx = SUM( gneq(1:) * ex(1:) * ex(1:) )
	sxy = SUM( gneq(1:) * ex(1:) * ey(1:) )
	syy = SUM( gneq(1:) * ey(1:) * ey(1:) )

END

!**********************************************************************

SUBROUTINE RESULTS_Output
	USE SHARE, ONLY: X, Y, Nx, Ny, Ux, Uy, C, T, L0, P, Rho
	IMPLICIT NONE

	WRITE(1,*) 'Zone T = "', T, '" F=Point, I=', Nx, ', J=', Ny

	DO Y = 1, Ny
	DO X = 1, Nx
		WRITE(1,'(2F8.3,4E14.6)') (X-0.5d0)/L0, (Y-0.5d0)/L0, Ux(X,Y), Uy(X,Y), C(X,Y)-0.5d0, P(X,Y)*Rho(X,Y)/3
	END DO
	END DO

END
