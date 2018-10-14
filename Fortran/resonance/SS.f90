! ############### MODEL ###################
MODULE MODEL
	INTEGER, PARAMETER :: LK = 60  ! lattice size (for fermions)
	INTEGER, PARAMETER :: LQ = 30   ! lattice size (for bosons)
	INTEGER, PARAMETER :: NW = 30
	REAL, PARAMETER  :: MAXW = 1.5
	REAL, PARAMETER  :: DW = MAXW/NW
	! model parameters
	REAL, PARAMETER :: AL = 0.001/3.
	REAL, PARAMETER :: MU = 1.
	REAL, PARAMETER :: T = 0.02
END MODULE MODEL
! ############## MATH ####################
MODULE MATH
	COMPLEX, PARAMETER :: Z0 = (0.,0.), Z1 = (1.,0.), ZI = (0.,1.)
	REAL, PARAMETER :: PI = 4*ATAN(1.)
	! reciprocal vectors and primitive vectors, Ai.Bj = delta_ij
	REAL, PARAMETER :: B1(2) = [1.,1./SQRT(3.)], B2(2) = [-1.,1./SQRT(3.)]
	REAL, PARAMETER :: A1(2) = [1.,SQRT(3.)]/2., A2(2) = [-1.,SQRT(3.)]/2.
CONTAINS
! norm of a momentum vector
PURE FUNCTION NORM(K)
	REAL, INTENT(IN) :: K(2)
	REAL :: NORM
	
	NORM = SQRT(K(1)**2 + K(2)**2)
END FUNCTION NORM
! end of module MATH
END MODULE MATH
! ############## GRID ####################
MODULE GRID
	USE MODEL
	USE MATH
	! parameters
	REAL, PARAMETER :: MAXK = MIN(2., 2./(3.*AL)) ! momentum cutoff (for fermions)
	REAL, PARAMETER :: MAXQ = 2.*MAXK    ! momentum cutoff (for bosons)
	REAL, PARAMETER :: DK = MAXK/LK ! momentum resolution (for fermions)
	REAL, PARAMETER :: DQ = MAXQ/LQ ! momentum resolution (for bosons)
	INTEGER, PARAMETER :: NK = 3*LK*(LK+1)+1 ! number of momentums (for fermions)
	INTEGER, PARAMETER :: NQ = 3*LQ*(LQ+1)+1 ! number of momentums (for bosons) 
	! grids
	REAL :: KS(NK,2), QS(NQ,2)
	! index maps
	INTEGER :: KMAP(-LK:LK, -LK:LK), QMAP(-LQ:LQ, -LQ:LQ)
CONTAINS
! initialization
SUBROUTINE GRID_INIT()
	! local
	INTEGER :: I, I1, I2
	REAL :: V(2)
	
	! clear index maps
	KMAP = 0
	QMAP = 0
	! construct momentum grid for fermions
	I = 1 ! initialize pointer
	DO I1 = -LK, LK
		DO I2 = -LK, LK
			V = DK * (I1 * B1 + I2 * B2) ! momentum vector
			IF (ABS(V(1)) < MAXK+DK/2) THEN
				KS(I, :) = V
				KMAP(I1, I2) = I
				I = I + 1
			END IF
		END DO ! I2
	END DO ! I1
	! construct momentum grid for bosons
	I = 1 ! initialize pointer
	DO I1 = -LQ, LQ
		DO I2 = -LQ, LQ
			V = DQ * (I1 * B1 + I2 * B2) ! momentum vector
			IF (ABS(V(1)) < MAXQ+DQ/2) THEN 
				QS(I, :) = V
				QMAP(I1, I2) = I
				I = I + 1
			END IF
		END DO ! I2
	END DO ! I1
END SUBROUTINE GRID_INIT
! K grid index function
FUNCTION KIND(K)
	REAL, INTENT(IN) :: K(2)
	INTEGER :: KIND
	! local
	INTEGER :: I1, I2
	
	I1 = NINT(DOT_PRODUCT(A1, K) / DK)
	I2 = NINT(DOT_PRODUCT(A2, K) / DK)
	IF (ABS(I1) <= LK .AND. ABS(I2) <= LK) THEN
		KIND = KMAP(I1, I2)
	ELSE
		KIND = 0
	END IF
END FUNCTION KIND
! Q grid index function
FUNCTION QIND(Q)
	REAL, INTENT(IN) :: Q(2)
	INTEGER :: QIND
	! local
	INTEGER :: I1, I2
	
	I1 = NINT(DOT_PRODUCT(A1, Q) / DQ)
	I2 = NINT(DOT_PRODUCT(A2, Q) / DQ)
	IF (ABS(I1) <= LQ .AND. ABS(I2) <= LQ) THEN
		QIND = QMAP(I1, I2)
	ELSE
		QIND = 0
	END IF
END FUNCTION QIND
! end of module GRID
END MODULE GRID
! ############# PHYSICS ###################
MODULE PHYSICS
	USE MODEL
	USE GRID
	! parameters
	INTEGER, PARAMETER :: NH = 4 ! single-particle Hilbert space dim
	REAL :: ES(NK, NH), US(NK, NH, NH)
	REAL :: CHI0(4, NQ, 0:NW)
CONTAINS
! band structure function beta(k)
PURE FUNCTION BETA(K)
	REAL :: BETA
	REAL, INTENT(IN) :: K(2) ! momentum vector
	
	BETA = K(1)**2 + K(2)**2 - MU
END FUNCTION BETA
! band structure function alpha(k)
PURE FUNCTION ALPH(K)
	REAL :: ALPH 
	REAL, INTENT(IN) :: K(2) ! momentum vector
	
	ALPH = AL * (K(1)**3 - 3*K(1)*K(2)**2)
END FUNCTION ALPH
! construct and diagonalize Hamiltonian
SUBROUTINE DIAG(D)
	REAL, INTENT(IN) :: D ! pairing gap
	! local
	INTEGER, PARAMETER :: LWORK = 64*NH
	REAL :: K(2), A, B, H(NH, NH), W(NH), WORK(LWORK)
	INTEGER :: IK, INFO
	
	DO IK = 1, NK
		K = KS(IK, :)
		A = ALPH(K)
		B = BETA(K)
		H = 0.
		H(1, 1) = A
		H(1, 2) = B
		H(1, 3) = D
		H(2, 2) = A
		H(2, 1) = B
		H(2, 4) = -D
		H(3, 3) = -A
		H(3, 4) = B
		H(3, 1) = D
		H(4, 4) = -A
		H(4, 3) = B
		H(4, 2) = -D
        ! call LAPACK to diagonalize
        CALL DSYEV('V', 'U', NH, H, NH, W, WORK, LWORK, INFO)
        IF (INFO /= 0) THEN ! error handling
			IF (INFO > 0) THEN
				WRITE (*,'(A,I12)') 'the QR algorithm failed.', INFO
			ELSE
				WRITE (*,'(A,I12)') 'DSYEV error ', INFO
			END IF
		END IF
        ES(IK, :) = W ! keep eigen values
        US(IK, :, :) = H ! keep eigen vectors
	END DO ! IK
END SUBROUTINE DIAG
! I0 expectation
FUNCTION I0(V1, V2)
	REAL, INTENT(IN) :: V1(NH), V2(NH)
	REAL :: I0
	
	I0 = V1(1)*V2(2) + V1(2)*V2(1) + V1(3)*V2(4) + V1(4)*V2(3)
END FUNCTION I0
! Ix expectation
FUNCTION IX(V1, V2)
	REAL, INTENT(IN) :: V1(NH), V2(NH)
	REAL :: IX
	
	IX = V1(1)*V2(4) + V1(2)*V2(3) + V1(3)*V2(2) + V1(4)*V2(1)
END FUNCTION IX
! Iz expectation
FUNCTION IZ(V1, V2)
	REAL, INTENT(IN) :: V1(NH), V2(NH)
	REAL :: IZ
	
	IZ = V1(1)*V2(2) + V1(2)*V2(1) - V1(3)*V2(4) - V1(4)*V2(3)
END FUNCTION IZ
! i Iy expectation
FUNCTION IY(V1, V2)
	REAL, INTENT(IN) :: V1(NH), V2(NH)
	REAL :: IY
	
	IY = V1(1)*V2(3) + V1(2)*V2(4) - V1(3)*V2(1) - V1(4)*V2(2)
END FUNCTION IY
! calculate CHI0
SUBROUTINE MK_CHI0()
	! local
	INTEGER :: IK1, IK2, IQ, IH1, IH2, IW
	REAL :: K1(2), K2(2), Q(2), W, FAC, I, X, Y, Z
	REAL :: E1(NH), U1(NH,NH), E2(NH), U2(NH,NH)
	
	CHI0 = 0. ! initialization
	DO IK1 = 1, NK
		K1 = KS(IK1, :)
		! get E1, U1
		E1 = ES(IK1, :)
		U1 = US(IK1, :, :)
		DO IK2 = 1, NK
			K2 = KS(IK2, :)
			! check momentum transfer in range
			Q = K2 - K1
			IQ = QIND(Q)
			IF (IQ == 0) CYCLE
			! get E2, U2
			E2 = ES(IK2, :)
			U2 = US(IK2, :, :)
			DO IH1 = 1, NH
				IF (E1(IH1) > 0.) CYCLE
				DO IH2 = 1, NH
					IF (E2(IH2) < 0.) CYCLE
					W = E2(IH2) - E1(IH1)
					IW = NINT(W/DW)
					IF (IW > NW) CYCLE
					I = I0(U1(:,IH1),U2(:,IH2))
					X = IX(U1(:,IH1),U2(:,IH2))
					Y = IY(U1(:,IH1),U2(:,IH2))
					Z = IZ(U1(:,IH1),U2(:,IH2))
					CHI0(1, IQ, IW) = CHI0(1, IQ, IW) + I * I
					CHI0(2, IQ, IW) = CHI0(2, IQ, IW) + X * X
					CHI0(3, IQ, IW) = CHI0(3, IQ, IW) + X * Y
					CHI0(4, IQ, IW) = CHI0(4, IQ, IW) + Z * Z
				END DO ! IH2
			END DO ! IH1			
		END DO ! IK2
	END DO ! IK1
	CHI0 = CHI0 / (NH*NK)
END SUBROUTINE MK_CHI0
! end of module PHYSICS
END MODULE PHYSICS
! ############### TASK ####################
MODULE TASK
CONTAINS
! ----------- data ------------------
! initialization
SUBROUTINE INIT()
	USE GRID
	
	CALL GRID_INIT() ! initialize momentum grids
END SUBROUTINE INIT
! collect CHI0
SUBROUTINE COLLECT_CHI0()
	USE PHYSICS
	USE MATHIO
	! local
	REAL :: WS(0:NW)
	INTEGER :: IW
	
	CALL DIAG(0.)
	CALL MK_CHI0()
	CALL EXPORT('CHI0', CHI0)
	CALL EXPORT('QS', QS)
	DO IW = 0, NW
		WS(IW) = IW * DW
	END DO
	CALL EXPORT('WS', WS)
END SUBROUTINE COLLECT_CHI0
! ----------- tests -----------------
! check dispersion
SUBROUTINE TEST_EPS()
	USE GRID
	USE PHYSICS
	USE MATHIO
	INTEGER :: IK
	REAL :: EPS(NK, 2), K(2)
	
	DO IK = 1, NK
		K = KS(IK,:)
		EPS(IK,1) = BETA(K) + ALPH(K)
		EPS(IK,2) = BETA(K) - ALPH(K)		
	END DO
	CALL EXPORT('EPS',EPS)
	CALL EXPORT('KS', KS)
END SUBROUTINE TEST_EPS
! test
SUBROUTINE TEST()
	USE MATH
	USE GRID
	USE PHYSICS
	USE MATHIO
	REAL :: K(2)
	
	CALL DIAG(0.)
	CALL EXPORT('ES', ES)
	CALL EXPORT('KS', KS)
END SUBROUTINE TEST
! end of module task
END MODULE TASK
! ############## PROGRAM ##################
PROGRAM MAIN
	USE TASK
	PRINT *, '------------ SS -------------'
	CALL INIT()
	
	CALL COLLECT_CHI0()
! 	CALL TEST_EPS()
! 	CALL TEST()
END PROGRAM