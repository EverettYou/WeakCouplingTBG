! ############### MODEL ###################
MODULE MODEL
	INTEGER, PARAMETER :: LK = 100  ! lattice size (for fermions)
	INTEGER, PARAMETER :: LQ = 21   ! lattice size (for bosons)
	! model parameters
	REAL, PARAMETER :: AL = 1./3.
	REAL, PARAMETER :: MU = 1.
	REAL, PARAMETER :: T = 0.02
END MODULE MODEL
! ############## MATH ####################
MODULE MATH
	COMPLEX, PARAMETER :: Z0 = (0.,0.), Z1 = (1.,0.), ZI = (0.,1.)
	REAL, PARAMETER :: PI = 4*ATAN(1.)
	REAL, PARAMETER :: B1(2) = [1.,1./SQRT(3.)], B2(2) = [-1.,1./SQRT(3.)]
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
	REAL, PARAMETER :: MAXK = MIN(3., 2./(3.*AL)) ! momentum cutoff (for fermions)
	REAL, PARAMETER :: MAXQ = 2.*MAXK    ! momentum cutoff (for bosons)
	REAL, PARAMETER :: DK = MAXK/LK ! momentum resolution (for fermions)
	REAL, PARAMETER :: DQ = MAXQ/LQ ! momentum resolution (for bosons)
	INTEGER, PARAMETER :: NK = 3*LK*(LK+1)+1 ! number of momentums (for fermions)
	INTEGER, PARAMETER :: NQ = 3*LQ*(LQ+1)+1 ! number of momentums (for bosons) 
	! grids
	REAL :: KS(NK,2), QS(NQ,2)
CONTAINS
! initialization
SUBROUTINE GRID_INIT()
	! local
	INTEGER :: I, I1, I2
	REAL :: V(2)
	
	! construct momentum grid for fermions
	I = 1 ! initialize pointer
	DO I1 = -LK, LK
		DO I2 = -LK, LK
			V = DK * (I1 * B1 + I2 * B2) ! momentum vector
			IF (ABS(V(1)) < MAXK+DK/2) THEN 
				KS(I, :) = V
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
				I = I + 1
			END IF
		END DO ! I2
	END DO ! I1
END SUBROUTINE GRID_INIT
! check if momentum is in range
PURE FUNCTION INQ(K)
	REAL, INTENT(IN) :: K(2)
	LOGICAL :: INQ
	
! 	INQ = .TRUE.
	INQ = ABS(K(1)) < MAXK+DK/2 .AND. (ABS(K(1))+SQRT(3.)*ABS(K(2)))/2 < MAXK+DK/2
END FUNCTION INQ
! end of module GRID
END MODULE GRID
! ############# PHYSICS ###################
MODULE PHYSICS
	USE MODEL
	USE GRID
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
! structural function f(q)
FUNCTION F(Q)
	USE MATH
	REAL :: F(2,2) ! f is constructed as 2x2 mat to contain s1, s2 = +/-1
	REAL, INTENT(IN) :: Q(2)
	! local
	REAL :: K0(2), K1(2), K2(2), ALPH1, BETA1, ALPH2, BETA2, EPS1, EPS2, W
	INTEGER :: IK, S1, S2, I1, I2
	
	F = 0. ! initialize F
	! momentum summation
	DO IK = 1, NK
		K0 = KS(IK,:) ! load from momentum grid
		! momentum K1, K2 are shifted from K0 by +-Q/2
		K1 = K0 - Q/2
		K2 = K0 + Q/2		
		IF (INQ(K1) .AND. INQ(K2)) THEN ! if momentums are in range
			! calculate band structure
			ALPH1 = ALPH(K1)
			BETA1 = BETA(K1)
			ALPH2 = ALPH(K2)
			BETA2 = BETA(K2)
			! loop over signs S1, S2
			DO S1 = -1, 1, 2
				I1 = (1 + S1)/2 + 1
				EPS1 = BETA1 + S1 * ALPH1 ! calculate energy
				DO S2 = -1, 1, 2
					I2 = (1 + S2)/2 + 1
					EPS2 = BETA2 + S2 * ALPH2 ! calculate energy
					IF (EPS1 /= EPS2) THEN
						W = (TANH(EPS2/(2*T))-TANH(EPS1/(2*T)))/(EPS2 - EPS1)
					ELSE
						W = 1./(2*T*COSH(EPS1/(2*T))**2)
					END IF
					! add to structural factor
					F(I1, I2) = F(I1, I2) + W
				END DO ! S2
			END DO ! S1
		END IF ! momentum in range
	END DO ! IK
	! normalize the momentum summation
	F = F / NK
END FUNCTION F
! density of state
FUNCTION DOS(W, DELTA)
	USE MATH
	REAL, INTENT(IN) :: W, DELTA
	REAL :: DOS
	! local
	INTEGER :: IK
	REAL :: H, K(2)
	COMPLEX :: ZW
	
	DOS = 0.
	ZW = W + ZI * T
	! momentum summation
	DO IK = 1, NK
		K = KS(IK, :) ! load momentum from grid
		H = BETA(K) + ALPH(K)
		DOS = DOS - 2*IMAG((ZW + H)/(ZW**2 - H**2 - DELTA**2))
	END DO ! IK
	DOS = DOS / NK
END FUNCTION DOS
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
! collect static susceptibility
SUBROUTINE COLLECT_SS()
	USE GRID
	USE PHYSICS
	USE MATHIO
	! local
	INTEGER :: IQ
	REAL :: SS(NQ,3), FMAT(2,2), Q(2)
	
	! collect data on bosonic momentum grid
	DO IQ = 1, NQ
		Q = QS(IQ,:) ! get momentum from grid
		FMAT = F(Q) ! calculate the spectral function
		SS(IQ,1) = FMAT(1,1)+FMAT(2,2)
		SS(IQ,2) = FMAT(1,2)+FMAT(2,1)
		SS(IQ,3) = FMAT(1,2)-FMAT(2,1)
	END DO ! IQ
	CALL EXPORT('SS',SS)
	CALL EXPORT('QS',QS)
END SUBROUTINE COLLECT_SS
! collect DOS
SUBROUTINE COLLECT_DOS()
	USE PHYSICS
	USE MATHIO
	REAL, PARAMETER :: W1 = -1., W2 = 0.6, DW = 0.01
	INTEGER, PARAMETER :: NW = INT((W2 - W1)/DW)
	REAL :: DAT(NW, 3), W
	INTEGER :: IW
	
	DO IW = 1, NW
		W = W1 + (IW - 1)*DW
		DAT(IW, 1) = W
		DAT(IW, 2) = DOS(W, 0.)
		DAT(IW, 3) = DOS(W, 0.2)
	END DO ! IW
	CALL EXPORT('DOS', DAT)
END SUBROUTINE COLLECT_DOS
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
	REAL :: FMAT(2,2)
	
! 	FMAT = F([0.,0.])
! 	PRINT *, FMAT
	PRINT *, DOS(-0.2,0.2)
END SUBROUTINE TEST
! end of module task
END MODULE TASK
! ############## PROGRAM ##################
PROGRAM MAIN
	USE TASK
	PRINT *, '------------ SS -------------'
	CALL INIT()
	
! 	CALL COLLECT_SS()
	CALL COLLECT_DOS() ! LK = 400, T = 0.01
! 	CALL TEST_EPS()
! 	CALL TEST()
END PROGRAM