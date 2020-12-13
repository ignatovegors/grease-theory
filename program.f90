PROGRAM greaseTheory
    IMPLICIT NONE
    INTEGER(2), PARAMETER :: io = 12
    INTEGER(4) :: ni, nj
    INTEGER(4) :: s_max
    REAL(8) :: omega, beta, phi, m, dx, dy, eps, k
    REAL(8), ALLOCATABLE :: x(:,:), y(:,:), p(:,:)


    CALL DataInput(io, omega, beta, phi, m, ni, nj, s_max, eps)

    ALLOCATE(x(0:ni, 0:nj))
    ALLOCATE(y(0:ni, 0:nj))
    ALLOCATE(p(0:ni, 0:nj))

    CALL MeshMaking(ni, nj, omega, dx, dy, x, y)

    CALL InitialConditions(ni, nj, p, x, y, omega)

    CALL Solver(ni, nj, beta, dx, dy, p, m, phi, omega, s_max, x, y, io, eps, k)

    CALL InitialConditions(ni, nj, p, x, y, omega)

    CALL SolverSecondOrder(ni, nj, beta, dx, dy, p, m, phi, omega, s_max, x, y, io, eps, k)

    DEALLOCATE(x, y, p)


END PROGRAM


SUBROUTINE DataInput(io, omega, beta, phi, m, ni, nj, s_max, eps)
    ! Takes input data from file input.txt
    IMPLICIT NONE
    INTEGER(2) :: io
    INTEGER(4) :: ni, nj, s_max
    REAL(8) :: omega, beta, phi, m, eps
    INTENT(IN) io
    INTENT(OUT) omega, beta, phi, m, ni, nj, s_max, eps

    WRITE(*,*) 'READING INPUT FILE'
    OPEN(io,FILE='INPUT.TXT')
    READ(io,*) omega
    READ(io,*) beta
    READ(io,*) phi
    READ(io,*) m
    READ(io,*) ni
    READ(io,*) nj
    READ(io,*) s_max
    READ(io,*) eps

    CLOSE(io)
    WRITE(*,*) 'SUCCESS'

    END SUBROUTINE


SUBROUTINE MeshMaking(ni, nj, omega, dx, dy, x, y)
    ! Makes mesh for numerical solution
    IMPLICIT NONE
    INTEGER(4) :: ni, nj, i, j
    REAL(8) :: omega, dx, dy
    REAL(8), DIMENSION(0:ni, 0:nj) :: x, y
    INTENT(IN) omega, ni, nj
    INTENT(OUT) dx, dy, x, y

    WRITE(*,*) 'MESH MAKING'

    dx = 1D0 / ni
    dy = omega / nj

    DO i = 0, ni
        DO j = 0, nj
            x(i,j) = i * dx
            y(i,j) = j * dy
        END DO
    END DO

    WRITE(*,*) 'SUCCESS'

    END SUBROUTINE


SUBROUTINE InitialConditions(ni, nj, p, x, y, omega)
    ! Initial conditions in the area
    IMPLICIT NONE
    INTEGER(4) :: ni, nj, i, j
    REAL(8), DIMENSION(0:ni,0:nj) :: p, x, y
    REAL(8) :: omega
    INTENT(IN) ni, nj, x, y, omega
    INTENT(OUT) p
   
    p = 0D0
    
    END SUBROUTINE


SUBROUTINE BoundaryZeroDerivative(ni, nj, p)
    ! Boundary zero derivative conditions in normal direction for 
    ! pressure on top and right sides of the area 
    IMPLICIT NONE
    INTEGER(4) :: ni, nj
    REAL(8), DIMENSION(0:ni,0:nj) :: p
    INTENT(IN) ni, nj
    INTENT(OUT) p
   
    p(1:ni, nj) = p(1:ni, nj - 1)
    p(ni, 1:nj) = p(ni - 1, 1:nj)
  
    END SUBROUTINE


SUBROUTINE BoundaryZeroDerivativeSecondOrder(ni, nj, p)
    ! Boundary zero derivative second order conditions 
    ! in normal direction for pressure on top and right 
    ! sides of the area 
    IMPLICIT NONE
    INTEGER(4) :: ni, nj
    REAL(8), DIMENSION(0:ni,0:nj) :: p
    INTENT(IN) ni, nj
    INTENT(OUT) p
   
    p(1:ni, nj) = (4D0 * p(1:ni, nj - 1) - p(1:ni, nj - 2)) / 3D0
    p(ni, 1:nj) = (4D0 * p(ni - 1, 1:nj) - p(ni - 2, 1:nj)) / 3D0
  
    END SUBROUTINE


SUBROUTINE Solver(ni, nj, beta, dx, dy, p, m, phi, omega, s_max, x, y, io, eps, k)
    ! Solver for Reynolds equation
    IMPLICIT NONE
    REAL(8), EXTERNAL :: PressureForce
    INTEGER(2) :: io
    INTEGER(4) :: s, ni, nj, s_max
    REAL(8) :: beta, dx, dy, m, phi, norm_coeff, omega, p_force, p_force_temp, p_force_a, p_force_b, k, eps
    REAL(8), DIMENSION(0:ni,0:nj) :: p, p_temp, x, y

    CALL InitialConditions(ni, nj, p_temp, x, y, omega)

    OPEN(io, FILE='RESIDUALS_FO.PLT')

    DO s = 1, s_max
        
        CALL SolverIteration(ni, nj, beta, dx, dy, p, m, phi, 1D0)

        CALL BoundaryZeroDerivative(ni, nj, p)

        p_force = PressureForce(p, ni, nj, dx, dy)
        p_force_temp = PressureForce(p_temp, ni, nj, dx, dy)

        IF (s == 1) THEN 
            norm_coeff = MAXVAL(ABS(p - p_temp))
            WRITE(*,*) 1, 'ITERATION MADE, RESIDUAL = ', 1D0, 'PRESSURE FORCE = ', p_force
        END IF

        IF (MOD(s,1000) == 0) THEN
            WRITE(*,*) s, 'ITERATIONS MADE, RESIDUAL = ', MAXVAL(ABS(p - p_temp)) / norm_coeff, 'PRESSURE FORCE = ', p_force
        END IF

        WRITE(io,*) s, MAXVAL(ABS(p - p_temp)) / norm_coeff, p_force

        IF (s == s_max) THEN
            WRITE(*,*) 'SOLUTION CONVERGED BY ITERATIONS BOUNDARY'
        END IF

        IF (MAXVAL(ABS(p - p_temp)) / norm_coeff < eps) THEN
            WRITE(*,*) 'SOLUTION CONVERGED BY DIFFERENCE BOUNDARY'
            WRITE(*,*) s, 'ITERATIONS MADE, RESIDUAL = ', MAXVAL(ABS(p - p_temp)) / norm_coeff, 'PRESSURE FORCE = ', p_force
            EXIT
        END IF        

        p_temp = p

    END DO

    CLOSE(io)

    CALL DataOutput(io, ni, nj, x, y, p, 'RES_FO.PLT')

    p_force_a = p_force

    CALL InitialConditions(ni, nj, p_temp, x, y, omega)

    DO s = 1, s_max

        CALL SolverIteration(ni, nj, beta, dx, dy, p, m, phi, 1.001D0)

        CALL BoundaryZeroDerivative(ni, nj, p)

        p_force = PressureForce(p, ni, nj, dx, dy)
        p_force_temp = PressureForce(p_temp, ni, nj, dx, dy)

        IF (MAXVAL(ABS(p - p_temp)) / norm_coeff < eps) THEN
            EXIT
        END IF        

        p_temp = p

    END DO

    p_force_b = p_force

    k = (p_force_b - p_force_a) / 1D-3

    WRITE(*,*) 'K = ', k
    WRITE(*,*)

    END SUBROUTINE


SUBROUTINE SolverSecondOrder(ni, nj, beta, dx, dy, p, m, phi, omega, s_max, x, y, io, eps, k)
    ! Second-order solver for Reynolds equation
    IMPLICIT NONE
    REAL(8), EXTERNAL :: PressureForce
    INTEGER(2) :: io
    INTEGER(4) :: s, ni, nj, s_max
    REAL(8) :: beta, dx, dy, m, phi, norm_coeff, omega, p_force, p_force_temp, p_force_a, p_force_b, k, eps
    REAL(8), DIMENSION(0:ni,0:nj) :: p, p_temp, x, y

    CALL InitialConditions(ni, nj, p_temp, x, y, omega)

    OPEN(io, FILE='RESIDUALS_SO.PLT')

    DO s = 1, s_max

        CALL SolverSecondOrderIteration(ni, nj, beta, dx, dy, p, m, phi, 1D0)

        CALL BoundaryZeroDerivativeSecondOrder(ni, nj, p)

        p_force = PressureForce(p, ni, nj, dx, dy)
        p_force_temp = PressureForce(p_temp, ni, nj, dx, dy)

        IF (s == 1) THEN 
            norm_coeff = MAXVAL(ABS(p - p_temp))
            WRITE(*,*) 1, 'ITERATION MADE, RESIDUAL = ', 1D0, 'PRESSURE FORCE = ', p_force
        END IF

        IF (MOD(s,1000) == 0) THEN
            WRITE(*,*) s, 'ITERATIONS MADE, RESIDUAL = ', MAXVAL(ABS(p - p_temp)) / norm_coeff, 'PRESSURE FORCE = ', p_force
        END IF

        WRITE(io,*) s, MAXVAL(ABS(p - p_temp)) / norm_coeff, p_force

        IF (s == s_max) THEN
            WRITE(*,*) 'SOLUTION CONVERGED BY ITERATIONS BOUNDARY'
        END IF

        IF (MAXVAL(ABS(p - p_temp)) / norm_coeff < eps) THEN
            WRITE(*,*) 'SOLUTION CONVERGED BY DIFFERENCE BOUNDARY'
            WRITE(*,*) s, 'ITERATIONS MADE, RESIDUAL = ', MAXVAL(ABS(p - p_temp)) / norm_coeff, 'PRESSURE FORCE = ', p_force
            EXIT
        END IF        

        p_temp = p

    END DO

    CLOSE(io)

    CALL DataOutput(io, ni, nj, x, y, p, 'RES_SO.PLT')

    p_force_a = p_force

    CALL InitialConditions(ni, nj, p_temp, x, y, omega)

    DO s = 1, s_max
        
        CALL SolverSecondOrderIteration(ni, nj, beta, dx, dy, p, m, phi, 1.001D0)

        CALL BoundaryZeroDerivativeSecondOrder(ni, nj, p)

        p_force = PressureForce(p, ni, nj, dx, dy)
        p_force_temp = PressureForce(p_temp, ni, nj, dx, dy)

        IF (MAXVAL(ABS(p - p_temp)) / norm_coeff < eps) THEN
            EXIT
        END IF        

        p_temp = p

    END DO

    p_force_b = p_force

    k = (p_force_b - p_force_a) / 1D-3

    WRITE(*,*) 'K = ', k
    WRITE(*,*)

    END SUBROUTINE


SUBROUTINE SolverIteration(ni, nj, beta, dx, dy, p, m, phi, h)
    ! 1 iteration of second-order solver for Reynolds equation
    IMPLICIT NONE
    LOGICAL(1), EXTERNAL :: IsInSource, IsInDitchX, IsInDitchY
    REAL(8), EXTERNAL :: SourceFunction
    INTEGER(4) :: i, j, ni, nj
    REAL(8) :: beta, dx, dy, m, phi, h
    REAL(8), DIMENSION(0:ni,0:nj) :: p
    INTENT(IN) ni, nj, beta, dx, dy, m, phi, h
    INTENT(INOUT) p

    DO i = 1, ni - 1
        DO j = 1, nj - 1
            IF (IsInSource(i, j, beta, dx, dy)) THEN 
                p(i,j) = (m * SourceFunction(p(i,j)) + phi * (p(i + 1, j) / dx + p(i, j + 1) / dy)) &
                    / (phi * (1D0 / dx + 1D0 / dy)) 
            ELSE IF (IsInDitchX(i, j, beta, dx, dy)) THEN
                p(i,j) = (phi * (p(i + 1, j) + p(i - 1, j)) / (dx * dx) + h**3 * (p(i, j - 1) + p(i, j + 1)) / dy) &
                    / (2D0 * phi / (dx * dx) + h**3 * 2D0 / dy)
            ELSE IF (IsInDitchY(i, j, beta, dx, dy)) THEN
                p(i,j) = (phi * (p(i, j + 1) + p(i, j - 1)) / (dy * dy) + h**3 * (p(i - 1, j) + p(i + 1, j)) / dx) &
                    / (2D0 * phi / (dy * dy) + h**3 * 2D0 / dx)
            ELSE 
                p(i,j) = ((p(i + 1, j) + p(i - 1, j)) / (dx * dx) + (p(i, j + 1) + p(i, j - 1)) / (dy * dy)) &
                    / (2D0 * (1D0 / (dx * dx) + 1D0 / (dy * dy)))
            END IF
        END DO
    END DO

    END SUBROUTINE


SUBROUTINE SolverSecondOrderIteration(ni, nj, beta, dx, dy, p, m, phi, h)
    ! 1 iteration of second-order solver for Reynolds equation
    IMPLICIT NONE
    LOGICAL(1), EXTERNAL :: IsInSource, IsInDitchX, IsInDitchY
    REAL(8), EXTERNAL :: SourceFunction
    INTEGER(4) :: i, j, ni, nj
    REAL(8) :: beta, dx, dy, m, phi, h
    REAL(8), DIMENSION(0:ni,0:nj) :: p
    INTENT(IN) ni, nj, beta, dx, dy, m, phi, h
    INTENT(INOUT) p

    DO i = 1, ni - 1
        DO j = 1, nj - 1
            IF (IsInSource(i, j, beta, dx, dy)) THEN 
                p(i,j) = (m * SourceFunction(p(i,j)) + phi * (4D0 * p(i + 1, j) - p(i + 2, j)) / (2D0 * dx) &
                    + phi * (4D0 * p(i, j + 1) - p(i, j + 2)) / (2D0 * dy)) / (phi * 1.5D0 * (1D0 / dx + 1D0 / dy)) 
            ELSE IF (IsInDitchX(i, j, beta, dx, dy)) THEN
                p(i,j) = (phi * (p(i + 1, j) + p(i - 1, j)) / (dx * dx) + h**3 * (4D0 * p(i, j - 1) - &
                    p(i, j - 2) + 4D0 * p(i, j + 1) - p(i, j + 2)) / (2D0 * dy)) / (2D0 * phi / (dx * dx) + h**3 * 3D0 / dy)
            ELSE IF (IsInDitchY(i, j, beta, dx, dy)) THEN
                p(i,j) = (phi * (p(i, j + 1) + p(i, j - 1)) / (dy * dy) + h**3 * (4D0 * p(i - 1, j) - &
                    p(i - 2, j) + 4D0 * p(i + 1, j) - p(i + 2, j)) / (2D0 * dx)) / (2D0 * phi / (dy * dy) + h**3 * 3D0 / dx)
            ELSE 
                p(i,j) = ((p(i + 1, j) + p(i - 1, j)) / (dx * dx) + (p(i, j + 1) + p(i, j - 1)) / (dy * dy)) &
                    / (2D0 * (1D0 / (dx * dx) + 1D0 / (dy * dy)))
            END IF
        END DO
    END DO 

    END SUBROUTINE


SUBROUTINE DataOutput(io, ni, nj, x, y, p, filename)
    ! Nodes-based results output
    IMPLICIT NONE
    INTEGER(2) :: io
    INTEGER(4) :: ni, nj
    REAL(8), DIMENSION(0:ni,0:nj) :: x, y, p
    CHARACTER(LEN=10) :: filename
    INTENT(IN) io, ni, nj, x, y, p, filename
    
    WRITE(*,*) 'RESULTS OUTPUT' 
    OPEN(io,FILE=filename)
    WRITE(io,*) 'VARIABLES = "X", "Y", "P"' 
    WRITE(io,*) 'ZONE I=', ni + 1, ', J=', nj + 1, ', DATAPACKING=BLOCK'
    WRITE(io,'(100E25.16)') x(0:ni, 0:nj)
    WRITE(io,*) 
    WRITE(io,'(100E25.16)') y(0:ni, 0:nj)
    WRITE(io,*)
    WRITE(io,'(100E25.16)') p(0:ni, 0:nj)
    CLOSE(io)
    WRITE(*,*) 'SUCCESS'

    END SUBROUTINE 


SUBROUTINE LubricantSourceLocation(beta, dx, dy, res_i, res_j)
    !   Calculation coordinates of lubricant source
    IMPLICIT NONE
    REAL(8) :: beta, dx, dy
    INTEGER(4) :: res_i, res_j
    INTENT(IN) beta, dx, dy
    INTENT(OUT) res_i, res_j

    res_i = NINT(beta / dx)
    res_j = NINT(beta / dy)

    END SUBROUTINE


LOGICAL(1) FUNCTION IsInDitchX(i, j, beta, dx, dy)
    ! Checks if node is in ditch that is along x-axis
    IMPLICIT NONE
    INTEGER(4) :: i, j, lub_src_i, lub_src_j
    REAL(8) :: beta, dx, dy

    CALL LubricantSourceLocation(beta, dx, dy, lub_src_i, lub_src_j)

    IF ((j == lub_src_j) .AND. (i > lub_src_i)) THEN
        IsInDitchX = .TRUE.
    ELSE
        IsInDitchX = .FALSE.
    END IF

    END FUNCTION


LOGICAL(1) FUNCTION IsInDitchY(i, j, beta, dx, dy)
    ! Checks if node is in ditch that is along y axis
    IMPLICIT NONE
    INTEGER(4) :: i, j, lub_src_i, lub_src_j
    REAL(8) :: beta, dx, dy

    CALL LubricantSourceLocation(beta, dx, dy, lub_src_i, lub_src_j)

    IF ((i == lub_src_i) .AND. (j > lub_src_j)) THEN
        IsInDitchY = .TRUE.
    ELSE
        IsInDitchY = .FALSE.
    END IF

    END FUNCTION


LOGICAL(1) FUNCTION IsInSource(i, j, beta, dx, dy)
    ! Checks if node is in lubricant source
    IMPLICIT NONE
    INTEGER(4) :: i, j, lub_src_i, lub_src_j
    REAL(8) :: beta, dx, dy

    CALL LubricantSourceLocation(beta, dx, dy, lub_src_i, lub_src_j)

    IF ((i == lub_src_i) .AND. (j == lub_src_j)) THEN
        IsInSource = .TRUE.
    ELSE
        IsInSource = .FALSE.
    END IF

    END FUNCTION


REAL(8) FUNCTION SourceFunction(x)
    ! Source function of lubricant source
    IMPLICIT NONE
    REAL(8) :: x

    SourceFunction = DSQRT(1D0 - x)

    END FUNCTION


REAL(8) FUNCTION PressureForce(a, ni, nj, dx, dy)
    IMPLICIT NONE
    REAL(8) :: dx, dy, sum
    REAL(8), DIMENSION(0:ni, 0:nj) :: a
    INTEGER(4) :: ni, nj, i, j

    sum = 0

    DO i = 0, ni - 1
        DO j = 0, nj - 1
            sum = sum + (a(i,j) + a(i + 1, j) + a(i, j + 1) + a(i + 1, j + 1)) / 4D0
        END DO
    END DO

    PressureForce = sum * dx * dy

    END FUNCTION
