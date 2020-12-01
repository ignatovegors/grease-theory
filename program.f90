PROGRAM greaseTheory
    IMPLICIT NONE
    LOGICAL(1), EXTERNAL :: IsInSource
    INTEGER(2), PARAMETER :: io = 12
    INTEGER(4) :: ni, nj
    INTEGER(4) :: i, j, s_max
    REAL(8) :: omega, beta, phi, m, dx, dy, eps, k
    REAL(8), ALLOCATABLE :: x(:,:), y(:,:), p(:,:)


    CALL DataInput(io, omega, beta, phi, m, ni, nj, s_max, eps)

    ALLOCATE(x(0:ni, 0:nj))
    ALLOCATE(y(0:ni, 0:nj))
    ALLOCATE(p(0:ni, 0:nj))

    CALL MeshMaking(ni, nj, omega, dx, dy, x, y)

    CALL InitialConditions(ni, nj, p, x, y, omega)

    CALL BoundaryZeroPressure(ni, nj, p)
    
    CALL BoundaryZeroDerivative(ni, nj, p) 

    ! CALL BoundaryZeroDerivativeSecondOrder(ni, nj, p) 

    CALL Solver(ni, nj, beta, dx, dy, p, m, phi, omega, s_max, x, y, io, k)

    CALL DataOutput(io, ni, nj, x, y, p)

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

    WRITE(*,*) 'INITIAL CONDITIONS APPLYING'
    
    DO i = 0, ni
        DO j = 0, nj
            ! IF (x(i,j) < y(i,j) / omega) THEN ! More accurate initial conditions
            !     p(i,j) = x(i,j)
            ! ELSE
            !     p(i,j) = y(i,j) / omega
            ! END IF

            p(i,j) = 0D0
        END DO
    END DO

    WRITE(*,*) 'SUCCESS'
    
    END SUBROUTINE


SUBROUTINE BoundaryZeroPressure(ni, nj, p)
    ! Boundary zero pressure conditions for pressure on bottom 
    ! and left sides of the area 
    IMPLICIT NONE
    INTEGER(4) :: ni, nj
    REAL(8), DIMENSION(0:ni,0:nj) :: p
    INTENT(IN) ni, nj
    INTENT(OUT) p
   
    p(0:ni, 0) = 0D0
    p(0, 0:nj) = 0D0
    
    END SUBROUTINE


SUBROUTINE BoundaryZeroDerivative(ni, nj, p)
    ! Boundary zero derivative conditions in normal direction for 
    ! pressure on top and right sides of the area 
    IMPLICIT NONE
    INTEGER(4) :: ni, nj
    REAL(8), DIMENSION(0:ni,0:nj) :: p
    INTENT(IN) ni, nj
    INTENT(OUT) p
   
    p(1:ni - 1, nj) = p(1:ni - 1, nj - 1)
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
   
    p(1:ni - 1, nj) = (4D0 * p(1:ni - 1, nj - 1) - p(1:ni - 1, nj - 2)) / 3D0
    p(ni, 1:nj) = (4D0 * p(ni - 1, 1:nj) - p(ni - 2, 1:nj)) / 3D0
  
    END SUBROUTINE


SUBROUTINE Solver(ni, nj, beta, dx, dy, p, m, phi, omega, s_max, x, y, io, k)
    ! Solver for Reynolds equation
    IMPLICIT NONE
    LOGICAL(1), EXTERNAL :: IsInSource, IsInDitchX, IsInDitchY
    REAL(8), EXTERNAL :: SourceFunction, PressureForce
    INTEGER(2) :: io
    INTEGER(4) :: i, j, s, ni, nj, s_max
    REAL(8) :: beta, dx, dy, m, phi, norm_coeff, omega, p_force, p_force_temp, p_force_a, p_force_b, k
    REAL(8), DIMENSION(0:ni,0:nj) :: p, p_temp, x, y

    CALL InitialConditions(ni, nj, p_temp, x, y, omega)

    OPEN(io, FILE='RESIDUALS.PLT')

    DO s = 1, s_max
        DO i = 1, ni - 1
            DO j = 1, nj - 1
                IF (IsInSource(i, j, beta, dx, dy)) THEN 
                    p(i,j) = (dx * dy * m * SourceFunction(p(i,j)) + phi * (dy * p(i + 1, j) + dx * p(i, j + 1))) &
                        / (phi * (dx + dy)) 
                ELSE IF (IsInDitchX(i, j, beta, dx, dy)) THEN
                    p(i,j) = (phi * dy * (p(i + 1, j) + p(i - 1, j)) + dx * dx * (p(i, j - 1) + p(i, j + 1))) &
                        / (2D0 * (phi * dy + dx * dx))
                ELSE IF (IsInDitchY(i, j, beta, dx, dy)) THEN
                    p(i,j) = (phi * dx * (p(i, j + 1) + p(i, j - 1)) + dy * dy * (p(i - 1, j) + p(i + 1, j))) &
                        / (2D0 * (phi * dx + dy * dy))
                ELSE 
                    p(i,j) = (dy * dy * (p(i + 1, j) + p(i - 1, j)) + dx * dx * (p(i, j + 1) + p(i, j - 1))) &
                        / (2D0 * (dx * dx + dy * dy))
                END IF
            END DO
        END DO

        CALL BoundaryZeroDerivative(ni, nj, p)

        p_force = PressureForce(p, ni, nj, dx, dy)
        p_force_temp = PressureForce(p_temp, ni, nj, dx, dy)

        IF (s == 1) THEN 
            norm_coeff = MAXVAL(ABS(p - p_temp))
            WRITE(*,*) 1, 'ITERATION MADE, RESIDUAL = ', 1D0, 'PRESSURE FORCE = ', p_force
        END IF

        IF (MOD(s,100) == 0) THEN
            WRITE(*,*) s, 'ITERATIONS MADE, RESIDUAL = ', MAXVAL(ABS(p - p_temp)) / norm_coeff, 'PRESSURE FORCE = ', p_force
        END IF

        WRITE(io,*) s, MAXVAL(ABS(p - p_temp)) / norm_coeff, p_force

        IF (s == s_max) THEN
            WRITE(*,*) 'SOLUTION CONVERGED BY ITERATIONS BOUNDARY'
        END IF

        IF (ABS(1D0 - p_force_temp / p_force) < 1D-6) THEN
            EXIT            
        END IF        

        p_temp = p

    END DO

    CLOSE(io)

    p_force_a = p_force

    CALL InitialConditions(ni, nj, p_temp, x, y, omega)

    DO s = 1, s_max
        DO i = 1, ni - 1
            DO j = 1, nj - 1
                IF (IsInSource(i, j, beta, dx, dy)) THEN 
                    p(i,j) = (dx * dy * m * SourceFunction(p(i,j)) + phi * (dy * p(i + 1, j) + dx * p(i, j + 1))) &
                        / (phi * (dx + dy)) 
                ELSE IF (IsInDitchX(i, j, beta, dx, dy)) THEN
                    p(i,j) = (phi * dy * (p(i + 1, j) + p(i - 1, j)) + dx * dx * 1.001D0**3 * (p(i, j - 1) + p(i, j + 1))) &
                        / (2D0 * (phi * dy + 1.001D0**3 * dx * dx))
                ELSE IF (IsInDitchY(i, j, beta, dx, dy)) THEN
                    p(i,j) = (phi * dx * (p(i, j + 1) + p(i, j - 1)) + dy * dy * 1.001D0**3 * (p(i - 1, j) + p(i + 1, j))) &
                        / (2D0 * (phi * dx + 1.001D0**3 * dy * dy))
                ELSE 
                    p(i,j) = (dy * dy * (p(i + 1, j) + p(i - 1, j)) + dx * dx * (p(i, j + 1) + p(i, j - 1))) &
                        / (2D0 * (dx * dx + dy * dy))
                END IF
            END DO
        END DO

        CALL BoundaryZeroDerivative(ni, nj, p)

        p_force = PressureForce(p, ni, nj, dx, dy)
        p_force_temp = PressureForce(p_temp, ni, nj, dx, dy)

        IF (ABS(1D0 - p_force_temp / p_force) < 1D-6) THEN
            EXIT            
        END IF        

        p_temp = p

    END DO

    p_force_b = p_force

    k = (p_force_b - p_force_a) / 1D-3

    WRITE(*,*) 'K = ', k

    END SUBROUTINE


SUBROUTINE DataOutput(io, ni, nj, x, y, p)
    ! Nodes-based results output
    IMPLICIT NONE
    INTEGER(2) :: io
    INTEGER(4) :: ni, nj
    REAL(8), DIMENSION(0:ni,0:nj) :: x, y, p
    INTENT(IN) io, ni, nj, x, y, p
    
    WRITE(*,*) 'RESULTS OUTPUT' 
    OPEN(io,FILE='RES.PLT')
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

    res_i = NINT(beta / dx)
    res_j = NINT(beta / dy)

    END SUBROUTINE


LOGICAL(1) FUNCTION IsInDitchX(i, j, beta, dx, dy)
    ! Checks if node is in ditch that is along x axis
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


LOGICAL(1) FUNCTION ConvergenceCheck(a, b, eps)
    ! Convergence check for iteration process 
    IMPLICIT NONE
    REAL(8) :: a, b, eps

    ConvergenceCheck = (ABS(a - b) < eps)

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
