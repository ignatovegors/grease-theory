PROGRAM greaseTheory
    IMPLICIT NONE
    INTEGER(2), PARAMETER :: io = 12
    INTEGER(4) :: ni, nj
    INTEGER(4) :: i, j, s_max
    REAL(8) :: omega, beta, phi, m, dx, dy, eps
    REAL(8), ALLOCATABLE :: x(:,:), y(:,:), p(:,:)


    CALL DataInput(io, omega, beta, phi, m, ni, nj, s_max, eps)

    ALLOCATE(x(0:ni, 0:nj))
    ALLOCATE(y(0:ni, 0:nj))
    ALLOCATE(p(0:ni, 0:nj))

    CALL MeshMaking(ni, nj, omega, dx, dy, x, y)

    CALL InitialConditions(ni, nj, p, x, y, omega)

    CALL BoundaryZeroPressure(ni, nj, p)
    
    CALL BoundaryZeroDerivative(ni, nj, p) 

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
            IF (x(i,j) < y(i,j) / omega) THEN
                p(i,j) = x(i,j)
            ELSE
                p(i,j) = y(i,j) / omega
            END IF
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


! SUBROUTINE SolverPrandtl(ni, nj, s_max, dx, dy, nu, eps, u_0, u, v)
!     ! Solver for Prandtl (Simplified Navier-Stokes) equations system
!     IMPLICIT NONE
!     LOGICAL(1), EXTERNAL :: ConvergenceCheckPrandtl
!     INTEGER(4) :: i, j, s, ni, nj, s_max
!     REAL(8) :: dx, dy, nu, eps, u_0
!     REAL(8), DIMENSION(ni,nj) :: u, v
!     REAL(8), DIMENSION(nj) :: u_temp, v_temp, a, b, c, d
!     INTENT(IN) :: ni, nj, s_max, dx, dy, nu, eps, u_0
!     INTENT(OUT) :: u, v

!     WRITE(*,*) 'SOLVING EQUATIONS (PRANDTL)'

!     DO i = 2, ni
            
!         u_temp = u(i - 1, :)
!         v_temp = v(i - 1, :)

!         DO s = 1, s_max

!             a(1) = 0D0
!             b(1) = 1D0
!             c(1) = 0D0
!             d(1) = 0D0

!             DO j = 2, nj - 1
!                 a(j) = - v_temp(j - 1) / (2D0 * dy) - nu / dy**2
!                 b(j) = u_temp(j) / dx + 2D0 * nu / dy**2
!                 c(j) = v_temp(j + 1) / (2D0 * dy) - nu / dy**2
!                 d(j) = u(i - 1, j)**2D0 / dx
!             END DO

!             a(nj) = 0D0
!             b(nj) = 1D0
!             c(nj) = 0D0
!             d(nj) = u_0

!             CALL ThomasAlgorithm(nj, a, b, c, d, u(i, :))

!             DO j = 2, nj
!                 v(i,j) = v(i, j - 1) - dy / (2 * dx) * (u(i,j) - u(i - 1, j) + u(i, j - 1) - u(i - 1, j - 1))
!             END DO    
            
!             IF ((ConvergenceCheckPrandtl(u(i, :), u_temp, nj, eps)) .AND. (ConvergenceCheckPrandtl(v(i, :), v_temp, nj, eps))) THEN
!                 WRITE(*,*) 'SOLUTION CONVERGED BY RESIDUALS, NODE №', I, ', s = ', s
!                 EXIT 
!             END IF

!             u_temp = u(i, :)
!             v_temp = v(i, :)

!             IF (s == s_max) THEN
!                 WRITE(*,*) 'SOLUTION CONVERGED BY ITERATIONS BOUNDARY, NODE №', I
!             END IF

!         END DO

!     END DO

!     WRITE(*,*) 'SUCCESS'

!     END SUBROUTINE


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


LOGICAL(1) FUNCTION IsInDitch(i, j, beta, dx, dy)
    ! Checks if node is in ditch
    IMPLICIT NONE
    INTEGER(4) :: i, j, lub_src_i, lub_src_j
    REAL(8) :: beta, dx, dy

    CALL LubricantSourceLocation(beta, dx, dy, lub_src_i, lub_src_j)

    IF ((i == lub_src_i) .AND. (j > lub_src_j)) THEN
        IsInDitch = .TRUE.
    ELSE IF ((j == lub_src_j) .AND. (i > lub_src_i)) THEN
        IsInDitch = .TRUE.
    ELSE
        IsInDitch = .FALSE.
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
