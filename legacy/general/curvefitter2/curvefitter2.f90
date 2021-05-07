
          PROGRAM curvefitter2
            INTEGER :: i
            CHARACTER(len=24) :: arg_k
            real(8) RMS
            character(254) file,file_out
            integer k,n_exp
            REAL(8), ALLOCATABLE, DIMENSION(:) :: P
            
           file_out='curvefitter2.temp'
            
            CALL getarg(1, file)
!             WRITE (*,*) file
            CALL getarg(2, arg_k)
            read(arg_k, '(i10)') k
!             WRITE (*,*) k
            allocate(P(0:k))
            
            call Curve_Fitting(file, k,n_exp,RMS, P, file_out)
            
            WRITE (0,*) 'RMS = ', RMS
            
          END PROGRAM

        !---------------------
    subroutine Curve_Fitting(file, k, n_exp,RMS, P, file_out)
    ! MAXIMUM DEGREE OF THE POLYNOM IS 6!!!!
    implicit none
    character(254) file, file_out
    integer k
    real(8) RMS, P(0:k),mult !output
    real(8) A(0:k,0:k),sumi(0:2*k),qq,rr,Y(0:k),AA(k+1,k+1),AA_inv(k+1,k+1), A_inv(0:k,0:k),prov(0:k,0:k),&
    Y_calc, X1a, X1b,aaa,bbb,d(0:6),PP(0:6)
    
    REAL(8), ALLOCATABLE, DIMENSION(:) :: X,B
    integer i,j,i_exp,n_exp,i9,j9,ind_err
    
    !reading from file
    !determination of number of experimental points
    
    if(k.gt.6)then
    write(*,*)'MAXIMUM DEGREE OF POLYNOM MUST BE 6'
    stop
    end if
    
    
    open(1,file = file)
    do i = 1,100000000
        read(1,*,end = 500)qq,rr
    end do
    
    open(13,file = 'Unit_matrix?.txt')
           
500 n_exp = i-1
    close(1)
    
    allocate(X(N_exp))
    allocate(B(N_exp))
    
     open(1,file = file)
    do i = 1,n_exp
        read(1,*)X(i),B(i)
        end do
        close(1)
        
        !goto 1512
        
        X1a = 0
        X1b = 1.
        
               
        bbb = (X1a-X1b)/(X(1)-X(n_exp))
        aaa = X1a - bbb*X(1)
        aaa = X1b-bbb*X(N_exp)
        
        !write(*,*)'aaa = ',aaa
        !write(*,*)'bbb = ',bbb
        !mult = X(n_exp/2)
    X = aaa+bbb*X
    !1512    mult = (X(n_exp) - x(1))
    !write(*,*)'mult_1 = ',mult
    !X = X/mult  !(X(n_exp) - x(1))
    
    !do i = 1,n_exp
    !write(13,*)x(i)
    !end do
    !forming matrix to be inverted
    
    sumi = 0.
    Y = 0.
    
    do i = 0,2*k     
        do j = 1,n_exp
      sumi(i) = sumi(i)+X(j)**(i)  
        end do      !j
        !write(*,*)sumi(i)
    end do
    
    
    
    do i = 0,k     
        do j = 1,n_exp
            Y(i) = Y(i)+B(j)*X(j)**(i)
        end do      !j
        end do
    !------
    
    do i = 0,k
        do j = 0,k
      A(i,j) = sumi(i+j) 
       AA(i+1,j+1) = A(i,j) 
        end do
    end do
    
   ! write(*,*)'matrix A'
    !write(*,'(1x,4(1x,f12.6))')((A(i9,j9),i9 = 0,k), j9 = 0,k)
        
    !call MATR_INVERT(AA, AA_inv, k+1, ind_err)
    !call inverse(aa,aa_inv,k+1)
    call matrixinv(aa,aa_inv,k+1)   !GOOD
    !write(*,*)'ind_err = ',ind_err
    
    !write(*,*)'matrix A_inv'
    !write(*,'(1x,4(1x,f12.6))')((AA_inv(i9,j9),i9 = 1,k+1), j9 = 1,k+1)
    
    do i = 0,k
        do j = 0,k
         A_inv(i,j) = AA_inv(i+1,j+1)   
            
        end do
    end do
    
    prov = 0.
    
    do i = 0,k
        do j = 0,k
        do i9 = 0,k
            prov(i,j) = prov(i,j)+A_inv(i,i9)*A(i9,j)
            
        end do
        end do
    end do
    
!     write(*,*)'Unit matrix?'
!     write(*,'(1x,7(1x,f12.6))')((prov(i9,j9),i9 = 0,k), j9 = 0,k)
    
    write(13,*)'Unit matrix?'
    write(13,'(1x,7(1x,e12.6))')((prov(i9,j9),i9 = 0,k), j9 = 0,k)
   
    !pause
    !write(*,*)'P'
!----------------------------- 
write(13,*)'-------------------' 
   ! X = (x-aaa)/bbb
     !do i = 1,n_exp
    !write(13,*)x(i)
    !end do
    !pause
    
    Y=0.
    
     do i = 0,k     
        do j = 1,n_exp
           ! Y(i) = Y(i)+B(j)*(aaa+bbb*X(j))**(i)
                        Y(i) = Y(i)+B(j)*X(j)**(i)

        end do      !j
        end do
    
    P = 0.
    do i = 0,k
        do j = 0,k
     p(i) = p(i)+A_inv(i,j)*Y(j)
        end do  
     !d0 := a^6*p6+a^5*p5+a^4*p4+a^3*p3+a^2*p2+a*p1   
      
        
       ! write(*,*)P(i)
    end do
    
    pp = 0.
    do i = 0,k
    pp(i) = p(i)
    end do
    
     !d0 := a^6*p6+a^5*p5+a^4*p4+a^3*p3+a^2*p2+a*p1+p0   
    d(0)=aaa**6*pp(6)+aaa**5*pp(5)+aaa**4*pp(4)+aaa**3*pp(3)+aaa**2*pp(2)+aaa*pp(1)+pp(0)
    
    !d1 := 6*a^5*b*p6+5*a^4*b*p5+4*a^3*b*p4+3*a^2*b*p3+2*a*b*p2+b*p1
    d(1)= 6*aaa**5*bbb*pp(6)+5*aaa**4*bbb*pp(5)+4*aaa**3*bbb*pp(4)+3*aaa**2*bbb*pp(3)+2*aaa*bbb*pp(2)+bbb*pp(1)
    
    !d2 := 15*a^4*b^2*p6+10*a^3*b^2*p5+6*a^2*b^2*p4+3*a*b^2*p3+b^2*p2
    d(2)= 15*aaa**4*bbb**2*pp(6)+10*aaa**3*bbb**2*pp(5)+6*aaa**2*bbb**2*pp(4)+3*aaa*bbb**2*pp(3)+bbb**2*pp(2)
    
    !d3 := 20*a^3*b^3*p6+10*a^2*b^3*p5+4*a*b^3*p4+b^3*p3
    d(3)= 20*aaa**3*bbb**3*pp(6)+10*aaa**2*bbb**3*pp(5)+4*aaa*bbb**3*pp(4)+bbb**3*pp(3)
    
    !d4 := 15*a^2*b^4*p6+5*a*b^4*p5+b^4*p4
    d(4)= 15*aaa**2*bbb**4*pp(6)+5*aaa*bbb**4*pp(5)+bbb**4*pp(4)
    
    !d5 := 6*a*b^5*p6+b^5*p5
    d(5)= 6*aaa*bbb**5*pp(6)+bbb**5*pp(5)
    
    !d6 := p6*b^6
    d(6)= pp(6)*bbb**6
    
    do i = 0,k
    p(i) = d(i)
    write(*,*)d(i)
    end do
    
    !pause
    open(300,file = file_out)
    !-----RMS-----
    
    X = (x-aaa)/bbb
    RMS = 0.
   
    do j = 1,n_exp
        Y_calc= 0.
        do i = 0,k
            Y_calc = Y_calc+d(i)*X(j)**(i)
                 end do  !i
        write(300,*)X(j),Y_calc
        RMS = RMS+(Y_calc-B(j))**2/n_exp
    end do
    
    RMS = dsqrt(RMS)
    write(300,*)'RMS = ',RMS
!     write(*,*)'RMS = ',RMS
    do i = 0,k
     write(300,*)p(i)
     end do
     
    close(300)
    
    return
    end
    !-----------------
    SUBROUTINE MATR_INVERT(matrix, inverse, n, ind_err)
	IMPLICIT NONE
	!Declarations
	INTEGER, INTENT(IN) :: n
	INTEGER, INTENT(OUT) :: ind_err  !Return error status. -1 for error, 0 for normal
	REAL(8), INTENT(IN), DIMENSION(n,n) :: matrix  !Input matrix
	REAL(8), INTENT(OUT), DIMENSION(n,n) :: inverse !Inverted matrix
	
	LOGICAL :: FLAG = .TRUE.
	INTEGER :: i, j, k, l,errorflag
	REAL(8) :: m
	REAL(8), DIMENSION(n,2*n) :: augmatrix !augmented matrix
	
	!Augment input matrix with an identity matrix
	DO i = 1, n
		DO j = 1, 2*n
			IF (j <= n ) THEN
				augmatrix(i,j) = matrix(i,j)
			ELSE IF ((i+n) == j) THEN
				augmatrix(i,j) = 1
			Else
				augmatrix(i,j) = 0
			ENDIF
		END DO
	END DO
	
	!Reduce augmented matrix to upper traingular form
	DO k =1, n-1
		IF (augmatrix(k,k) == 0) THEN
			FLAG = .FALSE.
			DO i = k+1, n
				IF (augmatrix(i,k) /= 0) THEN
					DO j = 1,2*n
						augmatrix(k,j) = augmatrix(k,j)+augmatrix(i,j)
					END DO
					FLAG = .TRUE.
					EXIT
				ENDIF
				IF (FLAG .EQV. .FALSE.) THEN
					PRINT*, "Matrix is non - invertible"
					inverse = 0
					errorflag = -1
					return
				ENDIF
			END DO
		ENDIF
		DO j = k+1, n			
			m = augmatrix(j,k)/augmatrix(k,k)
			DO i = k, 2*n
				augmatrix(j,i) = augmatrix(j,i) - m*augmatrix(k,i)
			END DO
		END DO
	END DO
	
	!Test for invertibility
	DO i = 1, n
		IF (augmatrix(i,i) == 0) THEN
			PRINT*, "Matrix is non - invertible"
			inverse = 0
			errorflag = -1
			return
		ENDIF
	END DO
	
	!Make diagonal elements as 1
	DO i = 1 , n
		m = augmatrix(i,i)
		DO j = i , (2 * n)				
			   augmatrix(i,j) = (augmatrix(i,j) / m)
		END DO
	END DO
	
	!Reduced right side half of augmented matrix to identity matrix
	DO k = n-1, 1, -1
		DO i =1, k
		m = augmatrix(i,k+1)
			DO j = k, (2*n)
				augmatrix(i,j) = augmatrix(i,j) -augmatrix(k+1,j) * m
			END DO
		END DO
	END DO				
	
	!store answer
	DO i =1, n
		DO j = 1, n
			inverse(i,j) = augmatrix(i,j+n)
		END DO
	END DO
	errorflag = 0
    ind_err = errorflag
END SUBROUTINE MATR_INVERT
!------------OTHER----------

  subroutine inverse(a,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
implicit none 
integer n
double precision a(n,n), c(n,n)
double precision L(n,n), U(n,n), b(n), d(n), x(n)
double precision coeff
integer i, j, k

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
L=0.0
U=0.0
b=0.0

! step 1: forward elimination
do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
end do

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,n
  L(i,i) = 1.0
end do
! U matrix is the upper triangular part of A
do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
end do

! Step 3: compute columns of the inverse matrix C
do k=1,n
  b(k)=1.0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
end do
end subroutine inverse
!-----------------
 subroutine matrixinv(a,b,n)
 ! subroutine to calculate the inverse of a matrix using Gauss-Jordan elimination
 ! the inverse of matrix a(n,n) is calculated and stored in the matrix b(n,n)
 integer :: i,j,k,l,m,n,irow
 real(8):: big,dum
 real(8) a(n,n),b(n,n)

 !build the identity matrix
!PSPad editor 4.5.2 (2241) www.pspad.com 9/11/2008 2:20:27 PM cad
!matrix-inverse.f90 Page:2/3
!D:\Gabriel\Classes\Fall08\Finite-Element-ME549\Homeworks\Matrix-inverse\ Last modification: 9/11/2008 2:20:22 PM
 do i = 1,n
do j = 1,n
 b(i,j) = 0.0
 end do
 b(i,i) = 1.0
 end do

 do i = 1,n ! this is the big loop over all the columns of a(n,n)
 ! in case the entry a(i,i) is zero, we need to find a good pivot; this pivot
 ! is chosen as the largest value on the column i from a(j,i) with j = 1,n
 big = a(i,i)
 do j = i,n
 if (a(j,i).gt.big) then
 big = a(j,i)
 irow = j
 end if
 end do
 ! interchange lines i with irow for both a() and b() matrices
 if (big.gt.a(i,i)) then
 do k = 1,n
 dum = a(i,k) ! matrix a()
 a(i,k) = a(irow,k)
 a(irow,k) = dum
 dum = b(i,k) ! matrix b()
 b(i,k) = b(irow,k)
 b(irow,k) = dum
 end do
 end if
 ! divide all entries in line i from a(i,j) by the value a(i,i);
 ! same operation for the identity matrix
dum = a(i,i)
 do j = 1,n
 a(i,j) = a(i,j)/dum
 b(i,j) = b(i,j)/dum
 end do
 ! make zero all entries in the column a(j,i); same operation for indent()
 do j = i+1,n
 dum = a(j,i)
 do k = 1,n
 a(j,k) = a(j,k) - dum*a(i,k)
 b(j,k) = b(j,k) - dum*b(i,k)
 end do
 end do
 end do

 ! substract appropiate multiple of row j from row j-1
 do i = 1,n-1
 do j = i+1,n
 dum = a(i,j)
 do l = 1,n
 a(i,l) = a(i,l)-dum*a(j,l)
 b(i,l) = b(i,l)-dum*b(j,l)
 end do
 end do
 end do

 end
 
 !-------------------------------
 !!*******************************************************
!*    LU decomposition routines used by test_lu.f90    *
!*                                                     *
!*                 F90 version by J-P Moreau, Paris    *
!*                        (www.jpmoreau.fr)            *
!* --------------------------------------------------- *
!* Reference:                                          *
!*                                                     *
!* "Numerical Recipes By W.H. Press, B. P. Flannery,   *
!*  S.A. Teukolsky and W.T. Vetterling, Cambridge      *
!*  University Press, 1986" [BIBLI 08].                *
!*                                                     * 
!*******************************************************
MODULE LU

CONTAINS

!  ***************************************************************
!  * Given an N x N matrix A, this routine replaces it by the LU *
!  * decomposition of a rowwise permutation of itself. A and N   *
!  * are input. INDX is an output vector which records the row   *
!  * permutation effected by the partial pivoting; D is output   *
!  * as -1 or 1, depending on whether the number of row inter-   *
!  * changes was even or odd, respectively. This routine is used *
!  * in combination with LUBKSB to solve linear equations or to  *
!  * invert a matrix. Return code is 1, if matrix is singular.   *
!  ***************************************************************
 Subroutine LUDCMP(A,N,INDX,D,CODE)
 PARAMETER(NMAX=100,TINY=1.5D-16)
 REAL*8  AMAX,DUM, SUM, A(N,N),VV(NMAX)
 INTEGER CODE, D, INDX(N)

 D=1; CODE=0

 DO I=1,N
   AMAX=0.d0
   DO J=1,N
     IF (DABS(A(I,J)).GT.AMAX) AMAX=DABS(A(I,J))
   END DO ! j loop
   IF(AMAX.LT.TINY) THEN
     CODE = 1
     RETURN
   END IF
   VV(I) = 1.d0 / AMAX
 END DO ! i loop

 DO J=1,N
   DO I=1,J-1
     SUM = A(I,J)
     DO K=1,I-1
       SUM = SUM - A(I,K)*A(K,J) 
     END DO ! k loop
     A(I,J) = SUM
   END DO ! i loop
   AMAX = 0.d0
   DO I=J,N
     SUM = A(I,J)
     DO K=1,J-1
       SUM = SUM - A(I,K)*A(K,J) 
     END DO ! k loop
     A(I,J) = SUM
     DUM = VV(I)*DABS(SUM)
     IF(DUM.GE.AMAX) THEN
       IMAX = I
       AMAX = DUM
     END IF
   END DO ! i loop  
   
   IF(J.NE.IMAX) THEN
     DO K=1,N
       DUM = A(IMAX,K)
       A(IMAX,K) = A(J,K)
       A(J,K) = DUM
     END DO ! k loop
     D = -D
     VV(IMAX) = VV(J)
   END IF

   INDX(J) = IMAX
   IF(DABS(A(J,J)) < TINY) A(J,J) = TINY

   IF(J.NE.N) THEN
     DUM = 1.d0 / A(J,J)
     DO I=J+1,N
       A(I,J) = A(I,J)*DUM
     END DO ! i loop
   END IF 
 END DO ! j loop

 RETURN
 END subroutine LUDCMP


!  ******************************************************************
!  * Solves the set of N linear equations A . X = B.  Here A is     *
!  * input, not as the matrix A but rather as its LU decomposition, *
!  * determined by the routine LUDCMP. INDX is input as the permuta-*
!  * tion vector returned by LUDCMP. B is input as the right-hand   *
!  * side vector B, and returns with the solution vector X. A, N and*
!  * INDX are not modified by this routine and can be used for suc- *
!  * cessive calls with different right-hand sides. This routine is *
!  * also efficient for plain matrix inversion.                     *
!  ******************************************************************
 Subroutine LUBKSB(A,N,INDX,B)
 REAL*8  SUM, A(N,N),B(N)
 INTEGER INDX(N)

 II = 0

 DO I=1,N
   LL = INDX(I)
   SUM = B(LL)
   B(LL) = B(I)
   IF(II.NE.0) THEN
     DO J=II,I-1
       SUM = SUM - A(I,J)*B(J)
     END DO ! j loop
   ELSE IF(SUM.NE.0.d0) THEN
     II = I
   END IF
   B(I) = SUM
 END DO ! i loop

 DO I=N,1,-1
   SUM = B(I)
   IF(I < N) THEN
     DO J=I+1,N
       SUM = SUM - A(I,J)*B(J)
     END DO ! j loop
   END IF
   B(I) = SUM / A(I,I)
 END DO ! i loop

 RETURN
 END subroutine LUBKSB

END MODULE LU

! end of file lu.f90
    
