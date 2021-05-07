    !---------------------
    subroutine Curve_Fitting(file$, k, RMS, P, file_out$)
    implicit none
    character(254) file$, file_out$
    integer k
    real(8) RMS, P(0:k),mult !output
    real(8) A(0:k,0:k),sumi(0:2*k),qq,rr,Y(0:k),AA(k+1,k+1),AA_inv(k+1,k+1), A_inv(0:k,0:k),prov(0:k,0:k),Y_calc
    REAL(8), ALLOCATABLE, DIMENSION(:) :: X,B
    integer i,j,i_exp,n_exp,i9,j9,ind_err
    
    !reading from file
    !determination of number of experimental points
    open(1,file = file$)
    do i = 1,100000000
        read(1,*,end = 500)qq,rr
    end do
           
500 n_exp = i-1
    close(1)
    
    allocate(X(N_exp))
    allocate(B(N_exp))
    
     open(1,file = file$)
    do i = 1,n_exp
        read(1,*)X(i),B(i)
        end do
        
        mult = X(n_exp)
    X = x/mult
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
        
    call MATR_INVERT(AA, AA_inv, k+1, ind_err)
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
    
    !write(*,*)'Prov'
    !write(*,'(1x,4(1x,f12.6))')((prov(i9,j9),i9 = 0,k), j9 = 0,k)
   
    
    !write(*,*)'P'
    
    X = x*mult
    
    P = 0.
    do i = 0,k
        do j = 0,k
     p(i) = p(i)+A_inv(i,j)*Y(j)/mult**(i)
        end do  
        
        write(*,*)P(i)
    end do
    open(300,file = file_out$)
    !-----RMS-----
    RMS = 0.
    do j = 1,n_exp
        Y_calc= 0.
        do i = 0,k
            Y_calc = Y_calc+P(i)*X(j)**(i)
                 end do  !i
        write(300,*)X(j),Y_calc
        RMS = RMS+(Y_calc-B(j))**2/n_exp
    end do
    
    RMS = dsqrt(RMS)
    write(300,*)'RMS = ',RMS
    
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
    