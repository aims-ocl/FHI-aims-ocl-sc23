

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Solving Quantum Chemistry Systems using Purification.
MODULE DensityMatrixSolversModule
  USE DataTypesModule, ONLY : NTREAL
  USE EigenBoundsModule, ONLY : GershgorinBounds
  USE LoadBalancerModule, ONLY : PermuteMatrix, UndoPermuteMatrix
  USE LoggingModule, ONLY : WriteElement, WriteListElement, WriteHeader, &
       & WriteCitation, EnterSubLog, ExitSubLog
  USE PMatrixMemoryPoolModule, ONLY : MatrixMemoryPool_p, &
       & DestructMatrixMemoryPool
  USE PSMatrixAlgebraModule, ONLY : IncrementMatrix, MatrixMultiply, &
       & DotMatrix, MatrixTrace, ScaleMatrix
  USE PSMatrixModule, ONLY : Matrix_ps, ConstructEmptyMatrix, DestructMatrix, &
       & CopyMatrix, PrintMatrixInformation, FillMatrixIdentity
  USE SolverParametersModule, ONLY : SolverParameters_t, PrintParameters, &
       & SolverParameters_init, DestructSolverParameters
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Solvers
  PUBLIC :: PM
  PUBLIC :: TRS2
  PUBLIC :: TRS2_cpscf
  PUBLIC :: TRS4
  PUBLIC :: HPCP
  PUBLIC :: ScaleAndFold
  PUBLIC :: EnergyDensityMatrix
  ! PUBLIC :: HPCPPlus
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the density matrix from a Hamiltonian using the PM method.
  !> Based on the PM algorithm presented in \cite palser1998canonical
  SUBROUTINE PM(Hamiltonian, InverseSquareRoot, nel, Density, &
       & energy_value_out, chemical_potential_out, solver_parameters_in)
    !> The matrix to compute the corresponding density from.
    TYPE(Matrix_ps), INTENT(IN) :: Hamiltonian
    !> The inverse square root of the overlap matrix.
    TYPE(Matrix_ps), INTENT(IN) :: InverseSquareRoot
    !> The number of electrons.
    INTEGER, INTENT(IN) :: nel
    !> The density matrix computed by this routine.
    TYPE(Matrix_ps), INTENT(INOUT) :: Density
    !> The energy of the system (optional).
    REAL(NTREAL), INTENT(OUT), OPTIONAL :: energy_value_out
    !> The chemical potential (optional).
    REAL(NTREAL), INTENT(OUT), OPTIONAL :: chemical_potential_out
    !> Parameters for the solver (optional).
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Optional Parameters
    TYPE(SolverParameters_t) :: solver_parameters
    !! Local Matrices
    TYPE(Matrix_ps) :: WorkingHamiltonian
    TYPE(Matrix_ps) :: Identity
    TYPE(Matrix_ps) :: X_k, X_k2, X_k3, TempMat
    !! Local Variables
    REAL(NTREAL) :: e_min, e_max
    REAL(NTREAL) :: factor
    REAL(NTREAL) :: lambda,alpha,alpha1,alpha2
    REAL(NTREAL) :: a1,a2,a3
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: sigma_array
    REAL(NTREAL) :: trace_value
    REAL(NTREAL) :: trace_value2
    REAL(NTREAL) :: norm_value
    REAL(NTREAL) :: energy_value, energy_value2
    !! For computing the chemical potential
    REAL(NTREAL) :: zero_value, midpoint, interval_a, interval_b
    !! Temporary Variables
    TYPE(MatrixMemoryPool_p) :: pool
    INTEGER :: outer_counter, inner_counter
    INTEGER :: total_iterations

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = SolverParameters_init()
    END IF

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Density Matrix Solver")
       CALL EnterSubLog
       CALL WriteElement(key="Method", value="PM")
       CALL WriteCitation("palser1998canonical")
       CALL PrintParameters(solver_parameters)
    END IF

    ALLOCATE(sigma_array(solver_parameters%max_iterations))

    !! Construct All The Necessary Matrices
    CALL ConstructEmptyMatrix(Density, Hamiltonian)
    CALL ConstructEmptyMatrix(WorkingHamiltonian, Hamiltonian)
    CALL ConstructEmptyMatrix(X_k, Hamiltonian)
    CALL ConstructEmptyMatrix(X_k2, Hamiltonian)
    CALL ConstructEmptyMatrix(X_k3, Hamiltonian)
    CALL ConstructEmptyMatrix(TempMat, Hamiltonian)
    CALL ConstructEmptyMatrix(Identity, Hamiltonian)
    CALL FillMatrixIdentity(Identity)

    !! Compute the working hamiltonian.
    CALL MatrixMultiply(InverseSquareRoot,Hamiltonian,TempMat, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
    CALL MatrixMultiply(TempMat,InverseSquareRoot,WorkingHamiltonian, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool)

    !! Load Balancing Step
    IF (solver_parameters%do_load_balancing) THEN
       CALL PermuteMatrix(WorkingHamiltonian, WorkingHamiltonian, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
       CALL PermuteMatrix(Identity, Identity, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
    END IF

    !! Compute the lambda scaling value.
    CALL GershgorinBounds(WorkingHamiltonian,e_min,e_max)

    !! Initialize
    CALL CopyMatrix(WorkingHamiltonian,X_k)

    !! Compute lambda
    CALL MatrixTrace(X_k, trace_value)
    lambda = trace_value/X_k%actual_matrix_dimension

    !! Compute alpha
    alpha1 = nel*0.5_NTREAL/(e_max-lambda)
    alpha2 = (X_k%actual_matrix_dimension-nel*0.5_NTREAL)/(lambda-e_min)
    alpha = MIN(alpha1,alpha2)

    factor = -alpha/X_k%actual_matrix_dimension

    CALL ScaleMatrix(X_k,factor)

    factor = (alpha*lambda+nel*0.5_NTREAL)/X_k%actual_matrix_dimension

    CALL IncrementMatrix(Identity,X_k,alpha_in=factor)

    !! Iterate
    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Iterations")
       CALL EnterSubLog
    END IF
    outer_counter = 1
    norm_value = solver_parameters%converge_diff + 1.0_NTREAL
    energy_value = 0.0_NTREAL
    DO outer_counter = 1,solver_parameters%max_iterations
       !! Compute X_k2
       CALL MatrixMultiply(X_k,X_k,X_k2, &
            & threshold_in=solver_parameters%threshold, &
            & memory_pool_in=pool)

       !! Compute X_k3
       CALL MatrixMultiply(X_k,X_k2,X_k3, &
            & threshold_in=solver_parameters%threshold, &
            & memory_pool_in=pool)

       !! Compute X_k - X_k2
       CALL CopyMatrix(X_k,TempMat)
       CALL IncrementMatrix(X_k2,TempMat, &
            & alpha_in=-1.0_NTREAL, &
            & threshold_in=solver_parameters%threshold)

       !! Compute Sigma
       CALL MatrixTrace(TempMat, trace_value)
       CALL DotMatrix(TempMat,X_k,trace_value2)
       sigma_array(outer_counter) = trace_value2/trace_value

       IF (sigma_array(outer_counter) .GT. 0.5_NTREAL) THEN
          a1 = 0.0_NTREAL
          a2 = 1.0_NTREAL+1.0_NTREAL/sigma_array(outer_counter)
          a3 = -1.0_NTREAL/sigma_array(outer_counter)
       ELSE
          a1 = (1.0_NTREAL-2.0_NTREAL*sigma_array(outer_counter)) &
               & / (1.0_NTREAL-sigma_array(outer_counter))
          a2 = (1.0_NTREAL+sigma_array(outer_counter)) &
               & / (1.0_NTREAL-sigma_array(outer_counter))
          a3 = -1.0_NTREAL/(1.0_NTREAL-sigma_array(outer_counter))
       END IF

       !! Update X_k
       CALL ScaleMatrix(X_k,a1)
       CALL IncrementMatrix(X_k2,X_k, &
            & alpha_in=a2, &
            & threshold_in=solver_parameters%threshold)
       CALL IncrementMatrix(X_k3,X_k, &
            & alpha_in=a3, &
            & threshold_in=solver_parameters%threshold)

       !! Energy value based convergence
       energy_value2 = energy_value
       CALL DotMatrix(X_k, WorkingHamiltonian, energy_value)
       energy_value = 2.0_NTREAL*energy_value
       norm_value = ABS(energy_value - energy_value2)

       IF (solver_parameters%be_verbose) THEN
          CALL WriteListElement(key="Round", value=outer_counter)
          CALL EnterSubLog
          CALL WriteElement(key="Convergence", value=norm_value)
          CALL WriteElement("Energy_Value", value=energy_value)
          CALL ExitSubLog
       END IF

       IF (norm_value .LE. solver_parameters%converge_diff) THEN
          EXIT
       END IF
    END DO
    total_iterations = outer_counter-1
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
       CALL WriteElement(key="Total_Iterations", value=outer_counter)
       CALL PrintMatrixInformation(X_k)
    END IF

    IF (PRESENT(energy_value_out)) THEN
       energy_value_out = energy_value
    END IF

    !! Undo Load Balancing Step
    IF (solver_parameters%do_load_balancing) THEN
       CALL UndoPermuteMatrix(X_k, X_k, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
    END IF

    !! Compute the density matrix in the non-orthogonalized basis
    CALL MatrixMultiply(InverseSquareRoot,X_k,TempMat, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
    CALL MatrixMultiply(TempMat,InverseSquareRoot,Density, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool)

    !! Cleanup
    CALL DestructMatrix(WorkingHamiltonian)
    CALL DestructMatrix(X_k)
    CALL DestructMatrix(X_k2)
    CALL DestructMatrix(X_k3)
    CALL DestructMatrix(TempMat)
    CALL DestructMatrix(Identity)
    CALL DestructMatrixMemoryPool(pool)

    !! Compute The Chemical Potential
    IF (PRESENT(chemical_potential_out)) THEN
       interval_a = 0.0_NTREAL
       interval_b = 1.0_NTREAL
       midpoint = 0.0_NTREAL
       midpoints: DO outer_counter = 1,solver_parameters%max_iterations
          midpoint = (interval_b - interval_a)/2.0_NTREAL + interval_a
          zero_value = midpoint
          !! Compute polynomial function at the guess point.
          polynomial: DO inner_counter=1,total_iterations
             IF (sigma_array(inner_counter) .GT. 0.5_NTREAL) THEN
                zero_value = ((1.0_NTREAL+sigma_array(inner_counter)) &
                     & *zero_value**2) - (zero_value**3)
                zero_value = zero_value/sigma_array(inner_counter)
             ELSE
                zero_value = ((1.0_NTREAL - 2.0_NTREAL* &
                     & sigma_array(inner_counter))*zero_value) &
                     & + ((1.0_NTREAL+sigma_array(inner_counter))* &
                     & zero_value**2) - (zero_value**3)
                zero_value = zero_value/(1.0_NTREAL-sigma_array(inner_counter))
             END IF
          END DO polynomial
          !! Change bracketing.
          IF (zero_value .LT. 0.5_NTREAL) THEN
             interval_a = midpoint
          ELSE
             interval_b = midpoint
          END IF
          !! Check convergence.
          IF (ABS(zero_value-0.5_NTREAL) .LT. &
               & solver_parameters%converge_diff) THEN
             EXIT
          END IF
       END DO midpoints
       !! Undo scaling.
       chemical_potential_out = lambda - &
            & (Hamiltonian%actual_matrix_dimension*midpoint - nel*0.5_NTREAL) &
            & /alpha
    END IF

    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
    END IF

    !! Cleanup
    DEALLOCATE(sigma_array)
    CALL DestructSolverParameters(solver_parameters)
  END SUBROUTINE PM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the density matrix from a Hamiltonian using the TRS2 method.
  !> Based on the TRS2 algorithm presented in \cite niklasson2002.
  SUBROUTINE TRS2(Hamiltonian, InverseSquareRoot, nel, Density, &
       & energy_value_out, chemical_potential_out, solver_parameters_in)
    !> The matrix to compute the corresponding density from
    TYPE(Matrix_ps), INTENT(IN) :: Hamiltonian
    !> The inverse square root of the overlap matrix.
    TYPE(Matrix_ps), INTENT(IN) :: InverseSquareRoot
    !> The number of electrons.
    INTEGER, INTENT(IN) :: nel
    !> The density matrix computed by this routine.
    TYPE(Matrix_ps), INTENT(INOUT) :: Density
    !> The energy of the system (optional).
    REAL(NTREAL), INTENT(OUT), OPTIONAL :: energy_value_out
    !> The chemical potential (optional).
    REAL(NTREAL), INTENT(OUT), OPTIONAL :: chemical_potential_out
    !> Parameters for the solver (optional).
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Optional Parameters
    TYPE(SolverParameters_t) :: solver_parameters
    !! Local Matrices
    TYPE(Matrix_ps) :: WorkingHamiltonian
    TYPE(Matrix_ps) :: Identity
    TYPE(Matrix_ps) :: X_k, X_k2, TempMat
    !! Local Variables
    REAL(NTREAL) :: e_min, e_max
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: sigma_array
    REAL(NTREAL) :: trace_value
    REAL(NTREAL) :: norm_value
    REAL(NTREAL) :: energy_value, energy_value2
    !! For computing the chemical potential
    REAL(NTREAL) :: zero_value, midpoint, interval_a, interval_b
    !! Temporary Variables
    TYPE(MatrixMemoryPool_p) :: pool
    INTEGER :: outer_counter, inner_counter
    INTEGER :: total_iterations

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = SolverParameters_init()
    END IF

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Density Matrix Solver")
       CALL EnterSubLog
       CALL WriteElement(key="Method", value="TRS2")
       CALL WriteCitation("niklasson2002expansion")
       CALL PrintParameters(solver_parameters)
    END IF

    ALLOCATE(sigma_array(solver_parameters%max_iterations))

    !! Construct All The Necessary Matrices
    CALL ConstructEmptyMatrix(Density, Hamiltonian)
    CALL ConstructEmptyMatrix(WorkingHamiltonian, Hamiltonian)
    CALL ConstructEmptyMatrix(X_k, Hamiltonian)
    CALL ConstructEmptyMatrix(X_k2, Hamiltonian)
    CALL ConstructEmptyMatrix(TempMat, Hamiltonian)
    CALL ConstructEmptyMatrix(Identity, Hamiltonian)
    CALL FillMatrixIdentity(Identity)

    !! Compute the working hamiltonian.
    CALL MatrixMultiply(InverseSquareRoot,Hamiltonian,TempMat, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
    CALL MatrixMultiply(TempMat,InverseSquareRoot,WorkingHamiltonian, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool)

    !! Load Balancing Step
    IF (solver_parameters%do_load_balancing) THEN
       CALL PermuteMatrix(WorkingHamiltonian, WorkingHamiltonian, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
       CALL PermuteMatrix(Identity, Identity, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
    END IF

    !! Compute the lambda scaling value.
    CALL GershgorinBounds(WorkingHamiltonian,e_min,e_max)

    !! Initialize
    CALL CopyMatrix(WorkingHamiltonian,X_k)
    CALL ScaleMatrix(X_k,-1.0_NTREAL)
    CALL IncrementMatrix(Identity,X_k,alpha_in=e_max)
    CALL ScaleMatrix(X_k,1.0_NTREAL/(e_max-e_min))

    !! Iterate
    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Iterations")
       CALL EnterSubLog
    END IF
    outer_counter = 1
    norm_value = solver_parameters%converge_diff + 1.0_NTREAL
    energy_value = 0.0_NTREAL

       !write(6,*) 'shanghui test in NTpoly(num, electron_num, trace):'

    DO outer_counter = 1,solver_parameters%max_iterations
       !! Compute Sigma
       CALL MatrixTrace(X_k, trace_value)

       IF (nel*0.5_NTREAL .LT. trace_value) THEN
          sigma_array(outer_counter) = -1.0_NTREAL
       ELSE
          sigma_array(outer_counter) = 1.0_NTREAL
       END IF

       !! Compute X_k2
       CALL MatrixMultiply(X_k,X_k,X_k2, &
            & threshold_in=solver_parameters%threshold, &
            & memory_pool_in=pool)

       !! Update X_k
       IF (sigma_array(outer_counter) .GT. 0.0_NTREAL) THEN
          CALL ScaleMatrix(X_k,2.0_NTREAL)
          CALL IncrementMatrix(X_k2,X_k, &
               & alpha_in=-1.0_NTREAL, &
               & threshold_in=solver_parameters%threshold)
       ELSE
          CALL CopyMatrix(X_k2,X_k)
       END IF

       !! Energy value based convergence
       energy_value2 = energy_value
       CALL DotMatrix(X_k, WorkingHamiltonian, energy_value)
       energy_value = 2.0_NTREAL*energy_value
       norm_value = ABS(energy_value - energy_value2)


       IF (solver_parameters%be_verbose) THEN
          CALL WriteListElement(key="Round", value=outer_counter)
          CALL EnterSubLog
          CALL WriteElement(key="Convergence", value=norm_value)
          CALL WriteElement("Energy_Value", value=energy_value)
          CALL ExitSubLog
       END IF

       IF (norm_value .LE. solver_parameters%converge_diff) THEN
          EXIT
       END IF
    END DO
    total_iterations = outer_counter-1
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
       CALL WriteElement(key="Total_Iterations", value=outer_counter)
       CALL PrintMatrixInformation(X_k)
    END IF

    IF (PRESENT(energy_value_out)) THEN
       energy_value_out = energy_value
    END IF

    !! Undo Load Balancing Step
    IF (solver_parameters%do_load_balancing) THEN
       CALL UndoPermuteMatrix(X_k, X_k, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
    END IF

    !! Compute the density matrix in the non-orthogonalized basis
    CALL MatrixMultiply(InverseSquareRoot,X_k,TempMat, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
    CALL MatrixMultiply(TempMat,InverseSquareRoot,Density, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool)

    !! Cleanup
    CALL DestructMatrix(WorkingHamiltonian)
    CALL DestructMatrix(X_k)
    CALL DestructMatrix(X_k2)
    CALL DestructMatrix(TempMat)
    CALL DestructMatrix(Identity)
    CALL DestructMatrixMemoryPool(pool)

    !! Compute The Chemical Potential
    IF (PRESENT(chemical_potential_out)) THEN
       interval_a = 0.0_NTREAL
       interval_b = 1.0_NTREAL
       midpoint = 0.0_NTREAL
       midpoints: DO outer_counter = 1,solver_parameters%max_iterations
          midpoint = (interval_b - interval_a)/2.0_NTREAL + interval_a
          zero_value = midpoint
          !! Compute polynomial function at the guess point.
          polynomial: DO inner_counter=1,total_iterations
             IF (sigma_array(inner_counter) .LT. 0.0_NTREAL) THEN
                zero_value = zero_value*zero_value
             ELSE
                zero_value = 2.0_NTREAL*zero_value - zero_value*zero_value
             END IF
          END DO polynomial
          !! Change bracketing.
          IF (zero_value .LT. 0.5_NTREAL) THEN
             interval_a = midpoint
          ELSE
             interval_b = midpoint
          END IF
          !! Check convergence.
          IF (ABS(zero_value-0.5_NTREAL) .LT. &
               & solver_parameters%converge_diff) THEN
             EXIT
          END IF
       END DO midpoints
       !! Undo scaling.
       chemical_potential_out = e_max + (e_min - e_max)*midpoint
    END IF

    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
    END IF

    !! Cleanup
    DEALLOCATE(sigma_array)
    CALL DestructSolverParameters(solver_parameters)
  END SUBROUTINE TRS2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the first_order density matrix from a first_order Hamiltonian using the TRS2_cpscf method.
  !> Based on the TRS2_cpscf algorithm.
  SUBROUTINE TRS2_cpscf(Hamiltonian, InverseSquareRoot, nel, Density, &
                        first_order_H, first_order_S, first_order_DM, &               
       & energy_value_out, chemical_potential_out, solver_parameters_in)
    !> The matrix to compute the corresponding density from
    TYPE(Matrix_ps), INTENT(IN) :: Hamiltonian
    !> The inverse square root of the overlap matrix.
    TYPE(Matrix_ps), INTENT(IN) :: InverseSquareRoot
    !> The number of electrons.
    INTEGER, INTENT(IN) :: nel
    !> The density matrix computed by this routine.
    TYPE(Matrix_ps), INTENT(INOUT) :: Density
    !> 
    TYPE(Matrix_ps), INTENT(IN) :: first_order_H
    !> this is zeros
    TYPE(Matrix_ps), INTENT(IN) :: first_order_S
    !> The first_order_density_matrix computed by this routine.
    TYPE(Matrix_ps), INTENT(INOUT) :: first_order_DM

    !> The energy of the system (optional).
    REAL(NTREAL), INTENT(OUT), OPTIONAL :: energy_value_out
    !> The chemical potential (optional).
    REAL(NTREAL), INTENT(OUT), OPTIONAL :: chemical_potential_out
    !> Parameters for the solver (optional).
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Optional Parameters
    TYPE(SolverParameters_t) :: solver_parameters
    !! Local Matrices
    TYPE(Matrix_ps) :: WorkingHamiltonian
    TYPE(Matrix_ps) :: Identity
    TYPE(Matrix_ps) :: X_k, X_k2, TempMat
    TYPE(Matrix_ps) :: first_order_WorkingHamiltonian
    TYPE(Matrix_ps) :: first_order_X_k, first_order_X_k2
    !! Local Variables
    REAL(NTREAL) :: e_min, e_max
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: sigma_array
    REAL(NTREAL) :: trace_value
    REAL(NTREAL) :: norm_value
    REAL(NTREAL) :: energy_value, energy_value2
    !! For computing the chemical potential
    REAL(NTREAL) :: zero_value, midpoint, interval_a, interval_b
    !! Temporary Variables
    TYPE(MatrixMemoryPool_p) :: pool
    INTEGER :: outer_counter, inner_counter
    INTEGER :: total_iterations

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = SolverParameters_init()
    END IF

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("first_order_Density Matrix Solver")
       CALL EnterSubLog
       CALL WriteElement(key="Method", value="TRS2_cpscf")
       CALL WriteCitation("niklasson2002expansion")
       CALL PrintParameters(solver_parameters)
    END IF

    ALLOCATE(sigma_array(solver_parameters%max_iterations))

    !! Construct All The Necessary Matrices
    CALL ConstructEmptyMatrix(Density, Hamiltonian)
    CALL ConstructEmptyMatrix(WorkingHamiltonian, Hamiltonian)
    CALL ConstructEmptyMatrix(X_k, Hamiltonian)
    CALL ConstructEmptyMatrix(X_k2, Hamiltonian)
    CALL ConstructEmptyMatrix(TempMat, Hamiltonian)
    CALL ConstructEmptyMatrix(Identity, Hamiltonian)
    CALL FillMatrixIdentity(Identity)

    CALL ConstructEmptyMatrix(first_order_DM, Hamiltonian)
    CALL ConstructEmptyMatrix(first_order_WorkingHamiltonian, Hamiltonian)
    CALL ConstructEmptyMatrix(first_order_X_k, Hamiltonian)
    CALL ConstructEmptyMatrix(first_order_X_k2, Hamiltonian)

    !! Compute the working hamiltonian.
      ! H(0) = S^{-1/2} H_ao(0) S^{-1/2}
    CALL MatrixMultiply(InverseSquareRoot,Hamiltonian,TempMat, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
    CALL MatrixMultiply(TempMat,InverseSquareRoot,WorkingHamiltonian, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool)

      ! H(1) = S^{-1/2} H_ao(1) S^{-1/2}
    CALL MatrixMultiply(InverseSquareRoot,first_order_H,TempMat, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
    CALL MatrixMultiply(TempMat,InverseSquareRoot, first_order_WorkingHamiltonian, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool)

    !! Load Balancing Step
    IF (solver_parameters%do_load_balancing) THEN
       CALL PermuteMatrix(WorkingHamiltonian, WorkingHamiltonian, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
       CALL PermuteMatrix(Identity, Identity, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
       CALL PermuteMatrix(first_order_WorkingHamiltonian, first_order_WorkingHamiltonian, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
    END IF

    !! Compute the lambda scaling value.
    CALL GershgorinBounds(WorkingHamiltonian,e_min,e_max)

    !! Initialize 
       ! D(0)= ( e_max-H(0) )/( e_max-e_min )
    CALL CopyMatrix(WorkingHamiltonian,X_k)
    CALL ScaleMatrix(X_k,-1.0_NTREAL)
    CALL IncrementMatrix(Identity,X_k,alpha_in=e_max)
    CALL ScaleMatrix(X_k,1.0_NTREAL/(e_max-e_min))
       ! D(1)= ( -H(1) )/( e_max-e_min )
    CALL CopyMatrix(first_order_WorkingHamiltonian, first_order_X_k)
    CALL ScaleMatrix(first_order_X_k, -1.0_NTREAL/(e_max-e_min))


    !! Iterate
    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Iterations")
       CALL EnterSubLog
    END IF
    outer_counter = 1
    norm_value = solver_parameters%converge_diff + 1.0_NTREAL
    energy_value = 0.0_NTREAL


    DO outer_counter = 1,solver_parameters%max_iterations
       !! Compute Sigma
       CALL MatrixTrace(X_k, trace_value)

       IF (nel*0.5_NTREAL .LT. trace_value) THEN
          sigma_array(outer_counter) = -1.0_NTREAL
       ELSE
          sigma_array(outer_counter) = 1.0_NTREAL
       END IF

       !! Compute X_k2  = D(0)*D(0)
       CALL MatrixMultiply(X_k,X_k,X_k2, &
            & threshold_in=solver_parameters%threshold, &
            & memory_pool_in=pool)
       !! Compute first_order_X_k2 = D(1)*D(0) + D(0)*D(1)
       CALL MatrixMultiply(first_order_X_k,X_k,first_order_X_k2, &
            & threshold_in=solver_parameters%threshold, &
            & memory_pool_in=pool)
       CALL MatrixMultiply(X_k,first_order_X_k, TempMat, &
            & threshold_in=solver_parameters%threshold, &
            & memory_pool_in=pool)
       CALL IncrementMatrix(TempMat, first_order_X_k2, &
            & alpha_in=1.0_NTREAL, &
            & threshold_in=solver_parameters%threshold)


       !! Update X_k
       IF (sigma_array(outer_counter) .GT. 0.0_NTREAL) THEN

          ! D(0)_{n+1} = 2*D(0)_{n} - D(0)_{n}*D(0)_{n} 
          CALL ScaleMatrix(X_k,2.0_NTREAL)
          CALL IncrementMatrix(X_k2,X_k, &
               & alpha_in=-1.0_NTREAL, &
               & threshold_in=solver_parameters%threshold)

          ! D(1)_{n+1} = 2*D(1)_{n} - D(1)_{n}*D(0)_{n} - D(0)_{n}*D(1)_{n}
          CALL ScaleMatrix(first_order_X_k,2.0_NTREAL)
          CALL IncrementMatrix(first_order_X_k2, first_order_X_k, &
               & alpha_in=-1.0_NTREAL, &
               & threshold_in=solver_parameters%threshold)

       ELSE
          ! D(0)_{n+1} = D(0)_{n}*D(0)_{n} 
          CALL CopyMatrix(X_k2,X_k)
    
          ! D(1)_{n+1} = D(1)_{n}*D(0)_{n} - D(0)_{n}*D(1)_{n}
          CALL CopyMatrix(first_order_X_k2, first_order_X_k)

       END IF

       !! Energy value based convergence
       energy_value2 = energy_value
       CALL DotMatrix(X_k, WorkingHamiltonian, energy_value)
       energy_value = 2.0_NTREAL*energy_value
       norm_value = ABS(energy_value - energy_value2)
       !write(6,*) outer_counter,energy_value 


       IF (solver_parameters%be_verbose) THEN
          CALL WriteListElement(key="Round", value=outer_counter)
          CALL EnterSubLog
          CALL WriteElement(key="Convergence", value=norm_value)
          CALL WriteElement("Energy_Value", value=energy_value)
          CALL ExitSubLog
       END IF

       IF (norm_value .LE. solver_parameters%converge_diff) THEN
          EXIT
       END IF
    END DO
    total_iterations = outer_counter-1
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
       CALL WriteElement(key="Total_Iterations", value=outer_counter)
       CALL PrintMatrixInformation(first_order_X_k)
    END IF

    IF (PRESENT(energy_value_out)) THEN
       energy_value_out = energy_value
    END IF

    !! Undo Load Balancing Step
    IF (solver_parameters%do_load_balancing) THEN
       CALL UndoPermuteMatrix(X_k, X_k, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
       CALL UndoPermuteMatrix(first_order_X_k, first_order_X_k, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
    END IF

    !! Compute the density matrix in the non-orthogonalized basis
      ! DM_ao(0) = S^{-1/2} DM(0) S^{-1/2}
    CALL MatrixMultiply(InverseSquareRoot,X_k,TempMat, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
    CALL MatrixMultiply(TempMat,InverseSquareRoot,Density, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool)

      ! DM_ao(1) = S^{-1/2} DM(1) S^{-1/2}
    CALL MatrixMultiply(InverseSquareRoot,first_order_X_k,TempMat, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
    CALL MatrixMultiply(TempMat,InverseSquareRoot,first_order_DM, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool)

    !! Cleanup
    CALL DestructMatrix(WorkingHamiltonian)
    CALL DestructMatrix(X_k)
    CALL DestructMatrix(X_k2)
    CALL DestructMatrix(TempMat)
    CALL DestructMatrix(Identity)
    CALL DestructMatrixMemoryPool(pool)

    CALL DestructMatrix(first_order_WorkingHamiltonian)
    CALL DestructMatrix(first_order_X_k)
    CALL DestructMatrix(first_order_X_k2)


    !! Compute The Chemical Potential
    IF (PRESENT(chemical_potential_out)) THEN
       interval_a = 0.0_NTREAL
       interval_b = 1.0_NTREAL
       midpoint = 0.0_NTREAL
       midpoints: DO outer_counter = 1,solver_parameters%max_iterations
          midpoint = (interval_b - interval_a)/2.0_NTREAL + interval_a
          zero_value = midpoint
          !! Compute polynomial function at the guess point.
          polynomial: DO inner_counter=1,total_iterations
             IF (sigma_array(inner_counter) .LT. 0.0_NTREAL) THEN
                zero_value = zero_value*zero_value
             ELSE
                zero_value = 2.0_NTREAL*zero_value - zero_value*zero_value
             END IF
          END DO polynomial
          !! Change bracketing.
          IF (zero_value .LT. 0.5_NTREAL) THEN
             interval_a = midpoint
          ELSE
             interval_b = midpoint
          END IF
          !! Check convergence.
          IF (ABS(zero_value-0.5_NTREAL) .LT. &
               & solver_parameters%converge_diff) THEN
             EXIT
          END IF
       END DO midpoints
       !! Undo scaling.
       chemical_potential_out = e_max + (e_min - e_max)*midpoint
    END IF

    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
    END IF

    !! Cleanup
    DEALLOCATE(sigma_array)
    CALL DestructSolverParameters(solver_parameters)
  END SUBROUTINE TRS2_cpscf






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the density matrix from a Hamiltonian using the TRS4 method.
  !> Based on the TRS4 algorithm presented in \cite niklasson2002
  SUBROUTINE TRS4(Hamiltonian, InverseSquareRoot, nel, Density, &
       & energy_value_out, chemical_potential_out, solver_parameters_in)
    !> The matrix to compute the corresponding density from.
    TYPE(Matrix_ps), INTENT(IN)  :: Hamiltonian
    !> The inverse square root of the overlap matrix.
    TYPE(Matrix_ps), INTENT(IN) :: InverseSquareRoot
    !> The number of electrons.
    INTEGER, INTENT(IN) :: nel
    !> The density matrix computed by this routine.
    TYPE(Matrix_ps), INTENT(INOUT) :: Density
    !> The energy of the system (optional).
    REAL(NTREAL), INTENT(OUT), OPTIONAL :: energy_value_out
    !> The chemical potential (optional).
    REAL(NTREAL), INTENT(OUT), OPTIONAL :: chemical_potential_out
    !> Parameters for the solver (optional).
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    REAL(NTREAL), PARAMETER :: sigma_min = 0.0_NTREAL
    REAL(NTREAL), PARAMETER :: sigma_max = 6.0_NTREAL
    !! Handling Optional Parameters
    TYPE(SolverParameters_t) :: solver_parameters
    !! Local Matrices
    TYPE(Matrix_ps) :: WorkingHamiltonian
    TYPE(Matrix_ps) :: Identity
    TYPE(Matrix_ps) :: X_k, X_k2, Fx_right, GX_right, TempMat
    !! Local Variables
    REAL(NTREAL) :: e_min, e_max
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: sigma_array
    REAL(NTREAL) :: norm_value
    REAL(NTREAL) :: energy_value, energy_value2
    !! For computing the chemical potential
    REAL(NTREAL) :: zero_value, midpoint, interval_a, interval_b
    REAL(NTREAL) :: tempfx,tempgx
    !! Temporary Variables
    TYPE(MatrixMemoryPool_p) :: pool
    INTEGER :: outer_counter, inner_counter
    INTEGER :: total_iterations
    REAL(NTREAL) :: trace_fx, trace_gx

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = SolverParameters_init()
    END IF

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Density Matrix Solver")
       CALL EnterSubLog
       CALL WriteElement(key="Method", value="TRS4")
       CALL WriteCitation("niklasson2002expansion")
       CALL PrintParameters(solver_parameters)
    END IF

    ALLOCATE(sigma_array(solver_parameters%max_iterations))

    !! Construct All The Necessary Matrices
    CALL ConstructEmptyMatrix(Density, Hamiltonian)
    CALL ConstructEmptyMatrix(WorkingHamiltonian, Hamiltonian)
    CALL ConstructEmptyMatrix(X_k, Hamiltonian)
    CALL ConstructEmptyMatrix(X_k2, Hamiltonian)
    CALL ConstructEmptyMatrix(TempMat, Hamiltonian)
    CALL ConstructEmptyMatrix(Fx_right, Hamiltonian)
    CALL ConstructEmptyMatrix(Gx_right, Hamiltonian)
    CALL ConstructEmptyMatrix(Identity, Hamiltonian)
    CALL FillMatrixIdentity(Identity)

    !! Compute the working hamiltonian.
    CALL MatrixMultiply(InverseSquareRoot,Hamiltonian,TempMat, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
    CALL MatrixMultiply(TempMat,InverseSquareRoot,WorkingHamiltonian, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool)

    !! Load Balancing Step
    IF (solver_parameters%do_load_balancing) THEN
       CALL PermuteMatrix(WorkingHamiltonian, WorkingHamiltonian, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
       CALL PermuteMatrix(Identity, Identity, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
    END IF

    !! Compute the lambda scaling value.
    CALL GershgorinBounds(WorkingHamiltonian,e_min,e_max)

    !! Initialize
    CALL CopyMatrix(WorkingHamiltonian,X_k)
    CALL ScaleMatrix(X_k,-1.0_NTREAL)
    CALL IncrementMatrix(Identity,X_k,alpha_in=e_max)
    CALL ScaleMatrix(X_k,1.0_NTREAL/(e_max-e_min))

    !! Iterate
    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Iterations")
       CALL EnterSubLog
    END IF
    outer_counter = 1
    norm_value = solver_parameters%converge_diff + 1.0_NTREAL
    energy_value = 0.0_NTREAL
    DO outer_counter = 1,solver_parameters%max_iterations
       !! Compute X_k2
       CALL MatrixMultiply(X_k, X_k, X_k2, &
            & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
       !! Compute Fx_right
       CALL CopyMatrix(X_k2,Fx_right)
       CALL ScaleMatrix(Fx_right,-3.0_NTREAL)
       CALL IncrementMatrix(X_k,Fx_right,alpha_in=4.0_NTREAL)
       !! Compute Gx_right
       CALL CopyMatrix(Identity,Gx_right)
       CALL IncrementMatrix(X_k,Gx_right,alpha_in=-2.0_NTREAL)
       CALL IncrementMatrix(X_k2,Gx_right)

       !! Compute Traces
       CALL DotMatrix(X_k2,Fx_right,trace_fx)
       CALL DotMatrix(X_k2,Gx_right,trace_gx)

       !! Avoid Overflow
       IF (ABS(trace_gx) .LT. 1.0e-14_NTREAL) THEN
          EXIT
       END IF

       !! Compute Sigma
       sigma_array(outer_counter) = (nel*0.5_NTREAL-trace_fx)/trace_gx

       !! Update The Matrix
       IF (sigma_array(outer_counter) .GT. sigma_max) THEN
          CALL CopyMatrix(X_k, TempMat)
          CALL ScaleMatrix(TempMat, 2.0_NTREAL)
          CALL IncrementMatrix(X_k2, TempMat, alpha_in=-1.0_NTREAL)
       ELSE IF (sigma_array(outer_counter) .LT. sigma_min) THEN
          CALL CopyMatrix(X_k2, TempMat)
       ELSE
          CALL ScaleMatrix(Gx_right,sigma_array(outer_counter))
          CALL IncrementMatrix(Fx_right,Gx_right)
          CALL MatrixMultiply(X_k2, Gx_right, TempMat, &
               & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
       END IF

       CALL IncrementMatrix(TempMat,X_k,alpha_in=-1.0_NTREAL)
       CALL CopyMatrix(TempMat,X_k)

       !! Energy value based convergence
       energy_value2 = energy_value
       CALL DotMatrix(X_k,WorkingHamiltonian,energy_value)
       energy_value = 2.0_NTREAL*energy_value
       norm_value = ABS(energy_value - energy_value2)

       IF (solver_parameters%be_verbose) THEN
          CALL WriteListElement(key="Round", value=outer_counter)
          CALL EnterSubLog
          CALL WriteElement(key="Convergence", value=norm_value)
          CALL WriteElement("Energy_Value", value=energy_value)
          CALL ExitSubLog
       END IF

       IF (norm_value .LE. solver_parameters%converge_diff) THEN
          EXIT
       END IF
    END DO
    total_iterations = outer_counter-1
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
       CALL WriteElement(key="Total_Iterations", value=outer_counter)
       CALL PrintMatrixInformation(X_k)
    END IF

    IF (PRESENT(energy_value_out)) THEN
       energy_value_out = energy_value
    END IF

    !! Undo Load Balancing Step
    IF (solver_parameters%do_load_balancing) THEN
       CALL UndoPermuteMatrix(X_k, X_k, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
    END IF

    !! Compute the density matrix in the non-orthogonalized basis
    CALL MatrixMultiply(InverseSquareRoot,X_k,TempMat, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
    CALL MatrixMultiply(TempMat,InverseSquareRoot,Density, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool)

    !! Cleanup
    CALL DestructMatrix(WorkingHamiltonian)
    CALL DestructMatrix(X_k)
    CALL DestructMatrix(X_k2)
    CALL DestructMatrix(Fx_right)
    CALL DestructMatrix(Gx_right)
    CALL DestructMatrix(TempMat)
    CALL DestructMatrix(Identity)
    CALL DestructMatrixMemoryPool(pool)

    !! Compute The Chemical Potential
    IF (PRESENT(chemical_potential_out)) THEN
       interval_a = 0.0_NTREAL
       interval_b = 1.0_NTREAL
       midpoint = 0.0_NTREAL
       midpoints: DO outer_counter = 1,solver_parameters%max_iterations
          midpoint = (interval_b - interval_a)/2.0_NTREAL + interval_a
          zero_value = midpoint
          !! Compute polynomial function at the guess point.
          polynomial: DO inner_counter=1,total_iterations
             IF (sigma_array(inner_counter) .GT. sigma_max) THEN
                zero_value = 2.0_NTREAL*zero_value - zero_value*zero_value
             ELSE IF (sigma_array(inner_counter) .LT. sigma_min) THEN
                zero_value = zero_value*zero_value
             ELSE
                tempfx = (zero_value*zero_value) * &
                     & (4.0_NTREAL*zero_value - &
                     & 3.0_NTREAL*zero_value*zero_value)
                tempgx = (zero_value*zero_value) * (1.0_NTREAL-zero_value) &
                     & * (1.0_NTREAL-zero_value)
                zero_value = tempfx + sigma_array(inner_counter)*tempgx
             END IF
          END DO polynomial
          !! Change bracketing.
          IF (zero_value .LT. 0.5_NTREAL) THEN
             interval_a = midpoint
          ELSE
             interval_b = midpoint
          END IF
          !! Check convergence.
          IF (ABS(zero_value-0.5_NTREAL) .LT. &
               & solver_parameters%converge_diff) THEN
             EXIT
          END IF
       END DO midpoints
       !! Undo scaling.
       chemical_potential_out = e_max + (e_min - e_max)*midpoint
    END IF

    !! Cleanup
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
    END IF

    DEALLOCATE(sigma_array)
    CALL DestructSolverParameters(solver_parameters)
  END SUBROUTINE TRS4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the density matrix from a Hamiltonian using the HPCP method.
  !> Based on the algorithm presented in \cite truflandier2016communication.
  SUBROUTINE HPCP(Hamiltonian, InverseSquareRoot, nel, Density, &
       & energy_value_out, chemical_potential_out, solver_parameters_in)
    !> The matrix to compute the corresponding density from.
    TYPE(Matrix_ps), INTENT(IN) :: Hamiltonian
    !> The inverse square root of the overlap matrix.
    TYPE(Matrix_ps), INTENT(IN) :: InverseSquareRoot
    !> The number of electrons.
    INTEGER, INTENT(IN) :: nel
    !> The density matrix computed by this routine.
    TYPE(Matrix_ps), INTENT(INOUT) :: Density
    !> The energy of the system (optional).
    REAL(NTREAL), INTENT(OUT), OPTIONAL :: energy_value_out
    !> The chemical potential (optional).
    REAL(NTREAL), INTENT(OUT), OPTIONAL :: chemical_potential_out
    !> Parameters for the solver (optional).
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Optional Parameters
    TYPE(SolverParameters_t) :: solver_parameters
    !! Local Matrices
    TYPE(Matrix_ps) :: WorkingHamiltonian
    TYPE(Matrix_ps) :: TempMat
    TYPE(Matrix_ps) :: Identity
    TYPE(Matrix_ps) :: D1, DH, DDH, D2DH
    !! Local Variables
    REAL(NTREAL) :: e_min, e_max
    REAL(NTREAL) :: beta_1, beta_2
    REAL(NTREAL) :: beta, beta_bar
    REAL(NTREAL) :: sigma, sigma_bar
    REAL(NTREAL) :: mu
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: sigma_array
    REAL(NTREAL) :: trace_value
    REAL(NTREAL) :: norm_value, norm_value2
    REAL(NTREAL) :: energy_value, energy_value2
    !! For computing the chemical potential
    REAL(NTREAL) :: zero_value, midpoint, interval_a, interval_b
    !! Temporary Variables
    TYPE(MatrixMemoryPool_p) :: pool
    INTEGER :: outer_counter, inner_counter
    INTEGER :: total_iterations
    INTEGER :: matrix_dimension

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = SolverParameters_init()
    END IF

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Density Matrix Solver")
       CALL EnterSubLog
       CALL WriteElement(key="Method", value="HPCP")
       CALL WriteCitation("truflandier2016communication")
       CALL PrintParameters(solver_parameters)
    END IF

    ALLOCATE(sigma_array(solver_parameters%max_iterations))

    matrix_dimension = Hamiltonian%actual_matrix_dimension

    !! Construct All The Necessary Matrices
    CALL ConstructEmptyMatrix(Density, Hamiltonian)
    CALL ConstructEmptyMatrix(WorkingHamiltonian, Hamiltonian)
    CALL ConstructEmptyMatrix(TempMat, Hamiltonian)
    CALL ConstructEmptyMatrix(D1, Hamiltonian)
    CALL ConstructEmptyMatrix(DH, Hamiltonian)
    CALL ConstructEmptyMatrix(DDH, Hamiltonian)
    CALL ConstructEmptyMatrix(D2DH, Hamiltonian)
    CALL ConstructEmptyMatrix(Identity, Hamiltonian)
    CALL FillMatrixIdentity(Identity)

    !! Compute the working hamiltonian.
    CALL MatrixMultiply(InverseSquareRoot,Hamiltonian,TempMat, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
    CALL MatrixMultiply(TempMat,InverseSquareRoot,WorkingHamiltonian, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool)

    !! Load Balancing Step
    IF (solver_parameters%do_load_balancing) THEN
       CALL PermuteMatrix(WorkingHamiltonian, WorkingHamiltonian, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
       CALL PermuteMatrix(Identity, Identity, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
    END IF

    !! Compute the initial matrix.
    CALL GershgorinBounds(WorkingHamiltonian,e_min,e_max)
    CALL MatrixTrace(WorkingHamiltonian, mu)
    mu = mu/matrix_dimension
    sigma_bar = (matrix_dimension - 0.5_NTREAL*nel)/matrix_dimension
    sigma = 1.0_NTREAL - sigma_bar
    beta = sigma/(e_max - mu)
    beta_bar = sigma_bar/(mu - e_min)
    beta_1 = sigma
    beta_2 = MIN(beta,beta_bar)

    !! Initialize
    CALL CopyMatrix(Identity,D1)
    CALL ScaleMatrix(D1,beta_1)
    CALL CopyMatrix(Identity,TempMat)
    CALL ScaleMatrix(TempMat,mu)
    CALL IncrementMatrix(WorkingHamiltonian, TempMat, -1.0_NTREAL)
    CALL ScaleMatrix(TempMat,beta_2)
    CALL IncrementMatrix(TempMat,D1)
    trace_value = 0.0_NTREAL

    !! Iterate
    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Iterations")
       CALL EnterSubLog
    END IF
    outer_counter = 1
    norm_value = solver_parameters%converge_diff + 1.0_NTREAL
    norm_value2 = norm_value
    energy_value = 0.0_NTREAL

    DO outer_counter = 1,solver_parameters%max_iterations
       !! Compute the hole matrix DH
       CALL CopyMatrix(D1,DH)
       CALL IncrementMatrix(Identity,DH,alpha_in=-1.0_NTREAL)
       CALL ScaleMatrix(DH,-1.0_NTREAL)

       !! Compute DDH, as well as convergence check
       CALL MatrixMultiply(D1,DH,DDH,threshold_in=solver_parameters%threshold,&
            & memory_pool_in=pool)
       CALL MatrixTrace(DDH, trace_value)
       norm_value = ABS(trace_value)

       !! Compute D2DH
       CALL MatrixMultiply(D1,DDH,D2DH, &
            & threshold_in=solver_parameters%threshold, &
            & memory_pool_in=pool)

       !! Compute Sigma
       CALL MatrixTrace(D2DH, sigma_array(outer_counter))
       sigma_array(outer_counter) = sigma_array(outer_counter)/trace_value

       CALL CopyMatrix(D1,TempMat)

       !! Compute D1 + 2*D2DH
       CALL IncrementMatrix(D2DH,D1,alpha_in=2.0_NTREAL)

       !! Compute D1 + 2*D2DH -2*Sigma*DDH
       CALL IncrementMatrix(DDH, D1, &
            & alpha_in=-1.0_NTREAL*2.0_NTREAL*sigma_array(outer_counter))

       !! Energy value based convergence
       energy_value2 = energy_value
       CALL DotMatrix(D1,WorkingHamiltonian,energy_value)
       energy_value = 2.0_NTREAL*energy_value
       norm_value = ABS(energy_value - energy_value2)

       IF (solver_parameters%be_verbose) THEN
          CALL WriteListElement(key="Round", value=outer_counter)
          CALL EnterSubLog
          CALL WriteElement(key="Convergence", value=norm_value)
          CALL WriteElement("Energy_Value", value=energy_value)
          CALL ExitSubLog
       END IF

       IF (norm_value .LE. solver_parameters%converge_diff) THEN
          EXIT
       END IF
    END DO
    total_iterations = outer_counter-1
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
       CALL WriteElement(key="Total_Iterations", value=outer_counter)
       CALL PrintMatrixInformation(D1)
    END IF

    IF (PRESENT(energy_value_out)) THEN
       energy_value_out = energy_value
    END IF

    !! Undo Load Balancing Step
    IF (solver_parameters%do_load_balancing) THEN
       CALL UndoPermuteMatrix(D1, D1, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
    END IF

    !! Compute the density matrix in the non-orthogonalized basis
    CALL MatrixMultiply(InverseSquareRoot,D1,TempMat, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
    CALL MatrixMultiply(TempMat,InverseSquareRoot,Density, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool)

    !! Cleanup
    CALL DestructMatrix(WorkingHamiltonian)
    CALL DestructMatrix(TempMat)
    CALL DestructMatrix(D1)
    CALL DestructMatrix(DH)
    CALL DestructMatrix(DDH)
    CALL DestructMatrix(D2DH)
    CALL DestructMatrix(Identity)
    CALL DestructMatrixMemoryPool(pool)

    !! Compute The Chemical Potential
    IF (PRESENT(chemical_potential_out)) THEN
       interval_a = 0.0_NTREAL
       interval_b = 1.0_NTREAL
       midpoint = 0.0_NTREAL
       midpoints: DO outer_counter = 1,solver_parameters%max_iterations
          midpoint = (interval_b - interval_a)/2.0_NTREAL + interval_a
          zero_value = midpoint
          !! Compute polynomial function at the guess point.
          polynomial: DO inner_counter=1,total_iterations
             zero_value = zero_value + &
                  & 2.0_NTREAL*((zero_value**2)*(1.0_NTREAL-zero_value) &
                  & - sigma_array(inner_counter)* &
                  & zero_value*(1.0_NTREAL-zero_value))
          END DO polynomial
          !! Change bracketing.
          IF (zero_value .LT. 0.5_NTREAL) THEN
             interval_a = midpoint
          ELSE
             interval_b = midpoint
          END IF
          !! Check convergence.
          IF (ABS(zero_value-0.5_NTREAL) .LT. &
               & solver_parameters%converge_diff) THEN
             EXIT
          END IF
       END DO midpoints
       !! Undo scaling.
       chemical_potential_out = mu + (beta_1 - midpoint)/beta_2
    END IF
    !! Cleanup
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
    END IF
    DEALLOCATE(sigma_array)
    CALL DestructSolverParameters(solver_parameters)
  END SUBROUTINE HPCP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the density matrix from a Hamiltonian using the HPCP+ method.
  SUBROUTINE HPCPPlus(Hamiltonian, InverseSquareRoot, nel, Density, &
       & energy_value_out, chemical_potential_out, solver_parameters_in)
    !! Parameters
    TYPE(Matrix_ps), INTENT(IN) :: Hamiltonian
    TYPE(Matrix_ps), INTENT(IN) :: InverseSquareRoot
    INTEGER, INTENT(IN) :: nel
    TYPE(Matrix_ps), INTENT(INOUT) :: Density
    REAL(NTREAL), INTENT(OUT), OPTIONAL :: energy_value_out
    REAL(NTREAL), INTENT(OUT), OPTIONAL :: chemical_potential_out
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: &
         & solver_parameters_in
    !! Handling Optional Parameters
    TYPE(SolverParameters_t) :: solver_parameters
    !! Local Matrices
    TYPE(Matrix_ps) :: WorkingHamiltonian
    TYPE(Matrix_ps) :: TempMat
    TYPE(Matrix_ps) :: Identity
    TYPE(Matrix_ps) :: D1, DH, DDH, D2DH
    !! Local Variables
    REAL(NTREAL) :: e_min, e_max
    REAL(NTREAL) :: beta_1, beta_2
    REAL(NTREAL) :: beta_1_h, beta_2_h
    REAL(NTREAL) :: a, b, c, d
    REAL(NTREAL) :: one_third
    REAL(NTREAL) :: mixing_interior, mixing_value
    REAL(NTREAL) :: beta, beta_bar
    REAL(NTREAL) :: sigma, sigma_bar
    REAL(NTREAL) :: mu
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: sigma_array
    REAL(NTREAL) :: trace_value
    REAL(NTREAL) :: norm_value
    REAL(NTREAL) :: energy_value, energy_value2
    !! For computing the chemical potential
    REAL(NTREAL) :: zero_value, midpoint, interval_a, interval_b
    !! Temporary Variables
    TYPE(MatrixMemoryPool_p) :: pool
    INTEGER :: outer_counter, inner_counter
    INTEGER :: total_iterations
    INTEGER :: matrix_dimension

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = SolverParameters_init()
    END IF

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Density Matrix Solver")
       CALL EnterSubLog
       CALL WriteElement(key="Method", value="HPCP+")
       CALL WriteCitation("truflandier2016communication")
       CALL PrintParameters(solver_parameters)
    END IF

    ALLOCATE(sigma_array(solver_parameters%max_iterations))

    matrix_dimension = Hamiltonian%actual_matrix_dimension

    !! Construct All The Necessary Matrices
    CALL ConstructEmptyMatrix(Density, Hamiltonian)
    CALL ConstructEmptyMatrix(WorkingHamiltonian, Hamiltonian)
    CALL ConstructEmptyMatrix(TempMat, Hamiltonian)
    CALL ConstructEmptyMatrix(D1, Hamiltonian)
    CALL ConstructEmptyMatrix(DH, Hamiltonian)
    CALL ConstructEmptyMatrix(DDH, Hamiltonian)
    CALL ConstructEmptyMatrix(D2DH, Hamiltonian)
    CALL ConstructEmptyMatrix(Identity, Hamiltonian)
    CALL FillMatrixIdentity(Identity)

    !! Compute the working hamiltonian.
    CALL MatrixMultiply(InverseSquareRoot,Hamiltonian,TempMat, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
    CALL MatrixMultiply(TempMat,InverseSquareRoot,WorkingHamiltonian, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool)

    !! Load Balancing Step
    IF (solver_parameters%do_load_balancing) THEN
       CALL PermuteMatrix(WorkingHamiltonian, WorkingHamiltonian, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
       CALL PermuteMatrix(Identity, Identity, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
    END IF

    !! Compute the initial matrix.
    CALL GershgorinBounds(WorkingHamiltonian,e_min,e_max)
    CALL MatrixTrace(WorkingHamiltonian, mu)
    mu = mu/matrix_dimension
    sigma_bar = (matrix_dimension - 0.5_NTREAL*nel)/matrix_dimension
    sigma = 1.0_NTREAL - sigma_bar
    beta = sigma/((e_max) - mu)
    beta_bar = sigma_bar/(mu - (e_min))
    beta_1 = sigma
    beta_1_h = sigma_bar
    beta_2 = MIN(beta,beta_bar)
    beta_2_h = -1.0_NTREAL*MAX(beta,beta_bar)

    !! Initialize
    CALL CopyMatrix(Identity,D1)
    CALL ScaleMatrix(D1,beta_1)
    CALL CopyMatrix(Identity,TempMat)
    CALL ScaleMatrix(TempMat,mu)
    CALL IncrementMatrix(WorkingHamiltonian, TempMat, &
         & -1.0_NTREAL)
    CALL ScaleMatrix(TempMat,beta_2)
    CALL IncrementMatrix(TempMat,D1)

    CALL CopyMatrix(Identity,DH)
    CALL ScaleMatrix(DH,beta_1_h)
    CALL CopyMatrix(Identity,TempMat)
    CALL ScaleMatrix(TempMat,mu)
    CALL IncrementMatrix(WorkingHamiltonian,TempMat,-1.0_NTREAL)
    CALL ScaleMatrix(TempMat,beta_2_h)
    CALL IncrementMatrix(TempMat,DH)

    CALL CopyMatrix(Identity,TempMat)
    CALL IncrementMatrix(DH,TempMat,-1.0_NTREAL)

    CALL DotMatrix(D1,D1,a)
    CALL DotMatrix(TempMat,TempMat,b)
    CALL DotMatrix(D1,TempMat,c)

    one_third = 1.0_NTREAL/3.0_NTREAL
    IF (sigma < one_third) THEN
       d = (nel*0.5_NTREAL) - 2.0_NTREAL*one_third*(nel*0.5_NTREAL)
    ELSE
       d = (nel*0.5_NTREAL) - &
            & 2.0_NTREAL*one_third*(matrix_dimension - nel*0.5_NTREAL)
    END IF

    mixing_interior = SQRT((2.0_NTREAL*c-2.0_NTREAL*b)**2 &
         & - 4.0_NTREAL*(a+b-2.0_NTREAL*c)*(b-d))
    mixing_value = ((2.0_NTREAL*b - 2.0_NTREAL*c) + mixing_interior) &
         & / (2.0_NTREAL*(a+b - 2.0_NTREAL*c))
    IF (.NOT. mixing_value .LE. 1.0_NTREAL &
         & .AND. mixing_value .GE. 0.0_NTREAL) THEN
       mixing_value = ((2.0_NTREAL*b - 2.0_NTREAL*c) - mixing_interior) &
            & / (2.0_NTREAL*(a+b - 2.0_NTREAL*c))
    ENDIF

    CALL ScaleMatrix(D1,mixing_value)
    CALL CopyMatrix(Identity,TempMat)
    CALL IncrementMatrix(DH,TempMat,-1.0_NTREAL)
    CALL IncrementMatrix(TempMat,D1,1.0_NTREAL-mixing_value)

    !! Iterate
    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Iterations")
       CALL EnterSubLog
    END IF
    outer_counter = 1
    norm_value = solver_parameters%converge_diff + 1.0_NTREAL
    energy_value = 0.0_NTREAL
    DO outer_counter = 1,solver_parameters%max_iterations
       !! Compute the hole matrix DH
       CALL CopyMatrix(D1,DH)
       CALL IncrementMatrix(Identity,DH,alpha_in=-1.0_NTREAL)
       CALL ScaleMatrix(DH,-1.0_NTREAL)

       !! Compute DDH, as well as convergence check
       CALL MatrixMultiply(D1,DH,DDH,threshold_in=solver_parameters%threshold,&
            & memory_pool_in=pool)
       CALL MatrixTrace(DDH, trace_value)
       norm_value = ABS(trace_value)

       !! Compute D2DH
       CALL MatrixMultiply(D1,DDH,D2DH, &
            & threshold_in=solver_parameters%threshold, memory_pool_in=pool)

       !! Compute Sigma
       CALL MatrixTrace(D2DH, sigma_array(outer_counter))
       sigma_array(outer_counter) = sigma_array(outer_counter)/trace_value

       CALL CopyMatrix(D1,TempMat)

       !! Compute D1 + 2*D2DH
       CALL IncrementMatrix(D2DH,D1,alpha_in=2.0_NTREAL)

       !! Compute D1 + 2*D2DH -2*Sigma*DDH
       CALL IncrementMatrix(DDH, D1, &
            & alpha_in=-1.0_NTREAL*2.0_NTREAL*sigma_array(outer_counter))

       !! Energy value based convergence
       energy_value2 = energy_value
       CALL DotMatrix(D1,WorkingHamiltonian,energy_value)
       energy_value = 2.0_NTREAL*energy_value
       norm_value = ABS(energy_value - energy_value2)

       IF (solver_parameters%be_verbose) THEN
          CALL WriteListElement(key="Round", value=outer_counter)
          CALL EnterSubLog
          CALL WriteElement(key="Convergence", value=norm_value)
          CALL WriteElement("Energy_Value", value=energy_value)
          CALL ExitSubLog
       END IF

       IF (norm_value .LE. solver_parameters%converge_diff) THEN
          EXIT
       END IF
    END DO
    total_iterations = outer_counter-1
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
       CALL WriteElement(key="Total_Iterations", value=outer_counter)
       CALL PrintMatrixInformation(D1)
    END IF

    IF (PRESENT(energy_value_out)) THEN
       energy_value_out = energy_value
    END IF

    !! Undo Load Balancing Step
    IF (solver_parameters%do_load_balancing) THEN
       CALL UndoPermuteMatrix(D1, D1, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
    END IF

    !! Compute the density matrix in the non-orthogonalized basis
    CALL MatrixMultiply(InverseSquareRoot,D1,TempMat, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
    CALL MatrixMultiply(TempMat,InverseSquareRoot,Density, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool)

    !! Cleanup
    CALL DestructMatrix(WorkingHamiltonian)
    CALL DestructMatrix(TempMat)
    CALL DestructMatrix(D1)
    CALL DestructMatrix(DH)
    CALL DestructMatrix(DDH)
    CALL DestructMatrix(D2DH)
    CALL DestructMatrix(Identity)
    CALL DestructMatrixMemoryPool(pool)

    !! Compute The Chemical Potential
    IF (PRESENT(chemical_potential_out)) THEN
       interval_a = 0.0_NTREAL
       interval_b = 1.0_NTREAL
       midpoint = 0.0_NTREAL
       midpoints: DO outer_counter = 1,solver_parameters%max_iterations
          midpoint = (interval_b - interval_a)/2.0_NTREAL + interval_a
          zero_value = midpoint
          !! Compute polynomial function at the guess point.
          polynomial: DO inner_counter=1,total_iterations
             zero_value = zero_value + &
                  & 2.0_NTREAL*((zero_value*zero_value)*(1.0_NTREAL-zero_value)&
                  & - sigma_array(inner_counter)* &
                  & (zero_value*(1.0_NTREAL-zero_value)))
          END DO polynomial
          !! Change bracketing.
          IF (zero_value .LT. 0.5_NTREAL) THEN
             interval_a = midpoint
          ELSE
             interval_b = midpoint
          END IF
          !! Check convergence.
          IF (ABS(zero_value-0.5_NTREAL) .LT. &
               & solver_parameters%converge_diff) THEN
             EXIT
          END IF
       END DO midpoints
       !! Undo scaling.
       chemical_potential_out = mu + (beta_1 - midpoint)/beta_2
    END IF
    !! Cleanup
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
    END IF
    DEALLOCATE(sigma_array)
    CALL DestructSolverParameters(solver_parameters)
  END SUBROUTINE HPCPPlus
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the density matrix from a Hamiltonian using the Scale and Fold
  !> method. Based on the method of \cite rubensson2011nonmonotonic .
  !> Note that for this method, you must provide the value of the homo and
  !> lumo gap. It is not necessary for these to be accurate, but give a
  !> conservative value.
  SUBROUTINE ScaleAndFold(Hamiltonian, InverseSquareRoot, nel, Density, &
       & homo, lumo, energy_value_out, solver_parameters_in)
    !> The matrix to compute the corresponding density from
    TYPE(Matrix_ps), INTENT(IN) :: Hamiltonian
    !> The inverse square root of the overlap matrix.
    TYPE(Matrix_ps), INTENT(IN) :: InverseSquareRoot
    !> The number of electrons.
    INTEGER, INTENT(IN) :: nel
    !> The density matrix computed by this routine.
    TYPE(Matrix_ps), INTENT(INOUT) :: Density
    !> The energy of the system (optional).
    REAL(NTREAL), INTENT(OUT), OPTIONAL :: energy_value_out
    !> A conservative estimate of the highest occupied eigenvalue.
    REAL(NTREAL), INTENT(IN) :: homo
    !> A conservative estimate of the lowest unoccupied eigenvalue.
    REAL(NTREAL), INTENT(IN) :: lumo
    !> Parameters for the solver (optional).
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Optional Parameters
    TYPE(SolverParameters_t) :: solver_parameters
    !! Local Matrices
    TYPE(Matrix_ps) :: WorkingHamiltonian
    TYPE(Matrix_ps) :: Identity
    TYPE(Matrix_ps) :: X_k, X_k2, TempMat
    !! Local Variables
    REAL(NTREAL) :: e_min, e_max
    REAL(NTREAL) :: Beta, BetaBar, alpha
    REAL(NTREAL) :: trace_value
    REAL(NTREAL) :: norm_value
    REAL(NTREAL) :: energy_value, energy_value2
    !! Temporary Variables
    TYPE(MatrixMemoryPool_p) :: pool
    INTEGER :: outer_counter
    INTEGER :: total_iterations

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = SolverParameters_init()
    END IF

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Density Matrix Solver")
       CALL EnterSubLog
       CALL WriteElement(key="Method", value="Scale and Fold")
       CALL WriteCitation("rubensson2011nonmonotonic")
       CALL PrintParameters(solver_parameters)
    END IF

    !! Construct All The Necessary Matrices
    CALL ConstructEmptyMatrix(Density, Hamiltonian)
    CALL ConstructEmptyMatrix(WorkingHamiltonian, Hamiltonian)
    CALL ConstructEmptyMatrix(X_k, Hamiltonian)
    CALL ConstructEmptyMatrix(X_k2, Hamiltonian)
    CALL ConstructEmptyMatrix(TempMat, Hamiltonian)
    CALL ConstructEmptyMatrix(Identity, Hamiltonian)
    CALL FillMatrixIdentity(Identity)

    !! Compute the working hamiltonian.
    CALL MatrixMultiply(InverseSquareRoot,Hamiltonian,TempMat, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
    CALL MatrixMultiply(TempMat,InverseSquareRoot,WorkingHamiltonian, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool)

    !! Load Balancing Step
    IF (solver_parameters%do_load_balancing) THEN
       CALL PermuteMatrix(WorkingHamiltonian, WorkingHamiltonian, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
       CALL PermuteMatrix(Identity, Identity, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
    END IF

    !! Compute the lambda scaling value.
    CALL GershgorinBounds(WorkingHamiltonian,e_min,e_max)

    !! Initialize
    CALL CopyMatrix(WorkingHamiltonian,X_k)
    CALL ScaleMatrix(X_k,-1.0_NTREAL)
    CALL IncrementMatrix(Identity,X_k,alpha_in=e_max)
    CALL ScaleMatrix(X_k,1.0_NTREAL/(e_max-e_min))
    Beta = (e_max - lumo) / (e_max - e_min)
    BetaBar = (e_max - homo) / (e_max - e_min)

    !! Iterate
    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Iterations")
       CALL EnterSubLog
    END IF
    outer_counter = 1
    norm_value = solver_parameters%converge_diff + 1.0_NTREAL
    energy_value = 0.0_NTREAL
    DO outer_counter = 1,solver_parameters%max_iterations
       !! Determine the path
       CALL MatrixTrace(X_k, trace_value)
       IF (trace_value .GT. 0.5_NTREAL*nel) THEN
          alpha = 2.0_NTREAL/(2.0_NTREAL - Beta)
          CALL ScaleMatrix(X_k, alpha)
          CALL IncrementMatrix(Identity, X_k, alpha_in=(1.0_NTREAL-alpha))
          CALL MatrixMultiply(X_k,X_k,X_k2, &
               & threshold_in=solver_parameters%threshold, &
               & memory_pool_in=pool)
          CALL CopyMatrix(X_k2, X_k)
          Beta = (alpha * Beta + 1 - alpha)**2
          BetaBar = (alpha * BetaBar + 1 - alpha)**2
       ELSE
          alpha = 2.0/(1.0 + BetaBar)
          CALL MatrixMultiply(X_k,X_k,X_k2, &
               & threshold_in=solver_parameters%threshold, &
               & memory_pool_in=pool)
          CALL ScaleMatrix(X_k, 2*alpha)
          CALL IncrementMatrix(X_k2, X_k, alpha_in=-1.0_NTREAL*alpha**2)
          Beta = 2.0_NTREAL * alpha * Beta - alpha**2 * Beta**2
          BetaBar = 2.0_NTREAL * alpha * BetaBar - alpha**2 * BetaBar ** 2
       END IF

       !! Energy value based convergence
       energy_value2 = energy_value
       CALL DotMatrix(X_k, WorkingHamiltonian, energy_value)
       energy_value = 2.0_NTREAL*energy_value
       norm_value = ABS(energy_value - energy_value2)

       IF (solver_parameters%be_verbose) THEN
          CALL WriteListElement(key="Round", value=outer_counter)
          CALL EnterSubLog
          CALL WriteElement(key="Convergence", value=norm_value)
          CALL WriteElement("Energy_Value", value=energy_value)
          CALL ExitSubLog
       END IF

       IF (norm_value .LE. solver_parameters%converge_diff) THEN
          EXIT
       END IF
    END DO
    total_iterations = outer_counter-1
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
       CALL WriteElement(key="Total_Iterations", value=outer_counter)
       CALL PrintMatrixInformation(X_k)
    END IF

    IF (PRESENT(energy_value_out)) THEN
       energy_value_out = energy_value
    END IF

    !! Undo Load Balancing Step
    IF (solver_parameters%do_load_balancing) THEN
       CALL UndoPermuteMatrix(X_k, X_k, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
    END IF

    !! Compute the density matrix in the non-orthogonalized basis
    CALL MatrixMultiply(InverseSquareRoot,X_k,TempMat, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
    CALL MatrixMultiply(TempMat,InverseSquareRoot,Density, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool)

    !! Cleanup
    CALL DestructMatrix(WorkingHamiltonian)
    CALL DestructMatrix(X_k)
    CALL DestructMatrix(X_k2)
    CALL DestructMatrix(TempMat)
    CALL DestructMatrix(Identity)
    CALL DestructMatrixMemoryPool(pool)

    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
    END IF

    !! Cleanup
    CALL DestructSolverParameters(solver_parameters)
  END SUBROUTINE SCALEANDFOLD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the energy-weighted density matrix.
  SUBROUTINE EnergyDensityMatrix(Hamiltonian, Density, EnergyDensity, &
       & threshold_in)
    !> The matrix to compute from.
    TYPE(Matrix_ps), INTENT(IN) :: Hamiltonian
    !> The density matrix.
    TYPE(Matrix_ps), INTENT(IN) :: Density
    !> The energy-weighted density matrix.
    TYPE(Matrix_ps), INTENT(INOUT) :: EnergyDensity
    !> Threshold for flushing small values (default = 0).
    REAL(NTREAL), INTENT(IN), OPTIONAL :: threshold_in
    !! Handling Optional Parameters
    REAL(NTREAL) :: threshold
    !! Local Matrices
    TYPE(Matrix_ps) :: TempMat
    !! Temporary Variables
    TYPE(MatrixMemoryPool_p) :: pool

    !! Optional Parameters
    IF (PRESENT(threshold_in)) THEN
       threshold = threshold_in
    ELSE
       threshold = 0.0_NTREAL
    END IF

    !! Construct All The Necessary Matrices
    CALL ConstructEmptyMatrix(EnergyDensity, Hamiltonian)
    CALL ConstructEmptyMatrix(TempMat, Hamiltonian)

    !! EDM = DM * H * DM
    CALL MatrixMultiply(Density, Hamiltonian, TempMat, &
         & threshold_in=threshold, memory_pool_in=pool)
    CALL MatrixMultiply(TempMat, Density, EnergyDensity, &
         & threshold_in=threshold, memory_pool_in=pool)

    !! Cleanup
    CALL DestructMatrix(TempMat)
    CALL DestructMatrixMemoryPool(pool)
  END SUBROUTINE EnergyDensityMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE DensityMatrixSolversModule
