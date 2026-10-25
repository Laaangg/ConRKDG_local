module Solver2D
  use PDE2D
  use Prob2D
  use Mesh2D
  use WENO
  use GaussQuadrature1D
  implicit none

  type :: Solver2D_t
    ! Time integration and flux settings
    character(16) :: typeTime    = 'RK3'
    character(16) :: typeNumFlux = 'LxF'
    character(16) :: typeCFL     = 'NumCFL'
    real(8)       :: CFL         = 0.5_8
    real(8)       :: t_final     = 0.0_8
    real(8)       :: dt          = 0.0_8
    real(8)       :: t           = 0.0_8
    integer       :: nt          = 0
    integer       :: step_count  = 0

    ! Solution state and spatial operator RHS (n_sub^2 x numCells)
    real(8), allocatable :: U(:,:)
    real(8), allocatable :: rhs(:,:)

    ! WENO limiter and troubled-cell indicator
    type(WENO_t)         :: we
    character(16)        :: typeWENO        = 'none'
    logical              :: isLimiterActive = .false.
    integer, allocatable :: isTroubled(:)
    integer              :: totalTroubled   = 0
    real(8)              :: M_TVB           = 50.0_8

    ! Error metrics
    real(8) :: L1err   = 0.0_8
    real(8) :: L2err   = 0.0_8
    real(8) :: Linferr = 0.0_8

    ! ANSYS-style per-step file output settings (Choice A)
    logical        :: enableOutput   = .true.
    character(16)  :: outputFormat   = 'BINARY'  ! 'BINARY' (.bin) or 'ASCII' (.txt)
    character(16)  :: saveMode       = 'FPS'     ! 'FPS' (frame interval) or 'STEP' (every N steps)
    real(8)        :: saveFPS        = 20.0_8    ! Target frames per second of simulation time (e.g. 20 FPS)
    real(8)        :: frameInterval  = 0.05_8    ! dt between saved frames = 1.0 / saveFPS
    real(8)        :: nextFrameTime  = 0.0_8     ! Next target time to save a frame
    integer        :: frameCount     = 0         ! Number of frames saved
    integer        :: saveEvery      = 1         ! Step stride if saveMode == 'STEP'
    logical        :: saveCurrent    = .true.    ! Overwrite sol_current at EACH time step
    character(256) :: outputDir      = ''
    integer        :: monitorUnit    = 0
  end type Solver2D_t

contains

  ! -------------------------------------------------------------------
  ! Constructor: initialize 2D ConRKDG Solver (Solver.m, Lines 43-100)
  ! -------------------------------------------------------------------
  subroutine init_Solver2D(solver, prob, mesh, typeTime, typeNumFlux, typeWENO, CFL, &
                           saveEvery, enableOutput, outputFormat, saveMode, saveFPS, &
                           frameInterval, saveCurrent)
    type(Solver2D_t), intent(inout) :: solver
    type(Prob2D_t),   intent(in)    :: prob
    type(Mesh2D_t),   intent(in)    :: mesh
    character(*),     intent(in)    :: typeTime, typeNumFlux, typeWENO
    real(8),          intent(in)    :: CFL
    integer, optional, intent(in)   :: saveEvery
    logical, optional, intent(in)   :: enableOutput
    character(*), optional, intent(in) :: outputFormat
    character(*), optional, intent(in) :: saveMode
    real(8), optional, intent(in)      :: saveFPS
    real(8), optional, intent(in)      :: frameInterval
    logical, optional, intent(in)      :: saveCurrent

    integer :: numSub2, numCells

    solver%typeTime    = trim(typeTime)
    solver%typeNumFlux = trim(typeNumFlux)
    solver%t_final     = prob%finalTime
    solver%t           = 0.0_8
    solver%step_count  = 0
    solver%nt          = 0
    solver%frameCount  = 0
    solver%nextFrameTime = 0.0_8

    if (CFL <= 0.0_8) then
      solver%typeCFL = 'Exact'
      solver%CFL     = 0.5_8
    else
      solver%typeCFL = 'NumCFL'
      solver%CFL     = CFL
    end if

    numSub2  = mesh%n_sub * mesh%n_sub
    numCells = mesh%numCells

    allocate(solver%U(numSub2, numCells))
    allocate(solver%rhs(numSub2, numCells))
    allocate(solver%isTroubled(numCells))

    solver%U             = 0.0_8
    solver%rhs           = 0.0_8
    solver%isTroubled    = 0
    solver%totalTroubled = 0

    if (trim(typeWENO) /= 'none' .and. trim(typeWENO) /= 'off' .and. len_trim(typeWENO) > 0) then
      solver%isLimiterActive = .true.
      solver%typeWENO        = trim(typeWENO)
      call init_WENO(solver%we, trim(typeWENO))
    else
      solver%isLimiterActive = .false.
      solver%typeWENO        = 'none'
    end if

    if (present(saveEvery)) then
      solver%saveEvery = max(1, saveEvery)
    else
      solver%saveEvery = 1
    end if

    if (present(enableOutput)) then
      solver%enableOutput = enableOutput
    else
      solver%enableOutput = .true.
    end if

    if (present(outputFormat)) then
      solver%outputFormat = trim(outputFormat)
    else
      solver%outputFormat = 'BINARY'
    end if

    if (present(saveMode)) then
      solver%saveMode = trim(saveMode)
    else
      solver%saveMode = 'FPS'
    end if

    if (present(saveFPS)) then
      solver%saveFPS = saveFPS
    else
      solver%saveFPS = 20.0_8
    end if

    if (present(frameInterval)) then
      if (frameInterval > 0.0_8) then
        solver%frameInterval = frameInterval
      else if (solver%saveFPS > 0.0_8) then
        solver%frameInterval = 1.0_8 / solver%saveFPS
      else
        solver%frameInterval = 0.05_8
      end if
    else if (solver%saveFPS > 0.0_8) then
      solver%frameInterval = 1.0_8 / solver%saveFPS
    else
      solver%frameInterval = 0.05_8
    end if

    if (present(saveCurrent)) then
      solver%saveCurrent = saveCurrent
    else
      solver%saveCurrent = .true.
    end if

    if (solver%enableOutput) then
      call initOutputDirectory(solver, prob, mesh)
    end if
  end subroutine init_Solver2D

  ! -------------------------------------------------------------------
  ! Evaluate Initial Conditions via 2D Sub-Cell Gauss Quadrature
  ! Matches Solver.setIC (Solver.m, Lines 103-166; Paper Sec. 2.3 Eq. 13)
  ! -------------------------------------------------------------------
  subroutine setIC(solver, prob, mesh, n_quad_pts)
    type(Solver2D_t),  intent(inout) :: solver
    type(Prob2D_t),    intent(in)    :: prob
    type(Mesh2D_t),    intent(inout) :: mesh
    integer, optional, intent(in)    :: n_quad_pts

    integer :: nG, c, sCell_x, sCell_y, sub_id, qx, qy
    real(8) :: x_start, y_start, x0, y0, u_sum, u_exact, weight_2D
    real(8), allocatable :: xq(:), yq(:)
    type(GaussQuadrature1D_t) :: GQ_IC

    if (present(n_quad_pts)) then
      nG = n_quad_pts
    else
      nG = 5
    end if

    call init_GaussQuadrature1D(GQ_IC, nG)
    allocate(xq(nG), yq(nG))

    do c = 1, mesh%numCells
      x_start = mesh%cells(c)%Centroid(1) - 0.5_8 * mesh%dx
      y_start = mesh%cells(c)%Centroid(2) - 0.5_8 * mesh%dy

      do sCell_y = 1, mesh%n_sub
        do sCell_x = 1, mesh%n_sub
          sub_id = (sCell_y - 1) * mesh%n_sub + sCell_x

          x0 = x_start + real(sCell_x - 1, 8) * mesh%dx_sub
          y0 = y_start + real(sCell_y - 1, 8) * mesh%dy_sub

          xq = 0.5_8 * mesh%dx_sub * GQ_IC%xGauss + (x0 + 0.5_8 * mesh%dx_sub)
          yq = 0.5_8 * mesh%dy_sub * GQ_IC%xGauss + (y0 + 0.5_8 * mesh%dy_sub)

          mesh%cells(c)%SubCellCentroids(sub_id, 1) = x_start + (real(sCell_x, 8) - 0.5_8) * mesh%dx_sub
          mesh%cells(c)%SubCellCentroids(sub_id, 2) = y_start + (real(sCell_y, 8) - 0.5_8) * mesh%dy_sub

          u_sum = 0.0_8
          do qy = 1, nG
            do qx = 1, nG
              u_exact   = eval_IC(prob, xq(qx), yq(qy))
              weight_2D = 0.25_8 * GQ_IC%wGauss(qx) * GQ_IC%wGauss(qy)
              u_sum     = u_sum + u_exact * weight_2D
            end do
          end do

          solver%U(sub_id, c) = u_sum
        end do
      end do
    end do

    deallocate(xq, yq, GQ_IC%xGauss, GQ_IC%wGauss)
  end subroutine setIC

  ! -------------------------------------------------------------------
  ! Compute Time Step dt from CFL Condition
  ! Matches Solver.computeTimeStep (Solver.m, Lines 169-195; Paper Eqs. 51-52)
  ! -------------------------------------------------------------------
  subroutine computeTimeStep(solver, prob, mesh, dt_out)
    type(Solver2D_t), intent(inout) :: solver
    type(Prob2D_t),   intent(in)    :: prob
    type(Mesh2D_t),   intent(in)    :: mesh
    real(8),          intent(out)   :: dt_out

    real(8) :: inv_dt, dt_cfl

    inv_dt = (prob%maxSpeedX / mesh%dx_sub) + (prob%maxSpeedY / mesh%dy_sub)
    dt_cfl = 1.0_8 / inv_dt

    if (trim(solver%typeCFL) == 'Exact' .and. prob%withExact == 1 .and. mesh%n_sub > 3) then
      dt_out = solver%CFL * (dt_cfl)**(real(mesh%n_sub, 8) / 3.0_8)
    else
      dt_out = solver%CFL * dt_cfl
    end if

    ! Clip dt if saving by frame interval to hit next frame time exactly
    if (trim(solver%saveMode) == 'FPS' .and. solver%enableOutput) then
      if (solver%nextFrameTime > solver%t + 1.0e-14_8 .and. &
          solver%t + dt_out > solver%nextFrameTime) then
        dt_out = solver%nextFrameTime - solver%t
      end if
    end if

    if (solver%t + dt_out > solver%t_final) then
      dt_out = solver%t_final - solver%t
    end if
    solver%dt = dt_out
  end subroutine computeTimeStep

  ! -------------------------------------------------------------------
  ! Global Lax-Friedrichs Numerical Flux in X and Y (Paper Eq. 11)
  ! -------------------------------------------------------------------
  pure elemental real(8) function computeNumFluxX(prob, uL, uR) result(numFluxX)
    type(Prob2D_t), intent(in) :: prob
    real(8),        intent(in) :: uL, uR
    real(8) :: fL, fR
    fL       = analyticalFluxX(prob%pde, uL)
    fR       = analyticalFluxX(prob%pde, uR)
    numFluxX = 0.5_8 * (fL + fR - prob%maxSpeedX * (uR - uL))
  end function computeNumFluxX

  pure elemental real(8) function computeNumFluxY(prob, uB, uT) result(numFluxY)
    type(Prob2D_t), intent(in) :: prob
    real(8),        intent(in) :: uB, uT
    real(8) :: gB, gT
    gB       = analyticalFluxY(prob%pde, uB)
    gT       = analyticalFluxY(prob%pde, uT)
    numFluxY = 0.5_8 * (gB + gT - prob%maxSpeedY * (uT - uB))
  end function computeNumFluxY

  ! -------------------------------------------------------------------
  ! Spatial DG Operator L(U) / computeRHS
  ! Matches Solver.computeRHS (Solver.m, Lines 216-363; Paper Eqs. 23-27)
  ! -------------------------------------------------------------------
  subroutine computeRHS(solver, U_state, prob, mesh, dUdt)
    type(Solver2D_t), intent(inout) :: solver
    real(8),          intent(in)    :: U_state(:,:)
    type(Prob2D_t),   intent(in)    :: prob
    type(Mesh2D_t),   intent(in)    :: mesh
    real(8),          intent(out)   :: dUdt(size(U_state, 1), size(U_state, 2))

    integer :: num_faces, num_cells, num_sub, n1D_sub, num_quad
    integer :: f, c, i, j, Cell_Minus_ID, Cell_Plus_ID
    real(8) :: nx, ny, dxi_dx, deta_dy, J_k
    real(8), allocatable :: edge_flux_integral(:,:), vol_flux_integral(:,:)
    real(8) :: W_1D_Y(mesh%n_quad), W_1D_X(mesh%n_quad)
    real(8) :: W_2D(mesh%n_quad, mesh%n_quad)
    real(8) :: U_Minus_mat(mesh%n_sub, mesh%n_sub), U_Plus_mat(mesh%n_sub, mesh%n_sub)
    real(8) :: U_mat(mesh%n_sub, mesh%n_sub)
    real(8) :: U_Minus_edge(mesh%n_quad), U_Plus_edge(mesh%n_quad)
    real(8) :: F_L(mesh%n_quad), F_R(mesh%n_quad), G_B(mesh%n_quad), G_T(mesh%n_quad)
    real(8) :: Flux_q(mesh%n_quad)
    real(8) :: Flux_Y_proj(mesh%n_sub), Flux_X_proj(mesh%n_sub)
    real(8) :: Flux_mat_L(mesh%n_sub, mesh%n_sub), Flux_mat_R(mesh%n_sub, mesh%n_sub)
    real(8) :: Flux_mat_B(mesh%n_sub, mesh%n_sub), Flux_mat_T(mesh%n_sub, mesh%n_sub)
    real(8) :: U_q(mesh%n_quad, mesh%n_quad), F_qw(mesh%n_quad, mesh%n_quad), G_qw(mesh%n_quad, mesh%n_quad)
    real(8) :: VolX_mat(mesh%n_sub, mesh%n_sub), VolY_mat(mesh%n_sub, mesh%n_sub)
    real(8) :: total_spatial_flux(mesh%n_sub * mesh%n_sub)

    num_faces = mesh%numFaces
    num_cells = mesh%numCells
    n1D_sub   = mesh%n_sub
    num_sub   = n1D_sub * n1D_sub
    num_quad  = mesh%n_quad

    allocate(edge_flux_integral(num_sub, num_cells))
    allocate(vol_flux_integral(num_sub, num_cells))

    dUdt               = 0.0_8
    edge_flux_integral = 0.0_8
    vol_flux_integral  = 0.0_8

    W_1D_Y = mesh%wq * mesh%dy_sub
    W_1D_X = mesh%wq * mesh%dx_sub

    ! 1. Loop over all faces (interface fluxes, Eq. 11 & 2nd term of Eq. 23)
    do f = 1, num_faces
      nx = mesh%faces(f)%Normal(1)
      ny = mesh%faces(f)%Normal(2)

      if (nx > 0.5_8 .or. ny > 0.5_8) then
        Cell_Minus_ID = mesh%faces(f)%LeftCellID
        Cell_Plus_ID  = mesh%faces(f)%RightCellID
      else
        Cell_Minus_ID = mesh%faces(f)%RightCellID
        Cell_Plus_ID  = mesh%faces(f)%LeftCellID
      end if

      if (abs(nx) > 0.5_8) then
        ! Vertical face (flux in +X direction)
        if (Cell_Minus_ID == 0) then
          U_Plus_mat   = reshape(U_state(:, Cell_Plus_ID), [n1D_sub, n1D_sub])
          U_Minus_edge = matmul(matmul(mesh%H_bound_L, U_Plus_mat), mesh%H_quad)
        else
          U_Minus_mat  = reshape(U_state(:, Cell_Minus_ID), [n1D_sub, n1D_sub])
          U_Minus_edge = matmul(matmul(mesh%H_bound_R, U_Minus_mat), mesh%H_quad)
        end if

        if (Cell_Plus_ID == 0) then
          U_Minus_mat = reshape(U_state(:, Cell_Minus_ID), [n1D_sub, n1D_sub])
          U_Plus_edge = matmul(matmul(mesh%H_bound_R, U_Minus_mat), mesh%H_quad)
        else
          U_Plus_mat  = reshape(U_state(:, Cell_Plus_ID), [n1D_sub, n1D_sub])
          U_Plus_edge = matmul(matmul(mesh%H_bound_L, U_Plus_mat), mesh%H_quad)
        end if

        ! Global Lax-Friedrichs Numerical Flux Eq. (11)
        F_L    = analyticalFluxX(prob%pde, U_Minus_edge)
        F_R    = analyticalFluxX(prob%pde, U_Plus_edge)
        Flux_q = 0.5_8 * (F_L + F_R - prob%maxSpeedX * (U_Plus_edge - U_Minus_edge))

        Flux_Y_proj = matmul(mesh%H_quad, Flux_q * W_1D_Y)

        if (Cell_Minus_ID /= 0) then
          do j = 1, n1D_sub
            do i = 1, n1D_sub
              Flux_mat_L(i, j) = mesh%H_bound_R(i) * Flux_Y_proj(j)
            end do
          end do
          edge_flux_integral(:, Cell_Minus_ID) = edge_flux_integral(:, Cell_Minus_ID) + &
                                                 reshape(Flux_mat_L, [num_sub])
        end if

        if (Cell_Plus_ID /= 0) then
          do j = 1, n1D_sub
            do i = 1, n1D_sub
              Flux_mat_R(i, j) = mesh%H_bound_L(i) * Flux_Y_proj(j)
            end do
          end do
          edge_flux_integral(:, Cell_Plus_ID) = edge_flux_integral(:, Cell_Plus_ID) - &
                                                reshape(Flux_mat_R, [num_sub])
        end if

      else if (abs(ny) > 0.5_8) then
        ! Horizontal face (flux in +Y direction)
        if (Cell_Minus_ID == 0) then
          U_Plus_mat   = reshape(U_state(:, Cell_Plus_ID), [n1D_sub, n1D_sub])
          U_Minus_edge = matmul(transpose(mesh%H_quad), matmul(U_Plus_mat, mesh%H_bound_L))
        else
          U_Minus_mat  = reshape(U_state(:, Cell_Minus_ID), [n1D_sub, n1D_sub])
          U_Minus_edge = matmul(transpose(mesh%H_quad), matmul(U_Minus_mat, mesh%H_bound_R))
        end if

        if (Cell_Plus_ID == 0) then
          U_Minus_mat = reshape(U_state(:, Cell_Minus_ID), [n1D_sub, n1D_sub])
          U_Plus_edge = matmul(transpose(mesh%H_quad), matmul(U_Minus_mat, mesh%H_bound_R))
        else
          U_Plus_mat  = reshape(U_state(:, Cell_Plus_ID), [n1D_sub, n1D_sub])
          U_Plus_edge = matmul(transpose(mesh%H_quad), matmul(U_Plus_mat, mesh%H_bound_L))
        end if

        ! Global Lax-Friedrichs Numerical Flux Eq. (11)
        G_B    = analyticalFluxY(prob%pde, U_Minus_edge)
        G_T    = analyticalFluxY(prob%pde, U_Plus_edge)
        Flux_q = 0.5_8 * (G_B + G_T - prob%maxSpeedY * (U_Plus_edge - U_Minus_edge))

        Flux_X_proj = matmul(mesh%H_quad, Flux_q * W_1D_X)

        if (Cell_Minus_ID /= 0) then
          do j = 1, n1D_sub
            do i = 1, n1D_sub
              Flux_mat_B(i, j) = Flux_X_proj(i) * mesh%H_bound_R(j)
            end do
          end do
          edge_flux_integral(:, Cell_Minus_ID) = edge_flux_integral(:, Cell_Minus_ID) + &
                                                 reshape(Flux_mat_B, [num_sub])
        end if

        if (Cell_Plus_ID /= 0) then
          do j = 1, n1D_sub
            do i = 1, n1D_sub
              Flux_mat_T(i, j) = Flux_X_proj(i) * mesh%H_bound_L(j)
            end do
          end do
          edge_flux_integral(:, Cell_Plus_ID) = edge_flux_integral(:, Cell_Plus_ID) - &
                                                reshape(Flux_mat_T, [num_sub])
        end if
      end if
    end do

    ! 2. Loop over all cells (volume integrals & final RHS, Eqs. 23-27)
    dxi_dx  = 1.0_8 / mesh%dx_sub
    deta_dy = 1.0_8 / mesh%dy_sub
    J_k     = mesh%dx_sub * mesh%dy_sub

    do j = 1, num_quad
      do i = 1, num_quad
        W_2D(i, j) = mesh%wq(i) * mesh%wq(j)
      end do
    end do

    do c = 1, num_cells
      U_mat = reshape(U_state(:, c), [n1D_sub, n1D_sub])

      U_q  = matmul(transpose(mesh%H_quad), matmul(U_mat, mesh%H_quad))
      F_qw = analyticalFluxX(prob%pde, U_q) * W_2D
      G_qw = analyticalFluxY(prob%pde, U_q) * W_2D

      VolX_mat = matmul(mesh%dH_quad, matmul(F_qw, transpose(mesh%H_quad)))  * (dxi_dx  * J_k)
      VolY_mat = matmul(mesh%H_quad,  matmul(G_qw, transpose(mesh%dH_quad))) * (deta_dy * J_k)

      vol_flux_integral(:, c) = reshape(VolX_mat + VolY_mat, [num_sub])

      total_spatial_flux = vol_flux_integral(:, c) - edge_flux_integral(:, c)
      dUdt(:, c)         = matmul(mesh%S_inv_2D, total_spatial_flux) / J_k
    end do

    solver%rhs = dUdt
    deallocate(edge_flux_integral, vol_flux_integral)
  end subroutine computeRHS

  ! -------------------------------------------------------------------
  ! Apply WENO Limiter to Troubled Cells
  ! Matches Solver.applyWENOLimiter (Solver.m, Lines 366-446)
  ! -------------------------------------------------------------------
  subroutine applyWENOLimiter(solver, U_state, mesh, U_limited)
    type(Solver2D_t), intent(inout) :: solver
    real(8),          intent(in)    :: U_state(:,:)
    type(Mesh2D_t),   intent(inout) :: mesh
    real(8),          intent(out)   :: U_limited(size(U_state, 1), size(U_state, 2))

    integer :: Nx, Ny, num_sub, stencil_size, left_nbr, right_nbr
    integer :: i, j, k, cID, idx_X, idx_Y
    integer, allocatable :: offsets(:)
    real(8), allocatable :: U_macro_avg(:), stencil_X(:), stencil_Y(:)
    real(8) :: avg_C, internal_diff, neighbor_diff_X, neighbor_diff_Y, threshold
    logical :: isTrouble

    U_limited = U_state
    if (.not. solver%isLimiterActive) return

    Nx           = mesh%Nx
    Ny           = mesh%Ny
    num_sub      = size(U_state, 1)
    stencil_size = solver%we%sizeBigStencil

    allocate(offsets(stencil_size))
    allocate(stencil_X(stencil_size), stencil_Y(stencil_size))
    allocate(U_macro_avg(mesh%numCells))

    if (stencil_size == 5) then
      offsets   = [-2, -1, 0, 1, 2]
      left_nbr  = 2
      right_nbr = 4
    else
      offsets   = [-2, -1, 0, 1, 2, 3]
      left_nbr  = 2
      right_nbr = 4
    end if

    ! 1. Calculate Macro-Cell Averages
    do cID = 1, mesh%numCells
      U_macro_avg(cID) = sum(U_state(:, cID)) / real(num_sub, 8)
    end do

    solver%isTroubled = 0
    threshold         = solver%M_TVB * (max(mesh%dx, mesh%dy)**2)

    ! 2. Loop over all cells to detect and limit shocks
    do j = 1, Ny
      do i = 1, Nx
        cID   = i + (j - 1) * Nx
        avg_C = U_macro_avg(cID)

        ! A. Fetch Stencil (Periodic Wrapping)
        do k = 1, stencil_size
          idx_X        = modulo(i + offsets(k) - 1, Nx) + 1
          stencil_X(k) = U_macro_avg(idx_X + (j - 1) * Nx)

          idx_Y        = modulo(j + offsets(k) - 1, Ny) + 1
          stencil_Y(k) = U_macro_avg(i + (idx_Y - 1) * Nx)
        end do

        ! B. Troubled-Cell Indicator (TCI)
        internal_diff   = maxval(U_state(:, cID)) - minval(U_state(:, cID))
        neighbor_diff_X = max(abs(avg_C - stencil_X(left_nbr)), abs(avg_C - stencil_X(right_nbr)))
        neighbor_diff_Y = max(abs(avg_C - stencil_Y(left_nbr)), abs(avg_C - stencil_Y(right_nbr)))

        isTrouble = (internal_diff > 1.5_8 * max(neighbor_diff_X, neighbor_diff_Y)) .and. &
                    (internal_diff > threshold)

        ! C. WENO Reconstruction / Fallback
        if (isTrouble) then
          mesh%cells(cID)%IsTrouble = .true.
          solver%isTroubled(cID)    = 1

          call computeNonlinearWeights(solver%we, stencil_X, mesh%dx)
          call computeNonlinearWeights(solver%we, stencil_Y, mesh%dy)

          U_limited(:, cID) = avg_C
        else
          mesh%cells(cID)%IsTrouble = .false.
        end if
      end do
    end do

    solver%totalTroubled = sum(solver%isTroubled)
    deallocate(offsets, stencil_X, stencil_Y, U_macro_avg)
  end subroutine applyWENOLimiter

  ! -------------------------------------------------------------------
  ! Main Solve Method: 3-Stage Third-Order SSP-RK3 Time Integration
  ! Matches Solver.solve (Solver.m, Lines 449-536; Paper Eq. 28)
  ! Also writes per-step output files (Choice A) when enableOutput=.true.
  ! -------------------------------------------------------------------
  subroutine solve(solver, prob, mesh)
    type(Solver2D_t), intent(inout) :: solver
    type(Prob2D_t),   intent(in)    :: prob
    type(Mesh2D_t),   intent(inout) :: mesh

    integer :: num_sub, num_cells
    real(8) :: dt_step
    real(8), allocatable :: U_curr(:,:), U1(:,:), U2(:,:), U_new(:,:)
    real(8), allocatable :: L1(:,:), L2(:,:), L3(:,:)

    num_sub   = mesh%n_sub * mesh%n_sub
    num_cells = mesh%numCells

    allocate(U_curr(num_sub, num_cells), U1(num_sub, num_cells))
    allocate(U2(num_sub, num_cells),     U_new(num_sub, num_cells))
    allocate(L1(num_sub, num_cells),     L2(num_sub, num_cells), L3(num_sub, num_cells))

    solver%t          = 0.0_8
    solver%step_count = 0
    solver%frameCount = 0
    U_curr            = solver%U

    write(*,'(A)') 'Starting 2D ConRKDG SSP-RK3 Integration...'
    if (solver%isLimiterActive) then
      write(*,'(A,A)') 'WENO Limiter enabled: ', trim(solver%typeWENO)
    else
      write(*,'(A)') 'WENO Limiter disabled.'
    end if

    if (solver%enableOutput) then
      if (trim(solver%saveMode) == 'FPS') then
        write(*,'(A,F6.2,A,F8.4,A)') 'Export Mode: FPS (Target FPS: ', solver%saveFPS, &
                                    ', Frame Interval dt: ', solver%frameInterval, ')'
      else
        write(*,'(A,I0,A)') 'Export Mode: STEP (Every ', solver%saveEvery, ' steps)'
      end if
      if (solver%saveCurrent) then
        write(*,'(A)') 'At-The-Moment Checkpoint: enabled (sol_current overwritten each step)'
      end if
    end if

    ! Compute initial dt preview
    call computeTimeStep(solver, prob, mesh, dt_step)

    ! Write initial snapshot at t = 0
    if (solver%enableOutput) then
      if (trim(solver%saveMode) == 'FPS') then
        call writeStepSolution(solver, mesh, 0, solver%t, isFrame=.true.)
        solver%frameCount    = 1
        solver%nextFrameTime = solver%frameInterval
      else
        call writeStepSolution(solver, mesh, 0, solver%t, isFrame=.false.)
      end if

      ! Write initial sol_current checkpoint
      if (solver%saveCurrent) then
        call writeCurrentState(solver, mesh)
      end if
    end if

    do while (solver%t < solver%t_final - 1.0e-14_8)
      ! 1. Compute time step dt (Eqs. 51-52)
      call computeTimeStep(solver, prob, mesh, dt_step)

      ! 2. SSP-RK3 Stages (Eq. 28)
      ! --- Stage 1 ---
      call computeRHS(solver, U_curr, prob, mesh, L1)
      U1 = U_curr + dt_step * L1
      if (solver%isLimiterActive) then
        call applyWENOLimiter(solver, U1, mesh, U_new)
        U1 = U_new
      end if

      ! --- Stage 2 ---
      call computeRHS(solver, U1, prob, mesh, L2)
      U2 = 0.75_8 * U_curr + 0.25_8 * U1 + 0.25_8 * dt_step * L2
      if (solver%isLimiterActive) then
        call applyWENOLimiter(solver, U2, mesh, U_new)
        U2 = U_new
      end if

      ! --- Stage 3 ---
      call computeRHS(solver, U2, prob, mesh, L3)
      U_new = (1.0_8 / 3.0_8) * U_curr + (2.0_8 / 3.0_8) * U2 + (2.0_8 / 3.0_8) * dt_step * L3
      if (solver%isLimiterActive) then
        call applyWENOLimiter(solver, U_new, mesh, U_curr)
      else
        U_curr = U_new
      end if

      ! 3. Update global state, time, and step counter
      solver%U          = U_curr
      solver%t          = solver%t + dt_step
      solver%step_count = solver%step_count + 1

      ! 4. At-The-Moment Checkpoint: Overwrite sol_current at EACH time step
      if (solver%enableOutput .and. solver%saveCurrent) then
        call writeCurrentState(solver, mesh)
      end if

      ! 5. Write frame snapshot or step snapshot
      if (solver%enableOutput) then
        if (trim(solver%saveMode) == 'FPS') then
          if (solver%t >= solver%nextFrameTime - 1.0e-12_8 .or. &
              abs(solver%t - solver%t_final) <= 1.0e-12_8) then
            call writeStepSolution(solver, mesh, solver%frameCount, solver%t, isFrame=.true.)
            solver%frameCount    = solver%frameCount + 1
            solver%nextFrameTime = min(solver%t_final, solver%nextFrameTime + solver%frameInterval)
          end if
        else
          if (mod(solver%step_count, solver%saveEvery) == 0 .or. &
              abs(solver%t - solver%t_final) <= 1.0e-12_8) then
            call writeStepSolution(solver, mesh, solver%step_count, solver%t, isFrame=.false.)
          end if
        end if
      end if

      ! 6. Print progress to terminal
      if (mod(solver%step_count, 10) == 0 .or. abs(solver%t - solver%t_final) <= 1.0e-12_8) then
        write(*,'(A,F10.4,A,F10.4,A,I0,A,I0,A)') &
          '  Time: ', solver%t, ' / ', solver%t_final, &
          ' (Step: ', solver%step_count, ', Troubled: ', solver%totalTroubled, ')'
      end if
    end do

    solver%nt = solver%step_count
    if (solver%enableOutput .and. solver%monitorUnit /= 0) then
      close(solver%monitorUnit)
      solver%monitorUnit = 0
    end if

    write(*,'(A,I0,A)') 'Simulation Complete! Total steps: ', solver%nt, '.'
    deallocate(U_curr, U1, U2, U_new, L1, L2, L3)
  end subroutine solve

  ! -------------------------------------------------------------------
  ! Map sub-cell state array U(n_sub^2, numCells) to 2D grid
  ! Matches Solver.mapStateToGrid (Solver.m, Lines 539-559)
  ! -------------------------------------------------------------------
  subroutine mapStateToGrid(U_state, mesh, U_grid)
    real(8),        intent(in)  :: U_state(:,:)
    type(Mesh2D_t), intent(in)  :: mesh
    real(8),        intent(out) :: U_grid(mesh%Ny * mesh%n_sub, mesh%Nx * mesh%n_sub)

    integer :: i, j, sx, sy, cID
    real(8) :: U_mat(mesh%n_sub, mesh%n_sub)

    do j = 1, mesh%Ny
      do i = 1, mesh%Nx
        cID   = i + (j - 1) * mesh%Nx
        U_mat = reshape(U_state(:, cID), [mesh%n_sub, mesh%n_sub])
        do sy = 1, mesh%n_sub
          do sx = 1, mesh%n_sub
            U_grid((j - 1) * mesh%n_sub + sy, (i - 1) * mesh%n_sub + sx) = U_mat(sx, sy)
          end do
        end do
      end do
    end do
  end subroutine mapStateToGrid

  ! -------------------------------------------------------------------
  ! Compute L1, L2 (Eq. 54), and Linf (Eq. 55) Errors at Current Time
  ! Matches Solver.computeErrors (Solver.m, Lines 562-594)
  ! -------------------------------------------------------------------
  subroutine computeErrors(solver, prob, mesh)
    type(Solver2D_t), intent(inout) :: solver
    type(Prob2D_t),   intent(in)    :: prob
    type(Mesh2D_t),   intent(in)    :: mesh

    integer :: ny_grid, nx_grid
    real(8) :: sub_area
    real(8), allocatable :: X_plot(:,:), Y_plot(:,:), U_plot(:,:), U_exact(:,:), diff_U(:,:)

    if (prob%withExact /= 1) return

    ny_grid = mesh%Ny * mesh%n_sub
    nx_grid = mesh%Nx * mesh%n_sub
    allocate(X_plot(ny_grid, nx_grid), Y_plot(ny_grid, nx_grid))
    allocate(U_plot(ny_grid, nx_grid), U_exact(ny_grid, nx_grid), diff_U(ny_grid, nx_grid))

    call getSpatialGrid(mesh, X_plot, Y_plot)
    call mapStateToGrid(solver%U, mesh, U_plot)

    U_exact  = eval_exact(prob, X_plot, Y_plot, solver%t)
    diff_U   = abs(U_plot - U_exact)
    sub_area = mesh%dx_sub * mesh%dy_sub

    solver%L1err   = sum(diff_U) * sub_area
    solver%L2err   = sqrt(sum(diff_U * diff_U) * sub_area)
    solver%Linferr = maxval(diff_U)

    deallocate(X_plot, Y_plot, U_plot, U_exact, diff_U)
  end subroutine computeErrors

  ! ===================================================================
  ! ANSYS-STYLE OUTPUT SUBROUTINES (CHOICE A: ELEMENT TABLE PER STEP)
  ! ===================================================================

  ! -------------------------------------------------------------------
  ! Initialize output subfolder, write static mesh.txt once, and open
  ! monitor.txt
  ! -------------------------------------------------------------------
  subroutine initOutputDirectory(solver, prob, mesh)
    type(Solver2D_t), intent(inout) :: solver
    type(Prob2D_t),   intent(in)    :: prob
    type(Mesh2D_t),   intent(in)    :: mesh

    character(256) :: folderCmd, meshFile, monitorFile
    integer        :: fid, i, j, sx, sy
    real(8), allocatable :: x_sub_1d(:), y_sub_1d(:)

    write(solver%outputDir, '(A,A,A,I0,A,I0,A,I0,A,A)') &
      'results/', trim(prob%shortName), '_nx', mesh%Nx, '_ny', mesh%Ny, &
      '_m', mesh%m_degree, '_', trim(solver%typeWENO)

    ! Create results/<run_dir> on Windows/Linux
    folderCmd = 'powershell -NoProfile -Command "New-Item -ItemType Directory -Force -Path ''' // &
                trim(solver%outputDir) // ''' | Out-Null"'
    call execute_command_line(trim(folderCmd))

    ! 1. Write static mesh.txt once at t = 0
    meshFile = trim(solver%outputDir) // '/mesh.txt'
    open(newunit=fid, file=trim(meshFile), status='replace', action='write')
    write(fid, '(A)') '# ConRKDG 2D Static Mesh & Sub-Cell Coordinates'
    write(fid, '(A)') '# Nx Ny m_degree n_sub N_sub_total xLeft xRight yLeft yRight dx dy dx_sub dy_sub'
    write(fid, '(5(I0,1X),8(ES22.15,1X))') &
      mesh%Nx, mesh%Ny, mesh%m_degree, mesh%n_sub, mesh%n_sub * mesh%n_sub, &
      mesh%xLeft, mesh%xRight, mesh%yLeft, mesh%yRight, &
      mesh%dx, mesh%dy, mesh%dx_sub, mesh%dy_sub

    allocate(x_sub_1d(mesh%Nx * mesh%n_sub), y_sub_1d(mesh%Ny * mesh%n_sub))
    do i = 1, mesh%Nx
      do sx = 1, mesh%n_sub
        x_sub_1d((i - 1) * mesh%n_sub + sx) = &
          mesh%xLeft + real(i - 1, 8) * mesh%dx + (real(sx, 8) - 0.5_8) * mesh%dx_sub
      end do
    end do
    do j = 1, mesh%Ny
      do sy = 1, mesh%n_sub
        y_sub_1d((j - 1) * mesh%n_sub + sy) = &
          mesh%yLeft + real(j - 1, 8) * mesh%dy + (real(sy, 8) - 0.5_8) * mesh%dy_sub
      end do
    end do

    write(fid, '(A,I0,A)') '# X_sub_centers (', mesh%Nx * mesh%n_sub, ' values)'
    write(fid, '(*(ES22.15,1X))') x_sub_1d
    write(fid, '(A,I0,A)') '# Y_sub_centers (', mesh%Ny * mesh%n_sub, ' values)'
    write(fid, '(*(ES22.15,1X))') y_sub_1d
    close(fid)
    deallocate(x_sub_1d, y_sub_1d)

    ! 2. Open monitor.txt for live per-step logging
    monitorFile = trim(solver%outputDir) // '/monitor.txt'
    open(newunit=solver%monitorUnit, file=trim(monitorFile), status='replace', action='write')
    write(solver%monitorUnit, '(A,A,A,I0,A,I0,A,I0,A,A,A,F6.3)') &
      '# Problem: ', trim(prob%typeProb), ' | Nx=', mesh%Nx, ' Ny=', mesh%Ny, &
      ' m=', mesh%m_degree, ' | WENO=', trim(solver%typeWENO), ' | CFL=', solver%CFL
    write(solver%monitorUnit, '(A)') &
      '# Step  Time                    dt                      TroubledCells  TotalCells  Troubled(%)             SnapshotFile'
    flush(solver%monitorUnit)
  end subroutine initOutputDirectory

  ! -------------------------------------------------------------------
  ! Write Snapshot (Frame or Step Mode)
  ! and append row to monitor.txt
  ! -------------------------------------------------------------------
  subroutine writeStepSolution(solver, mesh, snap_idx, current_time, isFrame)
    type(Solver2D_t), intent(in) :: solver
    type(Mesh2D_t),   intent(in) :: mesh
    integer,          intent(in) :: snap_idx
    real(8),          intent(in) :: current_time
    logical, optional, intent(in):: isFrame

    character(64)  :: stepFileShort
    character(256) :: stepFileFull
    integer        :: fid, c, num_sub2
    real(8)        :: troubled_pct
    logical        :: useFrame

    num_sub2 = mesh%n_sub * mesh%n_sub
    useFrame = .false.
    if (present(isFrame)) useFrame = isFrame

    if (trim(solver%outputFormat) == 'BINARY') then
      ! 1. Binary Stream format (.bin) — fastest I/O, minimal size
      if (useFrame) then
        write(stepFileShort, '(A,I5.5,A)') 'sol_frame_', snap_idx, '.bin'
      else
        write(stepFileShort, '(A,I5.5,A)') 'sol_step_', snap_idx, '.bin'
      end if
      stepFileFull = trim(solver%outputDir) // '/' // trim(stepFileShort)

      open(newunit=fid, file=trim(stepFileFull), access='stream', form='unformatted', &
           status='replace', action='write')
      write(fid) solver%isTroubled(1:mesh%numCells)
      write(fid) solver%U(1:num_sub2, 1:mesh%numCells)
      close(fid)
    else
      ! 2. ASCII text format (.txt) / (.dat) — human-readable Choice A table
      if (useFrame) then
        ! write(stepFileShort, '(A,I5.5,A)') 'sol_frame_', snap_idx, '.txt'
        write(stepFileShort, '(A,I5.5,A)') 'sol_frame_', snap_idx, '.dat'
      else
        ! write(stepFileShort, '(A,I5.5,A)') 'sol_step_', snap_idx, '.txt'
        write(stepFileShort, '(A,I5.5,A)') 'sol_step_', snap_idx, '.dat'
      end if
      stepFileFull = trim(solver%outputDir) // '/' // trim(stepFileShort)

      open(newunit=fid, file=trim(stepFileFull), status='replace', action='write')
      write(fid, '(A,I0,A,ES22.15,A,ES22.15,A,I0,A,I0,A,I0)') &
        '# Step: ', solver%step_count, ' | Time: ', current_time, ' | dt: ', solver%dt, &
        ' | Nx=', mesh%Nx, ' Ny=', mesh%Ny, ' n_sub=', mesh%n_sub
      write(fid, '(A,I0,A)') &
        '# CellID  ix  iy  isTroubled  U_sub(1..', num_sub2, ')'

      do c = 1, mesh%numCells
        write(fid, '(4(I6,1X),*(ES22.15,1X))') &
          mesh%cells(c)%ID, mesh%cells(c)%ix, mesh%cells(c)%iy, &
          solver%isTroubled(c), solver%U(1:num_sub2, c)
      end do
      close(fid)
    end if

    ! Update monitor.txt
    if (solver%monitorUnit /= 0) then
      troubled_pct = 100.0_8 * real(solver%totalTroubled, 8) / real(mesh%numCells, 8)
      write(solver%monitorUnit, '(I6,2X,ES22.15,2X,ES22.15,2X,I13,2X,I10,2X,ES22.15,2X,A)') &
        solver%step_count, current_time, solver%dt, solver%totalTroubled, mesh%numCells, &
        troubled_pct, trim(stepFileShort)
      flush(solver%monitorUnit)
    end if
  end subroutine writeStepSolution

  ! -------------------------------------------------------------------
  ! At-The-Moment Checkpoint: Overwrite sol_current at EACH time step
  ! Just in case the simulation is halted or interrupted
  ! -------------------------------------------------------------------
  subroutine writeCurrentState(solver, mesh)
    type(Solver2D_t), intent(in) :: solver
    type(Mesh2D_t),   intent(in) :: mesh

    character(256) :: curFile, curInfoFile
    integer        :: fid, c, num_sub2

    num_sub2 = mesh%n_sub * mesh%n_sub

    if (trim(solver%outputFormat) == 'BINARY') then
      curFile = trim(solver%outputDir) // '/sol_current.bin'
      open(newunit=fid, file=trim(curFile), access='stream', form='unformatted', &
           status='replace', action='write')
      write(fid) solver%isTroubled(1:mesh%numCells)
      write(fid) solver%U(1:num_sub2, 1:mesh%numCells)
      close(fid)
    else
      curFile = trim(solver%outputDir) // '/sol_current.txt'
      open(newunit=fid, file=trim(curFile), status='replace', action='write')
      write(fid, '(A,I0,A,ES22.15,A,ES22.15,A,I0,A,I0,A,I0)') &
        '# Step: ', solver%step_count, ' | Time: ', solver%t, ' | dt: ', solver%dt, &
        ' | Nx=', mesh%Nx, ' Ny=', mesh%Ny, ' n_sub=', mesh%n_sub
      write(fid, '(A,I0,A)') &
        '# CellID  ix  iy  isTroubled  U_sub(1..', num_sub2, ')'
      do c = 1, mesh%numCells
        write(fid, '(4(I6,1X),*(ES22.15,1X))') &
          mesh%cells(c)%ID, mesh%cells(c)%ix, mesh%cells(c)%iy, &
          solver%isTroubled(c), solver%U(1:num_sub2, c)
      end do
      close(fid)
    end if

    ! Write small human-readable metadata info file
    curInfoFile = trim(solver%outputDir) // '/current_info.txt'
    open(newunit=fid, file=trim(curInfoFile), status='replace', action='write')
    write(fid, '(A)') '# ConRKDG At-The-Moment Checkpoint Info'
    write(fid, '(A,I0)')      'step: ', solver%step_count
    write(fid, '(A,ES22.15)') 'time: ', solver%t
    write(fid, '(A,ES22.15)') 'dt: ', solver%dt
    write(fid, '(A,I0)')      'troubled_cells: ', solver%totalTroubled
    write(fid, '(A,I0)')      'total_cells: ', mesh%numCells
    write(fid, '(A,A)')       'file: ', trim(curFile)
    close(fid)
  end subroutine writeCurrentState

  ! -------------------------------------------------------------------
  ! Compute and write final L1, L2, Linf errors to errors.txt
  ! -------------------------------------------------------------------
  subroutine writeFinalErrors(solver, prob, mesh)
    type(Solver2D_t), intent(inout) :: solver
    type(Prob2D_t),   intent(in)    :: prob
    type(Mesh2D_t),   intent(in)    :: mesh

    character(256) :: errFile
    integer        :: fid

    call computeErrors(solver, prob, mesh)

    write(*,'(A)') '------------------------------------------------------------'
    write(*,'(A,ES16.8)') '  L1   Error : ', solver%L1err
    write(*,'(A,ES16.8)') '  L2   Error : ', solver%L2err
    write(*,'(A,ES16.8)') '  Linf Error : ', solver%Linferr
    write(*,'(A)') '------------------------------------------------------------'

    if (solver%enableOutput) then
      errFile = trim(solver%outputDir) // '/errors.txt'
      open(newunit=fid, file=trim(errFile), status='replace', action='write')
      write(fid, '(A)') '# Nx Ny m_degree n_sub nt L1_Error L2_Error Linf_Error'
      write(fid, '(5(I0,1X),3(ES22.15,1X))') &
        mesh%Nx, mesh%Ny, mesh%m_degree, mesh%n_sub, solver%nt, &
        solver%L1err, solver%L2err, solver%Linferr
      close(fid)
      write(*,'(A,A)') 'Results & step snapshots written to: ', trim(solver%outputDir)
    end if
  end subroutine writeFinalErrors

  ! -------------------------------------------------------------------
  ! Deallocate Solver2D arrays
  ! -------------------------------------------------------------------
  subroutine free_Solver2D(solver)
    type(Solver2D_t), intent(inout) :: solver
    if (allocated(solver%U))          deallocate(solver%U)
    if (allocated(solver%rhs))        deallocate(solver%rhs)
    if (allocated(solver%isTroubled)) deallocate(solver%isTroubled)
  end subroutine free_Solver2D

end module Solver2D
