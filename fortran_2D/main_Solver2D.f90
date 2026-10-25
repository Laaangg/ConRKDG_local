! ===================================================================
! 2D Conservative RKDG Method with WENO Limiter — Main Driver
! Translated from MATLAB_prototype_refactor/test_main.m
! ===================================================================
program main_Solver2D
  use Prob2D
  use Mesh2D
  use Solver2D
  implicit none

  type(Prob2D_t)   :: prob
  type(Mesh2D_t)   :: mesh
  type(Solver2D_t) :: solver

  integer       :: Nx, Ny, m_degree, save_every
  real(8)       :: CFL, t_start_cpu, t_end_cpu
  character(64) :: probName
  character(16) :: typeTime, typeNumFlux, typeWENO, outputFormat, save_mode
  real(8)       :: save_fps, frame_interval
  logical       :: enable_output, save_current


  integer, dimension(5) :: grid_array  = [10, 20, 40, 80, 160]
  integer :: curr_grid

  ! -----------------------------------------------------------------
  ! 1. Simulation Parameters (Matches test_main.m)
  ! -----------------------------------------------------------------
  probName = 'Scalar_Linear_Advection_Smooth'
  ! probName = 'Scalar_Linear_Advection_Smooth_2'
  ! probName = 'Scalar_Linear_Advection_WavePackage'
  ! probName = 'Scalar_Linear_Advection_Shock_Diagonal'
  ! probName = 'Scalar_Linear_Advection_Quadrant'
  ! probName = 'Scalar_Inviscid_Burgers'

  do curr_grid = lbound(grid_array, 1), ubound(grid_array, 1)
    Nx = grid_array(curr_grid)
    Ny = grid_array(curr_grid)

    write(*,*) 'Running simulation for grid size: ', Nx, ' x ', Ny

    ! Initialize problem and mesh
  call init_Prob2D(prob, trim(probName))

  ! Nx       = 20
  ! Ny       = 20
  m_degree = 4
  call init_Mesh2D(mesh, Nx, Ny, m_degree, prob)

  CFL         = 0.5_8
  typeTime    = 'RK3'
  typeNumFlux = 'LxF'
  ! Options: 'WENO6theta', 'WENO5z', 'WENO6z', 'WENO6cz', 'none'
  typeWENO    = 'WENO6z'

  ! -----------------------------------------------------------------
  ! 2. Output Settings (Frame Interval / FPS & At-The-Moment Checkpoint)
  ! -----------------------------------------------------------------
  enable_output  = .true.
  ! outputFormat   = 'BINARY'    ! 'BINARY' (.bin stream, ~70% smaller) or 'ASCII' (.txt)
  outputFormat   = 'ASCII'
  save_mode      = 'FPS'       ! 'FPS' (frame interval) or 'STEP' (every N steps)
  save_fps       = 50.0_8      ! When save_mode='FPS': frames per unit simulation time (dt_frame = 1.0 / save_fps)
  frame_interval = 0.0_8       ! If > 0, overrides save_fps with explicit dt_frame (e.g. 0.05)
  save_every     = 1           ! When save_mode='STEP': snapshot every N global steps
  save_current   = .true.      ! Overwrite sol_current at EACH time step (checkpoint if halted)

  ! -----------------------------------------------------------------
  ! 3. Initialize Solver and Evaluate Initial Conditions
  ! -----------------------------------------------------------------
  call cpu_time(t_start_cpu)

  call init_Solver2D(solver, prob, mesh, typeTime, typeNumFlux, typeWENO, &
                     CFL, save_every, enable_output, outputFormat, &
                     save_mode, save_fps, frame_interval, save_current)
  call setIC(solver, prob, mesh)

  ! -----------------------------------------------------------------
  ! 4. Solve for Solution at Final Time (SSP-RK3)
  ! -----------------------------------------------------------------
  call solve(solver, prob, mesh)

  ! -----------------------------------------------------------------
  ! 5. Verification & Error Calculation
  ! -----------------------------------------------------------------
  call writeFinalErrors(solver, prob, mesh)

  call cpu_time(t_end_cpu)
  write(*,'(A,F10.4,A)') 'Total Elapsed CPU Time: ', t_end_cpu - t_start_cpu, ' s'

  ! Clean up allocated memory
  call free_Solver2D(solver)
  call free_Mesh2D(mesh)
  end do

end program main_Solver2D
