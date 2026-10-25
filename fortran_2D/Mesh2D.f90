module Mesh2D
  use Prob2D
  use GaussQuadrature1D
  use LinAlg
  implicit none

  ! -------------------------------------------------------------------
  ! Macro-Cell Derived Type (Cell.m)
  ! -------------------------------------------------------------------
  type :: Cell_t
    integer              :: ID = 0
    integer              :: ix = 0
    integer              :: iy = 0
    real(8)              :: Area = 0.0_8
    real(8)              :: Centroid(2) = [0.0_8, 0.0_8]
    real(8), allocatable :: SubCellCentroids(:,:) ! (n_sub^2, 2)
    logical              :: IsGhost    = .false.
    logical              :: AtBoundary = .false.
    logical              :: IsTrouble  = .false.
  end type Cell_t

  ! -------------------------------------------------------------------
  ! Interface Face Derived Type (Face.m)
  ! LeftCellID / RightCellID = 0 indicates empty / boundary side
  ! -------------------------------------------------------------------
  type :: Face_t
    integer       :: ID          = 0
    real(8)       :: Normal(2)   = [0.0_8, 0.0_8]
    real(8)       :: Length      = 0.0_8
    integer       :: LeftCellID  = 0
    integer       :: RightCellID = 0
    character(16) :: BCType      = 'INTERIOR'
    real(8)       :: Flux(4)     = 0.0_8
  end type Face_t

  ! -------------------------------------------------------------------
  ! 2D Cartesian Mesh & Reference Element Math Derived Type (Mesh.m)
  ! -------------------------------------------------------------------
  type :: Mesh2D_t
    ! Domain bounds and lengths
    real(8) :: xLeft  = 0.0_8
    real(8) :: xRight = 0.0_8
    real(8) :: yLeft  = 0.0_8
    real(8) :: yRight = 0.0_8
    real(8) :: Lx     = 0.0_8
    real(8) :: Ly     = 0.0_8

    ! Macro-element grid parameters
    integer :: Nx       = 0
    integer :: Ny       = 0
    integer :: numCells = 0
    integer :: numFaces = 0
    real(8) :: dx       = 0.0_8
    real(8) :: dy       = 0.0_8

    ! Polynomial degree and sub-cell parameters
    integer :: m_degree = 0
    integer :: n_sub    = 0
    integer :: n_quad   = 5
    real(8) :: dx_sub   = 0.0_8
    real(8) :: dy_sub   = 0.0_8

    ! Topology arrays
    type(Cell_t), allocatable :: cells(:)
    type(Face_t), allocatable :: faces(:)

    ! Reference element math arrays (Sec. 2.3 Eqs. 16-22 & Eq. 26)
    real(8), allocatable :: alpha(:,:)     ! (n_sub, n_sub)
    real(8), allocatable :: S_1D(:,:)      ! (n_sub, n_sub)
    real(8), allocatable :: S_inv_2D(:,:)  ! (n_sub^2, n_sub^2)
    real(8), allocatable :: H_quad(:,:)    ! (n_sub, n_quad)
    real(8), allocatable :: dH_quad(:,:)   ! (n_sub, n_quad)
    real(8), allocatable :: H_bound_L(:)   ! (n_sub)
    real(8), allocatable :: H_bound_R(:)   ! (n_sub)
    real(8), allocatable :: xq(:)          ! (n_quad)
    real(8), allocatable :: wq(:)          ! (n_quad)
  end type Mesh2D_t

contains

  ! -------------------------------------------------------------------
  ! Constructor: initialize 2D Mesh & Reference Element Math
  ! -------------------------------------------------------------------
  subroutine init_Mesh2D(mesh, Nx, Ny, m_degree, prob)
    type(Mesh2D_t), intent(inout) :: mesh
    integer,        intent(in)    :: Nx, Ny, m_degree
    type(Prob2D_t), intent(in)    :: prob

    mesh%xLeft    = prob%xLeft
    mesh%xRight   = prob%xRight
    mesh%yLeft    = prob%yLeft
    mesh%yRight   = prob%yRight
    mesh%Lx       = mesh%xRight - mesh%xLeft
    mesh%Ly       = mesh%yRight - mesh%yLeft

    mesh%Nx       = Nx
    mesh%Ny       = Ny
    mesh%numCells = Nx * Ny
    mesh%dx       = mesh%Lx / real(Nx, 8)
    mesh%dy       = mesh%Ly / real(Ny, 8)

    mesh%m_degree = m_degree
    mesh%n_sub    = m_degree + 1
    mesh%n_quad   = 5
    mesh%dx_sub   = mesh%dx / real(mesh%n_sub, 8)
    mesh%dy_sub   = mesh%dy / real(mesh%n_sub, 8)

    call generateMesh(mesh, prob)
    call initReferenceMath(mesh)
  end subroutine init_Mesh2D

  ! -------------------------------------------------------------------
  ! Generate 2D cells, sub-cell centroids, and vertical/horizontal faces
  ! Matches Mesh.generateMesh (Mesh.m, Lines 75-204)
  ! -------------------------------------------------------------------
  subroutine generateMesh(mesh, prob)
    type(Mesh2D_t), intent(inout) :: mesh
    type(Prob2D_t), intent(in)    :: prob

    integer :: i, j, sCell_x, sCell_y, sub_id, ID_cell, fID
    integer :: numVFaces, numHFaces, cLeftID, cRightID
    real(8) :: x_start, y_start, len, normal(2)
    logical :: isPeriodic
    character(16) :: bc

    allocate(mesh%cells(mesh%numCells))

    ! 1. Create Cells and populate SubCellCentroids
    do j = 0, mesh%Ny - 1
      do i = 0, mesh%Nx - 1
        ID_cell = i + 1 + j * mesh%Nx
        mesh%cells(ID_cell)%ID          = ID_cell
        mesh%cells(ID_cell)%ix          = i + 1
        mesh%cells(ID_cell)%iy          = j + 1
        mesh%cells(ID_cell)%Area        = mesh%dx * mesh%dy
        mesh%cells(ID_cell)%Centroid(1) = 0.5_8 * mesh%dx + real(i, 8) * mesh%dx + mesh%xLeft
        mesh%cells(ID_cell)%Centroid(2) = 0.5_8 * mesh%dy + real(j, 8) * mesh%dy + mesh%yLeft
        mesh%cells(ID_cell)%IsTrouble   = .false.

        x_start = mesh%cells(ID_cell)%Centroid(1) - 0.5_8 * mesh%dx
        y_start = mesh%cells(ID_cell)%Centroid(2) - 0.5_8 * mesh%dy

        allocate(mesh%cells(ID_cell)%SubCellCentroids(mesh%n_sub * mesh%n_sub, 2))
        do sCell_y = 1, mesh%n_sub
          do sCell_x = 1, mesh%n_sub
            sub_id = (sCell_y - 1) * mesh%n_sub + sCell_x
            mesh%cells(ID_cell)%SubCellCentroids(sub_id, 1) = &
              x_start + (real(sCell_x, 8) - 0.5_8) * mesh%dx_sub
            mesh%cells(ID_cell)%SubCellCentroids(sub_id, 2) = &
              y_start + (real(sCell_y, 8) - 0.5_8) * mesh%dy_sub
          end do
        end do
      end do
    end do

    ! 2. Check Boundary Condition Type
    isPeriodic = (trim(prob%typeBC) == 'periodic' .or. trim(prob%typeBC) == 'PERIODIC')

    if (isPeriodic) then
      numVFaces = mesh%Ny * mesh%Nx
      numHFaces = mesh%Nx * mesh%Ny
    else
      numVFaces = mesh%Ny * (mesh%Nx + 1)
      numHFaces = mesh%Nx * (mesh%Ny + 1)
    end if

    mesh%numFaces = numVFaces + numHFaces
    allocate(mesh%faces(mesh%numFaces))
    fID = 1

    ! 3. Create Vertical Faces
    do j = 0, mesh%Ny - 1
      if (isPeriodic) then
        do i = 1, mesh%Nx
          len = mesh%dy
          if (i == mesh%Nx) then
            cLeftID  = mesh%Nx + j * mesh%Nx
            cRightID = 1 + j * mesh%Nx
            bc       = prob%typeBC
          else
            cLeftID  = i + j * mesh%Nx
            cRightID = i + 1 + j * mesh%Nx
            bc       = 'INTERIOR'
          end if
          mesh%faces(fID)%ID          = fID
          mesh%faces(fID)%Normal      = [1.0_8, 0.0_8]
          mesh%faces(fID)%Length      = len
          mesh%faces(fID)%LeftCellID  = cLeftID
          mesh%faces(fID)%RightCellID = cRightID
          mesh%faces(fID)%BCType      = bc
          fID = fID + 1
        end do
      else
        do i = 0, mesh%Nx
          len = mesh%dy
          if (i == 0) then
            cLeftID  = 1 + j * mesh%Nx
            cRightID = 0
            normal   = [-1.0_8, 0.0_8]
            bc       = prob%typeBC
          else if (i == mesh%Nx) then
            cLeftID  = mesh%Nx + j * mesh%Nx
            cRightID = 0
            normal   = [1.0_8, 0.0_8]
            bc       = prob%typeBC
          else
            cLeftID  = i + j * mesh%Nx
            cRightID = i + 1 + j * mesh%Nx
            normal   = [1.0_8, 0.0_8]
            bc       = 'INTERIOR'
          end if
          mesh%faces(fID)%ID          = fID
          mesh%faces(fID)%Normal      = normal
          mesh%faces(fID)%Length      = len
          mesh%faces(fID)%LeftCellID  = cLeftID
          mesh%faces(fID)%RightCellID = cRightID
          mesh%faces(fID)%BCType      = bc
          fID = fID + 1
        end do
      end if
    end do

    ! 4. Create Horizontal Faces
    if (isPeriodic) then
      do j = 1, mesh%Ny
        do i = 0, mesh%Nx - 1
          len = mesh%dx
          if (j == mesh%Ny) then
            cLeftID  = (i + 1) + (mesh%Ny - 1) * mesh%Nx
            cRightID = (i + 1)
            bc       = prob%typeBC
          else
            cLeftID  = (i + 1) + (j - 1) * mesh%Nx
            cRightID = (i + 1) + j * mesh%Nx
            bc       = 'INTERIOR'
          end if
          mesh%faces(fID)%ID          = fID
          mesh%faces(fID)%Normal      = [0.0_8, 1.0_8]
          mesh%faces(fID)%Length      = len
          mesh%faces(fID)%LeftCellID  = cLeftID
          mesh%faces(fID)%RightCellID = cRightID
          mesh%faces(fID)%BCType      = bc
          fID = fID + 1
        end do
      end do
    else
      do j = 0, mesh%Ny
        do i = 0, mesh%Nx - 1
          len = mesh%dx
          if (j == 0) then
            cLeftID  = (i + 1)
            cRightID = 0
            normal   = [0.0_8, -1.0_8]
            bc       = prob%typeBC
          else if (j == mesh%Ny) then
            cLeftID  = (i + 1) + (mesh%Ny - 1) * mesh%Nx
            cRightID = 0
            normal   = [0.0_8, 1.0_8]
            bc       = prob%typeBC
          else
            cLeftID  = (i + 1) + (j - 1) * mesh%Nx
            cRightID = (i + 1) + j * mesh%Nx
            normal   = [0.0_8, 1.0_8]
            bc       = 'INTERIOR'
          end if
          mesh%faces(fID)%ID          = fID
          mesh%faces(fID)%Normal      = normal
          mesh%faces(fID)%Length      = len
          mesh%faces(fID)%LeftCellID  = cLeftID
          mesh%faces(fID)%RightCellID = cRightID
          mesh%faces(fID)%BCType      = bc
          fID = fID + 1
        end do
      end do
    end if
  end subroutine generateMesh

  ! -------------------------------------------------------------------
  ! Initialize Reference Element Shape Functions & Mass Matrices
  ! Matches Mesh.initReferenceMath (Mesh.m, Lines 207-289)
  ! -------------------------------------------------------------------
  subroutine initReferenceMath(mesh)
    type(Mesh2D_t), intent(inout) :: mesh

    integer :: n1D_sub, num_quad, i, j, q, power
    real(8), allocatable :: A(:,:), S_2D(:,:), H_weighted(:,:)
    real(8) :: a_int, b_int, x_val, val, d_val, val_L, val_R, coef
    type(GaussQuadrature1D_t) :: GQ

    n1D_sub  = mesh%n_sub
    num_quad = mesh%n_quad

    allocate(A(n1D_sub, n1D_sub))
    allocate(mesh%alpha(n1D_sub, n1D_sub))
    allocate(mesh%S_1D(n1D_sub, n1D_sub))
    allocate(S_2D(n1D_sub * n1D_sub, n1D_sub * n1D_sub))
    allocate(mesh%S_inv_2D(n1D_sub * n1D_sub, n1D_sub * n1D_sub))
    allocate(mesh%H_quad(n1D_sub, num_quad))
    allocate(mesh%dH_quad(n1D_sub, num_quad))
    allocate(H_weighted(n1D_sub, num_quad))
    allocate(mesh%H_bound_L(n1D_sub))
    allocate(mesh%H_bound_R(n1D_sub))
    allocate(mesh%xq(num_quad))
    allocate(mesh%wq(num_quad))

    ! 1. Derive 1D Shape Functions on Reference Element [0, m+1] (Eqs. 16-21)
    do j = 1, n1D_sub
      do i = 1, n1D_sub
        A(i, j) = (real(i, 8)**j - real(i - 1, 8)**j) / real(j, 8)
      end do
    end do
    call invert_matrix(n1D_sub, A, mesh%alpha)

    ! 2. Setup Gauss-Legendre Quadrature on [0, m+1]
    call init_GaussQuadrature1D(GQ, num_quad)
    a_int   = 0.0_8
    b_int   = real(n1D_sub, 8)
    mesh%xq = 0.5_8 * (b_int - a_int) * GQ%xGauss + 0.5_8 * (b_int + a_int)
    mesh%wq = 0.5_8 * (b_int - a_int) * GQ%wGauss

    ! 3. Evaluate Shape Functions (H) and Derivatives (dH) at Quadrature Nodes
    do i = 1, n1D_sub
      do q = 1, num_quad
        x_val = mesh%xq(q)
        val   = 0.0_8
        d_val = 0.0_8
        do power = 0, mesh%m_degree
          coef = mesh%alpha(power + 1, i)
          if (power == 0) then
            val = val + coef
          else
            val   = val   + coef * (x_val**power)
            d_val = d_val + coef * real(power, 8) * (x_val**(power - 1))
          end if
        end do
        mesh%H_quad(i, q)  = val
        mesh%dH_quad(i, q) = d_val
        H_weighted(i, q)   = val * mesh%wq(q)
      end do
    end do

    ! 4. Construct 1D and 2D Mass/Stiffness Matrices (Eq. 22 & Eq. 26)
    mesh%S_1D = matmul(H_weighted, transpose(mesh%H_quad))
    call kron_square(n1D_sub, mesh%S_1D, mesh%S_1D, S_2D)
    call invert_matrix(n1D_sub * n1D_sub, S_2D, mesh%S_inv_2D)

    ! 5. Evaluate Shape Functions at Boundaries (xi = 0 and xi = m + 1)
    do i = 1, n1D_sub
      val_L = mesh%alpha(1, i)
      val_R = 0.0_8
      do power = 0, mesh%m_degree
        coef = mesh%alpha(power + 1, i)
        if (power == 0) then
          val_R = val_R + coef
        else
          val_R = val_R + coef * (real(n1D_sub, 8)**power)
        end if
      end do
      mesh%H_bound_L(i) = val_L
      mesh%H_bound_R(i) = val_R
    end do

    deallocate(A, S_2D, H_weighted, GQ%xGauss, GQ%wGauss)
  end subroutine initReferenceMath

  ! -------------------------------------------------------------------
  ! Compute sub-cell spatial grid coordinates (Ny*n_sub x Nx*n_sub)
  ! Matches Mesh.getSpatialGrid (Mesh.m, Lines 292-309)
  ! -------------------------------------------------------------------
  subroutine getSpatialGrid(mesh, X_plot, Y_plot)
    type(Mesh2D_t), intent(in)  :: mesh
    real(8),        intent(out) :: X_plot(mesh%Ny * mesh%n_sub, mesh%Nx * mesh%n_sub)
    real(8),        intent(out) :: Y_plot(mesh%Ny * mesh%n_sub, mesh%Nx * mesh%n_sub)

    integer :: i, j, sx, sy, cID
    real(8) :: x_start, y_start

    do j = 1, mesh%Ny
      do i = 1, mesh%Nx
        cID     = i + (j - 1) * mesh%Nx
        x_start = mesh%cells(cID)%Centroid(1) - 0.5_8 * mesh%dx
        y_start = mesh%cells(cID)%Centroid(2) - 0.5_8 * mesh%dy
        do sy = 1, mesh%n_sub
          do sx = 1, mesh%n_sub
            X_plot((j - 1) * mesh%n_sub + sy, (i - 1) * mesh%n_sub + sx) = &
              x_start + (real(sx, 8) - 0.5_8) * mesh%dx_sub
            Y_plot((j - 1) * mesh%n_sub + sy, (i - 1) * mesh%n_sub + sx) = &
              y_start + (real(sy, 8) - 0.5_8) * mesh%dy_sub
          end do
        end do
      end do
    end do
  end subroutine getSpatialGrid

  ! -------------------------------------------------------------------
  ! Deallocate Mesh2D arrays
  ! -------------------------------------------------------------------
  subroutine free_Mesh2D(mesh)
    type(Mesh2D_t), intent(inout) :: mesh
    integer :: c

    if (allocated(mesh%cells)) then
      do c = 1, size(mesh%cells)
        if (allocated(mesh%cells(c)%SubCellCentroids)) then
          deallocate(mesh%cells(c)%SubCellCentroids)
        end if
      end do
      deallocate(mesh%cells)
    end if
    if (allocated(mesh%faces))     deallocate(mesh%faces)
    if (allocated(mesh%alpha))     deallocate(mesh%alpha)
    if (allocated(mesh%S_1D))      deallocate(mesh%S_1D)
    if (allocated(mesh%S_inv_2D))  deallocate(mesh%S_inv_2D)
    if (allocated(mesh%H_quad))    deallocate(mesh%H_quad)
    if (allocated(mesh%dH_quad))   deallocate(mesh%dH_quad)
    if (allocated(mesh%H_bound_L)) deallocate(mesh%H_bound_L)
    if (allocated(mesh%H_bound_R)) deallocate(mesh%H_bound_R)
    if (allocated(mesh%xq))        deallocate(mesh%xq)
    if (allocated(mesh%wq))        deallocate(mesh%wq)
  end subroutine free_Mesh2D

end module Mesh2D
