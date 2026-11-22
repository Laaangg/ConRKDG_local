! ===================================================================
! GaussQuadrature2D Module: 2D Tensor-Product Gauss-Legendre Quadrature
! Pre-computes reference nodes and normalized weights on [-1, 1] x [-1, 1]
! ===================================================================
module GaussQuadrature2D
  use GaussQuadrature1D
  implicit none

  type :: GaussQuadrature2D_t
    integer              :: n1D    = 5     ! Number of 1D Gauss points per axis
    integer              :: nTotal = 25    ! Total 2D Gauss points = n1D^2
    real(8), allocatable :: xi(:)          ! Reference xi in [-1, 1], size (nTotal)
    real(8), allocatable :: eta(:)         ! Reference eta in [-1, 1], size (nTotal)
    real(8), allocatable :: w2D(:)         ! Normalized 2D weights, sum(w2D) = 1.0
  end type GaussQuadrature2D_t

contains

  ! -------------------------------------------------------------------
  ! Initialize 2D Tensor-Product Gauss-Legendre Quadrature
  ! -------------------------------------------------------------------
  subroutine init_GaussQuadrature2D(gq2D, n1D)
    type(GaussQuadrature2D_t), intent(inout) :: gq2D
    integer,                   intent(in)    :: n1D

    type(GaussQuadrature1D_t) :: gq1D
    integer :: i, j, idx

    gq2D%n1D    = n1D
    gq2D%nTotal = n1D * n1D

    if (allocated(gq2D%xi))  deallocate(gq2D%xi)
    if (allocated(gq2D%eta)) deallocate(gq2D%eta)
    if (allocated(gq2D%w2D)) deallocate(gq2D%w2D)

    allocate(gq2D%xi(gq2D%nTotal))
    allocate(gq2D%eta(gq2D%nTotal))
    allocate(gq2D%w2D(gq2D%nTotal))

    call init_GaussQuadrature1D(gq1D, n1D)

    ! Construct 2D tensor product.
    ! Normalization factor 0.25 ensures sum(w2D) = 1.0 on [-1, 1] x [-1, 1]
    idx = 0
    do j = 1, n1D
      do i = 1, n1D
        idx = idx + 1
        gq2D%xi(idx)  = gq1D%xGauss(i)
        gq2D%eta(idx) = gq1D%xGauss(j)
        gq2D%w2D(idx) = 0.25_8 * gq1D%wGauss(i) * gq1D%wGauss(j)
      end do
    end do

    deallocate(gq1D%xGauss, gq1D%wGauss)
  end subroutine init_GaussQuadrature2D

  ! -------------------------------------------------------------------
  ! Deallocate 2D Gauss Quadrature arrays
  ! -------------------------------------------------------------------
  subroutine free_GaussQuadrature2D(gq2D)
    type(GaussQuadrature2D_t), intent(inout) :: gq2D

    if (allocated(gq2D%xi))  deallocate(gq2D%xi)
    if (allocated(gq2D%eta)) deallocate(gq2D%eta)
    if (allocated(gq2D%w2D)) deallocate(gq2D%w2D)
    gq2D%n1D    = 0
    gq2D%nTotal = 0
  end subroutine free_GaussQuadrature2D

end module GaussQuadrature2D
