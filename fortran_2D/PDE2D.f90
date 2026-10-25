module PDE2D
  implicit none

  integer, parameter :: PDE_LINEAR_ADVECTION_2D = 1
  integer, parameter :: PDE_INVISCID_BURGERS_2D = 2

  type :: PDE2D_t
    character(64) :: typeName = ''
    integer       :: pdeID    = PDE_LINEAR_ADVECTION_2D
    real(8)       :: speedX   = 1.0_8
    real(8)       :: speedY   = 1.0_8
  end type PDE2D_t

contains

  ! -------------------------------------------------------------------
  ! Constructor: initialize 2D PDE
  ! -------------------------------------------------------------------
  subroutine init_PDE2D(eq, typeStr)
    type(PDE2D_t), intent(inout) :: eq
    character(*),  intent(in)    :: typeStr

    eq%typeName = trim(typeStr)

    select case (trim(typeStr))
    case ('Scalar_Linear_Advection_2D')
      eq%pdeID  = PDE_LINEAR_ADVECTION_2D
      eq%speedX = 1.0_8
      eq%speedY = 1.0_8

    case ('Scalar_Inviscid_Burgers_2D')
      eq%pdeID  = PDE_INVISCID_BURGERS_2D
      eq%speedX = 1.0_8
      eq%speedY = 1.0_8

    case default
      eq%pdeID  = PDE_LINEAR_ADVECTION_2D
      eq%speedX = 1.0_8
      eq%speedY = 1.0_8
    end select
  end subroutine init_PDE2D

  ! -------------------------------------------------------------------
  ! Dynamically update linear advection wave speeds
  ! -------------------------------------------------------------------
  subroutine setSpeed_PDE2D(eq, cx, cy)
    type(PDE2D_t), intent(inout) :: eq
    real(8),       intent(in)    :: cx, cy
    eq%speedX = cx
    eq%speedY = cy
  end subroutine setSpeed_PDE2D

  ! -------------------------------------------------------------------
  ! Elemental Analytical Flux F(u) in X direction
  ! Works on scalars, 1D arrays, and 2D matrices seamlessly
  ! -------------------------------------------------------------------
  pure elemental real(8) function analyticalFluxX(eq, u) result(F)
    type(PDE2D_t), intent(in) :: eq
    real(8),       intent(in) :: u

    select case (eq%pdeID)
    case (PDE_LINEAR_ADVECTION_2D)
      F = eq%speedX * u
    case (PDE_INVISCID_BURGERS_2D)
      F = 0.5_8 * (u * u)
    case default
      F = eq%speedX * u
    end select
  end function analyticalFluxX

  ! -------------------------------------------------------------------
  ! Elemental Analytical Flux G(u) in Y direction
  ! Works on scalars, 1D arrays, and 2D matrices seamlessly
  ! -------------------------------------------------------------------
  pure elemental real(8) function analyticalFluxY(eq, u) result(G)
    type(PDE2D_t), intent(in) :: eq
    real(8),       intent(in) :: u

    select case (eq%pdeID)
    case (PDE_LINEAR_ADVECTION_2D)
      G = eq%speedY * u
    case (PDE_INVISCID_BURGERS_2D)
      G = 0.5_8 * (u * u)
    case default
      G = eq%speedY * u
    end select
  end function analyticalFluxY

end module PDE2D
