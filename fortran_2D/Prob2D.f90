module Prob2D
  use PDE2D
  implicit none

  real(8), parameter :: PI_VAL = 4.0_8 * atan(1.0_8)

  integer, parameter :: PROB_ADV_SMOOTH          = 1
  integer, parameter :: PROB_ADV_SMOOTH_2        = 2
  integer, parameter :: PROB_ADV_WAVEPACKAGE     = 3
  integer, parameter :: PROB_ADV_SHOCK_DIAGONAL  = 4
  integer, parameter :: PROB_ADV_QUADRANT        = 5
  integer, parameter :: PROB_INVISCID_BURGERS_2D = 6

  type :: Prob2D_t
    character(64) :: typeProb  = ''
    character(32) :: shortName = ''
    integer       :: probID    = PROB_ADV_SMOOTH_2
    real(8)       :: xLeft     = 0.0_8
    real(8)       :: xRight    = 2.0_8
    real(8)       :: yLeft     = 0.0_8
    real(8)       :: yRight    = 2.0_8
    real(8)       :: finalTime = 2.0_8
    real(8)       :: maxSpeedX = 1.0_8
    real(8)       :: maxSpeedY = 1.0_8
    type(PDE2D_t) :: pde
    character(16) :: typeBC    = 'periodic'
    integer       :: withExact = 1
  end type Prob2D_t

contains

  ! -------------------------------------------------------------------
  ! Constructor: initialize 2D benchmark problem (eg_Prob_2D.m)
  ! -------------------------------------------------------------------
  subroutine init_Prob2D(prob, typeStr)
    type(Prob2D_t), intent(inout) :: prob
    character(*),   intent(in)    :: typeStr

    prob%typeProb = trim(typeStr)

    select case (trim(typeStr))

    case ('Scalar_Linear_Advection_Smooth')
      prob%probID    = PROB_ADV_SMOOTH
      prob%shortName = 'AdvSmooth'
      prob%xLeft     = 0.0_8
      prob%xRight    = 2.0_8
      prob%yLeft     = 0.0_8
      prob%yRight    = 2.0_8
      prob%finalTime = 2.0_8
      call init_PDE2D(prob%pde, 'Scalar_Linear_Advection_2D')
      prob%typeBC    = 'periodic'
      prob%maxSpeedX = abs(prob%pde%speedX)
      prob%maxSpeedY = abs(prob%pde%speedY)
      prob%withExact = 1

    case ('Scalar_Linear_Advection_Smooth_2')
      prob%probID    = PROB_ADV_SMOOTH_2
      prob%shortName = 'AdvSmooth2'
      prob%xLeft     = -1.0_8
      prob%xRight    =  1.0_8
      prob%yLeft     = -1.0_8
      prob%yRight    =  1.0_8
      prob%finalTime =  2.0_8
      call init_PDE2D(prob%pde, 'Scalar_Linear_Advection_2D')
      prob%typeBC    = 'periodic'
      prob%maxSpeedX = abs(prob%pde%speedX)
      prob%maxSpeedY = abs(prob%pde%speedY)
      prob%withExact = 1

    case ('Scalar_Linear_Advection_WavePackage')
      prob%probID    = PROB_ADV_WAVEPACKAGE
      prob%shortName = 'WavePackage'
      prob%xLeft     = -1.0_8
      prob%xRight    =  1.0_8
      prob%yLeft     = -1.0_8
      prob%yRight    =  1.0_8
      prob%finalTime =  6.3_8
      call init_PDE2D(prob%pde, 'Scalar_Linear_Advection_2D')
      prob%typeBC    = 'periodic'
      prob%maxSpeedX = abs(prob%pde%speedX)
      prob%maxSpeedY = abs(prob%pde%speedY)
      prob%withExact = 1

    case ('Scalar_Linear_Advection_Shock_Diagonal')
      prob%probID    = PROB_ADV_SHOCK_DIAGONAL
      prob%shortName = 'ShockDiag'
      prob%xLeft     = -1.0_8
      prob%xRight    =  1.0_8
      prob%yLeft     = -1.0_8
      prob%yRight    =  1.0_8
      prob%finalTime =  1.0_8
      call init_PDE2D(prob%pde, 'Scalar_Linear_Advection_2D')
      prob%typeBC    = 'fixed'
      prob%maxSpeedX = abs(prob%pde%speedX)
      prob%maxSpeedY = abs(prob%pde%speedY)
      prob%withExact = 1

    case ('Scalar_Linear_Advection_Quadrant')
      prob%probID    = PROB_ADV_QUADRANT
      prob%shortName = 'Quadrant'
      prob%xLeft     = -1.0_8
      prob%xRight    =  1.0_8
      prob%yLeft     = -1.0_8
      prob%yRight    =  1.0_8
      prob%finalTime =  1.0_8
      call init_PDE2D(prob%pde, 'Scalar_Linear_Advection_2D')
      prob%typeBC    = 'fixed'
      prob%maxSpeedX = abs(prob%pde%speedX)
      prob%maxSpeedY = abs(prob%pde%speedY)
      prob%withExact = 1

    case ('Scalar_Inviscid_Burgers', 'Scalar_Inviscid_Burgers_2D')
      prob%probID    = PROB_INVISCID_BURGERS_2D
      prob%shortName = 'Burgers2D'
      prob%xLeft     = 0.0_8
      prob%xRight    = 2.0_8
      prob%yLeft     = 0.0_8
      prob%yRight    = 2.0_8
      prob%finalTime = 1.0_8
      call init_PDE2D(prob%pde, 'Scalar_Inviscid_Burgers_2D')
      prob%typeBC    = 'periodic'
      prob%maxSpeedX = 1.0_8
      prob%maxSpeedY = 1.0_8
      prob%withExact = 1

    case default
      write(*,*) 'Error: Unknown 2D problem type: ', trim(typeStr)
      stop 1
    end select
  end subroutine init_Prob2D

  ! -------------------------------------------------------------------
  ! Mathematical modulo matching MATLAB mod(a, m) for positive m
  ! -------------------------------------------------------------------
  pure elemental real(8) function matlab_mod(a, m) result(res)
    real(8), intent(in) :: a, m
    res = a - floor(a / m) * m
  end function matlab_mod

  ! -------------------------------------------------------------------
  ! 1D Wave Package helper on [-1, 1]
  ! -------------------------------------------------------------------
  pure elemental real(8) function eval_wave_1D(x) result(q1)
    real(8), intent(in) :: x
    real(8), parameter  :: z     = -0.7_8
    real(8), parameter  :: delta = 0.005_8
    real(8), parameter  :: a     = 0.5_8
    real(8), parameter  :: alpha = 10.0_8
    real(8) :: beta

    beta = log(2.0_8) / (36.0_8 * delta * delta)
    q1   = 0.0_8

    if (x >= -0.8_8 .and. x <= -0.6_8) then
      q1 = ( exp(-beta * (x - (z - delta))**2) &
           + 4.0_8 * exp(-beta * (x - z)**2)   &
           + exp(-beta * (x - (z + delta))**2) ) / 6.0_8
    else if (x >= -0.4_8 .and. x <= -0.2_8) then
      q1 = 1.0_8
    else if (x >= 0.0_8 .and. x <= 0.2_8) then
      q1 = 1.0_8 - abs(10.0_8 * (x - 0.1_8))
    else if (x >= 0.4_8 .and. x <= 0.6_8) then
      q1 = ( sqrt(max(1.0_8 - alpha**2 * (x - (a - delta))**2, 0.0_8)) &
           + 4.0_8 * sqrt(max(1.0_8 - alpha**2 * (x - a)**2, 0.0_8))   &
           + sqrt(max(1.0_8 - alpha**2 * (x - (a + delta))**2, 0.0_8)) ) / 6.0_8
    end if
  end function eval_wave_1D

  ! -------------------------------------------------------------------
  ! 2D Wave Package Initial Condition & Exact Solution (WavePackage2D.m)
  ! -------------------------------------------------------------------
  pure elemental real(8) function WavePackage2D(x, y) result(q)
    real(8), intent(in) :: x, y
    real(8) :: xw, yw

    xw = matlab_mod(x + 1.0_8, 2.0_8) - 1.0_8
    yw = matlab_mod(y + 1.0_8, 2.0_8) - 1.0_8

    q = eval_wave_1D(xw) * eval_wave_1D(yw)
  end function WavePackage2D

  ! -------------------------------------------------------------------
  ! 2D 4-Quadrant Riemann Problem (Riemann2D_Quadrants.m)
  ! -------------------------------------------------------------------
  pure elemental real(8) function Riemann2D_Quadrants(x, y) result(q)
    real(8), intent(in) :: x, y
    real(8), parameter  :: xShock = 0.0_8, yShock = 0.0_8
    real(8), parameter  :: qLL = 0.5_8, qLR = 1.0_8, qRL = 1.5_8, qRR = 2.0_8

    if (x < xShock .and. y < yShock) then
      q = qLL
    else if (x < xShock .and. y >= yShock) then
      q = qLR
    else if (x >= xShock .and. y < yShock) then
      q = qRL
    else
      q = qRR
    end if
  end function Riemann2D_Quadrants

  ! -------------------------------------------------------------------
  ! 2D Inviscid Burgers Exact Solution via Newton-Raphson (exact_Burgers_2D.m)
  ! -------------------------------------------------------------------
  pure elemental real(8) function exact_Burgers_2D(x, y, t, xLeft, xRight, yLeft, yRight) result(U)
    real(8), intent(in) :: x, y, t, xLeft, xRight, yLeft, yRight
    real(8) :: Lx, Ly, xi, eta, F, dF, dF_safe, delta
    integer :: it
    integer, parameter :: maxIter = 40
    real(8), parameter :: tol = 1.0e-12_8

    Lx = xRight - xLeft
    Ly = yRight - yLeft

    U = sin(PI_VAL * x) * sin(PI_VAL * y)
    if (abs(t) < 1.0e-15_8) return

    do it = 1, maxIter
      xi  = matlab_mod(x - U * t - xLeft, Lx) + xLeft
      eta = matlab_mod(y - U * t - yLeft, Ly) + yLeft

      F  = U - sin(PI_VAL * xi) * sin(PI_VAL * eta)
      dF = 1.0_8 + PI_VAL * t * sin(PI_VAL * (xi + eta))

      if (abs(dF) < 1.0e-8_8) then
        if (dF >= 0.0_8) then
          dF_safe = 1.0e-8_8
        else
          dF_safe = -1.0e-8_8
        end if
      else
        dF_safe = dF
      end if

      delta = F / dF_safe
      delta = max(min(delta, 0.5_8), -0.5_8)
      U     = U - delta

      if (abs(delta) < tol) exit
    end do
  end function exact_Burgers_2D

  ! -------------------------------------------------------------------
  ! Evaluate Initial Condition u0(x, y)
  ! -------------------------------------------------------------------
  pure elemental real(8) function eval_IC(prob, x, y) result(u0)
    type(Prob2D_t), intent(in) :: prob
    real(8),        intent(in) :: x, y
    real(8), parameter         :: qLeft = 1.5_8, qRight = 0.5_8

    select case (prob%probID)
    case (PROB_ADV_SMOOTH)
      u0 = sin(PI_VAL * x) * sin(PI_VAL * y)

    case (PROB_ADV_SMOOTH_2)
      u0 = sin(PI_VAL * x - sin(PI_VAL * x) / PI_VAL) * &
           sin(PI_VAL * y - sin(PI_VAL * y) / PI_VAL)

    case (PROB_ADV_WAVEPACKAGE)
      u0 = WavePackage2D(x, y)

    case (PROB_ADV_SHOCK_DIAGONAL)
      if ((x + y) < 0.0_8) then
        u0 = qLeft
      else
        u0 = qRight
      end if

    case (PROB_ADV_QUADRANT)
      u0 = Riemann2D_Quadrants(x, y)

    case (PROB_INVISCID_BURGERS_2D)
      u0 = sin(PI_VAL * x) * sin(PI_VAL * y)

    case default
      u0 = 0.0_8
    end select
  end function eval_IC

  ! -------------------------------------------------------------------
  ! Evaluate Exact Analytical Solution u(x, y, t)
  ! -------------------------------------------------------------------
  pure elemental real(8) function eval_exact(prob, x, y, t) result(u_ex)
    type(Prob2D_t), intent(in) :: prob
    real(8),        intent(in) :: x, y, t
    real(8) :: xs, ys
    real(8), parameter :: qLeft = 1.5_8, qRight = 0.5_8

    xs = x - prob%pde%speedX * t
    ys = y - prob%pde%speedY * t

    select case (prob%probID)
    case (PROB_ADV_SMOOTH)
      u_ex = sin(PI_VAL * xs) * sin(PI_VAL * ys)

    case (PROB_ADV_SMOOTH_2)
      u_ex = sin(PI_VAL * xs - sin(PI_VAL * xs) / PI_VAL) * &
             sin(PI_VAL * ys - sin(PI_VAL * ys) / PI_VAL)

    case (PROB_ADV_WAVEPACKAGE)
      u_ex = WavePackage2D(xs, ys)

    case (PROB_ADV_SHOCK_DIAGONAL)
      if ((xs + ys) < 0.0_8) then
        u_ex = qLeft
      else
        u_ex = qRight
      end if

    case (PROB_ADV_QUADRANT)
      u_ex = Riemann2D_Quadrants(xs, ys)

    case (PROB_INVISCID_BURGERS_2D)
      u_ex = exact_Burgers_2D(x, y, t, prob%xLeft, prob%xRight, prob%yLeft, prob%yRight)

    case default
      u_ex = 0.0_8
    end select
  end function eval_exact

end module Prob2D
