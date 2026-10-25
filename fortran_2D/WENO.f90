module WENO
  implicit none

  ! WENO type identifiers
  integer, parameter :: WENO5Z    = 1
  integer, parameter :: WENO6THETA = 2
  integer, parameter :: WENO6Z    = 3
  integer, parameter :: WENO6CZ   = 4

  type :: WENO_t
    integer  :: typeWENO
    integer  :: sizeBigStencil
    integer  :: numberStencil
    real(8)  :: gamma(4)
    real(8)  :: beta(4)
    real(8)  :: alpha(4)
    real(8)  :: omega(4)
    real(8)  :: tau
    real(8)  :: theta
    real(8)  :: reconstructQ(4)
    real(8)  :: epsWENO
    real(8)  :: alphaR
    real(8)  :: Ctau
    real(8)  :: pWENO
    ! WENO6cz / WENO6z extras
    real(8)  :: gamma6
    real(8)  :: beta6
    real(8)  :: omega6
    real(8)  :: reconstructQ6
  end type WENO_t

contains

  ! -------------------------------------------------------
  ! Constructor
  ! -------------------------------------------------------
  subroutine init_WENO(we, typeWE)
    type(WENO_t),   intent(inout) :: we
    character(*), intent(in)    :: typeWE

    we%epsWENO = 1.0e-10_8
    we%alphaR  = 20.0_8
    we%Ctau    = 5.0_8
    we%pWENO   = 1.0_8

    select case (trim(typeWE))
    case ('WENO5z')
      we%typeWENO       = WENO5Z
      we%numberStencil  = 3
      we%sizeBigStencil = 5
    case ('WENO6theta')
      we%typeWENO       = WENO6THETA
      we%numberStencil  = 4
      we%sizeBigStencil = 6
    case ('WENO6z')
      we%typeWENO       = WENO6Z
      we%numberStencil  = 4
      we%sizeBigStencil = 6
    case ('WENO6cz')
      we%typeWENO       = WENO6CZ
      we%numberStencil  = 4
      we%sizeBigStencil = 6
    case default
      write(*,*) 'Error: unknown WENO type: ', trim(typeWE)
    end select

    we%gamma       = 0.0_8
    we%beta        = 0.0_8
    we%omega       = 0.0_8
    we%reconstructQ = 0.0_8

  end subroutine init_WENO

  ! -------------------------------------------------------
  ! Set optimal weights
  ! -------------------------------------------------------
  subroutine setOptimalWeights(we)
    type(WENO_t), intent(inout) :: we

    select case (we%typeWENO)
    case (WENO5Z)
      we%gamma(1) = 1.0_8/10.0_8
      we%gamma(2) = 6.0_8/10.0_8
      we%gamma(3) = 3.0_8/10.0_8
    case (WENO6THETA)
      we%gamma(1) = 1.0_8/20.0_8
      we%gamma(2) = 3.0_8/20.0_8
      we%gamma(3) = 3.0_8/20.0_8
      we%gamma(4) = 1.0_8/20.0_8
    case (WENO6Z)
      we%gamma(1) = 1.0_8/20.0_8
      we%gamma(2) = 9.0_8/20.0_8
      we%gamma(3) = 3.0_8/10.0_8
      we%gamma(4) = 1.0_8/5.0_8
    case (WENO6CZ)
      we%gamma(1) = 0.002_8
      we%gamma(2) = 0.088_8
      we%gamma(3) = 0.008_8
      we%gamma(4) = 0.002_8
      we%gamma6   = 0.9_8
    end select

  end subroutine setOptimalWeights


  ! -------------------------------------------------------
  ! Compute beta 1
  ! -------------------------------------------------------
  subroutine computeBeta1(we, s)
    type(WENO_t), intent(inout) :: we
    real(8),    intent(in)    :: s(3)
    real(8) :: t1, t2

    select case (we%typeWENO)
    case (WENO5Z, WENO6Z, WENO6CZ)
      t1 = s(1) - 2.0_8*s(2) + s(3)
      t2 = s(1) - 4.0_8*s(2) + 3.0_8*s(3)
      we%beta(1) = 13.0_8/12.0_8*t1**2 + t2**2/4.0_8
    case (WENO6THETA)
      t1 = s(1) - 2.0_8*s(2) + s(3)
      t2 = s(1) - 3.0_8*s(2) + 2.0_8*s(3)
      we%beta(1) = 13.0_8/12.0_8*t1**2 + t2**2
    end select
  end subroutine computeBeta1

  ! -------------------------------------------------------
  ! Compute beta 2
  ! -------------------------------------------------------
  subroutine computeBeta2(we, s)
    type(WENO_t), intent(inout) :: we
    real(8),    intent(in)    :: s(3)
    real(8) :: t1, t2

    select case (we%typeWENO)
    case (WENO5Z, WENO6Z, WENO6CZ)
      t1 = s(1) - 2.0_8*s(2) + s(3)
      t2 = s(1) - s(3)
      we%beta(2) = 13.0_8/12.0_8*t1**2 + t2**2/4.0_8
    case (WENO6THETA)
      t1 = s(1) - 2.0_8*s(2) + s(3)
      t2 = s(3) - s(2)
      we%beta(2) = 13.0_8/12.0_8*t1**2 + t2**2
    end select
  end subroutine computeBeta2

  ! -------------------------------------------------------
  ! Compute beta 3
  ! -------------------------------------------------------
  subroutine computeBeta3(we, s)
    type(WENO_t), intent(inout) :: we
    real(8),    intent(in)    :: s(3)
    real(8) :: t1, t2

    select case (we%typeWENO)
    case (WENO5Z, WENO6Z, WENO6CZ)
      t1 = s(1) - 2.0_8*s(2) + s(3)
      t2 = 3.0_8*s(1) - 4.0_8*s(2) + s(3)
      we%beta(3) = 13.0_8/12.0_8*t1**2 + t2**2/4.0_8
    case (WENO6THETA)
      t1 = s(1) - 2.0_8*s(2) + s(3)
      !t2 = s(2) - s(3)
	  t2 = s(1) - s(2)
      we%beta(3) = 13.0_8/12.0_8*t1**2 + t2**2
    end select
  end subroutine computeBeta3

  ! -------------------------------------------------------
  ! Compute beta 4
  ! -------------------------------------------------------
  subroutine computeBeta4(we, s)
    type(WENO_t), intent(inout) :: we
    real(8),    intent(in)    :: s(:)
    real(8) :: t1, t2, t3

    select case (we%typeWENO)
    case (WENO6THETA)
      t1 = 3.0_8*s(1) - 7.0_8*s(2) + 5.0_8*s(3) - s(4)
      t2 = 2.0_8*s(2) - 3.0_8*s(3) + s(4)
      we%beta(4) = 13.0_8/48.0_8*t1**2 + t2**2
    case (WENO6Z, WENO6CZ)
      t1 = 11.0_8*s(1) - 18.0_8*s(2) + 9.0_8*s(3) - 2.0_8*s(4)
      t2 =  2.0_8*s(1) -  5.0_8*s(2) + 4.0_8*s(3) -       s(4)
      t3 =        s(1) -  3.0_8*s(2) + 3.0_8*s(3) -       s(4)
      we%beta(4) = 1.0_8/36.0_8*t1**2 + 13.0_8/12.0_8*t2**2 + 781.0_8/720.0_8*t3**2
    end select
  end subroutine computeBeta4


  ! -------------------------------------------------------
  ! Compute smoothness indicators
  ! -------------------------------------------------------
  subroutine computeSmoothness(we, bs)
    type(WENO_t), intent(inout) :: we
    real(8),    intent(in)    :: bs(:)
    real(8) :: tauCZ, tmp

    select case (we%typeWENO)
    case (WENO5Z)
      call computeBeta1(we, bs(1:3))
      call computeBeta2(we, bs(2:4))
      call computeBeta3(we, bs(3:5))

    case (WENO6THETA, WENO6Z)
      call computeBeta1(we, bs(1:3))
      call computeBeta2(we, bs(2:4))
      call computeBeta3(we, bs(3:5))
      call computeBeta4(we, bs(3:6))

    case (WENO6CZ)
      call computeBeta1(we, bs(1:3))
      call computeBeta2(we, bs(2:4))
      call computeBeta3(we, bs(3:5))
      call computeBeta4(we, bs(3:6))
      ! beta6 for WENO-CZ6
      tmp   = we%beta(1) + 4.0_8*we%beta(2) + we%beta(3)
      tauCZ = abs(we%beta(4) - tmp/6.0_8)
      we%beta6 = (1.0_8 + 2.0_8*tauCZ/(1.0_8 + tauCZ)) * (we%beta(1) + we%beta(4)) / 2.0_8
    end select

  end subroutine computeSmoothness

  ! -------------------------------------------------------
  ! Compute nonlinear weights
  ! -------------------------------------------------------
  subroutine computeNonlinearWeights(we, bs, dx)
    type(WENO_t), intent(inout) :: we
    real(8),    intent(in)    :: bs(:)
    real(8),    intent(in)    :: dx

    integer :: i
    real(8) :: sumAlpha, tauZ, tau6Z, tauCZ, tmp
    real(8) :: tau5, tau6, term1, term2
    real(8) :: maxBeta, minBeta, rBeta, alpha6

    call setOptimalWeights(we)

    select case (we%typeWENO)

    case (WENO5Z)
      call computeSmoothness(we, bs)
      tauZ = abs(we%beta(1) - we%beta(3))
      sumAlpha = 0.0_8
      do i = 1, we%numberStencil
        tmp = (tauZ / (we%epsWENO + we%beta(i)))**we%pWENO
        we%alpha(i) = we%gamma(i) * (1.0_8 + tmp)
        sumAlpha = sumAlpha + we%alpha(i)
      end do
      do i = 1, we%numberStencil
        we%omega(i) = we%alpha(i) / sumAlpha
      end do

    case (WENO6THETA)
      ! tau5
      term1 = bs(1) - 4.0_8*bs(2) + 6.0_8*bs(3) - 4.0_8*bs(4) + bs(5)
      term2 = -bs(2) + 3.0_8*bs(3) - 3.0_8*bs(4) + bs(5)
      tau5  = 13.0_8/12.0_8*term1**2 + term2**2 + we%Ctau*dx**6
      ! tau6
      term1 = -bs(1) + 5.0_8*bs(2) - 10.0_8*bs(3) + 10.0_8*bs(4) - 5.0_8*bs(5) + bs(6)
      term2 =  bs(1) - 3.0_8*bs(2) +  2.0_8*bs(3) +  2.0_8*bs(4) - 3.0_8*bs(5) + bs(6)
      tau6  = 13.0_8/12.0_8*term1**2 + term2**2/4.0_8 + we%Ctau*dx**8
      ! theta
      if (tau6 < tau5) then
        we%tau   = tau6
        we%theta = 0.0_8
      else
        we%tau   = tau5
        we%theta = 1.0_8
      end if
      ! re-scale optimal weights
      we%gamma(1) = we%gamma(1) * (1.0_8 + we%theta)
      we%gamma(2) = we%gamma(2) * (3.0_8 + we%theta)
      we%gamma(3) = we%gamma(3) * (3.0_8 - we%theta)
      we%gamma(4) = we%gamma(4) * (1.0_8 - we%theta)
      ! smoothness
      call computeSmoothness(we, bs)
      maxBeta = maxval(we%beta(1:we%numberStencil))
      minBeta = minval(we%beta(1:we%numberStencil))
      rBeta   = maxBeta / (we%epsWENO + minBeta)
      if (rBeta <= we%alphaR) we%beta(1:we%numberStencil) = 0.0_8
      ! nonlinear weights
      sumAlpha = 0.0_8
      do i = 1, we%numberStencil
        tmp = (we%tau / (we%epsWENO + we%beta(i)))**we%pWENO
        we%alpha(i) = we%gamma(i) * (1.0_8 + tmp)
        sumAlpha = sumAlpha + we%alpha(i)
      end do
      do i = 1, we%numberStencil
        we%omega(i) = we%alpha(i) / sumAlpha
      end do

    case (WENO6Z)
      call computeSmoothness(we, bs)
      tmp   = we%beta(1) + 4.0_8*we%beta(2) + we%beta(3)
      tau6Z = abs(we%beta(4) - tmp/6.0_8)
      sumAlpha = 0.0_8
      do i = 1, we%numberStencil
        tmp = tau6Z / (we%epsWENO + we%beta(i))
        we%alpha(i) = we%gamma(i) * (1.0_8 + tmp)
        sumAlpha = sumAlpha + we%alpha(i)
      end do
      do i = 1, we%numberStencil
        we%omega(i) = we%alpha(i) / sumAlpha
      end do

    case (WENO6CZ)
      call computeSmoothness(we, bs)
      tmp   = we%beta(1) + 4.0_8*we%beta(2) + we%beta(3)
      tauCZ = abs(we%beta(4) - tmp/6.0_8)
      sumAlpha = 0.0_8
      do i = 1, we%numberStencil
        tmp = tauCZ / (we%epsWENO + we%beta(i))
        we%alpha(i) = we%gamma(i) * (1.0_8 + tmp)
        sumAlpha = sumAlpha + we%alpha(i)
      end do
      tmp    = tauCZ / (we%epsWENO + we%beta6)
      alpha6 = we%gamma6 * (1.0_8 + tmp)
      sumAlpha = sumAlpha + alpha6
      do i = 1, we%numberStencil
        we%omega(i) = we%alpha(i) / sumAlpha
      end do
      we%omega6 = alpha6 / sumAlpha

    end select

  end subroutine computeNonlinearWeights


  ! -------------------------------------------------------
  ! Sub-stencil reconstructions
  ! -------------------------------------------------------
  subroutine reconstructSub1(we, s)
    type(WENO_t), intent(inout) :: we
    real(8),    intent(in)    :: s(3)
    we%reconstructQ(1) = 2.0_8/6.0_8*s(1) - 7.0_8/6.0_8*s(2) + 11.0_8/6.0_8*s(3)
  end subroutine reconstructSub1

  subroutine reconstructSub2(we, s)
    type(WENO_t), intent(inout) :: we
    real(8),    intent(in)    :: s(3)
    we%reconstructQ(2) = -1.0_8/6.0_8*s(1) + 5.0_8/6.0_8*s(2) + 2.0_8/6.0_8*s(3)
  end subroutine reconstructSub2

  subroutine reconstructSub3(we, s)
    type(WENO_t), intent(inout) :: we
    real(8),    intent(in)    :: s(3)
    we%reconstructQ(3) = 2.0_8/6.0_8*s(1) + 5.0_8/6.0_8*s(2) - 1.0_8/6.0_8*s(3)
  end subroutine reconstructSub3

  subroutine reconstructSub4(we, s)
    type(WENO_t), intent(inout) :: we
    real(8),    intent(in)    :: s(:)

    select case (we%typeWENO)
    case (WENO6THETA)
      we%reconstructQ(4) = 11.0_8/6.0_8*s(2) - 7.0_8/6.0_8*s(3) + 2.0_8/6.0_8*s(4)
    case (WENO6Z, WENO6CZ)
      we%reconstructQ(4) = 1.0_8/4.0_8*s(1) + 13.0_8/12.0_8*s(2) &
                         - 5.0_8/12.0_8*s(3) + 1.0_8/12.0_8*s(4)
    end select
  end subroutine reconstructSub4

  subroutine reconstruct6(we, bs)
    type(WENO_t), intent(inout) :: we
    real(8),    intent(in)    :: bs(6)

    if (we%typeWENO == WENO6CZ) then
      we%reconstructQ6 = (bs(1) - 8.0_8*bs(2) + 37.0_8*bs(3) &
                        + 37.0_8*bs(4) - 8.0_8*bs(5) + bs(6)) / 60.0_8
    end if
  end subroutine reconstruct6

  ! -------------------------------------------------------
  ! Reconstruct at cell interface
  ! -------------------------------------------------------
  real(8) function reconstructAt(we, bs, dx)
    type(WENO_t), intent(inout) :: we
    real(8),    intent(in)    :: bs(:)
    real(8),    intent(in)    :: dx

    integer :: i
    real(8) :: tmp1, tmp2

    reconstructAt = 0.0_8

    select case (we%typeWENO)

    case (WENO5Z)
      call reconstructSub1(we, bs(1:3))
      call reconstructSub2(we, bs(2:4))
      call reconstructSub3(we, bs(3:5))
      call computeNonlinearWeights(we, bs, dx)
      do i = 1, we%numberStencil
        reconstructAt = reconstructAt + we%omega(i) * we%reconstructQ(i)
      end do

    case (WENO6THETA)
      call reconstructSub1(we, bs(1:3))
      call reconstructSub2(we, bs(2:4))
      call reconstructSub3(we, bs(3:5))
      call reconstructSub4(we, bs(3:6))
      call computeNonlinearWeights(we, bs, dx)
      do i = 1, we%numberStencil
        reconstructAt = reconstructAt + we%omega(i) * we%reconstructQ(i)
      end do

    case (WENO6Z)
      call reconstructSub1(we, bs(1:3))
      call reconstructSub2(we, bs(2:4))
      call reconstructSub3(we, bs(3:5))
      call reconstructSub4(we, bs(3:6))
      call computeNonlinearWeights(we, bs, dx)
      do i = 1, we%numberStencil
        reconstructAt = reconstructAt + we%omega(i) * we%reconstructQ(i)
      end do

    case (WENO6CZ)
      call reconstructSub1(we, bs(1:3))
      call reconstructSub2(we, bs(2:4))
      call reconstructSub3(we, bs(3:5))
      call reconstructSub4(we, bs(3:6))
      call reconstruct6(we, bs(1:6))
      call computeNonlinearWeights(we, bs, dx)
      tmp1 = 0.0_8
      tmp2 = 0.0_8
      do i = 1, we%numberStencil
        tmp1 = tmp1 + we%gamma(i) * we%reconstructQ(i)
        tmp2 = tmp2 + we%omega(i) * we%reconstructQ(i)
      end do
      reconstructAt = we%omega6/we%gamma6 * (we%reconstructQ6 - tmp1) + tmp2

    end select

  end function reconstructAt

end module WENO
