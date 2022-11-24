program cyl
  implicit none
  !Real units are in SI units
  !Underscore Variables are in LBM units
  logical, parameter:: &
    isDebug = .true.

  integer, parameter:: &
    chanH_ = 82, &
    dim = 2, &
    q = 9, &
    time_ = 100000, &
    noOfSnaps = 5, &
    dispFreq = 10, &
    noOfPtOnCircle = 400

  double precision, parameter:: &
    rhoF_ = 1.0d0, &
    rhoF = 1000.0d0, &
    chanL = 2.2d0, & !Length of channel
    chanH = 0.41d0, & !Width of channel
    barL = 0.35d0, &
    barH = 0.02d0, &
    uMean = 0.2d0, &
    nu = 0.001d0, &
    dia = 0.1d0, &
    xc = 0.2d0, &
    yc = 0.2d0

  double precision, parameter:: &
    d0 = 0.0d0, &
    d1 = 1.0d0, &
    d2 = 2.0d0, &
    haf = 0.5d0, &
    one6th = 1.0d0/6.0d0, &
    one3rd = 1.0d0/3.0d0, &
    pi = 4.0d0*atan(1.0d0), &
    one36th = 1.0d0/36.0d0

  type lbmTriplet_t
    double precision::r, u, v
  end type lbmTriplet_t
  type doublet_t
    double precision::x, y
  end type doublet_t

  type triplet_t
    double precision::x, y, f(0:q - 1), fneq(0:q - 1)
  end type triplet_t

  type, extends(triplet_t):: box_t
    ! double precision::x, y, z(0:q - 1)
    type(triplet_t)::fluidNode(2)
    logical:: isInside
  end type box_t

  type custom_t
    type(doublet_t)::unitVec, force
    type(triplet_t)::Pt
    type(box_t)::box(4)
  end type custom_t

  type cf_t
    double precision:: length, rho, time, velocity, nu, force
  end type

  type si_lb_t
    double precision:: channelHeight, channelLength, totalTime, diameter
    double precision:: xCentre, yCentre, barHeight, barLength, meanVel, nu, tau
  end type

  type(cf_t)::cf
  type(si_lb_t)::si, lb
  type(custom_t), allocatable, dimension(:)::ptOnCircle
  type(doublet_t)::totalForce
  type(lbmTriplet_t)::onSurf

  integer::nx, ny
  double precision::  uPara_, uParaRamp_
  double precision:: t, invTau, sigma(2, 2)
  integer:: i, j, k, a, a1, t_, ia, ja, solnumber
  integer, allocatable, dimension(:, :)::isn
  double precision:: tmp1, tmp2, tmp3, rhoSum, fx_t, fy_t, Cd, Cl, Cd2, Cl2
  double precision:: fx(2), fy(2), dudx, dudy, dvdx, dvdy, f_neq
  double precision, dimension(0:q - 1):: wi, feq, Q_xx, Q_yy, Q_xy
  integer, dimension(0:q - 1):: kb
  integer:: ci(0:q - 1, 2)
  double precision, dimension(0:q - 1):: tmpA
  double precision, allocatable, dimension(:, :, :):: f, ft
  type(lbmTriplet_t), allocatable, dimension(:, :):: lbmVar
  logical::isCyl, isBar
  character(len=30):: filename
  solnumber = 0

  invTau = 1.0d0/lb%tau

  write (*, *) '======================================================'
  write (*, *) 'Program started at :', dateTime()
!===Conversion Factors===
  call calcConvFactor(cf)
!===Other LBM parameters===
  call calcLBparams(si, lb)
  nx = int(lb%channelLength)
  ny = int(lb%channelHeight)

  write (*, *) 'lb.tau  = ', lb%tau
  write (*, *) 'lb.uMax = ', 1.5d0*lb%meanVel
  write (*, *) 'si.Re = ', si%meanVel*si%diameter/si%nu
  write (*, *) 'lb.Re = ', lb%meanVel*lb%diameter/lb%nu
  ! write (*, *) 'Aborted for Checking'; stop
!----------------------------------------------------------------------
  allocate (isn(nx + 2, ny + 2))
  allocate (lbmVar(nx + 2, ny + 2))
  ! allocate (uy(nx + 2, ny + 2))
  ! allocate (rho(nx + 2, ny + 2))
  allocate (f(0:q - 1, nx + 2, ny + 2))
  allocate (ft(0:q - 1, nx + 2, ny + 2))
!----------------------------------------------------------------------
  call setupCiWi(ci, wi, kb)
!----------------------------------------------------------------------
  call initDistFun(f)
!----------------------------------------------------------------------
  fx(1) = d0
  fy(1) = d0
  open (unit=10, file="tRhoCdCl.dat")
  write (10, *) "Variables=timeLBM,timeReal,rho,Cd,Cl,Cd2,Cl2"
  ! open (unit=11, file="VelSig.dat")
!----------------------------------------------------------------------
  do i = 1, nx + 1
    do j = 1, ny + 2

      ! ii = i - 1.5
      ! jj = j - 1.5

      isCyl = ((i - lb%xCentre)**2.0 + (j - lb%yCentre)**2.0)**0.5 .le. 0.5*lb%diameter

      isBar = (i .ge. lb%xCentre) .and. (i .le. lb%xCentre + lb%diameter/2 + lb%barLength) .and. &
              (j .ge. lb%xCentre - lb%barHeight/2) .and. (j .le. lb%xCentre + lb%barHeight/2)

      isBar = .false.

      if (isCyl .or. isBar) then
        isn(i, j) = 1
      elseif (j == 1 .or. j == ny + 2) then
        isn(i, j) = 2
      else
        isn(i, j) = 0
      end if

    end do
  end do

  call createDataStructOnCircle(noOfPtOnCircle, ptOnCircle)
!----------------------------------------------------------------------
  do t_ = 0, time_

    t = t_*cf%time

    rhoSum = d0
    do i = 2, nx + 1
      do j = 2, ny + 1
        lbmVar(i, j) = macroVar(f(:, i, j))
        rhoSum = rhoSum + lbmVar(i, j)%r
      end do
    end do

    rhoSum = rhoSum/(nx*ny)
    if (rhoSum .gt. 10.0d0) then
      write (*, *) 'Code Diverged. Aborted by User'
      stop
    end if
!----------------------------------------------------------------------
    do i = 2, nx + 1
      do j = 2, ny + 1
        feq = eqDist(lbmVar(i, j))
        ! do a = 0, q - 1
        ft(:, i, j) = f(:, i, j) - (f(:, i, j) - feq)/lb%tau !collision
        ! end do
      end do
    end do
!----------------------------------------------------------------------
    do i = 2, nx + 1 !Streaming post-collision
      do j = 2, ny + 1
        do a = 0, q - 1
          ia = i + ci(a, 1)
          ja = j + ci(a, 2)
          !if (ia<1 )        { ia = nx  }
          !if (ia>nx)        { ia = 1   }

          f(a, ia, ja) = ft(a, i, j)
        end do
      end do
    end do
!----------------------------------------------------------------------
    do j = 2, ny + 1
      !Inlet
      i = 2

      uPara_ = 6.0d0*lb%meanVel*(ny - (j - 1.5))*(j - 1.5)/ny**2.0d0; 
      if (t .lt. 2.0d0) then
        uParaRamp_ = uPara_*(1 - cos(pi*t/2.0d0))/2.0d0
      else
        uParaRamp_ = uPara_
      end if

      lbmVar(i, j)%u = uParaRamp_
      lbmVar(i, j)%v = 0.0
      lbmVar(i, j)%r = (f(0, i, j) + f(2, i, j) + f(4, i, j) + 2.0*(f(3, i, j) + f(6, i, j) + f(7, i, j)))/(1 - lbmVar(i, j)%u)
      tmp2 = lbmVar(i, j)%u*lbmVar(i, j)%u + lbmVar(i, j)%v*lbmVar(i, j)%v

      dudx = (-3.0*lbmVar(i, j)%u + 4.0*lbmVar(i + 1, j)%u - 1.0*lbmVar(i + 2, j)%u)/2.0 !dudx
      if (j == 2) then
        dudy = (-3.0*lbmVar(i, j)%u + 4.0*lbmVar(i, j + 1)%u - 1.0*lbmVar(i, j + 2)%u)/2.0
      else if (j == ny + 1) then
        dudy = (3.0*lbmVar(i, j)%u - 4.0*ux(i, j - 1)%u + 1.0*ux(i, j - 2)%u)/2.0
      else
        dudy = (ux(i, j + 1) - ux(i, j - 1))/2.0
      end if

      dvdx = (-3.0*lbmVar(i, j)%v + 4.0*uy(i + 1, j) - 1.0*uy(i + 2, j))/2.0
      if (j == 2) then
        dvdy = (-3.0*lbmVar(i, j)%v + 4.0*uy(i, j + 1) - 1.0*uy(i, j + 2))/2.0
      else if (j == ny + 1) then
        dvdy = (3.0*lbmVar(i, j)%v - 4.0*uy(i, j - 1) + 1.0*uy(i, j - 2))/2.0
      else
        dvdy = (uy(i, j + 1) - uy(i, j - 1))/2.0
      end if

      do a = 0, q - 1
        Q_xx(a) = ci(a, 1)*ci(a, 1) - 1.0/3.0
        Q_yy(a) = ci(a, 2)*ci(a, 2) - 1.0/3.0
        Q_xy(a) = ci(a, 1)*ci(a, 2)
        f_neq = -wi(a)*3.0*lb%tau*lbmVar(i, j)%r*(Q_xx(a)*dudx + Q_xy(a)*dudy + Q_xy(a)*dvdx + Q_yy(a)*dvdy)
        tmp1 = ci(a, 1)*lbmVar(i, j)%u + ci(a, 2)*lbmVar(i, j)%v
        feq(a) = wi(a)*lbmVar(i, j)%r*(1 + tmp1*3.0 + 0.5*tmp1*tmp1*3.0*3.0 - 0.5*3.0*tmp2)
        f(a, i, j) = feq(a) + f_neq
      end do

      !Outlet
      i = nx + 1
      lbmVar(i, j)%r = rhoF_
      lbmVar(i, j)%u = (f(0, i, j) + f(2, i, j) + f(4, i, j) + 2*(f(1, i, j) + f(5, i, j) + f(8, i, j)))/lbmVar(i, j)%r - 1
      lbmVar(i, j)%v = 0.0
      tmp2 = lbmVar(i, j)%u*lbmVar(i, j)%u + lbmVar(i, j)%v*lbmVar(i, j)%v

      dudx = (3.0*lbmVar(i, j)%u - 4.0*ux(i - 1, j) + 1.0*ux(i - 2, j))/2.0 !dudx
      if (j == 2) then
        dudy = (-3.0*lbmVar(i, j)%u + 4.0*ux(i, j + 1) - 1.0*ux(i, j + 2))/2.0
      else if (j == ny + 1) then
        dudy = (3.0*lbmVar(i, j)%u - 4.0*ux(i, j - 1) + 1.0*ux(i, j - 2))/2.0
      else
        dudy = (ux(i, j + 1) - ux(i, j - 1))/2.0
      end if

      dvdx = (3.0*lbmVar(i, j)%v - 4.0*uy(i - 1, j) + 1.0*uy(i - 2, j))/2.0 !dvdx
      if (j == 2) then
        dvdy = (-3.0*lbmVar(i, j)%v + 4.0*uy(i, j + 1) - 1.0*uy(i, j + 2))/2.0
      else if (j == ny + 1) then
        dvdy = (3.0*lbmVar(i, j)%v - 4.0*uy(i, j - 1) + 1.0*uy(i, j - 2))/2.0
      else
        dvdy = (uy(i, j + 1) - uy(i, j - 1))/2.0
      end if

      do a = 0, q - 1
        Q_xx(a) = ci(a, 1)*ci(a, 1) - 1.0/3.0
        Q_yy(a) = ci(a, 2)*ci(a, 2) - 1.0/3.0
        Q_xy(a) = ci(a, 1)*ci(a, 2)
        f_neq = -wi(a)*3.0*lb%tau*lbmVar(i, j)%r*(Q_xx(a)*dudx + Q_xy(a)*dudy + Q_xy(a)*dvdx + Q_yy(a)*dvdy)
        tmp1 = ci(a, 1)*lbmVar(i, j)%u + ci(a, 2)*lbmVar(i, j)%v
        feq(a) = wi(a)*lbmVar(i, j)%r*(1 + tmp1*3.0 + 0.5*tmp1*tmp1*3.0*3.0 - 0.5*3.0*tmp2)
        f(a, i, j) = feq(a) + f_neq
      end do
    end do
!----------------------------------------------------------------------
    ! do j = 2, ny + 1
    !    i = 2
    !    uPara_ = 6.0d0*lb%meanVel*(ny - (j - 1.5))*(j - 1.5)/ny**2.0d0;
    !    if (t .lt. 2.0d0) then
    !       uParaRamp_ = uPara_*(1 - cos(pi*t/2.0d0))/2.0d0
    !    else
    !       uParaRamp_ = uPara_
    !    end if
    !    lbmVar(i, j)%r = (f(0, i, j) + f(2, i, j) + f(4, i, j) + 2*(f(6, i, j) + f(3, i, j) + f(7, i, j)))/(1 - uParaRamp_)
    !    f(1, i, j) = f(3, i, j) + ((2.0/3.0)*lbmVar(i, j)%r*uParaRamp_)
    !    f(5, i, j) = f(7, i, j) - (0.5*(f(2, i, j) - f(4, i, j))) + ((1.0/6.0)*lbmVar(i, j)%r*uParaRamp_)
    !    f(8, i, j) = f(6, i, j) + (0.5*(f(2, i, j) - f(4, i, j))) + ((1.0/6.0)*lbmVar(i, j)%r*uParaRamp_)
    ! end do

    ! do j = 2, ny + 1
    !    i = nx + 1
    !    lbmVar(i, j)%u = (f(0, i, j) + f(2, i, j) + f(4, i, j) + 2*(f(1, i, j) + f(5, i, j) + f(8, i, j)))/rhoF_ - 1
    !    f(3, i, j) = f(1, i, j) - ((2.0/3.0)*rhoF_*lbmVar(i, j)%u)
    !    f(6, i, j) = f(8, i, j) - 0.5*(f(2, i, j) - f(4, i, j)) - ((1.0/6.0)*rhoF_*lbmVar(i, j)%u)
    !    f(7, i, j) = f(5, i, j) + 0.5*(f(2, i, j) - f(4, i, j)) - ((1.0/6.0)*rhoF_*lbmVar(i, j)%u)
    ! end do
!----------------------------------------------------------------------
!=================working================================
    do i = 1, size(ptOnCircle)
      do k = 1, 4

        associate (lf => ptOnCircle(i)%box(k)%fluidNode, lbox => ptOnCircle(i)%box(k))

          ! lf(2)%x = 45.0d0
          ! lf(2)%y = 44.0d0
          ! lf(2)%z = [20.0d0, 20.0d0, 20.0d0, 30.0d0, 20.0d0, 20.0d0, 20.0d0, 20.0d0, 200.0d0]

          ! lf(1)%x = 45.0d0
          ! lf(1)%y = 45.0d0
          ! lf(1)%z = [15.0d0, 20.0d0, 10.0d0, 10.0d0, 10.0d0, 10.0d0, 10.0d0, 10.0d0, 100.0d0]

          ! lbox%x = 45.0d0
          ! lbox%y = 46.0d0
          ! call linearExterp(lbox)
          ! write (*, *) "Hello", lbox%z
          ! stop

          if (lbox%isInside) then
            tmpA = f(:, int(lf(1)%x), int(lf(1)%y))
            lf(1)%f = tmpA
            lf(1)%fneq = tmpA - eqDist(macroVar(tmpA))

            tmpA = f(:, int(lf(2)%x), int(lf(2)%y))
            lf(2)%f = tmpA
            lf(2)%fneq = tmpA - eqDist(macroVar(tmpA))

            call linearExterp(lbox)
          else
            tmpA = f(:, int(lbox%x), int(lbox%y))
            lbox%f = tmpA
            lbox%fneq = tmpA - eqDist(macroVar(tmpA))
          end if

        end associate

      end do

      ! ptOnCircle(1)%box(1)%z = [10.0d0, 20.0d0, 10.0d0, 10.0d0, 10.0d0, 10.0d0, 10.0d0, 10.0d0, 100.0d0]
      ! ptOnCircle(1)%box(2)%z = [20.0d0, 20.0d0, 10.0d0, 10.0d0, 10.0d0, 10.0d0, 10.0d0, 10.0d0, 100.0d0]
      ! ptOnCircle(1)%box(3)%z = [30.0d0, 20.0d0, 10.0d0, 10.0d0, 10.0d0, 10.0d0, 10.0d0, 10.0d0, 100.0d0]
      ! ptOnCircle(1)%box(4)%z = [40.0d0, 20.0d0, 10.0d0, 10.0d0, 10.0d0, 10.0d0, 10.0d0, 10.0d0, 100.0d0]
      ! ! ptOnCircle(1)%Pt%x = 51
      ! ! ptOnCircle(1)%Pt%y = 42
      ! call bilinearInterp(ptOnCircle(1)%Pt, ptOnCircle(1)%box)
      ! write (*, *) ptOnCircle(1)%Pt%x
      ! write (*, *) ptOnCircle(1)%Pt%y
      ! write (*, *) ptOnCircle(1)%Pt%z
      ! stop

      call bilinearInterp(ptOnCircle(i)%Pt, ptOnCircle(i)%box)

      call calcStressTensor2NEQ(ptOnCircle(i)%Pt%f, ptOnCircle(i)%Pt%fneq, sigma)

      associate (uv => ptOnCircle(i)%unitVec)
        ptOnCircle(i)%force%x = sigma(1, 1)*uv%x + sigma(1, 2)*uv%y! - onSurf%r*onSurf%u*(onSurf%u*uv%x + onSurf%v*uv%y)
        ptOnCircle(i)%force%y = sigma(2, 1)*uv%x + sigma(2, 2)*uv%y! - onSurf%r*onSurf%v*(onSurf%u*uv%x + onSurf%v*uv%y)
      end associate
    end do

    ! ptOnCircle%force%x = 1.75
    ! write (*, *) lb%diameter
    ! write (*, *) trapIntegrate(ptOnCircle%force%x)
    ! stop

    totalForce%x = trapIntegrate(ptOnCircle%force%x)
    totalForce%y = trapIntegrate(ptOnCircle%force%y)

    ! write (*, *) '======================================'
    ! write (*, *) (ptOnCircle%Pt)
    ! stop

!=================working================================

    fx(2) = d0
    fy(2) = d0
    do i = 2, nx + 1 !BC
      do j = 2, ny + 1
        if (isn(i, j) .eq. 0) then

          do a = 0, q - 1

            ia = i + ci(a, 1)
            ja = j + ci(a, 2)

            if (isn(ia, ja) .eq. 1) then !structure
              f(kb(a), i, j) = ft(a, i, j)
              f(a, ia, ja) = ft(kb(a), ia, ja)
              fx(2) = fx(2) + ci(a, 1)*2.0*(-ft(kb(a), ia, ja) + ft(a, i, j))
              fy(2) = fy(2) + ci(a, 2)*2.0*(-ft(kb(a), ia, ja) + ft(a, i, j))
            end if

            if (isn(ia, ja) .eq. 2) then !wall
              f(kb(a), i, j) = ft(a, i, j)
            end if

          end do
        end if
      end do
    end do
!----------------------------------------------------------------------
    fx_t = 0.5*(fx(1) + fx(2))
    fy_t = 0.5*(fy(1) + fy(2))
    ! Cd = fx_t*cf%force!/(0.5*rhoF_*lb%meanVel*lb%meanVel*lb%diameter)
    ! Cl = fy_t*cf%force!/(0.5*rhoF_*lb%meanVel*lb%meanVel*lb%diameter)
    Cd = fx_t/(0.5*rhoF_*lb%meanVel*lb%meanVel*lb%diameter)
    Cd2 = totalForce%x/(0.5*rhoF_*lb%meanVel*lb%meanVel*lb%diameter)
    Cl = fy_t/(0.5*rhoF_*lb%meanVel*lb%meanVel*lb%diameter)
    Cl2 = totalForce%y/(0.5*rhoF_*lb%meanVel*lb%meanVel*lb%diameter)
    write (*, *) rhoF_
!----------------------------------------------------------------------
    if (mod(t_, dispFreq) .eq. 0) then
      write (10, '(E16.6,2X,I10,5(2X,E12.4))') t, t_, rhoSum, Cd, Cl, Cd2, Cl2
      ! write (*, '(E16.6,2X,I10,5(2X,E12.4))') t, t_, rhoSum, Cd, Cl, Cd2, Cl2
      !write (*, '(I8,4(3X,F10.6))') ts, ts*cf%time, rhoAvg, Cd, Cl

      !write (11, '(I8,6(3X,F10.6))') ts, ux(150, 201), uy(150, 201), ux(200, 250), uy(200, 250), ux(250, 201), uy(250, 201)
    end if
!----------------------------------------------------------------------
    fx(1) = fx(2)
    fy(1) = fy(2)
!----------------------------------------------------------------------
    if (t_ .le. time_ .and. mod(t_, (time_/(noOfSnaps - 1))) .eq. 0) then

      solnumber = solnumber + 1
      write (filename, '(a,i3.3,a)') "snap", solnumber, ".dat"
      open (unit=12, file=filename)

      write (12, *) "Variables=x,y,u,v,rho,region"
      write (12, '(2(a,I6))') "Zone I=", nx, ",J=", ny

      do j = 2, ny + 1
        do i = 2, nx + 1
          write (12, '(2(2X,I6),3(2X,E12.4),2X,I3)') i, j, lbmVar(i, j)%u, lbmVar(i, j)%v, lbmVar(i, j)%r, isn(i, j)
        end do
        write (12, *)
      end do
      close (12)
      write (*, '(a,I3,a,I10)') "snap", solnumber, " recorded at LBM time %d", t_
    end if
!----------------------------------------------------------------------
  end do!Time loop Ends

  close (10)
  ! close (11)
  write (*, *) 'Program ended at :', dateTime()
  write (*, *) '======================================================'
contains

  pure subroutine calcStressTensor(f, sigma, OnSurf)
    implicit none

    double precision, dimension(0:q - 1), intent(in) :: f
    double precision, intent(out) ::  sigma(2, 2)
    type(lbmTriplet_t), intent(out) :: onSurf
    double precision:: u(2), tmp(3), rho
    integer:: i, j, a

    tmp = d0
    do a = 0, q - 1
      tmp(1) = tmp(1) + f(a)
      tmp(2) = tmp(2) + f(a)*ci(a, 1)
      tmp(3) = tmp(3) + f(a)*ci(a, 2)
    end do

    rho = tmp(1)
    u(1) = tmp(2)/rho
    u(2) = tmp(3)/rho

    onSurf%r = rho
    onSurf%u = u(1)
    onSurf%v = u(2)

    do i = 1, 2
      do j = 1, 2

        tmp(1) = d0
        do a = 0, q - 1
          tmp(1) = tmp(1) + (ci(a, i) - u(i))*(ci(a, j) - u(j))*f(a)
        end do

        sigma(i, j) = -one6th*invTau*rho*kdf(i, j) - (d1 - haf*invTau)*tmp(1)

      end do
    end do

  end subroutine calcStressTensor

  pure subroutine calcStressTensor2(f, sigma, lbmVar)
    implicit none

    double precision, dimension(0:q - 1), intent(in) :: f
    double precision, intent(out) ::  sigma(2, 2)
    type(lbmTriplet_t), intent(out) :: lbmVar

    double precision, dimension(0:q - 1) :: feq, fneq
    double precision:: tmp(3)
    integer:: i, j, a

    lbmVar = macroVar(f)

    feq = eqDist(lbmVar)

    fneq = f - feq

    do i = 1, 2
      do j = 1, 2

        tmp(1) = d0
        do a = 0, q - 1
          tmp(1) = tmp(1) + fneq(a)*(ci(a, i)*ci(a, j) - haf*(ci(a, 1)*ci(a, 1) + ci(a, 2)*ci(a, 2))*kdf(i, j))
        end do

        sigma(i, j) = (d1 - haf*invTau)*tmp(1)

      end do
    end do

    sigma(1, 1) = sigma(1, 1) - one3rd*lbmVar%r
    sigma(2, 2) = sigma(2, 2) - one3rd*lbmVar%r

  end subroutine calcStressTensor2

  pure subroutine calcStressTensor2NEQ(f, fneq, sigma)
    implicit none

    double precision, dimension(0:q - 1), intent(in) :: f, fneq
    double precision, intent(out) ::  sigma(2, 2)
    type(lbmTriplet_t)::lbmVar
    double precision:: tmp(3)
    integer:: i, j, a

    lbmVar = macroVar(f)

    do i = 1, 2
      do j = 1, 2

        tmp(1) = d0
        do a = 0, q - 1
          tmp(1) = tmp(1) + fneq(a)*(ci(a, i)*ci(a, j) - haf*(ci(a, 1)*ci(a, 1) + ci(a, 2)*ci(a, 2))*kdf(i, j))
        end do

        sigma(i, j) = (d1 - haf*invTau)*tmp(1)

      end do
    end do

    sigma(1, 1) = sigma(1, 1) - one3rd*lbmVar%r
    sigma(2, 2) = sigma(2, 2) - one3rd*lbmVar%r

  end subroutine calcStressTensor2NEQ

  pure function eqDist(lbmVar) result(feq)
    ! double precision, dimension(0:q - 1), intent(in) :: f
    type(lbmTriplet_t), intent(in) :: lbmVar
    double precision, dimension(0:q - 1) :: feq
    integer:: a
    double precision:: tmp1, tmp2

    do a = 0, q - 1
      tmp1 = lbmVar%u*ci(a, 1) + lbmVar%v*ci(a, 2)
      tmp2 = lbmVar%u**d2 + lbmVar%v**d2
      feq(a) = wi(a)*lbmVar%r*(1.0 + 3.0*tmp1 + 4.5*tmp1*tmp1 - 1.5*tmp2)
      ! ft(a) = f(a) - (f(a) - feq)/lb%tau !collision
    end do

  end function eqDist

  pure function macroVar(f) result(lbmVar)
    implicit none

    double precision, dimension(0:q - 1), intent(in) :: f
    type(lbmTriplet_t) :: lbmVar

    double precision:: tmp(3)
    integer:: a

    tmp = d0
    do a = 0, q - 1
      tmp(1) = tmp(1) + f(a)
      tmp(2) = tmp(2) + f(a)*ci(a, 1)
      tmp(3) = tmp(3) + f(a)*ci(a, 2)
    end do

    lbmVar%r = tmp(1)
    lbmVar%u = tmp(2)/lbmVar%r
    lbmVar%v = tmp(3)/lbmVar%r

  end function macroVar

  pure function kdf(i, j) result(retval)
    integer, intent(in) :: i, j
    integer :: retval

    if (i == j) then
      retval = 1
    else
      retval = 0
    end if
  end function kdf

  pure subroutine createDataStructOnCircle(noOfPts, ptOnCircle)
    integer, intent(in) :: noOfPts
    type(custom_t), allocatable, dimension(:), intent(out)::ptOnCircle
    double precision::theta0, theta, dTheta, dirDotUnitVec(q - 1)
    integer::i, k, a, outDir, itmp(1)
    allocate (ptOnCircle(noOfPts))

    theta0 = d0
    dTheta = 2*pi/noOfPts

    do i = 0, noOfPts - 1
      theta = theta0 + i*dTheta
      ptOnCircle(i + 1)%Pt%x = lb%xCentre + 0.5*lb%diameter*cos(theta)
      ptOnCircle(i + 1)%Pt%y = lb%yCentre + 0.5*lb%diameter*sin(theta)

      ptOnCircle(i + 1)%unitVec%x = cos(theta)
      ptOnCircle(i + 1)%unitVec%y = sin(theta)

      associate (locBox => ptOnCircle(i + 1)%box)
        locBox(1)%x = ceiling(ptOnCircle(i + 1)%Pt%x)
        locBox(2)%x = ceiling(ptOnCircle(i + 1)%Pt%x)
        locBox(3)%x = floor(ptOnCircle(i + 1)%Pt%x)
        locBox(4)%x = floor(ptOnCircle(i + 1)%Pt%x)

        locBox(1)%y = floor(ptOnCircle(i + 1)%Pt%y)
        locBox(2)%y = ceiling(ptOnCircle(i + 1)%Pt%y)
        locBox(3)%y = ceiling(ptOnCircle(i + 1)%Pt%y)
        locBox(4)%y = floor(ptOnCircle(i + 1)%Pt%y)

        do a = 1, q - 1
          dirDotUnitVec(a) = (ci(a, 1)*ptOnCircle(i + 1)%unitVec%x + ci(a, 2)*ptOnCircle(i + 1)%unitVec%y) &
                             /(sqrt(ci(a, 1)**d2 + ci(a, 2)**d2))
        end do

        itmp = maxloc(dirDotUnitVec)
        outDir = itmp(1)
        ! write (*, *) outDir

        do k = 1, 4
          locBox(k)%isInside = ((locBox(k)%x - lb%xCentre)**d2 + (locBox(k)%y - lb%yCentre)**d2)**haf .le. haf*lb%diameter

          if (locBox(k)%isInside) then
            locBox(k)%fluidNode(1)%x = locBox(k)%x + ci(outDir, 1)
            locBox(k)%fluidNode(1)%y = locBox(k)%y + ci(outDir, 2)
            locBox(k)%fluidNode(2)%x = locBox(k)%x + 2*ci(outDir, 1)
            locBox(k)%fluidNode(2)%y = locBox(k)%y + 2*ci(outDir, 2)
          end if

        end do

      end associate

    end do

  end subroutine createDataStructOnCircle

  function dateTime()

    implicit none
    character(len=30)::dateTime
    character(len=10):: ampm
    integer:: d, h, m, n, s, y, mm, values(8)
    character(len=3), parameter, dimension(12) :: &
      month = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

    call date_and_time(values=values)

    y = values(1)
    m = values(2)
    d = values(3)
    h = values(5)
    n = values(6)
    s = values(7)
    mm = values(8)

    if (h < 12) then
      ampm = 'AM'
    elseif (h == 12) then
      if (n == 0 .and. s == 0) then
        ampm = 'Noon'
      else
        ampm = 'PM'
      end if
    else
      h = h - 12
      if (h < 12) then
        ampm = 'PM'
      elseif (h == 12) then
        if (n == 0 .and. s == 0) then
          ampm = 'Midnight'
        else
          ampm = 'AM'
        end if
      end if
    end if

    write (dateTime, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)') &
      d, trim(month(m)), y, h, ':', n, ':', s, '.', mm, trim(ampm)
  end function dateTime

  ! function mulMatVec(A, B)
  !   implicit none

  !   double precision, dimension(:, :), contiguous, intent(in)  :: A
  !   double precision, dimension(:), contiguous, intent(in)  :: B
  !   double precision, dimension(size(A, 1)) :: mulMatVec
  !   integer :: i, k

  !   if (size(A, 2) == size(B)) then
  !     mulMatVec = 0.0d0
  !     do k = 1, size(A, 2)
  !       do i = 1, size(A, 1)
  !         mulMatVec(i) = mulMatVec(i) + A(i, k)*B(k)
  !       end do
  !     end do
  !   else
  !     write (*, *) 'Incompatible Martix-Vector Multiplication'
  !     stop
  !   end if

  ! end function mulMatVec

  pure subroutine linearExterp(b)
    ! double precision, intent(in) :: x, y
    type(box_t), intent(inout) ::  b
    double precision:: x0x1, x2x1

    associate (f => b%fluidNode)
    if (b%x /= f(1)%x) then
      x0x1 = b%x - f(1)%x
      x2x1 = f(2)%x - f(1)%x
    else
      x0x1 = b%y - f(1)%y
      x2x1 = f(2)%y - f(1)%y
    end if

    b%f = f(1)%f + x0x1*(f(2)%f - f(1)%f)/x2x1
    b%fneq = f(1)%fneq + x0x1*(f(2)%fneq - f(1)%fneq)/x2x1
    end associate

  end subroutine linearExterp

  pure subroutine bilinearInterp(pt, box, err)
    type(triplet_t), intent(inout) :: pt
    type(box_t), intent(in) :: box(4)
    integer, optional, intent(out) :: err
    double precision::locX, locY, coeff(4)
    integer::i

    locX = pt%x - box(4)%x
    locY = pt%y - box(4)%y

    if (present(err)) then
      err = 0
      if (locX .gt. 1.0005 .or. locY .gt. 1.0005) then
        ! write (*, *) "Bilinear Interpolation: Ask-point outside box. Aborting"
        ! stop
        err = 1
      end if
    end if

    coeff = rectShape(locX, locY)

    pt%fneq = 0
    do i = 1, 4
      pt%f = pt%f + coeff(i)*box(i)%f
      pt%fneq = pt%fneq + coeff(i)*box(i)%fneq
    end do

  end subroutine bilinearInterp

  pure function rectShape(x, y) result(psi)
    implicit none

    double precision, intent(in) :: x, y
    double precision::psi(4)

    psi = [x*(1 - y), x*y, (1 - x)*y, (1 - x)*(1 - y)]

  end function rectShape

  pure function trapIntegrate(force) result(totalForce)
    implicit none
    double precision, dimension(noOfPtOnCircle), intent(in)::force
    double precision::totalForce, dx
    integer::i, ip1

    dx = pi*lb%diameter/noOfPtOnCircle

    totalForce = d0
    do i = 1, noOfPtOnCircle

      if (i == noOfPtOnCircle) then
        ip1 = 1
      else
        ip1 = i + 1
      end if

      totalForce = totalForce + haf*(force(i) + force(ip1))*dx
    end do

  end function trapIntegrate

  pure subroutine calcConvFactor(cf)
    type(cf_t), intent(out) :: cf

    cf%length = chanH/chanH_
    cf%rho = rhoF/rhoF_
    cf%time = dia/uMean*0.0025d0
    cf%nu = cf%length**2.0d0/cf%time
    cf%velocity = cf%length/cf%time
    ! cf%force = cf%rho*cf%length**4.0d0*cf%time**(-2.0d0)
    cf%force = cf%rho*cf%length**3.0d0*cf%time**(-2.0d0)

  end subroutine calcConvFactor

  pure subroutine calcLBparams(si, lb)
    type(si_lb_t), intent(out) :: si, lb

    si%channelLength = chanL
    si%channelHeight = chanH
    si%barLength = barL
    si%barHeight = barH
    si%diameter = dia
    si%xCentre = xc
    si%yCentre = yc
    si%meanVel = uMean
    si%nu = nu

    lb%channelLength = si%channelLength/cf%length
    lb%channelHeight = si%channelHeight/cf%length
    lb%diameter = si%diameter/cf%length
    lb%xCentre = si%xCentre/cf%length + 1.5d0
    lb%yCentre = si%yCentre/cf%length + 1.5d0
    lb%barLength = si%barLength/cf%length
    lb%barHeight = si%barHeight/cf%length
    lb%nu = si%nu/cf%nu
    lb%meanVel = si%meanVel/cf%velocity
    lb%tau = 3.0d0*lb%nu + 0.5d0

    lb%totalTime = time_
    si%totalTime = lb%totalTime*cf%time

  end subroutine calcLBparams

  pure subroutine setupCiWi(ci, wi, kb)
    integer, intent(out) :: ci(0:q - 1, 2), kb(0:q - 1)
    double precision, intent(out) :: wi(0:q - 1)
    integer:: ex(0:q - 1), ey(0:q - 1), a, a1

    ex(0) = 0; ey(0) = 0; 
    ex(1) = 1; ey(1) = 0; 
    ex(2) = 0; ey(2) = 1; 
    ex(3) = -1; ey(3) = 0; 
    ex(4) = 0; ey(4) = -1; 
    ex(5) = 1; ey(5) = 1; 
    ex(6) = -1; ey(6) = 1; 
    ex(7) = -1; ey(7) = -1; 
    ex(8) = 1; ey(8) = -1; 
    ci(:, 1) = ex
    ci(:, 2) = ey
!----------------------------------------------------------------------
    do a = 0, q - 1
      if (a .eq. 0) wi(a) = 16*one36th
      if (a .ge. 1 .and. a .le. 4) wi(a) = 4*one36th
      if (a .ge. 5 .and. a .le. 8) wi(a) = one36th
    end do

    do a = 0, q - 1
      do a1 = a, q - 1
        if (ci(a, 1) + ci(a1, 1) .eq. 0 .and. ci(a, 2) + ci(a1, 2) .eq. 0) then
          kb(a) = a1
          kb(a1) = a
        end if
      end do
    end do

  end subroutine setupCiWi

  pure subroutine initDistFun(f)
    double precision, dimension(:, :, :), intent(inout) :: f
    double precision:: uPara_, tmp(3)
    integer:: i, j, a

    do i = 1, nx + 2
      do j = 1, ny + 2
        uPara_ = d0! 6*lb%meanVel*(ny - (j - 1.5))*(j - 1.5)/ny**2;
        do a = 0, q - 1
          tmp(1) = uPara_*ci(a, 1)
          tmp(2) = uPara_*uPara_
          f(a, i, j) = wi(a)*rhoF_*(1.0 + 3.0*tmp(1) + 4.5*tmp(1)*tmp(1) - 1.5*tmp(2))
        end do
      end do
    end do

  end subroutine initDistFun

end program cyl
