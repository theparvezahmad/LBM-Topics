program cyl
  implicit none
  !Real units are in SI units
  !Underscore Variables are in LBM units

  integer, parameter:: &
    chanH_ = 82, &
    q = 9, &
    time_ = 100000, &
    noOfSnaps = 5, &
    dispFreq = 1

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
    pi = 4.0d0*datan(1.0d0), &
    one36th = 1.0d0/36.0d0

  integer::nx, ny
  double precision:: nu_, uMean_, uPara_, uParaRamp_, dia_, xc_, yc_, chanL_, barL_, barH_
  double precision:: Clen, Crho, Ct, Cnu, CVel, CFor, tau, t
  integer:: i, j, a, a1, t_, ia, ja, solnumber
  integer, allocatable, dimension(:, :)::isn
  double precision:: tmp1, tmp2, tmp3, rhoSum, feq, fx_t, fy_t, Cd, Cl
  double precision:: fx(2), fy(2), dudx, dudy, dvdx, dvdy, f_neq
  double precision, dimension(0:q - 1):: wt, Q_xx, Q_yy, Q_xy
  integer, dimension(0:q - 1):: ex, ey, kb
  double precision, allocatable, dimension(:, :, :):: f, ft
  double precision, allocatable, dimension(:, :):: ux, uy, rho
  logical::isCyl, isBar
  character(len=30):: filename
  solnumber = 0

  write (*, *) '======================================================'
  write (*, *) 'Program started at :', dateTime()
!===Conversion Factors===
  Clen = chanH/chanH_
  Crho = rhoF/rhoF_
  Ct = dia/uMean*0.0025d0
  Cnu = Clen**2.0d0/Ct
  CVel = Clen/Ct
  ! CFor = Crho*Clen**4.0d0*Ct**(-2.0d0)
  CFor = Crho*Clen**3.0d0*Ct**(-2.0d0)

!===Other LBM parameters===
  chanL_ = chanL/Clen
  dia_ = dia/Clen
  xc_ = xc/Clen + 1.5d0
  yc_ = yc/Clen + 1.5d0
  barL_ = barL/Clen
  barH_ = barH/Clen
  nx = int(chanL_)
  ny = chanH_
  nu_ = nu/Cnu
  uMean_ = uMean/CVel
  tau = 3*nu_ + 0.5d0
  sNu = 1.0d0/tau

  write (*, *) 'Clen = ', Clen
  write (*, *) 'Crho = ', Crho
  write (*, *) 'Ct   = ', Ct
  write (*, *) 'CVel = ', CVel
  write (*, *) 'tau  = ', tau
  write (*, *) 'uMax_ = ', 1.5d0*uMean_
  write (*, *) 'Re_SI = ', uMean*dia/nu
  write (*, *) 'Re_LBM = ', uMean_*dia_/nu_
  write (*, *) 'Aborted for Checking'; stop
!----------------------------------------------------------------------
  allocate (isn(nx + 2, ny + 2))
  allocate (ux(nx + 2, ny + 2))
  allocate (uy(nx + 2, ny + 2))
  allocate (rho(nx + 2, ny + 2))
  allocate (f(0:q - 1, nx + 2, ny + 2))
  allocate (ft(0:q - 1, nx + 2, ny + 2))
!----------------------------------------------------------------------
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
    if (a .eq. 0) wt(a) = 16*one36th
    if (a .ge. 1 .and. a .le. 4) wt(a) = 4*one36th
    if (a .ge. 5 .and. a .le. 8) wt(a) = one36th
  end do
!----------------------------------------------------------------------
  do a = 0, q - 1
    do a1 = a, q - 1
      if (ex(a) + ex(a1) .eq. 0 .and. ey(a) + ey(a1) .eq. 0) then
        kb(a) = a1
        kb(a1) = a
      end if
    end do
  end do
!----------------------------------------------------------------------
  do i = 1, nx + 2
    do j = 1, ny + 2
      uPara_ = d0! 6*uMean_*(ny - (j - 1.5))*(j - 1.5)/ny**2;
      do a = 0, q - 1
        tmp1 = uPara_*ex(a)
        tmp2 = uPara_*uPara_
        f(a, i, j) = wt(a)*rhoF_*(1.0 + 3.0*tmp1 + 4.5*tmp1*tmp1 - 1.5*tmp2)
      end do
    end do
  end do
!----------------------------------------------------------------------
  fx(1) = d0
  fy(1) = d0
  open (unit=10, file="tRhoCdCl.dat")
  write (10, *) "Variables=timeLBM,timeReal,rho,Cd,Cl"
  ! open (unit=11, file="VelSig.dat")
!----------------------------------------------------------------------
  do i = 1, nx + 1
    do j = 1, ny + 2

      ! ii = i - 1.5
      ! jj = j - 1.5

      isCyl = ((i - xc_)**2.0 + (j - yc_)**2.0)**0.5 .le. 0.5*dia_

      isBar = (i .ge. xc_) .and. (i .le. xc_ + dia_/2 + barL_) .and. &
              (j .ge. xc_ - barH_/2) .and. (j .le. xc_ + barH_/2)

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
!----------------------------------------------------------------------
  do t_ = 0, time_

    t = t_*Ct

    rhoSum = d0
    do i = 2, nx + 1
      do j = 2, ny + 1
        tmp1 = d0
        tmp2 = d0
        tmp3 = d0
        pxx = d0
        pxy = d0
        do a = 0, q - 1
          tmp1 = tmp1 + f(a, i, j)
          tmp2 = tmp2 + f(a, i, j)*ex(a)
          tmp3 = tmp3 + f(a, i, j)*ey(a)
          pxx(i, j) = pxx(i, j) + f(a, i, j)*(ex(a)**2.0d0 - ey(a)**2.0d0)
          pxy(i, j) = pxy(i, j) + f(a, i, j)*(ex(a)*ey(a))
        end do

        rho(i, j) = tmp1
        ux(i, j) = tmp2/tmp1
        uy(i, j) = tmp3/tmp1
        ! dudx(i, j) = 0.75d0*sNu*(ux(i, j)**2.0d0 - uy(i, j)**2.0d0 - pxx(i, j)/rho(i, j))
        ! dvdy(i, j) = 0.75d0*sNu*(uy(i, j)**2.0d0 - ux(i, j)**2.0d0 + pxx(i, j)/rho(i, j))
        ! dudyPdvdx(i, j) = 3.0d0*sNu*ux(i, j)*uy(i, j) - pxy(i, j)/rho(i, j)

        S11(i, j) = -1.0d0/3.0d0*rho(i, j) + 0.75d0*sNu*(ux(i, j)**2.0d0 - uy(i, j)**2.0d0 - pxx(i, j)/rho(i, j))
        ! S12(i, j) = 0.5d0*dudyPdvdx(i, j)
        S12(i, j) = 0.5d0*3.0d0*sNu*ux(i, j)*uy(i, j) - pxy(i, j)/rho(i, j)
        ! S21(i, j) = 0.5d0*dudyPdvdx(i, j)
        S22(i, j) = -1.0d0/3.0d0*rho(i, j) + 0.75d0*sNu*(uy(i, j)**2.0d0 - ux(i, j)**2.0d0 + pxx(i, j)/rho(i, j))

        rhoSum = rhoSum + tmp1
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
        do a = 0, q - 1
          tmp1 = ux(i, j)*ex(a) + uy(i, j)*ey(a)
          tmp2 = ux(i, j)*ux(i, j) + uy(i, j)*uy(i, j)
          feq = wt(a)*rho(i, j)*(1.0 + 3.0*tmp1 + 4.5*tmp1*tmp1 - 1.5*tmp2)
          ft(a, i, j) = f(a, i, j) - (f(a, i, j) - feq)/tau !collision
        end do
      end do
    end do
!----------------------------------------------------------------------
    do i = 2, nx + 1 !Streaming post-collision
      do j = 2, ny + 1
        do a = 0, q - 1
          ia = i + ex(a)
          ja = j + ey(a)
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

      uPara_ = 6.0d0*uMean_*(ny - (j - 1.5))*(j - 1.5)/ny**2.0d0; 
      if (t .lt. 2.0d0) then
        uParaRamp_ = uPara_*(1 - cos(pi*t/2.0d0))/2.0d0
      else
        uParaRamp_ = uPara_
      end if

      ux(i, j) = uParaRamp_
      uy(i, j) = 0.0
      rho(i, j) = (f(0, i, j) + f(2, i, j) + f(4, i, j) + 2.0*(f(3, i, j) + f(6, i, j) + f(7, i, j)))/(1 - ux(i, j))
      tmp2 = ux(i, j)*ux(i, j) + uy(i, j)*uy(i, j)

      dudx = (-3.0*ux(i, j) + 4.0*ux(i + 1, j) - 1.0*ux(i + 2, j))/2.0 !dudx
      if (j == 2) then
        dudy = (-3.0*ux(i, j) + 4.0*ux(i, j + 1) - 1.0*ux(i, j + 2))/2.0
      else if (j == ny + 1) then
        dudy = (3.0*ux(i, j) - 4.0*ux(i, j - 1) + 1.0*ux(i, j - 2))/2.0
      else
        dudy = (ux(i, j + 1) - ux(i, j - 1))/2.0
      end if

      dvdx = (-3.0*uy(i, j) + 4.0*uy(i + 1, j) - 1.0*uy(i + 2, j))/2.0
      if (j == 2) then
        dvdy = (-3.0*uy(i, j) + 4.0*uy(i, j + 1) - 1.0*uy(i, j + 2))/2.0
      else if (j == ny + 1) then
        dvdy = (3.0*uy(i, j) - 4.0*uy(i, j - 1) + 1.0*uy(i, j - 2))/2.0
      else
        dvdy = (uy(i, j + 1) - uy(i, j - 1))/2.0
      end if

      do a = 0, q - 1
        Q_xx(a) = ex(a)*ex(a) - 1.0/3.0
        Q_yy(a) = ey(a)*ey(a) - 1.0/3.0
        Q_xy(a) = ex(a)*ey(a)
        f_neq = -wt(a)*3.0*tau*rho(i, j)*(Q_xx(a)*dudx + Q_xy(a)*dudy + Q_xy(a)*dvdx + Q_yy(a)*dvdy)
        tmp1 = ex(a)*ux(i, j) + ey(a)*uy(i, j)
        feq = wt(a)*rho(i, j)*(1 + tmp1*3.0 + 0.5*tmp1*tmp1*3.0*3.0 - 0.5*3.0*tmp2)
        f(a, i, j) = feq + f_neq
      end do

      !Outlet
      i = nx + 1
      rho(i, j) = rhoF_
      ux(i, j) = (f(0, i, j) + f(2, i, j) + f(4, i, j) + 2*(f(1, i, j) + f(5, i, j) + f(8, i, j)))/rho(i, j) - 1
      uy(i, j) = 0.0
      tmp2 = ux(i, j)*ux(i, j) + uy(i, j)*uy(i, j)

      dudx = (3.0*ux(i, j) - 4.0*ux(i - 1, j) + 1.0*ux(i - 2, j))/2.0 !dudx
      if (j == 2) then
        dudy = (-3.0*ux(i, j) + 4.0*ux(i, j + 1) - 1.0*ux(i, j + 2))/2.0
      else if (j == ny + 1) then
        dudy = (3.0*ux(i, j) - 4.0*ux(i, j - 1) + 1.0*ux(i, j - 2))/2.0
      else
        dudy = (ux(i, j + 1) - ux(i, j - 1))/2.0
      end if

      dvdx = (3.0*uy(i, j) - 4.0*uy(i - 1, j) + 1.0*uy(i - 2, j))/2.0 !dvdx
      if (j == 2) then
        dvdy = (-3.0*uy(i, j) + 4.0*uy(i, j + 1) - 1.0*uy(i, j + 2))/2.0
      else if (j == ny + 1) then
        dvdy = (3.0*uy(i, j) - 4.0*uy(i, j - 1) + 1.0*uy(i, j - 2))/2.0
      else
        dvdy = (uy(i, j + 1) - uy(i, j - 1))/2.0
      end if

      do a = 0, q - 1
        Q_xx(a) = ex(a)*ex(a) - 1.0/3.0
        Q_yy(a) = ey(a)*ey(a) - 1.0/3.0
        Q_xy(a) = ex(a)*ey(a)
        f_neq = -wt(a)*3.0*tau*rho(i, j)*(Q_xx(a)*dudx + Q_xy(a)*dudy + Q_xy(a)*dvdx + Q_yy(a)*dvdy)
        tmp1 = ex(a)*ux(i, j) + ey(a)*uy(i, j)
        feq = wt(a)*rho(i, j)*(1 + tmp1*3.0 + 0.5*tmp1*tmp1*3.0*3.0 - 0.5*3.0*tmp2)
        f(a, i, j) = feq + f_neq
      end do
    end do
!----------------------------------------------------------------------
    ! do j = 2, ny + 1
    !    i = 2
    !    uPara_ = 6.0d0*uMean_*(ny - (j - 1.5))*(j - 1.5)/ny**2.0d0;
    !    if (t .lt. 2.0d0) then
    !       uParaRamp_ = uPara_*(1 - cos(pi*t/2.0d0))/2.0d0
    !    else
    !       uParaRamp_ = uPara_
    !    end if
    !    rho(i, j) = (f(0, i, j) + f(2, i, j) + f(4, i, j) + 2*(f(6, i, j) + f(3, i, j) + f(7, i, j)))/(1 - uParaRamp_)
    !    f(1, i, j) = f(3, i, j) + ((2.0/3.0)*rho(i, j)*uParaRamp_)
    !    f(5, i, j) = f(7, i, j) - (0.5*(f(2, i, j) - f(4, i, j))) + ((1.0/6.0)*rho(i, j)*uParaRamp_)
    !    f(8, i, j) = f(6, i, j) + (0.5*(f(2, i, j) - f(4, i, j))) + ((1.0/6.0)*rho(i, j)*uParaRamp_)
    ! end do

    ! do j = 2, ny + 1
    !    i = nx + 1
    !    ux(i, j) = (f(0, i, j) + f(2, i, j) + f(4, i, j) + 2*(f(1, i, j) + f(5, i, j) + f(8, i, j)))/rhoF_ - 1
    !    f(3, i, j) = f(1, i, j) - ((2.0/3.0)*rhoF_*ux(i, j))
    !    f(6, i, j) = f(8, i, j) - 0.5*(f(2, i, j) - f(4, i, j)) - ((1.0/6.0)*rhoF_*ux(i, j))
    !    f(7, i, j) = f(5, i, j) + 0.5*(f(2, i, j) - f(4, i, j)) - ((1.0/6.0)*rhoF_*ux(i, j))
    ! end do
!----------------------------------------------------------------------
    fx(2) = d0
    fy(2) = d0
    do i = 2, nx + 1 !BC
      do j = 2, ny + 1
        if (isn(i, j) .eq. 0) then

          do a = 0, q - 1

            ia = i + ex(a)
            ja = j + ey(a)

            if (isn(ia, ja) .eq. 1) then !structure
              f(kb(a), i, j) = ft(a, i, j)
              f(a, ia, ja) = ft(kb(a), ia, ja)
              fx(2) = fx(2) + ex(a)*2.0*(-ft(kb(a), ia, ja) + ft(a, i, j))
              fy(2) = fy(2) + ey(a)*2.0*(-ft(kb(a), ia, ja) + ft(a, i, j))
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
    ! Cd = fx_t*CFor!/(0.5*rhoF_*uMean_*uMean_*dia_)
    ! Cl = fy_t*CFor!/(0.5*rhoF_*uMean_*uMean_*dia_)
    Cd = fx_t/(0.5*rhoF_*uMean_*uMean_*dia_)
    Cl = fy_t/(0.5*rhoF_*uMean_*uMean_*dia_)
!----------------------------------------------------------------------
    if (mod(t_, dispFreq) .eq. 0) then
      write (10, '(E16.6,2X,I10,3(2X,E12.4))') t, t_, rhoSum, Cd, Cl
      !write (*, '(I8,4(3X,F10.6))') ts, ts*Ct, rhoAvg, Cd, Cl

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
          write (12, '(2(2X,I6),3(2X,E12.4),2X,I3)') i, j, ux(i, j), uy(i, j), rho(i, j), isn(i, j)
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

  subroutine calcStressTensor(f, sigma)
    implicit none

    double precision, dimension(:), intent(in) :: f
    double precision, dimension(2, 2), intent(out) ::  sigma
    double precision:: u(2), tmp(2)

    tmp = zero
    do a = 0, q - 1
      tmp(1) = tmp(1) + f(i)
      tmp(2) = tmp(2) + f(i)*ex(a)
      tmp(3) = tmp(3) + f(i)*ey(a)
    end do

    rho = tmp(1)
    u(1) = tmp(2)/rho
    u(2) = tmp(3)/rho

    do i = 1, 2
      do j = 1, 2

        tmp(1) = zero
        do a = 0, q - 1
          tmp(1) = tmp(1) + (ci(a, i) - u(i))*(ci(a, j) - u(j))*f(a)
        end do

        sigma(i, j) = -one6th*invTau*rho*kdf(i, j) - (one - haf*invTau)*tmp(1)

      end do
    end do

  end subroutine calcStressTensor

  function kdf(i, j) result(retval)
    integer, intent(in) :: i, j
    integer :: retval

    if (i == j) then
      retval = 1
    else
      retval = 0
    end if
  end function kdf

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

end program cyl
