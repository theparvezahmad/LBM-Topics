program cyl
  implicit none
  !Real units are in SI units
  !Underscore Variables are in LBM units

  integer, parameter:: &
    chanH_ = 82, &
    dim = 2, &
    q = 9, &
    time_ = 100000, &
    noOfSnaps = 5, &
    dispFreq = 100, &
    noOfPtOnCircle = 400, &
    noOfPtOnBar = 592

  double precision, parameter:: &
    rhoF_ = 1.0d0, &
    rhoF = 1000.0d0, &
    chanL = 2.5d0, & !Length of channel
    chanH = 0.41d0, & !Width of channel
    barL = 0.35d0, &
    barH = 0.02d0, &
    uMean = 0.2d0, &
    nu = 0.001d0, &
    dia = 0.1d0, &
    ! diaNew = 1.32d0*0.1d0, &
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
    double precision::x, y, z
  end type triplet_t
  type doubletInt_t
    integer::x, y
  end type doubletInt_t

  ! type triplet_t
  !   double precision::x, y, z(0:q - 1)
  ! end type triplet_t

  ! type, extends(triplet_t):: box_t
  !   ! double precision::x, y, z(0:q - 1)
  !   type(triplet_t)::fluidNode(2)
  !   logical:: isInside
  ! end type box_t

  type tensor2D_t
    double precision::t11, t12, t21, t22
  end type

  type ptWithTensorN0_t
    double precision::x, y
    ! type(tensor2D_t)::st
    double precision::st(2, 2)
  end type

  type ptWithTensorN1_t
    integer::x, y
    ! type(tensor2D_t)::st
    double precision::st(2, 2)
  end type

  type ptWithTensorN2_t
    double precision::x, y
    ! type(tensor2D_t)::st
    double precision::st(2, 2)
    type(doubletInt_t)::b(2)
  end type

  type custom_t
    integer::noOfRD
    type(doublet_t)::uv, force
    type(ptWithTensorN0_t), dimension(3)::n0
    type(ptWithTensorN1_t), dimension(3)::n1
    type(ptWithTensorN2_t), dimension(3)::n2
    type(ptWithTensorN2_t), dimension(3)::n3
  end type custom_t

  type(custom_t), allocatable, dimension(:)::ptOnCircle, ptOnBar
  type(doublet_t)::totalForceCir, totalForceBar
  type(triplet_t)::dataPt1, dataPt2, dataPt3, dataPt, dataPtArr(0:2)
  type(lbmTriplet_t)::lbmVarA(noOfPtOnBar)
  integer::nx, ny
  double precision, dimension(noOfPtOnCircle):: forceXCir, forceYCir
  double precision, dimension(noOfPtOnBar):: forceXBar, forceYBar
  double precision:: nu_, uMean_, uPara_, uParaRamp_, dia_, xc_, yc_, chanL_, barL_, barH_, arcSkip
  double precision:: Clen, Crho, Ct, Cnu, CVel, CFor, tau, t, invTau, sigma(2, 2), avgST(2, 2)
  integer:: i, iRD, j, k, p, a, a1, t_, ia, ja, solnumber, cnt
  double precision:: i_, j_, ia_, ja_, ii, jj, Delta(1000), cDelta, chi, ubfx, ubfy, uwx, uwy, iMid, jMid
  integer, allocatable, dimension(:, :)::isn, wb
  double precision:: tmp1, tmp2, tmp3, rhoSum, feq, fx_t, fy_t, Cd, Cl, Cd2, Cl2
  double precision:: fx(2), fy(2), dudx, dudy, dvdx, dvdy, f_neq
  double precision, dimension(0:q - 1):: wt, Q_xx, Q_yy, Q_xy, tmpA
  integer, dimension(0:q - 1):: ex, ey, kb, ci(0:q - 1, 2)
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
  ! diaNew_ = diaNew/Clen
  xc_ = xc/Clen + 1.5d0
  yc_ = yc/Clen + 1.5d0
  barL_ = barL/Clen
  barH_ = barH/Clen
  nx = int(chanL_)
  ny = chanH_
  nu_ = nu/Cnu
  uMean_ = uMean/CVel
  tau = 3*nu_ + 0.5d0
  invTau = 1.0d0/tau
  arcSkip = atan(barH_/dia_)

  write (*, *) 'Clen = ', Clen
  write (*, *) 'Crho = ', Crho
  write (*, *) 'Ct   = ', Ct
  write (*, *) 'CVel = ', CVel
  write (*, *) 'tau  = ', tau
  write (*, *) 'uMax_ = ', 1.5d0*uMean_
  write (*, *) 'Re_SI = ', uMean*dia/nu
  write (*, *) 'Re_LBM = ', uMean_*dia_/nu_
  ! write (*, *) 'Aborted for Checking'; stop
!----------------------------------------------------------------------
  allocate (isn(nx + 2, ny + 2))
  allocate (wb(nx + 2, ny + 2))
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
  write (10, *) "Variables=timeLBM,timeReal,rho,Cd,Cl,Cd2,Cl2"
  ! open (unit=11, file="VelSig.dat")
!----------------------------------------------------------------------
  do i = 1, nx + 1
    do j = 1, ny + 2

      ! ii = i - 1.5
      ! jj = j - 1.5

      ! isCyl = ((i - xc_)**2.0 + (j - yc_)**2.0)**0.5 .le. 0.5*dia_
      isCyl = .false.

      isBar = (i .ge. xc_ + haf*dia_) .and. (i .le. xc_ + dia_/2 + barL_) .and. &
              (j .ge. xc_ - barH_/2) .and. (j .le. xc_ + barH_/2)
      ! isBar = .false.

      if (isCyl) then
        isn(i, j) = 1
      elseif (isBar) then
        isn(i, j) = 3
      elseif (j == 1 .or. j == ny + 2) then
        isn(i, j) = 2
      else
        isn(i, j) = 0
      end if

    end do
  end do

  open (unit=12, file="region.dat")
  write (12, *) "Variables=x,y,region"
  write (12, '(2(a,I6))') "Zone I=", nx + 2, ",J=", ny + 2
  do j = 1, ny + 2
    do i = 1, nx + 2
      write (12, '(2(2X,I6),2X,I3)') i, j, isn(i, j)
    end do
    write (12, *)
  end do
  close (12)
  ! stop
  !------------------------------------------------------------------
  cnt = 0
  do i = 2, nx + 1
    do j = 2, ny + 1
      if (isn(i, j) == 0) then
        do a = 0, q - 1
          ia = i + ex(a)
          ja = j + ey(a)
          if (isn(ia, ja) == 1) then !Cylinder
            wb(ia, ja) = 1
            i_ = i
            j_ = j
            ia_ = ia
            ja_ = ja
            ii = (i_ + ia_)/2.0d0
            jj = (j_ + ja_)/2.0d0
            !printf("%s\n","A" )
            do while (dabs((ii - xc_)*(ii - xc_) + (jj - yc_)*(jj - yc_) - dia_*dia_/4.0d0) .gt. 0.001d0)
              !printf("%s\n","in" )
              if ((ii - xc_)*(ii - xc_) + (jj - yc_)*(jj - yc_) - dia_*dia_/4.0d0 .le. 0.0d0) then
                ia_ = ii
                ja_ = jj
              else
                i_ = ii
                j_ = jj
              end if

              ii = (i_ + ia_)/2.0d0
              jj = (j_ + ja_)/2.0d0
            end do

            cnt = cnt + 1
            Delta(cnt) = dsqrt((i - ii)*(i - ii) + dble(j - jj)*(j - jj))/dsqrt(dble(i - ia)*(i - ia) + dble(j - ja)*(j - ja))
            ! write(*,'(6(I4),3(F8.4))') cnt,a,i,j,ia,ja,ii,jj,Delta(cnt)
          end if
        end do
      end if
    end do
  end do
  ! Delta = 0.5d0;
  write (*, *) "Total number of interactions: ", cnt
  ! stop

  call createDataStructOnCircle(noOfPtOnCircle, ptOnCircle)
  call createDSonBar(noOfPtOnBar, ptOnBar)

!----------------------------------------------------------------------
  do t_ = 0, time_

    t = t_*Ct

    rhoSum = d0
    do i = 2, nx + 1
      do j = 2, ny + 1
        tmp1 = d0
        tmp2 = d0
        tmp3 = d0
        do a = 0, q - 1
          tmp1 = tmp1 + f(a, i, j)
          tmp2 = tmp2 + f(a, i, j)*ex(a)
          tmp3 = tmp3 + f(a, i, j)*ey(a)
        end do

        rho(i, j) = tmp1
        ux(i, j) = tmp2/tmp1
        uy(i, j) = tmp3/tmp1

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
    fx(2) = d0
    fy(2) = d0
    uwx = d0
    uwy = d0
    cnt = 0
    do i = 2, nx + 1 !BC
      do j = 2, ny + 1
        if (isn(i, j) .eq. 0) then

          do a = 0, q - 1

            ia = i + ex(a)
            ja = j + ey(a)

            iMid = haf*(i + ia)
            jMid = haf*(j + ja)

            if (isn(ia, ja) .eq. 1) then !structure
              cnt = cnt + 1
              if (Delta(cnt) .ge. d0 .and. Delta(cnt) .lt. 0.5d0) then
                chi = (2.0d0*Delta(cnt) - 1.0d0)/(tau - 2.0d0)
                ubfx = ux(i - ex(a), j - ey(a))
                ubfy = uy(i - ex(a), j - ey(a))
              end if

              if (Delta(cnt) .ge. 0.5d0 .and. Delta(cnt) .lt. 1.0d0) then
                chi = (2.0d0*Delta(cnt) - 1.0d0)/(tau + 0.5d0)
                ubfx = 0.5d0*(2.0d0*Delta(cnt) - 3.0d0)/Delta(cnt)*ux(i, j) + 1.5d0*uwx/Delta(cnt)
                ubfy = 0.5d0*(2.0d0*Delta(cnt) - 3.0d0)/Delta(cnt)*uy(i, j) + 1.5d0*uwy/Delta(cnt)
              end if

              tmp1 = ux(i, j)*ex(a) + uy(i, j)*ey(a)
              tmp2 = ux(i, j)*ux(i, j) + uy(i, j)*uy(i, j)
              feq = wt(a)*rho(i, j)*(1.0d0 + 3.0d0*tmp1 + 4.5d0*tmp1*tmp1 - 1.5d0*tmp2)
              ft(kb(a), ia, ja) = ft(a, i, j) - chi*(ft(a, i, j) - feq) + &
                                  wt(a)*rho(i, j)*3.0d0*(ex(a)*(chi*(ubfx - ux(i, j)) - 2.0d0*uwx) + &
                                                         ey(a)*(chi*(ubfy - uy(i, j)) - 2.0d0*uwy))

              fx(2) = fx(2) + ex(a)*(ft(a, i, j) + ft(kb(a), ia, ja))
              fy(2) = fy(2) + ey(a)*(ft(a, i, j) + ft(kb(a), ia, ja))
            end if

            if (isn(ia, ja) .eq. 3) then !bar
              ! if ((iMid .gt. xc_ + haf*dia_) .and. (iMid .lt. xc_ + haf*dia_ + barL_) .and. (jMid .eq. yc_)) then
              ! write (*, *) t_, a, iMid, jMid
              cDelta = 0.5d0

              chi = (2.0d0*cDelta - 1.0d0)/(tau + 0.5d0)
              ubfx = 0.5d0*(2.0d0*cDelta - 3.0d0)/cDelta*ux(i, j) + 1.5d0*uwx/cDelta
              ubfy = 0.5d0*(2.0d0*cDelta - 3.0d0)/cDelta*uy(i, j) + 1.5d0*uwy/cDelta

              tmp1 = ux(i, j)*ex(a) + uy(i, j)*ey(a)
              tmp2 = ux(i, j)*ux(i, j) + uy(i, j)*uy(i, j)
              feq = wt(a)*rho(i, j)*(1.0d0 + 3.0d0*tmp1 + 4.5d0*tmp1*tmp1 - 1.5d0*tmp2)
              ft(kb(a), ia, ja) = ft(a, i, j) - chi*(ft(a, i, j) - feq) + &
                                  wt(a)*rho(i, j)*3.0d0*(ex(a)*(chi*(ubfx - ux(i, j)) - 2.0d0*uwx) + &
                                                         ey(a)*(chi*(ubfy - uy(i, j)) - 2.0d0*uwy))

              fx(2) = fx(2) + ex(a)*(ft(a, i, j) + ft(kb(a), ia, ja))
              fy(2) = fy(2) + ey(a)*(ft(a, i, j) + ft(kb(a), ia, ja))
            end if

            if (isn(ia, ja) .eq. 2) then !wall
              f(kb(a), i, j) = ft(a, i, j)
            end if

          end do
        end if
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
    ! call calcStressTensor2(f(:, int(xc_), int(yc_ + dia_)), sigma)
    ! write (*, '(4(2X,E12.4))') sigma(1, 1), sigma(1, 2), sigma(2, 1), sigma(2, 2)

    ! dataPt2%x = 59.0d0
    ! dataPt2%y = 11.0d0
    ! dataPt2%z = 300.0d0

    ! dataPt1%x = 60.0d0
    ! dataPt1%y = 10.0d0
    ! dataPt1%z = 200.0d0

    ! dataPt%x = 60.5d0
    ! dataPt%y = 9.5d0
    ! call linearExtInt(dataPt1, dataPt2, dataPt)
    ! write (*, *) dataPt%z
    ! stop

    ! dataPt1%x = 0.0d0
    ! dataPt1%y = 3.0d0
    ! dataPt1%z = 9.0d0

    ! dataPt2%x = 0.0d0
    ! dataPt2%y = 2.0d0
    ! dataPt2%z = 4.0d0

    ! dataPt3%x = 0.0d0
    ! dataPt3%y = 4.0d0
    ! dataPt3%z = 16.0d0

    ! dataPt%x = 0.0d0
    ! dataPt%y = 1.5d0
    ! call quadExtInt(dataPt1, dataPt2, dataPt3, dataPt)
    ! write (*, *) dataPt%z
    ! stop
!=================working================================
    do i = 1, size(ptOnCircle)
      associate (poc => ptOnCircle(i))

        do iRD = 1, poc%noOfRD

          dataPt1%x = poc%n2(iRD)%b(1)%x
          dataPt1%y = poc%n2(iRD)%b(1)%y
          dataPt2%x = poc%n2(iRD)%b(2)%x
          dataPt2%y = poc%n2(iRD)%b(2)%y
          dataPt%x = poc%n2(iRD)%x
          dataPt%y = poc%n2(iRD)%y
          do a = 0, q - 1
            dataPt1%z = f(a, int(dataPt1%x), int(dataPt1%y))
            dataPt2%z = f(a, int(dataPt2%x), int(dataPt2%y))
            call linearExtInt(dataPt1, dataPt2, dataPt)
            tmpA(a) = dataPt%z
          end do

          call calcStressTensor2(tmpA, sigma)
          poc%n2(iRD)%st = sigma
          call calcStressTensor2(f(:, int(poc%n1(iRD)%x), int(poc%n1(iRD)%y)), sigma)
          poc%n1(iRD)%st = sigma

          dataPt1%x = poc%n2(iRD)%x
          dataPt1%y = poc%n2(iRD)%y
          dataPt2%x = poc%n1(iRD)%x
          dataPt2%y = poc%n1(iRD)%y
          dataPt%x = poc%n0(iRD)%x
          dataPt%y = poc%n0(iRD)%y
          do k = 1, 2
            do p = 1, 2
              dataPt1%z = poc%n2(iRD)%st(k, p)
              dataPt2%z = poc%n1(iRD)%st(k, p)
              call linearExtInt(dataPt1, dataPt2, dataPt)
              poc%n0(iRD)%st(k, p) = dataPt%z
            end do
          end do

        end do

        do k = 1, 2
          do p = 1, 2
            avgST(k, p) = sum(poc%n0(1:poc%noOfRD)%st(k, p))/poc%noOfRD
          end do
        end do
        ! do k = 1, 4

        !   associate (lf => poc%box(k)%fluidNode, lb => poc%box(k))

        !     ! lf(2)%x = 45.0d0
        !     ! lf(2)%y = 44.0d0
        !     ! lf(2)%z = [20.0d0, 20.0d0, 20.0d0, 30.0d0, 20.0d0, 20.0d0, 20.0d0, 20.0d0, 200.0d0]

        !     ! lf(1)%x = 45.0d0
        !     ! lf(1)%y = 45.0d0
        !     ! lf(1)%z = [15.0d0, 20.0d0, 10.0d0, 10.0d0, 10.0d0, 10.0d0, 10.0d0, 10.0d0, 100.0d0]

        !     ! lb%x = 45.0d0
        !     ! lb%y = 46.0d0
        !     ! call linearExterp(lb)
        !     ! write (*, *) "Hello", lb%z
        !     ! stop

        !     if (lb%isInside) then
        !       lf(1)%z = f(:, int(lf(1)%x), int(lf(1)%y))
        !       lf(2)%z = f(:, int(lf(2)%x), int(lf(2)%y))
        !       call linearExterp(lb)
        !     else
        !       lb%z = f(:, int(lb%x), int(lb%y))
        !     end if

        !   end associate

        ! end do

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

        ! call bilinearInterp(poc%Pt, poc%box)

        ! call calcStressTensor2(poc%Pt%z, sigma, onSurf)

        poc%force%x = avgST(1, 1)*poc%uv%x + avgST(1, 2)*poc%uv%y
        ! - onSurf%r*onSurf%u*(onSurf%u*poc%uv%x + onSurf%v*poc%uv%y)
        poc%force%y = avgST(2, 1)*poc%uv%x + avgST(2, 2)*poc%uv%y
        ! - onSurf%r*onSurf%v*(onSurf%u*poc%uv%x + onSurf%v*poc%uv%y)

        ! if (t_ == 100000) then
        !   write (*, *) i, poc%force
        ! end if
      end associate
    end do

    do i = 1, size(ptOnBar)
      associate (pob => ptOnBar(i))

        do iRD = 1, pob%noOfRD

          dataPt1%x = pob%n3(iRD)%b(1)%x
          dataPt1%y = pob%n3(iRD)%b(1)%y
          dataPt2%x = pob%n3(iRD)%b(2)%x
          dataPt2%y = pob%n3(iRD)%b(2)%y
          dataPt%x = pob%n3(iRD)%x
          dataPt%y = pob%n3(iRD)%y
          do a = 0, q - 1
            dataPt1%z = f(a, int(dataPt1%x), int(dataPt1%y))
            dataPt2%z = f(a, int(dataPt2%x), int(dataPt2%y))
            call linearExtInt(dataPt1, dataPt2, dataPt)
            tmpA(a) = dataPt%z
          end do
          call calcStressTensor2(tmpA, sigma)
          pob%n3(iRD)%st = sigma

          dataPt1%x = pob%n2(iRD)%b(1)%x
          dataPt1%y = pob%n2(iRD)%b(1)%y
          dataPt2%x = pob%n2(iRD)%b(2)%x
          dataPt2%y = pob%n2(iRD)%b(2)%y
          dataPt%x = pob%n2(iRD)%x
          dataPt%y = pob%n2(iRD)%y
          do a = 0, q - 1
            dataPt1%z = f(a, int(dataPt1%x), int(dataPt1%y))
            dataPt2%z = f(a, int(dataPt2%x), int(dataPt2%y))
            call linearExtInt(dataPt1, dataPt2, dataPt)
            tmpA(a) = dataPt%z
          end do
          call calcStressTensor2(tmpA, sigma)
          pob%n2(iRD)%st = sigma

          call calcStressTensor2(f(:, int(pob%n1(iRD)%x), int(pob%n1(iRD)%y)), sigma)
          pob%n1(iRD)%st = sigma

          dataPt3%x = pob%n3(iRD)%x
          dataPt3%y = pob%n3(iRD)%y
          dataPt2%x = pob%n2(iRD)%x
          dataPt2%y = pob%n2(iRD)%y
          dataPt1%x = pob%n1(iRD)%x
          dataPt1%y = pob%n1(iRD)%y
          dataPt%x = pob%n0(iRD)%x
          dataPt%y = pob%n0(iRD)%y
          do k = 1, 2
            do p = 1, 2
              dataPt1%z = pob%n1(iRD)%st(k, p)
              dataPt2%z = pob%n2(iRD)%st(k, p)
              dataPt3%z = pob%n3(iRD)%st(k, p)
              ! call linearExtInt(dataPt1, dataPt2, dataPt)
              call quadExtInt(dataPt1, dataPt2, dataPt3, dataPt)
              pob%n0(iRD)%st(k, p) = dataPt%z
            end do
          end do

        end do

        do k = 1, 2
          do p = 1, 2
            avgST(k, p) = sum(pob%n0(1:pob%noOfRD)%st(k, p))/pob%noOfRD
          end do
        end do

        pob%force%x = avgST(1, 1)*pob%uv%x + avgST(1, 2)*pob%uv%y
        pob%force%y = avgST(2, 1)*pob%uv%x + avgST(2, 2)*pob%uv%y

        ! if (t_ == 100000) then
        !   write (*, *) i, pob%force
        ! end if
      end associate
    end do

    ! ptOnCircle%force%x = 1.75
    ! write (*, *) dia_
    ! write (*, *) trapIntegrate(ptOnCircle%force%x)
    ! stop

    do k = 1, noOfPtOnCircle
      forceXCir(k) = ptOnCircle(k)%force%x
      forceYCir(k) = ptOnCircle(k)%force%y
    end do

    do k = 1, noOfPtOnBar
      forceXBar(k) = ptOnBar(k)%force%x
      forceYBar(k) = ptOnBar(k)%force%y
    end do
    ! write (*, *) forceX
    ! totalForceCir%x = trapIntegrate(forceXCir)
    ! totalForceCir%y = trapIntegrate(forceYCir)
    totalForceCir%x = sum(forceXCir)*(2*pi - 2*arcSkip)*haf*dia_/noOfPtOnCircle
    totalForceCir%y = sum(forceYCir)*(2*pi - 2*arcSkip)*haf*dia_/noOfPtOnCircle

    totalForceBar%x = sum(forceXBar)*(2*barL_ + 2*barH_)/noOfPtOnBar
    totalForceBar%y = sum(forceYBar)*(2*barL_ + 2*barH_)/noOfPtOnBar

    ! if (t_ == 100000) then
    !   write (*, *) '======================================'
    !   write (*, *) totalForce%x, sum(forceX)
    !   write (*, *) 0.5*rhoF_*uMean_*uMean_*dia_
    ! end if

    if (mod(t_, 500) == 0) then
      do i = 71, 74
        associate (pob => ptOnBar(i))
          iRD = 1
          dataPt1%x = pob%n2(iRD)%x
          dataPt1%y = pob%n2(iRD)%y
          dataPt2%x = pob%n1(iRD)%x
          dataPt2%y = pob%n1(iRD)%y
          dataPt%x = pob%n0(iRD)%x
          dataPt%y = pob%n0(iRD)%y
          do a = 0, q - 1
            dataPt1%z = f(a, int(dataPt1%x), int(dataPt1%y))
            dataPt2%z = f(a, int(dataPt2%x), int(dataPt2%y))
            call linearExtInt(dataPt1, dataPt2, dataPt)
            tmpA(a) = dataPt%z
          end do
          lbmVarA(i) = macroVar(tmpA)
        end associate
      end do
      write (*, '(I10,4(F8.4))') t_, lbmVarA(71)%r, lbmVarA(72)%r, lbmVarA(73)%r, lbmVarA(74)%r
    end if

!=================working================================
    fx_t = 0.5*(fx(1) + fx(2))
    fy_t = 0.5*(fy(1) + fy(2))
    Cd = fx_t*CFor!/(0.5*rhoF_*uMean_*uMean_*dia_)
    Cl = fy_t*CFor!/(0.5*rhoF_*uMean_*uMean_*dia_)
    ! Cd2 = (totalForceCir%x + totalForceBar%x)*CFor!/(0.5*rhoF_*uMean_*uMean_*dia_)
    Cd2 = (totalForceBar%x)*CFor!/(0.5*rhoF_*uMean_*uMean_*dia_)
    ! Cl2 = (totalForcecir%y + totalForceBar%y)*CFor!/(0.5*rhoF_*uMean_*uMean_*dia_)
    Cl2 = (totalForceBar%y)*CFor!/(0.5*rhoF_*uMean_*uMean_*dia_)

    ! Cd = fx_t/(0.5*rhoF_*uMean_*uMean_*dia_)
    ! Cd2 = totalForce%x/(0.5*rhoF_*uMean_*uMean_*dia_)
    ! Cl = fy_t/(0.5*rhoF_*uMean_*uMean_*dia_)
    ! Cl2 = totalForce%y/(0.5*rhoF_*uMean_*uMean_*dia_)
!----------------------------------------------------------------------
    if (mod(t_, dispFreq) .eq. 0) then
      write (10, '(E16.6,2X,I10,7(2X,E12.4))') t, t_, rhoSum, Cd, Cl, Cd2, Cl2
      ! write (*, '(E16.6,2X,I10,5(2X,E12.4))') t, t_, rhoSum, Cd, Cl, Cd2, Cl2
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

  pure subroutine calcStressTensor(f, sigma)
    implicit none

    double precision, dimension(0:q - 1), intent(in) :: f
    double precision, intent(out) ::  sigma(2, 2)
    ! type(lbmTriplet_t), intent(out) :: onSurf
    double precision:: u(2), tmp(3), rho
    integer:: i, j, a

    tmp = d0
    do a = 0, q - 1
      tmp(1) = tmp(1) + f(a)
      tmp(2) = tmp(2) + f(a)*ex(a)
      tmp(3) = tmp(3) + f(a)*ey(a)
    end do

    rho = tmp(1)
    u(1) = tmp(2)/rho
    u(2) = tmp(3)/rho

    ! onSurf%r = rho
    ! onSurf%u = u(1)
    ! onSurf%v = u(2)

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

  pure subroutine calcStressTensor2(f, sigma)
    implicit none

    double precision, dimension(0:q - 1), intent(in) :: f
    double precision, intent(out) ::  sigma(2, 2)
    type(lbmTriplet_t) :: lbmVar

    double precision, dimension(0:q - 1) :: feq, fneq
    double precision:: tmp(3)
    integer:: i, j, a

    lbmVar = macroVar(f)

    feq = eqDist(lbmVar)

    fneq = f - feq

    sigma = d0
    do i = 1, 2
      do j = 1, 2

        tmp(1) = d0
        do a = 0, q - 1
          tmp(1) = tmp(1) + fneq(a)*(ci(a, i)*ci(a, j) - haf*(ci(a, 1)*ci(a, 1) + ci(a, 2)*ci(a, 2))*kdf(i, j))
        end do

        sigma(i, j) = -(d1 - haf*invTau)*tmp(1)

      end do
    end do

    sigma(1, 1) = sigma(1, 1) - one3rd*lbmVar%r
    sigma(2, 2) = sigma(2, 2) - one3rd*lbmVar%r

  end subroutine calcStressTensor2

  pure function eqDist(lbmVar) result(feq)
    ! double precision, dimension(0:q - 1), intent(in) :: f
    type(lbmTriplet_t), intent(in) :: lbmVar
    double precision, dimension(0:q - 1) :: feq
    integer:: a
    double precision:: tmp1, tmp2

    do a = 0, q - 1
      tmp1 = lbmVar%u*ex(a) + lbmVar%v*ey(a)
      tmp2 = lbmVar%u**d2 + lbmVar%v**d2
      feq(a) = wt(a)*lbmVar%r*(1.0 + 3.0*tmp1 + 4.5*tmp1*tmp1 - 1.5*tmp2)
      ! ft(a) = f(a) - (f(a) - feq)/tau !collision
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
      tmp(2) = tmp(2) + f(a)*ex(a)
      tmp(3) = tmp(3) + f(a)*ey(a)
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
    type(doublet_t)::dir
    type(doubletInt_t)::locBox(4)
    double precision::theta0, theta, dTheta, dirDotUnitVec(4)
    double precision:: x0, x1, y0, y1, a1, b1, c1, a2, b2, c2, determinant
    double precision:: x_xConst, y_xConst, x_yConst, y_yConst
    double precision:: dist_xConst, dist_yConst
    integer::i, a, k, outDir(3), xConst, yConst!, itmp(1)
    allocate (ptOnCircle(noOfPts))

    theta0 = d0 + arcSkip
    dTheta = (2*pi - 2*arcSkip)/noOfPts

    ! theta0 = d0
    ! dTheta = (2*pi)/noOfPts

    do i = 0, noOfPts - 1
      associate (poc => ptOnCircle(i + 1))
        theta = theta0 + (i + haf)*dTheta
        poc%n0(:)%x = xc_ + 0.5*dia_*cos(theta)
        poc%n0(:)%y = yc_ + 0.5*dia_*sin(theta)

        poc%uv%x = cos(theta)
        poc%uv%y = sin(theta)

        locBox(1)%x = ceiling(poc%n0(1)%x)
        locBox(2)%x = ceiling(poc%n0(1)%x)
        locBox(3)%x = floor(poc%n0(1)%x)
        locBox(4)%x = floor(poc%n0(1)%x)

        locBox(1)%y = floor(poc%n0(1)%y)
        locBox(2)%y = ceiling(poc%n0(1)%y)
        locBox(3)%y = ceiling(poc%n0(1)%y)
        locBox(4)%y = floor(poc%n0(1)%y)

        do a = 1, 4
          dir%x = locBox(a)%x - poc%n0(1)%x
          dir%y = locBox(a)%y - poc%n0(1)%y
          dirDotUnitVec(a) = (dir%x*poc%uv%x + dir%y*poc%uv%y) &
                             /(sqrt(dir%x**d2 + dir%y**d2))
        end do

        poc%noOfRD = 0
        do a = 1, 4
          ! if (dirDotUnitVec(a) .eq. maxval(dirDotUnitVec)) then
          if (dirDotUnitVec(a) .gt. d0) then
            poc%noOfRD = poc%noOfRD + 1
            outDir(poc%noOfRD) = a
            ! exit
          end if
        end do

        ! itmp = maxloc(dirDotUnitVec)
        ! outDir = itmp(1)
        do k = 1, poc%noOfRD
          dir%x = locBox(outDir(k))%x - poc%n0(k)%x
          dir%y = locBox(outDir(k))%y - poc%n0(k)%y

          poc%n1(k)%x = locBox(outDir(k))%x
          poc%n1(k)%y = locBox(outDir(k))%y
          ! write (*, *) outDir, poc%n1(1)%x, poc%n1(1)%y

          xConst = int(poc%n1(k)%x) + int(sign(d1, dir%x))
          yconst = int(poc%n1(k)%y) + int(sign(d1, dir%y))

          x0 = poc%n0(k)%x
          y0 = poc%n0(k)%y
          x1 = poc%n1(k)%x
          y1 = poc%n1(k)%y

          a1 = y1 - y0
          b1 = -(x1 - x0)
          c1 = (y1 - y0)*x0 - (x1 - x0)*y0

          a2 = 1.0
          b2 = 0.0
          c2 = xConst

          determinant = a1*b2 - a2*b1

          if (abs(determinant) < 0.0001) then
            x_xConst = xConst
            y_xConst = 99999
          else
            x_xConst = (c1*b2 - c2*b1)/determinant
            y_xConst = (a1*c2 - a2*c1)/determinant
          end if

          a2 = 0.0
          b2 = 1.0
          c2 = yConst

          determinant = a1*b2 - a2*b1

          if (abs(determinant) < 0.0001) then
            x_yConst = 99999
            y_yConst = yConst
          else
            x_yConst = (c1*b2 - c2*b1)/determinant
            y_yConst = (a1*c2 - a2*c1)/determinant
          end if

          dist_xConst = sqrt((x_xConst - x0)**d2 + (y_xConst - y0)**d2)
          dist_yConst = sqrt((x_yConst - x0)**d2 + (y_yConst - y0)**d2)

          if (dist_xConst <= dist_yConst) then
            poc%n2(k)%x = x_xConst
            poc%n2(k)%y = y_xConst
            poc%n2(k)%b(1)%x = x_xConst
            poc%n2(k)%b(1)%y = floor(y_xConst)
            poc%n2(k)%b(2)%x = x_xConst
            poc%n2(k)%b(2)%y = ceiling(y_xConst)
          else
            poc%n2(k)%x = x_yConst
            poc%n2(k)%y = y_yConst
            poc%n2(k)%b(1)%x = floor(x_yConst)
            poc%n2(k)%b(1)%y = y_yConst
            poc%n2(k)%b(2)%x = ceiling(x_yConst)
            poc%n2(k)%b(2)%y = y_yConst
          end if
        end do
      end associate
    end do

  end subroutine createDataStructOnCircle

  pure subroutine createDSonBar(noOfPtOnBar, ptOnBar)
    integer, intent(in)::noOfPtOnBar
    type(custom_t), allocatable, dimension(:), intent(out) :: ptOnBar
    type(doublet_t):: startPtOnBar, dir
    type(doubletInt_t)::locBox(4)
    double precision::dirDotUnitVec(4)
    double precision:: x0, x1, y0, y1, a1, b1, c1, a2, b2, c2, determinant
    double precision:: x_xConst, y_xConst, x_yConst, y_yConst
    double precision:: x_xConst2, y_xConst2, x_yConst2, y_yConst2
    double precision:: dist_xConst, dist_yConst, dist_xConst2, dist_yConst2
    integer::i, a, k, outDir(3), xConst, yConst, xConst2, yConst2!, itmp(1)
    double precision:: ds, arcLen, leftOverDist, min2ndDist

    allocate (ptOnBar(noOfPtOnBar))

    startPtOnBar%x = xc_ + haf*dia_
    startPtOnBar%y = yc_ - haf*barH_
    ds = (d2*barL_ + d2*barH_)/noOfPtOnBar

    ! xFlag = 1
    ! yFlag = 0
    arcLen = d0

    do i = 0, noOfPtOnBar - 1
      associate (pob => ptOnBar(i + 1))
        arcLen = (i + haf)*ds

        if (arcLen .lt. barL_) then
          pob%n0(:)%x = startPtOnBar%x + arcLen
          pob%n0(:)%y = startPtOnBar%y
          pob%uv%x = 0.0d0
          pob%uv%y = -1.0d0
        elseif (arcLen .lt. (barL_ + barH_)) then
          pob%n0(:)%x = startPtOnBar%x + barL_
          pob%n0(:)%y = startPtOnBar%y + arcLen - barL_
          pob%uv%x = 1.0d0
          pob%uv%y = 0.0d0
        elseif (arcLen .lt. (d2*barL_ + barH_)) then
          pob%n0(:)%x = startPtOnBar%x + (2*barL_ + barH_ - arcLen)
          pob%n0(:)%y = startPtOnBar%y + barH_
          pob%uv%x = 0.0d0
          pob%uv%y = 1.0d0
        else
          pob%n0(:)%x = startPtOnBar%x
          pob%n0(:)%y = startPtOnBar%y + arcLen - (d2*barL_ + barH_)
          pob%uv%x = -1.0d0
          pob%uv%y = 0.0d0
        end if
        ! write (*, '(I4,5(F10.4))') (i + 1), arcLen, pob%n0(1)%x, pob%n0(1)%y, pob%uv%x, pob%uv%y

        locBox(1)%x = ceiling(pob%n0(1)%x)
        locBox(2)%x = ceiling(pob%n0(1)%x)
        locBox(3)%x = floor(pob%n0(1)%x)
        locBox(4)%x = floor(pob%n0(1)%x)

        locBox(1)%y = floor(pob%n0(1)%y)
        locBox(2)%y = ceiling(pob%n0(1)%y)
        locBox(3)%y = ceiling(pob%n0(1)%y)
        locBox(4)%y = floor(pob%n0(1)%y)

        do a = 1, 4
          dir%x = locBox(a)%x - pob%n0(1)%x
          dir%y = locBox(a)%y - pob%n0(1)%y
          dirDotUnitVec(a) = (dir%x*pob%uv%x + dir%y*pob%uv%y) &
                             /(sqrt(dir%x**d2 + dir%y**d2))
        end do

        pob%noOfRD = 0
        do a = 1, 4
          ! if (dirDotUnitVec(a) .eq. maxval(dirDotUnitVec)) then
          if (dirDotUnitVec(a) .gt. d0) then
            pob%noOfRD = pob%noOfRD + 1
            outDir(pob%noOfRD) = a
            ! exit
          end if
        end do

        ! itmp = maxloc(dirDotUnitVec)
        ! outDir = itmp(1)
        do k = 1, pob%noOfRD
          dir%x = locBox(outDir(k))%x - pob%n0(k)%x
          dir%y = locBox(outDir(k))%y - pob%n0(k)%y

          pob%n1(k)%x = locBox(outDir(k))%x
          pob%n1(k)%y = locBox(outDir(k))%y
          ! write (*, *) outDir, pob%n1(1)%x, pob%n1(1)%y

          xConst = int(pob%n1(k)%x) + int(sign(d1, dir%x))
          yconst = int(pob%n1(k)%y) + int(sign(d1, dir%y))

          xConst2 = int(pob%n1(k)%x) + 2*int(sign(d1, dir%x))
          yconst2 = int(pob%n1(k)%y) + 2*int(sign(d1, dir%y))

          x0 = pob%n0(k)%x
          y0 = pob%n0(k)%y
          x1 = pob%n1(k)%x
          y1 = pob%n1(k)%y

          a1 = y1 - y0
          b1 = -(x1 - x0)
          c1 = (y1 - y0)*x0 - (x1 - x0)*y0

          a2 = 1.0
          b2 = 0.0
          c2 = xConst

          determinant = a1*b2 - a2*b1

          if (abs(determinant) < 0.0001) then
            x_xConst = xConst
            y_xConst = 99999
          else
            x_xConst = (c1*b2 - c2*b1)/determinant
            y_xConst = (a1*c2 - a2*c1)/determinant
          end if

          a2 = 0.0
          b2 = 1.0
          c2 = yConst

          determinant = a1*b2 - a2*b1

          if (abs(determinant) < 0.0001) then
            x_yConst = 99999
            y_yConst = yConst
          else
            x_yConst = (c1*b2 - c2*b1)/determinant
            y_yConst = (a1*c2 - a2*c1)/determinant
          end if

          !---------------------
          a1 = y1 - y0
          b1 = -(x1 - x0)
          c1 = (y1 - y0)*x0 - (x1 - x0)*y0

          a2 = 1.0
          b2 = 0.0
          c2 = xConst2

          determinant = a1*b2 - a2*b1

          if (abs(determinant) < 0.0001) then
            x_xConst2 = xConst2
            y_xConst2 = 99999
          else
            x_xConst2 = (c1*b2 - c2*b1)/determinant
            y_xConst2 = (a1*c2 - a2*c1)/determinant
          end if

          a2 = 0.0
          b2 = 1.0
          c2 = yConst2

          determinant = a1*b2 - a2*b1

          if (abs(determinant) < 0.0001) then
            x_yConst2 = 99999
            y_yConst2 = yConst
          else
            x_yConst2 = (c1*b2 - c2*b1)/determinant
            y_yConst2 = (a1*c2 - a2*c1)/determinant
          end if

          dist_xConst = sqrt((x_xConst - x0)**d2 + (y_xConst - y0)**d2)
          dist_yConst = sqrt((x_yConst - x0)**d2 + (y_yConst - y0)**d2)
          dist_xConst2 = sqrt((x_xConst2 - x0)**d2 + (y_xConst2 - y0)**d2)
          dist_yConst2 = sqrt((x_yConst2 - x0)**d2 + (y_yConst2 - y0)**d2)
          !---------------------
          if (dist_xConst <= dist_yConst) then
            leftOverDist = dist_yConst
            pob%n2(k)%x = x_xConst
            pob%n2(k)%y = y_xConst
            pob%n2(k)%b(1)%x = x_xConst
            pob%n2(k)%b(1)%y = floor(y_xConst)
            pob%n2(k)%b(2)%x = x_xConst
            pob%n2(k)%b(2)%y = ceiling(y_xConst)
          else
            leftOverDist = dist_xConst
            pob%n2(k)%x = x_yConst
            pob%n2(k)%y = y_yConst
            pob%n2(k)%b(1)%x = floor(x_yConst)
            pob%n2(k)%b(1)%y = y_yConst
            pob%n2(k)%b(2)%x = ceiling(x_yConst)
            pob%n2(k)%b(2)%y = y_yConst
          end if

          min2ndDist = min(leftOverDist, dist_xConst2, dist_yConst2)
          if (min2ndDist .eq. dist_xConst) then
            pob%n3(k)%x = x_xConst
            pob%n3(k)%y = y_xConst
            pob%n3(k)%b(1)%x = x_xConst
            pob%n3(k)%b(1)%y = floor(y_xConst)
            pob%n3(k)%b(2)%x = x_xConst
            pob%n3(k)%b(2)%y = ceiling(y_xConst)
          elseif (min2ndDist .eq. dist_yConst) then
            pob%n3(k)%x = x_yConst
            pob%n3(k)%y = y_yConst
            pob%n3(k)%b(1)%x = floor(x_yConst)
            pob%n3(k)%b(1)%y = y_yConst
            pob%n3(k)%b(2)%x = ceiling(x_yConst)
            pob%n3(k)%b(2)%y = y_yConst
          elseif (min2ndDist .eq. dist_xConst2) then
            pob%n3(k)%x = x_xConst2
            pob%n3(k)%y = y_xConst2
            pob%n3(k)%b(1)%x = x_xConst2
            pob%n3(k)%b(1)%y = floor(y_xConst2)
            pob%n3(k)%b(2)%x = x_xConst2
            pob%n3(k)%b(2)%y = ceiling(y_xConst2)
          elseif (min2ndDist .eq. dist_yConst2) then
            pob%n3(k)%x = x_yConst2
            pob%n3(k)%y = y_yConst2
            pob%n3(k)%b(1)%x = floor(x_yConst2)
            pob%n3(k)%b(1)%y = y_yConst2
            pob%n3(k)%b(2)%x = ceiling(x_yConst2)
            pob%n3(k)%b(2)%y = y_yConst2
          end if

        end do

      end associate
    end do
    ! stop
  end subroutine createDSonBar

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

  pure subroutine linearExtInt(pt1, pt2, pt)
    ! double precision, intent(in) :: x, y
    type(triplet_t), intent(in) :: pt1, pt2
    type(triplet_t), intent(inout) :: pt
    double precision:: x0x1, x2x1

    if (pt%x /= pt1%x) then
      x0x1 = pt%x - pt1%x
      x2x1 = pt2%x - pt1%x
    elseif (pt%y /= pt1%y) then
      x0x1 = pt%y - pt1%y
      x2x1 = pt2%y - pt1%y
    else
      x0x1 = 0.0d0
      x2x1 = 1.0d0
    end if

    pt%z = pt1%z + x0x1*(pt2%z - pt1%z)/x2x1

  end subroutine linearExtInt

  pure subroutine quadExtInt(pt1, pt2, pt3, ptAsk)
    ! double precision, intent(in) :: x, y
    type(triplet_t), intent(in) :: pt1, pt2, pt3
    type(triplet_t), intent(inout) :: ptAsk
    type(triplet_t) :: pt(0:2)
    double precision:: xAsk, x(0:2), sum, prod(0:2)
    integer:: i, j, m

    pt(0) = pt1
    pt(1) = pt2
    pt(2) = pt3

    xAsk = sqrt(ptAsk%x**d2 + ptAsk%y**d2)
    do i = 0, 2
      x(i) = sqrt(pt(i)%x**d2 + pt(i)%y**d2)
    end do

    sum = d0
    do j = 0, 2
      prod(j) = d1
      do m = 0, 2
        if (m == j) cycle
        prod(j) = prod(j)*(xAsk - x(m))/(x(j) - x(m))
      end do
      sum = sum + pt(j)%z*prod(j)
    end do

    ptAsk%z = sum

  end subroutine quadExtInt

  ! pure subroutine bilinearInterp(pt, box, err)
  !   type(triplet_t), intent(inout) :: pt
  !   type(box_t), intent(in) :: box(4)
  !   integer, optional, intent(out) :: err
  !   double precision::locX, locY, coeff(4)
  !   integer::i

  !   locX = pt%x - box(4)%x
  !   locY = pt%y - box(4)%y

  !   if (present(err)) then
  !     err = 0
  !     if (locX .gt. 1.0005 .or. locY .gt. 1.0005) then
  !       ! write (*, *) "Bilinear Interpolation: Ask-point outside box. Aborting"
  !       ! stop
  !       err = 1
  !     end if
  !   end if

  !   coeff = rectShape(locX, locY)

  !   pt%z = 0
  !   do i = 1, 4
  !     pt%z = pt%z + coeff(i)*box(i)%z
  !   end do

  ! end subroutine bilinearInterp

  ! pure function rectShape(x, y) result(psi)
  !   implicit none

  !   double precision, intent(in) :: x, y
  !   double precision::psi(4)

  !   psi = [x*(1 - y), x*y, (1 - x)*y, (1 - x)*(1 - y)]

  ! end function rectShape

  pure function trapIntegrate(force) result(totalForce)
    implicit none
    double precision, dimension(noOfPtOnCircle), intent(in)::force
    double precision::totalForce, dx
    integer::i, ip1

    dx = pi*dia_/noOfPtOnCircle

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

end program cyl
