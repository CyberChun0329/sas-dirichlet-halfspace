module three_dimensional_analysis
   !! 三维工程积分模型的现代 Fortran 实现。
   !!
   !! 本模块把原来散落在固定格式 .for 文件中的过程集中到一个模块，目的有三点：
   !! 1. 所有过程都具有显式接口，编译器可以检查参数个数、类型和 intent；
   !! 2. 统一使用 real64，避免旧代码依赖隐式类型和编译器默认精度；
   !! 3. 同时保留原始全域 20 点积分和“按观测点分区”的改进积分，便于复现实验并
   !!    单独评估近奇异积分的误差。
   use, intrinsic :: iso_fortran_env, only : real64
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_quiet_nan, &
                                             ieee_value
   implicit none
   private

   integer, parameter, public :: dp = real64
   integer, parameter, public :: quadrature_global = 0
   integer, parameter, public :: quadrature_split_at_observer = 1

   real(dp), parameter :: pi = acos(-1.0_dp)
   real(dp), parameter :: exp_underflow_limit = sqrt(-log(tiny(1.0_dp)))
   integer, parameter :: adaptive_limit = 1000
   integer, parameter :: kernel_steady = 0
   integer, parameter :: kernel_transient = 1
   integer, parameter :: kernel_loss = 2

   ! 20 点 Gauss-Legendre 节点与权重。数组只读并由所有二维积分共用。
   real(dp), parameter :: gl20_nodes(20) = [ &
      -0.993128599185094924786_dp, -0.963971927277913791268_dp, &
      -0.912234428251325905868_dp, -0.839116971822218823395_dp, &
      -0.746331906460150792614_dp, -0.636053680726515025453_dp, &
      -0.510867001950827098004_dp, -0.373706088715419560673_dp, &
      -0.227785851141645078080_dp, -0.076526521133497333755_dp, &
       0.076526521133497333755_dp,  0.227785851141645078080_dp, &
       0.373706088715419560673_dp,  0.510867001950827098004_dp, &
       0.636053680726515025453_dp,  0.746331906460150792614_dp, &
       0.839116971822218823395_dp,  0.912234428251325905868_dp, &
       0.963971927277913791268_dp,  0.993128599185094924786_dp ]

   real(dp), parameter :: gl20_weights(20) = [ &
      0.017614007139152118312_dp, 0.040601429800386941331_dp, &
      0.062672048334109063570_dp, 0.083276741576704748725_dp, &
      0.101930119817240435037_dp, 0.118194531961518417312_dp, &
      0.131688638449176626898_dp, 0.142096109318382051329_dp, &
      0.149172986472603746788_dp, 0.152753387130725850698_dp, &
      0.152753387130725850698_dp, 0.149172986472603746788_dp, &
      0.142096109318382051329_dp, 0.131688638449176626898_dp, &
      0.118194531961518417312_dp, 0.101930119817240435037_dp, &
      0.083276741576704748725_dp, 0.062672048334109063570_dp, &
      0.040601429800386941331_dp, 0.017614007139152118312_dp ]

   abstract interface
      function scalar_integrand(argument) result(value)
         import :: dp
         real(dp), intent(in) :: argument
         real(dp) :: value
      end function scalar_integrand
   end interface

   public :: semianalytic_loss
   public :: direct_loss
   public :: steady_integral_direct
   public :: transient_integral_direct
   public :: loss_integral_direct
   public :: steady_integral_semianalytic
   public :: find_time_to_fraction
   public :: lower_incomplete_gamma_indexed
   public :: relative_error
   public :: integrate_gk15

contains

   function semianalytic_loss(m, n, x_bounds, y_bounds, x, y, z, time_limit, &
                              diffusivity, status) result(value)
      !! 利用积分等式计算三维温度亏损的半解析解，对应旧程序 THETA3/F3。
      !!
      !! 原 F3 在函数内部强制把 M、N 改为 0，实际只实现均匀面源情形。现代版
      !! 明确拒绝非零阶次，以免调用者误以为高阶输入已经被使用。
      integer, intent(in) :: m, n
      real(dp), intent(in) :: x_bounds(2), y_bounds(2)
      real(dp), intent(in) :: x, y, z, time_limit, diffusivity
      integer, intent(out), optional :: status

      real(dp) :: value
      real(dp) :: split_point, effective_split, upper_limit
      real(dp) :: first_part, second_part, first_error, second_error
      integer :: first_status, second_status, local_status

      local_status = validate_common_inputs(m, n, x_bounds, y_bounds, x, y, z, &
                                            time_limit, diffusivity)
      if (local_status == 0 .and. (m /= 0 .or. n /= 0)) local_status = 3

      if (local_status /= 0) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
         if (present(status)) status = local_status
         return
      end if

      split_point = 1.0_dp / sqrt(4.0_dp * diffusivity * time_limit)

      ! exp(-(z*u)^2) 在 z*u=6 时已小于 2.4e-16。保留旧程序至少取 100
      ! 的截断规则，以便结果可比；更远尾部对 real64 精度已没有实质贡献。
      upper_limit = max(6.0_dp / z, 100.0_dp)
      effective_split = min(split_point, upper_limit)

      call integrate_gk15(kernel, 0.0_dp, effective_split, 1.0e-12_dp, &
                          1.0e-10_dp, first_part, first_error, first_status)
      if (split_point < upper_limit) then
         call integrate_gk15(kernel, split_point, upper_limit, 1.0e-12_dp, &
                             1.0e-10_dp, second_part, second_error, second_status)
      else
         ! 极短时间下分段点可能位于截断上限之外；此时整个有效区间都用低段公式。
         second_part = 0.0_dp
         second_error = 0.0_dp
         second_status = 0
      end if

      value = first_part + second_part
      local_status = max(first_status, second_status)
      if (present(status)) status = local_status

   contains

      function kernel(u) result(kernel_value)
         !! 经过积分等式降维后的一维被积函数。
         real(dp), intent(in) :: u
         real(dp) :: kernel_value
         real(dp) :: corner_sum, exponential_factor
         integer :: i, j

         corner_sum = 0.0_dp
         do j = 1, 2
            do i = 1, 2
               corner_sum = corner_sum + (-1.0_dp)**(i + j) * &
                  signum(x_bounds(i) - x) * signum(y_bounds(j) - y) * &
                  erf(abs(x_bounds(i) - x) * u) * &
                  erf(abs(y_bounds(j) - y) * u)
            end do
         end do

         exponential_factor = exp_negative_square(z * u)
         if (u <= split_point) then
            kernel_value = diffusivity * z * time_limit * corner_sum * &
                           exponential_factor / (2.0_dp * sqrt(pi))
         else
            kernel_value = z * corner_sum * exponential_factor / &
                           (8.0_dp * sqrt(pi) * u**2)
         end if
      end function kernel

   end function semianalytic_loss


   function direct_loss(m, n, x_bounds, y_bounds, x, y, z, time_limit, &
                        diffusivity, spatial_mode, use_stable_loss_kernel, &
                        status) &
                        result(value)
      !! 对“稳态解 - 瞬态解”做时间积分，对应旧程序 THETA4/F4。
      !!
      !! use_stable_loss_kernel=.false. 可严格复现旧算法：每个时间积分点分别计算
      !! 稳态与瞬态二维积分。设为 .true. 时直接积分稳定亏损核，是结果等价的优化。
      integer, intent(in) :: m, n
      real(dp), intent(in) :: x_bounds(2), y_bounds(2)
      real(dp), intent(in) :: x, y, z, time_limit, diffusivity
      integer, intent(in), optional :: spatial_mode
      logical, intent(in), optional :: use_stable_loss_kernel
      integer, intent(out), optional :: status

      real(dp) :: value
      real(dp) :: absolute_error
      integer :: mode, local_status
      logical :: use_stable_kernel

      local_status = validate_common_inputs(m, n, x_bounds, y_bounds, x, y, z, &
                                            time_limit, diffusivity)
      mode = quadrature_global
      if (present(spatial_mode)) mode = spatial_mode
      if (mode /= quadrature_global .and. &
          mode /= quadrature_split_at_observer) local_status = 4

      use_stable_kernel = .true.
      if (present(use_stable_loss_kernel)) then
         use_stable_kernel = use_stable_loss_kernel
      end if

      if (local_status /= 0) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
         if (present(status)) status = local_status
         return
      end if

      call integrate_gk15(time_kernel, 0.0_dp, time_limit, 1.0e-12_dp, &
                          1.0e-10_dp, value, absolute_error, local_status)
      if (present(status)) status = local_status

   contains

      function time_kernel(time) result(kernel_value)
         real(dp), intent(in) :: time
         real(dp) :: kernel_value
         real(dp) :: current_steady

         if (use_stable_kernel) then
            ! 优化路径直接积“稳态-瞬态”的稳定亏损核，既省掉稳态重复计算，也
            ! 避免大时间下两个接近的二维积分结果相减。
            kernel_value = loss_integral_direct(m, n, x_bounds, y_bounds, x, y, &
                                                 z, time, diffusivity, mode)
            return
         else
            current_steady = steady_integral_direct(m, n, x_bounds, &
                                                    y_bounds, x, y, z, mode)
         end if
         kernel_value = current_steady - transient_integral_direct( &
                           m, n, x_bounds, y_bounds, x, y, z, time, &
                           diffusivity, mode)
      end function time_kernel

   end function direct_loss


   function steady_integral_direct(m, n, x_bounds, y_bounds, x, y, z, &
                                   spatial_mode, status) result(value)
      !! 用张量积 Gauss-Legendre 公式直接计算稳态二维面积分。
      integer, intent(in) :: m, n
      real(dp), intent(in) :: x_bounds(2), y_bounds(2)
      real(dp), intent(in) :: x, y, z
      integer, intent(in), optional :: spatial_mode
      integer, intent(out), optional :: status
      real(dp) :: value
      integer :: mode, local_status

      mode = quadrature_global
      if (present(spatial_mode)) mode = spatial_mode
      local_status = validate_spatial_inputs(m, n, x_bounds, y_bounds, x, y, &
                                             z, mode)
      if (local_status /= 0) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
      else
         value = spatial_integral(m, n, x_bounds, y_bounds, x, y, z, &
                                  0.0_dp, 1.0_dp, kernel_steady, mode)
      end if
      if (present(status)) status = local_status
   end function steady_integral_direct


   function transient_integral_direct(m, n, x_bounds, y_bounds, x, y, z, &
                                      time, diffusivity, spatial_mode, &
                                      status) result(value)
      !! 用张量积 Gauss-Legendre 公式直接计算给定时刻的瞬态二维面积分。
      integer, intent(in) :: m, n
      real(dp), intent(in) :: x_bounds(2), y_bounds(2)
      real(dp), intent(in) :: x, y, z, time, diffusivity
      integer, intent(in), optional :: spatial_mode
      integer, intent(out), optional :: status
      real(dp) :: value
      integer :: mode, local_status

      mode = quadrature_global
      if (present(spatial_mode)) mode = spatial_mode
      local_status = validate_spatial_inputs(m, n, x_bounds, y_bounds, x, y, &
                                             z, mode)
      ! Fortran 不保证 .or. 短路求值：先排除 NaN/Inf，再做大小比较。
      if (.not. ieee_is_finite(time) .or. &
          .not. ieee_is_finite(diffusivity)) then
         local_status = 2
      else if (time < 0.0_dp .or. diffusivity <= 0.0_dp) then
         local_status = 2
      end if

      if (local_status /= 0) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
      else if (time <= 0.0_dp) then
         ! t -> 0+ 时旧公式的括号整体趋于 0；显式处理可避免 1/sqrt(t)。
         value = 0.0_dp
      else
         value = spatial_integral(m, n, x_bounds, y_bounds, x, y, z, &
                                  time, diffusivity, kernel_transient, mode)
      end if
      if (present(status)) status = local_status
   end function transient_integral_direct


   function loss_integral_direct(m, n, x_bounds, y_bounds, x, y, z, time, &
                                 diffusivity, spatial_mode, status) result(value)
      !! 直接计算“稳态-瞬态”二维亏损积分，避免先算两个近似相等的大数再相减。
      integer, intent(in) :: m, n
      real(dp), intent(in) :: x_bounds(2), y_bounds(2)
      real(dp), intent(in) :: x, y, z, time, diffusivity
      integer, intent(in), optional :: spatial_mode
      integer, intent(out), optional :: status
      real(dp) :: value
      integer :: mode, local_status

      mode = quadrature_global
      if (present(spatial_mode)) mode = spatial_mode
      local_status = validate_spatial_inputs(m, n, x_bounds, y_bounds, x, y, &
                                             z, mode)
      if (.not. ieee_is_finite(time) .or. &
          .not. ieee_is_finite(diffusivity)) then
         local_status = 2
      else if (time < 0.0_dp .or. diffusivity <= 0.0_dp) then
         local_status = 2
      end if

      if (local_status /= 0) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
      else if (time <= 0.0_dp) then
         value = steady_integral_direct(m, n, x_bounds, y_bounds, x, y, z, mode)
      else
         value = spatial_integral(m, n, x_bounds, y_bounds, x, y, z, time, &
                                  diffusivity, kernel_loss, mode)
      end if
      if (present(status)) status = local_status
   end function loss_integral_direct


   function spatial_integral(m, n, x_bounds, y_bounds, x, y, z, time, &
                             diffusivity, kernel_kind, mode) result(value)
      !! 在一个或四个子矩形上进行 20x20 张量积积分。
      !!
      !! 当观测点位于积分域内部时，分区模式以 x、y 为分界形成四个子矩形。
      !! Gauss 节点会自然聚集到各子区端点，因而能更好地解析 z 很小时位于
      !! (x,y) 附近的尖峰；全域模式则与两份旧程序完全一致。
      integer, intent(in) :: m, n, kernel_kind, mode
      real(dp), intent(in) :: x_bounds(2), y_bounds(2)
      real(dp), intent(in) :: x, y, z, time, diffusivity
      real(dp) :: value

      real(dp) :: x_cuts(3), y_cuts(3)
      real(dp) :: x_mid, x_half, y_mid, y_half
      real(dp) :: xi, eta, radius, scaled_radius, kernel
      integer :: nx, ny, ix, iy, i, j

      call make_cuts(x_bounds, x, mode, x_cuts, nx)
      call make_cuts(y_bounds, y, mode, y_cuts, ny)

      value = 0.0_dp
      do iy = 1, ny
         y_mid = 0.5_dp * (y_cuts(iy) + y_cuts(iy + 1))
         y_half = 0.5_dp * (y_cuts(iy + 1) - y_cuts(iy))
         do ix = 1, nx
            x_mid = 0.5_dp * (x_cuts(ix) + x_cuts(ix + 1))
            x_half = 0.5_dp * (x_cuts(ix + 1) - x_cuts(ix))

            do j = 1, size(gl20_nodes)
               eta = y_mid + y_half * gl20_nodes(j)
               do i = 1, size(gl20_nodes)
                  xi = x_mid + x_half * gl20_nodes(i)
                  radius = sqrt((xi - x)**2 + (eta - y)**2 + z**2)

                  select case (kernel_kind)
                  case (kernel_transient)
                     scaled_radius = radius / sqrt(4.0_dp * diffusivity * time)
                     ! erfc(s) 比 1-erf(s) 在 s 较大时稳定，第二项与旧公式等价。
                     ! 超过阈值后两项均低于 real64 正规数范围，直接置零可避免
                     ! 大量无意义的指数下溢及 IEEE_UNDERFLOW_FLAG。
                     if (scaled_radius >= exp_underflow_limit) then
                        kernel = 0.0_dp
                     else
                        kernel = erfc(scaled_radius) + 2.0_dp * scaled_radius * &
                                 exp_negative_square(scaled_radius) / sqrt(pi)
                     end if
                  case (kernel_loss)
                     scaled_radius = radius / sqrt(4.0_dp * diffusivity * time)
                     if (scaled_radius >= exp_underflow_limit) then
                        kernel = 1.0_dp
                     else
                        kernel = stable_loss_factor(scaled_radius)
                     end if
                  case default
                     kernel = 1.0_dp
                  end select

                  value = value + x_half * y_half * gl20_weights(i) * &
                     gl20_weights(j) * z * xi**m * eta**n * kernel / &
                     (2.0_dp * pi * radius**3)
               end do
            end do
         end do
      end do
   end function spatial_integral


   subroutine make_cuts(bounds, observer, mode, cuts, part_count)
      real(dp), intent(in) :: bounds(2), observer
      integer, intent(in) :: mode
      real(dp), intent(out) :: cuts(3)
      integer, intent(out) :: part_count
      real(dp) :: scale, tolerance

      cuts = 0.0_dp
      scale = max(1.0_dp, abs(bounds(1)), abs(bounds(2)))
      tolerance = 32.0_dp * epsilon(1.0_dp) * scale
      if (mode == quadrature_split_at_observer .and. &
          observer > bounds(1) + tolerance .and. &
          observer < bounds(2) - tolerance) then
         cuts(1:3) = [bounds(1), observer, bounds(2)]
         part_count = 2
      else
         cuts(1:2) = bounds
         part_count = 1
      end if
   end subroutine make_cuts


   function steady_integral_semianalytic(m, n, x_bounds, y_bounds, x, y, z, &
                                         status) result(value)
      !! 旧程序 TSS/AJ 系列的现代化版本，用于稳态积分的独立交叉校验。
      integer, intent(in) :: m, n
      real(dp), intent(in) :: x_bounds(2), y_bounds(2), x, y, z
      integer, intent(out), optional :: status
      real(dp) :: value
      integer :: k, l, local_status

      local_status = validate_spatial_inputs(m, n, x_bounds, y_bounds, x, y, &
                                             z, quadrature_global)
      if (local_status /= 0) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
      else
         value = 0.0_dp
         do l = 0, n
            do k = 0, m
               value = value + binomial_coefficient(m, k) * &
                  binomial_coefficient(n, l) * x**(m - k) * y**(n - l) * &
                  auxiliary_j(k, l, x_bounds, y_bounds, x, y, z)
            end do
         end do
      end if
      if (present(status)) status = local_status
   end function steady_integral_semianalytic


   function binomial_coefficient(n, k) result(value)
      integer, intent(in) :: n, k
      real(dp) :: value
      integer :: i, k_min

      if (k < 0 .or. k > n) then
         value = 0.0_dp
         return
      end if
      k_min = min(k, n - k)
      value = 1.0_dp
      do i = 1, k_min
         value = value * real(n - k_min + i, dp) / real(i, dp)
      end do
   end function binomial_coefficient


   function auxiliary_j(m, n, x_bounds, y_bounds, x, y, z) result(value)
      integer, intent(in) :: m, n
      real(dp), intent(in) :: x_bounds(2), y_bounds(2), x, y, z
      real(dp) :: value
      integer :: i, j, jj, first_order

      if (n == 0) then
         value = auxiliary_j0(m, x_bounds, y_bounds, x, y, z)
         return
      else if (n == 1) then
         ! 旧代码误写为 M==1，导致许多混合阶次走错递推起点。
         value = auxiliary_j1(m, x_bounds, y_bounds, x, y, z)
         return
      else if (mod(n, 2) == 0) then
         first_order = 2
         value = auxiliary_j0(m, x_bounds, y_bounds, x, y, z)
      else
         first_order = 3
         value = auxiliary_j1(m, x_bounds, y_bounds, x, y, z)
      end if

      do jj = first_order, n, 2
         value = -real(jj - 1, dp) * z**2 * value / real(m + jj - 1, dp)
         do i = 1, 2
            ! 旧式在此多除以 (M+1)，会系统性破坏高阶递推。
            value = value + real(jj - 1, dp) * z * (-1.0_dp)**i * &
               (x_bounds(i) - x)**(m + 1) * &
               auxiliary_h(jj - 2, y_bounds, x, y, z, x_bounds(i)) / &
               (2.0_dp * pi * real(m + jj - 1, dp))
         end do
         do j = 1, 2
            value = value - real(m, dp) * z * (-1.0_dp)**j * &
               (y_bounds(j) - y)**(jj - 1) * &
               auxiliary_g(m, x_bounds, x, y, z, y_bounds(j)) / &
               (2.0_dp * pi * real(m + jj - 1, dp))
         end do
      end do
   end function auxiliary_j


   function auxiliary_j0(m, x_bounds, y_bounds, x, y, z) result(value)
      integer, intent(in) :: m
      real(dp), intent(in) :: x_bounds(2), y_bounds(2), x, y, z
      real(dp) :: value
      real(dp) :: base_even, base_odd, dx, dy, radius
      integer :: i, j, ii, first_order

      base_even = 0.0_dp
      base_odd = 0.0_dp
      do j = 1, 2
         dy = y_bounds(j) - y
         do i = 1, 2
            dx = x_bounds(i) - x
            radius = sqrt(dx**2 + dy**2 + z**2)
            base_even = base_even + (-1.0_dp)**(i + j) * &
               atan(dx * dy / (z * radius)) / (2.0_dp * pi)
            base_odd = base_odd - z * (-1.0_dp)**(i + j) * &
               stable_log_sum(dy, dx**2 + z**2) / (2.0_dp * pi)
         end do
      end do

      if (m == 0) then
         value = base_even
         return
      else if (m == 1) then
         value = base_odd
         return
      else if (mod(m, 2) == 0) then
         first_order = 2
         value = base_even
      else
         first_order = 3
         value = base_odd
      end if

      do ii = first_order, m, 2
         value = -z**2 * value
         do j = 1, 2
            value = value + z * (-1.0_dp)**j * (y_bounds(j) - y) * &
               auxiliary_g(ii - 2, x_bounds, x, y, z, y_bounds(j)) / &
               (2.0_dp * pi)
         end do
      end do
   end function auxiliary_j0


   function auxiliary_j1(m, x_bounds, y_bounds, x, y, z) result(value)
      integer, intent(in) :: m
      real(dp), intent(in) :: x_bounds(2), y_bounds(2), x, y, z
      real(dp) :: value
      integer :: j

      value = 0.0_dp
      do j = 1, 2
         value = value - z * (-1.0_dp)**j * &
            auxiliary_g(m, x_bounds, x, y, z, y_bounds(j)) / (2.0_dp * pi)
      end do
   end function auxiliary_j1


   function auxiliary_h(n, y_bounds, x, y, z, xi) result(value)
      integer, intent(in) :: n
      real(dp), intent(in) :: y_bounds(2), x, y, z, xi
      real(dp) :: value
      real(dp) :: base_even, base_odd, dy, radius, radial_square
      integer :: j, jj, first_order

      radial_square = (xi - x)**2 + z**2
      base_even = 0.0_dp
      base_odd = 0.0_dp
      do j = 1, 2
         dy = y_bounds(j) - y
         radius = sqrt(radial_square + dy**2)
         base_even = base_even + (-1.0_dp)**j * stable_log_sum(dy, radial_square)
         base_odd = base_odd + (-1.0_dp)**j * radius
      end do

      if (n == 0) then
         value = base_even
         return
      else if (n == 1) then
         value = base_odd
         return
      else if (mod(n, 2) == 0) then
         first_order = 2
         value = base_even
      else
         first_order = 3
         value = base_odd
      end if

      do jj = first_order, n, 2
         value = -real(jj - 1, dp) * radial_square * value / real(jj, dp)
         do j = 1, 2
            dy = y_bounds(j) - y
            radius = sqrt(radial_square + dy**2)
            value = value + (-1.0_dp)**j * dy**(jj - 1) * radius / real(jj, dp)
         end do
      end do
   end function auxiliary_h


   function auxiliary_g(m, x_bounds, x, y, z, eta) result(value)
      integer, intent(in) :: m
      real(dp), intent(in) :: x_bounds(2), x, y, z, eta
      real(dp) :: value
      real(dp) :: base_even, base_odd, dx, radius, radial_square
      integer :: i, ii, first_order

      radial_square = (eta - y)**2 + z**2
      base_even = 0.0_dp
      base_odd = 0.0_dp
      do i = 1, 2
         dx = x_bounds(i) - x
         radius = sqrt(dx**2 + radial_square)
         base_even = base_even + (-1.0_dp)**i * stable_log_sum(dx, radial_square)
         base_odd = base_odd + (-1.0_dp)**i * radius
      end do

      if (m == 0) then
         value = base_even
         return
      else if (m == 1) then
         value = base_odd
         return
      else if (mod(m, 2) == 0) then
         first_order = 2
         value = base_even
      else
         first_order = 3
         value = base_odd
      end if

      do ii = first_order, m, 2
         value = -real(ii - 1, dp) * radial_square * value / real(ii, dp)
         do i = 1, 2
            dx = x_bounds(i) - x
            radius = sqrt(dx**2 + radial_square)
            value = value + (-1.0_dp)**i * dx**(ii - 1) * radius / real(ii, dp)
         end do
      end do
   end function auxiliary_g


   pure function stable_log_sum(delta, positive_square) result(value)
      !! 稳定计算 log(sqrt(delta^2+q)+delta)，避免 delta<0 时相消。
      real(dp), intent(in) :: delta, positive_square
      real(dp) :: value, radius

      radius = sqrt(delta**2 + positive_square)
      if (delta >= 0.0_dp) then
         value = log(radius + delta)
      else
         value = log(positive_square) - log(radius - delta)
      end if
   end function stable_log_sum


   subroutine find_time_to_fraction(x_bounds, y_bounds, x, y, z, diffusivity, &
                                    target_fraction, time, status)
      !! 二分求解瞬态解达到给定稳态比例的时间，对应旧程序 EQUA/FK。
      real(dp), intent(in) :: x_bounds(2), y_bounds(2), x, y, z, diffusivity
      real(dp), intent(in) :: target_fraction
      real(dp), intent(out) :: time
      integer, intent(out) :: status

      real(dp) :: lower, upper, midpoint, f_lower, f_upper, f_midpoint
      real(dp) :: target_steady
      integer :: iteration

      status = validate_common_inputs(0, 0, x_bounds, y_bounds, x, y, z, &
                                      1.0_dp, diffusivity)
      if (.not. ieee_is_finite(target_fraction)) then
         status = 2
      else if (target_fraction <= 0.0_dp .or. target_fraction >= 1.0_dp) then
         status = 2
      end if
      if (status /= 0) then
         time = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if

      ! 稳态组合与搜索时间无关，只计算一次；旧程序在每个试探时间重复计算。
      target_steady = steady_integral_direct(0, 0, x_bounds, y_bounds, x, y, &
         z, quadrature_split_at_observer) &
         - steady_integral_direct(2, 0, x_bounds, y_bounds, x, y, z, &
            quadrature_split_at_observer) / (x_bounds(2) - x_bounds(1)) &
         - steady_integral_direct(0, 2, x_bounds, y_bounds, x, y, z, &
            quadrature_split_at_observer) / (y_bounds(2) - y_bounds(1))
      if (.not. ieee_is_finite(target_steady)) then
         status = 3
         time = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if

      ! 从 t=0 起步，不漏掉小于旧固定起点 1e-3 的根；上界按倍数扩张，
      ! 比每次线性增加 0.1 更适合跨越多个时间尺度。
      lower = 0.0_dp
      f_lower = fraction_residual(lower)
      upper = 1.0e-6_dp
      f_upper = fraction_residual(upper)
      if (.not. ieee_is_finite(f_lower) .or. &
          .not. ieee_is_finite(f_upper)) then
         status = 3
         time = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if

      do iteration = 1, 200
         if (brackets_zero(f_lower, f_upper)) exit
         lower = upper
         f_lower = f_upper
         upper = 2.0_dp * upper
         f_upper = fraction_residual(upper)
         if (.not. ieee_is_finite(f_upper)) then
            status = 3
            time = ieee_value(0.0_dp, ieee_quiet_nan)
            return
         end if
      end do
      if (.not. brackets_zero(f_lower, f_upper)) then
         status = 1
         time = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if

      do iteration = 1, 200
         midpoint = 0.5_dp * (lower + upper)
         f_midpoint = fraction_residual(midpoint)
         if (.not. ieee_is_finite(f_midpoint)) then
            status = 3
            time = ieee_value(0.0_dp, ieee_quiet_nan)
            return
         end if
         if (abs(upper - lower) <= 1.0e-6_dp * max(midpoint, 1.0_dp)) exit
         if (brackets_zero(f_lower, f_midpoint)) then
            upper = midpoint
            f_upper = f_midpoint
         else
            lower = midpoint
            f_lower = f_midpoint
         end if
      end do
      time = midpoint

   contains

      function fraction_residual(t) result(value)
         real(dp), intent(in) :: t
         real(dp) :: value, transient_value

         transient_value = transient_integral_direct(0, 0, x_bounds, y_bounds, &
            x, y, z, t, diffusivity, quadrature_split_at_observer) &
            - transient_integral_direct(2, 0, x_bounds, y_bounds, x, y, z, t, &
               diffusivity, quadrature_split_at_observer) / &
               (x_bounds(2) - x_bounds(1)) &
            - transient_integral_direct(0, 2, x_bounds, y_bounds, x, y, z, t, &
               diffusivity, quadrature_split_at_observer) / &
               (y_bounds(2) - y_bounds(1))

         value = transient_value - target_fraction * target_steady
      end function fraction_residual


      pure function brackets_zero(left_value, right_value) result(bracketed)
         !! 通过符号而不是乘积判断根是否被夹住，避免大残差乘法溢出。
         real(dp), intent(in) :: left_value, right_value
         logical :: bracketed

         bracketed = (left_value <= 0.0_dp .and. right_value >= 0.0_dp) .or. &
                     (right_value <= 0.0_dp .and. left_value >= 0.0_dp)
      end function brackets_zero

   end subroutine find_time_to_fraction


   function lower_incomplete_gamma_indexed(order_twice, x, status) result(value)
      !! 计算 gamma(a,x)，其中 a=order_twice/2；兼容旧函数 GAMMA(N1,X)。
      !! 小 x 区域直接使用正项级数，避免旧递推在高阶下发生灾难性相消。
      integer, intent(in) :: order_twice
      real(dp), intent(in) :: x
      integer, intent(out), optional :: status
      real(dp) :: value
      real(dp) :: a, ap, term, series_sum
      real(dp) :: b, c, d, h, an, q, log_complete
      integer :: iteration, local_status

      local_status = 0
      if (order_twice <= 0 .or. .not. ieee_is_finite(x)) then
         local_status = 2
      else if (x < 0.0_dp) then
         local_status = 2
      end if
      if (local_status /= 0) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
         if (present(status)) status = local_status
         return
      else if (x <= 0.0_dp) then
         value = 0.0_dp
         if (present(status)) status = 0
         return
      end if

      a = 0.5_dp * real(order_twice, dp)
      if (x < a + 1.0_dp) then
         ! gamma(a,x)=exp(-x+a*log(x))*[1/a + x/(a(a+1)) + ...]
         ap = a
         term = 1.0_dp / a
         series_sum = term
         do iteration = 1, 100000
            ap = ap + 1.0_dp
            term = term * x / ap
            series_sum = series_sum + term
            if (abs(term) <= abs(series_sum) * 8.0_dp * epsilon(1.0_dp)) exit
         end do
         if (iteration > 100000) local_status = 1
         value = series_sum * exp(-x + a * log(x))
      else
         ! Lentz 连分式先求正规化上不完整伽马 Q，再由 gamma(a)*(1-Q) 得下部。
         b = x + 1.0_dp - a
         c = 1.0_dp / tiny(1.0_dp)
         d = 1.0_dp / b
         h = d
         do iteration = 1, 100000
            an = -real(iteration, dp) * (real(iteration, dp) - a)
            b = b + 2.0_dp
            d = an * d + b
            if (abs(d) < tiny(1.0_dp)) d = tiny(1.0_dp)
            c = b + an / c
            if (abs(c) < tiny(1.0_dp)) c = tiny(1.0_dp)
            d = 1.0_dp / d
            term = d * c
            h = h * term
            if (abs(term - 1.0_dp) <= 8.0_dp * epsilon(1.0_dp)) exit
         end do
         if (iteration > 100000) local_status = 1
         log_complete = log_gamma(a)
         q = exp(-x + a * log(x) - log_complete) * h
         value = exp(log_complete) * max(0.0_dp, 1.0_dp - q)
      end if

      if (present(status)) status = local_status
   end function lower_incomplete_gamma_indexed


   pure function relative_error(approximation, reference) result(value)
      real(dp), intent(in) :: approximation, reference
      real(dp) :: value

      if (.not. ieee_is_finite(approximation) .or. &
          .not. ieee_is_finite(reference)) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
      else if (abs(reference) > tiny(1.0_dp)) then
         value = abs(approximation - reference) / abs(reference)
      else
         value = abs(approximation - reference)
      end if
   end function relative_error


   subroutine integrate_gk15(function_to_integrate, lower_bound, upper_bound, &
                             absolute_tolerance, relative_tolerance, result, &
                             absolute_error, status)
      !! 自适应 15 点 Kronrod / 7 点 Gauss 积分。
      !!
      !! 旧代码把 Gauss 节点错接到 XGK(1,3,5)。正确的 7 点 Gauss 子集是
      !! XGK(2,4,6) 加中心点；错误索引会夸大误差估计并造成大量无效细分。
      procedure(scalar_integrand) :: function_to_integrate
      real(dp), intent(in) :: lower_bound, upper_bound
      real(dp), intent(in) :: absolute_tolerance, relative_tolerance
      real(dp), intent(out) :: result, absolute_error
      integer, intent(out) :: status

      real(dp) :: left_list(adaptive_limit), right_list(adaptive_limit)
      real(dp) :: result_list(adaptive_limit), error_list(adaptive_limit)
      real(dp) :: area, error_sum, tolerance
      real(dp) :: left, middle, right, left_area, right_area
      real(dp) :: left_error, right_error, unused_resabs, unused_resasc
      integer :: interval_count, worst, rule_status

      status = 0
      if (.not. ieee_is_finite(lower_bound) .or. &
          .not. ieee_is_finite(upper_bound) .or. &
          .not. ieee_is_finite(absolute_tolerance) .or. &
          .not. ieee_is_finite(relative_tolerance)) then
         status = 2
         result = ieee_value(0.0_dp, ieee_quiet_nan)
         absolute_error = result
         return
      end if
      if (upper_bound < lower_bound .or. absolute_tolerance < 0.0_dp .or. &
          relative_tolerance < 0.0_dp) then
         status = 2
         result = ieee_value(0.0_dp, ieee_quiet_nan)
         absolute_error = result
         return
      else if (upper_bound <= lower_bound) then
         result = 0.0_dp
         absolute_error = 0.0_dp
         return
      end if

      call kronrod15(function_to_integrate, lower_bound, upper_bound, area, &
                     absolute_error, unused_resabs, unused_resasc, rule_status)
      if (rule_status /= 0 .or. .not. ieee_is_finite(area) .or. &
          .not. ieee_is_finite(absolute_error)) then
         status = 3
         result = ieee_value(0.0_dp, ieee_quiet_nan)
         absolute_error = result
         return
      end if
      interval_count = 1
      left_list(1) = lower_bound
      right_list(1) = upper_bound
      result_list(1) = area
      error_list(1) = absolute_error
      error_sum = absolute_error
      tolerance = max(absolute_tolerance, relative_tolerance * abs(area))

      do while (error_sum > tolerance)
         if (interval_count >= adaptive_limit) then
            status = 1
            exit
         end if

         worst = maxloc(error_list(1:interval_count), dim=1)
         left = left_list(worst)
         right = right_list(worst)
         middle = 0.5_dp * (left + right)

         call kronrod15(function_to_integrate, left, middle, left_area, &
                        left_error, unused_resabs, unused_resasc, rule_status)
         if (rule_status /= 0) then
            status = 3
            result = ieee_value(0.0_dp, ieee_quiet_nan)
            absolute_error = result
            return
         end if
         call kronrod15(function_to_integrate, middle, right, right_area, &
                        right_error, unused_resabs, unused_resasc, rule_status)
         if (rule_status /= 0 .or. .not. ieee_is_finite(left_area) .or. &
             .not. ieee_is_finite(right_area) .or. &
             .not. ieee_is_finite(left_error) .or. &
             .not. ieee_is_finite(right_error)) then
            status = 3
            result = ieee_value(0.0_dp, ieee_quiet_nan)
            absolute_error = result
            return
         end if

         area = area - result_list(worst) + left_area + right_area
         error_sum = max(0.0_dp, error_sum - error_list(worst) + &
                         left_error + right_error)

         left_list(worst) = left
         right_list(worst) = middle
         result_list(worst) = left_area
         error_list(worst) = left_error

         interval_count = interval_count + 1
         left_list(interval_count) = middle
         right_list(interval_count) = right
         result_list(interval_count) = right_area
         error_list(interval_count) = right_error

         tolerance = max(absolute_tolerance, relative_tolerance * abs(area))
      end do

      result = area
      absolute_error = error_sum
   end subroutine integrate_gk15


   subroutine kronrod15(function_to_integrate, lower_bound, upper_bound, result, &
                        absolute_error, result_abs, result_asc, status)
      procedure(scalar_integrand) :: function_to_integrate
      real(dp), intent(in) :: lower_bound, upper_bound
      real(dp), intent(out) :: result, absolute_error, result_abs, result_asc
      integer, intent(out) :: status

      real(dp), parameter :: xgk(8) = [ &
         0.991455371120812639207_dp, 0.949107912342758524526_dp, &
         0.864864423359769072790_dp, 0.741531185599394439864_dp, &
         0.586087235467691130294_dp, 0.405845151377397166907_dp, &
         0.207784955007898467601_dp, 0.0_dp ]
      real(dp), parameter :: wgk(8) = [ &
         0.022935322010529224964_dp, 0.063092092629978553291_dp, &
         0.104790010322250183840_dp, 0.140653259715525918745_dp, &
         0.169004726639267902827_dp, 0.190350578064785409913_dp, &
         0.204432940075298892414_dp, 0.209482141084727828013_dp ]
      real(dp), parameter :: wg(4) = [ &
         0.129484966168869693271_dp, 0.279705391489276667902_dp, &
         0.381830050505118944950_dp, 0.417959183673469387755_dp ]

      real(dp) :: center, half_length, center_value, result_gauss, result_kronrod
      real(dp) :: f_left(7), f_right(7), pair_sum, mean_value, raw_error
      integer :: j

      status = 0
      center = 0.5_dp * (lower_bound + upper_bound)
      half_length = 0.5_dp * (upper_bound - lower_bound)
      center_value = function_to_integrate(center)
      if (.not. ieee_is_finite(center_value)) then
         status = 3
         result = ieee_value(0.0_dp, ieee_quiet_nan)
         absolute_error = result
         result_abs = result
         result_asc = result
         return
      end if
      result_gauss = wg(4) * center_value
      result_kronrod = wgk(8) * center_value
      result_abs = wgk(8) * abs(center_value)

      do j = 1, 7
         f_left(j) = function_to_integrate(center - half_length * xgk(j))
         f_right(j) = function_to_integrate(center + half_length * xgk(j))
         if (.not. ieee_is_finite(f_left(j)) .or. &
             .not. ieee_is_finite(f_right(j))) then
            status = 3
            result = ieee_value(0.0_dp, ieee_quiet_nan)
            absolute_error = result
            result_abs = result
            result_asc = result
            return
         end if
         pair_sum = f_left(j) + f_right(j)
         result_kronrod = result_kronrod + wgk(j) * pair_sum
         result_abs = result_abs + wgk(j) * (abs(f_left(j)) + abs(f_right(j)))
         select case (j)
         case (2)
            result_gauss = result_gauss + wg(1) * pair_sum
         case (4)
            result_gauss = result_gauss + wg(2) * pair_sum
         case (6)
            result_gauss = result_gauss + wg(3) * pair_sum
         end select
      end do

      mean_value = 0.5_dp * result_kronrod
      result_asc = wgk(8) * abs(center_value - mean_value)
      do j = 1, 7
         result_asc = result_asc + wgk(j) * &
            (abs(f_left(j) - mean_value) + abs(f_right(j) - mean_value))
      end do

      result = result_kronrod * half_length
      result_abs = result_abs * abs(half_length)
      result_asc = result_asc * abs(half_length)
      raw_error = abs((result_kronrod - result_gauss) * half_length)

      absolute_error = raw_error
      if (result_asc > 0.0_dp .and. raw_error > 0.0_dp) then
         absolute_error = result_asc * &
            min(1.0_dp, (200.0_dp * raw_error / result_asc)**1.5_dp)
      end if
      if (result_abs > tiny(1.0_dp) / (50.0_dp * epsilon(1.0_dp))) then
         absolute_error = max(50.0_dp * epsilon(1.0_dp) * result_abs, &
                              absolute_error)
      end if
   end subroutine kronrod15


   pure function signum(value) result(sign_value)
      real(dp), intent(in) :: value
      real(dp) :: sign_value

      if (value > 0.0_dp) then
         sign_value = 1.0_dp
      else if (value < 0.0_dp) then
         sign_value = -1.0_dp
      else
         sign_value = 0.0_dp
      end if
   end function signum


   pure function exp_negative_square(value) result(exponential)
      !! 稳定计算 exp(-value^2)，并在结果必然下溢前直接返回零。
      real(dp), intent(in) :: value
      real(dp) :: exponential

      if (abs(value) >= exp_underflow_limit) then
         exponential = 0.0_dp
      else
         exponential = exp(-value * value)
      end if
   end function exp_negative_square


   pure function stable_loss_factor(scaled_radius) result(value)
      !! 稳定计算 erf(s)-2*s*exp(-s^2)/sqrt(pi)。
      !! s 很小时两项都约为 2s/sqrt(pi)，直接相减会丢失有效数字；此时使用
      !! 从 s^3 起始的收敛级数。s>=0.25 时直接式已经具有充足相对精度。
      real(dp), intent(in) :: scaled_radius
      real(dp) :: value
      real(dp) :: series_sum, term, power, factorial
      integer :: k

      if (scaled_radius < 0.25_dp) then
         series_sum = 0.0_dp
         power = scaled_radius**3
         factorial = 1.0_dp
         do k = 1, 40
            if (k > 1) then
               factorial = factorial * real(k, dp)
               power = power * scaled_radius**2
            end if
            term = (-1.0_dp)**(k + 1) * 2.0_dp * real(k, dp) * power / &
                   (factorial * real(2 * k + 1, dp))
            series_sum = series_sum + term
            if (abs(term) <= 4.0_dp * epsilon(1.0_dp) * &
                             max(abs(series_sum), tiny(1.0_dp))) exit
         end do
         value = 2.0_dp * series_sum / sqrt(pi)
      else
         value = erf(scaled_radius) - 2.0_dp * scaled_radius * &
                 exp_negative_square(scaled_radius) / sqrt(pi)
      end if
   end function stable_loss_factor


   pure function validate_common_inputs(m, n, x_bounds, y_bounds, x, y, z, &
                                        time_limit, diffusivity) result(status)
      integer, intent(in) :: m, n
      real(dp), intent(in) :: x_bounds(2), y_bounds(2)
      real(dp), intent(in) :: x, y, z, time_limit, diffusivity
      integer :: status

      status = validate_spatial_inputs(m, n, x_bounds, y_bounds, x, y, z, &
                                       quadrature_global)
      if (.not. ieee_is_finite(time_limit) .or. &
          .not. ieee_is_finite(diffusivity)) then
         status = 2
      else if (time_limit <= 0.0_dp .or. diffusivity <= 0.0_dp) then
         status = 2
      end if
   end function validate_common_inputs


   pure function validate_spatial_inputs(m, n, x_bounds, y_bounds, x, y, z, &
                                         mode) result(status)
      integer, intent(in) :: m, n, mode
      real(dp), intent(in) :: x_bounds(2), y_bounds(2), x, y, z
      integer :: status

      status = 0
      if (.not. all(ieee_is_finite(x_bounds)) .or. &
          .not. all(ieee_is_finite(y_bounds)) .or. &
          .not. ieee_is_finite(x) .or. .not. ieee_is_finite(y) .or. &
          .not. ieee_is_finite(z)) then
         status = 2
      else
         if (m < 0 .or. n < 0 .or. z <= 0.0_dp) status = 2
         if (x_bounds(2) <= x_bounds(1) .or. &
             y_bounds(2) <= y_bounds(1)) status = 2
      end if
      if (mode /= quadrature_global .and. &
          mode /= quadrature_split_at_observer) status = 4
   end function validate_spatial_inputs

end module three_dimensional_analysis
