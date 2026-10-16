module prepared_sas
   use temperature_common, only: rk, ok, invalid, numerical_failure, pi, sqrtpi, &
      patch_source_t, gauss_rule_t, query_work_t, init_source, valid_rule, &
      validate_query, mapped_rule, quiet_nan, stable_expm1
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   private
   integer, parameter :: moment_orders(25) = [0,1,2,3,4,5,6,7,8,9,10,11,12, &
                                            13,14,15,16,17,18,19,20,21,22,23,24]
   real(rk), parameter :: moments0(25) = merge(2.0_rk/(real(moment_orders,rk)+1.0_rk), &
                                             0.0_rk,mod(moment_orders,2)==0)
   real(rk), parameter :: moments1(25) = merge(2.0_rk/(real(moment_orders,rk)+2.0_rk), &
                                             0.0_rk,mod(moment_orders,2)==1)
   real(rk), parameter :: moments2(25) = merge(4.0_rk/(real(moment_orders,rk)+3.0_rk) &
                                             -2.0_rk/(real(moment_orders,rk)+1.0_rk), &
                                             0.0_rk,mod(moment_orders,2)==0)
   ! Input-independent library constants, not a hidden per-source/query cost.
   real(rk), parameter :: moment_library(0:2,0:24) = &
      transpose(reshape([moments0,moments1,moments2],[25,3]))

   ! Only source geometry and degree moments are prepared. All shifted Gaussian
   ! moments, panel geometry, and final responses are fresh query work.
   type, public :: prepared_sas_t
      private
      type(patch_source_t) :: source
      logical :: ready = .false.
   end type
   public :: prepare_sas, evaluate_sas
contains
   subroutine prepare_sas(source,prepared,status)
      type(patch_source_t), intent(in) :: source
      type(prepared_sas_t), intent(out) :: prepared
      integer, intent(out) :: status
      prepared = prepared_sas_t()
      status = invalid
      if (.not.source%ready) return
      ! Revalidate the public source fields and recompute effective degrees.
      call init_source(prepared%source,source%c,source%a,source%b,source%alpha, &
                       status,source%cx,source%cy)
      if (status /= ok) return
      prepared%ready = .true.
   end subroutine

   function evaluate_sas(prepared,query,rule,status,work) result(value)
      type(prepared_sas_t), intent(in) :: prepared
      real(rk), intent(in) :: query(4)
      type(gauss_rule_t), intent(in) :: rule
      integer, intent(out) :: status
      type(query_work_t), intent(out), optional :: work
      real(rk) :: value, x, y, z, dt, v0, h, span, u
      real(rk) :: ax(0:2), ay(0:2), contraction, contraction_correction
      real(rk) :: integral, correction, term, candidates(8)
      real(rk), allocatable :: v(:), w(:)
      integer :: i, p, q
      value = quiet_nan()
      status = invalid
      if (present(work)) work = query_work_t()
      if (.not.prepared%ready) return
      if (.not.valid_rule(rule)) return
      call validate_query(prepared%source,query,x,y,z,dt,status)
      if (status /= ok) return

      ! The omitted step tail is bounded absolutely by sum(abs(C))*erfc(12).
      ! This does not certify relative accuracy for weak or cancelling signals.
      ! Compare before dividing so a very large v0 cannot overflow.
      if (z >= 24.0_rk*sqrt(dt)) then
         value = 0.0_rk
         return
      end if
      v0 = z/(2.0_rk*sqrt(dt))
      span = max(prepared%source%a+abs(x),prepared%source%b+abs(y))
      h = z/span
      candidates = [h/4.0_rk,h,4.0_rk*h,16.0_rk*h,0.25_rk,1.0_rk,3.0_rk,6.0_rk]
      call mapped_rule(v0,12.0_rk,candidates,rule,v,w,status)
      if (status /= ok) return
      integral = 0.0_rk
      correction = 0.0_rk
      do i = 1, size(v)
         if (present(work)) then
            work%v_nodes = work%v_nodes+1
            work%kernel_evals = work%kernel_evals+1
         end if
         u = v(i)/z
         call axis_values(u,x,prepared%source%a,prepared%source%mx,moment_library,ax,status)
         if (status /= ok) return
         call axis_values(u,y,prepared%source%b,prepared%source%my,moment_library,ay,status)
         if (status /= ok) return
         ! u^2 Ax Ay is evaluated as (u Ax)(u Ay) to avoid overflowing u^2.
         ax = u*ax
         ay = u*ay
         contraction = 0.0_rk
         contraction_correction = 0.0_rk
         do p = 0, prepared%source%mx
            do q = 0, prepared%source%my
               term = prepared%source%c(q,p)*ax(p)*ay(q)
               if (.not.ieee_is_finite(term)) then
                  status = numerical_failure
                  return
               end if
               call compensated_add(contraction,contraction_correction,term)
            end do
         end do
         term = w(i)*exp(-v(i)*v(i))*(contraction+contraction_correction)
         if (.not.ieee_is_finite(term)) then
            status = numerical_failure
            return
         end if
         call compensated_add(integral,correction,term)
      end do
      value = (2.0_rk/(pi*sqrtpi))*(integral+correction)
      if (.not.ieee_is_finite(value)) then
         value = quiet_nan()
         status = numerical_failure
         return
      end if
      status = ok
   end function

   subroutine axis_values(u,position,width,degree,moments,values,status)
      real(rk), intent(in) :: u, position, width, moments(0:2,0:24)
      integer, intent(in) :: degree
      real(rk), intent(out) :: values(0:2)
      integer, intent(out) :: status
      real(rk) :: beta, delta, linear, quadratic, lo, hi, difference
      real(rk) :: k0, k1, k2, e0, e1, change, inverse_beta, factor, shift
      real(rk) :: previous, current, following, total(0:2), correction(0:2)
      real(rk) :: term, variation
      integer :: k, p
      logical :: small, huge_change
      values = 0.0_rk
      status = numerical_failure
      beta = width*u
      delta = position*u
      if (.not.all(ieee_is_finite([beta,delta]))) return
      ! The small region necessarily has beta<=1/2. This test also avoids
      ! squaring a large beta merely to decide which analytic branch to use.
      small = .false.
      if (beta <= 0.5_rk) then
         linear = (2.0_rk*beta)*delta
         quadratic = -beta*beta
         variation = abs(linear)+abs(quadratic)
         small = variation <= 0.25_rk
      end if
      if (degree == 0) then
         ! Match the oracle's direct-erf branch for the uniform modest case.
         small = small .and. (abs(delta) > 1.0_rk .or. &
                 abs(position)/width > 128.0_rk .or. u == 0.0_rk)
      end if
      if (small) then
         previous = 1.0_rk
         current = linear
         total = moments(:,0)+moments(:,1)*linear
         correction = 0.0_rk
         do k = 2, 24
            following = (linear*current+2.0_rk*quadratic*previous)/real(k,rk)
            do p = 0, degree
               term = moments(p,k)*following
               call compensated_add(total(p),correction(p),term)
            end do
            previous = current
            current = following
            ! Two successive coefficients cover the recurrence's parity zeros.
            if (k >= 4 .and. max(abs(previous),abs(current)) < 1.0e-20_rk) exit
         end do
         values(0:degree) = width*gaussian(delta)*(total(0:degree)+correction(0:degree))
      else
         if (beta == 0.0_rk) return
         if (beta > huge(1.0_rk)-abs(delta)) return
         lo = -beta-delta
         hi = beta-delta
         if (max(abs(lo),abs(hi)) <= 1.0_rk) then
            difference = erf(hi)-erf(lo)
         else if (lo >= 0.0_rk) then
            difference = erfc(lo)-erfc(hi)
         else if (hi <= 0.0_rk) then
            difference = erfc(-hi)-erfc(-lo)
         else
            difference = erf(hi)-erf(lo)
         end if
         k0 = (sqrtpi/2.0_rk)*difference
         factor = width/beta
         values(0) = factor*k0
         if (degree > 0) then
            e0 = gaussian(lo)
            e1 = gaussian(hi)
            ! Underflowed Gaussian tails contain no recoverable binary64
            ! shifted moments; do not attempt a hidden spatial-reference repair.
            if (k0 == 0.0_rk .and. e0 == 0.0_rk .and. e1 == 0.0_rk) then
               status = ok
               return
            end if
            ! change=-4 beta delta. A magnitude beyond binary64's range gives
            ! expm1(-abs(change))=-1, without forming an overflowing product.
            huge_change = .false.
            if (beta > 0.25_rk) then
               huge_change = abs(delta) > (huge(1.0_rk)/4.0_rk)/beta
            end if
            if (huge_change) then
               if (delta <= 0.0_rk) then
                  k1 = 0.5_rk*e0
               else
                  k1 = -0.5_rk*e1
               end if
            else
               change = -4.0_rk*(beta*delta)
               if (change >= 0.0_rk) then
                  k1 = -0.5_rk*e0*stable_expm1(-change)
               else
                  k1 = 0.5_rk*e1*stable_expm1(change)
               end if
            end if
            inverse_beta = 1.0_rk/beta
            shift = position/width
            values(1) = factor*(shift*k0+k1*inverse_beta)
            if (degree == 2) then
               if (abs(shift) > sqrt(huge(1.0_rk))/4.0_rk) return
               k2 = (k0-hi*e1+lo*e0)/2.0_rk
               values(2) = factor*((2.0_rk*shift*shift-1.0_rk)*k0 &
                  +4.0_rk*shift*k1*inverse_beta+2.0_rk*k2*inverse_beta*inverse_beta)
            end if
         end if
      end if
      if (.not.all(ieee_is_finite(values(0:degree)))) return
      status = ok
      ! Analytic binary64 moment differences and shifted contractions can lose
      ! relative accuracy outside the patch and for signed cancellations.
      ! A finite result is not a relative error certificate; refine the v rule
      ! and compare independent references where that accuracy is required.
   end subroutine

   pure function gaussian(argument) result(value)
      real(rk), intent(in) :: argument
      real(rk) :: value
      ! Very large finite arguments must not overflow when squared.
      if (abs(argument) > sqrt(huge(1.0_rk))) then
         value = 0.0_rk
      else
         value = exp(-argument*argument)
      end if
   end function

   pure subroutine compensated_add(total,correction,term)
      real(rk), intent(inout) :: total, correction
      real(rk), intent(in) :: term
      real(rk) :: following
      following = total+term
      if (abs(total) >= abs(term)) then
         correction = correction+((total-following)+term)
      else
         correction = correction+((term-following)+total)
      end if
      total = following
   end subroutine
end module
