! Periodic verification wrapper. This is separate from the retained timing code.
module periodic_validation
  use validation_sas, only: rk,qk,center_integrand
  use temperature_common, only: gauss_rule_t,init_rule
  use mp_moments, only: mp_steady
  use validation_reference, only: ref_periodic_context,ref_prepare_periodic, &
    ref_periodic_value,ref_periodic_derivative_value
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  implicit none
  private
  real(rk), parameter, public :: period=0.5_rk
  type, public :: periodic_context
    real(rk) :: z=0,L=0,steady=0
    character :: profile='U'
    type(gauss_rule_t) :: young,old
  end type
  type, public :: extrema_result
    real(rk) :: tmin=0,tmax=0,minimum=0,maximum=0,swing=0,sampled_swing=0
  end type
  type, public :: tail_result
    real(rk) :: temperature=0,swing=0,zeta_remainder=0,zeta_refinement=0
  end type
  ! Numeric CSV columns after j,L,profile,d,z. NaN means not evaluated.
  integer, parameter, public :: nmetrics=34
  integer, parameter, public :: m_tmin=1,m_tmax=2,m_min=3,m_max=4,m_swing=5, &
    m_sampled=6,m_steady=7,m_tailt=8,m_tails=9,m_zeta_rem=10,m_zeta_ref=11, &
    m_ref_min=12,m_ref_max=13,m_ref_swing=14,m_ref_value_diff=15,m_ref_swing_diff=16, &
    m_root_tmin=17,m_root_tmax=18,m_root_swing=19,m_root_swing_diff=20, &
    m_root_phase_diff=21,m_derivative_diff=22,m_derivative_residual=23, &
    m_phase_change=24,m_order_change=25,m_history_change=26,m_history_bound=27, &
    m_spatial48_diff=28,m_spatial48_change=29,m_order_swing=30,m_history_swing=31, &
    m_root_brackets=32,m_selected=33,m_endpoint_correction=34
  character(*), parameter, public :: periodic_header= &
    'j,L,profile,d,z,tmin,tmax,minimum,maximum,swing,sampled_swing,steady,'// &
    'tail_temperature,tail_swing,zeta_remainder,zeta_refinement,ref_min,ref_max,ref_swing,'// &
    'ref_temperature_difference,ref_swing_difference,root_tmin,root_tmax,root_swing,'// &
    'root_swing_difference,root_phase_difference,derivative_difference,derivative_residual,'// &
    'phase_change,order_change,history_change,history_bound,spatial48_difference,'// &
    'spatial48_change,order_swing,history_swing,root_brackets,selected,endpoint_correction'
  public :: periodic_prepare,periodic_deficit,periodic_value,periodic_derivative
  public :: periodic_extrema,tail_bounds,hurwitz_tail,periodic_depth
contains
  subroutine periodic_prepare(ctx,z,L,profile,ng)
    type(periodic_context), intent(out) :: ctx
    real(rk), intent(in) :: z,L
    character, intent(in) :: profile
    integer, intent(in), optional :: ng
    real(qk) :: c(0:2),lq,zq
    integer :: n,status
    if(z<=0.or.L<=0.or.(profile/='U'.and.profile/='Q')) error stop 'periodic: invalid source'
    n=64; if(present(ng)) n=ng
    if(n<4) error stop 'periodic: quadrature order too small'
    ctx%z=z;ctx%L=L;ctx%profile=profile
    call init_rule(ctx%young,n,status)
    if(status/=0) error stop 'periodic: young GL construction failed'
    call init_rule(ctx%old,max(12,n/4),status)
    if(status/=0) error stop 'periodic: old GL construction failed'
    lq=real(L,qk);zq=real(z,qk);c=0;c(0)=1
    if(profile=='Q') c(2)=-1/(lq*lq)
    ctx%steady=real(mp_steady(0.0_qk,0.0_qk,zq,[-lq,lq,-lq,lq],c,c,224),rk)
  end subroutine

  function periodic_deficit(ctx,age) result(v)
    ! Instantaneous S-T, integrated directly on [0,v0]; never the cumulative D.
    type(periodic_context), intent(in) :: ctx
    real(rk), intent(in) :: age
    real(rk) :: v,v0,x,s,c,y,t
    integer :: i
    if(age<=0) then
      v=ctx%steady;return
    end if
    v0=min(ctx%z/(2*sqrt(age)),10.0_rk);s=0;c=0
    if(age<2) then
      do i=1,ctx%young%n
        x=v0*(ctx%young%x(i)+1)/2
        y=ctx%young%w(i)*center_integrand(x,ctx%z,ctx%L,ctx%L,ctx%profile)-c
        t=s+y;c=(t-s)-y;s=t
      end do
    else
      do i=1,ctx%old%n
        x=v0*(ctx%old%x(i)+1)/2
        y=ctx%old%w(i)*center_integrand(x,ctx%z,ctx%L,ctx%L,ctx%profile)-c
        t=s+y;c=(t-s)-y;s=t
      end do
    end if
    v=s*v0/2
  end function

  function periodic_value(ctx,phase,duty,n_history) result(v)
    type(periodic_context), intent(in) :: ctx
    real(rk), intent(in) :: phase,duty
    integer, intent(in), optional :: n_history
    real(rk) :: v,a,b,c,y,t
    integer :: n,k
    n=400;if(present(n_history)) n=n_history
    if(duty<0.or.duty>1.or.n<1.or.phase<0.or.phase>period) error stop 'periodic: invalid waveform query'
    v=0;c=0
    if(duty==0) return
    do k=0,n-1
      a=phase+period*k;b=a-period*duty
      y=periodic_deficit(ctx,b)-periodic_deficit(ctx,a)-c
      t=v+y;c=(t-v)-y;v=t
    end do
  end function

  function step_derivative(ctx,age) result(v)
    type(periodic_context), intent(in) :: ctx
    real(rk), intent(in) :: age
    real(rk) :: v,x
    v=0
    if(age<=0) return
    x=ctx%z/(2*sqrt(age))
    if(x>=10) return
    v=center_integrand(x,ctx%z,ctx%L,ctx%L,ctx%profile)*x/(2*age)
  end function

  function periodic_derivative(ctx,phase,duty,n_history) result(v)
    type(periodic_context), intent(in) :: ctx
    real(rk), intent(in) :: phase,duty
    integer, intent(in), optional :: n_history
    real(rk) :: v,a,b,c,y,t
    integer :: n,k
    n=400;if(present(n_history)) n=n_history
    if(duty<0.or.duty>1.or.n<1) error stop 'periodic: invalid derivative query'
    v=0;c=0
    do k=0,n-1
      a=phase+period*k;b=a-period*duty
      y=step_derivative(ctx,a)-step_derivative(ctx,b)-c
      t=v+y;c=(t-v)-y;v=t
    end do
  end function

  function bounded_phase(ctx,duty,n,lo,hi,sign) result(t)
    ! Golden-section bounded minimization; both signs are searched on each cell.
    type(periodic_context), intent(in) :: ctx
    real(rk), intent(in) :: duty,lo,hi,sign
    integer, intent(in) :: n
    real(rk) :: t,a,b,x1,x2,f1,f2
    real(rk), parameter :: ratio=0.6180339887498948482045868343656381_rk
    integer :: k
    a=lo;b=hi;x1=b-ratio*(b-a);x2=a+ratio*(b-a)
    f1=sign*periodic_value(ctx,x1,duty,n);f2=sign*periodic_value(ctx,x2,duty,n)
    do k=1,100
      if(b-a<=2e-12_rk) exit
      if(f1<=f2) then
        b=x2;x2=x1;f2=f1;x1=b-ratio*(b-a);f1=sign*periodic_value(ctx,x1,duty,n)
      else
        a=x1;x1=x2;f1=f2;x2=a+ratio*(b-a);f2=sign*periodic_value(ctx,x2,duty,n)
      end if
    end do
    if(f1<=f2) then
      t=x1
    else
      t=x2
    end if
  end function

  subroutine periodic_extrema(ctx,duty,n,seeds,result,endpoint_search)
    type(periodic_context), intent(in) :: ctx
    real(rk), intent(in) :: duty
    integer, intent(in) :: n,seeds
    type(extrema_result), intent(out) :: result
    logical, intent(in), optional :: endpoint_search
    real(rk) :: p(0:seeds),v(0:seeds),lo(0:seeds+3),hi(0:seeds+3),h,t,f
    integer :: i,j,k,count
    logical :: endpoints
    if(seeds<4) error stop 'periodic: phase grid too small'
    endpoints=.true.;if(present(endpoint_search)) endpoints=endpoint_search
    h=period/seeds
    do i=0,seeds
      p(i)=i*h;v(i)=periodic_value(ctx,p(i),duty,n)
    end do
    result=extrema_result(0.0_rk,0.0_rk,v(0),v(0),0.0_rk,maxval(v)-minval(v))
    call accept(period,v(seeds));t=duty*period;call accept(t,periodic_value(ctx,t,duty,n))
    count=0
    do i=1,seeds-1
      if((v(i)<=v(i-1).and.v(i)<=v(i+1)).or.(v(i)>=v(i-1).and.v(i)>=v(i+1))) then
        count=count+1;lo(count)=p(i-1);hi(count)=p(i+1)
      end if
    end do
    if(endpoints) then
      count=count+1;lo(count)=0;hi(count)=h
      count=count+1;lo(count)=period-h;hi(count)=period
      count=count+1;lo(count)=max(0.0_rk,duty*period-h);hi(count)=min(period,duty*period+h)
    end if
    do j=1,count
      do k=1,2
        t=bounded_phase(ctx,duty,n,lo(j),hi(j),real(3-2*k,rk))
        f=periodic_value(ctx,t,duty,n);call accept(t,f)
      end do
    end do
    result%swing=result%maximum-result%minimum
  contains
    subroutine accept(phase,value)
      real(rk), intent(in) :: phase,value
      if(value<result%minimum) then
        result%minimum=value;result%tmin=phase
      end if
      if(value>result%maximum) then
        result%maximum=value;result%tmax=phase
      end if
    end subroutine
  end subroutine

  function zeta_em(s,a,m,remainder) result(v)
    real(qk), intent(in) :: s,a
    integer, intent(in) :: m
    real(qk), intent(out) :: remainder
    real(qk) :: v,Aend,rising
    real(qk), parameter :: coeff(6)=[1.0_qk/12,-1.0_qk/720,1.0_qk/30240, &
      -1.0_qk/1209600,1.0_qk/47900160,-691.0_qk/1307674368000.0_qk]
    integer :: i,k,j
    v=0
    do i=0,m-1
      v=v+(a+i)**(-s)
    end do
    Aend=a+m;v=v+Aend**(1-s)/(s-1)+Aend**(-s)/2
    do k=1,6
      rising=1
      do j=0,2*k-2
        rising=rising*(s+j)
      end do
      v=v+coeff(k)*rising*Aend**(-s-2*k+1)
    end do
    ! |B_12(x)| <= |B_12| bounds the periodic-Bernoulli remainder.
    remainder=abs(coeff(6))*rising*Aend**(-s-11)
  end function

  subroutine hurwitz_tail(s,a,value,remainder,refinement)
    real(qk), intent(in) :: s,a
    real(qk), intent(out) :: value,remainder,refinement
    real(qk) :: coarse,coarse_rem
    if(s<=1.or.a<=0) error stop 'periodic: invalid Hurwitz zeta arguments'
    coarse=zeta_em(s,a,16,coarse_rem)
    value=zeta_em(s,a,32,remainder)
    refinement=abs(value-coarse)
  end subroutine

  function tail_bounds(z,L,profile,duty,n) result(b)
    real(rk), intent(in) :: z,L,duty
    character, intent(in) :: profile
    integer, intent(in), optional :: n
    type(tail_result) :: b
    real(qk) :: zz(3),rem(3),change(3),area,K,geom,P,qpi,a
    integer :: history,i
    history=400;if(present(n)) history=n
    if(duty<0.or.duty>1.or.history<=duty) error stop 'periodic: invalid tail arguments'
    a=real(history,qk)-real(duty,qk);P=real(period,qk);qpi=acos(-1.0_qk)
    area=4*real(L,qk)**2
    if(profile=='Q') area=16*real(L,qk)**2/9
    if(profile/='U'.and.profile/='Q') error stop 'periodic: invalid tail profile'
    K=real(z,qk)*area/(8*qpi**1.5_qk);geom=(2*real(L,qk)**2+real(z,qk)**2)/4
    do i=1,3
      call hurwitz_tail(real(i,qk)+1.5_qk,a,zz(i),rem(i),change(i))
    end do
    b%zeta_remainder=real(maxval(rem),rk);b%zeta_refinement=real(maxval(change),rk)
    ! Include analytic remainders and a rounding allowance when forming upper bounds.
    zz=zz+rem+32*epsilon(1.0_qk)*abs(zz)
    b%temperature=real(K*real(duty,qk)*P**(-1.5_qk)*zz(1),rk)*(1+16*epsilon(1.0_rk))
    b%swing=real(K*real(duty,qk)*(2.5_qk*P**(-1.5_qk)*zz(2)+ &
      geom*P**(-2.5_qk)*zz(3)),rk)*(1+16*epsilon(1.0_rk))
  end function

  subroutine stationary_root(ref,phase,duty,n,root,bracket,residual)
    type(ref_periodic_context), intent(in) :: ref
    real(rk), intent(in) :: phase,duty
    integer, intent(in) :: n
    real(rk), intent(out) :: root,residual
    logical, intent(out) :: bracket
    real(rk) :: lo,hi,fl,fh,mid,fm
    integer :: k
    lo=max(1e-12_rk,phase-period/128);hi=min(period-1e-12_rk,phase+period/128)
    fl=ref_periodic_derivative_value(ref,lo,duty,n);fh=ref_periodic_derivative_value(ref,hi,duty,n)
    bracket=fl*fh<0
    root=ieee_value(0.0_rk,ieee_quiet_nan);residual=root
    if(.not.bracket) return
    do k=1,100
      mid=(lo+hi)/2;fm=ref_periodic_derivative_value(ref,mid,duty,n)
      if(fm==0) then
        lo=mid;hi=mid;exit
      end if
      if((fm>0).eqv.(fl>0)) then
        lo=mid;fl=fm
      else
        hi=mid;fh=fm
      end if
      if(hi-lo<=2e-12_rk+1e-12_rk*abs(mid)) exit
    end do
    root=(lo+hi)/2;residual=abs(ref_periodic_derivative_value(ref,root,duty,n))
  end subroutine

  subroutine periodic_depth(j,L,profile,duty,z,v)
    integer, intent(in) :: j
    real(rk), intent(in) :: L,duty,z
    character, intent(in) :: profile
    real(rk), intent(out) :: v(nmetrics)
    type(periodic_context) :: ctx,fine
    type(ref_periodic_context) :: ref24,ref32,ref48
    type(extrema_result) :: base,dense,order_fine,history_fine,naive
    type(tail_result) :: tails,tails800
    real(rk) :: roots(2),residual(2),values(2),deriv(2),lo,hi,lo32,hi32
    logical :: brackets(2),selected
    integer :: k
    v=ieee_value(0.0_rk,ieee_quiet_nan)
    call periodic_prepare(ctx,z,L,profile)
    call periodic_extrema(ctx,duty,400,128,base)
    call periodic_extrema(ctx,duty,400,512,dense)
    tails=tail_bounds(z,L,profile,duty)
    v(1:11)=[base%tmin,base%tmax,base%minimum,base%maximum,base%swing, &
      base%sampled_swing,ctx%steady,tails%temperature,tails%swing,tails%zeta_remainder,tails%zeta_refinement]
    call ref_prepare_periodic(ref32,z,L,profile,32)
    lo=ref_periodic_value(ref32,base%tmin,duty,400);hi=ref_periodic_value(ref32,base%tmax,duty,400)
    v(m_ref_min:m_ref_swing_diff)=[lo,hi,hi-lo,max(abs(lo-base%minimum),abs(hi-base%maximum)), &
      abs(hi-lo-base%swing)]
    call ref_prepare_periodic(ref24,z,L,profile,24)
    call stationary_root(ref24,base%tmin,duty,400,roots(1),brackets(1),residual(1))
    call stationary_root(ref24,base%tmax,duty,400,roots(2),brackets(2),residual(2))
    v(m_root_brackets)=real(count(brackets),rk)
    if(all(brackets)) then
      do k=1,2
        values(k)=ref_periodic_value(ref32,roots(k),duty,400)
        deriv(k)=abs(periodic_derivative(ctx,roots(k),duty,400)- &
          ref_periodic_derivative_value(ref24,roots(k),duty,400))
      end do
      v(m_root_tmin:m_derivative_residual)=[roots(1),roots(2),values(2)-values(1), &
        abs(values(2)-values(1)-base%swing),maxval(abs(roots-[base%tmin,base%tmax])), &
        maxval(deriv),maxval(residual)]
    end if
    v(m_phase_change)=abs(dense%swing-base%swing)
    ! Fixed controls: the 12 L=1 table depths and two L=.2,Q,d=.1 depths.
    selected=j==7.or.(j==0.and.profile=='Q'.and.abs(duty-0.1_rk)<1e-14_rk)
    v(m_selected)=0
    if(selected) then
      v(m_selected)=1
      call periodic_prepare(fine,z,L,profile,128)
      call periodic_extrema(fine,duty,400,256,order_fine)
      call periodic_extrema(fine,duty,800,512,history_fine)
      tails800=tail_bounds(z,L,profile,duty,800)
      call ref_prepare_periodic(ref48,z,L,profile,48)
      lo=ref_periodic_value(ref48,order_fine%tmin,duty,400)
      hi=ref_periodic_value(ref48,order_fine%tmax,duty,400)
      lo32=ref_periodic_value(ref32,order_fine%tmin,duty,400)
      hi32=ref_periodic_value(ref32,order_fine%tmax,duty,400)
      v(m_order_change:m_history_swing)=[abs(order_fine%swing-base%swing), &
        abs(history_fine%swing-order_fine%swing),tails%swing+tails800%swing, &
        abs(hi-lo-order_fine%swing),abs((hi-lo)-(hi32-lo32)),order_fine%swing,history_fine%swing]
    end if
    if(j==7.and.profile=='Q'.and.abs(duty-0.4_rk)<1e-14_rk.and.abs(z-0.1_rk)<1e-14_rk) then
      call periodic_extrema(ctx,duty,400,128,naive,.false.)
      v(m_endpoint_correction)=base%swing-naive%swing
    end if
  end subroutine
end module
