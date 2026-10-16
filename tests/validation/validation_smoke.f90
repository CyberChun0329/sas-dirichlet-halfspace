program validation_smoke
  use iso_fortran_env, only: real64,real128
  use validation_reference, only: ref_rule,ref_temperature,ref_moment,ref_deficit
  use validation_sas, only: sas_temperature,sas_uniform_deficit
  use mp_arithmetic, only: mp_real,mp_begin,mp_set_precision,mp_precision_of,mp_abi_bytes, &
    mp_version,mp,mp_integer,mp_compare,mp_string,mp_reset,operator(-),operator(/),operator(+),abs
  use mp_moments, only: mp_moment_table
  use three_dimensional_analysis, only: steady_integral_semianalytic
  implicit none
  real(real64), parameter :: bounds(4)=[-1.0_real64,1.0_real64,-1.0_real64,1.0_real64]
  real(real64) :: nodes(24),weights(24),a,b,closed,relative64
  type(mp_real) :: x,y,z,bd(4),lo(0:8,0:8),hi(0:8,0:8),change
  integer :: status
  if(precision(1.0_real128)<30) error stop 'binary128 unavailable'
  call ref_rule(24,nodes,weights)
  if(abs(sum(weights)-2.0_real64)>2e-15_real64) error stop 'independent GL mass'
  closed=2/acos(-1.0_real64)*atan(4.0_real64/3)
  a=ref_temperature(0.0_real64,0.0_real64,.5_real64,bounds,[1.0_real64],[1.0_real64], &
    -1.0_real64,40,.true.)
  if(abs(a-closed)>2e-14_real64) error stop 'independent uniform solid angle'
  a=ref_moment(0.0_real64,0.0_real64,.5_real64,bounds,3,2,32,.true.,.false.)
  if(abs(a)>1e-16_real64) error stop 'independent signed odd moment'
  a=ref_temperature(.17_real64,-.21_real64,.43_real64,bounds,[1.0_real64,.2_real64], &
    [1.0_real64,-.3_real64],.7_real64,40,.true.)
  b=sas_temperature(.17_real64,-.21_real64,.43_real64,.7_real64,bounds, &
    [1.0_real64,.2_real64],[1.0_real64,-.3_real64],1e-12_real64,status)
  if(status/=0.or.abs(a-b)>2e-12_real64) error stop 'independent transient/SAS interface'
  a=ref_deficit(.17_real64,-.21_real64,.43_real64,1.0_real64,1.0_real64,bounds,48,.true.)
  b=sas_uniform_deficit(.17_real64,-.21_real64,.43_real64,1.0_real64,1.0_real64,bounds, &
    1e-12_real64,status)
  if(status/=0.or.abs(a-b)>2e-12_real64) error stop 'independent accumulated deficit/SAS interface'
  call mp_begin(224)
  if(mp_abi_bytes()/=32) error stop 'unexpected MPFR ABI'
  x=mp_integer(0);y=x;z=mp_integer(10)
  bd=[mp_integer(-1),mp_integer(1),mp_integer(-1),mp_integer(1)]
  call mp_moment_table(x,y,z,bd,8,8,lo)
  if(mp_precision_of(lo(8,8))/=224) error stop 'actual MPFR224 precision'
  call mp_set_precision(320)
  x=mp_integer(0);y=x;z=mp_integer(10)
  bd=[mp_integer(-1),mp_integer(1),mp_integer(-1),mp_integer(1)]
  call mp_moment_table(x,y,z,bd,8,8,hi)
  if(mp_precision_of(hi(8,8))/=320) error stop 'actual MPFR320 precision'
  change=abs(lo(8,8)-hi(8,8))/abs(hi(8,8))
  if(mp_compare(change,mp('1e-35'))>=0) error stop 'unrounded MPFR224/320 difference'
  print '(A)', 'MPFR '//mp_version()//' 224/320 relative J88 change='//mp_string(change,30)
  call mp_reset()
  ! Retain the original binary64 failure diagnostic rather than silently fixing
  ! or routing it through MPFR. This is a single known diagnostic, not the bank.
  a=steady_integral_semianalytic(8,8,bounds(:2),bounds(3:),0.0_real64,0.0_real64,10.0_real64,status)
  if(status/=0) error stop 'binary64 recurrence status'
  b=ref_moment(0.0_real64,0.0_real64,10.0_real64,bounds,8,8,32,.false.,.false.)
  relative64=abs(a-b)/abs(b)
  if(relative64<.1_real64) error stop 'binary64 instability diagnostic disappeared'
  print '(A,ES16.7)', 'retained binary64 J88 z10 relative error=',relative64
  print '(A)', 'validation_smoke PASS: local interface checks'
end program
