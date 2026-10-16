program gaussian_coeff_check
  use iso_fortran_env, only: real128
  use ieee_arithmetic, only: ieee_is_finite
  use mp_moments, only: mp_gaussian_coeff
  implicit none
  integer,parameter :: q=real128,orders(6)=[2,4,6,8,10,12]
  real(q),parameter :: atol=1e-27_q,rtol=1e-25_q
  real(q) :: b224(0:12),b320(0:12),c512(0:12),c1024(0:12),p512(0:12),p1024(0:12)
  real(q) :: dcos,dbessel,dprec,limit,maxratio
  integer :: k,m,j,u,s,holds
  character(1024) :: outdir
  character(8) :: state
  call get_command_argument(1,outdir)
  if(len_trim(outdir)==0)error stop 'gaussian_coeff_check NEW_OUTPUT_DIR'
  open(newunit=u,file=trim(outdir)//'/gaussian_power_coefficients.csv',status='new',action='write')
  open(newunit=s,file=trim(outdir)//'/gaussian_coeff_summary.txt',status='new',action='write')
  write(u,'(a)')'M,power,Bessel224,Bessel320,cosine512_power,cosine1024_power,cosine_change,Bessel_cosine_difference,Bessel_precision_change,mixed_limit,state'
  holds=0;maxratio=0
  do k=1,size(orders)
    m=orders(k)
    call mp_gaussian_coeff(m,224,b224(:m))
    call mp_gaussian_coeff(m,320,b320(:m))
    call cosine_coefficients(m,512,c512(:m))
    call cosine_coefficients(m,1024,c1024(:m))
    call powers(m,c512(:m),p512(:m))
    call powers(m,c1024(:m),p1024(:m))
    do j=0,m
      dcos=abs(p512(j)-p1024(j));dbessel=abs(b320(j)-p1024(j));dprec=abs(b224(j)-b320(j))
      limit=atol+rtol*abs(p1024(j));state='PASS'
      if(.not.all(ieee_is_finite([dcos,dbessel,dprec])).or.max(dcos,dbessel,dprec)>limit)then
        holds=holds+1;state='HOLD'
      end if
      maxratio=max(maxratio,max(dcos,dbessel,dprec)/limit)
      write(u,'(i0,",",i0,8(",",es45.34e4),",",a)')m,j,b224(j),b320(j),p512(j),p1024(j),dcos,dbessel,dprec,limit,trim(state)
    end do
  end do
  close(u)
  write(s,'(a)')'target: exp(-4 X^2); finite even Chebyshev degrees M=2,4,6,8,10,12'
  write(s,'(a)')'reference: independent binary128 midpoint cosine sums at 512/1024 nodes, then binary128 Chebyshev-to-power conversion'
  write(s,'(a)')'comparison: existing MPFR Bessel construction at 224/320 bits, rounded once to binary128 physical-power coefficients'
  write(s,'(a)')'coefficient differences here are binary128 differences, not an unrounded MPFR precision certificate'
  write(s,'(a,2(1x,es45.34e4))')'fixed absolute/relative tolerances:',atol,rtol
  write(s,'(a,es45.34e4)')'maximum_difference_over_mixed_limit: ',maxratio
  write(s,'(a,i0)')'HOLD_coefficients: ',holds
  write(s,'(a)')'scope: coefficient cross-check at the degrees and tolerances listed above'
  if(holds==0)then
    write(s,'(a)')'status=PASS'
  else
    write(s,'(a)')'status=HOLD'
  end if
  close(s)
  if(holds/=0)error stop 'coefficient HOLD; retain detailed output'
contains
  subroutine cosine_coefficients(m,n,c)
    integer,intent(in) :: m,n
    real(q),intent(out) :: c(0:m)
    real(q) :: theta,total,pi
    integer :: i,k
    pi=acos(-1.0_q)
    do k=0,m
      if(mod(k,2)==1)then
        c(k)=0;cycle
      end if
      total=0
      do i=1,n
        theta=(real(i,q)-.5_q)*pi/n
        total=total+exp(-4.0_q*cos(theta)**2)*cos(k*theta)
      end do
      c(k)=total/n
      if(k>0)c(k)=2*c(k)
    end do
  end subroutine
  subroutine powers(m,c,p)
    integer,intent(in) :: m
    real(q),intent(in) :: c(0:m)
    real(q),intent(out) :: p(0:m)
    real(q) :: t(0:m,0:m)
    integer :: k,j
    t=0;t(0,0)=1
    if(m>=1)t(1,1)=1
    do k=2,m
      do j=1,k
        t(j,k)=2*t(j-1,k-1)
      end do
      t(:,k)=t(:,k)-t(:,k-2)
    end do
    p=matmul(t,c)
  end subroutine
end program
