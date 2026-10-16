program summarize_validation
   use iso_fortran_env, only: real64
   use ieee_arithmetic, only: ieee_is_finite
   implicit none
   integer,parameter :: rk=real64
   integer,parameter :: s2m(6)=[1,0,2,3,4,5],s2n(6)=[0,2,1,2,3,5]
   integer :: u,o,t,ios,i,j,k,n,row,part,mark(401,2),m,expected(8),unit2
   real(rk) :: maxsas,maxref,e(9,5),z(9),v,tableerr,tableref,lo,hi,mid,q90,pi,nodez
   character(4096) :: line
   character(40) :: name
   character(16) :: tag
   character(2) :: profile
   character(32),parameter :: banks(8)=[character(32)::'stability.csv','moment_regression.csv','gamma_regression.csv','signed_transient.csv','final_corner_checks.csv','final_operator_checks.csv','case1_profiles.csv','case2_profiles.csv']
   expected=[45,72,21,36,12,12,150,600]
   open(newunit=o,file='results/validation_summary.txt',status='new')
   call require_pass('core_summary.txt')
   call require_pass('case1_summary.txt')
   call require_pass('gaussian_summary.txt')
   call require_pass('precision_summary.txt')
   call require_pass('table_s2_summary.txt')
   call read_open('table_s2.csv',u);n=0
   do
      read(u,'(a)',iostat=ios)line
      if(ios/=0)exit
      n=n+1
      if(n>6)error stop 'extra Table S2 case'
      if(nint(number(line,1))/=s2m(n).or.nint(number(line,2))/=s2n(n))error stop 'Table S2 index identity'
      if(number(line,7)>=1e-11_rk.or.number(line,8)>=1e-11_rk)error stop 'Table S2 threshold'
   end do
   close(u);if(n/=6)error stop 'incomplete original six Table S2 cases'
   write(o,*)'table_s2_cases=',n
   do i=1,8
      call read_open(trim(banks(i)),u);n=0
      do
         read(u,'(a)',iostat=ios)line
         if(ios/=0)exit
         n=n+1
         do k=1,fields(line)
            v=number(line,k)
         end do
      end do
      close(u)
      write(o,'(a,i0)')trim(banks(i))//' rows=',n
      if(n/=expected(i))error stop 'bank count mismatch'
   end do
   call read_open('case2_convergence.csv',u);n=0
   do
      read(u,'(a)',iostat=ios)line
      if(ios/=0)exit
      n=n+1
      if(n>6)error stop 'extra Gaussian convergence record'
      if(nint(number(line,1))/=2*n)error stop 'Gaussian order identity'
      do k=2,10
         v=number(line,k)
      end do
      if(max(number(line,3),number(line,8))>=1e-11_rk)error stop 'Gaussian evaluation threshold'
   end do
   close(u);if(n/=6)error stop 'incomplete Gaussian convergence records'
   write(o,*)'Gaussian_convergence_orders=',n,' all_banks_and_precision_status=PASS'
   mark=0;maxsas=0;maxref=0;n=0
   do part=1,8
      write(tag,'(i0)')part
      call require_pass('roots_'//trim(tag)//'_summary.txt')
      call read_open('case3_grid_'//trim(tag)//'.csv',u)
      do
         read(u,'(a)',iostat=ios)line
         if(ios/=0)exit
         row=nint(number(line,1))+1;profile=field(line,3)
         j=1;if(profile=='Q')j=2
         if(row<1.or.row>401.or.(profile/='U'.and.profile/='Q'))error stop 'bad fresh node'
         if(mark(row,j)/=0)error stop 'duplicate original node'
         mark(row,j)=1;n=n+1
         nodez=.1_rk+.00975_rk*(row-1)
         if(abs(number(line,2)-nodez)>2e-15_rk)error stop 'fresh depth identity'
         do k=4,11
            v=number(line,k)
         end do
         maxsas=max(maxsas,number(line,7));maxref=max(maxref,number(line,8))
         if(abs(number(line,9)-log(nodez))>2e-15_rk)error stop 'fresh logarithmic plot coordinate'
         if(abs(number(line,10)-number(line,4)/nodez**2)>2e-14_rk)error stop 'fresh SAS plot coordinate'
         if(abs(number(line,11)-number(line,6)/nodez**2)>2e-14_rk)error stop 'fresh reference plot coordinate'
      end do
      close(u)
   end do
   if(n/=802.or.any(mark/=1))error stop 'incomplete original grid'
   write(o,*)'fresh_depth_nodes=',401,' fresh_roots=',n
   write(o,*)'fresh_max_SAS_relative=',maxsas,' reference_refinement=',maxref
   call read_open('case3_table.csv',u);n=0;tableerr=0;tableref=0
   do
      read(u,'(a)',iostat=ios)line
      if(ios/=0)exit
      n=n+1;tableerr=max(tableerr,number(line,7));tableref=max(tableref,number(line,8))
   end do
   close(u);if(n/=24)error stop 'table roots incomplete'
   write(o,*)'table_root_count=',n,' max_relative=',tableerr,' reference_refinement=',tableref
   if(maxsas>=1e-9_rk.or.tableerr>=1e-9_rk.or.maxref>=1e-11_rk.or.tableref>=1e-11_rk)error stop 'root validation threshold'
   call read_open('stability.csv',u);n=0
   do
      read(u,'(a)',iostat=ios)line
      if(ios/=0)exit
      n=n+1;i=(n-1)/5+1;j=mod(n-1,5)+1
      z(i)=number(line,1);e(i,j)=number(line,9)
   end do
   close(u)
   open(newunit=t,file='results/table_s3_fortran.md',status='new')
   write(t,'(a)')'| $z$ | $(0,0)$ | $(2,2)$ | $(4,4)$ | $(6,6)$ | $(8,8)$ |'
   write(t,'(a)')'|---:|---:|---:|---:|---:|---:|'
   do i=1,9
      write(t,'(a,f5.1,a)',advance='no')'| ',z(i),' |'
      do j=1,5
         if(e(i,j)<1e-12_rk)then
            write(t,'(a)',advance='no')' $<10^{-12}$ |'
         else
            write(t,'(1x,es10.2,a)',advance='no')e(i,j),' |'
         end if
      end do
      write(t,*)
   end do
   close(t)
   ! Recompute the scalar far-field constant, without a historical solver value.
   pi=acos(-1.0_rk);lo=0;hi=2
   do i=1,100
      mid=(lo+hi)/2
      v=erfc(sqrt(mid))+2*sqrt(mid/pi)*exp(-mid)
      if(v>.9_rk)then
         lo=mid
      else
         hi=mid
      end if
   end do
   q90=(lo+hi)/2
   write(o,*)'q90=',q90,' c90=',1/(4*q90)
   write(o,*)'status=PASS';close(o)
   print *,'non-periodic aggregation: PASS; fresh roots=',sum(mark)
contains
   subroutine require_pass(file)
      character(*),intent(in) :: file
      integer :: unit,status
      logical :: passed
      character(4096) :: s
      open(newunit=unit,file='results/'//file,status='old',action='read',iostat=status)
      if(status/=0)error stop 'missing bank completion record'
      passed=.false.
      do
         read(unit,'(a)',iostat=status)s
         if(status/=0)exit
         s=adjustl(s)
         if(trim(s)=='status=PASS'.or.trim(s)=='status: PASS')passed=.true.
         if(index(s,'FAIL')>0)error stop 'bank reports failure'
      end do
      close(unit)
      if(.not.passed)error stop 'bank has no successful completion'
   end subroutine
   subroutine read_open(file,unit)
      character(*),intent(in) :: file
      integer,intent(out) :: unit
      open(newunit=unit,file='results/'//file,status='old',action='read',iostat=ios)
      if(ios/=0)error stop 'missing result'
      read(unit,'(a)')line
   end subroutine
   function field(text,index) result(v)
      character(*),intent(in) :: text
      integer,intent(in) :: index
      character(256) :: v
      integer :: p,start,count
      start=1;count=1;v=''
      do p=1,len_trim(text)+1
         if(p<=len_trim(text))then
            if(text(p:p)/=',')cycle
         end if
         if(count==index)then
            v=adjustl(text(start:p-1));return
         end if
         start=p+1;count=count+1
      end do
      error stop 'missing CSV field'
   end function
   function number(text,index) result(v)
      character(*),intent(in) :: text
      integer,intent(in) :: index
      real(rk) :: v
      integer :: status
      character(256) :: s
      s=field(text,index);read(s,*,iostat=status)v
      if(status/=0)error stop 'invalid numerical CSV'
      if(.not.ieee_is_finite(v))error stop 'nonfinite numerical CSV'
   end function
   function fields(text)result(n)
      character(*),intent(in) :: text
      integer :: n,i
      n=1
      do i=1,len_trim(text)
         if(text(i:i)==',')n=n+1
      end do
   end function
end program
