module unified_poly_reference
   use replay_types
   use replay_case_bank, only: case_spec_t
   use spatial_routes, only: spatial_t,prepare_spatial,evaluate_spatial,evaluate_shared
   use unified_engine, only: candidate_errors
   implicit none
   private
   type,public :: reference_t
      real(rk),allocatable :: values(:),v48(:),v72(:),v96(:),shared(:)
      real(rk) :: EA=0,ER=0,crossEA=0,crossER=0,floor=0,afloor=0
      integer :: order=0
      logical :: certified=.false.
   end type
   public :: fresh_reference,write_reference,read_reference
contains
   subroutine fresh_reference(c,r,status)
      type(case_spec_t),intent(in) :: c
      type(reference_t),intent(out) :: r
      integer,intent(out) :: status
      type(spatial_t) :: truth
      type(gauss_rule_t) :: r48,r72,r96
      logical :: active
      integer :: i,st
      status=invalid
      allocate(r%values(size(c%q,2)),r%v48(size(c%q,2)),r%v72(size(c%q,2)), &
         r%v96(size(c%q,2)),r%shared(size(c%q,2)))
      r%values=quiet_nan();r%v48=r%values;r%v72=r%values;r%v96=r%values;r%shared=r%values
      call prepare_spatial(c%source,truth,status)
      if(status/=ok) return
      call init_rule(r48,48,status)
      if(status/=ok) return
      call init_rule(r72,72,status)
      if(status/=ok) return
      call init_rule(r96,96,status)
      if(status/=ok) return
      do i=1,size(r%values)
         r%v48(i)=evaluate_spatial(truth,c%q(:,i),r48,r48,c%observable,st)
         if(st/=ok) then
            status=st;return
         end if
         r%v72(i)=evaluate_spatial(truth,c%q(:,i),r72,r72,c%observable,st)
         if(st/=ok) then
            status=st;return
         end if
         r%v96(i)=evaluate_spatial(truth,c%q(:,i),r96,r96,c%observable,st)
         if(st/=ok) then
            status=st;return
         end if
         r%shared(i)=evaluate_shared(truth,c%q(:,i),r96,c%observable,st)
         if(st/=ok) then
            status=st;return
         end if
      end do
      r%values=r%v72;r%order=72
      call candidate_errors(r%v48,r%values,r%EA,r%ER,active)
      if(r%EA>1e-11_rk.or.r%ER>1e-9_rk) then
         r%values=r%v96;r%order=96
         call candidate_errors(r%v72,r%values,r%EA,r%ER,active)
      end if
      call candidate_errors(r%shared,r%values,r%crossEA,r%crossER,active)
      r%certified=all(ieee_is_finite(r%values)).and.r%EA<=1e-11_rk.and.r%ER<=1e-9_rk.and. &
         r%crossEA<=max(1e-11_rk,5*r%EA).and.r%crossER<=max(1e-9_rk,5*r%ER)
      r%floor=max(1e-12_rk,5*r%ER,5*r%crossER);r%afloor=max(1e-15_rk,5*r%EA,5*r%crossEA)
      status=ok
   end subroutine

   subroutine write_reference(u,c,r)
      integer,intent(in) :: u
      type(case_spec_t),intent(in) :: c
      type(reference_t),intent(in) :: r
      integer :: i
      write(u,'(a)')trim(c%id)
      write(u,*)size(c%q,2),c%observable,size(c%source%c,1),size(c%source%c,2)
      write(u,'(*(es26.17e3,1x))')c%source%a,c%source%b,c%source%alpha,c%source%cx,c%source%cy
      write(u,'(*(es26.17e3,1x))')c%source%c
      do i=1,size(c%q,2)
         write(u,'(*(es26.17e3,1x))')c%q(:,i)
      end do
      write(u,*)r%order,r%certified
      write(u,'(*(es26.17e3,1x))')r%EA,r%ER,r%crossEA,r%crossER,r%floor,r%afloor
      do i=1,size(c%q,2)
         write(u,'(*(es26.17e3,1x))')r%values(i),r%v48(i),r%v72(i),r%v96(i),r%shared(i)
      end do
   end subroutine

   subroutine read_reference(u,c,r,status)
      integer,intent(in) :: u
      type(case_spec_t),intent(in) :: c
      type(reference_t),intent(out) :: r
      integer,intent(out) :: status
      character(128) :: id
      integer :: np,obs,ny,nx,ios,i
      real(rk) :: geometry(5),q(4)
      real(rk),allocatable :: coeff(:,:)
      status=invalid
      read(u,'(a)',iostat=ios)id
      if(ios/=0.or.trim(id)/=trim(c%id)) return
      read(u,*,iostat=ios)np,obs,ny,nx
      if(ios/=0) return
      if(np/=size(c%q,2).or.obs/=c%observable.or.ny/=size(c%source%c,1).or.nx/=size(c%source%c,2)) return
      read(u,*,iostat=ios)geometry
      if(ios/=0) return
      if(any(geometry/=[c%source%a,c%source%b,c%source%alpha,c%source%cx,c%source%cy])) return
      allocate(coeff(ny,nx));read(u,*,iostat=ios)coeff
      if(ios/=0) return
      if(any(coeff/=c%source%c)) return
      do i=1,np
         read(u,*,iostat=ios)q
         if(ios/=0) return
         if(any(q/=c%q(:,i))) return
      end do
      allocate(r%values(np),r%v48(np),r%v72(np),r%v96(np),r%shared(np))
      read(u,*,iostat=ios)r%order,r%certified
      if(ios/=0) return
      read(u,*,iostat=ios)r%EA,r%ER,r%crossEA,r%crossER,r%floor,r%afloor
      if(ios/=0) return
      do i=1,np
         read(u,*,iostat=ios)r%values(i),r%v48(i),r%v72(i),r%v96(i),r%shared(i)
         if(ios/=0) return
      end do
      ! The release controller verifies file/binary hashes; this additionally
      ! prevents a different source, observable or query batch being paired.
      if(.not.all(ieee_is_finite(r%values))) return
      status=ok
   end subroutine
end module
