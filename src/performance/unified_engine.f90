! Mathematical adapter only. No old result or timing files are read here.
module unified_engine
   use replay_types
   use source_projection, only: project_source,truncate_poly_rank
   use sas_general, only: sas_general_t,prepare_sas_general,evaluate_sas_general
   use spatial_routes, only: spatial_t,prepare_spatial,evaluate_spatial,evaluate_shared
   use svd_factor_routes, only: poly_factor_t,prepare_poly_factor,evaluate_poly_factor
   implicit none
   private
   integer,parameter,public :: RT_SAS=1,RT_SOURCE_2D=2,RT_SOURCE_SHARED=3,RT_PROJECT_2D=4, &
      RT_PROJECT_SHARED=5,RT_POLY_FACTOR=6
   type,public :: descriptor_t
      integer :: route=0,m=-1,rank=-1,nx=0,ny=0
      character(128) :: key='',label=''
   end type
   type,public :: prepared_t
      type(descriptor_t) :: spec
      type(source_t) :: model
      type(sas_general_t) :: sas
      type(spatial_t) :: space
      type(poly_factor_t) :: factor
      type(gauss_rule_t) :: rx,ry
      logical :: ready=.false.
   end type
   public :: compact_source,prepare_candidate,evaluate_candidate,build_candidate_bank,make_descriptor
   public :: candidate_errors,preparation_fingerprint,same_candidate
contains
   subroutine compact_source(source,small,status)
      type(source_t),intent(in) :: source
      type(source_t),intent(out) :: small
      integer,intent(out) :: status
      if(source%kind==SRC_POLY.and.allocated(source%c)) then
         ! Remove only exact trailing zeros; M=0/1 truncations remain distinct.
         call init_poly(small,source%c(0:source%my,0:source%mx),source%a,source%b,source%alpha, &
            status,source%cx,source%cy)
      else
         small=source;status=invalid
         if(small%ready) status=ok
      end if
   end subroutine

   function make_descriptor(route,m,rank,nx,ny) result(d)
      integer,intent(in) :: route,m,rank,nx,ny
      type(descriptor_t) :: d
      d%route=route;d%m=m;d%rank=rank;d%nx=nx;d%ny=ny
      select case(route)
      case(RT_SAS);d%label='SAS_PROJECTED'
      case(RT_SOURCE_2D);d%label='SPACE_RECT_ERFC'
      case(RT_SOURCE_SHARED);d%label='SPACE_SHARED_SOURCE'
      case(RT_PROJECT_2D);d%label='SPACE_POLY'
      case(RT_PROJECT_SHARED);d%label='SPACE_SHARED_PROJECTED'
      case(RT_POLY_FACTOR);d%label='SPACE_POLY_FACTOR'
      end select
      write(d%key,'(a,"_M",i0,"_R",i0,"_nx",i0,"_ny",i0)')trim(d%label),m,rank,nx,ny
   end function

   pure logical function same_candidate(a,b)
      type(descriptor_t),intent(in) :: a,b
      same_candidate=a%route==b%route.and.a%m==b%m.and.a%rank==b%rank.and.a%nx==b%nx.and.a%ny==b%ny
   end function

   subroutine build_candidate_bank(source,nodes,allow_truncations,bank,status,degrees,independent_axes)
      type(source_t),intent(in) :: source
      integer,intent(in) :: nodes(:)
      logical,intent(in) :: allow_truncations
      type(descriptor_t),allocatable,intent(out) :: bank(:)
      integer,intent(out) :: status
      integer,optional,intent(in) :: degrees(:)
      logical,optional,intent(in) :: independent_axes
      type(descriptor_t),allocatable :: work(:)
      integer,allocatable :: ms(:)
      integer :: i,j,k,m,n,count,maximum
      logical :: axes
      status=invalid;axes=.false.
      if(present(independent_axes)) axes=independent_axes
      if(.not.source%ready.or..not.any(source%kind==[SRC_POLY,SRC_GAUSSIAN]).or.size(nodes)==0) return
      if(any(nodes<1).or.any(nodes>48)) return
      if(present(degrees)) then
         ms=degrees
      else if(source%kind==SRC_POLY) then
         m=max(source%mx,source%my)
         if(allow_truncations.and.m>1) then
            ms=[0,1,m]
         else if(allow_truncations.and.m==1) then
            ms=[0,1]
         else
            ms=[m]
         end if
      else
         ! Non-polynomial banks must supply their original degree range.
         return
      end if
      if(any(ms<0).or.any(ms>32)) return
      maximum=(4*size(ms)+2)*size(nodes)**2
      allocate(work(maximum));count=0
      do k=1,size(ms)
         m=ms(k)
         do i=1,size(nodes)
            n=nodes(i)
            call add(RT_SAS,m,-1,n,n)
            call add(RT_PROJECT_SHARED,m,-1,n,n)
            if(axes) then
               do j=1,size(nodes)
                  call add(RT_PROJECT_2D,m,-1,n,nodes(j))
                  call add(RT_POLY_FACTOR,m,-1,n,nodes(j))
               end do
            else
               call add(RT_PROJECT_2D,m,-1,n,n)
               call add(RT_POLY_FACTOR,m,-1,n,n)
            end if
         end do
      end do
      do i=1,size(nodes)
         n=nodes(i)
         call add(RT_SOURCE_SHARED,-1,-1,n,n)
         if(axes) then
            do j=1,size(nodes)
               call add(RT_SOURCE_2D,-1,-1,n,nodes(j))
            end do
         else
            call add(RT_SOURCE_2D,-1,-1,n,n)
         end if
      end do
      bank=work(:count);status=ok
   contains
      subroutine add(route,degree,rank,nx,ny)
         integer,intent(in) :: route,degree,rank,nx,ny
         type(descriptor_t) :: d
         integer :: t
         d=make_descriptor(route,degree,rank,nx,ny)
         do t=1,count
            if(same_candidate(work(t),d)) return
         end do
         count=count+1;work(count)=d
      end subroutine
   end subroutine

   subroutine prepare_candidate(source,spec,prepared,status)
      type(source_t),intent(in) :: source
      type(descriptor_t),intent(in) :: spec
      type(prepared_t),intent(out) :: prepared
      integer,intent(out) :: status
      type(source_t) :: projected,truncated
      prepared=prepared_t();prepared%spec=spec;status=invalid
      call init_rule(prepared%rx,spec%nx,status)
      if(status/=ok) return
      if(spec%nx==spec%ny) then
         prepared%ry=prepared%rx
      else
         call init_rule(prepared%ry,spec%ny,status)
         if(status/=ok) return
      end if
      select case(spec%route)
      case(RT_SAS,RT_PROJECT_2D,RT_PROJECT_SHARED,RT_POLY_FACTOR)
         call project_source(source,spec%m,projected,status)
         if(status/=ok) return
         if(spec%rank>=0) then
            call truncate_poly_rank(projected,spec%rank,truncated,status)
            if(status/=ok) return
            projected=truncated
         end if
         call compact_source(projected,prepared%model,status)
      case default
         call compact_source(source,prepared%model,status)
      end select
      if(status/=ok) return
      select case(spec%route)
      case(RT_SAS);call prepare_sas_general(prepared%model,prepared%sas,status)
      case(RT_SOURCE_2D,RT_SOURCE_SHARED,RT_PROJECT_2D,RT_PROJECT_SHARED)
         call prepare_spatial(prepared%model,prepared%space,status)
      case(RT_POLY_FACTOR);call prepare_poly_factor(prepared%model,prepared%factor,status)
      case default;status=invalid
      end select
      prepared%ready=status==ok
   end subroutine

   subroutine evaluate_candidate(prepared,queries,observable,values,status)
      type(prepared_t),intent(in) :: prepared
      real(rk),intent(in) :: queries(:,:)
      integer,intent(in) :: observable
      real(rk),intent(out) :: values(:)
      integer,intent(out) :: status
      integer :: i
      values=quiet_nan();status=invalid
      if(.not.prepared%ready.or.size(queries,1)/=4.or.size(queries,2)/=size(values)) return
      do i=1,size(values)
         select case(prepared%spec%route)
         case(RT_SAS)
            values(i)=evaluate_sas_general(prepared%sas,queries(:,i),prepared%rx,observable,status)
         case(RT_SOURCE_2D,RT_PROJECT_2D)
            values(i)=evaluate_spatial(prepared%space,queries(:,i),prepared%rx,prepared%ry,observable,status)
         case(RT_SOURCE_SHARED,RT_PROJECT_SHARED)
            values(i)=evaluate_shared(prepared%space,queries(:,i),prepared%rx,observable,status)
         case(RT_POLY_FACTOR)
            values(i)=evaluate_poly_factor(prepared%factor,queries(:,i),prepared%rx,prepared%ry,observable,status)
         case default;status=invalid
         end select
         if(status/=ok) return
      end do
      if(.not.all(ieee_is_finite(values))) status=numerical_failure
   end subroutine

   subroutine preparation_fingerprint(prepared,values,status)
      type(prepared_t),intent(in) :: prepared
      real(rk),intent(out) :: values(:)
      integer,intent(out) :: status
      ! Read the newly prepared representation, without evaluating a query.
      ! Keeping these responses separate prevents PREP being inferred by subtraction.
      values=0;status=invalid
      if(.not.prepared%ready.or.size(values)<1) return
      values(1)=real(prepared%model%mx+prepared%model%my+prepared%spec%nx+prepared%spec%ny,rk)
      if(allocated(prepared%model%c)) values(1)=values(1)+prepared%model%c(0,0)
      if(allocated(prepared%factor%px)) values(1)=values(1)+prepared%factor%px(0,1)+prepared%factor%py(0,1)
      status=ok
      if(.not.all(ieee_is_finite(values))) status=numerical_failure
   end subroutine

   subroutine candidate_errors(values,reference,EA,ER,active)
      real(rk),intent(in) :: values(:),reference(:)
      real(rk),intent(out) :: EA,ER
      logical,intent(out) :: active
      integer :: i
      EA=huge(1.0_rk);ER=EA;active=any(abs(reference)>=1e-6_rk)
      if(size(values)/=size(reference)) return
      if(.not.all(ieee_is_finite(values)).or..not.all(ieee_is_finite(reference))) return
      EA=maxval(abs(values-reference));ER=0
      do i=1,size(values)
         if(abs(reference(i))>=1e-6_rk) ER=max(ER,abs(values(i)-reference(i))/abs(reference(i)))
      end do
   end subroutine
end module
