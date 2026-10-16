program test_case_bank
   use replay_types, only: rk,ok,OBS_D
   use replay_case_bank
   implicit none
   type(case_spec_t),allocatable :: cases(:)
   integer,allocatable :: nodes(:),nx(:),ny(:)
   integer :: status,i,audits

   call family_count('temperature_scale',128,384,3)
   call family_count('paper',100,900,9)
   call family_count('deficit_original',8,106,-1)

   call build_cases('temperature_scale',cases,status)
   call require(status==ok,'temperature_scale build status')
   call require(cases(1)%profile=='U','temperature_scale first profile')
   call require(cases(1)%source%a==1.0_rk .and. cases(1)%source%b==2.0_rk,'temperature_scale first geometry')
   call require(cases(1)%Fo==.001_rk .and. abs(cases(1)%q(3,1)-.1_rk)<1e-15_rk,'temperature_scale first Fo/depth')
   call require(all(cases(1)%q(1:2,1)==[0.0_rk,0.0_rk]),'temperature_scale center query')
   call require(cases(128)%profile=='Q' .and. cases(128)%source%b==.75_rk,'temperature_scale last group')

   call build_cases('paper',cases,status)
   audits=count([(cases(i)%audit,i=1,size(cases))])
   call require(audits==20,'paper rho=.025 audit membership')
   call require(cases(1)%bank==3 .and. cases(1)%audit,'paper base bank and audit flag')
   call require(all(cases(1)%q(1:2,2)==[.9_rk,1.8_rk]),'paper near-corner point')
   call require(all(cases(1)%q(1:2,6)==[1.2_rk,0.0_rk]),'paper axial outside point')

   call build_cases('deficit_original',cases,status)
   call require(cases(1)%profile=='U' .and. cases(1)%observable==OBS_D,'original deficit first profile/observable')
   call require(size(cases(1)%q,2)==50,'original deficit 50-depth case')
   call require(abs(cases(1)%q(3,1)-.1_rk)<1e-15_rk .and. abs(cases(1)%q(3,50)-1.5_rk)<1e-15_rk, &
      'original deficit depth endpoints')
   call require(cases(2)%q(1,1)==.5_rk .and. cases(3)%q(1,1)==.9_rk .and. cases(4)%q(1,1)==2.0_rk, &
      'original deficit U offaxis/edge/outside positions')
   call require(cases(2)%q(3,1)==.05_rk .and. cases(2)%q(3,8)==10.0_rk, &
      'original deficit extension depths')
   call require(cases(7)%profile=='Q' .and. cases(7)%Fo==100.0_rk .and. &
      cases(8)%profile=='Q' .and. cases(8)%Fo==1.0_rk,'original deficit Q time cases')
   call require(cases(1)%rho==0,'mixed depth metadata sentinel')
   call build_cases('temperature_scale',cases,status)
   call node_list(cases(1),nodes)
   call require(size(nodes)==31 .and. nodes(1)==2 .and. nodes(31)==32,'dense node bank')
   call build_cases('deficit_original',cases,status)
   call node_list(cases(1),nodes)
   call require(size(nodes)==5 .and. all(nodes==[2,4,8,16,32]),'sparse node bank')

   call build_cases('paper',cases,status)
   call pair_list(cases(1),nx,ny)
   call require(size(nx)==121 .and. size(ny)==121,'rectangular base pair count')
   call require(unique_pairs(nx,ny),'rectangular base pairs unique')
   call require(has_pair(4,6,nx,ny) .and. has_pair(32,28,nx,ny),'rectangular base representative pairs')
   cases(1)%bank=4
   call pair_list(cases(1),nx,ny)
   call require(size(nx)==229 .and. size(ny)==229,'rectangular audit pair count')
   call require(unique_pairs(nx,ny),'rectangular audit pairs unique')
   call require(has_pair(48,44,nx,ny) .and. has_pair(36,40,nx,ny),'rectangular audit representative pairs')
   call require(.not.has_pair(46,48,nx,ny),'audit bank excludes non-grid off-diagonal extras')

   call build_cases('not_a_family',cases,status)
   call require(status/=ok .and. .not.allocated(cases),'unknown family rejected')
   write(*,'(A)') 'test_case_bank: PASS'
contains
   subroutine family_count(name,expected_cases,expected_points,points_per_case)
      character(len=*),intent(in) :: name
      integer,intent(in) :: expected_cases,expected_points,points_per_case
      integer :: st,j,npoints
      call build_cases(name,cases,st)
      call require(st==ok,'build status for '//trim(name))
      call require(size(cases)==expected_cases,'case count for '//trim(name))
      npoints=0
      do j=1,size(cases)
         npoints=npoints+size(cases(j)%q,2)
         if(points_per_case>0) call require(size(cases(j)%q,2)==points_per_case, &
            'point count per case for '//trim(name))
      end do
      call require(npoints==expected_points,'total point count for '//trim(name))
   end subroutine

   subroutine require(condition,label)
      logical,intent(in) :: condition
      character(len=*),intent(in) :: label
      if(.not.condition) then
         write(*,'(A)') 'FAIL: '//trim(label)
         error stop 1
      end if
   end subroutine

   logical function has_pair(x,y,xx,yy)
      integer,intent(in) :: x,y,xx(:),yy(:)
      integer :: j
      has_pair=.false.
      do j=1,size(xx)
         if(xx(j)==x .and. yy(j)==y) then
            has_pair=.true.;return
         end if
      end do
   end function

   logical function unique_pairs(xx,yy)
      integer,intent(in) :: xx(:),yy(:)
      integer :: j,k
      unique_pairs=.false.
      do j=1,size(xx)
         do k=j+1,size(xx)
            if(xx(j)==xx(k) .and. yy(j)==yy(k)) return
         end do
      end do
      unique_pairs=.true.
   end function
end program test_case_bank
