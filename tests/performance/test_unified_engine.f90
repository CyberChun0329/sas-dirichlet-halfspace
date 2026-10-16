program test_unified_engine
   use replay_types
   use unified_engine
   use unified_policy, only: policy_t
   use paired_measure, only: pair_order,measure_paired
   use replay_measure, only: measurement_t
   implicit none
   type(source_t) :: source,small
   type(prepared_t) :: a,b
   type(descriptor_t) :: da,db
   type(descriptor_t),allocatable :: bank(:)
   type(policy_t) :: p
   type(measurement_t) :: measured(2)
   real(rk) :: queries(4,3),va(3),vb(3),coeff(0:2,0:4)
   integer :: st,statuses(2),i,calls
   queries(:,1)=[.1_rk,.2_rk,.3_rk,.75_rk]
   queries(:,2)=[.9_rk,.8_rk,.5_rk,1.0_rk]
   queries(:,3)=[1.2_rk,-.1_rk,2.0_rk,1000.0_rk]
   call profile_source('N1',1.0_rk,2.0_rk,source,st);call require(st==ok,'N1')
   call compact_source(source,small,st)
   call require(st==ok.and.all(shape(small%c)==[2,2]),'minimal rectangular N1')
   da=make_descriptor(RT_POLY_FACTOR,8,-1,12,16)
   call prepare_candidate(source,da,a,st)
   call require(st==ok.and.a%factor%rank==2.and.all(shape(a%model%c)==[2,2]),'padded alias compact before SVD')
   coeff=0;coeff(0,0)=.7_rk;coeff(1,1)=.2_rk;coeff(2,4)=-.13_rk;coeff(0,3)=.1_rk
   call init_poly(source,coeff,1.0_rk,2.0_rk,1.0_rk,st);call require(st==ok,'signed polynomial')
   da=make_descriptor(RT_PROJECT_2D,4,-1,12,16);db=make_descriptor(RT_POLY_FACTOR,4,-1,12,16)
   call prepare_candidate(source,da,a,st);call require(st==ok,'space preparation')
   call prepare_candidate(source,db,b,st);call require(st==ok.and.b%factor%rank==3,'complete rectangular thin SVD')
   call evaluate_candidate(a,queries,OBS_T,va,st);call require(st==ok,'space query')
   call evaluate_candidate(b,queries,OBS_T,vb,st);call require(st==ok,'factor query')
   call require(maxval(abs(va-vb))<2e-12_rk,'factor Horner temperature fidelity')
   call evaluate_candidate(a,queries,OBS_D,va,st);call require(st==ok,'space deficit')
   call evaluate_candidate(b,queries,OBS_D,vb,st);call require(st==ok,'factor deficit')
   call require(maxval(abs(va-vb))<2e-11_rk,'factor Horner deficit fidelity')
   db=make_descriptor(RT_SAS,4,0,12,12)
   call prepare_candidate(source,db,b,st);call require(st==ok,'rank-zero preparation')
   call evaluate_candidate(b,queries,OBS_T,vb,st)
   call require(st==ok.and.all(vb==0),'explicit rank truncation actually consumed')
   call profile_source('Q',1.0_rk,2.0_rk,source,st)
   da=make_descriptor(RT_SAS,0,-1,12,12);db=make_descriptor(RT_SAS,2,-1,12,12)
   call prepare_candidate(source,da,a,st);call require(st==ok.and.size(a%model%c)==1,'M0 true truncation')
   call prepare_candidate(source,db,b,st);call require(st==ok.and.size(b%model%c)==9,'Q full degree')
   call evaluate_candidate(a,queries,OBS_T,va,st)
   call evaluate_candidate(b,queries,OBS_T,vb,st)
   call require(maxval(abs(va-vb))>1e-5_rk,'M0 not relabelled full Q')
   call build_candidate_bank(source,[2,4],.false.,bank,st,degrees=[0,1,2],independent_axes=.true.)
   call require(st==ok.and.all(bank%route>=1).and.all(bank%route<=6),'polynomial routes retained')
   call require(any(bank%route==RT_POLY_FACTOR.and.bank%nx/=bank%ny),'independent axes preserved')
   call require(all(pair_order(1,0)==[1,2]).and.all(pair_order(2,0)==[2,1]),'alternating window schedule')
   calls=0;p%formal_allowed=.false.
   call measure_paired(never_called,never_called,va,vb,measured,statuses,p,.true.,0)
   call require(calls==0.and.all(statuses/=ok),'default formal gate before callback or clock')
   print '(a)','PASS unified engine: finite math fixtures, source/rank semantics and non-executing formal gate'
contains
   subroutine require(test,label)
      logical,intent(in) :: test
      character(*),intent(in) :: label
      if(.not.test) then
         print *,trim(label);error stop 'engine test failed'
      end if
   end subroutine
   subroutine never_called(values,status)
      real(rk),intent(out) :: values(:)
      integer,intent(out) :: status
      calls=calls+1;values=0;status=ok
   end subroutine
end program
