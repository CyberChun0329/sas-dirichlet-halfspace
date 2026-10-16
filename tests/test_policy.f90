program test_policy
   use unified_policy
   implicit none
   type(policy_t) :: p,bad
   integer :: status
   integer,allocatable :: nodes(:)
   call load_policy('config/policy.nml',p,status)
   if(status/=0.or.p%formal_allowed)error stop 'policy read/current no-formal boundary'
   call candidate_nodes(p,.false.,.false.,nodes,status)
   if(status/=0.or.size(nodes)/=8)error stop 'base nodes'
   if(any(nodes/=[2,4,6,8,12,16,24,32]))error stop 'base identity'
   call candidate_nodes(p,.true.,.true.,nodes,status)
   if(status/=0.or.size(nodes)/=13)error stop 'audit with permitted one'
   if(nodes(1)/=1.or.nodes(13)/=48)error stop 'audit endpoints'
   bad=p;bad%coarse_nodes(3)=bad%coarse_nodes(2)
   call validate_policy(bad,status)
   if(status==0)error stop 'duplicate node accepted'
   bad=p;bad%max_workers=9
   call validate_policy(bad,status)
   if(status==0)error stop 'oversubscription accepted'
   bad=p;bad%max_refinement_rounds=4
   call validate_policy(bad,status)
   if(status==0)error stop 'unbounded refinement accepted'
   bad=p;bad%formal_selection_passes=4
   call validate_policy(bad,status)
   if(status==0)error stop 'unbounded formal closure accepted'
   bad=p;bad%callback_limit_seconds=31
   call validate_policy(bad,status)
   if(status==0)error stop 'callback limit expansion accepted'
   bad=p;bad%formal_min_window=.001
   call validate_policy(bad,status)
   if(status==0)error stop 'short formal window accepted'
   print '(a)','PASS policy: fixed nodes, bounded refinement/resources, formal disabled'
end program
