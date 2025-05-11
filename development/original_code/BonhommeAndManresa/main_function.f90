program main
  use heuristics
  use heuristics_nocovariates
  use heuristics_bootstrap
  use heuristics_nocovariates_bootstrap
  use imsl
  use linalg
  use utils


  implicit none

  ! Parameters of the data
  integer :: T !Time periods if the panel is balanced and maxtime span if the panel is unbalanced
  integer :: N !Number of countries
  integer :: K !Number of covariates
  integer :: G !Groups
  integer :: Aux(1,2)

  ! Parameters of the algorithms
  integer :: Nsim, Nsim_bootstrap !Number of simulations
  integer :: Neighbourmax, Neighbourmax_bootstrap !Number of maximum realocations
  integer :: step_max, step_max_bootstrap
  
  ! Data matrices
  integer :: flag_VNS, flag_unbalanced,flag_nocovariates
  integer :: flag_bootstrap, flag_VNS_bootstrap, Iter_bootstrap, dd
 

  open(unit = 5001, file = 'InputNT.txt')
  read(5001,*) Aux
  close(5001)

  N = Aux(1,1)
  T = Aux(1,2)



  WRITE(*,*) '------ 1. DATA & GROUPS------'
  WRITE(*,*)'1.1.TYPE THE NUMBER OF GROUPS'
  READ(*,*) G

  WRITE(*,*)'1.2.TYPE THE NUMBER OF COVARIATES'
  READ(*,*) K

  WRITE(*,*) '------ 2. HEURISTICS -------'
  
  WRITE(*,*)'2.1.TYPE 0 FOR ALGORITHM 1 (ITERATIVE) AND 1 FOR ALGORITHM 2 (VNS)'
  READ(*,*) flag_VNS

  WRITE(*,*)'2.2.Type the number of simulations Nsim'
  READ(*,*) Nsim

  If(flag_VNS == 1)then
	
	WRITE(*,*)'2.2.1.Type the number of Neighbours'
	READ(*,*) Neighbourmax
	
	WRITE(*,*)'2.2.2.Type the number of steps (time-max)'
	READ(*,*) step_max

  end if

  WRITE(*,*) '------ 3. STANDARD ERRORS -------'
 
  WRITE(*,*)'3.1.TYPE 1 FOR BOOTSTRAP STD. ERRORS CALCULATION'
  READ(*,*) flag_bootstrap

  if (flag_bootstrap == 1) then
	
	WRITE(*,*)'3.2.TYPE 0 FOR ALGORITHM 1 (ITERATIVE) AND 1 FOR ALGORITHM 2 (VNS)'
	READ(*,*) flag_VNS_bootstrap

	WRITE(*,*)'3.3.Type the number of simulations Nsim'
    READ(*,*) Nsim_bootstrap

  If(flag_VNS_bootstrap == 1)then
	
	WRITE(*,*)'3.3.1.Type the number of Neighbours'
	READ(*,*) Neighbourmax_bootstrap
	
	WRITE(*,*)'3.3.2.Type the number of steps (time-max)'
	READ(*,*) step_max_bootstrap

  end if

	WRITE(*,*)'3.3.TYPE THE NUMBER OF BOOTSTRAP REPLICATIONS'
	READ(*,*) Iter_bootstrap
  
  end if

  
  if (K>0) then
	call heuristics_complete(N,T,K,G,Nsim,Neighbourmax,step_max,flag_VNS)
  else
	call heuristics_complete_nocov(N,T,G,Nsim,Neighbourmax,step_max,flag_VNS)
  end if

  if (flag_bootstrap == 1) then
		if (K>0) then 
			call heuristics_complete_bootstrap(N,T,K,G,Nsim_bootstrap,Neighbourmax_bootstrap,step_max_bootstrap,flag_VNS_bootstrap,Iter_bootstrap)
		else
			call heuristics_complete_nocov_bootstrap(N,T,G,Nsim_bootstrap,Neighbourmax_bootstrap,step_max_bootstrap,flag_VNS_bootstrap,Iter_bootstrap)
		end if
  end if

end program