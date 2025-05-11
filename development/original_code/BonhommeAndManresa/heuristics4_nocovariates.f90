module heuristics_nocovariates
  implicit none
contains

!********************* HEURISTICS *************************************!
subroutine heuristics_complete_nocov(N,T,G,Nsim,Neighbourmax,step_max,flag_VNS)
!  use global_module
  use imsl
  use linalg
  use utils
  
  implicit none

  integer, intent(in) :: N,T,G,Nsim,Neighbourmax,step_max,flag_VNS
  
  real(8) :: Y(N*T,1) 
  real(8) :: D(N,T), auxdebug

  ! Parameters of the model
  real(8) :: deltapar_new
  real(8) :: agt_vect(G*T,1), agt_vect0(G*T,1)
  real(8) :: agt_tren(T,G), agt_tren0(T,G),agt_tren_opt(T,G), agt_tren_empty(T,G)
  real(8) :: countg(1,G),countg0(1,G), countg_opt(1,G),countg_prime(1,G), countgaux(1,G)
  real(8) :: countgT(T,G), countgT0(T,G),countgT_opt(T,G),countgT_prime(T,G), countgauxT(T,G)
  real(8) :: Obj, Obj0, Obj00, Objh0, Objk, Objk0,Obj_prime, Obj_matrix(Nsim),Obj_max

  ! Parameters of the heuristics
  ! Parameters of K-MEANS
  ! Parameters of classification
  integer :: i, tt, ii, j, control1, control2,nn, step, loopstep, flag, nolocal
  real(8) :: yaux(T,1)
  real(8) :: agt_aux(T,1)
  real(8) :: aux_u, auxaux_u 
  integer :: gi_num(N), gi_num0(N),gi_num_opt(N),gi_final(N), gi_num_prime(N),gi_final_matrix(N,Nsim), neighbour,flagg(G) ! contains the group at wich i belongs
  real(8) :: gi(N,G) !contains a 1 in the column gg where i is assigned
  real(8) :: gi_class(N),gi_classempty(N),gi_class_opt(N),gi_class_prime(N) !contains the distances to the centroid associated with classification gi
  
  ! Parameters of the local search
  integer :: phi, phi_max

  !Parameters of updating
  real(8) :: ygt(T,G), gi_classaux(N,G)
  real(8) :: Y_dm(N*T,1), aux_vec(G)


  ! Auxiliary parameters
  real(8) :: tmp2(G), tmp3(Neighbourmax), tmp4(Neighbourmax), deltapar
  real(8) :: V(G), W(N*T,1)
  integer :: gg, dd, jsim, ggg, ggi,ggaux,flagempty,flag_shorttrend
  integer :: time_array_0(8), time_array_1(8)
  real(8) :: start_time, finish_time, aux
  real(8) :: Auxmat(N*T,1)


  !To write down the output
  open(unit=13,file='assignment.txt')
  open(unit=14,file='outputobj.txt')

  !write(*,*) 'Type 1 if you want to do VNS, type 0 if only k-means'
  !read(*,*), flag_VNS

!call date_and_time(values=time_array_0)
!  start_time = time_array_0 (5) * 3600 + time_array_0 (6) * 60 &
!           + time_array_0 (7) + 0.001 * time_array_0 (8)

  
  ! We read the data 
  !open(unit=11,file='ldem_linc_unbalanced.txt')
  !open(unit=12,file='dem_unbalanced.txt')
  !open(unit=111,file='Ti_final.txt')
  
  !open(unit=11,file='inputX.txt')
  !open(unit=12,file='inputY.txt')
  open(unit = 5000, file = 'data.txt')
  open(unit=111,file='Ti_unbalanced.txt')

	!do dd = 1,N*T
	!	read(11,*) X(dd,:)
	!end do

	!read(12,*) Y
	
	do dd = 1,N*T
		read(5000,*) Auxmat(dd,:)
	end do
	
	do dd = 1,N*T
		Y(dd,1) = Auxmat(dd,1)
	end do

	do dd = 1,N
		read(111,*) D(dd,:)
	end do

  !close(11)
  !close(12)
  close(5000)
  close(111)


  Obj00 = 100000.d0
  Obj_max = 10000.d0
  
  ! Kmeans loop:
	do jsim = 1,Nsim
200   Obj00 = 100000.d0

	print*, jsim
	
	
	!Initialization:
	
	!We pick up an initial condition for beta

	!call drnnor(K,tmp1)
    !theta0 = reshape(tmp1, (/ K, 1 /))

	!We compute the residuals and select randomly G "centers"
	!W = Y - matmul(X,theta0)
	W = Y
	call drnun(G,tmp2)
	do gg = 1,G
		V(gg) = max(nint(tmp2(gg)*N),1)
		agt_tren0(:,gg) = W((V(gg)-1)*T+1:V(gg)*T,1)
	end do

	!Iterations
	loopstep = 1
	deltapar = 1.0d0
	do while (deltapar > 0)
		!print*, deltapar
		! Step 1: assignment
		do i = 1,N
			do tt = 1,T
				yaux(tt,1) = Y((i-1)*T+tt,1)
				!xaux(tt,:) = X((i-1)*T+tt,:)
			end do
			do gg = 1,G
				aux_u = 0.0d0
				auxaux_u = 0.0d0
				do tt = 1,T
					!do kk = 1,K
					!	auxaux_u = auxaux_u + xaux(tt,kk)*theta0(kk,1)
					!end do
					aux_u = aux_u + D(i,tt)*(yaux(tt,1) - auxaux_u - agt_tren0(tt,gg))**2
					auxaux_u = 0.0d0
				end do ! tt
				gi_classaux(i,gg) = aux_u
			end do !gg
		end do !i

		
		! Classification
		gi_num = minloc(gi_classaux,2)
		gi_class = minval(gi_classaux,2)
		
		! Initializing the matrices for counting the number of elements in each group and time
		! observations: countg and countgT
		do gg = 1,G
			countg(1,gg) = 0.d0
			flagg(gg) = 0
			do tt = 1,T
				countgT(tt,gg) = 0.d0
			end do
		end do

		! Counting the number of observations in each group and time period
		do gg = 1,G
			do i = 1,N
				if (gi_num(i) == gg) then
					flagg(gg) = 1
					countg(1,gg) = countg(1,gg)+1
					do tt = 1,T
						countgT(tt,gg) = countgT(tt,gg) + D(i,tt)
					end do
				end if
			end do
		end do !i
				

		!We reassign furthest individuals to empty groups
		!countgaux = countg
		
		do gg = 1,G
			if (countg(1,gg) == 0) then
			flagempty = 0
			do while (flagempty==0)
				ii = maxloc(abs(gi_class),1)
				ggaux = gi_num(ii)
				if (countg(1,ggaux) == 1.d0) then
					!gi_classempty(ii) = 0.0d0
					gi_class(ii) = 0.0d0
				else
					flagempty = 1
					gi_num(ii) = gg
					countg(1,gg) = 1.d0
					countg(1,ggaux) = countg(1,ggaux)-1.0d0
					do tt = 1,T
						countgT(tt,gg) = D(ii,tt)
						countgT(tt,ggaux) = countgT(tt,ggaux)-D(ii,tt)
					end do
					gi_class(ii) = 0.0d0
				end if
			end do
			end if
		end do

		! The initial conditions has no short trends
		!do gg = 1,G
		!	do tt = 1,T
		!		if (countgT(tt,gg) == 0.0d0) then
		!			!go to 200
		!			print*, 'sort trend'
		!			flag_shorttrend
		!		end if
		!	end do
		!end do

		if(minval(countg(1,1:G),1)<1.0d0) then
			!print*, 'kk'
			Obj_matrix(jsim) = 10000.d0
			go to 100
		end if

		do tt = 1,T
			do gg = 1,G
				if(countgT(tt,gg)<1.0d0) then
					!print*, 'short trend'
					flag_shorttrend = 1
					!go to 200
				end if
			end do
		end do

		! Step 2: update

		! We compute the within group means
		! We initialize to zero:
		do gg = 1,G
			do tt = 1,T
		!		do kk = 1,K
		!			xgt(tt,(kk-1)*G+gg) = 0.0d0
		!		end do
				ygt(tt,gg) = 0.0d0
			end do
		end do
		
		do i = 1,N
			do gg = 1,G
				if (gi_num(i) == gg) then
					do tt = 1,T
						if(countgT(tt,gg)>=1.0d0) then
							!do kk = 1,K
							!	xgt(tt,(kk-1)*G+gg) = xgt(tt,(kk-1)*G+gg) + X((i-1)*T+tt,kk)/countgT(tt,gg)
							!end do
							ygt(tt,gg) = ygt(tt,gg) + Y((i-1)*T+tt,1)/countgT(tt,gg)
						else
							!do kk = 1,K
							!	xgt(tt,(kk-1)*G+gg) = 0.0d0
							!end do
							ygt(tt,gg) = 0.0d0
						end if
					end do
				end if
			end do
		end do

		! We compute the demeaned vectors:
		do i = 1,N
			do gg = 1,G
				if (gi_num(i) == gg) then
					do tt = 1,T
						!do kk = 1,K
						!	X_dm((i-1)*T+tt,kk) = D(i,tt)*(X((i-1)*T+tt,kk) - xgt(tt,(kk-1)*G+gg))
						!end do
						Y_dm((i-1)*T+tt,1) = D(i,tt)*(Y((i-1)*T+tt,1) - ygt(tt,gg))
					end do !tt
				end if
			end do !gg
		end do !i
		
		!mat1 = matmul(transpose(X_dm),X_dm)
		!print*, mat1

		!call lu_inverse(mat1,imat1)
		!theta = matmul(matmul(imat1,transpose(X_dm)),Y_dm)

		Obj = 0.0d0
		do i = 1,N
			do tt = 1,T
				auxaux_u = 0.d0
				!do kk = 1,K
				!	auxaux_u = auxaux_u + X_dm((i-1)*T+tt,kk)*theta(kk,1)
				!end do
				Obj = Obj + D(i,tt)*(Y_dm((i-1)*T+tt,1) - auxaux_u)**2
			end do
		end do
		!print*, Obj

		!write(13,*) Obj, theta, ' jsim'
		
		! We compute the time-trends
		! We initialize them first:
		do gg = 1,G
			do tt = 1,T
				agt_tren(tt,gg) = 0.d0
			end do
		end do

		do i = 1,N
			do gg = 1,G
				if(gi_num(i) == gg) then
					do tt = 1,T
						if (countgT(tt,gg) >= 1.0d0) then
							auxaux_u = 0.d0
							!do kk = 1,K
							!	auxaux_u = auxaux_u + X((i-1)*T+tt,kk)*theta(kk,1)
							!end do
							agt_tren(tt,gg) = agt_tren(tt,gg) + (Y((i-1)*T+tt,1) - auxaux_u)/countgT(tt,gg)
						else
							agt_tren(tt,gg) = -10.000d0
						end if
					end do
				end if
			end do
		end do

		if (flag_shorttrend == 0) then
			deltapar_new = 0.0d0
			!do kk = 1,K
			!	deltapar_new = deltapar_new + (theta(kk,1)-theta0(kk,1))**2
			!end do
			do gg = 1,G
				do tt = 1,T
					deltapar_new = deltapar_new + (agt_tren(tt,gg) - agt_tren0(tt,gg))**2
				end do
			end do
			!auxdebug = minval(countgT)
			!print*, auxdebug

			deltapar = deltapar_new
			!print*, deltapar

			loopstep  = loopstep + 1
			if (loopstep > 10000) then
				print*, deltapar
			!	call drnnor(K,tmp1)
			!	theta = theta + reshape(0.1*tmp1, (/ K, 1 /))
			!	deltapar = 1
			!	loopstep = 1
				!agt_tren(tt,gg) = agt_tren(tt,gg) + ??
			end if
		else
			deltapar = -10.0d0
			!theta = theta0
			agt_tren = agt_tren0
		end if ! end if flagshort_trend == 0

		!theta0 = theta
		agt_tren0 = agt_tren


		!print*, Obj
		!print*, deltapar
				
		end do !deltapar
		
		if ((Obj < Obj00).or.(flag_shorttrend == 1)) then
			!theta_opt = theta
			agt_tren_opt = agt_tren
			Obj00 = Obj
			!print*, Obj00
			!print*, theta_opt
			gi_num_opt = gi_num
			countg_opt = countg
			countgT_opt = countgT
			flag_shorttrend = 0
		end if

	if (flag_VNS == 1) then
	!---------------------------Beginning of VNS--------------------------------------------!

	! Step 0: read any initial K-cluster partition:

	!theta0 = theta_opt
	agt_tren0 = agt_tren_opt
	Obj = Obj00
	
	gi_num = gi_num_opt
	countg = countg_opt
	countgT = countgT_opt
	step = 0
	
	control1 = 1 !WRITE WHAT IS THIS??
	do while (control1 == 1)
	
	!print*, 'VNS begins'

	!Step 1: Set neighbour = 1
	neighbour = 1

	!Step 2: Relocate neighbor randomly selected objets to randomly selected new clusters
	! to create a new partition gi_num_prime:
	
	control2 = 1 !WRITE WHAT IS THIS CONTROL??
	do while (control2 == 1)
	
	!print*, 'neighbour', neighbour , '        step', step
	!print*, Objk

	! We randomize the individuals and clusters
	call drnun(Neighbourmax,tmp3)
	call drnun(Neighbourmax,tmp4)
	
	gi_num_prime = gi_num
	countg_prime = countg
	countgT_prime = countgT

	do nn = 1,neighbour
		ii = min(max(nint(tmp3(nn)*N),1),N)
		ggaux = gi_num_prime(ii)
		! we choose another individual if it belongs to a cluster where she is alone:
		do while (countg_prime(1,ggaux)==1)
		!print*, 'loop1'
			call drnun(Neighbourmax,tmp3)
			ii = min(max(nint(tmp3(nn)*N),1),N)
			ggaux = gi_num_prime(ii)
		end do
		ggi = min(max(nint(tmp4(nn)*G),1),G)
		! if the random cluster is the one to which ii belongs or has only one element,
		! then, keep on drawing pontential clusters.
		do while(ggaux == ggi)
		!print*, 'loop2'
		!print*, 'loop there is too few elements in each cluster'
			call drnun(1,aux)
			ggi = max(nint(aux*G),1)
		end do
		gi_num_prime(ii) = ggi
		countg_prime(1,ggi) = countg_prime(1,ggi) + 1.0d0
		countg_prime(1,ggaux) = countg_prime(1,ggaux)-1.0d0
		do tt = 1,T
				countgT_prime(tt,ggi) = countgT_prime(tt,ggi) + D(ii,tt)
				countgT_prime(tt,ggaux) = countgT_prime(tt,ggaux) - D(ii,tt)
		end do
	end do

	if (minval(countg_prime(1,1:G),1) < 1.0d0) then
		!print*, 'countg_prime is 0'
		Obj_matrix(jsim) = 10000.d0
		go to 100
	end if

	!Step 3: Apply HK-means unsing gi_num_prime as the starting partition
	!print*, 'HK-means'
	deltapar = 1.0d0

	!do while (deltapar > 0 .and. loopstep < 10000)
	
	flag = 0
	loopstep = 1
	do while (deltapar > 0.0d0)
	!print*, deltapar
	

	if (flag == 1) then !IE this is not the first time I enter this loop
		! Step 1: assignment
		do i = 1,N
			do tt = 1,T
				yaux(tt,1) = Y((i-1)*T+tt,1)
				!xaux(tt,:) = X((i-1)*T+tt,:)
			end do
			do gg = 1,G
				aux_u = 0.0d0
				auxaux_u = 0.0d0
				do tt = 1,T
					!do kk = 1,K
					!	auxaux_u = auxaux_u + xaux(tt,kk)*theta0(kk,1)
					!end do
					aux_u = aux_u + D(i,tt)*(yaux(tt,1) - auxaux_u - agt_tren0(tt,gg))**2
					auxaux_u = 0.0d0
				end do ! tt
				gi_classaux(i,gg) = aux_u
			end do !gg
		end do !i

		! Classification
		gi_num_prime = minloc(gi_classaux,2)
		gi_class_prime = minval(gi_classaux,2)

	! Counting number of elements in each group
		do gg = 1,G
			countg_prime(1,gg) = 0.0d0
			do tt = 1,T
				countgT_prime(tt,gg) = 0.0d0
			end do
		end do

		! Counting the number of observations in each group and time period
		do gg = 1,G
			do i = 1,N
				if (gi_num_prime(i) == gg) then
					flagg(gg) = 1
					countg_prime(1,gg) = countg_prime(1,gg)+1
					do tt = 1,T
						countgT_prime(tt,gg) = countgT_prime(tt,gg) + D(i,tt)
					end do
				end if
			end do
		end do !i
		
		
		!We reassign furthest individuals to empty groups 
		do gg = 1,G
			if (countg_prime(1,gg) == 0.0d0) then
			flagempty = 0
				do while (flagempty == 0)
				!print*, 'flagempty'
					ii = maxloc(abs(gi_class_prime),1)
					ggaux = gi_num_prime(ii)
					if (countg_prime(1,ggaux) == 1.0d0) then
						gi_class_prime(ii) = 0.0d0
					else
						flagempty = 1
						gi_num_prime(ii) = gg
						countg_prime(1,gg) = 1.0d0
						countg_prime(1,ggaux) = countg_prime(1,ggaux)-1.0d0
						gi_class_prime(ii) = 0.0d0
						do tt = 1,T
							countgT_prime(tt,gg) = D(ii,tt)
							countgT_prime(tt,ggaux) = countgT_prime(tt,ggaux) - D(ii,tt)
						end do
						gi_class(ii) = 0.0d0
					end if
				end do
			end if
		end do

		end if !it's not the first round

		! Step 2: Update

		
		! We compute the within group means
		! We initialize to zero:

		flag = 1
		
		if (minval(countg_prime(1,1:G),1) == 0.0d0) then
			!print*, 'countg_prime is 0'
			Obj_matrix(jsim) = 10000.d0
			go to 100
		end if

		do gg = 1,G
			do tt = 1,T
				!do kk = 1,K
				!	xgt(tt,(kk-1)*G+gg) = 0.0d0
				!end do
				ygt(tt,gg) = 0.0d0
			end do
		end do

		do i = 1,N
			do gg = 1,G
				if (gi_num_prime(i) == gg) then
					do tt = 1,T
						if (countgT_prime(tt,gg) >= 1.0d0) then
							!do kk = 1,K
							!	xgt(tt,(kk-1)*G+gg) = xgt(tt,(kk-1)*G+gg) + X((i-1)*T+tt,kk)/countgT_prime(tt,gg)
							!end do
							ygt(tt,gg) = ygt(tt,gg) + Y((i-1)*T+tt,1)/countgT_prime(tt,gg)
						else
							flag_shorttrend = 1;
							!do kk = 1,K
							!	xgt(tt,(kk-1)*G+gg) = 0.0d0
							!end do
							ygt(tt,gg) = 0.0d0
						end if
					end do
				end if
			end do
		end do

		!aux = sum(countgT_prime)
		!print*, aux

		! We compute the demeaned vectors:
		do i = 1,N
			do gg = 1,G
				if (gi_num_prime(i) == gg) then
					do tt = 1,T
						!do kk = 1,K
						!	X_dm((i-1)*T+tt,kk) = D(i,tt)*(X((i-1)*T+tt,kk) - xgt(tt,(kk-1)*G+gg))
						!end do
						Y_dm((i-1)*T+tt,1) = D(i,tt)*(Y((i-1)*T+tt,1) - ygt(tt,gg))
					end do !tt
				end if
			end do !gg
		end do !i
		
		!mat1 = matmul(transpose(X_dm),X_dm)
		!print*, mat1

		!call lu_inverse(mat1,imat1)
		!theta = matmul(matmul(imat1,transpose(X_dm)),Y_dm)

		Objh0 = 0.0d0
		do i = 1,N
			do tt = 1,T
				auxaux_u = 0.d0
				!do kk = 1,K
				!	auxaux_u = auxaux_u + X_dm((i-1)*T+tt,kk)*theta(kk,1)
				!end do
				Objh0 = Objh0 + D(i,tt)*(Y_dm((i-1)*T+tt,1) - auxaux_u)**2
			end do
		end do
		!print*, Objh0
		
		!write(13,*), Objh0, theta,' N'
		
		! We compute the time-trends
		! We initialize them first:
		do gg = 1,G
			do tt = 1,T
				agt_tren(tt,gg) = 0.d0
			end do
		end do
		
		do i = 1,N
			do gg = 1,G
				if(gi_num_prime(i) == gg) then
					do tt = 1,T
					if (countgT_prime(tt,gg) >= 1.0d0) then
					auxaux_u = 0.d0
						!do kk = 1,K
						!	auxaux_u = auxaux_u + X((i-1)*T+tt,kk)*theta(kk,1)
						!end do
						agt_tren(tt,gg) = agt_tren(tt,gg) + (Y((i-1)*T+tt,1) - auxaux_u)/countgT_prime(tt,gg)
					else
						agt_tren(tt,gg) = -10.0d0
					end if
					end do
				end if
			end do
		end do

		if (flag_shorttrend == 0) then

		deltapar_new = 0.0d0
		!do kk = 1,K
		!	deltapar_new = deltapar_new + (theta(kk,1)-theta0(kk,1))**2
		!end do
		do gg = 1,G
			do tt = 1,T
				!deltapar_new = deltapar_new + countgT_prime(tt,gg)*(agt_tren(tt,gg) - agt_tren0(tt,gg))**2
				deltapar_new = deltapar_new + (agt_tren(tt,gg) - agt_tren0(tt,gg))**2
				!deltapar_new = deltapar_new + (agt_tren(tt,gg) - agt_tren0(tt,gg))**2
			end do
		end do

		deltapar = deltapar_new
		
		loopstep  = loopstep + 1
		!print*, deltapar
		if (loopstep > 1000) then
			!print*, deltapar
			!call drnnor(K,tmp1)
			!theta0 = theta0 + reshape(0.1*tmp1, (/ K, 1 /))
			deltapar = 1
			loopstep = 1
			!print*, X_dm
			!print*, countgT_prime
			!print*, theta
			!print*, agt_tren
		end if
		else
			deltapar = -10.0d0
			!theta = theta0
			agt_tren = agt_tren0
			flag_shorttrend = 0
		end if ! end if flagshort_trend == 0

		
		!theta0 = theta
		agt_tren0 = agt_tren

		!print*, deltapar

		end do !deltapar

		
		! Beginning of the local heuristic search
		phi_max = 1
		Objk = Objh0
		loopstep = 1
		!if (nolocal == 0) then
		do while (phi_max == 1)
		loopstep = loopstep + 1
		phi = 0
		do i = 1,N
			! We identify the group to which i belongs
			ggi = gi_num_prime (i)
			if (countg_prime(1,ggi) > 1.d0) then
				do gg = 1,G
					if (gg.ne.ggi) then
						! We need to compute the objective function if we change i from its current 
						! group to another one,hence we need to compute theta and the means across 
						! the new groups.

						! We first modify the allocation of group
						! We change i from group ggi to group gg
						gi_num0 = gi_num_prime
						gi_num0(i) = gg
						
						! Adjust the number of elements in each group
						countg0 = countg_prime
						countgT0 = countgT_prime

						countg0(1,ggi) = countg_prime(1,ggi)-1.0d0
						countg0(1,gg) = countg_prime(1,gg)+1.0d0
						
						do tt = 1,T
							countgT0(tt,ggi) = countgT0(tt,ggi)-D(i,tt)
							countgT0(tt,gg) = countgT0(tt,gg)+D(i,tt)
						end do

						do tt = 1,T
							do ggg = 1,G
								if(countgT0(tt,ggg)<1.0d0) then
								!print*, 'short trend'
								!flag_shorttrend = 1
								end if
							end do
						end do

						! We compute the within group means
						! We initialize to zero:
						do ggg = 1,G
							do tt = 1,T
								!do kk = 1,K
								!	xgt(tt,(kk-1)*G+ggg) = 0.0d0
								!end do
								ygt(tt,ggg) = 0.0d0
							end do
						end do

						do j = 1,N
							do ggg = 1,G
								if (gi_num0(j) == ggg) then
									do tt = 1,T
										if (countgT0(tt,ggg) >= 1.0d0) then
											!do kk = 1,K
											!	xgt(tt,(kk-1)*G+ggg) = xgt(tt,(kk-1)*G+ggg) + X((j-1)*T+tt,kk)/countgT0(tt,ggg)
											!end do
											ygt(tt,ggg) = ygt(tt,ggg) + Y((j-1)*T+tt,1)/countgT0(tt,ggg)
										else
											!do kk = 1,K
											!	xgt(tt,(kk-1)*G+ggg) = 0.0d0
											!end do
											ygt(tt,ggg) = 0.0d0
										end if
									end do
								end if
							end do
						end do
						
						! We compute the demeaned vectors:
						do j = 1,N
							do ggg = 1,G
								if (gi_num0(j) == ggg) then
									do tt = 1,T
										!do kk = 1,K
										!	X_dm((j-1)*T+tt,kk) = D(j,tt)*(X((j-1)*T+tt,kk) - xgt(tt,(kk-1)*G+ggg))
										!end do
										Y_dm((j-1)*T+tt,1) = D(j,tt)*(Y((j-1)*T+tt,1) - ygt(tt,ggg))
									end do !tt
								end if
							end do !ggg
						end do !j

						!mat1 = matmul(transpose(X_dm),X_dm)
						!print*, mat1

						!call lu_inverse(mat1,imat1)
						!theta0 = matmul(matmul(imat1,transpose(X_dm)),Y_dm)
						
						Objk0 = 0.d0
						do j = 1,N
							do tt = 1,T
								auxaux_u = 0.d0
								!do kk = 1,K
								!	auxaux_u = auxaux_u + X_dm((j-1)*T+tt,kk)*theta0(kk,1)
								!end do
								Objk0 = Objk0 +  D(j,tt)*(Y_dm((j-1)*T+tt,1) - auxaux_u)**2
							end do
						end do

						!write(13,*), Objk0, theta0,' L'

						if (Objk0 < Objk) then
							!print*, 'bingo'
							gi_num_prime = gi_num0
							countg_prime = countg0
							countgT_prime = countgT0
							Objk = Objk0
							!theta = theta0
							phi = 1
						end if
					end if !gg neq ggi
				end do !! gg
			end if !! there is more than 1 element in ggi
		end do ! i

		if (phi == 0) then
			phi_max = 0
		end if
		loopstep = loopstep + 1
		if (loopstep > 10000) then
			!print*, 'phi_loop'
		end if

		end do !while phi_max
		!end if

		! Step 4: If the Objetive has improved, we reset gi_num = gi_prime (after HK) and 
		! restart fron neighbor 1
		Obj_prime = Objk
		!print*, Obj
		if (Obj_prime < Obj) then
			gi_num = gi_num_prime
			gi_num_opt = gi_num
			!theta_opt = theta
			countg = countg_prime
			countgT = countgT_prime
			Obj = Obj_prime
			control2 = 0 !---> go to Step1
		else
		! Step 5: If neighbour is not yet max, set neighbour = neighbour + 1 and go to Step 2 
			if (neighbour < Neighbourmax) then
				neighbour = neighbour + 1
			else
				control2 = 0 !---> if not, go to Step 1: with the same partition, draw all neighbours
				step = step + 1
				print*, '.'
			end if
		end if
		 
		end do !control2
		if (step > step_max) then
			control1 = 0
		end if
	end do !control1

	end if ! end of flag_VNS

	!print*, theta_opt
	!print*, Obj
	
	Obj_matrix(jsim) = Obj
	!do kk = 1,K
	!	theta_matrix(jsim,kk) = theta_opt(kk,1)
	!end do
	gi_final_matrix(1:N,jsim) = gi_num_opt
	write(14,*), 'Iteration', jsim 
	write(14,*), 'Objective function', Obj
	write(14,*), ' '
	
	!if(Obj_max>Obj) then
	!	Obj_max = Obj
	!	gi_final = gi_num
	!end if

100	end do ! jsim

	call date_and_time(values=time_array_1)
      finish_time = time_array_1 (5) * 3600 + time_array_1 (6) * 60 &
           + time_array_1 (7) + 0.001 * time_array_1 (8)

	
	!write (13, '(8x, 1a, 1f16.6)') 'elapsed wall clock time:', &
     !      finish_time - start_time
	
	nn = minloc(Obj_matrix,1)
	write(14,*) ' '
	write(14,*) '--- The optimal solution is ---'
	write(14,*) ' '
	!write(14,*), theta_matrix(nn,:)
	write(14,*) 'Objective function:', minval(Obj_matrix)
	!write(14,*), minval(Obj_matrix)
	!write(14,*), 'assignment'
	write(13,822), gi_final_matrix(1:N,nn)
	!write(13,*), N
	
	!write(14,*), 'assignment'
	
	close(13)
	close(14)

	822 FORMAT (1I2)
	
	end subroutine heuristics_complete_nocov
	
	end module heuristics_nocovariates