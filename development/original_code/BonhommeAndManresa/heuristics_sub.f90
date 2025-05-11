!module global_module
!  implicit none

  ! Parameters of the data
!  integer, parameter :: T = 7				!Time 
!  integer, parameter :: N = 90				!Number of countries
!  integer, parameter :: K = 2				!Number of covariates
!  integer, parameter :: G = 10

  ! Parameters of the algorithms
!  integer, parameter :: Neighbourmax = 10	!Number of maximum realocations
!  integer, parameter :: step_max = 1

  ! Data matrices
!  real(8) :: X(N*T,K) !
!  real(8) :: Y(N*T,1) !
  !integer :: G
  
  !real(8), parameter :: pi = 3.14159265358979d0

!end module global_module


!********************* MAIN *******************************************!
!program main
!use global_module

!real(8) :: Obj, theta_opt(K,1), gi_num(N)
  

!call heuristics(G,N,T,K,Nsim,Neighbourmax,step_max,Obj,theta_opt,gi_num)


!end program main


!********************* SUBROUTINE *************************************!
module heuristics_sub
implicit none
contains

subroutine heuristics(G,N,T,K,Nsim,Neighbourmax,step_max,Obj,theta_opt,agt_tren_opt,gi_num_opt,gi_class)

  !use global_module
  use imsl
  use linalg
  use utils
  
  implicit none
  
  ! Number of Groups to estimate (should be loaded)
  integer, intent(in) :: G,N,T,K,Nsim,Neighbourmax,step_max
  integer, intent(out) :: gi_num_opt(N)
  real(8), intent(out) :: Obj, theta_opt(K,1),agt_tren_opt(T,G),gi_class(N)
  
  !storage of the data		
  real(8) :: X(N*T,K) !
  real(8) :: Y(N*T,1) !

  ! Parameters of the model
  real(8) :: theta(K,1), theta0(K,1),deltapar_new
  real(8) :: agt_vect(G*T,1), agt_vect0(G*T,1)
  real(8) :: agt_tren(T,G), agt_tren0(T,G)
  real(8) :: countg(1,G),countg0(1,G), countg_opt(1,G),countg_prime(1,G)
  real(8) :: Obj0, Obj00, Objh0, Objk, Objk0,Obj_prime

  ! Parameters of the heuristics
  ! Parameters of K-MEANS
  ! Parameters of classification
  integer :: i, tt, kk, ii, j, control1, control2,nn, step, loopstep, flag
  real(8) :: xaux(T,K), yaux(T,1)
  real(8) :: agt_aux(T,1)
  real(8) :: aux_u, auxaux_u 
  integer :: gi_num0(N),gi_num(N), gi_num_prime(N), neighbour ! contains the group at wich i belongs
  real(8) :: gi(N,G) !contains a 1 in the column gg where i is assigned
  real(8) :: gi_class_prime(N) !contains the distances to the centroid associated with classification gi
  real(8) :: mat1(K,K), imat1(K,K)

  ! Parameters of the local search
  integer :: phi, phi_max

  !Parameters of updating
  real(8) :: x1gt(T,G), x2gt(T,G), ygt(T,G),gi_classaux(N,G)
  real(8) :: X_dm(N*T,K), Y_dm(N*T,1), aux_vec(G)


  ! Auxiliary parameters
  real(8) :: tmp1(K), tmp2(G), tmp3(Neighbourmax), tmp4(Neighbourmax), deltapar
  real(8) :: V(G), W(N*T,1)
  integer :: gg, dd, jsim, ggg, ggi
  integer :: time_array_0(8), time_array_1(8),flagflag
  real(8) :: start_time, finish_time, aux


  call date_and_time(values=time_array_0)
  start_time = time_array_0 (5) * 3600 + time_array_0 (6) * 60 &
           + time_array_0 (7) + 0.001 * time_array_0 (8)

  

  !To write down the output
  open(unit=13,file='outputtheta.txt')
  
  ! We read the data
  open(unit=11,file='inputX.txt')
  open(unit=12,file='inputY.txt')

	do dd = 1,N*T
		read(11,*) X(dd,:)
	end do

	read(12,*) Y

  close(11)
  close(12)
  Obj00 = 100000.d0

  ! Kmeans loop:
	do jsim = 1,Nsim
	print*, jsim
	
	
	!Initialization:
	
	!We pick up an initial condition for beta
	call drnnor(K,tmp1)
    theta0 = reshape(tmp1, (/ K, 1 /))

	!We compute the residuals and select randomly G "centers"
	W = Y - matmul(X,theta0)
	call drnun(G,tmp2)
	do gg = 1,G
		V(gg) = max(nint(tmp2(gg)*N),1)
		agt_tren0(:,gg) = W((V(gg)-1)*T+1:V(gg)*T,1)
	end do
	
	

	!Iterations
	deltapar = 1.0d0
	do while (deltapar > 0)
		
		! Step 1: assignment
		do i = 1,N
			do tt = 1,T
				yaux(tt,1) = Y((i-1)*T+tt,1)
				xaux(tt,:) = X((i-1)*T+tt,:)
			end do
			do gg = 1,G
				aux_u = 0.0
				auxaux_u = 0.0
				do tt = 1,T
					do kk = 1,K
						auxaux_u = auxaux_u + xaux(tt,kk)*theta0(kk,1)
					end do
					aux_u = aux_u + (yaux(tt,1) - auxaux_u - agt_tren0(tt,gg))**2
					auxaux_u = 0.0
				end do ! tt
				gi_classaux(i,gg) = aux_u
			end do !gg
		end do !i

		! Classification
		gi_num = minloc(gi_classaux,2)
		gi_class = minval(gi_classaux,2)
		
		! Counting number of elements in each group

		do gg = 1,G
			countg(1,gg) = 0
		end do
		do i = 1,N
			do gg = 1,G
				if (gi_num(i) == gg) then
					countg(1,gg) = countg(1,gg)+1
				end if
			end do
		end do !i
		
		!We check for empty groups
		do gg = 1,G
			if (countg(1,gg) == 0) then
				ii = maxloc(abs(gi_class),1)
				gi_num(ii) = gg
				countg(1,gg) = 1
				gi_class(ii) = 0.0d0
			end if
		end do
		

		! Step 2: update

		! We compute the within group means
		! We initialize to zero:
		do gg = 1,G
			do tt = 1,T
				x1gt(tt,gg) = 0.0d0
				x2gt(tt,gg) = 0.0d0
				ygt(tt,gg) = 0.0d0
			end do
		end do

		do i = 1,N
			do gg = 1,G
				if (gi_num(i) == gg) then
					do tt = 1,T
					x1gt(tt,gg) = x1gt(tt,gg) + X((i-1)*T+tt,1)/countg(1,gg)
					x2gt(tt,gg) = x2gt(tt,gg) + X((i-1)*T+tt,2)/countg(1,gg)
					ygt(tt,gg) = ygt(tt,gg) + Y((i-1)*T+tt,1)/countg(1,gg)
					end do
				end if
			end do
		end do

		! We compute the demeaned vectors:
		do i = 1,N
			do gg = 1,G
			if (gi_num(i) == gg) then
				do tt = 1,T
				X_dm((i-1)*T+tt,1) = X((i-1)*T+tt,1) - x1gt(tt,gg)
				X_dm((i-1)*T+tt,2) = X((i-1)*T+tt,2) - x2gt(tt,gg)
				Y_dm((i-1)*T+tt,1) = Y((i-1)*T+tt,1) - ygt(tt,gg)
				end do !tt
			end if
			end do !gg
		end do !i
		
		mat1 = matmul(transpose(X_dm),X_dm)
		!print*, mat1

		call lu_inverse(mat1,imat1)
		theta = matmul(matmul(imat1,transpose(X_dm)),Y_dm)

		Obj = 0.0d0
		do i = 1,N
			do tt = 1,T
				Obj = Obj + (Y_dm((i-1)*T+tt,1) - X_dm((i-1)*T+tt,1)*theta(1,1) - X_dm((i-1)*T+tt,2)*theta(2,1))**2
			end do
		end do
		
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
						agt_tren(tt,gg) = agt_tren(tt,gg) + (Y((i-1)*T+tt,1) - X((i-1)*T+tt,1)*theta(1,1) - X((i-1)*T+tt,2)*theta(2,1))/countg(1,gg)
					end do
				end if
			end do
		end do

		deltapar_new = 0.0d0
		do kk = 1,K
			deltapar_new = deltapar_new + (theta(kk,1)-theta0(kk,1))**2
		end do
		do gg = 1,G
			do tt = 1,T
				deltapar_new = deltapar_new + (agt_tren(tt,gg) - agt_tren0(tt,gg))**2
			end do
		end do

		deltapar = deltapar_new

		theta0 = theta
		agt_tren0 = agt_tren
		
		end do !deltapar
		
		if (Obj < Obj00) then
			theta_opt = theta
			agt_tren_opt = agt_tren
			Obj00 = Obj
			gi_num_opt = gi_num
			countg_opt = countg
		end if

	!---------------------------Beginning of VNS--------------------------------------------!

	! Step 0: read any initial K-cluster partition:

	theta0 = theta_opt
	gi_num = gi_num_opt
	countg = countg_opt
	agt_tren0 = agt_tren_opt
	Obj = Obj00
	step = 0
	
	control1 = 1 !WRITE WHAT IS THIS??
	do while (control1 == 1)
	
	print*, 'VNS begins'

	!Step 1: Set neighbour = 1
	neighbour = 1

	!Step 2: Relocate neighbor randomly selected objets to randomly selected new clusters
	! to create a new partition gi_num_prime:
	
	control2 = 1 !WRITE WHAT IS THIS CONTROL??
	do while (control2 == 1)
	
	print*, 'neighbour', neighbour , '        step', step

	! We randomize the individuals and clusters
	call drnun(Neighbourmax,tmp3)
	call drnun(Neighbourmax,tmp4)
	

	gi_num_prime = gi_num
	countg_prime = countg
	
	
	do nn = 1,neighbour
		ii = min(max(nint(tmp3(nn)*N),1),N)
		ggi = min(max(nint(tmp4(nn)*G),1),G)
		do while (gi_num(ii) == ggi)
			call drnun(1,aux)
			ggi = max(nint(aux*G),1)
		end do
		gi_num_prime(ii) = ggi
		countg_prime(1,ggi) = countg(1,ggi) + 1
		countg_prime(1,gi_num(ii)) = countg(1,gi_num(ii))-1
	end do

	!Step 3: Apply HK-means unsing gi_num_prime as the starting soluiton
	print*, 'HK-means'
	deltapar = 1.0d0
	!do while (deltapar > 0 .and. loopstep < 10000)
	
	flag = 0
	do while (deltapar > 0)
	if (flag == 1) then !IE this is not the first time I enter this loop
		! Step 1: assignment
		do i = 1,N
			do tt = 1,T
				yaux(tt,1) = Y((i-1)*T+tt,1)
				xaux(tt,:) = X((i-1)*T+tt,:)
			end do
			do gg = 1,G
				aux_u = 0.0
				auxaux_u = 0.0
				do tt = 1,T
					do kk = 1,K
						auxaux_u = auxaux_u + xaux(tt,kk)*theta0(kk,1)
					end do
					aux_u = aux_u + (yaux(tt,1) - auxaux_u - agt_tren0(tt,gg))**2
					auxaux_u = 0.0
				end do ! tt
				gi_classaux(i,gg) = aux_u
			end do !gg
		end do !i

		! Classification
		gi_num_prime = minloc(gi_classaux,2)
		gi_class_prime = minval(gi_classaux,2)

	! Counting number of elements in each group	
		
		do gg = 1,G
			countg_prime(1,gg) = 0
		end do
		do i = 1,N
			do gg = 1,G
				if (gi_num_prime(i) == gg) then
					countg_prime(1,gg) = countg_prime(1,gg)+1
				end if
			end do
		end do !i
		
		!We check for empty groups
		do gg = 1,G
			if (countg_prime(1,gg) == 0) then
				ii = maxloc(abs(gi_class),1)
				gi_num(ii) = gg
				countg(1,gg) = 1
				gi_class(ii) = 0.0d0
			end if
		end do
		

		end if !it's not the first round
	
		! Step 2: Update

		! We compute the within group means
		! We initialize to zero:

		flag = 1
		do gg = 1,G
			do tt = 1,T
				x1gt(tt,gg) = 0.0d0
				x2gt(tt,gg) = 0.0d0
				ygt(tt,gg) = 0.0d0
			end do
		end do

		do i = 1,N
			do gg = 1,G
				if (gi_num_prime(i) == gg) then
					do tt = 1,T
						x1gt(tt,gg) = x1gt(tt,gg) + X((i-1)*T+tt,1)/countg_prime(1,gg)
						x2gt(tt,gg) = x2gt(tt,gg) + X((i-1)*T+tt,2)/countg_prime(1,gg)
						ygt(tt,gg) = ygt(tt,gg) + Y((i-1)*T+tt,1)/countg_prime(1,gg)
					end do
				end if
			end do
		end do

		! We compute the demeaned vectors:
		do i = 1,N
			do gg = 1,G
				if (gi_num_prime(i) == gg) then
					do tt = 1,T
						X_dm((i-1)*T+tt,1) = X((i-1)*T+tt,1) - x1gt(tt,gg)
						X_dm((i-1)*T+tt,2) = X((i-1)*T+tt,2) - x2gt(tt,gg)
						Y_dm((i-1)*T+tt,1) = Y((i-1)*T+tt,1) - ygt(tt,gg)
					end do !tt
				end if
			end do !gg
		end do !i
		
		mat1 = matmul(transpose(X_dm),X_dm)
		!print*, mat1

		call lu_inverse(mat1,imat1)
		theta = matmul(matmul(imat1,transpose(X_dm)),Y_dm)

		Objh0 = 0.0d0
		do i = 1,N
			do tt = 1,T
				Objh0 = Objh0 + (Y_dm((i-1)*T+tt,1) - X_dm((i-1)*T+tt,1)*theta(1,1) - X_dm((i-1)*T+tt,2)*theta(2,1))**2
			end do
		end do
		
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
						agt_tren(tt,gg) = agt_tren(tt,gg) + (Y((i-1)*T+tt,1) - X((i-1)*T+tt,1)*theta(1,1) - X((i-1)*T+tt,2)*theta(2,1))/countg_prime(1,gg)
					end do
				end if
			end do
		end do

		deltapar_new = 0.0d0
		do kk = 1,K
			deltapar_new = deltapar_new + (theta(kk,1)-theta0(kk,1))**2
		end do
		do gg = 1,G
			do tt = 1,T
				deltapar_new = deltapar_new + (agt_tren(tt,gg) - agt_tren0(tt,gg))**2
			end do
		end do

		deltapar = deltapar_new

		theta0 = theta
		agt_tren0 = agt_tren
		
		loopstep = loopstep+1
		if(loopstep == 10000) then
			print*, deltapar
		end if
		end do !deltapar

		! Beginning of the local heuristic search
		phi_max = 1
		Objk = Objh0
		do while (phi_max == 1)
		phi = 0
		do i = 1,N
			! We identify the group to which i belongs
			ggi = gi_num_prime (i)
			if (countg_prime(1,ggi) > 1) then
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
						countg0(1,ggi) = countg_prime(1,ggi)-1
						countg0(1,gg) = countg_prime(1,gg)+1
						
						! We compute the within group means
						! We initialize to zero:
						do ggg = 1,G
							do tt = 1,T
								x1gt(tt,ggg) = 0.0d0
								x2gt(tt,ggg) = 0.0d0
								ygt(tt,ggg) = 0.0d0
							end do
						end do

						do j = 1,N
							do ggg = 1,G
								if (gi_num0(j) == ggg) then
									do tt = 1,T
										x1gt(tt,ggg) = x1gt(tt,ggg) + X((j-1)*T+tt,1)/countg0(1,ggg)
										x2gt(tt,ggg) = x2gt(tt,ggg) + X((j-1)*T+tt,2)/countg0(1,ggg)
										ygt(tt,ggg) = ygt(tt,ggg) + Y((j-1)*T+tt,1)/countg0(1,ggg)
									end do
								end if
							end do
						end do
						
						! We compute the demeaned vectors:
						do j = 1,N
							do ggg = 1,G
								if (gi_num0(j) == ggg) then
									do tt = 1,T
										X_dm((j-1)*T+tt,1) = X((j-1)*T+tt,1) - x1gt(tt,ggg)
										X_dm((j-1)*T+tt,2) = X((j-1)*T+tt,2) - x2gt(tt,ggg)
										Y_dm((j-1)*T+tt,1) = Y((j-1)*T+tt,1) - ygt(tt,ggg)
									end do !tt
								end if
							end do !ggg
						end do !j

						mat1 = matmul(transpose(X_dm),X_dm)
						!print*, mat1

						call lu_inverse(mat1,imat1)
						theta0 = matmul(matmul(imat1,transpose(X_dm)),Y_dm)
						
						Objk0 = 0.d0
						do j = 1,N
							do tt = 1,T
								Objk0 = Objk0 +  (Y_dm((j-1)*T+tt,1) - X_dm((j-1)*T+tt,1)*theta0(1,1) - X_dm((j-1)*T+tt,2)*theta0(2,1))**2
							end do
						end do

						if (Objk0 < Objk) then
							gi_num_prime = gi_num0
							countg_prime = countg0
							Objk = Objk0
							theta = theta0
							phi = 1
						end if
					end if !gg neq ggi
				end do !! gg
			end if !! there is more than 1 element in ggi
		end do ! i

		if (phi == 0) then
			phi_max = 0
		end if

		end do !while phi_max

		! Step 4: If the Objetive has improved, we reset gi_num = gi_prime (after HK) and 
		! restart fron neighbor 1

		Obj_prime = Objk
		if (Obj_prime < Obj) then
			print*, 'Updating partition'
			gi_num = gi_num_prime
			gi_num_opt = gi_num_prime
			theta_opt = theta
			countg = countg_prime
			Obj = Obj_prime
			control2 = 0 !---> go to Step1
		else
		! Step 5: If neighbour is not yet max, set neighbour = neighbour + 1 and go to Step 2 
			if (neighbour < Neighbourmax) then
				neighbour = neighbour + 1
			else
				control2 = 0 !---> if not, go to Step 1: with the same partition, draw all neighbours
				step = step + 1
			end if
		end if
		
		write(13,*), Obj, theta
		 
		end do !control2
		if (step > step_max) then
			control1 = 0
		end if
	end do !control1
	
	!We compute the final trends
	do gg = 1,G
			do tt = 1,T
				agt_tren_opt(tt,gg) = 0.d0
			end do
		end do
		
		do i = 1,N
			do gg = 1,G
				if(gi_num(i) == gg) then
					do tt = 1,T
						agt_tren_opt(tt,gg) = agt_tren_opt(tt,gg) + (Y((i-1)*T+tt,1) - X((i-1)*T+tt,1)*theta_opt(1,1) - X((i-1)*T+tt,2)*theta_opt(2,1))/countg(1,gg)
					end do
				end if
			end do
		end do
	
	! And the distance of each element to the trend (centroid)
	do i = 1,N
		do tt = 1,T
			yaux(tt,1) = Y((i-1)*T+tt,1)
			xaux(tt,:) = X((i-1)*T+tt,:)
		end do
		aux_u = 0.0
		auxaux_u = 0.0
		do tt = 1,T
			do kk = 1,K
				auxaux_u = auxaux_u + xaux(tt,kk)*theta_opt(kk,1)
			end do
			aux_u = aux_u + (yaux(tt,1) - auxaux_u - agt_tren_opt(tt,gi_num_opt(i)))**2
			auxaux_u = 0.0
		end do ! tt
	gi_class(i) = aux_u
	end do !i
	


	call date_and_time(values=time_array_1)
      finish_time = time_array_1 (5) * 3600 + time_array_1 (6) * 60 &
           + time_array_1 (7) + 0.001 * time_array_1 (8)

	write (13, '(8x, 1a, 1f16.6)') 'elapsed wall clock time:', &
           finish_time - start_time
	
	close(13)

	print*, theta_opt
	print*, Obj

	end do ! jsim
	end subroutine heuristics

  subroutine heuristics_simple(G,N,T,K,Y,X,Nsim,Neighbourmax,step_max,Obj00,gi_num_opt)


  !use global_module
  use imsl
  use linalg
  use utils
  
  implicit none
  
  ! Number of Groups to estimate (should be loaded)
  integer, intent(in) :: G,N,T,K,Nsim,Neighbourmax,step_max
  real(8), intent(in) :: X(N*T,K),Y(N*T,1)
  real(8), intent(out) :: Obj00
  integer, intent(out) :: gi_num_opt(N)

  ! Parameters of the model
  real(8) :: theta(K,1), theta0(K,1),deltapar_new,theta_opt(K,1),agt_tren_opt(T,G),gi_class(N)
  real(8) :: agt_vect(G*T,1), agt_vect0(G*T,1)
  real(8) :: agt_tren(T,G), agt_tren0(T,G)
  real(8) :: countg(1,G),countg0(1,G), countgaux(1,G), countg_opt(1,G),countg_prime(1,G)
  real(8) :: Obj0, Obj, Objh0, Objk, Objk0,Obj_prime

  ! Parameters of the heuristics
  ! Parameters of K-MEANS
  ! Parameters of classification
  integer :: i, tt, kk, ii, j, control1, control2,nn, step, loopstep, flag,ggaux
  real(8) :: xaux(T,K), yaux(T,1)
  real(8) :: agt_aux(T,1)
  real(8) :: aux_u, auxaux_u 
  integer :: gi_num0(N),gi_num(N), gi_num_prime(N), neighbour ! contains the group at wich i belongs
  real(8) :: gi(N,G) !contains a 1 in the column gg where i is assigned
  real(8) :: gi_class_prime(N) !contains the distances to the centroid associated with classification gi
  real(8) :: mat1(K,K), imat1(K,K)

  ! Parameters of the local search
  integer :: phi, phi_max,flagflag

  !Parameters of updating
  real(8) :: x1gt(T,G), x2gt(T,G), ygt(T,G),gi_classaux(N,G)
  real(8) :: X_dm(N*T,K), Y_dm(N*T,1), aux_vec(G)


  ! Auxiliary parameters
  real(8) :: tmp1(K), tmp2(G), tmp3(Neighbourmax), tmp4(Neighbourmax), deltapar
  real(8) :: V(G), W(N*T,1)
  integer :: gg, dd, jsim, ggg, ggi
  integer :: time_array_0(8), time_array_1(8)
  real(8) :: start_time, finish_time, aux


  call date_and_time(values=time_array_0)
  start_time = time_array_0 (5) * 3600 + time_array_0 (6) * 60 &
           + time_array_0 (7) + 0.001 * time_array_0 (8)

  !To write down the output
  open(unit=13,file='outputtheta.txt')
  
  ! We read the data
 ! open(unit=11,file='inputX_aux.txt')
 ! open(unit=12,file='inputY_aux.txt')

!	do dd = 1,N*T
!		read(11,*) X(dd,:)
!	end do

!	read(12,*) Y

!  close(11)
!  close(12)
  Obj00 = 100000.d0
  flagflag = 0

  ! Kmeans loop:
	do jsim = 1,Nsim
	!print*, jsim
	
	
	!Initialization:
	
	!We pick up an initial condition for beta
	call drnnor(K,tmp1)
    theta0 = reshape(tmp1, (/ K, 1 /))

	!We compute the residuals and select randomly G "centers"
	W = Y - matmul(X,theta0)
	call drnun(G,tmp2)
	do gg = 1,G
		V(gg) = max(nint(tmp2(gg)*N),1)
		agt_tren0(:,gg) = W((V(gg)-1)*T+1:V(gg)*T,1)
	end do
	
	

	!Iterations
	deltapar = 1.0d0
	do while (deltapar > 0)
		
		! Step 1: assignment
		do i = 1,N
			do tt = 1,T
				yaux(tt,1) = Y((i-1)*T+tt,1)
				xaux(tt,:) = X((i-1)*T+tt,:)
			end do
			do gg = 1,G
				aux_u = 0.0
				auxaux_u = 0.0
				do tt = 1,T
					do kk = 1,K
						auxaux_u = auxaux_u + xaux(tt,kk)*theta0(kk,1)
					end do
					aux_u = aux_u + (yaux(tt,1) - auxaux_u - agt_tren0(tt,gg))**2
					auxaux_u = 0.0
				end do ! tt
				gi_classaux(i,gg) = aux_u
			end do !gg
		end do !i

		! Classification
		gi_num = minloc(gi_classaux,2)
		gi_class = minval(gi_classaux,2)
		
		! Counting number of elements in each group

		do gg = 1,G
			countg(1,gg) = 0
		end do
		do i = 1,N
			do gg = 1,G
				if (gi_num(i) == gg) then
					countg(1,gg) = countg(1,gg)+1
				end if
			end do
		end do !i
		
		!We check for empty groups
		countgaux = countg
		do gg = 1,G
			if (countgaux(1,gg) == 0) then
				ii = maxloc(abs(gi_class),1)
				ggaux = gi_num(ii)
				gi_num(ii) = gg
				countg(1,gg) = 1
				countg(1,ggaux) = countg(1,ggaux)-1
				gi_class(ii) = 0.0d0
			end if
		end do
		

		! Step 2: update

		! We compute the within group means
		! We initialize to zero:
		do gg = 1,G
			do tt = 1,T
				x1gt(tt,gg) = 0.0d0
				x2gt(tt,gg) = 0.0d0
				ygt(tt,gg) = 0.0d0
			end do
		end do

		do i = 1,N
			do gg = 1,G
				if (gi_num(i) == gg) then
					do tt = 1,T
					x1gt(tt,gg) = x1gt(tt,gg) + X((i-1)*T+tt,1)/countg(1,gg)
					x2gt(tt,gg) = x2gt(tt,gg) + X((i-1)*T+tt,2)/countg(1,gg)
					ygt(tt,gg) = ygt(tt,gg) + Y((i-1)*T+tt,1)/countg(1,gg)
					end do
				end if
			end do
		end do

		! We compute the demeaned vectors:
		do i = 1,N
			do gg = 1,G
			if (gi_num(i) == gg) then
				do tt = 1,T
				X_dm((i-1)*T+tt,1) = X((i-1)*T+tt,1) - x1gt(tt,gg)
				X_dm((i-1)*T+tt,2) = X((i-1)*T+tt,2) - x2gt(tt,gg)
				Y_dm((i-1)*T+tt,1) = Y((i-1)*T+tt,1) - ygt(tt,gg)
				end do !tt
			end if
			end do !gg
		end do !i
		
		mat1 = matmul(transpose(X_dm),X_dm)
		!print*, mat1

		call lu_inverse(mat1,imat1)
		theta = matmul(matmul(imat1,transpose(X_dm)),Y_dm)

		Obj = 0.0d0
		do i = 1,N
			do tt = 1,T
				Obj = Obj + (Y_dm((i-1)*T+tt,1) - X_dm((i-1)*T+tt,1)*theta(1,1) - X_dm((i-1)*T+tt,2)*theta(2,1))**2
			end do
		end do
		
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
						agt_tren(tt,gg) = agt_tren(tt,gg) + (Y((i-1)*T+tt,1) - X((i-1)*T+tt,1)*theta(1,1) - X((i-1)*T+tt,2)*theta(2,1))/countg(1,gg)
					end do
				end if
			end do
		end do

		deltapar_new = 0.0d0
		do kk = 1,K
			deltapar_new = deltapar_new + (theta(kk,1)-theta0(kk,1))**2
		end do
		do gg = 1,G
			do tt = 1,T
				deltapar_new = deltapar_new + (agt_tren(tt,gg) - agt_tren0(tt,gg))**2
			end do
		end do

		deltapar = deltapar_new

		theta0 = theta
		agt_tren0 = agt_tren
		
		end do !deltapar
		
		if (Obj < Obj00) then
			theta_opt = theta
			agt_tren_opt = agt_tren
			Obj00 = Obj
			gi_num_opt = gi_num
			countg_opt = countg
		end if

	!---------------------------Beginning of VNS--------------------------------------------!

	! Step 0: read any initial K-cluster partition:

	flagflag = 0
	if (flagflag == 1) then
	theta0 = theta_opt
	gi_num = gi_num_opt
	countg = countg_opt
	agt_tren0 = agt_tren_opt
	Obj = Obj00
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

	! We randomize the individuals and clusters
	call drnun(Neighbourmax,tmp3)
	call drnun(Neighbourmax,tmp4)
	

	gi_num_prime = gi_num
	countg_prime = countg
	
	
	do nn = 1,neighbour
		ii = min(max(nint(tmp3(nn)*N),1),N)
		ggi = min(max(nint(tmp4(nn)*G),1),G)
		do while (gi_num(ii) == ggi)
			call drnun(1,aux)
			ggi = max(nint(aux*G),1)
		end do
		gi_num_prime(ii) = ggi
		countg_prime(1,ggi) = countg(1,ggi) + 1
		countg_prime(1,gi_num(ii)) = countg(1,gi_num(ii))-1
	end do

	!Step 3: Apply HK-means unsing gi_num_prime as the starting soluiton
	!print*, 'HK-means'
	deltapar = 1.0d0
	!do while (deltapar > 0 .and. loopstep < 10000)
	
	flag = 0
	do while (deltapar > 0)
	if (flag == 1) then !IE this is not the first time I enter this loop
		! Step 1: assignment
		do i = 1,N
			do tt = 1,T
				yaux(tt,1) = Y((i-1)*T+tt,1)
				xaux(tt,:) = X((i-1)*T+tt,:)
			end do
			do gg = 1,G
				aux_u = 0.0
				auxaux_u = 0.0
				do tt = 1,T
					do kk = 1,K
						auxaux_u = auxaux_u + xaux(tt,kk)*theta0(kk,1)
					end do
					aux_u = aux_u + (yaux(tt,1) - auxaux_u - agt_tren0(tt,gg))**2
					auxaux_u = 0.0
				end do ! tt
				gi_classaux(i,gg) = aux_u
			end do !gg
		end do !i

		! Classification
		gi_num_prime = minloc(gi_classaux,2)
		gi_class_prime = minval(gi_classaux,2)

	! Counting number of elements in each group	
		
		do gg = 1,G
			countg_prime(1,gg) = 0
		end do
		do i = 1,N
			do gg = 1,G
				if (gi_num_prime(i) == gg) then
					countg_prime(1,gg) = countg_prime(1,gg)+1
				end if
			end do
		end do !i
		
		!We check for empty groups
		do gg = 1,G
			if (countg_prime(1,gg) == 0) then
				ii = maxloc(abs(gi_class),1)
				gi_num(ii) = gg
				countg(1,gg) = 1
				gi_class(ii) = 0.0d0
			end if
		end do
		

		end if !it's not the first round
	
		! Step 2: Update

		! We compute the within group means
		! We initialize to zero:

		flag = 1
		do gg = 1,G
			do tt = 1,T
				x1gt(tt,gg) = 0.0d0
				x2gt(tt,gg) = 0.0d0
				ygt(tt,gg) = 0.0d0
			end do
		end do

		do i = 1,N
			do gg = 1,G
				if (gi_num_prime(i) == gg) then
					do tt = 1,T
						x1gt(tt,gg) = x1gt(tt,gg) + X((i-1)*T+tt,1)/countg_prime(1,gg)
						x2gt(tt,gg) = x2gt(tt,gg) + X((i-1)*T+tt,2)/countg_prime(1,gg)
						ygt(tt,gg) = ygt(tt,gg) + Y((i-1)*T+tt,1)/countg_prime(1,gg)
					end do
				end if
			end do
		end do

		! We compute the demeaned vectors:
		do i = 1,N
			do gg = 1,G
				if (gi_num_prime(i) == gg) then
					do tt = 1,T
						X_dm((i-1)*T+tt,1) = X((i-1)*T+tt,1) - x1gt(tt,gg)
						X_dm((i-1)*T+tt,2) = X((i-1)*T+tt,2) - x2gt(tt,gg)
						Y_dm((i-1)*T+tt,1) = Y((i-1)*T+tt,1) - ygt(tt,gg)
					end do !tt
				end if
			end do !gg
		end do !i
		
		mat1 = matmul(transpose(X_dm),X_dm)
		!print*, mat1

		call lu_inverse(mat1,imat1)
		theta = matmul(matmul(imat1,transpose(X_dm)),Y_dm)

		Objh0 = 0.0d0
		do i = 1,N
			do tt = 1,T
				Objh0 = Objh0 + (Y_dm((i-1)*T+tt,1) - X_dm((i-1)*T+tt,1)*theta(1,1) - X_dm((i-1)*T+tt,2)*theta(2,1))**2
			end do
		end do
		
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
						agt_tren(tt,gg) = agt_tren(tt,gg) + (Y((i-1)*T+tt,1) - X((i-1)*T+tt,1)*theta(1,1) - X((i-1)*T+tt,2)*theta(2,1))/countg_prime(1,gg)
					end do
				end if
			end do
		end do

		deltapar_new = 0.0d0
		do kk = 1,K
			deltapar_new = deltapar_new + (theta(kk,1)-theta0(kk,1))**2
		end do
		do gg = 1,G
			do tt = 1,T
				deltapar_new = deltapar_new + (agt_tren(tt,gg) - agt_tren0(tt,gg))**2
			end do
		end do

		deltapar = deltapar_new

		theta0 = theta
		agt_tren0 = agt_tren
		
		loopstep = loopstep+1
		if(loopstep == 10000) then
			!print*, deltapar
		end if
		end do !deltapar

		! Beginning of the local heuristic search
		phi_max = 1
		Objk = Objh0
		do while (phi_max == 1)
		phi = 0
		do i = 1,N
			! We identify the group to which i belongs
			ggi = gi_num_prime (i)
			if (countg_prime(1,ggi) > 1) then
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
						countg0(1,ggi) = countg_prime(1,ggi)-1
						countg0(1,gg) = countg_prime(1,gg)+1
						
						! We compute the within group means
						! We initialize to zero:
						do ggg = 1,G
							do tt = 1,T
								x1gt(tt,ggg) = 0.0d0
								x2gt(tt,ggg) = 0.0d0
								ygt(tt,ggg) = 0.0d0
							end do
						end do

						do j = 1,N
							do ggg = 1,G
								if (gi_num0(j) == ggg) then
									do tt = 1,T
										x1gt(tt,ggg) = x1gt(tt,ggg) + X((j-1)*T+tt,1)/countg0(1,ggg)
										x2gt(tt,ggg) = x2gt(tt,ggg) + X((j-1)*T+tt,2)/countg0(1,ggg)
										ygt(tt,ggg) = ygt(tt,ggg) + Y((j-1)*T+tt,1)/countg0(1,ggg)
									end do
								end if
							end do
						end do
						
						! We compute the demeaned vectors:
						do j = 1,N
							do ggg = 1,G
								if (gi_num0(j) == ggg) then
									do tt = 1,T
										X_dm((j-1)*T+tt,1) = X((j-1)*T+tt,1) - x1gt(tt,ggg)
										X_dm((j-1)*T+tt,2) = X((j-1)*T+tt,2) - x2gt(tt,ggg)
										Y_dm((j-1)*T+tt,1) = Y((j-1)*T+tt,1) - ygt(tt,ggg)
									end do !tt
								end if
							end do !ggg
						end do !j

						mat1 = matmul(transpose(X_dm),X_dm)
						!print*, mat1

						call lu_inverse(mat1,imat1)
						theta0 = matmul(matmul(imat1,transpose(X_dm)),Y_dm)
						
						Objk0 = 0.d0
						do j = 1,N
							do tt = 1,T
								Objk0 = Objk0 +  (Y_dm((j-1)*T+tt,1) - X_dm((j-1)*T+tt,1)*theta0(1,1) - X_dm((j-1)*T+tt,2)*theta0(2,1))**2
							end do
						end do

						if (Objk0 < Objk) then
							gi_num_prime = gi_num0
							countg_prime = countg0
							Objk = Objk0
							theta = theta0
							phi = 1
						end if
					end if !gg neq ggi
				end do !! gg
			end if !! there is more than 1 element in ggi
		end do ! i

		if (phi == 0) then
			phi_max = 0
		end if

		end do !while phi_max

		! Step 4: If the Objetive has improved, we reset gi_num = gi_prime (after HK) and 
		! restart fron neighbor 1

		Obj_prime = Objk
		if (Obj_prime < Obj) then
			!print*, 'Updating partition'
			gi_num = gi_num_prime
			gi_num_opt = gi_num_prime
			theta_opt = theta
			countg = countg_prime
			Obj = Obj_prime
			control2 = 0 !---> go to Step1
		else
		! Step 5: If neighbour is not yet max, set neighbour = neighbour + 1 and go to Step 2 
			if (neighbour < Neighbourmax) then
				neighbour = neighbour + 1
			else
				control2 = 0 !---> if not, go to Step 1: with the same partition, draw all neighbours
				step = step + 1
			end if
		end if
		
		write(13,*), Obj, theta
		 
		end do !control2
		if (step > step_max) then
			control1 = 0
		end if
	end do !control1
	
	!We compute the final trends
!	do gg = 1,G
!			do tt = 1,T
!				agt_tren_opt(tt,gg) = 0.d0
!			end do
!		end do
!		
!		do i = 1,N
!			do gg = 1,G
!				if(gi_num(i) == gg) then
!					do tt = 1,T
!						agt_tren_opt(tt,gg) = agt_tren_opt(tt,gg) + (Y((i-1)*T+tt,1) - X((i-1)*T+tt,1)*theta_opt(1,1) - X((i-1)*T+tt,2)*theta_opt(2,1))/countg(1,gg)
!					end do
!				end if
!			end do
!		end do
	
	! And the distance of each element to the trend (centroid)
!	do i = 1,N
!		do tt = 1,T
!			yaux(tt,1) = Y((i-1)*T+tt,1)
!			xaux(tt,:) = X((i-1)*T+tt,:)
!		end do
!		aux_u = 0.0
!		auxaux_u = 0.0
!		do tt = 1,T
!			do kk = 1,K
!				auxaux_u = auxaux_u + xaux(tt,kk)*theta_opt(kk,1)
!			end do
!			aux_u = aux_u + (yaux(tt,1) - auxaux_u - agt_tren_opt(tt,gi_num_opt(i)))**2
!			auxaux_u = 0.0
!		end do ! tt
!	gi_class(i) = aux_u
!	end do !i
	


!	call date_and_time(values=time_array_1)
!      finish_time = time_array_1 (5) * 3600 + time_array_1 (6) * 60 &
!           + time_array_1 (7) + 0.001 * time_array_1 (8)

!	write (6, '(8x, 1a, 1f16.6)') 'elapsed wall clock time:', &
!           finish_time - start_time
	
	close(13)

	!print*, theta_opt
	!print*, Obj

end if

	end do ! jsim
	end subroutine heuristics_simple



	subroutine CALCFUNCOBJ(N,T,G,gi_num,Y,X,Obj)
	use imsl
	use linalg
	use utils

	real(8), intent(out) :: Obj
	integer, intent(in) :: N,T,G,gi_num(N)
	real(8), intent(in) :: Y(N*T,1),X(N*T,2)

	real(8) :: x1gt(T,G),x2gt(T,G),ygt(t,g)
	real(8) :: X_dm(N*T,2), Y_dm(N*T,1)
	integer :: countg(1,G),gg,i,tt,j
	real(8) :: mat1(2,2), imat1(2,2), theta(2,1),disc,Objy,Objy2,Objxy1,Objxy12,Objxy2,Objxy22

	do gg = 1,G
		countg(1,gg) = 0
	end do
	
	do i = 1,N
		do gg = 1,G
			if (gi_num(i) == gg) then
				countg(1,gg) = countg(1,gg) + 1
			end if
		end do
	end do
	

	! We compute the within group means
	! We initialize to zero:
	do gg = 1,G
		do tt = 1,T
			x1gt(tt,gg) = 0.0d0
			x2gt(tt,gg) = 0.0d0
			ygt(tt,gg) = 0.0d0
		end do
	end do

	do i = 1,N
		do gg = 1,G
			if (gi_num(i) == gg) then
				do tt = 1,T
				x1gt(tt,gg) = x1gt(tt,gg) + X((i-1)*T+tt,1)/dfloat(countg(1,gg))
				x2gt(tt,gg) = x2gt(tt,gg) + X((i-1)*T+tt,2)/dfloat(countg(1,gg))
				ygt(tt,gg) = ygt(tt,gg) + Y((i-1)*T+tt,1)/dfloat(countg(1,gg))
				end do
			end if
		end do
	end do

   ! We compute the demeaned vectors:
	do i = 1,N
		do gg = 1,G
		if (gi_num(i) == gg) then
			do tt = 1,T
			X_dm((i-1)*T+tt,1) = X((i-1)*T+tt,1) - x1gt(tt,gg)
			X_dm((i-1)*T+tt,2) = X((i-1)*T+tt,2) - x2gt(tt,gg)
			Y_dm((i-1)*T+tt,1) = Y((i-1)*T+tt,1) - ygt(tt,gg)
			end do !tt
		end if
		end do !gg
	end do !i

	mat1 = matmul(transpose(X_dm),X_dm)
	!print*, mat1

	!call lu_inverse(mat1,imat1)
	disc = mat1(1,1)*mat1(2,2) - mat1(1,2)*mat1(2,1)
	!IF (DISC.EQ.0.) THEN
	!	PRINT*, 'ERROR'
	!END IF
	imat1(1,1) = mat1(2,2)/disc
	imat1(1,2) = - mat1(1,2)/disc
	imat1(2,1) = - mat1(1,2)/disc
	imat1(2,2) = mat1(1,1)/disc
	theta = matmul(matmul(imat1,transpose(X_dm)),Y_dm)

	Obj = 0.0d0
	do i = 1,N
		do tt = 1,T
			Obj = Obj + (Y_dm((i-1)*T+tt,1) - X_dm((i-1)*T+tt,1)*theta(1,1) - X_dm((i-1)*T+tt,2)*theta(2,1))**2
		end do
	end do

	Objy = 0.0d0
	do i = 1,N
		do tt = 1,T
			Objy = Objy + (Y_dm((i-1)*T+tt,1))**2
		end do
	end do

	Objy2 = 0.0d0
	do gg = 1,G
	do i = 1,N
		do j = i+1,N
			if (gi_num(i) == gg) then
				if (gi_num(j) == gg) then
					do tt = 1,T
						Objy2 = Objy2 + (Y((i-1)*T+tt,1)-Y((j-1)*T+tt,1))**2/countg(1,gg)
					end do
				end if
			end if
		end do
	end do
	end do

	
	Objxy1 = 0.0d0
	do i = 1,N
		do tt = 1,T
			Objxy1 = Objxy1 + (X_dm((i-1)*T+tt,1)*X_dm((i-1)*T+tt,2))
		end do
	end do


	Objxy2 = 0.0d0
	do i = 1,N
		do tt = 1,T
			Objxy2 = Objxy2 + (X_dm((i-1)*T+tt,2))**2
		end do
	end do

end subroutine CALCFUNCOBJ

subroutine COMPUTE_MATRIX(N,G,T,Y,X,A,BX1,BX2,CX1,CX2,CX3)
	
INTEGER, intent(in) :: N,G,T
DOUBLE PRECISION, intent(in) :: Y(N*T,1),X(N*T,2)
DOUBLE PRECISION, intent(out) :: A(N,N),BX1(N,N),BX2(N,N),CX1(N,N),CX2(N,N),CX3(N,N)

integer i,j,tt

do i = 1,N
	do j = 1,N
		A(i,j) = 0.d0
		BX1(i,j) = 0.d0
		BX2(i,j) = 0.d0
		CX1(i,j) = 0.d0
		CX2(i,j) = 0.d0
		CX3(i,j) = 0.d0
	end do
end do


do i = 1,N
	do j = 1,N
		do tt = 1,T
			A(i,j) = A(i,j) + (Y((i-1)*T+tt,1)-Y((j-1)*T+tt,1))**2
			BX1(i,j) = BX1(i,j) + (Y((i-1)*T+tt,1)-Y((j-1)*T+tt,1))*(X((i-1)*T+tt,1)-X((j-1)*T+tt,1))
			BX2(i,j) = BX2(i,j) + (Y((i-1)*T+tt,1)-Y((j-1)*T+tt,1))*(X((i-1)*T+tt,2)-X((j-1)*T+tt,2))
			CX1(i,j) = CX1(i,j) + (X((i-1)*T+tt,1)-X((j-1)*T+tt,1))**2
			CX3(i,j) = CX3(i,j) + (X((i-1)*T+tt,2)-X((j-1)*T+tt,2))**2
			CX2(i,j) = CX2(i,j) + (X((i-1)*T+tt,1)-X((j-1)*T+tt,1))*(X((i-1)*T+tt,2)-X((j-1)*T+tt,2))
		end do
	end do
end do




end subroutine COMPUTE_MATRIX


end module heuristics_sub