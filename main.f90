      ! This is a small 3-D Lid-Driven Cavity incompressible flow solver in Fortran.
      ! Details:
      !      - Staggered grid with uniform spacing is used: dx=dy=dz=h
      !      - The PPE is solved, to enforce div(u)=0, with red-black Gauss-Seidel method
      !      - Time marching: Explicit Euler (1st order accurate)
      !      - Pressure is treated purely implicitly
      !      - One ghost cell surrounds all interior cells
      !      - OpenMP parallelization is implemented
      !      - Compile command: "gfortran -fopenmp -std=gnu -O3 -g main.f90 -o main"
      ! 
      ! Developed by Charlie Kawczynski.
      module PPE_solver
      implicit none
      integer,parameter :: cp = selected_real_kind(14)
      private; public :: red_black
      contains
      subroutine red_black(p,divU,fact,N,h2,odd)
      implicit none
      integer,intent(in) :: N
      integer,dimension(3),intent(in) :: odd
      real(cp),intent(in) :: fact,h2
      real(cp),dimension(N+2,N+2,N+2),intent(in) :: divU
      real(cp),dimension(N+2,N+2,N+2),intent(inout) :: p
      integer :: i,j,k
      !$OMP PARALLEL DO
      do k=2+odd(3),N+1,2; do j=2+odd(2),N+1,2; do i=2+odd(1),N+1,2
      p(i,j,k)=fact*(p(i+1,j,k)+p(i,j+1,k)+p(i,j,k+1)+p(i-1,j,k)+p(i,j-1,k)+p(i,j,k-1)-divU(i,j,k)*h2)
      enddo; enddo; enddo
      !$OMP END PARALLEL DO
      end subroutine
      end module

      program main
      use PPE_solver
      implicit none
      integer,parameter :: cp = selected_real_kind(14)
      integer,parameter :: N = 45
      real(cp),dimension(:,:,:),allocatable :: p,divU,u,ustar,v,vstar,w,wstar
      real(cp),dimension(:,:,:),allocatable :: E1_x,E2_x,E1_y,E2_y,E1_z,E2_z
      real(cp) :: h,dt,Re,hdt_inv,h_inv,dt_inv,Re_inv,h2,h2_inv,fact,KE,KE_old,KE_temp,dV,tol
      integer :: i,j,k,N_output,N_PPE,iter_PPE,iter,N_iter

      h = 1.0_cp/real(N,cp)            ! spatial step size (hard coded and uniform)
      Re = 400.0_cp                    ! Reynolds number
      N_iter = 10**5                   ! number of time steps
      dt = 1.0_cp*10.0_cp**(-3.0_cp)   ! time step
      N_PPE = 5                        ! number of PPE iterations
      N_output = 100                   ! output transient data every N_output time steps
      tol = 1.0_cp*10.0_cp**(-12.0_cp) ! stops simulation when |KE-KE_old|<tol

      ! Initialize data
      KE_old = 0.0_cp; KE = 0.0_cp; fact = 1.0_cp/6.0_cp; h_inv = 1.0_cp/h; h2 = h**2.0_cp; h2_inv = 1.0_cp/h2
      Re_inv = 1.0_cp/Re; dt_inv = 1.0_cp/dt; dV = h**3.0_cp; hdt_inv = h_inv*dt_inv
      write(*,*) '3-D Lid-driven cavity flow. Re,h,dt,N_iter,N_PPE:'
      write(*,*) Re,h,dt,N_iter,N_PPE
      write(*,*) ''

      allocate(u(N+3,N+2,N+2),ustar(N+3,N+2,N+2),E1_x(N+2,N+3,N+3),E2_x(N+2,N+3,N+3),divU(N+2,N+2,N+2))
      allocate(v(N+2,N+3,N+2),vstar(N+2,N+3,N+2),E1_y(N+3,N+2,N+3),E2_y(N+3,N+2,N+3),p(N+2,N+2,N+2))
      allocate(w(N+2,N+2,N+3),wstar(N+2,N+2,N+3),E1_z(N+3,N+3,N+2),E2_z(N+3,N+3,N+2))

      u = 0.0_cp; ustar = u; E1_x = 0.0_cp; E2_x = 0.0_cp; p = 0.0_cp
      v = 0.0_cp; vstar = v; E1_y = 0.0_cp; E2_y = 0.0_cp; divU = p
      w = 0.0_cp; wstar = w; E1_z = 0.0_cp; E2_z = 0.0_cp

      call system('mkdir output')
      open(1,file='output/KE.dat')
      write(1,*) 'TITLE="KINETIC ENERGY VS TIME"'
      write(1,*) 'VARIABLES = "Time","Kinetic Energy"'
      write(1,*) 'ZONE DATAPACKING = POINT'

      do iter=1,N_iter ! Momentum equation

        !$OMP PARALLEL DO
        do k=2,N+2; do j=2,N+2; do i=2,N+2 ! Advection term in form d/dxj (uj ui)
        ! d/dxj (uj ui) for i=j
        ustar(i,j,k) = -0.25_cp*( (u(i,j,k)+u(i+1,j,k))**2 - ((u(i-1,j,k)+u(i,j,k))**2) )*h_inv
        vstar(i,j,k) = -0.25_cp*( (v(i,j,k)+v(i,j+1,k))**2 - ((v(i,j-1,k)+v(i,j,k))**2) )*h_inv
        wstar(i,j,k) = -0.25_cp*( (w(i,j,k)+w(i,j,k+1))**2 - ((w(i,j,k-1)+w(i,j,k))**2) )*h_inv
        ! d/dxj (uj ui) for i≠j
        E1_y(i,j,k) = 0.25_cp*(u(i,j,k-1)+u(i,j,k))*(w(i-1,j,k)+w(i,j,k)) ! x (y edge)
        E1_z(i,j,k) = 0.25_cp*(u(i,j-1,k)+u(i,j,k))*(v(i-1,j,k)+v(i,j,k)) ! x (z edge)
        E1_x(i,j,k) = 0.25_cp*(v(i,j,k-1)+v(i,j,k))*(w(i,j-1,k)+w(i,j,k)) ! y (x edge)
        E2_z(i,j,k) = 0.25_cp*(v(i-1,j,k)+v(i,j,k))*(u(i,j-1,k)+u(i,j,k)) ! y (z edge)
        E2_x(i,j,k) = 0.25_cp*(w(i,j-1,k)+w(i,j,k))*(v(i,j,k-1)+v(i,j,k)) ! z (x edge)
        E2_y(i,j,k) = 0.25_cp*(w(i-1,j,k)+w(i,j,k))*(u(i,j,k-1)+u(i,j,k)) ! z (y edge)
        enddo; enddo; enddo
        !$OMP END PARALLEL DO
        !$OMP PARALLEL DO
        do k=2,N+1; do j=2,N+1; do i=2,N+1 ! Advection term d/dxj (uj ui) (continued)
        ustar(i,j,k)=ustar(i,j,k)-(E1_y(i,j,k+1)-E1_y(i,j,k)+E1_z(i,j+1,k)-E1_z(i,j,k))*h_inv
        vstar(i,j,k)=vstar(i,j,k)-(E1_x(i,j,k+1)-E1_x(i,j,k)+E2_z(i+1,j,k)-E2_z(i,j,k))*h_inv
        wstar(i,j,k)=wstar(i,j,k)-(E2_x(i,j+1,k)-E2_x(i,j,k)+E2_y(i+1,j,k)-E2_y(i,j,k))*h_inv
        enddo; enddo; enddo
        !$OMP END PARALLEL DO

        !$OMP PARALLEL DO
        do k=2,N+1; do j=2,N+1; do i=2,N+1 ! Diffusion term
        ustar(i,j,k)=u(i,j,k)+dt*(ustar(i,j,k)+Re_inv*h2_inv*((u(i+1,j,k)-2.0_cp*u(i,j,k)+u(i-1,j,k)) + &
                                                              (u(i,j+1,k)-2.0_cp*u(i,j,k)+u(i,j-1,k)) + &
                                                              (u(i,j,k+1)-2.0_cp*u(i,j,k)+u(i,j,k-1))))
        vstar(i,j,k)=v(i,j,k)+dt*(vstar(i,j,k)+Re_inv*h2_inv*((v(i+1,j,k)-2.0_cp*v(i,j,k)+v(i-1,j,k)) + &
                                                              (v(i,j+1,k)-2.0_cp*v(i,j,k)+v(i,j-1,k)) + &
                                                              (v(i,j,k+1)-2.0_cp*v(i,j,k)+v(i,j,k-1))))
        wstar(i,j,k)=w(i,j,k)+dt*(wstar(i,j,k)+Re_inv*h2_inv*((w(i+1,j,k)-2.0_cp*w(i,j,k)+w(i-1,j,k)) + &
                                                              (w(i,j+1,k)-2.0_cp*w(i,j,k)+w(i,j-1,k)) + &
                                                              (w(i,j,k+1)-2.0_cp*w(i,j,k)+w(i,j,k-1))))
        enddo; enddo; enddo
        !$OMP END PARALLEL DO

        ustar( 2 ,:,:) = 0.0_cp; vstar(:, 2 ,:) = 0.0_cp; wstar(:,:, 2 ) = 0.0_cp ! Remove forcing on wall (kinematic)
        ustar(N+2,:,:) = 0.0_cp; vstar(:,N+2,:) = 0.0_cp; wstar(:,:,N+2) = 0.0_cp ! Remove forcing on wall (kinematic)

        !$OMP PARALLEL DO
        do k=2,N+1; do j=2,N+1; do i=2,N+1 ! Compute PPE source
        divU(i,j,k) = hdt_inv*(ustar(i+1,j,k)-ustar(i,j,k)+vstar(i,j+1,k)-vstar(i,j,k)+wstar(i,j,k+1)-wstar(i,j,k))
        enddo; enddo; enddo
        !$OMP END PARALLEL DO

        do iter_PPE=1,N_PPE ! solve PPE: red-black Gauss-Seidel
          call red_black(p,divU,fact,N,h2,(/0,0,0/)); call red_black(p,divU,fact,N,h2,(/1,0,0/))
          call red_black(p,divU,fact,N,h2,(/0,1,0/)); call red_black(p,divU,fact,N,h2,(/0,0,1/))
          call red_black(p,divU,fact,N,h2,(/1,1,1/)); call red_black(p,divU,fact,N,h2,(/0,1,1/))
          call red_black(p,divU,fact,N,h2,(/1,0,1/)); call red_black(p,divU,fact,N,h2,(/1,1,0/))
          p( 1 ,:,:) = p( 2 ,:,:); p(:, 1 ,:) = p(:, 2 ,:); p(:,:, 1 ) = p(:,:, 2 ) ! Apply p BCs
          p(N+2,:,:) = p(N+1,:,:); p(:,N+2,:) = p(:,N+1,:); p(:,:,N+2) = p(:,:,N+1) ! Apply p BCs
        enddo

        !$OMP PARALLEL DO
        do k=2,N+1; do j=2,N+1; do i=2,N+1 ! Pressure correction
          u(i,j,k) = ustar(i,j,k) - dt*(h_inv*(p(i,j,k) - p(i-1,j,k)))
          v(i,j,k) = vstar(i,j,k) - dt*(h_inv*(p(i,j,k) - p(i,j-1,k)))
          w(i,j,k) = wstar(i,j,k) - dt*(h_inv*(p(i,j,k) - p(i,j,k-1)))
        enddo; enddo; enddo
        !$OMP END PARALLEL DO

        u( 2 ,:,:) = 0.0_cp; v(:, 2 ,:) = 0.0_cp; w(:,:, 2 ) = 0.0_cp ! Apply u BCs (wall coincident)
        u(N+2,:,:) = 0.0_cp; v(:,N+2,:) = 0.0_cp; w(:,:,N+2) = 0.0_cp ! Apply u BCs (wall coincident)

        ! Apply u BCs (ghost cells, including sliding lid)
        u(:,1,:)=-u(:,2,:); u(:,N+2,:)=2.0_cp-u(:,N+1,:);u(:,:,1)=-u(:,:,2); u(:,:,N+2)=-u(:,:,N+1)
        v(1,:,:)=-v(2,:,:); v(N+2,:,:)=-v(N+1,:,:);v(:,:,1)=-v(:,:,2); v(:,:,N+2)=-v(:,:,N+1)
        w(1,:,:)=-w(2,:,:); w(N+2,:,:)=-w(N+1,:,:);w(:,1,:)=-w(:,2,:); w(:,N+2,:)=-w(:,N+1,:)

        if (mod(iter+1,N_output).eq.1) then
          KE_temp = 0.0_cp
          !$OMP PARALLEL DO REDUCTION(+:KE_temp)
          do k=2,N+1; do j=2,N+1; do i=2,N+1 ! Pressure correction
          KE_temp = KE_temp + (u(i,j,k)+u(i+1,j,k))**2+(v(i,j,k)+v(i,j+1,k))**2+(w(i,j,k)+w(i,j,k+1))**2
          enddo; enddo; enddo
          !$OMP END PARALLEL DO
          KE_old = KE_temp; KE_old = 0.25_cp*KE_old*dV
          if ((iter.gt.1).and.(abs(KE-KE_old).lt.tol)) then; write(*,*) 'Exited early'; exit
          endif
          write(1,*) iter*dt,KE; flush(1)
          write(*,'(A35,F16.8,I13,2F16.8,E20.4E2)') 't,iter,percent complete,KE,|dKE| = ',&
          iter*dt,iter,real(iter,cp)/real(N_iter,cp)*100.0_cp,KE,abs(KE-KE_old)
          KE = KE_old
        endif
      enddo
      close(1) ! Close KE unit

      !$OMP PARALLEL DO
      do k=2,N+1; do j=2,N+1; do i=2,N+1 ! Compute divU
      divU(i,j,k) = h_inv*(u(i+1,j,k)-u(i,j,k)+v(i,j+1,k)-v(i,j,k)+w(i,j,k+1)-w(i,j,k))
      enddo; enddo; enddo
      !$OMP END PARALLEL DO

      open(2,file='output/slice.dat') ! Export solution at center plane
      write(2,*) 'TITLE="3-D CUBIC LDC AT CENTER PLANE"'
      write(2,*) 'VARIABLES = "x","y","u","v","p","divU"'
      write(2,*) 'ZONE, I=',N,',J=',N,' DATAPACKING = POINT'
      k=N/2; do j=2,N+1; do i=2,N+1
      write(2,*) (i-1)*h-0.5*h,(j-1)*h-0.5*h,&
      0.5_cp*(u(i+1,j,k)+u(i,j,k)),0.5_cp*(v(i,j+1,k)+v(i,j,k)),p(i,j,k),divU(i,j,k)
      enddo; enddo; close(2)

      open(3,file='output/solution.dat') ! Export solution at cell centers
      write(3,*) 'TITLE="3-D CUBIC LDC AT CELL CENTERS"'
      write(3,*) 'VARIABLES = "x","y","z","u","v","w","p","divU"'
      write(3,*) 'ZONE, I=',N,',J=',N,',K=',N,' DATAPACKING = POINT'
      do k=2,N+1; do j=2,N+1; do i=2,N+1
      write(3,*) (i-1)*h-0.5*h,(j-1)*h-0.5*h,(k-1)*h-0.5*h,&
      0.5_cp*(u(i+1,j,k)+u(i,j,k)),0.5_cp*(v(i,j+1,k)+v(i,j,k)),0.5_cp*(w(i,j,k+1)+w(i,j,k)),&
      p(i,j,k),divU(i,j,k)
      enddo; enddo; enddo; close(3)

      open(5,file='output/solution_n.dat') ! Export solution at nodes
      write(5,*) 'TITLE="3-D CUBIC LDC AT CELL CORNERS"'
      write(5,*) 'VARIABLES = "x","y","z","u","v","w"'
      write(5,*) 'ZONE, I=',N+1,',J=',N+1,',K=',N+1,' DATAPACKING = POINT'
      do k=2,N+2; do j=2,N+2; do i=2,N+2
      write(5,*) (i-2)*h,(j-2)*h,(k-2)*h,0.25_cp*(u(i,j-1,k)+u(i,j-1,k-1)+u(i,j,k-1)+u(i,j,k)),&
                                         0.25_cp*(v(i-1,j,k)+v(i-1,j,k-1)+v(i,j,k-1)+v(i,j,k)),&
                                         0.25_cp*(w(i-1,j,k)+w(i-1,j-1,k)+w(i,j-1,k)+w(i,j,k))
      enddo; enddo; enddo; close(5)

      deallocate(p,divU,u,v,w,ustar,vstar,wstar,E1_x,E1_y,E1_z,E2_x,E2_y,E2_z)

     end program