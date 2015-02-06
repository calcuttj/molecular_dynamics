! PHY 480 Project 2
! 2/3/15
Program lattice
  Implicit None
  real(8):: temp
  real(8):: mass
  integer ::cell_dim,N,counter
  real(8), allocatable :: pos(:,:)
  real(8), allocatable :: vel(:,:)

  temp = 1d0
  mass = 1d0
  cell_dim = 2
  N = 4*cell_dim**3
  allocate( pos(N,3) )
  allocate( vel(N,3) )

  call init_random_seed()
  call position_initializer(N,cell_dim,pos)
  do counter=1,N
     print *,counter, pos(counter,1),pos(counter,2),pos(counter,3)
  end do
  call velocity_initializer(temp,N,mass,vel)
  do counter=1,N
     print *,counter, vel(counter,1),vel(counter,2),vel(counter,3)
  end do
end program

subroutine init_random_seed()
  implicit none

  integer, allocatable :: seed(:)
  integer :: i, n, un, istat, dt(8), pid, t(2), s
  integer(8) :: count, tms
  call random_seed(size = n)
  allocate(seed(n))
  open(newunit=un, file="/dev/urandom", access="stream",&
       form="unformatted", action="read", status="old", &
       iostat=istat)
  if (istat == 0) then
     read(un) seed
     close(un)
  else
     call system_clock(count)
     if (count /= 0) then
        t = transfer(count, t)
     else
        call date_and_time(values=dt)
        tms = (dt(1) - 1970)*365_8 * 24 * 60 * 60 * 1000 &
             + dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
             + dt(3) * 24 * 60 * 60 * 60 * 1000 &
             + dt(5) * 60 * 60 * 1000 &
             + dt(6) * 60 * 1000 + dt(7) * 1000 &
             + dt(8)
        t = transfer(tms, t)
     end if
     s = ieor(t(1), t(2))
     pid = getpid() + 1099279 ! Add a prime
     s = ieor(s, pid)
     if (n >= 3) then
        seed(1) = t(1) + 36269
        seed(2) = t(2) + 72551
        seed(3) = pid
        if (n > 3) then
           seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
        end if
     else
        seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
     end if
  end if
  call random_seed(put=seed)
end subroutine init_random_seed

subroutine position_initializer(N,dim_cells,position)
  integer :: N,dim_cells
  real(8),intent(out), dimension(N,3) :: position
  integer :: counter1,x_counter,y_counter,z_counter,particle_counter
  real(8) :: cell_width
  cell_width = 1d0
  !--------------------------------------------------------------
  !initializing positions of particles in room as FCC lattice
  particle_counter =1
  do z_counter = 0, dim_cells-1
     do y_counter =0,dim_cells-1
        do x_counter = 0, dim_cells-1
           position(particle_counter,:) = (/ 0d0+ x_counter,0d0+ y_counter,0d0 + z_counter /)
           position(particle_counter + 1,:) = (/ 0.5d0 + x_counter, 0.5d0 + y_counter,0d0+ z_counter /)
           position(particle_counter + 2,:) = (/ 0.5d0 + x_counter,0d0+ y_counter, 0.5d0 + z_counter /)
           if (particle_counter+3 > N) then
              print *, "Out of bounds in position array"
              exit
           end if
           position(particle_counter + 3,:) = (/ 0d0+x_counter, 0.5d0 + y_counter, 0.5d0 + z_counter /)
           particle_counter = particle_counter + 4
        end do
     end do
  end do
  !---------------------------------------------------------------
  return
end subroutine position_initializer

subroutine velocity_initializer(temp, N, mass,velocity)
  real(8) :: temp, mass, stand_dev, kb
  integer :: N, counter
  real(8) :: p_x, p_y, p_z
  real(8), intent(out), dimension(N,3) :: velocity
  kb = 1d0
  p_x = 0d0
  p_y = 0d0
  p_z = 0d0
  stand_dev = sqrt(kb*mass*temp)
  !---------------------------------------------------------------
  !initializing velocities of particles according to a Gaussian distribution
  do counter =1,N
     call get_gaussian_velocity(stand_dev,velocity(counter,1))
     call get_gaussian_velocity(stand_dev,velocity(counter,2))
     call get_gaussian_velocity(stand_dev,velocity(counter,3))
  end do
  !---------------------------------------------------------------
  !---------------------------------------------------------------
  !finds total momentum in each direction and reduces to get closer to 0
  do counter=1,N
     p_x = p_x + velocity(counter,1)
     p_y = p_y + velocity(counter,2)
     p_z = p_z + velocity(counter,3)
  end do
  velocity(:,1) = velocity(:,1) - p_x/N
  velocity(:,2) = velocity(:,2) - p_y/N
  velocity(:,3) = velocity(:,3) - p_z/N
  p_x =0d0
  p_y = 0d0
  p_z = 0d0
  do counter=1,N
     p_x = p_x + velocity(counter,1)
     p_y = p_y + velocity(counter,2)
     p_z = p_z + velocity(counter,3)
  end do
  print *,p_x,p_y,p_z
  !---------------------------------------------------------------
  return
end subroutine velocity_initializer

subroutine get_gaussian_velocity(stand_dev,x)

  real(8), intent(in) :: stand_dev
  real(8), intent(out) :: x


  real(8) :: random_1,random_y, I, x_integ, dx, pi, prefactor, y, test
  integer :: counter, check
  x = 3.0d0  
  check = 0
  pi = 4d0*atan(1.0d0)
  dx = stand_dev/10d0
  x= 0d0
  prefactor = 1d0/(stand_dev)*sqrt(2d0*pi)
  do while (check == 0)
     call random_number(random_1)
     random_1 = random_1
     x_integ = -10d0*stand_dev
     I = 0d0
     do while (x_integ < 10d0*stand_dev)
        I = I + prefactor*exp(-x_integ*x_integ/(2d0*stand_dev*stand_dev))*dx
        if (I > random_1) then
           x = x_integ
           exit
        else
           x_integ = x_integ + dx
        end if
     end do
     call random_number(random_y)
     random_y = 10d0*random_y
     test = prefactor*exp(-x*x/(2d0*stand_dev*stand_dev))
     if (random_y < test) then
        return
        check = 1
     end if
  end do
end subroutine get_gaussian_velocity
