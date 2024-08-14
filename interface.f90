module interface_mod
  use datatypes
  use Variables
  use Geometry, only: NumNP, TempN, ThNew, coord
  use timedata, only: EpsTime, tMax
  implicit none
  
  integer, parameter :: interface_output_unit = 101
  character(len=*), parameter :: output_dir = 'interface'
  character(len=*), parameter :: output_file = 'interface_output.csv'
  integer, parameter :: output_frequency = 10000  ! Output every 10000th time step
  
  logical :: is_initialized = .false.
  
contains
  subroutine initialize_interface()
    character(len=256) :: command
    logical :: dir_exists
    integer :: iostat
    
    ! Check if directory exists
    inquire(file=output_dir, exist=dir_exists)
    
    if (.not. dir_exists) then
      ! Create directory
      command = "mkdir -p " // output_dir
      call system(command)
    end if
    
    ! Open file in the new directory
    open(unit=interface_output_unit, file=output_dir//'/'//output_file, status='replace', iostat=iostat)
    if (iostat /= 0) then
      print *, "Error opening file: ", output_dir//'/'//output_file
      return
    end if
    
    write(interface_output_unit, '(A)') "Time(hours),Time(days),Depth(cm),Temperature(C),SoilWaterContent"
    close(interface_output_unit)
    
    is_initialized = .true.
  end subroutine initialize_interface

  subroutine interface(t, timestep)
    real(dp), intent(in) :: t
    integer, intent(in) :: timestep
    
    ! Initialize if not already done
    if (.not. is_initialized) then
      call initialize_interface()
    end if
    
    ! Only output data every output_frequency time steps
    if (mod(timestep, output_frequency) == 0) then
      call output_data(t)
    end if
    
    ! Check if it's the last time step
    if (abs(t - tMax) <= EpsTime) then
      call output_data(t)  ! Ensure the last step is always output
    end if
  end subroutine interface

  subroutine output_data(t)
    real(dp), intent(in) :: t
    integer :: i, iostat
    real(dp) :: t_days
    
    t_days = t / 24.0_dp  ! Convert hours to days
    
    open(unit=interface_output_unit, file=output_dir//'/'//output_file, status='old', position='append', iostat=iostat)
    if (iostat /= 0) then
      print *, "Error opening file for appending: ", output_dir//'/'//output_file
      return
    end if
    
    do i = 1, NumNP
      write(interface_output_unit, '(F12.6,",",F12.6,",",F12.6,",",F12.6,",",F12.6)') &
        t, t_days, coord(i), TempN(i), ThNew(i)
    end do
    
    close(interface_output_unit)
  end subroutine output_data
end module interface_mod
