module io
    use parameters
    use particles
    implicit none
    
    private
    public :: write_trajectory, write_energies, write_rdf, write_initial_info
    
    character(len=*), parameter :: outdir = 'output'
    logical :: dir_exists = .false.
    
contains
    subroutine ensure_output_directory()
        ! Local variables for directory operations
        logical :: dir_check
        integer :: stat
        character(len=256) :: cmd
        
        ! Check if directory exists
        inquire(file=trim(outdir)//'/.', exist=dir_check)
        
        if (.not. dir_check) then
            ! Create directory using system command
#ifdef _WIN32
            cmd = 'if not exist "'//trim(outdir)//'" mkdir "'//trim(outdir)//'"'
#else
            cmd = 'mkdir '//trim(outdir)
#endif
            call execute_command_line(trim(cmd), wait=.true., exitstat=stat)
            if (stat /= 0) then
                write(*,*) 'Error: Failed to create output directory'
                stop
            end if
        end if
        
        dir_exists = .true.
    end subroutine ensure_output_directory

    subroutine write_initial_info()
        ! Ensure output directory exists
        if (.not. dir_exists) call ensure_output_directory()
        
        ! Write simulation parameters to a file
        open(unit=12, file=trim(outdir)//'/simulation_info.txt', status='replace')
        write(12,*) '=== Molecular Dynamics Simulation Parameters ==='
        write(12,*) 'Number of particles:', n_particles
        write(12,*) 'Box size:', box_size
        write(12,*) 'Temperature:', temperature
        write(12,*) 'Time step:', dt
        write(12,*) 'Number of steps:', n_steps
        write(12,*) 'Output frequency:', output_frequency
        write(12,*) 'Initial configuration:', initial_configuration
        write(12,*) ''
        write(12,*) '=== Lennard-Jones Parameters ==='
        write(12,*) 'Epsilon:', epsilon
        write(12,*) 'Sigma:', sigma
        write(12,*) 'Cutoff distance:', cutoff
        write(12,*) ''
        write(12,*) '=== Additional Features ==='
        write(12,*) 'Thermostat:', use_thermostat
        if (use_thermostat) then
            write(12,*) 'Thermostat coupling:', thermostat_coupling
        end if
        write(12,*) 'Neighbor list:', use_neighbor_list
        if (use_neighbor_list) then
            write(12,*) 'Neighbor skin distance:', neighbor_skin
            write(12,*) 'Neighbor update frequency:', neighbor_update_frequency
        end if
        write(12,*) 'RDF calculation:', calculate_rdf
        if (calculate_rdf) then
            write(12,*) 'RDF bins:', n_rdf_bins
            write(12,*) 'RDF max distance:', rdf_max_distance
        end if
        close(12)
        
        ! Create header for energies file
        open(unit=11, file=trim(outdir)//'/energies.dat', status='replace')
        write(11,'(A)') '# Step PotentialEnergy KineticEnergy TotalEnergy Temperature Pressure'
        close(11)
    end subroutine write_initial_info
    
    subroutine write_trajectory(system, step)
        type(particle_system), intent(in) :: system
        integer, intent(in) :: step
        integer :: i
        character(len=100) :: comment
        
        ! Ensure output directory exists
        if (.not. dir_exists) call ensure_output_directory()
        
        ! Write positions to trajectory file
        open(unit=10, file=trim(outdir)//'/trajectory.xyz', position='append')
        write(10,'(I0)') n_particles
        write(comment,'(A,I0,A,E16.8,A,F10.3,A,E16.8)') &
            'Step: ', step, ' E_tot: ', &
            system%potential_energy + system%kinetic_energy, &
            ' T: ', system%temperature, &
            ' P: ', system%pressure
        write(10,'(A)') trim(comment)
        
        do i = 1, n_particles
            write(10,'(A,3E16.8)') 'Ar', system%positions(:,i)
        end do
        
        close(10)
    end subroutine write_trajectory
    
    subroutine write_energies(system, step)
        type(particle_system), intent(in) :: system
        integer, intent(in) :: step
        
        ! Ensure output directory exists
        if (.not. dir_exists) call ensure_output_directory()
        
        ! Write energies to energy file
        open(unit=11, file=trim(outdir)//'/energies.dat', position='append')
        write(11,'(I10,5E16.8)') step, &
            system%potential_energy, &
            system%kinetic_energy, &
            system%potential_energy + system%kinetic_energy, &
            system%temperature, &
            system%pressure
        close(11)
    end subroutine write_energies
    
    subroutine write_rdf(system, step)
        type(particle_system), intent(in) :: system
        integer, intent(in) :: step
        integer :: i
        real(8) :: r, dr
        character(len=100) :: filename
        
        if (.not. calculate_rdf) return
        
        ! Ensure output directory exists
        if (.not. dir_exists) call ensure_output_directory()
        
        ! Write RDF to file
        write(filename,'(A,A,I0,A)') trim(outdir), '/rdf_', step, '.dat'
        open(unit=13, file=filename, status='replace')
        write(13,'(A)') '# r g(r)'
        
        dr = rdf_max_distance / dble(n_rdf_bins)
        do i = 1, n_rdf_bins
            r = (dble(i)-0.5d0) * dr
            write(13,'(2E16.8)') r, system%rdf(i)
        end do
        
        close(13)
    end subroutine write_rdf
    
end module io 