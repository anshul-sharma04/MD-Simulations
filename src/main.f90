program molecular_dynamics
    use parameters
    use particles
    use forces
    use integrator
    use io
    implicit none
    
    type(particle_system) :: system
    integer :: step
    integer :: start_time, end_time, count_rate, count_max
    real(8) :: elapsed_time
    real(8) :: initial_energy, energy_drift
    
    ! Get CPU time rate
    call system_clock(count_rate=count_rate, count_max=count_max)
    
    ! Record start time
    call system_clock(count=start_time)
    
    ! Write initial simulation information
    call write_initial_info()
    
    ! Initialize random number generator
    call random_seed()
    
    ! Initialize the system
    call system%init()
    
    ! Calculate initial forces
    call calculate_forces(system)
    
    ! Calculate initial kinetic energy and temperature
    call calculate_kinetic_energy(system)
    
    ! Store initial energy
    initial_energy = system%potential_energy + system%kinetic_energy
    
    ! Write initial configuration
    call write_trajectory(system, 0)
    call write_energies(system, 0)
    if (calculate_rdf) call write_rdf(system, 0)
    
    ! Print initial state
    print '(A)', '=== Molecular Dynamics Simulation Started ==='
    print '(A,I0)', 'Number of particles: ', n_particles
    print '(A,F10.3)', 'Initial temperature: ', system%temperature
    print '(A,F10.3)', 'Initial energy: ', initial_energy
    print '(A)', 'Step     Temperature     Energy     Energy Drift     Progress'
    
    ! Main simulation loop
    do step = 1, n_steps
        ! Integrate equations of motion
        call velocity_verlet_step(system)
        
        ! Calculate energy drift
        energy_drift = abs((system%potential_energy + system%kinetic_energy - initial_energy) / initial_energy)
        
        ! Write output at specified frequency
        if (mod(step, output_frequency) == 0) then
            call write_trajectory(system, step)
            call write_energies(system, step)
            if (calculate_rdf) call write_rdf(system, step)
            
            ! Print progress
            print '(I6,F14.3,F12.3,E16.3,F10.1,A)', &
                step, &
                system%temperature, &
                system%potential_energy + system%kinetic_energy, &
                energy_drift, &
                100.0 * step / n_steps, '%'
        end if
    end do
    
    ! Record end time
    call system_clock(count=end_time)
    
    ! Calculate elapsed time in seconds
    elapsed_time = real(end_time - start_time, kind=8) / real(count_rate, kind=8)
    
    ! Print final statistics
    print '(A)', ''
    print '(A)', '=== Simulation Complete ==='
    print '(A,F10.3,A)', 'Elapsed time: ', elapsed_time, ' seconds'
    print '(A,F10.3)', 'Final temperature: ', system%temperature
    print '(A,F10.3)', 'Final energy: ', system%potential_energy + system%kinetic_energy
    print '(A,E12.3)', 'Energy drift: ', energy_drift
    print '(A,F10.3)', 'Final pressure: ', system%pressure
    
    ! Print output files information
    print '(A)', ''
    print '(A)', 'Output files generated:'
    print '(A)', '- simulation_info.txt: Simulation parameters and settings'
    print '(A)', '- trajectory.xyz: Particle positions (can be visualized with VMD/OVITO)'
    print '(A)', '- energies.dat: Time evolution of energies, temperature, and pressure'
    if (calculate_rdf) then
        print '(A)', '- rdf_*.dat: Radial distribution function at output steps'
    end if
    
end program molecular_dynamics 