module integrator
    use parameters
    use particles
    use forces
    implicit none
    
contains
    subroutine velocity_verlet_step(system)
        type(particle_system), intent(inout) :: system
        real(8) :: old_forces(3,n_particles)
        
        ! Store old forces
        old_forces = system%forces
        
        ! First half of velocity update
        system%velocities = system%velocities + 0.5d0 * system%forces * dt / mass
        
        ! Update positions
        system%positions = system%positions + system%velocities * dt
        
        ! Apply periodic boundary conditions
        call system%apply_periodic_boundaries()
        
        ! Calculate new forces
        call calculate_forces(system)
        
        ! Second half of velocity update
        system%velocities = system%velocities + 0.5d0 * system%forces * dt / mass
        
        ! Calculate kinetic energy and temperature
        call calculate_kinetic_energy(system)
        
        ! Apply thermostat if enabled
        if (use_thermostat) then
            call apply_thermostat(system)
        end if
        
        ! Calculate RDF if enabled
        if (calculate_rdf) then
            call system%calculate_rdf()
        end if
    end subroutine velocity_verlet_step
    
    subroutine apply_thermostat(system)
        type(particle_system), intent(inout) :: system
        real(8) :: scaling_factor
        
        ! Berendsen thermostat
        scaling_factor = sqrt(1.0d0 + (dt/thermostat_coupling) * &
            (temperature/system%temperature - 1.0d0))
        
        system%velocities = system%velocities * scaling_factor
        
        ! Recalculate kinetic energy and temperature
        call calculate_kinetic_energy(system)
    end subroutine apply_thermostat
    
end module integrator 