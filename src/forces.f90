module forces
    use parameters
    use particles
    implicit none
    
contains
    subroutine calculate_forces(system)
        type(particle_system), intent(inout) :: system
        integer :: i, j, k
        real(8) :: rij(3), r2, r6, r12, force_magnitude
        real(8) :: cutoff2, virial
        
        ! Reset forces and potential energy
        system%forces = 0.0d0
        system%potential_energy = 0.0d0
        virial = 0.0d0
        cutoff2 = cutoff * cutoff
        
        if (use_neighbor_list) then
            ! Check if neighbor list needs updating
            if (maxval(sum((system%positions - system%neighbors%positions_at_update)**2, dim=1)) > &
                (0.5d0 * neighbor_skin)**2) then
                call system%update_neighbor_list()
            end if
            
            ! Calculate forces using neighbor list
            do i = 1, n_particles
                do k = 1, system%neighbors%n_neighbors(i)
                    j = system%neighbors%list(k,i)
                    if (j > i) then  ! Avoid double counting
                        rij = system%positions(:,i) - system%positions(:,j)
                        rij = rij - box_size * nint(rij/box_size)  ! Minimum image convention
                        r2 = sum(rij*rij)
                        
                        if (r2 < cutoff2) then
                            r6 = (sigma*sigma/r2)**3
                            r12 = r6*r6
                            
                            ! Calculate force magnitude
                            force_magnitude = 24.0d0 * epsilon * (2.0d0*r12 - r6) / r2
                            
                            ! Update forces
                            system%forces(:,i) = system%forces(:,i) + force_magnitude * rij
                            system%forces(:,j) = system%forces(:,j) - force_magnitude * rij
                            
                            ! Update potential energy and virial
                            system%potential_energy = system%potential_energy + &
                                4.0d0 * epsilon * (r12 - r6)
                            virial = virial + force_magnitude * sum(rij*rij)
                        end if
                    end if
                end do
            end do
        else
            ! Calculate forces between all pairs
            do i = 1, n_particles-1
                do j = i+1, n_particles
                    rij = system%positions(:,i) - system%positions(:,j)
                    rij = rij - box_size * nint(rij/box_size)
                    r2 = sum(rij*rij)
                    
                    if (r2 < cutoff2) then
                        r6 = (sigma*sigma/r2)**3
                        r12 = r6*r6
                        
                        force_magnitude = 24.0d0 * epsilon * (2.0d0*r12 - r6) / r2
                        
                        system%forces(:,i) = system%forces(:,i) + force_magnitude * rij
                        system%forces(:,j) = system%forces(:,j) - force_magnitude * rij
                        
                        system%potential_energy = system%potential_energy + &
                            4.0d0 * epsilon * (r12 - r6)
                        virial = virial + force_magnitude * sum(rij*rij)
                    end if
                end do
            end do
        end if
        
        ! Add long-range corrections
        call add_long_range_corrections(system, virial)
        
        ! Calculate pressure
        system%pressure = (n_particles * kb * system%temperature / (box_size**3)) + &
                         virial / (3.0d0 * box_size**3)
    end subroutine calculate_forces
    
    subroutine add_long_range_corrections(system, virial)
        type(particle_system), intent(inout) :: system
        real(8), intent(inout) :: virial
        real(8) :: r_cut3, r_cut9, density
        
        ! Calculate density
        density = n_particles / (box_size**3)
        
        ! Calculate corrections for potential energy and pressure
        r_cut3 = (sigma/cutoff)**3
        r_cut9 = r_cut3**3
        
        ! Potential energy correction
        system%potential_energy = system%potential_energy + &
            (8.0d0/3.0d0) * pi * density * n_particles * epsilon * sigma**3 * &
            (r_cut9/3.0d0 - r_cut3)
        
        ! Virial correction
        virial = virial + &
            (32.0d0/3.0d0) * pi * density * n_particles * epsilon * sigma**3 * &
            (r_cut9 - r_cut3)
    end subroutine add_long_range_corrections
    
    subroutine calculate_kinetic_energy(system)
        type(particle_system), intent(inout) :: system
        
        system%kinetic_energy = 0.5d0 * mass * sum(system%velocities * system%velocities)
        call system%calculate_temperature()
    end subroutine calculate_kinetic_energy
    
end module forces 