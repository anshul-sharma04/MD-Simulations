module particles
    use parameters
    implicit none
    
    private
    public :: particle_system
    
    type :: neighbor_list_type
        integer, allocatable :: list(:,:)  ! List of neighbors for each particle
        integer, allocatable :: n_neighbors(:)  ! Number of neighbors for each particle
        real(8), allocatable :: positions_at_update(:,:)  ! Positions when list was last updated
        integer :: last_update = 0  ! Step of last update
    end type neighbor_list_type
    
    type :: particle_system
        real(8), allocatable :: positions(:,:)  ! Particle positions (3,n_particles)
        real(8), allocatable :: velocities(:,:) ! Particle velocities (3,n_particles)
        real(8), allocatable :: forces(:,:)     ! Forces on particles (3,n_particles)
        real(8) :: potential_energy             ! Total potential energy
        real(8) :: kinetic_energy              ! Total kinetic energy
        real(8) :: temperature                 ! Current temperature
        real(8) :: pressure                    ! Current pressure
        type(neighbor_list_type) :: neighbors  ! Neighbor list
        real(8), allocatable :: rdf(:)         ! Radial distribution function
    contains
        procedure :: init => initialize_system
        procedure :: apply_periodic_boundaries
        procedure :: calculate_temperature
        procedure :: calculate_pressure
        procedure :: update_neighbor_list
        procedure :: calculate_rdf => calculate_rdf_impl  ! Renamed to avoid conflict
    end type particle_system
    
contains
    
    subroutine initialize_system(this)
        class(particle_system), intent(inout) :: this
        
        ! Allocate arrays
        allocate(this%positions(3,n_particles))
        allocate(this%velocities(3,n_particles))
        allocate(this%forces(3,n_particles))
        if (calculate_rdf) allocate(this%rdf(n_rdf_bins))
        
        ! Initialize neighbor list if used
        if (use_neighbor_list) then
            allocate(this%neighbors%list(n_particles, n_particles))
            allocate(this%neighbors%n_neighbors(n_particles))
            allocate(this%neighbors%positions_at_update(3,n_particles))
            this%neighbors%n_neighbors = 0
            this%neighbors%last_update = 0
        end if
        
        ! Initialize positions based on selected configuration
        select case (initial_configuration)
        case ('fcc')
            call initialize_fcc_lattice(this)
        case ('random')
            call initialize_random_positions(this)
        case ('file')
            call read_positions_from_file(this)
        end select
        
        ! Initialize velocities from Maxwell-Boltzmann distribution
        call initialize_velocities(this)
        
        ! Initialize forces to zero
        this%forces = 0.0d0
        this%potential_energy = 0.0d0
        this%kinetic_energy = 0.0d0
        this%temperature = 0.0d0
        this%pressure = 0.0d0
        if (calculate_rdf) this%rdf = 0.0d0
    end subroutine initialize_system
    
    subroutine initialize_fcc_lattice(this)
        class(particle_system), intent(inout) :: this
        integer :: n, ix, iy, iz, i
        real(8) :: lattice_constant, pos(3)
        
        ! Calculate lattice constant for FCC structure
        n = nint((real(n_particles)/4.0d0)**(1.0d0/3.0d0))  ! Number of unit cells in each dimension
        lattice_constant = box_size / dble(n)
        
        i = 1
        do ix = 0, n-1
            do iy = 0, n-1
                do iz = 0, n-1
                    if (i > n_particles) exit
                    ! Add particle at lattice point
                    pos = [dble(ix), dble(iy), dble(iz)] * lattice_constant
                    this%positions(:,i) = pos
                    i = i + 1
                    
                    if (i > n_particles) exit
                    ! Add particle at face center (xy face)
                    pos = [dble(ix)+0.5d0, dble(iy)+0.5d0, dble(iz)] * lattice_constant
                    this%positions(:,i) = pos
                    i = i + 1
                    
                    if (i > n_particles) exit
                    ! Add particle at face center (xz face)
                    pos = [dble(ix)+0.5d0, dble(iy), dble(iz)+0.5d0] * lattice_constant
                    this%positions(:,i) = pos
                    i = i + 1
                    
                    if (i > n_particles) exit
                    ! Add particle at face center (yz face)
                    pos = [dble(ix), dble(iy)+0.5d0, dble(iz)+0.5d0] * lattice_constant
                    this%positions(:,i) = pos
                    i = i + 1
                end do
                if (i > n_particles) exit
            end do
            if (i > n_particles) exit
        end do
    end subroutine initialize_fcc_lattice
    
    subroutine initialize_random_positions(this)
        class(particle_system), intent(inout) :: this
        integer :: i, j
        real(8) :: pos(3), min_distance
        logical :: overlap
        
        do i = 1, n_particles
            do
                overlap = .false.
                call random_number(pos)
                pos = pos * box_size
                
                ! Check for overlap with existing particles
                do j = 1, i-1
                    min_distance = sqrt(sum((pos - this%positions(:,j))**2))
                    if (min_distance < sigma) then
                        overlap = .true.
                        exit
                    end if
                end do
                
                if (.not. overlap) then
                    this%positions(:,i) = pos
                    exit
                end if
            end do
        end do
    end subroutine initialize_random_positions
    
    subroutine read_positions_from_file(this)
        class(particle_system), intent(inout) :: this
        integer :: i, io_stat, n_atoms
        character(len=80) :: comment
        
        open(unit=20, file=input_file, status='old', action='read')
        read(20,*) n_atoms
        if (n_atoms /= n_particles) then
            print *, "Error: Number of particles in file does not match n_particles"
            stop
        end if
        read(20,*) comment
        
        do i = 1, n_particles
            read(20,*) comment, this%positions(:,i)
        end do
        
        close(20)
    end subroutine read_positions_from_file
    
    subroutine initialize_velocities(this)
        class(particle_system), intent(inout) :: this
        integer :: i, j
        real(8) :: velocity_scale, total_momentum(3)
        
        ! Initialize random velocities
        do i = 1, n_particles
            do j = 1, 3
                call random_number(this%velocities(j,i))
                this%velocities(j,i) = this%velocities(j,i) - 0.5d0
            end do
        end do
        
        ! Remove center of mass motion
        total_momentum = sum(this%velocities, dim=2)
        this%velocities = this%velocities - spread(total_momentum/n_particles, 2, n_particles)
        
        ! Scale velocities to match desired temperature
        call this%calculate_temperature()
        velocity_scale = sqrt(temperature/this%temperature)
        this%velocities = this%velocities * velocity_scale
    end subroutine initialize_velocities
    
    subroutine apply_periodic_boundaries(this)
        class(particle_system), intent(inout) :: this
        integer :: i
        
        do i = 1, n_particles
            this%positions(:,i) = modulo(this%positions(:,i), box_size)
        end do
    end subroutine apply_periodic_boundaries
    
    subroutine calculate_temperature(this)
        class(particle_system), intent(inout) :: this
        this%temperature = sum(this%velocities**2) * mass / (3.0d0 * n_particles * kb)
    end subroutine calculate_temperature
    
    subroutine calculate_pressure(this)
        class(particle_system), intent(inout) :: this
        this%pressure = (n_particles * kb * this%temperature / (box_size**3)) + &
                       (sum(this%positions * this%forces) / (3.0d0 * box_size**3))
    end subroutine calculate_pressure
    
    subroutine update_neighbor_list(this)
        class(particle_system), intent(inout) :: this
        integer :: i, j
        real(8) :: rij(3), r2
        real(8) :: cutoff_with_skin
        
        cutoff_with_skin = cutoff + neighbor_skin
        this%neighbors%n_neighbors = 0
        
        do i = 1, n_particles-1
            do j = i+1, n_particles
                rij = this%positions(:,i) - this%positions(:,j)
                rij = rij - box_size * nint(rij/box_size)  ! Minimum image convention
                r2 = sum(rij*rij)
                
                if (r2 < cutoff_with_skin**2) then
                    this%neighbors%n_neighbors(i) = this%neighbors%n_neighbors(i) + 1
                    this%neighbors%n_neighbors(j) = this%neighbors%n_neighbors(j) + 1
                    this%neighbors%list(this%neighbors%n_neighbors(i),i) = j
                    this%neighbors%list(this%neighbors%n_neighbors(j),j) = i
                end if
            end do
        end do
        
        this%neighbors%positions_at_update = this%positions
        this%neighbors%last_update = this%neighbors%last_update + 1
    end subroutine update_neighbor_list
    
    subroutine calculate_rdf_impl(this)  ! Renamed from calculate_rdf to avoid conflict
        class(particle_system), intent(inout) :: this
        integer :: i, j, bin
        real(8) :: rij(3), r, dr, volume, norm
        
        if (.not. calculate_rdf) return
        
        dr = rdf_max_distance / dble(n_rdf_bins)
        this%rdf = 0.0d0
        
        do i = 1, n_particles-1
            do j = i+1, n_particles
                rij = this%positions(:,i) - this%positions(:,j)
                rij = rij - box_size * nint(rij/box_size)
                r = sqrt(sum(rij*rij))
                
                if (r < rdf_max_distance) then
                    bin = int(r/dr) + 1
                    if (bin <= n_rdf_bins) then
                        this%rdf(bin) = this%rdf(bin) + 2.0d0  ! Count each pair twice
                    end if
                end if
            end do
        end do
        
        ! Normalize RDF
        do i = 1, n_rdf_bins
            r = (dble(i)-0.5d0) * dr
            volume = 4.0d0 * pi * (r**2) * dr
            norm = (dble(n_particles) * volume) / (box_size**3)
            this%rdf(i) = this%rdf(i) / (dble(n_particles) * norm)
        end do
    end subroutine calculate_rdf_impl
    
end module particles 