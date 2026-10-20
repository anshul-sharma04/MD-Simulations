module parameters
    implicit none
    
    ! Physical constants
    real(8), parameter :: kb = 1.380649d-23  ! Boltzmann constant (J/K)
    real(8), parameter :: avogadro = 6.02214076d23  ! Avogadro's number
    real(8), parameter :: atomic_mass_unit = 1.66053906660d-27  ! kg
    real(8), parameter :: pi = 3.14159265359d0  ! Pi constant
    
    ! Simulation parameters
    integer, parameter :: n_particles = 100   ! Number of particles
    real(8), parameter :: box_size = 10.0d0   ! Simulation box size
    real(8), parameter :: dt = 0.001d0        ! Time step
    real(8), parameter :: temperature = 300.0d0 ! Temperature in Kelvin
    integer, parameter :: n_steps = 1000      ! Number of simulation steps
    integer, parameter :: output_frequency = 10 ! Steps between outputs
    
    ! Lennard-Jones parameters (for Argon)
    real(8), parameter :: epsilon = 1.65d-21   ! Energy parameter (J)
    real(8), parameter :: sigma = 3.4d-10      ! Distance parameter (m)
    real(8), parameter :: cutoff = 2.5d0 * sigma ! Cutoff distance
    real(8), parameter :: mass = 39.95d0 * atomic_mass_unit  ! Argon mass
    
    ! Thermostat parameters
    logical, parameter :: use_thermostat = .true.
    real(8), parameter :: thermostat_coupling = 0.1d0  ! Thermostat coupling strength
    
    ! Neighbor list parameters
    logical, parameter :: use_neighbor_list = .true.
    real(8), parameter :: neighbor_skin = 0.3d0 * sigma  ! Skin distance for neighbor list
    integer, parameter :: neighbor_update_frequency = 10  ! Steps between neighbor list updates
    
    ! Analysis parameters
    logical, parameter :: calculate_rdf = .true.  ! Radial distribution function
    integer, parameter :: n_rdf_bins = 100        ! Number of RDF histogram bins
    real(8), parameter :: rdf_max_distance = 5.0d0 * sigma  ! Maximum distance for RDF
    
    ! System initialization
    character(len=20), parameter :: initial_configuration = 'fcc'  ! 'fcc', 'random', or 'file'
    character(len=100), parameter :: input_file = 'initial.xyz'    ! Input file for 'file' configuration
    
contains
    ! Function to convert reduced units to SI units
    function reduced_to_si_energy(reduced_energy) result(si_energy)
        real(8), intent(in) :: reduced_energy
        real(8) :: si_energy
        si_energy = reduced_energy * epsilon
    end function reduced_to_si_energy
    
    function reduced_to_si_distance(reduced_distance) result(si_distance)
        real(8), intent(in) :: reduced_distance
        real(8) :: si_distance
        si_distance = reduced_distance * sigma
    end function reduced_to_si_distance
    
    function reduced_to_si_time(reduced_time) result(si_time)
        real(8), intent(in) :: reduced_time
        real(8) :: si_time
        si_time = reduced_time * sqrt(mass * sigma * sigma / epsilon)
    end function reduced_to_si_time
    
end module parameters 