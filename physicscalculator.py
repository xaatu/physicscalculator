import math


# ROUGHLY TESTED FOR ERRORS AS I WAS MAKING IT BUT MAY NOT BE 100% PERFECT YET


# PHYSICS CALC FUNCTIONS
def calculate_force(mass, acceleration):
    return mass * acceleration

def calculate_work(force, distance):
    return force * distance

def calculate_kinetic_energy(mass, velocity):
    return 0.5 * mass * velocity ** 2

def calculate_potential_energy(mass, height, gravity=9.81):
    return mass * gravity * height

def calculate_projectile_range(velocity, angle, gravity=9.81):
    angle_rad = math.radians(angle)
    return (velocity ** 2 * math.sin(2 * angle_rad)) / gravity

def calculate_projectile_time_of_flight(velocity, angle, gravity=9.81):
    angle_rad = math.radians(angle)
    return (2 * velocity * math.sin(angle_rad)) / gravity

def calculate_centripetal_force(mass, velocity, radius):
    return (mass * velocity ** 2) / radius

def calculate_thermodynamic_work(pressure, volume_change):
    return pressure * volume_change

def calculate_lens_magnification(image_distance, object_distance):
    return image_distance / object_distance

def calculate_electric_force(charge1, charge2, distance):
    k = 8.9875e9  # COULOMBS CONSTANT IN N m²/C²
    return k * (charge1 * charge2) / distance ** 2

def calculate_magnetic_force(charge, velocity, magnetic_field):
    return charge * velocity * magnetic_field

def calculate_electric_field(charge, distance):
    k = 8.9875e9  # COULOMBS CONSTANT IN N m²/C²
    return k * charge / distance ** 2

def calculate_relativistic_energy(mass, velocity):
    c = 3.0e8  # SOL IN m/s
    gamma = 1 / math.sqrt(1 - (velocity ** 2 / c ** 2))
    return mass * c ** 2 * gamma

def calculate_relativistic_momentum(mass, velocity):
    c = 3.0e8  # SOL IN m/s
    gamma = 1 / math.sqrt(1 - (velocity ** 2 / c ** 2))
    return mass * velocity * gamma

def calculate_photon_energy(wavelength):
    h = 6.626e-34  # PLANCK'S CONSTANT IN J·s
    c = 3.0e8  # SOL IN m/s
    return h * c / wavelength

def calculate_de_broglie_wavelength(mass, velocity):
    h = 6.626e-34  # PLANCK'S CONSTANT IN J·s
    return h / (mass * velocity)

def calculate_fluid_pressure(force, area):
    return force / area

def calculate_buoyant_force(density, volume, gravity=9.81):
    return density * volume * gravity

def calculate_schwarzschild_radius(mass):
    G = 6.67430e-11  # GRAVITATIONAL CONSTANT m³/kg/s²
    c = 3.0e8  # SOL IN m/s
    return 2 * G * mass / c ** 2

def calculate_luminosity(apparent_magnitude, distance):
    # USING DISTANCE MODULUS FORMULA
    M_sun = 4.83  # ABS MAGNITUDE OF THE SUN
    d_pc = distance * 3.086e16  # DISTANCE FROM LIGHT YEARS TO M
    return 3.846e26 * 10 ** ((apparent_magnitude - M_sun) / -2.5)

def calculate_stellar_distance(apparent_magnitude, luminosity):
    # INVERSE SQUARE LAW FOR LUMINOSITY
    M_sun = 4.83  # ABS MAGNITUDE OF SUN
    return 10 ** ((apparent_magnitude - M_sun) / -5) * 3.086e16  # CONVERT TO M

def calculate_decay_constant(half_life):
    return math.log(2) / half_life

def calculate_radioactive_decay(initial_activity, decay_constant, time):
    return initial_activity * math.exp(-decay_constant * time)

def calculate_stress(force, area):
    return force / area

def calculate_strain(original_length, deformed_length):
    return (deformed_length - original_length) / original_length

def calculate_youngs_modulus(stress, strain):
    return stress / strain

def calculate_wave_speed(frequency, wavelength):
    return frequency * wavelength

def calculate_wave_frequency(speed, wavelength):
    return speed / wavelength

def calculate_wave_wavelength(speed, frequency):
    return speed / frequency

def display_menu():
    print("\nAdvanced Physics Calculator")
    print("1. Calculate Force (F = m * a)")
    print("2. Calculate Work (W = F * d)")
    print("3. Calculate Kinetic Energy (KE = 0.5 * m * v^2)")
    print("4. Calculate Potential Energy (PE = m * g * h)")
    print("5. Calculate Projectile Range")
    print("6. Calculate Projectile Time of Flight")
    print("7. Calculate Centripetal Force (F_c = m * v^2 / r)")
    print("8. Calculate Thermodynamic Work (W = P * ΔV)")
    print("9. Calculate Lens Magnification (M = i / o)")
    print("10. Calculate Electric Force (F = k * q1 * q2 / r^2)")
    print("11. Calculate Magnetic Force (F = q * v * B)")
    print("12. Calculate Electric Field (E = k * q / r^2)")
    print("13. Calculate Relativistic Energy (E = m * c^2 / sqrt(1 - v^2 / c^2))")
    print("14. Calculate Relativistic Momentum (p = m * v / sqrt(1 - v^2 / c^2))")
    print("15. Calculate Photon Energy (E = h * c / λ)")
    print("16. Calculate de Broglie Wavelength (λ = h / (m * v))")
    print("17. Calculate Fluid Pressure (P = F / A)")
    print("18. Calculate Buoyant Force (F_b = ρ * V * g)")
    print("19. Calculate Schwarzschild Radius (R_s = 2 * G * m / c^2)")
    print("20. Calculate Luminosity")
    print("21. Calculate Stellar Distance")
    print("22. Calculate Decay Constant (λ = ln(2) / t_1/2)")
    print("23. Calculate Radioactive Decay")
    print("24. Calculate Stress (σ = F / A)")
    print("25. Calculate Strain (ε = ΔL / L₀)")
    print("26. Calculate Young's Modulus (E = σ / ε)")
    print("27. Calculate Wave Speed (v = f * λ)")
    print("28. Calculate Wave Frequency (f = v / λ)")
    print("29. Calculate Wave Wavelength (λ = v / f)")
    print("30. Exit")

def main():
    while True:
        display_menu()
        choice = input("Choose an option (1-30): ")

        if choice == '1':
            mass = float(input("Enter mass (kg): "))
            acceleration = float(input("Enter acceleration (m/s^2): "))
            force = calculate_force(mass, acceleration)
            print(f"Force = {force:.2f} N")

        elif choice == '2':
            force = float(input("Enter force (N): "))
            distance = float(input("Enter distance (m): "))
            work = calculate_work(force, distance)
            print(f"Work = {work:.2f} J")

        elif choice == '3':
            mass = float(input("Enter mass (kg): "))
            velocity = float(input("Enter velocity (m/s): "))
            kinetic_energy = calculate_kinetic_energy(mass, velocity)
            print(f"Kinetic Energy = {kinetic_energy:.2f} J")

        elif choice == '4':
            mass = float(input("Enter mass (kg): "))
            height = float(input("Enter height (m): "))
            gravity = float(input("Enter gravity (m/s^2) [default 9.81]: ") or 9.81)
            potential_energy = calculate_potential_energy(mass, height, gravity)
            print(f"Potential Energy = {potential_energy:.2f} J")

        elif choice == '5':
            velocity = float(input("Enter initial velocity (m/s): "))
            angle = float(input("Enter launch angle (degrees): "))
            gravity = float(input("Enter gravity (m/s^2) [default 9.81]: ") or 9.81)
            range_ = calculate_projectile_range(velocity, angle, gravity)
            print(f"Projectile Range = {range_:.2f} m")

        elif choice == '6':
            velocity = float(input("Enter initial velocity (m/s): "))
            angle = float(input("Enter launch angle (degrees): "))
            gravity = float(input("Enter gravity (m/s^2) [default 9.81]: ") or 9.81)
            time_of_flight = calculate_projectile_time_of_flight(velocity, angle, gravity)
            print(f"Time of Flight = {time_of_flight:.2f} s")

        elif choice == '7':
            mass = float(input("Enter mass (kg): "))
            velocity = float(input("Enter velocity (m/s): "))
            radius = float(input("Enter radius (m): "))
            centripetal_force = calculate_centripetal_force(mass, velocity, radius)
            print(f"Centripetal Force = {centripetal_force:.2f} N")

        elif choice == '8':
            pressure = float(input("Enter pressure (Pa): "))
            volume_change = float(input("Enter volume change (m^3): "))
            work = calculate_thermodynamic_work(pressure, volume_change)
            print(f"Thermodynamic Work = {work:.2f} J")

        elif choice == '9':
            image_distance = float(input("Enter image distance (m): "))
            object_distance = float(input("Enter object distance (m): "))
            magnification = calculate_lens_magnification(image_distance, object_distance)
            print(f"Lens Magnification = {magnification:.2f}")

        elif choice == '10':
            charge1 = float(input("Enter charge 1 (C): "))
            charge2 = float(input("Enter charge 2 (C): "))
            distance = float(input("Enter distance between charges (m): "))
            electric_force = calculate_electric_force(charge1, charge2, distance)
            print(f"Electric Force = {electric_force:.2f} N")

        elif choice == '11':
            charge = float(input("Enter charge (C): "))
            velocity = float(input("Enter velocity (m/s): "))
            magnetic_field = float(input("Enter magnetic field strength (T): "))
            magnetic_force = calculate_magnetic_force(charge, velocity, magnetic_field)
            print(f"Magnetic Force = {magnetic_force:.2f} N")

        elif choice == '12':
            charge = float(input("Enter charge (C): "))
            distance = float(input("Enter distance from the charge (m): "))
            electric_field = calculate_electric_field(charge, distance)
            print(f"Electric Field = {electric_field:.2f} N/C")

        elif choice == '13':
            mass = float(input("Enter mass (kg): "))
            velocity = float(input("Enter velocity (m/s): "))
            energy = calculate_relativistic_energy(mass, velocity)
            print(f"Relativistic Energy = {energy:.2e} J")

        elif choice == '14':
            mass = float(input("Enter mass (kg): "))
            velocity = float(input("Enter velocity (m/s): "))
            momentum = calculate_relativistic_momentum(mass, velocity)
            print(f"Relativistic Momentum = {momentum:.2e} kg·m/s")

        elif choice == '15':
            wavelength = float(input("Enter wavelength (m): "))
            photon_energy = calculate_photon_energy(wavelength)
            print(f"Photon Energy = {photon_energy:.2e} J")

        elif choice == '16':
            mass = float(input("Enter mass (kg): "))
            velocity = float(input("Enter velocity (m/s): "))
            de_broglie_wavelength = calculate_de_broglie_wavelength(mass, velocity)
            print(f"de Broglie Wavelength = {de_broglie_wavelength:.2e} m")

        elif choice == '17':
            force = float(input("Enter force (N): "))
            area = float(input("Enter area (m^2): "))
            fluid_pressure = calculate_fluid_pressure(force, area)
            print(f"Fluid Pressure = {fluid_pressure:.2f} Pa")

        elif choice == '18':
            density = float(input("Enter fluid density (kg/m^3): "))
            volume = float(input("Enter volume (m^3): "))
            gravity = float(input("Enter gravity (m/s^2) [default 9.81]: ") or 9.81)
            buoyant_force = calculate_buoyant_force(density, volume, gravity)
            print(f"Buoyant Force = {buoyant_force:.2f} N")

        elif choice == '19':
            mass = float(input("Enter mass of the black hole (kg): "))
            schwarzschild_radius = calculate_schwarzschild_radius(mass)
            print(f"Schwarzschild Radius = {schwarzschild_radius:.2e} m")

        elif choice == '20':
            apparent_magnitude = float(input("Enter apparent magnitude: "))
            distance = float(input("Enter distance (light-years): "))
            luminosity = calculate_luminosity(apparent_magnitude, distance)
            print(f"Luminosity = {luminosity:.2e} W")

        elif choice == '21':
            apparent_magnitude = float(input("Enter apparent magnitude: "))
            luminosity = float(input("Enter luminosity (W): "))
            distance = calculate_stellar_distance(apparent_magnitude, luminosity)
            print(f"Stellar Distance = {distance:.2e} m")

        elif choice == '22':
            half_life = float(input("Enter half-life (s): "))
            decay_constant = calculate_decay_constant(half_life)
            print(f"Decay Constant = {decay_constant:.2e} s^-1")

        elif choice == '23':
            initial_activity = float(input("Enter initial activity (Bq): "))
            decay_constant = float(input("Enter decay constant (s^-1): "))
            time = float(input("Enter time (s): "))
            activity = calculate_radioactive_decay(initial_activity, decay_constant, time)
            print(f"Activity after {time:.2f} s = {activity:.2e} Bq")

        elif choice == '24':
            force = float(input("Enter force (N): "))
            area = float(input("Enter area (m^2): "))
            stress = calculate_stress(force, area)
            print(f"Stress = {stress:.2f} Pa")

        elif choice == '25':
            original_length = float(input("Enter original length (m): "))
            deformed_length = float(input("Enter deformed length (m): "))
            strain = calculate_strain(original_length, deformed_length)
            print(f"Strain = {strain:.2f}")

        elif choice == '26':
            stress = float(input("Enter stress (Pa): "))
            strain = float(input("Enter strain: "))
            youngs_modulus = calculate_youngs_modulus(stress, strain)
            print(f"Young's Modulus = {youngs_modulus:.2e} Pa")

        elif choice == '27':
            frequency = float(input("Enter frequency (Hz): "))
            wavelength = float(input("Enter wavelength (m): "))
            wave_speed = calculate_wave_speed(frequency, wavelength)
            print(f"Wave Speed = {wave_speed:.2e} m/s")

        elif choice == '28':
            speed = float(input("Enter wave speed (m/s): "))
            wavelength = float(input("Enter wavelength (m): "))
            frequency = calculate_wave_frequency(speed, wavelength)
            print(f"Frequency = {frequency:.2e} Hz")

        elif choice == '29':
            speed = float(input("Enter wave speed (m/s): "))
            frequency = float(input("Enter frequency (Hz): "))
            wavelength = calculate_wave_wavelength(speed, frequency)
            print(f"Wavelength = {wavelength:.2e} m")

        elif choice == '30':
            print("Exiting...")
            break

        else:
            print("Invalid choice, please try again.")

if __name__ == "__main__":
    main()
