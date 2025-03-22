--[[
RebelCALC Engineering Modules - Lua Script Examples

This script demonstrates how to use the RebelCALC engineering modules from Lua scripts.
It includes examples for the CAD, Physics, and Electrical modules.
]]

-- Print a section header
function print_header(title)
    print(string.rep("=", 80))
    print("  " .. title)
    print(string.rep("=", 80))
end

-- Print a subsection header
function print_subheader(title)
    print("\n" .. title)
    print(string.rep("-", #title))
end

-- Format a 3D point for display
function format_point(point)
    return string.format("(%.2f, %.2f, %.2f)", point.x, point.y, point.z)
end

-- Format a complex number for display
function format_complex(complex)
    if complex.imag >= 0 then
        return string.format("%.2f + %.2fj", complex.real, complex.imag)
    else
        return string.format("%.2f - %.2fj", complex.real, -complex.imag)
    end
end

-- CAD Module Examples
function cad_examples()
    print_header("CAD Module Examples")
    
    -- Create points
    print_subheader("Creating and manipulating 3D points")
    local p1 = cad.Point3D.new(1, 2, 3)
    local p2 = cad.Point3D.new(4, 5, 6)
    local p3 = cad.Point3D.new(7, 8, 9)
    
    print("Point 1: " .. format_point(p1))
    print("Point 2: " .. format_point(p2))
    print("Point 3: " .. format_point(p3))
    
    -- Point operations
    print("\nPoint operations:")
    local sum = p1:add(p2)
    local diff = p2:subtract(p1)
    local scaled = p1:multiply(2)
    
    print("Point 1 + Point 2 = " .. format_point(sum))
    print("Point 2 - Point 1 = " .. format_point(diff))
    print("Point 1 * 2 = " .. format_point(scaled))
    
    -- Vector operations
    print("\nVector operations:")
    local dot = p1:dot(p2)
    local cross = p1:cross(p2)
    local magnitude = p1:magnitude()
    local normalized = p1:normalize()
    
    print("Dot product of Point 1 and Point 2: " .. string.format("%.2f", dot))
    print("Cross product of Point 1 and Point 2: " .. format_point(cross))
    print("Magnitude of Point 1: " .. string.format("%.2f", magnitude))
    print("Normalized Point 1: " .. format_point(normalized))
    
    -- Create a triangle
    print_subheader("Triangle operations")
    local triangle = cad.Triangle.new(p1, p2, p3)
    
    local area = triangle:area()
    local normal = triangle:normal()
    local perimeter = triangle:perimeter()
    
    print("Triangle area: " .. string.format("%.2f", area))
    print("Triangle perimeter: " .. string.format("%.2f", perimeter))
    print("Triangle normal: " .. format_point(normal))
    
    -- Create a circle
    print_subheader("Circle operations")
    local center = cad.Point3D.new(0, 0, 0)
    local normal = cad.Point3D.new(0, 0, 1)
    local radius = 5.0
    local circle = cad.Circle.new(center, normal, radius)
    
    local circle_area = circle:area()
    local circumference = circle:circumference()
    
    print("Circle area: " .. string.format("%.2f", circle_area))
    print("Circle circumference: " .. string.format("%.2f", circumference))
    
    -- Create points on the circle
    print("\nPoints on the circle:")
    for i = 0, 3 do
        local angle = i * math.pi / 2
        local point = circle:pointAtAngle(angle)
        print(string.format("Point at angle %.2f radians: %s", angle, format_point(point)))
    end
    
    -- Intersection detection
    print_subheader("Intersection detection")
    local line = cad.LineSegment.new(center, p1)
    local intersection = cad.lineTriangleIntersection(line, triangle)
    
    if intersection then
        print("Line intersects triangle at: " .. format_point(intersection))
    else
        print("Line does not intersect triangle")
    end
end

-- Physics Module Examples
function physics_examples()
    print_header("Physics Module Examples")
    
    -- Projectile motion
    print_subheader("Projectile motion")
    local initial_pos = cad.Point3D.new(0, 0, 0)
    local initial_vel = cad.Point3D.new(20, 15, 0)
    local gravity = cad.Point3D.new(0, -9.81, 0)
    
    print("Initial position: " .. format_point(initial_pos))
    print("Initial velocity: " .. format_point(initial_vel))
    print("Gravity: " .. format_point(gravity))
    
    -- Calculate trajectory
    print("\nTrajectory:")
    print(string.format("%10s %15s %15s %15s", "Time (s)", "X Position (m)", "Y Position (m)", "Z Position (m)"))
    print(string.rep("-", 55))
    
    for t = 0, 3, 0.5 do
        local pos = physics.projectileTrajectory(initial_pos, initial_vel, gravity, t)
        print(string.format("%10.2f %15.2f %15.2f %15.2f", t, pos.x, pos.y, pos.z))
    end
    
    -- Calculate key parameters
    local time_of_flight = physics.projectileTimeOfFlight(initial_pos.y, initial_vel.y)
    local range = physics.projectileRange(initial_pos.y, initial_vel)
    local max_height = physics.projectileMaxHeight(initial_pos.y, initial_vel.y)
    
    print("\nProjectile parameters:")
    print("Time of flight: " .. string.format("%.2f", time_of_flight) .. " s")
    print("Range: " .. string.format("%.2f", range) .. " m")
    print("Maximum height: " .. string.format("%.2f", max_height) .. " m")
    
    -- Energy calculations
    print_subheader("Energy calculations")
    local mass = 1.0 -- kg
    local kinetic_energy = physics.kineticEnergy(mass, initial_vel)
    local potential_energy = physics.potentialEnergy(mass, max_height)
    
    print("Mass: " .. mass .. " kg")
    print("Initial kinetic energy: " .. string.format("%.2f", kinetic_energy) .. " J")
    print("Maximum potential energy: " .. string.format("%.2f", potential_energy) .. " J")
    print("Total energy: " .. string.format("%.2f", kinetic_energy + potential_energy) .. " J")
    
    -- Spring force
    print_subheader("Spring force")
    local spring_pos1 = cad.Point3D.new(0, 0, 0)
    local spring_pos2 = cad.Point3D.new(1.5, 0, 0)
    local spring_constant = 10.0 -- N/m
    local rest_length = 1.0 -- m
    
    local spring_force = physics.springForce(spring_constant, rest_length, spring_pos1, spring_pos2)
    
    print("Spring constant: " .. spring_constant .. " N/m")
    print("Rest length: " .. rest_length .. " m")
    print("Current length: " .. string.format("%.2f", spring_pos1:distance(spring_pos2)) .. " m")
    print("Spring force: " .. format_point(spring_force) .. " N")
    
    -- Gravitational force
    print_subheader("Gravitational force")
    local pos1 = cad.Point3D.new(0, 0, 0)
    local pos2 = cad.Point3D.new(10, 0, 0)
    local mass1 = 5.97e24 -- Earth's mass in kg
    local mass2 = 7.35e22 -- Moon's mass in kg
    
    local grav_force = physics.gravitationalForce(mass1, mass2, pos1, pos2)
    
    print("Mass 1: " .. string.format("%.2e", mass1) .. " kg")
    print("Mass 2: " .. string.format("%.2e", mass2) .. " kg")
    print("Distance: " .. string.format("%.2f", pos1:distance(pos2)) .. " m")
    print("Gravitational force: " .. format_point(grav_force) .. " N")
    print("Force magnitude: " .. string.format("%.2e", grav_force:magnitude()) .. " N")
end

-- Electrical Module Examples
function electrical_examples()
    print_header("Electrical Module Examples")
    
    -- Component values
    print_subheader("Component values")
    local resistance = 100.0 -- ohms
    local capacitance = 1.0e-6 -- farads
    local inductance = 0.1 -- henries
    
    print("Resistance: " .. resistance .. " Ω")
    print("Capacitance: " .. capacitance * 1e6 .. " μF")
    print("Inductance: " .. inductance * 1e3 .. " mH")
    
    -- Impedance calculations
    print_subheader("Impedance calculations")
    print(string.format("%15s %20s %20s %20s", "Frequency (Hz)", "Resistor (Ω)", "Capacitor (Ω)", "Inductor (Ω)"))
    print(string.rep("-", 75))
    
    for i = 0, 5 do
        local frequency = 10 ^ i
        
        local z_r = electrical.impedanceResistor(resistance)
        local z_c = electrical.impedanceCapacitor(capacitance, frequency)
        local z_l = electrical.impedanceInductor(inductance, frequency)
        
        print(string.format("%15.2f %20s %20s %20s", 
            frequency, 
            format_complex(z_r), 
            format_complex(z_c), 
            format_complex(z_l)))
    end
    
    -- Series and parallel RLC circuits
    print_subheader("Series and parallel RLC circuits")
    local frequency = 1000.0 -- Hz
    
    local z_series = electrical.impedanceSeriesRLC(resistance, inductance, capacitance, frequency)
    local z_parallel = electrical.impedanceParallelRLC(resistance, inductance, capacitance, frequency)
    
    print("Frequency: " .. frequency .. " Hz")
    print("Series RLC impedance: " .. format_complex(z_series) .. " Ω")
    print("Parallel RLC impedance: " .. format_complex(z_parallel) .. " Ω")
    
    -- Resonant frequency
    print_subheader("Resonant frequency")
    local resonant_freq = electrical.resonantFrequency(inductance, capacitance)
    local q_factor_series = electrical.qualityFactorSeries(resistance, inductance, capacitance)
    local bandwidth = electrical.bandwidth(resonant_freq, q_factor_series)
    
    print("Resonant frequency: " .. string.format("%.2f", resonant_freq) .. " Hz")
    print("Quality factor (series): " .. string.format("%.2f", q_factor_series))
    print("Bandwidth: " .. string.format("%.2f", bandwidth) .. " Hz")
    
    -- Filter calculations
    print_subheader("Filter calculations")
    local cutoff_freq_lp = electrical.cutoffFrequencyLowPassRC(resistance, capacitance)
    local cutoff_freq_hp = electrical.cutoffFrequencyHighPassRC(resistance, capacitance)
    
    print("Low-pass RC cutoff frequency: " .. string.format("%.2f", cutoff_freq_lp) .. " Hz")
    print("High-pass RC cutoff frequency: " .. string.format("%.2f", cutoff_freq_hp) .. " Hz")
    
    -- Power calculations
    print_subheader("Power calculations")
    local voltage = 120.0 -- V
    local current = 2.0 -- A
    local power_factor = 0.8 -- lagging
    
    local power_dc = electrical.powerDC(voltage, current)
    local power_ac = electrical.powerAC(voltage, current, power_factor)
    
    print("Voltage: " .. voltage .. " V")
    print("Current: " .. current .. " A")
    print("Power factor: " .. power_factor)
    print("DC power: " .. power_dc .. " W")
    print("AC power: " .. power_ac .. " W")
    
    -- Voltage divider
    print_subheader("Voltage divider")
    local r1 = 1000.0 -- 1 kOhm
    local r2 = 2000.0 -- 2 kOhm
    local input_voltage = 12.0 -- 12 V
    
    local output_voltage = electrical.voltageDivider(input_voltage, r1, r2)
    
    print("Input voltage: " .. input_voltage .. " V")
    print("R1: " .. r1 .. " Ω")
    print("R2: " .. r2 .. " Ω")
    print("Output voltage: " .. output_voltage .. " V")
end

-- Run all examples
print("RebelCALC Engineering Modules - Lua Script Examples")
print(string.rep("*", 80))

cad_examples()
physics_examples()
electrical_examples()

print(string.rep("*", 80))
print("End of Engineering Modules Examples")
