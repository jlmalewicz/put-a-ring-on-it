#inputs: initial x position [SI units], initial y position [SI units], lens mass [in solar mass units]
def rk4(x_init, y_init, z_init, lens_mass = 1e3, N = 1000, detector_loc = -75): 
    
    ##! Import root finding method to solve for the turning point.
    ##! TODO: implement root finding by hand?
    from scipy.optimize import brentq 
    
    ##! Define necessary physical constants
    G = constants.G
    c = constants.c
    solar_mass = 1.98847e30
    
    ##! Lens Properties
    M = solar_mass*lens_mass # convert lens mass to kg
    R = 2*G*M/(c**2) # Schwarzschild Radius, to be used as as length normalization
    
    detector_location = detector_loc*R
    
    ##! In the 'inbound' portion of the trajectory, r must decrease, therefore
    ##! dr/dphi must be negative (phi is always increasing). However, once the particle passes the turning point
    ##! r1 (closest approach), r must INCREASE, therefore dr/dphi becomes positive.
    ##! To model this, the sign of dr/dphi is initially set to -1. Once the distance
    ##! r is within a certain distance from the turning point distance, the sign is changed
    ##! to +1 to reflect the outbound orbit. This tolerance must be finite, as at r1 itself,
    ##! the equation for dr/dphi diverges. The result is that the orbit "jumps" over the turning
    ##! point. Epsilon is the tolerance, i.e., the difference between r and r1 at which
    ##! point the code knows to 'jump' past the turning point, and reverse the sign of dr/dphi
    epsilon = 0.01*R 
    
    ##! Transform input coordinates (unprimed system) into planar (primed) coordinates
    
    delta = np.arctan2(z_init,y_init)
    
    x0_prime = x_init
    y0_prime = np.linalg.norm([y_init, z_init])
    
    ##! initialize photon locations, in both cartesian and polar coodinates
    x0 = np.array([x0_prime])
    y0 = np.array([y0_prime])
    r0 = np.sqrt(x0**2 + y0**2)
    phi_0 = np.arctan2(y0,x0)
    
    ##! Set the impact parameter of the particle.
    ##! RIGHT NOW, IMPACT PARAMETER IS THE INITIAL Y COORDINATE.
    ##! THIS ASSUMES THE PARTICLE ENTERS THE REGION WITH A
    ##! VELOCITY PARALLEL TO THE X AXIS.
    ##! TODO: Generalize this.
    b = y0
    
    ##! Check whether the photon orbit is bounded or unbounded.
    ##! If the impact parameter b is below this threshold,
    ##! the photon will be captured by the mass. In this case
    ##! the function returns the string 'Captured'
    if b**2 < 27*(R/2)**2:
        return -1, -1
    
    ##! Function to determine the turning point r1.
    ##! Finding the root of this function returns the 
    ##! turning point
    def turning_point(r):
        return r**3 - r*b**2+R*b**2
    
    ##! Calculate the turning point via the Brent method.
    ##! First input: function of which the root is desired
    ##! Second input: lower bound of the root.
    ##! Third input: upper bound of the root. 
    ##! THERE SHOULD BE TWO POSITIVE SOLUTIONS FOR THE TURNING POINT,
    ##! ONE VERY NEAR THE SCHWARZCHILD RADIUS, AND ONE FURTHER AWAY.
    ##! THIS REFLECTS THE TWO POINTS AT WHICH 1/b^2 INTERSECTS Weff.
    ##! WHAT 
    r1 = brentq(turning_point, 2*R, b)
    
    ##! Define the step size in phi (azimuthal angle), as the difference
    ##! between the initial angle phi_0 and 2*pi,
    ##! divided by the number of steps N
    dphi = (2*np.pi-phi_0)/N
    
    ##! Set the sign of the derivative to negative initially, reflecting
    ##! the inbound portion of the trajectory
    sign = -1
    
    ##! Create an array of r values that range from phi_0 to 2*pi with step dphi
    phi = np.arange(phi_0,2*np.pi,dphi)

    
    ##! Initialize an array to hold the values of r calculated for each respective phi
    ##! via the ODE solver. Make the first value the initial value r0
    r = np.array([r0])
    
    ##! Define the ODE for r(phi), in terms of the effective potential weff.
    ##! Input: radial distance, in SI units
    def f(r):
        weff = (1-R/r)/(r**2)
        return 1/(sign*(1/(b**2) - weff)**(-1/2)*(1/(r**2)))
    
    ##! RK4 SOLVER
    ##! Calculate the respective distance r for the given value of phi.
    ##! If the particle exceeds the initial distance r0, terminate
    ##! the calculation. Later, this condition will be replaced with
    ##! the detector location. 
    ##! If the particle is within epsilon distance of the turning
    ##! point, jump over the turning point and change the sign of 
    ##! dr/dphi
    for i in range(1,len(phi)):
        k1 = dphi*f(r[i-1])
        k2 = dphi*f(r[i-1] + 0.5*k1)
        k3 = dphi*f(r[i-1] + 0.5*k2)
        k4 = dphi*f(r[i-1] + k3)
        new_r = r[i-1]+(1/6)*(k1+2*k2+2*k3+k4)
        r = np.append(r, new_r)
        if np.linalg.norm(r[i])*np.cos(phi[i]) < detector_location:
            phi = phi[0:len(r)]
            break
        elif np.linalg.norm(r[i]) - r1 < epsilon:
            sign = 1
        
    ##! Convert the output from polar to cartesian coordinates.
    ##! Y axis: Axis connecting the lens location to the turning point.
    ##! X axis: Axis perpendicular to Y, and contained in the plane
    ##! defined by the lens location and any two points along the 
    ##! photon trajectory (Conservation of angular momentum requires
    ##! that the photon trajectory be confined to a plane containing 
    ##! the lens.)
    
    x = r*np.cos(phi)
    y = r*np.sin(phi)
    
    ##! We now need to backtrace, and calculate the y value exactly when the 
    ##! particle crossed the detector
    
    y2, y1 = y[-1], y[-2]
    x2, x1 = x[-1], x[-2]
    
    m = (y2-y1)/(x2-x1)
    
    y_detector = m*detector_location + y2-m*x2
    
    x[-1] = detector_location
    y[-1] = y_detector
    
    ##! now we must convert the final coordinate, y_detector
    ##! which is in the primed system, back into the unprimed system
    
    
    ##! Particles that have an initial y position y_init that is negative are simply
    ##! flipped, treated as if the y_init was positive, and then flipped back. This 
    ##! avoids having to redefine the angle phi, which would change the treatment of the
    ##! sign of the derivative. 
    #if y_init < 0:
    #    y = -1*y
    
    ##! now we must convert the final coordinate, y_detector
    ##! which is in the primed system, back into the unprimed system
    y_out = y[-1]*np.cos(delta)
    z_out = y[-1]*np.sin(delta)
    
    ##! Return the array of x and y coordinates. For efficiency, this will eventually
    ##! be replaced with returning only the final coordinates x[-1] and y[-1].
    ##! STILL NEED TO ADD THE THRID DIMENSION Z TO THIS. Plan to do this by rotating
    ##! each plane into the xy plane, modeling, and then rotating back, and calculating the 
    ##! final (x,y,z) position from that.
    return y_out, z_out