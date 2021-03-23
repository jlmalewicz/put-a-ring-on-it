#inputs: initial x position [SI units], initial y position [SI units], lens mass [in solar mass units]
def photon_orbit5(x_init, y_init, lens_mass = 1e3, N = 1000): 
    
    ##! Import root finding method to solve for the turning point.
    ##! TODO: implement root finding by hand?
    from scipy.optimize import brentq 
    
    ##! Set font size for plots
    plt.rcParams.update({'font.size': 22})
    
    ##! Define necessary physical constants
    G = constants.G
    c = constants.c
    solar_mass = 1.98847e30
    
    ##! Lens Properties
    M = solar_mass*lens_mass # convert lens mass to kg
    R = 2*G*M/(c**2) # Schwarzschild Radius, to be used as as length normalization
    
    ##! initialize photon locations, in both cartesian and polar coodinates
    x0 = np.array([x_init])
    y0 = np.array([y_init])
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
        print('Captured')
        return 'Captured'
    
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
    
    ##! Define the step size in r (radial distance), as the distance 
    ##! from the turning point (r1) to the initial location (r0),
    ##! divided by the number of steps N
    dr = (r0-r1)/N
    
    ##! Create an array of r values that range from r1 to r0 with step dr
    r = np.arange(r1,r0,dr)
    
    ##! Initialize an array to hold the values of phi calculated for each respective r
    ##! via the ODE solver
    phi = np.zeros(len(r))
    
    ##! Initialize the first location as the turning point, with distance r1,
    ##! and angle (by definition) of pi/2.
    ##! RIGHT NOW THE SOLVER TRACKS THE PARTICLE OUTWARD, FROM THE 
    ##! TURNING POINT TO THE INITIAL POINT.
    ##! TODO: REVERSE THIS, SO THAT IT TRACKS THE PARTICLE INWARD.
    ##! THIS WILL ALLOW EASIER TRANSLATION OF THE COORDINATE SYSTEM,
    ##! AS UNDER THE CURRENT FRAMEWORK THE Y AXIS IS DEFINED AS THE 
    ##! AXIS THAT INTERSECTS THE TURNING POINT (phi = pi/2)
    r[0] = r1
    phi[0] = np.pi/2
    
    ##! Define the ODE for phi(r), in terms of the effective potential weff.
    ##! Input: radial distance, in SI units
    def f(r):
        weff = (1-R/r)/(r**2)
        return (1/(b**2) - weff)**(-1/2)*(1/(r**2))
    
    ##! RK4 SOLVER
    ##! Calculate the respective angle phi for the given value of r.
    ##! Again, this only calculates from the turning point r1 to the 
    ##! initial point r0
    for i in range(1,len(phi)):
        k1 = dr*f(r[i-1]+ dr)
        k2 = dr*f(r[i-1]+ 0.5*dr)
        k3 = dr*f(r[i-1]+ 0.5*dr)
        k4 = dr*f(r[i-1]+ dr)
        phi[i] = phi[i-1]+(1/6)*(k1+2*k2+2*k3+k4)
        
    ##! Convert the output from polar to cartesian coordinates.
    ##! Y axis: Axis connecting the lens location to the turning point.
    ##! X axis: Axis perpendicular to Y, and contained in the plane
    ##! defined by the lens location and any two points along the 
    ##! photon trajectory (Conservation of angular momentum requires
    ##! that the photon trajectory be confined to a plane containing 
    ##! the lens.)
    x = r*np.cos(phi)
    y = r*np.sin(phi)
    
    ##! Initialize axes for trajectory plot.
    fig, ax = plt.subplots(1,1, figsize = (10,8))
    
    ##! Plot trajectory from origin to turning point. Map radial
    ##! distance into a colormap defining the color of each point.
    fig1 = ax.scatter(x/R,y/R, c = r/R, cmap = 'jet', s = 10)
    ##! Plot the trajectory "mirrored" around the turning point.
    ##! Both origin->turning point and turning point-> exit trajectories
    ##! should be symmetric around the turning point.
    fig2 = ax.scatter(-x/R, y/R, c = r/R, cmap = 'jet', s = 10)
    
    eh_theta = np.linspace(0,2*np.pi,1000)
    eh_x = R*np.cos(eh_theta)
    eh_y = R*np.sin(eh_theta)
    
    capture_theta = np.linspace(0,2*np.pi,1000)
    capture_x = (3/2)*np.sqrt(3)*R*np.cos(eh_theta)
    capture_y = (3/2)*np.sqrt(3)*R*np.sin(eh_theta)
    
    ##! Draw the event horizon (circle radius R, black)
    ax.plot(eh_x/R,eh_y/R, color = 'k')
    
    ##! Draw the "capture radius" (circle radius sqrt(27/4)*R, orange)
    ax.plot(capture_x/R,capture_y/R, color = 'tab:orange')
    
    ##! Draw the axes
    ax.axvline(0, linestyle = 'dashed')
    ax.axhline(0, linestyle = 'dashed')
    
    ##! Set the axes limits
    ##! TODO: Generalize this
    ax.set_xlim(-x0/R, x0/R)
    ax.set_ylim(-x0/R, x0/R)
    
    ##! Set colorbar and respective labels
    c1 = plt.colorbar(fig1, ax = ax)
    c1.set_label('Distance r [R$_S$]')
    ax.set_xlabel('X ($R_S$)', fontsize = 20)
    ax.set_ylabel('Y ($R_S$)', fontsize = 20)
    ax.set_title('Impact Parameter = %s R' % (b/R))
    
    plt.tight_layout()
    
    return r, phi
