import numpy as np

def point_plane_distance(A, B, C, D, p):
    '''
    '''
    px, py, pz = p

    dist = np.abs(A*px+B*py+C*pz+D)/np.sqrt(A**2+B**2+C**2)

    return dist

def proj_plane_distance(x, y, z, A, B, C, D, p):
    '''
    Plane equation form: A*x + B*y +C*z + D = 0
    '''
    #X,Y = np.meshgrid(np.linspace(min(x), max(x), 100),
    #                  np.linspace(min(y), max(y), 100))

    #Z = np.zeros(X.shape)

    px, py, pz = p[0], p[1], p[2]
    dist = np.abs(A*px+B*py+C*pz+D)/np.sqrt(A**2+B**2+C**2)

    # normal to the plane
    v = np.asarray([A, B, C])/np.linalg.norm([A, B, C])
    line = []

    # Create a carpet of points, measure point-wise dist
    #for r in range(X.shape[0]):
    #    for c in range(X.shape[1]):
    #        Z[r,c] = (A * X[r,c] + B * Y[r,c] + D)/-C

    # Projection of the centroid onto plane
    n = p + dist*v
    nx, ny, nz = n[0], n[1], n[2]

    # Inner Check
    checkdist = np.abs(A*nx+B*ny+C*nz+D)/np.sqrt(A**2+B**2+C**2)
    if checkdist > 0.1:
        n = p - dist*v

    nx, ny, nz = n[0], n[1], n[2]
    checkdist2 = np.abs(A*nx+B*ny+C*nz+D)/np.sqrt(A**2+B**2+C**2)

    return n, v

def project_atoms_onto_slice(x, y, z, a, b, c, d, out=False):
    '''
    Project atom 3D coordinates onto a 2D slice (however oriented)
    '''
    proj = []
    for j in range(len(x)):
        ref = x[j], y[j], z[j]
        proj.append(proj_plane_distance(x, y, z, a, b, c, d, ref)[0])

    if out:
        with open(out, 'w') as o:
            o.write(str(len(proj))+'\n\n')
            np.savetxt(o, np.c_[np.ones(len(proj)), proj], fmt='%3d %10.6f %10.6f %10.6f')

    return np.c_[proj]

def surf_plane_intersection(A, B, C, D, surf_points):
    '''
    '''

    intersec = []
    for p in surf_points:
        dist = point_plane_distance(A, B, C, D, p)
        if dist < 0.02:
            intersec.append(p)

    return np.array(intersec)
