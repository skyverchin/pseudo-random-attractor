import math
import random as rnd
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from mpl_toolkits.mplot3d import Axes3D
#from stl import mesh
#from scipy.spatial import Delaunay
import surf2stl

MAXITERATIONS = 100000
NEXAMPLES = 1000

def CreatingAttractor():
    for n in range(NEXAMPLES):        
        lyapunov = 0
        xmin= 1e32
        xmax=-1e32
        ymin= 1e32
        ymax=-1e32
        zmin= 1e32
        zmax=-1e32
        fx, fy, fz, x, y, z = [], [], [], [], [], []
        
        # Initialize coefficients for this attractor
        for i in range(6):
            fx.append(rnd.uniform(-2, 2))
            fy.append(rnd.uniform(-2, 2))
            fz.append(rnd.uniform(-2, 2))
    
        # Calculate the attractor
        classify = True
        x.append(rnd.uniform(-0.5, 0.5))
        y.append(rnd.uniform(-0.5, 0.5))
        z.append(rnd.uniform(-0.5, 0.5))
        
        d0 = -1
        while d0 <= 0:
            xe = x[0] + rnd.uniform(-0.5, 0.5) / 1000.0
            ye = y[0] + rnd.uniform(-0.5, 0.5) / 1000.0
            ze = z[0] + rnd.uniform(-0.5, 0.5) / 1000.0
            dx = x[0] - xe
            dy = y[0] - ye
            dz = z[0] - ze
            d0 = math.sqrt(dx * dx + dy * dy + dz * dz)

        for i in range(MAXITERATIONS):
            # Calculate next term  
            x.append(fx[0] + fx[1]*x[i-1] + fx[2]*x[i-1]*x[i-1] + fx[3]*x[i-1]*y[i-1] + fx[4]*y[i-1] + fx[5]*y[i-1]*y[i-1])
            y.append(fy[0] + fy[1]*x[i-1] + fy[2]*x[i-1]*x[i-1] + fy[3]*x[i-1]*y[i-1] + fy[4]*y[i-1] + fy[5]*y[i-1]*y[i-1])
            z.append(fz[0] + fz[1]*x[i-1] + fz[2]*x[i-1]*x[i-1] + fz[3]*x[i-1]*y[i-1] + fz[4]*y[i-1] + fz[5]*y[i-1]*y[i-1])

            xenew = fx[0] + fx[1]*xe + fx[2]*xe*xe + fx[3]*xe*ye + fx[4]*ye + fx[5]*ye*ye
            yenew = fy[0] + fy[1]*xe + fy[2]*xe*xe + fy[3]*xe*ye + fy[4]*ye + fy[5]*ye*ye
            zenew = fz[0] + fz[1]*xe + fz[2]*xe*xe + fz[3]*xe*ye + fz[4]*ye + fz[5]*ye*ye

            # Update the bounds 
            xmin = min(xmin,x[i])
            xmax = max(xmax,x[i])

            ymin = min(ymin,y[i])
            ymax = max(ymax,y[i])
            
            zmin = min(zmin,z[i])
            zmax = max(zmax,z[i])

            # Does the series tend to infinity
            if xmin < -1e10 or xmax > 1e10 or ymin < -1e10 or ymax > 1e10 or zmin < -1e10 or zmax > 1e10:
                classify = False
                print ("Infinite attractor")
                break

            # Does the series tend to a point
            dx = x[i] - x[i-1]
            dy = y[i] - y[i-1]
            dz = z[i] - z[i-1]
            if abs(dx) < 1e-10 and abs(dy) < 1e-10 and abs(dz) < 1e-10:
                classify = False
                print ("Point attractor")
                break
            
            # Calculate the lyapunov exponents
            if i > 1000:
                dx = x[i] - xenew
                dy = y[i] - yenew
                dz = z[i] - zenew
                dd = math.sqrt(dx * dx + dy * dy + dz * dz)
                lyapunov += math.log(math.fabs(dd / d0))
                xe = x[i] + d0 * dx / dd
                ye = y[i] + d0 * dy / dd
                ze = z[i] + d0 * dy / dd
            
            # Classify the series according to lyapunov
            if classify:
                if abs(lyapunov) < 10:
                    print ("Neutrally stable")
                elif lyapunov < 0:
                    print ("Periodic {} ".format(lyapunov))
                else:
                    print ("Chaotic {} ".format(lyapunov))
                    PlotAtractor(x,y,z)

def PlotAtractor(x,y,z):
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')    

    for i in range(MAXITERATIONS):
        if i > 100:
            #ax.scatter(x,y,z)
            plt.scatter(x,y,z)

            #ax.plot(x,y,z)
            #plt.plot(x, y, z)
            X, Y = np.meshgrid(x, y)
            Z = np.meshgrid(z)                       
            plt.show()
            surf2stl.write('result.stl', X, Y, Z)
            
            
    #Something go wrong here ¯\_(ツ)_/¯, while use this code we use surf2stl(modified)
    #for i in range(MAXITERATIONS):
    #    if i > 100:
    #        points = np.column_stack((np.array(x), np.array(y), np.array(z)))
    #        triangles = Delaunay(points)

    #        normal = np.cross(points[triangles.simplices[:, 1]] - points[triangles.simplices[:, 0]], points[triangles.simplices[:, 2]] - points[triangles.simplices[:, 0]])
    #        normal /= np.linalg.norm(normal, axis=1).reshape(-1, 1)

    #        surface = mesh.Mesh(np.zeros(len(triangles.simplices), dtype=mesh.Mesh.dtype))
    #        for j, triangle in enumerate(triangles.simplices):
    #            for k, point in enumerate(triangle):
    #                surface.vectors[j][k] = points[point]
    #            surface.normals[j] = normal[j]

    #        surface.save('example.stl')

CreatingAttractor()
