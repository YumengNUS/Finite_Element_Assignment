# =============================================================================
# Program to solve 2D Helmholtz equation using Galerkin FE
# 
# by   :XU Yumeng
# date :Nov/2018
# =============================================================================

import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from dxiIdxJ_Yumeng import dxiIdxJ
from psi_Yumeng import psi

# =============================================================================
# Set up the mesh
#
# numNodes = total number of global nodes
# numElems = total number of global elements
#
# nodePos(n,0) = x position of node n
# nodePos(n,1) = y position of node n
#
# elemNode(e,0) = first global node in element e (local node 1) 
# elemNode(e,1) = second global node in element e (local node 2) 
# elemNode(e,2) = third global node in element e (local node 3) 
# elemNode(e,3) = fourth global node in element e (local node 4) 
#
# BC(n,0) = type of boundary condition at node n
#           BC(n,0) = 1 => essential (potential)
#           BC(n,0) = 2 => natural (flux)
# BC(n,1) = value of the boundary condition at node n
# =============================================================================

# Change resolution using this variable
elemsPerSide = 32


# Constants in the governing equation
mu = 1000.0    # Shear modulus, Pa
omega =1000.0  # Frequency, Hz
rho = 1060.0   # Density, kg/m3
kappa = rho*omega*omega/mu  # 1/m^2
#kappa = 0
# Size of the domain (square, remains constant)
sideLength = 0.08  # m
# Number of nodes in each element (bilinear, remains constant)
nodesPerElement = 4

# =============================================================================
# Generate nodes
# =============================================================================
nodesPerSide = elemsPerSide + 1
numNodes = nodesPerSide*nodesPerSide
nodePos = np.zeros([numNodes, 2])
incr = sideLength/elemsPerSide
n = 0
for j in range (0, nodesPerSide):
    for i in range (0, nodesPerSide):
        nodePos[n, 0] = i*incr
        nodePos[n, 1] = j*incr
        n += 1

# =============================================================================
# Generate elements
# =============================================================================
numElems = elemsPerSide*elemsPerSide
elemNode = np.zeros([numElems, nodesPerElement], dtype=int)
e = 0
for j in range (0, elemsPerSide):
    for i in range (0, elemsPerSide):
        elemNode[e, 0] = j*nodesPerSide + i
        elemNode[e, 1] = elemNode[e, 0] + 1
        elemNode[e, 2] = elemNode[e, 0] + nodesPerSide
        elemNode[e, 3] = elemNode[e, 0] + nodesPerSide + 1
        e += 1

# =============================================================================
# Generate boundary conditions
# =============================================================================
BCs = np.zeros([numNodes, 2])
n = 0
for j in range (0, nodesPerSide):
    for i in range (0, nodesPerSide):
        if (i == 0):
            BCs[n, 0] = 1
            BCs[n, 1] = 1.0e-6
        else:
            BCs[n, 0] = 2
            BCs[n, 1] = 0
        n += 1

# =============================================================================
# Set up the numerical integration points 'gaussPos' and weights 'gaussWei' 
# =============================================================================
numGaussPoints = 4  # 2x2 Gauss points

gp1 = 0.5 - (1.0/(2.0*math.sqrt(3.0))) 
gp2 = 0.5 + (1.0/(2.0*math.sqrt(3.0)))

gaussPos = np.array([[gp1, gp1], [gp2, gp1], [gp1, gp2], [gp2, gp2]])
gaussWei = np.array([0.25, 0.25, 0.25, 0.25]) 

# =============================================================================
# Assemble the global stiffness matrix, K
# =============================================================================
K = np.zeros([numNodes, numNodes])
f = np.zeros([numNodes])
u = np.zeros([numNodes])

# Loop over elements
for elem in range (0, numElems):
    
    # Initialise element stiffness matrix (EK) 
    EK = np.zeros([nodesPerElement, nodesPerElement])
    
    # =========================================================================
    # Create the element stiffness matrix
    # =========================================================================
    # Q. Write your code to assemble EK here
        
    # loop over local nodes of each element       
    for m in range (0,4):
        for n in range (0,4):
            for iGP in range (0, numGaussPoints):
                # Get the values of gauss point
                xi1, xi2 = gaussPos[iGP,:]
                dpsimIdxi1 = psi (m, 1, xi1, xi2) #dPsim/dxi1
                dpsimIdxi2 = psi (m, 2, xi1, xi2) #dPsim/dxi2
                dpsinIdxi1 = psi (n, 1, xi1, xi2) #dPsin/dxi1
                dpsinIdxi2 = psi (n, 2, xi1, xi2) #dPsin/dxi2
                #Jacobian and dxi/dx
                dxiIdx, J = dxiIdxJ(elem, nodesPerElement, xi1, xi2, nodePos, elemNode)
                Psim = psi (m, 0, xi1, xi2)
                # Calculate the intergration using 4 gauss points
                EK[m,n] = EK[m,n] + gaussWei[iGP]*(dpsimIdxi1*dxiIdx[0,0]*dpsinIdxi1*dxiIdx[0,0]+dpsimIdxi2*dxiIdx[1,1]*dpsinIdxi2*dxiIdx[1,1]-Psim*kappa)*J
                
  
    # =========================================================================
    # Assemble EK into the global stiffness matrix
    # =========================================================================
    
    # Q. Write your code to add EK to K here
            # Transform local node number to global node number
            K[elemNode[elem,m],elemNode[elem,n]] = K[elemNode[elem,m],elemNode[elem,n]] + EK[m,n]
# =============================================================================
# Apply boundary & initial conditions
# =============================================================================

# Q. Write your code to apply the boundary conditions here

# Boundary condition
u = BCs[:,1]
f = BCs[:,1]
for n in range(0,numNodes):  
   if BCs[n,0] == 1:
       # Overwrite row with zeros and then put 1 at the position with essential BC
       K[n,:] = 0
       K[n,n] = 1
       # Boundary condition for RHS vector f
       f[n] = BCs[n,1]

    
# =============================================================================
# Solve the linear system Ku=f for u
# =============================================================================
u = np.linalg.solve(K, f)

# =============================================================================
# Surface plot of the solution
# =============================================================================
x = nodePos[0:nodesPerSide, 0]
x, y = np.meshgrid(x, x)
u = u.reshape(nodesPerSide, nodesPerSide)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(x, y, u, cmap=plt.cm.coolwarm, linewidth=0, alpha=0.5)
