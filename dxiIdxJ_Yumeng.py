# =============================================================================
# Function to calculate dxi_i/dx_j and the element Jacobian, J
#
# elem = global element number
# nodesPerElement = number of nodes per element 
# xi1 = xi1 coordinate (range 0-1)
# xi2 = xi2 coordinate (range 0-1)
# nodePos = matrix of node coordinates
# elemNode = list of nodes for each element
# =============================================================================

import numpy as np
from psi_Yumeng import psi

def dxiIdxJ (elem, nodesPerElement, xi1, xi2, nodePos, elemNode):
    
    # Q: Implement your code to calculate dxi_i/dx_j and the Jacobian here
    #Create a 2x2 matrix, dxi/dx
    dxi_iIdx_j = np.zeros([2,2])
    nodenum1, nodenum2, nodenum3, nodenum4 = elemNode[elem, :]
    x1, y1 = nodePos[nodenum1, :] 
    x2, y2 = nodePos[nodenum2, :]
    x3, y3 = nodePos[nodenum3, :]
    dxIdxi = np.zeros([2,2])
    
    
    
    
    
    dxIdxi[0,0] = (x2-x1)/(1-0)
    dxIdxi[0,1] = (x3-x1)/(1-0)
    dxIdxi[1,0]= (y2-y1)/(1-0)
    dxIdxi[1,1] = (y3-y1)/(1-0)
    det = np.linalg.det(dxIdxi)
    
    
    dxi_iIdx_j[0,0] = (1/det)*(dxIdxi[1,1])
    dxi_iIdx_j[0,1] = -(1/det)*(dxIdxi[1,0])
    dxi_iIdx_j[1,0] = -(1/det)*(dxIdxi[0,1])
    dxi_iIdx_j[1,1] = (1/det)*(dxIdxi[0,0])
    J = abs(det)
        # absolute value, the Jacobian must be positive
    return dxi_iIdx_j, J
    # returen dx_i/dx_j, Jacobian

