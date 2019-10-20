# =============================================================================
# Function to evaluate the bilinear Lagrange basis functions and derivatives
#
# num = local node number
# der = derivative: 0 = function; 1 = d/dxi_1; 2 = d/dxi_2
# xi1 = xi1 coordinate (range 0-1)
# xi2 = xi2 coordinate (range 0-1)
# =============================================================================

def psi (num, der, xi1, xi2):
    
    if (der == 0):
        if (num == 0):
            return (1.0-xi1)*(1.0-xi2)
        elif (num == 1):
            return xi1*(1.0-xi2)
        elif (num == 2):
            return (1.0-xi1)*xi2            
        elif (num == 3):
            return xi1*xi2            
        else:
            return 0
    elif (der == 1):
        if (num == 0):
            return -(1.0-xi2)
        elif (num == 1):
            return 1.0-xi2
        elif (num == 2):
            return -xi2            
        elif (num == 3):
            return xi2            
        else:
            return 0      
    elif (der == 2):
        if (num == 0):
            return -(1.0-xi1)
        elif (num == 1):
            return -xi1
        elif (num == 2):
            return 1.0-xi1            
        elif (num == 3):
            return xi1            
        else:
            return 0       
    else:
        return 0
