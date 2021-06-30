# -*- coding: utf-8 -*-
"""
Created on Fri Mar 19 16:02:31 2021

@author: nicol
"""

# This function is used to interpolate the refrigerant properties obtained using thermeprop on Matlab.
# The input 'mat' is the property matrix that is uploaded in the code at beginning of the 'Main module.
# In order to compute a property the user must specify which property has to be computed ('prop') and then specify
# the pair of known properties on which the computation is based.
def propfunction(mat, prop, T = False, P = False, rho = False, s = False, h = False, Q = False  ):
    """
    Parameters
    ----------
    prop : TYPE, mandatory
        properties wanted: can be T, P, h, s, rho, vmisc, k, T_sat, Q.
    T : TYPE, optional
        Temperature in K. The default is False.
    P : TYPE, optional
        Pressure in Pa. The default is False.
    rho : TYPE, optional
        density in kg/m^3. The default is False.
    s : TYPE, optional
        entropy in J/(kgK). The default is False.
    h : TYPE, optional
        enthalpy in J/kg. The default is False.
    Q : TYPE, optional
        enthalpy in J/kg. The default is False.
    Returns
    -------
    None.
    
    """
    # 'prop_value' is the output of the function. It is initialized here as 'None'
    prop_value = None
    
    # Procedure to compute the property value given temperature and vapor quality.
    if T != False and Q != False: 
        
        # The property matrix is sliced accordingly with the required property and the known properties.
        mat = mat['properties']['fluid_properties_TQ']
        T_vect = mat[0,0][0,0][1][0]
        Q_vect = mat[0,0][0,0][0][0]
        
        idx1 = 0
        idx2 = 0
        
        # Interpolation nodes are identified.
        while T_vect[idx1] < T:
            idx1 += 1
        lever1 = (T - T_vect[idx1 - 1])/(T_vect[idx1] - T_vect[idx1 - 1])   
        
        while Q_vect[idx2] < Q:
            idx2 += 1
        lever2 = (Q - Q_vect[idx2 - 1])/(Q_vect[idx2] - Q_vect[idx2 - 1]) 
        
        # Each property is characterised by an index 'N'
        if prop == 'rho':
            N = 2     
        elif prop == 'h':
            N = 3
        elif prop == 's':
            N = 4
        elif prop == 'P':
            N = 5
        else:
            raise Exception('Invalid inputs')
        
        # Interpolation is performed.
        prop_int = mat[0,0][0,0][N][idx1-1][idx2-1:(idx2+1)] + lever1*(mat[0,0][0,0][N][idx1][idx2-1:(idx2+1)] - mat[0,0][0,0][N][idx1-1][idx2-1:(idx2+1)])
        prop_value = prop_int[0] + lever2*(prop_int[1] - prop_int[0])    
    
    # Procedure to compute the property value given temperature and pressure.
    elif P != False and T != False: 
        
        # The property matrix is sliced accordingly with the required property and the known properties. 
        mat = mat['properties']['fluid_properties_PT']
        P_vect = mat[0,0][0,0][0][0]
        T_vect = mat[0,0][0,0][1][0]
        
        idx1 = 0
        idx2 = 0
        
        # Interpolation nodes are identified.
        while P_vect[idx1] < P:
            idx1 += 1
        lever1 = (P - P_vect[idx1 - 1])/(P_vect[idx1] - P_vect[idx1 - 1])   
        
        while T_vect[idx2] < T:
            idx2 += 1
        lever2 = (T - T_vect[idx2 - 1])/(T_vect[idx2] - T_vect[idx2 - 1])        
        
        # Each property is characterised by an index 'N'
        if prop == 'H':
            N = 2    
        elif prop == 'rho':
            N = 3
        else:
            raise Exception('Invalid inputs')
        
        # Interpolation is performed.
        prop_int = mat[0,0][0,0][N][idx1-1][idx2-1:(idx2+1)] + lever1*(mat[0,0][0,0][N][idx1][idx2-1:(idx2+1)] - mat[0,0][0,0][N][idx1-1][idx2-1:(idx2+1)])
        prop_value = prop_int[0] + lever2*(prop_int[1] - prop_int[0])
         
    # Procedure to compute the property value given temperature and density.
    elif T != False and rho != False: 
        
        # The property matrix is sliced accordingly with the required property and the known properties. 
        mat = mat['properties']['fluid_properties_Trho']
        T_vect = mat[0,0][0,0][0][0]
        rho_vect = mat[0,0][0,0][1][0]
        
        idx1 = 0
        idx2 = 0
        
        # Interpolation nodes are identified.
        while T_vect[idx1] < T:
            idx1 += 1
        lever1 = (T - T_vect[idx1 - 1])/(T_vect[idx1] - T_vect[idx1 - 1])   
        while rho_vect[idx2] < rho:
            idx2 += 1
        lever2 = (rho - rho_vect[idx2 - 1])/(rho_vect[idx2] - rho_vect[idx2 - 1])        
        
        # Each property is characterised by an index 'N'
        if prop == 'P':
            N = 2    
        elif prop == 's':
            N = 3
        elif prop == 'Q':
            N = 4
        elif prop == 'h': 
            N = 5
        elif prop == 'u':
            N = 6
        elif prop == 'cv':
            N = 7
        elif prop == 'dPdT':
            N = 8
        elif prop == 'dudT':
            N = 9
        elif prop == 'k_star':
            N = 11
        else:
            raise Exception('Invalid inputs')
        
        # Interpolation is performed.
        prop_int = mat[0,0][0,0][N][idx1-1][idx2-1:(idx2+1)] + lever1*(mat[0,0][0,0][N][idx1][idx2-1:(idx2+1)] - mat[0,0][0,0][N][idx1-1][idx2-1:(idx2+1)])
        prop_value = prop_int[0] + lever2*(prop_int[1] - prop_int[0])
    
    # Procedure to compute the property value given pressure and density.
    elif P != False and rho != False: 
        
        # The property matrix is sliced accordingly with the required property and the known properties. 
        mat = mat['properties']['fluid_properties_Prho']
        P_vect = mat[0,0][0,0][0][0]
        rho_vect = mat[0,0][0,0][1][0]
        
        idx1 = 0
        idx2 = 0
        
        # Interpolation nodes are identified.
        while P_vect[idx1] < P:
            idx1 += 1
        lever1 = (P - P_vect[idx1 - 1])/(P_vect[idx1] - P_vect[idx1 - 1])   
    
        while rho_vect[idx2] < rho:
            idx2 += 1
        lever2 = (rho - rho_vect[idx2 - 1])/(rho_vect[idx2] - rho_vect[idx2 - 1])        
        
        # Each property is characterised by an index 'N'
        if prop == 'T':
            N = 2    
        else:
            raise Exception('Invalid inputs')
        
        # Interpolation is performed.
        prop_int = mat[0,0][0,0][N][idx1-1][idx2-1:(idx2+1)] + lever1*(mat[0,0][0,0][N][idx1][idx2-1:(idx2+1)] - mat[0,0][0,0][N][idx1-1][idx2-1:(idx2+1)])
        prop_value = prop_int[0] + lever2*(prop_int[1] - prop_int[0])
            
    # Procedure to compute the property value given pressure and enthalpy.
    elif h != False and P != False:
        
        # The property matrix is sliced accordingly with the required property and the known properties.
        mat = mat['properties']['fluid_properties_PH']
        P_vect = mat[0,0][0,0][0][0]
        h_vect = mat[0,0][0,0][1][0]
        
        idx1 = 0
        idx2 = 0
        
        # Interpolation nodes are identified.
        while P_vect[idx1] < P:
            idx1 += 1
        lever1 = (P - P_vect[idx1 - 1])/(P_vect[idx1] - P_vect[idx1 - 1])   
        
        while h_vect[idx2] < h:
            idx2 += 1
        lever2 = (h - h_vect[idx2 - 1])/(h_vect[idx2] - h_vect[idx2 - 1])        
        
        # Each property is characterised by an index 'N'
        if prop == 'rho': 
            N = 2
        elif prop == 'T':
            N = 3
        elif prop == 's':
            N = 4
        elif prop == 'Q':
            N = 5
        elif prop == 'mvisc':
            N = 6
        elif prop == 'k':
            N = 7
        else:
            raise Exception('Invalid inputs')
        
        # Interpolation is performed.
        prop_int = mat[0,0][0,0][N][idx1-1][idx2-1:(idx2+1)] + lever1*(mat[0,0][0,0][N][idx1][idx2-1:(idx2+1)] - mat[0,0][0,0][N][idx1-1][idx2-1:(idx2+1)])
        prop_value = prop_int[0] + lever2*(prop_int[1] - prop_int[0])    
    
    # Procedure to compute the property value given entropy and density.    
    elif s != False and rho != False:
        
        # The property matrix is sliced accordingly with the required property and the known properties.
        mat = mat['properties']['fluid_properties_Srho']
        s_vect = mat[0,0][0,0][0][0]
        rho_vect = mat[0,0][0,0][1][0]
        
        idx1 = 0
        idx2 = 0
        
        # Interpolation nodes are identified.
        while s_vect[idx1] < s:
            idx1 += 1
        lever1 = (s - s_vect[idx1 - 1])/(s_vect[idx1] - s_vect[idx1 - 1])   
        
        while rho_vect[idx2] < rho:
            idx2 += 1
        lever2 = (rho - rho_vect[idx2 - 1])/(rho_vect[idx2] - rho_vect[idx2 - 1])        
    
        # Each property is characterised by an index 'N'
        if prop == 'P': 
            N = 2
        elif prop == 'T':
            N = 3
        else:
            raise Exception('Invalid inputs')
        
        # Interpolation is performed.
        prop_int = mat[0,0][0,0][N][idx1-1][idx2-1:(idx2+1)] + lever1*(mat[0,0][0,0][N][idx1][idx2-1:(idx2+1)] - mat[0,0][0,0][N][idx1-1][idx2-1:(idx2+1)])
        prop_value = prop_int[0] + lever2*(prop_int[1] - prop_int[0])    
    
    # Procedure to compute the property value given entropy and density.   
    elif s != False and P != False:
        
        # The property matrix is sliced accordingly with the required property and the known properties.
        mat = mat['properties']['fluid_properties_PS']
        P_vect = mat[0,0][0,0][0][0]
        s_vect = mat[0,0][0,0][1][0]
        
        idx1 = 0
        idx2 = 0
        
        # Interpolation nodes are identified.
        while P_vect[idx1] < P:
            idx1 += 1
        lever1 = (P - P_vect[idx1 - 1])/(P_vect[idx1] - P_vect[idx1 - 1])   
    
        while s_vect[idx2] < s:
            idx2 += 1
        lever2 = (s - s_vect[idx2 - 1])/(s_vect[idx2] - s_vect[idx2 - 1])        
    
        # Each property is characterised by an index 'N'
        if prop == 'h': 
            N = 2
        elif prop == 'rho': 
            N = 3
        elif prop == 'Q': 
            N = 4
        elif prop == 'T': 
            N = 5     
        else:
            raise Exception('Invalid inputs')
        
        # Interpolation is performed.
        prop_int = mat[0,0][0,0][N][idx1-1][idx2-1:(idx2+1)] + lever1*(mat[0,0][0,0][N][idx1][idx2-1:(idx2+1)] - mat[0,0][0,0][N][idx1-1][idx2-1:(idx2+1)])
        prop_value = prop_int[0] + lever2*(prop_int[1] - prop_int[0])        
        
    # The following properties can be obtained using providing only pressure as a known quantity:
        # - Saturation temperature
        # - Saturated liquid enthalpy
        # - Saturated vapor enthalpy
        # - Enthalpy of vaporizaion
    elif P != False: 
    
        # The property matrix is sliced accordingly with the required property and the known properties.
        mat = mat['properties']['fluid_properties']
        P_vect = mat[0,0][0,0][0][0]*1e5
        
        # Interpolation node is identified.
        index = 0
        while P_vect[index] < P:
            index += 1
        lever = (P - P_vect[index-1])/(P_vect[index] - P_vect[index-1])    
        
        # Interpolation is performed
        if prop == 'T_sat': 
            prop_value = mat[0,0][0,0][1][0][index-1] + lever*(mat[0,0][0,0][1][0][index] - mat[0,0][0,0][1][0][index-1])
        if prop == 'h_sat_l': 
            prop_value = mat[0,0][0,0][3][0][index-1] + lever*(mat[0,0][0,0][3][0][index] - mat[0,0][0,0][3][0][index-1])
        if prop == 'h_sat_v': 
            prop_value = mat[0,0][0,0][4][0][index-1] + lever*(mat[0,0][0,0][4][0][index] - mat[0,0][0,0][4][0][index-1])
        if prop == 'h_fg': 
            prop_value = mat[0,0][0,0][5][0][index-1] + lever*(mat[0,0][0,0][5][0][index] - mat[0,0][0,0][5][0][index-1])


    else:
        raise Exception('Invalid inputs')

    return prop_value
  
    
#OIl function ---> see ian bell thesis 
    
    
# TO BE IMPLEMENTED
def propfunction_mixture(mat, oil_rate, prop, T = False, P = False, rho = False, s = False, h = False, Q = False ): 
    oil_rate = 0
    
    prop_value = propfunction(mat, prop, T, P, rho, s, h, Q )
    
    return prop_value
    
    