#########################
### LICENSE
#########################

"""
This file is part of TopOpt-DMS2019. TopOpt-DMS2019 is licensed 
under the terms of GNU General Public License as published by the 
Free Software Foundation. For more information and the LICENSE file, 
see <https://github.com/arjendeetman/TopOpt-DMS2019>.
"""

#########################
### LOADING MODULUES
#########################

from __future__ import division
from scipy.sparse.linalg import spsolve
from scipy.sparse import coo_matrix
from math import pi
import numpy as np
import logging
import time
import sys
import os
import gc


#########################
### DEFINITIONS FEM
#########################

# Global stiffness matrix
def MaStiffnMatTruss(_KE, _iK, _jK, _ndof, _free, _Emax, _n, _x, _Emin, _pen):
    
    # Modified Young's modules per element
    Ex=_Emin+_x**_pen*(_Emax-_Emin)

    # Apply geometry and material related properties of the finite element
    Ecal=(np.ones(_n)*Ex)[np.newaxis]
    
    # Setup global stiffness matrix K as an aparse matrix
    sK=(Ecal*_KE.T).flatten(order='F')
    K=coo_matrix((sK,(_iK,_jK)),shape=(_ndof,_ndof)).tocsc()
    K=(K+K.T)*0.5

    # Preconditioner for iterative solver (incomplete LU decomposition)
    Mit=None

    # Return values
    return K, Mit, Ex

# Global stiffness matrix
def MaStiffnMatBeam(_KE, _iK, _jK, _ndof, _free, _Emax, _n, _den, _Emin, _pen, _fil, _con):

    # Modified Young's modules per element
    ## All (power law approach)
    Ex=_Emin+_den**_pen*(_Emax-_Emin)
    ## Overwrite connection
    Ex[_con]=_Emin[_con]+_den[_con]**_pen*(_Emax[_con]-_Emin[_con])

    # Apply geometry and material related properties of the finite element
    Ecal=(np.ones(_n)*Ex)[np.newaxis]
    
    # Setup global stiffness matrix K as an aparse matrix
    sK=(Ecal*_KE.T).flatten(order='F')
    K=coo_matrix((sK,(_iK,_jK)),shape=(_ndof,_ndof)).tocsc()
    K=(K+K.T)*0.5

    # Preconditioner for iterative solver (incomplete LU decomposition)
    Mit=None

    # Return values
    return K, Mit, Ex

# Construct element stiffness matrix 3D truss
def ElStiffnMatBeam(_cx, _cy, _cz, _A, _Iy, _Iz, _J, _L, _v, _n):
    
    # Construct local element stiffness matrix
    ## Initialize array
    KEl=np.zeros((_n,144))
    ## Unit value for Young's modulus
    E=1.0
    ## Shear modulus
    G=(E)/(2*(1+_v))
    ## k-values
    k1=(12.0*E*_Iy)/(_L**3)
    k2=(4.0*E*_Iy)/(_L)
    k3=(6.0*E*_Iy)/(_L**2)
    k4=(2.0*E*_Iy)/(_L)
    k5=(12.0*E*_Iz)/(_L**3)
    k6=(4.0*E*_Iz)/(_L)
    k7=(6.0*E*_Iz)/(_L**2)
    k8=(2.0*E*_Iz)/(_L)
    k9=(E*_A)/_L
    k10=(G*_J)/_L
    ## 1st row
    KEl[:,0]=k9;        KEl[:,6]=-k9
    ## 2nd row
    KEl[:,13]=k5;       KEl[:,17]=k7;       KEl[:,19]=-k5;      KEl[:,23]=k7
    ## 3rd row
    KEl[:,26]=k1;       KEl[:,28]=-k3;      KEl[:,32]=-k1;      KEl[:,34]=-k3
    ## 4th row
    KEl[:,39]=k10;      KEl[:,45]=-k10
    ## 5th row
    KEl[:,50]=-k3;      KEl[:,52]=k2;       KEl[:,56]=k3;       KEl[:,58]=k4
    ## 6th row
    KEl[:,61]=k7;       KEl[:,65]=k6;       KEl[:,67]=-k7;      KEl[:,71]=k8
    ## 7th row
    KEl[:,72]=-k9;      KEl[:,78]=k9
    ## 8th row
    KEl[:,85]=-k5;      KEl[:,89]=-k7;      KEl[:,91]=k5;       KEl[:,95]=-k7
    ## 9th row
    KEl[:,98]=-k1;      KEl[:,100]=k3;      KEl[:,104]=k1;      KEl[:,106]=k3
    ## 10th row
    KEl[:,111]=-k10;    KEl[:,117]=k10;
    ## 11th row
    KEl[:,122]=-k3;     KEl[:,124]=k4;      KEl[:,128]=k3;      KEl[:,130]=k2
    ## 12th row
    KEl[:,133]=k7;      KEl[:,137]=k8;      KEl[:,139]=-k7;     KEl[:,143]=k6

    # KE-transposed
    ## Initialize
    KElT=np.zeros((_n,144))
    ## Tranposed matrix of KEl
    for i in range(0,12):
        j=i*12
        KElT[:,j:j+12]=KEl[:,i::12].copy()

    # Symmetric
    KEl=0.5*(KEl+KElT)

    # Construct element transformation matrix
    ## Initialize array
    T=np.zeros((_n,144))
    ## R-matrix (only valid if element is in xy plane)
    ### Initialize
    R11=np.zeros(_n)
    R12=np.zeros(_n)
    R13=np.zeros(_n)
    R21=np.zeros(_n)
    R22=np.zeros(_n)
    R23=np.zeros(_n)
    R31=np.zeros(_n)
    R32=np.zeros(_n)
    R33=np.zeros(_n)
    ### Loop up
    ind1=np.where(_cz==-1)
    ind2=np.where(_cz==1)
    ind3=np.where(np.abs(_cz)!=1)
    ### Fill ind1
    R13[ind1]=-1.0
    R22[ind1]=1.0
    R31[ind1]=-1.0
    ### Fill ind2
    R13[ind2]=1.0
    R22[ind2]=1.0
    R31[ind2]=-1.0
    ### Fill ind3
    R11[ind3]=_cx[ind3]
    R12[ind3]=_cy[ind3]
    R13[ind3]=_cz[ind3]
    R21[ind3]=-1*_cy[ind3]/(np.sqrt(_cx[ind3]**2+_cy[ind3]**2))
    R22[ind3]=_cx[ind3]/(np.sqrt(_cx[ind3]**2+_cy[ind3]**2))
    R23[ind3]=0.0
    R31[ind3]=-_cx[ind3]*_cz[ind3]/(np.sqrt(_cx[ind3]**2+_cy[ind3]**2))
    R32[ind3]=-_cy[ind3]*_cz[ind3]/(np.sqrt(_cx[ind3]**2+_cy[ind3]**2))
    R33[ind3]=_cx[ind3]**2/(np.sqrt(_cx[ind3]**2+_cy[ind3]**2))+_cy[ind3]**2/(np.sqrt(_cx[ind3]**2+_cy[ind3]**2))
    ## 1st row
    T[:,0]=R11;         T[:,1]=R12;         T[:,2]=R13
    ## 2nd row
    T[:,12]=R21;        T[:,13]=R22;        T[:,14]=R23
    ## 3rd row
    T[:,24]=R31;        T[:,25]=R32;        T[:,26]=R33
    ## 4th row
    T[:,39]=R11;        T[:,40]=R12;        T[:,41]=R13
    ## 5th row  
    T[:,51]=R21;        T[:,52]=R22;        T[:,53]=R23
    ## 6th row
    T[:,63]=R31;        T[:,64]=R32;        T[:,65]=R33
    ## 7th row  
    T[:,78]=R11;        T[:,79]=R12;        T[:,80]=R13
    ## 8th row  
    T[:,90]=R21;        T[:,91]=R22;        T[:,92]=R23
    ## 9th row
    T[:,102]=R31;       T[:,103]=R32;       T[:,104]=R33
    ## 10th     
    T[:,117]=R11;       T[:,118]=R12;       T[:,119]=R13
    ## 11th row
    T[:,129]=R21;       T[:,130]=R22;       T[:,131]=R23
    ## 12th row
    T[:,141]=R31;       T[:,142]=R32;       T[:,143]=R33

    # T-transposed
    ## Initialize
    TT=np.zeros((_n,144))
    ## Tranposed matrix of T (TT=T.T)
    for i in range(0,12):
        j=i*12
        TT[:,j:j+12]=T[:,i::12].copy()

    # Transform local element stiffness matrix to global system
    ## Initialize
    dum=np.zeros((_n,144))
    KEg=np.zeros((_n,144))
    ## Matrix multiplication of TT and KEl
    ### i = row index, j = column index
    for i in range(0,12):
        for j in range(0,12):
            dum[:,i*12+j]=(TT[:,i*12:i*12+12]*KEl[:,j::12]).sum(1)
    ## Matrix multiplication of dum and T
    ### i = row index, j = column index
    for i in range(0,12):
        for j in range(0,12):
            KEg[:,i*12+j]=(dum[:,i*12:i*12+12]*T[:,j::12]).sum(1)
    
    # KE-transposed
    ## Initialize
    KEgT=np.zeros((_n,144))
    ## Tranposed matrix of KEg
    for i in range(0,12):
        j=i*12
        KEgT[:,j:j+12]=KEg[:,i::12].copy()

    # Symmetric
    KEg=0.5*(KEg+KEgT)

    # Return values
    return KEg

# Construct element stiffness matrix 3D truss
def ElStiffnMatTruss(_cx, _cy, _cz, _n, _A, _L):
    
    # Unit value for Young's modulus, cross section area and length
    E=1.0

    # Apply geometry and material related properties of the finite element
    EAL=(np.ones(_n)*E*_A/_L)[np.newaxis].T

    # Initialize array
    KE=np.empty((_n,36))

    # Fill array
    ## First row of element stiffness matrix
    KE[:,0]=_cx**2
    KE[:,1]=_cx*_cy
    KE[:,2]=_cx*_cz
    KE[:,3]=-_cx**2
    KE[:,4]=-_cx*_cy
    KE[:,5]=-_cx*_cz
    ## Second row of element stiffness matrix
    KE[:,6]=_cx*_cy
    KE[:,7]=_cy**2
    KE[:,8]=_cy*_cz
    KE[:,9]=-_cx*_cy
    KE[:,10]=-_cy**2
    KE[:,11]=-_cy*_cz
    ## Third row of element stiffness matrix
    KE[:,12]=_cx*_cz
    KE[:,13]=_cy*_cz
    KE[:,14]=_cz**2
    KE[:,15]=-_cx*_cz
    KE[:,16]=-_cy*_cz
    KE[:,17]=-_cz**2
    ## Fourth row of element stiffness matrix
    KE[:,18]=-_cx**2
    KE[:,19]=-_cx*_cy
    KE[:,20]=-_cx*_cz
    KE[:,21]=_cx**2
    KE[:,22]=_cx*_cy
    KE[:,23]=_cx*_cz
    ## Fifth row of element stiffness matrix
    KE[:,24]=-_cx*_cy
    KE[:,25]=-_cy**2
    KE[:,26]=-_cy*_cz
    KE[:,27]=_cx*_cy
    KE[:,28]=_cy**2
    KE[:,29]=_cy*_cz
    ## Sixth row of element stiffness matrix
    KE[:,30]=-_cx*_cz
    KE[:,31]=-_cy*_cz
    KE[:,32]=-_cz**2
    KE[:,33]=_cx*_cz
    KE[:,34]=_cy*_cz
    KE[:,35]=_cz**2
    ## EA/L
    KE=EAL*KE

    # KE-transposed
    ## Initialize
    KET=np.zeros((_n,36))
    ## Tranposed matrix of KE
    for i in range(0,6):
        j=i*6
        KET[:,j:j+6]=KE[:,i::6].copy()

    # Symmetric
    KE=0.5*(KE+KET)

    # Return values
    return KE


#########################
### DEFINITIONS TOP OPT                                    
#########################

# Function for calculation of objective and constraint function values and sensitivities
def evalfBeam(_K, _U, _F, _ndof, _free, _Mit, _n, _edofMat, _KE, _Ex, _den, _penalty, _Emin, _Emax, _dvdx, \
        _m, _L, _A, _volfrac, _H, _Hs, _lmin, _edof, _fil, _con, _nf, _nc, _x, _mu0, _mu1, _xPhys, _ft):
    
    # Solve displacement field with FEA
    _U[_free,0]=spsolve(_K[_free,:][:,_free],  _F[_free,0])

    # Objective function value and partial derivative
    ## Function value per element and apply scale factor
    f0i=_mu0*(((_U[_edofMat].reshape(_n,_edof,1)*_KE.reshape(_n,_edof,_edof)).sum(1)* _U[_edofMat].reshape(_n,_edof)).sum(1))
    ## Function value f0
    f0val=(_Ex*f0i).sum()
    ## Partial derivates of Young's modulus
    dEdx=_penalty*_den**(_penalty-1)*(_Emax-_Emin)*np.ones(_n)
    ## Partial derivate of objective function f0
    ### Filaments
    df0dx=-dEdx*f0i

    # Constraint function value and partial derivative
    ## Initialize
    fval=np.zeros((_m,1))
    dfdx=np.zeros((_m,_n))
    ## Function values of volume constraints
    fval[0]=_mu1*((_x*_L*_A).sum()/(_L*_A).sum()-_volfrac)
    ## Partial derivates of volume constraints
    dfdx[0]=_mu1*_dvdx.copy()

    # Filtering
    if _lmin>0.0:
        ## Avoid division by zero with gamma
        gamma=0.001
        ## Sensetivity filteirng
        if  _ft=="SENS":
            ## Objective function
            df0dx=np.asarray((_H*(_x*df0dx))[np.newaxis].T/_Hs).flatten()/np.maximum(gamma,_x)
            ## Constraint functions
            ## @@ dfdx[0]=np.asarray((H[fil][:,fil]*(x*dfdx[0]))[np.newaxis].T/Hs[fil]).flatten()/np.maximum(gamma,x)
        ## Density filtering
        elif  _ft=="DENS":
            ## Objective function
            df0dx=np.asarray(_H*(df0dx[np.newaxis].T/_Hs)).flatten()
            ## Costraint function
            dfdx[0]=np.asarray(_H*(dfdx[0][np.newaxis].T/_Hs)).flatten()

    # Return values and gradients of the objective and constraints functions
    return f0val,df0dx,fval,dfdx,_U

# Function for calculation of objective and constraint function values and sensitivities
def evalfTruss(_K, _U, _F, _ndof, _free, _Mit, _n, _edofMat, _KE, _Ex, _x, _penalty, _Emin, \
        _Emax, _dvdx, _m, _L, _A, _volfrac, _H, _Hs, _lmin, _edof, _mu0, _mu1, _xPhys, _ft):
    
    # Solve displacement field with FEA
    _U[_free,0]=spsolve(_K[_free,:][:,_free], _F[_free,0])

    # Objective function value and partial derivative
    ## Function value per element
    f0i=_mu0*(((_U[_edofMat].reshape(_n,_edof,1)*_KE.reshape(_n,_edof,_edof)).sum(1)*_U[_edofMat].reshape(_n,_edof)).sum(1))
    ## Function value f0
    f0val=(_Ex*f0i).sum()
    ## Partial derivates of Young's modulus
    dEdx=_penalty*_x**(_penalty-1)*(_Emax-_Emin)*np.ones(_n)
    ## Partial derivate of objective function f0
    df0dx=-dEdx*f0i

    # Constraint function value and partial derivative
    ## Initialize
    fval=np.zeros((_m,1))
    dfdx=np.zeros((_m,_n))
    ## Function values of volume constraints
    fval[0]=_mu1*((_x*_L*_A).sum()/(_L*_A).sum()-_volfrac)
    ## Partial derivates of volume constraints
    dfdx[0]=_mu1*_dvdx.copy()

    # Sensetivity filtering
    if _lmin>0.0:
        ## Avoid division by zero with gamma
        gamma=0.001
        ## Sensetivity filtering
        if _ft=="SENS":
            ## Objective function
            df0dx=np.asarray((_H*(_x*df0dx))[np.newaxis].T/_Hs).flatten()/np.maximum(gamma,_x)
            ## Constraint functions
            ## @@ dfdx[0]=np.asarray((H*(x*dfdx[0]))[np.newaxis].T/Hs).flatten()/np.maximum(gamma,x)
        ## Density filtering
        elif _ft=="DENS":
            ## Objective function
            df0dx=np.asarray(_H*(df0dx[np.newaxis].T/_Hs)).flatten()
            ## Costraint function
            dfdx[0]=np.asarray(_H*(dfdx[0][np.newaxis].T/_Hs)).flatten()

    # Return values and gradients of the objective and constraints functions
    return f0val,df0dx,fval,dfdx,_U

# Optimality criterion
def oc(_n, _x, _volfrac, _df0dx, _dvdx, _g, _move, _lam, _logger, _passive):

    ## Check input
    count=len(np.where(_df0dx>0)[0])
    if count!=0:
        _logger.warning("Negative sensitivities found in OC method.")
        _logger.warning("{} negative sensitivities are replaced with zero".format(count))
        _df0dx[np.where(_df0dx>0)]=0.0

    ## Initial upper and lower limit of the bi-section method
    l1=0.0
    l2=_lam*4
    
    ## Settings for design variables
    xnew=np.zeros(_n)
    xmin=0.0
    xmax=1.0
    
    ## Find values with bi-section method
    while (l2-l1)/(l1+l2)>1e-4:
        lmid=0.5*(l2+l1)
        B=np.sqrt(np.divide(-_df0dx, (_dvdx*lmid)))
        xnew=np.maximum(xmin,np.maximum(_x-_move, np.minimum(xmax, np.minimum(_x+_move,_x*B))))
        if len(_passive)>0:
            xnew[_passive]=1.0
        _g = 0.0
        gt=_g+np.sum((_dvdx*(xnew-_volfrac)))
        if gt>0: 
            l1=lmid
        else: 
            l2=lmid

    ## Return values
    return xnew,gt,lmid

# Filtering with filter length
def filterLength(_centX, _centY, _centZ, _groups, _lmin, _n, _cx, _cy, _cz, _logger, _remove=True):
 
    # Initialize arrays for constructing the sparse matrix
    iH=np.array([]).astype(int)
    jH=np.array([]).astype(int)
    sH=np.array([]).astype(float)
    # Find all unique groups
    uni=np.array(list(set(_groups))).flatten()
    # Remove passive elements
    try:
        uni.remove(0) # passive elements
    except:
        pass
    # Remove values
    if _remove == True:
        uniIter1=uni[uni>0].astype(int) ## Pos
        uniIter2=uni[uni<0].astype(int) ## Neg
    else:
        uniIter1=uni
    # Filter
    ## Groups >0
    for i in uniIter1:
        ## Find index of group members
        ind1=np.array(np.where(_groups==i)).flatten()
        ## Get coordinates of elements in group
        xg=_centX[ind1][np.newaxis]
        yg=_centY[ind1][np.newaxis]
        zg=_centZ[ind1][np.newaxis]
        ## Calculate distance
        ### 2D array
        dist=_lmin-np.sqrt((xg.T-xg)**2+(yg.T-yg)**2+(zg.T-zg)**2)
        ### 1D array
        val=dist.flatten()[np.where(dist.flatten()>0.0)]
        ## Add values to array
        ind=np.where(dist>0.0)
        ind2a=ind[0]
        ind2b=ind[1]
        ind3a=ind1[ind2a]
        ind3b=ind1[ind2b]
        iH=np.append(iH, ind3a).astype(int)
        jH=np.append(jH, ind3b).astype(int)
        sH=np.append(sH, val).astype(float)
    ## Passive elements: Place a one on the diagonal
    try:
        ind=np.array(np.where(_groups==0)).flatten()
        val=np.ones(len(ind))*_lmin
        iH=np.append(iH, ind).astype(int)
        jH=np.append(jH, ind).astype(int)
        sH=np.append(sH, val).astype(float)
    except:
        pass
    ## Connector elements: Place a one on the diagonal
    if _remove==True:
        for i in uniIter2:
            ind=np.array(np.where(_groups==i)).flatten()
            val=np.ones(len(ind))*_lmin
            iH=np.append(iH, ind).astype(int)
            jH=np.append(jH, ind).astype(int)
            sH=np.append(sH, val).astype(float)

    # Finalize assembly and convert to csc format
    H=coo_matrix((sH,(iH,jH)),shape=(_n,_n)).tocsc()
    Hs=H.sum(1) 

    # Return
    return (H,Hs)

# Calculate scaling factor for objective function
def scalingFactorObjTruss(_KE, _iK, _jK, _ndof, _free, _Emax, _ne, _Emin, _volfrac, _U, _F):
    ## Set density of elements equal to volume fraction
    den=np.ones(_ne)*_volfrac
    ## Construct the global stiffness matrix with p=1
    K,Mit,Ex=MaStiffnMatTruss(_KE,_iK,_jK,_ndof,_free,_Emax,_ne,den,_Emin,1.0)
    ## Solve displacement field
    _U[_free,0]=spsolve(K[_free,:][:,_free], _F[_free,0])
    ## Calculate function value f0
    f0val=(_F*_U).sum()
    ## Calculate scale factor
    mu0=(float(50)/f0val)*1.50
    ## Return value
    return mu0

# Calculate scaling factor for objective function
def scalingFactorObjBeam(_KE, _iK, _jK, _ndof, _free, _Emax, _ne, _Emin, _fil, _con, _volfrac, _U, _F):
    ## Set density of elements equal to volume fraction
    den=np.ones(_ne)*_volfrac
    ## Construct the global stiffness matrix with p=1
    K,Mit,Ex=MaStiffnMatBeam(_KE,_iK,_jK,_ndof,_free,_Emax,_ne,den,_Emin,1.0,_fil,_con)
    ## Solve displacement field
    _U[_free,0]=spsolve(K[_free,:][:,_free], _F[_free,0])
    ## Calculate function value f0
    f0val=(_F*_U).sum()
    ## Calculate scale factor
    mu0=(float(50)/f0val)*1.50
    ## Return value
    return mu0


#########################
### LOGGING              
#########################

# Setup logger
def setupLogger(_file):

    # Create logger
    logger=logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
    
    # Create console handler and set level to debug
    ch=logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    
    # Create file handler and set level to debug
    fh=logging.FileHandler(_file)
    fh.setLevel(logging.DEBUG)
    
    # Add formatter to console and file handler
    formatter=logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    fh.setFormatter(formatter)
    
    # Add ch and fh to logger
    logger.addHandler(ch)
    logger.addHandler(fh)
   
    # Open logfile and reset (clear file if it exsists)
    with open(_file, 'w'): pass
    
    # Return logger
    return logger

# Output densities
def savedens(_x, _path, _k):
    file=os.path.join(_path, r"DEN", r"ITERATION-{}.csv".format(_k))
    np.savetxt(file,_x,fmt='%.3f')


#########################
### MAIN BEAM                                                                                  
#########################

# Main program
def mainBeam(_volfrac, _penalty, _move, _lmin, _tolx, _kmax, _width, _height, _filtering, _freeze, _save, _continuation):

    # Timer
    starttimer=time.process_time()

    #####################
    #####################
    #####################

    # Input
    nameLog="BEAM_Optimization"
    nameGeo="element.csv"
    nameBc="bc.csv"
    nameLoad="load.csv"
    nameFrozen="passive.csv"

    # Material and cross section properties
    ## Filament
    a=_height ## Height     = 1.0 mm
    b=_width ## Width       = 6.0 mm
    E1max=float(1.0) ## Maximum Young's modulus
    E1min=E1max*float(1e-9) ## Minimum Young's modulus
    I1y=1.0/12.0*b*a**3 ## Moment of inertia
    I1z=1.0/12.0*a*b**3 ## Moment of inertia
    J1=a*b**3*(1.0/3.0-0.21*(b/a)*(1-b**4/(12*a**4))) ## Torsion constant
    A1=a*b ## Cross section area
    v1=0.3 ## Poisson's ratio
    ## Glue (connection elements between main layer)
    r=(float(a*b)/pi)**0.5 ## Radius such that the area is the same
    E2max=float(1.0) ## Maximum Young's modulus
    E2min=E2max*float(1e-9) ## Minimum Young's modulus
    I2y=1.0/4.0*r**4 ## Moment of inertia
    I2z=1.0/4.0*r**4 ## Moment of inertia
    J2=pi*r**4/2.0 ## Torsion constant
    A2=r**2*pi ## Cross section area
    v2=0.3 ## Poisson's ratio
    ## Delete
    del a,b,r

    # Loads (unit load magnitude of 1.0 is multiplied with the loadfactor)
    loadfactor=float(1.0)

    # Filtering
    if _filtering == 1:
        ft="SENS"
    if _filtering == 2:
        ft="DENS"


    #####################
    #####################
    #####################

    # Setup logger and folders
    ##  Find path of script
    path=os.path.dirname(os.path.realpath(__file__))
    ## Check if folder exist
    for folder in ["LOG", "DEN"]:
        directory=os.path.join(path, r"{}".format(folder))
        try: os.stat(directory)
        except: os.mkdir(directory)
        del folder, directory
    ## Setup logger
    logFile=os.path.join(path, r"LOG", "{}_{}.log".format(time.strftime("%Y-%m-%d_%H-%M-%S"), nameLog))
    logger=setupLogger(logFile)
    ## First log
    logger.info("START OF THE MAIN PROGRAM\n")

    # Import geometry
    ## File
    inputFile=os.path.join(path, r"INPUT", "{}".format(nameGeo))
    ## Import data
    data=np.genfromtxt(inputFile, delimiter=',')
    ## Set data
    ### Connectivity (nid)
    con1=data[:,0].astype(int)
    con2=data[:,1].astype(int)
    ### Element orientation in x, y and y-directionW
    cx=data[:,2].astype(float)
    cy=data[:,3].astype(float)
    cz=data[:,4].astype(float)
    ### Length
    L=data[:,5].astype(float)
    ### Total number of elements and nodes
    ne=int(len(cx))
    nNod=int(np.max([con1, con2]))+1
    ### Coordinates of the center of the element
    centX=data[:,6].astype(float)
    centY=data[:,7].astype(float)
    centZ=data[:,8].astype(float)
    ### Layers and groups
    layers=data[:,9].astype(int)
    groups=data[:,10].astype(int)
    groups[groups>=0]+=1
    ## Delete
    del data, inputFile

    # Import boundary conditions
    ## File
    inputFile=os.path.join(path, r"INPUT", "{}".format(nameBc))
    ## Import data
    data=np.genfromtxt(inputFile, delimiter=',')
    ## Degrees of freedom
    bcdof=np.array(data).flatten().astype(int)
    ## Delete
    del data,inputFile

    # Import loads
    ## File
    inputFile=os.path.join(path, r"INPUT", "{}".format(nameLoad))
    ## Import data
    data=np.genfromtxt(inputFile, delimiter=',')
    ## Degrees of freedom and load magnitude
    if data.size == 2:
        loaddof=int(data[0])
        load=data[1]*loadfactor
    else:
        loaddof=np.array(data[:,0]).flatten().astype(int)
        load=np.array(data[:,1]).flatten().astype(float)*loadfactor
    ## Delete
    del data, inputFile

    # Read frozen
    if _freeze==1:
        try:
            ## File
            inputFile=os.path.join(path, r"INPUT", "{}".format(nameFrozen))
            ## Import data
            frozen=np.genfromtxt(inputFile, delimiter=',').flatten().astype(int)
            ## Delete
            del inputFile
            ## Set group to 0
            groups[frozen]=0
        except:
            frozen=np.array([])
    else:
        frozen=np.array([])

    # Initialize FEA
    ## Total number of degrees of freedom
    ndof=int(nNod*6) ## Structure
    edof=int(2*6) ## Element
    ## Initialze vectors for forces and displacements
    F=np.zeros((ndof,1))
    U=np.zeros((ndof,1))
    ## Set load
    F[loaddof,0]=load
    ## Free degrees of freedom
    free=np.setdiff1d(np.arange(0,ndof),bcdof)

    # Set array's with material and section properties
    ### Index of filament elements
    fil=np.where(groups>=0)[0].flatten()
    ### Set property
    Emax=np.ones(ne)*E1max
    Emin=np.ones(ne)*E1min
    Iy=np.ones(ne)*I1y
    Iz=np.ones(ne)*I1z
    J=np.ones(ne)*J1
    A=np.ones(ne)*A1
    v=np.ones(ne)*v1
    ## Glue (connection elements)
    ### Index of connection elements
    con=np.where(groups<0)[0].flatten()
    ### Set property
    Emax[con]=E2max
    Emin[con]=E2min
    Iy[con]=I2y
    Iz[con]=I2z
    J[con]=J2
    A[con]=A2
    v[con]=v2

    # Meaning of variables
    nf=len(fil)
    nc=len(con)

    # Log input
    logger.info("Total number of finite elements            :   {}".format(ne))
    logger.info("Number of finite elements for filament     :   {}".format(nf))
    logger.info("Number of finite elements for connection   :   {}".format(nc))
    logger.info("Total number of design variables           :   {}".format(nf))
    logger.info("Total number of nodes                      :   {}".format(nNod))
    logger.info("Total degrees of freedom                   :   {}".format(ndof))
    logger.info("Number of free degrees of freedom          :   {}".format(len(free)))
    logger.info("Volume fraction                            :   {}".format(_volfrac))
    logger.info("Penalty factor                             :   {}".format(_penalty))
    logger.info("Continuation step                          :   {}".format(_continuation))
    logger.info("Filter length                              :   {}".format(_lmin))
    logger.info("Filter approach                            :   {}".format(ft))
    logger.info("Maximum number of iterations               :   {}".format(_kmax))
    logger.info("Tolerance (max. density change)            :   {}".format(_tolx))
    logger.info("Maximum Young's modulus for fil (if x_i=1) :   {}".format(E1max))
    logger.info("Minimum Young's modulus for fil (if x_i=0) :   {}".format(E1min))
    logger.info("Maximum Young's modulus for con (if x_i=1) :   {}".format(E2max))
    logger.info("Minimum Young's modulus for con (if x_i=0) :   {}\n".format(E2min))

    # Prepare data for constructing the global stiffness matrix
    ## Element stiffness matrix
    KE=ElStiffnMatBeam(cx,cy,cz,A,Iy,Iz,J,L,v,ne)
    ## Connectivity matrix
    edofMat=np.zeros((ne,edof),dtype=int)
    ### Element node 1
    edofMat[:,0]=con1*6
    edofMat[:,1]=con1*6+1
    edofMat[:,2]=con1*6+2
    edofMat[:,3]=con1*6+3
    edofMat[:,4]=con1*6+4
    edofMat[:,5]=con1*6+5
    ### Element node 2
    edofMat[:,6]=con2*6
    edofMat[:,7]=con2*6+1
    edofMat[:,8]=con2*6+2
    edofMat[:,9]=con2*6+3
    edofMat[:,10]=con2*6+4
    edofMat[:,11]=con2*6+5
    ## Construct the index pointers for the coo format
    iK=np.kron(edofMat,np.ones((edof,1))).flatten()
    jK=np.kron(edofMat,np.ones((1,edof))).flatten()

    # Calculate scaling parameters
    ## 1 < f0 < 100 
    mu0=scalingFactorObjBeam(KE,iK,jK,ndof,free,Emax,ne,Emin,fil,con,_volfrac,U,F)
    ## dfi/dx ~= 1
    mu1=float(1)/(A*L/np.sum(A*L)).mean()

    ## Log scaling values
    logger.info("Scaling factor for objective function f0   :   {:.3f}".format(mu0))
    logger.info("Scaling factor for constraint function f1  :   {:.3f}\n".format(mu1))
  
    # Preparation for filtering
    ## Initialize
    H=None
    Hs=None
    ## Filtering the length of the
    if _lmin>0.0: 
       H,Hs=filterLength(centX,centY,centZ,groups,_lmin,ne,cx,cy,cz,logger)

    # Initialize iteration
    ## Design variables x
    if len(frozen)>0:
        xlower=0.0
        xupper=1.0
        while (xupper-xlower)/(xlower+xupper)>1e-4:
            xmid=0.5*(xlower+xupper)
            x=np.ones(ne)*xmid
            x[frozen]=1.0
            vol=(x*L*A).sum()/(L*A).sum()
            if vol<_volfrac:
                xlower = xmid
            else: 
                xupper = xmid
        del xlower, xupper, xmid
    else:
        x=_volfrac*np.ones(ne)
    ## Physical densities
    xPhys=x.copy()
    ## If filtering is applied
    if _lmin>0.0: 
        if ft=='SENS': 
            xPhys=x.copy()
        elif ft=='DENS': 
            xPhys=np.asarray(H*(x[np.newaxis].T)/Hs)[:,0]
	## Partial derivates of element volume
    dvdx=(A*L/np.sum(A*L))
    ## Initialize iteration number and maximum density change
    k=0 ## Iteration number
    change=1.0 ## Maximum density changes between two iterations
    ## Construct the global stiffness matrix by calling the function
    K,Mit,Ex=MaStiffnMatBeam(KE,iK,jK,ndof,free,Emax,ne,xPhys,Emin,_penalty,fil,con)

    # Set solver for design variables x
    m=1 ## Number of constraint functions (dummy value, OC can handle one contraint function)
    g=0.0 ## Must be initialized to use the NGuyen/Paulino OC approach
    xold1=x.copy() ## Design variables one iteration ago
    lam=1e9 ## Initial upper bound for the bi sectio method

    # Continuation step
    if _continuation==True:
        penMin=1.0 ## Minimum penalty factor, default value = 1.0
        penMax=_penalty ## Maximum penalty factor, default value = penalty
        switch=20 ## Number of iterations before increase of penalty factor, default value = 20
        step=1.02 ## Multiplication factor, default value = 1.02
        # Settings for first iteration
        if k<=switch: _penalty=penMin
        else: _penalty=min(penMax,step*pen)
    else:
        penMax=_penalty

    # Iteration k=0
    ## Function values and gradients of the objective and constraint functions
    f0val,df0dx,fval,dfdx,U=evalfBeam(K,U,F,ndof,free,Mit,ne,edofMat,KE,Ex,x,_penalty,Emin,Emax,dvdx,m, \
        L,A,_volfrac,H,Hs,_lmin,edof,fil,con,nf,nc,x,mu0,mu1,xPhys,ft)
    ## Save csv file with initial densities
    if _save == True:
        savedens(xPhys,path,k)
    ## Calculate the current volume fraction
    vol=(xPhys*L*A).sum()/(L*A).sum()
    ## Log
    logger.info("Ite.: {:03}  Vol.: {:.3f}  Ch.: {:.4f}  Obj.: {:.1f}".format(k,vol,0.0,f0val))
    ## Update iteration number
    k+=1

    # Start optimization process	
    while (change>_tolx or _penalty!=penMax) and k<=_kmax:

        # Calculalte new design variables
        ## Save old design variables x
        xold1=x.copy()
        ## Calculate new design variables with OC method
        x,g,lam=oc(nf,x,_volfrac,df0dx,dfdx[0].flatten(),g,_move,lam,logger,frozen)
        ## Update physical densities
        xPhys=x.copy()
        ## If filtering is applied
        if _lmin>0.0:
            ## Sensetivy based filtering
            if ft=='SENS': xPhys=x.copy()
            ## Density based filtering
            elif ft=='DENS': xPhys=np.asarray(H*(x[np.newaxis].T)/Hs).flatten()

        ## Update master stiffness matrix K
        K,Mit,Ex=MaStiffnMatBeam(KE,iK,jK,ndof,free,Emax,ne,xPhys,Emin,_penalty,fil,con)
        ## Re-calculate function values and gradients of the objective and constraint functions
        f0val,df0dx,fval,dfdx,U=evalfBeam(K,U,F,ndof,free,Mit,ne,edofMat,KE,Ex,x,_penalty,Emin,Emax,dvdx,\
            m,L,A,_volfrac,H,Hs,_lmin,edof,fil,con,nf,nc,x,mu0,mu1,xPhys,ft)
        ## Compute maximum density change
        change=np.max(np.abs(x.flatten()-xold1.flatten()))
        ## Save densities
        if _save == True:
            savedens(xPhys,path,k)
        ## Calculate the current volume fraction
        vol=(xPhys*L*A).sum()/(L*A).sum()
        ## Log and print iteration history
        logger.info("Ite.: {:03}  Vol.: {:.3f}  Ch.: {:.4f}  Obj.: {:.1f}".format(k,vol,change,f0val))
        ## Update iteration number
        k+=1
        ## Continuation step
        if _continuation==True and _penalty!=penMax:
            if k<=switch: _penalty=penMin
            else: _penalty=min(penMax,step*_penalty)
        # End of while loop

    # Save and log final results
    ## Save final densities
    savedens(xPhys,path,"FINAL")
    ## Timer
    endtimer=time.process_time()
    timer=np.ceil(endtimer-starttimer)
    ## Final log
    if change<=_tolx: logger.info("Ite.: {:03}  Convergence criterion reached\n".format(k))
    elif k>=_kmax: logger.info("Ite.: {:03}  Maximum number of iterations reached\n".format(k))
    logger.info("FINAL OBJECTIVE FUNCTION VALUE     :   {}\n".format(f0val/mu0))

    # Check values black and white filtering
    logger.info("Black and white filtering data for different threshold values")
    for threshold in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
        ## Copy densitys
        xbw = xPhys.copy()
        ## Apply black and white filter
        ind1=np.where(xPhys<threshold)
        ind2=np.where(xPhys>=threshold)
        xbw[ind1]=0.0
        xbw[ind2]=1.0 
        ## Update master stiffness matrix K
        K,Mit,Ex=MaStiffnMatBeam(KE,iK,jK,ndof,free,Emax,ne,xbw,Emin,_penalty,fil,con)
        ## Get objective function value
        f0val,df0dx,fval,dfdx,U=evalfBeam(K,U,F,ndof,free,Mit,ne,edofMat,KE,Ex,x,_penalty,Emin,Emax,dvdx,\
            m,L,A,_volfrac,H,Hs,_lmin,edof,fil,con,nf,nc,x,mu0,mu1,xPhys,ft)
        ## Calculate volume fraction
        volbw=(xbw*L*A).sum()/(L*A).sum()
        ## Log data
        logger.info("Threshold value            :   {}".format(threshold))
        logger.info("Volume fraction            :   {}".format(volbw))
        logger.info("Objective function value   :   {}\n".format(f0val/mu0))

    # End log
    logger.info("OPTIMIZATION FINISHED IN {:.0f} SEC. /// {:.0f} MIN.\n".format(timer,float(timer)/60))
    logger.info("END OF THE MAIN PROGRAM\n")
    logger.info("YOU CAN CLOSE THIS WINDOW AND CONTINUE IN GRASSHOPPER")


#########################
### MAN TRUSS   
#########################

# Main program
def mainTruss(_volfrac, _penalty, _move, _lmin, _tolx, _kmax, _width, _height, _filtering, _freeze, _save, _continuation):

    # Timer 
    starttimer=time.process_time()

    #####################
    #####################
    #####################

    # Input
    nameLog="TRUSS_Optimization"
    nameGeo="element.csv"
    nameBc="bc.csv"
    nameLoad="load.csv"
    nameFrozen="passive.csv"
    
    # Material and cross section properties
    A=float(_width*_height) ## Cross section area 
    Emax=float(1.0) ## Maximum Young's modulus
    Emin=Emax*float(1e-9) ## Minimum Young's modulus
    
    # Loads (unit loads magnitude of 1.0 is multiplied with the loadfactor)
    loadfactor=float(1.0)

    # Filtering
    if _filtering == 1:
        ft="SENS"
    if _filtering == 2:
        ft="DENS"


    #####################
    #####################
    #####################

    # Setup logger and folders
    ##  Find path of script
    path=os.path.dirname(os.path.realpath(__file__))
    ## Check if folder exist
    for folder in ["LOG", "DEN"]:
        directory=os.path.join(path, r"{}".format(folder))
        try: os.stat(directory)
        except: os.mkdir(directory)
        del folder, directory
    ## Setup logger
    logFile=os.path.join(path, r"LOG", "{}_{}.log".format(time.strftime("%Y-%m-%d_%H-%M-%S"), nameLog))
    logger=setupLogger(logFile)
    ## First log
    logger.info("START OF THE MAIN PROGRAM\n")

    # Import geometry
    ## File
    inputFile=os.path.join(path, r"INPUT", "{}".format(nameGeo))
    ## Import data
    data=np.genfromtxt(inputFile, delimiter=',')
    ## Set data
    ### Connectivity (nid)
    con1=data[:,0].astype(int)
    con2=data[:,1].astype(int)
    ### Element orientation in x, y and y-direction
    cx=data[:,2].astype(float)
    cy=data[:,3].astype(float)
    cz=data[:,4].astype(float)
    ### Length
    L=data[:,5].astype(float)
    ### Total number of elements and nodes
    nEle=int(len(cx))
    nNod=int(np.max([con1, con2]))+1
    ### Coordinates of the center of the element
    centX=data[:,6].astype(float)
    centY=data[:,7].astype(float)
    centZ=data[:,8].astype(float)
    ### Layers and groups
    layers=data[:,9].astype(int)
    groups=data[:,10].astype(int)
    groups[groups>=0]+=1
    ## Delete
    del data, inputFile

    # Import boundary conditions
    ## File
    inputFile=os.path.join(path, r"INPUT", "{}".format(nameBc))
    ## Import data
    data=np.genfromtxt(inputFile, delimiter=',')
    ## Degrees of freedom
    bcdof=np.array(data).flatten().astype(int)
    ## Delete
    del data,inputFile

    # Import loads
    ## File
    inputFile=os.path.join(path, r"INPUT", "{}".format(nameLoad))
    ## Import data
    data=np.genfromtxt(inputFile, delimiter=',')
    ## Degrees of freedom and load magnitude
    if data.size == 2:
        loaddof=int(data[0])
        load=data[1]*loadfactor
    else:
        loaddof=np.array(data[:,0]).flatten().astype(int)
        load=np.array(data[:,1]).flatten().astype(float)*loadfactor
    ## Delete
    del data, inputFile

    # Read frozen
    if _freeze==1:
        try:
            ## File
            inputFile=os.path.join(path, r"INPUT", "{}".format(nameFrozen))
            ## Import data
            frozen=np.genfromtxt(inputFile, delimiter=',').flatten().astype(int)
            ## Delete
            del inputFile
            ## Set group to 0
            groups[frozen]=0
        except:
            frozen=np.array([])
    else:
        frozen=np.array([])

    # Initialize FEA
    ## Total number of degrees of freedom
    ndof=int(nNod*3) ## Structure
    edof=int(2*3) ## Element
    ## Initialze vectors for forces and displacements
    F=np.zeros((ndof,1))
    U=np.zeros((ndof,1))   
    ## Set load
    F[loaddof,0]=load
    ## Free degrees of freedom
    free=np.setdiff1d(np.arange(0,ndof),bcdof)

    # Log input
    logger.info("Number of finite elements                  :   {}".format(nEle))
    logger.info("Number of nodes                            :   {}".format(nNod))
    logger.info("Total degrees of freedom                   :   {}".format(ndof))
    logger.info("Number of free degrees of freedom          :   {}".format(len(free)))
    logger.info("Volume fraction                            :   {}".format(_volfrac))
    logger.info("Penalty factor                             :   {}".format(_penalty))
    logger.info("Continuation step                          :   {}".format(_continuation))
    logger.info("Filter length                              :   {}".format(_lmin))
    logger.info("Filter approach                            :   {}".format(ft))
    logger.info("Maximum number of iterations               :   {}".format(_kmax))
    logger.info("Tolerance (max. density change)            :   {}".format(_tolx))
    logger.info("Maximum Young's modulus (if x_i=1)         :   {}".format(Emax))
    logger.info("Minimum Young's modulus (if x_i=0)         :   {}\n".format(Emin))

    # Construct the global stiffness matrix
    ## Element stiffness matrix
    KE=ElStiffnMatTruss(cx,cy,cz,nEle,A,L)
    ## Connectivity matrix
    edofMat=np.zeros((nEle,edof),dtype=int)
    edofMat[:,0]=con1*3
    edofMat[:,1]=con1*3+1
    edofMat[:,2]=con1*3+2
    edofMat[:,3]=con2*3
    edofMat[:,4]=con2*3+1
    edofMat[:,5]=con2*3+2
    ## Construct the index pointers for the coo format
    iK=np.kron(edofMat,np.ones((edof,1))).flatten()
    jK=np.kron(edofMat,np.ones((1,edof))).flatten()

    # Calculate scaling parameters
    ## 1 < f0 < 100 
    mu0=scalingFactorObjTruss(KE,iK,jK,ndof,free,Emax,nEle,Emin,_volfrac,U,F)
    ## dfi/dx ~= 1
    mu1=float(1)/(A*L/np.sum(A*L)).mean()

    ## Log scaling values
    logger.info("Scaling factor for objective function f0   :   {:.3f}".format(mu0))
    logger.info("Scaling factor for constraint function f1  :   {:.3f}\n".format(mu1))

    # Preparation for filtering
    ## Initialize
    H=None
    Hs=None
    ## Filtering the length of the elements
    if _lmin>0.0: 
       H,Hs=filterLength(centX,centY,centZ,groups,_lmin,nEle,cz,cy,cz,logger)

    # Initialize iteration
    ## Design variables x
    if len(frozen)>0:
        xlower=0.0
        xupper=1.0
        while (xupper-xlower)/(xlower+xupper)>1e-4:
            xmid=0.5*(xlower+xupper)
            x=np.ones(nEle)*xmid
            x[_frozen]=1.0
            vol=(x*L*A).sum()/(L*A).sum()
            if vol<_volfrac:
                xlower = xmid
            else: 
                xupper = xmid
        del xlower, xupper, xmid
    else:
        x=_volfrac*np.ones(nEle)
    ## Physical densities
    xPhys=x.copy()
    ## If filtering is applied
    if _lmin>0.0: 
        if ft=='SENS': 
            xPhys=x.copy()
        elif ft=='DENS':
            xPhys=np.asarray(H*(x[np.newaxis].T)/Hs)[:,0]
    ## Partial derivates of element volume
    dvdx=A*L/np.sum(A*L)
    ## Initialize iteration number and maximum density change
    k=0 ## Iteration number
    change=1.0 ## Maximum density changes between two iterations
    ## Construct the global stiffness matrix by calling the function
    K,Mit,Ex=MaStiffnMatTruss(KE,iK,jK,ndof,free,Emax,nEle,xPhys,Emin,_penalty)

    # Set solver for design variables x
    m=1 ## Number of constraint functions (dummy value, OC can handle one constraint function)
    g=0.0 ## Must be initialized to use the NGuyen/Paulino OC approach
    xold1=x.copy() ## Design variables one iteration ago
    lam=1e9 ## Initial upper bound for the bi sectio method

    # Continuation step
    if _continuation==True:
        penMin=1.0 ## Minimum penalty factor, default value = 1.0
        penMax=_penalty ## Maximum penalty factor, default value = penalty
        switch=20 ## Number of iterations before increase of penalty factor, default value = 20
        step=1.02 ## Multiplication factor, default value = 1.02
        # Settings for first iteration
        if k<=switch: _penalty=penMin
        else: _penalty=min(penMax,step*pen)
    else:
        penMax=_penalty

    # Iteration k=0
    ## Function values and gradients of the objective and constraint functions
    f0val,df0dx,fval,dfdx,U=evalfTruss(K,U,F,ndof,free,Mit,nEle,edofMat,KE,Ex,x,_penalty,Emin,Emax,dvdx,m, \
        L,A,_volfrac,H,Hs,_lmin,edof,mu0,mu1,xPhys,ft)
    ## Save csv file with initial densities
    if _save == True:
        savedens(xPhys,path,k)
    ## Calculate the current volume fraction
    vol=(xPhys*L*A).sum()/(L*A).sum()
    ## Log
    logger.info("Ite.: {:03}  Vol.: {:.3f}  Ch.: {:.4f}  Obj.: {:.1f}".format(k,vol,0.0,f0val))
    ## Update iteration number
    k+=1

    # Start optimization process	
    while (change>_tolx or _penalty!=penMax) and k<=_kmax:

        # Calculate new design variables
        ## Save old design variables x
        xold1=x.copy()
        ## Calculate new design variables with OC method
        x,g,lam=oc(nEle,x,_volfrac,df0dx,dfdx[0].flatten(),g,_move,lam,logger,frozen)
        ## Update physical densities
        xPhys=x.copy()
        ## If filtering is applied
        if _lmin>0.0:
            ## Sensetivy based filtering
            if ft=='SENS': xPhys=x.copy()
            ## Density based filtering
            elif ft=='DENS': xPhys=np.asarray(H*(x[np.newaxis].T)/Hs).flatten()

        ## Update master stiffness matrix K
        K,Mit,Ex=MaStiffnMatTruss(KE,iK,jK,ndof,free,Emax,nEle,xPhys,Emin,_penalty)
        ## Re-calculate function values and gradients of the objective and constraint functions
        f0val,df0dx,fval,dfdx,U=evalfTruss(K,U,F,ndof,free,Mit,nEle,edofMat,KE,Ex,x,_penalty,Emin,Emax,dvdx,m, \
            L,A,_volfrac,H,Hs,_lmin,edof,mu0,mu1,xPhys,ft)
        ## Compute maximum density change
        change=np.max(np.abs(x.flatten()-xold1.flatten()))
        ## Save densities
        if _save == True:
            savedens(xPhys,path,k)
        ## Calculate the current volume fraction
        vol=(xPhys*L*A).sum()/(L*A).sum()
        ## Log and print iteration history
        logger.info("Ite.: {:03}  Vol.: {:.3f}  Ch.: {:.4f}  Obj.: {:.1f}".format(k,vol,change,f0val))
        ## Update iteration number
        k+=1
        ## Continuation step
        if _continuation==True and _penalty!=penMax:
            if k<=switch: _penalty=penMin
            else: _penalty=min(penMax,step*_penalty)
        # End of while loop

    # Save and log final results
    ## Save final densities
    savedens(xPhys,path,"FINAL")
    ## Timer
    endtimer=time.process_time()
    timer=np.ceil(endtimer-starttimer)
    ## Final log
    if change<=_tolx: logger.info("Ite.: {:03}  Convergence criterion reached\n".format(k))
    elif k>=_kmax: logger.info("Ite.: {:03}  Maximum number of iterations reached\n".format(k))
    logger.info("FINAL OBJECTIVE FUNCTION VALUE     :   {}\n".format(f0val/mu0))

    # Check values black and white filtering
    logger.info("Black and white filtering data for different threshold values")
    for threshold in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
        ## Copy densitys
        xbw = xPhys.copy()
        ## Apply black and white filter
        ind1=np.where(xPhys<threshold)
        ind2=np.where(xPhys>=threshold)
        xbw[ind1]=0.0
        xbw[ind2]=1.0 
        ## Update master stiffness matrix K
        K,Mit,Ex=MaStiffnMatTruss(KE,iK,jK,ndof,free,Emax,nEle,xbw,Emin,_penalty)
        ## Get objective function value
        f0val,df0dx,fval,dfdx,U=evalfTruss(K,U,F,ndof,free,Mit,nEle,edofMat,KE,Ex,x,_penalty,Emin,Emax,dvdx,m, \
            L,A,_volfrac,H,Hs,_lmin,edof,mu0,mu1,xbw,ft)
        ## Calculate volume fraction
        volbw=(xbw*L*A).sum()/(L*A).sum()
        ## Log data
        logger.info("Threshold value            :   {}".format(threshold))
        logger.info("Volume fraction            :   {}".format(volbw))
        logger.info("Objective function value   :   {}\n".format(f0val/mu0))

    # Log end 
    logger.info("OPTIMIZATION FINISHED IN {:.0f} SEC. /// {:.0f} MIN.\n".format(timer,float(timer)/60))
    logger.info("END OF THE MAIN PROGRAM\n")
    logger.info("YOU CAN CLOSE THIS WINDOW AND CONTINUE IN GRASSHOPPER")
