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
### MAIN PROGRAM
#########################

# Main program
def main():

    # Loading modules
    import numpy as np
    import os

    # Read parameters
    ## Find path of script
    path = os.path.dirname(os.path.realpath(__file__))
    ## File
    inputFile = os.path.join(path, r"INPUT", "parameters.csv")
    ## Import parameters
    parameters = np.genfromtxt(inputFile, delimiter=',')
    ## Set parameters
    volfrac = parameters[0]
    penalty = parameters[1]
    move = parameters[2]
    lmin = parameters[3]
    tolx = parameters[4]
    kmax = parameters[5]
    width = parameters[6]
    height = parameters[7]
    edof = parameters[8]
    filtering = parameters[9]
    freeze = parameters[10]
    if parameters[11] == 1: 
        save = True
    else: 
        save = False
    if parameters[12] == 1: 
        continuation = True
    else: 
        continuation = False

    # Run optimization process
    ## Truss
    if edof == 3:
        from TOPOPT import mainTruss
        mainTruss(volfrac, penalty, move, lmin, tolx, kmax, width, height, filtering, freeze, save, continuation)
    ## Beam
    if edof == 6:
        from TOPOPT import mainBeam
        mainBeam(volfrac, penalty, move, lmin, tolx, kmax, width, height, filtering, freeze, save, continuation)


#########################
### RUN THE MAIN PROGRAM
#########################

# Run the main program
if __name__ == "__main__":
    main()