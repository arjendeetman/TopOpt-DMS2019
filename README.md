# TopOpt-DMS2019
![TopologyOptimization](/Utility/TopologyOptimization.gif)

TopOpt-DMS2019 is a topology optimization plug-in for Rhinoceros Grasshopper. This plug-in was used in the 'Robotic Wood Printing' workshop at the Design Modeling Symposium in Berlin. The plug-in is made open source with as purpose that the results that are published in the related conference paper can be reproduced. The data and log files of the results as prestened in this paper can be found in the folder with [example files](https://github.com/arjendeetman/TopOpt-DMS2019/tree/master/ExampleFiles). The conference paper can be found [here](https://link.springer.com/chapter/10.1007/978-3-030-29829-6_36).

TopOpt2019-DMS is software as it is. The plug-in can be used, however, I do not garuantee that the software will be maintained and bugs will be fixed. I only made it available to make it easier for other people to reproduce my work. 

## TopOpt DMS2019 installation
#### Grasshopper installation
- Download the latest release file from the [release page](https://github.com/arjendeetman/TopOpt-DMS2019/releases). 
- Unzip the release file.
- Place the folder "DMS2019" in the Grasshopper Libraries folder (File -> Special Folder -> Components)
    - Keep the folder structure intact in order that the following paths, folders and files exists
	```
	/Grasshopper/Libraries/DMS2019/DEN 
	/Grasshopper/Libraries/DMS2019/INPUT
	/Grasshopper/Libraries/DMS2019/LOG
	/Grasshopper/Libraries/DMS2019/MAIN.py
	/Grasshopper/Libraries/DMS2019/TOPOPT.py
	/Grasshopper/Libraries/DMS2019/TopOptDMS2019.gha
	```
- Unblock the files: right-click on the file and select properties from the menu. Click unblock on the general tab. 
- Restart Rhino

#### Python 3.x installation and libraries
Inside Grasshopper an external python program is launched. Python 3.x needs to be installed. Make sure python is installed and can be launched from the windows command line. Make sure that the libraries `numpy` and `scipy` are installed. After you installed python, you can install `numpy` and `scipy` from the command line with the following commands: 

```
pip install numpy
pip install scipy
```

By default the Grasshopper component runs the topology optimization program with the command `python` from the windows command line. If you launch your python installation with a different command, e.g. `python36`, you can change this by setting `python36` as last input for the component that runs the topolgy optimization program. An example is shown in the figure below. 

<img src="/Utility/gh_python_call.png" width="50%" height="50%" />

#### Additional tools (optional)
For visualization of the results the plug-in Human can be helpful. This plug-in can be downloaded from food4rhino: https://www.food4rhino.com/app/human. Components of this plug-in are used in the [example files](https://github.com/arjendeetman/TopOpt-DMS2019/tree/master/ExampleFiles).

## License
TopOpt-DMS2019

Copyright (c) 2019-2020 Arjen Deetman

TopOpt-DMS219 is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

TODO

## Credits
Significant parts of TopOpt-DMS2019 have been developed by Arjen Deetman as part of a research project on continuous timber fibre placement at [EDEK Uni Kassel](https://edek.uni-kassel.de/).

## References
Aage,  N.,  Johansen,  V.E.  (2013).  A  165 LINE  TOPOLOGY  OPTIMIZATION  CODE.  Retrieved  November  2,  2019  from 
http://www.topopt.mek.dtu.dk/Apps-and-software/Topology-optimization-codes-written-in-Python 

[Andreassen, E., Clausen, A., Schevenels, M., Lazarov, B.S., Sigmund, O. (2011). Efficient topology optimization in MATLAB 
using 88 lines of code. Structural and Multidisciplinary Optimization 43. 1-16. doi:10.1007/s00158-010-0594-7](https://link.springer.com/article/10.1007/s00158-010-0594-7)

[Dawod M. et al. (2020) Continuous Timber Fibre Placement. In: Gengnagel C., Baverel O., Burry J., Ramsgaard Thomsen M., Weinzierl S. (eds) Impact: Design With All Senses. DMSB 2019. Springer, Cham](https://link.springer.com/chapter/10.1007/978-3-030-29829-6_36)
