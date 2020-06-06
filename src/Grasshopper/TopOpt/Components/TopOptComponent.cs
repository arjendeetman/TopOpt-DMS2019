// This file is part of TopOpt-DMS2019. TopOpt-DMS2019 is licensed 
// under the terms of GNU General Public License as published by the 
// Free Software Foundation. For more information and the LICENSE file, 
// see <https://github.com/arjendeetman/TopOpt-DMS2019>.

using System;
using System.Linq;
using System.Text;
using System.IO;
using System.Collections.Generic;
using Grasshopper.Kernel;
using Rhino.Geometry;
using TopOptDMS2019.Utils;

namespace TopOptDMS2019.Components
{
    /// <summary>
    /// Component tor prepare the topolgy optimization data
    /// </summary>
    public class TopOptComponent : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the component.
        /// </summary>
        public TopOptComponent()
          : base("Preperation for Topology Optimization",
              "PREPERATION",
              "This component prepares the main data (e.g. geometry) for the Topology Optimization program.",
              "DMS2019",
              "Topology Optimization")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddCurveParameter("Elements", "E", "Straight curves that represent the finite elements", GH_ParamAccess.list);
            pManager.AddCurveParameter("Boundary conditions", "BC", "Straight curves that represent the supports of the structure", GH_ParamAccess.list);
            pManager.AddCurveParameter("Forces", "F", "Straight curves that represent the forces acting on the structure", GH_ParamAccess.list);
            pManager.AddGenericParameter("Element type", "ET", "Choose the element type (input a integer: 3 for truss elements or 6 for beam elements.", GH_ParamAccess.item);
            pManager.AddGenericParameter("Working plane", "WP", "Working plane: XY, XZ or YZ-plane.", GH_ParamAccess.item);
            pManager.AddGenericParameter("Dimension", "D", "Dimension: 2D or 3D.", GH_ParamAccess.item);
            pManager.AddTextParameter("Path", "P", "Path for writing the .csv files.", GH_ParamAccess.item, "DEFAULT");
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddPointParameter("Nodes", "N", "Sorted list with Finite Element Nodes", GH_ParamAccess.list);
            pManager.AddCurveParameter("Elements", "E", "Sorted list with Finite Elements", GH_ParamAccess.list);
            pManager.AddNumberParameter("Orientations", "OR", "Index number of the different orientations", GH_ParamAccess.list);
            pManager.AddNumberParameter("Element set", "SE", "The index number of the different element set that form a straight line", GH_ParamAccess.list);
            pManager.AddNumberParameter("Layers", "LA", "The index number of the different layers.", GH_ParamAccess.list);
            pManager.AddTextParameter("Output files", "F", "The external csv files that this component wrote.", GH_ParamAccess.list);
        }

        // Fields
        private bool _createElementTypeList = false;
        private bool _createWorkingPlaneList = false;
        private bool _createDimensionList = false;
        private List<int> _myBcDof = new List<int>();
        private List<int> _myLoadDof = new List<int>();
        private List<double> _myLoadMag = new List<double>();
        private List<int> _nid1 = new List<int>();
        private List<int> _nid2 = new List<int>();
        private List<double> _length = new List<double>();
        private List<double> _coordX = new List<double>();
        private List<double> _coordY = new List<double>();
        private List<double> _coordZ = new List<double>();
        private List<double> _cx = new List<double>();
        private List<double> _cy = new List<double>();
        private List<double> _cz = new List<double>();
        private List<int> _orientation = new List<int>();
        private List<int> _set = new List<int>();
        private List<int> _layer = new List<int>();

        /// <summary>
        /// This is the method that actually does the work (main program).
        /// </summary>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // Declaring variables for input parameters
            List<Curve> myElements = new List<Curve>();
            List<Curve> myBc = new List<Curve>();
            List<Curve> myLoad = new List<Curve>();
            int eDof = new int();
            int workingPlane = new int();
            int dimension = new int();
            string path = String.Empty;
            string pathInput = String.Empty;
            string GhCompFolder = String.Empty;

            // Create value lists
            // Element type
            GetElementTypeList();
            if (this.Params.Input[3].SourceCount != 0) _createElementTypeList = false;
            // Working plane
            GetWorkingPlaneList();
            if (this.Params.Input[4].SourceCount != 0) _createWorkingPlaneList = false;
            // Working plane
            GetDimensionList();
            if (this.Params.Input[5].SourceCount != 0) _createDimensionList = false;

            // Access the input parameters
            if (!DA.GetDataList(0, myElements)) return;
            if (!DA.GetDataList(1, myBc)) return;
            if (!DA.GetDataList(2, myLoad)) return;
            if (!DA.GetData(3, ref eDof)) return;
            if (!DA.GetData(4, ref workingPlane)) return;
            if (!DA.GetData(5, ref dimension)) return;
            if (!DA.GetData(6, ref pathInput)) return;

            // Get path
            if (pathInput == "DEFAULT")
            {
                GhCompFolder = Grasshopper.Folders.DefaultAssemblyFolder;
                path = Path.Combine(GhCompFolder, "DMS2019");
            }
            else
            {
                path = pathInput;
            }
            
            // Count number of elements
            int n = myElements.Count;

            // Declaring variables for output parameter
            List<string> outputFiles = new List<string>();

            // Check input data of element type (should be 3 or 6)
            if (eDof != 3 && eDof != 6)
            {
                _createElementTypeList = true;
                ExpireSolution(true);
            }

            // Check input data of working plane (should be 12 or 13 or 23)
            if (workingPlane != 12 && workingPlane != 13 && workingPlane != 23)
            {
                _createWorkingPlaneList = true;
                ExpireSolution(true);
            }

            // Check input data of dimension should be 2 or 3)
            if (dimension != 2 && dimension != 3)
            {
                _createDimensionList = true;
                ExpireSolution(true);
            }

            // Reset the lists
            _myBcDof.Clear();
            _myLoadDof.Clear();
            _myLoadMag.Clear();
            _nid1.Clear();
            _nid2.Clear();
            _length.Clear();
            _coordX.Clear();
            _coordY.Clear();
            _coordZ.Clear();
            _cx.Clear();
            _cy.Clear();
            _cz.Clear();
            _orientation.Clear();
            _set.Clear();
            _layer.Clear();

            // Sort start and end points of elements
            for (int i = 0; i < n; i++)
            {
                myElements[i] = HelperMethods.ReDrawCurve(myElements[i]);
            }

            // Sort list with elements based on coordinates
            myElements = myElements.OrderBy(i => i.PointAtStart.X).
                ThenBy(i => i.PointAtEnd.X).
                ThenBy(i => i.PointAtStart.Y).
                ThenBy(i => i.PointAtEnd.Y).
                ThenBy(i => i.PointAtStart.Z).
                ThenBy(i => i.PointAtEnd.Z).ToList();

            // Main program
            // Find all unique Finite Element nodes
            List<Point3d> myNodes = FindNodes(myElements);
            // Find degrees of freedom of boundary conditions
            FindBcDof(myNodes, myBc, eDof, workingPlane, dimension);
            // Find degrees of freedom and magninute of applied loads
            FindLoadValues(myNodes, myLoad, eDof);
            // Get element data
            GetElementData(myNodes, myElements);
            // Get sets, orientations and layers
            FindSet(myElements, workingPlane);
            if (dimension == 2) _layer = _orientation.ToList();
            else if (dimension == 3) FindLayer(workingPlane);

            // Write element data to file
            StringBuilder csv1 = new StringBuilder();
            for (int i = 0; i < n; i++)
            {
                string c0 = _nid1[i].ToString();
                string c1 = _nid2[i].ToString();
                string c2 = _cx[i].ToString();
                string c3 = _cy[i].ToString();
                string c4 = _cz[i].ToString();
                string c5 = _length[i].ToString();
                string c6 = _coordX[i].ToString();
                string c7 = _coordY[i].ToString();
                string c8 = _coordZ[i].ToString();
                string c9 = _orientation[i].ToString();
                string c10 = _set[i].ToString();
                string c11 = _layer[i].ToString();
                string newLine = string.Format("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11}", 
                    c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10,c11);
                csv1.AppendLine(newLine);
            }
            string filePath1 = Path.Combine(path, "INPUT", "element.csv");
            File.WriteAllText(filePath1, csv1.ToString());
            outputFiles.Add(filePath1);

            // Write loads data to file
            StringBuilder csv2 = new StringBuilder();
            for (int i = 0; i < _myLoadDof.Count; i++)
            {
                string c0 = _myLoadDof[i].ToString();
                string c1 = _myLoadMag[i].ToString();
                string newLine = string.Format("{0},{1}", c0, c1);
                csv2.AppendLine(newLine);
            }
            string filePath2 = Path.Combine(path, "INPUT", "load.csv");
            File.WriteAllText(filePath2, csv2.ToString());
            outputFiles.Add(filePath2);

            // Write bc data to file
            StringBuilder csv3 = new StringBuilder();
            for (int i = 0; i < _myBcDof.Count; i++)
            {
                string c0 = _myBcDof[i].ToString();
                string newLine = string.Format("{0}", c0);
                csv3.AppendLine(newLine);
            }
            string filePath3 = Path.Combine(path, "INPUT", "bc.csv");
            File.WriteAllText(filePath3, csv3.ToString());
            outputFiles.Add(filePath3);

            // Output
            DA.SetDataList(0, myNodes);
            DA.SetDataList(1, myElements);
            DA.SetDataList(2, _orientation);
            DA.SetDataList(3, _set);
            DA.SetDataList(4, _layer);
            DA.SetDataList(5, outputFiles);
        }

        #region value lists
        /// <summary>
        /// Function generation of value list for element types
        /// </summary>
        void GetElementTypeList()
        {
            if (_createElementTypeList == true)
            {
                var par = this.Params.Input[3];

                foreach (var activeObj in OnPingDocument().ActiveObjects())
                {
                    if (par.DependsOn(activeObj) == true)
                    {
                        var obj = activeObj as Grasshopper.Kernel.Special.GH_ValueList;

                        obj.Name = "Element types";
                        obj.NickName = "Element type";
                        obj.ListMode = Grasshopper.Kernel.Special.GH_ValueListMode.DropDown;

                        obj.ListItems.Clear();

                        var item1 = new Grasshopper.Kernel.Special.GH_ValueListItem("TRUSS", "3");
                        var item2 = new Grasshopper.Kernel.Special.GH_ValueListItem("BEAM", "6");

                        obj.ListItems.Add(item1);
                        obj.ListItems.Add(item2);

                        obj.ExpirePreview(true);
                        obj.ExpireSolution(true);

                        _createElementTypeList = false;
                    }
                }
            }
        }

        /// <summary>
        /// Function generation of value list for element types
        /// </summary>
        private void GetWorkingPlaneList()
        {
            if (_createWorkingPlaneList == true)
            {
                var par = this.Params.Input[4];

                foreach (var activeObj in OnPingDocument().ActiveObjects())
                {
                    if (par.DependsOn(activeObj) == true)
                    {
                        var obj = activeObj as Grasshopper.Kernel.Special.GH_ValueList;

                        obj.Name = "Working planes";
                        obj.NickName = "Working plane";
                        obj.ListMode = Grasshopper.Kernel.Special.GH_ValueListMode.DropDown;

                        obj.ListItems.Clear();

                        var item1 = new Grasshopper.Kernel.Special.GH_ValueListItem("XY-PLANE", "12");
                        var item2 = new Grasshopper.Kernel.Special.GH_ValueListItem("XZ-PLANE", "13");
                        var item3 = new Grasshopper.Kernel.Special.GH_ValueListItem("YZ-PLANE", "23");

                        obj.ListItems.Add(item1);
                        obj.ListItems.Add(item2);
                        obj.ListItems.Add(item3);

                        obj.ExpirePreview(true);
                        obj.ExpireSolution(true);

                        _createWorkingPlaneList = false;
                    }
                }
            }
        }

        /// <summary>
        /// Function generation of value list for element types
        /// </summary>
        private void GetDimensionList()
        {
            if (_createDimensionList == true)
            {
                var par = this.Params.Input[5];

                foreach (var activeObj in OnPingDocument().ActiveObjects())
                {
                    if (par.DependsOn(activeObj) == true)
                    {
                        var obj = activeObj as Grasshopper.Kernel.Special.GH_ValueList;

                        obj.Name = "Dimension";
                        obj.NickName = "Dimension";
                        obj.ListMode = Grasshopper.Kernel.Special.GH_ValueListMode.DropDown;

                        obj.ListItems.Clear();

                        var item1 = new Grasshopper.Kernel.Special.GH_ValueListItem("2D", "2");
                        var item2 = new Grasshopper.Kernel.Special.GH_ValueListItem("3D", "3");

                        obj.ListItems.Add(item1);
                        obj.ListItems.Add(item2);

                        obj.ExpirePreview(true);
                        obj.ExpireSolution(true);

                        _createDimensionList = false;
                    }
                }
            }
        }
        #endregion

        /// <summary>
        /// Find all unique nodes 
        /// </summary>
        private List<Point3d> FindNodes(List<Curve> curves)
        {
            // Initialize list for finite element nodes
            List<Point3d> Node1 = new List<Point3d>();
            List<Point3d> Node2 = new List<Point3d>();

            // Obtaond all the start and end points of the lines
            for (int i = 0; i < curves.Count; i++)
            {
                Node1.Add(curves[i].PointAtStart);
                Node2.Add(curves[i].PointAtEnd);
            }

            // Make one list and obtain all unique nodes
            List<Point3d> allNodes = new List<Point3d>(Node1.Count + Node2.Count);
            allNodes.AddRange(Node1);
            allNodes.AddRange(Node2);
            List<Point3d> uniList = new HashSet<Point3d>(allNodes).ToList();

            // Sort list based on coordinates
            uniList = uniList.OrderBy(i => i.X).ThenBy(i => i.Y).ThenBy(i => i.Z).ToList();

            // Return values
            return uniList;
        }

        /// <summary>
        /// Gives the numbers of the degrees of freedom of the boundary conditions
        /// </summary>
        private void FindBcDof(List<Point3d> allNodes, List<Curve> allBcs, int eDof, int basePlane, int dim)
        {
            // Loop over all curves
            for (int i = 0; i < allBcs.Count; i++)
            {
                // Get points / nodes
                Point3d pt1 = HelperMethods.RoundPoint3d(allBcs[i].PointAtStart);
                Point3d pt2 = HelperMethods.RoundPoint3d(allBcs[i].PointAtEnd);

                // Get direction / vector of the boundary conditions
                Vector3d vector = new Vector3d(pt2 - pt1);

                // Get index of node in list with all nodes
                int index1 = allNodes.FindIndex(j => j == pt1);
                int index2 = allNodes.FindIndex(j => j == pt2);

                // Get degree of freedom
                if (index1 != -1)
                {
                    if (vector[0] != 0) _myBcDof.Add(index1 * eDof + 0);
                    if (vector[1] != 0) _myBcDof.Add(index1 * eDof + 1);
                    if (vector[2] != 0) _myBcDof.Add(index1 * eDof + 2);
                }
                else if (index2 != -1)
                {
                    if (vector[0] != 0) _myBcDof.Add(index2 * eDof + 0);
                    if (vector[1] != 0) _myBcDof.Add(index2 * eDof + 1);
                    if (vector[2] != 0) _myBcDof.Add(index2 * eDof + 2);
                }
                else
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Applied boundary condition is not connected to one of the finite element nodes.");
                }
            }

            // Add  translational degrees of freedom if dimension is 2D
            if (dim == 2)
            {
                if (basePlane == 12) _myBcDof.AddRange(Enumerable.Range(0, allNodes.Count).Select(x => x * eDof + 2));
                else if (basePlane == 13) _myBcDof.AddRange(Enumerable.Range(0, allNodes.Count).Select(x => x * eDof + 1));
                else if (basePlane == 23) _myBcDof.AddRange(Enumerable.Range(0, allNodes.Count).Select(x => x * eDof + 0));
                else AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Wrong input for working plane.");
            }

            // Add rotational degress of freedom if dimension is 2D and element type is beam
            if (dim == 2 && eDof == 6)
            {
                if (basePlane == 12)
                {
                    _myBcDof.AddRange(Enumerable.Range(0, allNodes.Count).Select(x => x * eDof + 3));
                    _myBcDof.AddRange(Enumerable.Range(0, allNodes.Count).Select(x => x * eDof + 4));
                }
                else if (basePlane == 13)
                {
                    _myBcDof.AddRange(Enumerable.Range(0, allNodes.Count).Select(x => x * eDof + 3));
                    _myBcDof.AddRange(Enumerable.Range(0, allNodes.Count).Select(x => x * eDof + 5));
                }
                else if (basePlane == 23)
                {
                    _myBcDof.AddRange(Enumerable.Range(0, allNodes.Count).Select(x => x * eDof + 4));
                    _myBcDof.AddRange(Enumerable.Range(0, allNodes.Count).Select(x => x * eDof + 5));
                }
                else
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Wrong input for working plane and / or element type");
                }

            }

            // Sort list and check for duplicates
            _myBcDof = new HashSet<int>(_myBcDof).ToList();
            _myBcDof.Sort();
        }

        /// <summary>
        /// Gives the numbers of the degrees of freedom of the applied loads and the magnitude of the load
        /// </summary>
        private void FindLoadValues(List<Point3d> allNodes, List<Curve> allLoads, int eDof)
        {
            // Get load magnitude based on line length and degrees of freedom
            for (int i = 0; i < allLoads.Count; i++)
            {
                // Get points / nodes
                Point3d pt1 = HelperMethods.RoundPoint3d(allLoads[i].PointAtStart);
                Point3d pt2 = HelperMethods.RoundPoint3d(allLoads[i].PointAtEnd);

                // Get direction / vector of the applied loads
                Point3d direction = new Point3d(pt2 - pt1);

                // Get index of node in list with all nodes
                int index1 = allNodes.FindIndex(j => j == pt1);
                int index2 = allNodes.FindIndex(j => j == pt2);

                // Get degree of freedom and magnitude
                if (index1 != -1)
                {
                    if (direction[0] != 0)
                    {
                        _myLoadDof.Add(index1 * eDof + 0);
                        _myLoadMag.Add(direction[0]);
                    }
                    if (direction[1] != 0)
                    {
                        _myLoadDof.Add(index1 * eDof + 1);
                        _myLoadMag.Add(direction[1]);
                    }
                    if (direction[2] != 0)
                    {
                        _myLoadDof.Add(index1 * eDof + 2);
                        _myLoadMag.Add(direction[2]);
                    }
                }
                else if (index2 != -1)
                {
                    if (direction[0] != 0)
                    {
                        _myLoadDof.Add(index2 * eDof + 0);
                        _myLoadMag.Add(direction[0]);
                    }
                    if (direction[1] != 0)
                    {
                        _myLoadDof.Add(index2 * eDof + 1);
                        _myLoadMag.Add(direction[1]);
                    }
                    if (direction[2] != 0)
                    {
                        _myLoadDof.Add(index2 * eDof + 2);
                        _myLoadMag.Add(direction[2]);
                    }
                }
                else
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Applied load is not connected to one of the finite element nodes.");
                }
            }

            // Scale load magnitude
            var absoluteMax = _myLoadMag.Select(x => Math.Abs(x)).Max();
            for (int i = 0; i < _myLoadMag.Count; i++)
            {
                _myLoadMag[i] = _myLoadMag[i] / absoluteMax;
            }
        }

        /// <summary>
        /// Gets the geometric data of the finite elements
        /// </summary>
        private void GetElementData(List<Point3d> allNodes, List<Curve> allElements)
        {     
            // Initialize variables
            int n = allElements.Count;
            double x1 = new double();
            double y1 = new double();
            double z1 = new double();
            double x2 = new double();
            double y2 = new double();
            double z2 = new double();
            double L = new double();

            // Get the node number of the start and the end node
            for (int i = 0; i < n; i++)
            {
                // Get points / nodes
                Point3d pt1 = allElements[i].PointAtStart;
                Point3d pt2 = allElements[i].PointAtEnd;

                // Get index
                int index1 = allNodes.FindIndex(j => j == pt1);
                int index2 = allNodes.FindIndex(j => j == pt2);

                // Add index to list
                _nid1.Add(index1);
                _nid2.Add(index2);

                // Add length to list
                L = allElements[i].GetLength();
                _length.Add(L);

                // Add coordinates to list
                x1 = pt1.X;
                y1 = pt1.Y;
                z1 = pt1.Z;
                x2 = pt2.X;
                y2 = pt2.Y;
                z2 = pt2.Z;

                // The center coordinates of the element
                _coordX.Add(0.5 * (x1 + x2));
                _coordY.Add(0.5 * (y1 + y2));
                _coordZ.Add(0.5 * (z1 + z2));

                // Orientatien of the element in the three different directions
                _cx.Add((x2 - x1) / L);
                _cy.Add((y2 - y1) / L);
                _cz.Add((z2 - z1) / L);
            }
        }

        /// <summary>
        /// Gives the set number and the orientation
        /// </summary>
        private void FindSet(List<Curve> allElements, int basePlane)
        {
            // Initialize variables
            int n = allElements.Count;
            List<Vector3d> values = new List<Vector3d>(n);
            _orientation = Enumerable.Repeat(-1, n).ToList();
            _set = Enumerable.Repeat(-1, n).ToList();

            // Get unique layers
            for (int i = 0; i < n; i++)
            {
                values.Add(new Vector3d(_cx[i], _cy[i], _cz[i]));
            }
            List<Vector3d> uniValues = new HashSet<Vector3d>(values).ToList();
            uniValues.Sort();

            // Save positive and negative orientation numbers in a seperate list
            List<int> orientationPos = new List<int>();
            List<int> orientationNeg = new List<int>();

            // Get orientation number
            for (int i = 0; i < n; i++)
            {

                // Get orientation index number
                var val = uniValues.FindIndex(j => j == values[i]);
                val += 1; // Avoid zero values

                // Make the number negative if it is a in between layer
                // Remove connection elements in 3D
                if (basePlane == 12 && _cz[i] != 0)
                {
                    _orientation[i] = -1 * val;
                }
                else if (basePlane == 13 && _cy[i] != 0)
                {
                    _orientation[i] = -1 * val;
                }
                else if (basePlane == 23 && _cx[i] != 0)
                {
                    _orientation[i] = -1 * val;
                }
                else
                {
                    _orientation[i] = val;
                }

                // Save in list for positive and negative values
                if (_orientation[i] > 0)
                {
                    orientationPos.Add(val);
                }
                else
                {
                    orientationNeg.Add(-val);
                }
            }
            
            // Renumber orientation 
            // Get all the unique values
            List<int> uniOriPos = new HashSet<int>(orientationPos).ToList();
            List<int> uniOriNeg = new HashSet<int>(orientationNeg).ToList();
            // Sort all the unique values
            uniOriPos.Sort();
            uniOriNeg.Sort();
            // Change numbering
            for (int i = 0; i < n; i++)
            {
                if (_orientation[i] > 0)
                {
                    _orientation[i] = uniOriPos.FindIndex(j => j == _orientation[i]);
                }
                else
                {
                    _orientation[i] = uniOriNeg.FindIndex(j => j == _orientation[i]);
                    _orientation[i] += 1;
                    _orientation[i] *= -1;
                }
            }

            // Join curves
            // Unique layer numbers
            List<int> uniLayerIndex = new HashSet<int>(_orientation).ToList();
            // Initial list with curves
            List<Curve> joinedCurves = new List<Curve>();
            List<int> joinedCurvesLayer = new List<int>();
            // Loop over layers (orientation groups)
            for (int i = 0; i < uniLayerIndex.Count; i++)
            {
                // Get layer number
                int layerNumber = uniLayerIndex[i];

                // Get index numbers of elements
                List<int> index = Enumerable.Range(0, n).Where(x=> _orientation[x] == layerNumber).ToList();

                // Create empty list
                List<Curve> curves = new List<Curve>();

                // Add elements to list (from index)
                for (int j = 0; j < index.Count; j++)
                {
                    var ind = index[j];
                    curves.Add(allElements[ind]);
                }
                // Join curves
                List<Curve> joined = Curve.JoinCurves(curves).ToList();

                // Add joined curves to list and layer number
                joinedCurves.AddRange(joined);
                joinedCurvesLayer.AddRange(Enumerable.Repeat(layerNumber, joined.Count).ToList());
            }

            // Sort list with joined curves based on coordinates
            var joinedCurvesZip = joinedCurves.Zip(joinedCurvesLayer, (x, y) => new { x, y })
                .OrderBy(pair => pair.x.PointAtStart.X).
                ThenBy(pair => pair.x.PointAtStart.Y).
                ThenBy(pair => pair.x.PointAtStart.Z).
                ThenBy(pair => pair.x.PointAtEnd.X).
                ThenBy(pair => pair.x.PointAtEnd.Y).
                ThenBy(pair => pair.x.PointAtEnd.Z).ToList();
            joinedCurves = joinedCurvesZip.Select(pair => pair.x).ToList();
            joinedCurvesLayer = joinedCurvesZip.Select(pair => pair.y).ToList();

            // Make list for positive and negative set numbers
            List<int> setPos = new List<int>();
            List<int> setNeg = new List<int>();

            // Get group numbers based on joined curves
            for (int i = 0; i < allElements.Count; i++)
            {
                // Loop over the joined curves
                for (int j = 0; j < joinedCurves.Count; j++)
                {
                    if (_orientation[i] == joinedCurvesLayer[j])
                    { 
                        // Check for an intersection
                        var crv1 = allElements[i];
                        var crv2 = joinedCurves[j];
                        var intersect = Rhino.Geometry.Intersect.Intersection.CurveCurve(crv1, crv2, 0.0, 0.0);
            
                        // If there is an intersection the right curve is found: change group number and break loop
                        if (intersect.Count != 0)
                        {   
                            if (_orientation[i] >= 0)
                            {
                                _set[i] = j+1;
                                setPos.Add(_set[i]);
                            }
                            else
                            {
                                _set[i] = -(j+1);
                                setNeg.Add(_set[i]);
                            }
                        break;
                        }
                    }
                }                    
            }

            // Renumber set numbers
            // Get all the unique values
            List<int> uniSetPos = new HashSet<int>(setPos).ToList();
            List<int> uniSetNeg = new HashSet<int>(setNeg).ToList();
            // Sort all the unique values
            uniSetPos.Sort();
            uniSetNeg.Sort();
            // Change numbering
            for (int i = 0; i < n; i++)
            {
                if (_set[i] > 0)
                {
                    _set[i] = uniSetPos.FindIndex(j => j == _set[i]);
                }
                else
                {
                    _set[i] = uniSetNeg.FindIndex(j => j == _set[i]);
                    _set[i] += 1;
                    _set[i] *= -1;
                }
            }

        }

        /// <summary>
        /// Gives the layer number
        /// </summary>
        private void FindLayer(int basePlane)
        {
            // Create intial list
            List<int> uniLayer = new List<int>();
            List<double> coord = new List<double>();
            List<double> uniCoord = new List<double>();

            // Get right coordinate 
            if (basePlane == 12) coord = _coordZ;
            else if (basePlane == 13) coord = _coordY;
            else if (basePlane == 23) coord = _coordX;

            // Unique coords
            uniCoord = new HashSet<double>(coord).ToList();
            uniCoord.Sort();

            // Get layer number
            for (int i = 0; i < _set.Count; i++)
            {
                if (_set[i] >= 0)
                {
                    int index = uniCoord.FindIndex(a => a == coord[i]);
                    _layer.Add(index);
                }
                else
                {
                    _layer.Add(-1);
                }
            }

            // Unique coords
            uniLayer = new HashSet<int>(_layer).ToList();
            uniLayer.Sort();

            // Check if there are connection elemenets (-1)
            bool condition = false;
            if (uniLayer[0] == -1) condition = true;

            // Update list such that the step between two indices is always 1
            if (condition)
                for (int i = 0; i < _layer.Count; i++)
                {
                    int index = uniLayer.FindIndex(a => a == _layer[i]);
                    _layer[i] = index-1;
                }
            else
                for (int i = 0; i < _layer.Count; i++)
                {
                    int index = uniLayer.FindIndex(a => a == _layer[i]);
                    _layer[i] = index;
                }
        }

        /// <summary>
        /// The Exposure property controls where in the panel a component icon 
        /// will appear. There are seven possible locations (primary to septenary), 
        /// each of which can be combined with the GH_Exposure.obscure flag, which 
        /// ensures the component will only be visible on panel dropdowns.
        /// </summary>
        public override GH_Exposure Exposure
        {
            get { return GH_Exposure.primary; }
        }

        /// <summary>
        /// Provides an Icon for every component that will be visible in the User Interface.
        /// Icons need to be 24x24 pixels.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get { return Properties.Resources.prep_icon; }
        }

        /// <summary>
        /// Each component must have a unique Guid to identify it. 
        /// It is vital this Guid doesn't change otherwise old ghx files 
        /// that use the old ID will partially fail during loading.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("2c165f43-9584-43e6-8813-346964d94dd3"); }
        }
    }
}
