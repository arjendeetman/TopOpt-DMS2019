// This file is part of TopOpt-DMS2019. TopOpt-DMS2019 is licensed 
// under the terms of GNU General Public License as published by the 
// Free Software Foundation. For more information and the LICENSE file, 
// see <https://github.com/arjendeetman/TopOpt-DMS2019>.

using System;
using System.Text;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using Grasshopper.Kernel;
using Rhino.Geometry;
using TopOptDMS2019.Utils;

namespace TopOptDMS2019.Components
{
    /// <summary>
    /// Component that identifies which elements are marked as passive.
    /// </summary>
    public class PassiveElements : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the component.
        /// </summary>
        public PassiveElements()
          : base("Passive elements",
              "PASSIVE",
              "This component gets the index numbers of the passive elements.",
              "DMS2019",
              "Topology Optimization")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddCurveParameter("Passive elements", "PA", "Straight curves that represent the finite elements and are passive during the optimization process.", GH_ParamAccess.list);
            pManager.AddCurveParameter("Other elements", "OE", "Other straight curves that represent the finite elements.", GH_ParamAccess.list);
            pManager.AddTextParameter("Path", "P", "Path for writing the .csv file.", GH_ParamAccess.item, "DEFAULT");
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddCurveParameter("Elements", "E", "Sorted list with Finite Elements", GH_ParamAccess.list);
            pManager.AddTextParameter("Index", "I", "The index number of the passive elements (from the sorted list with elements).", GH_ParamAccess.item);
            pManager.AddTextParameter("Output file", "F", "The external csv file that this component wrote.", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // Declaring variables for input parameters
            List<Curve> elements1 = new List<Curve>();
            List<Curve> elements2 = new List<Curve>();
            string path = String.Empty;
            string pathInput = String.Empty;
            string GhCompFolder = String.Empty;

            // Access the input parameters
            if (!DA.GetDataList(0, elements2)) return;
            if (!DA.GetDataList(1, elements1)) return;
            if (!DA.GetData(2, ref pathInput)) return;

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

            // Declaring variables for output parameter
            List<string> outputFiles = new List<string>();

            // Number of elements
            int n1 = elements1.Count;
            int n2 = elements2.Count;
            int n = n1 + n2;

            // Redraw curves
            for (int i = 0; i < n1; i++) elements1[i] = HelperMethods.ReDrawCurve(elements1[i]);
            for (int i = 0; i < n2; i++) elements2[i] = HelperMethods.ReDrawCurve(elements2[i]);

            // Make list with fixed elements value (0 or 1, 1 = fixed)
            List<int> fix = new List<int>();
            for (int i = 0; i < n1; i++) fix.Add(0);
            for (int i = 0; i < n2; i++) fix.Add(1);

            // Make list with all elements
            List<Curve> myElements = new List<Curve>();
            for (int i = 0; i < n1; i++) myElements.Add(elements1[i]);
            for (int i = 0; i < n2; i++) myElements.Add(elements2[i]);

            // Sort the list with elements
            List<Curve> myElementsSorted = new List<Curve>();
            myElementsSorted = myElements.OrderBy(i => i.PointAtStart.X).
                ThenBy(i => i.PointAtEnd.X).
                ThenBy(i => i.PointAtStart.Y).
                ThenBy(i => i.PointAtEnd.Y).
                ThenBy(i => i.PointAtStart.Z).
                ThenBy(i => i.PointAtEnd.Z).ToList();

            // Compare index of items of new list
            List<int> index = new List<int>();
            for (int i = 0; i < n; i++) index.Add(myElements.FindIndex(j => j == myElementsSorted[i]));

            // Sort list with fixed elements (list with 0/1 values)
            List<int> fixSorted = new List<int>();
            for (int i = 0; i < n; i++) fixSorted.Add(fix[index[i]]);

            // Get index number of fixed elements
            List<int> indexFixed = new List<int>();
            for (int i = 0; i < n; i++) if (fixSorted[i] == 1) indexFixed.Add(i);

            // Write to file
            StringBuilder csv = new StringBuilder();
            for (int i = 0; i < indexFixed.Count; i++)
            {
                string c0 = indexFixed[i].ToString();
                string newLine = string.Format("{0}", c0);
                csv.AppendLine(newLine);
            }
            string filePath = Path.Combine(path, "INPUT", "passive.csv");
            File.WriteAllText(filePath, csv.ToString());
            outputFiles.Add(filePath);

            // Assign the output parameters
            DA.SetDataList(0, myElementsSorted);
            DA.SetDataList(1, indexFixed);
            DA.SetDataList(2, outputFiles);
        }

        /// <summary>
        /// Provides an Icon for the component.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get { return Properties.Resources.passive; }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("f8651d2c-5785-4ddf-9604-5e8d90df9028"); }
        }
    }
}