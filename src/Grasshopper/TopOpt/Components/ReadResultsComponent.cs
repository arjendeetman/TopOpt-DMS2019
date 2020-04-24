using System;
using System.Collections.Generic;
using Grasshopper.Kernel;
using System.IO;

namespace TopOptDMS2019.Components
{
    public class ReadResultsComponent : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the component.
        /// </summary>
        public ReadResultsComponent()
          : base("Read the Topology Optimization results", 
              "READ",
              "This component reads the result from the Topology Optimization simulation.",
              "DMS2019",
              "Topology Optimization")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddTextParameter("Path", "P", "Path for reading the .csv file.", GH_ParamAccess.item, "DEFAULT");
            pManager.AddBooleanParameter("Run", "R", "Toggle or button to read the external .csv file.", GH_ParamAccess.item, false);
            pManager.AddIntegerParameter("Iteration", "K", "The iteration number. By default the results of the last iteration will be read.", GH_ParamAccess.item, -1);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddNumberParameter("Densities", "D", "The densities of the finite elements.", GH_ParamAccess.list);
            pManager.AddTextParameter("Readed file", "F", "The external csv file that this component was reading.", GH_ParamAccess.item);
        }

        //Fields
        private List<double> _densOut = new List<double>();
        private string _filePath = "";

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // Declaring the variables
            string path;
            string pathInput = String.Empty;
            string GhCompFolder;
            bool run = new bool();
            int iter = new int();

            // Access the input parameters
            if (!DA.GetData(0, ref pathInput)) return;
            if (!DA.GetData(1, ref run)) return;
            if (!DA.GetData(2, ref iter)) return;

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

            // Read the data
            if (run)
            {
                // Clear list
                _densOut.Clear();

                // Read last iteration
                if (iter == -1)
                {
                    _filePath = Path.Combine(path, "DEN", "ITERATION-FINAL.csv");
                    if (File.Exists(_filePath))
                    {
                        var lines = File.ReadAllLines(_filePath);
                        foreach (string i in lines)
                        {
                            _densOut.Add(Convert.ToDouble(i));
                        }
                    }
                }
                else
                {
                    _filePath = Path.Combine(path, "DEN", string.Format("ITERATION-{0}.csv", iter));
                    if (File.Exists(_filePath))
                    {
                        var lines = File.ReadAllLines(_filePath);
                        foreach (string i in lines)
                        {
                            _densOut.Add(Convert.ToDouble(i));
                        }
                    }
                }
            }

            // Assign the output parameters
            DA.SetDataList(0, _densOut);
            DA.SetData(1, _filePath);
        }

        /// <summary>
        /// Provides an Icon for the component.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get { return Properties.Resources.read_icon; }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("305964ce-6a86-448f-9530-acfc3954b9d4"); }
        }
    }
}