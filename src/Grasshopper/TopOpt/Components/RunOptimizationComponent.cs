using System;
using System.Collections.Generic;
using Grasshopper.Kernel;
using System.Text;
using System.IO;
using System.Diagnostics;

namespace TopOptDMS2019.Components
{
    /// <summary>
    /// Component to run the topology optimization program (calls python in the command line)
    /// </summary>
    public class RunOptimizationComponent : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the component.
        /// </summary>
        public RunOptimizationComponent()
          : base("Run the Topology Optimization program",
              "RUN",
              "This component writes the input parameters and runs the external Topology Optimization program.",
              "DMS2019",
              "Topology Optimization")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddNumberParameter("Volume fraction", "VF", "The volume fraction of the structure.", GH_ParamAccess.item, 0.50);
            pManager.AddNumberParameter("Minimum length", "LMIN", "The minumum length.", GH_ParamAccess.item, 0.0);
            pManager.AddNumberParameter("Cross section width", "W", "The widht of the cross section (in mm).", GH_ParamAccess.item, 4.0);
            pManager.AddNumberParameter("Cross section height", "H", "The height of the cross section (in mm).", GH_ParamAccess.item, 0.6);
            pManager.AddGenericParameter("Element type", "ET", "Choose the element type (input a integer: 3 for truss elements or 6 for beam elements.", GH_ParamAccess.item);
            pManager.AddBooleanParameter("Passive", "PA", "Toggle of the passive elements should be used.", GH_ParamAccess.item, false);
            pManager.AddNumberParameter("Advanced settings", "SET", "Advanced settings for the topology optimization simulation", GH_ParamAccess.list, new List<double>{3, 0.20, 0.001, 1000, 1, 0});
            pManager.AddBooleanParameter("Save densities", "S", "Save the densities of all iterations.", GH_ParamAccess.item, false);
            pManager.AddTextParameter("Path", "P", "Path for writing the .csv files.", GH_ParamAccess.item, "DEFAULT");
            pManager.AddBooleanParameter("Run", "R", "Toggle or button to call the external Topology Optimization program", GH_ParamAccess.item, false);
            pManager.AddTextParameter("Python caller", "PY", "Name for calling the python installation", GH_ParamAccess.item, "python");
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddTextParameter("Output files", "F", "The external csv file that this component wrote.", GH_ParamAccess.item);
        }

        // Fields
        private bool _createElementTypeList = false;

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // Declare variables
            double volumeFraction = new double();
            double penalty;
            double move;
            double lmin = new double();
            double tol;
            double kmax;
            double width = new double();
            double height = new double();
            int eDof = new int();
            double filtering;
            List<double> settings = new List<double>();
            string path;
            string pathInput = String.Empty;
            string GhCompFolder;
            bool run = new bool();
            bool save = new bool();
            int saveOut;
            bool freeze = new bool();
            int freezeOut;
            double continuation;
            string pyCaller = "python";

            // Create value lists
            GetElementTypeList();
            if (this.Params.Input[8].SourceCount != 0) _createElementTypeList = false;

            // Access the input parameters
            if (!DA.GetData(0, ref volumeFraction)) return;
            if (!DA.GetData(1, ref lmin)) return;
            if (!DA.GetData(2, ref width)) return;
            if (!DA.GetData(3, ref height)) return;
            if (!DA.GetData(4, ref eDof)) return;
            if (!DA.GetData(5, ref freeze)) return;
            if (!DA.GetDataList(6, settings)) return;
            if (!DA.GetData(7, ref save)) return;
            if (!DA.GetData(8, ref pathInput)) return;
            if (!DA.GetData(9, ref run)) return;
            if (!DA.GetData(10, ref pyCaller)) return;

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

            // Check input data of element type (should be 3 or 6)
            if (eDof != 3 && eDof != 6)
            {
                _createElementTypeList = true;
                ExpireSolution(true);
            }

            // freeze value
            if (freeze == true) { freezeOut = 1; }
            else { freezeOut = 0; }

            // save value
            if (save == true) { saveOut = 1; }
            else {saveOut = 0; }

            penalty = settings[0];
            move = settings[1];
            tol = settings[2];
            kmax = Math.Round(settings[3], 0);
            filtering = settings[4];
            continuation = settings[5];

            // Check input parameters
            /// volumeFraction
            if (volumeFraction <= 0 || volumeFraction >= 1)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "The volume fraction should be larger than 0 and smaller than 1.");
            }
            /// lmin
            if (lmin < 0)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Negative values for the minimum filtering length are not allowed");
            }
            /// width
            if (width <= 0)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "The width should be a positive value and be unequal to zero.");
            }
            /// height
            if (height <= 0)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "The height should be a positive value and be unequal to zero.");
            }

            // Check input parameters: advanced settings
            /// penalty
            if (penalty < 1)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "The penalization factor should be at least equal to 1");
            }
            else if (penalty > 5)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "An unreaslistic high value for the penalization factor is used. Realistic values are in the range from 1 till 5.");
            }
            /// move
            if (move > 1)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "The move limit should be equal or less then 1. The move limit is changed to 1. The recommended value is 0.2");
                move = 1;
            }
            else if (move < 0)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Negative values for the move limit are not allowed. The move limit should be in the range of 0 till 1. The recommended value is 0.2");
            }
            else if (move == 0)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "A move limit equal to zero is not allowed. The move limit should be in the range of >0 till 1. The recommended value is 0.2");
            }
            // tolerance
            if (tol < 0)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "A negative value for the tolerance is used. The optimization process will stop when the maximum number of iterations kmax is reached.");
            }
            else if (tol < 0.000000001)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "An unrealistic small value for the tolerance is used. Recommended values are 0.01 and 0.001");
            }
            else if (tol > 0.5)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "An unrealistic high value for the tolerance is used. Recommended values are 0.01 and 0.001");
            }
            // kmax
            if (kmax < 0)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "A negative value for the maximum number of iterations is not allowed.");
            }
            // kmax and continuation
            if (kmax < 70 && continuation == 1)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "The continuation strategy is enabled, which needs a mininum of 70 iterations. The maximum number of iterations is set to 70.");
                kmax = 70;
            }
            // move limit and tolerance
            if (move < tol)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "It is not allowed to have a move limit set smaller than the tolerance.");
            }

            // Write parameters
            // CSV builder
            var csv = new StringBuilder();
            // Add lines
            csv.AppendLine(volumeFraction.ToString());
            csv.AppendLine(penalty.ToString());
            csv.AppendLine(move.ToString());
            csv.AppendLine(lmin.ToString());
            csv.AppendLine(tol.ToString());
            csv.AppendLine(kmax.ToString());
            csv.AppendLine(width.ToString());
            csv.AppendLine(height.ToString());
            csv.AppendLine(eDof.ToString());
            csv.AppendLine(filtering.ToString());
            csv.AppendLine(freezeOut.ToString());
            csv.AppendLine(saveOut.ToString());
            csv.AppendLine(continuation.ToString());
            // Write file to path
            var filePath = Path.Combine(path, "INPUT", "parameters.csv");
            File.WriteAllText(filePath, csv.ToString());

            // Run the Python program
            if (run)
            {
                Process proc = new Process();
                proc.StartInfo.FileName = "cmd.exe";
                proc.StartInfo.WorkingDirectory = path;
                proc.StartInfo.UseShellExecute = false;
                proc.StartInfo.RedirectStandardOutput = true;
                proc.StartInfo.RedirectStandardInput = true;
                proc.Start();
                proc.StandardInput.WriteLine(pyCaller + " MAIN.py");
                proc.WaitForExit();
            }

            // Assign the output parameters
            DA.SetData(0, filePath);
        }

        /// <summary>
        /// Function generation of value list for element types
        /// </summary>
        private void GetElementTypeList()
        {
            if (_createElementTypeList == true)
            {
                var par = this.Params.Input[4];

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
        /// Provides an Icon for the component.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get { return Properties.Resources.run_icon; }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("822a55fc-0aba-4412-90d0-ed5556865b98"); }
        }
    }
}