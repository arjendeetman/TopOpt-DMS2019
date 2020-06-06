// This file is part of TopOpt-DMS2019. TopOpt-DMS2019 is licensed 
// under the terms of GNU General Public License as published by the 
// Free Software Foundation. For more information and the LICENSE file, 
// see <https://github.com/arjendeetman/TopOpt-DMS2019>.

using System;
using System.Collections.Generic;
using Grasshopper.Kernel;

namespace TopOptDMS2019.Components
{
    /// <summary>
    /// Advanced settings component. 
    /// </summary>
    public class AdvancedSettingsComponent : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the component.
        /// </summary>
        public AdvancedSettingsComponent()
          : base("Advanced settings component",
              "SET",
              "This component sets advanced settings for the simulation",
              "DMS2019",
              "Topology Optimization")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddNumberParameter("Penalization factor", "P", "The penalization factor for the SIMP law.", GH_ParamAccess.item, 3);
            pManager.AddNumberParameter("Move limit", "M", "The maximum density change between two iterations.", GH_ParamAccess.item, 0.20);
            pManager.AddNumberParameter("Tolerance", "TOL", "Stop criteria for the optimization process.", GH_ParamAccess.item, 0.001);
            pManager.AddIntegerParameter("Maximum iterations", "KMAX", "The maximum number of iterations.", GH_ParamAccess.item, 1000);
            pManager.AddGenericParameter("Filtering approach", "FI", "Choose the filtering approach (sensitivity or density filtering).", GH_ParamAccess.item);
            pManager.AddBooleanParameter("Continuation", "C", "Set True value to enable continuation strategy", GH_ParamAccess.item, false);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddNumberParameter("Advanced settings", "SET", "Advanced settings for the topology optimization simulation", GH_ParamAccess.list);
        }

        // Fields
        private bool _createFilteringList = false;

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // Declare variables
            double penalty = new double();
            double move = new double();
            double tol = new double();
            int kmax = new int();
            int filtering = new int();
            bool continuation = new bool();
            List<double> settings = new List<double>();

            // Create value lists
            GetFilteringList();
            if (this.Params.Input[4].SourceCount != 0) _createFilteringList = false;

            // Access the input parameters
            if (!DA.GetData(0, ref penalty)) return;
            if (!DA.GetData(1, ref move)) return;
            if (!DA.GetData(2, ref tol)) return;
            if (!DA.GetData(3, ref kmax)) return;
            if (!DA.GetData(4, ref filtering)) return;
            if (!DA.GetData(5, ref continuation)) return;

            // Filtering type list
            if (filtering != 1 && filtering != 2)
            {
                _createFilteringList = true;
                ExpireSolution(true);
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
                AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "An unrealistic small value for the tolerance is used. Recommend values are 0.01 and 0.001");
            }
            else if (tol > 0.5)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "An unrealistic high value for the tolerance is used. Recommend values are 0.01 and 0.001");
            }
            // kmax
            if (kmax < 0)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "A negative value for the maximum number of iterations is not allowed.");
            }
            // kmax and continuation
            if (kmax < 70 && continuation == true)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "The continuation strategy is enabled, which needs a mininum of 70 iterations. The maximum number of iterations is set to 70. ");
                kmax = 70;
            }
            // move limit and tolerance
            if (move < tol)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "It is not allowed to have a move limit set smaller than the tolerance.");
            }

            // Add values to list
            settings.Add(penalty);
            settings.Add(move);
            settings.Add(tol);
            settings.Add(kmax);
            settings.Add(filtering);
            if (continuation == true)
            {
                settings.Add(1);
            }
            else
            {
                settings.Add(0);
            }

            // Assign the output parameters
            DA.SetDataList(0, settings);
        }


        /// <summary>
        /// Function generation of value list for filtering types
        /// </summary>
        private void GetFilteringList()
        {
            if (_createFilteringList == true)
            {
                var par = this.Params.Input[4];

                foreach (var activeObj in OnPingDocument().ActiveObjects())
                {
                    if (par.DependsOn(activeObj) == true)
                    {
                        var obj = activeObj as Grasshopper.Kernel.Special.GH_ValueList;

                        obj.Name = "Filtering approach";
                        obj.NickName = "Filtering approach";
                        obj.ListMode = Grasshopper.Kernel.Special.GH_ValueListMode.DropDown;

                        obj.ListItems.Clear();

                        var item1 = new Grasshopper.Kernel.Special.GH_ValueListItem("SENS", "1");
                        var item2 = new Grasshopper.Kernel.Special.GH_ValueListItem("DENS", "2");

                        obj.ListItems.Add(item1);
                        obj.ListItems.Add(item2);

                        obj.ExpirePreview(true);
                        obj.ExpireSolution(true);

                        _createFilteringList = false;
                    }
                }
            }
        }

        /// <summary>
        /// Provides an Icon for the component.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get { return Properties.Resources.professor_icon; }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("b84d5dc9-0dfc-4f5f-88e7-f4d34534212f"); }
        }
    }
}