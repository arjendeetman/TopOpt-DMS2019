﻿// This file is part of TopOpt-DMS2019. TopOpt-DMS2019 is licensed 
// under the terms of GNU General Public License as published by the 
// Free Software Foundation. For more information and the LICENSE file, 
// see <https://github.com/arjendeetman/TopOpt-DMS2019>.

using System;
using System.IO;
using System.Collections.Generic;
using Grasshopper.Kernel;

namespace TopOptDMS2019.Components
{
    /// <summary>
    /// Delete saved local files generated by other components
    /// </summary>
    public class DeleteFiles : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the component.
        /// </summary>
        public DeleteFiles()
          : base("Delete",
              "DEL",
              "This component deletes the selected Topology Optimization files.",
              "DMS2019",
              "Utility")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddTextParameter("Files", "F", "Files that should be deleted.", GH_ParamAccess.list);
            pManager.AddBooleanParameter("Delete", "D", "Toggle to delete the selected files.", GH_ParamAccess.item, false);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddTextParameter("Deleted files", "D", "List with files that was deleted", GH_ParamAccess.list);
            pManager.AddTextParameter("Failed", "F", "List with files that could not be deleted", GH_ParamAccess.list);
        }

        //Fields
        private List<string> _deletedFiles = new List<string>();
        private List<string> _failedFiles = new List<string>();

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // Declaring the variables
            bool delete = new bool();
            List<string> files = new List<string>();

            // Access the input parameters
            if (!DA.GetDataList(0, files)) return;
            if (!DA.GetData(1, ref delete)) return;

            // Delete files if the input is toggled
            if (delete == true)
            {
                _deletedFiles.Clear();
                _failedFiles.Clear();

                foreach (string file in files)
                {
                    try
                    {
                        File.Delete(file);
                        _deletedFiles.Add(file);
                    }
                    catch
                    {
                        _failedFiles.Add(file);
                    }
                }
            }

            // Assign the output parameters
            DA.SetDataList(0, _deletedFiles);
            DA.SetDataList(1, _failedFiles);
        }

        /// <summary>
        /// Provides an Icon for the component.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get { return Properties.Resources.red_cross; }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("d6aaeeea-727a-43f6-bac8-a926b983f0c0"); }
        }
    }
}