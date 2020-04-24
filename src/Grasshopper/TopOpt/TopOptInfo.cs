using System;
using System.Drawing;
using Grasshopper.Kernel;

namespace TopOptDMS2019
{
    public class TopOptInfo : GH_AssemblyInfo
    {
        public override string Name
        {
            get
            {
                return "TopOpt-DMS2019";
            }
        }
        public override Bitmap Icon
        {
            get
            {
                //Return a 24x24 pixel bitmap to represent this GHA library.
                return null;
            }
        }
        public override string Description
        {
            get
            {
                //Return a short string describing the purpose of this GHA library.
                return "";
            }
        }
        public override Guid Id
        {
            get
            {
                return new Guid("1774b68a-c40a-4687-9bd9-e7bcd2c33be4");
            }
        }

        public override string AuthorName
        {
            get
            {
                //Return a string identifying you or your company.
                return "Arjen Deetman";
            }
        }
        public override string AuthorContact
        {
            get
            {
                //Return a string representing your preferred contact details.
                return "info@arjendeetman.nl";
            }
        }
    }
}
