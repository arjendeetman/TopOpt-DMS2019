// This file is part of TopOpt-DMS2019. TopOpt-DMS2019 is licensed 
// under the terms of GNU General Public License as published by the 
// Free Software Foundation. For more information and the LICENSE file, 
// see <https://github.com/arjendeetman/TopOpt-DMS2019>.

using System;
using System.Collections.Generic;
using System.Linq;
using Rhino.Geometry;

namespace TopOptDMS2019.Utils
{
    /// <summary>
    /// Helpermethods
    /// </summary>
    public static class HelperMethods
    {
        /// <summary>
        /// Redraws the curve to give the curves the same orientation
        /// </summary>
        /// <param name="crv"> The curve to re draw</param>
        /// <returns> Returns the reconstructed curve. </returns>
        public static Curve ReDrawCurve(Curve crv)
        {
            // Get points
            List<Point3d> points = new List<Point3d>() { };
            points.Add(RoundPoint3d(crv.PointAtStart));
            points.Add(RoundPoint3d(crv.PointAtEnd));

            // Sort points
            points = points.OrderBy(j => j.X).ThenBy(j => j.Y).ThenBy(j => j.Z).ToList();

            // Draw new curve
            Curve newCurve = new Line(points[0], points[1]).ToNurbsCurve();

            // Return values
            return newCurve;
        }

        /// <summary>
        /// Round of point coordinates and return new points and arrays with coordinates
        /// </summary>
        /// <param name="point"> The point </param>
        /// <param name="tolerance"> Tolerance for round off </param>
        /// <returns> Returns the new point. </returns>
        public static Point3d RoundPoint3d(Point3d point, int tolerance = 6)
        {
            // Get coords and round oof
            double x = Math.Round(point.X, tolerance);
            double y = Math.Round(point.Y, tolerance);
            double z = Math.Round(point.Z, tolerance);

            // Create new point
            Point3d newPoint = new Point3d(x, y, z);

            // Return values
            return newPoint;
        }
    }
}
