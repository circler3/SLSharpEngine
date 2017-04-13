using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using TriangleNet;
using TriangleNet.Geometry;
using System.Threading.Tasks;
using TriangleNet.Topology;

namespace IIS.SLSharp.Examples.GeoClipmap
{
    class Triangulation
    {
        private float[,] matrix;

        public Triangulation(ref float[,] _matrix, List<Vertex> list)
        {
            matrix = _matrix;
            Polygon poly = new Polygon();
            //List<Vertex> list = new List<Vertex>();
            //for (int i = 0; i < matrix.GetLength(0); i++)
            //{
            //    for (int j = 0; j < matrix.GetLength(1); j++)
            //    {
            //        list.Add(new Vertex(i, j, matrix[i, j]));
            //    }
            //}
            poly.Add(new Contour(list));
            var mesh = (Mesh)poly.Triangulate(null, null);
            Parallel.ForEach(mesh.Triangles, ProjectionFix);
        }

        Func<double, double, double> max = (a, b) => a > b ? a : b;
        Func<double, double, double> min = (a, b) => a < b ? a : b;

        private void ProjectionFix(Triangle tri)
        {
            // 先计算三角形形成的平面，以便后边计算栅格点的投影高度
            var x1 = tri.GetVertex(0).X;
            var y1 = tri.GetVertex(0).Y;
            var z1 = tri.GetVertex(0).Z;
            var x2 = tri.GetVertex(1).X;
            var y2 = tri.GetVertex(1).Y;
            var z2 = tri.GetVertex(1).Z;
            var x3 = tri.GetVertex(2).X;
            var y3 = tri.GetVertex(2).Y;
            var z3 = tri.GetVertex(2).Z;

            var A = (y2 - y1) * (z3 - z1) - (z2 - z1) * (y3 - y1);
            var B = (z2 - z1) * (x3 - x1) - (x2 - x1) * (z3 - z1);
            var C = (x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1);
            var D = -(A * x1 + B * y1 + C * z1);

            double tri_left = min(min(x1, x2), x3);
            double tri_right = max(max(x1, x2), x3);
            double tri_bottom = min(min(y1, y2), y3);
            double tri_top = max(max(y1, y2), y3);

            //int tri_col_index_start = (int)((tri_left - left) / cell_size);
            //int tri_col_index_end = (int)((tri_right - left) / cell_size);
            //int tri_row_index_start = (int)((tri_bottom - bottom) / cell_size);
            //int tri_row_index_end = (int)((tri_top - bottom) / cell_size);

            int tri_col_index_start = (int)tri_left;
            int tri_col_index_end = (int)tri_right;
            int tri_row_index_start = (int)tri_bottom;
            int tri_row_index_end = (int)tri_top;

            //tri_col_index_start = (tri_col_index_start - 10) >= 0 ? (tri_col_index_start - 10) : 0;
            //tri_col_index_end = (tri_col_index_start + 10) < cell_col ? (tri_col_index_end + 10) : cell_col - 1;
            //tri_row_index_start = (tri_row_index_start - 10) >= 0 ? (tri_row_index_start - 10) : 0;
            //tri_row_index_end = (tri_row_index_end + 10) < cell_row ? (tri_row_index_end + 10) : cell_row - 1;
            for (int i = tri_row_index_start; i <= tri_row_index_end; i+=1)
            {
                for (int j = tri_col_index_start; j <= tri_col_index_end; j+=1)
                {
                    double x = j;
                    double y = i;

                    // 判断点是否落在三角形内
                    double[] plgx = { x1, x2, x3 };
                    double[] plgy = { y1, y2, y3 };
                    int ii = 0;
                    int jj = 2;
                    bool oddNodes = false;

                    for (ii = 0; ii < 3; ii++)
                    {
                        double xi = plgx[ii];
                        double yi = plgy[ii];
                        double xj = plgx[jj];
                        double yj = plgy[jj];
                        if ((yi < y && yj >= y || yj < y && yi >= y) && (xi <= x || xj <= x))
                        {
                            if (xi + (y - yi) / (yj - yi) * (xj - xi) < x)
                            {
                                oddNodes = !oddNodes;
                            }
                        }
                        jj = ii;
                    }

                    if (oddNodes)
                    {
                        double result = -(A * x + B * y + D) / C;
                        matrix[i,j] = (float)result;
                    }

                }
            }

        }
    }
}
