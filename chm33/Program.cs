using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics.LinearAlgebra;

namespace chm33
{
    class Program
    {
        
        static void Main(string[] args)
        {
            Console.WriteLine("xd");

            Console.WriteLine("kek");
            /*   2 1 2   
            *    1 2 0
            *    2 0 3
            */
            double[,] array = { { 2, 1, 2 }, { 1, 2, 0 }, { 2, 0, 3 } };
            Matrix<double> matrix = Matrix<double>.Build.DenseOfArray(array);
            double[] x00 = { 1, 1, 1 };
            Vector<double> x0 = Vector<double>.Build.DenseOfArray(x00);
            Console.WriteLine("***CHARACTERISTIC VALUE***");
            int k = 1;
            double eps = 0.001;
            double maxvalue(Matrix<double> A, Vector<double> x)
            {
                Vector<double> x1 = A * x;
                double m0 = 0;
                double m1 = x[k - 1] / x[k - 1];
                Vector<double> x2;
                double m2 = 0;
                do
                {
                    x2 = A * x1;
                    m2 = x2[k - 1] / x1[k - 1];
                    x1 = x2;
                    m0 = m1;
                    m1 = m2;
                }
                while (Math.Abs(m1 - m0) > eps);
                return m2;
            }
            Console.WriteLine("Maximum: " + maxvalue(matrix, x0));
            double mnorm = matrix.L1Norm();
            double[,] array1 = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
            Matrix<double> E = Matrix<double>.Build.DenseOfArray(array1);
            Matrix<double> B = mnorm * E - matrix;
            Console.WriteLine("Minimum: " +(mnorm - maxvalue(B, x0)));


            Console.WriteLine("***NEUTON`S***");

            /*
              sin(x+y) - 1.2x = 0
              x*x + y*y = 1
            */

            double f1(double x, double y)
            {
                return Math.Sin(x + y) - 1.2 * x;
            }

            double f2(double x, double y)
            {
                return x * x + y * y - 1;
            }

            double f1derx(double x, double y)
            {
                return Math.Cos(x + y) - 1.2;
            }

            double f1dery(double x, double y)
            {
                return Math.Cos(x + y);
            }

            double f2derx(double x, double y)
            {
                return 2 * x;
            }

            double f2dery(double x, double y)
            {
                return 2 * y;
            }
        
            double[] mz = { 0, 0 };           
            void jacobs(Matrix<double> m, Vector<double> x)
            {
                Vector<double> xx0 = Vector<double>.Build.DenseOfArray(mz);
                double sum;
                int i = 0, j = 0;
                do
                {
                    for (i = 0; i < xx0.Count; i++)
                    {
                        xx0[i] = x[i];
                    }
                    for (i = 0; i < xx0.Count; i++)
                    {
                        sum = 0;
                        for (j = 0; j < xx0.Count; j++)
                        {
                            if (j != i)
                            {
                                sum += m[i, j] / m[i, i] * xx0[j];
                            }
                        }
                        x[i] = m[i, m.ColumnCount - 1] / m[i, i] - sum;
                    }
                } while (Math.Max(Math.Abs((x - xx0)[0]), Math.Abs((x - xx0)[1])) > eps);
            }
            double[,] ma = { { 0, 0, 0 }, { 0, 0, 0 } };          
            Matrix<double> a = Matrix<double>.Build.DenseOfArray(ma);
            Vector<double> z = Vector<double>.Build.DenseOfArray(mz);
            void method(double x, double y)
            {
                int i = 0;               
                do
                {
                    double[,] farray = { { f1derx(x, y), f1dery(x, y), f1(x,y) }, { f2derx(x, y), f2dery(x, y), f2(x, y) } };
                    a = Matrix<double>.Build.DenseOfArray(farray);
                    jacobs(a, z);
                    x = x - z[0];
                    y = y - z[1];                  
                    i++;
                }
                while (Math.Max(Math.Abs(z[0]), Math.Abs(z[1])) >= eps);
                Console.WriteLine("\nx = " + x);
                Console.WriteLine("y = " + y);
                Console.WriteLine("\niterations: {0}", i);
                Console.WriteLine("\nf1(x, y) = {0}", Math.Round(f1(x, y), 4));
                Console.WriteLine("f2(x, y) = {0}", Math.Round(f2(x, y), 4) + 1);
                Console.WriteLine("\nprescision: {0}", eps);
            }
            method(1.25, 0.75);
            Console.ReadKey();         
        }
    }
}
