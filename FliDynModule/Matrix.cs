using System;
using System.Text;
//using FliDynModule;

namespace FliDynModule
{
    public class Matrix
    {
        private double[][] a;
        private int m, n;

        #region Constructors
        public Matrix(int m, int n)
        {
            this.m = m;
            this.n = n;
            a = new double[m][];
            for (int i = 0; i < m; i++)
            {
                a[i] = new double[n];
            }
        }
        public Matrix(int m, int n, double s)
        {
            this.m = m;
            this.n = n;
            a = new double[m][];
            for (int i = 0; i < m; i++)
            {
                a[i] = new double[n];
            }
            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    a[i][j] = s;
                }
            }
        }
        public Matrix(double[][] A)
        {
            if (A == null)
                throw new ArgumentNullException("A");
            m = A.Length;
            n = A[0].Length;
            for (int i = 0; i < m; i++)
            {
                if (A[i].Length != n)
                {
                    throw new System.ArgumentException("All rows must have the same length.");
                }
            }
            this.a = A;
        }
        public Matrix(double[,] arr)
        {
            if (arr == null)
                throw new ArgumentNullException("A");
            m = arr.GetLength(0);
            n = arr.GetLength(1);
            this.a = new double[m][];
            for (int i = 0; i < m; i++)
            {
                this.a[i] = new double[n];
                for (int j = 0; j < n; j++)
                {
                    this.a[i][j] = arr[i, j];
                }
            }
        }		
        public Matrix(double[][] A, int m, int n)
        {
            this.a = A;
            this.m = m;
            this.n = n;
        }


        #endregion

        #region Getters, setters and accessors
        public int RowDimension
        {
            get
            {
                return m;
            }
        }
        public int ColumnDimension
        {
            get
            {
                return n;
            }
        }
        public static explicit operator double[,] (Matrix a)
        {
            if (a == null)
                throw new ArgumentNullException("A");
            double[,] X = new double[a.m, a.n];
            for (int i = 0; i < a.m; i++)
            {
                for (int j = 0; j < a.n; j++)
                {
                    X[i, j] = a.a[i][j];
                }
            }
            return X;
        }
        public static explicit operator double[][] (Matrix a)
        {
            if (a == null)
                throw new ArgumentNullException("A");
            return a.a;
        }	
        public Matrix Clone()
        {
            Matrix X = new Matrix(m, n);
            double[][] C = X.a;
            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    C[i][j] = a[i][j];
                }
            }
            return X;
        }
        public Vector CloneColumn(int columnNum)
        {
            if (0 > columnNum || columnNum > ColumnDimension - 1)
                throw new IndexOutOfRangeException("Column index is out of range");
            Vector v = Vector.Zeros(RowDimension);
            for (int i = 0; i < RowDimension; i++)
            {
                v[i] = a[i][columnNum];
            }
            return v;
        }
        public double this[int i, int j]
        {
            get { return a[i][j]; }
            set { a[i][j] = value; }
        }
        public double[] this[int i]
        {
            get { return a[i]; }
            set { a[i] = value; }
        }
        public Matrix Submatrix(int i0, int i1, int j0, int j1)
        {
            Matrix X = new Matrix(i1 - i0 + 1, j1 - j0 + 1);
            double[][] B = X.a;
            try
            {
                for (int i = i0; i <= i1; i++)
                {
                    for (int j = j0; j <= j1; j++)
                    {
                        B[i - i0][j - j0] = a[i][j];
                    }
                }
            }
            catch (System.IndexOutOfRangeException e)
            {
                throw new System.IndexOutOfRangeException("Submatrix indices", e);
            }
            return X;
        }
        public Matrix Submatrix(int[] r, int[] c)
        {
            if (r == null)
                throw new ArgumentNullException("r");
            if (c == null)
                throw new ArgumentNullException("c");
            Matrix X = new Matrix(r.Length, c.Length);
            double[][] B = X.a;
            try
            {
                for (int i = 0; i < r.Length; i++)
                {
                    for (int j = 0; j < c.Length; j++)
                    {
                        B[i][j] = a[r[i]][c[j]];
                    }
                }
            }
            catch (System.IndexOutOfRangeException e)
            {
                throw new System.IndexOutOfRangeException("Submatrix indices", e);
            }
            return X;
        }
        public Matrix Submatrix(int i0, int i1, int[] c)
        {
            if (c == null)
                throw new ArgumentNullException("c");
            Matrix X = new Matrix(i1 - i0 + 1, c.Length);
            double[][] B = X.a;
            try
            {
                for (int i = i0; i <= i1; i++)
                {
                    for (int j = 0; j < c.Length; j++)
                    {
                        B[i - i0][j] = a[i][c[j]];
                    }
                }
            }
            catch (System.IndexOutOfRangeException e)
            {
                throw new System.IndexOutOfRangeException("Submatrix indices", e);
            }
            return X;
        }		
        public Matrix Submatrix(int[] r, int j0, int j1)
        {
            if (r == null)
                throw new ArgumentNullException("r");
            Matrix X = new Matrix(r.Length, j1 - j0 + 1);
            double[][] B = X.a;
            try
            {
                for (int i = 0; i < r.Length; i++)
                {
                    for (int j = j0; j <= j1; j++)
                    {
                        B[i][j - j0] = a[r[i]][j];
                    }
                }
            }
            catch (System.IndexOutOfRangeException e)
            {
                throw new System.IndexOutOfRangeException("Submatrix indices", e);
            }
            return X;
        }
        #endregion

        #region Arithmetic operators

        public static Matrix operator +(Matrix A, Matrix B)
        {
            if (A == null)
                throw new ArgumentNullException("A");
            if (B == null)
                throw new ArgumentNullException("B");
            A.CheckMatrixDimensions(B);
            Matrix X = new Matrix(A.m, A.n);
            double[][] C = X.a;
            for (int i = 0; i < A.m; i++)
            {
                for (int j = 0; j < A.n; j++)
                {
                    C[i][j] = A.a[i][j] + B.a[i][j];
                }
            }
            return X;
        }

        public static Matrix operator -(Matrix A, Matrix B)
        {
            if (A == null)
                throw new ArgumentNullException("A");
            if (B == null)
                throw new ArgumentNullException("B");
            A.CheckMatrixDimensions(B);
            Matrix X = new Matrix(A.m, A.n);
            double[][] C = X.a;
            for (int i = 0; i < A.m; i++)
            {
                for (int j = 0; j < A.n; j++)
                {
                    C[i][j] = A.a[i][j] - B.a[i][j];
                }
            }
            return X;
        }
        public Matrix times(double s)
        {
            Matrix X = new Matrix(m, n);
            double[][] C = X.a;
            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    C[i][j] = s * a[i][j];
                }
            }
            return X;
        }

        public static Matrix operator *(Matrix A, double s)
        {
            if (A == null)
                throw new ArgumentNullException("A");
            Matrix X = new Matrix(A.m, A.n);
            double[][] C = X.a;
            for (int i = 0; i < A.m; i++)
            {
                for (int j = 0; j < A.n; j++)
                {
                    C[i][j] = s * A.a[i][j];
                }
            }
            return X;
        }

        public static Matrix operator *(double s, Matrix A)
        {
            return A * s;
        }
        public Matrix Mul(double s)
        {
            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    a[i][j] = s * a[i][j];
                }
            }
            return this;
        }

        public static Matrix operator *(Matrix A, Matrix B)
        {
            if (B == null)
                throw new ArgumentNullException("B");
            if (B.m != A.n)
            {
                throw new System.ArgumentException("Matrix inner dimensions must agree.");
            }
            Matrix X = new Matrix(A.m, B.n);
            double[][] C = X.a;
            double[] Bcolj = new double[A.n];
            for (int j = 0; j < B.n; j++)
            {
                for (int k = 0; k < A.n; k++)
                {
                    Bcolj[k] = B.a[k][j];
                }
                for (int i = 0; i < A.m; i++)
                {
                    double[] Arowi = A.a[i];
                    double s = 0;
                    for (int k = 0; k < A.n; k++)
                    {
                        s += Arowi[k] * Bcolj[k];
                    }
                    C[i][j] = s;
                }
            }
            return X;
        }

        #endregion

        #region Matrix operations
        public Vector SolveGE(Vector b)
        {
            return Gauss.Solve(this, b);
        }
        
        public static Matrix Identity(int m, int n)
        {
            Matrix A = new Matrix(m, n);
            double[][] X = A.a;
            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    X[i][j] = (i == j ? 1.0 : 0.0);
                }
            }
            return A;
        }
        
        public Matrix Transpose()
        {
            Matrix result = new Matrix(n, m);
            double[][] T = result.a;
            for (int i = 0; i < m; i++)
                for (int j = 0; j < n; j++)
                    T[j][i] = a[i][j];
            return result;
        }
        public Matrix Cholesky()
        {
            Matrix result = new Matrix(m, m);
            var Li = result.a;

            // Main loop
            for (int i = 0; i < n; i++)
            {
                var Lrowi = Li[i];
                for (int j = 0; j < (i + 1); j++)
                {
                    var Lrowj = Li[j];
                    double s = 0;
                    for (int k = 0; k < j; k++)
                        s += Lrowi[k] * Lrowj[k];
                    if (i == j)
                        Lrowi[j] = Math.Sqrt(a[i][i] - s);
                    else
                        Lrowi[j] = (a[i][j] - s) / Lrowj[j];
                }
            }

            return result;
        }
        public Matrix InverseLower()
        {
            int n = this.ColumnDimension;
            var I = Matrix.Identity(n, n);
            var invLtr = new double[n][];
            for (int col = 0; col < n; col++)
            {
                Vector x = Vector.Zeros(n);
                x[col] = 1;
                invLtr[col] = this.SolveLower(x);
            }
            var invL = new Matrix(invLtr).Transpose();

            return invL;
        }

        public Vector SolveLower(Vector b)
        {
            double[] x = new double[m];

            for (int i = 0; i < m; i++)
            {
                x[i] = b[i];
                for (int j = 0; j < i; j++)
                    x[i] -= a[i][j] * x[j];
                x[i] /= a[i][i];
            }

            return new Vector(x);
        }

        public Vector SolveUpper(Vector b)
        {
            double[] x = new double[m];

            for (int i = m - 1; i >= 0; i--)
            {
                x[i] = b[i];
                for (int j = i + 1; j < n; j++)
                    x[i] -= a[i][j] * x[j];
                x[i] /= a[i][i];
            }

            return new Vector(x);
        }
	
        private void CheckMatrixDimensions(Matrix B)
        {
            if (B.m != m || B.n != n)
            {
                throw new System.ArgumentException("Matrix dimensions must agree.");
            }
        }

        #endregion
        public Matrix Inverse3()
        {
            Matrix a = this;
            double[,] inv = new double[3, 3];
            if ((a.ColumnDimension == a.RowDimension && a.ColumnDimension == 3) && ((a[0, 0] * a[1, 1] * a[2, 2] - a[0, 0] * a[1, 2] * a[2, 1] - a[0, 1] * a[1, 0] * a[2, 2] + a[0, 1] * a[1, 2] * a[2, 0] + a[0, 2] * a[1, 0] * a[2, 1] - a[0, 2] * a[1, 1] * a[2, 0]) != 0))
            {
                double[,] inv_temp = {{  (a[1,1]*a[2,2] - a[1,2]*a[2,1])/(a[0,0]*a[1,1]*a[2,2] - a[0,0]*a[1,2]*a[2,1] - a[0,1]*a[1,0]*a[2,2] + a[0,1]*a[1,2]*a[2,0] + a[0,2]*a[1,0]*a[2,1] - a[0,2]*a[1,1]*a[2,0]), -(a[0,1]*a[2,2] - a[0,2]*a[2,1])/(a[0,0]*a[1,1]*a[2,2] - a[0,0]*a[1,2]*a[2,1] - a[0,1]*a[1,0]*a[2,2] + a[0,1]*a[1,2]*a[2,0] + a[0,2]*a[1,0]*a[2,1] - a[0,2]*a[1,1]*a[2,0]),  (a[0,1]*a[1,2] - a[0,2]*a[1,1])/(a[0,0]*a[1,1]*a[2,2] - a[0,0]*a[1,2]*a[2,1] - a[0,1]*a[1,0]*a[2,2] + a[0,1]*a[1,2]*a[2,0] + a[0,2]*a[1,0]*a[2,1] - a[0,2]*a[1,1]*a[2,0])},

                { -(a[1,0]*a[2,2] - a[1,2]*a[2,0])/(a[0,0]*a[1,1]*a[2,2] - a[0,0]*a[1,2]*a[2,1] - a[0,1]*a[1,0]*a[2,2] + a[0,1]*a[1,2]*a[2,0] + a[0,2]*a[1,0]*a[2,1] - a[0,2]*a[1,1]*a[2,0]),  (a[0,0]*a[2,2] - a[0,2]*a[2,0])/(a[0,0]*a[1,1]*a[2,2] - a[0,0]*a[1,2]*a[2,1] - a[0,1]*a[1,0]*a[2,2] + a[0,1]*a[1,2]*a[2,0] + a[0,2]*a[1,0]*a[2,1] - a[0,2]*a[1,1]*a[2,0]), -(a[0,0]*a[1,2] - a[0,2]*a[1,0])/(a[0,0]*a[1,1]*a[2,2] - a[0,0]*a[1,2]*a[2,1] - a[0,1]*a[1,0]*a[2,2] + a[0,1]*a[1,2]*a[2,0] + a[0,2]*a[1,0]*a[2,1] - a[0,2]*a[1,1]*a[2,0])},

                {  (a[1,0]*a[2,1] - a[1,1]*a[2,0])/(a[0,0]*a[1,1]*a[2,2] - a[0,0]*a[1,2]*a[2,1] - a[0,1]*a[1,0]*a[2,2] + a[0,1]*a[1,2]*a[2,0] + a[0,2]*a[1,0]*a[2,1] - a[0,2]*a[1,1]*a[2,0]), -(a[0,0]*a[2,1] - a[0,1]*a[2,0])/(a[0,0]*a[1,1]*a[2,2] - a[0,0]*a[1,2]*a[2,1] - a[0,1]*a[1,0]*a[2,2] + a[0,1]*a[1,2]*a[2,0] + a[0,2]*a[1,0]*a[2,1] - a[0,2]*a[1,1]*a[2,0]),  (a[0,0]*a[1,1] - a[0,1]*a[1,0])/(a[0,0]*a[1,1]*a[2,2] - a[0,0]*a[1,2]*a[2,1] - a[0,1]*a[1,0]*a[2,2] + a[0,1]*a[1,2]*a[2,0] + a[0,2]*a[1,0]*a[2,1] - a[0,2]*a[1,1]*a[2,0])}};
                inv = inv_temp;
            }
            else
                throw new System.ArgumentException(" The Matrix is not Invertible.");
            return new Matrix(inv);
        }
        public double[] Size()
        {
            double[] temp = new double[2];
            temp[0] = this.m; temp[1] = this.n;
            return temp;
        }
        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();
            sb.Append("[");
            if (a != null)
                for (int i = 0; i < this.m; i++)
                {
                    for (int j=0; j<this.n; j++)
                    {
                        if (j > 0) sb.Append("\t");
                        sb.Append(a[i][j]);
                    }
                    if (i != this.m-1) sb.Append("]\n[");
                }
            sb.Append("]");
            return sb.ToString();
        }
    }
}
