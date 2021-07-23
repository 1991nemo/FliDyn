﻿using System;
using System.Linq;
using System.Globalization;

namespace FliDynModule
{
    public class SparseMatrix
    {
        const int Delta = 1;

        private int n, m;
        private double[][] items;
        private int[][] indices;
        private int[] count;
        public SparseMatrix(int m, int n)
        {
            this.m = m;
            this.n = n;

            items = new double[m][];
            indices = new int[m][];
            count = new int[m];
        }

        public SparseMatrix(int m, int n, double[][] matrixItems, int[][] matrixIndices, int[] _count)
        {
            this.m = m;
            this.n = n;
            this.items = matrixItems;
            this.indices = matrixIndices;
            this.count = _count;
        }

        public int[] Count
        {
            get { return this.count; }
        }
        
        public int RowDimension
        {
            get { return m; }
        }
        
        public int ColumnDimension
        {
            get { return n; }
        }

        public SparseMatrix Copy()
        {
            var A = new SparseMatrix(m, n);
            for (int i = 0; i < m; i++)
            {
                A.indices[i] = (int[])indices[i].Clone();
                A.items[i] = (double[])items[i].Clone();
                A.count = (int[])count.Clone();
            }
            return A;
        }
        
        public Matrix DenseMatrix()
        {
            var DM = new Matrix(m, n, 0.0);
            for (int i = 0; i < m; i++)
                for (int j = 0; j < count[i]; j++)
                    DM[i, indices[i][j]] = items[i][j];
            return DM;
        }
        public SparseMatrix plus(SparseMatrix B)
        {
            if (B == null)
                throw new ArgumentNullException("B");
            var C = new SparseMatrix(m, n);

            for (int i = 0; i < m; i++)
            {
                if (indices[i] != null)
                {
                    C.indices[i] = new int[count[i]];
                    C.items[i] = new double[count[i]];
                    for (int j = 0; j < count[i]; j++)
                    {
                        C.indices[i][j] = this.indices[i][j];
                        C.items[i][j] = this.items[i][j] + B.items[i][j];
                    }
                }
                C.count[i] = this.count[i];
            }

            return C;
        }
        public static SparseMatrix operator +(SparseMatrix A, SparseMatrix B)
        {
            if (A == null)
                throw new ArgumentNullException("A");
            if (B == null)
                throw new ArgumentNullException("B");
            return A.plus(B);
        }
        public SparseMatrix minus(SparseMatrix B)
        {
            if (B == null)
                throw new ArgumentNullException("B");
            var C = new SparseMatrix(m, n);

            for (int i = 0; i < m; i++)
            {
                if (indices[i] != null)
                {
                    C.indices[i] = new int[count[i]];
                    C.items[i] = new double[count[i]];
                    for (int j = 0; j < count[i]; j++)
                    {
                        C.indices[i][j] = this.indices[i][j];
                        C.items[i][j] = this.items[i][j] - B.items[i][j];
                    }
                }
                C.count[i] = this.count[i];
            }

            return C;
        }
        public static SparseMatrix operator -(SparseMatrix A, SparseMatrix B)
        {
            if (A == null)
                throw new ArgumentNullException("A");
            if (B == null)
                throw new ArgumentNullException("B");
            return A.minus(B);
        }
        public Vector times(Vector v)
        {
            double[] vv = v;
            Vector result = Vector.Zeros(m);
            unchecked 
            {
                for (int i = 0; i < m; i++)
                    if (indices[i] != null)
                    {
                        double s = 0;
                        for (int k = 0; k < count[i]; k++)
                            s += items[i][k] * vv[indices[i][k]];
                        result[i] = s;
                    }
            }
            return result;
        }
        public static Vector operator *(SparseMatrix m, Vector v)
        {
            if (m == null)
                throw new ArgumentNullException("m");
            if ((double[])v == null)
                throw new ArgumentNullException("v");
            return m.times(v);
        }
        public Vector timesRight(Vector v)
        {
            double[] vv = v;
            Vector result = Vector.Zeros(n);
            unchecked
            {
                for (int i = 0; i < n; i++)
                    if (indices[i] != null)
                    {
                        for (int k = 0; k < count[i]; k++)
                            result[indices[i][k]] += vv[i] * items[i][k];
                    }
            }
            return result;
        }
        public static Vector operator *(Vector v, SparseMatrix m)
        {
            if (m == null)
                throw new ArgumentNullException("m");
            if ((double[])v == null)
                throw new ArgumentNullException("v");
            return m.timesRight(v);
        }
        public SparseMatrix times(double s)
        {
            var B = new SparseMatrix(m, n);

            for (int i = 0; i < m; i++)
            {
                if (indices[i] != null)
                {
                    B.indices[i] = new int[count[i]];
                    B.items[i] = new double[count[i]];
                    for (int j = 0; j < count[i]; j++)
                    {
                        B.indices[i][j] = this.indices[i][j];
                        B.items[i][j] = s * this.items[i][j];
                    }
                }
                B.count[i] = this.count[i];
            }
            return B;
        }
        public static SparseMatrix operator *(SparseMatrix A, double s)
        {
            if (A == null)
                throw new ArgumentNullException("A");
            return A.times(s);
        }
        public SparseMatrix Mul(double s)
        {
            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < count[i]; j++)
                {
                    this.items[i][j] *= s;
                }
            }
            return this;
        }
        public SparseMatrix times(SparseMatrix B)
        {
            if (B == null)
                throw new ArgumentNullException("B");
            if (B.m != n)
            {
                throw new System.ArgumentException("Sparse natrix inner dimensions must agree.");
            }

            var C = new SparseMatrix(m, B.n);
            int idx, ii;
            for (int i = 0; i < m; i++)
            {
                if (indices[i] != null)
                    for (int j = 0; j < B.n; j++)
                    {
                        for (int jj = 0; jj < count[i]; jj++)
                        {
                            ii = indices[i][jj];
                            if (B.indices[ii] != null)
                            {
                                idx = Array.BinarySearch(B.indices[ii], 0, B.count[ii], j);
                                if (idx >= 0)
                                    C[i, j] += items[i][jj] * B.items[ii][idx];
                            }
                        }
                    }
            }

            return C;
        }

        public static SparseMatrix operator *(SparseMatrix A, SparseMatrix B)
        {
            if (A == null)
                throw new ArgumentNullException("A");
            if (B == null)
                throw new ArgumentNullException("B");
            return A.times(B);
        }
        public double this[int i, int j]
        {
            get
            {
#if DEBUG
                if (i < 0 || j < 0 || i >= m || j >= n)
                    throw new IndexOutOfRangeException(String.Format(CultureInfo.InvariantCulture,
                        "Element index ({0},{1}) is out of range", i, j));
#endif
                if (indices[i] == null)
                    return 0;
                int jidx = Array.BinarySearch(indices[i], 0, count[i], j);
                if (jidx < 0)
                    return 0;
                else
                    return items[i][jidx];
            }
            set
            {
#if DEBUG
                if (i < 0 || j < 0 || i >= m || j >= n)
                    throw new IndexOutOfRangeException(String.Format(CultureInfo.InvariantCulture,
                        "Element index ({0},{1}) is out of range", i, j));
#endif
                if (indices[i] == null)
                {
                    indices[i] = new int[Delta];
                    items[i] = new double[Delta];
                    indices[i][0] = j;
                    items[i][0] = value;
                    count[i] = 1;
                }
                else
                {
                    int jidx = Array.BinarySearch(indices[i], 0, count[i], j);
                    if (jidx >= 0)
                        items[i][jidx] = value;
                    else
                    {
                        int indexToAdd = ~jidx;
                        if (count[i] >= items[i].Length)
                        {
                            int delta = Math.Min(Delta, n - items[i].Length);
                            int[] newIndices = new int[indices[i].Length + delta];
                            double[] newItems = new double[items[i].Length + delta];
                            Array.Copy(indices[i], newIndices, indexToAdd);
                            Array.Copy(items[i], newItems, indexToAdd);
                            Array.Copy(indices[i], indexToAdd, newIndices, indexToAdd + 1, count[i] - indexToAdd);
                            Array.Copy(items[i], indexToAdd, newItems, indexToAdd + 1, count[i] - indexToAdd);
                            items[i] = newItems;
                            indices[i] = newIndices;
                        }
                        else
                        {
                            Array.Copy(indices[i], indexToAdd, indices[i], indexToAdd + 1, count[i] - indexToAdd);
                            Array.Copy(items[i], indexToAdd, items[i], indexToAdd + 1, count[i] - indexToAdd);
                        }
                        count[i]++;
                        indices[i][indexToAdd] = j;
                        items[i][indexToAdd] = value;
                    }
                }
            }
        }
        public SparseMatrix Transpose()
        {
            var At = new SparseMatrix(this.ColumnDimension, this.RowDimension);

            for (int i = 0; i < this.RowDimension; i++)
                for (int j = 0; j < this[i].count; j++)
                    At[this[i].indices[j], i] = this[i].items[j];

            return At;
        }
        public SparseVector this[int i]
        {
            get
            {
                return new SparseVector(items[i], indices[i], n);
            }
            set
            {
                indices[i] = value.indices;
                items[i] = value.items;
                count[i] = value.Length;
            }
        }
        public void ScaleRow(int i, int j1, int j2, double sf)
        {
            for (int k = 0; k < count[i]; k++)
                if ((indices[i][k] >= j1) && (indices[i][k] <= j2))
                    items[i][k] *= sf;
        }
        
        public void SwitchRows(int i, int j)
        {
            var tempItems = items[i];
            var tempIndices = indices[i];
            var tempCount = count[i];

            items[i] = items[j];
            indices[i] = indices[j];
            count[i] = count[j];

            items[j] = tempItems;
            indices[j] = tempIndices;
            count[j] = tempCount;
        }
        public Vector SolveGE(Vector b) { return Gauss.SolveCore(Copy(), b); }
        public Vector SolveLower(Vector b)
        {
            Vector x = Vector.Zeros(m);

            for (int i = 0; i < m; i++)
            {
                x[i] = b[i];
                var idx = indices[i];
                var its = items[i];
                for (int k = 0; k < count[i]; k++)
                {
                    int j = idx[k];
                    if (j < i)
                        x[i] -= its[k] * x[j];
                }
                x[i] /= this[i][i];
            }

            return x;
        }
        public static SparseMatrix Identity(int m, int n)
        {
            int o = Math.Min(m, n);
            var A = new SparseMatrix(m, n);
            for (int i = 0; i < o; i++)
                A[i, i] = 1.0;

            return A;
        }
        public bool IsLowerTriangular()
        {
            for (int i = 0; i < m; i++)
                if (indices[i].Last() > i)
                    return false;

            return true;
        }	
        private void CheckMatrixDimensions(Matrix B)
        {
            if (B.RowDimension != m || B.ColumnDimension != n)
            {
                throw new System.ArgumentException("Sparse matrix dimensions must agree.");
            }
        }
    }
}
