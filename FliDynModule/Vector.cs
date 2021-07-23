using System;
using System.Text;

namespace FliDynModule
{
    public struct Vector
    {
        private double[] v;
        public int Length
        {
            get { return v == null ? 0 : v.Length; }
        }
        public Vector(params double[] elts)
        {
            if (elts == null) throw new ArgumentNullException("elts");
            v = elts;
        }
        public static Vector Zeros(int n)
        {
            Vector v = new Vector();
            v.v = new double[n];
            return v;
        }
        public Vector Clone()
        {
            return v == null ? new Vector() : new Vector((double[])v.Clone());
        }
        public double[] ToArray()
        {
            return (double[])v.Clone();
        }
        public static void Copy(Vector src, Vector dst)
        {
            if (src.v == null)
                throw new ArgumentNullException("src");
            if (dst.v == null)
                throw new ArgumentNullException("dst");
            int n = src.v.Length;
            if (dst.v.Length != n)
                dst.v = new double[n];
            Array.Copy(src.v, dst.v, n);
        }
        public double LInfinityNorm
        {
            get
            {
                double max = 0;

                for (int i = 0; i < v.Length; i++)
                {
                    if (Math.Abs(v[i]) > max)
                        max = v[i];
                }

                return System.Math.Abs(max);
            }
        }
        public double EuclideanNorm
        {
            get
            {
                double lsq = 0;

                for (int i = 0; i < v.Length; i++)
                {
                    lsq += v[i] * v[i];

                }

                return Math.Sqrt(lsq);
            }
        }
        public double Sum
        {
            get
            {
                double sum = 0;

                for (int i = 0; i < v.Length; i++)
                {
                    sum += v[i];
                }

                return sum;
            }
        }
        public static double GetEuclideanNorm(Vector v1, Vector v2)
        {
            double[] av1 = v1.v;
            double[] av2 = v2.v;
            if (av1 == null)
                throw new ArgumentNullException("v1");
            if (av2 == null)
                throw new ArgumentNullException("v2");
            if (av1.Length != av2.Length)
                throw new ArgumentException("Vector lenghtes do not match");
            double norm = 0;
            for (int i = 0; i < av1.Length; i++)
                norm += (av1[i] - av2[i]) * (av1[i] - av2[i]);
            return Math.Sqrt(norm);
        }
        public static double GetLInfinityNorm(Vector v1, Vector v2)
        {
            double[] av1 = v1.v;
            double[] av2 = v2.v;
            if (av1 == null)
                throw new ArgumentNullException("v1");
            if (av2 == null)
                throw new ArgumentNullException("v2");
            if (av1.Length != av2.Length)
                throw new ArgumentException("Vector lenghtes do not match");
            double norm = 0;
            for (int i = 0; i < av1.Length; i++)
                norm = Math.Max(norm, Math.Abs(av1[i] - av2[i]));
            return norm;
        }
        public static Vector Lerp(double t, double t0, Vector v0, double t1, Vector v1)
        {
            return (v0 * (t1 - t) + v1 * (t - t0)) / (t1 - t0);
        }
        public double this[int idx]
        {
            get { return v[idx]; }
            set { v[idx] = value; }
        }
        public static implicit operator double[] (Vector v)
        {
            return v.v;
        }
        public static implicit operator double (Vector v)
        {
            double[] av = v;
            if (av == null)
                throw new ArgumentNullException("v");
            if (av.Length != 1)
                throw new InvalidOperationException("Cannot convert multi-element vector to scalar");
            return av[0];
        }
        public static implicit operator Vector(double[] v)
        {
            return new Vector(v);
        }
        public static implicit operator Vector(double v)
        {
            return new Vector(v);
        }
        public void MulAdd(Vector v1, double factor)
        {
            double[] av1 = v1.v;
            if (av1 == null)
                throw new ArgumentNullException("v1");
            if (this.Length != av1.Length)
                throw new InvalidOperationException("Cannot add vectors of different length");

            for (int i = 0; i < v.Length; i++)
                v[i] = v[i] + factor * av1[i];
        }
        public static Vector operator +(Vector v1, Vector v2)
        {
            double[] av1 = v1;
            double[] av2 = v2;
            if (av1.Length != av2.Length)
                throw new InvalidOperationException("Cannot add vectors of different length");
            double[] result = new double[av1.Length];
            for (int i = 0; i < av1.Length; i++)
                result[i] = av1[i] + av2[i];
            return new Vector(result);
        }
        public static Vector operator +(Vector v, double c)
        {
            double[] av = v;
            double[] result = new double[av.Length];
            for (int i = 0; i < av.Length; i++)
                result[i] = av[i] + c;
            return new Vector(result);
        }
        public static Vector operator -(Vector v1, Vector v2)
        {
            double[] av1 = v1;
            double[] av2 = v2;
            if (av1.Length != av2.Length)
                throw new InvalidOperationException("Cannot subtract vectors of different length");
            double[] result = new double[av1.Length];
            for (int i = 0; i < av1.Length; i++)
                result[i] = av1[i] - av2[i];
            return new Vector(result);
        }
        public static Vector operator *(Matrix m, Vector v)
        {
            if (m == null)
                throw new ArgumentNullException("m");
            if (v.Length != m.ColumnDimension)
                throw new ArgumentException("Dimensions of vector and matrix do not match");

            double[] av = v;
            int rowDimension = m.RowDimension;
            int columnDimension = m.ColumnDimension;

            double[] result = new double[rowDimension];

            for (int i = 0; i < rowDimension; i++)
            {
                var acc = 0.0;
                var column = m[i];

                for (int j = 0; j < columnDimension; j++)
                {
                    acc += column[j] * av[j];
                }
                result[i] = acc;
            }
            return new Vector(result);
        }
        public static Vector operator *(Vector v, Matrix m)
        {
            double[] av = v;
            if (m == null)
                throw new ArgumentNullException("m");
            if (v.Length != m.RowDimension)
                throw new ArgumentException("Dimensions of matrix and vector do not match");
            double[] result = new double[m.ColumnDimension];
            for (int i = 0; i < m.RowDimension; i++)
            {

                for (int j = 0; j < m.ColumnDimension; j++)
                {
                    result[j] = result[j] + av[i] * m[i, j];
                }
            }
            return new Vector(result);
        }
        public static Vector operator *(Vector v, double a)
        {
            double[] av = v;
            double[] result = new double[av.Length];
            for (int i = 0; i < av.Length; i++)
                result[i] = a * av[i];
            return new Vector(result);
        }
        public static Vector operator *(double a, Vector v)
        {
            double[] av = v;
            double[] result = new double[av.Length];
            for (int i = 0; i < av.Length; i++)
                result[i] = a * av[i];
            return new Vector(result);
        }
        public static double operator *(Vector a, Vector b)
        {
            double res = 0;
            if (a.Length != a.Length)
                throw new InvalidOperationException("Cannot multiply vectors of different length");

            for (int i = 0; i < a.Length; i++)
            {
                res = res + a[i] * b[i];
            }
            return res;
        }
        public static Matrix operator &(Vector a, Vector b)
        {
            int m = a.Length, n = b.Length;
            Matrix res = new Matrix(m, n);
            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    res[i, j] = a[i] * b[j];
                }
            }
            return res;
        }
        public static Vector operator /(Vector v, double a)
        {
            double[] av = v;

            if (a == 0.0d) throw new DivideByZeroException("Cannot divide by zero");

            double[] result = new double[av.Length];
            for (int i = 0; i < av.Length; i++)
                result[i] = av[i] / a;
            return new Vector(result);
        }
        public static Vector operator /(Vector a, Vector b)
        {
            if (a.Length != b.Length)
                throw new InvalidOperationException("Cannot element-wise divide vectors of different length");
            double[] res = Vector.Zeros(a.Length);

            for (int i = 0; i < a.Length; i++)
            {
                if (b[i] == 0.0d) throw new DivideByZeroException("Cannot divide by zero");
                res[i] = a[i] / b[i];
            }

            return res;
        }
        public static Vector Max(Vector v1, Vector v2)
        {
            double[] av1 = v1.v;
            double[] av2 = v2.v;

            if (av1 == null)
                throw new ArgumentNullException("v1");
            if (av2 == null)
                throw new ArgumentNullException("v2");

            if (av1.Length != av2.Length)
                throw new ArgumentException("Vector lengths do not match");
            Vector y = Vector.Zeros(av1.Length);
            for (int i = 0; i < av1.Length; i++)
                y[i] = Math.Max(av1[i], av2[i]);

            return y;
        }
        public Vector Abs()
        {
            if (v == null)
                return new Vector();

            int n = v.Length;
            Vector y = Vector.Zeros(n);
            for (int i = 0; i < n; i++)
                y[i] = Math.Abs(v[i]);
            return y;
        }
        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();
            sb.Append("[");
            if (v != null)
                for (int i = 0; i < v.Length; i++)
                {
                    if (i > 0) sb.Append(", ");
                    sb.Append(v[i]);
                }
            sb.Append("]");
            return sb.ToString();
        }

        public override bool Equals(object obj)
        {
            if (obj is Vector)
            {
                var v2 = (Vector)obj;
                if (v2.Length != Length)
                    return false;
                var av2 = v2.v;
                for (var i = 0; i < v.Length; i++)
                    if (v[i] != av2[i])
                        return false;
                return true;
            }
            else
                return false;
        }

        public override int GetHashCode()
        {
            return v == null ? base.GetHashCode() : v.GetHashCode();
        }
    }
}
