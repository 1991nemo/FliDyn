using System;

namespace FliDynModule
{
    public struct SparseVector
    {
        const int IncrementSize = 16;
        private int n;
        public double[] items;
        public int[] indices;
        public int count;
        public SparseVector(int n)
        {
            this.n = n;
            this.items = new double[IncrementSize];
            this.indices = new int[IncrementSize];
            this.count = 0;
        }
        public SparseVector(double[] items, int[] indices, int n)
        {
            if (items == null)
                throw new ArgumentNullException("items");
            if (indices == null)
                throw new ArgumentNullException("indices");
            this.items = items;
            this.indices = indices;
            this.count = items.Length;
            this.n = n;
        }
        public int Length
        {
            get { return n; }
        }

        public SparseVector Clone()
        {
            return n == 0 ? new SparseVector() : new SparseVector((double[])items.Clone(), (int[])indices.Clone(), n);
        }
        public double this[int i]
        {
            get
            {
                if (i < 0 || i >= n)
                    throw new IndexOutOfRangeException();
                int idx = Array.BinarySearch(indices, 0, count, i);
                if (idx < 0)
                    return 0;
                else
                    return items[idx];
            }
            set
            {
                if (i < 0 || i >= n)
                    throw new IndexOutOfRangeException();
                int idx = Array.BinarySearch(indices, 0, count, i);
                if (idx >= 0)
                    items[idx] = value;
                else
                {
                    int indexToAdd = ~idx;
                    if (count >= items.Length)
                    {
                        int delta = Math.Min(IncrementSize, n - items.Length);
                        int[] newIndices = new int[indices.Length + delta];
                        double[] newItems = new double[items.Length + delta];
                        Array.Copy(indices, newIndices, indices.Length);
                        Array.Copy(items, newItems, items.Length);
                        items = newItems;
                        indices = newIndices;

                    }
                    Array.Copy(indices, indexToAdd, indices, indexToAdd + 1, count - indexToAdd);
                    Array.Copy(items, indexToAdd, items, indexToAdd + 1, count - indexToAdd);
                    count++;
                    indices[indexToAdd] = i;
                    items[indexToAdd] = value;
                }
            }
        }
    }
}
