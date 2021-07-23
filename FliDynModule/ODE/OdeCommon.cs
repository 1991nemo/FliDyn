using System;

namespace FliDynModule.ODE
{
    public class Options
    {
        public double InitialStep { get; set; }
        public double AbsoluteTolerance { get; set; }
        public double RelativeTolerance { get; set; }
        public double OutputStep { get; set; }
        public double MaxStep { get; set; }
        public double MinStep { get; set; }
        public double MaxScale { get; set; }
        public double MinScale { get; set; }
        public int NumberOfIterations { get; set; }
        public Matrix Jacobian { get; set; }
        public SparseMatrix SparseJacobian { get; set; }
        public Options()
        {
            InitialStep = 0.0d;
            AbsoluteTolerance = 1e-6;
            RelativeTolerance = 1e-3;
            MaxStep = Double.MaxValue;
            MinStep = 0.0d;
            MaxScale = 1.1d;
            MinScale = 0.9d;
            OutputStep = 0.0d;
            NumberOfIterations = 5;
        }

        private static readonly Options defaultOpts = new Options();
        public static Options Default
        {
            get { return defaultOpts; }
        }
    }
    public struct SolPoint
    {
        private Vector x; //Problem's phase variables
        private double t; //Current time
        public Vector X
        {
            get { return x; }
        }
        
        public double T
        {
            get { return t; }
        }
        internal SolPoint(double t, Vector x)
        {
            this.x = x;
            this.t = t;
        }
    }
}
