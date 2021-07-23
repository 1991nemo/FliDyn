using System;
using System.Collections.Generic;
using System.Linq;
//using FliDynModule;

namespace FliDynModule.ODE
{
    public static class OdeHelpers
    {
        public static IEnumerable<SolPoint> SolveTo(this IEnumerable<SolPoint> solution, double to)
        {
            return solution.TakeWhile(p => p.T <= to);
        }
        public static IEnumerable<SolPoint> SolveFromTo(this IEnumerable<SolPoint> solution, double from, double to)
        {
            return solution.SkipWhile(p => p.T < from).TakeWhile(p => p.T <= to);
        }
        public static IEnumerable<SolPoint> SolveFromToStep(this IEnumerable<SolPoint> solution, double from, double to, double step)
        {
            return solution.WithStep(step).SkipWhile(p => p.T < from).TakeWhile(p => p.T <= to);
        }
        public static IEnumerable<SolPoint> WithStep(this IEnumerable<SolPoint> solution, double delta)
        {
            var en = solution.GetEnumerator();
            if (!en.MoveNext())
                yield break;
            SolPoint prev = en.Current;
            double n = Math.Ceiling(prev.T / delta);
            double tout = n * delta;
            if (tout == prev.T)
                yield return prev;
            n++;
            tout = n * delta;
            while (true)
            {
                if (!en.MoveNext())
                    yield break;
                SolPoint current = en.Current;
                while (current.T >= tout)
                {
                    yield return new SolPoint(tout, Vector.Lerp(tout, prev.T, prev.X, current.T, current.X));
                    n++;
                    tout = n * delta;
                }
                prev = current;
            }
        }
        public static IEnumerable<SolPoint> AppendStep(this IEnumerable<SolPoint> solution)
        {
            var en = solution.GetEnumerator();
            if (!en.MoveNext())
                yield break;
            SolPoint prev = en.Current;
            double[] temp = new double[prev.X.Length + 1];
            Array.Copy(prev.X, temp, prev.X.Length);
            yield return new SolPoint(prev.T, new Vector((double[])temp.Clone()));
            while (en.MoveNext())
            {
                SolPoint curr = en.Current;
                Array.Copy(curr.X, temp, curr.X.Length);
                temp[curr.X.Length] = curr.T - prev.T;
                yield return new SolPoint(curr.T, new Vector((double[])temp.Clone()));
                prev = curr;
            }
            yield break;
        }
    }
}
