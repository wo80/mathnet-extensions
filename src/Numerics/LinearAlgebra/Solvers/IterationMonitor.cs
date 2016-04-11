
namespace MathNet.Numerics.LinearAlgebra.Solvers
{
    using System;
    using System.Collections.Generic;

    public class IterationMonitor<T> : IIterationStopCriterion<T>
        where T : struct, IEquatable<T>, IFormattable
    {
        // Iteration progress reporter.
        Action<int, double> _progress;

        /// <summary>
        /// The iteration number of the last iteration.
        /// </summary>
        int _lastIteration = -1;

        /// <summary>
        /// The iteration number before the solver stopped.
        /// </summary>
        public int LastIteration
        {
            get { return _lastIteration; }
        }

        /// <summary>
        /// The array that holds the tracking information.
        /// </summary>
        List<double> _residualNormHistory;

        /// <summary>
        /// Gets the residual norm history.
        /// </summary>
        public List<double> ResidualNormHistory
        {
            get { return _residualNormHistory; }
        }

        public IterationMonitor()
        {
            _residualNormHistory = new List<double>();
        }

        public IterationMonitor(Action<int, double> progress)
        {
            _residualNormHistory = new List<double>();
            _progress = progress;
        }

        public IterationStatus DetermineStatus(int iterationNumber, Vector<T> solutionVector, Vector<T> sourceVector, Vector<T> residualVector)
        {
            if (iterationNumber < 0)
            {
                throw new ArgumentOutOfRangeException("iterationNumber");
            }

            if (iterationNumber > _lastIteration)
            {
                double residual = residualVector.InfinityNorm();

                if (_progress != null)
                {
                    _progress(iterationNumber, residual);
                }

                _residualNormHistory.Add(residual);

                _lastIteration = iterationNumber;
            }

            return IterationStatus.Continue;
        }

        public IIterationStopCriterion<T> Clone()
        {
            return new IterationMonitor<T>();
        }

        public void Reset()
        {
            _lastIteration = -1;
            _residualNormHistory.Clear();
        }

        public IterationStatus Status
        {
            get { return IterationStatus.Continue; }
        }
    }
}
