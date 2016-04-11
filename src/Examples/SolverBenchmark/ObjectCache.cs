using MathNet.Numerics.LinearAlgebra.Solvers;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace SolverBenchmark
{
    /// <summary>
    /// Cache solver and preconditioner objects.
    /// </summary>
    class ObjectCache
    {
        private Dictionary<Type, IIterativeSolver<double>>  _doubleSolvers;
        private Dictionary<Type, IIterativeSolver<Complex>> _complexSolvers;

        private Dictionary<Type, IPreconditioner<double>> _doublePreconditioners;
        private Dictionary<Type, IPreconditioner<Complex>> _complexPreconditioners;

        public ObjectCache()
        {
            _doubleSolvers = new Dictionary<Type, IIterativeSolver<double>>();
            _complexSolvers = new Dictionary<Type, IIterativeSolver<Complex>>();

            _doublePreconditioners = new Dictionary<Type, IPreconditioner<double>>();
            _complexPreconditioners = new Dictionary<Type, IPreconditioner<Complex>>();
        }

        public void Refresh()
        {
            _doubleSolvers.Clear();
            _complexSolvers.Clear();

            _doublePreconditioners.Clear();
            _complexPreconditioners.Clear();
        }

        public object GetSolver(string typeName, Type solverType)
        {
            if (typeName.Equals("Double", StringComparison.OrdinalIgnoreCase))
            {
                IIterativeSolver<double> solver;

                if (!_doubleSolvers.TryGetValue(solverType, out solver))
                {
                    solver = (IIterativeSolver<double>)ReflectionHelper.CreateInstance(solverType);
                    _doubleSolvers[solverType] = solver;
                }

                return solver;
            }
            else if (typeName.Equals("Complex", StringComparison.OrdinalIgnoreCase))
            {
                IIterativeSolver<Complex> solver;

                if (!_complexSolvers.TryGetValue(solverType, out solver))
                {
                    solver = (IIterativeSolver<Complex>)ReflectionHelper.CreateInstance(solverType);
                    _complexSolvers[solverType] = solver;
                }

                return solver;
            }

            throw new Exception();
        }

        public object GetPreconditioner(string typeName, Type precondType)
        {
            if (typeName.Equals("Double", StringComparison.OrdinalIgnoreCase))
            {
                IPreconditioner<double> solver;

                if (!_doublePreconditioners.TryGetValue(precondType, out solver))
                {
                    solver = (IPreconditioner<double>)ReflectionHelper.CreateInstance(precondType);
                    _doublePreconditioners[precondType] = solver;
                }

                return solver;
            }
            else if (typeName.Equals("Complex", StringComparison.OrdinalIgnoreCase))
            {
                IPreconditioner<Complex> solver;

                if (!_complexPreconditioners.TryGetValue(precondType, out solver))
                {
                    solver = (IPreconditioner<Complex>)ReflectionHelper.CreateInstance(precondType);
                    _complexPreconditioners[precondType] = solver;
                }

                return solver;
            }

            throw new Exception();
        }
    }
}
