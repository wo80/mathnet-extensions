
namespace SolverBenchmark
{
    using MathNet.Numerics.LinearAlgebra.Solvers;
    using System;
    using System.Collections.Generic;
    using System.Globalization;
    using System.Linq;
    using System.Numerics;
    using System.Reflection;
    using System.Threading.Tasks;

    public static class ReflectionHelper
    {
        static List<Type> _doubleSolvers;
        static List<Type> _complexSolvers;

        static List<Type> _doublePreconditioners;
        static List<Type> _complexPreconditioners;

        public static bool HasProperties(object obj)
        {
            return obj.GetType().GetProperties().Length > 0;
        }

        public static object CreateInstance(Type type)
        {
            if (type.GetConstructors().Any(z => z.GetParameters().Count() == 0))
            {
                return Activator.CreateInstance(type);
            }

            // This is a hack for class with constructors with one optional parameter.
            return Activator.CreateInstance(type,
                BindingFlags.CreateInstance |
                BindingFlags.Public |
                BindingFlags.Instance |
                BindingFlags.OptionalParamBinding,
                null, new object[] { Type.Missing },
                CultureInfo.CurrentCulture);
        }

        public static async Task<List<Type>> GetSolvers<T>()
            where T : struct, IEquatable<T>, IFormattable
        {
            if (typeof(T) == typeof(double))
            {
                if (_doubleSolvers == null)
                {
                    await Load();
                }

                return _doubleSolvers;
            }
            else
            {
                if (_complexSolvers == null)
                {
                    await Load();
                }

                return _complexSolvers;
            }
        }

        public static async Task<List<Type>> GetPreconditioners<T>()
            where T : struct, IEquatable<T>, IFormattable
        {
            if (typeof(T) == typeof(double))
            {
                if (_doublePreconditioners == null)
                {
                    await Load();
                }

                return _doublePreconditioners;
            }
            else
            {
                if (_complexPreconditioners == null)
                {
                    await Load();
                }

                return _complexPreconditioners;
            }
        }

        private static async Task Load()
        {
            await Task.Run(() =>
            {
                // Make sure extensions assembly gets loaded.
                var s = new MathNet.Numerics.LinearAlgebra.Double.Solvers.SSORPreconditioner();

                _doubleSolvers = LoadSolvers<double>().ToList();
                _complexSolvers = LoadSolvers<Complex>().ToList();

                _doublePreconditioners = LoadPreconditioners<double>().ToList();
                _complexPreconditioners = LoadPreconditioners<Complex>().ToList();

                // Manually add identity preconditioner.
                _doublePreconditioners.Insert(0, typeof(UnitPreconditioner<double>));
                _complexPreconditioners.Insert(0, typeof(UnitPreconditioner<Complex>));
            });
        }

        private static IEnumerable<Type> LoadSolvers<T>()
            where T : struct, IEquatable<T>, IFormattable
        {
            var type = typeof(IIterativeSolver<T>);

            return AppDomain.CurrentDomain.GetAssemblies()
                .SelectMany(s => s.GetTypes())
                .Where(p => type.IsAssignableFrom(p));
        }

        private static IEnumerable<Type> LoadPreconditioners<T>()
            where T : struct, IEquatable<T>, IFormattable
        {
            var type = typeof(IPreconditioner<T>);

            return AppDomain.CurrentDomain.GetAssemblies()
                .SelectMany(s => s.GetTypes())
                .Where(p => type.IsAssignableFrom(p));
        }
    }
}
