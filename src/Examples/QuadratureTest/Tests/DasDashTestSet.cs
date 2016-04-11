
namespace QuadratureTest.Tests
{
    using System;
    using System.Collections.Generic;

    /// <summary>
    /// From "Application of Mixed Quadrature Rules in the Adaptive Quadrature Routine"
    /// Debasish Das, Rajani B. Dash
    /// </summary>
    public class DasDashTestSet : ITestSet
    {
        public string Name
        {
            get { return "Das-Dash"; }
        }

        public override string ToString()
        {
            return Name;
        }

        public IEnumerable<ITestFunction> GetTestFunctions()
        {
            yield return new Function1();
            yield return new Function2();
            yield return new Function3();
            yield return new Function4();
            yield return new Function5();
            yield return new Function6();
            yield return new Function7();
            yield return new Function8();
            yield return new Function9();
            yield return new Function10();
            yield return new Function11();
            yield return new Function12();
            yield return new Function13();
            yield return new Function14();
            yield return new Function15();
            yield return new Function16();
        }

        public class Function1 : TestFunction
        {
            public Function1()
                : base(0.0, 10.0 * Math.PI, -21.921477854236)
            {
            }

            public override double Eval(double x)
            {
                count++;

                return Math.Sin(x) * Math.Exp(x / 10);
            }
        }

        public class Function2 : TestFunction
        {
            public Function2()
                : base(0.0, 4.0, -1.54878837252795)
            {
            }

            public override double Eval(double x)
            {
                count++;

                return 13 * (x - Math.Pow(x, 2)) * Math.Exp(-3 * x / 2);
            }
        }

        public class Function3 : TestFunction
        {
            public Function3()
                : base(0.0, 2.0 * Math.PI, -0.209672479661)
            {
            }

            public override double Eval(double x)
            {
                count++;

                return x * Math.Sin(30 * x) * Math.Cos(x);
            }
        }

        public class Function4 : TestFunction
        {
            public Function4()
                : base(0.0, 1.0, 1.154700538379)
            {
            }

            public override double Eval(double x)
            {
                count++;

                return 2 / (2 + Math.Sin(10 * Math.PI * x));
            }
        }

        public class Function5 : TestFunction
        {
            public Function5()
                : base(0.0, 1.0, 0.04912172951763)
            {
            }

            public override double Eval(double x)
            {
                count++;

                return Math.Pow(x, 16) * Math.Cos(Math.Pow(x, 16));
            }
        }

        public class Function6 : TestFunction
        {
            public Function6()
                : base(0.0, 1.0, 0.666666666667)
            {
            }

            public override double Eval(double x)
            {
                count++;

                return Math.Sqrt(x);
            }
        }

        public class Function7 : TestFunction
        {
            public Function7()
                : base(0.0, 1.0, 0.849726325420)
            {
            }

            public override double Eval(double x)
            {
                count++;

                return Math.Sin(Math.Sqrt(Math.PI * x));
            }
        }

        public class Function8 : TestFunction
        {
            public Function8()
                : base(0.0, 2.0, 1.141592653589)
            {
            }

            public override double Eval(double x)
            {
                count++;

                return Math.Asin(Math.Sqrt(x / (2 + x)));
            }
        }

        public class Function9 : TestFunction
        {
            public Function9()
                : base(0.0, 1.0, 0.130996908077)
            {
            }

            public override double Eval(double x)
            {
                count++;

                return Math.Pow(Sech(10 * x - 0.2), 2)
                    + Math.Pow(Sech(100 * x - 0.4), 4)
                    + Math.Pow(Sech(1000 * x - 0.6), 6);
            }

            private static double Sech(double x)
            {
                return (1.0 / Math.Cosh(x));
            }
        }

        public class Function10 : TestFunction
        {
            public Function10()
                : base(0.0, 5.0, 0.49872676724581)
            {
            }

            public override double Eval(double x)
            {
                count++;

                return 50 / (Math.PI * (1 + 2500 * Math.Pow(x, 2)));
            }
        }

        public class Function11 : TestFunction
        {
            public Function11()
                : base(0.0, 2.0, -1.115957990932)
            {
            }

            public override double Eval(double x)
            {
                count++;

                return Math.Exp(x) * Math.Sin(Math.Pow(x, 2) * Math.Cos(Math.Exp(x)));
            }
        }

        public class Function12 : TestFunction
        {
            public Function12()
                : base(0.0, 1.0, -0.704379707168)
            {
            }

            public override double Eval(double x)
            {
                count++;

                return (30 * Math.Pow(x, 9) * (Math.Cos(Math.Pow(x, 6)) - 1) / (1 + Math.Pow(x, 10))) * Math.Exp(Math.Pow(x, 15));
            }
        }

        public class Function13 : TestFunction
        {
            public Function13()
                : base(0.0, 1.0, 0.866972987339)
            {
            }

            public override double Eval(double x)
            {
                count++;

                return 1 / (Math.Pow(x, 4) + 1);
            }
        }

        public class Function14 : TestFunction
        {
            public Function14()
                : base(-1.0, 1.0, 1.582232963729)
            {
            }

            public override double Eval(double x)
            {
                count++;

                return 1 / (Math.Pow(x, 4) + Math.Pow(x, 2) + 0.9);
            }
        }

        public class Function15 : TestFunction
        {
            public Function15()
                : base(0.0, 4.0, 0.966440320387)
            {
            }

            public override double Eval(double x)
            {
                count++;

                return Math.Cos(Math.Cos(x) + 3 * Math.Sin(x) + 2 * Math.Cos(2 * x) + 3 * Math.Sin(2 * x) + 3 * Math.Cos(3 * x));
            }
        }

        public class Function16 : TestFunction
        {
            public Function16()
                : base(0.0, 2 * Math.PI, 0.002514279834)
            {
            }

            public override double Eval(double x)
            {
                count++;

                return x * Math.Cos(50 * x) * Math.Sin(x);
            }
        }
    }
}
