
namespace QuadratureTest.Tests
{
    using System;
    using System.Collections.Generic;

    /// <summary>
    /// From http://people.sc.fsu.edu/~jburkardt/c_src/test_int/test_int.html
    /// </summary>
    public class EngelsTestSet : ITestSet
    {
        public string Name
        {
            get { return "Engels"; }
        }

        public override string ToString()
        {
            return Name;
        }

        public IEnumerable<ITestFunction> GetTestFunctions()
        {
            for (int i = 1; i <= 20; i++)
            {
                yield return new FunctionK(i);
            }

            yield return new Function21();
            yield return new Function22();
            yield return new Function23();
            yield return new Function24();
            yield return new Function25();
            yield return new Function26();
            yield return new Function27();
            yield return new Function28();
            yield return new Function29();
            yield return new Function30();
            yield return new Function31();
            yield return new Function32();
            yield return new Function33();
            yield return new Function34();
            yield return new Function35();
            yield return new Function36();
            yield return new Function37();
            yield return new Function38();
            yield return new Function39();
            yield return new Function40();
            yield return new Function41();
            yield return new Function42();
            yield return new Function43();
            yield return new Function44();
            yield return new Function45();
            yield return new Function46();
            yield return new Function47();
            yield return new Function48();
            yield return new Function49();
            yield return new Function50();
            yield return new Function51();
            yield return new Function52();
            yield return new Function53();
            yield return new Function54();
            yield return new Function55();
            yield return new Function56();
            yield return new Function57();
            yield return new Function58();
            yield return new Function59();
            yield return new Function60();
            yield return new Function61();
            yield return new Function62();
            yield return new Function63();
            yield return new Function64();
            yield return new Function65();
            yield return new Function66();
            yield return new Function67();
            yield return new Function68();
            yield return new Function69();
            yield return new Function70();
        }

        #region Functions 1-29

        /// <summary>
        /// F_k(x) = sum(j = 0..k-1) (-1)^j (j+1) x^(k-j-1), k = 1..n
        /// </summary>
        class FunctionK : TestFunction
        {
            int k;

            static double[] values =
            {
                1.0,
                -3.0 / 2.0,
                7.0 / 3.0,
                -35.0 / 12.0,
                37.0 / 10.0,
                -4.316666666666667, // =-259.0 / 60.0,
                 5.076190476190476, // = 533.0 / 105.0,
                -5.710714285714285, // =-1599.0 / 280.0,
                 6.456349206349207, // = 1627.0 / 252.0,
                -7.101984126984127, // =-17897.0 / 2520.0,
                 7.838528138528138, // = 18107.0 / 2310.0,
                -8.49173881673882, // =-235391.0 / 27720.0,
                 9.22187257187257, // = 237371.0 / 25740.0,
                -9.88057775557776, // != -1237371.0 / 24024.0
                10.6059496059496, // = 95549.0 / 9009.0,
               -11.26882145632146, // =-1624333.0 / 144144.0,
                11.99051683610507, // = 1632341.0 / 136136.0,
               -12.65665666033313, // !=-31014479.0 / 24504485.0
                13.37542806350856, // = 155685007.0 / 11639628.0,
               -14.04419946668399, // =-155685007.0 / 11085360.0
            };

            public FunctionK(int k)
                : base(values[k - 1])
            {
                this.k = k;
            }

            public override double Eval(double x)
            {
                count++;

                double p, y = 0.0;

                for (int j = 0; j < k; j++)
                {
                    p = (j % 2 == 0) ? 1 : -1;

                    y += p * (j + 1) * Math.Pow(x, k - j - 1);
                }

                return y;
            }
        }

        class Function21 : TestFunction
        {
            public Function21()
                : base(1085.252666666667)
            {
                // != 1627279.0 / 1500.0
            }

            public override double Eval(double x)
            {
                count++;

                return (10 * x - 1.0) * (10 * x - 1.1) * (10 * x - 1.2) * (10 * x - 1.3);
            }
        }

        class Function22 : TestFunction
        {
            public Function22()
                : base(0.01, 1.1, 4999.586776859504)
            {
            }

            public override double Eval(double x)
            {
                count++;

                return 1.0 / (x * x * x);
            }
        }

        class Function23 : TestFunction
        {
            public Function23()
                : base(0.01, 1.1, 333333.0828950663)
            {
            }

            public override double Eval(double x)
            {
                count++;

                return 1.0 / (x * x * x * x);
            }
        }

        class Function24 : TestFunction
        {
            public Function24()
                : base(0.01, 1.1, 25000000 - 1.0 / 5.8564)
            {
                // != 610.1808091361587
            }

            public override double Eval(double x)
            {
                count++;

                return 1.0 / (x * x * x * x * x);
            }
        }

        class Function25 : TestFunction
        {
            public Function25()
                : base(0.6931471805599453)
            {
                // ln2
            }

            public override double Eval(double x)
            {
                count++;

                return 1.0 / (1.0 + x);
            }
        }

        class Function26 : TestFunction
        {
            public Function26()
                : base(0.7853981633974483)
            {
                // pi / 4
            }

            public override double Eval(double x)
            {
                count++;

                return 1.0 / (1.0 + x * x);
            }
        }

        class Function27 : TestFunction
        {
            public Function27()
                : base(0.866972987339911)
            {
                // (ln((2 + sqrt(2)) / (2 - sqrt(2))) + pi) / sqrt(32)
            }

            public override double Eval(double x)
            {
                count++;

                return 1.0 / (1.0 + x * x * x * x);
            }
        }

        class Function28 : TestFunction
        {
            public Function28()
                : base(-1.0, 1.0, 1.456205826451164)
            {
                // ln(9) / 4 + pi / sqrt(12)
            }

            public override double Eval(double x)
            {
                count++;

                return 1.0 / (1.0 + x * x + x * x * x * x);
            }
        }

        class Function29 : TestFunction
        {
            public Function29()
                : base(-1.0, 1.0, 1.56439644406905)
            {
                // 2 / sqrt(1.005) * atan(1.0 / sqrt(1.005))
            }

            public override double Eval(double x)
            {
                count++;

                return 1.0 / (1.005 + x * x);
            }
        }

        #endregion

        #region Functions 30-39

        class Function30 : TestFunction
        {
            public Function30()
                : base(0.0, 10.0, 0.4993633810764567)
            {
                // atan(500 / pi)
            }

            public override double Eval(double x)
            {
                count++;

                return 50  / ((2500 * x * x + 1) * Math.PI);
            }
        }

        class Function31 : TestFunction
        {
            const int k = 2;

            static double[] values =
            {
                29.422553486074, // k=2
                97.3465489, // k=3
                312.159332, // k=4
                991.458833, // k=5
                3139.5926 // k=6
            };

            public Function31()
                : base(-1.0, 1.0, values[k - 2])
            {
            }

            public override double Eval(double x)
            {
                count++;

                return 1 / (x * x + Math.Pow(10, -k));
            }
        }

        class Function32 : TestFunction
        {
            public Function32()
                : base(0.5144128009905458)
            {
            }

            public override double Eval(double x)
            {
                count++;

                return 1 / (1 + 5 * x * x);
            }
        }

        class Function33 : TestFunction
        {
            public Function33()
                : base(0.3998760050557662)
            {
            }

            public override double Eval(double x)
            {
                count++;

                return 1 / (1 + 10 * x * x);
            }
        }

        class Function34 : TestFunction
        {
            public Function34()
                : base(1.246450480280461)
            {
            }

            public override double Eval(double x)
            {
                count++;

                return 1 / (1 - 0.5 * x * x);
            }
        }

        class Function35 : TestFunction
        {
            public Function35()
                : base(2.670965314886704)
            {
            }

            public override double Eval(double x)
            {
                count++;

                return 1 / (1 - 0.98 * x * x);
            }
        }

        class Function36 : TestFunction
        {
            public Function36()
                : base(3.803756514651015)
            {
            }

            public override double Eval(double x)
            {
                count++;

                return 1 / (1 - 0.998 * x * x);
            }
        }

        class Function37 : TestFunction
        {
            const double r = 2;
            const double alpha = 0.5;

            public Function37()
                : base(Math.PI)
            {
            }

            public override double Eval(double x)
            {
                count++;

                return Math.Pow(2, r) / (1 + Math.Pow(2, r) * (x - alpha) * (x - alpha));
            }
        }

        class Function38 : TestFunction
        {
            public Function38()
                : base(0.4)
            {
            }

            public override double Eval(double x)
            {
                count++;

                return Math.Pow(x, 3.0 / 2.0);
            }
        }

        class Function39 : TestFunction
        {
            public Function39()
                : base(2.0 / 7.0)
            {
            }

            public override double Eval(double x)
            {
                count++;

                return Math.Pow(x, 5.0 / 2.0);
            }
        }

        #endregion

        #region Functions 40-49

        class Function40 : TestFunction
        {
            const int k = 1;

            static double[] values =
            {
                0.666666666666666, // k=1
                0.8, // k=2
                0.8888888888888899, // k=3
                0.941176470588235 // k=4
            };

            public Function40()
                : base(values[k - 1])
            {
            }

            public override double Eval(double x)
            {
                count++;

                return Math.Pow(x, 1 / (2.0 * k));
            }
        }

        class Function41 : TestFunction
        {
            const int k = 3;

            static double[] values =
            {
                0.4647525054, // k=1
                0.0, // k=2
                0.1488716212, // k=3
                0.0, // k=4
                0.06551476837 // k=5
            };

            public Function41()
                : base(values[k - 1])
            {
            }

            public override double Eval(double x)
            {
                count++;

                return Math.Sqrt(Math.Pow(Math.Abs(x * x - 0.25), k));
            }
        }

        class Function42 : TestFunction
        {
            public Function42()
                : base(-9.0, 100.0, 26.0)
            {
            }

            public override double Eval(double x)
            {
                count++;

                return 1 / Math.Sqrt(Math.Abs(x));
            }
        }

        class Function43 : TestFunction
        {
            public Function43()
                : base(0.833333333333333)
            {
            }

            public override double Eval(double x)
            {
                count++;

                return Math.Pow(x, 1.0 / 5.0);
            }
        }

        class Function44 : TestFunction
        {
            public Function44()
                : base(0.909090909090909)
            {
            }

            public override double Eval(double x)
            {
                count++;

                return Math.Pow(x, 1.0 / 10.0);
            }
        }

        class Function45 : TestFunction
        {
            public Function45()
                : base(-1.0, 1.0, 1.460447131787105)
            {
            }

            public override double Eval(double x)
            {
                count++;

                return Math.Sqrt(Math.Abs(x + 0.5));
            }
        }

        class Function46 : TestFunction
        {
            const double alpha = 0.1;
            const double beta = 0.5;

            public Function46()
                : base((Math.Pow(beta, alpha + 1) + Math.Pow(1 - beta, alpha + 1)) / (alpha + 1))
            {
            }

            public override double Eval(double x)
            {
                count++;

                return Math.Pow(Math.Abs(x - beta), alpha);
            }
        }

        class Function47 : TestFunction
        {
            const double alpha = -0.55;
            const double beta = 0.5;

            public Function47()
                : base(Math.Pow(beta, alpha + 1) / Math.Pow(alpha + 1, 2) * ((alpha + 1) * Math.Log(beta) - 1)
                + Math.Pow(1 - beta, alpha + 1) / Math.Pow(alpha + 1, 2) * ((alpha + 1) * Math.Log(1 - beta) - 1))
            {
            }

            public override double Eval(double x)
            {
                count++;

                return Math.Pow(Math.Abs(x - beta), alpha) * Math.Log(Math.Abs(x - beta));
            }
        }

        class Function48 : TestFunction
        {
            public Function48()
                : base(1.718281828459045)
            {
            }

            public override double Eval(double x)
            {
                count++;

                return Math.Exp(x);
            }
        }

        class Function49 : TestFunction
        {
            public Function49()
                : base(0.3798854930417224)
            {
            }

            public override double Eval(double x)
            {
                count++;

                return 1 / (Math.Exp(x) + 1);
            }
        }

        #endregion

        #region Functions 50-59

        class Function50 : TestFunction
        {
            public Function50()
                : base(0.1705573495024382)
            {
            }

            public override double Eval(double x)
            {
                count++;

                return x / (Math.Exp(x) + 1);
            }
        }

        class Function51 : TestFunction
        {
            public Function51()
                : base(0.0, 10.0, 0.5000002111661001)
            {
            }

            public override double Eval(double x)
            {
                count++;

                return Math.Sqrt(50) * Math.Exp(-50 * 3.14159 * x * x);
            }
        }

        class Function52 : TestFunction
        {
            public Function52()
                : base(0.0, 10.0, 1.0)
            {
            }

            public override double Eval(double x)
            {
                count++;

                return 25 * Math.Exp(-25 * x);
            }
        }

        class Function53 : TestFunction
        {
            public Function53()
                : base(0.5321250988215339)
            {
            }

            public override double Eval(double x)
            {
                count++;

                return Math.Exp(-x) - Math.Exp(-10 * x);
            }
        }

        class Function54 : TestFunction
        {
            static double c = 2 / Math.Sqrt(Math.PI);

            public Function54()
                : base(0.3958259698343338)
            {
            }

            public override double Eval(double x)
            {
                count++;

                return c * (Math.Exp(-9 * x * x) + Math.Exp(-1024 * (x - 0.25) * (x - 0.25)));
            }
        }

        class Function55 : TestFunction
        {
            public Function55()
                : base(0.0, 10.0, 14.02585092994046)
            {
            }

            public override double Eval(double x)
            {
                count++;

                return Math.Log(x);
            }
        }

        class Function56 : TestFunction
        {
            public Function56()
                : base(0.0001, 7.0, -1.967546385989968)
            {
            }

            public override double Eval(double x)
            {
                count++;

                return Math.Log(x) * Math.Sin(x);
            }
        }

        class Function57 : TestFunction
        {
            public Function57()
                : base(1.711857371268652)
            {
            }

            public override double Eval(double x)
            {
                count++;

                return Math.Pow(2, x) * x / (Math.Pow(2, x) - 1);
            }
        }

        class Function58 : TestFunction
        {
            public Function58()
                : base(0.6366197723675814)
            {
            }

            public override double Eval(double x)
            {
                count++;

                return Math.Sin(Math.PI * x);
            }
        }

        class Function59 : TestFunction
        {
            public Function59()
                : base(0.841470984807897)
            {
            }

            public override double Eval(double x)
            {
                count++;

                return Math.Cos(x);
            }
        }

        #endregion

        #region Functions 60-69

        class Function60 : TestFunction
        {
            public Function60()
                : base(1.154700538379252)
            {
            }

            public override double Eval(double x)
            {
                count++;

                return 2 / (2 + Math.Sin(10 * Math.PI * x));
            }
        }

        class Function61 : TestFunction
        {
            public Function61()
                : base(0.0, 2 * Math.PI, -0.2096724796611653)
            {
            }

            public override double Eval(double x)
            {
                count++;

                return x * Math.Sin(30 * x) * Math.Cos(x);
            }
        }

        class Function62 : TestFunction
        {
            public Function62()
                : base(0.0, 2 * Math.PI, 0.1178097245096172)
            {
            }

            public override double Eval(double x)
            {
                count++;

                return x * Math.Sin(30 * x) * Math.Cos(50 * x);
            }
        }

        class Function63 : TestFunction
        {
            public Function63()
                : base(0.0, 2 * Math.PI, -1.27162980944677)
            {
            }

            public override double Eval(double x)
            {
                count++;

                return Math.PI * x * Math.Sin(30 * x) / Math.Sqrt(4 * Math.PI * Math.PI - x * x);
            }
        }

        class Function64 : TestFunction
        {
            public Function64()
                : base(-1.0, 1.0, 0.4794282266888015)
            {
            }

            public override double Eval(double x)
            {
                count++;

                return 0.92 * Math.Cosh(x) - Math.Cos(x);
            }
        }

        class Function65 : TestFunction
        {
            public Function65()
                : base(0.0)
            {
            }

            public override double Eval(double x)
            {
                count++;

                return Math.Sin(100 * Math.PI * x);
            }
        }

        class Function66 : TestFunction
        {
            public Function66()
                : base(1.0 / 20.0, 1.0 / 3.0, -0.3004108269560286)
            {
            }

            public override double Eval(double x)
            {
                count++;

                return 1 / x * Math.Sin(1 / x);
            }
        }

        class Function67 : TestFunction
        {
            const int m = 1;

            public Function67()
                : base(1.0)
            {
            }

            public override double Eval(double x)
            {
                count++;

                return Math.Tan(Math.PI * x / 2) * Math.Sin(m * Math.PI * x);
            }
        }

        class Function68 : TestFunction
        {
            const int m = 2;

            public Function68()
                : base(-1.0)
            {
            }

            public override double Eval(double x)
            {
                count++;

                return Math.Tan(Math.PI * x / 2) * Math.Sin(m * Math.PI * x);
            }
        }

        class Function69 : TestFunction
        {
            const int m = 3;

            public Function69()
                : base(1.0)
            {
            }

            public override double Eval(double x)
            {
                count++;

                return Math.Tan(Math.PI * x / 2) * Math.Sin(m * Math.PI * x);
            }
        }

        #endregion

        #region Functions 70-79

        class Function70 : TestFunction
        {
            const int m = 4;

            public Function70()
                : base(-1.0)
            {
            }

            public override double Eval(double x)
            {
                count++;

                return Math.Tan(Math.PI * x / 2) * Math.Sin(m * Math.PI * x);
            }
        }

        #endregion

        #region Functions 80-89

        class Function80 : TestFunction
        {
            public Function80()
                : base(0.0, 0.0, 0.0)
            {
            }

            public override double Eval(double x)
            {
                count++;

                return 0.0;
            }
        }

        #endregion

        #region Functions 90-96

        class Function90 : TestFunction
        {
            public Function90()
                : base(0.0, 0.0, 0.0)
            {
            }

            public override double Eval(double x)
            {
                count++;

                return 0.0;
            }
        }

        #endregion
    }
}
