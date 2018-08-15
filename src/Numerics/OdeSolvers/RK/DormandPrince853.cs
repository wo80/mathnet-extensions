// Based on Fortran code DOP853
// Copyright (c) 2004, Ernst Hairer
// License: Simplified BSD License (https://www.unige.ch/~hairer/software.html)

namespace MathNet.Numerics.OdeSolvers.RK
{
    using System;

    /// <summary>
    /// Numerical solution of a system of first order ordinary differential equations y'=f(x,y).
    /// This is an explicit runge-kutta method of order 8(5,3) due to Dormand & Prince (with stepsize
    /// control and dense output)
    ///
    /// Authors: E. Hairer and G. Wanner
    ///
    /// This code is described in:
    ///         E. Hairer, S.P. Norsett and G. Wanner
    ///         Solving Ordinary Differential Equations I. Nonstiff Problems (2nd edition)
    ///         Springer-Verlag (1993)
    /// </summary>
    public class DormandPrince853 : IRungeKuttaStepper
    {
        #region Runge-Kutta coefficients

        private const double c2 = 0.0526001519587677318785587544488;
        private const double c3 = 0.0789002279381515978178381316732;
        private const double c4 = 0.11835034190722739672675719751;
        private const double c5 = 0.28164965809277260327324280249;
        private const double c6 = 0.333333333333333333333333333333;
        private const double c7 = 0.25;
        private const double c8 = 0.307692307692307692307692307692;
        private const double c9 = 0.651282051282051282051282051282;
        private const double c10 = 0.6;
        private const double c11 = 0.857142857142857142857142857142;
        private const double c14 = 0.1;
        private const double c15 = 0.2;
        private const double c16 = 0.777777777777777777777777777778;

        private const double a21 = 0.0526001519587677318785587544488;

        private const double a31 = 0.0197250569845378994544595329183;
        private const double a32 = 0.0591751709536136983633785987549;

        private const double a41 = 0.0295875854768068491816892993775;
        private const double a43 = 0.0887627564304205475450678981324;

        private const double a51 = 0.241365134159266685502369798665;
        private const double a53 = -0.884549479328286085344864962717;
        private const double a54 = 0.924834003261792003115737966543;

        private const double a61 = 0.037037037037037037037037037037;
        private const double a64 = 0.170828608729473871279604482173;
        private const double a65 = 0.125467687566822425016691814123;

        private const double a71 = 0.037109375;
        private const double a74 = 0.170252211019544039314978060272;
        private const double a75 = 0.0602165389804559606850219397283;
        private const double a76 = -0.017578125;

        private const double a81 = 0.0370920001185047927108779319836;
        private const double a84 = 0.170383925712239993810214054705;
        private const double a85 = 0.107262030446373284651809199168;
        private const double a86 = -0.0153194377486244017527936158236;
        private const double a87 = 0.00827378916381402288758473766002;

        private const double a91 = 0.624110958716075717114429577812;
        private const double a94 = -3.36089262944694129406857109825;
        private const double a95 = -0.868219346841726006818189891453;
        private const double a96 = 27.5920996994467083049415600797;
        private const double a97 = 20.1540675504778934086186788979;
        private const double a98 = -43.4898841810699588477366255144;

        private const double a101 = 0.477662536438264365890433908527;
        private const double a104 = -2.48811461997166764192642586468;
        private const double a105 = -0.590290826836842996371446475743;
        private const double a106 = 21.2300514481811942347288949897;
        private const double a107 = 15.2792336328824235832596922938;
        private const double a108 = -33.2882109689848629194453265587;
        private const double a109 = -0.0203312017085086261358222928593;

        private const double a111 = -0.93714243008598732571704021658;
        private const double a114 = 5.18637242884406370830023853209;
        private const double a115 = 1.09143734899672957818500254654;
        private const double a116 = -8.14978701074692612513997267357;
        private const double a117 = -18.5200656599969598641566180701;
        private const double a118 = 22.7394870993505042818970056734;
        private const double a119 = 2.49360555267965238987089396762;
        private const double a1110 = -3.0467644718982195003823669022;

        private const double a121 = 2.27331014751653820792359768449;
        private const double a124 = -10.5344954667372501984066689879;
        private const double a125 = -2.00087205822486249909675718444;
        private const double a126 = -17.9589318631187989172765950534;
        private const double a127 = 27.9488845294199600508499808837;
        private const double a128 = -2.85899827713502369474065508674;
        private const double a129 = -8.87285693353062954433549289258;
        private const double a1210 = 12.3605671757943030647266201528;
        private const double a122 = 0.643392746015763530355970484046;

        private const double a141 = 0.0561675022830479523392909219681;
        private const double a147 = 0.253500210216624811088794765333;
        private const double a148 = -0.246239037470802489917441475441;
        private const double a149 = -0.124191423263816360469010140626;
        private const double a1410 = 0.15329179827876569731206322685;
        private const double a1411 = 0.00820105229563468988491666602057;
        private const double a1412 = 0.00756789766054569976138603589584;
        private const double a1413 = -0.008298;

        private const double a151 = 0.0318346481635021405060768473261;
        private const double a156 = 0.0283009096723667755288322961402;
        private const double a157 = 0.0535419883074385676223797384372;
        private const double a158 = -0.0549237485713909884646569340306;
        private const double a1511 = -1.08347328697249322858509316994e-4;
        private const double a1512 = 3.82571090835658412954920192323e-4;
        private const double a1513 = -3.40465008687404560802977114492e-4;
        private const double a1514 = 0.141312443674632500278074618366;

        private const double a161 = -0.428896301583791923408573538692;
        private const double a166 = -4.69762141536116384314449447206;
        private const double a167 = 7.68342119606259904184240953878;
        private const double a168 = 4.06898981839711007970213554331;
        private const double a169 = 0.356727187455281109270669543021;
        private const double a1613 = -0.00139902416515901462129418009734;
        private const double a1614 = 2.9475147891527723389556272149;
        private const double a1615 = -9.15095847217987001081870187138;

        private const double b1 = 0.0542937341165687622380535766363;
        private const double b6 = 4.45031289275240888144113950566;
        private const double b7 = 1.89151789931450038304281599044;
        private const double b8 = -5.8012039600105847814672114227;
        private const double b9 = 0.31116436695781989440891606237;
        private const double b10 = -0.152160949662516078556178806805;
        private const double b11 = 0.201365400804030348374776537501;
        private const double b12 = 0.0447106157277725905176885569043;

        private const double bh1 = 0.244094488188976377952755905512;
        private const double bh2 = 0.733846688281611857341361741547;
        private const double bh3 = 0.0220588235294117647058823529412;

        private const double e1 = 0.01312004499419488073250102996;
        private const double e6 = -1.225156446376204440720569753;
        private const double e7 = -0.4957589496572501915214079952;
        private const double e8 = 1.664377182454986536961530415;
        private const double e9 = -0.350328848749973681688648729;
        private const double e10 = 0.3341791187130174790297318841;
        private const double e11 = 0.08192320648511571246570742613;
        private const double e12 = -0.02235530786388629525884427845;

        private const double d41 = -8.4289382761090128651353491142;
        private const double d46 = 0.5667149535193777696253178359;
        private const double d47 = -3.0689499459498916912797304727;
        private const double d48 = 2.384667656512069828772814968;
        private const double d49 = 2.1170345824450282767155149946;
        private const double d410 = -0.8713915837779729920678990749;
        private const double d411 = 2.240437430260788275854177165;
        private const double d412 = 0.6315787787694688181557024929;

        private const double d413 = -0.0889903364513333108206981174;
        private const double d414 = 18.148505520854727256656404962;
        private const double d415 = -9.1946323924783554000451984436;
        private const double d416 = -4.4360363875948939664310572;

        private const double d51 = 10.427508642579134603413151009;
        private const double d56 = 242.28349177525818288430175319;
        private const double d57 = 165.20045171727028198505394887;
        private const double d58 = -374.54675472269020279518312152;
        private const double d59 = -22.113666853125306036270938578;
        private const double d510 = 7.7334326684722638389603898808;
        private const double d511 = -30.674084731089398182061213626;
        private const double d512 = -9.3321305264302278729567221706;

        private const double d513 = 15.697238121770843886131091075;
        private const double d514 = -31.139403219565177677282850411;
        private const double d515 = -9.3529243588444783865713862664;
        private const double d516 = 35.81684148639408375246589854;

        private const double d61 = 19.985053242002433820987653617;
        private const double d66 = -387.03730874935176555105901742;
        private const double d67 = -189.17813819516756882830838328;
        private const double d68 = 527.80815920542364900561016686;
        private const double d69 = -11.573902539959630126141871134;
        private const double d610 = 6.8812326946963000169666922661;
        private const double d611 = -1.000605096691083840318386098;
        private const double d612 = 0.7777137798053443209286926574;

        private const double d613 = -2.7782057523535084065932004339;
        private const double d614 = -60.196695231264120758267380846;
        private const double d615 = 84.320405506677161018159903784;
        private const double d616 = 11.99229113618278932803513003;

        private const double d71 = -25.693933462703749003312586129;
        private const double d76 = -154.18974869023643374053993627;
        private const double d77 = -231.52937917604549567536039109;
        private const double d78 = 357.6391179106141237828534991;
        private const double d79 = 93.405324183624310003907691704;
        private const double d710 = -37.458323136451633156875139351;
        private const double d711 = 104.09964950896230045147246184;
        private const double d712 = 29.84029342666050312334436357;

        private const double d713 = -43.533456590011143754432175058;
        private const double d714 = 96.3245539591882829483949506;
        private const double d715 = -39.177261675615439165231486172;
        private const double d716 = -149.72683625798562581422125276;

        #endregion

        IErrorController controller;
        StiffnessChecker stiff;

        int n;

        Action<double, double[], double[]> fcn;

        double told, dtold;
        double[] dxdt, xout, xtemp, xerr, xerr2;
        double[] k2, k3, k4, k5, k6, k7, k8, k9, k10;
        double[] r;

        double rtol, atol;

        public int Order { get; private set; }

        /// <summary>
        /// Initializes a new instance of the <see cref="DormandPrince853"/> class.
        /// </summary>
        /// <param name="n">Dimension of the system.</param>
        /// <param name="controller">The error controller.</param>
        /// <param name="stiff">The stiffness detector.</param>
        public DormandPrince853(int n, Action<double, double[], double[]> fcn, double rtol, double atol,
            IErrorController controller, StiffnessChecker stiff)
        {
            this.n = n;
            this.fcn = fcn;
            this.rtol = rtol;
            this.atol = atol;

            this.controller = controller;
            this.stiff = stiff;

            this.Order = 8;

            xout = new double[n];
            xtemp = new double[n];
            dxdt = new double[n];
            xerr = new double[n];
            xerr2 = new double[n];

            k2 = new double[n];
            k3 = new double[n];
            k4 = new double[n];
            k5 = new double[n];
            k6 = new double[n];
            k7 = new double[n];
            k8 = new double[n];
            k9 = new double[n];
            k10 = new double[n];

            r = new double[8 * n];
        }

        double sign(double a, double b)
        {
            double x = Math.Abs(a);
            return b >= 0 ? x : -x;
        }

        // Core integrator for DOP853
        public double Integrate(ref double t, ref double dt, double[] x, ref int step, int nmax, double posneg, bool dense)
        {
            double th = t, err = 0;

            // Basic integration step
            while (step < nmax)
            {
                Step(t, dt, x);

                step++;

                th = t + dt;

                // Error estimation
                err = Error(dt, x);

                dtold = dt;
                told = t;

                // Computation of HNEW
                if (controller.Success(err, posneg, ref dt))
                {
                    break;
                }
            }

            fcn(th, xout, k4);

            // Stiffness detection
            if (!stiff.Check(controller.Accepted, dtold, k4, k3, xout, xtemp))
            {
                throw new Exception(" The problem seems to become stiff at t = " + t);
            }

            if (dense)
            {
                PrepareInterpolation(t, dtold, x);
            }

            for (int i = 0; i < n; ++i)
            {
                dxdt[i] = k4[i];
                x[i] = xout[i];
            }

            t = th;

            return err;
        }

        private void PrepareInterpolation(double t, double dt, double[] x)
        {
            int i, n = this.n;

            for (i = 0; i < n; ++i)
            {
                double dx = xout[i] - x[i];
                double bspl = dt * dxdt[i] - dx;

                r[i] = x[i];
                r[i + n] = dx;
                r[i + n * 2] = bspl;
                r[i + n * 3] = dx - dt * k4[i] - bspl;
                r[i + n * 4] = dxdt[i] * d41 + k6[i] * d46 + k7[i] * d47 + k8[i] * d48 + k9[i] * d49 + k10[i] * d410 + k2[i] * d411 + k3[i] * d412;
                r[i + n * 5] = dxdt[i] * d51 + k6[i] * d56 + k7[i] * d57 + k8[i] * d58 + k9[i] * d59 + k10[i] * d510 + k2[i] * d511 + k3[i] * d512;
                r[i + n * 6] = dxdt[i] * d61 + k6[i] * d66 + k7[i] * d67 + k8[i] * d68 + k9[i] * d69 + k10[i] * d610 + k2[i] * d611 + k3[i] * d612;
                r[i + n * 7] = dxdt[i] * d71 + k6[i] * d76 + k7[i] * d77 + k8[i] * d78 + k9[i] * d79 + k10[i] * d710 + k2[i] * d711 + k3[i] * d712;
            }

            // The next three function evaluations

            for (i = 0; i < n; ++i)
            {
                xtemp[i] = x[i] + dt * (dxdt[i] * a141 + k7[i] * a147 + k8[i] * a148 + k9[i] * a149 + k10[i] * a1410 + k2[i] * a1411 + k3[i] * a1412 + k4[i] * a1413);
            }

            fcn(t + dt * c14, xtemp, k10);

            for (i = 0; i < n; ++i)
            {
                xtemp[i] = x[i] + dt * (dxdt[i] * a151 + k6[i] * a156 + k7[i] * a157 + k8[i] * a158 + k2[i] * a1511 + k3[i] * a1512 + k4[i] * a1513 + k10[i] * a1514);
            }

            fcn(t + dt * c15, xtemp, k2);

            for (i = 0; i < n; ++i)
            {
                xtemp[i] = x[i] + dt * (dxdt[i] * a161 + k6[i] * a166 + k7[i] * a167 + k8[i] * a168 + k9[i] * a169 + k4[i] * a1613 + k10[i] * a1614 + k2[i] * a1615);
            }

            fcn(t + dt * c16, xtemp, k3);

            // Final preparation
            for (i = 0; i < n; ++i)
            {
                r[i + n * 4] = dt * (r[i + n * 4] + k4[i] * d413 + k10[i] * d414 + k2[i] * d415 + k3[i] * d416);
                r[i + n * 5] = dt * (r[i + n * 5] + k4[i] * d513 + k10[i] * d514 + k2[i] * d515 + k3[i] * d516);
                r[i + n * 6] = dt * (r[i + n * 6] + k4[i] * d613 + k10[i] * d614 + k2[i] * d615 + k3[i] * d616);
                r[i + n * 7] = dt * (r[i + n * 7] + k4[i] * d713 + k10[i] * d714 + k2[i] * d715 + k3[i] * d716);
            }
        }

        private void Step(double t, double dt, double[] x)
        {
            int i;

            for (i = 0; i < n; ++i)
            {
                xtemp[i] = x[i] + dt * a21 * dxdt[i];
            }

            fcn(t + c2 * dt, xtemp, k2);
            for (i = 0; i < n; ++i)
            {
                xtemp[i] = x[i] + dt * (a31 * dxdt[i] + a32 * k2[i]);
            }

            fcn(t + c3 * dt, xtemp, k3);
            for (i = 0; i < n; ++i)
            {
                xtemp[i] = x[i] + dt * (a41 * dxdt[i] + a43 * k3[i]);
            }

            fcn(t + dt * c4, xtemp, k4);
            for (i = 0; i < n; ++i)
            {
                xtemp[i] = x[i] + dt * (dxdt[i] * a51 + k3[i] * a53 + k4[i] * a54);
            }

            fcn(t + dt * c5, xtemp, k5);
            for (i = 0; i < n; ++i)
            {
                xtemp[i] = x[i] + dt * (dxdt[i] * a61 + k4[i] * a64 + k5[i] * a65);
            }

            fcn(t + dt * c6, xtemp, k6);
            for (i = 0; i < n; ++i)
            {
                xtemp[i] = x[i] + dt * (dxdt[i] * a71 + k4[i] * a74 + k5[i] * a75 + k6[i] * a76);
            }

            fcn(t + dt * c7, xtemp, k7);
            for (i = 0; i < n; ++i)
            {
                xtemp[i] = x[i] + dt * (dxdt[i] * a81 + k4[i] * a84 + k5[i] * a85 + k6[i] * a86 + k7[i] * a87);
            }

            fcn(t + dt * c8, xtemp, k8);
            for (i = 0; i < n; ++i)
            {
                xtemp[i] = x[i] + dt * (dxdt[i] * a91 + k4[i] * a94 + k5[i] * a95 + k6[i] * a96 + k7[i] * a97 + k8[i] * a98);
            }

            fcn(t + dt * c9, xtemp, k9);
            for (i = 0; i < n; ++i)
            {
                xtemp[i] = x[i] + dt * (dxdt[i] * a101 + k4[i] * a104 + k5[i] * a105 + k6[i] * a106 + k7[i] * a107 + k8[i] * a108 + k9[i] * a109);
            }

            fcn(t + dt * c10, xtemp, k10);
            for (i = 0; i < n; ++i)
            {
                xtemp[i] = x[i] + dt * (dxdt[i] * a111 + k4[i] * a114 + k5[i] * a115 + k6[i] * a116 + k7[i] * a117 + k8[i] * a118 + k9[i] * a119 + k10[i] * a1110);
            }

            fcn(t + dt * c11, xtemp, k2);
            for (i = 0; i < n; ++i)
            {
                xtemp[i] = x[i] + dt * (dxdt[i] * a121 + k4[i] * a124 + k5[i] * a125 + k6[i] * a126 + k7[i] * a127 + k8[i] * a128 + k9[i] * a129 + k10[i] * a1210 + k2[i] * a122);
            }

            fcn(t + dt, xtemp, k3);
            for (i = 0; i < n; ++i)
            {
                k4[i] = dxdt[i] * b1 + k6[i] * b6 + k7[i] * b7 + k8[i] * b8 + k9[i] * b9 + k10[i] * b10 + k2[i] * b11 + k3[i] * b12;
                xout[i] = x[i] + dt * k4[i];
            }

            for (i = 0; i < n; ++i)
            {
                xerr[i] = dxdt[i] * e1 + k6[i] * e6 + k7[i] * e7 + k8[i] * e8 + k9[i] * e9 + k10[i] * e10 + k2[i] * e11 + k3[i] * e12;
                xerr2[i] = k4[i] - dxdt[i] * bh1 - k9[i] * bh2 - k3[i] * bh3;
            }
        }

        public double Error(double dt, double[] x)
        {
            double err = 0.0, err2 = 0.0, sk, deno, temp;

            for (int i = 0; i < n; ++i)
            {
                sk = atol + rtol * Math.Max(Math.Abs(x[i]), Math.Abs(xout[i]));

                temp = xerr2[i] / sk;
                err2 += temp * temp;

                temp = xerr[i] / sk;
                err += temp * temp;
            }

            deno = err + err2 * 0.01;
            if (deno <= 0.0)
            {
                deno = 1.0;
            }

            return Math.Abs(dt) * err * Math.Sqrt(1.0 / (n * deno));
        }
        
        /// <summary>
        /// Provides an approximation to the i-th component of the solution at t.
        /// </summary>
        /// <param name="i"></param>
        /// <param name="t"></param>
        /// <param name="con"></param>
        /// <returns></returns>
        public double Interpolate(int i, double t, double[] con)
        {
            int n = this.n;

            double s = (t - told) / dtold;
            double s1 = 1.0 - s;
            double conpar = con[4 * n + i] + s * (con[5 * n + i] + s1 * (con[6 * n + i] + s * con[7 * n + i]));

            return con[i] + s * (con[n + i] + s1 * (con[2 * n + i] + s * (con[3 * n + i] + s1 * conpar)));
        }
    }
}
