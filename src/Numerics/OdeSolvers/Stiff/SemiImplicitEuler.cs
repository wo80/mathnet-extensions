// Based on Fortran code SEULEX
// Copyright (c) 2004, Ernst Hairer
// License: Simplified BSD License (https://www.unige.ch/~hairer/software.html)

namespace MathNet.Numerics.OdeSolvers.Stiff
{
    using MathNet.Numerics.LinearAlgebra.Double;
    using MathNet.Numerics.LinearAlgebra.Double.Factorization;
    using System;
    using S_fp = System.Action<int, double, double[], double[]>;
    using J_fp = System.Action<int, double, double[], double[], int>;
    using M_fp = System.Action<int, double[], int>;

    /// <summary>
    /// 
    /// Authors: E. Hairer and G. Wanner
    ///
    /// This code is part of the book:
    ///         E. Hairer and G. Wanner
    ///         Solving Ordinary Differential Equations II.
    ///         Stiff and differential-algebraic problems. (2nd edition)
    ///         Springer-Verlag (1996)
    /// </summary>
    public class SemiImplicitEuler
    {

        double d_sign(double a, double b)
        {
            double x = Math.Abs(a);
            return b >= 0 ? x : -x;
        }

        struct coseu__1
        {
            public double xoldd, hhh;
            public int nnrd, kright;
        }
        struct coseu__2
        {
            public double xold, h;
            public int nrd, ir;
        }

        coseu__1 coseu_1;
        coseu__2 coseu_2;

        /* Table of constant values */

        /* Subroutine */
        public int seulex_(int n, S_fp fcn, int ifcn, double x,
            double[] y, double xend, double h, double[] rtol, double[] atol, int itol,
            J_fp jac, int ijac, M_fp mas, int imas, S_fp solout, int iout,
            double[] work, int lwork, int[] iwork, int liwork, int idid)
        {
            /* Local variables */
            int i, m1, m2, km, km2, nm1, lde, nrd;
            double fac1, fac2, fac3, fac4;
            int ndec, ijob;
            double hmax;
            int nmax;
            double thet;
            int nsol;
            double safe1, safe2;
            int ldjac;
            bool jband;
            double wkdec;
            double wkjac;
            int ldmas;
            double wkfcn;
            bool arret;
            int nstep, nsequ;
            double wksol, wkrow;
            int ldmas2, lambda, naccpt, nrejct;
            bool implct;
            int nrdens;

            bool autnms;
            double uround;


            /**
             *     NUMERICAL SOLUTION OF A STIFF (OR DIFFERENTIAL ALGEBRAIC)
             *     SYSTEM OF FIRST 0RDER ORDINARY DIFFERENTIAL EQUATIONS  MY'=F(X,Y).
             *     THIS IS AN EXTRAPOLATION-ALGORITHM, BASED ON THE
             *     LINEARLY IMPLICIT EULER METHOD (WITH STEP SIZE CONTROL
             *     AND ORDER SELECTION).
             *
             *     AUTHORS: E. HAIRER AND G. WANNER
             *              UNIVERSITE DE GENEVE, DEPT. DE MATHEMATIQUES
             *              CH-1211 GENEVE 24, SWITZERLAND
             *              E-MAIL:  Ernst.Hairer@math.unige.ch
             *                       Gerhard.Wanner@math.unige.ch
             *              INCLUSION OF DENSE OUTPUT BY E. HAIRER AND A. OSTERMANN
             *
             *     THIS CODE IS PART OF THE BOOK:
             *         E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
             *         EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
             *         SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS 14,
             *         SPRINGER-VERLAG 1991, SECOND EDITION 1996.
             *
             *     VERSION OF SEPTEMBER 30, 1995
             *         SMALL CORRECTIONS ON JUNE 11, 1999
             *
             *     INPUT PARAMETERS
             *     ----------------
             *     N           DIMENSION OF THE SYSTEM
             *
             *     FCN         NAME (EXTERNAL) OF SUBROUTINE COMPUTING THE
             *                 VALUE OF F(X,Y):
             *                    SUBROUTINE FCN(N,X,Y,F,RPAR,IPAR)
             *                    DOUBLE PRECISION X,Y(N),F(N)
             *                    F(1)=...   ETC.
             *                 RPAR, IPAR (SEE BELOW)
             *
             *     IFCN        GIVES INFORMATION ON FCN:
             *                    IFCN=0: F(X,Y) INDEPENDENT OF X (AUTONOMOUS)
             *                    IFCN=1: F(X,Y) MAY DEPEND ON X (NON-AUTONOMOUS)
             *
             *     X           INITIAL X-VALUE
             *
             *     Y(N)        INITIAL VALUES FOR Y
             *
             *     XEND        FINAL X-VALUE (XEND-X MAY BE POSITIVE OR NEGATIVE)
             *
             *     H           INITIAL STEP SIZE GUESS;
             *                 FOR STIFF EQUATIONS WITH INITIAL TRANSIENT,
             *                 H=1.D0/(NORM OF F'), USUALLY 1.D-2 OR 1.D-3, IS GOOD.
             *                 THIS CHOICE IS NOT VERY IMPORTANT, THE CODE QUICKLY
             *                 ADAPTS ITS STEP SIZE (IF H=0.D0, THE CODE PUTS H=1.D-6).
             *
             *     RTOL,ATOL   RELATIVE AND ABSOLUTE ERROR TOLERANCES. THEY
             *                 CAN BE BOTH SCALARS OR ELSE BOTH VECTORS OF LENGTH N.
             *
             *     ITOL        SWITCH FOR RTOL AND ATOL:
             *                   ITOL=0: BOTH RTOL AND ATOL ARE SCALARS.
             *                     THE CODE KEEPS, ROUGHLY, THE LOCAL ERROR OF
             *                     Y(I) BELOW RTOLaBS(Y(I))+ATOL
             *                   ITOL=1: BOTH RTOL AND ATOL ARE VECTORS.
             *                     THE CODE KEEPS THE LOCAL ERROR OF Y(I) BELOW
             *                     RTOL(I)aBS(Y(I))+ATOL(I).
             *
             *     JAC         NAME (EXTERNAL) OF THE SUBROUTINE WHICH COMPUTES
             *                 THE PARTIAL DERIVATIVES OF F(X,Y) WITH RESPECT TO Y
             *                 (THIS ROUTINE IS ONLY CALLED IF IJAC=1; SUPPLY
             *                 A DUMMY SUBROUTINE IN THE CASE IJAC=0).
             *                 FOR IJAC=1, THIS SUBROUTINE MUST HAVE THE FORM
             *                    SUBROUTINE JAC(N,X,Y,DFY,LDFY,RPAR,IPAR)
             *                    DOUBLE PRECISION X,Y(N),DFY(LDFY,N)
             *                    DFY(1,1)= 0...
             *                 LDFY, THE COLOMN-LENGTH OF THE ARRAY, IS
             *                 FURNISHED BY THE CALLING PROGRAM.
             *                 IF (MLJAC.EQ.N) THE JACOBIAN IS SUPPOSED TO
             *                    BE FULL AND THE PARTIAL DERIVATIVES ARE
             *                    STORED IN DFY AS
             *                       DFY(I,J) = PARTIAL F(I) / PARTIAL Y(J)
             *                 ELSE, THE JACOBIAN IS TAKEN AS BANDED AND
             *                    THE PARTIAL DERIVATIVES ARE STORED
             *                    DIAGONAL-WISE AS
             *                       DFY(I-J+MUJAC+1,J) = PARTIAL F(I) / PARTIAL Y(J).
             *
             *     IJAC        SWITCH FOR THE COMPUTATION OF THE JACOBIAN:
             *                    IJAC=0: JACOBIAN IS COMPUTED INTERNALLY BY FINITE
             *                       DIFFERENCES, SUBROUTINE "JAC" IS NEVER CALLED.
             *                    IJAC=1: JACOBIAN IS SUPPLIED BY SUBROUTINE JAC.
             *
             *     MLJAC       SWITCH FOR THE BANDED STRUCTURE OF THE JACOBIAN:
             *                    MLJAC=N: JACOBIAN IS A FULL MATRIX. THE LINEAR
             *                       ALGEBRA IS DONE BY FULL-MATRIX GAUSS-ELIMINATION.
             *                    0<=MLJAC<N: MLJAC IS THE LOWER BANDWITH OF JACOBIAN
             *                       MATRIX (>= NUMBER OF NON-ZERO DIAGONALS BELOW
             *                       THE MAIN DIAGONAL).
             *
             *     MUJAC       UPPER BANDWITH OF JACOBIAN  MATRIX (>= NUMBER OF NON-
             *                 ZERO DIAGONALS ABOVE THE MAIN DIAGONAL).
             *                 NEED NOT BE DEFINED IF MLJAC=N.
             *
             *     ----   MAS,IMAS,MLMAS, AND MUMAS HAVE ANALOG MEANINGS      -----
             *     ----   FOR THE "MASS MATRIX" (THE MATRIX "M" OF SECTION IV.8): -
             *
             *     MAS         NAME (EXTERNAL) OF SUBROUTINE COMPUTING THE MASS-
             *                 MATRIX M.
             *                 IF IMAS=0, THIS MATRIX IS ASSUMED TO BE THE IDENTITY
             *                 MATRIX AND NEEDS NOT TO BE DEFINED;
             *                 SUPPLY A DUMMY SUBROUTINE IN THIS CASE.
             *                 IF IMAS=1, THE SUBROUTINE MAS IS OF THE FORM
             *                    SUBROUTINE MAS(N,AM,LMAS,RPAR,IPAR)
             *                    DOUBLE PRECISION AM(LMAS,N)
             *                    AM(1,1)= 0....
             *                    IF (MLMAS.EQ.N) THE MASS-MATRIX IS STORED
             *                    AS FULL MATRIX LIKE
             *                         AM(I,J) = M(I,J)
             *                    ELSE, THE MATRIX IS TAKEN AS BANDED AND STORED
             *                    DIAGONAL-WISE AS
             *                         AM(I-J+MUMAS+1,J) = M(I,J).
             *
             *     IMAS       GIVES INFORMATION ON THE MASS-MATRIX:
             *                    IMAS=0: M IS SUPPOSED TO BE THE IDENTITY
             *                       MATRIX, MAS IS NEVER CALLED.
             *                    IMAS=1: MASS-MATRIX  IS SUPPLIED.
             *
             *     MLMAS       SWITCH FOR THE BANDED STRUCTURE OF THE MASS-MATRIX:
             *                    MLMAS=N: THE FULL MATRIX CASE. THE LINEAR
             *                       ALGEBRA IS DONE BY FULL-MATRIX GAUSS-ELIMINATION.
             *                    0<=MLMAS<N: MLMAS IS THE LOWER BANDWITH OF THE
             *                       MATRIX (>= NUMBER OF NON-ZERO DIAGONALS BELOW
             *                       THE MAIN DIAGONAL).
             *                 MLMAS IS SUPPOSED TO BE .LE. MLJAC.
             *
             *     MUMAS       UPPER BANDWITH OF MASS-MATRIX (>= NUMBER OF NON-
             *                 ZERO DIAGONALS ABOVE THE MAIN DIAGONAL).
             *                 NEED NOT BE DEFINED IF MLMAS=N.
             *                 MUMAS IS SUPPOSED TO BE .LE. MUJAC.
             *
             *     SOLOUT      NAME (EXTERNAL) OF SUBROUTINE PROVIDING THE
             *                 NUMERICAL SOLUTION DURING INTEGRATION.
             *                 IF IOUT>=1, IT IS CALLED AFTER EVERY SUCCESSFUL STEP.
             *                 SUPPLY A DUMMY SUBROUTINE IF IOUT=0.
             *                 IT MUST HAVE THE FORM
             *                    SUBROUTINE SOLOUT (NR,XOLD,X,Y,RC,LRC,IC,LIC,N,
             *                                       RPAR,IPAR,IRTRN)
             *                    DOUBLE PRECISION X,Y(N),RC(LRC),IC(LIC)
             *                    ....
             *                 SOLOUT FURNISHES THE SOLUTION "Y" AT THE NR-TH
             *                    GRID-POINT "X" (THEREBY THE INITIAL VALUE IS
             *                    THE FIRST GRID-POINT).
             *                 "XOLD" IS THE PRECEEDING GRID-POINT.
             *                 "IRTRN" SERVES TO INTERRUPT THE INTEGRATION. IF IRTRN
             *                    IS SET <0, SEULEX RETURNS TO THE CALLING PROGRAM.
             *                 DO NOT CHANGE THE ENTRIES OF RC(LRC),IC(LIC)!
             *
             *          -----  CONTINUOUS OUTPUT (IF IOUT=2): -----
             *                 DURING CALLS TO "SOLOUT", A CONTINUOUS SOLUTION
             *                 FOR THE INTERVAL [XOLD,X] IS AVAILABLE THROUGH
             *                 THE DOUBLE PRECISION FUNCTION
             *                        >>>   CONTEX(I,S,RC,LRC,IC,LIC)   <<<
             *                 WHICH PROVIDES AN APPROXIMATION TO THE I-TH
             *                 COMPONENT OF THE SOLUTION AT THE POINT S. THE VALUE
             *                 S SHOULD LIE IN THE INTERVAL [XOLD,X].
             *
             *     IOUT        GIVES INFORMATION ON THE SUBROUTINE SOLOUT:
             *                    IOUT=0: SUBROUTINE IS NEVER CALLED
             *                    IOUT=1: SUBROUTINE IS USED FOR OUTPUT
             *                    IOUT=2: DENSE OUTPUT IS PERFORMED IN SOLOUT
             *
             *     WORK        ARRAY OF WORKING SPACE OF LENGTH "LWORK".
             *                 SERVES AS WORKING SPACE FOR ALL VECTORS AND MATRICES.
             *                 "LWORK" MUST BE AT LEAST
             *                        N*(LJAC+LMAS+LE1+KM+8)+4kM+20+KM2nRDENS
             *                 WHERE
             *                    KM2=2+KM*(KM+3)/2  AND  NRDENS=IWORK(6) (SEE BELOW)
             *                 AND
             *                    LJAC=N              IF MLJAC=N (FULL JACOBIAN)
             *                    LJAC=MLJAC+MUJAC+1  IF MLJAC<N (BANDED JAC.0)
             *                 AND
             *                    LMAS=0              IF IMAS=0
             *                    LMAS=N              IF IMAS=1 AND MLMAS=N (FULL)
             *                    LMAS=MLMAS+MUMAS+1  IF MLMAS<N (BANDED MASS-M.0)
             *                 AND
             *                    LE1=N               IF MLJAC=N (FULL JACOBIAN)
             *                    LE1=2mLJAC+MUJAC+1 IF MLJAC<N (BANDED JAC.0).
             *                 AND
             *                    KM=12               IF IWORK(3)=0
             *                    KM=IWORK(3)         IF IWORK(3).GT.0
             *
             *                 IN THE USUAL CASE WHERE THE JACOBIAN IS FULL AND THE
             *                 MASS-MATRIX IS THE INDENTITY (IMAS=0), THE MINIMUM
             *                 STORAGE REQUIREMENT IS
             *                         LWORK = 2nn+(KM+8)n+4kM+13+KM2nRDENS.
             *                 IF IWORK(9)=M1>0 THEN "LWORK" MUST BE AT LEAST
             *                    N*(LJAC+KM+8)+(N-M1)*(LMAS+LE1)+4kM+20+KM2nRDENS
             *                 WHERE IN THE DEFINITIONS OF LJAC, LMAS AND LE1 THE
             *                 NUMBER N CAN BE REPLACED BY N-M1.
             *
             *     LWORK       DECLARED LENGTH OF ARRAY "WORK".
             *
             *     IWORK       int WORKING SPACE OF LENGTH "LIWORK".
             *                 "LIWORK" MUST BE AT LEAST  2n+KM+20+NRDENS.
             *
             *     LIWORK      DECLARED LENGTH OF ARRAY "IWORK".
             *
             *     RPAR, IPAR  REAL AND int PARAMETERS (OR PARAMETER ARRAYS) WHICH
             *                 CAN BE USED FOR COMMUNICATION BETWEEN YOUR CALLING
             *                 PROGRAM AND THE FCN, JAC, MAS, SOLOUT SUBROUTINES.
             *
             * ----------------------------------------------------------------------
             *
             *     SOPHISTICATED SETTING OF PARAMETERS
             *     -----------------------------------
             *              SEVERAL PARAMETERS OF THE CODE ARE TUNED TO MAKE IT WORK
             *              WELL. THEY MAY BE DEFINED BY SETTING WORK(1),..,WORK(13)
             *              AS WELL AS IWORK(1),..,IWORK(4) DIFFERENT FROM ZERO.
             *              FOR ZERO INPUT, THE CODE CHOOSES DEFAULT VALUES:
             *
             *    IWORK(1)  IF IWORK(1).NE.0, THE CODE TRANSFORMS THE JACOBIAN
             *              MATRIX TO HESSENBERG FORM. THIS IS PARTICULARLY
             *              ADVANTAGEOUS FOR LARGE SYSTEMS WITH FULL JACOBIAN.
             *              IT DOES NOT WORK FOR BANDED JACOBIAN (MLJAC<N)
             *              AND NOT FOR IMPLICIT SYSTEMS (IMAS=1).
             *
             *    IWORK(2)  THIS IS THE MAXIMAL NUMBER OF ALLOWED STEPS.
             *              THE DEFAULT VALUE (FOR IWORK(2)=0) IS 100000.
             *
             *    IWORK(3)  THE MAXIMUM NUMBER OF COLUMNS IN THE EXTRAPOLATION
             *              TABLE. THE DEFAULT VALUE (FOR IWORK(3)=0) IS 12.
             *              IF IWORK(3).NE.0 THEN IWORK(3) SHOULD BE .GE.3.
             *
             *    IWORK(4)  SWITCH FOR THE STEP SIZE SEQUENCE
             *              IF IWORK(4).EQ.1 THEN 1,2,3,4,6,8,12,16,24,32,48,...
             *              IF IWORK(4).EQ.2 THEN 2,3,4,6,8,12,16,24,32,48,64,...
             *              IF IWORK(4).EQ.3 THEN 1,2,3,4,5,6,7,8,9,10,...
             *              IF IWORK(4).EQ.4 THEN 2,3,4,5,6,7,8,9,10,11,...
             *              THE DEFAULT VALUE (FOR IWORK(4)=0) IS IWORK(4)=2.
             *
             *    IWORK(5)  PARAMETER "LAMBDA" OF DENSE OUTPUT; POSSIBLE VALUES
             *              ARE 0 AND 1; DEFAULT IWORK(5)=0.
             *
             *    IWORK(6)  = NRDENS = NUMBER OF COMPONENTS, FOR WHICH DENSE OUTPUT
             *              IS REQUIRED
             *
             *    IWORK(21),...,IWORK(NRDENS+20) INDICATE THE COMPONENTS, FOR WHICH
             *              DENSE OUTPUT IS REQUIRED
             *
             *       IF THE DIFFERENTIAL SYSTEM HAS THE SPECIAL STRUCTURE THAT
             *            Y(I)' = Y(I+M2)   FOR  I=1,...,M1,
             *       WITH M1 A MULTIPLE OF M2, A SUBSTANTIAL GAIN IN COMPUTERTIME
             *       CAN BE ACHIEVED BY SETTING THE PARAMETERS IWORK(9) AND IWORK(10).
             *       E.G., FOR SECOND ORDER SYSTEMS P'=V, V'=G(P,V), WHERE P AND V ARE
             *       VECTORS OF DIMENSION N/2, ONE HAS TO PUT M1=M2=N/2.
             *       FOR M1>0 SOME OF THE INPUT PARAMETERS HAVE DIFFERENT MEANINGS:
             *       - JAC: ONLY THE ELEMENTS OF THE NON-TRIVIAL PART OF THE
             *              JACOBIAN HAVE TO BE STORED
             *              IF (MLJAC.EQ.N-M1) THE JACOBIAN IS SUPPOSED TO BE FULL
             *                 DFY(I,J) = PARTIAL F(I+M1) / PARTIAL Y(J)
             *                FOR I=1,N-M1 AND J=1,N.
             *              ELSE, THE JACOBIAN IS BANDED ( M1 = M2 * MM )
             *                 DFY(I-J+MUJAC+1,J+Km2) = PARTIAL F(I+M1) / PARTIAL Y(J+Km2)
             *                FOR I=1,MLJAC+MUJAC+1 AND J=1,M2 AND K=0,MM.
             *       - MLJAC: MLJAC=N-M1: IF THE NON-TRIVIAL PART OF THE JACOBIAN IS FULL
             *                0<=MLJAC<N-M1: IF THE (MM+1) SUBMATRICES (FOR K=0,MM)
             *                     PARTIAL F(I+M1) / PARTIAL Y(J+Km2),  I,J=1,M2
             *                    ARE BANDED, MLJAC IS THE MAXIMAL LOWER BANDWIDTH
             *                    OF THESE MM+1 SUBMATRICES
             *       - MUJAC: MAXIMAL UPPER BANDWIDTH OF THESE MM+1 SUBMATRICES
             *                NEED NOT BE DEFINED IF MLJAC=N-M1
             *       - MAS: IF IMAS=0 THIS MATRIX IS ASSUMED TO BE THE IDENTITY AND
             *              NEED NOT BE DEFINED. SUPPLY A DUMMY SUBROUTINE IN THIS CASE.
             *              IT IS ASSUMED THAT ONLY THE ELEMENTS OF RIGHT LOWER BLOCK OF
             *              DIMENSION N-M1 DIFFER FROM THAT OF THE IDENTITY MATRIX.
             *              IF (MLMAS.EQ.N-M1) THIS SUBMATRIX IS SUPPOSED TO BE FULL
             *                 AM(I,J) = M(I+M1,J+M1)     FOR I=1,N-M1 AND J=1,N-M1.
             *              ELSE, THE MASS MATRIX IS BANDED
             *                 AM(I-J+MUMAS+1,J) = M(I+M1,J+M1)
             *       - MLMAS: MLMAS=N-M1: IF THE NON-TRIVIAL PART OF M IS FULL
             *                0<=MLMAS<N-M1: LOWER BANDWIDTH OF THE MASS MATRIX
             *       - MUMAS: UPPER BANDWIDTH OF THE MASS MATRIX
             *                NEED NOT BE DEFINED IF MLMAS=N-M1
             *
             *    IWORK(9)  THE VALUE OF M1.  DEFAULT M1=0.
             *
             *    IWORK(10) THE VALUE OF M2.  DEFAULT M2=M1.
             *
             *    WORK(1)   UROUND, THE ROUNDING UNIT, DEFAULT 1.D-16.
             *
             *    WORK(2)   MAXIMAL STEP SIZE, DEFAULT XEND-X.
             *
             *    WORK(3)   DECIDES WHETHER THE JACOBIAN SHOULD BE RECOMPUTED;
             *              INCREASE WORK(3), TO 0.01 SAY, WHEN JACOBIAN EVALUATIONS
             *              ARE COSTLY. FOR SMALL SYSTEMS WORK(3) SHOULD BE SMALLER.
             *              DEFAULT MIN(1.0D-4,RTOL(1))
             *
             *    WORK(4), WORK(5)   PARAMETERS FOR STEP SIZE SELECTION
             *              THE NEW STEP SIZE FOR THE J-TH DIAGONAL ENTRY IS
             *              CHOSEN SUBJECT TO THE RESTRICTION
             *                 FACMIN/WORK(5) <= HNEW(J)/HOLD <= 1/FACMIN
             *              WHERE FACMIN=WORK(4)**(1/(J-1))
             *              DEFAULT VALUES: WORK(4)=0.1D0, WORK(5)=4.D0
             *
             *    WORK(6), WORK(7)   PARAMETERS FOR THE ORDER SELECTION
             *              ORDER IS DECREASED IF    W(K-1) <= W(K)*WORK(6)
             *              ORDER IS INCREASED IF    W(K) <= W(K-1)*WORK(7)
             *              DEFAULT VALUES: WORK(6)=0.7D0, WORK(7)=0.9D0
             *
             *    WORK(8), WORK(9)   SAFETY FACTORS FOR STEP CONTROL ALGORITHM
             *             HNEW=H*WORK(9)*(WORK(8)tOL/ERR)**(1/(J-1))
             *              DEFAULT VALUES: WORK(8)=0.8D0, WORK(9)=0.93D0
             *
             *    WORK(10), WORK(11), WORK(12), WORK(13)   ESTIMATED WORKS FOR
             *             A CALL TO  FCN, JAC, DEC, SOL, RESPECTIVELY.
             *             DEFAULT VALUES ARE: WORK(10)=1.D0, WORK(11)=5.D0,
             *             WORK(12)=1.D0, WORK(13)=1.D0.
             *
             * -----------------------------------------------------------------------
             *
             *     OUTPUT PARAMETERS
             *     -----------------
             *     X           X-VALUE WHERE THE SOLUTION IS COMPUTED
             *                 (AFTER SUCCESSFUL RETURN X=XEND)
             *
             *     Y(N)        SOLUTION AT X
             *
             *     H           PREDICTED STEP SIZE OF THE LAST ACCEPTED STEP
             *
             *     IDID        REPORTS ON SUCCESSFULNESS UPON RETURN:
             *                   IDID=1  COMPUTATION SUCCESSFUL,
             *                   IDID=-1 COMPUTATION UNSUCCESSFUL.
             *
             *   IWORK(14)  NFCN    NUMBER OF FUNCTION EVALUATIONS (THOSE FOR NUMERICAL
             *                      EVALUATION OF THE JACOBIAN ARE NOT COUNTED)
             *   IWORK(15)  NJAC    NUMBER OF JACOBIAN EVALUATIONS (EITHER ANALYTICALLY
             *                      OR NUMERICALLY)
             *   IWORK(16)  NSTEP   NUMBER OF COMPUTED STEPS
             *   IWORK(17)  NACCPT  NUMBER OF ACCEPTED STEPS
             *   IWORK(18)  NREJCT  NUMBER OF REJECTED STEPS (AFTER AT LEAST ONE STEP
             *                      HAS BEEN ACCEPTED)
             *   IWORK(19)  NDEC    NUMBER OF LU-DECOMPOSITIONS (N-DIMENSIONAL MATRIX)
             *   IWORK(20)  NSOL    NUMBER OF FORWARD-BACKWARD SUBSTITUTIONS
             */

            /* Parameter adjustments */
            //--y;
            //--rtol;
            //--atol;
            //--work;
            //--iwork;
            //--rpar;
            //--ipar;

            /* Function Body */
            nstep = 0;
            naccpt = 0;
            nrejct = 0;
            ndec = 0;
            nsol = 0;
            arret = false;
            /* -------- NMAX , THE MAXIMAL NUMBER OF STEPS ----- */
            if (iwork[1] == 0)
            {
                nmax = 100000;
            }
            else
            {
                nmax = iwork[1];
                if (nmax <= 0)
                {

                    Console.WriteLine(" WRONG INPUT IWORK(2)=", iwork[1]);

                    arret = true;
                }
            }
            /* -------- KM     MAXIMUM NUMBER OF COLUMNS IN THE EXTRAPOLATION */
            if (iwork[2] == 0)
            {
                km = 12;
            }
            else
            {
                km = iwork[2];
                if (km <= 2)
                {

                    Console.WriteLine(" CURIOUS INPUT IWORK(3)=", iwork[2]);

                    arret = true;
                }
            }
            /* -------- NSEQU     CHOICE OF STEP SIZE SEQUENCE */
            nsequ = iwork[3];
            if (iwork[3] == 0)
            {
                nsequ = 2;
            }
            if (nsequ <= 0 || nsequ >= 5)
            {

                Console.WriteLine(" CURIOUS INPUT IWORK(4)=", iwork[3]);

                arret = true;
            }
            /* -------- LAMBDA   PARAMETER FOR DENSE OUTPUT */
            lambda = iwork[4];
            if (lambda < 0 || lambda >= 2)
            {

                Console.WriteLine(" CURIOUS INPUT IWORK(5)=", iwork[4]);

                arret = true;
            }
            /* -------- NRDENS   NUMBER OF DENSE OUTPUT COMPONENTS */
            nrdens = iwork[5];
            if (nrdens < 0 || nrdens > n)
            {

                Console.WriteLine(" CURIOUS INPUT IWORK(6)=", iwork[5]);

                arret = true;
            }
            /* -------- PARAMETER FOR SECOND ORDER EQUATIONS */
            m1 = iwork[8];
            m2 = iwork[9];
            nm1 = n - m1;
            if (m1 == 0)
            {
                m2 = n;
            }
            if (m2 == 0)
            {
                m2 = m1;
            }
            if (m1 < 0 || m2 < 0 || m1 + m2 > n)
            {

                Console.WriteLine(" CURIOUS INPUT FOR IWORK(9,10)=", m1, m2);

                arret = true;
            }
            /* -------- UROUND   SMALLEST NUMBER SATISFYING 1.D0+UROUND>1.D0 */
            if (work[0] == 0.0)
            {
                uround = 1e-16;
            }
            else
            {
                uround = work[0];
                if (uround <= 0.0 || uround >= 1.0)
                {

                    Console.WriteLine("  UROUND=", work);

                    arret = true;
                }
            }
            /* -------- MAXIMAL STEP SIZE */
            if (work[1] == 0.0)
            {
                hmax = xend - x;
            }
            else
            {
                hmax = work[1];
            }
            /* ------ THET     DECIDES WHETHER THE JACOBIAN SHOULD BE RECOMPUTED; */
            if (work[2] == 0.0)
            {
                thet = Math.Min(1e-4, rtol[0]);
            }
            else
            {
                thet = work[2];
            }
            /* -------  FAC1,FAC2     PARAMETERS FOR STEP SIZE SELECTION */
            if (work[3] == 0.0)
            {
                fac1 = 0.1;
            }
            else
            {
                fac1 = work[3];
            }
            if (work[4] == 0.0)
            {
                fac2 = 4.0;
            }
            else
            {
                fac2 = work[4];
            }
            /* -------  FAC3, FAC4   PARAMETERS FOR THE ORDER SELECTION */
            if (work[5] == 0.0)
            {
                fac3 = 0.7;
            }
            else
            {
                fac3 = work[5];
            }
            if (work[6] == 0.0)
            {
                fac4 = 0.9;
            }
            else
            {
                fac4 = work[6];
            }
            /* ------- SAFE1, SAFE2 SAFETY FACTORS FOR STEP SIZE PREDICTION */
            if (work[7] == 0.0)
            {
                safe1 = 0.6;
            }
            else
            {
                safe1 = work[7];
            }
            if (work[8] == 0.0)
            {
                safe2 = 0.93;
            }
            else
            {
                safe2 = work[8];
            }
            /* ------- WKFCN,WKJAC,WKDEC,WKSOL  ESTIMATED WORK FOR  FCN,JAC,DEC,SOL */
            if (work[9] == 0.0)
            {
                wkfcn = 1.0;
            }
            else
            {
                wkfcn = work[9];
            }
            if (work[10] == 0.0)
            {
                wkjac = 5.0;
            }
            else
            {
                wkjac = work[10];
            }
            if (work[11] == 0.0)
            {
                wkdec = 1.0;
            }
            else
            {
                wkdec = work[11];
            }
            if (work[12] == 0.0)
            {
                wksol = 1.0;
            }
            else
            {
                wksol = work[12];
            }
            wkrow = wkfcn + wksol;
            /* --------- CHECK IF TOLERANCES ARE O.K. */
            if (itol == 0)
            {
                if (atol[0] <= 0.0 || rtol[0] <= uround * 10.0)
                {

                    Console.WriteLine(" TOLERANCES ARE TOO SMALL");

                    arret = true;
                }
            }
            else
            {
                for (i = 0; i < n; ++i)
                {
                    if (atol[i] <= 0.0 || rtol[i] <= uround * 10.0)
                    {

                        Console.WriteLine(" TOLERANCES(" + i + ") ARE TOO SMALL");

                        arret = true;
                    }
                }
            }
            /* *** *** *** *** *** *** *** *** *** *** *** *** *** */
            /*         COMPUTATION OF ARRAY ENTRIES */
            /* *** *** *** *** *** *** *** *** *** *** *** *** *** */
            /* ---- AUTONOMOUS, IMPLICIT, BANDED OR NOT ? */
            autnms = ifcn == 0;
            implct = imas != 0;

            /* -------- COMPUTATION OF THE ROW-DIMENSIONS OF THE 2-ARRAYS --- */
            /* -- JACOBIAN AND MATRIX E */
            {
                ldjac = nm1;
                lde = nm1;
            }
            /* -- MASS MATRIX */
            if (implct)
            {
                {
                    ldmas = nm1;
                    ijob = 5;
                }
            }
            else
            {
                ldmas = 0;
                {
                    ijob = 1;
                    if (n > 2 && iwork[0] != 0)
                    {
                        ijob = 7;
                    }
                }
            }
            ldmas2 = Math.Max(1, ldmas);
            /* ------ HESSENBERG OPTION ONLY FOR EXPLICIT EQU. WITH FULL JACOBIAN */
            if (implct && ijob == 7)
            {

                Console.WriteLine(" HESSENBERG OPTION ONLY FOR EXPLICIT EQUATIONS WITH FULL JACOBIAN");

                arret = true;
            }
            /* ------- PREPARE THE ENTRY-POINTS FOR THE ARRAYS IN WORK ----- */
            km2 = km * (km + 1) / 2;
            var yh = new double[n];// 21;
            var dy = new double[n];// yh + n;
            var fx = new double[n];// dy + n;
            var yhh = new double[n];// fx + n;
            var dyh = new double[n];// yhh + n;
            var del = new double[n];// dyh + n;
            var wh = new double[n];// del + n;
            var scal = new double[n];// wh + n;
            var hh = new double[km];// scal + n;
            var w = new double[km];// hh + km;
            var a = new double[km];// w + km;
            var _jac = new double[n * ldjac];// a + km;
            var e = new double[nm1 * lde];// jac + n * ldjac;
            var _mas = new double[nm1 * ldmas];// e + nm1 * lde;
            var t = new double[n * km];// mas + nm1 * ldmas;
            var ifac = new double[km];// t + n * km;
            var de = new double[(km + 2) * nrdens];// ifac + km;
            var fsafe = new double[km2 * nrdens];// de + (km + 2) * nrdens;

            /* ------ TOTAL STORAGE REQUIREMENT ----------- */
            //istore = fsafe + km2 * nrdens - 1;
            //if (istore > lwork)
            //{
            //    Console.WriteLine(" INSUFFICIENT STORAGE FOR WORK, MIN. LWORK=", istore);
            //    arret = true;
            //}

            /* ------- ENTRY POINTS FOR int WORKSPACE ----- */
            var co = new int[nrdens];// 21;
            var ip = new int[n];// nrdens + 21;
            var nj = new int[km];// ip + n;
            var iph = new int[km];// nj + km;

            /* --------- TOTAL REQUIREMENT --------------- */
            //istore = co + km - 1;
            //if (istore > liwork)
            //{
            //    Console.WriteLine(" INSUFF. STORAGE FOR IWORK, MIN. LIWORK=", istore);
            //    arret = true;
            //}

            /* ------ WHEN A FAIL HAS OCCURED, WE RETURN WITH IDID=-1 */
            if (arret)
            {
                idid = -1;
                return 0;
            }
            nrd = Math.Max(1, nrdens);
            /* -------- CALL TO CORE INTEGRATOR ------------ */
            seucor_(n, fcn, x, y, xend, hmax, h, km, rtol, atol,
                 itol, jac, ijac, mas,
                 solout, iout, idid, ijob, m1, m2, nm1, nmax, uround,
                nsequ, autnms, implct, ldjac, lde, ldmas2, yh,
                dy, fx, yhh, dyh, del,
                 wh, scal, hh, w,
                a, _jac, e, _mas, t, ip,
                 nj, iph, fac1, fac2, fac3, fac4,
                 thet, safe1, safe2, wkjac, wkdec, wkrow, km2, nrd, ifac,
                 fsafe, lambda, ref nstep, ref naccpt,
                 ref nrejct, ref ndec, ref nsol, de, co);

            iwork[16] = nstep;
            iwork[17] = naccpt;
            iwork[18] = nrejct;
            iwork[19] = ndec;
            iwork[20] = nsol;
            /* ----------- RETURN ----------- */
            return 0;
        } /* seulex_ */



        /*  ----- ... AND HERE IS THE CORE INTEGRATOR  ---------- */

        /* Subroutine */
        int seucor_(int n, S_fp fcn, double x, double[]
y, double xend, double hmax, double h, int km,
double[] rtol, double[] atol, int itol, J_fp jac, int
ijac, M_fp mas,S_fp solout, int iout, int idid, int ijob,
int m1, int m2, int nm1, int nmax, double
uround, int nsequ, bool autnms, bool implct,
int lfjac, int le, int ldmas, double[] yh,
double[] dy, double[] fx, double[] yhh, double[] dyh,
double[] del, double[] wh, double[] scal, double[] hh,
double[] w, double[] a, double[] fjac, double[] e,
double[] fmas, double[] t, int[] ip, int[] nj, int[]
iphes, double fac1, double fac2, double fac3,
double fac4, double thet, double safe1, double
safe2, double wkjac, double wkdec, double wkrow,
int km2, int nrd, double[] facul, double[] fsafe,
int lambda, ref int nstep, ref int naccpt, ref int nrejct, ref int ndec, ref int nsol,
double[] dens, int[] icomp)
        {
            /* Format strings */
            //static char fmt_979[] = "(\002 EXIT OF SEULEX AT X=\002,d14.7,\002   H=\002,d14.7)";

            /* System generated locals */
            int i1, i2, i3;
            double d1;

            /* Builtin functions */
            //double d_sign(double *, double *), d_lg10(double *), Math.Sqrt(double), pow_di(double *, int *);

            /* Local variables */
            int i, j, k, l, kc = 0, ii, kk, i_n, mm, kx;
            double t1i, fac = 0, err;
            int ipt, klr, krn, lbeg, lend, lrde;
            double delt;
            bool last;
            double xold;
            bool atov;
            double hopt = 0;

            int kopt;
            double xsol, facnj, theta, ysafe, hmaxn;
            int nrsol, irtrn, nsolu;
            bool caljac;
            double dblenj;
            bool calhes = false;
            bool reject;
            double factor;
            double errold = 0, posneg;


            /* ---------------------------------------------------------- */
            /*     CORE INTEGRATOR FOR SEULEX */
            /*     PARAMETERS SAME AS IN SEULEX WITH WORKSPACE ADDED */
            /* ---------------------------------------------------------- */
            /*         DECLARATIONS */
            /* ---------------------------------------------------------- */
            /* --- COMPUTE COEFFICIENTS FOR DENSE OUTPUT */
            /* Parameter adjustments */
            //--iphes;
            //--scal;
            //--wh;
            //--del;
            //--dyh;
            //--yhh;
            //--fx;
            //--dy;
            //--yh;
            //--y;
            //--facul;
            //--nj;
            //t -= 1 + km;
            //--a;
            //--w;
            //--hh;
            //--rtol;
            //--atol;
            //--ip;
            //fjac -= 1 + lfjac;
            //1 + le = 1 + le;
            //e -= 1 + le;
            //1 + ldmas = 1 + ldmas;
            //fmas -= 1 + ldmas;
            //--icomp;
            //--dens;
            //fsafe -= 1 + km2;
            //--rpar;
            //--ipar;

            /* Function Body */
            if (iout == 2)
            {
                coseu_1.nnrd = nrd;
                /* --- COMPUTE THE FACTORIALS -------- */
                facul[0] = 1.0;
                for (i = 0; i < km - 1; ++i)
                {
                    facul[i + 1] = i * facul[i];
                }
            }
            /* ------- COMPUTE MASS MATRIX FOR IMPLICIT CASE ---------- */
            if (implct)
            {
                mas(nm1, fmas, ldmas);
            }
            /* *** *** *** *** *** *** *** */
            /*  INITIALISATIONS */
            /* *** *** *** *** *** *** *** */
            lrde = (km + 2) * nrd;

            if (m1 > 0)
            {
                ijob += 10;
            }
            /* --- DEFINE THE STEP SIZE SEQUENCE */
            if (nsequ == 1)
            {
                nj[0] = 1;
                nj[1] = 2;
                nj[2] = 3;
                for (i = 3; i < km; ++i)
                {
                    nj[i] = nj[i - 2] << 1;
                }
            }
            if (nsequ == 2)
            {
                nj[0] = 2;
                nj[1] = 3;
                for (i = 2; i < km; ++i)
                {
                    nj[i] = nj[i - 2] << 1;
                }
            }
            for (i = 0; i < km; ++i)
            {
                if (nsequ == 3)
                {
                    nj[i] = i;
                }
                if (nsequ == 4)
                {
                    nj[i] = i + 1;
                }
            }
            a[0] = wkjac + nj[0] * wkrow + wkdec;
            for (i = 1; i < km; ++i)
            {
                a[i] = a[i - 1] + (nj[i] - 1) * wkrow + wkdec;
            }
            d1 = xend - x;
            posneg = d_sign(1.0, d1);
            k = Math.Max(2, Math.Min(km - 2, (int)(-Math.Log10(rtol[0] + atol[0]) * 0.6 + 1.5)));
            hmaxn = Math.Min(Math.Abs(hmax), Math.Abs(xend - x));
            h = Math.Max(Math.Abs(h), 1e-6);
            h = posneg * Math.Min(h, hmaxn);
            theta = Math.Abs(thet) * 2;
            if (iout != 0)
            {
                irtrn = 1;
                nrsol = 1;
                xsol = x;
                for (i = 0; i < n; ++i)
                {
                    yh[i] = y[i];
                }
                nsolu = n;
                xold = x;
                //solout(nrsol, xold, xsol, yh, dens, lrde, icomp, nrd, nsolu, irtrn);
                if (irtrn < 0)
                {
                    goto L120;
                }
                ipt = 0;
            }
            err = 0.0;
            w[0] = 1e30;
            for (i = 0; i < n; ++i)
            {
                if (itol == 0)
                {
                    scal[i] = atol[0] + rtol[0] * Math.Abs(y[i]);
                }
                else
                {
                    scal[i] = atol[i] + rtol[i] * Math.Abs(y[i]);
                }
            }
            caljac = false;
            reject = false;
            last = false;
            L10:
            if (reject)
            {
                theta = Math.Abs(thet) * 2;
            }
            atov = false;
            /* *** *** *** *** *** *** *** */
            /* --- IS XEND REACHED IN THE NEXT STEP? */
            /* *** *** *** *** *** *** *** */
            if (Math.Abs(xend - x) * 0.1 <= Math.Abs(x) * uround)
            {
                goto L110;
            }
            hopt = h;
            h = posneg * Math.Min(Math.Min(Math.Abs(h), Math.Abs(xend - x)), hmaxn);
            if ((x + h * 1.01 - xend) * posneg > 0.0)
            {
                h = xend - x;
                last = true;
            }
            if (autnms)
            {
                fcn(n, x, y, dy);
            }
            if (theta > thet && !caljac)
            {
                /* *** *** *** *** *** *** *** */
                /*  COMPUTATION OF THE JACOBIAN */
                /* *** *** *** *** *** *** *** */
                if (ijac == 0)
                {
                    /* --- COMPUTE JACOBIAN MATRIX NUMERICALLY */
                    if (!(autnms))
                    {
                        fcn(n, x, y, dy);
                    }

                    /* --- JACOBIAN IS FULL */
                    for (i = 0; i < n; ++i)
                    {
                        ysafe = y[i];
                        delt = Math.Sqrt(uround * Math.Max(1e-5, Math.Abs(ysafe)));
                        y[i] = ysafe + delt;
                        fcn(n, x, y, yh);
                        for (j = m1; j < n; ++j)
                        {
                            fjac[j - m1 + i * lfjac] = (yh[j] - dy[j]) / delt;
                        }
                        y[i] = ysafe;
                    }
                }
                else
                {
                    /* --- COMPUTE JACOBIAN MATRIX ANALYTICALLY */
                    jac(n, x, y, fjac, lfjac);
                }
                caljac = true;
                calhes = true;
            }
            /* *** *** *** *** *** *** *** */
            /* --- THE FIRST AND LAST STEP */
            /* *** *** *** *** *** *** *** */
            if (nstep == 0 || last)
            {
                ipt = 0;
                ++(nstep);
                for (j = 0; j < k; ++j)
                {
                    kc = j;
                    seul_(j, n, fcn, x, y, dy, fx, fjac,
                        lfjac, fmas, ldmas, e,
                         le, ip, h, km, hmaxn, t, scal,
                         nj, hh, w, a, yhh, dyh, del,
                        wh, err, safe1, fac, fac1, fac2, safe2, theta,
                        ref ndec, ref nsol, errold,
                        iphes, icomp, autnms, implct, reject,
                        atov, fsafe, km2, nrd, iout, ipt, m1, m2,
                        nm1, ijob, calhes);
                    if (atov)
                    {
                        goto L10;
                    }
                    if (j > 1 && err <= 1.0)
                    {
                        goto L60;
                    }
                }
                goto L55;
            }
            /* --- BASIC INTEGRATION STEP */
            L30:
            ipt = 0;
            ++(nstep);
            if (nstep >= nmax)
            {
                goto L120;
            }
            kc = k - 1;
            for (j = 0; j < kc; ++j)
            {
                seul_(j, n, fcn, x, y, dy, fx, fjac,
                    lfjac, fmas, ldmas, e, le, ip,
                    h, km, hmaxn, t, scal, nj, hh, w
                    , a, yhh, dyh, del, wh, err, safe1, fac,
                     fac1, fac2, safe2, theta, ref ndec, ref nsol,
                    errold, iphes, icomp, autnms, implct,
                    reject, atov, fsafe, km2, nrd, iout,
                    ipt, m1, m2, nm1, ijob, calhes);
                if (atov)
                {
                    goto L10;
                }
            }
            /* *** *** *** *** *** *** *** */
            /* --- CONVERGENCE MONITOR */
            /* *** *** *** *** *** *** *** */
            if (k == 2 || reject)
            {
                goto L50;
            }
            if (err <= 1.0)
            {
                goto L60;
            }
            if (err > (double)(nj[k + 1] * nj[k]) * 4.0)
            {
                goto L100;
            }
            L50:
            seul_(k, n, fcn, x, y, dy, fx, fjac,
                lfjac, fmas, ldmas, e, le, ip, h,
                km, hmaxn, t, scal, nj, hh, w, a,
                yhh, dyh, del, wh, err, safe1, fac, fac1, fac2,
                safe2, theta, ref ndec, ref nsol, errold,
                iphes, icomp, autnms, implct, reject, atov,
                fsafe, km2, nrd, iout, ipt, m1, m2, nm1, ijob,
                calhes);
            if (atov)
            {
                goto L10;
            }
            kc = k;
            if (err <= 1.0)
            {
                goto L60;
            }
            /* --- HOPE FOR CONVERGENCE IN LINE K+1 */
            L55:
            if (err > (double)nj[k + 1] * 2.0)
            {
                goto L100;
            }
            kc = k + 1;
            seul_(kc, n, fcn, x, y, dy, fx, fjac,
                lfjac, fmas, ldmas, e, le, ip, h,
                km, hmaxn, t, scal, nj, hh, w, a,
                yhh, dyh, del, wh, err, safe1, fac, fac1, fac2,
                safe2, theta, ref ndec, ref nsol, errold,
                iphes, icomp, autnms, implct, reject, atov,
                fsafe, km2, nrd, iout, ipt, m1, m2, nm1, ijob,
                calhes);
            if (atov)
            {
                goto L10;
            }
            if (err > 1.0)
            {
                goto L100;
            }
            /* *** *** *** *** *** *** *** */
            /* --- STEP IS ACCEPTED */
            /* *** *** *** *** *** *** *** */
            L60:
            xold = x;
            x += h;
            if (iout == 2)
            {
                coseu_1.kright = kc;
                for (i = 0; i < nrd; ++i)
                {
                    dens[i] = y[icomp[i]];
                }
            }
            for (i = 0; i < n; ++i)
            {
                t1i = t[i * km + 1];
                if (itol == 0)
                {
                    scal[i] = atol[0] + rtol[0] * Math.Abs(t1i);
                }
                else
                {
                    scal[i] = atol[i] + rtol[i] * Math.Abs(t1i);
                }
                y[i] = t1i;
            }
            ++(naccpt);
            caljac = false;
            if (iout == 2)
            {
                coseu_1.xoldd = xold;
                coseu_1.hhh = h;
                for (i = 0; i < nrd; ++i)
                {
                    dens[nrd + i] = y[icomp[i]];
                }
                i1 = coseu_1.kright - 1;
                for (klr = 1; klr <= i1; ++klr)
                {
                    /* --- COMPUTE DIFFERENCES */
                    if (klr >= 2)
                    {
                        i2 = kc;
                        for (kk = klr; kk <= i2; ++kk)
                        {
                            lbeg = (kk + 1) * kk / 2;
                            lend = lbeg - kk + 2;
                            for (l = lbeg; l >= lend; --l)
                            {
                                for (i = 0; i < nrd; ++i)
                                {
                                    fsafe[l + i * km2] -= fsafe[l - 1 + i * km2];
                                }
                            }
                        }
                    }
                    /* --- COMPUTE DERIVATIVES AT RIGHT END ---- */
                    i2 = kc;
                    for (kk = klr + lambda; kk <= i2; ++kk)
                    {
                        facnj = (double)nj[kk];
                        facnj = Math.Pow(facnj, klr) / facul[klr + 1]; // TODO: pow_di()
                        ipt = (kk + 1) * kk / 2;
                        for (i = 0; i < nrd; ++i)
                        {
                            krn = (kk - lambda + 1) * nrd;
                            dens[krn + i] = fsafe[ipt + i * km2] * facnj;
                        }
                    }
                    for (j = klr + lambda; j < kc; ++j)
                    {
                        dblenj = (double)nj[j];
                        i3 = klr + lambda + 1;
                        for (l = j; l >= i3; --l)
                        {
                            factor = dblenj / nj[l - 1] - 1.0;
                            for (i = 0; i < nrd; ++i)
                            {
                                krn = (l - lambda + 1) * nrd + i;
                                dens[krn - nrd] = dens[krn] + (dens[krn] - dens[krn - nrd]) / factor;
                            }
                        }
                    }
                }
                /* ---  COMPUTE THE COEFFICIENTS OF THE INTERPOLATION POLYNOMIAL */
                for (i_n = 0; i_n < nrd; ++i_n)
                {
                    for (j = 0; j < coseu_1.kright; ++j)
                    {
                        ii = nrd * j + i_n;
                        dens[ii] -= dens[ii - nrd];
                    }
                }
            }
            if (iout != 0)
            {
                irtrn = 1;
                nrsol = naccpt + 1;
                xsol = x;
                for (i = 0; i < n; ++i)
                {
                    yh[i] = y[i];
                }
                nsolu = n;
                //solout(nrsol, xold, xsol, yh, dens, lrde, icomp, nrd, nsolu, irtrn);
                if (irtrn < 0)
                {
                    goto L120;
                }
            }
            /* --- COMPUTE OPTIMAL ORDER */
            if (kc == 2)
            {
                kopt = Math.Min(3, km - 1);
                if (reject)
                {
                    kopt = 2;
                }
                goto L80;
            }
            if (kc <= k)
            {
                kopt = kc;
                if (w[kc - 1] < w[kc] * fac3)
                {
                    kopt = kc - 1;
                }
                if (w[kc] < w[kc - 1] * fac4)
                {
                    kopt = Math.Min(kc + 1, km - 1);
                }
            }
            else
            {
                kopt = kc - 1;
                if (kc > 3 && w[kc - 2] < w[kc - 1] * fac3)
                {
                    kopt = kc - 2;
                }
                if (w[kc] < w[kopt] * fac4)
                {
                    kopt = Math.Min(kc, km - 1);
                }
            }
            /* --- AFTER A REJECTED STEP */
            L80:
            if (reject)
            {
                k = Math.Min(kopt, kc);
                h = posneg * Math.Min(Math.Abs(h), Math.Abs(hh[k]));
                reject = false;
                goto L10;
            }
            /* --- COMPUTE STEP SIZE FOR NEXT STEP */
            if (kopt <= kc)
            {
                h = hh[kopt];
            }
            else
            {
                if (kc < k && w[kc] < w[kc - 1] * fac4)
                {
                    h = hh[kc] * a[kopt + 1] / a[kc];
                }
                else
                {
                    h = hh[kc] * a[kopt] / a[kc];
                }
            }
            k = kopt;
            h = posneg * Math.Abs(h);
            goto L10;
            /* *** *** *** *** *** *** *** */
            /* --- STEP IS REJECTED */
            /* *** *** *** *** *** *** *** */
            L100:
            k = Math.Min(Math.Min(k, kc), km - 1);
            if (k > 2 && w[k - 1] < w[k] * fac3)
            {
                --k;
            }
            ++(nrejct);
            h = posneg * hh[k];
            last = false;
            reject = true;
            if (caljac)
            {
                goto L30;
            }
            goto L10;
            /* --- SOLUTION EXIT */
            L110:
            h = hopt;
            idid = 1;
            return 0;
            /* --- FAIL EXIT */
            L120:
            //do_fio(x);
            //do_fio(h);
            idid = -1;
            return 0;
        } /* seucor_ */



        /* *** *** *** *** *** *** *** */
        /*     S U B R O U T I N E    S E U L */
        /* *** *** *** *** *** *** *** */

        /* Subroutine */
        int seul_(int jj, int n, S_fp fcn, double x,
double[] y, double[] dy, double[] fx, double[] fjac,
int lfjac, double[] fmas, int ldmas, double[] e,
int le, int[] ip, double h, int km, double
hmaxn, double[] t, double[] scal, int[] nj, double[] hh,
double[] w, double[] a, double[] yh, double[] dyh,
double[] del, double[] wh, double err, double safe1,
double fac, double fac1, double fac2, double
safe2, double theta,
ref int ndec, ref int nsol,
double errold, int[] iphes, int[] icomp, bool autnms,
bool implct, bool reject, bool atov,
double[] fsafe, int km2, int nrd, int iout, int
 ipt, int m1, int m2, int nm1, int ijob, bool
calhes)
        {
            /* System generated locals */
            double d1, d2;

            /* Builtin functions */
            //double Math.Sqrt(double), Math.Pow(double *, double *);

            /* Local variables */
            int i, j, l, m;
            double hj;
            int mm;
            double hji;
            int ier = 0;
            double sum, del1, del2, expo;
            double facmin;

            /* --- THIS SUBROUTINE COMPUTES THE J-TH LINE OF THE */
            /* --- EXTRAPOLATION TABLE AND PROVIDES AN ESTIMATE */
            /* --- OF THE OPTIMAL STEP SIZE */
            /* *** *** *** *** *** *** *** */
            /*  COMPUTE THE MATRIX E AND ITS DECOMPOSITION */
            /* *** *** *** *** *** *** *** */
            /* Parameter adjustments */
            //--iphes;
            //--wh;
            //--del;
            //--dyh;
            //--yh;
            //--scal;
            //--ip;
            //--fx;
            //--dy;
            //--y;
            //lfjac = lfjac;
            //1 + lfjac = 1 + lfjac;
            //fjac -= 1 + lfjac;
            //ldmas = ldmas;
            //1 + ldmas = 1 + ldmas;
            //fmas -= 1 + ldmas;
            //le = le;
            //1 + le = 1 + le;
            //e -= 1 + le;
            //--a;
            //--w;
            //--hh;
            //--nj;
            //km = km;
            //1 + km = 1 + km;
            //t -= 1 + km;
            //km2 = km2;
            //1 + km2 = 1 + km2;
            //fsafe -= 1 + km2;
            //--icomp;
            //--rpar;
            //--ipar;

            /* Function Body */
            hj = h / nj[jj];
            hji = 1.0 / hj;
            dc_decsol.decomr_(n, fjac, lfjac, fmas, ldmas, m1, m2, nm1, hji, e, le, ip, ref ier, ijob, calhes, iphes);
            if (ier != 0)
            {
                goto L79;
            }
            ++(ndec);
            /* *** *** *** *** *** *** *** */
            /* --- STARTING PROCEDURE */
            /* *** *** *** *** *** *** *** */
            if (!(autnms))
            {
                d1 = x + hj;
                fcn(n, d1, y, dy);
            }
            for (i = 0; i < n; ++i)
            {
                yh[i] = y[i];
                del[i] = dy[i];
            }
            dc_decsol.slvseu_(n, fjac, lfjac, fmas, ldmas, m1, m2, nm1, hji, e, le, ip, iphes, del, ijob);
            ++(nsol);
            m = nj[jj];
            if (iout == 2 && m == jj)
            {
                ++(ipt);
                for (i = 0; i < nrd; ++i)
                {
                    fsafe[ipt + i * km2] = del[icomp[i]];
                }
            }
            /* *** *** *** *** *** *** *** */
            /* --- SEMI-IMPLICIT EULER METHOD */
            /* *** *** *** *** *** *** *** */
            if (m > 1)
            {
                for (mm = 0; mm < m - 1; ++mm)
                {
                    for (i = 0; i < n; ++i)
                    {
                        yh[i] += del[i];
                    }
                    if (autnms)
                    {
                        fcn(n, x + hj * mm, yh, dyh);
                    }
                    else
                    {
                        fcn(n, x + hj * (mm + 1), yh, dyh);
                    }
                    if (mm == 1 && jj <= 2)
                    {
                        /* --- STABILITY CHECK */
                        del1 = 0.0;
                        for (i = 0; i < n; ++i)
                        {
                            /* Computing 2nd power */
                            d1 = del[i] / scal[i];
                            del1 += d1 * d1;
                        }
                        del1 = Math.Sqrt(del1);
                        if (implct)
                        {
                            for (i = 0; i < nm1; ++i)
                            {
                                wh[i] = del[i + m1];
                            }
                            //if (mlb == nm1)
                            {
                                for (i = 0; i < nm1; ++i)
                                {
                                    sum = 0.0;
                                    for (j = 0; j < nm1; ++j)
                                    {
                                        sum += fmas[i + j * ldmas] * wh[j];
                                    }
                                    del[i + m1] = sum;
                                }
                            }
                        }
                        if (!(autnms))
                        {
                            d1 = x + hj;
                            fcn(n, d1, yh, wh);
                            for (i = 0; i < n; ++i)
                            {
                                del[i] = wh[i] - del[i] * hji;
                            }
                        }
                        else
                        {
                            for (i = 0; i < n; ++i)
                            {
                                del[i] = dyh[i] - del[i] * hji;
                            }
                        }
                        dc_decsol.slvseu_(n, fjac, lfjac, fmas, ldmas, m1, m2, nm1, hji, e, le, ip, iphes, del, ijob);
                        ++(nsol);
                        del2 = 0.0;
                        for (i = 0; i < n; ++i)
                        {
                            /* Computing 2nd power */
                            d1 = del[i] / scal[i];
                            del2 += d1 * d1;
                        }
                        del2 = Math.Sqrt(del2);
                        theta = del2 / Math.Max(1.0, del1);
                        if (theta > 1.0)
                        {
                            goto L79;
                        }
                    }
                    dc_decsol.slvseu_(n, fjac, lfjac, fmas, ldmas, m1, m2, nm1, hji, e, le, ip, iphes, dyh, ijob);
                    ++(nsol);
                    for (i = 0; i < n; ++i)
                    {
                        del[i] = dyh[i];
                    }
                    if (iout == 2 && mm >= m - jj)
                    {
                        ++(ipt);
                        for (i = 0; i < nrd; ++i)
                        {
                            fsafe[ipt + i * km2] = del[icomp[i]];
                        }
                    }
                }
            }
            for (i = 0; i < n; ++i)
            {
                t[jj + i * km] = yh[i] + del[i];
            }
            /* *** *** *** *** *** *** *** */
            /* --- POLYNOMIAL EXTRAPOLATION */
            /* *** *** *** *** *** *** *** */
            if (jj == 1)
            {
                return 0;
            }
            for (l = jj; l >= 2; --l)
            {
                fac = nj[jj] / (double)nj[l - 1] - 1.0;
                for (i = 0; i < n; ++i)
                {
                    t[l - 1 + i * km] = t[l + i * km] + (t[l + i * km] - t[l - 1 + i * km]) / fac;
                }
            }
            err = 0.0;
            for (i = 0; i < n; ++i)
            {
                /* Computing 2nd power */
                d2 = Math.Min(Math.Abs(t[i * km + 1] - t[i * km + 2]) / scal[i], 1e15);
                err += d2 * d2;
            }
            if (err >= 1e30)
            {
                goto L79;
            }
            err = Math.Sqrt(err / (double)(n));
            if (jj > 2 && err >= errold)
            {
                goto L79;
            }
            /* Computing MAX */
            d1 = err * 4;
            errold = Math.Max(d1, 1.0);
            /* --- COMPUTE OPTIMAL STEP SIZES */
            expo = 1.0 / jj;
            facmin = Math.Pow(fac1, expo);
            fac = Math.Min(fac2 / facmin, Math.Max(facmin, Math.Pow(err / safe1, expo) / safe2));
            fac = 1.0 / fac;

            hh[jj] = Math.Min(Math.Abs(h) * fac, hmaxn);
            w[jj] = a[jj] / hh[jj];
            return 0;
            L79:
            atov = true;
            h *= 0.5;
            reject = true;
            return 0;
        } /* seul_ */


        double contex_(int ii, double x, double[] rc, int lrc,
            int[] ic, int lic)
        {
            /* System generated locals */
            double ret_val = 0;

            /* Local variables */
            int i, j;
            double theta;

            /* ---------------------------------------------------------- */
            /*     THIS FUNCTION CAN BE USED FOR CONINUOUS OUTPUT IN CONECTION */
            /*     WITH THE OUTPUT-SUBROUTINE FOR SEULEX. IT PROVIDES AN */
            /*     APPROXIMATION TO THE II-TH COMPONENT OF THE SOLUTION AT X. */
            /* ---------------------------------------------------------- */
            /* ----- COMPUTE PLACE OF II-TH COMPONENT */
            /* Parameter adjustments */
            //--rc;
            //--ic;

            /* Function Body */
            i = 0;
            for (j = 0; j < coseu_2.nrd; ++j)
            {
                if (ic[j] == ii)
                {
                    i = j;
                }
            }
            if (i == 0)
            {
                Console.WriteLine(" NO DENSE OUTPUT AVAILABLE FOR COMP.", ii);
                return ret_val;
            }
            /* ----- COMPUTE THE INTERPOLATED VALUE */
            theta = (x - coseu_2.xold) / coseu_2.h;
            ret_val = rc[coseu_2.ir * coseu_2.nrd + i];
            for (j = 1; j < coseu_2.ir; ++j)
            {
                ret_val = rc[(coseu_2.ir + 1 - j) * coseu_2.nrd + i] + ret_val * (theta - 1.0);
            }
            ret_val = rc[i] + ret_val * theta;
            return ret_val;
        } /* contex_ */
    }
}
