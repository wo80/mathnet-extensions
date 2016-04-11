using System;
using System.Collections.Generic;
using System.Text;

namespace ExtensionsBenchmark.Benchmark
{
    public interface IReporter
    {
        void Write(StringBuilder sb, EnvironmentHelper env);
        void Write(StringBuilder sb, IEnumerable<BenchmarkContext> results);
    }

    public class MarkdownReporter : IReporter
    {
        public void Write(StringBuilder sb, EnvironmentHelper env)
        {
            sb.AppendLine(env.ToFormattedString());
        }

        public void Write(StringBuilder sb, IEnumerable<BenchmarkContext> results)
        {
            sb.AppendLine("|      Method      | Math.NET | Extension |     |");
            sb.AppendLine("| ---------------- |---------:| ---------:| --- |");

            foreach (var context in results)
            {
                var time = context.TotalTime();

                sb.AppendFormat("| {0,40} ", context.Name);
                sb.AppendFormat("| {0,10} ms ", time.Item1);
                sb.AppendFormat("| {0,10} ms ", time.Item2);
                sb.AppendFormat("| {0,-5} ", Report.CompareResult((double)time.Item1 / time.Item2));
                sb.AppendLine("|");
            }
        }
    }

    public class Report
    {
        IReporter reporter;
        StringBuilder sb;

        public Report()
            : this(new MarkdownReporter())
        {
        }

        public Report(IReporter reporter)
        {
            this.reporter = reporter;
            this.sb = new StringBuilder();
        }

        public void Append(EnvironmentHelper env)
        {
            reporter.Write(sb, env);
        }

        public void Append(IEnumerable<BenchmarkContext> results)
        {
            reporter.Write(sb, results);
        }

        public override string ToString()
        {
            return sb.ToString();
        }

        public static string CompareResult(double ratio)
        {
            if (ratio > 1)
            {
                return GetCompareResultString(ratio, '+');
            }

            if (ratio < 1)
            {
                return GetCompareResultString(1d / ratio, '-');
            }

            return string.Empty;
        }

        private static string GetCompareResultString(double ratio, char symbol)
        {
            if (ratio > 1000)
            {
                return new String(symbol, 5);
            }

            if (ratio > 500)
            {
                return new String(symbol, 4);
            }

            if (ratio > 100)
            {
                return new String(symbol, 3);
            }

            if (ratio > 10)
            {
                return new String(symbol, 2);
            }

            if (ratio > 1.5)
            {
                return new String(symbol, 1);
            }

            return string.Empty;
        }
    }
}
