// From BenchmarkDotNet, https://github.com/PerfDotNet/BenchmarkDotNet/
// Copyright (c) 2013–2015 Andrey Akinshin, Jon Skeet, Matt Warren
// MIT License

namespace ExtensionsBenchmark
{
    using System;
    using System.Diagnostics;
    using System.Globalization;
    using System.Linq;
    using System.Management;
    using System.Reflection;
    using System.Text;

    public sealed class EnvironmentHelper
    {
        static EnvironmentHelper()
        {
            MainCultureInfo = (CultureInfo)CultureInfo.InvariantCulture.Clone();
            MainCultureInfo.NumberFormat.NumberDecimalSeparator = ".";
        }

        public static readonly CultureInfo MainCultureInfo;

        public string OsVersion { get; set; }
        public string ProcessorName { get; set; }
        public int ProcessorCount { get; set; }
        public string ClrVersion { get; set; }
        public string Architecture { get; set; }
        public bool HasAttachedDebugger { get; set; }
        public string Configuration { get; set; }

        /// <summary>
        /// The frequency of the timer as the number of ticks per second.
        /// </summary>
        public long ChronometerFrequency { get; set; }

        public static EnvironmentHelper GetCurrentInfo()
        {
            return new EnvironmentHelper
        {
            OsVersion = GetOsVersion(),
            ProcessorName = GetProcessorName(),
            ProcessorCount = GetProcessorCount(),
            ClrVersion = GetClrVersion(),
            Architecture = GetArchitecture(),
            HasAttachedDebugger = GetHasAttachedDebugger(),
            Configuration = GetConfiguration()
        };
        }

        public string ToFormattedString(string clrHint = "")
        {
            var sb = new StringBuilder();

            sb.AppendFormat("OS={0}", OsVersion);
            sb.AppendLine();
            sb.AppendFormat("Processor={0}, ProcessorCount={1}", ProcessorName, ProcessorCount);
            sb.AppendLine();
            sb.AppendFormat("{0}CLR={1}, Arch={2} {3}{4}", clrHint, ClrVersion, Architecture, Configuration, GetDebuggerFlag());
            sb.AppendLine();

            return sb.ToString();
        }

        private string GetDebuggerFlag()
        {
            return HasAttachedDebugger ? " [AttachedDebugger]" : "";
        }

        private static string GetOsVersion()
        {
            return Environment.OSVersion.ToString();
        }

        private static string GetProcessorName()
        {
            string info = string.Empty;
            if (IsWindows() && !IsMono())
            {
                try
                {
                    var mosProcessor = new ManagementObjectSearcher("SELECT * FROM Win32_Processor");
                    foreach (var moProcessor in mosProcessor.Get().Cast<ManagementObject>())
                    {
                        var name = moProcessor["name"];
                        if (name != null)
                        {
                            info += name.ToString();
                        }
                    }
                }
                catch (Exception)
                {
                }
            }
            else
                info = "?";
            return info;
        }

        private static int GetProcessorCount()
        {
            return Environment.ProcessorCount;
        }

        private static string GetClrVersion()
        {
            if (IsMono())
            {
                var monoRuntimeType = Type.GetType("Mono.Runtime");
                if (monoRuntimeType != null)
                {
                    var monoDisplayName = monoRuntimeType.GetMethod("GetDisplayName", BindingFlags.NonPublic | BindingFlags.Static);
                    if (monoDisplayName != null)
                        return "Mono " + monoDisplayName.Invoke(null, null);
                }
            }
            return "MS.NET " + Environment.Version;
        }

        private static string GetArchitecture()
        {
            return IntPtr.Size == 4 ? "32-bit" : "64-bit";
        }

        private static bool GetHasAttachedDebugger()
        {
            return Debugger.IsAttached;
        }

        private static string GetConfiguration()
        {
            string configuration = "RELEASE";
#if DEBUG
            configuration = "DEBUG";
#endif
            return configuration;
        }

        public static bool IsMono()
        {
            return Type.GetType("Mono.Runtime") != null;
        }

        public static bool IsWindows()
        {
            return Environment.OSVersion.Platform.ToString().Contains("Win");
        }
    }
}