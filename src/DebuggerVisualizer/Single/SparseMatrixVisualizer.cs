using MathNet.MatrixDebuggerVisualizer.UI;
using MathNet.Numerics.LinearAlgebra.Single;
using Microsoft.VisualStudio.DebuggerVisualizers;
using System;

namespace MathNet.MatrixDebuggerVisualizer.Single
{
	public class SparseMatrixVisualizer : DialogDebuggerVisualizer
	{
		protected override void Show(IDialogVisualizerService windowService, IVisualizerObjectProvider objectProvider)
		{
			if (windowService == null)
			{
				throw new ArgumentNullException("windowService");
			}

			if (objectProvider == null)
			{
				throw new ArgumentNullException("objectProvider");
			}

			var obj = (SparseMatrix)objectProvider.GetObject();

			using (var form = new SparseMatrixVisualizerForm())
			{
                form.SetStorageAdapter(new SparseStorageAdapter(obj));
				windowService.ShowDialog(form);
			}
		}
	}
}
