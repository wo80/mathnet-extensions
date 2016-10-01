using MathNet.MatrixDebuggerVisualizer.UI;
using MathNet.Numerics.LinearAlgebra.Complex32;
using Microsoft.VisualStudio.DebuggerVisualizers;
using System;

namespace MathNet.MatrixDebuggerVisualizer.Complex32
{
	public class SparseMatrixVisualizer : DialogDebuggerVisualizer
	{
		protected override void Show(IDialogVisualizerService windowService, IVisualizerObjectProvider objectProvider)
		{
			if (windowService == null)
			{
				throw new System.ArgumentNullException("windowService");
			}

			if (objectProvider == null)
			{
				throw new System.ArgumentNullException("objectProvider");
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
