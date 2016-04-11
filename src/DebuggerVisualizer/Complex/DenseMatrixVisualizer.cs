using MathNet.MatrixDebuggerVisualizer.UI;
using MathNet.Numerics.LinearAlgebra.Complex;
using Microsoft.VisualStudio.DebuggerVisualizers;
using System;

namespace MathNet.MatrixDebuggerVisualizer.Complex
{
	public class DenseMatrixVisualizer : DialogDebuggerVisualizer
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

			var obj = (DenseMatrix)objectProvider.GetObject();

			using (var form = new DenseMatrixVisualizerForm())
            {
                form.SetStorageAdapter(new DenseStorageAdapter(obj));
				windowService.ShowDialog(form);
			}
		}
	}
}
