from .base_plot import BasePlot
from typing import Sequence


class CoveragePlot(BasePlot):

    def add_xaxis(self, xaxis_data: Sequence):
        self.options["xAxis"][0].update(data=xaxis_data)
        return self

    def add_yaxis(
            self,
            series_name: str,
            y_axis: Sequence,
            is_large: bool = False,
    ):
        self.options.get("series").append(
            {
                "type": "bar",
                "name": series_name,
                "data": y_axis,
                "large": is_large,
            }
        )
        return self
