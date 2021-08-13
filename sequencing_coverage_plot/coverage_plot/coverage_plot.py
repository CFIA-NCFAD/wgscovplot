from .base_plot import BasePlot
from typing import Sequence, Optional
import pandas as pd


class CoveragePlot(BasePlot):

    def prepare_dataset(self, df_metadata: pd.DataFrame, header: list = []):
        self.options.get("dataset").append(
            {
                "dimensions": header,  # header of table
                "source": df_metadata.values.tolist()
            }
        )
        return self

    def add_xaxis(self, grid_index: int = 0):
        if Sequence is not None:
            self.options.get("xAxis").append(
                {
                    "type": "category",
                    "gridIndex": grid_index,
                    "name": "Pos",
                    "axisLabel": {
                        "interval": "auto"
                    }
                }
            )
        return self

    def add_yaxis(
            self,
            series_name: str,
            scale: str = "log",
            is_large: bool = True,
            grid_index: int = 0,
            x_axis_index: int = 0,
            y_axis_index: int = 0,
            x_data_name: str = "",
            y_data_name: str = "",
            dataset_index: Optional[int] = 0,
            tooltip_data: list = [],
            range_min: Optional[int] = 10,
            range_max: Optional[int] = 100000
    ):
        self.options["title"].update(text="Coverage Plot")
        self.options.get("yAxis").append(
            {
                "type": scale,
                "gridIndex": grid_index,
                "name": series_name,
                "nameTextStyle": {
                    "fontStyle": "normal",
                    "fontWeight": "bolder",
                },
                "nameLocation": "end",
                "min": range_min,
                "max": range_max
            }
        )
        self.options.get("series").append(
            {
                "type": "bar",
                "xAxisIndex": x_axis_index,
                "yAxisIndex": y_axis_index,
                "encode": {
                    "x": x_data_name,
                    "y": y_data_name,
                    "tooltip": tooltip_data
                },
                "datasetIndex": dataset_index,
                "large": is_large,
                "itemStyle": {
                    "color": '#1e6793'
                }
            }
        )
        return self

    def data_zoom(
            self,
            zoom_style: Optional[str] = "slider",
            orient: Optional[str] = "horizontal",
            is_show: bool = True,
            x_axis_index: list = [],
            y_axis_index: list = []
    ):
        self.options.get("dataZoom")[0].update(xAxisIndex=x_axis_index,
                                               yAxisIndex=y_axis_index)
        if zoom_style == "slider":
            self.options.get("dataZoom")[1].update(show=is_show,
                                                   xAxisIndex=x_axis_index,
                                                   yAxisIndex=y_axis_index,
                                                   orient=orient)
        return self

    def add_grid(
            self,
            pos_top: str,
            pos_bottom: str,
            height: str = "auto"
    ):
        self.options.get("grid").append(
            {
                "top": pos_top,
                "bottom": pos_bottom,
                "height": height
            }
        )
        return self
