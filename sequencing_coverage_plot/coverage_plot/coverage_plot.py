from .base_plot import BasePlot
from typing import Sequence, Union, Optional
import pandas as pd


class CoveragePlot(BasePlot):

    def prepare_dataset(self, df_metadata: pd.DataFrame):
        columns = ["SampleName", "Reference", "Pos", "CoverageDepth"]
        source = [columns]
        for data in df_metadata.values.tolist():
            source.append(data)  # source : [[],[],...,[]]
        self.options.get("dataset").append(
            {
                "source": source
            }
        )

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
            is_large: bool = True,
            grid_index: int = 0,
            x_axis_index: int = 0,
            y_axis_index: int = 0,
            dataset_index: Optional[int] = 0
    ):
        self.options["title"].update(text="Coverage Plot")
        self.options.get("yAxis").append(
            {
                "type": "log",
                "gridIndex": grid_index,
                "name": series_name,
                "nameTextStyle": {
                    "fontStyle": "normal",
                    "fontWeight": "bolder",
                },
                "nameLocation": "end",
                "min": 10,
                "max": 100000
            }
        )
        self.options.get("series").append(
            {
                "type": "bar",
                "xAxisIndex": x_axis_index,
                "yAxisIndex": y_axis_index,
                "encode": {
                    "x": "Pos",
                    "y": "CoverageDepth",
                    "tooltip": [0, 1, 2, 3]
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
        data = self.options.get("dataZoom")
        data[0]["xAxisIndex"] = x_axis_index
        data[0]["yAxisIndex"] = y_axis_index
        if zoom_style == "slider":
            data[1]["show"] = is_show
            data[1]["xAxisIndex"] = x_axis_index
            data[1]["yAxisIndex"] = y_axis_index
            data[1]["orient"] = orient
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
