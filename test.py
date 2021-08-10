from sequencing_coverage_plot import CoveragePlot
import pandas as pd


def read_depths(fpath) -> pd.DataFrame:
    df = pd.read_table(fpath,
                       names=['reference', 'pos', 'depth'],
                       header=None)
    return df


cov_file1 = '/home/hnguyen/Documents/LabData/results/mapping/sample-01.fastq/sample-01.fastq-MN908947.3-depths.tsv'
df_depth = read_depths(cov_file1)

x = df_depth['pos'].to_list()
y = df_depth['depth'].to_list()

cov_plot = CoveragePlot()
cov_plot.add_xaxis(x)
cov_plot.add_yaxis("Sample_01_CFIA", y)
cov_plot.render_html("first_coverage_plot.html")
