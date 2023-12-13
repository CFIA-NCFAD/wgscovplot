import {WgsCovPlotDB} from "../db";
import {getDepthSeries} from "./depthSeries";
import {isEmpty, isNil, some} from "lodash";
import {getGeneFeatureSeries, getRegionAmpliconDepthSeries} from "./geneFeatures";
import {getVariantsSeries} from "./variantSeries";
import {getSegmentVariantSeries} from "./segmented/variantSeries";
import {getSegmentPrimerSeries} from "./segmented/pcr";
import {setState} from "../state";

export const getSeries = (db: WgsCovPlotDB) => {
  let series: any[] = [];
  series.push(...getDepthSeries(db));
  if (!isEmpty(db.amplicon_depths)) {
    series.push(...getRegionAmpliconDepthSeries(db));
  }
  if (!isEmpty(db.variants) && db.chartOptions.showVariants) {
    if (isNil(db.segments)) {
      series.push(...getVariantsSeries(db));
    } else {
      series.push(...getSegmentVariantSeries(db));
    }
  }
  if (!isEmpty(db.primer_matches) && db.show_primer_matches) {
    series.push(...getSegmentPrimerSeries(db));
  }
  let isNegativeStrand: boolean = some(db.echart_features, {value: {strand: -1}})
  let isPositiveStrand: boolean = some(db.echart_features, {value: {strand: 1}})
  if (isNegativeStrand && isPositiveStrand){
      setState("doubleStrand", true)
  }
  if (db.chartOptions.showFeatures && (db.show_genes || db.show_amplicons)) {
    let geneFeatureSeries: any = getGeneFeatureSeries(db);
    series.push(geneFeatureSeries);
  }
  return series;
}