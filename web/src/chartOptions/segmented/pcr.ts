import {graphic} from "echarts/core";
import {get, has, isEmpty, isNil} from "lodash";
import {WgsCovPlotDB} from "../../db";

export function getSegmentPrimerData(sample: string, db: WgsCovPlotDB) {
  if (isNil(db.segCoords)) {
    return;
  }
  const primerFeatures = [];
  const segments = Object.keys(db.segCoords); // trigger only segCoords is updated from get Dataset
  for (const segment of segments) {
    let matchesData;
    if (has(db.primer_matches, [sample, segment])) {
      matchesData = get(db.primer_matches, [sample, segment]);
    }
    const coords = get(db.segCoords, [segment]);
    if (!isNil(matchesData)) {
      for (const match of matchesData) {
        const {
          start,
          end,
          name,
          query_aligned,
          target_aligned,
          matched_aligned,
          cigar,
          edit_distance,
          other_locations,
        } = match;
        // TODO: push object instead of array? does echarts need value to be an array instead of object?
        primerFeatures.push({
          value: [
            start + coords.start,
            end + coords.start,
            10, // height for primer annotation in log scale
            name,
            query_aligned,
            target_aligned,
            matched_aligned,
            cigar,
            edit_distance,
            other_locations
          ],
          itemStyle: {color: "#b71ae3"}
        });
      }
    }
  }
  return primerFeatures;
}

export function primerMatchRenderer(db: WgsCovPlotDB) {
  // eslint-disable-next-line
  return function ({coordSys} : {coordSys: any}, api: any) {
    const [startX, startY] = api.coord([api.value(0), api.value(2)]);
    const [endX, endY] = api.coord([api.value(1), 1]);
    const rectShape = graphic.clipRectByRect(
      {
        x: startX,
        y: startY,
        width: endX - startX,
        height: (endY - startY) * 0.3
      },
      coordSys
    );
    return rectShape && {
      type: "rect",
      shape: rectShape,
      style: api.style(),
      invisible: !db.show_primer_matches
    };
  };
}

export function getSegmentPrimerSeries(db: WgsCovPlotDB) {
  // eslint-disable-next-line
  const primerSeries: any[] = [];
  if (isEmpty(db.primer_matches)) {
    return primerSeries;
  }
  for (const [i, sample] of db.chartOptions.selectedSamples.entries()) {
    primerSeries.push({
      type: "custom",
      xAxisIndex: i,
      yAxisIndex: i,
      renderItem: primerMatchRenderer(db),
      label: {
        show: db.show_primer_matches,
        position: "top",
        distance: 20,
        rotate: 45
      },
      encode: {
        x: [0, 1],
        y: 3,
      },
      silent: true,
      data: getSegmentPrimerData(sample, db),
    });
  }
  return primerSeries;
}

export function getPrimerInfo(sample: string, position: number, segment: string, db: WgsCovPlotDB) {
  const primerInfos = get(db.primer_matches, [sample, segment], []);
  let primerInfoRows: string[][] = [];
  for (const primerInfo of primerInfos) {
    const {
      start,
      end,
      name,
      query_aligned,
      matched_aligned,
      target_aligned,
      edit_distance,
      cigar,
      other_locations,
    } = primerInfo;
    if (position >= start && position <= end) {
      primerInfoRows = [
        ["Primer Name", name],
        ["Primer Sequence", query_aligned],
        ["Match Aligned", matched_aligned],
        ["Ref Sequence", target_aligned],
        ["Cigar", cigar],
        ["Start", (start + 1).toLocaleString()],
        ["End", (end + 1).toLocaleString()],
        ["Edit Distance", edit_distance.toString()],
        ["Other Locations", other_locations]
      ];
    }
  }
  return primerInfoRows;
}