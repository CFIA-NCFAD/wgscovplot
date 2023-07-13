import {get, isNil, find, uniqBy, flatMap, keys} from "lodash";
import {genomeCoverage, meanCoverage, medianCoverage} from "../stats";
import {SampleDepths, VariantCall, WgsCovPlotDB} from "../db";
import {unwrap} from "solid-js/store";

/**
 * Variant calling result comparison table at a given position
 */
export const getVariantComparison = (db: WgsCovPlotDB, position: number): Array<Array<string>> => {
    if (isNil(db.variants)) {
        return [];
    }
    let rows: any[][] = [];
    let variantArr: VariantCall[] = [];
    for (let sample of db.chartOptions.selectedSamples) {
        let sampleVariants = db.variants[sample];
        let variant = find(sampleVariants, {POS: position + ""}) as VariantCall;
        if (!isNil(variant)) {
            variantArr.push(variant);
        } else {
            variantArr.push({sample, POS: position + ""});
        }
    }
    // @ts-ignore
    let unionKeys = uniqBy(flatMap(variantArr, keys));
    let depths: SampleDepths = unwrap(db.depths) as SampleDepths;
    unionKeys.push("Coverage Depth");
    unionKeys.forEach((key: string) => {
        let row = [];
        row.push(key);
        if (key === "Coverage Depth") {
            for (let sample of db.chartOptions.selectedSamples) {
                row.push(depths[sample][position - 1].toLocaleString());
            }
        } else {
            for (let variant of variantArr) {
                let val = get(variant, key, "");
                row.push(val);
            }
        }
        rows.push(...[row]);
    });
    return rows;
}

export const getCoverageStatComparison = (db: WgsCovPlotDB, start: number, end: number, position: number) => {
    const selectedSamples = db.chartOptions.selectedSamples;
    const depths = unwrap(db.depths) as SampleDepths;
    const low_coverage_threshold = db.chartOptions.low_coverage_threshold;
    let rows = [];
    let tableHeader = [
        "",
        `Depth at position ${position.toLocaleString()}`,
        "Range",
        "Mean Coverage (X)",
        "Median Coverage (X)",
        `Genome Coverage (>=${low_coverage_threshold}X) (%)`,
    ];

    for (let sample of selectedSamples) {
        let sampleDepths = depths[sample];
        let meanCov = meanCoverage(sampleDepths, start, end).toFixed(2);
        let medianCov = medianCoverage(sampleDepths, start, end).toFixed(2);
        let genomeCov = genomeCoverage(sampleDepths, start, end, low_coverage_threshold).toFixed(2);
        let coverageDepth = sampleDepths[position - 1];
        let row = [
            sample,
            coverageDepth.toLocaleString(),
            `${start.toLocaleString()}-${end.toLocaleString()}`,
            meanCov,
            medianCov,
            genomeCov,
        ];
        rows.push(...[row]);
    }
    return {headers: tableHeader, rows: rows};
}
