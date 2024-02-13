import {createStore} from "solid-js/store";
import {defaultDB, WgsCovPlotDB} from "./db";

// global state DB taking values from window.db if available
// @ts-expect-error db should exist in window
export const [state, setState] = createStore<WgsCovPlotDB>({...defaultDB, ...db} ?? defaultDB);
