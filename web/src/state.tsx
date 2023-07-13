import {createStore} from "solid-js/store";
import {defaultDB, WgsCovPlotDB} from "./db";

// global state DB taking values from window.db if available
// @ts-ignore
export const [state, setState] = createStore<WgsCovPlotDB>({...defaultDB, ...window.db} ?? defaultDB);
