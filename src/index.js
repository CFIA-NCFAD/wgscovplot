import "bootstrap/dist/css/bootstrap.css";
import "bootstrap/dist/js/bootstrap.js";
import $ from "jquery";
import JQuery from "jquery";
import {popper} from "@popperjs/core";
import select2 from "select2";
import "select2/dist/css/select2.css";

import {echarts, initEventHandlers, initWgscovplotRenderEnv} from "./wgscovplot";

export {
    echarts,
    initEventHandlers,
    initWgscovplotRenderEnv,
    $, JQuery,
    select2,
    popper,
};
