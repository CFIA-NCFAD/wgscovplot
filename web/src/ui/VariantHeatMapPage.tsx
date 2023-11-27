import {Component} from "solid-js";

export const VariantHeatMapPage: Component = () => {
  return (
    <div class="w-full">
      <div class="container mx-auto px-4">
        <div class="flex flex-wrap">
          <div class="w-full lg:w-12/12 px-4">
            <div class="relative flex flex-col min-w-0 break-words bg-white w-full mb-6 shadow-lg rounded">
              <div class="px-4 py-5 flex-auto">
                <h4 class="text-2xl font-semibold">Variant Heatmap</h4>
              </div>
            </div>
          </div>
        </div>
      </div>
    </div>
  );
}