import {Component} from "solid-js";

export const AboutPage: Component = () => {
  return (
    <div class="w-full">
      <div class="container mx-auto px-4">
        <div class="flex flex-wrap">
          <div class="w-full lg:w-12/12 px-4">
            <div class="relative flex flex-col min-w-0 break-words bg-white w-full mb-6 shadow-lg rounded">
              <div class="px-4 py-5 flex-auto">
                <h4 class="text-2xl font-semibold">About wgscovplot</h4>
                <p class="mt-2 mb-4 text-gray-600">
                  wgscovplot is a tool for visualising coverage data from WGS sequencing experiments.
                  It is best suited for viral sequencing data, but could be used for other microbial organisms.
                </p>
                <p class="mt-2 mb-4 text-gray-600">
                  More information can be found on
                  the <a class="text-blue-800 hover:text-blue-500 hover:ring hover:rounded"
                         target="_blank" rel="noopener noreferrer"
                         href="https://github.com/CFIA-NCFAD/wgscovplot">
                  GitHub repository
                </a>.
                </p>
              </div>
            </div>
          </div>
        </div>
      </div>
    </div>
  );
}