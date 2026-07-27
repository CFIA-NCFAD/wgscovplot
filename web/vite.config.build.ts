import {defineConfig} from 'vite';
import solidPlugin from 'vite-plugin-solid';
import cssInjectedByJsPlugin from "vite-plugin-css-injected-by-js";


export default defineConfig({
    plugins: [
        solidPlugin(),
        cssInjectedByJsPlugin({topExecutionPriority: false}),
    ],
    define: {
        'process.env.NODE_ENV': JSON.stringify('production'),
        'process.env': {}
    },
    build: {
        lib: {
            entry: './src/index.tsx',
            name: 'wgscovplot',
            fileName: 'wgscovplot'
        },
        outDir: './build',
        rollupOptions: {
            output: {
                entryFileNames: `wgscovplot.js`,
            }
        },
        watch: null
    },
});