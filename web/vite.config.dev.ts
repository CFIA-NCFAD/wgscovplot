import {defineConfig} from 'vite';
import solidPlugin from 'vite-plugin-solid';
import handlebars from 'vite-plugin-handlebars';
import cssInjectedByJsPlugin from "vite-plugin-css-injected-by-js";


export default defineConfig({
    plugins: [
        solidPlugin(),
        cssInjectedByJsPlugin({topExecutionPriority: false}),
    ],
    define: {
        'process.env': {}
    },
    build: {
        lib: {
            entry: './src/index.tsx',
            name: 'wgscovplot',
            fileName: 'wgscovplot'
        },
        outDir: './tmpl',
        rollupOptions: {
            output: {
                entryFileNames: `assets/[name].js`,
            }
        },
        watch: {}
    },
});