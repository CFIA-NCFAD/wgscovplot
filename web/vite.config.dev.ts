import {defineConfig} from 'vite';
import solidPlugin from 'vite-plugin-solid';
import handlebars from 'vite-plugin-handlebars';


export default defineConfig({
    plugins: [
        solidPlugin(),
        // @ts-ignore
        handlebars({
            context: {
            }
        }),
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
                chunkFileNames: `assets/[name].js`,
                assetFileNames: `assets/[name].[ext]`
            }
        },
        watch: {}
    },
});