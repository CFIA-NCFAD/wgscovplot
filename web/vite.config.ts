import { defineConfig } from 'vite';
import solidPlugin from 'vite-plugin-solid';
import handlebars from 'vite-plugin-handlebars';


export default defineConfig({
  plugins: [
    solidPlugin(),
    // @ts-ignore
    handlebars({
      context: {
        db_js: 'data/sars-cov-2.js',
      }
    }),
  ],
  server: {
    port: 3000,
    open: true
  },
  build: {
    target: 'esnext',
  },
});
