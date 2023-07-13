import { defineConfig } from 'vite';
import solidPlugin from 'vite-plugin-solid';
import handlebars from 'vite-plugin-handlebars';

export default defineConfig({
  plugins: [
    solidPlugin(),
    // @ts-ignore
    handlebars({
      context: {
        db_js: 'data/flu-db.js',
      }
    }),
  ],
  server: {
    port: 3000,
  },
  build: {
    target: 'esnext',
  },
});
