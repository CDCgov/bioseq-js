import { defineConfig } from "tsup";

export default defineConfig([
  // npm: ESM + CommonJS + type declarations
  {
    entry: { index: "src/index.ts" },
    format: ["esm", "cjs"],
    dts: true,
    sourcemap: true,
    clean: true,
    target: "es2020",
  },
  // <script> tag builds exposing a global `bioseq`
  {
    entry: { bioseq: "src/index.ts" },
    format: ["iife"],
    globalName: "bioseq",
    target: "es2020",
    outExtension: () => ({ js: ".js" }),
  },
  {
    entry: { "bioseq.min": "src/index.ts" },
    format: ["iife"],
    globalName: "bioseq",
    minify: true,
    target: "es2020",
    outExtension: () => ({ js: ".js" }),
  },
]);
