# MolFrame documentation site

Fumadocs + Next.js static export for GitHub Pages.

```bash
nix shell nixpkgs#nodejs_22 -c npm ci
nix shell nixpkgs#nodejs_22 -c npm run dev
```

`npm run build` writes the deployable static artifact to `out/`. GitHub Actions sets the `/molframe` base path for the project page.
