import { createMDX } from 'fumadocs-mdx/next';

const withMDX = createMDX();

/** @type {import('next').NextConfig} */
const config = {
  output: 'export',
  basePath: process.env.GITHUB_ACTIONS ? '/molframe' : '',
  assetPrefix: process.env.GITHUB_ACTIONS ? '/molframe/' : undefined,
  trailingSlash: true,
  reactStrictMode: true,
};

export default withMDX(config);