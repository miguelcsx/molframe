import type { BaseLayoutProps } from 'fumadocs-ui/layouts/shared';

export function baseOptions(): BaseLayoutProps {
  return {
    nav: { title: 'MolFrame', url: '/docs' },
    githubUrl: 'https://github.com/miguelcsx/molframe',
    links: [{ text: 'MolGFX', url: 'https://miguelcsx.github.io/molgfx/', type: 'main' }],
  };
}