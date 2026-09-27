import { RootProvider } from 'fumadocs-ui/provider/next';
import './global.css';
import type { Metadata } from 'next';

export const metadata: Metadata = {
  title: { default: 'MolFrame documentation', template: '%s | MolFrame' },
  description: 'Practical structural bioinformatics for Python and Rust.',
  metadataBase: new URL('https://miguelcsx.github.io/molframe/'),
};

export default function Layout({ children }: LayoutProps<'/'>) {
  return <html lang="en" suppressHydrationWarning><body className="min-h-screen"><RootProvider search={{ options: { type: 'static' } }}>{children}</RootProvider></body></html>;
}