import { source } from '@/lib/source';
import { DocsBody, DocsDescription, DocsPage, DocsTitle } from 'fumadocs-ui/layouts/docs/page';
import { notFound } from 'next/navigation';
import { getMDXComponents } from '@/components/mdx';
import type { Metadata } from 'next';
export function generateStaticParams() { return source.generateParams(); }
export default async function Page(props: PageProps<'/docs/[[...slug]]'>) { const { slug } = await props.params; const page = source.getPage(slug); if (!page) notFound(); const MDX = page.data.body; return <DocsPage toc={page.data.toc} full={page.data.full}><DocsTitle>{page.data.title}</DocsTitle><DocsDescription>{page.data.description}</DocsDescription><DocsBody><MDX components={getMDXComponents()} /></DocsBody></DocsPage>; }
export async function generateMetadata(props: PageProps<'/docs/[[...slug]]'>): Promise<Metadata> { const { slug } = await props.params; const page = source.getPage(slug); if (!page) notFound(); return { title: page.data.title, description: page.data.description }; }