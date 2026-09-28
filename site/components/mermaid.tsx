'use client';

import { useEffect, useId, useState } from 'react';
import { useTheme } from 'next-themes';

export function Mermaid({ chart }: { chart: string }) {
  const id = useId().replace(/[^a-zA-Z0-9_-]/g, '');
  const { resolvedTheme } = useTheme();
  const [svg, setSvg] = useState<string>();
  const [error, setError] = useState<string>();

  useEffect(() => {
    let cancelled = false;

    void import('mermaid').then(async ({ default: mermaid }) => {
      mermaid.initialize({
        startOnLoad: false,
        securityLevel: 'strict',
        theme: resolvedTheme === 'dark' ? 'dark' : 'default',
      });
      const rendered = await mermaid.render(`molframe-mermaid-${id}`, chart);
      if (!cancelled) setSvg(rendered.svg);
    }).catch((cause: unknown) => {
      if (!cancelled) {
        setError(cause instanceof Error ? cause.message : 'Unable to render diagram.');
      }
    });

    return () => { cancelled = true; };
  }, [chart, id, resolvedTheme]);

  if (error) return <pre className="not-prose overflow-x-auto rounded-lg border p-4 text-sm">{chart}</pre>;
  if (!svg) return <div aria-busy="true" className="my-6 min-h-24 rounded-lg border bg-fd-muted/30" />;

  return <div className="not-prose my-6 overflow-x-auto rounded-lg border bg-fd-card p-4" dangerouslySetInnerHTML={{ __html: svg }} />;
}
