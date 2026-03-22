#!/usr/bin/env node

import fs from 'node:fs/promises';
import path from 'node:path';

const USAGE = [
  'Usage: node scripts/molstar_render_single.mjs --structure PATH --format mmcif|pdb|bcif --out-dir DIR [options]',
  '',
  'Options:',
  '  --width N         Output image width (default: 320)',
  '  --height N        Output image height (default: 240)',
  '  --style-json JSON Optional style payload passed from Python queue',
  '  -h, --help        Show this help text',
].join('\n');

function parseArgs(argv) {
  if (argv.includes('--help') || argv.includes('-h')) {
    return { help: true };
  }

  const result = {
    help: false,
    structure: null,
    format: 'mmcif',
    outDir: null,
    width: 320,
    height: 240,
    styleJson: '{}',
  };

  for (let index = 0; index < argv.length; index += 1) {
    const token = argv[index];
    if (!token.startsWith('--')) {
      throw new Error(`Unexpected argument: ${token}`);
    }

    const key = token.slice(2);
    const value = argv[index + 1];
    if (value == null || value.startsWith('--')) {
      throw new Error(`Missing value for --${key}`);
    }

    if (key === 'structure') result.structure = value;
    else if (key === 'format') result.format = value;
    else if (key === 'out-dir') result.outDir = value;
    else if (key === 'width') result.width = Number.parseInt(value, 10);
    else if (key === 'height') result.height = Number.parseInt(value, 10);
    else if (key === 'style-json') result.styleJson = value;
    else throw new Error(`Unknown argument: --${key}`);

    index += 1;
  }

  if (!result.structure) throw new Error('Missing required arg: --structure');
  if (!result.outDir) throw new Error('Missing required arg: --out-dir');
  if (!Number.isFinite(result.width) || result.width <= 0) {
    throw new Error(`Invalid --width value: ${result.width}`);
  }
  if (!Number.isFinite(result.height) || result.height <= 0) {
    throw new Error(`Invalid --height value: ${result.height}`);
  }

  return result;
}

async function loadPlaywright() {
  try {
    const playwright = await import('playwright');
    if (!playwright.chromium) {
      throw new Error('Playwright chromium export missing');
    }
    return playwright.chromium;
  } catch (error) {
    throw new Error(
      `Playwright is required for thumbnail rendering. Install with: npm install --save-dev playwright. Original error: ${error.message}`,
    );
  }
}

async function readStructure(structurePath, format) {
  const normalized = String(format).toLowerCase();
  const binary = normalized === 'bcif';

  if (binary) {
    const bytes = await fs.readFile(structurePath);
    return {
      isBinary: true,
      format: normalized,
      binaryBase64: bytes.toString('base64'),
      textData: null,
    };
  }

  const textData = await fs.readFile(structurePath, 'utf8');
  return {
    isBinary: false,
    format: normalized,
    binaryBase64: null,
    textData,
  };
}

async function setupViewer(page, width, height, backgroundHex) {
  const molstarJsUrl = process.env.MOLSTAR_VIEWER_JS_URL || 'https://unpkg.com/molstar/build/viewer/molstar.js';
  const molstarCssUrl = process.env.MOLSTAR_VIEWER_CSS_URL || 'https://unpkg.com/molstar/build/viewer/molstar.css';

  const html = [
    '<!doctype html>',
    '<html>',
    '<head>',
    '<meta charset="utf-8" />',
    `<link rel="stylesheet" href="${molstarCssUrl}" />`,
    '<style>',
    'html, body, #app {',
    '  margin: 0;',
    '  width: 100%;',
    '  height: 100%;',
    `  background: ${backgroundHex};`,
    '  overflow: hidden;',
    '}',
    '</style>',
    '</head>',
    '<body>',
    '<div id="app"></div>',
    `<script src="${molstarJsUrl}"></script>`,
    '</body>',
    '</html>',
  ].join('\n');

  await page.setViewportSize({ width, height });
  await page.setContent(html, { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => Boolean(globalThis.molstar && globalThis.molstar.Viewer), {
    timeout: 45_000,
  });

  await page.evaluate(async () => {
    const viewer = await globalThis.molstar.Viewer.create('app', {
      layoutIsExpanded: false,
      layoutShowControls: false,
      layoutShowLog: false,
      layoutShowLeftPanel: false,
      viewportShowExpand: false,
      viewportShowSelectionMode: false,
      viewportShowAnimation: false,
      collapseLeftPanel: true,
      pluginStateServer: false,
    });
    globalThis.__volumizerViewer = viewer;
  });
}

async function waitForViewerReady(page, timeoutMs = 15000) {
  await page.waitForFunction(
    () => {
      const viewer = globalThis.__volumizerViewer;
      if (!viewer || !viewer.plugin) return false;
      const canvas = viewer.plugin.canvas3d || viewer.plugin.canvas3dContext?.canvas3d;
      return Boolean(canvas);
    },
    { timeout: timeoutMs },
  );
}

async function loadStructureIntoViewer(page, structureData) {
  await page.evaluate(async (payload) => {
    const viewer = globalThis.__volumizerViewer;
    if (!viewer) throw new Error('Mol* viewer not initialized');

    const parts = [];
    if (payload.isBinary) {
      const bytes = Uint8Array.from(atob(payload.binaryBase64), (char) => char.charCodeAt(0));
      parts.push(bytes);
    } else {
      parts.push(payload.textData);
    }

    const blob = new Blob(parts, {
      type: payload.isBinary ? 'application/octet-stream' : 'text/plain',
    });
    const blobUrl = URL.createObjectURL(blob);
    try {
      await viewer.loadStructureFromUrl(blobUrl, payload.format, payload.isBinary);
      await new Promise((resolve) => requestAnimationFrame(() => requestAnimationFrame(resolve)));
    } finally {
      URL.revokeObjectURL(blobUrl);
    }

    // Render volumes as gaussian surface instead of ball-and-stick
    try {
      const plugin = viewer.plugin;
      if (plugin && plugin.managers && plugin.managers.structure) {
        const structures = plugin.managers.structure.hierarchy.current.structures;
        if (structures && structures.length > 0) {
          const structRef = structures[0];

          const build = plugin.state.data.build();
          for (const comp of structRef.components) {
            build.delete(comp.cell.transform.ref);
          }
          await build.commit();

          const polymer = await plugin.builders.structure.tryCreateComponentStatic(
            structRef.cell, 'polymer', { label: 'Protein' },
          );
          if (polymer) {
            await plugin.builders.structure.representation.addRepresentation(
              polymer, { type: 'cartoon' },
            );
          }

          for (const componentType of ['ligand', 'non-standard']) {
            const comp = await plugin.builders.structure.tryCreateComponentStatic(
              structRef.cell, componentType, { label: 'Volumes' },
            );
            if (comp) {
              await plugin.builders.structure.representation.addRepresentation(
                comp, { type: 'gaussian-surface' },
              );
            }
          }

          await new Promise((resolve) => requestAnimationFrame(() => requestAnimationFrame(resolve)));
        }
      }
    } catch (styleError) {
      console.warn('Volume surface styling failed:', styleError);
    }
  }, structureData);
}

async function setAxisView(page, axis) {
  return await page.evaluate(async (axisName) => {
    const viewer = globalThis.__volumizerViewer;
    const plugin = viewer?.plugin;
    const canvas3d = plugin?.canvas3d || plugin?.canvas3dContext?.canvas3d;
    if (!canvas3d || typeof canvas3d.requestCameraReset !== 'function') {
      return false;
    }

    const sphere = canvas3d.boundingSphereVisible || canvas3d.boundingSphere;
    const center = sphere?.center || [0, 0, 0];
    const radius = Math.max(Number(sphere?.radius) || 0, 40);
    const distance = radius * 2.5;

    let direction = [1, 0, 0];
    let up = [0, 1, 0];

    if (axisName === 'y') {
      direction = [0, 1, 0];
      up = [0, 0, 1];
    } else if (axisName === 'z') {
      direction = [0, 0, 1];
      up = [0, 1, 0];
    }

    const snapshot = canvas3d.camera.getSnapshot();
    snapshot.target = [center[0], center[1], center[2]];
    snapshot.position = [
      center[0] + direction[0] * distance,
      center[1] + direction[1] * distance,
      center[2] + direction[2] * distance,
    ];
    snapshot.up = up;

    canvas3d.requestCameraReset({ snapshot, durationMs: 0 });
    await new Promise((resolve) => setTimeout(resolve, 120));
    return true;
  }, axis);
}

async function renderThumbnails(args) {
  const style = JSON.parse(args.styleJson || '{}');
  const structurePath = path.resolve(args.structure);
  const outDir = path.resolve(args.outDir);

  await fs.mkdir(outDir, { recursive: true });

  const chromium = await loadPlaywright();
  const structureData = await readStructure(structurePath, args.format);

  const browser = await chromium.launch({
    headless: true,
    args: ['--use-angle=swiftshader', '--disable-dev-shm-usage'],
  });

  try {
    const context = await browser.newContext({
      viewport: { width: args.width, height: args.height },
      deviceScaleFactor: 1,
    });
    const page = await context.newPage();

    await setupViewer(page, args.width, args.height, style.background_hex || '#ffffff');
    await loadStructureIntoViewer(page, structureData);

    try {
      await waitForViewerReady(page, 15000);
    } catch {
      // Continue and still emit screenshots; camera orientation may not be available.
    }

    const axes = ['x', 'y', 'z'];
    for (const axis of axes) {
      await setAxisView(page, axis);
      const outPath = path.join(outDir, `${axis}.png`);
      await page.screenshot({ path: outPath, type: 'png' });
    }

    await context.close();
  } finally {
    await browser.close();
  }
}

async function main() {
  const args = parseArgs(process.argv.slice(2));
  if (args.help) {
    process.stdout.write(`${USAGE}\n`);
    return;
  }

  await renderThumbnails(args);
}

main().catch((error) => {
  const message = error instanceof Error ? error.stack || error.message : String(error);
  process.stderr.write(`${message}\n`);
  process.exit(1);
});
