#!/usr/bin/env node

import fs from 'node:fs/promises';
import path from 'node:path';
import { performance } from 'node:perf_hooks';
import { fileURLToPath } from 'node:url';

const SCRIPT_PATH = fileURLToPath(import.meta.url);
const SCRIPT_DIR = path.dirname(SCRIPT_PATH);
const REPO_ROOT = path.resolve(SCRIPT_DIR, '..');
const NODE_TIMING_PREFIX = '__VOLUMIZER_TIMING__ ';
const VALID_RENDER_BACKENDS = new Set(['software', 'hardware', 'auto']);
const VALID_AXIS_RENDER_MODES = new Set(['compatibility', 'fast']);

const USAGE = [
  'Usage: node scripts/molstar_render_single.mjs --structure PATH --format mmcif|pdb|bcif --out-dir DIR [options]',
  '',
  'Options:',
  '  --width N                  Output image width (default: 320)',
  '  --height N                 Output image height (default: 240)',
  '  --style-json JSON          Optional style payload passed from Python queue',
  '  --render-backend MODE      software|hardware|auto (default: software)',
  '  --axis-render-mode MODE    compatibility|fast (default: compatibility)',
  '  --timing-jsonl PATH        Optional JSONL path for renderer timing output',
  '  -h, --help                 Show this help text',
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
    renderBackend: 'software',
    axisRenderMode: 'compatibility',
    timingJsonl: null,
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
    else if (key === 'render-backend') result.renderBackend = String(value).toLowerCase();
    else if (key === 'axis-render-mode') result.axisRenderMode = String(value).toLowerCase();
    else if (key === 'timing-jsonl') result.timingJsonl = value;
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
  if (!VALID_RENDER_BACKENDS.has(result.renderBackend)) {
    throw new Error('Invalid --render-backend value');
  }
  if (!VALID_AXIS_RENDER_MODES.has(result.axisRenderMode)) {
    throw new Error('Invalid --axis-render-mode value');
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
      `Playwright is required for thumbnail rendering. Install with: npm install --save-dev playwright molstar. Original error: ${error.message}`,
    );
  }
}

async function pathExists(candidatePath) {
  try {
    await fs.access(candidatePath);
    return true;
  } catch {
    return false;
  }
}

async function appendJsonl(targetPath, payload) {
  const resolved = path.resolve(targetPath);
  await fs.mkdir(path.dirname(resolved), { recursive: true });
  await fs.appendFile(resolved, `${JSON.stringify(payload)}\n`, 'utf8');
}

async function resolveMolstarAssets() {
  const jsOverride = process.env.MOLSTAR_VIEWER_JS_URL || null;
  const cssOverride = process.env.MOLSTAR_VIEWER_CSS_URL || null;
  const assetRoot = process.env.MOLSTAR_ASSET_ROOT
    ? path.resolve(process.env.MOLSTAR_ASSET_ROOT)
    : path.join(REPO_ROOT, 'node_modules', 'molstar', 'build', 'viewer');
  const localJsPath = path.join(assetRoot, 'molstar.js');
  const localCssPath = path.join(assetRoot, 'molstar.css');

  const jsAsset = jsOverride
    ? { url: jsOverride }
    : ((await pathExists(localJsPath)) ? { path: localJsPath } : null);
  const cssAsset = cssOverride
    ? { url: cssOverride }
    : ((await pathExists(localCssPath)) ? { path: localCssPath } : null);

  if (!jsAsset || !cssAsset) {
    throw new Error(
      `Mol* viewer assets not found. Expected local assets under ${assetRoot}. `
      + 'Install npm dependencies or set MOLSTAR_ASSET_ROOT / '
      + 'MOLSTAR_VIEWER_JS_URL / MOLSTAR_VIEWER_CSS_URL.',
    );
  }

  return { jsAsset, cssAsset };
}

async function addScriptAsset(page, asset) {
  if (asset.path) {
    await page.addScriptTag({ path: asset.path });
    return;
  }
  await page.addScriptTag({ url: asset.url });
}

async function addStyleAsset(page, asset) {
  if (asset.path) {
    await page.addStyleTag({ path: asset.path });
    return;
  }
  await page.addStyleTag({ url: asset.url });
}

async function waitForFrames(page, frameCount = 1) {
  await page.evaluate(async (requestedCount) => {
    const count = Math.max(1, Number(requestedCount) || 1);
    await new Promise((resolve) => {
      let remaining = count;
      const tick = () => {
        remaining -= 1;
        if (remaining <= 0) {
          resolve(true);
          return;
        }
        requestAnimationFrame(tick);
      };
      requestAnimationFrame(tick);
    });
  }, frameCount);
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

function clipCifByAxis(cifText, axisIdx) {
  const axisField = ['Cartn_x', 'Cartn_y', 'Cartn_z'][axisIdx];
  const lines = cifText.split('\n');
  const result = [];
  let inAtomSite = false;
  let readingFields = false;
  const fieldNames = [];
  let coordColIdx = -1;

  for (let i = 0; i < lines.length; i += 1) {
    const line = lines[i];
    const trimmed = line.trim();

    if (inAtomSite) {
      if (readingFields && trimmed.startsWith('_atom_site.')) {
        const fieldName = trimmed.split(/\s+/)[0].replace('_atom_site.', '');
        fieldNames.push(fieldName);
        if (fieldName === axisField) {
          coordColIdx = fieldNames.length - 1;
        }
        result.push(line);
        continue;
      }

      if (readingFields) {
        readingFields = false;
      }

      if (
        trimmed === ''
        || trimmed.startsWith('#')
        || trimmed.startsWith('_')
        || trimmed === 'loop_'
      ) {
        inAtomSite = false;
        result.push(line);
      } else if (coordColIdx >= 0) {
        const cols = trimmed.split(/\s+/);
        const coord = Number.parseFloat(cols[coordColIdx]);
        if (!Number.isNaN(coord) && coord <= 0) {
          result.push(line);
        }
      } else {
        result.push(line);
      }
      continue;
    }

    if (trimmed === 'loop_') {
      const nextLine = (lines[i + 1] || '').trim();
      if (nextLine.startsWith('_atom_site.')) {
        inAtomSite = true;
        readingFields = true;
        fieldNames.length = 0;
        coordColIdx = -1;
      }
    }
    result.push(line);
  }

  return result.join('\n');
}

function buildAxisStructureData(structureData, axisRenderMode) {
  if (axisRenderMode === 'fast' || structureData.isBinary || !structureData.textData) {
    return {
      x: structureData,
      y: structureData,
      z: structureData,
    };
  }

  return {
    x: { ...structureData, textData: clipCifByAxis(structureData.textData, 0) },
    y: { ...structureData, textData: clipCifByAxis(structureData.textData, 1) },
    z: { ...structureData, textData: clipCifByAxis(structureData.textData, 2) },
  };
}

function getBrowserLaunchArgs(renderBackend) {
  const args = ['--disable-dev-shm-usage'];
  if (renderBackend === 'software') {
    args.unshift('--use-angle=swiftshader');
  }
  return args;
}

async function setupViewer(page, width, height, backgroundHex, molstarAssets) {
  const html = [
    '<!doctype html>',
    '<html>',
    '<head>',
    '<meta charset="utf-8" />',
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
    '</body>',
    '</html>',
  ].join('\n');

  await page.setViewportSize({ width, height });
  await page.setContent(html, { waitUntil: 'domcontentloaded' });
  await addStyleAsset(page, molstarAssets.cssAsset);
  await addScriptAsset(page, molstarAssets.jsAsset);
  await page.waitForFunction(
    () => Boolean(globalThis.molstar && globalThis.molstar.Viewer),
    { timeout: 45_000 },
  );

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

async function clearViewer(page) {
  await page.evaluate(async () => {
    const viewer = globalThis.__volumizerViewer;
    if (!viewer) return;
    await viewer.plugin.clear();
  });
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
    } finally {
      URL.revokeObjectURL(blobUrl);
    }

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

          const proteinPalette = [
            0x2E8B57, 0x808080, 0x3CB371, 0xA0A0A0,
            0x228B22, 0x6B6B6B, 0x006400, 0xB5B5B5,
            0x4CAF50, 0x959595, 0x66BB6A, 0x7A7A7A,
          ];
          const volumeIds = ['HUB', 'POR', 'POK', 'CAV', 'OCC'];
          const notVolumeExpr = volumeIds
            .map((id) => `(= atom.label_comp_id ${id})`)
            .join(' ');
          const polymer = await plugin.builders.structure.tryCreateComponent(
            structRef.cell,
            {
              type: {
                name: 'script',
                params: {
                  language: 'mol-script',
                  expression: `(sel.atom.atom-groups :residue-test (not (or ${notVolumeExpr})))`,
                },
              },
              nullIfEmpty: true,
              label: 'Protein',
            },
            'protein',
          );
          if (polymer) {
            const polymerRepr = await plugin.builders.structure.representation.addRepresentation(
              polymer,
              {
                type: 'cartoon',
                color: 'chain-id',
                colorParams: {
                  palette: {
                    name: 'colors',
                    params: {
                      list: { kind: 'set', colors: proteinPalette },
                    },
                  },
                },
              },
            );
          }

          const volumeTypes = [
            { id: 'HUB', color: 0xCC3333, label: 'Hubs' },
            { id: 'POR', color: 0xDD8833, label: 'Pores' },
            { id: 'POK', color: 0x3366CC, label: 'Pockets' },
            { id: 'CAV', color: 0xCC33CC, label: 'Cavities' },
          ];

          for (const volumeType of volumeTypes) {
            try {
              const comp = await plugin.builders.structure.tryCreateComponent(
                structRef.cell,
                {
                  type: {
                    name: 'script',
                    params: {
                      language: 'mol-script',
                      expression: `(sel.atom.atom-groups :residue-test (= atom.label_comp_id ${volumeType.id}))`,
                    },
                  },
                  nullIfEmpty: true,
                  label: volumeType.label,
                },
                `volume-${volumeType.id.toLowerCase()}`,
              );
              if (comp) {
                const volumeSurfaceTypeParams = {
                  quality: 'custom',
                  doubleSided: false,
                  interior: {
                    color: volumeType.color,
                    colorStrength: 1,
                    substance: { metalness: 0, roughness: 1, bumpiness: 0 },
                    substanceStrength: 0,
                  },
                };
                const volumeRepr = await plugin.builders.structure.representation.addRepresentation(
                  comp,
                  {
                    type: 'gaussian-surface',
                    typeParams: volumeSurfaceTypeParams,
                    color: 'uniform',
                    colorParams: { value: volumeType.color },
                  },
                );
              }
            } catch (err) {
              console.warn(`Failed to create ${volumeType.id} volume component:`, err);
            }
          }
        }
      }
    } catch (styleError) {
      console.warn('Volume surface styling failed:', styleError);
    }

    await new Promise((resolve) => requestAnimationFrame(() => resolve(true)));
  }, structureData);
}

async function setAxisView(page, axis, timeoutMs = 1000) {
  const snapshot = await page.evaluate((axisName) => {
    const viewer = globalThis.__volumizerViewer;
    const plugin = viewer?.plugin;
    const canvas3d = plugin?.canvas3d || plugin?.canvas3dContext?.canvas3d;
    if (!canvas3d || typeof canvas3d.requestCameraReset !== 'function') {
      return null;
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

    const snapshotValue = canvas3d.camera.getSnapshot();
    snapshotValue.target = [center[0], center[1], center[2]];
    snapshotValue.position = [
      center[0] + direction[0] * distance,
      center[1] + direction[1] * distance,
      center[2] + direction[2] * distance,
    ];
    snapshotValue.up = up;

    canvas3d.requestCameraReset({ snapshot: snapshotValue, durationMs: 0 });
    return {
      target: snapshotValue.target,
      position: snapshotValue.position,
      up: snapshotValue.up,
    };
  }, axis);

  if (!snapshot) {
    return false;
  }

  try {
    await page.waitForFunction(
      (expectedSnapshot) => {
        const viewer = globalThis.__volumizerViewer;
        const plugin = viewer?.plugin;
        const canvas3d = plugin?.canvas3d || plugin?.canvas3dContext?.canvas3d;
        const current = canvas3d?.camera?.getSnapshot?.();
        if (!current) return false;

        const nearlyEqual = (left, right) => Math.abs(Number(left) - Number(right)) < 0.01;
        const sameVector = (left, right) =>
          Array.isArray(left)
          && Array.isArray(right)
          && left.length === right.length
          && left.every((value, index) => nearlyEqual(value, right[index]));

        return (
          sameVector(current.position, expectedSnapshot.position)
          && sameVector(current.target, expectedSnapshot.target)
          && sameVector(current.up, expectedSnapshot.up)
        );
      },
      snapshot,
      { timeout: timeoutMs },
    );
  } catch {
    // Fall through to a single frame wait below.
  }

  await waitForFrames(page, 1);
  return true;
}

async function renderThumbnailsOnce({
  chromium,
  args,
  molstarAssets,
  structureData,
  style,
  renderBackend,
}) {
  const timing = {
    requested_backend: args.renderBackend,
    backend_used: renderBackend,
    fallback_used: false,
    axis_render_mode: args.axisRenderMode,
    browser_launch_ms: 0,
    viewer_setup_ms: 0,
    first_structure_load_ms: null,
    clip_prep_ms: 0,
    axes: {},
  };

  const browserLaunchStartedAt = performance.now();
  const browser = await chromium.launch({
    headless: true,
    args: getBrowserLaunchArgs(renderBackend),
  });
  timing.browser_launch_ms = performance.now() - browserLaunchStartedAt;

  try {
    const context = await browser.newContext({
      viewport: { width: args.width, height: args.height },
      deviceScaleFactor: 1,
    });

    try {
      const page = await context.newPage();
      const viewerSetupStartedAt = performance.now();
      await setupViewer(page, args.width, args.height, style.background_hex || '#ffffff', molstarAssets);

      try {
        await waitForViewerReady(page, 15000);
      } catch {
        // Continue; camera orientation may not be available yet.
      }

      timing.viewer_setup_ms = performance.now() - viewerSetupStartedAt;

      const clipPrepStartedAt = performance.now();
      const axisStructureData = buildAxisStructureData(structureData, args.axisRenderMode);
      timing.clip_prep_ms = performance.now() - clipPrepStartedAt;

      const axes = ['x', 'y', 'z'];
      if (args.axisRenderMode === 'fast') {
        const loadStartedAt = performance.now();
        await loadStructureIntoViewer(page, structureData);
        timing.first_structure_load_ms = performance.now() - loadStartedAt;

        for (const axis of axes) {
          const axisStartedAt = performance.now();
          await setAxisView(page, axis);

          const outPath = path.join(path.resolve(args.outDir), `${axis}.png`);
          await page.screenshot({ path: outPath, type: 'png' });

          timing.axes[axis] = {
            render_ms: performance.now() - axisStartedAt,
            reloaded_structure: false,
            structure_load_ms: null,
            clip_update_ms: 0,
          };
        }
      } else {
        for (let axisIndex = 0; axisIndex < axes.length; axisIndex += 1) {
          const axis = axes[axisIndex];
          const axisStartedAt = performance.now();

          if (axisIndex > 0) {
            await clearViewer(page);
          }

          const loadStartedAt = performance.now();
          await loadStructureIntoViewer(page, axisStructureData[axis]);
          const structureLoadMs = performance.now() - loadStartedAt;
          if (timing.first_structure_load_ms === null) {
            timing.first_structure_load_ms = structureLoadMs;
          }

          await setAxisView(page, axis);

          const outPath = path.join(path.resolve(args.outDir), `${axis}.png`);
          await page.screenshot({ path: outPath, type: 'png' });

          timing.axes[axis] = {
            render_ms: performance.now() - axisStartedAt,
            reloaded_structure: axisIndex > 0,
            structure_load_ms: axisIndex > 0 ? structureLoadMs : null,
            clip_update_ms: 0,
          };
        }
      }
    } finally {
      await context.close();
    }
  } finally {
    await browser.close();
  }

  return timing;
}

async function renderThumbnails(args) {
  const style = JSON.parse(args.styleJson || '{}');
  const structurePath = path.resolve(args.structure);
  const outDir = path.resolve(args.outDir);
  const timingJsonlPath = args.timingJsonl ? path.resolve(args.timingJsonl) : null;

  await fs.mkdir(outDir, { recursive: true });

  const [chromium, molstarAssets, structureData] = await Promise.all([
    loadPlaywright(),
    resolveMolstarAssets(),
    readStructure(structurePath, args.format),
  ]);

  const renderStartedAt = performance.now();
  let timing = null;

  if (args.renderBackend === 'auto') {
    try {
      timing = await renderThumbnailsOnce({
        chromium,
        args,
        molstarAssets,
        structureData,
        style,
        renderBackend: 'hardware',
      });
      timing.requested_backend = 'auto';
      timing.backend_used = 'hardware';
    } catch (hardwareError) {
      const hardwareMessage =
        hardwareError instanceof Error ? hardwareError.message : String(hardwareError);
      process.stderr.write(
        `[molstar-render] hardware backend failed, retrying with software: ${hardwareMessage}\n`,
      );
      try {
        timing = await renderThumbnailsOnce({
          chromium,
          args,
          molstarAssets,
          structureData,
          style,
          renderBackend: 'software',
        });
      } catch (softwareError) {
        const softwareMessage =
          softwareError instanceof Error ? softwareError.message : String(softwareError);
        throw new Error(
          `hardware backend failed (${hardwareMessage}); software fallback failed (${softwareMessage})`,
        );
      }
      timing.requested_backend = 'auto';
      timing.backend_used = 'software';
      timing.fallback_used = true;
      timing.fallback_error = hardwareMessage;
    }
  } else {
    timing = await renderThumbnailsOnce({
      chromium,
      args,
      molstarAssets,
      structureData,
      style,
      renderBackend: args.renderBackend,
    });
  }

  timing.total_ms = performance.now() - renderStartedAt;
  timing.structure = structurePath;
  timing.out_dir = outDir;
  timing.format = args.format;

  if (timingJsonlPath) {
    await appendJsonl(timingJsonlPath, {
      event: 'render',
      ...timing,
    });
  }

  process.stdout.write(`${NODE_TIMING_PREFIX}${JSON.stringify(timing)}\n`);
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
