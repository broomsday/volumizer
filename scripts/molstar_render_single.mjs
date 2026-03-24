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

/**
 * Clip a text CIF file by removing atom_site records on the camera side
 * of the clipping plane (coordinate > 0 along the given axis).
 *
 * @param {string} cifText  Full CIF file text
 * @param {number} axisIdx  0 = x, 1 = y, 2 = z
 * @returns {string} CIF text with atoms filtered
 */
function clipCifByAxis(cifText, axisIdx) {
  const axisField = ['Cartn_x', 'Cartn_y', 'Cartn_z'][axisIdx];
  const lines = cifText.split('\n');
  const result = [];
  let inAtomSite = false;
  let readingFields = false;
  const fieldNames = [];
  let coordColIdx = -1;

  for (let i = 0; i < lines.length; i++) {
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
      } else if (readingFields) {
        // Transition from field names to data rows
        readingFields = false;
        // Fall through to data handling below
      }

      if (!readingFields) {
        if (
          trimmed === '' ||
          trimmed.startsWith('#') ||
          trimmed.startsWith('_') ||
          trimmed === 'loop_'
        ) {
          // End of atom_site data block
          inAtomSite = false;
          result.push(line);
        } else if (coordColIdx >= 0) {
          const cols = trimmed.split(/\s+/);
          const coord = parseFloat(cols[coordColIdx]);
          if (!isNaN(coord) && coord <= 0) {
            result.push(line);
          }
          // Skip atoms with coord > 0 (the camera-facing half)
        } else {
          result.push(line);
        }
      }
    } else {
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
  }
  return result.join('\n');
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

async function clearViewer(page) {
  await page.evaluate(async () => {
    const viewer = globalThis.__volumizerViewer;
    if (!viewer) return;
    const state = viewer.plugin.state.data;
    const root = state.tree.root.ref;
    const build = state.build();
    for (const child of state.tree.children.get(root)?.toArray() || []) {
      build.delete(child);
    }
    await build.commit();
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
      await new Promise((resolve) => requestAnimationFrame(() => requestAnimationFrame(resolve)));
    } finally {
      URL.revokeObjectURL(blobUrl);
    }

    // Apply colored volume and protein styling
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

          // Protein in green/grey shades — use explicit "not volume"
          // MolScript selection instead of 'polymer' preset because
          // biotite CIF lacks _entity_poly / _struct_asym tables.
          const proteinPalette = [
            0x2E8B57, 0x808080, 0x3CB371, 0xA0A0A0,
            0x228B22, 0x6B6B6B, 0x006400, 0xB5B5B5,
            0x4CAF50, 0x959595, 0x66BB6A, 0x7A7A7A,
          ];
          const volumeIds = ['HUB', 'POR', 'POK', 'CAV', 'OCC'];
          const notVolumeExpr = volumeIds
            .map(id => `(= atom.label_comp_id ${id})`)
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
            await plugin.builders.structure.representation.addRepresentation(
              polymer, {
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

          // Per-type volume components with distinct colors
          const volumeTypes = [
            { id: 'HUB', color: 0xCC3333, label: 'Hubs' },
            { id: 'POR', color: 0xDD8833, label: 'Pores' },
            { id: 'POK', color: 0x3366CC, label: 'Pockets' },
            { id: 'CAV', color: 0xCC33CC, label: 'Cavities' },
          ];

          for (const vt of volumeTypes) {
            try {
              const comp = await plugin.builders.structure.tryCreateComponent(
                structRef.cell,
                {
                  type: {
                    name: 'script',
                    params: {
                      language: 'mol-script',
                      expression: `(sel.atom.atom-groups :residue-test (= atom.label_comp_id ${vt.id}))`,
                    },
                  },
                  nullIfEmpty: true,
                  label: vt.label,
                },
                `volume-${vt.id.toLowerCase()}`,
              );
              if (comp) {
                await plugin.builders.structure.representation.addRepresentation(
                  comp, {
                    type: 'gaussian-surface',
                    color: 'uniform',
                    colorParams: { value: vt.color },
                  },
                );
              }
            } catch (err) {
              console.warn(`Failed to create ${vt.id} volume component:`, err);
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

    try {
      await waitForViewerReady(page, 15000);
    } catch {
      // Continue; camera orientation may not be available.
    }

    const axes = ['x', 'y', 'z'];
    const axisCoordIndex = { x: 0, y: 1, z: 2 };

    for (const axis of axes) {
      // Clip the CIF data: remove atoms on the camera side of the center plane
      let axisData = structureData;
      if (!structureData.isBinary && structureData.textData) {
        const clippedText = clipCifByAxis(structureData.textData, axisCoordIndex[axis]);
        axisData = { ...structureData, textData: clippedText };
      }

      await clearViewer(page);
      await loadStructureIntoViewer(page, axisData);
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
