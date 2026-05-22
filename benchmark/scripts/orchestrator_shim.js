/**
 *
 * Thin bridge between the benchmark runner and orchestrator.js.
 * Receives a base64-encoded JSON payload, calls orchestrator(), writes result to stdout.
  *
 * The orchestrator.js is the backend of the agent mode — it contextualizes the query,
 * routes to the right agent (pubmed / code / geo / vision / base), and returns text.
 * This shim calls that same pipeline, bypassing the Electron window.* APIs entirely,
 * which is exactly what we want for headless benchmarking.
 *
 * The shim also injects a <working_directory_instruction> block that mirrors what
 * mainChat.js injects when agent mode is active, so the orchestrator sees the same
 * context the real app would send.
 */

import "dotenv/config";
import path from "path";
import { fileURLToPath } from "url";

const __dirname = path.dirname(fileURLToPath(import.meta.url));

//  1. Replicate the agentConfig that orchestrator.js reads from global
// In the Electron main process this is set before any IPC call reaches the orchestrator.
// baseAgent.js checks: if (typeof window === "undefined" && global.agentConfig)
// orchestrator.js does the same. We replicate it here so the orchestrator can read the same config in this shim context.
global.agentConfig = {
  apiUrl:      process.env.ROLV_API_URL      || "http://localhost:3000",
  getJwtToken: () => process.env.ROLV_JWT_TOKEN || "",
};

//  2. Import orchestrator
import { orchestrator } from "../../lib/core/orchestrator.js";

//  3. Parse payload
const raw    = process.argv[2];
if (!raw) {
  process.stderr.write("[orchestrator_shim] No payload argument provided\n");
  process.exit(1);
}
const parsed = JSON.parse(Buffer.from(raw, "base64").toString("utf-8"));
const { prompt, history = [], dataset_dir = "." } = parsed;

//  4. Inject the same working directory hint mainChat.js adds in agent mode — this ensures the orchestrator gets the same crucial context whether it's called from the real app or this shim.
const workingDirHint = `

<working_directory_instruction>
A working directory is already selected in the app.
The files needed for this task are already available there.
Do NOT generate placeholder paths like "path/to/your/data" or "path/to/your/data/directory".
Do NOT generate absolute paths unless the user explicitly typed one.
Use relative paths only.
If the task involves barcodes.tsv, genes.tsv/features.tsv, and matrix.mtx in the selected directory, prefer Read10X(data.dir = ".").
</working_directory_instruction>`;

const fullPrompt = typeof prompt === "string"
  ? prompt + workingDirHint
  : prompt; // leave multimodal arrays as-is

//  5. Silence ALL console.log/warn/info so orchestrator's debug boxes don't pollute stdout and corrupt JSON parsing.
// We redirect them to stderr for debugging.
const realStdoutWrite = process.stdout.write.bind(process.stdout);
console.log  = (...a) => process.stderr.write("[orch-log] " + a.join(" ") + "\n");
console.warn = (...a) => process.stderr.write("[orch-warn] " + a.join(" ") + "\n");
console.info = (...a) => process.stderr.write("[orch-info] " + a.join(" ") + "\n");

// Also intercept any direct process.stdout.write calls from the orchestrator
process.stdout.write = (chunk, ...rest) => {
  // Allow only  final sentinel line through
  if (typeof chunk === "string" && chunk.startsWith('{"result"')) {
    return realStdoutWrite(chunk, ...rest);
  }
  process.stderr.write(chunk, ...rest);
  return true;
};

const noop = () => {};

try {
  const raw = await orchestrator(fullPrompt, history, noop);
  // orchestrator returns { response: string, telemetry: {...} }
  const result = (raw && typeof raw === "object" && raw.response !== undefined)
    ? raw.response
    : raw;
  // Write ONLY this JSON line to real stdout — benchmark_runner.js parses it
  realStdoutWrite(JSON.stringify({ result }) + "\n");
} catch (err) {
  process.stderr.write(`[orchestrator_shim] ${err.stack}\n`);
  realStdoutWrite(JSON.stringify({ result: null, error: err.message }) + "\n");
  process.exit(1);
}
