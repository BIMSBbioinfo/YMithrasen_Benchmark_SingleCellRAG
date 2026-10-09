// orchestrator_shim_nopreexec.js
// Agentic RAG without the pre-execution planner (model key "rolv-nopreexec").
// Same as orchestrator_shim.js except for preExecEnabled: false.
// Calls the ROLV orchestrator directly, without the Electron app, and prints one JSON
// line: { result, error? }. Orchestrator logging is redirected to stderr.
//
// Usage: node orchestrator_shim_nopreexec.js <base64 JSON { prompt, history? }>
// Needs ROLV_API_URL, ROLV_JWT_TOKEN and ROLV_ORCHESTRATOR_PATH in .env.

import "dotenv/config";
import path from "path";
import { fileURLToPath, pathToFileURL } from "url";

const __dirname = path.dirname(fileURLToPath(import.meta.url));

// 1. agentConfig that orchestrator.js reads from global (set by the Electron main process in the app)
global.agentConfig = {
  preExecEnabled: false,
  apiUrl:      process.env.ROLV_API_URL      || "http://localhost:3000",
  getJwtToken: () => process.env.ROLV_JWT_TOKEN || "",
};

// 2. Import the orchestrator. The ROLV app is a separate codebase;
// ROLV_ORCHESTRATOR_PATH points at its lib/core/orchestrator.js.
const orchestratorPath = process.env.ROLV_ORCHESTRATOR_PATH;
if (!orchestratorPath) {
  const message = "ROLV_ORCHESTRATOR_PATH is not set — point it at the ROLV app's lib/core/orchestrator.js";
  process.stderr.write(`[orchestrator_shim] ${message}\n`);
  process.stdout.write(JSON.stringify({ result: null, error: message }) + "\n");
  process.exit(1);
}
const { orchestrator } = await import(pathToFileURL(path.resolve(orchestratorPath)).href);

//  3. Parse payload
const raw    = process.argv[2];
if (!raw) {
  process.stderr.write("[orchestrator_shim] No payload argument provided\n");
  process.exit(1);
}
const parsed = JSON.parse(Buffer.from(raw, "base64").toString("utf-8"));
const { prompt, history = [], dataset_dir = "." } = parsed;

// 4. Same working-directory hint that mainChat.js adds in agent mode
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

// 5. Send console output to stderr so it cannot corrupt the JSON on stdout
const realStdoutWrite = process.stdout.write.bind(process.stdout);
console.log  = (...a) => process.stderr.write("[orch-log] " + a.join(" ") + "\n");
console.warn = (...a) => process.stderr.write("[orch-warn] " + a.join(" ") + "\n");
console.info = (...a) => process.stderr.write("[orch-info] " + a.join(" ") + "\n");

// also catch direct process.stdout.write calls
process.stdout.write = (chunk, ...rest) => {
  // only the final result line goes through
  if (typeof chunk === "string" && chunk.startsWith('{"result"')) {
    return realStdoutWrite(chunk, ...rest);
  }
  process.stderr.write(chunk, ...rest);
  return true;
};

const noop = () => {};

try {
  const raw = await orchestrator(fullPrompt, history, noop);
  // orchestrator returns { response, telemetry }
  const result = (raw && typeof raw === "object" && raw.response !== undefined)
    ? raw.response
    : raw;
  // the one line benchmark_runner.js parses
  realStdoutWrite(JSON.stringify({ result }) + "\n");
} catch (err) {
  process.stderr.write(`[orchestrator_shim] ${err.stack}\n`);
  realStdoutWrite(JSON.stringify({ result: null, error: err.message }) + "\n");
  process.exit(1);
}
