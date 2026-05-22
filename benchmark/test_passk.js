import { execSync, exec, execFile } from "child_process";
import { promisify } from "util";
import fs from "fs";
import path from "path";
import { fileURLToPath } from "url";
import { config } from "dotenv";

config(); // load .env

const execAsync = promisify(exec);
const __dirname = path.dirname(fileURLToPath(import.meta.url));

const DATASET_DIR = "C:/MDC_internship/rolvdesktop/rolvdesktop_electron/benchmark/datasets/filtered_gene_bc_matrices/hg19";
const R_EXEC = process.env.R_EXEC || "Rscript";

// Get ROLV's response to a simple PBMC prompt
console.log("Step 1: Calling orchestrator...");
const prompt = "I have three files: barcodes.tsv, genes.tsv, and matrix.mtx in the current directory. Write R code to load this data using Read10X and create a Seurat object called pbmc.";
const payload = Buffer.from(JSON.stringify({ prompt, history: [] })).toString("base64");

let response;
try {
  const stdout = execSync(
    `node "${path.join(__dirname, "scripts", "orchestrator_shim.js")}" "${payload}"`,
    { timeout: 120000, encoding: "utf8" }
  );
  const jsonLine = stdout.split("\n").map(l => l.trim()).filter(l => l.startsWith("{")).pop();
  const parsed = JSON.parse(jsonLine);
  response = parsed.result?.response ?? parsed.result ?? "";
  console.log("\n✅ Got response. First 600 chars:");
  console.log(String(response).slice(0, 600));
} catch (err) {
  console.error("❌ Orchestrator call failed:", err.message.slice(0, 300));
  process.exit(1);
}

// Extract R code
const match = String(response).match(/```(?:\{r[^}]*\}|r|R)\s*\n([\s\S]*?)```/i);
if (!match) {
  console.log("\n⚠️  No R code block found in response.");
  console.log("Full response:\n", String(response).slice(0, 1000));GIT 
  process.exit(1);
}

const rCode = match[1].trim();
console.log("\n✅ Extracted R code:");
console.log("─".repeat(50));
console.log(rCode.slice(0, 800));
console.log("─".repeat(50));

// execute it
console.log("\nStep 3: Executing R code...");
const tmpDir = fs.mkdtempSync(path.join(__dirname, "results", "test_"));
const scriptPath = path.join(tmpDir, "test.R");

// Apply same normalization as benchmark_runner
let normalized = rCode
  .replace(/data_dir\s*<-\s*["']path\/to\/your\/data["']/gi, 'data_dir <- "."')
  .replace(/["']path\/to\/your\/data["']/gi, '"."')
  .replace(/Read10X\s*\(\s*data\.dir\s*=\s*["']path\/to\/your\/data["']\s*\)/gi, 'Read10X(data.dir = ".")');

const fullScript = [
  `setwd("${DATASET_DIR}")`,
  normalized,
  `cat("BENCH_EXEC_OK\\n")`
].join("\n");

fs.writeFileSync(scriptPath, fullScript);
console.log("R script written to:", scriptPath);

const { stdout: rOut, stderr: rErr } = await new Promise((resolve) => {
  execFile(R_EXEC, ["--vanilla", scriptPath], { timeout: 60000 }, (error, stdout, stderr) => {
    resolve({ stdout: stdout || "", stderr: stderr || (error?.message || "") });
  });
});

const success = rOut.includes("BENCH_EXEC_OK");
console.log(success ? "\n✅ R execution PASSED" : "\n❌ R execution FAILED");
if (!success) {
  console.log("\nStdout:", rOut.slice(0, 500));
  console.log("\nStderr (first error):", rErr.slice(0, 800));
}

// Cleanup
try { fs.rmSync(tmpDir, { recursive: true }); } catch {}