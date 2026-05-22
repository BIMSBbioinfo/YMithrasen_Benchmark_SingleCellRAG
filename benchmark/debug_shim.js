/**
 * debug_shim.js — run this directly to see the real shim error
 * 
 * Place in: C:\MDC_internship\rolvdesktop\rolvdesktop_electron\benchmark\
 * Run with: node debug_shim.js
 */
import "dotenv/config";

global.agentConfig = {
  apiUrl:      process.env.ROLV_API_URL || "http://localhost:3000",
  getJwtToken: () => process.env.ROLV_JWT_TOKEN || "",
};

console.log("Step 1: agentConfig set");
console.log("  apiUrl:", global.agentConfig.apiUrl);
console.log("  hasToken:", !!global.agentConfig.getJwtToken());

console.log("\nStep 2: importing orchestrator...");
let orchestrator;
try {
  const mod = await import("../lib/core/orchestrator.js");
  orchestrator = mod.orchestrator;
  console.log("  ✅ orchestrator imported successfully");
  console.log("  type:", typeof orchestrator);
} catch (err) {
  console.error("  ❌ Import failed:");
  console.error("  ", err.message);
  console.error("\nFull error:");
  console.error(err.stack);
  process.exit(1);
}

console.log("\nStep 3: calling orchestrator with test prompt...");
try {
  const noop = () => {};
  const result = await orchestrator("What does UMAP stand for?", [], noop);
  console.log("  ✅ Got response:");
  console.log("  ", (result || "").slice(0, 200));
} catch (err) {
  console.error("  ❌ Call failed:");
  console.error("  ", err.message);
  console.error(err.stack);
}
