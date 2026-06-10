const state = {
  files: [],
  jobs: [],
  lastJobsSignature: "",
};

const elements = {
  form: document.getElementById("uploadForm"),
  fileInput: document.getElementById("files"),
  dropzone: document.getElementById("dropzone"),
  selectedFiles: document.getElementById("selectedFiles"),
  threads: document.getElementById("threads"),
  maxTargetSeqs: document.getElementById("maxTargetSeqs"),
  submitButton: document.getElementById("submitButton"),
  formMessage: document.getElementById("formMessage"),
  jobsEmpty: document.getElementById("jobsEmpty"),
  jobsList: document.getElementById("jobsList"),
  healthText: document.getElementById("healthText"),
  totalJobs: document.getElementById("totalJobs"),
  runningJobs: document.getElementById("runningJobs"),
  completedJobs: document.getElementById("completedJobs"),
};

function setFormMessage(message, tone = "") {
  elements.formMessage.textContent = message;
  elements.formMessage.className = `form-message ${tone}`.trim();
}

function escapeHtml(value) {
  return String(value ?? "")
    .replaceAll("&", "&amp;")
    .replaceAll("<", "&lt;")
    .replaceAll(">", "&gt;")
    .replaceAll('"', "&quot;")
    .replaceAll("'", "&#39;");
}

function formatTime(value) {
  if (!value) return "";
  const date = new Date(value);
  if (Number.isNaN(date.getTime())) return value;
  return date.toLocaleString();
}

function toDisplayInputKind(inputKind) {
  if (inputKind === "fastq") return "FASTQ -> consensus -> BLAST";
  if (inputKind === "fasta") return "FASTA -> BLAST";
  return inputKind || "";
}

function toDisplayStage(stage) {
  if (stage === "consensus") return "Building consensus";
  if (stage === "blast") return "Running BLAST";
  if (stage === "queued") return "Queued";
  if (stage === "completed") return "Completed";
  if (stage === "failed") return "Failed";
  return stage || "";
}

function updateSelectedFiles(files) {
  state.files = Array.from(files);
  if (!state.files.length) {
    elements.selectedFiles.innerHTML = "";
    return;
  }

  elements.selectedFiles.innerHTML = state.files
    .map((file) => {
      const sizeKb = Math.ceil(file.size / 1024);
      return `<li><span>${escapeHtml(file.name)}</span><span class="table-note">${sizeKb} KB</span></li>`;
    })
    .join("");
}

function updateOverview(jobs) {
  const runningCount = jobs.filter((job) => job.status === "queued" || job.status === "running").length;
  const completedCount = jobs.filter((job) => job.status === "completed").length;

  elements.totalJobs.textContent = String(jobs.length);
  elements.runningJobs.textContent = String(runningCount);
  elements.completedJobs.textContent = String(completedCount);
}

function stableStringify(value) {
  if (Array.isArray(value)) {
    return `[${value.map(stableStringify).join(",")}]`;
  }
  if (value && typeof value === "object") {
    const keys = Object.keys(value).sort();
    return `{${keys.map((key) => `${JSON.stringify(key)}:${stableStringify(value[key])}`).join(",")}}`;
  }
  return JSON.stringify(value);
}

function collectSummaryMetrics(job) {
  const counts = job.confidence_counts || {};
  const metrics = [
    `mode ${toDisplayInputKind(job.input_kind) || "-"}`,
    `stage ${toDisplayStage(job.stage) || "-"}`,
    `rows ${job.total_rows || 0}`,
    `files ${job.file_count || 0}`,
  ];

  Object.entries(counts).forEach(([key, value]) => {
    metrics.push(`${key} ${value}`);
  });

  return metrics.map((metric) => `<span class="metric">${escapeHtml(metric)}</span>`).join("");
}

function renderPreviewTable(rows) {
  if (!rows || !rows.length) {
    return "";
  }

  const configuredHeaders = window.SANGER16S_CONFIG.csvFields || [];
  const headers = configuredHeaders.length ? configuredHeaders : Object.keys(rows[0]);
  const thead = headers.map((header) => `<th>${escapeHtml(header)}</th>`).join("");
  const tbody = rows
    .map((row) => {
      const cells = headers.map((header) => `<td>${escapeHtml(row[header] ?? "")}</td>`).join("");
      return `<tr>${cells}</tr>`;
    })
    .join("");

  return `
    <div class="result-table-wrap">
      <table>
        <thead><tr>${thead}</tr></thead>
        <tbody>${tbody}</tbody>
      </table>
    </div>
  `;
}

function renderJobs(jobs) {
  state.jobs = jobs;
  updateOverview(jobs);
  elements.jobsEmpty.style.display = jobs.length ? "none" : "grid";

  elements.jobsList.innerHTML = jobs
    .map((job) => {
      const statusClass = `status-${job.status}`;
      const downloadLink = job.download_url ? `<a href="${escapeHtml(job.download_url)}">Download CSV</a>` : "";
      const preview = job.status === "completed" ? renderPreviewTable(job.preview_rows || []) : "";
      const completedAt = job.completed_at ? `Completed ${formatTime(job.completed_at)}` : "";
      const error = job.error ? `<div class="job-error">${escapeHtml(job.error)}</div>` : "";

      return `
        <article class="job-card">
          <div class="job-topbar">
            <div class="job-title">
              <div class="job-id">${escapeHtml(job.job_id)}</div>
              <div class="job-meta">Created ${escapeHtml(formatTime(job.created_at))}</div>
            </div>
            <span class="status-badge ${statusClass}">${escapeHtml(job.status)}</span>
          </div>
          <div class="job-body">
            <div class="job-files">${escapeHtml((job.files || []).join(", "))}</div>
            <div class="job-summary">${collectSummaryMetrics(job)}</div>
            ${error}
            <div class="job-actions">
              <span class="table-note">${escapeHtml(completedAt)}</span>
              ${downloadLink}
            </div>
            ${preview}
          </div>
        </article>
      `;
    })
    .join("");
}

async function refreshHealth() {
  try {
    const response = await fetch("/api/health");
    if (!response.ok) {
      throw new Error("Health check failed.");
    }

    const health = await response.json();
    const ok = health.ok && health.consensus_pipeline_ready;
    elements.healthText.textContent = ok ? "Ready" : "Check server";
    elements.healthText.className = `overview-value ${ok ? "status-ok" : "status-error"}`;
  } catch (error) {
    elements.healthText.textContent = "Unavailable";
    elements.healthText.className = "overview-value status-error";
  }
}

async function refreshJobs() {
  try {
    const response = await fetch("/api/jobs");
    if (!response.ok) {
      throw new Error("Job list request failed.");
    }

    const jobs = await response.json();
    const normalizedJobs = Array.isArray(jobs) ? jobs : [];
    const nextSignature = stableStringify(normalizedJobs);

    if (nextSignature === state.lastJobsSignature) {
      return;
    }

    state.lastJobsSignature = nextSignature;
    renderJobs(normalizedJobs);
  } catch (error) {
    setFormMessage("Failed to refresh job status.", "error");
  }
}

async function submitJob(event) {
  event.preventDefault();

  if (!state.files.length) {
    setFormMessage("Select FASTA or FASTQ files first.", "error");
    return;
  }

  const formData = new FormData();
  state.files.forEach((file) => formData.append("files", file));
  formData.append("threads", elements.threads.value || window.SANGER16S_CONFIG.defaultThreads);
  formData.append("max_target_seqs", elements.maxTargetSeqs.value || window.SANGER16S_CONFIG.defaultMaxTargetSeqs);

  elements.submitButton.disabled = true;
  setFormMessage("Job submitted.", "success");

  try {
    const response = await fetch("/api/jobs", {
      method: "POST",
      body: formData,
    });
    const payload = await response.json();

    if (!response.ok) {
      throw new Error(payload.error || "Job submission failed.");
    }

    state.lastJobsSignature = "";
    setFormMessage(`Job ${payload.job_id} queued.`, "success");
    elements.form.reset();
    updateSelectedFiles([]);
    await refreshJobs();
  } catch (error) {
    setFormMessage(error.message || "Job submission failed.", "error");
  } finally {
    elements.submitButton.disabled = false;
  }
}

function bindDropzone() {
  ["dragenter", "dragover"].forEach((eventName) => {
    elements.dropzone.addEventListener(eventName, (event) => {
      event.preventDefault();
      elements.dropzone.classList.add("dragover");
    });
  });

  ["dragleave", "drop"].forEach((eventName) => {
    elements.dropzone.addEventListener(eventName, (event) => {
      event.preventDefault();
      elements.dropzone.classList.remove("dragover");
    });
  });

  elements.dropzone.addEventListener("drop", (event) => {
    const droppedFiles = event.dataTransfer?.files;
    if (droppedFiles) {
      updateSelectedFiles(droppedFiles);
    }
  });
}

elements.fileInput.addEventListener("change", (event) => {
  updateSelectedFiles(event.target.files || []);
});

elements.form.addEventListener("submit", submitJob);

bindDropzone();
refreshHealth();
refreshJobs();
window.setInterval(refreshHealth, 15000);
window.setInterval(refreshJobs, 3000);
