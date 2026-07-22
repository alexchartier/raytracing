const form = document.getElementById("trace-form");
const presetButton = document.getElementById("preset-button");
const traceButton = document.getElementById("trace-button");
const statusBanner = document.getElementById("status-banner");
const summaryCopy = document.getElementById("summary-copy");
const rawOutput = document.getElementById("raw-output");

const metricReachable = document.getElementById("metric-reachable");
const metricAbsorption = document.getElementById("metric-absorption");
const metricError = document.getElementById("metric-error");
const metricLaunch = document.getElementById("metric-launch");

function setStatus(kind, text) {
  statusBanner.className = `status-banner status-${kind}`;
  statusBanner.textContent = text;
}

function formatNumber(value, digits = 2) {
  if (value === null || value === undefined || Number.isNaN(value)) {
    return "--";
  }
  return Number(value).toFixed(digits);
}

function formPayload(formElement) {
  const data = new FormData(formElement);
  return Object.fromEntries(data.entries());
}

function applyPreset() {
  form.elements.date.value = "2020-01-15";
  form.elements.time.value = "12:00";
  form.elements.frequency_mhz.value = "4.1";
  form.elements.tx_lat.value = "-77.8";
  form.elements.tx_lon.value = "166.4";
  form.elements.tx_alt_km.value = "0.0";
  form.elements.rx_lat.value = "-89.9";
  form.elements.rx_lon.value = "166.4";
  form.elements.rx_alt_km.value = "1.0";
  form.elements.f107.value = "";
}

async function runTrace(event) {
  event.preventDefault();
  traceButton.disabled = true;
  setStatus("running", "Tracing candidate paths through the ionosphere...");
  summaryCopy.textContent = "Running a 2-hop search through a 60-500 km model volume.";

  try {
    const response = await fetch("/api/trace", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(formPayload(form)),
    });
    const payload = await response.json();
    rawOutput.textContent = JSON.stringify(payload, null, 2);

    if (!response.ok) {
      setStatus("failure", payload.error || "Trace failed.");
      summaryCopy.textContent = "The server rejected the request or the solver failed before a result was produced.";
      metricReachable.textContent = "Unavailable";
      metricAbsorption.textContent = "--";
      metricError.textContent = "--";
      metricLaunch.textContent = "--";
      return;
    }

    const result = payload.result;
    metricReachable.textContent = result.reachable ? "Reachable" : "Not closed";
    metricAbsorption.textContent = `${formatNumber(result.absorption_db, 2)} dB`;
    metricError.textContent = `${formatNumber(result.error_m, 0)} m`;
    metricLaunch.textContent = `${formatNumber(result.launch_elevation_deg, 1)}° / ${formatNumber(result.launch_bearing_deg, 1)}°`;
    summaryCopy.textContent = result.message;
    setStatus(result.reachable ? "success" : "failure", result.message);
  } catch (error) {
    rawOutput.textContent = JSON.stringify({ error: String(error) }, null, 2);
    setStatus("failure", "Request failed before the raytrace completed.");
    summaryCopy.textContent = "The browser could not retrieve a result from the local raytrace service.";
  } finally {
    traceButton.disabled = false;
  }
}

presetButton.addEventListener("click", applyPreset);
form.addEventListener("submit", runTrace);
applyPreset();
