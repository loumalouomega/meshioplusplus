<script setup>
import { computed, onMounted, ref, watchEffect } from 'vue'
import { withBase } from 'vitepress'

const runs = ref([])
const state = ref('Loading benchmark history…')
const kind = ref('operations')
const selected = ref('')
const environment = ref('')
const dimensionKeys = {
  formats: ['size', 'format', 'input', 'cells'],
  operations: ['op', 'cells', 'backend', 'threads'],
}
const dimensions = row => dimensionKeys[kind.value].map(key => row[key]).join(' / ')
const matching = computed(() => runs.value.filter(run => run.kind === kind.value))
const environments = computed(() => [...new Set(matching.value.map(run => run.environment_id))])
const activeEnvironment = computed(() => environments.value.includes(environment.value)
  ? environment.value : environments.value.at(-1))
const activeRuns = computed(() => matching.value.filter(run => run.environment_id === activeEnvironment.value))
const cases = computed(() => [...new Set(activeRuns.value.flatMap(run => run.rows.map(dimensions)))].sort())
const activeCase = computed(() => cases.value.includes(selected.value) ? selected.value : cases.value[0])
watchEffect(() => { environment.value = activeEnvironment.value || '' })
watchEffect(() => { selected.value = activeCase.value || '' })
const measurements = computed(() => activeRuns.value.map(run => ({
  run, row: run.rows.find(row => dimensions(row) === activeCase.value),
})))
const metrics = computed(() => kind.value === 'formats' ? ['read_s', 'write_s'] : ['median_s'])
const maxTime = computed(() => Math.max(0.000001, ...measurements.value.flatMap(({ row }) =>
  row ? metrics.value.map(key => Number(row[key])) : [])))
const points = key => measurements.value.map(({ row }, index) => row ?
  `${40 + index * 680 / Math.max(1, measurements.value.length - 1)},${220 - Number(row[key]) / maxTime.value * 190}` : null)
// Break lines at missing rows: a skipped benchmark isn't a zero or interpolation.
const segments = key => {
  const result = []; let current = []
  for (const point of points(key)) {
    if (point) current.push(point)
    else if (current.length) { result.push(current.join(' ')); current = [] }
  }
  if (current.length) result.push(current.join(' '))
  return result
}
const raw = (run, extension) => withBase(`/benchmark-trends/runs/${run.run_id}-${run.attempt}/${run.kind}.${extension}`)
const workflow = run => `https://github.com/loumalouomega/meshioplusplus/actions/runs/${run.run_id}/attempts/${run.attempt}`

onMounted(async () => {
  try {
    const response = await fetch(withBase('/benchmark-trends/index.json'))
    if (!response.ok) throw new Error('History is not published yet.')
    const data = await response.json()
    if (data.schema !== 1 || !Array.isArray(data.runs)) throw new Error('Unsupported history schema.')
    runs.value = data.runs
    state.value = data.runs.length ? '' : 'No benchmark runs have been published yet.'
  } catch (error) { state.value = error.message }
})
</script>

<template>
  <p v-if="state" role="status">{{ state }}</p>
  <section v-else aria-label="Benchmark trends">
    <div class="trend-filters">
      <label>Suite <select v-model="kind"><option>operations</option><option>formats</option></select></label>
      <label>Environment <select v-model="environment">
        <option v-for="id in environments" :key="id" :value="id">{{ id }}</option>
      </select></label>
      <label>Case (size / backend / threads) <select v-model="selected">
        <option v-for="item in cases" :key="item">{{ item }}</option>
      </select></label>
    </div>
    <p v-if="!activeRuns.length">No records for this suite.</p>
    <template v-else>
      <details><summary>Machine, compiler and library versions</summary>
        <pre>{{ JSON.stringify(activeRuns.at(-1).environment, null, 2) }}</pre>
      </details>
      <p>Seconds (lower is better), maximum {{ maxTime.toPrecision(4) }} s. Points are successive runs; gaps are missing measurements.</p>
      <svg viewBox="0 0 760 250" role="img" :aria-label="`Timing history for ${activeCase}`">
        <line x1="40" y1="30" x2="40" y2="220" stroke="currentColor" />
        <line x1="40" y1="220" x2="720" y2="220" stroke="currentColor" />
        <g v-for="(key, color) in metrics" :key="key" :stroke="color ? '#e07b39' : '#42b883'">
          <polyline v-for="(segment, i) in segments(key)" :key="i" :points="segment" fill="none" stroke-width="2" />
          <template v-for="(point, i) in points(key)" :key="i">
            <circle v-if="point" :cx="point.split(',')[0]" :cy="point.split(',')[1]" r="3" :fill="color ? '#e07b39' : '#42b883'">
              <title>{{ key }}: {{ measurements[i].row[key] }} s, {{ measurements[i].run.timestamp }}</title>
            </circle>
          </template>
        </g>
      </svg>
      <p>{{ metrics.join(' (green), ') }}{{ metrics.length > 1 ? ' (orange)' : ' (green)' }}</p>
      <table>
        <thead><tr><th>Run</th><th>Commit</th><th v-for="key in metrics" :key="key">{{ key }}</th><th>Raw</th></tr></thead>
        <tbody><tr v-for="({ run, row }) in measurements" :key="`${run.run_id}-${run.attempt}`">
          <td><a :href="workflow(run)">{{ run.timestamp.slice(0, 10) }} / {{ run.attempt }}</a></td>
          <td><a :href="`https://github.com/loumalouomega/meshioplusplus/commit/${run.sha}`">{{ run.sha.slice(0, 8) }}</a></td>
          <td v-for="key in metrics" :key="key">{{ row ? row[key] : 'missing' }}</td>
          <td><a :href="raw(run, 'csv')">CSV</a> · <a :href="raw(run, 'json')">metadata</a></td>
        </tr></tbody>
      </table>
    </template>
  </section>
</template>

<style scoped>
.trend-filters { display: flex; flex-wrap: wrap; gap: 1rem; margin-bottom: 1rem; }
label { display: flex; flex-direction: column; max-width: 100%; }
select { border: 1px solid var(--vp-c-divider); border-radius: 4px; padding: .4rem; max-width: 100%; }
svg { width: 100%; }
pre { overflow: auto; font-size: .8rem; }
</style>
