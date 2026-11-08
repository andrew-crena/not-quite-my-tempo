/* Not Quite My Tempo: signal processing.
   Pure functions and small stateful analyzers, no DOM. Loaded by index.html
   and by tests/dsp.test.js (Node), so every number on screen can be checked
   against synthetic audio. */
(function (root) {
'use strict';

/* ---------- radix-2 complex FFT (in place) ---------- */
function makeFFT(n) {
  const levels = Math.round(Math.log2(n));
  const cos = new Float32Array(n / 2), sin = new Float32Array(n / 2);
  for (let i = 0; i < n / 2; i++) {
    cos[i] = Math.cos(2 * Math.PI * i / n);
    sin[i] = Math.sin(2 * Math.PI * i / n);
  }
  const rev = new Uint32Array(n);
  for (let i = 0; i < n; i++) {
    let r = 0;
    for (let j = 0; j < levels; j++) r = (r << 1) | ((i >>> j) & 1);
    rev[i] = r;
  }
  return function fft(re, im) {
    for (let i = 0; i < n; i++) {
      const j = rev[i];
      if (j > i) {
        let t = re[i]; re[i] = re[j]; re[j] = t;
        t = im[i]; im[i] = im[j]; im[j] = t;
      }
    }
    for (let size = 2; size <= n; size <<= 1) {
      const half = size >> 1, step = n / size;
      for (let i = 0; i < n; i += size) {
        for (let j = i, k = 0; j < i + half; j++, k += step) {
          const l = j + half;
          const tre = re[l] * cos[k] + im[l] * sin[k];
          const tim = im[l] * cos[k] - re[l] * sin[k];
          re[l] = re[j] - tre; im[l] = im[j] - tim;
          re[j] += tre; im[j] += tim;
        }
      }
    }
  };
}

function median(a) {
  const s = [...a].sort((x, y) => x - y);
  const m = s.length >> 1;
  return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2;
}

function hann(n) {
  const w = new Float32Array(n);
  for (let i = 0; i < n; i++) w[i] = 0.5 - 0.5 * Math.cos(2 * Math.PI * i / n);
  return w;
}

/* ---------- biquads ---------- */
function Biquad(b0, b1, b2, a0, a1, a2) {
  this.b0 = b0 / a0; this.b1 = b1 / a0; this.b2 = b2 / a0;
  this.a1 = a1 / a0; this.a2 = a2 / a0;
  this.z1 = 0; this.z2 = 0;
}
Biquad.prototype.process = function (x) {   // transposed direct form II
  const y = this.b0 * x + this.z1;
  this.z1 = this.b1 * x - this.a1 * y + this.z2;
  this.z2 = this.b2 * x - this.a2 * y;
  return y;
};
Biquad.prototype.reset = function () { this.z1 = 0; this.z2 = 0; };

function lowpass(fc, q, fs) {             // RBJ cookbook
  const w = 2 * Math.PI * fc / fs, c = Math.cos(w), al = Math.sin(w) / (2 * q);
  return new Biquad((1 - c) / 2, 1 - c, (1 - c) / 2, 1 + al, -2 * c, 1 - al);
}

/* ITU-R BS.1770-4 K-weighting, designed for any sample rate. Parameters are
   the analog prototype that reproduces the standard's 48 kHz coefficients. */
function kWeighting(fs) {
  let G = 3.999843853973347, Q = 0.7071752369554196, fc = 1681.974450955533;
  let K = Math.tan(Math.PI * fc / fs);
  const Vh = Math.pow(10, G / 20), Vb = Math.pow(Vh, 0.4996667741545416);
  const a0 = 1 + K / Q + K * K;
  const shelf = new Biquad(
    Vh + Vb * K / Q + K * K, 2 * (K * K - Vh), Vh - Vb * K / Q + K * K,
    a0, 2 * (K * K - 1), 1 - K / Q + K * K);
  Q = 0.5003270373238773; fc = 38.13547087602444;
  K = Math.tan(Math.PI * fc / fs);
  const h0 = 1 + K / Q + K * K;            // b stays [1,-2,1] as in the standard
  const hp = new Biquad(1, -2, 1, 1, 2 * (K * K - 1) / h0, (1 - K / Q + K * K) / h0);
  return [shelf, hp];
}

/* Anti-alias + decimate. Two Butterworth sections (4th order) at 0.42 of the
   output rate, replacing a plain box average that let 4-6 kHz alias down. */
function Decimator(fs, factor) {
  this.factor = factor;
  this.n = 0;
  const out = fs / factor, fc = 0.42 * out;
  this.f = factor > 1 ? [lowpass(fc, 0.5411961, fs), lowpass(fc, 1.3065630, fs)] : [];
}
Decimator.prototype.push = function (x) {       // returns NaN unless a sample is due
  for (let i = 0; i < this.f.length; i++) x = this.f[i].process(x);
  if (++this.n >= this.factor) { this.n = 0; return x; }
  return NaN;
};

/* ---------- loudness (BS.1770 / EBU R128, single channel) ----------
   Mono input is metered as one channel: a full-scale 997 Hz sine reads
   -3.01 LUFS. Momentary = 400 ms, short-term = 3 s, integrated = gated
   (absolute -70 LUFS, relative -10 LU) over 400 ms blocks with 75 % overlap. */
function LoudnessMeter(fs) {
  this.fs = fs;
  this.k = kWeighting(fs);
  this.sub = Math.round(fs / 10);            // 100 ms sub-blocks
  this.acc = 0; this.n = 0;
  this.subs = new Float64Array(30); this.w = 0;  // last 3 s of sub-block mean squares
  this.blocks = [];                         // 400 ms block mean squares (for integrated)
  this.maxBlocks = 36000;                   // 1 hour
}
LoudnessMeter.prototype.push = function (x) {
  const y = this.k[1].process(this.k[0].process(x));
  this.acc += y * y;
  if (++this.n >= this.sub) {
    this.subs[this.w % 30] = this.acc / this.n; this.w++;
    this.acc = 0; this.n = 0;
    if (this.w >= 4) {
      this.blocks.push(this.mean(4));
      if (this.blocks.length > this.maxBlocks) this.blocks.splice(0, 6000);
    }
  }
};
LoudnessMeter.prototype.mean = function (k) {
  const m = Math.min(k, this.w);
  if (!m) return 0;
  let s = 0;
  for (let i = 0; i < m; i++) s += this.subs[(this.w - 1 - i) % 30];
  return s / m;
};
function lufs(ms) { return ms > 1e-12 ? -0.691 + 10 * Math.log10(ms) : -Infinity; }
LoudnessMeter.prototype.momentary = function () { return this.w >= 4 ? lufs(this.mean(4)) : -Infinity; };
LoudnessMeter.prototype.shortTerm = function () { return this.w >= 30 ? lufs(this.mean(30)) : (this.w >= 4 ? lufs(this.mean(this.w)) : -Infinity); };
LoudnessMeter.prototype.integrated = function () {
  const abs = Math.pow(10, (-70 + 0.691) / 10);
  let s = 0, c = 0;
  for (const b of this.blocks) if (b > abs) { s += b; c++; }
  if (!c) return -Infinity;
  const rel = (s / c) * 0.1;                 // -10 LU
  s = 0; c = 0;
  for (const b of this.blocks) if (b > abs && b > rel) { s += b; c++; }
  return c ? lufs(s / c) : -Infinity;
};
LoudnessMeter.prototype.reset = function () {
  this.k.forEach(f => f.reset());
  this.acc = 0; this.n = 0; this.subs.fill(0); this.w = 0; this.blocks = [];
};

/* ---------- onset envelope: log-magnitude spectral flux ---------- */
function FluxTracker(fs, N, hop) {
  this.fs = fs; this.N = N; this.hop = hop;
  this.rate = fs / hop;
  this.fft = makeFFT(N); this.win = hann(N);
  this.re = new Float32Array(N); this.im = new Float32Array(N);
  this.prev = new Float32Array(N / 2 + 1);
  this.buf = new Float32Array(N); this.fill = 0;
  this.env = [];
  this.maxLen = Math.floor(this.rate * 20);
}
FluxTracker.prototype.push = function (x) {
  this.buf[this.fill++] = x;
  if (this.fill < this.N) return;
  const N = this.N, re = this.re, im = this.im;
  for (let i = 0; i < N; i++) { re[i] = this.buf[i] * this.win[i]; im[i] = 0; }
  this.fft(re, im);
  let flux = 0;
  for (let b = 1; b <= N / 2; b++) {
    const m = Math.log1p(100 * Math.sqrt(re[b] * re[b] + im[b] * im[b]));
    const d = m - this.prev[b];
    if (d > 0) flux += d;
    this.prev[b] = m;
  }
  this.env.push(flux);
  if (this.env.length > this.maxLen + 256) this.env.splice(0, this.env.length - this.maxLen);
  this.buf.copyWithin(0, this.hop); this.fill = N - this.hop;
};
FluxTracker.prototype.reset = function () {
  this.prev.fill(0); this.fill = 0; this.env = [];
};

/* Adaptive-threshold the envelope: subtract a ~0.25 s moving average and
   half-wave rectify, so slow loudness swells do not masquerade as tempo. */
function sharpenEnvelope(env, rate) {
  const n = env.length, h = Math.max(2, Math.round(rate * 0.125));
  const out = new Float32Array(n);
  const pre = new Float64Array(n + 1);
  for (let i = 0; i < n; i++) pre[i + 1] = pre[i] + env[i];
  for (let i = 0; i < n; i++) {
    const a = Math.max(0, i - h), b = Math.min(n, i + h + 1);
    const v = env[i] - (pre[b] - pre[a]) / (b - a);
    out[i] = v > 0 ? v : 0;
  }
  return out;
}

/* ---------- tempo: autocorrelation of the onset envelope ---------- */
const BPM_MIN = 60, BPM_MAX = 200;
function tempoFromEnv(envRaw, envRate, anchorBpm) {
  if (envRaw.length < Math.floor(envRate * 4)) return null;   // need >= 4 s
  const A = Math.min(envRaw.length, Math.floor(envRate * 16)); // last 16 s
  const sh = sharpenEnvelope(envRaw.slice(envRaw.length - A), envRate);
  // Light smoothing widens razor-thin onset peaks so fractional lags are not
  // missed between integer ACF samples.
  const env = new Float32Array(A);
  for (let i = 0; i < A; i++) env[i] = 0.5 * sh[i] + 0.25 * ((sh[i - 1] || 0) + (sh[i + 1] || 0));
  const seg = new Float32Array(A), wts = new Float32Array(A);
  const tau = envRate * 7;
  let mean = 0, wsum = 0;
  for (let i = 0; i < A; i++) {
    const w = Math.exp(-(A - 1 - i) / tau);
    wts[i] = w; seg[i] = env[i];
    mean += seg[i] * w; wsum += w;
  }
  mean /= wsum;
  let e0 = 0;
  for (let i = 0; i < A; i++) { seg[i] -= mean; e0 += seg[i] * seg[i] * wts[i]; }
  if (e0 < 1e-9) return null;

  const Lmin = Math.max(4, Math.floor(envRate * 60 / BPM_MAX));
  const Lmax = Math.min(A - 4, Math.ceil(envRate * 60 / BPM_MIN));
  if (Lmax <= Lmin + 2) return null;
  const maxLag = Math.min(A - 1, Lmax * 4 + 2);

  const acf = new Float32Array(maxLag + 2);
  for (let L = 1; L <= maxLag; L++) {
    let s = 0;
    for (let i = 0; i < A - L; i++) s += seg[i] * seg[i + L] * wts[i + L];
    acf[L] = (s / e0) * (A / (A - L));
  }
  // Interpolated ACF read so multiples of a fractional lag land correctly.
  const at = x => {
    const i = Math.floor(x), f = x - i;
    return i + 1 <= maxLag ? acf[i] * (1 - f) + acf[i + 1] * f : 0;
  };

  // Comb score against a wide, log-normal tempo prior centred at 120 BPM.
  const score = new Float32Array(Lmax + 2);
  for (let L = Lmin; L <= Lmax; L++) {
    const bpm = 60 * envRate / L;
    const prior = Math.exp(-0.5 * Math.pow(Math.log2(bpm / 120) / 1.0, 2));
    let s = acf[L];
    if (2 * L <= maxLag) s += 0.5 * at(2 * L);
    if (3 * L <= maxLag) s += 0.33 * at(3 * L);
    if (4 * L <= maxLag) s += 0.25 * at(4 * L);
    // Sub-multiple guard: if the half lag (a faster pulse still in range)
    // correlates almost as well, this lag is probably two beats, not one.
    let sub = 1;
    if (L / 2 >= Lmin && acf[L] > 0) {
      const h = at(L / 2) / acf[L];
      if (h > 0.75) sub = 1 - 0.35 * Math.min(1, (h - 0.75) / 0.2);
    }
    score[L] = s * prior * sub;
  }

  const peaks = [];
  for (let L = Lmin + 1; L < Lmax; L++) {
    if (score[L] >= score[L - 1] && score[L] >= score[L + 1] && score[L] > 0) peaks.push({ L, s: score[L] });
  }
  if (!peaks.length) return null;
  peaks.sort((a, b) => b.s - a.s);
  const top = peaks.slice(0, 6);
  const topMax = top[0].s;

  // Anchor: bounded tie-breaker toward the current tempo (or a simple ratio),
  // capped so real evidence still overrides a stale lock.
  let bestL = top[0].L, bestAdj = topMax;
  if (anchorBpm) {
    const anchorL = 60 * envRate / anchorBpm;
    for (const p of top) {
      let d = Infinity;
      for (const r of [1, 2, 0.5, 1.5, 1 / 1.5]) d = Math.min(d, Math.abs(Math.log2(p.L / (anchorL * r))));
      const adj = p.s * (1 + 0.6 * Math.exp(-0.5 * Math.pow(d / 0.03, 2)));
      if (adj > bestAdj) { bestAdj = adj; bestL = p.L; }
    }
  }

  // Sub-lag refinement: parabolic peak, then confirm with the 4th multiple,
  // which carries 4x the timing resolution.
  let Lstar = bestL;
  if (bestL > Lmin && bestL < Lmax) {
    const a = score[bestL - 1], b = score[bestL], c = score[bestL + 1];
    const den = a - 2 * b + c;
    if (den < 0) { const d = 0.5 * (a - c) / den; if (d > -1 && d < 1) Lstar = bestL + d; }
  }
  const M = 4 * Lstar;
  if (M + 2 <= maxLag) {
    const i = Math.round(M);
    let bi = i;
    for (let j = i - 2; j <= i + 2; j++) if (acf[j] > acf[bi]) bi = j;
    if (bi > 0 && bi < maxLag) {
      const a = acf[bi - 1], b = acf[bi], c = acf[bi + 1], den = a - 2 * b + c;
      const d = den < 0 ? 0.5 * (a - c) / den : 0;
      const refined = (bi + Math.max(-1, Math.min(1, d))) / 4;
      if (Math.abs(refined - Lstar) < 0.6) Lstar = refined;
    }
  }
  // Confidence: comb peak height, discounted when the runner-up (excluding
  // octave relatives) is nearly as strong.
  let rival = 0;
  for (const p of top) {
    const r = p.L / bestL;
    const rel = [0.5, 1, 2, 1.5, 2 / 3, 3, 1 / 3].some(k => Math.abs(Math.log2(r / k)) < 0.04);
    if (!rel && p.s > rival) rival = p.s;
  }
  const peak = score[bestL];
  const conf = Math.max(0, Math.min(1, (peak / 1.1) * (1 - 0.5 * Math.min(1, rival / Math.max(peak, 1e-9)))));
  return { bpm: 60 * envRate / Lstar, conf };
}

/* Beat phase: where the most recent beat sits in the envelope. Returns the
   number of envelope frames from the newest frame back to that beat. */
function beatPhase(envRaw, envRate, bpm, beats) {
  beats = beats || 8;
  const period = envRate * 60 / bpm;
  const W = Math.floor(period * beats);
  if (period < 2 || envRaw.length < W + 2) return null;
  const env = sharpenEnvelope(envRaw.slice(envRaw.length - W - 2), envRate);
  const n = env.length;
  const P = Math.floor(period);
  let bestOff = 0, best = -1;
  for (let off = 0; off < P; off++) {           // off = frames back from newest
    let s = 0;
    for (let k = 0; k < beats; k++) {
      const x = n - 1 - off - k * period;
      const i = Math.floor(x), f = x - i;
      if (i < 0) break;
      s += (env[i] * (1 - f) + (env[i + 1] || 0) * f) * (1 - 0.06 * k);  // favour recent beats
    }
    if (s > best) { best = s; bestOff = off; }
  }
  return { framesBack: bestOff, period, strength: best };
}

/* ---------- pitch: YIN (de Cheveigné & Kawahara 2002) ---------- */
function yinPitch(buf, sr, fmin, fmax, threshold) {
  const TH = threshold || 0.12;
  const W = Math.floor(buf.length / 2);
  const tauMin = Math.max(2, Math.floor(sr / fmax));
  const tauMax = Math.min(Math.ceil(sr / fmin), W - 1);
  if (tauMax <= tauMin + 2) return null;
  const d = new Float32Array(tauMax + 2);
  for (let t = 1; t <= tauMax + 1; t++) {
    let s = 0;
    for (let i = 0; i < W; i++) { const df = buf[i] - buf[i + t]; s += df * df; }
    d[t] = s;
  }
  const cm = new Float32Array(tauMax + 2);
  cm[0] = 1;
  let run = 0;
  for (let t = 1; t <= tauMax + 1; t++) { run += d[t]; cm[t] = run > 0 ? d[t] * t / run : 1; }

  let tau = -1;
  for (let t = tauMin; t <= tauMax; t++) {
    if (cm[t] < TH) { while (t + 1 <= tauMax && cm[t + 1] < cm[t]) t++; tau = t; break; }
  }
  if (tau < 0) {
    // No dip under the threshold (noisy or decaying note). The global minimum
    // is often a multiple of the period, so take the earliest dip that comes
    // close to it instead.
    let mv = Infinity;
    for (let t = tauMin; t <= tauMax; t++) if (cm[t] < mv) mv = cm[t];
    if (mv >= 0.3) return null;
    const lim = mv * 1.2 + 0.03;
    for (let t = tauMin; t <= tauMax; t++) {
      if (cm[t] < lim) { while (t + 1 <= tauMax && cm[t + 1] < cm[t]) t++; tau = t; break; }
    }
  }
  // A period at the long end of the search range is the range edge, not a
  // note (noise, rumble or hum residue); report nothing rather than a guess.
  if (tau >= tauMax - 1) return null;
  // Octave guard: if half the period is also a deep dip, the fundamental is there.
  const half = Math.round(tau / 2);
  if (half >= tauMin) {
    let h = half;
    for (let j = half - 1; j <= half + 1; j++) if (j >= tauMin && cm[j] < cm[h]) h = j;
    if (cm[h] < Math.max(TH, cm[tau] * 1.25) && cm[h] < 0.3) tau = h;
  }
  // Parabolic interpolation on the raw difference function (less biased than on cm).
  let better = tau;
  if (tau > 1 && tau <= tauMax) {
    const a = d[tau - 1], b = d[tau], c = d[tau + 1], den = a - 2 * b + c;
    if (den > 0) { const dd = 0.5 * (a - c) / den; if (dd > -1 && dd < 1) better = tau + dd; }
  }
  return { freq: sr / better, clarity: 1 - cm[tau] };
}

/* Refine a coarse frequency on a full-rate buffer: evaluate the difference
   function only around the expected period and interpolate. Restores the
   sub-cent precision that decimation costs on high notes. */
function refinePitch(buf, sr, fCoarse, lpFloor) {
  // Low-pass a copy a few harmonics above the note; skip the filter's start-up.
  const fc = Math.min(0.45 * sr, Math.max(lpFloor || 2500, 6 * fCoarse));
  const l1 = lowpass(fc, 0.5411961, sr), l2 = lowpass(fc, 1.3065630, sr);
  const SKIP = 256;
  const n = buf.length - SKIP;
  if (n <= 0) return fCoarse;
  const x = new Float32Array(n);
  for (let i = 0; i < buf.length; i++) { const y = l2.process(l1.process(buf[i])); if (i >= SKIP) x[i - SKIP] = y; }
  buf = x;
  const T = sr / fCoarse, W = Math.min(Math.floor(buf.length - T - 8), 4096);
  if (W < 256) return fCoarse;
  const t0 = Math.max(2, Math.floor(T) - 3), t1 = Math.ceil(T) + 3;
  let bi = -1, bv = Infinity;
  const d = {};
  for (let t = t0; t <= t1; t++) {
    let s = 0;
    for (let i = 0; i < W; i++) { const df = buf[i] - buf[i + t]; s += df * df; }
    d[t] = s;
    if (s < bv) { bv = s; bi = t; }
  }
  if (bi <= t0 || bi >= t1) return fCoarse;
  const a = d[bi - 1], b = d[bi], c = d[bi + 1], den = a - 2 * b + c;
  const dd = den > 0 ? 0.5 * (a - c) / den : 0;
  const f = sr / (bi + dd);
  return Math.abs(1200 * Math.log2(f / fCoarse)) < 40 ? f : fCoarse;
}

/* Tuner pre-filter, run continuously on the decimated stream so the narrow
   notches are settled (applying them per analysis window would ring).
   Removes DC/rumble below the lowest target and 50/60 Hz hum plus their
   second harmonics, which otherwise beat with low strings and capture the
   pitch search once a note decays. */
function notch(f0, q, fs) {               // RBJ cookbook
  const w = 2 * Math.PI * f0 / fs, c = Math.cos(w), al = Math.sin(w) / (2 * q);
  return new Biquad(1, -2 * c, 1, 1 + al, -2 * c, 1 - al);
}
function highpass(fc, q, fs) {            // RBJ cookbook
  const w = 2 * Math.PI * fc / fs, c = Math.cos(w), al = Math.sin(w) / (2 * q);
  return new Biquad((1 + c) / 2, -(1 + c), (1 + c) / 2, 1 + al, -2 * c, 1 - al);
}
function TunerPrefilter(fs) {
  this.fs = fs;
  this.notches = [50, 60, 100, 120].map(f => notch(f, 25, fs));
  this.setHighpass(20);
  this.setLowpass(0);
}
/* Optional low-pass a few harmonics above the highest target: stiff strings
   (bass especially) have upper partials that run sharp, and a period
   detector weighs them all. 0 turns it off. */
TunerPrefilter.prototype.setLowpass = function (fc) {
  if (this.lc === fc) return;
  this.lc = fc;
  this.lp = fc > 0 && fc < 0.45 * this.fs ? [lowpass(fc, 0.5411961, this.fs), lowpass(fc, 1.3065630, this.fs)] : [];
};
TunerPrefilter.prototype.setHighpass = function (fc) {
  if (this.fc === fc) return;
  this.fc = fc;
  this.hp = [highpass(fc, 0.5411961, this.fs), highpass(fc, 1.3065630, this.fs)];
};
TunerPrefilter.prototype.push = function (x) {
  for (let i = 0; i < 4; i++) x = this.notches[i].process(x);
  for (let i = 0; i < this.lp.length; i++) x = this.lp[i].process(x);
  return this.hp[1].process(this.hp[0].process(x));
};
TunerPrefilter.prototype.reset = function () {
  for (const b of this.notches) b.reset();
  for (const b of this.hp) b.reset();
  for (const b of this.lp) b.reset();
};

/* Combine recent pitch readings: the median rejects octave slips, the mean
   of readings near it smooths jitter. Returns null unless most readings
   agree. (With an even split the median falls between two clusters and
   nothing is near it; averaging an empty set used to give NaN.) */
function consensusPitch(freqs, tolCents) {
  if (!freqs.length) return null;
  const med = median(freqs), tol = tolCents || 15;
  let s = 0, n = 0;
  for (const f of freqs) if (Math.abs(1200 * Math.log2(f / med)) < tol) { s += f; n++; }
  return n * 2 > freqs.length ? s / n : null;
}

/* ---------- notes, tunings ---------- */
const NOTE_NAMES = ['C','C#','D','D#','E','F','F#','G','G#','A','A#','B'];
const FLAT_NAMES = ['C','Db','D','Eb','E','F','Gb','G','Ab','A','Bb','B'];
function midiToFreq(m, a4) { return (a4 || 440) * Math.pow(2, (m - 69) / 12); }
function freqToMidi(f, a4) { return 69 + 12 * Math.log2(f / (a4 || 440)); }
function noteName(m, flats) { const pc = ((m % 12) + 12) % 12; return (flats ? FLAT_NAMES : NOTE_NAMES)[pc]; }
function noteOctave(m) { return Math.floor(m / 12) - 1; }

const TUNINGS = {
  standard:  { label: 'Guitar · Standard (E A D G B E)', midi: [40, 45, 50, 55, 59, 64] },
  dropd:     { label: 'Guitar · Drop D (D A D G B E)',  midi: [38, 45, 50, 55, 59, 64] },
  halfdown:  { label: 'Guitar · Half step down (E♭)',   midi: [39, 44, 49, 54, 58, 63], flats: true },
  fulldown:  { label: 'Guitar · Whole step down (D)',   midi: [38, 43, 48, 53, 57, 62] },
  dadgad:    { label: 'Guitar · DADGAD',                midi: [38, 45, 50, 55, 57, 62] },
  openg:     { label: 'Guitar · Open G (D G D G B D)',  midi: [38, 43, 50, 55, 59, 62] },
  opend:     { label: 'Guitar · Open D (D A D F♯ A D)', midi: [38, 45, 50, 54, 57, 62] },
  bass4:     { label: 'Bass · 4-string (E A D G)',      midi: [28, 33, 38, 43] },
  bass5:     { label: 'Bass · 5-string (B E A D G)',    midi: [23, 28, 33, 38, 43] },
  ukulele:   { label: 'Ukulele · Standard (G C E A)',   midi: [67, 60, 64, 69] },
  chromatic: { label: 'Chromatic · any note',           midi: null },
};

/* Pitch search range for a tuning: a little below the lowest string and
   high enough for the top string's fundamental. */
function tuningRange(t, a4) {
  if (!t.midi) return { fmin: 27, fmax: 1400 };
  const lo = Math.min(...t.midi), hi = Math.max(...t.midi);
  return { fmin: midiToFreq(lo, a4) * 0.75, fmax: midiToFreq(hi, a4) * 1.45 };
}

/* Target for a detected frequency. Pinned string wins; otherwise nearest
   string (guitar) or nearest semitone (chromatic). */
function tunerTarget(f, t, a4, pinned) {
  const m = freqToMidi(f, a4);
  if (!t.midi) {
    const r = Math.round(m);
    return { idx: -1, midi: r, f: midiToFreq(r, a4), cents: 100 * (m - r) };
  }
  let best = -1, bd = Infinity;
  for (let i = 0; i < t.midi.length; i++) {
    const dd = Math.abs(m - t.midi[i]);
    if (dd < bd) { bd = dd; best = i; }
  }
  const idx = pinned != null && pinned >= 0 && pinned < t.midi.length ? pinned : best;
  const tm = t.midi[idx];
  return { idx, midi: tm, f: midiToFreq(tm, a4), cents: 100 * (m - tm) };
}

/* ---------- key: Krumhansl-Kessler profiles + Camelot ---------- */
const MAJ_PROF = [6.35,2.23,3.48,2.33,4.38,4.09,2.52,5.19,2.39,3.66,2.29,2.88];
const MIN_PROF = [6.33,2.68,3.52,5.38,2.60,3.53,2.54,4.75,3.98,2.69,3.34,3.17];
const MAJ_CAM = ['8B','3B','10B','5B','12B','7B','2B','9B','4B','11B','6B','1B'];
const MIN_CAM = ['5A','12A','7A','2A','9A','4A','11A','6A','1A','8A','3A','10A'];
const MAJ_NAMES = ['C','D♭','D','E♭','E','F','F♯','G','A♭','A','B♭','B'];
const MIN_NAMES = ['C','C♯','D','E♭','E','F','F♯','G','G♯','A','B♭','B'];

function pearson(x, y) {
  const n = x.length;
  let mx = 0, my = 0;
  for (let i = 0; i < n; i++) { mx += x[i]; my += y[i]; }
  mx /= n; my /= n;
  let sxy = 0, sxx = 0, syy = 0;
  for (let i = 0; i < n; i++) { const a = x[i] - mx, b = y[i] - my; sxy += a * b; sxx += a * a; syy += b * b; }
  const d = Math.sqrt(sxx * syy);
  return d > 0 ? sxy / d : 0;
}

function keyFromChroma(c) {
  let sum = 0;
  for (let i = 0; i < 12; i++) sum += c[i];
  if (sum < 1e-6) return null;
  const results = [];
  for (let root = 0; root < 12; root++) {
    const maj = new Array(12), min = new Array(12);
    for (let p = 0; p < 12; p++) { maj[p] = MAJ_PROF[(p - root + 12) % 12]; min[p] = MIN_PROF[(p - root + 12) % 12]; }
    results.push({ root, minor: false, r: pearson(c, maj) });
    results.push({ root, minor: true, r: pearson(c, min) });
  }
  results.sort((a, b) => b.r - a.r);
  const best = results[0];
  return {
    id: best.root + (best.minor ? 12 : 0),
    name: (best.minor ? MIN_NAMES : MAJ_NAMES)[best.root] + (best.minor ? ' minor' : ' major'),
    camelot: best.minor ? MIN_CAM[best.root] : MAJ_CAM[best.root],
    r: best.r,
    margin: best.r - results[1].r,
  };
}

/* Camelot neighbours: one step down, relative major/minor, one step up. */
function compatibleKeys(camelot) {
  const num = parseInt(camelot, 10);
  const letter = camelot.slice(-1);
  if (isNaN(num)) return [];
  const up = (num % 12) + 1;
  const dn = ((num + 10) % 12) + 1;
  return [dn + letter, num + (letter === 'A' ? 'B' : 'A'), up + letter];
}

/* Chroma from one FFT frame: spectral peaks folded to pitch classes, with a
   running tuning-offset estimate so a detuned record still lands correctly. */
function ChromaTracker(fs, N) {
  this.fs = fs; this.N = N;
  this.fft = makeFFT(N); this.win = hann(N);
  this.re = new Float32Array(N); this.im = new Float32Array(N);
  this.acc = new Float32Array(12);
  this.tx = 0; this.ty = 0;
}
ChromaTracker.prototype.frame = function (samples) {
  const N = this.N, re = this.re, im = this.im;
  for (let i = 0; i < N; i++) { re[i] = samples[i] * this.win[i]; im[i] = 0; }
  this.fft(re, im);
  const mag = new Float32Array(N / 2);
  for (let b = 0; b < N / 2; b++) mag[b] = Math.sqrt(re[b] * re[b] + im[b] * im[b]);
  const peaks = [];
  for (let b = 2; b < N / 2 - 1; b++) {
    const m = mag[b];
    if (m > mag[b - 1] && m >= mag[b + 1] && m > 1e-4) {
      const a = mag[b - 1], c = mag[b + 1], den = a - 2 * m + c;
      const d = den < 0 ? 0.5 * (a - c) / den : 0;      // sub-bin peak
      const f = (b + d) * this.fs / N;
      if (f >= 55 && f <= 4200) peaks.push({ p: 12 * Math.log2(f / 440) + 69, w: Math.log1p(m) });
    }
  }
  this.tx *= 0.99; this.ty *= 0.99;
  for (const pk of peaks) {
    const fr = pk.p - Math.round(pk.p);
    this.tx += pk.w * Math.cos(2 * Math.PI * fr);
    this.ty += pk.w * Math.sin(2 * Math.PI * fr);
  }
  const tuning = (this.tx || this.ty) ? Math.atan2(this.ty, this.tx) / (2 * Math.PI) : 0;
  for (let i = 0; i < 12; i++) this.acc[i] *= 0.994;
  for (const pk of peaks) {
    const pc = ((Math.round(pk.p - tuning)) % 12 + 12) % 12;
    this.acc[pc] += pk.w;
  }
};
ChromaTracker.prototype.reset = function () { this.acc.fill(0); this.tx = 0; this.ty = 0; };

/* ---------- tap tempo ---------- */
function TapTempo() { this.t = []; }
TapTempo.prototype.tap = function (now) {
  if (this.t.length && now - this.t[this.t.length - 1] > 2.0) this.t = [];
  this.t.push(now);
  if (this.t.length > 12) this.t.shift();
  if (this.t.length < 3) return null;
  const iv = [];
  for (let i = 1; i < this.t.length; i++) iv.push(this.t[i] - this.t[i - 1]);
  const m = median(iv);
  const kept = iv.filter(x => Math.abs(x - m) < m * 0.2);
  const avg = kept.reduce((a, b) => a + b, 0) / kept.length;
  return { bpm: 60 / avg, taps: this.t.length };
};

root.NQMT_DSP = {
  makeFFT, median, hann, Biquad, lowpass, kWeighting, Decimator,
  LoudnessMeter, FluxTracker, sharpenEnvelope, tempoFromEnv, beatPhase,
  yinPitch, refinePitch, TunerPrefilter, consensusPitch, NOTE_NAMES, midiToFreq, freqToMidi, noteName, noteOctave,
  TUNINGS, tuningRange, tunerTarget, keyFromChroma, compatibleKeys,
  ChromaTracker, TapTempo, BPM_MIN, BPM_MAX,
};
})(typeof globalThis !== 'undefined' ? globalThis : this);
