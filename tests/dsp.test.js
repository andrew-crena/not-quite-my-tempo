// Run: node tests/dsp.test.js
// Synthetic-signal checks for every analyzer in dsp.js.
'use strict';
require('../dsp.js');
const D = globalThis.NQMT_DSP;

let fails = 0, passes = 0;
function check(name, ok, detail) {
  if (ok) { passes++; console.log('  ok   ' + name + (detail ? '  (' + detail + ')' : '')); }
  else { fails++; console.log('  FAIL ' + name + (detail ? '  (' + detail + ')' : '')); }
}
let seed = 12345;
const rnd = () => ((seed = (seed * 1103515245 + 12345) & 0x7fffffff) / 0x7fffffff);
const noise = () => rnd() * 2 - 1;

/* ---- K-weighting matches the BS.1770 48 kHz coefficients ---- */
console.log('loudness');
{
  const [s, h] = D.kWeighting(48000);
  const near = (a, b) => Math.abs(a - b) < 1e-6;
  check('shelf b0', near(s.b0, 1.53512485958697), s.b0.toFixed(8));
  check('shelf a1', near(s.a1, -1.69065929318241), s.a1.toFixed(8));
  check('hp a1', near(h.a1, -1.99004745483398), h.a1.toFixed(8));
  check('hp a2', near(h.a2, 0.99007225036621), h.a2.toFixed(8));
  for (const fs of [44100, 48000]) {
    const m = new D.LoudnessMeter(fs);
    const A = Math.pow(10, -20 / 20);
    for (let i = 0; i < fs * 5; i++) m.push(A * Math.sin(2 * Math.PI * 997 * i / fs));
    const st = m.shortTerm(), it = m.integrated();
    check('997 Hz @ -20 dBFS sine reads -23.0 LUFS short-term @' + fs, Math.abs(st + 23.01) < 0.1, st.toFixed(2));
    check('integrated agrees @' + fs, Math.abs(it + 23.01) < 0.1, it.toFixed(2));
  }
  // Relative gate: 10 s at -20 dBFS then 10 s of near silence should still read ~-23.
  const m = new D.LoudnessMeter(48000);
  for (let i = 0; i < 480000; i++) m.push(0.1 * Math.sin(2 * Math.PI * 997 * i / 48000));
  for (let i = 0; i < 480000; i++) m.push(0.0001 * noise());
  check('gating ignores silence', Math.abs(m.integrated() + 23.01) < 0.2, m.integrated().toFixed(2));
}

/* ---- tempo on click tracks through the real pipeline ---- */
console.log('tempo');
function runTempo(bpm, opts) {
  opts = opts || {};
  const fs = 44100, DEC = 4, dsr = fs / DEC;
  const dec = new D.Decimator(fs, DEC);
  const fl = new D.FluxTracker(dsr, 512, 64);
  const secs = opts.secs || 14, period = 60 / bpm;
  let est = null, anchor = null;
  const total = fs * secs;
  for (let i = 0; i < total; i++) {
    const t = i / fs;
    let x = 0;
    const ph = t % period, beatN = Math.floor(t / period);
    if (ph < 0.03) {                                   // kick: decaying 60 Hz + click
      const e = Math.exp(-ph * 90);
      x += 0.6 * e * Math.sin(2 * Math.PI * 60 * ph) + (ph < 0.002 ? 0.5 * noise() : 0);
    }
    if (opts.hats) {                                   // off-beat hats, softer
      const hp = (t + period / 2) % period;
      if (hp < 0.02) x += 0.15 * Math.exp(-hp * 200) * noise();
    }
    if (opts.swingAccent && beatN % 4 === 0 && ph < 0.03) x *= 1.5;
    x += (opts.noise || 0.01) * noise();
    const y = dec.push(x);
    if (y === y) fl.push(y);
    if (i % (fs / 2) === 0 && i > 0) {
      const e = D.tempoFromEnv(fl.env, fl.rate, anchor);
      if (e) { est = e; anchor = e.bpm; }
    }
  }
  return { est, fl };
}
for (const [bpm, hats] of [[72, false], [90, true], [120, true], [128, true], [140, true], [174, false], [160, false]]) {
  const { est } = runTempo(bpm, { hats });
  const err = est ? Math.abs(est.bpm - bpm) : Infinity;
  check('click track ' + bpm + ' BPM within 0.3', err < 0.3, est ? est.bpm.toFixed(2) + ' conf ' + est.conf.toFixed(2) : 'none');
}
{
  const { est } = runTempo(123.5, { noise: 0.08 });
  check('fractional 123.5 BPM in noise within 0.3', est && Math.abs(est.bpm - 123.5) < 0.3, est && est.bpm.toFixed(2));
}
{
  const { est } = runTempo(120);
  const silent = D.tempoFromEnv(new Array(3000).fill(0), 172, null);
  check('no estimate from silence', silent === null);
  check('confidence high on clean clicks', est && est.conf > 0.4, est && est.conf.toFixed(2));
}
{
  const { fl } = runTempo(120, { secs: 10 });
  const ph = D.beatPhase(fl.env, fl.rate, 120);
  // 10 s at 120 BPM: last beat at t=9.5 s. Envelope frame lag is N/2 hop-frames.
  const newestT = (fl.env.length - 1) * 64 / (44100 / 4) + 256 / (44100 / 4);
  const beatT = newestT - ph.framesBack * 64 / (44100 / 4);
  const err = Math.min(Math.abs(beatT - 9.5), Math.abs(beatT - 10), Math.abs(beatT - 9.0));
  check('beat phase within 25 ms', err < 0.025, (err * 1000).toFixed(1) + ' ms');
}

/* ---- tap tempo ---- */
console.log('tap');
{
  const tt = new D.TapTempo();
  let r;
  for (let i = 0; i < 8; i++) r = tt.tap(i * 0.5 + (i === 4 ? 0.03 : 0));
  check('8 taps at 0.5 s -> 120', r && Math.abs(r.bpm - 120) < 1.5, r && r.bpm.toFixed(2));
  const r2 = tt.tap(10);
  check('long gap restarts', r2 === null);
}

/* ---- YIN on plucked-string-like tones ---- */
console.log('pitch');
function pluck(f, sr, n, opts) {
  opts = opts || {};
  const b = new Float32Array(n);
  for (let i = 0; i < n; i++) {
    const t = i / sr;
    let x = 0;
    const hs = opts.harm || [0.6, 1.0, 0.5, 0.35, 0.2, 0.12];   // strong 2nd harmonic, as on guitar
    hs.forEach((a, k) => { x += a * Math.sin(2 * Math.PI * f * (k + 1) * t + k); });
    b[i] = x * 0.2 + (opts.noise || 0) * noise() + (opts.hum || 0) * Math.sin(2 * Math.PI * 60 * t);
  }
  return b;
}
const sr = 12000;
const std = D.TUNINGS.standard;
for (const m of std.midi) {
  const f = D.midiToFreq(m, 440) * Math.pow(2, 7 / 1200);           // 7 cents sharp
  const r = D.tuningRange(std, 440);
  const p = D.yinPitch(pluck(f, sr, 4096, { noise: 0.01 }), sr, r.fmin, r.fmax);
  const tg = p && D.tunerTarget(p.freq, std, 440, null);
  check('string ' + D.noteName(m) + D.noteOctave(m) + ' +7c', tg && tg.midi === m && Math.abs(tg.cents - 7) < 1,
    tg ? tg.cents.toFixed(2) + 'c' : 'none');
}
{
  const f = D.midiToFreq(28, 440);                                    // bass low E, 41.2 Hz
  const r = D.tuningRange(D.TUNINGS.bass4, 440);
  const p = D.yinPitch(pluck(f, sr, 4096, { harm: [0.3, 1, 0.7, 0.5] }), sr, r.fmin, r.fmax);
  const tg = p && D.tunerTarget(p.freq, D.TUNINGS.bass4, 440, null);
  check('bass E1 weak fundamental', tg && tg.midi === 28 && Math.abs(tg.cents) < 1.5, tg && tg.cents.toFixed(2));
}
{
  const f = 110 * Math.pow(2, -3 / 1200) * 432 / 440;
  const p = D.yinPitch(pluck(f, sr, 4096, { hum: 0.05 }), sr, 60, 500);
  const tg = D.tunerTarget(p.freq, std, 432, null);
  check('A2 at A4=432 with mains hum', tg.midi === 45 && Math.abs(tg.cents + 3) < 1.5, tg.cents.toFixed(2));
}
for (const m of [64, 76, 81, 88]) {             // chromatic, coarse on 12 kHz then refine at 48 kHz
  const f = D.midiToFreq(m, 440) * Math.pow(2, 4 / 1200);
  const full = pluck(f, 48000, 16384), dec = new Float32Array(4096);
  for (let i = 0; i < 4096; i++) dec[i] = full[i * 4];
  const p = D.yinPitch(dec, 12000, 27, 1400);
  const fr = D.refinePitch(full, 48000, p.freq);
  const tg = D.tunerTarget(fr, D.TUNINGS.chromatic, 440, null);
  check('chromatic ' + D.noteName(m) + D.noteOctave(m) + ' +4c refined', tg.midi === m && Math.abs(tg.cents - 4) < 0.5, tg.cents.toFixed(2));
}
{
  // Pinned: low E tuned a full 3 semitones flat should still target E, not C#.
  const f = D.midiToFreq(37, 440);
  const tg = D.tunerTarget(f, std, 440, 0);
  check('pinned string reports -300c', Math.abs(tg.cents + 300) < 0.01 && tg.idx === 0);
}
check('silence -> no pitch', D.yinPitch(new Float32Array(4096), sr, 60, 500) === null);

/* ---- key ---- */
console.log('key');
function chordSeq(midis, secs) {
  const fs = 12000, n = fs * secs, ct = new D.ChromaTracker(fs, 4096);
  const buf = new Float32Array(n);
  for (let i = 0; i < n; i++) {
    const ch = midis[Math.floor(i / fs) % midis.length];
    let x = 0;
    for (const m of ch) { const f = D.midiToFreq(m, 440); x += 0.1 * Math.sin(2 * Math.PI * f * i / fs) + 0.04 * Math.sin(4 * Math.PI * f * i / fs); }
    buf[i] = x;
  }
  for (let p = 0; p + 4096 <= n; p += 2048) ct.frame(buf.subarray(p, p + 4096));
  return D.keyFromChroma(ct.acc);
}
{
  const k = chordSeq([[60, 64, 67], [65, 69, 72], [67, 71, 74], [60, 64, 67]], 16);   // I IV V I in C
  check('C major progression', k && k.camelot === '8B', k && k.name + ' ' + k.camelot);
  const k2 = chordSeq([[57, 60, 64], [62, 65, 69], [64, 68, 71], [57, 60, 64]], 16);  // i iv V i in A minor
  check('A minor progression', k2 && k2.camelot === '8A', k2 && k2.name + ' ' + k2.camelot);
}
check('Camelot neighbours of 12B', D.compatibleKeys('12B').join(',') === '11B,12A,1B', D.compatibleKeys('12B').join(','));
check('Camelot neighbours of 1A', D.compatibleKeys('1A').join(',') === '12A,1B,2A', D.compatibleKeys('1A').join(','));

console.log('\n' + passes + ' passed, ' + fails + ' failed');
process.exit(fails ? 1 : 0);
