# Not Quite My Tempo

Live tempo, tuner, key and loudness in the browser: <https://andrew-crena.github.io/not-quite-my-tempo/>

- **Tempo**: BPM from the onset envelope (log spectral flux, adaptive threshold, autocorrelation comb with sub-lag refinement), lock hysteresis, a phase-locked beat indicator, tap tempo, half/double and reset.
- **Tuner**: YIN pitch on a 12 kHz stream that is continuously notched at 50/60 Hz mains hum (and 100/120 Hz) and band-limited to the tuning's range, then refined against a band-limited copy of the full-rate signal. Readings at the edge of the search range are rejected, and a reading is shown only when most recent readings agree. Guitar standard and alternate tunings, 4- and 5-string bass, ukulele, chromatic. Adjustable A4 (415–466 Hz), select a string to lock onto it, in-tune within ±3 cents (with hysteresis).
- **Key**: chroma folded from spectral peaks with tuning compensation, matched against Krumhansl–Kessler profiles, with Camelot codes and compatible keys.
- **Loudness**: ITU-R BS.1770 K-weighted momentary, short-term and gated integrated LUFS.
- **Signal check**: peak, clips per minute, DC offset, noise floor, SNR, mains hum.

Everything runs locally; no audio leaves the device. No dependencies or build step.

## Files

- `index.html`: the interface and audio plumbing.
- `dsp.js`: all signal processing, free of DOM code so it can be tested.
- `tests/dsp.test.js`: synthetic-signal checks (`node tests/dsp.test.js`).

## Accessibility

Keyboard shortcuts (Space, 1, 2, K, T, R, H, D), screen-reader announcements on lock, key and tuning changes (not every frame), text equivalents for every colour cue, honours reduced motion and forced colours, and 4.5:1 text contrast.
