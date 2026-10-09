#!/usr/bin/env python3
"""Turns the game's event log into a soundtrack.

    python3 soundtrack.py events.txt out.wav seconds

Each log line is "time kind strength". Impacts come straight from the contact solver (the impulse of
a pair of bodies that has just started touching), so what you hear is what the simulation did: a
harder hit is louder and lower. Everything is synthesised here; there are no sample files.
"""
import sys
import wave

import numpy as np

RATE = 44100
rng = np.random.default_rng(7)


def envelope(n, attack, decay):
    t = np.arange(n) / RATE
    return np.minimum(t / attack, 1.0) * np.exp(-t / decay)


def tone(freq, seconds, decay, attack=0.004, harmonics=(1.0,)):
    n = int(seconds * RATE)
    t = np.arange(n) / RATE
    wave_ = sum(a * np.sin(2 * np.pi * freq * (k + 1) * t) for k, a in enumerate(harmonics))
    return wave_ * envelope(n, attack, decay)


def knock(strength):
    """A wooden knock: a short pitched thump plus a click of filtered noise."""
    level = min(1.0, 0.16 + 0.3 * np.log1p(strength))
    freq = 250.0 / (1.0 + 0.25 * strength) + rng.uniform(-18, 18)
    n = int(0.16 * RATE)
    t = np.arange(n) / RATE
    body = np.sin(2 * np.pi * freq * t * (1 + 0.6 * np.exp(-t / 0.012))) * np.exp(-t / 0.045)
    noise = rng.standard_normal(n) * np.exp(-t / 0.006)
    noise = np.convolve(noise, np.ones(6) / 6, mode="same")
    return level * (0.8 * body + 0.5 * noise)


def whoosh():
    n = int(0.22 * RATE)
    t = np.arange(n) / RATE
    sweep = np.sin(2 * np.pi * (520 * t - 1050 * t * t))
    out = 0.55 * sweep * envelope(n, 0.003, 0.06)
    thump = 0.25 * knock(1.2)
    out[:len(thump)] += thump
    return out


def chord(freqs, seconds, decay, gap=0.0):
    out = np.zeros(int((seconds + gap * len(freqs)) * RATE))
    for i, f in enumerate(freqs):
        note = tone(f, seconds, decay, harmonics=(1.0, 0.35, 0.12))
        start = int(i * gap * RATE)
        out[start:start + len(note)] += note
    return out


SOUNDS = {
    "fire": lambda s: whoosh(),
    "pick": lambda s: 0.22 * tone(660, 0.06, 0.02),
    "place": lambda s: 0.35 * tone(440, 0.16, 0.03, harmonics=(1.0, 0.5)) + 0.15 * knock(0.5),
    "select": lambda s: 0.3 * chord([523.25, 783.99], 0.25, 0.09, gap=0.07),
    "run": lambda s: 0.32 * chord([392.0, 523.25, 659.25], 0.22, 0.08, gap=0.06),
    "stop": lambda s: 0.28 * chord([523.25, 392.0], 0.2, 0.07, gap=0.07),
    "win": lambda s: 0.5 * chord([523.25, 659.25, 783.99, 1046.5], 1.4, 0.5, gap=0.09),
    "fail": lambda s: 0.4 * chord([311.13, 233.08], 0.7, 0.28, gap=0.16),
}


def music(seconds):
    """A quiet bed under the effects: a slow four-chord pad with a plucked line on top."""
    n = int(seconds * RATE)
    out = np.zeros(n)
    bar = 2.4
    chords = [(220.0, 261.63, 329.63), (174.61, 220.0, 261.63), (130.81, 196.0, 261.63 * 1.2599), (196.0, 246.94, 293.66)]
    pattern = [0, 1, 2, 1, 2, 1, 0, 1]
    for b in range(int(seconds / bar) + 1):
        notes = chords[b % 4]
        start = int(b * bar * RATE)
        m = min(int(bar * 1.15 * RATE), n - start)
        if m <= 0:
            break
        t = np.arange(m) / RATE
        swell = np.sin(np.pi * np.minimum(t / (bar * 1.15), 1.0)) ** 1.5
        for f in notes:
            out[start:start + m] += 0.035 * swell * (np.sin(2 * np.pi * f * t) + 0.3 * np.sin(2 * np.pi * 2 * f * t + 0.5))
        out[start:start + m] += 0.05 * swell * np.sin(2 * np.pi * notes[0] / 2 * t)
        for i, which in enumerate(pattern):
            pluck = 0.05 * tone(notes[which] * 2, 0.5, 0.16, harmonics=(1.0, 0.4, 0.15))
            at = start + int(i * bar / len(pattern) * RATE)
            k = min(len(pluck), n - at)
            if k > 0:
                out[at:at + k] += pluck[:k]
    fade = int(1.5 * RATE)
    out[:fade] *= np.linspace(0, 1, fade)
    out[-fade:] *= np.linspace(1, 0, fade)
    return out


def main():
    log, wav_path, seconds = sys.argv[1], sys.argv[2], float(sys.argv[3])
    n = int(seconds * RATE)
    left, right = np.zeros(n + RATE * 3), np.zeros(n + RATE * 3)
    last_knock, counts = -1.0, {}
    for line in open(log):
        when, kind, strength = line.split()
        when, strength = float(when), float(strength)
        if kind == "impact":
            # A pile settling makes dozens of tiny contacts a second: let the big ones through and thin the rest.
            if strength < 0.12 and when - last_knock < 0.03:
                continue
            last_knock = when
            sound = knock(strength)
        elif kind in SOUNDS:
            sound = SOUNDS[kind](strength)
        else:
            continue
        counts[kind] = counts.get(kind, 0) + 1
        pan = rng.uniform(0.35, 0.65) if kind == "impact" else 0.5
        at = int(when * RATE)
        left[at:at + len(sound)] += sound * (1 - pan) * 2
        right[at:at + len(sound)] += sound * pan * 2
    bed = music(seconds)
    left[:n] += bed
    right[:n] += bed
    mix = np.stack([left[:n], right[:n]], axis=1)
    mix = np.tanh(mix * 1.1) * 0.89  # soft limit rather than clip
    with wave.open(wav_path, "wb") as f:
        f.setnchannels(2)
        f.setsampwidth(2)
        f.setframerate(RATE)
        f.writeframes((mix * 32767).astype("<i2").tobytes())
    print("sounds:", ", ".join(f"{k} {v}" for k, v in sorted(counts.items())))


if __name__ == "__main__":
    main()
