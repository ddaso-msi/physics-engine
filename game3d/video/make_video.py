#!/usr/bin/env python3
"""Records a play-through of all five 3D levels and gives it a soundtrack.

    python3 game3d/video/make_video.py [output.mp4]

Needs a release build (build/release/game3d/puzzle_lab_3d), ffmpeg and numpy. The game is driven through its
own scripted mode with a dummy video driver, so no window opens: it plays for real, writes raw frames to a
pipe that ffmpeg encodes, and logs every audible event (impacts read from the contact solver, clicks,
results). soundtrack.py turns that log into audio and ffmpeg joins the two.

The cursor is steered onto points in the world (--at / --glide-at), so the script does not depend on where
each level's camera puts things on screen.
"""
import math
import os
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
GAME = os.path.join(ROOT, "build", "release", "game3d", "puzzle_lab_3d")
W, H = 1280, 720
args = []


def wait(seconds):
    args.extend(["--steps", str(round(seconds * 120))])


def to(p, seconds=0.5):
    args.extend(["--glide-at", *(f"{v:.3f}" for v in p), str(max(1, round(seconds * 120)))])


def click_at(p, travel=0.5):
    to(p, travel)
    wait(0.12)
    args.extend(["--down", "--up"])


def drag(a, b, travel=0.45, pull=0.6, keys=(), settle=0.15):
    to(a, travel)
    wait(0.1)
    args.append("--down")
    to(b, pull)
    for k in keys:
        args.extend(["--key", k])
        wait(0.12)
    wait(settle)
    args.append("--up")
    wait(0.1)


def key(name):
    args.extend(["--key", name])


def sweep(dyaw, dpitch, seconds):
    args.extend(["--sweep", str(dyaw), str(dpitch), str(round(seconds * 120))])


def run_button(travel=0.5):  # the Space Run / Stop button in the top bar
    args.extend(["--glide", "540", "28", str(round(travel * 120))])
    wait(0.12)
    args.extend(["--down", "--up"])


# ---- the play-through ----
args.extend(["--mouse", "1000", "600"])
wait(1.4)
args.extend(["--glide", "640", "282", "84"])  # level 1 on the menu
wait(0.15)
args.extend(["--down", "--up"])
wait(0.6)

# 1. Knock It Down: two balls. These two screen positions (through this level's starting camera) were found
# by playing every candidate headless, firing the first 1.9 s in and the second 1.1 s after it, as here.
args.extend(["--glide", "797", "136", "120"])
wait(0.3)
args.extend(["--down", "--up"])
args.extend(["--glide", "773", "192", "60"])
wait(0.6)
args.extend(["--down", "--up"])
sweep(-0.35, 0.05, 2.5)
wait(2.5)
key("next")

# 2. Bridge Builder: a bare deck drops the car; the same deck braced into triangles carries it.
def node(c, r):
    return (-4.0 + 2.0 * c, 3.9 if r == 1 else 1.9, 0.0)


def beam(a, b):
    drag(node(*a), node(*b), travel=0.28, pull=0.32, settle=0.05)


wait(0.8)
for c in range(4):
    beam((c, 1), (c + 1, 1))
run_button()
wait(2.6)
run_button(0.2)
wait(0.4)
for a, b in [((0, 0), (1, 1)), ((4, 0), (4, 1)), ((4, 0), (3, 1)), ((0, 0), (1, 0)), ((1, 0), (1, 1)), ((1, 0), (2, 1))]:
    beam(a, b)
run_button()
sweep(0.3, 0.0, 3.0)
wait(4.2)
key("next")

# 3. Pendulum Smash: pull the weight back 60 degrees, straight back from the bin.
wait(0.7)
pivot, rod = (0.0, 7.5, 0.0), 5.0
def bob(tilt, heading=0.0):
    out = (-math.cos(heading), 0.0, math.sin(heading))
    return (pivot[0] + out[0] * rod * math.sin(tilt), pivot[1] - rod * math.cos(tilt), pivot[2] + out[2] * rod * math.sin(tilt))
drag(bob(0.35), bob(math.radians(60)), travel=0.7, pull=1.4)
run_button(0.6)
wait(6.6)
key("next")

# 4. Chain Reaction: three spares into the gap, each turned to follow the curve.
wait(0.7)
def slot(k):
    a = k * 0.65 / 3.0 + 0.05
    return (3.0 * math.cos(a), 0.5, 3.0 * math.sin(a))
for i, (k, turns) in enumerate([(2, 4), (3, 3), (4, 3)]):
    drag((1.0 + 0.8 * i, 0.5, 1.2), slot(k), pull=0.7, keys=["left"] * turns)
run_button()
sweep(0.45, -0.1, 4.0)
wait(4.6)
key("next")

# 5. Impossible Tower: two blocks, and a plank stood on end on top.
wait(0.7)
drag((-3.4, 0.6, -1.6), (0, 0, 0), pull=0.8)
drag((-3.4, 0.6, 0.0), (0, 1.2, 0), pull=0.8)
drag((-6.0, 0.2, -2.0), (0, 2.4, 0), pull=0.9, keys=["tip"], settle=0.4)
run_button()
sweep(1.2, 0.05, 6.0)
wait(2.0)
key("escape")
wait(1.8)

# ---- record, score, join ----
out = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "puzzle-lab-3d-video.mp4"))
with tempfile.TemporaryDirectory() as tmp:
    fifo, silent, wav = os.path.join(tmp, "frames"), os.path.join(tmp, "silent.mp4"), os.path.join(tmp, "sound.wav")
    os.mkfifo(fifo)
    encoder = subprocess.Popen(["ffmpeg", "-loglevel", "error", "-y", "-f", "rawvideo", "-pixel_format", "rgb24", "-video_size", f"{W}x{H}",
                                "-framerate", "60", "-i", fifo, "-c:v", "libx264", "-preset", "slow", "-crf", "17", "-pix_fmt", "yuv420p", silent])
    game = subprocess.run([GAME, "--record", fifo] + args, env=dict(os.environ, SDL_VIDEODRIVER="dummy"))
    if game.returncode != 0 or encoder.wait() != 0:
        sys.exit("recording failed")
    frames = os.path.getsize(silent)  # (only to check something was written)
    if frames == 0:
        sys.exit("no frames")
    seconds = float(subprocess.run(["ffprobe", "-v", "error", "-show_entries", "format=duration", "-of", "csv=p=0", silent],
                                   capture_output=True, text=True, check=True).stdout)
    subprocess.run([sys.executable, os.path.join(os.path.dirname(__file__), "soundtrack.py"), fifo + ".events", wav, f"{seconds:.3f}"], check=True)
    subprocess.run(["ffmpeg", "-loglevel", "error", "-y", "-i", silent, "-i", wav, "-c:v", "copy", "-c:a", "aac", "-b:a", "192k", "-shortest",
                    "-movflags", "+faststart", out], check=True)
    os.replace(fifo + ".events", os.path.splitext(out)[0] + ".events.txt")
print(out)
