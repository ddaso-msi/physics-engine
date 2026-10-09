#!/usr/bin/env python3
"""Records a play-through of all five levels and gives it a soundtrack.

    python3 game/video/make_video.py [output.mp4]

Needs a release build (build/release/game/puzzle_lab), ffmpeg and numpy. The game is driven through its
own scripted mode with a dummy video driver, so no window opens: it plays for real, writes raw frames
to a pipe that ffmpeg encodes, and logs every audible event (impacts read from the contact solver,
clicks, results). soundtrack.py then turns that log into audio and ffmpeg joins the two.
"""
import math
import os
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
GAME = os.path.join(ROOT, "build", "release", "game", "puzzle_lab")
W, H, PPM = 1280, 720, 64.0

args = []


def world(px, py):  # screen pixels -> world metres
    return px / PPM, H / PPM - py / PPM


def wait(seconds):
    args.extend(["--steps", str(round(seconds * 120))])


def glide(x, y, seconds=0.5):
    args.extend(["--glide", f"{x:.3f}", f"{y:.3f}", str(round(seconds * 120))])


def click(x, y, travel=0.5):
    glide(x, y, travel)
    wait(0.12)
    args.extend(["--click", f"{x:.3f}", f"{y:.3f}"])


def drag(x0, y0, x1, y1, travel=0.45, pull=0.6, before_release=()):
    glide(x0, y0, travel)
    args.extend(["--press", f"{x0:.3f}", f"{y0:.3f}"])
    glide(x1, y1, pull)
    args.extend(before_release)
    wait(0.1)
    args.extend(["--release", f"{x1:.3f}", f"{y1:.3f}"])
    wait(0.1)


def key(name):
    args.extend(["--key", name])


RUN_BUTTON = world(540, 28)


def node(column, row):  # bridge connection points
    return 6.0 + 2.0 * column, 3.9 + 2.0 * (row - 1)


def beam(a, b):
    drag(*node(*a), *node(*b), travel=0.3, pull=0.35)


# ---- the play-through ----
args.extend(["--mouse", "16", "3"])
wait(1.3)
click(*world(640, 282), travel=0.7)  # level 1 on the menu
wait(0.6)

# 1. Knock It Down: three balls along one arc.
glide(14.5, 6.0, 1.0)
wait(0.4)
for _ in range(3):
    args.extend(["--click", "14.5", "6.0"])
    wait(0.5)
wait(4.6)
key("next")

# 2. Bridge Builder: a bare roadway fails, then the same roadway braced into triangles holds.
wait(0.8)
for c in range(4):
    beam((c, 1), (c + 1, 1))
click(*RUN_BUTTON)
wait(2.9)
click(*RUN_BUTTON, travel=0.2)
wait(0.4)
for a, b in [((0, 0), (1, 1)), ((4, 0), (4, 1)), ((4, 0), (3, 1)), ((0, 0), (1, 0)), ((1, 0), (1, 1)), ((1, 0), (2, 1))]:
    beam(a, b)
click(*RUN_BUTTON)
wait(6.6)
key("next")

# 3. Pendulum Smash: pull the weight back to 62 degrees.
wait(0.7)
pivot, rod = (6.6, 9.6), 5.2
start, angle = -0.35, math.radians(-62)
drag(pivot[0] + rod * math.sin(start), pivot[1] - rod * math.cos(start),
     pivot[0] + rod * math.sin(angle), pivot[1] - rod * math.cos(angle), travel=0.7, pull=1.3)
click(*RUN_BUTTON, travel=0.6)
wait(6.6)
key("next")

# 4. Chain Reaction: three spare dominoes into the gap.
wait(0.7)
for home, x in [(7.3, 7.25), (8.0, 8.0), (8.7, 8.75)]:
    drag(home, 6.6, x, 4.0)
click(*RUN_BUTTON)
wait(8.2)
key("next")

# 5. Impossible Tower: two blocks and a plank stood on end.
wait(0.7)
drag(2.1, 6.5, 10.0, 2.6, pull=0.8)
drag(3.7, 6.5, 10.0, 4.6, pull=0.8)
drag(2.9, 8.9, 10.0, 7.2, pull=0.9, before_release=["--key", "left", "--steps", "40"])
click(*RUN_BUTTON)
wait(7.6)
key("escape")
wait(1.8)

# ---- record, score, join ----
out = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "puzzle-lab-video.mp4"))
with tempfile.TemporaryDirectory() as tmp:
    fifo, silent, wav = os.path.join(tmp, "frames"), os.path.join(tmp, "silent.mp4"), os.path.join(tmp, "sound.wav")
    os.mkfifo(fifo)
    encoder = subprocess.Popen(["ffmpeg", "-loglevel", "error", "-y", "-f", "rawvideo", "-pixel_format", "rgb24", "-video_size", f"{W}x{H}",
                                "-framerate", "60", "-i", fifo, "-c:v", "libx264", "-preset", "slow", "-crf", "17", "-pix_fmt", "yuv420p", silent])
    game = subprocess.run([GAME, "--record", fifo] + args, env=dict(os.environ, SDL_VIDEODRIVER="dummy"))
    if game.returncode != 0 or encoder.wait() != 0:
        sys.exit("recording failed")
    steps = sum(int(args[i + 1]) for i, a in enumerate(args) if a == "--steps") + sum(int(args[i + 3]) for i, a in enumerate(args) if a == "--glide")
    subprocess.run([sys.executable, os.path.join(os.path.dirname(__file__), "soundtrack.py"), fifo + ".events", wav, f"{steps / 120:.3f}"], check=True)
    subprocess.run(["ffmpeg", "-loglevel", "error", "-y", "-i", silent, "-i", wav, "-c:v", "copy", "-c:a", "aac", "-b:a", "192k", "-shortest",
                    "-movflags", "+faststart", out], check=True)
    os.replace(fifo + ".events", os.path.splitext(out)[0] + ".events.txt")
print(out)
